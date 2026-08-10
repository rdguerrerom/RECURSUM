#!/usr/bin/env python3
"""Demonstration of the generic value-DAG front-end (recursum_value_dag) on a
determinant-CI p^dagger q (E_pq) contraction body.

Target: the `ci_rdm12` inner reduction from SHRIKE's include/shrike/ci.hpp,
    rdm1[pq]    = sum_K B[pq,K] cbra[K]
    rdm2[pq,rs] = sum_K A[qp,K] B[rs,K] - delta_qr rdm1[ps]
where A = E_ab cbra and B = E_ab cket are the resolved E-blocks (nact^2, ndet).

For a FIXED small active space (CAS(2,2), CAS(4,4)) we express this block as a
value-DAG and emit CSE'd C via emit_value_dag, then:
  1. compile & run the emitted C on random seeded active-space vectors;
  2. validate rdm1/rdm2 against PySCF make_rdm12 / trans_rdm12 to <= 1e-12
     (PySCF is the reference SHRIKE's ci.hpp is itself gated against, bit-for-bit,
     in python/test_shrike_ci.py);
  3. report the CSE statistics and a CSE-straight-line-C vs cblas_dgemm timing.

The A/B resolved blocks are produced by a faithful numpy port of ci.hpp's
ci_strings + ci_resolve_e (verified because the resulting rdm1/rdm2 match PySCF).
"""
import os, sys, subprocess, struct, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from recursum_value_dag import ValueDAG, emit_value_dag, cse_stats

OUT = os.path.join(HERE, "ci_valuedag_out")
os.makedirs(OUT, exist_ok=True)


# ── faithful numpy port of ci.hpp ci_strings + ci_resolve_e ──────────────────
def ci_strings(nact, nelec):
    """Return (occ, exc) where exc[J] is a list of (p,q,sign,addr) single
    excitations a^+_p a_q |J> = sign |addr> (incl. diagonal p==q). Colex rank."""
    from math import comb
    nstr = comb(nact, nelec) if 0 <= nelec <= nact else 0
    if nstr <= 0:
        nstr = 1
    # enumerate k-subsets in colex order -> occ[addr]
    def rank(occ):
        return sum(comb(occ[j], j + 1) for j in range(len(occ)))
    occ = [None] * nstr
    comb0 = list(range(nelec))
    while True:
        occ[rank(comb0)] = tuple(comb0)
        if nelec == 0:
            break
        i = nelec - 1
        while i >= 0 and comb0[i] == nact - nelec + i:
            i -= 1
        if i < 0:
            break
        comb0[i] += 1
        for j in range(i + 1, nelec):
            comb0[j] = comb0[j - 1] + 1
    if nelec == 0:
        occ[0] = ()
    exc = [[] for _ in range(nstr)]
    for J in range(nstr):
        o = occ[J]
        for iq in range(nelec):
            q = o[iq]
            for p in range(nact):
                if p == q:
                    exc[J].append((p, q, 1, J))
                    continue
                if p in o:
                    continue
                rest = [o[e] for e in range(nelec) if e != iq]        # k-1 entries
                ip = 0
                while ip < len(rest) and rest[ip] < p:
                    ip += 1
                sign = -1 if ((iq + ip) & 1) else 1
                newocc = rest[:ip] + [p] + rest[ip:]
                exc[J].append((p, q, sign, rank(newocc)))
    return occ, exc


def resolve_e(nact, sa_exc, sb_exc, cin):
    """out[p*nact+q, Ia, Ib] = (E_pq cin) — port of ci_resolve_e. cin:(ndeta,ndetb)."""
    ndeta, ndetb = cin.shape
    out = np.zeros((nact * nact, ndeta, ndetb))
    for Ja in range(ndeta):
        for (p, q, sg, addr) in sa_exc[Ja]:
            out[p * nact + q, addr, :] += sg * cin[Ja, :]
    for Jb in range(ndetb):
        for (p, q, sg, addr) in sb_exc[Jb]:
            out[p * nact + q, :, addr] += sg * cin[:, Jb]
    return out.reshape(nact * nact, ndeta * ndetb)


def build_blocks(nact, na, nb, cbra, cket):
    _, ea = ci_strings(nact, na)
    _, eb = ci_strings(nact, nb)
    A = resolve_e(nact, ea, eb, cbra)          # E cbra  (nact^2, ndet)
    B = resolve_e(nact, ea, eb, cket)          # E cket
    return A, B


# ── value-DAG for the rdm1/rdm2 reduction ────────────────────────────────────
def build_rdm_dag(nact, ndet):
    """rdm1[pq]=sum_K B[pq,K] cbra[K];  rdm2[pq,rs]=sum_K A[qp,K] B[rs,K] - d_qr rdm1[ps].
    Inputs A,B row-major (nact^2, ndet) flattened; cbra (ndet). Outputs rdm1
    (nact^2), rdm2 (nact^4)."""
    g = ValueDAG()
    n2 = nact * nact
    A = [[g.input("A", pq * ndet + K) for K in range(ndet)] for pq in range(n2)]
    B = [[g.input("B", pq * ndet + K) for K in range(ndet)] for pq in range(n2)]
    cb = [g.input("cbra", K) for K in range(ndet)]

    rdm1 = [g.dot(B[pq], cb) for pq in range(n2)]          # nact^2 nodes
    for pq in range(n2):
        g.output("rdm1", pq, rdm1[pq])

    for p in range(nact):
        for q in range(nact):
            pq = p * nact + q
            qp = q * nact + p
            for r in range(nact):
                for s in range(nact):
                    rs = r * nact + s
                    val = g.dot(A[qp], B[rs])              # sum_K A_qp[K] B_rs[K]
                    if q == r:                             # -delta_qr rdm1[ps]
                        val = g.sub(val, rdm1[p * nact + s])
                    g.output("rdm2", pq * n2 + rs, val)
    return g


# ── C driver + compile + run ─────────────────────────────────────────────────
DRIVER = r"""
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cblas.h>
#define NACT %(nact)d
#define NDET %(ndet)d
#define N2  (NACT*NACT)
#define N4  (N2*N2)
%(kernel)s
int main(void){
    static double A[N2*NDET], B[N2*NDET], cbra[NDET];
    static double rdm1[N2], rdm2[N4];
    if(fread(A,sizeof(double),N2*NDET,stdin)!=(size_t)N2*NDET) return 2;
    if(fread(B,sizeof(double),N2*NDET,stdin)!=(size_t)N2*NDET) return 2;
    if(fread(cbra,sizeof(double),NDET,stdin)!=(size_t)NDET) return 2;
    %(fname)s(A,B,cbra,rdm1,rdm2);
    fwrite(rdm1,sizeof(double),N2,stdout);
    fwrite(rdm2,sizeof(double),N4,stdout);

    /* --- timing: CSE straight-line kernel vs cblas_dgemm for the rdm2 block --- */
    const int REP = %(rep)d;
    struct timespec t0,t1;
    /* Ap = A with (q,p)->(p,q) pair transpose (the GEMM left operand, as ci.hpp) */
    static double Ap[N2*NDET], rd2[N4];
    for(int p=0;p<NACT;p++)for(int q=0;q<NACT;q++)
        for(int K=0;K<NDET;K++) Ap[(p*NACT+q)*NDET+K]=A[(q*NACT+p)*NDET+K];
    clock_gettime(CLOCK_MONOTONIC,&t0);
    for(int it=0;it<REP;it++) %(fname)s(A,B,cbra,rdm1,rdm2);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    double cse_ns=((t1.tv_sec-t0.tv_sec)*1e9+(t1.tv_nsec-t0.tv_nsec))/REP;
    clock_gettime(CLOCK_MONOTONIC,&t0);
    for(int it=0;it<REP;it++)
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,N2,N2,NDET,
                    1.0,Ap,NDET,B,NDET,0.0,rd2,N2);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    double gemm_ns=((t1.tv_sec-t0.tv_sec)*1e9+(t1.tv_nsec-t0.tv_nsec))/REP;
    fprintf(stderr,"TIMING nact=%%d ndet=%%d  CSE=%%.1f ns  GEMM=%%.1f ns  speedup=%%.2fx\n",
            NACT,NDET,cse_ns,gemm_ns,gemm_ns/cse_ns);
    return 0;
}
"""


def run_case(nact, na, nb, seed_bra=1, seed_ket=2, state=True):
    from math import comb
    ndeta, ndetb = comb(nact, na), comb(nact, nb)
    ndet = ndeta * ndetb
    rb = np.random.default_rng(seed_bra).standard_normal((ndeta, ndetb)); rb /= np.linalg.norm(rb)
    if state:
        cket = cbra = rb
    else:
        rk = np.random.default_rng(seed_ket).standard_normal((ndeta, ndetb)); rk /= np.linalg.norm(rk)
        cbra, cket = rb, rk

    A, B = build_blocks(nact, na, nb, cbra, cket)
    cbra_flat = cbra.reshape(-1)

    g = build_rdm_dag(nact, ndet)
    fname = f"ci_rdm_cas{na+nb}_{nact}"
    src = emit_value_dag(g, inputs=["A", "B", "cbra"], outputs=["rdm1", "rdm2"],
                         fname=fname, style="ssa")
    st = cse_stats(g)

    base = os.path.join(OUT, fname)
    with open(base + ".c", "w") as f:
        f.write(DRIVER % dict(nact=nact, ndet=ndet, kernel=src, fname=fname, rep=200000))
    subprocess.run(["gcc", "-O2", "-march=native", base + ".c", "-o", base,
                    "-lopenblas"], check=True)

    inp = A.reshape(-1).tobytes() + B.reshape(-1).tobytes() + cbra_flat.tobytes()
    p = subprocess.run([base], input=inp, capture_output=True, check=True)
    timing = p.stderr.decode().strip()
    buf = p.stdout
    n2, n4 = nact * nact, nact ** 4
    rdm1_c = np.frombuffer(buf[:n2 * 8], dtype=np.float64).reshape(nact, nact)
    rdm2_c = np.frombuffer(buf[n2 * 8:n2 * 8 + n4 * 8], dtype=np.float64).reshape(nact, nact, nact, nact)

    # numpy reference from the SAME resolved blocks (validates the DAG algebra)
    Ap = A.reshape(nact, nact, ndet).transpose(1, 0, 2).reshape(n2, ndet)   # (q,p)->(p,q)
    rdm1_np = (B @ cbra_flat).reshape(nact, nact)
    rdm2_np = (Ap @ B.T).reshape(nact, nact, nact, nact)
    for p_ in range(nact):
        for q_ in range(nact):
            for s_ in range(nact):
                rdm2_np[p_, q_, q_, s_] -= rdm1_np[p_, s_]

    # PySCF reference (the accepted oracle)
    from pyscf import fci
    if state:
        p1, p2 = fci.direct_spin1.make_rdm12(cbra, nact, (na, nb))
    else:
        p1, p2 = fci.direct_spin1.trans_rdm12(cbra, cket, nact, (na, nb))
        p1 = p1.T
    err_np1 = np.max(np.abs(rdm1_c - rdm1_np))
    err_np2 = np.max(np.abs(rdm2_c - rdm2_np))
    err_ps1 = np.max(np.abs(rdm1_c - p1))
    err_ps2 = np.max(np.abs(rdm2_c - p2))

    tag = "state" if state else "trans"
    print(f"\n=== CAS({na+nb},{nact}) [{tag}]  nact={nact} ndet={ndet} ===")
    print(f"  CSE stats: {st['nodes']} nodes, {st['ops_cse']} ops (CSE) vs "
          f"{st['ops_naive']} naive  ->  x{st['cse_factor']:.2f} fewer ops; "
          f"{st['n_outputs']} outputs")
    print(f"  emitted C: {base}.c  ({len(src.splitlines())} lines)")
    print(f"  rdm1 vs numpy = {err_np1:.2e}   rdm2 vs numpy = {err_np2:.2e}")
    print(f"  rdm1 vs PySCF = {err_ps1:.2e}   rdm2 vs PySCF = {err_ps2:.2e}"
          f"   -> {'PASS' if max(err_ps1,err_ps2)<=1e-12 else 'FAIL'} (tol 1e-12)")
    print(f"  {timing}")
    return max(err_ps1, err_ps2) <= 1e-12


if __name__ == "__main__":
    ok = True
    for nact, na, nb in [(2, 1, 1), (4, 2, 2)]:
        ok &= run_case(nact, na, nb, state=True)
        ok &= run_case(nact, na, nb, state=False)
    print("\n" + "=" * 62)
    print("value-DAG CI demo:", "ALL PASS (<=1e-12 vs PySCF)" if ok else "FAIL")
    sys.exit(0 if ok else 1)
