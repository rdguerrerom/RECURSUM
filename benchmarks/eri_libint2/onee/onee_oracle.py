"""Validate the ported one-electron OS DAGs (overlap S, kinetic T, nuclear V)
against PySCF int1e_*_cart, by building the full contracted-Cartesian matrices
*from the DAG* (onee_dag) rather than from a separate recursion.

For each shell pair and primitive pair we compute the OS scalars, evaluate the
DAG numerically (onee_dag.evaluate), contract, and apply the libcint per-shell
Cartesian normalization. This proves the ported recurrence DAGs reproduce PySCF.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import gamma as _gamma, gammainc as _gammainc_reg

from onee_dag import OneE, build_dag, evaluate, kinetic_terms
from onee_dag import all_ang_tuples

_XYZ = "xyz"


# ---- Boys F_0..F_nmax(T) (reference accuracy) -----------------------------
def boys_array(nmax, T):
    F = np.empty(nmax + 1)
    if T < 1e-13:
        for n in range(nmax + 1):
            F[n] = 1.0 / (2 * n + 1)
        return F
    for n in range(nmax + 1):
        a = n + 0.5
        F[n] = _gamma(a) * _gammainc_reg(a, T) / (2.0 * T ** a)
    return F


# ---- libcint per-shell Cartesian normalization ----------------------------
def _dfact(n):
    r = 1
    for i in range(1, n + 1):
        r *= 2 * i - 1
    return r if n > 0 else 1


def norm_l(l, a):
    return (2.0 * a / math.pi) ** 0.75 * (4.0 * a) ** (l / 2.0) / math.sqrt(_dfact(l))


def gto_norm(l, a):
    return math.sqrt(2.0 ** (2 * l + 3) * math.factorial(l + 1) * (2.0 * a) ** (l + 1.5)
                     / (math.factorial(2 * l + 2) * math.sqrt(math.pi)))


def cart_norm(l, a):
    return norm_l(l, a) if l <= 1 else gto_norm(l, a)


# ---- shells from a PySCF mol ----------------------------------------------
def get_shells(mol):
    sh = []
    off = 0
    for s in range(mol.nbas):
        l = int(mol.bas_angular(s))
        c = np.asarray(mol.bas_coord(s))
        exps = np.asarray(mol.bas_exp(s))
        ctr = np.asarray(mol.bas_ctr_coeff(s))  # (nprim, nctr)
        for col in range(ctr.shape[1]):
            ncart = (l + 1) * (l + 2) // 2
            sh.append(dict(l=l, c=c, e=exps, d=ctr[:, col].copy(), off=off, ncart=ncart))
            off += ncart
    return sh, off


def _scalars_AB(al, A, be, B):
    zeta = al + be
    P = (al * A + be * B) / zeta
    AB = A - B
    K = math.exp(-al * be / zeta * float(np.dot(AB, AB)))
    PA = P - A
    inv2z = 1.0 / (2.0 * zeta)
    base = {f"PA{_XYZ[i]}": PA[i] for i in range(3)}
    base.update({f"AB{_XYZ[i]}": AB[i] for i in range(3)})
    base["inv_2zeta"] = inv2z
    base["one"] = 1.0
    return P, zeta, K, base


def overlap_block(shA, shB):
    """Contracted Cartesian overlap block (ncA, ncB) via the DAG."""
    la, lb = shA["l"], shB["l"]
    A, B = shA["c"], shB["c"]
    ca, cb = all_ang_tuples(la), all_ang_tuples(lb)
    outs = [OneE(a=a, b=b, m=0) for a in ca for b in cb]
    dag = build_dag(outs, op="overlap")
    blk = np.zeros((len(ca), len(cb)))
    for ip, al in enumerate(shA["e"]):
        wa = shA["d"][ip] * cart_norm(la, al)
        for jp, be in enumerate(shB["e"]):
            wb = shB["d"][jp] * cart_norm(lb, be)
            P, zeta, K, sc = _scalars_AB(al, A, be, B)
            S00 = K * (math.pi / zeta) ** 1.5
            vals = evaluate(dag, sc, base_eval=lambda n: S00)
            w = wa * wb
            for ia, a in enumerate(ca):
                for ib, b in enumerate(cb):
                    blk[ia, ib] += w * vals[OneE(a=a, b=b, m=0)]
    return blk


def nuclear_block(shA, shB, atoms):
    """Contracted Cartesian nuclear-attraction block via the DAG (sum over nuclei)."""
    la, lb = shA["l"], shB["l"]
    A, B = shA["c"], shB["c"]
    ca, cb = all_ang_tuples(la), all_ang_tuples(lb)
    outs = [OneE(a=a, b=b, m=0) for a in ca for b in cb]
    dag = build_dag(outs, op="nuclear")
    nmax = la + lb
    blk = np.zeros((len(ca), len(cb)))
    for ip, al in enumerate(shA["e"]):
        wa = shA["d"][ip] * cart_norm(la, al)
        for jp, be in enumerate(shB["e"]):
            wb = shB["d"][jp] * cart_norm(lb, be)
            zeta = al + be
            P = (al * A + be * B) / zeta
            AB = A - B
            K = math.exp(-al * be / zeta * float(np.dot(AB, AB)))
            PA = P - A
            inv2z = 1.0 / (2.0 * zeta)
            w = wa * wb
            for (Z, C) in atoms:
                PC = P - C
                T = zeta * float(np.dot(PC, PC))
                F = boys_array(nmax, T)
                base_pref = -2.0 * math.pi / zeta * K * Z
                kf = base_pref * F
                sc = {f"PA{_XYZ[i]}": PA[i] for i in range(3)}
                sc.update({f"PC{_XYZ[i]}": PC[i] for i in range(3)})
                sc.update({f"AB{_XYZ[i]}": AB[i] for i in range(3)})
                sc["inv_2zeta"] = inv2z
                sc["one"] = 1.0
                vals = evaluate(dag, sc, base_eval=lambda n: kf[n.m])
                for ia, a in enumerate(ca):
                    for ib, b in enumerate(cb):
                        blk[ia, ib] += w * vals[OneE(a=a, b=b, m=0)]
    return blk


def kinetic_block(shA, shB):
    """Contracted Cartesian kinetic block via the overlap DAG + KinTerm combo."""
    la, lb = shA["l"], shB["l"]
    A, B = shA["c"], shB["c"]
    ca, cb = all_ang_tuples(la), all_ang_tuples(lb)
    # overlap outputs needed: (a,b), (a,b+2_i), (a,b-2_i) for every output+axis
    need = set()
    for a in ca:
        for b in cb:
            for kt in kinetic_terms(a, b):
                need.add(kt.source)          # OneE overlap nodes
    dag = build_dag(list(need), op="overlap")
    blk = np.zeros((len(ca), len(cb)))
    for ip, al in enumerate(shA["e"]):
        wa = shA["d"][ip] * cart_norm(la, al)
        for jp, be in enumerate(shB["e"]):
            wb = shB["d"][jp] * cart_norm(lb, be)
            P, zeta, K, sc = _scalars_AB(al, A, be, B)
            S00 = K * (math.pi / zeta) ** 1.5
            S = evaluate(dag, sc, base_eval=lambda n: S00)
            w = wa * wb
            beta = be
            for ia, a in enumerate(ca):
                for ib, b in enumerate(cb):
                    t = 0.0
                    for kt in kinetic_terms(a, b):
                        sval = S[kt.source]
                        if kt.kind == "beta_2b1":
                            coef = beta * (2 * kt.bi + 1)
                        elif kt.kind == "beta2":
                            coef = 2.0 * beta * beta
                        else:  # half_bb1
                            coef = 0.5 * kt.bi * (kt.bi - 1)
                        t += kt.sign * coef * sval
                    blk[ia, ib] += w * t
    return blk


def build_matrix(shells, nao, op, atoms=None):
    M = np.zeros((nao, nao))
    fn = {"overlap": overlap_block, "kinetic": kinetic_block, "nuclear": nuclear_block}[op]
    for i, sa in enumerate(shells):
        for j, sb in enumerate(shells):
            if j > i:
                continue
            blk = fn(sa, sb, atoms) if op == "nuclear" else fn(sa, sb)
            oa, ob = sa["off"], sb["off"]
            M[oa:oa + sa["ncart"], ob:ob + sb["ncart"]] = blk
            if i != j:
                M[ob:ob + sb["ncart"], oa:oa + sa["ncart"]] = blk.T
    return M


if __name__ == "__main__":
    from pyscf import gto
    for atom, basis, tag in [
        ("O 0 0 0; H 0 0 1.81; H 1.72 0 -0.55", "6-31g", "H2O/6-31G (s,p)"),
        ("N 0 0 0; H 0 0 1.9", "cc-pvdz", "NH/cc-pVDZ (s,p,d)"),
        ("C 0 0 0; O 0 0 2.1", "cc-pvtz", "CO/cc-pVTZ (s,p,d,f)"),
    ]:
        mol = gto.M(atom=atom, basis=basis, unit="Bohr", cart=True, verbose=0)
        sh, nao = get_shells(mol)
        maxL = max(s["l"] for s in sh)
        atoms = [(float(mol.atom_charge(k)), np.asarray(mol.atom_coord(k)))
                 for k in range(mol.natm)]
        S = build_matrix(sh, nao, "overlap")
        T = build_matrix(sh, nao, "kinetic")
        V = build_matrix(sh, nao, "nuclear", atoms)
        Sref = mol.intor("int1e_ovlp_cart")
        Tref = mol.intor("int1e_kin_cart")
        Vref = mol.intor("int1e_nuc_cart")
        def relerr(X, R):
            return np.max(np.abs(X - R)) / max(np.max(np.abs(R)), 1e-30)
        eS, eT, eV = relerr(S, Sref), relerr(T, Tref), relerr(V, Vref)
        ok = max(eS, eT, eV) < 1e-11
        print(f"{tag:24s} nao={nao:3d} maxL={maxL}  "
              f"S={eS:.1e} T={eT:.1e} V={eV:.1e}  {'PASS' if ok else 'FAIL'}")
