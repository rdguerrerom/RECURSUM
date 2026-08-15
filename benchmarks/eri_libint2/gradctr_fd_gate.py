"""FD-oracle correctness gate for the gradctr (augmented value+gradient)
codegen restored in recursum_dag_emit.py on 2026-08-15 (see
SHRIKE_Research/plan.md L76 for the full audit trail of why this needed
restoring: it existed and was validated in commits
fd52486->3b2cf49->8cd8438 on this repo's history, then silently dropped by
c39d4c6's wholesale file sync).

Methodology (mirrors the ORIGINAL validation quoted in those commit
messages, since no script was ever committed for it -- see plan.md L76(e)):
for a well-separated, non-degenerate single primitive quartet, compare the
gradctr kernel's analytic reconstruction

    d/dX_c = up(a+1_c) - a_c * dn(a-1_c)        (up is ALREADY 2*alpha-weighted
                                                  by the kernel at contraction,
                                                  spec SS8.4 -- do not re-weight)

against CENTRAL FINITE DIFFERENCE of the value kernel (recursum_call,
eri_classes.h) -- the SAME already-validated value kernels the shift
identity is defined in terms of (RECURSUM_GRADIENT_RECOMMENDATIONS.md SS4:
"validate it against finite difference of *those* blocks").

Usage:
    python3 gradctr_fd_gate.py [--classes ssss,ssps,pppp,...] [--workdir DIR]

Exits 0 and prints ALL PASS iff every class's worst |analytic-FD| < TOL.
"""
from __future__ import annotations

import argparse
import copy
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recursum_dag_emit as R
from oracle import Quartet, compute_scalars, boys
from eri_dag.integral import output_set, Integral
from eri_dag.dag import build_dag

L = "spdfg"
TOL = 2e-9  # matches the ~3.2e-11 originally-reported precision with margin

SCALAR_ORDER = ["PAx", "PAy", "PAz", "QCx", "QCy", "QCz", "WPx", "WPy", "WPz",
                "WQx", "WQy", "WQz", "ABx", "ABy", "ABz", "CDx", "CDy", "CDz",
                "inv_2zp", "inv_2zq", "inv_2zpq", "frac_q_over_pq", "frac_p_over_pq"]


def cls_name(c):
    return "".join(L[x] for x in c)


def default_quartet():
    # Well-separated, non-degenerate, no special symmetry -- avoids masking
    # sign/index bugs the way a symmetric geometry could.
    return Quartet(a=1.2, A=(0.05, 0.02, 0.01),
                   b=0.8, B=(0.5, 0.1, -0.2),
                   g=1.5, C=(0.0, 0.7, 0.15),
                   d=0.9, D=(-0.1, 0.05, 0.9))


class Harness:
    def __init__(self, classes, workdir):
        self.workdir = workdir
        os.makedirs(workdir, exist_ok=True)
        gtus, skipped = R.generate_gradctr_project(classes, outdir=workdir)
        if skipped:
            raise SystemExit(f"gradctr_fd_gate: classes skipped (>CHUNK_THRESHOLD, "
                              f"chunked gradctr not restored): {skipped}")
        vtus = R.generate_project(classes, outdir=workdir)
        import shutil
        shutil.copy(os.path.join(os.path.dirname(__file__), "recursum_eri_scalars.h"), workdir)
        shutil.copy(os.path.join(os.path.dirname(__file__), "validate_kernels.cpp"), workdir)
        shutil.copy(os.path.join(os.path.dirname(__file__), "gradctr_fd_gate.cpp"), workdir)
        value_srcs = [f for f in vtus if "_naive" not in f]
        self._run(["c++", "-O2", "-std=c++17", "-I", workdir,
                    os.path.join(workdir, "validate_kernels.cpp"),
                    *[os.path.join(workdir, f) for f in value_srcs],
                    "-o", os.path.join(workdir, "value_harness")])
        self._run(["c++", "-O2", "-std=c++17", "-I", workdir,
                    os.path.join(workdir, "gradctr_fd_gate.cpp"),
                    *[os.path.join(workdir, f) for f in gtus],
                    "-o", os.path.join(workdir, "gradctr_harness")])

    def _run(self, cmd):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit(f"gradctr_fd_gate: compile failed: {' '.join(cmd)}\n{r.stderr}")

    def _exec(self, name, inp):
        r = subprocess.run([os.path.join(self.workdir, name)], input=inp,
                            capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"{name} failed rc={r.returncode}: {r.stderr}")
        return r.stdout.split()

    def value_call(self, cls, q):
        dag = build_dag(output_set(*cls))
        mmax = dag.max_m()
        s = compute_scalars(q)
        F = boys(mmax, s["_T"]); K = s["_K"]
        kf = [K * f for f in F]
        n_out = len(output_set(*cls))
        inp = f"{cls_name(cls)} {mmax} {n_out}\n"
        inp += " ".join(f"{s[k]:.17e}" for k in SCALAR_ORDER) + "\n"
        inp += " ".join(f"{v:.17e}" for v in kf) + "\n"
        return [float(x) for x in self._exec("value_harness", inp)]

    def gradctr_sections(self, cls, q):
        outs_list, sections = R.grad_output_layout(*cls)
        dag = build_dag(outs_list)
        mmax = dag.max_m()
        nout = len(outs_list)
        s = compute_scalars(q)
        F = boys(mmax, s["_T"]); K = s["_K"]
        kf = [K * f for f in F]
        w2 = [2 * q.a, 2 * q.b, 2 * q.g]
        inp = f"{cls_name(cls)} {nout} {mmax}\n"
        inp += " ".join(f"{s[k]:.17e}" for k in SCALAR_ORDER) + "\n"
        inp += " ".join(f"{v:.17e}" for v in kf) + "\n"
        inp += " ".join(f"{v:.17e}" for v in w2) + "\n"
        out = self._exec("gradctr_harness", inp)
        if out and out[0] == "NOFN":
            raise RuntimeError(f"no gradctr fn dispatched for {cls_name(cls)}")
        vals = [float(x) for x in out]
        return {k: vals[st:st + cnt] for k, (st, cnt, _c) in sections.items()}


def shift_q(q, centre, comp, h):
    qq = copy.deepcopy(q)
    pts = [list(qq.A), list(qq.B), list(qq.C), list(qq.D)]
    pts[centre][comp] += h
    qq.A, qq.B, qq.C, qq.D = tuple(pts[0]), tuple(pts[1]), tuple(pts[2]), tuple(pts[3])
    return qq


def check_class(hz, cls, q, h=1e-5):
    Ls = list(cls)
    secs = hz.gradctr_sections(cls, q)
    outs_val = output_set(*cls)
    worst, worst_desc = 0.0, None
    for ci, cn in enumerate(["A", "B", "C"]):
        up_list = output_set(*[Ls[k] + (1 if k == ci else 0) for k in range(4)])
        up_idx = {t: i for i, t in enumerate(up_list)}
        has_dn = Ls[ci] >= 1
        if has_dn:
            dn_list = output_set(*[Ls[k] - (1 if k == ci else 0) for k in range(4)])
            dn_idx = {t: i for i, t in enumerate(dn_list)}
        for comp in range(3):
            fd = [(a - b) / (2 * h) for a, b in zip(
                hz.value_call(cls, shift_q(q, ci, comp, h)),
                hz.value_call(cls, shift_q(q, ci, comp, -h)))]
            for o, t in enumerate(outs_val):
                fields = [t.a, t.b, t.c, t.d]
                cen_tuple = fields[ci]
                up_t = list(cen_tuple); up_t[comp] += 1
                up_key = Integral(*[tuple(up_t) if k == ci else fields[k] for k in range(4)], m=0)
                # NOTE: "up" is ALREADY 2*alpha-weighted by the kernel at
                # contraction (spec SS8.4) -- do not re-apply the factor here.
                analytic = secs["up" + cn][up_idx[up_key]]
                if has_dn and cen_tuple[comp] > 0:
                    dn_t = list(cen_tuple); dn_t[comp] -= 1
                    dn_key = Integral(*[tuple(dn_t) if k == ci else fields[k] for k in range(4)], m=0)
                    analytic -= cen_tuple[comp] * secs["dn" + cn][dn_idx[dn_key]]
                diff = abs(analytic - fd[o])
                if diff > worst:
                    worst, worst_desc = diff, (cn, comp, o)
    return worst, worst_desc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", default=None,
                     help="comma-separated class names, e.g. ssss,ssps,pppp,dsds,ssff")
    ap.add_argument("--workdir", default=None)
    args = ap.parse_args()

    if args.classes:
        name_to_cls = {}
        for cls in R.CANON_LADDER:
            name_to_cls[cls_name(cls)] = cls
        # allow classes outside CANON_LADDER too, by letter-decoding
        def decode(nm):
            return tuple(L.index(ch) for ch in nm)
        classes = [name_to_cls.get(n, decode(n)) for n in args.classes.split(",")]
    else:
        classes = [(0, 0, 0, 0), (0, 0, 1, 0), (1, 0, 1, 0), (0, 0, 1, 1),
                   (1, 0, 1, 1), (1, 1, 1, 1), (0, 0, 2, 0), (2, 0, 2, 0),
                   (0, 0, 2, 2), (0, 0, 3, 0), (3, 0, 3, 0), (0, 0, 3, 3)]

    workdir = args.workdir or tempfile.mkdtemp(prefix="gradctr_fd_gate_")
    hz = Harness(classes, workdir)
    q = default_quartet()
    allok = True
    for cls in classes:
        w, desc = check_class(hz, cls, q)
        status = "PASS" if w < TOL else "FAIL"
        if status == "FAIL":
            allok = False
        print(f"{cls_name(cls):6s} {status}  max|analytic-FD|={w:.3e}  worst={desc}")
    print("ALL PASS" if allok else "SOME FAILED")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
