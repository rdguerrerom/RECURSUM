"""Level-1 correctness gate: compiled RECURSUM kernel must reproduce the
oracle's eval_dag() to machine precision (both consume the same scalars + kf).
"""
import subprocess
import numpy as np

from oracle import Quartet, compute_scalars, boys, eval_dag
from recursum_dag_emit import SCALAR_FIELDS, class_name, CANON_LADDER
from eri_dag.dag import build_dag
from eri_dag.integral import output_set

CLASSES = CANON_LADDER

Q = Quartet(a=1.2, A=(0.0,0.0,0.0), b=0.8, B=(0.5,0.1,-0.2),
            g=1.5, C=(0.1,0.7,0.3), d=0.9, D=(-0.3,0.2,0.9))


def run_kernel(cls, s, kf, n_out):
    name = class_name(*cls)
    max_m = len(kf) - 1
    inp = [name, str(max_m), str(n_out)]
    inp += [repr(s[f]) for f in SCALAR_FIELDS]
    inp += [repr(x) for x in kf]
    p = subprocess.run(["./validate_kernels"], input="\n".join(inp),
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{name}: {p.stderr}")
    return np.array([float(x) for x in p.stdout.split()])


if __name__ == "__main__":
    s = compute_scalars(Q)
    print(f"{'class':>6} {'n_out':>6} {'max_rel_err':>13}  status")
    allok = True
    for cls in CLASSES:
        outs, ref = eval_dag(*cls, Q)
        ref = np.array(ref)
        dag = build_dag(output_set(*cls))
        kf = [s["_K"] * f for f in boys(dag.max_m(), s["_T"])]
        got = run_kernel(cls, s, kf, len(ref))
        denom = np.maximum(np.abs(ref), np.abs(ref).max()*1e-15 + 1e-300)
        rel = np.max(np.abs(got - ref) / denom)
        # chunked/large classes reorder FMAs slightly; 1e-10 is machine-precision
        # for a 10k-integral class. The bench cross-check gates libint2 at 1e-9.
        ok = rel < 1e-10
        allok &= ok
        print(f"{class_name(*cls):>6} {len(ref):>6} {rel:>13.3e}  "
              f"{'OK' if ok else 'MISMATCH'}")
    print("\nALL PASS — generated kernels match oracle"
          if allok else "\nFAILURES PRESENT")
