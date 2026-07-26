"""HGP-OS primitive-quartet ERI oracle.

Two responsibilities:
  1. compute_scalars(...)  -> the exact per-quartet scalars the osx emitter uses
     (emit.py:215-250 + base-case K/T at emit.py:315-331). These are the SAME
     names/definitions the RECURSUM-generated kernel will consume, so a kernel
     fed these scalars + this Boys array must reproduce eval_dag() bit-for-bit
     (level-1 correctness gate, independent of physics conventions).
  2. eval_dag(...)         -> evaluate the osx HGP-OS DAG for a class at a
     primitive quartet, returning the Cartesian output integrals in canonical
     order. This is the reference the generated kernel is validated against.

A separate script (pyscf_reference.py) validates that eval_dag() matches
PySCF/libcint up to primitive normalization -- the physics gate that also
guarantees libint2 computes identical numbers.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple

from scipy.special import gammainc, gamma  # regularized lower incomplete gamma

from eri_dag.integral import Integral, output_set, all_ang_tuples
from eri_dag.dag import build_dag
from eri_dag.recurrence import Term

Vec3 = Tuple[float, float, float]

TWO_PI_5_2 = 34.98683665524972497057  # 2 * pi^(5/2)


def boys(m_max: int, T: float) -> List[float]:
    """F_m(T) for m=0..m_max. F_m(T) = P(m+0.5,T)*Gamma(m+0.5)/(2 T^{m+0.5})."""
    out = []
    if T < 1e-13:
        return [1.0 / (2 * m + 1) for m in range(m_max + 1)]
    for m in range(m_max + 1):
        a = m + 0.5
        # P(a,T) is the regularized lower incomplete gamma; multiply back Gamma(a)
        val = gammainc(a, T) * gamma(a) / (2.0 * T ** a)
        out.append(float(val))
    return out


@dataclass
class Quartet:
    a: float; A: Vec3          # exponent + center on A
    b: float; B: Vec3
    g: float; C: Vec3
    d: float; D: Vec3


def compute_scalars(q: Quartet) -> Dict[str, float]:
    """Exact transcription of emit.py:215-250 + K/T (emit.py:326-331)."""
    a, b, g, d = q.a, q.b, q.g, q.d
    Ax, Ay, Az = q.A; Bx, By, Bz = q.B; Cx, Cy, Cz = q.C; Dx, Dy, Dz = q.D

    zeta_p = a + b
    zeta_q = g + d
    zeta_pq = zeta_p + zeta_q
    rho = zeta_p * zeta_q / zeta_pq
    pref_AB = a * b / zeta_p
    pref_CD = g * d / zeta_q

    dABx, dABy, dABz = Ax - Bx, Ay - By, Az - Bz
    dCDx, dCDy, dCDz = Cx - Dx, Cy - Dy, Cz - Dz
    R2_AB = dABx * dABx + dABy * dABy + dABz * dABz
    R2_CD = dCDx * dCDx + dCDy * dCDy + dCDz * dCDz

    Px = (a * Ax + b * Bx) / zeta_p
    Py = (a * Ay + b * By) / zeta_p
    Pz = (a * Az + b * Bz) / zeta_p
    Qx = (g * Cx + d * Dx) / zeta_q
    Qy = (g * Cy + d * Dy) / zeta_q
    Qz = (g * Cz + d * Dz) / zeta_q
    Wx = (zeta_p * Px + zeta_q * Qx) / zeta_pq
    Wy = (zeta_p * Py + zeta_q * Qy) / zeta_pq
    Wz = (zeta_p * Pz + zeta_q * Qz) / zeta_pq

    dPQx, dPQy, dPQz = Px - Qx, Py - Qy, Pz - Qz
    R2_PQ = dPQx * dPQx + dPQy * dPQy + dPQz * dPQz

    s = {
        "PAx": Px - Ax, "PAy": Py - Ay, "PAz": Pz - Az,
        "QCx": Qx - Cx, "QCy": Qy - Cy, "QCz": Qz - Cz,
        "WPx": Wx - Px, "WPy": Wy - Py, "WPz": Wz - Pz,
        "WQx": Wx - Qx, "WQy": Wy - Qy, "WQz": Wz - Qz,
        # NB: emit.py:244-245 defines these as (B-A) and (D-C), not (A-B)/(C-D).
        "ABx": Bx - Ax, "ABy": By - Ay, "ABz": Bz - Az,
        "CDx": Dx - Cx, "CDy": Dy - Cy, "CDz": Dz - Cz,
        "inv_2zp": 0.5 / zeta_p, "inv_2zq": 0.5 / zeta_q, "inv_2zpq": 0.5 / zeta_pq,
        "frac_q_over_pq": zeta_q / zeta_pq, "frac_p_over_pq": zeta_p / zeta_pq,
        "one": 1.0,
    }
    # Base-case K prefactor and Boys argument (emit.py:326-331)
    geom = TWO_PI_5_2 / (zeta_p * zeta_q * math.sqrt(zeta_pq))
    K = geom * math.exp(-pref_AB * R2_AB) * math.exp(-pref_CD * R2_CD)
    T = rho * R2_PQ
    s["_K"] = K
    s["_T"] = T
    return s


# --- coefficient-string resolver ------------------------------------------
# recurrence.py emits both plain scalar names and composite placeholder names
# like "ai2_inv_2zp", "ai2_inv_2zp_xq", "ci1_inv_2zpq". Resolve them to floats.
_PATTERNS = [
    (re.compile(r"^ai(\d+)_inv_2zp_xq$"),  lambda n, s: n * s["inv_2zp"] * s["frac_q_over_pq"]),
    (re.compile(r"^ai(\d+)_inv_2zp$"),     lambda n, s: n * s["inv_2zp"]),
    (re.compile(r"^ci(\d+)_inv_2zq_xp$"),  lambda n, s: n * s["inv_2zq"] * s["frac_p_over_pq"]),
    (re.compile(r"^ci(\d+)_inv_2zq$"),     lambda n, s: n * s["inv_2zq"]),
    (re.compile(r"^ci(\d+)_inv_2zpq$"),    lambda n, s: n * s["inv_2zpq"]),
    (re.compile(r"^ai(\d+)_inv_2zpq$"),    lambda n, s: n * s["inv_2zpq"]),
]


def resolve_coef(coef: str, s: Dict[str, float]) -> float:
    if coef in s:
        return s[coef]
    for pat, fn in _PATTERNS:
        m = pat.match(coef)
        if m:
            return fn(int(m.group(1)), s)
    raise KeyError(f"unresolved coefficient string: {coef!r}")


def eval_dag(la: int, lb: int, lc: int, ld: int, q: Quartet) -> Tuple[List[Integral], List[float]]:
    """Evaluate the HGP-OS DAG for class (la lb|lc ld) at quartet q.

    Returns (ordered_output_integrals, values) with values in canonical
    Cartesian order (all_ang_tuples x all_ang_tuples x ...)."""
    s = compute_scalars(q)
    outs = output_set(la, lb, lc, ld)
    dag = build_dag(outs)
    m_max = dag.max_m()
    F = boys(m_max, s["_T"])
    K = s["_K"]

    val: Dict[Integral, float] = {}
    for node in dag.topological_order():
        terms = dag.nodes[node]
        if not terms:  # base integral (0000)^(m)
            val[node] = K * F[node.m]
        else:
            acc = 0.0
            for t in terms:
                acc += t.sign * resolve_coef(t.coef, s) * val[t.source]
            val[node] = acc
    return outs, [val[o] for o in outs]


if __name__ == "__main__":
    # smoke test: a well-separated, non-degenerate quartet
    q = Quartet(a=1.2, A=(0.0, 0.0, 0.0),
                b=0.8, B=(0.5, 0.0, 0.0),
                g=1.5, C=(0.0, 0.7, 0.0),
                d=0.9, D=(0.0, 0.0, 0.9))
    for cls in [(0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 1), (2, 0, 2, 0)]:
        outs, vals = eval_dag(*cls, q)
        print(f"class {cls}: n_out={len(vals):3d}  first={vals[0]:.10e}  "
              f"sum={sum(vals):.10e}")
