"""Validate the fused multipole DAG against an independent analytic reference.

The reference does not use any recurrence: it expands the Cartesian polynomial
about the Gaussian-product centre and integrates term by term, so agreement is
evidence that the recurrence and its scalar mapping are right rather than that
two implementations of the same recursion agree with each other.
"""
from __future__ import annotations

import itertools
import math
import sys

sys.path.insert(0, ".")

from onee_dag import MomE, all_ang_tuples
from onee_dag.dag import evaluate as _unused  # noqa: F401  (documents the sibling API)
import recursum_moment_emit as M


def poly_mul(p, q):
    r = [0.0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0.0:
            continue
        for j, b in enumerate(q):
            r[i + j] += a * b
    return r


def axis_integral(a, b, e, PA, PB, PC, zeta):
    """int (u+PA)^a (u+PB)^b (u+PC)^e exp(-zeta u^2) du."""
    poly = [1.0]
    for n, shift in ((a, PA), (b, PB), (e, PC)):
        for _ in range(n):
            poly = poly_mul(poly, [shift, 1.0])   # (shift + u)
    acc = 0.0
    for n, c in enumerate(poly):
        if c == 0.0 or n % 2:
            continue
        df = 1.0
        for k in range(n - 1, 1, -2):
            df *= k
        acc += c * df / (2.0 * zeta) ** (n // 2) * math.sqrt(math.pi / zeta)
    return acc


def reference_moment(a, b, e, alpha, beta, A, B, C):
    zeta = alpha + beta
    mu = alpha * beta / zeta
    val = 1.0
    for i in range(3):
        P = (alpha * A[i] + beta * B[i]) / zeta
        K = math.exp(-mu * (A[i] - B[i]) ** 2)
        val *= K * axis_integral(a[i], b[i], e[i],
                                 P - A[i], P - B[i], P - C[i], zeta)
    return val


def dag_values(la, lb, max_e, alpha, beta, A, B, C):
    es = M.moment_set(max_e)
    outs, _ = M._outputs_for(la, lb, es)
    dag = M.build_fused_dag(outs)
    zeta = alpha + beta
    P = [(alpha * A[i] + beta * B[i]) / zeta for i in range(3)]
    sc = {"one": 1.0, "inv_2zeta": 1.0 / (2.0 * zeta)}
    for i, ax in enumerate("xyz"):
        sc["PA" + ax] = P[i] - A[i]
        sc["PC" + ax] = P[i] - C[i]
        sc["AB" + ax] = A[i] - B[i]
    mu = alpha * beta / zeta
    S00 = 1.0
    for i in range(3):
        S00 *= math.exp(-mu * (A[i] - B[i]) ** 2) * math.sqrt(math.pi / zeta)

    vals = {}
    for n in dag.topological_order():
        terms = dag.nodes[n]
        if not terms:
            vals[n] = S00
            continue
        acc = 0.0
        for t in terms:
            acc += t.sign * t.imul * sc[t.coef] * vals[t.source]
        vals[n] = acc
    return vals, outs


def main():
    A = (0.11, -0.37, 0.52)
    B = (-0.63, 0.24, -0.18)
    C = (0.31, 0.17, -0.44)      # moment origin, deliberately not on A or B
    alpha, beta = 1.7, 0.83

    # The base overlap sets the scale of everything the recurrence builds, so
    # errors are judged against it. Pure relative error is the wrong criterion
    # for integrals that nearly cancel: a d-d quadrupole component can be 1e-7
    # while the base overlap is 1e-1, and an absolute error at the rounding
    # floor then reads as 1e-12 "relative". Both are reported; the gate is on
    # the scaled absolute error, which is what actually propagates.
    zeta = alpha + beta
    mu = alpha * beta / zeta
    s00 = 1.0
    for i in range(3):
        s00 *= math.exp(-mu * (A[i] - B[i]) ** 2) * math.sqrt(math.pi / zeta)

    worst_rel = 0.0
    worst_scaled = 0.0
    worst_where = None
    checked = 0
    for la, lb in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2),
                   (2, 1), (1, 2), (2, 2)]:
        vals, outs = dag_values(la, lb, 2, alpha, beta, A, B, C)
        for o in outs:
            ref = reference_moment(o.a, o.b, o.e, alpha, beta, A, B, C)
            err = abs(vals[o] - ref)
            worst_rel = max(worst_rel, err / max(abs(ref), 1e-30))
            scaled = err / s00
            if scaled > worst_scaled:
                worst_scaled = scaled
                worst_where = (la, lb, o.a, o.b, o.e)
            checked += 1
    print("checked %d moment integrals over 9 classes" % checked)
    print("worst error / |S00|          : %.3e   at %s" % (worst_scaled, worst_where))
    print("worst pure relative error    : %.3e   (near-cancelling components)"
          % worst_rel)
    ok = worst_scaled < 1e-14
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
