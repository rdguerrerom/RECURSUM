"""HGP-OS recurrence rules.

Each rule maps an Integral with L_? > 0 to a list of (coefficient,
source integral) terms whose linear combination equals the target.
Coefficients are symbolic strings referring to scalar variables that
the emitter is responsible for computing once per quartet (PA_i, WP_i,
inv_2zp, etc.) — see `emit.py::SCALARS`.

Two recurrence families:
  - VRR (Vertical, Obara-Saika): bumps angular momentum on A or C
    while keeping b = d = 0. Increases L_a or L_c by 1 and may
    increase m by 1.
  - HRR (Horizontal, Head-Gordon-Pople): transfers angular momentum
    A→B or C→D. Does NOT change m.

Convention: VRR is applied first (build full (a, 0 | c, 0) tower),
then HRR transfers. This matches the canonical HGP-OS DAG.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

from .integral import Integral, add_at, first_nonzero


@dataclass(frozen=True)
class Term:
    """One additive term in a recurrence expansion."""
    coef: str          # scalar variable name, e.g. "PAx", "inv_2zp"
    source: Integral   # the source integral being multiplied by coef
    sign: int = 1      # +1 or -1; multiplied into the emitted expression

    def __str__(self):
        s = ("-" if self.sign < 0 else "+")
        return f"{s} {self.coef} * {self.source.name()}"


# ----- VRR-A ---------------------------------------------------------------
#
# (a+1_i, 0 | c, 0)^(m) =
#      PA_i · (a, 0 | c, 0)^(m)
#    + WP_i · (a, 0 | c, 0)^(m+1)
#    + (a_i / 2ζ_p) · [(a-1_i, 0 | c, 0)^(m) − (ρ/ζ_p) · (a-1_i, 0 | c, 0)^(m+1)]
#    + (c_i / 2ζ_pq) · (a, 0 | c-1_i, 0)^(m+1)
#
# where ρ = ζ_p ζ_q / (ζ_p + ζ_q). Note (ρ/ζ_p) = ζ_q / (ζ_p+ζ_q) — we
# pre-compute this as `frac_q_over_pq` and a_i / 2ζ_p as `inv_2zp` so
# the rule emits a fixed pattern with up to 5 terms.


def vrr_a(target: Integral) -> List[Term]:
    """Apply VRR on A to a target with L_a > 0 and L_b = 0.

    Returns a list of Term entries summing to `target`."""
    assert target.L_a > 0 and target.L_b == 0, (
        f"vrr_a requires L_a>0 L_b==0: {target}"
    )
    a = target.a
    # Pick the index i to bump down. Convention: first nonzero (matches
    # gen_eri_kernels.py descent order).
    i = first_nonzero(a)
    coord = "xyz"[i]
    a_minus = add_at(a, i, -1)
    src1 = Integral(a=a_minus, b=target.b, c=target.c, d=target.d, m=target.m)
    src2 = Integral(a=a_minus, b=target.b, c=target.c, d=target.d, m=target.m + 1)

    terms: List[Term] = [
        Term(coef=f"PA{coord}", source=src1),
        Term(coef=f"WP{coord}", source=src2),
    ]
    # a_i lower term (a_i = a_minus[i] after bump-down).
    ai_lower = a_minus[i]
    if ai_lower > 0:
        a_minus2 = add_at(a_minus, i, -1)
        s3 = Integral(a=a_minus2, b=target.b, c=target.c, d=target.d, m=target.m)
        s4 = Integral(a=a_minus2, b=target.b, c=target.c, d=target.d, m=target.m + 1)
        # ai_lower · inv_2zp · (s3 - frac_q_over_pq · s4)
        cf = f"ai{ai_lower}_inv_2zp"  # placeholder; emit.py inlines as ai_lower * inv_2zp
        terms.append(Term(coef=cf, source=s3, sign=+1))
        terms.append(Term(coef=cf + "_xq", source=s4, sign=-1))
    # c_i cross term (c is target.c, NOT bumped).
    ci = target.c[i]
    if ci > 0:
        c_minus = add_at(target.c, i, -1)
        s5 = Integral(a=a_minus, b=target.b, c=c_minus, d=target.d, m=target.m + 1)
        cf = f"ci{ci}_inv_2zpq"
        terms.append(Term(coef=cf, source=s5))
    return terms


# ----- VRR-C ---------------------------------------------------------------
#
# Symmetric to VRR-A under (a↔c, b↔d, ζ_p↔ζ_q, PA↔QC, WP↔WQ).


def vrr_c(target: Integral) -> List[Term]:
    """Apply VRR on C to a target with L_c > 0 and L_d = 0."""
    assert target.L_c > 0 and target.L_d == 0, (
        f"vrr_c requires L_c>0 L_d==0: {target}"
    )
    c = target.c
    i = first_nonzero(c)
    coord = "xyz"[i]
    c_minus = add_at(c, i, -1)
    src1 = Integral(a=target.a, b=target.b, c=c_minus, d=target.d, m=target.m)
    src2 = Integral(a=target.a, b=target.b, c=c_minus, d=target.d, m=target.m + 1)

    terms: List[Term] = [
        Term(coef=f"QC{coord}", source=src1),
        Term(coef=f"WQ{coord}", source=src2),
    ]
    ci_lower = c_minus[i]
    if ci_lower > 0:
        c_minus2 = add_at(c_minus, i, -1)
        s3 = Integral(a=target.a, b=target.b, c=c_minus2, d=target.d, m=target.m)
        s4 = Integral(a=target.a, b=target.b, c=c_minus2, d=target.d, m=target.m + 1)
        cf = f"ci{ci_lower}_inv_2zq"
        terms.append(Term(coef=cf, source=s3, sign=+1))
        terms.append(Term(coef=cf + "_xp", source=s4, sign=-1))
    ai = target.a[i]
    if ai > 0:
        a_minus = add_at(target.a, i, -1)
        s5 = Integral(a=a_minus, b=target.b, c=c_minus, d=target.d, m=target.m + 1)
        cf = f"ai{ai}_inv_2zpq"
        terms.append(Term(coef=cf, source=s5))
    return terms


# ----- HRR-AB --------------------------------------------------------------
#
# (a, b+1_i | c, d) = (a+1_i, b | c, d) − (B-A)_i · (a, b | c, d)
#
# Used to transfer L from A to B once VRR has produced (a, 0 | c, d).


def hrr_ab(target: Integral) -> List[Term]:
    """Apply HRR on B to a target with L_b > 0."""
    assert target.L_b > 0, f"hrr_ab requires L_b>0: {target}"
    b = target.b
    i = first_nonzero(b)
    coord = "xyz"[i]
    b_minus = add_at(b, i, -1)
    a_plus = add_at(target.a, i, +1)
    src1 = Integral(a=a_plus, b=b_minus, c=target.c, d=target.d, m=target.m)
    src2 = Integral(a=target.a, b=b_minus, c=target.c, d=target.d, m=target.m)
    return [
        Term(coef="one", source=src1, sign=+1),
        Term(coef=f"AB{coord}", source=src2, sign=-1),
    ]


# ----- HRR-CD --------------------------------------------------------------


def hrr_cd(target: Integral) -> List[Term]:
    """Apply HRR on D to a target with L_d > 0."""
    assert target.L_d > 0, f"hrr_cd requires L_d>0: {target}"
    d = target.d
    i = first_nonzero(d)
    coord = "xyz"[i]
    d_minus = add_at(d, i, -1)
    c_plus = add_at(target.c, i, +1)
    src1 = Integral(a=target.a, b=target.b, c=c_plus, d=d_minus, m=target.m)
    src2 = Integral(a=target.a, b=target.b, c=target.c, d=d_minus, m=target.m)
    return [
        Term(coef="one", source=src1, sign=+1),
        Term(coef=f"CD{coord}", source=src2, sign=-1),
    ]


# ----- Dispatcher ----------------------------------------------------------
#
# Decision tree: HRR before VRR so that we transfer angular momentum
# back to B/D before building higher A/C towers. This matches the
# canonical HGP-OS scheme and minimizes the number of intermediates.


def expand(target: Integral) -> List[Term]:
    """Apply the right recurrence to expand `target` into a sum of
    lower-L sources. Returns [] if target is a base integral."""
    if target.is_base():
        return []
    # HRR first: transfer L from A to B / C to D when possible.
    if target.L_b > 0:
        return hrr_ab(target)
    if target.L_d > 0:
        return hrr_cd(target)
    # Then VRR on A or C.
    if target.L_a > 0:
        return vrr_a(target)
    if target.L_c > 0:
        return vrr_c(target)
    return []
