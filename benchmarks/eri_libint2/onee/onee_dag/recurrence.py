"""Obara-Saika recurrence rules for one-electron integrals (HGP-OS only).

Each rule maps a OneE target to a list of Terms whose signed linear
combination equals it. Coefficients are symbolic scalar names the emitter
(and the A2 numerical-reconstruction test) resolve per primitive:

    PAx/PAy/PAz   = P_i - A_i
    PCx/PCy/PCz   = P_i - C_i        (nuclear only)
    ABx/ABy/ABz   = A_i - B_i        (HRR)
    inv_2zeta     = 1 / (2ζ),  ζ = α + β
    one           = 1

A Term also carries an integer `imul` (the angular-momentum prefactor a_i)
so the down-step coefficient is `imul * inv_2zeta`.

Scheme (mirrors the ERI codegen + reference.py): build the (a, 0) tower by
VRR on A, then transfer A→B by HRR. Overlap/nuclear share the same HRR;
their VRR differs (nuclear carries the P-C term and the Boys index m).
Kinetic is a fixed linear combination of overlap nodes (no own recurrence).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List

from .node import OneE, MomE, add_at, first_nonzero

_XYZ = "xyz"


@dataclass(frozen=True)
class Term:
    coef: str            # symbolic scalar name (see module docstring)
    source: object       # OneE (S/T/V) or MomE (multipole moment) DAG node
    sign: int = 1        # +1 / -1
    imul: int = 1        # integer angular prefactor (multiplies coef)

    def __str__(self):
        s = "-" if self.sign < 0 else "+"
        mul = "" if self.imul == 1 else f"{self.imul}*"
        return f"{s} {mul}{self.coef} * {self.source.name()}"


# ----- overlap VRR on A:  (a+1_i,0) = PA_i (a,0) + a_i/2ζ (a-1_i,0) ---------

def vrr_a_overlap(t: OneE) -> List[Term]:
    assert t.L_a > 0 and t.L_b == 0
    i = first_nonzero(t.a)
    am = add_at(t.a, i, -1)
    out = [Term(f"PA{_XYZ[i]}", OneE(a=am, b=t.b, m=t.m))]
    if am[i] > 0:
        out.append(Term("inv_2zeta", OneE(a=add_at(am, i, -1), b=t.b, m=t.m),
                        sign=+1, imul=am[i]))
    return out


# ----- nuclear VRR on A ----------------------------------------------------
#  (a+1_i,0)^m = PA_i (a,0)^m - PC_i (a,0)^{m+1}
#               + a_i/2ζ [ (a-1_i,0)^m - (a-1_i,0)^{m+1} ]

def vrr_a_nuclear(t: OneE) -> List[Term]:
    assert t.L_a > 0 and t.L_b == 0
    i = first_nonzero(t.a)
    am = add_at(t.a, i, -1)
    out = [
        Term(f"PA{_XYZ[i]}", OneE(a=am, b=t.b, m=t.m), sign=+1),
        Term(f"PC{_XYZ[i]}", OneE(a=am, b=t.b, m=t.m + 1), sign=-1),
    ]
    if am[i] > 0:
        amm = add_at(am, i, -1)
        out.append(Term("inv_2zeta", OneE(a=amm, b=t.b, m=t.m),
                        sign=+1, imul=am[i]))
        out.append(Term("inv_2zeta", OneE(a=amm, b=t.b, m=t.m + 1),
                        sign=-1, imul=am[i]))
    return out


# ----- HRR transfer A→B (shared, no m change) ------------------------------
#  (a,b+1_i) = (a+1_i,b) + (A_i-B_i)(a,b)

def hrr_b(t: OneE) -> List[Term]:
    assert t.L_b > 0
    i = first_nonzero(t.b)
    bm = add_at(t.b, i, -1)
    return [
        Term("one", OneE(a=add_at(t.a, i, +1), b=bm, m=t.m), sign=+1),
        Term(f"AB{_XYZ[i]}", OneE(a=t.a, b=bm, m=t.m), sign=+1),
    ]


def expand(t: OneE, op: str) -> List[Term]:
    """Expand a target into lower sources. op in {'overlap','nuclear'}.
    Returns [] for a base node (a==b==0)."""
    if t.is_base():
        return []
    if t.L_b > 0:
        return hrr_b(t)
    if op == "overlap":
        return vrr_a_overlap(t)
    if op == "nuclear":
        return vrr_a_nuclear(t)
    raise ValueError(f"unknown op {op!r}")


# ----- multipole-moment VRR over the moment index e (HJO Eq. 9.5.43) --------
#  M(a,b, e+1_i) = PC_i M(a,b,e)
#                + 1/2ζ [ a_i M(a-1_i,b,e) + b_i M(a,b-1_i,e) + e_i M(a,b,e-1_i) ]
#  PC_i = P_i - C_i (C = moment origin). Base M(a,b,0) ≡ overlap S(a,b).
#  CRITICAL (theorist): the e-down-step prefactor is the SOURCE moment count
#  em[i] = e_i - 1 (so 0,1,2,3 climbing 1->2->3->4), mirroring vrr_a_overlap's
#  am[i] post-decrement convention.

def vrr_e_moment(t: MomE) -> List[Term]:
    assert t.n_mom > 0
    i = first_nonzero(t.e)
    em = add_at(t.e, i, -1)
    out = [Term(f"PC{_XYZ[i]}", MomE(a=t.a, b=t.b, e=em), sign=+1)]
    if t.a[i] > 0:
        out.append(Term("inv_2zeta", MomE(a=add_at(t.a, i, -1), b=t.b, e=em),
                        sign=+1, imul=t.a[i]))
    if t.b[i] > 0:
        out.append(Term("inv_2zeta", MomE(a=t.a, b=add_at(t.b, i, -1), e=em),
                        sign=+1, imul=t.b[i]))
    if em[i] > 0:
        out.append(Term("inv_2zeta", MomE(a=t.a, b=t.b, e=add_at(em, i, -1)),
                        sign=+1, imul=em[i]))
    return out


def expand_moment(t: MomE):
    """Expand a moment target. Returns [] for an e==0 base node (which the
    emitter bridges to the overlap value at (a,b))."""
    if t.is_base():
        return []
    return vrr_e_moment(t)


# ----- kinetic as a combination of overlap nodes ---------------------------
#  T(a,b) = Σ_i [ β(2b_i+1) S(a,b) - 2β² S(a,b+2_i) - ½ b_i(b_i-1) S(a,b-2_i) ]
# Coefficients depend on the runtime exponent β; we return them symbolically
# as (kind, factor_int, source) where the emitter forms the β-dependent value.

@dataclass(frozen=True)
class KinTerm:
    kind: str            # 'beta_2b1' | 'beta2' | 'half_bb1'
    bi: int              # b_i (for beta_2b1 / half_bb1 prefactor)
    source: OneE
    sign: int


def kinetic_terms(a, b) -> List[KinTerm]:
    out: List[KinTerm] = []
    for i in range(3):
        bi = b[i]
        out.append(KinTerm("beta_2b1", bi, OneE(a=a, b=b), sign=+1))
        out.append(KinTerm("beta2", 0, OneE(a=a, b=add_at(b, i, +2)), sign=-1))
        if bi >= 2:
            out.append(KinTerm("half_bb1", bi,
                               OneE(a=a, b=add_at(b, i, -2)), sign=-1))
    return out
