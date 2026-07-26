"""Integral node algebra for the HGP-OS ERI DAG.

An `Integral` represents a single intermediate value of the form
`(a b | c d)^(m)` where each of a, b, c, d is a 3-tuple of nonneg
ints (Cartesian angular momentum components) and m is the auxiliary
Boys index. Nodes are immutable + hashable so they can be DAG keys.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Tuple

# Angular momentum tuple: (lx, ly, lz)
AngTup = Tuple[int, int, int]


@dataclass(frozen=True, order=True)
class Integral:
    """One integral (a b | c d)^(m).

    Sort order is (L_total descending, m descending, a, b, c, d) so
    that topological order naturally emits the highest-aux first and
    bottoms out at (0000)^(m) base integrals."""
    # Field order matters for `order=True`; we want highest m and L
    # first in the natural sort, so we use negative for those keys.
    _neg_L_total: int = field(init=False, repr=False)
    _neg_m:       int = field(init=False, repr=False)
    a: AngTup = (0, 0, 0)
    b: AngTup = (0, 0, 0)
    c: AngTup = (0, 0, 0)
    d: AngTup = (0, 0, 0)
    m: int = 0

    def __post_init__(self):
        L = sum(self.a) + sum(self.b) + sum(self.c) + sum(self.d)
        object.__setattr__(self, "_neg_L_total", -L)
        object.__setattr__(self, "_neg_m", -self.m)

    @property
    def L_a(self) -> int: return sum(self.a)
    @property
    def L_b(self) -> int: return sum(self.b)
    @property
    def L_c(self) -> int: return sum(self.c)
    @property
    def L_d(self) -> int: return sum(self.d)
    @property
    def L_total(self) -> int: return -self._neg_L_total

    def is_base(self) -> bool:
        return self.L_a == 0 and self.L_b == 0 and self.L_c == 0 and self.L_d == 0

    def name(self) -> str:
        """C variable name. Format keeps the integral uniquely identified
        while staying within 60 chars. Example: `s_a000_b010_c100_d000_m1`."""
        def f(t): return f"{t[0]}{t[1]}{t[2]}"
        return f"s_a{f(self.a)}_b{f(self.b)}_c{f(self.c)}_d{f(self.d)}_m{self.m}"


def all_ang_tuples(L: int) -> list[AngTup]:
    """All Cartesian (lx, ly, lz) with lx+ly+lz == L, in canonical
    `osx_g_index` order (lz fastest, lx slowest).

    This matches PySCF / libcint Cartesian ordering — e.g. d-shell:
    [xx, xy, xz, yy, yz, zz]. Verified via mol.ao_labels()."""
    out = []
    for lx in range(L, -1, -1):
        for ly in range(L - lx, -1, -1):
            lz = L - lx - ly
            out.append((lx, ly, lz))
    return out


def _dfact_2nm1(n: int) -> int:
    """(2n−1)!! = 1·3·5···(2n−1) for n ≥ 1; by convention 1 for n ≤ 0.

    Used by the per-Cartesian-component normalization factor below."""
    if n <= 0:
        return 1
    result = 1
    for i in range(1, n + 1):
        result *= 2 * i - 1
    return result


def cart_norm_factor(lx: int, ly: int, lz: int) -> float:
    """Per-Cartesian-component normalization absorbed by PySCF / libcint.

    Each shell carries the shared radial norm N_l(α) = (2α/π)^{3/4} ·
    (4α)^{L/2} / sqrt((2L−1)!!), and the per-component factor is

        N_comp(lx, ly, lz) = sqrt( (2L−1)!! / ((2lx−1)!!·(2ly−1)!!·(2lz−1)!!) )

    For L ≤ 1 all components have N_comp = 1 (so (ss|ss), (sp|ss),
    (pp|ss), (sp|sp) test parity element-wise without this factor).
    For L ≥ 2 the xy / xz / yz mixed components get a sqrt(3); for
    L = 3 the fxyz gets a sqrt(15); etc. Missing this factor
    is what produced the 44% mismatch on the (dd|ss) test before
    this commit."""
    import math
    L = lx + ly + lz
    num = _dfact_2nm1(L)
    den = _dfact_2nm1(lx) * _dfact_2nm1(ly) * _dfact_2nm1(lz)
    return math.sqrt(num / den)


def output_set(la: int, lb: int, lc: int, ld: int) -> list[Integral]:
    """The set of output integrals for a bucket (la, lb, lc, ld) at
    m = 0 — flat list in canonical (a, b, c, d) order."""
    out = []
    for a in all_ang_tuples(la):
        for b in all_ang_tuples(lb):
            for c in all_ang_tuples(lc):
                for d in all_ang_tuples(ld):
                    out.append(Integral(a=a, b=b, c=c, d=d, m=0))
    return out


# --- vector arithmetic on AngTups -----------------------------------------

def add_at(t: AngTup, i: int, k: int = 1) -> AngTup:
    """Return t with t[i] += k. Index i ∈ {0, 1, 2}."""
    return tuple(v + (k if j == i else 0) for j, v in enumerate(t))


def first_nonzero(t: AngTup) -> int:
    """Index of first nonzero component, or -1 if all zero."""
    for i, v in enumerate(t):
        if v > 0:
            return i
    return -1
