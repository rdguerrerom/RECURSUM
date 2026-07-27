"""Node algebra for the one-electron (S, T, V_ne) Obara-Saika DAG.

A `OneE` node is one intermediate `(a, b)^(m)` where a, b are Cartesian
angular-momentum 3-tuples on centers A, B and m is the Boys auxiliary
index (overlap/kinetic ignore m; nuclear uses it). Nodes are immutable +
hashable so they can be DAG keys. Mirrors gpu/codegen/eri/integral.py.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

AngTup = Tuple[int, int, int]
MomTup = Tuple[int, int, int]


@dataclass(frozen=True, order=True)
class OneE:
    """One intermediate (a, b)^(m). Sort: (L_total desc, m desc, a, b)."""
    _neg_L: int = field(init=False, repr=False)
    _neg_m: int = field(init=False, repr=False)
    a: AngTup = (0, 0, 0)
    b: AngTup = (0, 0, 0)
    m: int = 0

    def __post_init__(self):
        object.__setattr__(self, "_neg_L", -(sum(self.a) + sum(self.b)))
        object.__setattr__(self, "_neg_m", -self.m)

    @property
    def L_a(self) -> int: return sum(self.a)
    @property
    def L_b(self) -> int: return sum(self.b)
    @property
    def L_total(self) -> int: return -self._neg_L

    def is_base(self) -> bool:
        return self.L_a == 0 and self.L_b == 0

    def name(self) -> str:
        def f(t): return f"{t[0]}{t[1]}{t[2]}"
        return f"o_a{f(self.a)}_b{f(self.b)}_m{self.m}"


@dataclass(frozen=True, order=True)
class MomE:
    """One multipole-moment intermediate (a, b; e). a, b are Cartesian
    angular-momentum 3-tuples on centers A, B; e=(ex,ey,ez) is the moment
    multi-index over the common origin C. NO Boys index (moments are
    separable like overlap; theorist HJO §9.5.2). Sort key
    (L_total desc, |e| asc, a, b, e): the overlap floor (|e|=0) and lower-|e|
    layers are emitted before higher |e| so the moment VRR only ever climbs e.
    A base node (n_mom==0) bridges to the overlap value at (a,b)."""
    _neg_L: int = field(init=False, repr=False)
    _e_ord: int = field(init=False, repr=False)
    a: AngTup = (0, 0, 0)
    b: AngTup = (0, 0, 0)
    e: MomTup = (0, 0, 0)

    def __post_init__(self):
        object.__setattr__(self, "_neg_L", -(sum(self.a) + sum(self.b)))
        object.__setattr__(self, "_e_ord", sum(self.e))

    @property
    def L_a(self) -> int: return sum(self.a)
    @property
    def L_b(self) -> int: return sum(self.b)
    @property
    def L_total(self) -> int: return -self._neg_L
    @property
    def n_mom(self) -> int: return self._e_ord

    def is_base(self) -> bool:
        """e==0 -> bottoms out into the overlap value at (a,b)."""
        return self.n_mom == 0

    def name(self) -> str:
        def f(t): return f"{t[0]}{t[1]}{t[2]}"
        return f"m_a{f(self.a)}_b{f(self.b)}_e{f(self.e)}"


def add_at(t: AngTup, i: int, k: int = 1) -> AngTup:
    return tuple(v + (k if j == i else 0) for j, v in enumerate(t))


def first_nonzero(t: AngTup) -> int:
    for i, v in enumerate(t):
        if v > 0:
            return i
    return -1


def all_ang_tuples(L):
    """All Cartesian (lx,ly,lz) with lx+ly+lz==L in libcint/PySCF order
    (lz fastest, lx slowest) — e.g. d: [xx,xy,xz,yy,yz,zz]."""
    out = []
    for lx in range(L, -1, -1):
        for ly in range(L - lx, -1, -1):
            out.append((lx, ly, L - lx - ly))
    return out
