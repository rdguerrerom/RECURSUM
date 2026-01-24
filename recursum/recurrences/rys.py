"""Rys quadrature recurrence relations for electron repulsion integrals."""

from ..codegen import Recurrence


def rys_2d_integral() -> Recurrence:
    """Rys 2D integral vertical recurrence I_x(n, 0, m, 0)."""
    rec = Recurrence("Rys2D", ["n", "m"],
                     ["B00", "B10", "B01", "C00", "C00p"],
                     namespace="rys_quadrature",
                     max_indices={"n": 3, "m": 3})  # Max L=3 (f-orbitals)
    rec.validity("n >= 0", "m >= 0")
    rec.base(n=0, m=0, value=1.0)
    # Augspurger 1990 Eq. 15: I(n,m) = C00*I(n-1,m) + (n-1)*B10*I(n-2,m) + m*B00*I(n-1,m-1)
    rec.rule("n > 0 && m == 0",
             "C00 * E[n-1, m] + (n-1) * B10 * E[n-2, m]",
             name="Bra VRR (m=0)")
    # Augspurger 1990 Eq. 16: I(n,m) = C00'*I(n,m-1) + (m-1)*B01*I(n,m-2) + n*B00*I(n-1,m-1)
    rec.rule("m > 0 && n == 0",
             "C00p * E[n, m-1] + (m-1) * B01 * E[n, m-2]",
             name="Ket VRR (n=0)")
    # General case: use bra-priority recurrence (Augspurger Eq. 15)
    rec.rule("n > 0 && m > 0",
             "C00 * E[n-1, m] + (n-1) * B10 * E[n-2, m] + m * B00 * E[n-1, m-1]",
             name="General (bra priority)")
    return rec


def rys_horizontal_transfer() -> Recurrence:
    """Rys horizontal recurrence (HRR) for angular momentum transfer."""
    rec = Recurrence("RysHRR", ["ix", "jx", "kx", "lx"],
                     ["Axi_Bxi", "Cxi_Dxi"],
                     namespace="rys_quadrature",
                     max_indices={"ix": 3, "jx": 3, "kx": 3, "lx": 3})  # Max L=3 (f-orbitals)
    rec.validity("ix >= 0", "jx >= 0", "kx >= 0", "lx >= 0")
    rec.base(ix=0, jx=0, kx=0, lx=0, value=1.0)
    rec.rule("jx > 0",
             "E[ix+1, jx-1, kx, lx] + Axi_Bxi * E[ix, jx-1, kx, lx]",
             name="Bra transfer")
    rec.rule("lx > 0 && jx == 0",
             "E[ix, jx, kx+1, lx-1] + Cxi_Dxi * E[ix, jx, kx, lx-1]",
             name="Ket transfer")
    return rec


def rys_vrr_full() -> Recurrence:
    """Full 4-index Rys VRR (Head-Gordon & Pople)."""
    rec = Recurrence("RysVRRFull", ["i", "j", "k", "l"],
                     ["PA", "PB", "QC", "QD", "alpha", "beta"],
                     namespace="rys_quadrature",
                     max_indices={"i": 3, "j": 3, "k": 3, "l": 3})  # Max L=3 (f-orbitals)
    rec.validity("i >= 0", "j >= 0", "k >= 0", "l >= 0")
    rec.base(i=0, j=0, k=0, l=0, value=1.0)
    rec.rule("i > 0",
             "PA * E[i-1, j, k, l] + (i-1) * alpha * E[i-2, j, k, l]",
             name="i-increment")
    return rec


def rys_polynomial_recursion() -> Recurrence:
    """Rys polynomial recursion for quadrature roots."""
    rec = Recurrence("RysPoly", ["n"], ["t", "X", "beta"],
                     namespace="rys_poly",
                     max_indices={"n": 6})  # Max order for f-orbitals (2*3)
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.rule("n > 0",
             "(t * t + (-beta)) * E[n-1] + (-beta) * E[n-2]",
             name="Three-term recurrence")
    return rec


__all__ = [
    "rys_2d_integral",
    "rys_horizontal_transfer",
    "rys_vrr_full",
    "rys_polynomial_recursion",
]
