"""Quantum chemistry recurrence relations."""

from ..codegen import Recurrence


def sto_auxiliary_B() -> Recurrence:
    """Filter-Steinborn auxiliary function B_{n,l}(x) for STO integrals."""
    rec = Recurrence("STOAuxB", ["n", "l"], ["x", "inv_x", "B00"],
                     namespace="sto", max_indices={"n": 6, "l": 3})  # Max L=3 (f-orbitals)
    rec.validity("n >= 0", "l >= 0", "n >= l")
    rec.base(n=0, l=0, value="B00")
    rec.rule("l == 0 && n > 0",
             "E[n-1, l] + (-(2*n-1)) * inv_x * E[n-1, l]",
             name="l=0 vertical")
    rec.rule("n == l && l > 0",
             "(-(2*l-1)) * inv_x * E[n-1, l-1]",
             name="Diagonal")
    rec.rule("n > l && l > 0",
             "E[n-1, l] + (-(2*n-1)) * inv_x * E[n-1, l] + E[n-1, l-1]",
             name="General")
    return rec


def boys_function() -> Recurrence:
    """Boys function F_n(T) for Gaussian integrals.

    Related to incomplete gamma function - DLMF section 8.2
    https://dlmf.nist.gov/8.2
    """
    rec = Recurrence("Boys", ["n"], ["T", "F0"],
                     namespace="boys_func", max_indices={"n": 6})  # Max order for f-orbitals (2*3)
    rec.validity("n >= 0")
    rec.base(n=0, value="F0")
    # Simplified: F_n = ((2n-1)*F_{n-1} - T*F_{n-1}) / (2n)
    rec.rule("n > 0",
             "(2*n-1) * E[n-1] + (-T) * E[n-1]",
             scale="1/(2*n)",
             name="Downward recurrence")
    return rec


def gaunt_coefficients() -> Recurrence:
    """
    Gaunt coefficients (integrals of three spherical harmonics).

    G(l1,l2,L; m1,m2,M) = ∫ Y_{l1}^{m1} Y_{l2}^{m2} Y_L^M dΩ

    These appear in the angular part of multi-center integrals in quantum chemistry
    (e.g., STO integrals, multipole expansions) and atomic/molecular physics.

    Selection rules:
        - |l1-l2| ≤ L ≤ l1+l2 (triangle inequality)
        - m1 + m2 + M = 0 (magnetic quantum number conservation)
        - l1 + l2 + L must be even (parity conservation)

    This is a simplified recurrence over angular momenta using Clebsch-Gordan-like
    coefficients. Full implementation requires handling of all 6 quantum numbers.

    Args:
        l1, l2, L: Angular momentum quantum numbers
        c1, c2: Clebsch-Gordan-like coupling coefficients (precomputed)

    Reference:
        Gaunt, Trans. Roy. Soc. London A228, 151 (1929)
        Xu, J. Comp. Phys. 139, 137 (1998) - efficient recursion schemes
    """
    rec = Recurrence("Gaunt", ["l1", "l2", "L"], ["c1", "c2"],
                     namespace="angular", max_indices={"l1": 3, "l2": 3, "L": 6})  # Max L=3 (f-orbitals)
    rec.validity("l1 >= 0", "l2 >= 0", "L >= 0")

    # Base case: G(0,0,0) = √(1/4π) integral of three Y_0^0
    rec.base(l1=0, l2=0, L=0, value=1.0)

    # l2=0 reduction (one angular momentum is zero)
    rec.rule("l1 > 0 && l2 == 0",
             "c1 * E[l1-1, l2, L-1] + c2 * E[l1-1, l2, L+1]",
             name="l2=0 reduction")

    # General reduction (l2 > 0)
    rec.rule("l2 > 0",
             "c1 * E[l1, l2-1, L-1] + c2 * E[l1, l2-1, L+1]",
             name="General reduction")

    return rec


__all__ = ["sto_auxiliary_B", "boys_function", "gaunt_coefficients"]
