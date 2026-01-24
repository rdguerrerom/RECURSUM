"""Orthogonal polynomial recurrence relations."""

from ..codegen import Recurrence


def hermite_coefficients() -> Recurrence:
    """McMurchie-Davidson Hermite coefficients E^{i,j}_t."""
    rec = Recurrence("Hermite", ["nA", "nB", "N"], ["PA", "PB", "aAB"],
                     namespace="hermite", max_indices={"nA": 3, "nB": 3, "N": 6})  # Max L=3 (f-orbitals)
    rec.validity("nA >= 0", "nB >= 0", "N >= 0", "nA + nB >= N")
    rec.base(nA=0, nB=0, N=0, value=1.0)
    rec.rule("nA == 0 && nB > 0",
             "aAB * E[nA, nB-1, N-1] + PB * E[nA, nB-1, N] + (N+1) * E[nA, nB-1, N+1]",
             name="B-side reduction")
    rec.rule("nB == 0 && nA > 0",
             "aAB * E[nA-1, nB, N-1] + PA * E[nA-1, nB, N] + (N+1) * E[nA-1, nB, N+1]",
             name="A-side reduction")
    rec.branch_average(
        "nA > 0 && nB > 0",
        ["aAB * E[nA, nB-1, N-1] + PB * E[nA, nB-1, N] + (N+1) * E[nA, nB-1, N+1]",
         "aAB * E[nA-1, nB, N-1] + PA * E[nA-1, nB, N] + (N+1) * E[nA-1, nB, N+1]"],
        name="Two-branch average")
    return rec


def legendre_polynomials() -> Recurrence:
    """Legendre polynomials P_n(x).

    DLMF section 18.9 - https://dlmf.nist.gov/18.9
    """
    rec = Recurrence("Legendre", ["n"], ["x"], namespace="legendre",
                     max_indices={"n": 15}, scipy_reference="scipy.special.eval_legendre")
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.base(n=1, value="x")
    rec.rule("n > 1", "(2*n-1) * x * E[n-1] + (-(n-1)) * E[n-2]",
             scale="1/n", name="Three-term recurrence")
    return rec


def chebyshev_T() -> Recurrence:
    """Chebyshev polynomials of the first kind T_n(x).

    DLMF section 18.3 - https://dlmf.nist.gov/18.3
    """
    rec = Recurrence("ChebyshevT", ["n"], ["x"], namespace="chebyshev",
                     max_indices={"n": 15}, scipy_reference="scipy.special.eval_chebyt")
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.base(n=1, value="x")
    rec.rule("n > 1", "(2) * x * E[n-1] + (-1) * E[n-2]", name="Three-term")
    return rec


def chebyshev_U() -> Recurrence:
    """Chebyshev polynomials of the second kind U_n(x).

    DLMF section 18.3 - https://dlmf.nist.gov/18.3
    """
    rec = Recurrence("ChebyshevU", ["n"], ["x", "two_x"], namespace="chebyshev",
                     max_indices={"n": 15}, scipy_reference="scipy.special.eval_chebyu")
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.base(n=1, value="two_x")
    rec.rule("n > 1", "two_x * E[n-1] + (-1) * E[n-2]", name="Three-term")
    return rec


def hermite_He() -> Recurrence:
    """Probabilist's Hermite polynomials He_n(x).

    DLMF section 18.3 - https://dlmf.nist.gov/18.3
    """
    rec = Recurrence("HermiteHe", ["n"], ["x"], namespace="hermite_poly",
                     max_indices={"n": 15}, scipy_reference="scipy.special.eval_hermite")
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.base(n=1, value="x")
    rec.rule("n > 1", "x * E[n-1] + (-(n-1)) * E[n-2]", name="Three-term")
    return rec


def hermite_H() -> Recurrence:
    """Physicist's Hermite polynomials H_n(x).

    DLMF section 18.3 - https://dlmf.nist.gov/18.3
    """
    rec = Recurrence("HermiteH", ["n"], ["x", "two_x"], namespace="hermite_poly",
                     max_indices={"n": 15})
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.base(n=1, value="two_x")
    rec.rule("n > 1", "two_x * E[n-1] + (-2*(n-1)) * E[n-2]", name="Three-term")
    return rec


def laguerre() -> Recurrence:
    """Laguerre polynomials L_n(x).

    DLMF section 18.3 - https://dlmf.nist.gov/18.3
    """
    rec = Recurrence("Laguerre", ["n"], ["x", "one_minus_x"], namespace="laguerre",
                     max_indices={"n": 15}, scipy_reference="scipy.special.eval_laguerre")
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.base(n=1, value="one_minus_x")
    rec.rule("n > 1", "(2*n-1-x) * E[n-1] + (-(n-1)) * E[n-2]", scale="1/n", name="Three-term")
    return rec


def associated_legendre() -> Recurrence:
    """Associated Legendre functions P_l^m(x).

    DLMF section 14.7 - https://dlmf.nist.gov/14.7
    """
    rec = Recurrence("AssocLegendre", ["l", "m"], ["x", "sqrt1mx2"], namespace="legendre",
                     max_indices={"l": 10, "m": 10})
    rec.validity("l >= 0", "m >= 0", "l >= m")
    rec.base(l=0, m=0, value=1.0)
    rec.rule("l == m && m > 0", "(-(2*m-1)) * sqrt1mx2 * E[l-1, m-1]", name="Diagonal")
    rec.rule("l == m + 1", "(2*m+1) * x * E[l-1, m]", name="First off-diagonal")
    rec.rule("l > m + 1", "(2*l-1) * x * E[l-1, m] + (-(l+m-1)) * E[l-2, m]",
             scale="1/(l-m)", name="General")
    return rec


__all__ = [
    "hermite_coefficients",
    "legendre_polynomials",
    "chebyshev_T",
    "chebyshev_U",
    "hermite_He",
    "hermite_H",
    "laguerre",
    "associated_legendre",
]
