"""Special function recurrence relations.

This module contains recurrence relations for various special functions
including Jacobi, Gegenbauer, Airy, and Bessel functions.
"""

from ..codegen import Recurrence


def jacobi_polynomials() -> Recurrence:
    """Jacobi polynomials P_n^(alpha,beta)(x).

    DLMF section 18.9 - https://dlmf.nist.gov/18.9

    Three-term recurrence relation for Jacobi polynomials with parameters alpha, beta.
    The Jacobi polynomials generalize many classical orthogonal polynomials.

    Note: P_0 = 1, P_1(x) depends on parameters - handled in dispatcher.
    """
    rec = Recurrence("Jacobi", ["n"], ["x", "alpha", "beta", "P0", "P1"],
                     namespace="jacobi", max_indices={"n": 12},
                     scipy_reference="scipy.special.eval_jacobi")
    rec.validity("n >= 0")
    rec.base(n=0, value="P0")
    rec.base(n=1, value="P1")
    # DLMF 18.9.2: (2n+α+β)(2n+α+β-1)(1-x²)P_n = ...
    # Simplified three-term recurrence: DLMF 18.9.3
    rec.rule("n > 1",
             "(2*n + alpha + beta - 1) * ((2*n + alpha + beta) * (2*n + alpha + beta - 2) * x + alpha*alpha - beta*beta) * E[n-1] + "
             "(-2) * (n + alpha - 1) * (n + beta - 1) * (2*n + alpha + beta) * E[n-2]",
             scale="1/(2*n*(n + alpha + beta)*(2*n + alpha + beta - 2))",
             name="Three-term recurrence")
    return rec


def gegenbauer_polynomials() -> Recurrence:
    """Gegenbauer (ultraspherical) polynomials C_n^(lambda)(x).

    DLMF section 18.3 - https://dlmf.nist.gov/18.3

    Gegenbauer polynomials are a generalization of Legendre polynomials.
    Also called ultraspherical polynomials.

    Note: C_0 = 1, C_1(x) = 2*lambda*x depends on parameter - handled in dispatcher.
    """
    rec = Recurrence("Gegenbauer", ["n"], ["x", "lambda", "C0", "C1"],
                     namespace="gegenbauer", max_indices={"n": 12},
                     scipy_reference="scipy.special.eval_gegenbauer")
    rec.validity("n >= 0")
    rec.base(n=0, value="C0")
    rec.base(n=1, value="C1")
    # DLMF 18.9.2 specialized to Gegenbauer: (n+1)C_{n+1} = 2(n+λ)x C_n - (n+2λ-1)C_{n-1}
    rec.rule("n > 1",
             "2 * (n + lambda - 1) * x * E[n-1] + (-(n + 2*lambda - 2)) * E[n-2]",
             scale="1/n",
             name="Three-term recurrence")
    return rec


def associated_laguerre_polynomials() -> Recurrence:
    """Associated Laguerre polynomials L_n^(alpha)(x).

    DLMF section 18.3 - https://dlmf.nist.gov/18.3

    Generalization of Laguerre polynomials with parameter alpha.

    Note: L_0 = 1, L_1(x) = 1 + alpha - x depends on parameter - handled in dispatcher.
    """
    rec = Recurrence("AssocLaguerre", ["n"], ["x", "alpha", "L0", "L1"],
                     namespace="laguerre", max_indices={"n": 12},
                     scipy_reference="scipy.special.eval_genlaguerre")
    rec.validity("n >= 0")
    rec.base(n=0, value="L0")
    rec.base(n=1, value="L1")
    # DLMF 18.9.2: (n+1)L_{n+1}^α = (2n+α+1-x)L_n^α - (n+α)L_{n-1}^α
    rec.rule("n > 1",
             "(2*n + alpha - 1 - x) * E[n-1] + (-(n + alpha - 1)) * E[n-2]",
             scale="1/n",
             name="Three-term recurrence")
    return rec


def airy_ai() -> Recurrence:
    """Airy function Ai(x) via recurrence for derivatives.

    DLMF section 9.2 - https://dlmf.nist.gov/9.2

    The Airy function satisfies Ai''(x) = x*Ai(x).
    This recurrence computes derivatives Ai^(n)(x).
    """
    rec = Recurrence("AiryAi", ["n"], ["x", "Ai0", "Ai1"],
                     namespace="airy", max_indices={"n": 10})
    rec.validity("n >= 0")
    rec.base(n=0, value="Ai0")
    rec.base(n=1, value="Ai1")
    # Differentiation of Ai''(x) = x*Ai(x) gives: Ai^(n+2) = x*Ai^(n) + n*Ai^(n-1)
    rec.rule("n > 1",
             "x * E[n-2] + (n - 1) * E[n-3]",
             name="Airy differential recurrence")
    return rec


def airy_bi() -> Recurrence:
    """Airy function Bi(x) via recurrence for derivatives.

    DLMF section 9.2 - https://dlmf.nist.gov/9.2

    The Airy function satisfies Bi''(x) = x*Bi(x).
    This recurrence computes derivatives Bi^(n)(x).
    """
    rec = Recurrence("AiryBi", ["n"], ["x", "Bi0", "Bi1"],
                     namespace="airy", max_indices={"n": 10})
    rec.validity("n >= 0")
    rec.base(n=0, value="Bi0")
    rec.base(n=1, value="Bi1")
    # Same recurrence as Ai but for Bi
    rec.rule("n > 1",
             "x * E[n-2] + (n - 1) * E[n-3]",
             name="Airy differential recurrence")
    return rec


def bessel_j() -> Recurrence:
    """Bessel function of the first kind J_n(x).

    DLMF section 10.2 - https://dlmf.nist.gov/10.2

    Standard Bessel function of the first kind.
    """
    rec = Recurrence("BesselJ", ["n"], ["x", "inv_x", "J0", "J1"],
                     namespace="bessel", max_indices={"n": 15},
                     scipy_reference="scipy.special.jv")
    rec.validity("n >= 0")
    rec.base(n=0, value="J0")
    rec.base(n=1, value="J1")
    # DLMF 10.2.19: J_{n+1}(x) = (2n/x)J_n(x) - J_{n-1}(x)
    rec.rule("n > 1",
             "(2 * (n - 1)) * inv_x * E[n-1] + (-1) * E[n-2]",
             name="Bessel upward recurrence")
    return rec


def bessel_y() -> Recurrence:
    """Bessel function of the second kind Y_n(x).

    DLMF section 10.2 - https://dlmf.nist.gov/10.2

    Standard Bessel function of the second kind (Neumann function).
    """
    rec = Recurrence("BesselY", ["n"], ["x", "inv_x", "Y0", "Y1"],
                     namespace="bessel", max_indices={"n": 15},
                     scipy_reference="scipy.special.yv")
    rec.validity("n >= 0")
    rec.base(n=0, value="Y0")
    rec.base(n=1, value="Y1")
    # Same recurrence as J_n
    rec.rule("n > 1",
             "(2 * (n - 1)) * inv_x * E[n-1] + (-1) * E[n-2]",
             name="Bessel upward recurrence")
    return rec


def spherical_bessel_j() -> Recurrence:
    """Spherical Bessel function of the first kind j_n(x).

    DLMF section 10.47 - https://dlmf.nist.gov/10.47

    Spherical Bessel functions: j_n(x) = sqrt(π/(2x)) J_{n+1/2}(x).
    """
    rec = Recurrence("SphericalBesselJ", ["n"], ["inv_x", "j0", "j1"],
                     namespace="bessel", max_indices={"n": 15},
                     scipy_reference="scipy.special.spherical_jn")
    rec.validity("n >= 0")
    rec.base(n=0, value="j0")
    rec.base(n=1, value="j1")
    # DLMF 10.51.2: j_{n+1}(x) = (2n+1)/x j_n(x) - j_{n-1}(x)
    rec.rule("n > 1",
             "(2 * (n - 1) + 1) * inv_x * E[n-1] + (-1) * E[n-2]",
             name="Spherical Bessel upward")
    return rec


def spherical_bessel_y() -> Recurrence:
    """Spherical Bessel function of the second kind y_n(x).

    DLMF section 10.47 - https://dlmf.nist.gov/10.47

    Spherical Bessel functions: y_n(x) = sqrt(π/(2x)) Y_{n+1/2}(x).
    """
    rec = Recurrence("SphericalBesselY", ["n"], ["inv_x", "y0", "y1"],
                     namespace="bessel", max_indices={"n": 15},
                     scipy_reference="scipy.special.spherical_yn")
    rec.validity("n >= 0")
    rec.base(n=0, value="y0")
    rec.base(n=1, value="y1")
    # Same recurrence as j_n
    rec.rule("n > 1",
             "(2 * (n - 1) + 1) * inv_x * E[n-1] + (-1) * E[n-2]",
             name="Spherical Bessel upward")
    return rec


def modified_bessel_i() -> Recurrence:
    """Modified Bessel function of the first kind I_n(x).

    DLMF section 10.29 - https://dlmf.nist.gov/10.29

    Modified Bessel function I_n(x) = i^(-n) J_n(ix).
    """
    rec = Recurrence("ModifiedBesselI", ["n"], ["inv_x", "I0", "I1"],
                     namespace="bessel", max_indices={"n": 15},
                     scipy_reference="scipy.special.iv")
    rec.validity("n >= 0")
    rec.base(n=0, value="I0")
    rec.base(n=1, value="I1")
    # DLMF 10.29.2: I_{n+1}(x) = I_{n-1}(x) - (2n/x)I_n(x)
    rec.rule("n > 1",
             "E[n-2] + (-(2 * (n - 1))) * inv_x * E[n-1]",
             name="Modified Bessel I upward")
    return rec


def modified_bessel_k() -> Recurrence:
    """Modified Bessel function of the second kind K_n(x).

    DLMF section 10.29 - https://dlmf.nist.gov/10.29

    Modified Bessel function K_n(x).
    """
    rec = Recurrence("ModifiedBesselK", ["n"], ["inv_x", "K0", "K1"],
                     namespace="bessel", max_indices={"n": 15},
                     scipy_reference="scipy.special.kv")
    rec.validity("n >= 0")
    rec.base(n=0, value="K0")
    rec.base(n=1, value="K1")
    # DLMF 10.29.2: K_{n+1}(x) = K_{n-1}(x) + (2n/x)K_n(x)
    rec.rule("n > 1",
             "E[n-2] + (2 * (n - 1)) * inv_x * E[n-1]",
             name="Modified Bessel K upward")
    return rec


__all__ = [
    "jacobi_polynomials",
    "gegenbauer_polynomials",
    "associated_laguerre_polynomials",
    "airy_ai",
    "airy_bi",
    "bessel_j",
    "bessel_y",
    "spherical_bessel_j",
    "spherical_bessel_y",
    "modified_bessel_i",
    "modified_bessel_k",
]
