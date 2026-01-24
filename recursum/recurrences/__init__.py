"""
Built-in recurrence relation definitions.

This module provides pre-defined recurrence relations for:
- Orthogonal polynomials (Legendre, Chebyshev, Hermite, Laguerre, Jacobi, Gegenbauer)
- Modified spherical Bessel functions (for STO integrals)
- Special functions (Bessel J/Y, Airy, Euler, Bernoulli)
- Quantum chemistry (McMurchie-Davidson, Boys function)
- Rys quadrature (for electron repulsion integrals)
- Combinatorics (binomial coefficients, Fibonacci)
"""

from .orthogonal import *
from .bessel import *
from .quantum import *
from .rys import *
from .combinatorics import *
from .mcmd import *
from .special import *

__all__ = [
    # Orthogonal polynomials
    "hermite_coefficients",
    "legendre_polynomials",
    "chebyshev_T",
    "chebyshev_U",
    "hermite_He",
    "hermite_H",
    "laguerre",
    "associated_legendre",
    # Special functions (Jacobi, Gegenbauer, Airy, Bessel, etc.)
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
    "euler_polynomials",
    "bernoulli_polynomials",
    # Bessel functions (STO-specific)
    "modified_spherical_bessel_i",
    "modified_spherical_bessel_k",
    "reduced_bessel_b",
    "reduced_bessel_a",
    # Quantum chemistry
    "sto_auxiliary_B",
    "boys_function",
    "gaunt_coefficients",
    # Rys quadrature
    "rys_2d_integral",
    "rys_horizontal_transfer",
    "rys_vrr_full",
    "rys_polynomial_recursion",
    # Combinatorics
    "binomial_coefficients",
    "fibonacci",
    # McMurchie-Davidson Hermite E coefficient
    "hermite_e_coefficient",
]


def get_all_recurrences():
    """Get all built-in recurrences organized by module."""
    return {
        "orthogonal": get_orthogonal_recurrences(),
        "bessel": get_bessel_recurrences(),
        "special": get_special_recurrences(),
        "quantum": get_quantum_recurrences(),
        "rys": get_rys_recurrences(),
        "combinatorics": get_combinatorics_recurrences(),
        "mcmd": get_mcmd_recurrences(),
    }


def get_orthogonal_recurrences():
    """Get all orthogonal polynomial recurrences."""
    from .orthogonal import (
        hermite_coefficients, legendre_polynomials, chebyshev_T,
        chebyshev_U, hermite_He, hermite_H, laguerre, associated_legendre
    )
    return [
        hermite_coefficients(),
        legendre_polynomials(),
        chebyshev_T(),
        chebyshev_U(),
        hermite_He(),
        hermite_H(),
        laguerre(),
        associated_legendre(),
    ]


def get_bessel_recurrences():
    """Get all Bessel function recurrences."""
    from .bessel import (
        modified_spherical_bessel_i, modified_spherical_bessel_k,
        reduced_bessel_b, reduced_bessel_a
    )
    return [
        modified_spherical_bessel_i(),
        modified_spherical_bessel_k(),
        reduced_bessel_b(),
        reduced_bessel_a(),
    ]


def get_quantum_recurrences():
    """Get all quantum chemistry recurrences."""
    from .quantum import sto_auxiliary_B, boys_function, gaunt_coefficients
    return [
        sto_auxiliary_B(),
        boys_function(),
        gaunt_coefficients(),
    ]


def get_rys_recurrences():
    """Get all Rys quadrature recurrences."""
    from .rys import (
        rys_2d_integral, rys_horizontal_transfer,
        rys_vrr_full, rys_polynomial_recursion
    )
    return [
        rys_2d_integral(),
        rys_horizontal_transfer(),
        rys_vrr_full(),
        rys_polynomial_recursion(),
    ]


def get_combinatorics_recurrences():
    """Get all combinatorics recurrences."""
    from .combinatorics import binomial_coefficients, fibonacci
    return [
        binomial_coefficients(),
        fibonacci(),
    ]


def get_mcmd_recurrences():
    """Get all McMurchie-Davidson Hermite recurrences."""
    from .mcmd import hermite_e_coefficient
    return [
        hermite_e_coefficient(),
    ]


def get_special_recurrences():
    """Get all special function recurrences."""
    from .special import (
        jacobi_polynomials, gegenbauer_polynomials, associated_laguerre_polynomials,
        airy_ai, airy_bi, bessel_j, bessel_y,
        spherical_bessel_j, spherical_bessel_y,
        modified_bessel_i, modified_bessel_k,
        euler_polynomials, bernoulli_polynomials
    )
    return [
        jacobi_polynomials(),
        gegenbauer_polynomials(),
        associated_laguerre_polynomials(),
        airy_ai(),
        airy_bi(),
        bessel_j(),
        bessel_y(),
        spherical_bessel_j(),
        spherical_bessel_y(),
        modified_bessel_i(),
        modified_bessel_k(),
        euler_polynomials(),
        bernoulli_polynomials(),
    ]
