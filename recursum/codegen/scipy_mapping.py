"""
SciPy reference function mappings for validation.

Maps RECURSUM recurrence names to corresponding SciPy special functions
for accuracy testing and validation.
"""

import numpy as np
from typing import Callable, Optional, Dict, Any

try:
    from scipy import special
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    special = None

try:
    from sympy.physics.wigner import gaunt as sympy_gaunt
    from sympy import bernoulli as sympy_bernoulli
    from sympy import euler as sympy_euler
    SYMPY_AVAILABLE = True
except ImportError:
    SYMPY_AVAILABLE = False
    sympy_gaunt = None
    sympy_bernoulli = None
    sympy_euler = None


class ScipyReference:
    """Wrapper for SciPy reference implementation."""

    def __init__(self, scipy_func: Callable, param_map: Optional[Dict[str, str]] = None,
                 preprocess: Optional[Callable] = None):
        """
        Initialize SciPy reference.

        Args:
            scipy_func: SciPy special function
            param_map: Mapping from RECURSUM param names to SciPy names
            preprocess: Optional function to preprocess input arguments
        """
        self.scipy_func = scipy_func
        self.param_map = param_map or {}
        self.preprocess = preprocess

    def __call__(self, n: int, **kwargs) -> np.ndarray:
        """
        Call SciPy function with mapped parameters.

        Args:
            n: Order/index
            **kwargs: Runtime parameters

        Returns:
            Result array from SciPy
        """
        if self.preprocess:
            kwargs = self.preprocess(n, kwargs)

        # Map parameter names
        scipy_kwargs = {
            self.param_map.get(k, k): v
            for k, v in kwargs.items()
            if k != 'n'
        }

        # Extract x and pass as positional argument (SciPy functions expect positional args)
        x = scipy_kwargs.pop('x', None)
        if x is not None:
            return self.scipy_func(n, x, **scipy_kwargs)
        else:
            return self.scipy_func(n, **scipy_kwargs)


def compute_base_cases(recurrence_name: str, x: np.ndarray) -> Dict[str, np.ndarray]:
    """
    Compute base cases needed for recurrence evaluation.

    Args:
        recurrence_name: Name of the recurrence
        x: Input points

    Returns:
        Dictionary of base case values
    """
    if not SCIPY_AVAILABLE:
        return {}

    if recurrence_name == "ModSphBesselI":
        return {
            "inv_x": 1.0 / x,
            "i0": np.sinh(x) / x,
            "i1": np.cosh(x) / x - np.sinh(x) / (x**2),
        }
    elif recurrence_name == "ModSphBesselK":
        return {
            "inv_x": 1.0 / x,
            "k0": (np.pi / 2) * np.exp(-x) / x,
            "k1": (np.pi / 2) * np.exp(-x) / x * (1 + 1/x),
        }
    elif recurrence_name == "ReducedBesselB":
        return {
            "inv_x": 1.0 / x,
            "b0": (1 - np.exp(-2*x)) / (2*x),
            "b1": (1 + np.exp(-2*x)) / (2*x) - (1 - np.exp(-2*x)) / (2*x**2),
        }
    elif recurrence_name == "ReducedBesselA":
        return {
            "inv_x": 1.0 / x,
            "a0": (np.pi / 2) / x,
            "a1": (np.pi / 2) / x * (1 + 1/x),
        }
    elif recurrence_name == "ChebyshevU":
        return {
            "two_x": 2 * x,
        }
    elif recurrence_name == "HermiteH":
        return {
            "two_x": 2 * x,
        }
    elif recurrence_name == "Laguerre":
        return {
            "one_minus_x": 1 - x,
        }
    elif recurrence_name == "AssocLegendre":
        # For associated Legendre, sqrt(1-x^2) is needed
        return {
            "sqrt1mx2": np.sqrt(1 - x**2),
        }
    # New special functions base cases
    elif recurrence_name == "BesselJ":
        return {
            "inv_x": 1.0 / x,
            "J0": special.j0(x) if SCIPY_AVAILABLE else np.zeros_like(x),
            "J1": special.j1(x) if SCIPY_AVAILABLE else np.zeros_like(x),
        }
    elif recurrence_name == "BesselY":
        return {
            "inv_x": 1.0 / x,
            "Y0": special.y0(x) if SCIPY_AVAILABLE else np.zeros_like(x),
            "Y1": special.y1(x) if SCIPY_AVAILABLE else np.zeros_like(x),
        }
    elif recurrence_name == "SphericalBesselJ":
        return {
            "inv_x": 1.0 / x,
            "j0": np.sin(x) / x,
            "j1": np.sin(x) / (x**2) - np.cos(x) / x,
        }
    elif recurrence_name == "SphericalBesselY":
        return {
            "inv_x": 1.0 / x,
            "y0": -np.cos(x) / x,
            "y1": -np.cos(x) / (x**2) - np.sin(x) / x,
        }
    elif recurrence_name == "ModifiedBesselI":
        return {
            "inv_x": 1.0 / x,
            "I0": special.i0(x) if SCIPY_AVAILABLE else np.zeros_like(x),
            "I1": special.i1(x) if SCIPY_AVAILABLE else np.zeros_like(x),
        }
    elif recurrence_name == "ModifiedBesselK":
        return {
            "inv_x": 1.0 / x,
            "K0": special.k0(x) if SCIPY_AVAILABLE else np.zeros_like(x),
            "K1": special.k1(x) if SCIPY_AVAILABLE else np.zeros_like(x),
        }
    elif recurrence_name == "AiryAi":
        ai, aip, _, _ = special.airy(x) if SCIPY_AVAILABLE else (np.zeros_like(x), np.zeros_like(x), None, None)
        return {
            "Ai0": ai,
            "Ai1": aip,
        }
    elif recurrence_name == "AiryBi":
        _, _, bi, bip = special.airy(x) if SCIPY_AVAILABLE else (None, None, np.zeros_like(x), np.zeros_like(x))
        return {
            "Bi0": bi,
            "Bi1": bip,
        }
    elif recurrence_name == "Jacobi":
        # Fixed parameters alpha=0, beta=0 for testing (reduces to Legendre)
        alpha, beta = 0.0, 0.0
        return {
            "x": x,
            "alpha": np.full_like(x, alpha),
            "beta": np.full_like(x, beta),
            "P0": np.ones_like(x),  # P^(alpha,beta)_0 = 1
            "P1": 0.5 * (alpha - beta + (alpha + beta + 2) * x),  # P^(alpha,beta)_1
        }
    elif recurrence_name == "Gegenbauer":
        # Fixed parameter lambda=0.5 for testing (Chebyshev second kind)
        lam = 0.5
        return {
            "x": x,
            "lambda": np.full_like(x, lam),
            "C0": np.ones_like(x),  # C^(lambda)_0 = 1
            "C1": 2 * lam * x,  # C^(lambda)_1 = 2*lambda*x
        }
    elif recurrence_name == "AssocLaguerre":
        # Fixed parameter alpha=0.0 for testing (reduces to ordinary Laguerre)
        alpha = 0.0
        return {
            "x": x,
            "alpha": np.full_like(x, alpha),
            "L0": np.ones_like(x),  # L^(alpha)_0 = 1
            "L1": 1 + alpha - x,  # L^(alpha)_1 = 1 + alpha - x
        }
    return {}


def _airy_ai_derivative(n: int, x: np.ndarray) -> np.ndarray:
    """Compute nth derivative of Airy Ai function using recurrence relation.

    Ai^(n+2)(x) = x * Ai^(n)(x) + n * Ai^(n-1)(x)
    """
    if n == 0:
        return special.airy(x)[0]  # Ai(x)
    elif n == 1:
        return special.airy(x)[1]  # Ai'(x)

    # Use recurrence: Ai^(n+2) = x * Ai^(n) + n * Ai^(n-1)
    ai_nm2 = special.airy(x)[0]  # Ai(x)
    ai_nm1 = special.airy(x)[1]  # Ai'(x)

    for k in range(2, n + 1):
        ai_n = x * ai_nm2 + (k - 1) * ai_nm1
        ai_nm2, ai_nm1 = ai_nm1, ai_n

    return ai_nm1


def _airy_bi_derivative(n: int, x: np.ndarray) -> np.ndarray:
    """Compute nth derivative of Airy Bi function using recurrence relation.

    Bi^(n+2)(x) = x * Bi^(n)(x) + n * Bi^(n-1)(x)
    """
    if n == 0:
        return special.airy(x)[2]  # Bi(x)
    elif n == 1:
        return special.airy(x)[3]  # Bi'(x)

    # Use recurrence: Bi^(n+2) = x * Bi^(n) + n * Bi^(n-1)
    bi_nm2 = special.airy(x)[2]  # Bi(x)
    bi_nm1 = special.airy(x)[3]  # Bi'(x)

    for k in range(2, n + 1):
        bi_n = x * bi_nm2 + (k - 1) * bi_nm1
        bi_nm2, bi_nm1 = bi_nm1, bi_n

    return bi_nm1


# Mapping registry
def _create_scipy_references() -> Dict[str, Optional[ScipyReference]]:
    """Create the mapping of recurrence names to SciPy functions."""
    if not SCIPY_AVAILABLE:
        return {}

    return {
        # Orthogonal polynomials
        "Legendre": ScipyReference(special.eval_legendre),
        "ChebyshevT": ScipyReference(special.eval_chebyt),
        "ChebyshevU": ScipyReference(
            special.eval_chebyu,
            preprocess=lambda n, kw: {"x": kw["x"]} if "x" in kw else {"x": kw.get("two_x", 0) / 2}
        ),
        "HermiteHe": ScipyReference(
            special.eval_hermitenorm,
            preprocess=lambda n, kw: {"x": kw["x"]} if "x" in kw else kw
        ),
        "HermiteH": ScipyReference(
            special.eval_hermite,
            preprocess=lambda n, kw: {"x": kw["x"]} if "x" in kw else {"x": kw.get("two_x", 0) / 2}
        ),
        "Laguerre": ScipyReference(
            special.eval_laguerre,
            preprocess=lambda n, kw: {"x": kw["x"]} if "x" in kw else {**kw, "x": 1 - kw.get("one_minus_x", 0)}
        ),

        # New orthogonal polynomials
        "Jacobi": ScipyReference(
            lambda n, x, alpha=0.0, beta=0.0: special.eval_jacobi(n, alpha, beta, x),
            preprocess=lambda n, kw: {"x": kw["x"], "alpha": kw.get("alpha", 0.0), "beta": kw.get("beta", 0.0)}
        ),
        "Gegenbauer": ScipyReference(
            lambda n, x, lam=0.5: special.eval_gegenbauer(n, lam, x),
            preprocess=lambda n, kw: {"x": kw["x"], "lam": kw.get("lambda", 0.5)}
        ),
        "AssocLaguerre": ScipyReference(
            lambda n, x, alpha=0.0: special.eval_genlaguerre(n, alpha, x),
            preprocess=lambda n, kw: {"x": kw["x"], "alpha": kw.get("alpha", 0.0)}
        ),

        # Modified spherical Bessel functions
        "ModSphBesselI": ScipyReference(
            lambda n, x: np.sqrt(np.pi / (2*x)) * special.iv(n + 0.5, x),
            preprocess=lambda n, kw: {"x": kw.get("x", 1.0 / kw.get("inv_x", 1.0))}
        ),
        "ModSphBesselK": ScipyReference(
            lambda n, x: np.sqrt(np.pi / (2*x)) * special.kv(n + 0.5, x),
            preprocess=lambda n, kw: {"x": kw.get("x", 1.0 / kw.get("inv_x", 1.0))}
        ),

        # Bessel functions
        "BesselJ": ScipyReference(
            special.jv,
            preprocess=lambda n, kw: {"x": kw.get("x", 1.0 / kw.get("inv_x", 1.0))}
        ),
        "BesselY": ScipyReference(
            special.yv,
            preprocess=lambda n, kw: {"x": kw.get("x", 1.0 / kw.get("inv_x", 1.0))}
        ),
        "SphericalBesselJ": ScipyReference(
            special.spherical_jn,
            preprocess=lambda n, kw: {"x": kw.get("x", 1.0 / kw.get("inv_x", 1.0))}
        ),
        "SphericalBesselY": ScipyReference(
            special.spherical_yn,
            preprocess=lambda n, kw: {"x": kw.get("x", 1.0 / kw.get("inv_x", 1.0))}
        ),
        "ModifiedBesselI": ScipyReference(
            special.iv,
            preprocess=lambda n, kw: {"x": kw.get("x", 1.0 / kw.get("inv_x", 1.0))}
        ),
        "ModifiedBesselK": ScipyReference(
            special.kv,
            preprocess=lambda n, kw: {"x": kw.get("x", 1.0 / kw.get("inv_x", 1.0))}
        ),

        # Airy functions - compute nth derivatives using recurrence relation
        # Ai^(n+2)(x) = x * Ai^(n)(x) + n * Ai^(n-1)(x)
        "AiryAi": ScipyReference(
            lambda n, x: _airy_ai_derivative(n, x),
            preprocess=lambda n, kw: {"x": kw["x"]}
        ),
        "AiryBi": ScipyReference(
            lambda n, x: _airy_bi_derivative(n, x),
            preprocess=lambda n, kw: {"x": kw["x"]}
        ),

        # No direct SciPy equivalent for these
        "Hermite": None,  # McMurchie-Davidson coefficients (no SciPy equivalent)
        "AssocLegendre": None,  # Complex; scipy.special.lpmv exists but different interface
        "Binomial": None,  # scipy.special.comb is not a recurrence
        "Fibonacci": None,  # No SciPy equivalent
        "ReducedBesselB": None,  # Scaled version, no direct SciPy equivalent
        "ReducedBesselA": None,  # Scaled version, no direct SciPy equivalent
        "STOAuxB": None,  # STO-specific, no SciPy equivalent
        "Boys": None,  # Boys function not in SciPy
        "Rys2D": None,  # Rys quadrature specific
        "RysHRR": None,
        "RysVRRFull": None,
        "RysPoly": None,
    }


SCIPY_REFERENCES = _create_scipy_references()


def get_scipy_reference(recurrence_name: str) -> Optional[ScipyReference]:
    """
    Get SciPy reference function for a recurrence.

    Args:
        recurrence_name: Name of the recurrence

    Returns:
        ScipyReference wrapper or None if no SciPy equivalent exists
    """
    return SCIPY_REFERENCES.get(recurrence_name)


def has_scipy_reference(recurrence_name: str) -> bool:
    """
    Check if a recurrence has a SciPy reference implementation.

    Args:
        recurrence_name: Name of the recurrence

    Returns:
        True if SciPy reference is available
    """
    return SCIPY_AVAILABLE and get_scipy_reference(recurrence_name) is not None
