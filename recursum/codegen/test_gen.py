"""
pytest test generator for validation against SciPy.

Generates comprehensive test suites comparing RECURSUM implementations
to SciPy reference functions.
"""

from typing import TYPE_CHECKING

from .scipy_mapping import has_scipy_reference

if TYPE_CHECKING:
    from .recurrence import Recurrence


class TestGenerator:
    """Generate pytest tests for recurrence validation."""

    def __init__(self, rec: "Recurrence"):
        """
        Initialize test generator.

        Args:
            rec: Recurrence definition to generate tests for
        """
        self.rec = rec

    def generate(self) -> str:
        """
        Generate pytest test file.

        Returns:
            String containing complete pytest test module
        """
        if has_scipy_reference(self.rec.name):
            return self._generate_scipy_comparison_test()
        else:
            return self._generate_consistency_test()

    def _generate_scipy_comparison_test(self) -> str:
        """Generate test comparing to SciPy reference."""
        module_name = self._infer_module_name()

        return f'''"""Auto-generated tests for {self.rec.name}."""
import pytest
import numpy as np
from numpy.testing import assert_allclose

try:
    import recursum._{module_name} as recursum_module
    RECURSUM_AVAILABLE = True
except ImportError:
    RECURSUM_AVAILABLE = False
    recursum_module = None

from recursum.codegen.scipy_mapping import (
    get_scipy_reference,
    compute_base_cases,
    SCIPY_AVAILABLE
)


@pytest.fixture
def test_points():
    """Generate test points."""
    return np.linspace({self._test_range()})


@pytest.mark.skipif(not RECURSUM_AVAILABLE or not SCIPY_AVAILABLE,
                   reason="RECURSUM C++ extension or SciPy not available")
def test_{self.rec.name.lower()}_vs_scipy(test_points):
    """Compare {self.rec.name} against SciPy reference."""
    scipy_ref = get_scipy_reference("{self.rec.name}")
    assert scipy_ref is not None, "SciPy reference not found"

    # Test multiple orders
    orders = {self._test_orders()}

    for n in orders:
        # Compute base cases
        base_cases = compute_base_cases("{self.rec.name}", test_points)

        # RECURSUM result (using positional arguments for pybind11 compatibility)
        params = [base_cases[k] if k in base_cases else test_points
                  for k in {self.rec.runtime_vars}]
        recursum_result = recursum_module.{self.rec.name.lower()}(n, *params)

        # SciPy reference
        scipy_result = scipy_ref(n, x=test_points)

        # Compare
        assert_allclose(
            recursum_result, scipy_result,
            rtol={self._relative_tolerance()},
            atol={self._absolute_tolerance()},
            err_msg=f"Mismatch at n={{n}}"
        )


@pytest.mark.skipif(not RECURSUM_AVAILABLE,
                   reason="RECURSUM C++ extension not available")
def test_{self.rec.name.lower()}_vectorization(test_points):
    """Test SIMD vectorization consistency."""
    if not RECURSUM_AVAILABLE:
        pytest.skip("RECURSUM not available")

    # Test with array sizes that are/aren't multiples of 8
    sizes = [7, 8, 15, 16, 32, 100]

    base_cases = compute_base_cases("{self.rec.name}", test_points)

    for size in sizes:
        x = test_points[:size]
        base_subset = {{k: v[:size] if isinstance(v, np.ndarray) and len(v) > 1 else v
                       for k, v in base_cases.items()}}

        params_vec = [base_subset[k] if k in base_subset else x
                      for k in {self.rec.runtime_vars}]
        result = recursum_module.{self.rec.name.lower()}(5, *params_vec)
        assert len(result) == size, f"Expected size {{size}}, got {{len(result)}}"


@pytest.mark.skipif(not RECURSUM_AVAILABLE,
                   reason="RECURSUM C++ extension not available")
def test_{self.rec.name.lower()}_boundary_cases():
    """Test boundary conditions."""
    x = np.array([1.0, 2.0, 3.0])
    base_cases = compute_base_cases("{self.rec.name}", x)

    # Base case (n=0)
    params_bound = [base_cases[k] if k in base_cases else x
                    for k in {self.rec.runtime_vars}]
    result_0 = recursum_module.{self.rec.name.lower()}(0, *params_bound)
    assert result_0.shape == x.shape

    # Maximum index
    max_n = {self.rec.max_indices[self.rec.indices[0]]}
    result_max = recursum_module.{self.rec.name.lower()}(max_n, *params_bound)
    assert result_max.shape == x.shape

    # Out of range should return zeros
    result_oor = recursum_module.{self.rec.name.lower()}(max_n + 1, *params_bound)
    assert_allclose(result_oor, 0.0)


@pytest.mark.skipif(not RECURSUM_AVAILABLE or not SCIPY_AVAILABLE,
                   reason="RECURSUM C++ extension or SciPy not available")
@pytest.mark.parametrize("n", {self._test_orders()})
def test_{self.rec.name.lower()}_accuracy(n, test_points):
    """Parameterized accuracy test."""
    scipy_ref = get_scipy_reference("{self.rec.name}")
    if scipy_ref is None:
        pytest.skip("No SciPy reference available")

    base_cases = compute_base_cases("{self.rec.name}", test_points)

    params_acc = [base_cases[k] if k in base_cases else test_points
                  for k in {self.rec.runtime_vars}]
    recursum_result = recursum_module.{self.rec.name.lower()}(n, *params_acc)
    scipy_result = scipy_ref(n, x=test_points)

    assert_allclose(recursum_result, scipy_result, rtol=1e-12, atol=1e-14)
'''

    def _generate_consistency_test(self) -> str:
        """Generate consistency test (no SciPy reference)."""
        module_name = self._infer_module_name()

        return f'''"""Auto-generated consistency tests for {self.rec.name}."""
import pytest
import numpy as np
from numpy.testing import assert_allclose

try:
    import recursum._{module_name} as recursum_module
    RECURSUM_AVAILABLE = True
except ImportError:
    RECURSUM_AVAILABLE = False
    recursum_module = None


@pytest.mark.skipif(not RECURSUM_AVAILABLE,
                   reason="RECURSUM C++ extension not available")
def test_{self.rec.name.lower()}_consistency():
    """Test internal consistency of {self.rec.name}."""
    # Basic smoke test - ensure function is callable
    # TODO: Add specific consistency checks for {self.rec.name}
    assert RECURSUM_AVAILABLE, "Extension not loaded"


@pytest.mark.skipif(not RECURSUM_AVAILABLE,
                   reason="RECURSUM C++ extension not available")
def test_{self.rec.name.lower()}_base_cases():
    """Verify base cases."""
    # TODO: Implement base case verification for {self.rec.name}
    pass
'''

    def _infer_module_name(self) -> str:
        """Infer which Python module this recurrence belongs to."""
        # Map namespaces to module names
        namespace_map = {
            "legendre": "orthogonal",
            "chebyshev": "orthogonal",
            "hermite_poly": "orthogonal",
            "hermite": "orthogonal",
            "laguerre": "orthogonal",
            "bessel_sto": "bessel",
            "sto": "quantum",
            "boys_func": "quantum",
            "angular": "quantum",
            "rys_quadrature": "rys",
            "rys_poly": "rys",
            "combinatorics": "combinatorics",
            "sequences": "combinatorics",
            # New special functions (all go to _special module)
            "jacobi": "special",
            "gegenbauer": "special",
            "assoc_laguerre": "special",
            "airy": "special",
            "bessel": "special",
            "euler": "special",
            "bernoulli": "special",
            "number_theory": "special",
        }
        return namespace_map.get(self.rec.namespace, "core")

    def _test_range(self) -> str:
        """Generate appropriate test range for this recurrence."""
        if "Bessel" in self.rec.name:
            return "0.1, 10.0, 100"  # Avoid x=0
        elif any(p in self.rec.name for p in ["Legendre", "Chebyshev"]):
            return "-1.0, 1.0, 50"  # Standard domain
        else:
            return "-5.0, 5.0, 100"

    def _test_orders(self) -> str:
        """Generate test orders."""
        max_n = self.rec.max_indices.get(self.rec.indices[0], 10)
        return f"[0, 1, 2, 5, {max_n//2}, {min(max_n, 15)}]"

    def _relative_tolerance(self) -> str:
        """Get relative tolerance for comparison."""
        return "1e-12"

    def _absolute_tolerance(self) -> str:
        """Get absolute tolerance for comparison."""
        return "1e-14"
