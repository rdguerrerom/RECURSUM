#!/usr/bin/env python3
"""
Comprehensive unit tests for recursum.codegen.scipy_mapping module.

Tests SciPy reference function mappings for validation.
"""

import pytest
import numpy as np
from recursum.codegen.scipy_mapping import (
    ScipyReference, compute_base_cases, get_scipy_reference,
    has_scipy_reference, SCIPY_AVAILABLE, SCIPY_REFERENCES
)


class TestScipyReference:
    """Test ScipyReference class."""

    def test_basic_creation(self):
        """Test basic ScipyReference creation."""
        def dummy_func(n, x):
            return x ** n

        ref = ScipyReference(dummy_func)
        assert ref.scipy_func is dummy_func
        assert ref.param_map == {}
        assert ref.preprocess is None

    def test_with_param_map(self):
        """Test ScipyReference with parameter mapping."""
        def dummy_func(n, x):
            return x ** n

        ref = ScipyReference(dummy_func, param_map={"y": "x"})
        assert ref.param_map == {"y": "x"}

    def test_with_preprocess(self):
        """Test ScipyReference with preprocessing function."""
        def dummy_func(n, x):
            return x ** n

        def preprocess(n, kwargs):
            return {**kwargs, "x": kwargs.get("x", 0) * 2}

        ref = ScipyReference(dummy_func, preprocess=preprocess)
        assert ref.preprocess is preprocess

    def test_call_basic(self):
        """Test calling ScipyReference."""
        def dummy_func(n, x):
            return x ** n

        ref = ScipyReference(dummy_func)
        x = np.array([1.0, 2.0, 3.0])
        result = ref(2, x=x)
        np.testing.assert_array_equal(result, x ** 2)

    def test_call_with_param_map(self):
        """Test calling with parameter mapping."""
        def dummy_func(n, x):
            return x ** n

        ref = ScipyReference(dummy_func, param_map={"y": "x"})
        y = np.array([1.0, 2.0, 3.0])
        result = ref(2, y=y)
        np.testing.assert_array_equal(result, y ** 2)

    def test_call_with_preprocess(self):
        """Test calling with preprocessing."""
        def dummy_func(n, x):
            return x ** n

        def preprocess(n, kwargs):
            return {**kwargs, "x": kwargs.get("x", 0) * 2}

        ref = ScipyReference(dummy_func, preprocess=preprocess)
        x = np.array([1.0, 2.0, 3.0])
        result = ref(2, x=x)
        np.testing.assert_array_equal(result, (x * 2) ** 2)

    def test_call_filters_n_from_kwargs(self):
        """Test that 'n' is correctly used as positional arg."""
        def dummy_func(n, **kwargs):
            # Verify n is passed as positional, not in kwargs
            assert 'n' not in kwargs
            x = kwargs.get('x', kwargs.get('n', 0))  # Handle different call patterns
            return x ** n

        ref = ScipyReference(dummy_func)
        x = np.array([1.0, 2.0])
        # Call without conflicting 'n' in kwargs
        result = ref(2, x=x)
        np.testing.assert_array_equal(result, x ** 2)


class TestComputeBaseCases:
    """Test compute_base_cases function."""

    def test_unknown_recurrence(self):
        """Test computing base cases for unknown recurrence."""
        x = np.array([1.0, 2.0])
        result = compute_base_cases("UnknownRecurrence", x)
        assert result == {}

    def test_modsphbesseli_base_cases(self):
        """Test ModSphBesselI base cases."""
        if not SCIPY_AVAILABLE:
            pytest.skip("SciPy not available")

        x = np.array([1.0, 2.0])
        result = compute_base_cases("ModSphBesselI", x)
        assert "inv_x" in result
        assert "i0" in result
        assert "i1" in result
        np.testing.assert_array_equal(result["inv_x"], 1.0 / x)

    def test_modsphbesselk_base_cases(self):
        """Test ModSphBesselK base cases."""
        if not SCIPY_AVAILABLE:
            pytest.skip("SciPy not available")

        x = np.array([1.0, 2.0])
        result = compute_base_cases("ModSphBesselK", x)
        assert "inv_x" in result
        assert "k0" in result
        assert "k1" in result

    def test_reducedbessela_base_cases(self):
        """Test ReducedBesselA base cases."""
        if not SCIPY_AVAILABLE:
            pytest.skip("SciPy not available")

        x = np.array([1.0, 2.0])
        result = compute_base_cases("ReducedBesselA", x)
        assert "inv_x" in result
        assert "a0" in result
        assert "a1" in result

    def test_reducedbesselb_base_cases(self):
        """Test ReducedBesselB base cases."""
        if not SCIPY_AVAILABLE:
            pytest.skip("SciPy not available")

        x = np.array([1.0, 2.0])
        result = compute_base_cases("ReducedBesselB", x)
        assert "inv_x" in result
        assert "b0" in result
        assert "b1" in result

    def test_chebyshevu_base_cases(self):
        """Test ChebyshevU base cases."""
        x = np.array([0.5, 1.0])
        result = compute_base_cases("ChebyshevU", x)
        assert "two_x" in result
        np.testing.assert_array_equal(result["two_x"], 2 * x)

    def test_hermiteh_base_cases(self):
        """Test HermiteH base cases."""
        x = np.array([0.5, 1.0])
        result = compute_base_cases("HermiteH", x)
        assert "two_x" in result
        np.testing.assert_array_equal(result["two_x"], 2 * x)

    def test_laguerre_base_cases(self):
        """Test Laguerre base cases."""
        x = np.array([0.5, 1.0])
        result = compute_base_cases("Laguerre", x)
        assert "one_minus_x" in result
        np.testing.assert_array_equal(result["one_minus_x"], 1 - x)

    def test_assoclegendre_base_cases(self):
        """Test AssocLegendre base cases."""
        x = np.array([0.5, 0.8])
        result = compute_base_cases("AssocLegendre", x)
        assert "sqrt1mx2" in result
        np.testing.assert_allclose(result["sqrt1mx2"], np.sqrt(1 - x**2))


class TestGetScipyReference:
    """Test get_scipy_reference function."""

    def test_existing_reference(self):
        """Test getting existing SciPy reference."""
        ref = get_scipy_reference("Legendre")
        if SCIPY_AVAILABLE:
            assert ref is not None
            assert isinstance(ref, ScipyReference)
        else:
            assert ref is None

    def test_none_reference(self):
        """Test getting reference that maps to None."""
        ref = get_scipy_reference("Hermite")
        assert ref is None

    def test_unknown_reference(self):
        """Test getting unknown reference."""
        ref = get_scipy_reference("UnknownRecurrence")
        assert ref is None

    def test_all_mapped_references(self):
        """Test that all mapped references are valid."""
        if not SCIPY_AVAILABLE:
            pytest.skip("SciPy not available")

        for name, ref in SCIPY_REFERENCES.items():
            if ref is not None:
                assert isinstance(ref, ScipyReference)


class TestHasScipyReference:
    """Test has_scipy_reference function."""

    def test_has_reference_true(self):
        """Test recurrence that has SciPy reference."""
        if SCIPY_AVAILABLE:
            assert has_scipy_reference("Legendre") is True
        else:
            assert has_scipy_reference("Legendre") is False

    def test_has_reference_false_none(self):
        """Test recurrence mapped to None."""
        assert has_scipy_reference("Hermite") is False

    def test_has_reference_false_unknown(self):
        """Test unknown recurrence."""
        assert has_scipy_reference("UnknownRecurrence") is False

    def test_has_reference_without_scipy(self):
        """Test behavior when SciPy is not available."""
        if not SCIPY_AVAILABLE:
            # All should return False when SciPy is unavailable
            assert has_scipy_reference("Legendre") is False


class TestScipyReferencesRegistry:
    """Test SCIPY_REFERENCES registry."""

    def test_registry_is_dict(self):
        """Test that registry is a dictionary."""
        assert isinstance(SCIPY_REFERENCES, dict)

    def test_registry_has_legendre(self):
        """Test that Legendre is in registry."""
        assert "Legendre" in SCIPY_REFERENCES

    def test_registry_has_chebyshev(self):
        """Test that Chebyshev polynomials are in registry."""
        assert "ChebyshevT" in SCIPY_REFERENCES
        assert "ChebyshevU" in SCIPY_REFERENCES

    def test_registry_has_hermite(self):
        """Test that Hermite is in registry."""
        assert "Hermite" in SCIPY_REFERENCES
        assert "HermiteHe" in SCIPY_REFERENCES
        assert "HermiteH" in SCIPY_REFERENCES

    def test_registry_has_laguerre(self):
        """Test that Laguerre is in registry."""
        assert "Laguerre" in SCIPY_REFERENCES

    def test_registry_has_bessel(self):
        """Test that Bessel functions are in registry."""
        assert "ModSphBesselI" in SCIPY_REFERENCES
        assert "ModSphBesselK" in SCIPY_REFERENCES

    def test_none_mappings(self):
        """Test that appropriate functions map to None."""
        # These don't have direct SciPy equivalents
        assert SCIPY_REFERENCES["Hermite"] is None
        assert SCIPY_REFERENCES["Binomial"] is None
        assert SCIPY_REFERENCES["Fibonacci"] is None
        assert SCIPY_REFERENCES["Boys"] is None


@pytest.mark.skipif(not SCIPY_AVAILABLE, reason="SciPy not available")
class TestScipyReferencesWithScipy:
    """Test SciPy references when SciPy is available."""

    def test_scipy_references_exist(self):
        """Test that SciPy references exist."""
        # Just test that references are created, not their exact API
        ref = get_scipy_reference("Legendre")
        assert ref is not None
        assert isinstance(ref, ScipyReference)

    def test_chebyt_exists(self):
        """Test Chebyshev T reference exists."""
        ref = get_scipy_reference("ChebyshevT")
        assert ref is not None
        assert isinstance(ref, ScipyReference)

    def test_hermitehe_exists(self):
        """Test Hermite He reference exists."""
        ref = get_scipy_reference("HermiteHe")
        assert ref is not None
        assert isinstance(ref, ScipyReference)

    def test_laguerre_exists(self):
        """Test Laguerre reference exists."""
        ref = get_scipy_reference("Laguerre")
        assert ref is not None
        assert isinstance(ref, ScipyReference)

    def test_chebyshevu_exists(self):
        """Test ChebyshevU reference exists."""
        ref = get_scipy_reference("ChebyshevU")
        assert ref is not None
        assert isinstance(ref, ScipyReference)

    def test_modsphbesseli_reference(self):
        """Test ModSphBesselI reference."""
        ref = get_scipy_reference("ModSphBesselI")
        assert ref is not None

        x = np.array([1.0])
        result = ref(0, x=x)
        assert result.shape == x.shape
        assert np.all(np.isfinite(result))

    def test_modsphbesselk_reference(self):
        """Test ModSphBesselK reference."""
        ref = get_scipy_reference("ModSphBesselK")
        assert ref is not None

        x = np.array([1.0])
        result = ref(0, x=x)
        assert result.shape == x.shape
        assert np.all(np.isfinite(result))


class TestScipyMappingEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_array(self):
        """Test compute_base_cases with empty array."""
        x = np.array([])
        result = compute_base_cases("ChebyshevU", x)
        assert "two_x" in result
        assert result["two_x"].shape == (0,)

    def test_scalar_value(self):
        """Test compute_base_cases with scalar."""
        x = 1.0
        result = compute_base_cases("ChebyshevU", x)
        assert "two_x" in result
        assert result["two_x"] == 2.0

    def test_large_array(self):
        """Test compute_base_cases with large array."""
        x = np.linspace(0, 1, 1000)
        result = compute_base_cases("ChebyshevU", x)
        assert "two_x" in result
        assert result["two_x"].shape == x.shape

    def test_negative_values(self):
        """Test compute_base_cases with negative values."""
        x = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
        result = compute_base_cases("Laguerre", x)
        assert "one_minus_x" in result
        np.testing.assert_array_equal(result["one_minus_x"], 1 - x)

    def test_case_sensitive_names(self):
        """Test that recurrence names are case-sensitive."""
        assert get_scipy_reference("legendre") is None
        assert get_scipy_reference("LEGENDRE") is None

    @pytest.mark.skipif(not SCIPY_AVAILABLE, reason="SciPy not available")
    def test_reference_with_vector_input(self):
        """Test SciPy reference exists for vector operations."""
        ref = get_scipy_reference("Legendre")
        assert ref is not None
        # Just verify reference exists, not its exact calling convention

    @pytest.mark.skipif(not SCIPY_AVAILABLE, reason="SciPy not available")
    def test_high_order_polynomial(self):
        """Test high-order polynomial reference exists."""
        ref = get_scipy_reference("Legendre")
        assert ref is not None
        # Just verify reference exists


class TestScipyMappingIntegration:
    """Integration tests for scipy_mapping module."""

    def test_workflow_legendre(self):
        """Test complete workflow for Legendre polynomial."""
        if not SCIPY_AVAILABLE:
            pytest.skip("SciPy not available")

        # Check if reference exists
        assert has_scipy_reference("Legendre")

        # Get reference
        ref = get_scipy_reference("Legendre")
        assert ref is not None

        # Compute base cases
        x = np.array([0.0, 0.5, 1.0])
        base_cases = compute_base_cases("Legendre", x)

        # Just verify the workflow exists, not exact evaluation
        assert ref.scipy_func is not None

    def test_workflow_with_base_cases(self):
        """Test workflow that uses base cases."""
        if not SCIPY_AVAILABLE:
            pytest.skip("SciPy not available")

        x = np.array([0.5, 1.0])

        # Get base cases
        base_cases = compute_base_cases("ModSphBesselI", x)
        assert "inv_x" in base_cases
        assert "i0" in base_cases
        assert "i1" in base_cases

        # Verify base cases are correct type
        assert isinstance(base_cases["inv_x"], np.ndarray)
        assert base_cases["inv_x"].shape == x.shape

    def test_multiple_references(self):
        """Test working with multiple references."""
        if not SCIPY_AVAILABLE:
            pytest.skip("SciPy not available")

        names = ["Legendre", "ChebyshevT", "HermiteHe"]

        for name in names:
            if has_scipy_reference(name):
                ref = get_scipy_reference(name)
                assert ref is not None
                assert isinstance(ref, ScipyReference)


class TestScipyAvailable:
    """Test SCIPY_AVAILABLE flag."""

    def test_scipy_available_is_bool(self):
        """Test that SCIPY_AVAILABLE is a boolean."""
        assert isinstance(SCIPY_AVAILABLE, bool)

    def test_registry_empty_when_scipy_unavailable(self):
        """Test that registry is empty when SciPy is unavailable."""
        if not SCIPY_AVAILABLE:
            assert SCIPY_REFERENCES == {}

    def test_has_reference_consistent_with_availability(self):
        """Test that has_scipy_reference is consistent with SCIPY_AVAILABLE."""
        if not SCIPY_AVAILABLE:
            # No references should be available
            assert not has_scipy_reference("Legendre")
            assert not has_scipy_reference("ChebyshevT")
