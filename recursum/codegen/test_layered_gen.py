"""
Tests for LayeredCppGenerator.

Tests:
1. Basic generation for Hermite E layered code
2. Verify structure (HermiteECoeffLayer templates)
3. Verify output parameter signature (void compute(Vec8d* out, ...))
4. Verify RECURSUM_FORCEINLINE is present
5. Verify no return-by-value (no std::array returns)
6. Verify base cases are correct
7. Compare structure to hand-written version
"""

try:
    import pytest
except ImportError:
    pytest = None

from recursum.codegen import Recurrence, LayeredCppGenerator


def test_hermite_e_layered_basic():
    """Test basic Hermite E layered code generation."""
    rec = Recurrence("HermiteE", ["nA", "nB", "t"], ["PA", "PB", "p"],
                     namespace="mcmd_hermite", max_indices={"nA": 5, "nB": 5, "t": 10})
    rec.validity("nA >= 0", "nB >= 0", "t >= 0", "t <= nA + nB")
    rec.base(nA=0, nB=0, t=0, value=1.0)

    # A-side recurrence
    rec.rule("nA > 0 && nB == 0 && t == 0",
             "PA * E[nA-1, 0, t] + (t + 1) * E[nA-1, 0, t+1]",
             name="A-side t=0")
    rec.rule("nA > 0 && nB == 0 && t > 0",
             "(0.5 / p) * E[nA-1, 0, t-1] + PA * E[nA-1, 0, t] + (t + 1) * E[nA-1, 0, t+1]",
             name="A-side t>0")

    # B-side recurrence
    rec.rule("nA == 0 && nB > 0 && t == 0",
             "PB * E[0, nB-1, t] + (t + 1) * E[0, nB-1, t+1]",
             name="B-side t=0")
    rec.rule("nA == 0 && nB > 0 && t > 0",
             "(0.5 / p) * E[0, nB-1, t-1] + PB * E[0, nB-1, t] + (t + 1) * E[0, nB-1, t+1]",
             name="B-side t>0")

    # General recurrence
    rec.rule("nA > 0 && nB > 0 && t == 0",
             "PA * E[nA-1, nB, t] + (t + 1) * E[nA-1, nB, t+1]",
             name="General t=0")
    rec.rule("nA > 0 && nB > 0 && t > 0",
             "(0.5 / p) * E[nA-1, nB, t-1] + PA * E[nA-1, nB, t] + (t + 1) * E[nA-1, nB, t+1]",
             name="General t>0")

    # Generate layered code
    code = rec.generate_layered()

    # Verify structure
    assert "HermiteECoeffLayer" in code, "Should have HermiteECoeffLayer struct"
    assert "static RECURSUM_FORCEINLINE void compute(Vec8d* out," in code, \
        "Should use output parameters"
    assert "RECURSUM_FORCEINLINE" in code, "Should have RECURSUM_FORCEINLINE macro"

    # Verify no return-by-value
    assert "std::array" not in code, "Should not use std::array return-by-value"
    assert "return std::array" not in code, "Should not return arrays"

    # Verify includes
    assert "#pragma once" in code
    assert "#include <recursum/vectorclass.h>" in code

    # Verify namespace
    assert "namespace mcmd_hermite" in code

    # Print code for inspection
    print("\n" + "="*80)
    print("GENERATED HERMITE E LAYERED CODE:")
    print("="*80)
    print(code)
    print("="*80)

    return code


def test_layered_has_base_case():
    """Test that base cases are generated correctly."""
    rec = Recurrence("Test", ["nA", "nB", "t"], ["x"])
    rec.validity("nA >= 0", "nB >= 0", "t >= 0")
    rec.base(nA=0, nB=0, t=0, value=1.0)

    code = rec.generate_layered()

    # Should have base case specialization
    assert "template<>" in code, "Should have template specialization"
    assert "Vec8d(1.0)" in code or "Vec8d(1)" in code, "Should have base case value"


def test_layered_has_accessor():
    """Test that accessor template is generated for API compatibility."""
    rec = Recurrence("Test", ["nA", "nB", "t"], ["x"])
    rec.validity("nA >= 0", "nB >= 0", "t >= 0")
    rec.base(nA=0, nB=0, t=0, value=1.0)

    code = rec.generate_layered()

    # Should have accessor template for API compatibility
    assert "API compatibility" in code or "accessor" in code.lower(), \
        "Should document API compatibility"


def test_layered_no_runtime_loops_in_small_cases():
    """
    Test that small cases don't have runtime loops.

    Note: For template-parameter-dependent sizes, we may still need
    loops, but they should be force-inlined for compiler optimization.
    """
    rec = Recurrence("Test", ["n"], ["x"])
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.rule("n > 0", "x * E[n-1]")

    code = rec.generate_layered()

    # If there are loops, they should be in force-inlined functions
    if "for (" in code:
        # Find the function containing the loop
        loop_idx = code.find("for (")
        # Look backwards to find the function signature
        func_start = code.rfind("static", 0, loop_idx)
        func_section = code[func_start:loop_idx]
        assert "RECURSUM_FORCEINLINE" in func_section, \
            "Loops should be in force-inlined functions"


def test_layered_exact_sized_buffers():
    """Test that exact-sized buffers are used, not MAX-sized."""
    rec = Recurrence("Test", ["nA", "nB", "t"], ["x"])
    rec.validity("nA >= 0", "nB >= 0", "t >= 0", "t <= nA + nB")
    rec.base(nA=0, nB=0, t=0, value=1.0)

    code = rec.generate_layered()

    # Should not have MAX_ constants
    assert "MAX_L" not in code, "Should not use MAX_L constant"
    assert "MAX_T" not in code, "Should not use MAX_T constant"

    # Should have template-parameter-dependent sizes
    # Like: Vec8d prev[nA + nB]  or  Vec8d layer[N_VALUES]
    assert "Vec8d" in code


def test_generate_via_recurrence_method():
    """Test that Recurrence.generate_layered() method works."""
    rec = Recurrence("Test", ["n"], ["x"])
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.rule("n > 0", "x * E[n-1]")

    # Should work through Recurrence method
    code = rec.generate_layered()
    assert code is not None
    assert len(code) > 0
    assert "Layer" in code


def test_compare_tmp_vs_layered_structure():
    """Compare TMP and Layered generator outputs."""
    rec = Recurrence("Test", ["n"], ["x"])
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.rule("n > 0", "x * E[n-1]")

    # Generate both
    tmp_code = rec.generate()
    layered_code = rec.generate_layered()

    print("\n" + "="*80)
    print("TMP CODE:")
    print("="*80)
    print(tmp_code)
    print("\n" + "="*80)
    print("LAYERED CODE:")
    print("="*80)
    print(layered_code)
    print("="*80)

    # Both should have RECURSUM_FORCEINLINE
    assert "RECURSUM_FORCEINLINE" in tmp_code
    assert "RECURSUM_FORCEINLINE" in layered_code

    # TMP returns values, Layered uses output params
    assert "return Vec8d" in tmp_code or "return (" in tmp_code
    assert "Vec8d* out" in layered_code

    # TMP has single-value struct, Layered has Layer struct
    assert "struct TestCoeff" in tmp_code
    assert "struct TestCoeffLayer" in layered_code


if __name__ == "__main__":
    # Run tests manually for debugging
    print("Testing Hermite E Layered Generation...")
    test_hermite_e_layered_basic()

    print("\n\nTesting base case generation...")
    test_layered_has_base_case()

    print("\n\nTesting accessor generation...")
    test_layered_has_accessor()

    print("\n\nTesting exact-sized buffers...")
    test_layered_exact_sized_buffers()

    print("\n\nTesting Recurrence.generate_layered() method...")
    test_generate_via_recurrence_method()

    print("\n\nComparing TMP vs Layered...")
    test_compare_tmp_vs_layered_structure()

    print("\n\n" + "="*80)
    print("ALL TESTS PASSED!")
    print("="*80)
