#!/usr/bin/env python3
"""
Demonstration of LayeredCppGenerator.

This script shows how to generate layer-by-layer C++ code for recurrences
using the new LayeredCppGenerator class.

The layered approach computes ALL auxiliary index values at once (e.g., all t
values for a given nA, nB in Hermite E), enabling true compile-time CSE.

Performance advantages over hand-written layered implementations:
- Output parameters (no return-by-value copies)
- RECURSUM_FORCEINLINE everywhere
- Exact-sized buffers (no MAX-sized arrays)
- Specialized t==0 handling (no branching in inner loop)
"""

from recursum.recurrences.mcmd import hermite_e_coefficient


def main():
    """Generate and display Hermite E layered code."""

    print("="*80)
    print("LAYERED CODE GENERATOR DEMONSTRATION")
    print("="*80)
    print()

    # Get the Hermite E recurrence definition
    rec = hermite_e_coefficient()

    print("Generating layered C++ code for Hermite E coefficients...")
    print(f"Recurrence: {rec.name}")
    print(f"Indices: {rec.indices}")
    print(f"Runtime vars: {rec.runtime_vars}")
    print()

    # Generate layered code
    layered_code = rec.generate_layered()

    print("="*80)
    print("GENERATED LAYERED CODE:")
    print("="*80)
    print(layered_code)
    print("="*80)
    print()

    # Show key features
    print("KEY FEATURES OF GENERATED CODE:")
    print("="*80)
    features = [
        ("Output parameters", "void compute(Vec8d* out, ...)"),
        ("Force-inline", "RECURSUM_FORCEINLINE everywhere"),
        ("Exact-sized buffers", "Vec8d prev[nA + nB + 1]"),
        ("Previous layer computed once", "HermiteECoeffLayer<nA-1, nB>::compute(prev, ...)"),
        ("Special t==0 case", "out[0] = PA * prev[0] + Vec8d(1) * prev[1]"),
        ("General t>0 loop", "for (int t = 1; t < N_VALUES; ++t)"),
        ("Using prev array", "prev[t-1], prev[t], prev[t+1]"),
        ("No return-by-value", "No std::array returns"),
        ("API compatibility", "Single-value accessor template"),
    ]

    for feature, example in features:
        # Case-sensitive search for exact matches
        status = "✓" if example in layered_code else "✗"
        print(f"  {status} {feature}")
        if example and status == "✓":
            print(f"      Example: {example}")

    print()
    print("="*80)
    print("PERFORMANCE COMPARISON:")
    print("="*80)
    print()
    print("Hand-written layered (current benchmarks):")
    print("  Hermite E L=0: 2.272 ns (6× slower than TMP)")
    print("  Problems:")
    print("    - Runtime loops with branch overhead")
    print("    - Return-by-value (std::array)")
    print("    - MAX-sized arrays (cache pollution)")
    print("    - Missing RECURSUM_FORCEINLINE")
    print()
    print("Generated layered (expected with this code):")
    print("  Hermite E L=0: ~0.4-0.5 ns (comparable to TMP)")
    print("  Improvements:")
    print("    - Output parameters (no copies)")
    print("    - RECURSUM_FORCEINLINE (aggressive inlining)")
    print("    - Exact-sized buffers (better cache)")
    print("    - Compiler-optimized loops")
    print()
    print("="*80)
    print("USAGE:")
    print("="*80)
    print("""
    from recursum.recurrences.mcmd import hermite_e_coefficient

    # Get recurrence
    rec = hermite_e_coefficient()

    # Generate TMP code (traditional approach)
    tmp_code = rec.generate()

    # Generate layered code (new approach)
    layered_code = rec.generate_layered()

    # Save to file
    with open("hermite_e_layered.hpp", "w") as f:
        f.write(layered_code)
    """)
    print("="*80)


if __name__ == "__main__":
    main()
