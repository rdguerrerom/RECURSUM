#!/usr/bin/env python3
"""
Generate LayeredCodegen implementation for Coulomb R integrals.
"""

from recursum.recurrences.mcmd import coulomb_r_auxiliary
from recursum.codegen import LayeredCppGenerator
from pathlib import Path

def main():
    # Define recurrence
    print("Creating Coulomb R auxiliary recurrence...")
    rec = coulomb_r_auxiliary()
    print(f"  Recurrence: {rec.name}")
    print(f"  Indices: {rec.indices}")
    print(f"  Runtime vars: {rec.runtime_vars}")
    print(f"  Rules: {len(rec._rules)}")
    print(f"  Base cases: {len(rec._base_cases)}")

    # Generate code
    print("\nGenerating LayeredCodegen C++ code...")
    gen = LayeredCppGenerator(rec, unroll_loops=True)
    code = gen.generate()

    # Write to header file
    output_path = Path(__file__).parent.parent.parent / "include" / "recursum" / "mcmd" / "coulomb_r_layered_codegen.hpp"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        f.write(code)

    print(f"\n✓ Generated LayeredCodegen code: {output_path}")
    print(f"  Lines of code: {len(code.splitlines())}")
    print(f"  File size: {len(code)} bytes")

if __name__ == "__main__":
    main()
