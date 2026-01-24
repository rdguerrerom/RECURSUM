"""
Master code generation orchestrator.

Coordinates generation of all C++ code, dispatchers, bindings, tests, and notebooks.
"""

from pathlib import Path
from typing import Dict, List

from .recurrence import Recurrence
from .cpp_generator import CppGenerator
from .dispatcher_gen import DispatcherGenerator
from .binding_gen import BindingGenerator
from .test_gen import TestGenerator
from .notebook_gen import NotebookGenerator


def generate_essential(output_dir: Path = None):
    """
    Generate only essential C++ code (headers, dispatchers, bindings).

    This is used during pip install to avoid requiring numpy/scipy.
    Tests and notebooks are skipped.

    Args:
        output_dir: Root directory for output (defaults to project root)
    """
    if output_dir is None:
        output_dir = Path(__file__).parent.parent.parent  # RECURSUM root
    elif isinstance(output_dir, str):
        output_dir = Path(output_dir)

    print(f"Output directory: {output_dir}")

    # Import recurrence definitions
    from ..recurrences import (
        get_orthogonal_recurrences,
        get_bessel_recurrences,
        get_special_recurrences,
        get_quantum_recurrences,
        get_rys_recurrences,
        get_combinatorics_recurrences,
        get_mcmd_recurrences,
    )

    # Collect all recurrences by module
    all_recurrences = {
        "orthogonal": get_orthogonal_recurrences(),
        "bessel": get_bessel_recurrences(),
        "special": get_special_recurrences(),
        "quantum": get_quantum_recurrences(),
        "rys": get_rys_recurrences(),
        "combinatorics": get_combinatorics_recurrences(),
        "mcmd": get_mcmd_recurrences(),
    }

    # 1. Generate C++ recurrence headers
    print("\n" + "=" * 70)
    print("Generating C++ recurrence headers...")
    print("=" * 70)
    generate_recurrence_headers(all_recurrences, output_dir)

    # 2. Generate dispatchers
    print("\n" + "=" * 70)
    print("Generating runtime dispatchers...")
    print("=" * 70)
    generate_dispatchers(all_recurrences, output_dir)

    # 3. Generate pybind11 bindings
    print("\n" + "=" * 70)
    print("Generating pybind11 bindings...")
    print("=" * 70)
    generate_bindings(all_recurrences, output_dir)

    print("\n" + "=" * 70)
    print("✓ Essential code generation complete!")
    print("=" * 70)


def generate_all(output_dir: Path = None):
    """
    Generate all C++ code, bindings, tests, and notebooks.

    Args:
        output_dir: Root directory for output (defaults to project root)
    """
    if output_dir is None:
        output_dir = Path(__file__).parent.parent.parent  # RECURSUM root
    elif isinstance(output_dir, str):
        output_dir = Path(output_dir)

    # First generate essential C++ code
    generate_essential(output_dir)

    # Import recurrence definitions again
    from ..recurrences import (
        get_orthogonal_recurrences,
        get_bessel_recurrences,
        get_special_recurrences,
        get_quantum_recurrences,
        get_rys_recurrences,
        get_combinatorics_recurrences,
        get_mcmd_recurrences,
    )

    # Collect all recurrences by module
    all_recurrences = {
        "orthogonal": get_orthogonal_recurrences(),
        "bessel": get_bessel_recurrences(),
        "special": get_special_recurrences(),
        "quantum": get_quantum_recurrences(),
        "rys": get_rys_recurrences(),
        "combinatorics": get_combinatorics_recurrences(),
        "mcmd": get_mcmd_recurrences(),
    }

    # Import numpy-dependent generators (lazy import)
    try:
        from .test_gen import TestGenerator
        from .notebook_gen import NotebookGenerator

        # 4. Generate pytest tests
        print("\n" + "=" * 70)
        print("Generating pytest tests...")
        print("=" * 70)
        generate_tests(all_recurrences, output_dir)

        # 5. Generate Jupyter notebooks
        print("\n" + "=" * 70)
        print("Generating validation notebooks...")
        print("=" * 70)
        generate_notebooks(all_recurrences, output_dir)

        print("\n" + "=" * 70)
        print("✓ Complete code generation finished!")
        print("=" * 70)
    except ImportError as e:
        print("\n" + "=" * 70)
        print(f"⚠️  Skipping tests/notebooks (numpy/scipy not available): {e}")
        print("✓ Essential code generation complete!")
        print("=" * 70)


def generate_recurrence_headers(recurrences_by_module: Dict[str, List[Recurrence]],
                                 output_dir: Path):
    """Generate C++ template headers."""
    header_dir = output_dir / "src" / "generated"
    header_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for module_name, recurrences in recurrences_by_module.items():
        for rec in recurrences:
            filename = f"{rec.name.lower()}_coeff.hpp"
            filepath = header_dir / filename

            generator = CppGenerator(rec)
            code = generator.generate()

            with open(filepath, 'w') as f:
                f.write(code)

            print(f"  ✓ {filename}")
            count += 1

    print(f"Generated {count} header files")


def generate_dispatchers(recurrences_by_module: Dict[str, List[Recurrence]],
                         output_dir: Path):
    """Generate runtime dispatchers."""
    dispatcher_dir = output_dir / "src" / "generated" / "dispatchers"
    dispatcher_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for module_name, recurrences in recurrences_by_module.items():
        for rec in recurrences:
            filename = f"{rec.name.lower()}_dispatcher.hpp"
            filepath = dispatcher_dir / filename

            generator = DispatcherGenerator(rec)
            code = generator.generate_header()

            with open(filepath, 'w') as f:
                f.write(code)

            print(f"  ✓ {filename}")
            count += 1

    print(f"Generated {count} dispatcher files")


def generate_bindings(recurrences_by_module: Dict[str, List[Recurrence]],
                      output_dir: Path):
    """Generate pybind11 bindings - single combined module."""
    bindings_dir = output_dir / "src" / "bindings"
    bindings_dir.mkdir(parents=True, exist_ok=True)

    # Flatten all recurrences into a single list
    all_recurrences = []
    for recurrences in recurrences_by_module.values():
        all_recurrences.extend(recurrences)

    # Generate single combined binding file
    filename = "recursum_bindings.cpp"
    filepath = bindings_dir / filename

    generator = BindingGenerator(all_recurrences, "recursum")
    code = generator.generate()

    with open(filepath, 'w') as f:
        f.write(code)

    print(f"  ✓ {filename} (combined module with {len(all_recurrences)} recurrences)")
    print(f"Generated 1 combined binding file")


def generate_tests(recurrences_by_module: Dict[str, List[Recurrence]],
                   output_dir: Path):
    """Generate pytest tests."""
    test_dir = output_dir / "tests" / "generated"
    test_dir.mkdir(parents=True, exist_ok=True)

    # Create __init__.py to make it a package
    (test_dir / "__init__.py").touch()

    count = 0
    for module_name, recurrences in recurrences_by_module.items():
        for rec in recurrences:
            filename = f"test_{rec.name.lower()}.py"
            filepath = test_dir / filename

            generator = TestGenerator(rec)
            code = generator.generate()

            with open(filepath, 'w') as f:
                f.write(code)

            print(f"  ✓ {filename}")
            count += 1

    print(f"Generated {count} test files")


def generate_notebooks(recurrences_by_module: Dict[str, List[Recurrence]],
                        output_dir: Path):
    """Generate Jupyter notebooks."""
    notebook_dir = output_dir / "notebooks" / "generated"
    notebook_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for module_name, recurrences in recurrences_by_module.items():
        for rec in recurrences:
            filename = f"validate_{rec.name.lower()}.ipynb"
            filepath = notebook_dir / filename

            generator = NotebookGenerator(rec)
            notebook_json = generator.generate()

            with open(filepath, 'w') as f:
                f.write(notebook_json)

            print(f"  ✓ {filename}")
            count += 1

    print(f"Generated {count} notebook files")


if __name__ == "__main__":
    # Allow running as script
    import sys
    output_path = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    generate_all(output_path)
