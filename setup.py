"""
Setup script for RECURSUM with CMake integration.

This script handles:
1. Auto-generation of C++ code from recurrence definitions
2. CMake build of pybind11 extensions
3. Installation of Python package
"""

import os
import subprocess
import sys
from pathlib import Path
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CMakeBuild(build_ext):
    """Custom build extension that runs CMake."""

    def run(self):
        """Run CMake build process."""
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build RECURSUM extensions. "
                "Install with: pip install cmake"
            )

        # Step 1: Initialize vectorclass dependency
        self.init_vectorclass()

        # Step 2: Generate C++ code from recurrence definitions
        self.generate_code()

        # Step 3: Run CMake for each extension
        for ext in self.extensions:
            self.build_cmake(ext)

        super().run()

    def init_vectorclass(self):
        """Initialize vectorclass library (download if needed)."""
        vectorclass_dir = Path(__file__).parent / "external" / "vectorclass"

        if vectorclass_dir.exists() and (vectorclass_dir / "vectorclass.h").exists():
            print("✓ Vectorclass already initialized")
            return

        print("=" * 70)
        print("Initializing vectorclass library...")
        print("=" * 70)

        init_script = Path(__file__).parent / "init_vectorclass.sh"
        if init_script.exists():
            try:
                subprocess.run(
                    ["bash", str(init_script)],
                    cwd=Path(__file__).parent,
                    check=True,
                    capture_output=True,
                    text=True
                )
                print("✓ Vectorclass initialized successfully")
            except subprocess.CalledProcessError as e:
                print(f"Error initializing vectorclass: {e.stderr}")
                raise RuntimeError(
                    "Failed to initialize vectorclass library. "
                    "Please run ./init_vectorclass.sh manually."
                )
        else:
            raise RuntimeError(
                "init_vectorclass.sh not found. "
                "Please ensure the script exists in the project root."
            )

    def generate_code(self):
        """Generate C++ headers, dispatchers, and bindings."""
        print("=" * 70)
        print("Generating C++ code from recurrence definitions...")
        print("=" * 70)

        try:
            # Import the generate_essential function (doesn't require numpy/scipy)
            sys.path.insert(0, str(Path(__file__).parent))
            from recursum.codegen.orchestrator import generate_essential

            # Generate essential code (headers, dispatchers, bindings)
            # Tests and notebooks are skipped during pip install
            output_dir = Path(__file__).parent
            generate_essential(output_dir)

            print("✓ Code generation complete")
        except ImportError as e:
            print(f"Warning: Could not generate code: {e}")
            print("Continuing with existing generated files...")
        except Exception as e:
            print(f"Error during code generation: {e}")
            raise

    def build_cmake(self, ext):
        """Run CMake build."""
        ext_dir = Path(self.get_ext_fullpath(ext.name)).parent.absolute()
        build_temp = Path(self.build_temp).absolute()
        build_temp.mkdir(parents=True, exist_ok=True)

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={ext_dir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE=Release",
            "-DRECURSUM_USE_NATIVE_ARCH=ON",
        ]

        build_args = ["--config", "Release"]

        # Detect number of CPUs for parallel build
        # Limit to 8 jobs to avoid OOM during template-heavy C++ compilation
        # Can be overridden with RECURSUM_BUILD_JOBS environment variable
        try:
            import multiprocessing
            max_jobs = int(os.environ.get("RECURSUM_BUILD_JOBS", "8"))
            num_jobs = min(max_jobs, multiprocessing.cpu_count())
            build_args.extend(["--", f"-j{num_jobs}"])
        except:
            build_args.extend(["--", "-j4"])

        print("=" * 70)
        print("Running CMake configuration...")
        print("=" * 70)

        subprocess.run(
            ["cmake", str(Path(__file__).parent)] + cmake_args,
            cwd=build_temp,
            check=True
        )

        print("=" * 70)
        print("Building C++ extensions...")
        print("=" * 70)

        subprocess.run(
            ["cmake", "--build", "."] + build_args,
            cwd=build_temp,
            check=True
        )

        print("✓ Build complete")


# Dummy extension to trigger CMake build
ext_modules = [
    Extension(
        "recursum._core",
        sources=[],
    )
]

if __name__ == "__main__":
    setup(
        ext_modules=ext_modules,
        cmdclass={
            "build_ext": CMakeBuild,
        },
        zip_safe=False,
    )
