"""
Jupyter notebook generator for interactive validation.

Generates notebooks with visualization and comparison against SciPy.
"""

import json
from typing import TYPE_CHECKING

from .scipy_mapping import has_scipy_reference

if TYPE_CHECKING:
    from .recurrence import Recurrence


class NotebookGenerator:
    """Generate Jupyter notebooks for recurrence validation."""

    def __init__(self, rec: "Recurrence"):
        """
        Initialize notebook generator.

        Args:
            rec: Recurrence definition to generate notebook for
        """
        self.rec = rec

    def generate(self) -> str:
        """
        Generate Jupyter notebook JSON.

        Returns:
            JSON string of notebook
        """
        cells = [
            self._title_cell(),
            self._imports_cell(),
            self._test_data_cell(),
        ]

        if has_scipy_reference(self.rec.name):
            cells.extend([
                self._recursum_computation_cell(),
                self._scipy_reference_cell(),
                self._comparison_plot_cell(),
                self._error_analysis_cell(),
            ])
        else:
            cells.append(self._basic_computation_cell())

        notebook = {
            "cells": cells,
            "metadata": {
                "kernelspec": {
                    "display_name": "Python 3",
                    "language": "python",
                    "name": "python3"
                },
                "language_info": {
                    "name": "python",
                    "version": "3.8.0"
                }
            },
            "nbformat": 4,
            "nbformat_minor": 4
        }

        return json.dumps(notebook, indent=2)

    def _create_cell(self, cell_type: str, source: list):
        """Helper to create a cell."""
        if cell_type == "markdown":
            return {
                "cell_type": "markdown",
                "metadata": {},
                "source": source
            }
        else:  # code
            return {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": source
            }

    def _title_cell(self):
        return self._create_cell("markdown", [
            f"# {self.rec.name} Validation\\n",
            "\\n",
            "Comparison of RECURSUM implementation vs SciPy reference.\\n",
            "\\n",
            f"**Recurrence**: {self.rec.name}\\n",
            f"**Namespace**: `{self.rec.namespace}`\\n",
            f"**Indices**: {', '.join(self.rec.indices)}\\n",
            f"**Runtime vars**: {', '.join(self.rec.runtime_vars)}\\n"
        ])

    def _imports_cell(self):
        module = self._infer_module_name()
        return self._create_cell("code", [
            "import numpy as np\\n",
            "import matplotlib.pyplot as plt\\n",
            "from scipy import special\\n",
            f"import recursum._{module} as recursum\\n",
            "from recursum.codegen.scipy_mapping import get_scipy_reference, compute_base_cases\\n",
            "\\n",
            "plt.style.use('seaborn-v0_8')\\n",
            "%matplotlib inline"
        ])

    def _test_data_cell(self):
        return self._create_cell("code", [
            "# Generate test data\\n",
            f"x = np.linspace({self._test_range()})\\n",
            f"base_cases = compute_base_cases('{self.rec.name}', x)\\n",
            f"print(f'Test points: {{len(x)}}')"
        ])

    def _recursum_computation_cell(self):
        return self._create_cell("code", [
            "# Compute using RECURSUM\\n",
            "orders = [0, 1, 2, 5, 10]\\n",
            "recursum_results = {}\\n",
            "\\n",
            "for n in orders:\\n",
            f"    recursum_results[n] = recursum.{self.rec.name.lower()}(n, x, **base_cases)"
        ])

    def _scipy_reference_cell(self):
        return self._create_cell("code", [
            "# Compute using SciPy\\n",
            f"scipy_ref = get_scipy_reference('{self.rec.name}')\\n",
            "scipy_results = {}\\n",
            "\\n",
            "for n in orders:\\n",
            "    scipy_results[n] = scipy_ref(n, x=x)"
        ])

    def _comparison_plot_cell(self):
        return self._create_cell("code", [
            "# Plot comparison\\n",
            "fig, axes = plt.subplots(2, 3, figsize=(15, 10))\\n",
            "axes = axes.flatten()\\n",
            "\\n",
            "for i, n in enumerate(orders + [15]):\\n",
            "    if i >= len(axes):\\n",
            "        break\\n",
            "    ax = axes[i]\\n",
            "    \\n",
            f"    rec_res = recursum.{self.rec.name.lower()}(n, x, **base_cases)\\n",
            "    sci_res = scipy_ref(n, x=x)\\n",
            "    \\n",
            "    ax.plot(x, rec_res, 'b-', label='RECURSUM', linewidth=2)\\n",
            "    ax.plot(x, sci_res, 'r--', label='SciPy', linewidth=1)\\n",
            "    ax.set_title(f'n={n}')\\n",
            "    ax.legend()\\n",
            "    ax.grid(True, alpha=0.3)\\n",
            "\\n",
            "plt.tight_layout()\\n",
            "plt.show()"
        ])

    def _error_analysis_cell(self):
        return self._create_cell("code", [
            "# Error analysis\\n",
            "orders_full = list(range(0, 20))\\n",
            "max_errors = []\\n",
            "\\n",
            "for n in orders_full:\\n",
            f"    rec_res = recursum.{self.rec.name.lower()}(n, x, **base_cases)\\n",
            "    sci_res = scipy_ref(n, x=x)\\n",
            "    rel_error = np.abs((rec_res - sci_res) / (sci_res + 1e-15))\\n",
            "    max_errors.append(np.max(rel_error))\\n",
            "\\n",
            "plt.figure(figsize=(10, 6))\\n",
            "plt.semilogy(orders_full, max_errors, 'o-')\\n",
            "plt.axhline(1e-12, color='r', linestyle='--', label='Target accuracy')\\n",
            "plt.xlabel('Order n')\\n",
            "plt.ylabel('Max Relative Error')\\n",
            "plt.title('Accuracy vs Order')\\n",
            "plt.legend()\\n",
            "plt.grid(True, alpha=0.3)\\n",
            "plt.show()\\n",
            "\\n",
            f"print(f'Maximum error: {{max(max_errors):.2e}}')"
        ])

    def _basic_computation_cell(self):
        return self._create_cell("code", [
            f"# Compute {self.rec.name} (no SciPy reference available)\\n",
            "# TODO: Add custom validation for this recurrence"
        ])

    def _infer_module_name(self) -> str:
        """Infer module name from namespace."""
        namespace_map = {
            "legendre": "orthogonal",
            "chebyshev": "orthogonal",
            "hermite_poly": "orthogonal",
            "hermite": "orthogonal",
            "laguerre": "orthogonal",
            "bessel_sto": "bessel",
            "sto": "quantum",
            "boys_func": "quantum",
            "rys_quadrature": "rys",
            "rys_poly": "rys",
            "combinatorics": "combinatorics",
            "sequences": "combinatorics",
        }
        return namespace_map.get(self.rec.namespace, "core")

    def _test_range(self) -> str:
        """Generate test range."""
        if "Bessel" in self.rec.name:
            return "0.1, 10.0, 100"
        elif any(p in self.rec.name for p in ["Legendre", "Chebyshev"]):
            return "-1.0, 1.0, 50"
        else:
            return "-5.0, 5.0, 100"
