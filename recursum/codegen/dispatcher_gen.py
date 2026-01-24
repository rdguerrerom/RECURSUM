"""
Dispatcher generator for runtime-to-compile-time index mapping.

Generates C++ switch statements that map runtime integer indices to
compile-time template parameter instantiations.
"""

from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from .recurrence import Recurrence


class DispatcherGenerator:
    """Generate runtime dispatchers for recurrence template instantiations."""

    def __init__(self, rec: "Recurrence"):
        """
        Initialize dispatcher generator.

        Args:
            rec: Recurrence definition to generate dispatcher for
        """
        self.rec = rec
        self.ctx = rec._ctx()

    def generate_header(self) -> str:
        """
        Generate complete dispatcher header file.

        Returns:
            String containing C++ dispatcher code with switch statements
        """
        num_indices = len(self.rec.indices)

        if num_indices == 1:
            return self._generate_1d_dispatcher()
        elif num_indices == 2:
            return self._generate_2d_dispatcher()
        elif num_indices == 3:
            return self._generate_3d_dispatcher()
        else:
            return self._generate_nd_dispatcher()

    def _generate_1d_dispatcher(self) -> str:
        """Generate dispatcher for single-index recurrence."""
        idx_name = self.rec.indices[0]
        max_idx = self.rec.max_indices[idx_name]
        func_name = f"dispatch_{self.rec.name}"

        # Generate switch cases
        cases = []
        for i in range(max_idx + 1):
            runtime_args = ", ".join(self.rec.runtime_vars)
            cases.append(
                f"        case {i}: return {self.ctx.struct_name}<{i}>::compute({runtime_args});"
            )

        # Format parameter list
        params = self._format_params()
        # Add comma before params only if params exist
        param_clause = f",\n    {params}" if params else ""

        ns_open = f"namespace {self.rec.namespace} {{\n" if self.rec.namespace else ""
        ns_close = f"\n}} // namespace {self.rec.namespace}" if self.rec.namespace else ""

        return f"""#pragma once
#include "{self.rec.name.lower()}_coeff.hpp"

{ns_open}
inline {self.rec.vec_type} {func_name}(
    int {idx_name}{param_clause}
) {{
    if ({idx_name} < 0 || {idx_name} > {max_idx}) {{
        return {self.rec.vec_type}(0.0);
    }}

    switch({idx_name}) {{
{chr(10).join(cases)}
        default: return {self.rec.vec_type}(0.0);
    }}
}}{ns_close}
"""

    def _generate_2d_dispatcher(self) -> str:
        """Generate dispatcher for two-index recurrence using nested switches."""
        idx1, idx2 = self.rec.indices[0], self.rec.indices[1]
        max_idx1 = self.rec.max_indices[idx1]
        max_idx2 = self.rec.max_indices[idx2]
        func_name = f"dispatch_{self.rec.name}"

        # Generate nested switch structure
        outer_cases = []
        for i in range(max_idx1 + 1):
            inner_cases = []
            for j in range(max_idx2 + 1):
                runtime_args = ", ".join(self.rec.runtime_vars)
                inner_cases.append(
                    f"            case {j}: return {self.ctx.struct_name}<{i}, {j}>::compute({runtime_args});"
                )

            inner_switch = f"""        case {i}:
            switch({idx2}) {{
{chr(10).join(inner_cases)}
                default: return {self.rec.vec_type}(0.0);
            }}"""
            outer_cases.append(inner_switch)

        params = self._format_params()
        # Add comma before params only if params exist
        param_clause = f",\n    {params}" if params else ""

        ns_open = f"namespace {self.rec.namespace} {{\n" if self.rec.namespace else ""
        ns_close = f"\n}} // namespace {self.rec.namespace}" if self.rec.namespace else ""

        return f"""#pragma once
#include "{self.rec.name.lower()}_coeff.hpp"

{ns_open}
inline {self.rec.vec_type} {func_name}(
    int {idx1},
    int {idx2}{param_clause}
) {{
    if ({idx1} < 0 || {idx1} > {max_idx1} ||
        {idx2} < 0 || {idx2} > {max_idx2}) {{
        return {self.rec.vec_type}(0.0);
    }}

    switch({idx1}) {{
{chr(10).join(outer_cases)}
        default: return {self.rec.vec_type}(0.0);
    }}
}}{ns_close}
"""

    def _generate_3d_dispatcher(self) -> str:
        """Generate dispatcher for three-index recurrence using triple nested switches."""
        idx1, idx2, idx3 = self.rec.indices
        max_idx1 = self.rec.max_indices[idx1]
        max_idx2 = self.rec.max_indices[idx2]
        max_idx3 = self.rec.max_indices[idx3]
        func_name = f"dispatch_{self.rec.name}"

        # Generate triple nested switch - this can be large but is still manageable
        # for typical quantum chemistry cases (e.g., 10x10x10 = 1000 cases)
        outer_cases = []
        for i in range(max_idx1 + 1):
            middle_cases = []
            for j in range(max_idx2 + 1):
                inner_cases = []
                for k in range(max_idx3 + 1):
                    runtime_args = ", ".join(self.rec.runtime_vars)
                    inner_cases.append(
                        f"                case {k}: return {self.ctx.struct_name}<{i}, {j}, {k}>::compute({runtime_args});"
                    )

                inner_switch = f"""            case {j}:
                switch({idx3}) {{
{chr(10).join(inner_cases)}
                    default: return {self.rec.vec_type}(0.0);
                }}"""
                middle_cases.append(inner_switch)

            middle_switch = f"""        case {i}:
            switch({idx2}) {{
{chr(10).join(middle_cases)}
                default: return {self.rec.vec_type}(0.0);
            }}"""
            outer_cases.append(middle_switch)

        params = self._format_params()
        # Add comma before params only if params exist
        param_clause = f",\n    {params}" if params else ""

        ns_open = f"namespace {self.rec.namespace} {{\n" if self.rec.namespace else ""
        ns_close = f"\n}} // namespace {self.rec.namespace}" if self.rec.namespace else ""

        return f"""#pragma once
#include "{self.rec.name.lower()}_coeff.hpp"

{ns_open}
inline {self.rec.vec_type} {func_name}(
    int {idx1},
    int {idx2},
    int {idx3}{param_clause}
) {{
    if ({idx1} < 0 || {idx1} > {max_idx1} ||
        {idx2} < 0 || {idx2} > {max_idx2} ||
        {idx3} < 0 || {idx3} > {max_idx3}) {{
        return {self.rec.vec_type}(0.0);
    }}

    switch({idx1}) {{
{chr(10).join(outer_cases)}
        default: return {self.rec.vec_type}(0.0);
    }}
}}{ns_close}
"""

    def _generate_nd_dispatcher(self) -> str:
        """Generate dispatcher for N-dimensional recurrence (N >= 4) using recursive dispatch."""
        func_name = f"dispatch_{self.rec.name}"
        params = self._format_params()
        ns_open = f"namespace {self.rec.namespace} {{\n" if self.rec.namespace else ""
        ns_close = f"\n}} // namespace {self.rec.namespace}" if self.rec.namespace else ""
        idx_params = ", ".join(f"int {idx}" for idx in self.rec.indices)

        # For 4+ dimensions, generate nested switches recursively
        # This handles cases like Rys VRR (4 indices)
        def generate_nested_switch(depth: int, current_indices: List[int]) -> List[str]:
            """Recursively generate nested switch statements."""
            if depth == len(self.rec.indices):
                # Base case: generate the actual function call
                template_args = ", ".join(str(i) for i in current_indices)
                runtime_args = ", ".join(self.rec.runtime_vars)
                indent = "    " * (depth + 2)
                return [f"{indent}return {self.ctx.struct_name}<{template_args}>::compute({runtime_args});"]

            idx_name = self.rec.indices[depth]
            max_idx = self.rec.max_indices[idx_name]
            indent = "    " * (depth + 1)

            cases = []
            for i in range(max_idx + 1):
                nested = generate_nested_switch(depth + 1, current_indices + [i])
                if depth + 1 == len(self.rec.indices):
                    # Last level - direct return
                    cases.append(f"{indent}case {i}: {nested[0].strip()}")
                else:
                    # Intermediate level - nested switch
                    next_idx = self.rec.indices[depth + 1]
                    cases.append(f"{indent}case {i}:")
                    cases.append(f"{indent}    switch({next_idx}) {{")
                    cases.extend(nested)
                    cases.append(f"{indent}        default: return {self.rec.vec_type}(0.0);")
                    cases.append(f"{indent}    }}")

            return cases

        # Generate bounds check
        bounds_check = " ||\n        ".join(
            f"{idx} < 0 || {idx} > {self.rec.max_indices[idx]}"
            for idx in self.rec.indices
        )

        # Generate the main switch
        first_idx = self.rec.indices[0]
        switch_body = generate_nested_switch(0, [])

        # Add comma before params only if params exist
        param_clause = f",\n    {params}" if params else ""

        return f"""#pragma once
#include "{self.rec.name.lower()}_coeff.hpp"

{ns_open}
inline {self.rec.vec_type} {func_name}(
    {idx_params}{param_clause}
) {{
    if ({bounds_check}) {{
        return {self.rec.vec_type}(0.0);
    }}

    switch({first_idx}) {{
{chr(10).join(switch_body)}
        default: return {self.rec.vec_type}(0.0);
    }}
}}{ns_close}
"""

    def _format_params(self) -> str:
        """Format runtime parameter list."""
        return ",\n    ".join(f"{self.rec.vec_type} {v}" for v in self.rec.runtime_vars)
