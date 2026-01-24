"""
pybind11 binding generator for Python-C++ interface.

Generates pybind11 wrapper code that handles numpy array conversion
and calls the C++ dispatchers.
"""

from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from .recurrence import Recurrence


class BindingGenerator:
    """Generate pybind11 bindings for recurrence functions."""

    def __init__(self, recurrences: List["Recurrence"], module_name: str):
        """
        Initialize binding generator.

        Args:
            recurrences: List of recurrence definitions for this module
            module_name: Name of the Python module (e.g., "orthogonal")
        """
        self.recurrences = recurrences
        self.module_name = module_name

    def generate(self) -> str:
        """
        Generate complete pybind11 binding file.

        Returns:
            String containing C++ pybind11 code
        """
        includes = self._generate_includes()
        wrappers = self._generate_wrapper_functions()
        module_def = self._generate_module_definition()

        return f"""// Auto-generated pybind11 bindings for {self.module_name}
{includes}

{wrappers}

{module_def}
"""

    def _generate_includes(self) -> str:
        """Generate include statements."""
        # Use dispatchers/ since src/generated is in include path
        headers = [f'"dispatchers/{r.name.lower()}_dispatcher.hpp"'
                   for r in self.recurrences]

        return "\n".join([
            "#include <pybind11/pybind11.h>",
            "#include <pybind11/numpy.h>",
            "#include <pybind11/stl.h>",
            *[f"#include {h}" for h in headers],
        ])

    def _generate_wrapper_functions(self) -> str:
        """Generate numpy-compatible wrapper functions."""
        wrappers = []
        for rec in self.recurrences:
            wrappers.append(self._generate_single_wrapper(rec))
        return "\n\n".join(wrappers)

    def _generate_single_wrapper(self, rec: "Recurrence") -> str:
        """Generate wrapper for one recurrence."""
        num_indices = len(rec.indices)

        if num_indices == 1:
            return self._generate_1d_wrapper(rec)
        elif num_indices == 2:
            return self._generate_2d_wrapper(rec)
        else:
            # For higher dimensions, use simpler implementation
            return self._generate_nd_wrapper(rec)

    def _generate_1d_wrapper(self, rec: "Recurrence") -> str:
        """Generate wrapper for single-index recurrence."""
        idx_name = rec.indices[0]
        func_name = rec.name.lower()
        dispatcher = f"dispatch_{rec.name}"

        # Primary input is usually the first runtime var
        if rec.runtime_vars:
            primary_var = rec.runtime_vars[0]
            other_vars = rec.runtime_vars[1:]
        else:
            # Handle recurrences with no runtime vars (like Binomial)
            primary_var = None
            other_vars = []

        ns = f"{rec.namespace}::" if rec.namespace else ""

        if primary_var is None:
            # No runtime parameters - pure index-based
            return f"""namespace py = pybind11;

py::array_t<double> {func_name}_wrapper(int {idx_name}) {{
    py::array_t<double> result(1);
    auto buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);

    Vec8d vec_result = {ns}{dispatcher}({idx_name});
    ptr[0] = vec_result[0];  // Extract first element only

    return result;
}}
"""

        # Generate wrapper with numpy array handling
        other_params = ", ".join(f"py::array_t<double> {v}" for v in other_vars)
        param_list = f"py::array_t<double> {primary_var}" + (f", {other_params}" if other_params else "")

        return f"""namespace py = pybind11;

py::array_t<double> {func_name}_wrapper(
    int {idx_name},
    {param_list}
) {{
    auto buf = {primary_var}.request();
    if (buf.ndim != 1) {{
        throw std::runtime_error("Input must be 1D array");
    }}

    size_t arr_size = buf.shape[0];
    py::array_t<double> result(arr_size);
    auto result_buf = result.request();
    double* result_ptr = static_cast<double*>(result_buf.ptr);
    double* input_ptr = static_cast<double*>(buf.ptr);

    // Load other parameters (assume scalar or array)
{self._generate_param_loads(other_vars, indent="    ")}

    // Process in Vec8d chunks
    size_t vec_idx = 0;
    for (; vec_idx + 8 <= arr_size; vec_idx += 8) {{
        Vec8d vec_input;
        vec_input.load(input_ptr + vec_idx);

{self._generate_dispatcher_call(rec, ns, dispatcher, idx_name, other_vars, indent="        ")}

        result_vec.store(result_ptr + vec_idx);
    }}

    // Handle remaining elements
    if (vec_idx < arr_size) {{
        Vec8d vec_input;
        vec_input.load_partial(static_cast<int>(arr_size - vec_idx), input_ptr + vec_idx);

{self._generate_dispatcher_call(rec, ns, dispatcher, idx_name, other_vars, indent="        ")}

        result_vec.store_partial(static_cast<int>(arr_size - vec_idx), result_ptr + vec_idx);
    }}

    return result;
}}
"""

    def _generate_2d_wrapper(self, rec: "Recurrence") -> str:
        """Generate wrapper for two-index recurrence."""
        idx1, idx2 = rec.indices[0], rec.indices[1]
        func_name = rec.name.lower()
        dispatcher = f"dispatch_{rec.name}"

        if rec.runtime_vars:
            primary_var = rec.runtime_vars[0]
            other_vars = rec.runtime_vars[1:]
        else:
            primary_var = None
            other_vars = []

        ns = f"{rec.namespace}::" if rec.namespace else ""

        if primary_var is None:
            # No runtime parameters
            return f"""namespace py = pybind11;

py::array_t<double> {func_name}_wrapper(int {idx1}, int {idx2}) {{
    py::array_t<double> result(1);
    auto buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);

    Vec8d vec_result = {ns}{dispatcher}({idx1}, {idx2});
    ptr[0] = vec_result[0];  // Extract first element only

    return result;
}}
"""

        other_params = ", ".join(f"py::array_t<double> {v}" for v in other_vars)
        param_list = f"py::array_t<double> {primary_var}" + (f", {other_params}" if other_params else "")

        return f"""namespace py = pybind11;

py::array_t<double> {func_name}_wrapper(
    int {idx1},
    int {idx2},
    {param_list}
) {{
    auto buf = {primary_var}.request();
    if (buf.ndim != 1) {{
        throw std::runtime_error("Input must be 1D array");
    }}

    size_t arr_size = buf.shape[0];
    py::array_t<double> result(arr_size);
    auto result_buf = result.request();
    double* result_ptr = static_cast<double*>(result_buf.ptr);
    double* input_ptr = static_cast<double*>(buf.ptr);

{self._generate_param_loads(other_vars, indent="    ")}

    // Process in Vec8d chunks
    size_t vec_idx = 0;
    for (; vec_idx + 8 <= arr_size; vec_idx += 8) {{
        Vec8d vec_input;
        vec_input.load(input_ptr + vec_idx);

{self._generate_dispatcher_call(rec, ns, dispatcher, f"{idx1}, {idx2}", other_vars, indent="        ")}

        result_vec.store(result_ptr + vec_idx);
    }}

    // Handle remaining elements
    if (vec_idx < arr_size) {{
        Vec8d vec_input;
        vec_input.load_partial(static_cast<int>(arr_size - vec_idx), input_ptr + vec_idx);

{self._generate_dispatcher_call(rec, ns, dispatcher, f"{idx1}, {idx2}", other_vars, indent="        ")}

        result_vec.store_partial(static_cast<int>(arr_size - vec_idx), result_ptr + vec_idx);
    }}

    return result;
}}
"""

    def _generate_nd_wrapper(self, rec: "Recurrence") -> str:
        """Generate wrapper for N-dimensional recurrence (N >= 3)."""
        func_name = rec.name.lower()
        idx_params = ", ".join(f"int {idx}" for idx in rec.indices)
        dispatcher = f"dispatch_{rec.name}"

        if rec.runtime_vars:
            primary_var = rec.runtime_vars[0]
            other_vars = rec.runtime_vars[1:]
        else:
            primary_var = None
            other_vars = []

        ns = f"{rec.namespace}::" if rec.namespace else ""

        if primary_var is None:
            # No runtime parameters (like Binomial)
            return f"""namespace py = pybind11;

py::array_t<double> {func_name}_wrapper({idx_params}) {{
    py::array_t<double> result(1);
    auto buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);

    Vec8d vec_result = {ns}{dispatcher}({', '.join(rec.indices)});
    ptr[0] = vec_result[0];  // Extract first element only

    return result;
}}
"""

        # With runtime parameters
        other_params = ", ".join(f"py::array_t<double> {v}" for v in other_vars)
        param_list = f"py::array_t<double> {primary_var}" + (f", {other_params}" if other_params else "")

        return f"""namespace py = pybind11;

py::array_t<double> {func_name}_wrapper(
    {idx_params},
    {param_list}
) {{
    auto buf = {primary_var}.request();
    if (buf.ndim != 1) {{
        throw std::runtime_error("Input must be 1D array");
    }}

    size_t arr_size = buf.shape[0];
    py::array_t<double> result(arr_size);
    auto result_buf = result.request();
    double* result_ptr = static_cast<double*>(result_buf.ptr);
    double* input_ptr = static_cast<double*>(buf.ptr);

{self._generate_param_loads(other_vars, indent="    ")}

    // Process in Vec8d chunks
    size_t vec_idx = 0;
    for (; vec_idx + 8 <= arr_size; vec_idx += 8) {{
        Vec8d vec_input;
        vec_input.load(input_ptr + vec_idx);

{self._generate_dispatcher_call(rec, ns, dispatcher, ', '.join(rec.indices), other_vars, indent="        ")}

        result_vec.store(result_ptr + vec_idx);
    }}

    // Handle remaining elements
    if (vec_idx < arr_size) {{
        Vec8d vec_input;
        vec_input.load_partial(static_cast<int>(arr_size - vec_idx), input_ptr + vec_idx);

{self._generate_dispatcher_call(rec, ns, dispatcher, ', '.join(rec.indices), other_vars, indent="        ")}

        result_vec.store_partial(static_cast<int>(arr_size - vec_idx), result_ptr + vec_idx);
    }}

    return result;
}}
"""

    def _generate_param_loads(self, other_vars: List[str], indent: str = "") -> str:
        """Generate code to load other runtime parameters."""
        if not other_vars:
            return ""

        lines = []
        for var in other_vars:
            lines.append(f"{indent}// Load {var} parameter pointers")
            lines.append(f"{indent}auto {var}_buf = {var}.request();")
            lines.append(f"{indent}double* {var}_ptr = static_cast<double*>({var}_buf.ptr);")
            lines.append(f"{indent}bool {var}_is_scalar = ({var}_buf.size == 1);")

        return "\n".join(lines)

    def _generate_dispatcher_call(self, rec: "Recurrence", ns: str,
                                   dispatcher: str, indices: str,
                                   other_vars: List[str], indent: str = "") -> str:
        """Generate the dispatcher function call."""
        if rec.runtime_vars:
            # Load auxiliary parameters for current chunk
            lines = []
            for var in other_vars:
                lines.append(f"{indent}Vec8d vec_{var};")
                lines.append(f"{indent}if ({var}_is_scalar) {{")
                lines.append(f"{indent}    vec_{var} = Vec8d({var}_ptr[0]);")
                lines.append(f"{indent}}} else {{")
                lines.append(f"{indent}    vec_{var}.load({var}_ptr + vec_idx);")
                lines.append(f"{indent}}}")

            # First runtime var is loaded as vec_input, others as vec_{varname}
            args = ["vec_input"] + [f"vec_{v}" for v in other_vars]
            call_args = ", ".join(args)
            lines.append(f"{indent}Vec8d result_vec = {ns}{dispatcher}({indices}, {call_args});")
            return "\n".join(lines)
        else:
            return f"{indent}Vec8d result_vec = {ns}{dispatcher}({indices});"

    def _generate_module_definition(self) -> str:
        """Generate PYBIND11_MODULE definition."""
        bindings = []
        for rec in self.recurrences:
            func_name = rec.name.lower()
            wrapper_name = f"{func_name}_wrapper"

            # Build parameter string for py::arg declarations
            idx_args = ", ".join(f'py::arg("{idx}")' for idx in rec.indices)
            var_args = ", ".join(f'py::arg("{v}")' for v in rec.runtime_vars)
            all_args = ", ".join(filter(None, [idx_args, var_args]))

            # Use shorter docstring to avoid template issues
            if all_args:
                bindings.append(
                    f'    m.def("{func_name}", &{wrapper_name}, "{rec.name}",\n'
                    f'          {all_args});'
                )
            else:
                bindings.append(
                    f'    m.def("{func_name}", &{wrapper_name}, "{rec.name}");'
                )

        return f"""
PYBIND11_MODULE(_{self.module_name}, m) {{
    m.doc() = "RECURSUM {self.module_name} recurrence relations";

{chr(10).join(bindings)}
}}
"""
