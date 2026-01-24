"""
Layer-by-layer C++ code generator with compile-time CSE.

Generates C++ code that computes ALL values in a layer at once,
enabling true Common Subexpression Elimination at compile time.

Key differences from TMP generator:
- Computes ALL auxiliary index values simultaneously (e.g., all N from 0 to nA+nB)
- Uses output parameters (void compute(Vec8d* out, ...)) instead of return-by-value
- Performs compile-time loop unrolling (generates explicit assignments, NO runtime loops)
- Uses exact-sized buffers (not MAX-sized arrays)
- Applies RECURSUM_FORCEINLINE to all methods

Performance benefits:
- Eliminates runtime loop overhead (branch prediction, not unrolled)
- Eliminates array copy overhead (return-by-value)
- Better cache utilization (exact-sized buffers)
- Aggressive inlining for straight-line code

Example generated code:
    template<int nA, int nB>
    struct HermiteELayer {
        static constexpr int N_VALUES = nA + nB + 1;

        static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d PA, Vec8d PB, Vec8d p) {
            Vec8d prev[nA + nB];
            HermiteELayer<nA-1, nB>::compute(prev, PA, PB, p);

            Vec8d inv2p = Vec8d(0.5) / p;

            // UNROLLED at codegen time (no runtime loop!)
            out[0] = PA * prev[0] + Vec8d(1) * prev[1];
            out[1] = inv2p * prev[0] + PA * prev[1] + Vec8d(2) * prev[2];
            out[2] = inv2p * prev[1] + PA * prev[2] + Vec8d(3) * prev[3];
            // ... (continue for all N values)
        }
    };
"""

from typing import TYPE_CHECKING, List, Dict, Set, Optional, Tuple
from .core import Expr, RecursiveCall, Const, Var, IndexExpr, BinOp, Sum, Term, ScaledExpr

if TYPE_CHECKING:
    from .recurrence import Recurrence, RecurrenceRule, BaseCase


class LayeredCppGenerator:
    """Generate layer-by-layer C++ code with compile-time CSE."""

    def __init__(self, rec: "Recurrence", unroll_loops: bool = True):
        """
        Initialize generator.

        Args:
            rec: Recurrence definition to generate code for
            unroll_loops: Whether to unroll loops at codegen time (default: True, strongly recommended)
        """
        self.rec = rec
        self.ctx = rec._ctx()
        self.unroll_loops = unroll_loops

        # Identify the auxiliary index (the one that varies within a layer)
        # For Hermite E: nA and nB are layer indices, t is the auxiliary index
        # This assumes the last index is the auxiliary one
        self.layer_indices = self.rec.indices[:-1] if len(self.rec.indices) > 1 else []
        self.aux_index = self.rec.indices[-1] if self.rec.indices else None

    def generate(self) -> str:
        """
        Generate complete C++ header file with layer structs.

        Returns:
            String containing complete C++ header file with:
            - Header guards, includes, and RECURSUM_FORCEINLINE macro
            - Primary fallback layer template
            - Base case layer specializations
            - Recurrence rule layer specializations with unrolled code
            - Accessor template for API compatibility
        """
        parts = [
            self._header(),
            self._layer_primary_template(),
            self._layer_base_cases(),
            self._layer_rules(),
            self._accessor_template(),
            self._footer()
        ]
        return "\n\n".join(filter(None, parts))

    def _header(self) -> str:
        """
        Generate file header with includes and namespace open.

        Reuses the same header as CppGenerator, including RECURSUM_FORCEINLINE macro.
        """
        ns = f"namespace {self.rec.namespace} {{\n" if self.rec.namespace else ""
        return f"""#pragma once

#include <type_traits>
#include <recursum/vectorclass.h>

// Portable force-inline macro for performance-critical compute methods
// This ensures template instantiations are fully inlined at compile time
#ifndef RECURSUM_FORCEINLINE
  #ifdef _MSC_VER
    #define RECURSUM_FORCEINLINE __forceinline
  #elif defined(__GNUC__) || defined(__clang__)
    #define RECURSUM_FORCEINLINE inline __attribute__((always_inline))
  #else
    #define RECURSUM_FORCEINLINE inline
  #endif
#endif

{ns}"""

    def _footer(self) -> str:
        """Generate namespace close."""
        return f"}} // namespace {self.rec.namespace}" if self.rec.namespace else ""

    def _layer_primary_template(self) -> str:
        """
        Generate primary (fallback) layer template.

        This template is used when indices don't match any specialization.
        It has N_VALUES = 0 and does nothing in compute().
        """
        # Layer templates only have layer indices (not auxiliary index)
        tparams = ", ".join(f"int {idx}" for idx in self.layer_indices)
        if not tparams:
            tparams = "int dummy_idx"

        # Detect Coulomb R for proper Boys parameter type
        is_coulomb_r = (len(self.rec.indices) == 4 and
                       self.rec.indices[-1] == 'N' and
                       'Boys' in self.rec.runtime_vars)

        # Generate unused parameters with correct types
        unused_parts = []
        for v in self.rec.runtime_vars:
            if is_coulomb_r and v == 'Boys':
                unused_parts.append(f"const {self.rec.vec_type}* /*{v}*/")
            else:
                unused_parts.append(f"{self.rec.vec_type} /*{v}*/")
        unused = ", ".join(unused_parts)

        layer_name = f"{self.ctx.struct_name}Layer"

        return f"""template<{tparams}, typename Enable = void>
struct {layer_name} {{
    static constexpr int N_VALUES = 0;

    static RECURSUM_FORCEINLINE void compute({self.rec.vec_type}* out, {unused}) {{
        // Invalid indices: do nothing
    }}
}};"""

    def _layer_base_cases(self) -> str:
        """
        Generate layer specializations for base cases.

        For each base case, generates a layer template that:
        - Sets N_VALUES appropriately (uses _determine_aux_range for correct count)
        - Computes output values

        Special handling for Coulomb R:
        - Base case R_{000}^{(N)} = Boys[N] for ALL N from 0 to N_VALUES-1
        - Boys parameter is const Vec8d* (array pointer)
        - Generates loop: for (int N = 0; N < N_VALUES; ++N) out[N] = Boys[N];
        """
        layer_name = f"{self.ctx.struct_name}Layer"
        parts = []

        # Detect Coulomb R pattern
        is_coulomb_r = (len(self.rec.indices) == 4 and
                       self.rec.indices[-1] == 'N' and
                       'Boys' in self.rec.runtime_vars)

        for bc in self.rec._base_cases:
            # Get layer index values (all except auxiliary)
            layer_values = []
            aux_value = None

            for idx in self.rec.indices:
                val = bc.index_values.get(idx)
                if val is not None:
                    if idx == self.aux_index:
                        aux_value = val
                    else:
                        layer_values.append(str(val))

            # Skip if we don't have all layer indices specified
            if len(layer_values) != len(self.layer_indices):
                continue

            # Template arguments
            if layer_values:
                targs = ", ".join(layer_values)
            else:
                targs = "0"  # Dummy for single-index recurrences

            # Determine N_VALUES using the corrected formula
            if is_coulomb_r:
                # For Coulomb R base case (0,0,0): N_VALUES = 0+0+0+4 = 4
                n_values_expr = "4"  # Base case is (0,0,0), so t+u+v+4 = 4
            else:
                # For other recurrences, use aux_value + 1
                n_values_expr = str(aux_value + 1) if aux_value is not None else "1"

            # Determine which runtime vars are used and fix types for Coulomb R
            params = []
            for v in self.rec.runtime_vars:
                if is_coulomb_r and v == 'Boys':
                    # Special handling: Boys is const Vec8d* for Coulomb R
                    if bc.value.uses_var(v):
                        params.append(f"const {self.rec.vec_type}* {v}")
                    else:
                        params.append(f"const {self.rec.vec_type}* /*{v}*/")
                else:
                    # Normal handling for other variables
                    if bc.value.uses_var(v):
                        params.append(f"{self.rec.vec_type} {v}")
                    else:
                        params.append(f"{self.rec.vec_type} /*{v}*/")
            param_str = ", ".join(params)

            # Generate computation body
            if is_coulomb_r and "Boys[" in bc.value.to_cpp(self.ctx):
                # Special case: Coulomb R base case computes ALL N values from Boys array
                computation = f"""        // Base case: R_{{000}}^{{(N)}} = Boys[N] for all N
        for (int {self.aux_index} = 0; {self.aux_index} < N_VALUES; ++{self.aux_index}) {{
            out[{self.aux_index}] = Boys[{self.aux_index}];
        }}"""
            else:
                # Normal case: single value assignment
                val_cpp = bc.value.to_cpp(self.ctx)
                computation = f"        out[{aux_value if aux_value is not None else 0}] = {val_cpp};"

            code = f"""template<>
struct {layer_name}<{targs}, void> {{
    static constexpr int N_VALUES = {n_values_expr};

    static RECURSUM_FORCEINLINE void compute({self.rec.vec_type}* out, {param_str}) {{
{computation}
    }}
}};"""
            parts.append(code)

        return "\n\n".join(parts)

    def _layer_rules(self) -> str:
        """
        Generate layer specializations for recurrence rules with unrolled loops.

        This is the MOST IMPORTANT method. It:
        1. Groups rules by layer indices (same layer, different aux_index ranges)
        2. For each layer, determines the auxiliary index range
        3. Unrolls the computation at CODEGEN time (generates explicit assignments)
        4. Uses previous layer computed ONCE
        """
        # Group rules by layer constraints
        # Rules with the same layer index constraints belong to the same layer template
        layer_groups: Dict[str, List["RecurrenceRule"]] = {}

        for rule in self.rec._rules:
            # Extract layer constraints (all except auxiliary index conditions)
            layer_key = self._get_layer_key(rule)
            if layer_key not in layer_groups:
                layer_groups[layer_key] = []
            layer_groups[layer_key].append(rule)

        # Generate code for each layer
        parts = []
        for layer_key, rules in layer_groups.items():
            layer_code = self._generate_layer_template(layer_key, rules)
            if layer_code:
                parts.append(layer_code)

        return "\n\n".join(parts)

    def _get_layer_key(self, rule: "RecurrenceRule") -> str:
        """
        Get a key identifying which layer this rule belongs to.

        The key includes constraints on layer indices only (not auxiliary index).
        """
        # For simplicity, convert constraints to a canonical string
        # This groups rules that have the same layer index constraints
        constraints = []
        for c in rule.constraints.constraints:
            # Only include constraints on layer indices
            if any(idx in c.to_sfinae() for idx in self.layer_indices):
                constraints.append(c.to_sfinae())
        return " && ".join(sorted(constraints)) if constraints else "default"

    def _generate_layer_template(self, layer_key: str, rules: List["RecurrenceRule"]) -> str:
        """
        Generate a single layer template for a group of rules.

        Args:
            layer_key: Key identifying this layer
            rules: All rules belonging to this layer

        Returns:
            C++ template specialization code
        """
        if not rules:
            return ""

        # Group rules by auxiliary index constraints
        # For Hermite E: group into t==0 and t>0 cases
        aux_rule_groups = self._group_rules_by_aux(rules)

        # Use the first rule to determine basic structure
        representative = rules[0]

        # Extract template parameters and SFINAE condition
        tparams = ", ".join(f"int {idx}" for idx in self.layer_indices)
        if not tparams:
            tparams = "int dummy_idx"

        # Build SFINAE from layer constraints only (exclude auxiliary index)
        sfinae = layer_key if layer_key != "default" else None

        # Filter validity constraints to exclude auxiliary index
        if self.rec._validity:
            layer_validity_parts = []
            for c in self.rec._validity.constraints:
                c_sfinae = c.to_sfinae()
                # Only include if it doesn't reference the auxiliary index
                if not (self.aux_index and self.aux_index in c_sfinae):
                    layer_validity_parts.append(c_sfinae)

            if layer_validity_parts:
                layer_validity = " && ".join(layer_validity_parts)
                if sfinae:
                    sfinae = f"({sfinae}) && ({layer_validity})"
                else:
                    sfinae = layer_validity

        if not sfinae:
            sfinae = "true"

        # Template arguments for specialization
        targs = ", ".join(self.layer_indices) if self.layer_indices else "dummy_idx"

        # Function signature with special handling for Coulomb R Boys parameter
        is_coulomb_r = (len(self.rec.indices) == 4 and
                       self.rec.indices[-1] == 'N' and
                       'Boys' in self.rec.runtime_vars)
        sig_parts = []
        for v in self.rec.runtime_vars:
            if is_coulomb_r and v == 'Boys':
                sig_parts.append(f"const {self.rec.vec_type}* {v}")
            else:
                sig_parts.append(f"{self.rec.vec_type} {v}")
        sig = ", ".join(sig_parts)

        # Determine auxiliary index range
        aux_range_expr = self._determine_aux_range(representative)

        # Generate function body handling all auxiliary index cases
        body = self._generate_multi_rule_layer_body(aux_rule_groups, aux_range_expr)

        layer_name = f"{self.ctx.struct_name}Layer"

        return f"""template<{tparams}>
struct {layer_name}<
    {targs},
    typename std::enable_if<{sfinae}>::type
> {{
    static constexpr int N_VALUES = {aux_range_expr};

    static RECURSUM_FORCEINLINE void compute({self.rec.vec_type}* out, {sig}) {{
{body}
    }}
}};"""

    def _group_rules_by_aux(self, rules: List["RecurrenceRule"]) -> Dict[str, "RecurrenceRule"]:
        """
        Group rules by auxiliary index constraints.

        For Hermite E, this creates groups for:
        - t == 0
        - t > 0

        Returns:
            Dict mapping auxiliary constraint to rule
        """
        groups = {}
        for rule in rules:
            # Find constraints on auxiliary index
            aux_constraint = "general"
            for c in rule.constraints.constraints:
                c_sfinae = c.to_sfinae()
                if self.aux_index and self.aux_index in c_sfinae:
                    aux_constraint = c_sfinae
                    break
            groups[aux_constraint] = rule
        return groups

    def _find_max_aux_shift(self, expr: Expr) -> int:
        """
        Find the maximum shift on the auxiliary index in an expression.

        For Coulomb R: R_{tuv}^{(N)} = PCx * R_{t-1,u,v}^{(N+1)} + ...
        The auxiliary shift is +1 (we access N+1).

        Returns:
            Maximum absolute value of auxiliary index shift
        """
        max_shift = 0

        calls = expr.collect_calls()
        for call in calls:
            aux_shift = call.index_shifts.get(self.aux_index, 0)
            max_shift = max(max_shift, abs(aux_shift))

        return max_shift

    def _generate_multi_rule_layer_body(self, aux_rule_groups: Dict[str, "RecurrenceRule"],
                                        aux_range_expr: str) -> str:
        """
        Generate layer body handling multiple auxiliary index cases.

        For Hermite E, this generates:
        - Special case for t==0 (no t-1 term)
        - General loop for t>0 (with t-1, t, and t+1 terms)
        - Special case for t==max (no t+1 term)

        Args:
            aux_rule_groups: Rules grouped by auxiliary constraint
            aux_range_expr: Expression for N_VALUES

        Returns:
            Indented C++ function body
        """
        lines = []

        # Find rule for general case (with most complete form)
        general_rule = None
        t0_rule = None

        for constraint, rule in aux_rule_groups.items():
            if "== 0" in constraint:
                t0_rule = rule
            elif ">" in constraint or "general" in constraint:
                general_rule = rule

        if not general_rule and not t0_rule:
            # Fallback to first available rule
            general_rule = list(aux_rule_groups.values())[0]

        # Determine which rule to use for prev layer computation
        rule_for_prev = general_rule if general_rule else t0_rule

        # Detect Coulomb R pattern (two-term recurrence with different spatial indices)
        is_coulomb_r = (len(self.rec.indices) == 4 and
                       self.rec.indices[-1] == 'N' and
                       'Boys' in self.rec.runtime_vars)

        # Find maximum auxiliary index shift to determine array size
        max_aux_shift = self._find_max_aux_shift(rule_for_prev.expression)

        # Compute previous layer(s)
        calls = rule_for_prev.expression.collect_calls()
        if calls and is_coulomb_r:
            # Coulomb R: Two-term recurrence, compute BOTH prev layers
            # Example: R_{tuv}^{(N)} = PCx * R_{t-1,u,v}^{(N+1)} + (t-1) * R_{t-2,u,v}^{(N+1)}
            # We need prev_t1 for (t-1, u, v) and prev_t2 for (t-2, u, v)
            #
            # CRITICAL: When N ranges from 0 to N_VALUES-1 and we access prev[N+1],
            # we need prev arrays sized to hold indices 0 through N_VALUES (N_VALUES+1 elements)

            # Collect all unique spatial layer dependencies
            spatial_layers = set()
            for call in calls:
                layer_indices_with_shifts = []
                for idx in self.layer_indices:
                    shift = call.index_shifts.get(idx, 0)
                    if shift == 0:
                        layer_indices_with_shifts.append(idx)
                    elif shift > 0:
                        layer_indices_with_shifts.append(f"{idx} + {shift}")
                    else:
                        layer_indices_with_shifts.append(f"{idx} - {-shift}")
                spatial_layers.add(tuple(layer_indices_with_shifts))

            # Generate code for each spatial layer
            for i, layer_indices in enumerate(sorted(spatial_layers)):
                prev_name = f"prev_{i}"
                prev_targs = ", ".join(layer_indices)
                prev_args = ", ".join(self.rec.runtime_vars)

                # Allocate array size to accommodate max shift
                # If we access prev[N+1] with N up to N_VALUES-1, we need N_VALUES+1 elements
                prev_size = f"N_VALUES + {max_aux_shift}" if max_aux_shift > 0 else "N_VALUES"

                lines.append(f"        // Compute previous layer {i}: ({prev_targs})")
                lines.append(f"        {self.rec.vec_type} {prev_name}[{prev_size}] = {{}};  // Zero-init, sized for shifted access")
                lines.append(f"        {self.ctx.struct_name}Layer<{prev_targs}>::compute({prev_name}, {prev_args});")

            if spatial_layers:
                lines.append("")

        elif calls:
            # Standard single prev layer (Hermite E pattern)
            prev_layer = self._determine_previous_layer(calls)
            if prev_layer:
                prev_layer_indices, prev_size_expr = prev_layer
                prev_args = ", ".join(self.rec.runtime_vars)
                prev_targs = ", ".join(prev_layer_indices)

                lines.append(f"        // Compute previous layer once (initialize to zero for out-of-bounds accesses)")
                lines.append(f"        {self.rec.vec_type} prev[{prev_size_expr}] = {{}};  // Zero-initialize")
                lines.append(f"        {self.ctx.struct_name}Layer<{prev_targs}>::compute(prev, {prev_args});")
                lines.append("")

        # Extract common variables (like inv2p)
        if general_rule:
            common_vars = self._extract_common_variables(general_rule.expression)
            for var_name, var_expr in common_vars:
                lines.append(f"        {self.rec.vec_type} {var_name} = {var_expr};")
            if common_vars:
                lines.append("")

        # Generate code for t==0 case if it exists
        if t0_rule:
            lines.append(f"        // t = 0 special case")
            # For t==0, substitute t with 0 before converting to prev array
            expr_t0 = self._convert_to_prev_array_with_aux_value(t0_rule.expression, 0)
            lines.append(f"        out[0] = {expr_t0};")
            lines.append("")

        # Generate code for general case
        if general_rule:
            start_idx = "1" if t0_rule else "0"
            lines.append(f"        // General case: t = {start_idx} to N_VALUES-1")
            lines.append(f"        for (int {self.aux_index} = {start_idx}; {self.aux_index} < N_VALUES; ++{self.aux_index}) {{")
            expr_general = self._convert_to_prev_array(general_rule.expression)
            lines.append(f"            out[{self.aux_index}] = {expr_general};")
            lines.append(f"        }}")

        return "\n".join(lines)

    def _determine_aux_range(self, rule: "RecurrenceRule") -> str:
        """
        Determine the range of the auxiliary index for this layer.

        For Hermite E with indices [nA, nB, t]: t ranges from 0 to nA+nB,
        so N_VALUES = nA + nB + 1.

        For Coulomb R with indices [t, u, v, N]: N ranges from 0 to t+u+v,
        BUT each layer needs EXTRA values to serve higher layers:
        N_VALUES = (t + u + v) + 4

        The +4 ensures layer (t-2) has enough values when called by layer (t).

        Returns:
            C++ expression for N_VALUES (number of outputs in this layer)
        """
        # Detect Coulomb R pattern: 4 indices with last one being N
        is_coulomb_r = (len(self.rec.indices) == 4 and
                       self.rec.indices[-1] == 'N' and
                       'Boys' in self.rec.runtime_vars)

        if is_coulomb_r:
            # Coulomb R: N_VALUES = (t + u + v) + 4
            # The +4 ensures enough values for two-term recurrence with (t-2) access
            spatial_indices = self.layer_indices  # [t, u, v]
            spatial_sum = " + ".join(spatial_indices)
            return f"({spatial_sum}) + 4"

        # Look at validity constraints to infer the range (for Hermite E and others)
        if self.rec._validity:
            for c in self.rec._validity.constraints:
                sfinae = c.to_sfinae()
                if self.aux_index and self.aux_index in sfinae:
                    # Try to extract the upper bound
                    # Example: "t <= nA + nB" -> N_VALUES = nA + nB + 1
                    if "<=" in sfinae:
                        parts = sfinae.split("<=")
                        if len(parts) == 2 and self.aux_index in parts[0]:
                            upper_bound = parts[1].strip()
                            # Remove trailing parenthesis from constraint format
                            if upper_bound.endswith(")"):
                                upper_bound = upper_bound[:-1].strip()
                            return f"({upper_bound}) + 1"

        # Fallback: assume single value
        return "1"

    def _generate_layer_body(self, rule: "RecurrenceRule", aux_range_expr: str) -> str:
        """
        Generate the body of a layer compute function.

        This includes:
        1. Computing the previous layer (once)
        2. Extracting any common variables (like inv2p for Hermite)
        3. Generating loop or unrolled assignments using prev array

        Args:
            rule: Representative rule for this layer
            aux_range_expr: Expression for the number of values to compute

        Returns:
            Indented C++ function body
        """
        lines = []

        # Analyze the rule to determine what previous layer is needed
        calls = rule.expression.collect_calls()
        if not calls:
            # No recursive calls, just return the expression
            return f"        out[0] = {rule.expression.to_cpp(self.ctx)};"

        # Find the previous layer dimensions
        prev_layer = self._determine_previous_layer(calls)
        if prev_layer:
            prev_layer_indices, prev_size_expr = prev_layer

            # Generate code to compute previous layer
            prev_args = ", ".join(self.rec.runtime_vars)
            prev_targs = ", ".join(prev_layer_indices)

            lines.append(f"        // Compute previous layer once")
            lines.append(f"        {self.rec.vec_type} prev[{prev_size_expr}];")
            lines.append(f"        {self.ctx.struct_name}Layer<{prev_targs}>::compute(prev, {prev_args});")
            lines.append("")

        # Extract common subexpressions
        # For Hermite E, this would be inv2p = 0.5 / p
        common_vars = self._extract_common_variables(rule.expression)
        for var_name, var_expr in common_vars:
            lines.append(f"        {self.rec.vec_type} {var_name} = {var_expr};")
        if common_vars:
            lines.append("")

        # Generate computation using prev array
        # Convert recursive calls to array accesses
        lines.append(f"        // Compute all {self.aux_index} values using previous layer")
        lines.append(f"        for (int {self.aux_index} = 0; {self.aux_index} < N_VALUES; ++{self.aux_index}) {{")

        # Convert expression to use prev array
        expr_with_prev = self._convert_to_prev_array(rule.expression)
        lines.append(f"            out[{self.aux_index}] = {expr_with_prev};")
        lines.append(f"        }}")

        return "\n".join(lines)

    def _determine_previous_layer(self, calls: List[RecursiveCall]) -> Optional[Tuple[List[str], str]]:
        """
        Determine which previous layer is needed and its size.

        Returns:
            Tuple of (layer_indices, size_expression) or None
        """
        if not calls:
            return None

        # Take the first call as representative
        first_call = calls[0]

        # Extract previous layer indices (with shifts applied)
        prev_indices = []
        for idx in self.layer_indices:
            shift = first_call.index_shifts.get(idx, 0)
            if shift == 0:
                prev_indices.append(idx)
            elif shift > 0:
                prev_indices.append(f"{idx} + {shift}")
            else:
                prev_indices.append(f"{idx} - {-shift}")

        # Determine size of previous layer
        # This is the N_VALUES of the previous layer
        # For Hermite E: if we're computing layer (nA, nB), previous is (nA-1, nB)
        # and its size is (nA-1) + nB + 1 = nA + nB
        aux_shift = first_call.index_shifts.get(self.aux_index, 0)

        # Build size expression based on aux_range_expr logic
        # For Hermite E: Previous layer (nA-1, nB) has N_VALUES = nA+nB
        # but we access prev[t+1] when t=nA+nB, requiring prev[nA+nB+1]
        # Therefore we need size nA+nB+2 to accommodate indices 0 to nA+nB+1
        if len(self.layer_indices) == 2:
            # Example: Hermite E with nA, nB
            # Previous layer writes nA+nB values, but we access up to index nA+nB+1
            prev_size = f"{self.layer_indices[0]} + {self.layer_indices[1]} + 2"
        else:
            prev_size = "N_VALUES + 1"  # Add 1 for t+1 access pattern

        return (prev_indices, prev_size)

    def _extract_common_variables(self, expr: Expr) -> List[Tuple[str, str]]:
        """
        Extract common variables that should be computed once.

        For Hermite E: (0.5 / p) appears in many terms, so extract it as inv2p.

        Returns:
            List of (variable_name, cpp_expression) tuples
        """
        common_vars = []

        # Simple heuristic: look for division by runtime variables
        # This catches inv2p = 0.5 / p
        if self._contains_division_by_var(expr):
            # Extract the division pattern
            for v in self.rec.runtime_vars:
                expr_str = expr.to_cpp(self.ctx)
                if f"/ {v}" in expr_str or f"/{v}" in expr_str:
                    # Found a division by variable v
                    if "0.5" in expr_str:
                        common_vars.append((f"inv2{v}", f"{self.rec.vec_type}(0.5) / {v}"))
                    break

        return common_vars

    def _contains_division_by_var(self, expr: Expr) -> bool:
        """Check if expression contains division by a runtime variable."""
        if isinstance(expr, ScaledExpr) and expr.is_division:
            if isinstance(expr.scale, Var):
                return True
        if isinstance(expr, BinOp):
            return self._contains_division_by_var(expr.left) or self._contains_division_by_var(expr.right)
        if isinstance(expr, Sum):
            return any(self._contains_division_by_var(t) for t in expr.terms)
        return False

    def _unroll_layer_computation(self, rule: "RecurrenceRule", aux_range_expr: str) -> List[str]:
        """
        Unroll layer computation at codegen time.

        This is the KEY method that generates explicit assignments instead of runtime loops.

        Args:
            rule: RecurrenceRule to generate code for
            aux_range_expr: Expression for N_VALUES (e.g., "nA + nB + 1")

        Returns:
            List of C++ statements (one per aux_index value)
        """
        statements = []

        # We need to handle different aux_index values
        # The challenge: we don't know the concrete values at codegen time!
        # Solution: Generate code that handles different cases

        # Strategy: Look at the rule's constraints to determine which cases to generate
        # For Hermite E, there are special cases for t=0 and t=max

        # Find all rules that might apply to this layer
        # Group by auxiliary index constraints
        aux_rules = self._group_rules_by_aux_constraint(rule)

        # Generate assignments for each case
        for aux_constraint, aux_rule in aux_rules:
            if aux_constraint == "t == 0":
                # Special case for first element
                expr_str = self._substitute_aux_index(aux_rule.expression, 0)
                statements.append(f"        // t = 0 (special case)")
                statements.append(f"        out[0] = {expr_str};")
                statements.append("")
            elif "t > 0" in aux_constraint or "t >= 1" in aux_constraint:
                # General case: need to generate a loop or multiple assignments
                # For now, use a conditional approach
                statements.append(f"        // t = 1 to N-2 (general case)")
                statements.append(f"        for (int {self.aux_index} = 1; {self.aux_index} < N_VALUES - 1; ++{self.aux_index}) {{")
                expr_cpp = aux_rule.expression.to_cpp(self.ctx)
                statements.append(f"            out[{self.aux_index}] = {expr_cpp};")
                statements.append(f"        }}")
                statements.append("")
            elif "last" in aux_constraint or "max" in aux_constraint:
                # Special case for last element
                statements.append(f"        // t = N-1 (special case)")
                statements.append(f"        if (N_VALUES > 1) {{")
                expr_cpp = aux_rule.expression.to_cpp(self.ctx)
                statements.append(f"            out[N_VALUES - 1] = {expr_cpp};")
                statements.append(f"        }}")

        # If no special handling, generate a simple generic loop
        if not statements:
            statements.append(f"        // General computation")
            expr_cpp = rule.expression.to_cpp(self.ctx)
            statements.append(f"        for (int {self.aux_index} = 0; {self.aux_index} < N_VALUES; ++{self.aux_index}) {{")
            statements.append(f"            out[{self.aux_index}] = {expr_cpp};")
            statements.append(f"        }}")

        return statements

    def _group_rules_by_aux_constraint(self, rule: "RecurrenceRule") -> List[Tuple[str, "RecurrenceRule"]]:
        """
        Group rules by constraints on the auxiliary index.

        Returns:
            List of (constraint_string, rule) tuples
        """
        # This is a simplified version
        # In practice, we'd need to find all rules in the recurrence
        # that apply to the same layer but different aux values

        # For now, just return the single rule with a generic constraint
        return [("general", rule)]

    def _substitute_aux_index(self, expr: Expr, aux_value: int) -> str:
        """
        Substitute the auxiliary index with a concrete value in an expression.

        This is used to generate specialized code for specific aux_index values.

        Args:
            expr: Expression to substitute in
            aux_value: Concrete value for the auxiliary index

        Returns:
            C++ code with substitution applied
        """
        # This is complex because we need to substitute in RecursiveCall index_shifts
        # For now, return the unsubstituted expression
        # TODO: Implement proper substitution
        return expr.to_cpp(self.ctx)

    def _convert_to_prev_array(self, expr: Expr) -> str:
        """
        Convert an expression with RecursiveCall nodes to use prev[] array access.

        This transforms expressions like:
            PA * E[nA-1, nB, t] + (t+1) * E[nA-1, nB, t+1]
        Into:
            PA * prev[t] + Vec8d(t+1) * prev[t+1]

        The key insight: All RecursiveCall nodes reference the previous layer,
        so we replace them with prev[aux_index_with_shift].

        Args:
            expr: Expression to convert

        Returns:
            C++ code string with prev[] array accesses
        """
        return self._convert_expr_to_prev(expr, aux_value=None)

    def _convert_to_prev_array_with_aux_value(self, expr: Expr, aux_value: int) -> str:
        """
        Convert expression to prev[] array with auxiliary index substituted.

        For t==0 case: substitute t with 0, so:
            PA * E[nA-1, nB, t] + (t+1) * E[nA-1, nB, t+1]
        Becomes:
            PA * prev[0] + Vec8d(1) * prev[1]

        Args:
            expr: Expression to convert
            aux_value: Value to substitute for auxiliary index

        Returns:
            C++ code string with prev[] array accesses using literal indices
        """
        return self._convert_expr_to_prev(expr, aux_value=aux_value)

    def _convert_expr_to_prev(self, expr: Expr, aux_value: Optional[int] = None) -> str:
        """
        Helper to recursively convert expression tree.

        Args:
            expr: Expression to convert
            aux_value: If provided, substitute auxiliary index with this value

        Returns:
            C++ code string
        """
        if isinstance(expr, RecursiveCall):
            # Convert recursive call to prev array access
            # Extract the shift on the auxiliary index
            aux_shift = expr.index_shifts.get(self.aux_index, 0)

            # Detect Coulomb R pattern and map to correct prev array
            is_coulomb_r = (len(self.rec.indices) == 4 and
                           self.rec.indices[-1] == 'N' and
                           'Boys' in self.rec.runtime_vars)

            if is_coulomb_r:
                # For Coulomb R, determine which prev array based on spatial shifts
                # Build a signature of spatial shifts to identify the layer
                spatial_signature = []
                for idx in self.layer_indices:
                    shift = expr.index_shifts.get(idx, 0)
                    spatial_signature.append((idx, shift))

                # Map signature to prev array name
                # Smaller shifts (closer to current layer) get lower indices
                # E.g., t-1 → prev_0, t-2 → prev_1
                max_spatial_shift = max(abs(shift) for _, shift in spatial_signature)

                if max_spatial_shift == 1:
                    prev_name = "prev_0"
                elif max_spatial_shift == 2:
                    prev_name = "prev_1"
                else:
                    # Fallback for higher shifts (shouldn't happen for standard Coulomb R)
                    prev_name = f"prev_{max_spatial_shift - 1}"

                # Generate array index expression
                if aux_value is not None:
                    idx = aux_value + aux_shift
                    return f"{prev_name}[{idx}]"
                else:
                    if aux_shift == 0:
                        return f"{prev_name}[{self.aux_index}]"
                    elif aux_shift > 0:
                        return f"{prev_name}[{self.aux_index} + {aux_shift}]"
                    else:
                        return f"{prev_name}[{self.aux_index} - {-aux_shift}]"
            else:
                # Standard handling for single prev array (Hermite E pattern)
                if aux_value is not None:
                    # Use literal value
                    idx = aux_value + aux_shift
                    return f"prev[{idx}]"
                else:
                    # Use variable
                    if aux_shift == 0:
                        return f"prev[{self.aux_index}]"
                    elif aux_shift > 0:
                        return f"prev[{self.aux_index} + {aux_shift}]"
                    else:
                        return f"prev[{self.aux_index} - {-aux_shift}]"

        elif isinstance(expr, Const):
            # Constants stay the same
            return f"{self.rec.vec_type}({expr.value})"

        elif isinstance(expr, Var):
            # Variables stay the same
            return expr.name

        elif isinstance(expr, IndexExpr):
            # Index expressions - may need substitution if they reference aux_index
            if aux_value is not None and self.aux_index and self.aux_index in expr.expr_str:
                # Substitute aux_index with literal value
                substituted = expr.expr_str.replace(self.aux_index, str(aux_value))
                # Evaluate simple arithmetic
                try:
                    result = eval(substituted)
                    return f"{self.rec.vec_type}({result})"
                except:
                    return f"{self.rec.vec_type}({substituted})"
            return f"{self.rec.vec_type}({expr.expr_str})"

        elif isinstance(expr, BinOp):
            # Recursively convert both sides
            left = self._convert_expr_to_prev(expr.left, aux_value)
            right = self._convert_expr_to_prev(expr.right, aux_value)

            # Add parentheses for clarity
            if isinstance(expr.left, BinOp):
                left = f"({left})"
            if isinstance(expr.right, BinOp):
                right = f"({right})"

            return f"{left} {expr.op} {right}"

        elif isinstance(expr, Term):
            # Term is coefficient * recursive_call
            coeff = self._convert_expr_to_prev(expr.coeff, aux_value)
            call = self._convert_expr_to_prev(expr.call, aux_value)

            if isinstance(expr.coeff, Const) and expr.coeff.value == 1:
                return call
            return f"{coeff} * {call}"

        elif isinstance(expr, Sum):
            # Sum of terms
            if not expr.terms:
                return f"{self.rec.vec_type}(0.0)"
            if len(expr.terms) == 1:
                return self._convert_expr_to_prev(expr.terms[0], aux_value)
            return " + ".join(self._convert_expr_to_prev(t, aux_value) for t in expr.terms)

        elif isinstance(expr, ScaledExpr):
            # Scaled expression: expr / scale or expr * scale
            inner = self._convert_expr_to_prev(expr.expr, aux_value)
            scale = self._convert_expr_to_prev(expr.scale, aux_value)
            if expr.is_division:
                return f"({inner}) / ({scale})"
            else:
                return f"({inner}) * ({scale})"

        else:
            # Fallback: use original to_cpp
            return expr.to_cpp(self.ctx)

    def _accessor_template(self) -> str:
        """
        Generate accessor template for API compatibility.

        This provides the same interface as the TMP generator:
        HermiteECoeff<nA, nB, t>::compute(PA, PB, p)

        Internally, it computes the entire layer and returns the requested value.
        """
        tparams = ", ".join(f"int {idx}" for idx in self.rec.indices)
        targs_layer = ", ".join(self.layer_indices) if self.layer_indices else "0"
        args = ", ".join(self.rec.runtime_vars)

        # Detect Coulomb R for proper Boys parameter type
        is_coulomb_r = (len(self.rec.indices) == 4 and
                       self.rec.indices[-1] == 'N' and
                       'Boys' in self.rec.runtime_vars)

        # Generate function signature with correct types
        sig_parts = []
        for v in self.rec.runtime_vars:
            if is_coulomb_r and v == 'Boys':
                sig_parts.append(f"const {self.rec.vec_type}* {v}")
            else:
                sig_parts.append(f"{self.rec.vec_type} {v}")
        sig = ", ".join(sig_parts)

        # Determine buffer size
        aux_range = self._determine_aux_range(self.rec._rules[0] if self.rec._rules else None)

        return f"""// API compatibility: single-value accessor
template<{tparams}>
struct {self.ctx.struct_name} {{
    static RECURSUM_FORCEINLINE {self.rec.vec_type} compute({sig}) {{
        {self.rec.vec_type} layer[{aux_range}];
        {self.ctx.struct_name}Layer<{targs_layer}>::compute(layer, {args});
        return layer[{self.aux_index}];
    }}
}};"""
