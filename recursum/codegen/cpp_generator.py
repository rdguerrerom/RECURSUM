"""
C++ template code generator for recurrence relations.

Generates SFINAE-enabled C++ template specializations that compile
recurrence relations into zero-overhead metaprograms.

Supports optimizations:
- CSE (Common Subexpression Elimination): Caches recursive calls in local variables
- Horner's method: Factors polynomial expressions to minimize multiplications
- Memoization patterns: Generates code with compile-time caching
"""

from typing import TYPE_CHECKING, Optional

from .core import Sum, ScaledExpr, BinOp, Expr, RecursiveCall, Term
from .optimizer import (
    ExpressionOptimizer, OptimizedCodeGenerator, OptimizedExpr,
    CachedVar, OptimizedSum, should_apply_cse, count_operations
)

if TYPE_CHECKING:
    from .recurrence import Recurrence, RecurrenceRule, BaseCase


class CppGenerator:
    """Generate C++ template code from Recurrence definitions.

    The generator supports several optimization modes:
    - 'none': No optimizations, direct code generation
    - 'cse': Common Subexpression Elimination only
    - 'full': CSE + Horner's method + strength reduction

    Optimizations are applied at the AST level before code generation,
    ensuring they work for ALL recurrence relations defined in RECURSUM.
    """

    def __init__(self, rec: "Recurrence", optimization: str = "cse"):
        """
        Initialize generator.

        Args:
            rec: Recurrence definition to generate code for
            optimization: Optimization level ('none', 'cse', 'full')
        """
        self.rec = rec
        self.ctx = rec._ctx()
        self.optimization = optimization

        # Initialize optimizer based on optimization level
        if optimization == "none":
            self.optimizer = None
        else:
            self.optimizer = ExpressionOptimizer(
                self.ctx,
                enable_cse=(optimization in ("cse", "full")),
                enable_horner=(optimization == "full"),
                cse_threshold=2
            )

    def generate(self) -> str:
        """
        Generate complete C++ header file.

        Returns:
            String containing full C++ template code with:
            - Header guards and includes
            - Primary template (fallback for invalid indices)
            - Base case specializations
            - Recurrence rule specializations (sorted by priority)
        """
        parts = [self._header(), self._primary_template()]

        # Add base cases
        for bc in self.rec._base_cases:
            parts.append(self._base_case(bc))

        # Add recurrence rules (sorted by priority)
        for rule in sorted(self.rec._rules, key=lambda r: r.priority_key()):
            parts.append(self._rule(rule))

        parts.append(self._footer())
        return "\n\n".join(filter(None, parts))

    def _header(self) -> str:
        """Generate file header with includes and namespace open.

        Includes portable force-inline macro for cross-platform optimization.
        The RECURSUM_FORCEINLINE macro ensures compile-time inlining on all
        major compilers (GCC, Clang, MSVC).
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

    def _primary_template(self) -> str:
        """
        Generate primary template (default/fallback).

        This template returns 0 for any indices that don't match
        a base case or recurrence rule specialization.
        All compute methods are force-inlined for optimal performance.
        """
        tparams = ", ".join(f"int {idx}" for idx in self.rec.indices)
        unused = ", ".join(f"{self.rec.vec_type} /*{v}*/" for v in self.rec.runtime_vars)

        return f"""template<{tparams}, typename Enable = void>
struct {self.ctx.struct_name} {{
    static RECURSUM_FORCEINLINE {self.rec.vec_type} compute({unused}) {{
        return {self.rec.vec_type}(0.0);
    }}
}};"""

    def _base_case(self, bc: "BaseCase") -> str:
        """
        Generate template specialization for a base case.

        Args:
            bc: Base case with specific index values

        Returns:
            C++ template specialization code with force-inlined compute method
        """
        targs = ", ".join(str(bc.index_values.get(idx, idx)) for idx in self.rec.indices)
        val = bc.value.to_cpp(self.ctx)

        # Determine which runtime vars are actually used in the value expression
        used_vars = set()
        for v in self.rec.runtime_vars:
            if bc.value.uses_var(v):
                used_vars.add(v)

        # Generate parameter list with unused markers for truly unused params
        params = []
        for v in self.rec.runtime_vars:
            if v in used_vars:
                params.append(f"{self.rec.vec_type} {v}")
            else:
                params.append(f"{self.rec.vec_type} /*{v}*/")
        param_str = ", ".join(params)

        return f"""template<>
struct {self.ctx.struct_name}<{targs}, void> {{
    static RECURSUM_FORCEINLINE {self.rec.vec_type} compute({param_str}) {{
        return {val};
    }}
}};"""

    def _rule(self, rule: "RecurrenceRule") -> str:
        """
        Generate SFINAE-enabled template specialization for a recurrence rule.

        Args:
            rule: Recurrence rule with constraints and expression

        Returns:
            C++ template specialization with std::enable_if and force-inlined compute
        """
        tparams = ", ".join(f"int {idx}" for idx in self.rec.indices)
        targs = ", ".join(self.rec.indices)
        sig = ", ".join(f"{self.rec.vec_type} {v}" for v in self.rec.runtime_vars)

        # Combine rule constraints with validity constraints
        sfinae = rule.constraints.to_sfinae()
        if self.rec._validity:
            sfinae = f"{sfinae} && {self.rec._validity.to_sfinae()}"

        body = self._body(rule)
        comment = f"        // {rule.name}\n" if rule.name else ""

        return f"""template<{tparams}>
struct {self.ctx.struct_name}<
    {targs},
    typename std::enable_if<{sfinae}>::type
> {{
    static RECURSUM_FORCEINLINE {self.rec.vec_type} compute({sig}) {{
{comment}{body}
    }}
}};"""

    def _body(self, rule: "RecurrenceRule") -> str:
        """
        Generate function body for a recurrence rule.

        Uses different strategies based on expression complexity:
        - Simple expressions (<=3 calls): inline single return
        - With CSE optimization: cache unique recursive calls in local variables
        - Sum expressions: intermediate variables for each term
        - Scaled sums: intermediate variables + final scaling
        - Branch averaging: separate variables for each branch

        Args:
            rule: Recurrence rule to generate body for

        Returns:
            C++ function body (indented, ready for template)
        """
        expr = rule.expression
        calls = expr.collect_calls()

        # Use optimizer if available and beneficial
        if self.optimizer and should_apply_cse(expr):
            return self._optimized_body(expr)

        # Simple case: inline return
        if len(calls) <= 3:
            return f"        return {expr.to_cpp(self.ctx)};"

        lines = []

        # Strategy for Sum expressions
        if isinstance(expr, Sum):
            for i, term in enumerate(expr.terms):
                lines.append(f"        {self.rec.vec_type} t{i+1} = {term.to_cpp(self.ctx)};")
            vars_str = " + ".join(f"t{i+1}" for i in range(len(expr.terms)))
            lines.append(f"        return {vars_str};")

        # Strategy for ScaledExpr with Sum inside
        elif isinstance(expr, ScaledExpr) and isinstance(expr.expr, Sum):
            inner = expr.expr
            for i, term in enumerate(inner.terms):
                lines.append(f"        {self.rec.vec_type} t{i+1} = {term.to_cpp(self.ctx)};")
            vars_str = " + ".join(f"t{i+1}" for i in range(len(inner.terms)))
            s = expr.scale.to_cpp(self.ctx)
            op = "/" if expr.is_division else "*"
            lines.append(f"        return ({vars_str}) {op} {s};")

        # Strategy for branch averaging (complex BinOp pattern)
        elif isinstance(expr, BinOp) and expr.op == '*':
            if isinstance(expr.left, BinOp) and expr.left.op == '+':
                if isinstance(expr.left.left, Sum) and isinstance(expr.left.right, Sum):
                    s1, s2 = expr.left.left, expr.left.right
                    lines.append("        // Branch A")
                    for i, t in enumerate(s1.terms):
                        lines.append(f"        {self.rec.vec_type} a{i+1} = {t.to_cpp(self.ctx)};")
                    lines.append("        // Branch B")
                    for i, t in enumerate(s2.terms):
                        lines.append(f"        {self.rec.vec_type} b{i+1} = {t.to_cpp(self.ctx)};")
                    a_str = " + ".join(f"a{i+1}" for i in range(len(s1.terms)))
                    b_str = " + ".join(f"b{i+1}" for i in range(len(s2.terms)))
                    scale = expr.right.to_cpp(self.ctx)
                    lines.append(f"        return ({a_str} + {b_str}) * {scale};")
                    return "\n".join(lines)
            lines.append(f"        return {expr.to_cpp(self.ctx)};")

        # Fallback: inline return
        else:
            lines.append(f"        return {expr.to_cpp(self.ctx)};")

        return "\n".join(lines)

    def _optimized_body(self, expr: Expr) -> str:
        """
        Generate an optimized function body using CSE.

        This method:
        1. Identifies all unique recursive calls in the expression
        2. Caches each unique call in a local variable
        3. Reconstructs the expression using cached variables
        4. Returns the optimized code

        Args:
            expr: Expression to optimize

        Returns:
            Optimized C++ function body
        """
        opt_result = self.optimizer.optimize_expression(expr)

        lines = []

        # Add CSE comment if we have intermediates
        if opt_result.intermediates:
            lines.append("        // CSE: Cache recursive calls")

        # Generate intermediate variable declarations
        for var_name, intermediate in opt_result.intermediates:
            if isinstance(intermediate, RecursiveCall):
                call_cpp = intermediate.to_cpp(self.ctx)
                lines.append(f"        {self.rec.vec_type} {var_name} = {call_cpp};")
            else:
                lines.append(f"        {self.rec.vec_type} {var_name} = {intermediate.to_cpp(self.ctx)};")

        # Handle the result expression
        if opt_result.result_expr:
            result = opt_result.result_expr

            # Check if we need further decomposition for scaled expressions
            if isinstance(result, ScaledExpr):
                inner_cpp = result.expr.to_cpp(self.ctx)
                scale_cpp = result.scale.to_cpp(self.ctx)
                op = "/" if result.is_division else "*"
                lines.append(f"        return ({inner_cpp}) {op} {scale_cpp};")
            else:
                lines.append(f"        return {result.to_cpp(self.ctx)};")

        return "\n".join(lines)
