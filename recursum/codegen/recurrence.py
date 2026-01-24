"""
Main Recurrence class for defining recurrence relations using a fluent API.

Example:
    rec = Recurrence("Hermite", ["i", "j", "t"], ["PA", "PB", "p2"])
    rec.validity("i >= 0", "j >= 0", "i + j >= t")
    rec.base(i=0, j=0, t=0, value=1.0)
    rec.rule("i == 0 && j > 0", "p2 * E[i,j-1,t-1] + PB * E[i,j-1,t]")
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Union

from .core import (
    Expr, Const, Var, BinOp, ScaledExpr,
    CodegenContext
)
from .constraints import ConstraintSet, ConstraintOp
from .parser import RecurrenceParser


# =============================================================================
# Rule and Base Case Classes
# =============================================================================

@dataclass
class BaseCase:
    """A base case for the recurrence with fixed index values."""
    index_values: Dict[str, int]
    value: Expr


@dataclass
class RecurrenceRule:
    """A recurrence rule with constraints and expression."""
    constraints: ConstraintSet
    expression: Expr
    scale: Optional[Expr] = None
    name: str = ""

    def priority_key(self) -> Tuple[int, int]:
        """
        Compute priority for sorting rules.

        Rules with more == constraints have higher priority (lower key),
        followed by rules with more total constraints.

        Returns:
            Tuple of (negative eq_count, negative total_count)
        """
        eq_count = sum(1 for c in self.constraints.constraints
                       if c.op == ConstraintOp.EQ)
        return (-eq_count, -len(self.constraints.constraints))


# =============================================================================
# Main Recurrence Class
# =============================================================================

@dataclass
class Recurrence:
    """
    Fluent API for defining recurrence relations.

    This class provides a domain-specific language for expressing recurrence
    relations that will be compiled to C++ template metaprogramming code.

    Attributes:
        name: Name of the recurrence (e.g., "Hermite")
        indices: List of template parameter names (compile-time integers)
        runtime_vars: List of runtime parameter names (Vec8d values)
        vec_type: SIMD vector type (default: "Vec8d")
        namespace: C++ namespace for generated code
        max_indices: Maximum values for each index (for dispatcher generation)
        scipy_reference: SciPy function name for validation

    Example:
        rec = Recurrence("Legendre", ["n"], ["x"])
        rec.validity("n >= 0")
        rec.base(n=0, value=1.0)
        rec.base(n=1, value="x")
        rec.rule("n > 1", "(2*n-1) * x * E[n-1] + (-(n-1)) * E[n-2]",
                 scale="1/n")
    """
    name: str
    indices: List[str]
    runtime_vars: List[str]
    vec_type: str = "Vec8d"
    namespace: str = ""
    max_indices: Optional[Dict[str, int]] = None
    scipy_reference: Optional[str] = None
    _base_cases: List[BaseCase] = field(default_factory=list)
    _rules: List[RecurrenceRule] = field(default_factory=list)
    _validity: Optional[ConstraintSet] = None

    def __post_init__(self):
        """Initialize max_indices if not provided."""
        if self.max_indices is None:
            # Default: max index of 20 for all indices
            self.max_indices = {idx: 20 for idx in self.indices}

    def validity(self, *constraints: str) -> "Recurrence":
        """
        Set validity constraints for the recurrence.

        These constraints define the domain where the recurrence produces
        non-zero results.

        Args:
            *constraints: Constraint expressions like "n >= 0", "n >= m"

        Returns:
            self for method chaining

        Example:
            rec.validity("n >= 0", "m >= 0", "n >= m")
        """
        self._validity = ConstraintSet.parse(*constraints)
        return self

    def base(self, value: Union[str, int, float, Expr], **index_values: int) -> "Recurrence":
        """
        Add a base case with specific index values.

        Args:
            value: The value for this base case (constant, variable name, or Expr)
            **index_values: Index assignments like n=0, m=1

        Returns:
            self for method chaining

        Example:
            rec.base(n=0, value=1.0)
            rec.base(n=1, value="x")
        """
        # Convert value to Expr if it's not already
        if not isinstance(value, Expr):
            if isinstance(value, str):
                if value in self.runtime_vars:
                    value = Var(value)
                else:
                    try:
                        value = Const(float(value))
                    except:
                        value = Var(value)
            else:
                value = Const(value)

        self._base_cases.append(BaseCase(index_values, value))
        return self

    def rule(self, constraints: Union[str, List[str], ConstraintSet],
             expression: Union[str, Expr],
             scale: Optional[str] = None,
             name: str = "") -> "Recurrence":
        """
        Add a recurrence rule.

        Args:
            constraints: When this rule applies (e.g., "n > 1")
            expression: The recurrence expression (e.g., "x * E[n-1] + E[n-2]")
            scale: Optional scaling factor (e.g., "1/n" for division)
            name: Optional descriptive name for the rule

        Returns:
            self for method chaining

        Example:
            rec.rule("n > 1", "x * E[n-1] + E[n-2]")
            rec.rule("n > 1", "(2*n-1) * x * E[n-1] + (-(n-1)) * E[n-2]",
                     scale="1/n", name="Three-term recurrence")
        """
        # Parse constraints
        if isinstance(constraints, str):
            constraints = ConstraintSet.parse(constraints)
        elif isinstance(constraints, list):
            constraints = ConstraintSet.parse(*constraints)

        # Parse expression
        parser = RecurrenceParser(self.indices, self.runtime_vars)

        if isinstance(expression, str):
            expression = parser.parse_expression(expression)

        # Handle scaling
        scale_expr = None
        if scale:
            scale_expr = parser.parse_scale(scale)
            if scale.startswith("1/"):
                expression = ScaledExpr(expression, scale_expr, is_division=True)
            else:
                expression = ScaledExpr(expression, scale_expr, is_division=False)

        self._rules.append(RecurrenceRule(constraints, expression, scale_expr, name))
        return self

    def branch_average(self, constraints: Union[str, List[str], ConstraintSet],
                       branches: List[str], name: str = "") -> "Recurrence":
        """
        Add a rule that averages multiple equivalent branches.

        This is useful for numerical stability when multiple recurrence
        directions give equivalent results (e.g., McMurchie-Davidson).

        Args:
            constraints: When this rule applies
            branches: List of equivalent expression strings
            name: Optional descriptive name

        Returns:
            self for method chaining

        Example:
            rec.branch_average(
                "nA > 0 && nB > 0",
                ["aAB * E[nA, nB-1, N-1] + PB * E[nA, nB-1, N]",
                 "aAB * E[nA-1, nB, N-1] + PA * E[nA-1, nB, N]"],
                name="Two-branch average"
            )
        """
        # Parse constraints
        if isinstance(constraints, str):
            constraints = ConstraintSet.parse(constraints)
        elif isinstance(constraints, list):
            constraints = ConstraintSet.parse(*constraints)

        # Parse all branches
        parser = RecurrenceParser(self.indices, self.runtime_vars)
        exprs = [parser.parse_expression(b) for b in branches]

        # Combine with averaging
        combined = exprs[0]
        for e in exprs[1:]:
            combined = BinOp('+', combined, e)

        if len(branches) > 1:
            combined = BinOp('*', combined, Const(1.0 / len(branches)))

        self._rules.append(RecurrenceRule(constraints, combined, name=name))
        return self

    def _ctx(self) -> CodegenContext:
        """Create a CodegenContext for this recurrence."""
        return CodegenContext(
            f"{self.name}Coeff",
            self.indices,
            self.runtime_vars,
            self.vec_type
        )

    def generate(self, optimization: str = "cse") -> str:
        """
        Generate C++ template code for this recurrence.

        Args:
            optimization: Optimization level for code generation:
                - 'none': No optimizations, direct code generation
                - 'cse': Common Subexpression Elimination (default)
                - 'full': CSE + Horner's method + strength reduction

        Returns:
            String containing complete C++ header file

        Note:
            This method imports CppGenerator to avoid circular dependencies.

        Example:
            # Generate with CSE optimization (default)
            code = rec.generate()

            # Generate without optimizations
            code_unopt = rec.generate(optimization='none')

            # Generate with full optimizations
            code_full = rec.generate(optimization='full')
        """
        from .cpp_generator import CppGenerator
        return CppGenerator(self, optimization=optimization).generate()

    def generate_layered(self, unroll: bool = True) -> str:
        """
        Generate C++ layer-by-layer code for this recurrence.

        This generates code that computes ALL auxiliary index values at once
        (e.g., all t values for a given nA, nB in Hermite E), enabling true
        compile-time Common Subexpression Elimination.

        Key differences from generate():
        - Uses output parameters: void compute(Vec8d* out, ...) instead of return-by-value
        - Applies RECURSUM_FORCEINLINE to all methods
        - Computes entire layers at once (not individual values)
        - Uses exact-sized buffers (not MAX-sized arrays)

        Performance benefits:
        - Eliminates array copy overhead
        - Better cache utilization
        - Enables aggressive compiler inlining and optimization
        - Each layer computed once and reused (true CSE)

        Args:
            unroll: Whether to structure code for compile-time unrolling (default: True)

        Returns:
            String containing complete C++ header file with layer templates

        Example:
            # Generate Hermite E with layered approach
            rec = hermite_e_coefficient()
            code = rec.generate_layered()

            # This generates HermiteECoeffLayer templates that compute
            # all E^{nA,nB}_t values (t=0 to nA+nB) simultaneously
        """
        from .layered_generator import LayeredCppGenerator
        return LayeredCppGenerator(self, unroll_loops=unroll).generate()
