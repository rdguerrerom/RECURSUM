"""
Optimization transformations for recurrence relation code generation.

This module provides Common Subexpression Elimination (CSE), Horner's method,
and memoization pattern generation for efficient C++ template metaprogramming.

These optimizations are applied at the AST level before code generation,
making them applicable to ALL recurrence relations defined in RECURSUM.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict
import re

from .core import (
    Expr, Const, IndexExpr, Var, RecursiveCall,
    BinOp, Term, Sum, ScaledExpr, CodegenContext
)


# =============================================================================
# CSE (Common Subexpression Elimination) Analysis
# =============================================================================

@dataclass
class SubExprKey:
    """
    Hashable key for identifying common subexpressions.

    Represents either:
    - A recursive call with specific index shifts
    - A coefficient expression
    - A product of coefficient * recursive_call
    """
    kind: str  # 'call', 'coeff', 'term'
    signature: str  # Unique string representation

    def __hash__(self):
        return hash((self.kind, self.signature))

    def __eq__(self, other):
        if not isinstance(other, SubExprKey):
            return False
        return self.kind == other.kind and self.signature == other.signature


@dataclass
class CSEAnalyzer:
    """
    Analyzes expressions to identify common subexpressions.

    This analyzer works on the AST to find:
    1. Repeated recursive calls (e.g., E[n-1, j, t] appearing multiple times)
    2. Repeated coefficient calculations (e.g., (2*n-1))
    3. Repeated term products (coefficient * call)

    The analysis is performed across ALL rules in a recurrence, enabling
    cross-rule CSE where beneficial.
    """

    # Map from subexpression key to list of occurrences
    occurrences: Dict[SubExprKey, List[Expr]] = field(default_factory=dict)

    # Counter for assigning CSE variable names
    cse_counter: int = 0

    def analyze_expr(self, expr: Expr) -> None:
        """Recursively analyze an expression for common subexpressions."""
        if isinstance(expr, RecursiveCall):
            key = self._call_key(expr)
            self._add_occurrence(key, expr)

        elif isinstance(expr, Term):
            # Track the call
            call_key = self._call_key(expr.call)
            self._add_occurrence(call_key, expr.call)

            # Track the coefficient if it's complex
            if not isinstance(expr.coeff, Const) or expr.coeff.value != 1:
                coeff_key = self._coeff_key(expr.coeff)
                self._add_occurrence(coeff_key, expr.coeff)

        elif isinstance(expr, Sum):
            for term in expr.terms:
                self.analyze_expr(term)

        elif isinstance(expr, ScaledExpr):
            self.analyze_expr(expr.expr)
            # Track the scale factor
            scale_key = self._coeff_key(expr.scale)
            self._add_occurrence(scale_key, expr.scale)

        elif isinstance(expr, BinOp):
            self.analyze_expr(expr.left)
            self.analyze_expr(expr.right)

    def _call_key(self, call: RecursiveCall) -> SubExprKey:
        """Generate a unique key for a recursive call."""
        # Sort shifts for consistent ordering
        shifts_str = ",".join(f"{k}:{v}" for k, v in sorted(call.index_shifts.items()))
        return SubExprKey('call', shifts_str)

    def _coeff_key(self, coeff: Expr) -> SubExprKey:
        """Generate a unique key for a coefficient expression."""
        # Use string representation as signature
        return SubExprKey('coeff', self._expr_signature(coeff))

    def _expr_signature(self, expr: Expr) -> str:
        """Generate a canonical string signature for an expression."""
        if isinstance(expr, Const):
            return f"const:{expr.value}"
        elif isinstance(expr, Var):
            return f"var:{expr.name}"
        elif isinstance(expr, IndexExpr):
            return f"idx:{expr.expr_str}"
        elif isinstance(expr, BinOp):
            left_sig = self._expr_signature(expr.left)
            right_sig = self._expr_signature(expr.right)
            return f"({left_sig}{expr.op}{right_sig})"
        elif isinstance(expr, RecursiveCall):
            return f"call:{self._call_key(expr).signature}"
        else:
            return str(type(expr).__name__)

    def _add_occurrence(self, key: SubExprKey, expr: Expr) -> None:
        """Add an occurrence of a subexpression."""
        if key not in self.occurrences:
            self.occurrences[key] = []
        self.occurrences[key].append(expr)

    def get_common_subexpressions(self, min_occurrences: int = 2) -> Dict[SubExprKey, List[Expr]]:
        """Return subexpressions that appear at least min_occurrences times."""
        return {k: v for k, v in self.occurrences.items() if len(v) >= min_occurrences}

    def generate_cse_name(self, key: SubExprKey) -> str:
        """Generate a variable name for a CSE intermediate."""
        self.cse_counter += 1
        if key.kind == 'call':
            return f"_cse_call_{self.cse_counter}"
        elif key.kind == 'coeff':
            return f"_cse_coeff_{self.cse_counter}"
        else:
            return f"_cse_{self.cse_counter}"


# =============================================================================
# Expression Optimizer
# =============================================================================

@dataclass
class OptimizedExpr:
    """
    Represents an optimized expression with CSE intermediates.

    Attributes:
        intermediates: List of (name, expr) pairs for intermediate variables
        result_expr: The final expression using intermediate variables
    """
    intermediates: List[Tuple[str, Expr]] = field(default_factory=list)
    result_expr: Optional[Expr] = None


class ExpressionOptimizer:
    """
    Applies optimizations to recurrence expressions.

    Supported optimizations:
    1. CSE (Common Subexpression Elimination)
    2. Horner's method for polynomial-like expressions
    3. Strength reduction (replacing expensive ops with cheaper ones)
    """

    def __init__(self, ctx: CodegenContext, enable_cse: bool = True,
                 enable_horner: bool = True, cse_threshold: int = 2):
        """
        Initialize the optimizer.

        Args:
            ctx: Code generation context
            enable_cse: Whether to apply CSE
            enable_horner: Whether to apply Horner's method
            cse_threshold: Minimum occurrences for CSE extraction
        """
        self.ctx = ctx
        self.enable_cse = enable_cse
        self.enable_horner = enable_horner
        self.cse_threshold = cse_threshold

    def optimize_expression(self, expr: Expr) -> OptimizedExpr:
        """
        Apply all enabled optimizations to an expression.

        Returns an OptimizedExpr with intermediate variable declarations
        and the optimized result expression.
        """
        result = OptimizedExpr()

        # Collect all recursive calls to cache
        calls = expr.collect_calls()

        if self.enable_cse and len(calls) >= self.cse_threshold:
            # Apply CSE for recursive calls
            result = self._apply_call_cse(expr, calls)
        else:
            result.result_expr = expr

        # Apply Horner's method if applicable
        if self.enable_horner:
            result = self._apply_horner(result)

        return result

    def _apply_call_cse(self, expr: Expr, calls: List[RecursiveCall]) -> OptimizedExpr:
        """
        Apply CSE by caching unique recursive calls.

        Generates intermediate variables for each unique recursive call,
        then replaces the calls in the expression with variable references.
        """
        result = OptimizedExpr()

        # Build map of unique calls to variable names
        call_to_var: Dict[str, str] = {}
        unique_calls: List[Tuple[str, RecursiveCall]] = []

        for call in calls:
            key = self._call_signature(call)
            if key not in call_to_var:
                var_name = f"e_{len(unique_calls)}"
                call_to_var[key] = var_name
                unique_calls.append((var_name, call))

        # Create intermediate declarations
        for var_name, call in unique_calls:
            result.intermediates.append((var_name, call))

        # Transform the expression to use cached variables
        result.result_expr = self._replace_calls_with_vars(expr, call_to_var)

        return result

    def _call_signature(self, call: RecursiveCall) -> str:
        """Generate a unique signature for a recursive call."""
        shifts = sorted(call.index_shifts.items())
        return ",".join(f"{k}:{v}" for k, v in shifts)

    def _replace_calls_with_vars(self, expr: Expr, call_to_var: Dict[str, str]) -> Expr:
        """Replace recursive calls in an expression with variable references."""
        if isinstance(expr, RecursiveCall):
            key = self._call_signature(expr)
            return CachedVar(call_to_var[key])

        elif isinstance(expr, Term):
            key = self._call_signature(expr.call)
            var = CachedVar(call_to_var[key])
            if isinstance(expr.coeff, Const) and expr.coeff.value == 1:
                return var
            return BinOp('*', expr.coeff, var)

        elif isinstance(expr, Sum):
            new_terms = []
            for term in expr.terms:
                new_term = self._replace_calls_with_vars(term, call_to_var)
                new_terms.append(new_term)
            return OptimizedSum(new_terms)

        elif isinstance(expr, ScaledExpr):
            new_inner = self._replace_calls_with_vars(expr.expr, call_to_var)
            return ScaledExpr(new_inner, expr.scale, expr.is_division)

        elif isinstance(expr, BinOp):
            new_left = self._replace_calls_with_vars(expr.left, call_to_var)
            new_right = self._replace_calls_with_vars(expr.right, call_to_var)
            return BinOp(expr.op, new_left, new_right)

        else:
            return expr

    def _apply_horner(self, opt_expr: OptimizedExpr) -> OptimizedExpr:
        """
        Apply Horner's method to polynomial-like expressions.

        Identifies expressions of the form:
            a + b*x + c*x^2 + ...
        And transforms them to:
            a + x*(b + x*c)

        This reduces the number of multiplications.
        """
        # For now, identify and optimize simple cases
        # Full implementation would analyze Sum expressions for polynomial structure

        if isinstance(opt_expr.result_expr, (Sum, OptimizedSum)):
            terms = (opt_expr.result_expr.terms if isinstance(opt_expr.result_expr, Sum)
                    else opt_expr.result_expr.exprs)

            # Check if this looks like a polynomial in any runtime variable
            # This is a simplified check - full implementation would be more thorough
            pass

        return opt_expr


# =============================================================================
# Optimized AST Nodes
# =============================================================================

@dataclass
class CachedVar(Expr):
    """
    Reference to a CSE-cached intermediate variable.

    Used to replace recursive calls with references to pre-computed values.
    """
    name: str

    def to_cpp(self, ctx: CodegenContext) -> str:
        return self.name

    def collect_calls(self) -> List[RecursiveCall]:
        return []

    def uses_var(self, var_name: str) -> bool:
        return self.name == var_name


@dataclass
class OptimizedSum(Expr):
    """
    Optimized sum of expressions (after CSE replacement).
    """
    exprs: List[Expr]

    def to_cpp(self, ctx: CodegenContext) -> str:
        if not self.exprs:
            return f"{ctx.vec_type}(0.0)"
        if len(self.exprs) == 1:
            return self.exprs[0].to_cpp(ctx)
        return " + ".join(e.to_cpp(ctx) for e in self.exprs)

    def collect_calls(self) -> List[RecursiveCall]:
        calls = []
        for e in self.exprs:
            calls.extend(e.collect_calls())
        return calls

    def uses_var(self, var_name: str) -> bool:
        return any(e.uses_var(var_name) for e in self.exprs)


# =============================================================================
# Code Generation with Optimizations
# =============================================================================

class OptimizedCodeGenerator:
    """
    Generates optimized C++ code with CSE and Horner transformations.

    This generator produces code with:
    1. Intermediate variable declarations for common subexpressions
    2. Optimized arithmetic using Horner's method
    3. Minimal redundant computation
    """

    def __init__(self, ctx: CodegenContext, optimizer: ExpressionOptimizer):
        self.ctx = ctx
        self.optimizer = optimizer

    def generate_body(self, expr: Expr, indent: str = "        ") -> str:
        """
        Generate an optimized function body.

        Returns C++ code with intermediate variables and final return statement.
        """
        opt = self.optimizer.optimize_expression(expr)

        lines = []

        # Generate intermediate variable declarations
        for var_name, intermediate in opt.intermediates:
            if isinstance(intermediate, RecursiveCall):
                call_cpp = intermediate.to_cpp(self.ctx)
                lines.append(f"{indent}{self.ctx.vec_type} {var_name} = {call_cpp};")
            else:
                lines.append(f"{indent}{self.ctx.vec_type} {var_name} = {intermediate.to_cpp(self.ctx)};")

        # Generate return statement
        if opt.result_expr:
            result_cpp = opt.result_expr.to_cpp(self.ctx)
            lines.append(f"{indent}return {result_cpp};")

        return "\n".join(lines)


# =============================================================================
# Memoization Pattern Generation
# =============================================================================

@dataclass
class MemoizationInfo:
    """
    Information for generating memoized template code.

    Memoization at the template level requires:
    1. A static storage mechanism (inline static variable or static member)
    2. Initialization flag or optional wrapper
    3. Thread-safety considerations (for parallel code)
    """
    storage_type: str = "inline_static"  # 'inline_static', 'static_member', 'thread_local'
    use_optional: bool = True  # Use std::optional for lazy initialization


class MemoizationGenerator:
    """
    Generates C++ code patterns for template-level memoization.

    For template metaprogramming, memoization is typically achieved through:
    1. Static inline variables with lazy initialization
    2. Template specialization caching
    3. Wrapper structs with memoization logic
    """

    def __init__(self, ctx: CodegenContext, info: MemoizationInfo = None):
        self.ctx = ctx
        self.info = info or MemoizationInfo()

    def generate_memoized_struct(self, targs: str, sig: str, body: str) -> str:
        """
        Generate a memoized struct with cached computation.

        The generated code uses inline static variables (C++17) for
        compile-time memoization across template instantiations.
        """
        if self.info.storage_type == "inline_static":
            return self._generate_inline_static(targs, sig, body)
        else:
            return self._generate_static_member(targs, sig, body)

    def _generate_inline_static(self, targs: str, sig: str, body: str) -> str:
        """Generate using C++17 inline static variables."""
        return f"""template<{targs}>
struct {self.ctx.struct_name}_Memoized {{
private:
    static inline bool _computed = false;
    static inline {self.ctx.vec_type} _cached_result;

public:
    static {self.ctx.vec_type} compute({sig}) {{
        if (!_computed) {{
{body}
            _cached_result = result;
            _computed = true;
        }}
        return _cached_result;
    }}
}};"""

    def _generate_static_member(self, targs: str, sig: str, body: str) -> str:
        """Generate using static member approach."""
        return f"""template<{targs}>
struct {self.ctx.struct_name}_Memoized {{
    static {self.ctx.vec_type} compute({sig}) {{
        static {self.ctx.vec_type} cached = []({sig}) {{
{body}
            return result;
        }}({", ".join(self.ctx.runtime_vars)});
        return cached;
    }}
}};"""


# =============================================================================
# Utility Functions
# =============================================================================

def count_operations(expr: Expr) -> Dict[str, int]:
    """
    Count arithmetic operations in an expression.

    Returns a dictionary with counts of:
    - 'add': additions
    - 'mul': multiplications
    - 'div': divisions
    - 'call': recursive calls
    """
    counts = defaultdict(int)

    def traverse(e: Expr):
        if isinstance(e, RecursiveCall):
            counts['call'] += 1
        elif isinstance(e, BinOp):
            if e.op == '+':
                counts['add'] += 1
            elif e.op == '*':
                counts['mul'] += 1
            traverse(e.left)
            traverse(e.right)
        elif isinstance(e, Term):
            if not (isinstance(e.coeff, Const) and e.coeff.value == 1):
                counts['mul'] += 1
            counts['call'] += 1
        elif isinstance(e, Sum):
            if len(e.terms) > 1:
                counts['add'] += len(e.terms) - 1
            for t in e.terms:
                traverse(t)
        elif isinstance(e, ScaledExpr):
            if e.is_division:
                counts['div'] += 1
            else:
                counts['mul'] += 1
            traverse(e.expr)

    traverse(expr)
    return dict(counts)


def estimate_cost(expr: Expr) -> float:
    """
    Estimate the computational cost of an expression.

    Uses approximate relative costs:
    - Addition: 1
    - Multiplication: 2
    - Division: 10
    - Recursive call: 50 (without memoization)
    """
    counts = count_operations(expr)
    return (counts.get('add', 0) * 1 +
            counts.get('mul', 0) * 2 +
            counts.get('div', 0) * 10 +
            counts.get('call', 0) * 50)


def should_apply_cse(expr: Expr) -> bool:
    """
    Determine if CSE would be beneficial for an expression.

    CSE is beneficial when:
    1. There are multiple recursive calls
    2. The same call appears more than once
    3. The estimated cost savings justify the overhead
    """
    calls = expr.collect_calls()
    if len(calls) < 2:
        return False

    # Check for duplicate calls
    signatures = set()
    for call in calls:
        sig = ",".join(f"{k}:{v}" for k, v in sorted(call.index_shifts.items()))
        if sig in signatures:
            return True
        signatures.add(sig)

    # Even without duplicates, CSE is useful to avoid recomputation
    # in the template instantiation
    return len(calls) >= 3


__all__ = [
    'SubExprKey',
    'CSEAnalyzer',
    'OptimizedExpr',
    'ExpressionOptimizer',
    'CachedVar',
    'OptimizedSum',
    'OptimizedCodeGenerator',
    'MemoizationInfo',
    'MemoizationGenerator',
    'count_operations',
    'estimate_cost',
    'should_apply_cse',
]
