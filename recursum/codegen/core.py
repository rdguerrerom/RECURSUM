"""
Core AST (Abstract Syntax Tree) classes for recurrence relation expressions.

These classes represent the structure of recurrence relation expressions
in a form that can be compiled to C++ template code.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Union
from abc import ABC, abstractmethod


# =============================================================================
# AST Expression Nodes
# =============================================================================

class Expr(ABC):
    """Base class for all expression nodes in the AST."""

    @abstractmethod
    def to_cpp(self, ctx: "CodegenContext") -> str:
        """Convert this expression to C++ code."""
        pass

    @abstractmethod
    def collect_calls(self) -> List["RecursiveCall"]:
        """Collect all recursive calls in this expression."""
        pass

    def uses_var(self, var_name: str) -> bool:
        """Check if this expression uses a given variable name."""
        return False

    def __add__(self, other) -> "Expr":
        return BinOp('+', self, other if isinstance(other, Expr) else Const(other))

    def __mul__(self, other) -> "Expr":
        return BinOp('*', self, other if isinstance(other, Expr) else Const(other))


@dataclass
class Const(Expr):
    """Constant value (numeric or string literal)."""
    value: Union[str, int, float]

    def to_cpp(self, ctx: "CodegenContext") -> str:
        if isinstance(self.value, (int, float)):
            return f"{ctx.vec_type}({self.value})"
        return str(self.value)

    def collect_calls(self) -> List["RecursiveCall"]:
        return []

    def uses_var(self, var_name: str) -> bool:
        return False


@dataclass
class IndexExpr(Expr):
    """Expression involving template indices (e.g., 2*n-1)."""
    expr_str: str

    def to_cpp(self, ctx: "CodegenContext") -> str:
        return f"{ctx.vec_type}({self.expr_str})"

    def collect_calls(self) -> List["RecursiveCall"]:
        return []

    def uses_var(self, var_name: str) -> bool:
        return var_name in self.expr_str


@dataclass
class Var(Expr):
    """Runtime variable reference (e.g., PA, x)."""
    name: str

    def to_cpp(self, ctx: "CodegenContext") -> str:
        return self.name

    def collect_calls(self) -> List["RecursiveCall"]:
        return []

    def uses_var(self, var_name: str) -> bool:
        return self.name == var_name


@dataclass
class RecursiveCall(Expr):
    """Recursive call to the recurrence with index shifts (e.g., E[i-1, j+1])."""
    index_shifts: Dict[str, int]

    def to_cpp(self, ctx: "CodegenContext") -> str:
        args = []
        for idx in ctx.indices:
            shift = self.index_shifts.get(idx, 0)
            if shift == 0:
                args.append(idx)
            elif shift > 0:
                args.append(f"{idx} + {shift}")
            else:
                args.append(f"{idx} - {-shift}")

        template_args = ", ".join(args)
        runtime_args = ", ".join(ctx.runtime_vars)
        return f"{ctx.struct_name}<{template_args}>::compute({runtime_args})"

    def collect_calls(self) -> List["RecursiveCall"]:
        return [self]

    def uses_var(self, var_name: str) -> bool:
        return False


@dataclass
class BinOp(Expr):
    """Binary operation (+ or *)."""
    op: str  # '+' or '*'
    left: Expr
    right: Expr

    def to_cpp(self, ctx: "CodegenContext") -> str:
        left_cpp = self.left.to_cpp(ctx)
        right_cpp = self.right.to_cpp(ctx)

        # Add parentheses for clarity
        if isinstance(self.left, BinOp):
            left_cpp = f"({left_cpp})"
        if isinstance(self.right, BinOp):
            right_cpp = f"({right_cpp})"

        return f"{left_cpp} {self.op} {right_cpp}"

    def collect_calls(self) -> List["RecursiveCall"]:
        return self.left.collect_calls() + self.right.collect_calls()

    def uses_var(self, var_name: str) -> bool:
        return self.left.uses_var(var_name) or self.right.uses_var(var_name)


@dataclass
class Term(Expr):
    """A term in a sum: coefficient * recursive_call."""
    coeff: Expr
    call: RecursiveCall

    def to_cpp(self, ctx: "CodegenContext") -> str:
        coeff_cpp = self.coeff.to_cpp(ctx)
        call_cpp = self.call.to_cpp(ctx)

        if isinstance(self.coeff, Const) and self.coeff.value == 1:
            return call_cpp
        return f"{coeff_cpp} * {call_cpp}"

    def collect_calls(self) -> List["RecursiveCall"]:
        return [self.call]

    def uses_var(self, var_name: str) -> bool:
        return self.coeff.uses_var(var_name)


@dataclass
class Sum(Expr):
    """Sum of terms."""
    terms: List[Term]

    def to_cpp(self, ctx: "CodegenContext") -> str:
        if not self.terms:
            return f"{ctx.vec_type}(0.0)"
        if len(self.terms) == 1:
            return self.terms[0].to_cpp(ctx)
        return " + ".join(t.to_cpp(ctx) for t in self.terms)

    def collect_calls(self) -> List["RecursiveCall"]:
        calls = []
        for term in self.terms:
            calls.extend(term.collect_calls())
        return calls

    def uses_var(self, var_name: str) -> bool:
        return any(t.uses_var(var_name) for t in self.terms)


@dataclass
class ScaledExpr(Expr):
    """Expression scaled by a factor: expr / scale or expr * scale."""
    expr: Expr
    scale: Expr
    is_division: bool = True

    def to_cpp(self, ctx: "CodegenContext") -> str:
        expr_cpp = self.expr.to_cpp(ctx)
        scale_cpp = self.scale.to_cpp(ctx)
        if self.is_division:
            return f"({expr_cpp}) / ({scale_cpp})"
        else:
            return f"({expr_cpp}) * ({scale_cpp})"

    def collect_calls(self) -> List["RecursiveCall"]:
        return self.expr.collect_calls()

    def uses_var(self, var_name: str) -> bool:
        return self.expr.uses_var(var_name) or self.scale.uses_var(var_name)


@dataclass
class CodegenContext:
    """Context for code generation containing naming and type information."""
    struct_name: str
    indices: List[str]
    runtime_vars: List[str]
    vec_type: str = "Vec8d"
