"""
RECURSUM Code Generation Module

This package provides a DSL and code generator for creating C++ template
metaprogramming implementations of recurrence relations.

Key classes:
    - Recurrence: Main fluent API for defining recurrence relations
    - RecurrenceParser: Parses DSL expressions
    - CppGenerator: Generates C++ template code
    - ConstraintSet: Manages SFINAE constraints

Example:
    >>> from recursum.codegen import Recurrence
    >>> rec = Recurrence("Legendre", ["n"], ["x"])
    >>> rec.validity("n >= 0")
    >>> rec.base(n=0, value=1.0)
    >>> rec.base(n=1, value="x")
    >>> rec.rule("n > 1", "(2*n-1) * x * E[n-1] + (-(n-1)) * E[n-2]",
    ...          scale="1/n")
    >>> cpp_code = rec.generate()
"""

# Core AST and context
from .core import (
    Expr,
    Const,
    IndexExpr,
    Var,
    RecursiveCall,
    BinOp,
    Term,
    Sum,
    ScaledExpr,
    CodegenContext,
)

# Constraints
from .constraints import (
    ConstraintOp,
    Constraint,
    ConstraintSet,
)

# Parser
from .parser import RecurrenceParser

# Main API
from .recurrence import (
    BaseCase,
    RecurrenceRule,
    Recurrence,
)

# C++ Generators
from .cpp_generator import CppGenerator
from .layered_generator import LayeredCppGenerator

# Optimizer
from .optimizer import (
    ExpressionOptimizer,
    OptimizedExpr,
    CachedVar,
    OptimizedSum,
    CSEAnalyzer,
    MemoizationGenerator,
    MemoizationInfo,
    count_operations,
    estimate_cost,
    should_apply_cse,
)

__all__ = [
    # Core AST
    "Expr",
    "Const",
    "IndexExpr",
    "Var",
    "RecursiveCall",
    "BinOp",
    "Term",
    "Sum",
    "ScaledExpr",
    "CodegenContext",
    # Constraints
    "ConstraintOp",
    "Constraint",
    "ConstraintSet",
    # Parser
    "RecurrenceParser",
    # Main API
    "BaseCase",
    "RecurrenceRule",
    "Recurrence",
    # Generators
    "CppGenerator",
    "LayeredCppGenerator",
    # Optimizer
    "ExpressionOptimizer",
    "OptimizedExpr",
    "CachedVar",
    "OptimizedSum",
    "CSEAnalyzer",
    "MemoizationGenerator",
    "MemoizationInfo",
    "count_operations",
    "estimate_cost",
    "should_apply_cse",
]
