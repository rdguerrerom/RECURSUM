"""
RECURSUM: Recurrence Relation Code Generator

A Python DSL and C++ code generator for zero-overhead recurrence relation
evaluation using template metaprogramming.

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

__version__ = "0.1.0"

# Main API
from .codegen import Recurrence

__all__ = ["Recurrence"]
