"""
Parser for recurrence relation DSL.

Parses einsum-inspired notation for recurrence relations:
    E[i,j,t]           - Recurrence at indices i,j,t
    E[i-1,j,t+1]       - With shifts
    coeff * E[...]     - Scaled term
    term1 + term2      - Sum of terms
"""

import re
from typing import Dict, List
from .core import (
    Expr, Const, IndexExpr, Var, RecursiveCall,
    BinOp, Term, Sum
)


class RecurrenceParser:
    """Parser for recurrence relation expressions."""

    def __init__(self, indices: List[str], runtime_vars: List[str]):
        """
        Initialize parser.

        Args:
            indices: Template parameter names (compile-time)
            runtime_vars: Runtime variable names
        """
        self.indices = indices
        self.runtime_vars = runtime_vars

    def parse_coefficient(self, s: str) -> Expr:
        """Parse a coefficient (constant, variable, or index expression)."""
        s = s.strip()
        if not s or s == "1":
            return Const(1)

        # Handle compound: (index_expr) * runtime_var  e.g., (2*n-1) * x
        if '*' in s and not s.startswith('('):
            parts = s.split('*')
            if len(parts) == 2:
                left = self.parse_coefficient(parts[0].strip())
                right = self.parse_coefficient(parts[1].strip())
                return BinOp('*', left, right)

        if s.startswith('(') and s.endswith(')'):
            inner = s[1:-1].strip()
            # Check if inner contains indices or variables
            has_idx = any(idx in inner for idx in self.indices)
            has_var = any(v in inner for v in self.runtime_vars)

            if has_idx and has_var:
                # Complex: treat as compound index expression
                return IndexExpr(inner)
            elif has_idx:
                return IndexExpr(inner)
            elif has_var:
                return Var(inner)
            try:
                return Const(eval(inner))
            except:
                return IndexExpr(inner)

        # Check if it's a runtime variable
        if s in self.runtime_vars:
            return Var(s)

        # Check if it's an index
        if s in self.indices:
            return IndexExpr(s)

        # Try to parse as numeric constant
        try:
            v = float(s)
            return Const(int(v) if v == int(v) else v)
        except:
            pass

        # Check if it contains any index (likely an index expression)
        for idx in self.indices:
            if idx in s:
                return IndexExpr(s)

        # Default: treat as variable
        return Var(s)

    def parse_index_shift(self, s: str) -> Dict[str, int]:
        """
        Parse index shifts from E[...] notation.

        Args:
            s: Index expression like "i-1, j, t+1"

        Returns:
            Dictionary mapping index names to shift amounts
        """
        shifts = {idx: 0 for idx in self.indices}
        parts = [p.strip() for p in s.split(',')]

        for i, part in enumerate(parts):
            if i >= len(self.indices):
                break
            idx = self.indices[i]

            if not part or part == idx:
                continue

            # Try to match patterns like "n-1" or "n+2"
            for known in self.indices:
                m = re.match(rf'^({known})\s*([+-])\s*(\d+)$', part)
                if m:
                    sign = 1 if m.group(2) == '+' else -1
                    shifts[m.group(1)] = sign * int(m.group(3))
                    break

        return shifts

    def parse_term(self, s: str) -> Term:
        """
        Parse a single term: coefficient * E[...].

        Args:
            s: Term string like "(2*n-1) * x * E[n-1]"

        Returns:
            Term AST node
        """
        s = s.strip()

        # Find the E[...] part
        m = re.search(r'E\[([^\]]+)\]', s)
        if not m:
            raise ValueError(f"No E[...] found in term: {s}")

        # Parse index shifts
        shifts = self.parse_index_shift(m.group(1))
        call = RecursiveCall(shifts)

        # Extract coefficient (everything before E[...])
        coeff_part = s[:m.start()].strip()
        if coeff_part.endswith('*'):
            coeff_part = coeff_part[:-1].strip()

        if not coeff_part or coeff_part == '1':
            return Term(Const(1), call)

        # Handle chained multiplication: (2*n-1) * x * E[...]
        # Split by * but respect parentheses
        parts = self._split_by_mult(coeff_part)
        if len(parts) == 1:
            coeff = self.parse_coefficient(parts[0])
        else:
            coeff = self.parse_coefficient(parts[0])
            for p in parts[1:]:
                coeff = BinOp('*', coeff, self.parse_coefficient(p))

        return Term(coeff, call)

    def _split_by_mult(self, s: str) -> List[str]:
        """Split by * respecting parentheses."""
        parts = []
        current = ""
        depth = 0

        for c in s:
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
            elif c == '*' and depth == 0:
                if current.strip():
                    parts.append(current.strip())
                current = ""
                continue
            current += c

        if current.strip():
            parts.append(current.strip())

        return parts if parts else ["1"]

    def parse_expression(self, s: str) -> Expr:
        """
        Parse a full expression (sum of terms).

        Args:
            s: Expression like "coeff1 * E[i-1] + coeff2 * E[i-2]"

        Returns:
            Expression AST (typically a Sum)
        """
        terms = []
        current = ""
        depth = 0

        for c in s:
            if c in '([':
                depth += 1
            elif c in ')]':
                depth -= 1
            elif c == '+' and depth == 0:
                if current.strip():
                    terms.append(self.parse_term(current))
                current = ""
                continue
            current += c

        if current.strip():
            terms.append(self.parse_term(current))

        return Sum(terms)

    def parse_scale(self, s: str) -> Expr:
        """
        Parse a scaling factor (typically from 1/n notation).

        Args:
            s: Scale expression like "1/n" or "1/(2*n)"

        Returns:
            Expression for the scale factor
        """
        s = s.strip()
        if s.startswith("1/"):
            d = s[2:].strip()
            if d.startswith('(') and d.endswith(')'):
                d = d[1:-1]

            # Check if it's an index expression
            for idx in self.indices:
                if idx in d:
                    return IndexExpr(d)

            # Try to parse as constant
            try:
                return Const(float(d))
            except:
                return Var(d)

        return self.parse_coefficient(s)
