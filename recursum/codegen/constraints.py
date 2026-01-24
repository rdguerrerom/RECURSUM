"""
Constraint system for SFINAE-based template specialization.

Constraints define when a particular recurrence rule applies,
and are compiled to C++ std::enable_if conditions.
"""

from dataclasses import dataclass
from typing import List
from enum import Enum


class ConstraintOp(Enum):
    """Comparison operators for constraints."""
    EQ = "=="
    NE = "!="
    LT = "<"
    LE = "<="
    GT = ">"
    GE = ">="


@dataclass
class Constraint:
    """Single constraint (e.g., n > 0)."""
    left: str
    op: ConstraintOp
    right: str

    def to_sfinae(self) -> str:
        """Convert to SFINAE condition for std::enable_if."""
        return f"({self.left} {self.op.value} {self.right})"

    @classmethod
    def parse(cls, expr: str) -> "Constraint":
        """Parse constraint from string (e.g., 'n > 0')."""
        # Try operators in order of precedence (longer first to avoid substring matches)
        for op in [ConstraintOp.EQ, ConstraintOp.NE, ConstraintOp.GE,
                   ConstraintOp.LE, ConstraintOp.GT, ConstraintOp.LT]:
            if op.value in expr:
                parts = expr.split(op.value)
                if len(parts) == 2:
                    return cls(parts[0].strip(), op, parts[1].strip())
        raise ValueError(f"Cannot parse constraint: {expr}")


@dataclass
class ConstraintSet:
    """Set of constraints combined with AND logic."""
    constraints: List[Constraint]

    def to_sfinae(self) -> str:
        """Convert to SFINAE condition."""
        if not self.constraints:
            return "true"
        return " && ".join(c.to_sfinae() for c in self.constraints)

    @classmethod
    def parse(cls, *exprs: str) -> "ConstraintSet":
        """Parse multiple constraint expressions.

        Args:
            *exprs: Constraint expressions like "n > 0", "m >= 0 && n >= m"

        Returns:
            ConstraintSet with all parsed constraints
        """
        all_constraints = []
        for expr in exprs:
            # Split by && to handle combined constraints
            for part in expr.split("&&"):
                part = part.strip()
                if part:
                    all_constraints.append(Constraint.parse(part))
        return cls(all_constraints)

    def merge(self, other: "ConstraintSet") -> "ConstraintSet":
        """Merge with another constraint set (AND logic)."""
        return ConstraintSet(self.constraints + other.constraints)
