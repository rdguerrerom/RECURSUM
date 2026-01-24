#!/usr/bin/env python3
"""
Comprehensive unit tests for recursum.codegen.constraints module.

Tests constraint parsing, SFINAE generation, and constraint set operations.
"""

import pytest
from recursum.codegen.constraints import Constraint, ConstraintOp, ConstraintSet


class TestConstraintOp:
    """Test ConstraintOp enumeration."""

    def test_all_operators_defined(self):
        """Test that all standard comparison operators are defined."""
        assert ConstraintOp.EQ.value == "=="
        assert ConstraintOp.NE.value == "!="
        assert ConstraintOp.LT.value == "<"
        assert ConstraintOp.LE.value == "<="
        assert ConstraintOp.GT.value == ">"
        assert ConstraintOp.GE.value == ">="


class TestConstraint:
    """Test Constraint class."""

    def test_constraint_creation(self):
        """Test direct constraint creation."""
        c = Constraint("n", ConstraintOp.GT, "0")
        assert c.left == "n"
        assert c.op == ConstraintOp.GT
        assert c.right == "0"

    def test_to_sfinae_equal(self):
        """Test SFINAE generation for equality."""
        c = Constraint("n", ConstraintOp.EQ, "0")
        assert c.to_sfinae() == "(n == 0)"

    def test_to_sfinae_not_equal(self):
        """Test SFINAE generation for inequality."""
        c = Constraint("n", ConstraintOp.NE, "0")
        assert c.to_sfinae() == "(n != 0)"

    def test_to_sfinae_greater_than(self):
        """Test SFINAE generation for greater than."""
        c = Constraint("n", ConstraintOp.GT, "1")
        assert c.to_sfinae() == "(n > 1)"

    def test_to_sfinae_greater_equal(self):
        """Test SFINAE generation for greater or equal."""
        c = Constraint("n", ConstraintOp.GE, "0")
        assert c.to_sfinae() == "(n >= 0)"

    def test_to_sfinae_less_than(self):
        """Test SFINAE generation for less than."""
        c = Constraint("n", ConstraintOp.LT, "10")
        assert c.to_sfinae() == "(n < 10)"

    def test_to_sfinae_less_equal(self):
        """Test SFINAE generation for less or equal."""
        c = Constraint("m", ConstraintOp.LE, "n")
        assert c.to_sfinae() == "(m <= n)"

    def test_parse_simple_gt(self):
        """Test parsing simple > constraint."""
        c = Constraint.parse("n > 0")
        assert c.left == "n"
        assert c.op == ConstraintOp.GT
        assert c.right == "0"

    def test_parse_simple_ge(self):
        """Test parsing simple >= constraint."""
        c = Constraint.parse("n >= 0")
        assert c.left == "n"
        assert c.op == ConstraintOp.GE
        assert c.right == "0"

    def test_parse_simple_lt(self):
        """Test parsing simple < constraint."""
        c = Constraint.parse("n < 10")
        assert c.left == "n"
        assert c.op == ConstraintOp.LT
        assert c.right == "10"

    def test_parse_simple_le(self):
        """Test parsing simple <= constraint."""
        c = Constraint.parse("m <= n")
        assert c.left == "m"
        assert c.op == ConstraintOp.LE
        assert c.right == "n"

    def test_parse_simple_eq(self):
        """Test parsing simple == constraint."""
        c = Constraint.parse("i == 0")
        assert c.left == "i"
        assert c.op == ConstraintOp.EQ
        assert c.right == "0"

    def test_parse_simple_ne(self):
        """Test parsing simple != constraint."""
        c = Constraint.parse("j != 0")
        assert c.left == "j"
        assert c.op == ConstraintOp.NE
        assert c.right == "0"

    def test_parse_with_spaces(self):
        """Test parsing with various whitespace."""
        c = Constraint.parse("  n   >   0  ")
        assert c.left == "n"
        assert c.op == ConstraintOp.GT
        assert c.right == "0"

    def test_parse_expression_left(self):
        """Test parsing with expression on left side."""
        c = Constraint.parse("n + m >= 0")
        assert c.left == "n + m"
        assert c.op == ConstraintOp.GE
        assert c.right == "0"

    def test_parse_expression_right(self):
        """Test parsing with expression on right side."""
        c = Constraint.parse("n > m + 1")
        assert c.left == "n"
        assert c.op == ConstraintOp.GT
        assert c.right == "m + 1"

    def test_parse_expression_both_sides(self):
        """Test parsing with expressions on both sides."""
        c = Constraint.parse("n + m >= i + j")
        assert c.left == "n + m"
        assert c.op == ConstraintOp.GE
        assert c.right == "i + j"

    def test_parse_invalid_no_operator(self):
        """Test parsing invalid constraint without operator."""
        with pytest.raises(ValueError, match="Cannot parse constraint"):
            Constraint.parse("n 0")

    def test_parse_invalid_empty(self):
        """Test parsing empty constraint."""
        with pytest.raises(ValueError, match="Cannot parse constraint"):
            Constraint.parse("")


class TestConstraintSet:
    """Test ConstraintSet class."""

    def test_empty_constraint_set(self):
        """Test empty constraint set."""
        cs = ConstraintSet([])
        assert cs.constraints == []
        assert cs.to_sfinae() == "true"

    def test_single_constraint(self):
        """Test constraint set with single constraint."""
        c = Constraint("n", ConstraintOp.GT, "0")
        cs = ConstraintSet([c])
        assert len(cs.constraints) == 1
        assert cs.to_sfinae() == "(n > 0)"

    def test_multiple_constraints(self):
        """Test constraint set with multiple constraints."""
        c1 = Constraint("n", ConstraintOp.GE, "0")
        c2 = Constraint("m", ConstraintOp.GE, "0")
        c3 = Constraint("n", ConstraintOp.GE, "m")
        cs = ConstraintSet([c1, c2, c3])
        assert len(cs.constraints) == 3
        assert cs.to_sfinae() == "(n >= 0) && (m >= 0) && (n >= m)"

    def test_parse_single_constraint(self):
        """Test parsing single constraint string."""
        cs = ConstraintSet.parse("n > 0")
        assert len(cs.constraints) == 1
        assert cs.constraints[0].left == "n"
        assert cs.constraints[0].op == ConstraintOp.GT
        assert cs.constraints[0].right == "0"

    def test_parse_multiple_constraints(self):
        """Test parsing multiple constraint strings."""
        cs = ConstraintSet.parse("n >= 0", "m >= 0", "n >= m")
        assert len(cs.constraints) == 3
        assert cs.to_sfinae() == "(n >= 0) && (m >= 0) && (n >= m)"

    def test_parse_combined_constraints(self):
        """Test parsing combined constraint string with &&."""
        cs = ConstraintSet.parse("n > 0 && m > 0")
        assert len(cs.constraints) == 2
        assert cs.constraints[0].left == "n"
        assert cs.constraints[0].op == ConstraintOp.GT
        assert cs.constraints[1].left == "m"
        assert cs.constraints[1].op == ConstraintOp.GT

    def test_parse_mixed_format(self):
        """Test parsing mix of combined and separate constraints."""
        cs = ConstraintSet.parse("n > 0 && m > 0", "n >= m")
        assert len(cs.constraints) == 3
        assert cs.constraints[0].left == "n"
        assert cs.constraints[1].left == "m"
        assert cs.constraints[2].left == "n"
        assert cs.constraints[2].right == "m"

    def test_parse_empty_string(self):
        """Test parsing empty string."""
        cs = ConstraintSet.parse("")
        assert len(cs.constraints) == 0
        assert cs.to_sfinae() == "true"

    def test_parse_whitespace_only(self):
        """Test parsing whitespace-only string."""
        cs = ConstraintSet.parse("   ")
        assert len(cs.constraints) == 0

    def test_merge_empty_sets(self):
        """Test merging two empty constraint sets."""
        cs1 = ConstraintSet([])
        cs2 = ConstraintSet([])
        merged = cs1.merge(cs2)
        assert len(merged.constraints) == 0

    def test_merge_with_empty(self):
        """Test merging non-empty set with empty set."""
        c = Constraint("n", ConstraintOp.GT, "0")
        cs1 = ConstraintSet([c])
        cs2 = ConstraintSet([])
        merged = cs1.merge(cs2)
        assert len(merged.constraints) == 1
        assert merged.constraints[0] == c

    def test_merge_two_sets(self):
        """Test merging two non-empty constraint sets."""
        c1 = Constraint("n", ConstraintOp.GE, "0")
        c2 = Constraint("m", ConstraintOp.GE, "0")
        cs1 = ConstraintSet([c1])
        cs2 = ConstraintSet([c2])
        merged = cs1.merge(cs2)
        assert len(merged.constraints) == 2
        assert merged.constraints[0] == c1
        assert merged.constraints[1] == c2

    def test_merge_preserves_order(self):
        """Test that merge preserves constraint order."""
        c1 = Constraint("n", ConstraintOp.GT, "0")
        c2 = Constraint("m", ConstraintOp.GT, "0")
        c3 = Constraint("i", ConstraintOp.EQ, "0")
        cs1 = ConstraintSet([c1, c2])
        cs2 = ConstraintSet([c3])
        merged = cs1.merge(cs2)
        assert merged.constraints[0] == c1
        assert merged.constraints[1] == c2
        assert merged.constraints[2] == c3


class TestConstraintIntegration:
    """Integration tests for constraint system."""

    def test_typical_recurrence_validity(self):
        """Test typical validity constraint for recurrence."""
        cs = ConstraintSet.parse("n >= 0")
        sfinae = cs.to_sfinae()
        assert sfinae == "(n >= 0)"

    def test_two_index_recurrence(self):
        """Test two-index recurrence constraints."""
        cs = ConstraintSet.parse("i >= 0", "j >= 0", "i + j >= t")
        sfinae = cs.to_sfinae()
        assert "(i >= 0)" in sfinae
        assert "(j >= 0)" in sfinae
        assert "(i + j >= t)" in sfinae

    def test_associated_legendre_constraints(self):
        """Test constraints for associated Legendre polynomials."""
        cs = ConstraintSet.parse("n >= 0", "m >= 0", "n >= m")
        sfinae = cs.to_sfinae()
        assert len(cs.constraints) == 3

    def test_equality_constraint_for_base_case(self):
        """Test equality constraint for base case."""
        cs = ConstraintSet.parse("n == 0")
        sfinae = cs.to_sfinae()
        assert sfinae == "(n == 0)"

    def test_compound_condition(self):
        """Test compound condition for recurrence rule."""
        cs = ConstraintSet.parse("n > 1 && n < 10")
        assert len(cs.constraints) == 2
        sfinae = cs.to_sfinae()
        assert "(n > 1)" in sfinae
        assert "(n < 10)" in sfinae


class TestConstraintEdgeCases:
    """Test edge cases and error handling."""

    def test_constraint_with_parentheses(self):
        """Test constraint with parentheses in expression."""
        c = Constraint.parse("(n+1) > 0")
        assert c.left == "(n+1)"
        assert c.op == ConstraintOp.GT
        assert c.right == "0"

    def test_constraint_with_multiplication(self):
        """Test constraint with multiplication."""
        c = Constraint.parse("2*n >= m")
        assert c.left == "2*n"
        assert c.op == ConstraintOp.GE
        assert c.right == "m"

    def test_multiple_equals_signs(self):
        """Test that == is parsed correctly (not as two =)."""
        c = Constraint.parse("n == 0")
        assert c.op == ConstraintOp.EQ
        # Should not parse as two separate = operators

    def test_constraint_with_negative_number(self):
        """Test constraint with negative number."""
        c = Constraint.parse("n > -1")
        assert c.left == "n"
        assert c.op == ConstraintOp.GT
        assert c.right == "-1"

    def test_operator_precedence_in_parsing(self):
        """Test that >= is parsed before > (longer operator first)."""
        c = Constraint.parse("n >= 0")
        assert c.op == ConstraintOp.GE
        # Should not parse as > with "=" left over

    def test_ne_parsed_correctly(self):
        """Test that != is parsed correctly."""
        c = Constraint.parse("n != 0")
        assert c.op == ConstraintOp.NE

    def test_le_parsed_correctly(self):
        """Test that <= is parsed correctly."""
        c = Constraint.parse("n <= 10")
        assert c.op == ConstraintOp.LE
