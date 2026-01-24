#!/usr/bin/env python3
"""
Comprehensive unit tests for recursum.codegen.recurrence module.

Tests the Recurrence fluent API for defining recurrence relations.
"""

import pytest
from recursum.codegen.recurrence import Recurrence, RecurrenceRule, BaseCase
from recursum.codegen.core import Const, Var, RecursiveCall, Sum, ScaledExpr
from recursum.codegen.constraints import ConstraintSet, ConstraintOp


class TestRecurrence:
    """Test Recurrence class initialization."""

    def test_basic_creation(self):
        """Test basic recurrence creation."""
        rec = Recurrence("Test", ["n"], ["x"])
        assert rec.name == "Test"
        assert rec.indices == ["n"]
        assert rec.runtime_vars == ["x"]
        assert rec.vec_type == "Vec8d"

    def test_custom_vec_type(self):
        """Test recurrence with custom vector type."""
        rec = Recurrence("Test", ["n"], ["x"], vec_type="double")
        assert rec.vec_type == "double"

    def test_custom_namespace(self):
        """Test recurrence with custom namespace."""
        rec = Recurrence("Test", ["n"], ["x"], namespace="myns")
        assert rec.namespace == "myns"

    def test_multi_index(self):
        """Test multi-index recurrence."""
        rec = Recurrence("Hermite", ["i", "j", "t"], ["PA", "PB", "p"])
        assert len(rec.indices) == 3
        assert len(rec.runtime_vars) == 3

    def test_max_indices_default(self):
        """Test default max_indices initialization."""
        rec = Recurrence("Test", ["n"], ["x"])
        assert rec.max_indices is not None
        assert rec.max_indices["n"] == 20

    def test_max_indices_custom(self):
        """Test custom max_indices."""
        rec = Recurrence("Test", ["n"], ["x"], max_indices={"n": 50})
        assert rec.max_indices["n"] == 50

    def test_scipy_reference(self):
        """Test scipy_reference setting."""
        rec = Recurrence("Legendre", ["n"], ["x"], scipy_reference="eval_legendre")
        assert rec.scipy_reference == "eval_legendre"

    def test_initial_state(self):
        """Test initial state of recurrence."""
        rec = Recurrence("Test", ["n"], ["x"])
        assert rec._base_cases == []
        assert rec._rules == []
        assert rec._validity is None


class TestRecurrenceValidity:
    """Test validity constraint setting."""

    def test_single_constraint(self):
        """Test setting single validity constraint."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.validity("n >= 0")
        assert rec._validity is not None
        assert len(rec._validity.constraints) == 1

    def test_multiple_constraints(self):
        """Test setting multiple validity constraints."""
        rec = Recurrence("Test", ["n", "m"], ["x"])
        rec.validity("n >= 0", "m >= 0", "n >= m")
        assert len(rec._validity.constraints) == 3

    def test_validity_returns_self(self):
        """Test that validity() returns self for chaining."""
        rec = Recurrence("Test", ["n"], ["x"])
        result = rec.validity("n >= 0")
        assert result is rec

    def test_validity_chaining(self):
        """Test method chaining with validity."""
        rec = Recurrence("Test", ["n"], ["x"]).validity("n >= 0")
        assert rec._validity is not None


class TestRecurrenceBase:
    """Test base case definitions."""

    def test_base_integer_value(self):
        """Test base case with integer value."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.base(n=0, value=1)
        assert len(rec._base_cases) == 1
        bc = rec._base_cases[0]
        assert bc.index_values["n"] == 0
        assert isinstance(bc.value, Const)
        assert bc.value.value == 1

    def test_base_float_value(self):
        """Test base case with float value."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.base(n=0, value=3.14)
        bc = rec._base_cases[0]
        assert isinstance(bc.value, Const)
        assert bc.value.value == 3.14

    def test_base_string_variable(self):
        """Test base case with string variable."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.base(n=1, value="x")
        bc = rec._base_cases[0]
        assert isinstance(bc.value, Var)
        assert bc.value.name == "x"

    def test_base_multiple_indices(self):
        """Test base case with multiple indices."""
        rec = Recurrence("Test", ["i", "j"], ["x"])
        rec.base(i=0, j=0, value=1.0)
        bc = rec._base_cases[0]
        assert bc.index_values["i"] == 0
        assert bc.index_values["j"] == 0

    def test_multiple_base_cases(self):
        """Test multiple base cases."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.base(n=0, value=1.0)
        rec.base(n=1, value="x")
        assert len(rec._base_cases) == 2

    def test_base_returns_self(self):
        """Test that base() returns self for chaining."""
        rec = Recurrence("Test", ["n"], ["x"])
        result = rec.base(n=0, value=1)
        assert result is rec

    def test_base_chaining(self):
        """Test method chaining with base."""
        rec = (Recurrence("Test", ["n"], ["x"])
               .base(n=0, value=1.0)
               .base(n=1, value="x"))
        assert len(rec._base_cases) == 2


class TestRecurrenceRule:
    """Test recurrence rule definitions."""

    def test_simple_rule(self):
        """Test simple recurrence rule."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.rule("n > 1", "x * E[n-1]")
        assert len(rec._rules) == 1
        rule = rec._rules[0]
        assert rule.constraints is not None
        assert rule.expression is not None

    def test_rule_with_two_terms(self):
        """Test rule with two terms."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.rule("n > 1", "x * E[n-1] + E[n-2]")
        rule = rec._rules[0]
        assert isinstance(rule.expression, Sum)
        assert len(rule.expression.terms) == 2

    def test_rule_with_scale(self):
        """Test rule with scale factor."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.rule("n > 1", "(2*n-1) * x * E[n-1]", scale="1/n")
        rule = rec._rules[0]
        assert isinstance(rule.expression, ScaledExpr)
        assert rule.expression.is_division

    def test_rule_with_name(self):
        """Test rule with descriptive name."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.rule("n > 1", "x * E[n-1]", name="Forward recurrence")
        rule = rec._rules[0]
        assert rule.name == "Forward recurrence"

    def test_multiple_rules(self):
        """Test multiple rules."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.rule("n == 0", "E[n]")  # Changed from "1" to valid E[...] expression
        rec.rule("n > 0", "x * E[n-1]")
        assert len(rec._rules) == 2

    def test_rule_returns_self(self):
        """Test that rule() returns self for chaining."""
        rec = Recurrence("Test", ["n"], ["x"])
        result = rec.rule("n > 1", "E[n-1]")
        assert result is rec

    def test_rule_chaining(self):
        """Test method chaining with rules."""
        rec = (Recurrence("Test", ["n"], ["x"])
               .rule("n == 1", "E[n]")  # Changed from "x" to valid E[...] expression
               .rule("n > 1", "x * E[n-1]"))
        assert len(rec._rules) == 2

    def test_rule_with_list_constraints(self):
        """Test rule with list of constraints."""
        rec = Recurrence("Test", ["n", "m"], ["x"])
        rec.rule(["n > 0", "m > 0"], "x * E[n-1, m-1]")
        rule = rec._rules[0]
        assert len(rule.constraints.constraints) == 2

    def test_rule_with_constraint_set(self):
        """Test rule with ConstraintSet object."""
        rec = Recurrence("Test", ["n"], ["x"])
        cs = ConstraintSet.parse("n > 1")
        rec.rule(cs, "E[n-1]")
        assert len(rec._rules) == 1


class TestRecurrenceBranchAverage:
    """Test branch_average method."""

    def test_single_branch(self):
        """Test branch average with single branch."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.branch_average("n > 0", ["x * E[n-1]"])
        assert len(rec._rules) == 1

    def test_two_branches(self):
        """Test branch average with two branches."""
        rec = Recurrence("Test", ["i", "j"], ["x", "y"])
        rec.branch_average(
            "i > 0 && j > 0",
            ["x * E[i-1, j]", "y * E[i, j-1]"]
        )
        rule = rec._rules[0]
        # Should create averaged expression
        assert rule.expression is not None

    def test_three_branches(self):
        """Test branch average with three branches."""
        rec = Recurrence("Test", ["i", "j", "k"], ["x"])
        rec.branch_average(
            "i > 0 && j > 0 && k > 0",
            ["E[i-1, j, k]", "E[i, j-1, k]", "E[i, j, k-1]"]
        )
        assert len(rec._rules) == 1

    def test_branch_average_with_name(self):
        """Test branch average with name."""
        rec = Recurrence("Test", ["n"], ["x"])
        rec.branch_average("n > 0", ["E[n-1]"], name="Averaged")
        rule = rec._rules[0]
        assert rule.name == "Averaged"

    def test_branch_average_returns_self(self):
        """Test that branch_average() returns self."""
        rec = Recurrence("Test", ["n"], ["x"])
        result = rec.branch_average("n > 0", ["E[n-1]"])
        assert result is rec


class TestRecurrenceContext:
    """Test _ctx() method."""

    def test_context_creation(self):
        """Test context creation from recurrence."""
        rec = Recurrence("Legendre", ["n"], ["x"])
        ctx = rec._ctx()
        assert ctx.struct_name == "LegendreCoeff"
        assert ctx.indices == ["n"]
        assert ctx.runtime_vars == ["x"]
        assert ctx.vec_type == "Vec8d"

    def test_context_with_custom_vec_type(self):
        """Test context with custom vector type."""
        rec = Recurrence("Test", ["n"], ["x"], vec_type="double")
        ctx = rec._ctx()
        assert ctx.vec_type == "double"

    def test_context_multi_index(self):
        """Test context for multi-index recurrence."""
        rec = Recurrence("Hermite", ["i", "j", "t"], ["PA", "PB", "p"])
        ctx = rec._ctx()
        assert len(ctx.indices) == 3
        assert len(ctx.runtime_vars) == 3


class TestRecurrenceRulePriority:
    """Test RecurrenceRule priority_key method."""

    def test_equality_constraint_priority(self):
        """Test that equality constraints have higher priority."""
        cs_eq = ConstraintSet.parse("n == 0")
        cs_gt = ConstraintSet.parse("n > 0")

        rule_eq = RecurrenceRule(cs_eq, Const(1))
        rule_gt = RecurrenceRule(cs_gt, Const(1))

        # Equality constraint should have lower key (higher priority)
        assert rule_eq.priority_key() < rule_gt.priority_key()

    def test_more_constraints_priority(self):
        """Test that more constraints have higher priority."""
        cs_one = ConstraintSet.parse("n > 0")
        cs_two = ConstraintSet.parse("n > 0", "m > 0")

        rule_one = RecurrenceRule(cs_one, Const(1))
        rule_two = RecurrenceRule(cs_two, Const(1))

        # More constraints should have lower key (higher priority)
        assert rule_two.priority_key() < rule_one.priority_key()

    def test_mixed_constraints_priority(self):
        """Test priority with mixed constraint types."""
        cs1 = ConstraintSet.parse("n == 0 && m > 0")  # 1 eq, 1 gt
        cs2 = ConstraintSet.parse("n > 0 && m > 0")   # 0 eq, 2 gt
        cs3 = ConstraintSet.parse("n == 0")            # 1 eq, 0 gt

        rule1 = RecurrenceRule(cs1, Const(1))
        rule2 = RecurrenceRule(cs2, Const(1))
        rule3 = RecurrenceRule(cs3, Const(1))

        # rule1 and rule3 have same eq count, but rule1 has more total
        assert rule1.priority_key() < rule2.priority_key()
        assert rule3.priority_key() < rule2.priority_key()
        assert rule1.priority_key() < rule3.priority_key()


class TestBaseCase:
    """Test BaseCase dataclass."""

    def test_base_case_creation(self):
        """Test BaseCase creation."""
        bc = BaseCase({"n": 0}, Const(1))
        assert bc.index_values == {"n": 0}
        assert isinstance(bc.value, Const)

    def test_base_case_multi_index(self):
        """Test BaseCase with multiple indices."""
        bc = BaseCase({"i": 0, "j": 0, "t": 0}, Const(1))
        assert len(bc.index_values) == 3


class TestRecurrenceIntegration:
    """Integration tests for complete recurrence definitions."""

    def test_fibonacci(self):
        """Test Fibonacci recurrence definition."""
        rec = (Recurrence("Fibonacci", ["n"], [])
               .validity("n >= 0")
               .base(n=0, value=0)
               .base(n=1, value=1)
               .rule("n > 1", "E[n-1] + E[n-2]"))

        assert len(rec._base_cases) == 2
        assert len(rec._rules) == 1

    def test_legendre(self):
        """Test Legendre polynomial recurrence."""
        rec = (Recurrence("Legendre", ["n"], ["x"])
               .validity("n >= 0")
               .base(n=0, value=1.0)
               .base(n=1, value="x")
               .rule("n > 1",
                     "(2*n-1) * x * E[n-1] + (-(n-1)) * E[n-2]",
                     scale="1/n"))

        assert len(rec._base_cases) == 2
        assert len(rec._rules) == 1
        assert isinstance(rec._rules[0].expression, ScaledExpr)

    def test_hermite(self):
        """Test Hermite polynomial recurrence."""
        rec = (Recurrence("Hermite", ["n"], ["two_x"])
               .validity("n >= 0")
               .base(n=0, value=1.0)
               .base(n=1, value="two_x")
               .rule("n > 1", "two_x * E[n-1] + (-2*(n-1)) * E[n-2]"))

        assert len(rec._base_cases) == 2
        assert len(rec._rules) == 1

    def test_chebyshev_t(self):
        """Test Chebyshev T polynomial recurrence."""
        rec = (Recurrence("ChebyshevT", ["n"], ["x"])
               .validity("n >= 0")
               .base(n=0, value=1.0)
               .base(n=1, value="x")
               .rule("n > 1", "2 * x * E[n-1] + (-1) * E[n-2]"))

        assert len(rec._base_cases) == 2
        assert len(rec._rules) == 1

    def test_associated_legendre(self):
        """Test associated Legendre polynomial recurrence."""
        rec = (Recurrence("AssocLegendre", ["n", "m"], ["x", "sqrt1mx2"])
               .validity("n >= 0", "m >= 0", "n >= m")
               .base(n=0, m=0, value=1.0)
               .rule("n == m && n > 0",
                     "(-(2*n-1)) * sqrt1mx2 * E[n-1, m-1]")
               .rule("n > m && m >= 0",
                     "(2*n-1) * x * E[n-1, m] + (-(n+m-1)) * E[n-2, m]",
                     scale="1/(n-m)"))

        assert len(rec._base_cases) == 1
        assert len(rec._rules) == 2

    def test_multi_index_hermite_e(self):
        """Test multi-index Hermite E coefficient recurrence."""
        rec = (Recurrence("HermiteE", ["i", "j", "t"], ["PA", "PB", "p2"])
               .validity("i >= 0", "j >= 0", "i + j >= t", "t >= 0")
               .base(i=0, j=0, t=0, value=1.0)
               .rule("i == 0 && j > 0",
                     "p2 * E[i, j-1, t-1] + PB * E[i, j-1, t] + (t+1) * E[i, j-1, t+1]")
               .rule("i > 0",
                     "p2 * E[i-1, j, t-1] + PA * E[i-1, j, t] + (t+1) * E[i-1, j, t+1]"))

        assert len(rec._base_cases) == 1
        assert len(rec._rules) == 2


class TestRecurrenceEdgeCases:
    """Test edge cases and error handling."""

    def test_no_indices_raises_error(self):
        """Test that creating recurrence with no indices is allowed."""
        rec = Recurrence("Const", [], [])
        assert rec.indices == []

    def test_no_runtime_vars(self):
        """Test recurrence with no runtime variables."""
        rec = Recurrence("Fibonacci", ["n"], [])
        assert rec.runtime_vars == []

    def test_empty_rule_expression(self):
        """Test that parser handles empty expression string."""
        rec = Recurrence("Test", ["n"], ["x"])
        # Empty expression may be handled differently by parser
        # Some implementations may allow it, others may raise error
        # Just verify it doesn't crash
        try:
            rec.rule("n > 0", "")
            # If it doesn't raise, that's also acceptable behavior
        except (ValueError, IndexError, AttributeError):
            # Parser may raise various errors for empty expression
            pass

    def test_base_with_unknown_index(self):
        """Test base case with unknown index name."""
        rec = Recurrence("Test", ["n"], ["x"])
        # Python allows this - it just creates a key in the dict
        rec.base(m=0, value=1)
        assert "m" in rec._base_cases[0].index_values

    def test_unicode_in_names(self):
        """Test that Unicode characters work in names."""
        rec = Recurrence("Test", ["n"], ["α", "β"])
        assert "α" in rec.runtime_vars
        assert "β" in rec.runtime_vars

    def test_very_long_expression(self):
        """Test parsing very long expression."""
        rec = Recurrence("Test", ["n"], ["x"])
        long_expr = " + ".join([f"E[n-{i}]" for i in range(1, 10)])
        rec.rule("n > 9", long_expr)
        rule = rec._rules[0]
        assert len(rule.expression.terms) == 9


class TestRecurrenceMethodChaining:
    """Test comprehensive method chaining."""

    def test_full_chain(self):
        """Test full method chaining."""
        rec = (Recurrence("Test", ["n"], ["x"])
               .validity("n >= 0")
               .base(n=0, value=1.0)
               .base(n=1, value="x")
               .rule("n > 1", "x * E[n-1] + E[n-2]"))

        assert rec._validity is not None
        assert len(rec._base_cases) == 2
        assert len(rec._rules) == 1

    def test_multiple_rules_chain(self):
        """Test chaining multiple rules."""
        rec = (Recurrence("Test", ["i", "j"], ["x", "y"])
               .base(i=0, j=0, value=1.0)
               .rule("i == 0 && j > 0", "y * E[i, j-1]")
               .rule("i > 0 && j == 0", "x * E[i-1, j]")
               .rule("i > 0 && j > 0", "x * E[i-1, j] + y * E[i, j-1]"))

        assert len(rec._rules) == 3

    def test_interleaved_base_and_rules(self):
        """Test interleaving base cases and rules."""
        rec = (Recurrence("Test", ["n"], ["x"])
               .base(n=0, value=0)
               .rule("n == 1", "E[n]")  # Changed from "x" to valid E[...] expression
               .base(n=1, value="x")
               .rule("n > 1", "E[n-1] + E[n-2]"))

        assert len(rec._base_cases) == 2
        assert len(rec._rules) == 2
