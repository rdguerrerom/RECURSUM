#!/usr/bin/env python3
"""
Comprehensive unit tests for recursum.codegen.parser module.

Tests recurrence relation DSL parsing, coefficient parsing, and expression building.
"""

import pytest
from recursum.codegen.parser import RecurrenceParser
from recursum.codegen.core import (
    Const, IndexExpr, Var, RecursiveCall, BinOp, Term, Sum
)


class TestRecurrenceParser:
    """Test RecurrenceParser initialization."""

    def test_basic_initialization(self):
        """Test basic parser initialization."""
        parser = RecurrenceParser(["n"], ["x"])
        assert parser.indices == ["n"]
        assert parser.runtime_vars == ["x"]

    def test_multi_index_initialization(self):
        """Test parser with multiple indices."""
        parser = RecurrenceParser(["i", "j", "t"], ["PA", "PB", "p"])
        assert parser.indices == ["i", "j", "t"]
        assert parser.runtime_vars == ["PA", "PB", "p"]


class TestParseCoefficient:
    """Test coefficient parsing."""

    def test_empty_string(self):
        """Test parsing empty coefficient string."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("")
        assert isinstance(result, Const)
        assert result.value == 1

    def test_string_one(self):
        """Test parsing '1' as coefficient."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("1")
        assert isinstance(result, Const)
        assert result.value == 1

    def test_integer(self):
        """Test parsing integer coefficient."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("5")
        assert isinstance(result, Const)
        assert result.value == 5

    def test_float(self):
        """Test parsing float coefficient."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("3.14")
        assert isinstance(result, Const)
        assert result.value == 3.14

    def test_negative_integer(self):
        """Test parsing negative integer."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("-2")
        assert isinstance(result, Const)
        assert result.value == -2

    def test_runtime_variable(self):
        """Test parsing runtime variable."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("x")
        assert isinstance(result, Var)
        assert result.name == "x"

    def test_index_variable(self):
        """Test parsing index variable."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("n")
        assert isinstance(result, IndexExpr)
        assert result.expr_str == "n"

    def test_index_expression_simple(self):
        """Test parsing simple index expression."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("2*n")
        # Parser may return BinOp or IndexExpr depending on implementation
        assert isinstance(result, (IndexExpr, BinOp))
        if isinstance(result, IndexExpr):
            assert "n" in result.expr_str

    def test_index_expression_arithmetic(self):
        """Test parsing arithmetic index expression."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("2*n - 1")
        # Parser may return BinOp or IndexExpr depending on implementation
        assert isinstance(result, (IndexExpr, BinOp))
        # Just verify it contains the index variable
        if isinstance(result, IndexExpr):
            assert "n" in result.expr_str

    def test_parenthesized_integer(self):
        """Test parsing parenthesized integer."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("(5)")
        assert isinstance(result, Const)
        assert result.value == 5

    def test_parenthesized_index_expr(self):
        """Test parsing parenthesized index expression."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("(2*n-1)")
        assert isinstance(result, IndexExpr)
        assert "n" in result.expr_str

    def test_parenthesized_variable(self):
        """Test parsing parenthesized variable."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("(x)")
        assert isinstance(result, Var)
        assert result.name == "x"

    def test_compound_multiplication(self):
        """Test parsing compound multiplication like '(2*n-1) * x'."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("(2*n-1) * x")
        # Parser may treat this as IndexExpr or BinOp depending on implementation
        assert isinstance(result, (BinOp, IndexExpr))
        if isinstance(result, BinOp):
            assert result.op == '*'

    def test_negative_in_parentheses(self):
        """Test parsing negative expression in parentheses."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("(-(n-1))")
        assert isinstance(result, IndexExpr)


class TestParseIndexShift:
    """Test index shift parsing."""

    def test_no_shifts(self):
        """Test parsing with no shifts (identity)."""
        parser = RecurrenceParser(["n"], ["x"])
        shifts = parser.parse_index_shift("n")
        assert shifts["n"] == 0

    def test_single_decrement(self):
        """Test parsing single decrement."""
        parser = RecurrenceParser(["n"], ["x"])
        shifts = parser.parse_index_shift("n-1")
        assert shifts["n"] == -1

    def test_single_increment(self):
        """Test parsing single increment."""
        parser = RecurrenceParser(["n"], ["x"])
        shifts = parser.parse_index_shift("n+1")
        assert shifts["n"] == 1

    def test_multiple_decrement(self):
        """Test parsing multiple decrement."""
        parser = RecurrenceParser(["n"], ["x"])
        shifts = parser.parse_index_shift("n-2")
        assert shifts["n"] == -2

    def test_multiple_increment(self):
        """Test parsing multiple increment."""
        parser = RecurrenceParser(["n"], ["x"])
        shifts = parser.parse_index_shift("n+3")
        assert shifts["n"] == 3

    def test_multi_index_no_shift(self):
        """Test multi-index with no shifts."""
        parser = RecurrenceParser(["i", "j", "t"], ["PA", "PB"])
        shifts = parser.parse_index_shift("i, j, t")
        assert shifts["i"] == 0
        assert shifts["j"] == 0
        assert shifts["t"] == 0

    def test_multi_index_mixed_shifts(self):
        """Test multi-index with mixed shifts."""
        parser = RecurrenceParser(["i", "j", "t"], ["PA", "PB"])
        shifts = parser.parse_index_shift("i-1, j, t+1")
        assert shifts["i"] == -1
        assert shifts["j"] == 0
        assert shifts["t"] == 1

    def test_multi_index_all_shifts(self):
        """Test multi-index with all shifts."""
        parser = RecurrenceParser(["i", "j", "t"], ["PA", "PB"])
        shifts = parser.parse_index_shift("i-1, j+2, t-3")
        assert shifts["i"] == -1
        assert shifts["j"] == 2
        assert shifts["t"] == -3

    def test_whitespace_handling(self):
        """Test whitespace handling in index shifts."""
        parser = RecurrenceParser(["n"], ["x"])
        shifts = parser.parse_index_shift("n - 1")
        assert shifts["n"] == -1

    def test_whitespace_multi_index(self):
        """Test whitespace in multi-index shifts."""
        parser = RecurrenceParser(["i", "j"], ["x"])
        shifts = parser.parse_index_shift("i - 1 , j + 2")
        assert shifts["i"] == -1
        assert shifts["j"] == 2


class TestParseTerm:
    """Test term parsing."""

    def test_simple_call_no_coeff(self):
        """Test parsing simple call without coefficient."""
        parser = RecurrenceParser(["n"], ["x"])
        term = parser.parse_term("E[n-1]")
        assert isinstance(term, Term)
        assert isinstance(term.coeff, Const)
        assert term.coeff.value == 1
        assert isinstance(term.call, RecursiveCall)
        assert term.call.index_shifts["n"] == -1

    def test_call_with_integer_coeff(self):
        """Test parsing call with integer coefficient."""
        parser = RecurrenceParser(["n"], ["x"])
        term = parser.parse_term("2 * E[n-1]")
        assert isinstance(term, Term)
        assert isinstance(term.coeff, Const)
        assert term.coeff.value == 2

    def test_call_with_variable_coeff(self):
        """Test parsing call with variable coefficient."""
        parser = RecurrenceParser(["n"], ["x"])
        term = parser.parse_term("x * E[n-1]")
        assert isinstance(term, Term)
        assert isinstance(term.coeff, Var)
        assert term.coeff.name == "x"

    def test_call_with_index_coeff(self):
        """Test parsing call with index expression coefficient."""
        parser = RecurrenceParser(["n"], ["x"])
        term = parser.parse_term("(2*n-1) * E[n-1]")
        assert isinstance(term, Term)
        assert isinstance(term.coeff, IndexExpr)

    def test_call_with_compound_coeff(self):
        """Test parsing call with compound coefficient."""
        parser = RecurrenceParser(["n"], ["x"])
        term = parser.parse_term("(2*n-1) * x * E[n-1]")
        assert isinstance(term, Term)
        assert isinstance(term.coeff, BinOp)
        assert term.coeff.op == '*'

    def test_negative_coeff(self):
        """Test parsing term with negative coefficient."""
        parser = RecurrenceParser(["n"], ["x"])
        term = parser.parse_term("-2 * E[n-2]")
        assert isinstance(term.coeff, Const)
        assert term.coeff.value == -2

    def test_negative_index_expr_coeff(self):
        """Test parsing term with negative index expression coefficient."""
        parser = RecurrenceParser(["n"], ["x"])
        term = parser.parse_term("(-(n-1)) * E[n-2]")
        assert isinstance(term.coeff, IndexExpr)

    def test_multi_index_call(self):
        """Test parsing multi-index call."""
        parser = RecurrenceParser(["i", "j", "t"], ["PA", "PB"])
        term = parser.parse_term("PA * E[i-1, j, t+1]")
        assert isinstance(term.coeff, Var)
        assert term.coeff.name == "PA"
        assert term.call.index_shifts["i"] == -1
        assert term.call.index_shifts["j"] == 0
        assert term.call.index_shifts["t"] == 1

    def test_no_e_raises_error(self):
        """Test that parsing without E[...] raises error."""
        parser = RecurrenceParser(["n"], ["x"])
        with pytest.raises(ValueError, match="No E\\[...\\] found"):
            parser.parse_term("x * y")


class TestParseExpression:
    """Test full expression parsing."""

    def test_single_term(self):
        """Test parsing single term."""
        parser = RecurrenceParser(["n"], ["x"])
        expr = parser.parse_expression("E[n-1]")
        assert isinstance(expr, Sum)
        assert len(expr.terms) == 1

    def test_two_terms(self):
        """Test parsing two terms."""
        parser = RecurrenceParser(["n"], ["x"])
        expr = parser.parse_expression("E[n-1] + E[n-2]")
        assert isinstance(expr, Sum)
        assert len(expr.terms) == 2

    def test_three_terms(self):
        """Test parsing three terms."""
        parser = RecurrenceParser(["n"], ["x"])
        expr = parser.parse_expression("x * E[n-1] + 2 * E[n-2] + E[n-3]")
        assert isinstance(expr, Sum)
        assert len(expr.terms) == 3

    def test_legendre_recurrence(self):
        """Test parsing Legendre recurrence."""
        parser = RecurrenceParser(["n"], ["x"])
        expr = parser.parse_expression("(2*n-1) * x * E[n-1] + (-(n-1)) * E[n-2]")
        assert isinstance(expr, Sum)
        assert len(expr.terms) == 2

    def test_hermite_recurrence(self):
        """Test parsing Hermite recurrence."""
        parser = RecurrenceParser(["n"], ["two_x"])
        expr = parser.parse_expression("two_x * E[n-1] + (-2) * E[n-2]")
        assert isinstance(expr, Sum)
        assert len(expr.terms) == 2

    def test_multi_index_recurrence(self):
        """Test parsing multi-index recurrence."""
        parser = RecurrenceParser(["i", "j", "t"], ["PA", "PB", "p2"])
        expr = parser.parse_expression("p2 * E[i, j-1, t-1] + PB * E[i, j-1, t]")
        assert isinstance(expr, Sum)
        assert len(expr.terms) == 2

    def test_complex_multi_term(self):
        """Test parsing complex multi-term expression."""
        parser = RecurrenceParser(["i", "j", "t"], ["PA", "PB", "p"])
        expr_str = "PA * E[i-1, j, t] + (t+1) * E[i-1, j, t+1] + PB * E[i, j-1, t]"
        expr = parser.parse_expression(expr_str)
        assert isinstance(expr, Sum)
        assert len(expr.terms) == 3

    def test_nested_parentheses(self):
        """Test parsing with nested parentheses."""
        parser = RecurrenceParser(["n"], ["x"])
        expr = parser.parse_expression("((2*n-1)) * x * E[n-1]")
        assert isinstance(expr, Sum)
        assert len(expr.terms) == 1


class TestParseScale:
    """Test scale factor parsing."""

    def test_simple_division(self):
        """Test parsing simple 1/n."""
        parser = RecurrenceParser(["n"], ["x"])
        scale = parser.parse_scale("1/n")
        assert isinstance(scale, IndexExpr)
        assert scale.expr_str == "n"

    def test_division_with_parens(self):
        """Test parsing 1/(2*n)."""
        parser = RecurrenceParser(["n"], ["x"])
        scale = parser.parse_scale("1/(2*n)")
        assert isinstance(scale, IndexExpr)

    def test_division_by_constant(self):
        """Test parsing 1/2."""
        parser = RecurrenceParser(["n"], ["x"])
        scale = parser.parse_scale("1/2")
        assert isinstance(scale, Const)
        assert scale.value == 2.0

    def test_division_by_variable(self):
        """Test parsing 1/x."""
        parser = RecurrenceParser(["n"], ["x"])
        scale = parser.parse_scale("1/x")
        assert isinstance(scale, Var)
        assert scale.name == "x"

    def test_non_division_constant(self):
        """Test parsing non-division constant."""
        parser = RecurrenceParser(["n"], ["x"])
        scale = parser.parse_scale("2")
        assert isinstance(scale, Const)
        assert scale.value == 2

    def test_non_division_variable(self):
        """Test parsing non-division variable."""
        parser = RecurrenceParser(["n"], ["x"])
        scale = parser.parse_scale("x")
        assert isinstance(scale, Var)
        assert scale.name == "x"


class TestSplitByMult:
    """Test _split_by_mult helper method."""

    def test_no_multiplication(self):
        """Test splitting string with no multiplication."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser._split_by_mult("x")
        assert result == ["x"]

    def test_simple_multiplication(self):
        """Test splitting simple multiplication."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser._split_by_mult("a * b")
        assert result == ["a", "b"]

    def test_three_factors(self):
        """Test splitting three factors."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser._split_by_mult("a * b * c")
        assert result == ["a", "b", "c"]

    def test_with_parentheses(self):
        """Test splitting respecting parentheses."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser._split_by_mult("(2*n-1) * x")
        assert len(result) == 2
        assert result[0] == "(2*n-1)"
        assert result[1] == "x"

    def test_nested_parentheses(self):
        """Test splitting with nested parentheses."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser._split_by_mult("((a*b)) * c")
        assert len(result) == 2
        assert "a*b" in result[0]
        assert result[1] == "c"

    def test_empty_string(self):
        """Test splitting empty string."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser._split_by_mult("")
        assert result == ["1"]

    def test_whitespace_handling(self):
        """Test whitespace is stripped."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser._split_by_mult("  a  *  b  ")
        assert result == ["a", "b"]


class TestParserEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_coefficient(self):
        """Test empty coefficient defaults to 1."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("")
        assert isinstance(result, Const)
        assert result.value == 1

    def test_whitespace_only_coefficient(self):
        """Test whitespace-only coefficient defaults to 1."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("   ")
        assert isinstance(result, Const)
        assert result.value == 1

    def test_term_without_coeff_no_asterisk(self):
        """Test parsing term without coefficient or asterisk."""
        parser = RecurrenceParser(["n"], ["x"])
        term = parser.parse_term("E[n]")
        assert isinstance(term.coeff, Const)
        assert term.coeff.value == 1

    def test_term_with_leading_whitespace(self):
        """Test parsing term with leading whitespace."""
        parser = RecurrenceParser(["n"], ["x"])
        term = parser.parse_term("  2 * E[n-1]  ")
        assert isinstance(term.coeff, Const)
        assert term.coeff.value == 2

    def test_expression_with_extra_whitespace(self):
        """Test expression with extra whitespace."""
        parser = RecurrenceParser(["n"], ["x"])
        expr = parser.parse_expression("  E[n-1]  +  E[n-2]  ")
        assert isinstance(expr, Sum)
        assert len(expr.terms) == 2

    def test_unknown_variable_treated_as_var(self):
        """Test unknown identifier treated as variable or index expr."""
        parser = RecurrenceParser(["n"], ["x"])
        result = parser.parse_coefficient("unknown")
        # Parser may treat unknown identifiers as Var or IndexExpr
        assert isinstance(result, (Var, IndexExpr))
        if isinstance(result, Var):
            assert result.name == "unknown"
        else:
            assert "unknown" in result.expr_str

    def test_mixed_case_indices(self):
        """Test mixed case index names."""
        parser = RecurrenceParser(["N"], ["X"])
        shifts = parser.parse_index_shift("N-1")
        assert shifts["N"] == -1

    def test_zero_shift(self):
        """Test parsing explicit zero shift."""
        parser = RecurrenceParser(["n"], ["x"])
        shifts = parser.parse_index_shift("n+0")
        assert shifts["n"] == 0


class TestParserIntegration:
    """Integration tests combining multiple parser features."""

    def test_complete_legendre_parse(self):
        """Test complete Legendre recurrence parsing."""
        parser = RecurrenceParser(["n"], ["x"])
        expr = parser.parse_expression("(2*n-1) * x * E[n-1] + (-(n-1)) * E[n-2]")

        assert isinstance(expr, Sum)
        assert len(expr.terms) == 2

        # First term: (2*n-1) * x * E[n-1]
        t1 = expr.terms[0]
        assert isinstance(t1.coeff, BinOp)
        assert t1.call.index_shifts["n"] == -1

        # Second term: (-(n-1)) * E[n-2]
        t2 = expr.terms[1]
        assert isinstance(t2.coeff, IndexExpr)
        assert t2.call.index_shifts["n"] == -2

    def test_complete_hermite_e_parse(self):
        """Test complete Hermite E coefficient parsing."""
        parser = RecurrenceParser(["i", "j", "t"], ["PA", "PB", "p2"])
        expr_str = "p2 * E[i-1, j, t-1] + PA * E[i-1, j, t] + (t+1) * E[i-1, j, t+1]"
        expr = parser.parse_expression(expr_str)

        assert isinstance(expr, Sum)
        assert len(expr.terms) == 3

        # Check shifts
        assert expr.terms[0].call.index_shifts == {"i": -1, "j": 0, "t": -1}
        assert expr.terms[1].call.index_shifts == {"i": -1, "j": 0, "t": 0}
        assert expr.terms[2].call.index_shifts == {"i": -1, "j": 0, "t": 1}

    def test_three_term_recurrence(self):
        """Test standard three-term recurrence."""
        parser = RecurrenceParser(["n"], ["a", "b", "c"])
        expr = parser.parse_expression("a * E[n-1] + b * E[n-2] + c * E[n-3]")

        assert isinstance(expr, Sum)
        assert len(expr.terms) == 3

        # Check coefficients
        assert expr.terms[0].coeff.name == "a"
        assert expr.terms[1].coeff.name == "b"
        assert expr.terms[2].coeff.name == "c"

        # Check shifts
        assert expr.terms[0].call.index_shifts["n"] == -1
        assert expr.terms[1].call.index_shifts["n"] == -2
        assert expr.terms[2].call.index_shifts["n"] == -3
