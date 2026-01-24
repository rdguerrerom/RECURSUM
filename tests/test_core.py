#!/usr/bin/env python3
"""
Comprehensive unit tests for recursum.codegen.core module.

Tests AST expression nodes, code generation context, and expression operations.
"""

import pytest
from recursum.codegen.core import (
    Expr, Const, IndexExpr, Var, RecursiveCall, BinOp, Term, Sum, ScaledExpr,
    CodegenContext
)


class TestCodegenContext:
    """Test CodegenContext class."""

    def test_basic_creation(self):
        """Test basic context creation."""
        ctx = CodegenContext("TestStruct", ["n"], ["x"])
        assert ctx.struct_name == "TestStruct"
        assert ctx.indices == ["n"]
        assert ctx.runtime_vars == ["x"]
        assert ctx.vec_type == "Vec8d"

    def test_custom_vec_type(self):
        """Test context with custom vector type."""
        ctx = CodegenContext("TestStruct", ["n"], ["x"], vec_type="double")
        assert ctx.vec_type == "double"

    def test_multiple_indices(self):
        """Test context with multiple indices."""
        ctx = CodegenContext("HermiteCoeff", ["i", "j", "t"], ["PA", "PB", "p"])
        assert len(ctx.indices) == 3
        assert len(ctx.runtime_vars) == 3


class TestConst:
    """Test Const expression node."""

    def test_integer_const(self):
        """Test integer constant."""
        c = Const(5)
        assert c.value == 5

    def test_float_const(self):
        """Test float constant."""
        c = Const(3.14)
        assert c.value == 3.14

    def test_string_const(self):
        """Test string constant."""
        c = Const("some_const")
        assert c.value == "some_const"

    def test_to_cpp_integer(self):
        """Test C++ code generation for integer."""
        c = Const(42)
        ctx = CodegenContext("Test", ["n"], ["x"])
        assert c.to_cpp(ctx) == "Vec8d(42)"

    def test_to_cpp_float(self):
        """Test C++ code generation for float."""
        c = Const(2.5)
        ctx = CodegenContext("Test", ["n"], ["x"])
        assert c.to_cpp(ctx) == "Vec8d(2.5)"

    def test_to_cpp_string(self):
        """Test C++ code generation for string literal."""
        c = Const("M_PI")
        ctx = CodegenContext("Test", ["n"], ["x"])
        assert c.to_cpp(ctx) == "M_PI"

    def test_collect_calls(self):
        """Test that Const has no recursive calls."""
        c = Const(10)
        assert c.collect_calls() == []

    def test_uses_var(self):
        """Test that Const doesn't use any variables."""
        c = Const(5)
        assert not c.uses_var("x")
        assert not c.uses_var("n")


class TestVar:
    """Test Var expression node."""

    def test_var_creation(self):
        """Test variable creation."""
        v = Var("x")
        assert v.name == "x"

    def test_to_cpp(self):
        """Test C++ code generation."""
        v = Var("alpha")
        ctx = CodegenContext("Test", ["n"], ["alpha"])
        assert v.to_cpp(ctx) == "alpha"

    def test_collect_calls(self):
        """Test that Var has no recursive calls."""
        v = Var("x")
        assert v.collect_calls() == []

    def test_uses_var_match(self):
        """Test uses_var with matching name."""
        v = Var("x")
        assert v.uses_var("x")

    def test_uses_var_no_match(self):
        """Test uses_var with non-matching name."""
        v = Var("x")
        assert not v.uses_var("y")


class TestIndexExpr:
    """Test IndexExpr expression node."""

    def test_simple_index(self):
        """Test simple index expression."""
        ie = IndexExpr("n")
        assert ie.expr_str == "n"

    def test_arithmetic_index(self):
        """Test arithmetic index expression."""
        ie = IndexExpr("2*n - 1")
        assert ie.expr_str == "2*n - 1"

    def test_to_cpp(self):
        """Test C++ code generation."""
        ie = IndexExpr("2*n + 1")
        ctx = CodegenContext("Test", ["n"], ["x"])
        assert ie.to_cpp(ctx) == "Vec8d(2*n + 1)"

    def test_collect_calls(self):
        """Test that IndexExpr has no recursive calls."""
        ie = IndexExpr("n + m")
        assert ie.collect_calls() == []

    def test_uses_var_match(self):
        """Test uses_var with variable in expression."""
        ie = IndexExpr("2*n - 1")
        assert ie.uses_var("n")

    def test_uses_var_no_match(self):
        """Test uses_var with variable not in expression."""
        ie = IndexExpr("2*n - 1")
        assert not ie.uses_var("m")


class TestRecursiveCall:
    """Test RecursiveCall expression node."""

    def test_no_shifts(self):
        """Test recursive call with no index shifts."""
        rc = RecursiveCall({"n": 0})
        assert rc.index_shifts == {"n": 0}

    def test_single_shift(self):
        """Test recursive call with single shift."""
        rc = RecursiveCall({"n": -1})
        assert rc.index_shifts["n"] == -1

    def test_multiple_shifts(self):
        """Test recursive call with multiple shifts."""
        rc = RecursiveCall({"i": -1, "j": 0, "t": 1})
        assert rc.index_shifts["i"] == -1
        assert rc.index_shifts["j"] == 0
        assert rc.index_shifts["t"] == 1

    def test_to_cpp_no_shift(self):
        """Test C++ code generation with no shifts."""
        rc = RecursiveCall({"n": 0})
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        cpp = rc.to_cpp(ctx)
        assert "TestCoeff<n>::compute(x)" == cpp

    def test_to_cpp_negative_shift(self):
        """Test C++ code generation with negative shift."""
        rc = RecursiveCall({"n": -1})
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        cpp = rc.to_cpp(ctx)
        assert "TestCoeff<n - 1>::compute(x)" == cpp

    def test_to_cpp_positive_shift(self):
        """Test C++ code generation with positive shift."""
        rc = RecursiveCall({"n": 2})
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        cpp = rc.to_cpp(ctx)
        assert "TestCoeff<n + 2>::compute(x)" == cpp

    def test_to_cpp_multiple_indices(self):
        """Test C++ code generation with multiple indices."""
        rc = RecursiveCall({"i": -1, "j": 0, "t": 1})
        ctx = CodegenContext("HermiteCoeff", ["i", "j", "t"], ["PA", "PB", "p"])
        cpp = rc.to_cpp(ctx)
        assert "HermiteCoeff<i - 1, j, t + 1>::compute(PA, PB, p)" == cpp

    def test_collect_calls(self):
        """Test that RecursiveCall returns itself."""
        rc = RecursiveCall({"n": -1})
        calls = rc.collect_calls()
        assert len(calls) == 1
        assert calls[0] is rc

    def test_uses_var(self):
        """Test that RecursiveCall doesn't use runtime variables."""
        rc = RecursiveCall({"n": -1})
        assert not rc.uses_var("x")


class TestBinOp:
    """Test BinOp expression node."""

    def test_addition(self):
        """Test binary addition."""
        left = Const(2)
        right = Const(3)
        binop = BinOp('+', left, right)
        assert binop.op == '+'
        assert binop.left is left
        assert binop.right is right

    def test_multiplication(self):
        """Test binary multiplication."""
        left = Var("x")
        right = Const(2)
        binop = BinOp('*', left, right)
        assert binop.op == '*'

    def test_to_cpp_simple(self):
        """Test C++ generation for simple binary op."""
        left = Const(2)
        right = Const(3)
        binop = BinOp('+', left, right)
        ctx = CodegenContext("Test", ["n"], ["x"])
        cpp = binop.to_cpp(ctx)
        assert "Vec8d(2) + Vec8d(3)" == cpp

    def test_to_cpp_with_var(self):
        """Test C++ generation with variable."""
        left = Var("x")
        right = Const(2)
        binop = BinOp('*', left, right)
        ctx = CodegenContext("Test", ["n"], ["x"])
        cpp = binop.to_cpp(ctx)
        assert "x * Vec8d(2)" == cpp

    def test_to_cpp_nested(self):
        """Test C++ generation for nested binary ops."""
        inner = BinOp('+', Const(1), Const(2))
        outer = BinOp('*', inner, Const(3))
        ctx = CodegenContext("Test", ["n"], ["x"])
        cpp = outer.to_cpp(ctx)
        assert "(Vec8d(1) + Vec8d(2)) * Vec8d(3)" == cpp

    def test_collect_calls_no_calls(self):
        """Test collect_calls with no recursive calls."""
        binop = BinOp('+', Const(1), Var("x"))
        assert binop.collect_calls() == []

    def test_collect_calls_with_calls(self):
        """Test collect_calls with recursive calls."""
        rc1 = RecursiveCall({"n": -1})
        rc2 = RecursiveCall({"n": -2})
        binop = BinOp('+', rc1, rc2)
        calls = binop.collect_calls()
        assert len(calls) == 2
        assert rc1 in calls
        assert rc2 in calls

    def test_uses_var_left(self):
        """Test uses_var when variable is in left operand."""
        binop = BinOp('+', Var("x"), Const(1))
        assert binop.uses_var("x")

    def test_uses_var_right(self):
        """Test uses_var when variable is in right operand."""
        binop = BinOp('*', Const(2), Var("y"))
        assert binop.uses_var("y")

    def test_uses_var_neither(self):
        """Test uses_var when variable is in neither operand."""
        binop = BinOp('+', Const(1), Const(2))
        assert not binop.uses_var("x")


class TestTerm:
    """Test Term expression node."""

    def test_unit_coefficient(self):
        """Test term with coefficient of 1."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(1), call)
        assert isinstance(term.coeff, Const)
        assert term.coeff.value == 1

    def test_constant_coefficient(self):
        """Test term with constant coefficient."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(2), call)
        assert term.coeff.value == 2

    def test_variable_coefficient(self):
        """Test term with variable coefficient."""
        call = RecursiveCall({"n": -1})
        term = Term(Var("x"), call)
        assert isinstance(term.coeff, Var)
        assert term.coeff.name == "x"

    def test_to_cpp_unit_coeff(self):
        """Test C++ generation with unit coefficient."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(1), call)
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        cpp = term.to_cpp(ctx)
        assert "TestCoeff<n - 1>::compute(x)" == cpp
        # Should not include multiplication by 1

    def test_to_cpp_const_coeff(self):
        """Test C++ generation with constant coefficient."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(2), call)
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        cpp = term.to_cpp(ctx)
        assert "Vec8d(2) * TestCoeff<n - 1>::compute(x)" == cpp

    def test_to_cpp_var_coeff(self):
        """Test C++ generation with variable coefficient."""
        call = RecursiveCall({"n": -1})
        term = Term(Var("x"), call)
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        cpp = term.to_cpp(ctx)
        assert "x * TestCoeff<n - 1>::compute(x)" == cpp

    def test_collect_calls(self):
        """Test that Term returns its recursive call."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(2), call)
        calls = term.collect_calls()
        assert len(calls) == 1
        assert calls[0] is call

    def test_uses_var_in_coeff(self):
        """Test uses_var when variable is in coefficient."""
        call = RecursiveCall({"n": -1})
        term = Term(Var("x"), call)
        assert term.uses_var("x")

    def test_uses_var_not_present(self):
        """Test uses_var when variable is not present."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(2), call)
        assert not term.uses_var("x")


class TestSum:
    """Test Sum expression node."""

    def test_empty_sum(self):
        """Test empty sum."""
        s = Sum([])
        assert s.terms == []

    def test_single_term(self):
        """Test sum with single term."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(1), call)
        s = Sum([term])
        assert len(s.terms) == 1

    def test_multiple_terms(self):
        """Test sum with multiple terms."""
        t1 = Term(Const(1), RecursiveCall({"n": -1}))
        t2 = Term(Const(2), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])
        assert len(s.terms) == 2

    def test_to_cpp_empty(self):
        """Test C++ generation for empty sum."""
        s = Sum([])
        ctx = CodegenContext("Test", ["n"], ["x"])
        cpp = s.to_cpp(ctx)
        assert cpp == "Vec8d(0.0)"

    def test_to_cpp_single_term(self):
        """Test C++ generation for single term."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(1), call)
        s = Sum([term])
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        cpp = s.to_cpp(ctx)
        assert "TestCoeff<n - 1>::compute(x)" == cpp

    def test_to_cpp_multiple_terms(self):
        """Test C++ generation for multiple terms."""
        t1 = Term(Const(1), RecursiveCall({"n": -1}))
        t2 = Term(Const(2), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        cpp = s.to_cpp(ctx)
        assert " + " in cpp
        assert "TestCoeff<n - 1>" in cpp
        assert "TestCoeff<n - 2>" in cpp

    def test_collect_calls(self):
        """Test collect_calls with multiple terms."""
        call1 = RecursiveCall({"n": -1})
        call2 = RecursiveCall({"n": -2})
        t1 = Term(Const(1), call1)
        t2 = Term(Const(2), call2)
        s = Sum([t1, t2])
        calls = s.collect_calls()
        assert len(calls) == 2
        assert call1 in calls
        assert call2 in calls

    def test_uses_var(self):
        """Test uses_var across multiple terms."""
        t1 = Term(Var("x"), RecursiveCall({"n": -1}))
        t2 = Term(Const(2), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])
        assert s.uses_var("x")
        assert not s.uses_var("y")


class TestScaledExpr:
    """Test ScaledExpr expression node."""

    def test_division(self):
        """Test scaled expression with division."""
        expr = Const(10)
        scale = Const(2)
        scaled = ScaledExpr(expr, scale, is_division=True)
        assert scaled.is_division is True

    def test_multiplication(self):
        """Test scaled expression with multiplication."""
        expr = Const(10)
        scale = Const(2)
        scaled = ScaledExpr(expr, scale, is_division=False)
        assert scaled.is_division is False

    def test_to_cpp_division(self):
        """Test C++ generation for division."""
        expr = Const(10)
        scale = Const(2)
        scaled = ScaledExpr(expr, scale, is_division=True)
        ctx = CodegenContext("Test", ["n"], ["x"])
        cpp = scaled.to_cpp(ctx)
        assert cpp == "(Vec8d(10)) / (Vec8d(2))"

    def test_to_cpp_multiplication(self):
        """Test C++ generation for multiplication."""
        expr = Const(10)
        scale = Const(2)
        scaled = ScaledExpr(expr, scale, is_division=False)
        ctx = CodegenContext("Test", ["n"], ["x"])
        cpp = scaled.to_cpp(ctx)
        assert cpp == "(Vec8d(10)) * (Vec8d(2))"

    def test_to_cpp_with_index_scale(self):
        """Test C++ generation with index expression scale."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(1), call)
        s = Sum([term])
        scale = IndexExpr("n")
        scaled = ScaledExpr(s, scale, is_division=True)
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        cpp = scaled.to_cpp(ctx)
        assert "/" in cpp
        assert "Vec8d(n)" in cpp

    def test_collect_calls(self):
        """Test collect_calls from inner expression."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(1), call)
        s = Sum([term])
        scaled = ScaledExpr(s, Const(2), is_division=True)
        calls = scaled.collect_calls()
        assert len(calls) == 1
        assert calls[0] is call

    def test_uses_var_in_expr(self):
        """Test uses_var when variable is in expression."""
        t = Term(Var("x"), RecursiveCall({"n": -1}))
        s = Sum([t])
        scaled = ScaledExpr(s, Const(2), is_division=True)
        assert scaled.uses_var("x")

    def test_uses_var_in_scale(self):
        """Test uses_var when variable is in scale."""
        call = RecursiveCall({"n": -1})
        t = Term(Const(1), call)
        s = Sum([t])
        scaled = ScaledExpr(s, Var("y"), is_division=True)
        assert scaled.uses_var("y")


class TestExprOperatorOverloads:
    """Test operator overloads on Expr base class."""

    def test_add_const(self):
        """Test adding constant to expression."""
        v = Var("x")
        result = v + 5
        assert isinstance(result, BinOp)
        assert result.op == '+'
        assert result.left is v
        assert isinstance(result.right, Const)
        assert result.right.value == 5

    def test_add_expr(self):
        """Test adding two expressions."""
        v1 = Var("x")
        v2 = Var("y")
        result = v1 + v2
        assert isinstance(result, BinOp)
        assert result.op == '+'
        assert result.left is v1
        assert result.right is v2

    def test_mul_const(self):
        """Test multiplying expression by constant."""
        v = Var("x")
        result = v * 3
        assert isinstance(result, BinOp)
        assert result.op == '*'
        assert result.left is v
        assert isinstance(result.right, Const)
        assert result.right.value == 3

    def test_mul_expr(self):
        """Test multiplying two expressions."""
        v1 = Var("x")
        v2 = Var("y")
        result = v1 * v2
        assert isinstance(result, BinOp)
        assert result.op == '*'
        assert result.left is v1
        assert result.right is v2

    def test_chained_operations(self):
        """Test chained operator overloads."""
        v = Var("x")
        result = (v + 1) * 2
        assert isinstance(result, BinOp)
        assert result.op == '*'
        assert isinstance(result.left, BinOp)
        assert result.left.op == '+'


class TestComplexExpressions:
    """Test complex expression trees."""

    def test_legendre_recurrence(self):
        """Test expression for Legendre recurrence."""
        # (2n-1)*x*P[n-1] - (n-1)*P[n-2]
        t1 = Term(BinOp('*', IndexExpr("2*n-1"), Var("x")), RecursiveCall({"n": -1}))
        t2 = Term(IndexExpr("-(n-1)"), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])
        scale = IndexExpr("n")
        expr = ScaledExpr(s, scale, is_division=True)

        ctx = CodegenContext("LegendreCoeff", ["n"], ["x"])
        cpp = expr.to_cpp(ctx)

        assert "n - 1" in cpp
        assert "n - 2" in cpp
        assert "/" in cpp

    def test_hermite_recurrence(self):
        """Test expression for Hermite recurrence."""
        # 2*x*H[n-1] - 2*(n-1)*H[n-2]
        t1 = Term(BinOp('*', Const(2), Var("x")), RecursiveCall({"n": -1}))
        t2 = Term(BinOp('*', Const(-2), IndexExpr("n-1")), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])

        ctx = CodegenContext("HermiteCoeff", ["n"], ["x"])
        cpp = s.to_cpp(ctx)

        assert "HermiteCoeff<n - 1>" in cpp
        assert "HermiteCoeff<n - 2>" in cpp

    def test_multi_index_recurrence(self):
        """Test multi-index recurrence expression."""
        # E[i-1, j] + E[i, j-1]
        t1 = Term(Const(1), RecursiveCall({"i": -1, "j": 0, "t": 0}))
        t2 = Term(Const(1), RecursiveCall({"i": 0, "j": -1, "t": 0}))
        s = Sum([t1, t2])

        ctx = CodegenContext("MultiCoeff", ["i", "j", "t"], ["x", "y"])
        cpp = s.to_cpp(ctx)

        assert "i - 1, j, t" in cpp
        assert "i, j - 1, t" in cpp
