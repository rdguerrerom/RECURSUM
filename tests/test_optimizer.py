#!/usr/bin/env python3
"""
Comprehensive unit tests for recursum.codegen.optimizer module.

Tests CSE analysis, expression optimization, memoization generation, and utility functions.
"""

import pytest
from recursum.codegen.optimizer import (
    SubExprKey, CSEAnalyzer, OptimizedExpr, ExpressionOptimizer,
    CachedVar, OptimizedSum, OptimizedCodeGenerator,
    MemoizationInfo, MemoizationGenerator,
    count_operations, estimate_cost, should_apply_cse
)
from recursum.codegen.core import (
    Const, Var, IndexExpr, RecursiveCall, BinOp, Term, Sum, ScaledExpr,
    CodegenContext
)


class TestSubExprKey:
    """Test SubExprKey class."""

    def test_creation(self):
        """Test SubExprKey creation."""
        key = SubExprKey('call', 'n:-1')
        assert key.kind == 'call'
        assert key.signature == 'n:-1'

    def test_equality(self):
        """Test SubExprKey equality."""
        key1 = SubExprKey('call', 'n:-1')
        key2 = SubExprKey('call', 'n:-1')
        assert key1 == key2

    def test_inequality_different_kind(self):
        """Test inequality with different kind."""
        key1 = SubExprKey('call', 'n:-1')
        key2 = SubExprKey('coeff', 'n:-1')
        assert key1 != key2

    def test_inequality_different_signature(self):
        """Test inequality with different signature."""
        key1 = SubExprKey('call', 'n:-1')
        key2 = SubExprKey('call', 'n:-2')
        assert key1 != key2

    def test_hashable(self):
        """Test that SubExprKey is hashable."""
        key1 = SubExprKey('call', 'n:-1')
        key2 = SubExprKey('call', 'n:-1')
        d = {key1: "value"}
        assert d[key2] == "value"

    def test_hash_inequality(self):
        """Test that different keys have different hashes."""
        key1 = SubExprKey('call', 'n:-1')
        key2 = SubExprKey('call', 'n:-2')
        assert hash(key1) != hash(key2)


class TestCSEAnalyzer:
    """Test CSEAnalyzer class."""

    def test_initialization(self):
        """Test CSEAnalyzer initialization."""
        analyzer = CSEAnalyzer()
        assert analyzer.occurrences == {}
        assert analyzer.cse_counter == 0

    def test_analyze_const(self):
        """Test analyzing constant expression."""
        analyzer = CSEAnalyzer()
        analyzer.analyze_expr(Const(5))
        # Constants don't create occurrences
        assert len(analyzer.occurrences) == 0

    def test_analyze_var(self):
        """Test analyzing variable expression."""
        analyzer = CSEAnalyzer()
        analyzer.analyze_expr(Var("x"))
        # Variables don't create occurrences
        assert len(analyzer.occurrences) == 0

    def test_analyze_recursive_call(self):
        """Test analyzing recursive call."""
        analyzer = CSEAnalyzer()
        call = RecursiveCall({"n": -1})
        analyzer.analyze_expr(call)
        assert len(analyzer.occurrences) == 1

    def test_analyze_term(self):
        """Test analyzing term."""
        analyzer = CSEAnalyzer()
        call = RecursiveCall({"n": -1})
        term = Term(Const(2), call)
        analyzer.analyze_expr(term)
        # Should track the call
        assert len(analyzer.occurrences) >= 1

    def test_analyze_term_with_complex_coeff(self):
        """Test analyzing term with complex coefficient."""
        analyzer = CSEAnalyzer()
        call = RecursiveCall({"n": -1})
        coeff = IndexExpr("2*n-1")
        term = Term(coeff, call)
        analyzer.analyze_expr(term)
        # Should track both call and coefficient
        assert len(analyzer.occurrences) >= 2

    def test_analyze_sum(self):
        """Test analyzing sum expression."""
        analyzer = CSEAnalyzer()
        t1 = Term(Const(1), RecursiveCall({"n": -1}))
        t2 = Term(Const(2), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])
        analyzer.analyze_expr(s)
        assert len(analyzer.occurrences) >= 2

    def test_analyze_duplicate_calls(self):
        """Test analyzing duplicate recursive calls."""
        analyzer = CSEAnalyzer()
        call1 = RecursiveCall({"n": -1})
        call2 = RecursiveCall({"n": -1})
        t1 = Term(Const(1), call1)
        t2 = Term(Var("x"), call2)
        s = Sum([t1, t2])
        analyzer.analyze_expr(s)

        # Same call should be detected
        common = analyzer.get_common_subexpressions(min_occurrences=2)
        assert len(common) >= 1

    def test_analyze_scaled_expr(self):
        """Test analyzing scaled expression."""
        analyzer = CSEAnalyzer()
        call = RecursiveCall({"n": -1})
        term = Term(Const(1), call)
        s = Sum([term])
        scale = IndexExpr("n")
        scaled = ScaledExpr(s, scale, is_division=True)
        analyzer.analyze_expr(scaled)
        assert len(analyzer.occurrences) >= 2  # call + scale

    def test_get_common_subexpressions(self):
        """Test getting common subexpressions."""
        analyzer = CSEAnalyzer()
        call = RecursiveCall({"n": -1})
        analyzer.analyze_expr(call)
        analyzer.analyze_expr(call)

        common = analyzer.get_common_subexpressions(min_occurrences=2)
        assert len(common) >= 1

    def test_get_common_subexpressions_threshold(self):
        """Test common subexpressions with different thresholds."""
        analyzer = CSEAnalyzer()
        call = RecursiveCall({"n": -1})
        analyzer.analyze_expr(call)

        # With threshold=2, should be empty
        common = analyzer.get_common_subexpressions(min_occurrences=2)
        assert len(common) == 0

        # With threshold=1, should have the call
        common = analyzer.get_common_subexpressions(min_occurrences=1)
        assert len(common) >= 1

    def test_generate_cse_name_call(self):
        """Test CSE name generation for calls."""
        analyzer = CSEAnalyzer()
        key = SubExprKey('call', 'n:-1')
        name = analyzer.generate_cse_name(key)
        assert '_cse_call_' in name

    def test_generate_cse_name_coeff(self):
        """Test CSE name generation for coefficients."""
        analyzer = CSEAnalyzer()
        key = SubExprKey('coeff', '2*n-1')
        name = analyzer.generate_cse_name(key)
        assert '_cse_coeff_' in name

    def test_generate_cse_name_increments(self):
        """Test that CSE name generation increments counter."""
        analyzer = CSEAnalyzer()
        key1 = SubExprKey('call', 'n:-1')
        key2 = SubExprKey('call', 'n:-2')
        name1 = analyzer.generate_cse_name(key1)
        name2 = analyzer.generate_cse_name(key2)
        assert name1 != name2


class TestExpressionOptimizer:
    """Test ExpressionOptimizer class."""

    def test_initialization(self):
        """Test optimizer initialization."""
        ctx = CodegenContext("Test", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx)
        assert opt.ctx is ctx
        assert opt.enable_cse is True
        assert opt.enable_horner is True

    def test_initialization_with_options(self):
        """Test optimizer with custom options."""
        ctx = CodegenContext("Test", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx, enable_cse=False, enable_horner=False)
        assert opt.enable_cse is False
        assert opt.enable_horner is False

    def test_optimize_simple_expr(self):
        """Test optimizing simple expression."""
        ctx = CodegenContext("Test", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx)
        call = RecursiveCall({"n": -1})
        result = opt.optimize_expression(call)
        assert result.result_expr is not None

    def test_optimize_with_cse_disabled(self):
        """Test optimization with CSE disabled."""
        ctx = CodegenContext("Test", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx, enable_cse=False)
        t1 = Term(Const(1), RecursiveCall({"n": -1}))
        t2 = Term(Const(2), RecursiveCall({"n": -1}))
        s = Sum([t1, t2])
        result = opt.optimize_expression(s)
        # Should not create intermediates when CSE disabled
        assert result.result_expr is s

    def test_optimize_with_duplicate_calls(self):
        """Test optimization extracts duplicate calls."""
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx, enable_cse=True, cse_threshold=2)
        call1 = RecursiveCall({"n": -1})
        call2 = RecursiveCall({"n": -1})
        t1 = Term(Var("x"), call1)
        t2 = Term(Const(2), call2)
        s = Sum([t1, t2])

        result = opt.optimize_expression(s)
        # Should create intermediates for duplicate call
        assert len(result.intermediates) >= 1

    def test_call_signature(self):
        """Test _call_signature method."""
        ctx = CodegenContext("Test", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx)
        call = RecursiveCall({"n": -1})
        sig = opt._call_signature(call)
        assert "n:-1" in sig

    def test_call_signature_multi_index(self):
        """Test _call_signature with multiple indices."""
        ctx = CodegenContext("Test", ["i", "j"], ["x"])
        opt = ExpressionOptimizer(ctx)
        call = RecursiveCall({"i": -1, "j": 2})
        sig = opt._call_signature(call)
        # Should be sorted
        assert "i:-1" in sig
        assert "j:2" in sig


class TestCachedVar:
    """Test CachedVar class."""

    def test_creation(self):
        """Test CachedVar creation."""
        cv = CachedVar("e_0")
        assert cv.name == "e_0"

    def test_to_cpp(self):
        """Test C++ code generation."""
        cv = CachedVar("e_0")
        ctx = CodegenContext("Test", ["n"], ["x"])
        assert cv.to_cpp(ctx) == "e_0"

    def test_collect_calls(self):
        """Test that CachedVar has no recursive calls."""
        cv = CachedVar("e_0")
        assert cv.collect_calls() == []

    def test_uses_var_match(self):
        """Test uses_var with matching name."""
        cv = CachedVar("e_0")
        assert cv.uses_var("e_0")

    def test_uses_var_no_match(self):
        """Test uses_var with non-matching name."""
        cv = CachedVar("e_0")
        assert not cv.uses_var("e_1")


class TestOptimizedSum:
    """Test OptimizedSum class."""

    def test_empty_sum(self):
        """Test empty optimized sum."""
        os = OptimizedSum([])
        assert os.exprs == []

    def test_single_expr(self):
        """Test optimized sum with single expression."""
        cv = CachedVar("e_0")
        os = OptimizedSum([cv])
        assert len(os.exprs) == 1

    def test_to_cpp_empty(self):
        """Test C++ generation for empty sum."""
        os = OptimizedSum([])
        ctx = CodegenContext("Test", ["n"], ["x"])
        assert os.to_cpp(ctx) == "Vec8d(0.0)"

    def test_to_cpp_single(self):
        """Test C++ generation for single expression."""
        cv = CachedVar("e_0")
        os = OptimizedSum([cv])
        ctx = CodegenContext("Test", ["n"], ["x"])
        assert os.to_cpp(ctx) == "e_0"

    def test_to_cpp_multiple(self):
        """Test C++ generation for multiple expressions."""
        cv1 = CachedVar("e_0")
        cv2 = CachedVar("e_1")
        os = OptimizedSum([cv1, cv2])
        ctx = CodegenContext("Test", ["n"], ["x"])
        cpp = os.to_cpp(ctx)
        assert " + " in cpp
        assert "e_0" in cpp
        assert "e_1" in cpp

    def test_collect_calls(self):
        """Test collect_calls for OptimizedSum."""
        cv = CachedVar("e_0")
        os = OptimizedSum([cv])
        assert os.collect_calls() == []

    def test_uses_var(self):
        """Test uses_var for OptimizedSum."""
        v1 = Var("x")
        v2 = Var("y")
        os = OptimizedSum([v1, v2])
        assert os.uses_var("x")
        assert os.uses_var("y")
        assert not os.uses_var("z")


class TestOptimizedCodeGenerator:
    """Test OptimizedCodeGenerator class."""

    def test_initialization(self):
        """Test OptimizedCodeGenerator initialization."""
        ctx = CodegenContext("Test", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx)
        gen = OptimizedCodeGenerator(ctx, opt)
        assert gen.ctx is ctx
        assert gen.optimizer is opt

    def test_generate_body_simple(self):
        """Test generating body for simple expression."""
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx, enable_cse=False)
        gen = OptimizedCodeGenerator(ctx, opt)

        call = RecursiveCall({"n": -1})
        body = gen.generate_body(call)

        assert "return" in body
        assert "TestCoeff<n - 1>" in body

    def test_generate_body_with_intermediates(self):
        """Test generating body with intermediate variables."""
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx, enable_cse=True, cse_threshold=2)
        gen = OptimizedCodeGenerator(ctx, opt)

        # Create expression with duplicate calls
        call1 = RecursiveCall({"n": -1})
        call2 = RecursiveCall({"n": -1})
        t1 = Term(Var("x"), call1)
        t2 = Term(Const(2), call2)
        s = Sum([t1, t2])

        body = gen.generate_body(s)

        assert "Vec8d" in body
        assert "return" in body


class TestMemoizationInfo:
    """Test MemoizationInfo class."""

    def test_default_values(self):
        """Test default MemoizationInfo values."""
        info = MemoizationInfo()
        assert info.storage_type == "inline_static"
        assert info.use_optional is True

    def test_custom_values(self):
        """Test custom MemoizationInfo values."""
        info = MemoizationInfo(storage_type="static_member", use_optional=False)
        assert info.storage_type == "static_member"
        assert info.use_optional is False


class TestMemoizationGenerator:
    """Test MemoizationGenerator class."""

    def test_initialization(self):
        """Test MemoizationGenerator initialization."""
        ctx = CodegenContext("Test", ["n"], ["x"])
        gen = MemoizationGenerator(ctx)
        assert gen.ctx is ctx
        assert gen.info is not None

    def test_initialization_with_info(self):
        """Test initialization with custom MemoizationInfo."""
        ctx = CodegenContext("Test", ["n"], ["x"])
        info = MemoizationInfo(storage_type="static_member")
        gen = MemoizationGenerator(ctx, info)
        assert gen.info.storage_type == "static_member"

    def test_generate_inline_static(self):
        """Test generating inline static memoization."""
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        gen = MemoizationGenerator(ctx)
        code = gen.generate_memoized_struct("int n", "double x", "            result = 1.0;\n")
        assert "inline static" in code or "static" in code
        assert "Memoized" in code

    def test_generate_static_member(self):
        """Test generating static member memoization."""
        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        info = MemoizationInfo(storage_type="static_member")
        gen = MemoizationGenerator(ctx, info)
        code = gen.generate_memoized_struct("int n", "double x", "    return 1.0;")
        assert "static" in code


class TestCountOperations:
    """Test count_operations utility function."""

    def test_const_no_ops(self):
        """Test counting operations in constant."""
        c = Const(5)
        counts = count_operations(c)
        assert counts.get('add', 0) == 0
        assert counts.get('mul', 0) == 0

    def test_recursive_call(self):
        """Test counting recursive call."""
        call = RecursiveCall({"n": -1})
        counts = count_operations(call)
        assert counts['call'] == 1

    def test_binop_add(self):
        """Test counting addition."""
        binop = BinOp('+', Const(1), Const(2))
        counts = count_operations(binop)
        assert counts['add'] == 1

    def test_binop_mul(self):
        """Test counting multiplication."""
        binop = BinOp('*', Var("x"), Const(2))
        counts = count_operations(binop)
        assert counts['mul'] == 1

    def test_term_with_coeff(self):
        """Test counting term with non-unit coefficient."""
        term = Term(Const(2), RecursiveCall({"n": -1}))
        counts = count_operations(term)
        assert counts['mul'] == 1
        assert counts['call'] == 1

    def test_term_with_unit_coeff(self):
        """Test counting term with unit coefficient."""
        term = Term(Const(1), RecursiveCall({"n": -1}))
        counts = count_operations(term)
        assert counts.get('mul', 0) == 0
        assert counts['call'] == 1

    def test_sum_of_terms(self):
        """Test counting sum of terms."""
        t1 = Term(Const(1), RecursiveCall({"n": -1}))
        t2 = Term(Const(2), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])
        counts = count_operations(s)
        assert counts['add'] == 1  # One addition between terms
        assert counts['call'] == 2  # Two calls
        assert counts['mul'] == 1   # One multiplication (from second term)

    def test_scaled_division(self):
        """Test counting scaled expression with division."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(1), call)
        s = Sum([term])
        scaled = ScaledExpr(s, Const(2), is_division=True)
        counts = count_operations(scaled)
        assert counts['div'] == 1
        assert counts['call'] == 1

    def test_scaled_multiplication(self):
        """Test counting scaled expression with multiplication."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(1), call)
        s = Sum([term])
        scaled = ScaledExpr(s, Const(2), is_division=False)
        counts = count_operations(scaled)
        assert counts['mul'] == 1
        assert counts['call'] == 1


class TestEstimateCost:
    """Test estimate_cost utility function."""

    def test_const_zero_cost(self):
        """Test that constant has zero cost."""
        c = Const(5)
        cost = estimate_cost(c)
        assert cost == 0

    def test_recursive_call_cost(self):
        """Test recursive call cost."""
        call = RecursiveCall({"n": -1})
        cost = estimate_cost(call)
        assert cost == 50  # Default cost for recursive call

    def test_addition_cost(self):
        """Test addition cost."""
        binop = BinOp('+', Const(1), Const(2))
        cost = estimate_cost(binop)
        assert cost == 1

    def test_multiplication_cost(self):
        """Test multiplication cost."""
        binop = BinOp('*', Var("x"), Const(2))
        cost = estimate_cost(binop)
        assert cost == 2

    def test_division_cost(self):
        """Test division cost."""
        call = RecursiveCall({"n": -1})
        term = Term(Const(1), call)
        s = Sum([term])
        scaled = ScaledExpr(s, Const(2), is_division=True)
        cost = estimate_cost(scaled)
        assert cost >= 10  # Division is expensive

    def test_complex_expression_cost(self):
        """Test cost of complex expression."""
        # x * E[n-1] + 2 * E[n-2]
        t1 = Term(Var("x"), RecursiveCall({"n": -1}))
        t2 = Term(Const(2), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])
        cost = estimate_cost(s)
        # 2 calls (100) + 1 mul (2) + 1 add (1) = 103
        assert cost > 100


class TestShouldApplyCSE:
    """Test should_apply_cse utility function."""

    def test_single_call_no_cse(self):
        """Test that single call doesn't need CSE."""
        call = RecursiveCall({"n": -1})
        assert should_apply_cse(call) is False

    def test_two_different_calls_needs_cse(self):
        """Test that multiple different calls benefit from CSE."""
        t1 = Term(Const(1), RecursiveCall({"n": -1}))
        t2 = Term(Const(2), RecursiveCall({"n": -2}))
        t3 = Term(Const(3), RecursiveCall({"n": -3}))
        s = Sum([t1, t2, t3])
        assert should_apply_cse(s) is True

    def test_duplicate_calls_needs_cse(self):
        """Test that duplicate calls need CSE."""
        call1 = RecursiveCall({"n": -1})
        call2 = RecursiveCall({"n": -1})
        t1 = Term(Var("x"), call1)
        t2 = Term(Const(2), call2)
        s = Sum([t1, t2])
        assert should_apply_cse(s) is True

    def test_two_calls_no_duplicates(self):
        """Test that two different calls don't trigger duplicate check."""
        t1 = Term(Const(1), RecursiveCall({"n": -1}))
        t2 = Term(Const(2), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])
        # Less than 3 calls and no duplicates
        assert should_apply_cse(s) is False


class TestOptimizerIntegration:
    """Integration tests for optimizer module."""

    def test_optimize_legendre_recurrence(self):
        """Test optimizing Legendre recurrence expression."""
        # (2n-1)*x*E[n-1] - (n-1)*E[n-2]
        t1 = Term(BinOp('*', IndexExpr("2*n-1"), Var("x")), RecursiveCall({"n": -1}))
        t2 = Term(IndexExpr("-(n-1)"), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])
        scale = IndexExpr("n")
        expr = ScaledExpr(s, scale, is_division=True)

        ctx = CodegenContext("LegendreCoeff", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx, enable_cse=True)
        result = opt.optimize_expression(expr)

        assert result.result_expr is not None

    def test_optimize_hermite_e_recurrence(self):
        """Test optimizing Hermite E coefficient recurrence."""
        # p2 * E[i-1, j, t-1] + PA * E[i-1, j, t] + (t+1) * E[i-1, j, t+1]
        t1 = Term(Var("p2"), RecursiveCall({"i": -1, "j": 0, "t": -1}))
        t2 = Term(Var("PA"), RecursiveCall({"i": -1, "j": 0, "t": 0}))
        t3 = Term(IndexExpr("t+1"), RecursiveCall({"i": -1, "j": 0, "t": 1}))
        s = Sum([t1, t2, t3])

        ctx = CodegenContext("HermiteECoeff", ["i", "j", "t"], ["PA", "PB", "p2"])
        opt = ExpressionOptimizer(ctx, enable_cse=True, cse_threshold=3)
        result = opt.optimize_expression(s)

        # With 3 calls, should apply CSE
        assert result.result_expr is not None

    def test_full_code_generation_pipeline(self):
        """Test complete code generation with optimization."""
        # Simple two-term recurrence: x * E[n-1] + E[n-2]
        t1 = Term(Var("x"), RecursiveCall({"n": -1}))
        t2 = Term(Const(1), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])

        ctx = CodegenContext("TestCoeff", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx, enable_cse=True, cse_threshold=2)
        gen = OptimizedCodeGenerator(ctx, opt)

        body = gen.generate_body(s)

        # Should generate valid C++ code
        assert "Vec8d" in body
        assert "return" in body
        assert "TestCoeff" in body


class TestOptimizerEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_empty_sum_optimization(self):
        """Test optimizing empty sum."""
        s = Sum([])
        ctx = CodegenContext("Test", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx)
        result = opt.optimize_expression(s)
        assert result.result_expr is not None

    def test_single_term_optimization(self):
        """Test optimizing single term."""
        term = Term(Var("x"), RecursiveCall({"n": -1}))
        s = Sum([term])
        ctx = CodegenContext("Test", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx)
        result = opt.optimize_expression(s)
        assert result.result_expr is not None

    def test_zero_threshold(self):
        """Test CSE with zero threshold."""
        ctx = CodegenContext("Test", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx, cse_threshold=0)
        call = RecursiveCall({"n": -1})
        result = opt.optimize_expression(call)
        assert result.result_expr is not None

    def test_high_threshold(self):
        """Test CSE with very high threshold."""
        ctx = CodegenContext("Test", ["n"], ["x"])
        opt = ExpressionOptimizer(ctx, cse_threshold=100)
        t1 = Term(Const(1), RecursiveCall({"n": -1}))
        t2 = Term(Const(2), RecursiveCall({"n": -2}))
        s = Sum([t1, t2])
        result = opt.optimize_expression(s)
        # With high threshold, should not apply CSE
        assert result.result_expr is s
