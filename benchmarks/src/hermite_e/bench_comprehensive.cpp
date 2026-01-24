/**
 * @file bench_comprehensive.cpp
 * @brief Comprehensive benchmarks: ALL shell pairs from (ss) to (gg)
 *
 * Coverage:
 * - Hermite E coefficients: E^{0,0}_0 to E^{4,4}_8 (125 coefficients)
 * - Hermite gradients: ∂E/∂PA (125 gradients)
 * - Shell pairs: ss, sp, pp, sd, pd, dd, sf, pf, df, ff, sg, pg, dg, fg, gg
 *
 * All measurements are from real execution.
 */

#include <benchmark/benchmark.h>
#include <iostream>
#include <vector>
#include <tuple>
#include "benchmark_dispatcher.hpp"
#include "benchmark_common.hpp"

using namespace recursum::benchmark;

// =============================================================================
// Shell Pair Benchmarks - E Coefficients
// =============================================================================

// Macro to define benchmark for a specific shell pair
#define DEFINE_SHELL_PAIR_BENCH(LA, LB, SHELL_NAME) \
static void BM_ShellPair_TMP_##SHELL_NAME(benchmark::State& state) { \
    GeometryParams params(42); \
    const auto& dispatcher = TMPDispatcher::instance(); \
    Vec8d result = Vec8d(0.0); \
    for (auto _ : state) { \
        for (int nA = 0; nA <= LA; ++nA) { \
            for (int nB = 0; nB <= LB; ++nB) { \
                for (int t = 0; t <= nA + nB; ++t) { \
                    result += dispatcher.compute(nA, nB, t, params.PA[0], params.PB[0], params.p); \
                } \
            } \
        } \
        benchmark::DoNotOptimize(result); \
    } \
    state.SetLabel("TMP"); \
    state.counters["impl"] = 0; \
    state.counters["LA"] = LA; \
    state.counters["LB"] = LB; \
} \
BENCHMARK(BM_ShellPair_TMP_##SHELL_NAME)->Unit(benchmark::kNanosecond)->MinTime(2.0); \
\
static void BM_ShellPair_Symbolic_##SHELL_NAME(benchmark::State& state) { \
    GeometryParams params(42); \
    const auto& dispatcher = SymbolicDispatcher::instance(); \
    Vec8d result = Vec8d(0.0); \
    for (auto _ : state) { \
        for (int nA = 0; nA <= LA; ++nA) { \
            for (int nB = 0; nB <= LB; ++nB) { \
                for (int t = 0; t <= nA + nB; ++t) { \
                    result += dispatcher.compute(nA, nB, t, params.PA[0], params.PB[0], params.one_over_2p); \
                } \
            } \
        } \
        benchmark::DoNotOptimize(result); \
    } \
    state.SetLabel("Symbolic"); \
    state.counters["impl"] = 1; \
    state.counters["LA"] = LA; \
    state.counters["LB"] = LB; \
} \
BENCHMARK(BM_ShellPair_Symbolic_##SHELL_NAME)->Unit(benchmark::kNanosecond)->MinTime(2.0);

// Define benchmarks for all shell pairs
DEFINE_SHELL_PAIR_BENCH(0, 0, ss)
DEFINE_SHELL_PAIR_BENCH(0, 1, sp)
DEFINE_SHELL_PAIR_BENCH(1, 1, pp)
DEFINE_SHELL_PAIR_BENCH(0, 2, sd)
DEFINE_SHELL_PAIR_BENCH(1, 2, pd)
DEFINE_SHELL_PAIR_BENCH(2, 2, dd)
DEFINE_SHELL_PAIR_BENCH(0, 3, sf)
DEFINE_SHELL_PAIR_BENCH(1, 3, pf)
DEFINE_SHELL_PAIR_BENCH(2, 3, df)
DEFINE_SHELL_PAIR_BENCH(3, 3, ff)
DEFINE_SHELL_PAIR_BENCH(0, 4, sg)
DEFINE_SHELL_PAIR_BENCH(1, 4, pg)
DEFINE_SHELL_PAIR_BENCH(2, 4, dg)
DEFINE_SHELL_PAIR_BENCH(3, 4, fg)
DEFINE_SHELL_PAIR_BENCH(4, 4, gg)

// =============================================================================
// Gradient Benchmarks
// =============================================================================

#define DEFINE_GRAD_BENCH(LA, LB, SHELL_NAME) \
static void BM_GradPA_TMP_##SHELL_NAME(benchmark::State& state) { \
    GeometryParams params(42); \
    const auto& dispatcher = TMPGradPADispatcher::instance(); \
    Vec8d result = Vec8d(0.0); \
    for (auto _ : state) { \
        for (int nA = 0; nA <= LA; ++nA) { \
            for (int nB = 0; nB <= LB; ++nB) { \
                for (int t = 0; t <= nA + nB; ++t) { \
                    result += dispatcher.compute(nA, nB, t, params.PA[0], params.PB[0], params.p); \
                } \
            } \
        } \
        benchmark::DoNotOptimize(result); \
    } \
    state.SetLabel("TMP dE/dPA"); \
    state.counters["impl"] = 0; \
    state.counters["LA"] = LA; \
    state.counters["LB"] = LB; \
} \
BENCHMARK(BM_GradPA_TMP_##SHELL_NAME)->Unit(benchmark::kNanosecond)->MinTime(2.0); \
\
static void BM_GradPA_Symbolic_##SHELL_NAME(benchmark::State& state) { \
    GeometryParams params(42); \
    const auto& dispatcher = SymbolicGradPADispatcher::instance(); \
    Vec8d result = Vec8d(0.0); \
    for (auto _ : state) { \
        for (int nA = 0; nA <= LA; ++nA) { \
            for (int nB = 0; nB <= LB; ++nB) { \
                for (int t = 0; t <= nA + nB; ++t) { \
                    result += dispatcher.compute(nA, nB, t, params.PA[0], params.PB[0], params.one_over_2p); \
                } \
            } \
        } \
        benchmark::DoNotOptimize(result); \
    } \
    state.SetLabel("Symbolic dE/dPA"); \
    state.counters["impl"] = 1; \
    state.counters["LA"] = LA; \
    state.counters["LB"] = LB; \
} \
BENCHMARK(BM_GradPA_Symbolic_##SHELL_NAME)->Unit(benchmark::kNanosecond)->MinTime(2.0);

// Define gradient benchmarks for key shell pairs
DEFINE_GRAD_BENCH(1, 1, pp)
DEFINE_GRAD_BENCH(2, 2, dd)
DEFINE_GRAD_BENCH(3, 3, ff)
DEFINE_GRAD_BENCH(4, 4, gg)

// =============================================================================
// Individual Coefficient Benchmarks (for detailed analysis)
// =============================================================================

// Generate benchmarks for representative coefficients at each shell pair level
#define DEFINE_COEFF_BENCH(NA, NB, T) \
static void BM_Coeff_TMP_E_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    const auto& dispatcher = TMPDispatcher::instance(); \
    Vec8d result; \
    for (auto _ : state) { \
        result = dispatcher.compute(NA, NB, T, params.PA[0], params.PB[0], params.p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.SetLabel("TMP"); \
    state.counters["impl"] = 0; \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
} \
BENCHMARK(BM_Coeff_TMP_E_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(2.0); \
\
static void BM_Coeff_Symbolic_E_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    const auto& dispatcher = SymbolicDispatcher::instance(); \
    Vec8d result; \
    for (auto _ : state) { \
        result = dispatcher.compute(NA, NB, T, params.PA[0], params.PB[0], params.one_over_2p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.SetLabel("Symbolic"); \
    state.counters["impl"] = 1; \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
} \
BENCHMARK(BM_Coeff_Symbolic_E_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(2.0);

// Key coefficients for gg shell pair
DEFINE_COEFF_BENCH(4, 4, 0)
DEFINE_COEFF_BENCH(4, 4, 1)
DEFINE_COEFF_BENCH(4, 4, 2)
DEFINE_COEFF_BENCH(4, 4, 3)
DEFINE_COEFF_BENCH(4, 4, 4)
DEFINE_COEFF_BENCH(4, 4, 5)
DEFINE_COEFF_BENCH(4, 4, 6)
DEFINE_COEFF_BENCH(4, 4, 7)
DEFINE_COEFF_BENCH(4, 4, 8)

// Key coefficients for fg shell pair
DEFINE_COEFF_BENCH(3, 4, 0)
DEFINE_COEFF_BENCH(3, 4, 3)
DEFINE_COEFF_BENCH(3, 4, 7)

// Key coefficients for dg shell pair
DEFINE_COEFF_BENCH(2, 4, 0)
DEFINE_COEFF_BENCH(2, 4, 3)
DEFINE_COEFF_BENCH(2, 4, 6)

// =============================================================================
// Main
// =============================================================================

int main(int argc, char** argv) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "RECURSUM COMPREHENSIVE BENCHMARK\n";
    std::cout << "TMP vs Symbolic: All Shell Pairs (ss to gg)\n";
    std::cout << std::string(80, '=') << "\n\n";

    std::cout << "Coverage:\n";
    std::cout << "  - Hermite E coefficients: E^{0,0}_0 to E^{4,4}_8\n";
    std::cout << "  - Hermite gradients: dE/dPA\n";
    std::cout << "  - Shell pairs: ss, sp, pp, sd, pd, dd, sf, pf, df, ff, sg, pg, dg, fg, gg\n";
    std::cout << std::string(80, '=') << "\n\n";

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "COMPREHENSIVE BENCHMARK COMPLETE\n";
    std::cout << std::string(80, '=') << "\n";

    return 0;
}
