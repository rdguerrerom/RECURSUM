/**
 * @file bench_layered_comparison.cpp
 * @brief Compare Original TMP vs Layered TMP vs Symbolic
 *
 * Tests whether the layered (CSE-enabled) approach closes the gap with symbolic.
 */

#include <benchmark/benchmark.h>
#include <iostream>

// Original TMP
#include <recursum/mcmd/hermite_e.hpp>

// Layered TMP with true CSE
#include <recursum/mcmd/hermite_e_layered.hpp>

// Symbolic
#include "hermite_e_symbolic.hpp"

#include "benchmark_common.hpp"

using namespace recursum::benchmark;

// =============================================================================
// Benchmark Macros
// =============================================================================

#define BENCH_ALL(NA, NB, T) \
/* Original TMP */ \
static void BM_Orig_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    Vec8d result; \
    for (auto _ : state) { \
        result = recursum::mcmd::HermiteE<NA, NB, T>::compute(params.PA[0], params.PB[0], params.p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.counters["impl"] = 0; \
} \
BENCHMARK(BM_Orig_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0); \
\
/* Layered TMP with CSE */ \
static void BM_Layer_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    Vec8d result; \
    for (auto _ : state) { \
        result = recursum::mcmd::layered::HermiteEOpt<NA, NB, T>::compute(params.PA[0], params.PB[0], params.p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.counters["impl"] = 2; \
} \
BENCHMARK(BM_Layer_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0); \
\
/* Symbolic */ \
static void BM_Sym_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    Vec8d result; \
    for (auto _ : state) { \
        result = recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_##T(params.PA[0], params.PB[0], params.one_over_2p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.counters["impl"] = 1; \
} \
BENCHMARK(BM_Sym_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// Key Coefficients to Compare
// =============================================================================

// Low angular momentum (should be similar)
BENCH_ALL(0, 0, 0)
BENCH_ALL(1, 0, 0)
BENCH_ALL(1, 0, 1)
BENCH_ALL(0, 1, 0)
BENCH_ALL(0, 1, 1)

// Medium angular momentum
BENCH_ALL(1, 1, 0)
BENCH_ALL(1, 1, 1)
BENCH_ALL(1, 1, 2)
BENCH_ALL(2, 0, 0)
BENCH_ALL(2, 0, 2)
BENCH_ALL(0, 2, 0)
BENCH_ALL(0, 2, 2)

// Higher angular momentum (dd)
BENCH_ALL(2, 2, 0)
BENCH_ALL(2, 2, 2)
BENCH_ALL(2, 2, 4)

// ff shell pairs
BENCH_ALL(3, 3, 0)
BENCH_ALL(3, 3, 3)
BENCH_ALL(3, 3, 6)

// gg shell pairs - where the difference should be largest
BENCH_ALL(4, 4, 0)
BENCH_ALL(4, 4, 4)
BENCH_ALL(4, 4, 8)

// =============================================================================
// Main
// =============================================================================

int main(int argc, char** argv) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "RECURSUM: LAYERED CSE COMPARISON\n";
    std::cout << "Original TMP vs Layered TMP (CSE) vs Symbolic\n";
    std::cout << std::string(80, '=') << "\n\n";
    std::cout << "impl=0: Original TMP (recursive)\n";
    std::cout << "impl=2: Layered TMP (CSE, compute all t together)\n";
    std::cout << "impl=1: Symbolic (SymPy-generated)\n";
    std::cout << std::string(80, '=') << "\n\n";

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
