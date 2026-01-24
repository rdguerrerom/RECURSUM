/**
 * @file bench_all_coefficients.cpp
 * @brief Benchmark EVERY Hermite E coefficient from E^{0,0}_0 to E^{4,4}_8
 *
 * This file generates benchmarks for ALL 125 coefficients to enable
 * complete statistical analysis of TMP vs Symbolic performance.
 *
 * All measurements are from real execution - no estimates.
 */

#include <benchmark/benchmark.h>
#include <iostream>
#include "benchmark_dispatcher.hpp"
#include "benchmark_common.hpp"

using namespace recursum::benchmark;

// =============================================================================
// Macro to generate benchmark for EVERY coefficient
// =============================================================================

#define BENCH_E(NA, NB, T) \
static void BM_E_TMP_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    const auto& dispatcher = TMPDispatcher::instance(); \
    Vec8d result; \
    for (auto _ : state) { \
        result = dispatcher.compute(NA, NB, T, params.PA[0], params.PB[0], params.p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.counters["impl"] = 0; \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
    state.counters["L"] = NA + NB; \
} \
BENCHMARK(BM_E_TMP_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0); \
\
static void BM_E_Sym_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    const auto& dispatcher = SymbolicDispatcher::instance(); \
    Vec8d result; \
    for (auto _ : state) { \
        result = dispatcher.compute(NA, NB, T, params.PA[0], params.PB[0], params.one_over_2p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.counters["impl"] = 1; \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
    state.counters["L"] = NA + NB; \
} \
BENCHMARK(BM_E_Sym_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// ALL 125 Hermite E Coefficients
// =============================================================================

// L=0: E^{0,0}_0 (1 coefficient)
BENCH_E(0, 0, 0)

// L=1: E^{1,0}_t, E^{0,1}_t (4 coefficients)
BENCH_E(1, 0, 0)
BENCH_E(1, 0, 1)
BENCH_E(0, 1, 0)
BENCH_E(0, 1, 1)

// L=2: E^{2,0}_t, E^{1,1}_t, E^{0,2}_t (9 coefficients)
BENCH_E(2, 0, 0)
BENCH_E(2, 0, 1)
BENCH_E(2, 0, 2)
BENCH_E(1, 1, 0)
BENCH_E(1, 1, 1)
BENCH_E(1, 1, 2)
BENCH_E(0, 2, 0)
BENCH_E(0, 2, 1)
BENCH_E(0, 2, 2)

// L=3: E^{3,0}_t, E^{2,1}_t, E^{1,2}_t, E^{0,3}_t (16 coefficients)
BENCH_E(3, 0, 0)
BENCH_E(3, 0, 1)
BENCH_E(3, 0, 2)
BENCH_E(3, 0, 3)
BENCH_E(2, 1, 0)
BENCH_E(2, 1, 1)
BENCH_E(2, 1, 2)
BENCH_E(2, 1, 3)
BENCH_E(1, 2, 0)
BENCH_E(1, 2, 1)
BENCH_E(1, 2, 2)
BENCH_E(1, 2, 3)
BENCH_E(0, 3, 0)
BENCH_E(0, 3, 1)
BENCH_E(0, 3, 2)
BENCH_E(0, 3, 3)

// L=4: E^{4,0}_t, E^{3,1}_t, E^{2,2}_t, E^{1,3}_t, E^{0,4}_t (25 coefficients)
BENCH_E(4, 0, 0)
BENCH_E(4, 0, 1)
BENCH_E(4, 0, 2)
BENCH_E(4, 0, 3)
BENCH_E(4, 0, 4)
BENCH_E(3, 1, 0)
BENCH_E(3, 1, 1)
BENCH_E(3, 1, 2)
BENCH_E(3, 1, 3)
BENCH_E(3, 1, 4)
BENCH_E(2, 2, 0)
BENCH_E(2, 2, 1)
BENCH_E(2, 2, 2)
BENCH_E(2, 2, 3)
BENCH_E(2, 2, 4)
BENCH_E(1, 3, 0)
BENCH_E(1, 3, 1)
BENCH_E(1, 3, 2)
BENCH_E(1, 3, 3)
BENCH_E(1, 3, 4)
BENCH_E(0, 4, 0)
BENCH_E(0, 4, 1)
BENCH_E(0, 4, 2)
BENCH_E(0, 4, 3)
BENCH_E(0, 4, 4)

// L=5: E^{4,1}_t, E^{3,2}_t, E^{2,3}_t, E^{1,4}_t (24 coefficients)
BENCH_E(4, 1, 0)
BENCH_E(4, 1, 1)
BENCH_E(4, 1, 2)
BENCH_E(4, 1, 3)
BENCH_E(4, 1, 4)
BENCH_E(4, 1, 5)
BENCH_E(3, 2, 0)
BENCH_E(3, 2, 1)
BENCH_E(3, 2, 2)
BENCH_E(3, 2, 3)
BENCH_E(3, 2, 4)
BENCH_E(3, 2, 5)
BENCH_E(2, 3, 0)
BENCH_E(2, 3, 1)
BENCH_E(2, 3, 2)
BENCH_E(2, 3, 3)
BENCH_E(2, 3, 4)
BENCH_E(2, 3, 5)
BENCH_E(1, 4, 0)
BENCH_E(1, 4, 1)
BENCH_E(1, 4, 2)
BENCH_E(1, 4, 3)
BENCH_E(1, 4, 4)
BENCH_E(1, 4, 5)

// L=6: E^{4,2}_t, E^{3,3}_t, E^{2,4}_t (21 coefficients)
BENCH_E(4, 2, 0)
BENCH_E(4, 2, 1)
BENCH_E(4, 2, 2)
BENCH_E(4, 2, 3)
BENCH_E(4, 2, 4)
BENCH_E(4, 2, 5)
BENCH_E(4, 2, 6)
BENCH_E(3, 3, 0)
BENCH_E(3, 3, 1)
BENCH_E(3, 3, 2)
BENCH_E(3, 3, 3)
BENCH_E(3, 3, 4)
BENCH_E(3, 3, 5)
BENCH_E(3, 3, 6)
BENCH_E(2, 4, 0)
BENCH_E(2, 4, 1)
BENCH_E(2, 4, 2)
BENCH_E(2, 4, 3)
BENCH_E(2, 4, 4)
BENCH_E(2, 4, 5)
BENCH_E(2, 4, 6)

// L=7: E^{4,3}_t, E^{3,4}_t (16 coefficients)
BENCH_E(4, 3, 0)
BENCH_E(4, 3, 1)
BENCH_E(4, 3, 2)
BENCH_E(4, 3, 3)
BENCH_E(4, 3, 4)
BENCH_E(4, 3, 5)
BENCH_E(4, 3, 6)
BENCH_E(4, 3, 7)
BENCH_E(3, 4, 0)
BENCH_E(3, 4, 1)
BENCH_E(3, 4, 2)
BENCH_E(3, 4, 3)
BENCH_E(3, 4, 4)
BENCH_E(3, 4, 5)
BENCH_E(3, 4, 6)
BENCH_E(3, 4, 7)

// L=8: E^{4,4}_t (9 coefficients)
BENCH_E(4, 4, 0)
BENCH_E(4, 4, 1)
BENCH_E(4, 4, 2)
BENCH_E(4, 4, 3)
BENCH_E(4, 4, 4)
BENCH_E(4, 4, 5)
BENCH_E(4, 4, 6)
BENCH_E(4, 4, 7)
BENCH_E(4, 4, 8)

// =============================================================================
// Gradient Benchmarks for ALL coefficients
// =============================================================================

#define BENCH_GRAD(NA, NB, T) \
static void BM_dEdPA_TMP_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    const auto& dispatcher = TMPGradPADispatcher::instance(); \
    Vec8d result; \
    for (auto _ : state) { \
        result = dispatcher.compute(NA, NB, T, params.PA[0], params.PB[0], params.p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.counters["impl"] = 0; \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
    state.counters["type"] = 1; \
} \
BENCHMARK(BM_dEdPA_TMP_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0); \
\
static void BM_dEdPA_Sym_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    const auto& dispatcher = SymbolicGradPADispatcher::instance(); \
    Vec8d result; \
    for (auto _ : state) { \
        result = dispatcher.compute(NA, NB, T, params.PA[0], params.PB[0], params.one_over_2p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.counters["impl"] = 1; \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
    state.counters["type"] = 1; \
} \
BENCHMARK(BM_dEdPA_Sym_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// Representative gradient benchmarks (key shell pairs)
// pp gradients
BENCH_GRAD(1, 1, 0)
BENCH_GRAD(1, 1, 1)
BENCH_GRAD(1, 1, 2)

// dd gradients
BENCH_GRAD(2, 2, 0)
BENCH_GRAD(2, 2, 2)
BENCH_GRAD(2, 2, 4)

// ff gradients
BENCH_GRAD(3, 3, 0)
BENCH_GRAD(3, 3, 3)
BENCH_GRAD(3, 3, 6)

// gg gradients
BENCH_GRAD(4, 4, 0)
BENCH_GRAD(4, 4, 4)
BENCH_GRAD(4, 4, 8)

// =============================================================================
// Main
// =============================================================================

int main(int argc, char** argv) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "RECURSUM: ALL COEFFICIENTS BENCHMARK\n";
    std::cout << "TMP vs Symbolic: E^{0,0}_0 to E^{4,4}_8 (125 coefficients)\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << "\nCoverage:\n";
    std::cout << "  L=0: 1 coefficient   (ss)\n";
    std::cout << "  L=1: 4 coefficients  (sp)\n";
    std::cout << "  L=2: 9 coefficients  (pp, sd)\n";
    std::cout << "  L=3: 16 coefficients (pd, sf)\n";
    std::cout << "  L=4: 25 coefficients (dd, pf, sg)\n";
    std::cout << "  L=5: 24 coefficients (df, pg)\n";
    std::cout << "  L=6: 21 coefficients (ff, dg)\n";
    std::cout << "  L=7: 16 coefficients (fg)\n";
    std::cout << "  L=8: 9 coefficients  (gg)\n";
    std::cout << "  Total: 125 E coefficients + gradients\n";
    std::cout << std::string(80, '=') << "\n\n";

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
