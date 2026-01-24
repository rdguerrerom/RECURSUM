/**
 * @file bench_direct_comparison.cpp
 * @brief Direct benchmark WITHOUT dispatcher overhead
 *
 * This benchmark calls TMP and Symbolic functions directly (no std::function)
 * to allow full compiler optimization including CSE, inlining, etc.
 */

#include <benchmark/benchmark.h>
#include <iostream>

// Include TMP implementation directly
#include <recursum/mcmd/hermite_e.hpp>

// Include symbolic implementation directly
#include "hermite_e_symbolic.hpp"

#include "benchmark_common.hpp"

using namespace recursum::benchmark;
using namespace recursum::mcmd;

// =============================================================================
// DIRECT BENCHMARK MACROS (no dispatcher, no std::function)
// =============================================================================

#define DIRECT_BENCH_E(NA, NB, T) \
static void BM_Direct_TMP_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    Vec8d result; \
    for (auto _ : state) { \
        result = HermiteE<NA, NB, T>::compute(params.PA[0], params.PB[0], params.p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.counters["impl"] = 0; \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
    state.counters["direct"] = 1; \
} \
BENCHMARK(BM_Direct_TMP_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0); \
\
static void BM_Direct_Sym_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    Vec8d result; \
    for (auto _ : state) { \
        result = recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_##T(params.PA[0], params.PB[0], params.one_over_2p); \
        benchmark::DoNotOptimize(result); \
    } \
    state.counters["impl"] = 1; \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
    state.counters["direct"] = 1; \
} \
BENCHMARK(BM_Direct_Sym_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// Key coefficients for direct comparison
// =============================================================================

// L=0: ss
DIRECT_BENCH_E(0, 0, 0)

// L=1: sp
DIRECT_BENCH_E(1, 0, 0)
DIRECT_BENCH_E(1, 0, 1)
DIRECT_BENCH_E(0, 1, 0)
DIRECT_BENCH_E(0, 1, 1)

// L=2: pp, sd
DIRECT_BENCH_E(1, 1, 0)
DIRECT_BENCH_E(1, 1, 1)
DIRECT_BENCH_E(1, 1, 2)
DIRECT_BENCH_E(2, 0, 0)
DIRECT_BENCH_E(2, 0, 1)
DIRECT_BENCH_E(2, 0, 2)

// L=3: pd, sf
DIRECT_BENCH_E(2, 1, 0)
DIRECT_BENCH_E(2, 1, 1)
DIRECT_BENCH_E(2, 1, 2)
DIRECT_BENCH_E(2, 1, 3)
DIRECT_BENCH_E(3, 0, 0)
DIRECT_BENCH_E(3, 0, 1)
DIRECT_BENCH_E(3, 0, 2)
DIRECT_BENCH_E(3, 0, 3)

// L=4: dd, pf, sg
DIRECT_BENCH_E(2, 2, 0)
DIRECT_BENCH_E(2, 2, 1)
DIRECT_BENCH_E(2, 2, 2)
DIRECT_BENCH_E(2, 2, 3)
DIRECT_BENCH_E(2, 2, 4)
DIRECT_BENCH_E(3, 1, 0)
DIRECT_BENCH_E(3, 1, 2)
DIRECT_BENCH_E(3, 1, 4)
DIRECT_BENCH_E(4, 0, 0)
DIRECT_BENCH_E(4, 0, 2)
DIRECT_BENCH_E(4, 0, 4)

// L=5: df, pg
DIRECT_BENCH_E(3, 2, 0)
DIRECT_BENCH_E(3, 2, 2)
DIRECT_BENCH_E(3, 2, 5)
DIRECT_BENCH_E(4, 1, 0)
DIRECT_BENCH_E(4, 1, 2)
DIRECT_BENCH_E(4, 1, 5)

// L=6: ff, dg
DIRECT_BENCH_E(3, 3, 0)
DIRECT_BENCH_E(3, 3, 3)
DIRECT_BENCH_E(3, 3, 6)
DIRECT_BENCH_E(4, 2, 0)
DIRECT_BENCH_E(4, 2, 3)
DIRECT_BENCH_E(4, 2, 6)

// L=7: fg
DIRECT_BENCH_E(4, 3, 0)
DIRECT_BENCH_E(4, 3, 3)
DIRECT_BENCH_E(4, 3, 7)
DIRECT_BENCH_E(3, 4, 0)
DIRECT_BENCH_E(3, 4, 3)
DIRECT_BENCH_E(3, 4, 7)

// L=8: gg (FULL coverage)
DIRECT_BENCH_E(4, 4, 0)
DIRECT_BENCH_E(4, 4, 1)
DIRECT_BENCH_E(4, 4, 2)
DIRECT_BENCH_E(4, 4, 3)
DIRECT_BENCH_E(4, 4, 4)
DIRECT_BENCH_E(4, 4, 5)
DIRECT_BENCH_E(4, 4, 6)
DIRECT_BENCH_E(4, 4, 7)
DIRECT_BENCH_E(4, 4, 8)

// =============================================================================
// Main
// =============================================================================

int main(int argc, char** argv) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "RECURSUM: DIRECT COMPARISON BENCHMARK (NO DISPATCHER OVERHEAD)\n";
    std::cout << "TMP vs Symbolic with FULL compiler optimization\n";
    std::cout << std::string(80, '=') << "\n\n";
    std::cout << "This benchmark bypasses std::function to allow:\n";
    std::cout << "  - Full inlining\n";
    std::cout << "  - Cross-boundary CSE\n";
    std::cout << "  - Loop unrolling\n";
    std::cout << "  - Register allocation optimization\n";
    std::cout << std::string(80, '=') << "\n\n";

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
