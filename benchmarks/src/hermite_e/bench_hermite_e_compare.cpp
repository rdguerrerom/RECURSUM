/**
 * @file bench_hermite_e_compare.cpp
 * @brief Side-by-side comparison benchmarks for TMP vs Symbolic vs Naive
 *
 * Uses dispatcher pattern (following McMD/hermite_dispatcher.hpp) to prevent
 * compiler from over-optimizing benchmark loops. This ensures fair comparison
 * between implementations by preventing constant folding and dead code elimination.
 *
 * IMPORTANT: All measurements are from real execution - no estimates or theoretical values.
 */

#include <benchmark/benchmark.h>
#include <iostream>
#include "benchmark_dispatcher.hpp"
#include "benchmark_common.hpp"

using namespace recursum::benchmark;

// =============================================================================
// Dispatched Comparison: TMP vs Symbolic (fair comparison via dispatcher)
// =============================================================================

// E^{3,3}_6 - highest complexity coefficient for ff shell pair
static void BM_Dispatched_TMP_3_3_6(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = TMPDispatcher::instance();
    Vec8d result;
    for (auto _ : state) {
        result = dispatcher.compute(3, 3, 6, params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel("TMP");
    state.counters["impl"] = 0;
    state.counters["nA"] = 3;
    state.counters["nB"] = 3;
    state.counters["t"] = 6;
}
BENCHMARK(BM_Dispatched_TMP_3_3_6)->Unit(benchmark::kNanosecond)->MinTime(2.0);

static void BM_Dispatched_Symbolic_3_3_6(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = SymbolicDispatcher::instance();
    Vec8d result;
    for (auto _ : state) {
        result = dispatcher.compute(3, 3, 6, params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel("Symbolic");
    state.counters["impl"] = 1;
    state.counters["nA"] = 3;
    state.counters["nB"] = 3;
    state.counters["t"] = 6;
}
BENCHMARK(BM_Dispatched_Symbolic_3_3_6)->Unit(benchmark::kNanosecond)->MinTime(2.0);

// E^{2,2}_4 - dd shell pair
static void BM_Dispatched_TMP_2_2_4(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = TMPDispatcher::instance();
    Vec8d result;
    for (auto _ : state) {
        result = dispatcher.compute(2, 2, 4, params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel("TMP");
    state.counters["impl"] = 0;
    state.counters["nA"] = 2;
    state.counters["nB"] = 2;
    state.counters["t"] = 4;
}
BENCHMARK(BM_Dispatched_TMP_2_2_4)->Unit(benchmark::kNanosecond)->MinTime(2.0);

static void BM_Dispatched_Symbolic_2_2_4(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = SymbolicDispatcher::instance();
    Vec8d result;
    for (auto _ : state) {
        result = dispatcher.compute(2, 2, 4, params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel("Symbolic");
    state.counters["impl"] = 1;
    state.counters["nA"] = 2;
    state.counters["nB"] = 2;
    state.counters["t"] = 4;
}
BENCHMARK(BM_Dispatched_Symbolic_2_2_4)->Unit(benchmark::kNanosecond)->MinTime(2.0);

// E^{1,1}_2 - pp shell pair
static void BM_Dispatched_TMP_1_1_2(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = TMPDispatcher::instance();
    Vec8d result;
    for (auto _ : state) {
        result = dispatcher.compute(1, 1, 2, params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel("TMP");
    state.counters["impl"] = 0;
    state.counters["nA"] = 1;
    state.counters["nB"] = 1;
    state.counters["t"] = 2;
}
BENCHMARK(BM_Dispatched_TMP_1_1_2)->Unit(benchmark::kNanosecond)->MinTime(2.0);

static void BM_Dispatched_Symbolic_1_1_2(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = SymbolicDispatcher::instance();
    Vec8d result;
    for (auto _ : state) {
        result = dispatcher.compute(1, 1, 2, params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel("Symbolic");
    state.counters["impl"] = 1;
    state.counters["nA"] = 1;
    state.counters["nB"] = 1;
    state.counters["t"] = 2;
}
BENCHMARK(BM_Dispatched_Symbolic_1_1_2)->Unit(benchmark::kNanosecond)->MinTime(2.0);

// E^{0,0}_0 - base case (ss shell pair)
static void BM_Dispatched_TMP_0_0_0(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = TMPDispatcher::instance();
    Vec8d result;
    for (auto _ : state) {
        result = dispatcher.compute(0, 0, 0, params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel("TMP");
    state.counters["impl"] = 0;
    state.counters["nA"] = 0;
    state.counters["nB"] = 0;
    state.counters["t"] = 0;
}
BENCHMARK(BM_Dispatched_TMP_0_0_0)->Unit(benchmark::kNanosecond)->MinTime(2.0);

static void BM_Dispatched_Symbolic_0_0_0(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = SymbolicDispatcher::instance();
    Vec8d result;
    for (auto _ : state) {
        result = dispatcher.compute(0, 0, 0, params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel("Symbolic");
    state.counters["impl"] = 1;
    state.counters["nA"] = 0;
    state.counters["nB"] = 0;
    state.counters["t"] = 0;
}
BENCHMARK(BM_Dispatched_Symbolic_0_0_0)->Unit(benchmark::kNanosecond)->MinTime(2.0);

// =============================================================================
// Scaling Analysis: Measure how performance changes with angular momentum
// Using dispatchers for fair comparison
// =============================================================================

static void BM_Scaling_TMP_L0(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = TMPDispatcher::instance();
    for (auto _ : state) {
        auto r = dispatcher.compute(0, 0, 0, params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(r);
    }
    state.counters["L_total"] = 0;
    state.counters["impl"] = 0;
}
BENCHMARK(BM_Scaling_TMP_L0)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Scaling_TMP_L2(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = TMPDispatcher::instance();
    for (auto _ : state) {
        auto r = dispatcher.compute(1, 1, 2, params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(r);
    }
    state.counters["L_total"] = 2;
    state.counters["impl"] = 0;
}
BENCHMARK(BM_Scaling_TMP_L2)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Scaling_TMP_L4(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = TMPDispatcher::instance();
    for (auto _ : state) {
        auto r = dispatcher.compute(2, 2, 4, params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(r);
    }
    state.counters["L_total"] = 4;
    state.counters["impl"] = 0;
}
BENCHMARK(BM_Scaling_TMP_L4)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Scaling_TMP_L6(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = TMPDispatcher::instance();
    for (auto _ : state) {
        auto r = dispatcher.compute(3, 3, 6, params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(r);
    }
    state.counters["L_total"] = 6;
    state.counters["impl"] = 0;
}
BENCHMARK(BM_Scaling_TMP_L6)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Scaling_Symbolic_L0(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = SymbolicDispatcher::instance();
    for (auto _ : state) {
        auto r = dispatcher.compute(0, 0, 0, params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(r);
    }
    state.counters["L_total"] = 0;
    state.counters["impl"] = 1;
}
BENCHMARK(BM_Scaling_Symbolic_L0)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Scaling_Symbolic_L2(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = SymbolicDispatcher::instance();
    for (auto _ : state) {
        auto r = dispatcher.compute(1, 1, 2, params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(r);
    }
    state.counters["L_total"] = 2;
    state.counters["impl"] = 1;
}
BENCHMARK(BM_Scaling_Symbolic_L2)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Scaling_Symbolic_L4(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = SymbolicDispatcher::instance();
    for (auto _ : state) {
        auto r = dispatcher.compute(2, 2, 4, params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(r);
    }
    state.counters["L_total"] = 4;
    state.counters["impl"] = 1;
}
BENCHMARK(BM_Scaling_Symbolic_L4)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Scaling_Symbolic_L6(benchmark::State& state) {
    GeometryParams params(42);
    const auto& dispatcher = SymbolicDispatcher::instance();
    for (auto _ : state) {
        auto r = dispatcher.compute(3, 3, 6, params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(r);
    }
    state.counters["L_total"] = 6;
    state.counters["impl"] = 1;
}
BENCHMARK(BM_Scaling_Symbolic_L6)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// Comprehensive coefficient comparison (all coefficients)
// =============================================================================

// Macro to generate dispatched benchmarks for any coefficient
#define DEFINE_DISPATCHED_BENCHMARK(NA, NB, T) \
static void BM_Compare_TMP_##NA##_##NB##_##T(benchmark::State& state) { \
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
BENCHMARK(BM_Compare_TMP_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(2.0); \
\
static void BM_Compare_Symbolic_##NA##_##NB##_##T(benchmark::State& state) { \
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
BENCHMARK(BM_Compare_Symbolic_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(2.0);

// Additional comparison points
DEFINE_DISPATCHED_BENCHMARK(1, 0, 0)
DEFINE_DISPATCHED_BENCHMARK(1, 0, 1)
DEFINE_DISPATCHED_BENCHMARK(0, 1, 0)
DEFINE_DISPATCHED_BENCHMARK(0, 1, 1)
DEFINE_DISPATCHED_BENCHMARK(1, 1, 0)
DEFINE_DISPATCHED_BENCHMARK(1, 1, 1)
DEFINE_DISPATCHED_BENCHMARK(2, 0, 0)
DEFINE_DISPATCHED_BENCHMARK(2, 0, 1)
DEFINE_DISPATCHED_BENCHMARK(2, 0, 2)
DEFINE_DISPATCHED_BENCHMARK(0, 2, 0)
DEFINE_DISPATCHED_BENCHMARK(0, 2, 1)
DEFINE_DISPATCHED_BENCHMARK(0, 2, 2)
DEFINE_DISPATCHED_BENCHMARK(2, 2, 0)
DEFINE_DISPATCHED_BENCHMARK(2, 2, 1)
DEFINE_DISPATCHED_BENCHMARK(2, 2, 2)
DEFINE_DISPATCHED_BENCHMARK(2, 2, 3)
DEFINE_DISPATCHED_BENCHMARK(3, 0, 0)
DEFINE_DISPATCHED_BENCHMARK(3, 0, 1)
DEFINE_DISPATCHED_BENCHMARK(3, 0, 2)
DEFINE_DISPATCHED_BENCHMARK(3, 0, 3)
DEFINE_DISPATCHED_BENCHMARK(0, 3, 0)
DEFINE_DISPATCHED_BENCHMARK(0, 3, 1)
DEFINE_DISPATCHED_BENCHMARK(0, 3, 2)
DEFINE_DISPATCHED_BENCHMARK(0, 3, 3)
DEFINE_DISPATCHED_BENCHMARK(3, 3, 0)
DEFINE_DISPATCHED_BENCHMARK(3, 3, 1)
DEFINE_DISPATCHED_BENCHMARK(3, 3, 2)
DEFINE_DISPATCHED_BENCHMARK(3, 3, 3)
DEFINE_DISPATCHED_BENCHMARK(3, 3, 4)
DEFINE_DISPATCHED_BENCHMARK(3, 3, 5)

// =============================================================================
// Main function with system info reporting
// =============================================================================

int main(int argc, char** argv) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "RECURSUM BENCHMARK: TMP vs Symbolic Hermite E Coefficients\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << "Implementation Types:\n";
    std::cout << "  impl=0: TMP (Template Metaprogramming via dispatcher)\n";
    std::cout << "  impl=1: Symbolic (SymPy-generated via dispatcher)\n";
    std::cout << "\n";
    std::cout << "Methodology:\n";
    std::cout << "  - Both implementations use std::function dispatchers\n";
    std::cout << "  - This prevents compiler over-optimization (fair comparison)\n";
    std::cout << "  - Pattern follows McMD/hermite_dispatcher.hpp\n";
    std::cout << std::string(80, '=') << "\n\n";

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "BENCHMARK COMPLETE\n";
    std::cout << "Results saved to JSON for analysis.\n";
    std::cout << std::string(80, '=') << "\n";

    return 0;
}
