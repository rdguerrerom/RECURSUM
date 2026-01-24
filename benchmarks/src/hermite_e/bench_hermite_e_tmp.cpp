/**
 * @file bench_hermite_e_tmp.cpp
 * @brief Benchmarks for RECURSUM TMP (Template Metaprogramming) Hermite E coefficients
 *
 * Tests the compile-time evaluated Hermite expansion coefficients using
 * SFINAE-enabled template specializations from include/recursum/mcmd/hermite_e.hpp
 */

#include <benchmark/benchmark.h>
#include <recursum/mcmd/hermite_e.hpp>
#include "benchmark_common.hpp"

using namespace recursum::mcmd;
using namespace recursum::benchmark;

// =============================================================================
// Individual Coefficient Benchmarks (TMP Implementation)
// =============================================================================

/**
 * Macro to generate benchmark for a specific (nA, nB, t) coefficient
 * Uses TMP HermiteE<nA,nB,t>::compute()
 */
#define BENCH_HERMITE_E_TMP(NA, NB, T) \
static void BM_HermiteE_TMP_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    Vec8d result; \
    int64_t total_geometries = 0; \
    for (auto _ : state) { \
        result = HermiteE<NA, NB, T>::compute( \
            params.PA[0], params.PB[0], params.p); \
        benchmark::DoNotOptimize(result); \
        total_geometries += 8; \
    } \
    state.counters["Geometries/s"] = benchmark::Counter( \
        total_geometries, benchmark::Counter::kIsRate); \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
    state.counters["impl"] = 0; \
} \
BENCHMARK(BM_HermiteE_TMP_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0)

// =============================================================================
// Base Cases (s-type)
// =============================================================================
BENCH_HERMITE_E_TMP(0, 0, 0);  // E^{0,0}_0 = 1

// =============================================================================
// p-type coefficients (L=1)
// =============================================================================
BENCH_HERMITE_E_TMP(1, 0, 0);
BENCH_HERMITE_E_TMP(1, 0, 1);
BENCH_HERMITE_E_TMP(0, 1, 0);
BENCH_HERMITE_E_TMP(0, 1, 1);
BENCH_HERMITE_E_TMP(1, 1, 0);
BENCH_HERMITE_E_TMP(1, 1, 1);
BENCH_HERMITE_E_TMP(1, 1, 2);

// =============================================================================
// d-type coefficients (L=2)
// =============================================================================
BENCH_HERMITE_E_TMP(2, 0, 0);
BENCH_HERMITE_E_TMP(2, 0, 1);
BENCH_HERMITE_E_TMP(2, 0, 2);
BENCH_HERMITE_E_TMP(0, 2, 0);
BENCH_HERMITE_E_TMP(0, 2, 1);
BENCH_HERMITE_E_TMP(0, 2, 2);
BENCH_HERMITE_E_TMP(2, 2, 0);
BENCH_HERMITE_E_TMP(2, 2, 2);
BENCH_HERMITE_E_TMP(2, 2, 4);

// =============================================================================
// f-type coefficients (L=3)
// =============================================================================
BENCH_HERMITE_E_TMP(3, 0, 0);
BENCH_HERMITE_E_TMP(3, 0, 1);
BENCH_HERMITE_E_TMP(3, 0, 2);
BENCH_HERMITE_E_TMP(3, 0, 3);
BENCH_HERMITE_E_TMP(0, 3, 0);
BENCH_HERMITE_E_TMP(0, 3, 3);
BENCH_HERMITE_E_TMP(3, 3, 0);
BENCH_HERMITE_E_TMP(3, 3, 3);
BENCH_HERMITE_E_TMP(3, 3, 6);

// =============================================================================
// Shell Pair Benchmarks (all coefficients for a shell pair)
// =============================================================================

/**
 * @brief Benchmark all Hermite E coefficients for ss shell pair
 */
static void BM_HermiteE_TMP_ShellPair_ss(benchmark::State& state) {
    GeometryParams params(42);
    int64_t total_geometries = 0;

    for (auto _ : state) {
        // ss: Only E^{0,0}_0
        auto e_000 = HermiteE<0, 0, 0>::compute(params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(e_000);
        total_geometries += 8;
    }

    state.counters["Geometries/s"] = benchmark::Counter(total_geometries, benchmark::Counter::kIsRate);
    state.counters["LA"] = 0;
    state.counters["LB"] = 0;
    state.counters["num_coeffs"] = 1;
}
BENCHMARK(BM_HermiteE_TMP_ShellPair_ss)->Unit(benchmark::kNanosecond)->MinTime(1.0);

/**
 * @brief Benchmark all Hermite E coefficients for pp shell pair
 */
static void BM_HermiteE_TMP_ShellPair_pp(benchmark::State& state) {
    GeometryParams params(42);
    int64_t total_geometries = 0;

    for (auto _ : state) {
        // pp: All E^{nA,nB}_t for nA,nB in {0,1}, t in {0..nA+nB}
        auto e_000 = HermiteE<0, 0, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_100 = HermiteE<1, 0, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_101 = HermiteE<1, 0, 1>::compute(params.PA[0], params.PB[0], params.p);
        auto e_010 = HermiteE<0, 1, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_011 = HermiteE<0, 1, 1>::compute(params.PA[0], params.PB[0], params.p);
        auto e_110 = HermiteE<1, 1, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_111 = HermiteE<1, 1, 1>::compute(params.PA[0], params.PB[0], params.p);
        auto e_112 = HermiteE<1, 1, 2>::compute(params.PA[0], params.PB[0], params.p);

        benchmark::DoNotOptimize(e_000);
        benchmark::DoNotOptimize(e_100);
        benchmark::DoNotOptimize(e_101);
        benchmark::DoNotOptimize(e_010);
        benchmark::DoNotOptimize(e_011);
        benchmark::DoNotOptimize(e_110);
        benchmark::DoNotOptimize(e_111);
        benchmark::DoNotOptimize(e_112);
        total_geometries += 8;
    }

    state.counters["Geometries/s"] = benchmark::Counter(total_geometries, benchmark::Counter::kIsRate);
    state.counters["LA"] = 1;
    state.counters["LB"] = 1;
    state.counters["num_coeffs"] = 8;
}
BENCHMARK(BM_HermiteE_TMP_ShellPair_pp)->Unit(benchmark::kNanosecond)->MinTime(1.0);

/**
 * @brief Benchmark all Hermite E coefficients for dd shell pair
 */
static void BM_HermiteE_TMP_ShellPair_dd(benchmark::State& state) {
    GeometryParams params(42);
    int64_t total_geometries = 0;

    for (auto _ : state) {
        // dd: E^{nA,nB}_t for nA,nB in {0,1,2}, t in {0..nA+nB}
        // Base
        auto e_000 = HermiteE<0, 0, 0>::compute(params.PA[0], params.PB[0], params.p);
        // nA=1
        auto e_100 = HermiteE<1, 0, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_101 = HermiteE<1, 0, 1>::compute(params.PA[0], params.PB[0], params.p);
        // nA=2
        auto e_200 = HermiteE<2, 0, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_201 = HermiteE<2, 0, 1>::compute(params.PA[0], params.PB[0], params.p);
        auto e_202 = HermiteE<2, 0, 2>::compute(params.PA[0], params.PB[0], params.p);
        // nB=1
        auto e_010 = HermiteE<0, 1, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_011 = HermiteE<0, 1, 1>::compute(params.PA[0], params.PB[0], params.p);
        // nB=2
        auto e_020 = HermiteE<0, 2, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_021 = HermiteE<0, 2, 1>::compute(params.PA[0], params.PB[0], params.p);
        auto e_022 = HermiteE<0, 2, 2>::compute(params.PA[0], params.PB[0], params.p);
        // Mixed
        auto e_110 = HermiteE<1, 1, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_111 = HermiteE<1, 1, 1>::compute(params.PA[0], params.PB[0], params.p);
        auto e_112 = HermiteE<1, 1, 2>::compute(params.PA[0], params.PB[0], params.p);
        auto e_120 = HermiteE<1, 2, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_121 = HermiteE<1, 2, 1>::compute(params.PA[0], params.PB[0], params.p);
        auto e_122 = HermiteE<1, 2, 2>::compute(params.PA[0], params.PB[0], params.p);
        auto e_123 = HermiteE<1, 2, 3>::compute(params.PA[0], params.PB[0], params.p);
        auto e_210 = HermiteE<2, 1, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_211 = HermiteE<2, 1, 1>::compute(params.PA[0], params.PB[0], params.p);
        auto e_212 = HermiteE<2, 1, 2>::compute(params.PA[0], params.PB[0], params.p);
        auto e_213 = HermiteE<2, 1, 3>::compute(params.PA[0], params.PB[0], params.p);
        auto e_220 = HermiteE<2, 2, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_221 = HermiteE<2, 2, 1>::compute(params.PA[0], params.PB[0], params.p);
        auto e_222 = HermiteE<2, 2, 2>::compute(params.PA[0], params.PB[0], params.p);
        auto e_223 = HermiteE<2, 2, 3>::compute(params.PA[0], params.PB[0], params.p);
        auto e_224 = HermiteE<2, 2, 4>::compute(params.PA[0], params.PB[0], params.p);

        benchmark::DoNotOptimize(e_000);
        benchmark::DoNotOptimize(e_100);
        benchmark::DoNotOptimize(e_200);
        benchmark::DoNotOptimize(e_220);
        benchmark::DoNotOptimize(e_224);
        total_geometries += 8;
    }

    state.counters["Geometries/s"] = benchmark::Counter(total_geometries, benchmark::Counter::kIsRate);
    state.counters["LA"] = 2;
    state.counters["LB"] = 2;
    state.counters["num_coeffs"] = 27;
}
BENCHMARK(BM_HermiteE_TMP_ShellPair_dd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

/**
 * @brief Benchmark all Hermite E coefficients for ff shell pair (highest complexity)
 */
static void BM_HermiteE_TMP_ShellPair_ff(benchmark::State& state) {
    GeometryParams params(42);
    int64_t total_geometries = 0;

    for (auto _ : state) {
        // ff: E^{nA,nB}_t for nA,nB in {0,1,2,3}, t in {0..nA+nB}
        // Selected representative coefficients for timing
        auto e_000 = HermiteE<0, 0, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_300 = HermiteE<3, 0, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_303 = HermiteE<3, 0, 3>::compute(params.PA[0], params.PB[0], params.p);
        auto e_030 = HermiteE<0, 3, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_033 = HermiteE<0, 3, 3>::compute(params.PA[0], params.PB[0], params.p);
        auto e_330 = HermiteE<3, 3, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_333 = HermiteE<3, 3, 3>::compute(params.PA[0], params.PB[0], params.p);
        auto e_336 = HermiteE<3, 3, 6>::compute(params.PA[0], params.PB[0], params.p);

        // Additional coefficients for complete coverage
        auto e_110 = HermiteE<1, 1, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_220 = HermiteE<2, 2, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_130 = HermiteE<1, 3, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_230 = HermiteE<2, 3, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_310 = HermiteE<3, 1, 0>::compute(params.PA[0], params.PB[0], params.p);
        auto e_320 = HermiteE<3, 2, 0>::compute(params.PA[0], params.PB[0], params.p);

        benchmark::DoNotOptimize(e_000);
        benchmark::DoNotOptimize(e_300);
        benchmark::DoNotOptimize(e_330);
        benchmark::DoNotOptimize(e_336);
        total_geometries += 8;
    }

    state.counters["Geometries/s"] = benchmark::Counter(total_geometries, benchmark::Counter::kIsRate);
    state.counters["LA"] = 3;
    state.counters["LB"] = 3;
    state.counters["num_coeffs"] = 64;  // (3+1)^2 * (6+1) / 2 approx
}
BENCHMARK(BM_HermiteE_TMP_ShellPair_ff)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// Memory Bandwidth Benchmarks
// =============================================================================

/**
 * @brief Memory bandwidth test for TMP implementation
 *
 * Tests cache behavior by computing many batches sequentially.
 */
static void BM_HermiteE_TMP_Memory(benchmark::State& state) {
    const int num_batches = state.range(0);
    std::vector<GeometryParams> params(num_batches);
    std::vector<Vec8d> results(num_batches);

    // Initialize with different seeds
    for (int i = 0; i < num_batches; ++i) {
        params[i] = GeometryParams(42 + i);
    }

    for (auto _ : state) {
        for (int i = 0; i < num_batches; ++i) {
            results[i] = HermiteE<3, 3, 6>::compute(
                params[i].PA[0], params[i].PB[0], params[i].p);
        }
        benchmark::DoNotOptimize(results.data());
        benchmark::ClobberMemory();
    }

    // Calculate bytes processed (input parameters + output)
    int64_t bytes_per_batch = 4 * sizeof(Vec8d);  // PA, PB, p, result
    state.SetBytesProcessed(state.iterations() * num_batches * bytes_per_batch);
    state.counters["Batches"] = num_batches;
}
BENCHMARK(BM_HermiteE_TMP_Memory)->RangeMultiplier(2)->Range(1, 1024)->Unit(benchmark::kMicrosecond);
