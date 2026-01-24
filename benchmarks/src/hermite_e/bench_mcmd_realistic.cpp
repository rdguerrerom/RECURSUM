/**
 * @file bench_mcmd_realistic.cpp
 * @brief REALISTIC McMurchie-Davidson Benchmark: Full Layer Computation
 *
 * In actual McMurchie-Davidson implementations, you ALWAYS need ALL t values
 * for a given (nA, nB) pair, not individual coefficients. This benchmark
 * tests the realistic use case.
 *
 * Compares:
 * - Original TMP: 9 separate recursive calls for gg
 * - Layered TMP: Compute full layer once with CSE
 * - Symbolic: 9 separate polynomial evaluations
 *
 * Coverage: All shell pairs from ss to gg
 */

#include <benchmark/benchmark.h>
#include <iostream>
#include <recursum/mcmd/hermite_e.hpp>
#include <recursum/mcmd/hermite_e_layered.hpp>
#include "hermite_e_symbolic.hpp"
#include "benchmark_common.hpp"

using namespace recursum::benchmark;
using namespace recursum::mcmd;

// =============================================================================
// Macro for Full Layer Benchmarks
// =============================================================================

// Helper to expand all t values for symbolic
#define SYM_SUM_0(NA, NB) \
    recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_0(params.PA[0], params.PB[0], params.one_over_2p)

#define SYM_SUM_1(NA, NB) \
    SYM_SUM_0(NA, NB) + recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_1(params.PA[0], params.PB[0], params.one_over_2p)

#define SYM_SUM_2(NA, NB) \
    SYM_SUM_1(NA, NB) + recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_2(params.PA[0], params.PB[0], params.one_over_2p)

#define SYM_SUM_3(NA, NB) \
    SYM_SUM_2(NA, NB) + recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_3(params.PA[0], params.PB[0], params.one_over_2p)

#define SYM_SUM_4(NA, NB) \
    SYM_SUM_3(NA, NB) + recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_4(params.PA[0], params.PB[0], params.one_over_2p)

#define SYM_SUM_5(NA, NB) \
    SYM_SUM_4(NA, NB) + recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_5(params.PA[0], params.PB[0], params.one_over_2p)

#define SYM_SUM_6(NA, NB) \
    SYM_SUM_5(NA, NB) + recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_6(params.PA[0], params.PB[0], params.one_over_2p)

#define SYM_SUM_7(NA, NB) \
    SYM_SUM_6(NA, NB) + recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_7(params.PA[0], params.PB[0], params.one_over_2p)

#define SYM_SUM_8(NA, NB) \
    SYM_SUM_7(NA, NB) + recursum::symbolic::hermite_e_symbolic_##NA##_##NB##_8(params.PA[0], params.PB[0], params.one_over_2p)

// Helper to expand all t values for original TMP
#define TMP_SUM_0(NA, NB) \
    HermiteE<NA, NB, 0>::compute(params.PA[0], params.PB[0], params.p)

#define TMP_SUM_1(NA, NB) \
    TMP_SUM_0(NA, NB) + HermiteE<NA, NB, 1>::compute(params.PA[0], params.PB[0], params.p)

#define TMP_SUM_2(NA, NB) \
    TMP_SUM_1(NA, NB) + HermiteE<NA, NB, 2>::compute(params.PA[0], params.PB[0], params.p)

#define TMP_SUM_3(NA, NB) \
    TMP_SUM_2(NA, NB) + HermiteE<NA, NB, 3>::compute(params.PA[0], params.PB[0], params.p)

#define TMP_SUM_4(NA, NB) \
    TMP_SUM_3(NA, NB) + HermiteE<NA, NB, 4>::compute(params.PA[0], params.PB[0], params.p)

#define TMP_SUM_5(NA, NB) \
    TMP_SUM_4(NA, NB) + HermiteE<NA, NB, 5>::compute(params.PA[0], params.PB[0], params.p)

#define TMP_SUM_6(NA, NB) \
    TMP_SUM_5(NA, NB) + HermiteE<NA, NB, 6>::compute(params.PA[0], params.PB[0], params.p)

#define TMP_SUM_7(NA, NB) \
    TMP_SUM_6(NA, NB) + HermiteE<NA, NB, 7>::compute(params.PA[0], params.PB[0], params.p)

#define TMP_SUM_8(NA, NB) \
    TMP_SUM_7(NA, NB) + HermiteE<NA, NB, 8>::compute(params.PA[0], params.PB[0], params.p)

// Helper to extract sum from layered
#define LAYER_SUM(NA, NB, TMAX) \
    [&]() { \
        auto layer = layered::HermiteELayer<NA, NB>::compute(params.PA[0], params.PB[0], params.p); \
        Vec8d sum = layer[0]; \
        for (int t = 1; t <= TMAX; ++t) sum += layer[t]; \
        return sum; \
    }()

// =============================================================================
// ss Shell Pair (L=0, 1 coefficient)
// =============================================================================

static void BM_Layer_Orig_ss(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_0(0, 0);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0;
    state.counters["L"] = 0;
    state.counters["n_coeffs"] = 1;
}
BENCHMARK(BM_Layer_Orig_ss)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_ss(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(0, 0, 0);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2;
    state.counters["L"] = 0;
    state.counters["n_coeffs"] = 1;
}
BENCHMARK(BM_Layer_CSE_ss)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_ss(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_0(0, 0);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1;
    state.counters["L"] = 0;
    state.counters["n_coeffs"] = 1;
}
BENCHMARK(BM_Layer_Sym_ss)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// sp Shell Pair (L=1, 2 coefficients)
// =============================================================================

static void BM_Layer_Orig_sp(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_1(0, 1);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0;
    state.counters["L"] = 1;
    state.counters["n_coeffs"] = 2;
}
BENCHMARK(BM_Layer_Orig_sp)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_sp(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(0, 1, 1);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2;
    state.counters["L"] = 1;
    state.counters["n_coeffs"] = 2;
}
BENCHMARK(BM_Layer_CSE_sp)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_sp(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_1(0, 1);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1;
    state.counters["L"] = 1;
    state.counters["n_coeffs"] = 2;
}
BENCHMARK(BM_Layer_Sym_sp)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// pp Shell Pair (L=2, 3 coefficients)
// =============================================================================

static void BM_Layer_Orig_pp(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_2(1, 1);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0;
    state.counters["L"] = 2;
    state.counters["n_coeffs"] = 3;
}
BENCHMARK(BM_Layer_Orig_pp)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_pp(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(1, 1, 2);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2;
    state.counters["L"] = 2;
    state.counters["n_coeffs"] = 3;
}
BENCHMARK(BM_Layer_CSE_pp)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_pp(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_2(1, 1);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1;
    state.counters["L"] = 2;
    state.counters["n_coeffs"] = 3;
}
BENCHMARK(BM_Layer_Sym_pp)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// dd Shell Pair (L=4, 5 coefficients)
// =============================================================================

static void BM_Layer_Orig_dd(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_4(2, 2);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0;
    state.counters["L"] = 4;
    state.counters["n_coeffs"] = 5;
}
BENCHMARK(BM_Layer_Orig_dd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_dd(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(2, 2, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2;
    state.counters["L"] = 4;
    state.counters["n_coeffs"] = 5;
}
BENCHMARK(BM_Layer_CSE_dd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_dd(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_4(2, 2);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1;
    state.counters["L"] = 4;
    state.counters["n_coeffs"] = 5;
}
BENCHMARK(BM_Layer_Sym_dd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// ff Shell Pair (L=6, 7 coefficients)
// =============================================================================

static void BM_Layer_Orig_ff(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_6(3, 3);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0;
    state.counters["L"] = 6;
    state.counters["n_coeffs"] = 7;
}
BENCHMARK(BM_Layer_Orig_ff)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_ff(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(3, 3, 6);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2;
    state.counters["L"] = 6;
    state.counters["n_coeffs"] = 7;
}
BENCHMARK(BM_Layer_CSE_ff)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_ff(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_6(3, 3);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1;
    state.counters["L"] = 6;
    state.counters["n_coeffs"] = 7;
}
BENCHMARK(BM_Layer_Sym_ff)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// gg Shell Pair (L=8, 9 coefficients)
// =============================================================================

static void BM_Layer_Orig_gg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_8(4, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0;
    state.counters["L"] = 8;
    state.counters["n_coeffs"] = 9;
}
BENCHMARK(BM_Layer_Orig_gg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_gg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(4, 4, 8);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2;
    state.counters["L"] = 8;
    state.counters["n_coeffs"] = 9;
}
BENCHMARK(BM_Layer_CSE_gg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_gg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_8(4, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1;
    state.counters["L"] = 8;
    state.counters["n_coeffs"] = 9;
}
BENCHMARK(BM_Layer_Sym_gg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// Additional Shell Pairs: sd, pd, pf, df, sg, pg, dg, fg
// =============================================================================

// sd (L=2, 3 coefficients)
static void BM_Layer_Orig_sd(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_2(0, 2);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0; state.counters["L"] = 2; state.counters["n_coeffs"] = 3;
}
BENCHMARK(BM_Layer_Orig_sd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_sd(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(0, 2, 2);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2; state.counters["L"] = 2; state.counters["n_coeffs"] = 3;
}
BENCHMARK(BM_Layer_CSE_sd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_sd(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_2(0, 2);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1; state.counters["L"] = 2; state.counters["n_coeffs"] = 3;
}
BENCHMARK(BM_Layer_Sym_sd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// pd (L=3, 4 coefficients)
static void BM_Layer_Orig_pd(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_3(1, 2);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0; state.counters["L"] = 3; state.counters["n_coeffs"] = 4;
}
BENCHMARK(BM_Layer_Orig_pd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_pd(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(1, 2, 3);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2; state.counters["L"] = 3; state.counters["n_coeffs"] = 4;
}
BENCHMARK(BM_Layer_CSE_pd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_pd(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_3(1, 2);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1; state.counters["L"] = 3; state.counters["n_coeffs"] = 4;
}
BENCHMARK(BM_Layer_Sym_pd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// sf (L=3, 4 coefficients)
static void BM_Layer_Orig_sf(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_3(0, 3);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0; state.counters["L"] = 3; state.counters["n_coeffs"] = 4;
}
BENCHMARK(BM_Layer_Orig_sf)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_sf(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(0, 3, 3);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2; state.counters["L"] = 3; state.counters["n_coeffs"] = 4;
}
BENCHMARK(BM_Layer_CSE_sf)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_sf(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_3(0, 3);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1; state.counters["L"] = 3; state.counters["n_coeffs"] = 4;
}
BENCHMARK(BM_Layer_Sym_sf)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// pf (L=4, 5 coefficients)
static void BM_Layer_Orig_pf(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_4(1, 3);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0; state.counters["L"] = 4; state.counters["n_coeffs"] = 5;
}
BENCHMARK(BM_Layer_Orig_pf)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_pf(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(1, 3, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2; state.counters["L"] = 4; state.counters["n_coeffs"] = 5;
}
BENCHMARK(BM_Layer_CSE_pf)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_pf(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_4(1, 3);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1; state.counters["L"] = 4; state.counters["n_coeffs"] = 5;
}
BENCHMARK(BM_Layer_Sym_pf)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// df (L=5, 6 coefficients)
static void BM_Layer_Orig_df(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_5(2, 3);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0; state.counters["L"] = 5; state.counters["n_coeffs"] = 6;
}
BENCHMARK(BM_Layer_Orig_df)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_df(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(2, 3, 5);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2; state.counters["L"] = 5; state.counters["n_coeffs"] = 6;
}
BENCHMARK(BM_Layer_CSE_df)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_df(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_5(2, 3);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1; state.counters["L"] = 5; state.counters["n_coeffs"] = 6;
}
BENCHMARK(BM_Layer_Sym_df)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// sg (L=4, 5 coefficients)
static void BM_Layer_Orig_sg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_4(0, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0; state.counters["L"] = 4; state.counters["n_coeffs"] = 5;
}
BENCHMARK(BM_Layer_Orig_sg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_sg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(0, 4, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2; state.counters["L"] = 4; state.counters["n_coeffs"] = 5;
}
BENCHMARK(BM_Layer_CSE_sg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_sg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_4(0, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1; state.counters["L"] = 4; state.counters["n_coeffs"] = 5;
}
BENCHMARK(BM_Layer_Sym_sg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// pg (L=5, 6 coefficients)
static void BM_Layer_Orig_pg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_5(1, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0; state.counters["L"] = 5; state.counters["n_coeffs"] = 6;
}
BENCHMARK(BM_Layer_Orig_pg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_pg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(1, 4, 5);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2; state.counters["L"] = 5; state.counters["n_coeffs"] = 6;
}
BENCHMARK(BM_Layer_CSE_pg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_pg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_5(1, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1; state.counters["L"] = 5; state.counters["n_coeffs"] = 6;
}
BENCHMARK(BM_Layer_Sym_pg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// dg (L=6, 7 coefficients)
static void BM_Layer_Orig_dg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_6(2, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0; state.counters["L"] = 6; state.counters["n_coeffs"] = 7;
}
BENCHMARK(BM_Layer_Orig_dg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_dg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(2, 4, 6);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2; state.counters["L"] = 6; state.counters["n_coeffs"] = 7;
}
BENCHMARK(BM_Layer_CSE_dg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_dg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_6(2, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1; state.counters["L"] = 6; state.counters["n_coeffs"] = 7;
}
BENCHMARK(BM_Layer_Sym_dg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// fg (L=7, 8 coefficients)
static void BM_Layer_Orig_fg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = TMP_SUM_7(3, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0; state.counters["L"] = 7; state.counters["n_coeffs"] = 8;
}
BENCHMARK(BM_Layer_Orig_fg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_CSE_fg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = LAYER_SUM(3, 4, 7);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2; state.counters["L"] = 7; state.counters["n_coeffs"] = 8;
}
BENCHMARK(BM_Layer_CSE_fg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_Layer_Sym_fg(benchmark::State& state) {
    GeometryParams params(42);
    for (auto _ : state) {
        Vec8d result = SYM_SUM_7(3, 4);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1; state.counters["L"] = 7; state.counters["n_coeffs"] = 8;
}
BENCHMARK(BM_Layer_Sym_fg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// Main
// =============================================================================

int main(int argc, char** argv) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "REALISTIC McMURCHIE-DAVIDSON BENCHMARK\n";
    std::cout << "Full Layer Computation (ALL t values for each shell pair)\n";
    std::cout << std::string(80, '=') << "\n\n";
    std::cout << "impl=0: Original TMP (separate recursive calls)\n";
    std::cout << "impl=2: Layered TMP with CSE (compute layer once)\n";
    std::cout << "impl=1: Symbolic (separate polynomial evaluations)\n";
    std::cout << std::string(80, '=') << "\n\n";
    std::cout << "Shell Pairs: ss, sp, pp, sd, pd, sf, pf, dd, df, sg, pg, ff, dg, fg, gg\n";
    std::cout << std::string(80, '=') << "\n\n";

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
