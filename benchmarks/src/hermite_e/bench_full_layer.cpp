/**
 * @file bench_full_layer.cpp
 * @brief Compare computing ALL t values (full layer) for gg shell pairs
 *
 * This tests whether the layered CSE approach wins when we need ALL t values,
 * not just a single one.
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
// Full Layer: Compute ALL t values for (4,4)
// =============================================================================

static void BM_FullLayer_Original_gg(benchmark::State& state) {
    GeometryParams params(42);
    Vec8d result;
    for (auto _ : state) {
        result = HermiteE<4, 4, 0>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<4, 4, 1>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<4, 4, 2>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<4, 4, 3>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<4, 4, 4>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<4, 4, 5>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<4, 4, 6>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<4, 4, 7>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<4, 4, 8>::compute(params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0;
}
BENCHMARK(BM_FullLayer_Original_gg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_FullLayer_Layered_gg(benchmark::State& state) {
    GeometryParams params(42);
    Vec8d result;
    for (auto _ : state) {
        auto layer = layered::HermiteELayer<4, 4>::compute(params.PA[0], params.PB[0], params.p);
        result = layer[0] + layer[1] + layer[2] + layer[3] + layer[4]
               + layer[5] + layer[6] + layer[7] + layer[8];
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2;
}
BENCHMARK(BM_FullLayer_Layered_gg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_FullLayer_Symbolic_gg(benchmark::State& state) {
    GeometryParams params(42);
    Vec8d result;
    for (auto _ : state) {
        result = recursum::symbolic::hermite_e_symbolic_4_4_0(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_4_4_1(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_4_4_2(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_4_4_3(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_4_4_4(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_4_4_5(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_4_4_6(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_4_4_7(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_4_4_8(params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1;
}
BENCHMARK(BM_FullLayer_Symbolic_gg)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// Full Layer: ff shell pairs
// =============================================================================

static void BM_FullLayer_Original_ff(benchmark::State& state) {
    GeometryParams params(42);
    Vec8d result;
    for (auto _ : state) {
        result = HermiteE<3, 3, 0>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<3, 3, 1>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<3, 3, 2>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<3, 3, 3>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<3, 3, 4>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<3, 3, 5>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<3, 3, 6>::compute(params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0;
}
BENCHMARK(BM_FullLayer_Original_ff)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_FullLayer_Layered_ff(benchmark::State& state) {
    GeometryParams params(42);
    Vec8d result;
    for (auto _ : state) {
        auto layer = layered::HermiteELayer<3, 3>::compute(params.PA[0], params.PB[0], params.p);
        result = layer[0] + layer[1] + layer[2] + layer[3] + layer[4] + layer[5] + layer[6];
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2;
}
BENCHMARK(BM_FullLayer_Layered_ff)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_FullLayer_Symbolic_ff(benchmark::State& state) {
    GeometryParams params(42);
    Vec8d result;
    for (auto _ : state) {
        result = recursum::symbolic::hermite_e_symbolic_3_3_0(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_3_3_1(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_3_3_2(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_3_3_3(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_3_3_4(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_3_3_5(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_3_3_6(params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1;
}
BENCHMARK(BM_FullLayer_Symbolic_ff)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// Full Layer: dd shell pairs
// =============================================================================

static void BM_FullLayer_Original_dd(benchmark::State& state) {
    GeometryParams params(42);
    Vec8d result;
    for (auto _ : state) {
        result = HermiteE<2, 2, 0>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<2, 2, 1>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<2, 2, 2>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<2, 2, 3>::compute(params.PA[0], params.PB[0], params.p);
        result += HermiteE<2, 2, 4>::compute(params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 0;
}
BENCHMARK(BM_FullLayer_Original_dd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_FullLayer_Layered_dd(benchmark::State& state) {
    GeometryParams params(42);
    Vec8d result;
    for (auto _ : state) {
        auto layer = layered::HermiteELayer<2, 2>::compute(params.PA[0], params.PB[0], params.p);
        result = layer[0] + layer[1] + layer[2] + layer[3] + layer[4];
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 2;
}
BENCHMARK(BM_FullLayer_Layered_dd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

static void BM_FullLayer_Symbolic_dd(benchmark::State& state) {
    GeometryParams params(42);
    Vec8d result;
    for (auto _ : state) {
        result = recursum::symbolic::hermite_e_symbolic_2_2_0(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_2_2_1(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_2_2_2(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_2_2_3(params.PA[0], params.PB[0], params.one_over_2p);
        result += recursum::symbolic::hermite_e_symbolic_2_2_4(params.PA[0], params.PB[0], params.one_over_2p);
        benchmark::DoNotOptimize(result);
    }
    state.counters["impl"] = 1;
}
BENCHMARK(BM_FullLayer_Symbolic_dd)->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// Main
// =============================================================================

int main(int argc, char** argv) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "FULL LAYER BENCHMARK: Computing ALL t values at once\n";
    std::cout << "Original TMP vs Layered CSE vs Symbolic\n";
    std::cout << std::string(80, '=') << "\n\n";
    std::cout << "impl=0: Original TMP (9 separate recursive calls)\n";
    std::cout << "impl=2: Layered TMP (compute layer once, extract all)\n";
    std::cout << "impl=1: Symbolic (9 separate polynomial evaluations)\n";
    std::cout << std::string(80, '=') << "\n\n";

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
