/**
 * @file bench_coulomb_hermite.cpp
 * @brief Benchmark: Coulomb Auxiliary Integrals R_{tuv}^{(N)} in Hermite Basis
 *
 * In the McMurchie-Davidson algorithm, ERIs require Coulomb auxiliary integrals
 * R_{tuv}^{(0)} for all (t,u,v) with t+u+v <= L_total.
 *
 * The number of R integrals needed is tetrahedral: (L+1)(L+2)(L+3)/6
 *
 * This benchmark measures the time to compute and STORE all R_{tuv} values
 * given pre-computed Boys function values.
 *
 * Recurrence relation (McMurchie-Davidson 1978):
 *   R_{tuv}^{(N)} = X_PC * R_{t-1,u,v}^{(N+1)} + (t-1) * R_{t-2,u,v}^{(N+1)}
 *   (with analogous relations for u and v)
 *
 * Base case:
 *   R_{000}^{(N)} = Boys_N(T)
 *
 * References:
 *   McMurchie & Davidson, J. Comput. Phys. 26 (1978) 218-231
 *   See benchmarks/McMD_METHODOLOGY.md Section 4 for detailed algorithm
 */

#include <benchmark/benchmark.h>
#include <array>
#include <cmath>
#include <iostream>
#include <recursum/mcmd/coulomb_aux.hpp>
#include <recursum/mcmd/coulomb_r_layered.hpp>
#include <recursum/mcmd/coulomb_r_symbolic.hpp>
#include "benchmark_common.hpp"

using namespace recursum::benchmark;
using namespace recursum::mcmd;

// Maximum L_total supported
constexpr int MAX_L = 8;
constexpr int MAX_R_SIZE = layered::tetrahedral(MAX_L);

// =============================================================================
// Boys Function Computation (simplified for benchmarking)
// =============================================================================

// Simple Boys function approximation for benchmark purposes
// In production, use quadbox or similar library
inline void compute_boys(double T, int n_max, Vec8d* Boys) {
    if (T < 1e-10) {
        for (int n = 0; n <= n_max; ++n) {
            Boys[n] = Vec8d(1.0 / (2 * n + 1));
        }
        return;
    }

    // Downward recursion from asymptotic
    double sqrtT = std::sqrt(T);
    double expT = std::exp(-T);

    // Start from large n with asymptotic approximation
    Boys[n_max] = Vec8d(0.5 * std::sqrt(M_PI / T) * expT);

    // Downward recursion: F_{n-1}(T) = (2T * F_n(T) + exp(-T)) / (2n - 1)
    for (int n = n_max; n > 0; --n) {
        double Fn = Boys[n][0];
        Boys[n - 1] = Vec8d((2.0 * T * Fn + expT) / (2.0 * n - 1.0));
    }
}

// =============================================================================
// Original TMP: Compute each R_{tuv} separately
// =============================================================================

template<int L_total>
static void BM_CoulombR_TMP(benchmark::State& state) {
    CoulombGeometryParams params(42);

    // Pre-compute Boys functions
    Vec8d Boys[MAX_L + 1];
    compute_boys(params.T[0], L_total, Boys);  // Use first element of T

    // Storage for R integrals
    std::array<Vec8d, MAX_R_SIZE> R_array;

    for (auto _ : state) {
        int idx = 0;
        // Compute all R_{tuv}^{(0)} for t+u+v <= L_total
        for (int t = 0; t <= L_total; ++t) {
            for (int u = 0; u <= L_total - t; ++u) {
                for (int v = 0; v <= L_total - t - u; ++v) {
                    // Use runtime dispatch for TMP (simulates realistic usage)
                    R_array[idx++] = CoulombR<0, 0, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
                }
            }
        }
        benchmark::DoNotOptimize(R_array.data());
        benchmark::ClobberMemory();
    }

    state.counters["impl"] = 0;  // Original TMP
    state.counters["L_total"] = L_total;
    state.counters["n_integrals"] = layered::tetrahedral(L_total);
}

// Specialized TMP benchmarks that actually use the template parameters
template<int L>
static void BM_CoulombR_TMP_Specialized(benchmark::State& state) {
    CoulombGeometryParams params(42);
    Vec8d Boys[L + 1];
    compute_boys(params.T[0], L, Boys);
    std::array<Vec8d, MAX_R_SIZE> R_array;

    for (auto _ : state) {
        int idx = 0;

        // L=0: just R_000
        if constexpr (L >= 0) {
            R_array[idx++] = CoulombR<0, 0, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
        }

        // L=1: R_000, R_100, R_010, R_001
        if constexpr (L >= 1) {
            R_array[idx++] = CoulombR<1, 0, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 1, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 0, 1, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
        }

        // L=2: add 6 more integrals
        if constexpr (L >= 2) {
            R_array[idx++] = CoulombR<2, 0, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<1, 1, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<1, 0, 1, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 2, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 1, 1, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 0, 2, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
        }

        // L=3: add 10 more integrals
        if constexpr (L >= 3) {
            R_array[idx++] = CoulombR<3, 0, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<2, 1, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<2, 0, 1, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<1, 2, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<1, 1, 1, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<1, 0, 2, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 3, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 2, 1, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 1, 2, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 0, 3, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
        }

        // L=4: add 15 more integrals
        if constexpr (L >= 4) {
            R_array[idx++] = CoulombR<4, 0, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<3, 1, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<3, 0, 1, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<2, 2, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<2, 1, 1, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<2, 0, 2, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<1, 3, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<1, 2, 1, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<1, 1, 2, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<1, 0, 3, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 4, 0, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 3, 1, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 2, 2, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 1, 3, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
            R_array[idx++] = CoulombR<0, 0, 4, 0>::compute(params.PC[0], params.PC[1], params.PC[2], Boys);
        }

        benchmark::DoNotOptimize(R_array.data());
        benchmark::ClobberMemory();
    }

    state.counters["impl"] = 0;
    state.counters["L_total"] = L;
    state.counters["n_integrals"] = layered::tetrahedral(L);
}

// =============================================================================
// Layered CSE: Compute all R_{tuv} at once with layer-by-layer CSE
// =============================================================================

template<int L_total>
static void BM_CoulombR_Layered(benchmark::State& state) {
    CoulombGeometryParams params(42);

    // Pre-compute Boys functions
    Vec8d Boys[MAX_L + 1];
    compute_boys(params.T[0], L_total, Boys);

    for (auto _ : state) {
        // Compute all R_{tuv}^{(0)} in one layer computation
        auto R_array = layered::CoulombRLayer<L_total>::compute(
            params.PC[0], params.PC[1], params.PC[2], Boys);

        benchmark::DoNotOptimize(R_array.data());
        benchmark::ClobberMemory();
    }

    state.counters["impl"] = 1;  // Layered CSE
    state.counters["L_total"] = L_total;
    state.counters["n_integrals"] = layered::tetrahedral(L_total);
}

// =============================================================================
// Symbolic/LayeredCodegen: Auto-generated optimized code with CSE
// =============================================================================

template<int L_total>
static void BM_CoulombR_Symbolic(benchmark::State& state) {
    CoulombGeometryParams params(42);

    // Pre-compute Boys functions
    Vec8d Boys[MAX_L + 1];
    compute_boys(params.T[0], L_total, Boys);

    for (auto _ : state) {
        // Compute all R_{tuv}^{(0)} using symbolic/codegen approach
        auto R_array = std::array<Vec8d, layered::tetrahedral(L_total)>();
        CoulombRSymbolic<L_total>::compute(R_array.data(),
            params.PC[0], params.PC[1], params.PC[2], Boys);

        benchmark::DoNotOptimize(R_array.data());
        benchmark::ClobberMemory();
    }

    state.counters["impl"] = 3;  // Symbolic/LayeredCodegen
    state.counters["L_total"] = L_total;
    state.counters["n_integrals"] = layered::tetrahedral(L_total);
}

// =============================================================================
// Register Benchmarks for Different L_total Values
// =============================================================================

// L=0: 1 integral (R_000)
BENCHMARK(BM_CoulombR_TMP_Specialized<0>)->Name("CoulombR/TMP/L0")->Unit(benchmark::kNanosecond)->MinTime(1.0);
BENCHMARK(BM_CoulombR_Layered<0>)->Name("CoulombR/Layered/L0")->Unit(benchmark::kNanosecond)->MinTime(1.0);
BENCHMARK(BM_CoulombR_Symbolic<0>)->Name("CoulombR/LayeredCodegen/L0")->Unit(benchmark::kNanosecond)->MinTime(1.0);

// L=1: 4 integrals
BENCHMARK(BM_CoulombR_TMP_Specialized<1>)->Name("CoulombR/TMP/L1")->Unit(benchmark::kNanosecond)->MinTime(1.0);
BENCHMARK(BM_CoulombR_Layered<1>)->Name("CoulombR/Layered/L1")->Unit(benchmark::kNanosecond)->MinTime(1.0);
BENCHMARK(BM_CoulombR_Symbolic<1>)->Name("CoulombR/LayeredCodegen/L1")->Unit(benchmark::kNanosecond)->MinTime(1.0);

// L=2: 10 integrals
BENCHMARK(BM_CoulombR_TMP_Specialized<2>)->Name("CoulombR/TMP/L2")->Unit(benchmark::kNanosecond)->MinTime(1.0);
BENCHMARK(BM_CoulombR_Layered<2>)->Name("CoulombR/Layered/L2")->Unit(benchmark::kNanosecond)->MinTime(1.0);
BENCHMARK(BM_CoulombR_Symbolic<2>)->Name("CoulombR/LayeredCodegen/L2")->Unit(benchmark::kNanosecond)->MinTime(1.0);

// L=3: 20 integrals
BENCHMARK(BM_CoulombR_TMP_Specialized<3>)->Name("CoulombR/TMP/L3")->Unit(benchmark::kNanosecond)->MinTime(1.0);
BENCHMARK(BM_CoulombR_Layered<3>)->Name("CoulombR/Layered/L3")->Unit(benchmark::kNanosecond)->MinTime(1.0);
BENCHMARK(BM_CoulombR_Symbolic<3>)->Name("CoulombR/LayeredCodegen/L3")->Unit(benchmark::kNanosecond)->MinTime(1.0);

// L=4: 35 integrals
BENCHMARK(BM_CoulombR_TMP_Specialized<4>)->Name("CoulombR/TMP/L4")->Unit(benchmark::kNanosecond)->MinTime(1.0);
BENCHMARK(BM_CoulombR_Layered<4>)->Name("CoulombR/Layered/L4")->Unit(benchmark::kNanosecond)->MinTime(1.0);
BENCHMARK(BM_CoulombR_Symbolic<4>)->Name("CoulombR/LayeredCodegen/L4")->Unit(benchmark::kNanosecond)->MinTime(1.0);

// L=5: 56 integrals
BENCHMARK(BM_CoulombR_Layered<5>)->Name("CoulombR/Layered/L5")->Unit(benchmark::kNanosecond)->MinTime(1.0);

// L=6: 84 integrals
BENCHMARK(BM_CoulombR_Layered<6>)->Name("CoulombR/Layered/L6")->Unit(benchmark::kNanosecond)->MinTime(1.0);

// L=7: 120 integrals
BENCHMARK(BM_CoulombR_Layered<7>)->Name("CoulombR/Layered/L7")->Unit(benchmark::kNanosecond)->MinTime(1.0);

// L=8: 165 integrals
BENCHMARK(BM_CoulombR_Layered<8>)->Name("CoulombR/Layered/L8")->Unit(benchmark::kNanosecond)->MinTime(1.0);

// =============================================================================
// Main
// =============================================================================

int main(int argc, char** argv) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "COULOMB AUXILIARY INTEGRALS BENCHMARK\n";
    std::cout << "McMurchie-Davidson R_{tuv}^{(0)} in Hermite Basis\n";
    std::cout << std::string(80, '=') << "\n\n";

    std::cout << "Timed operation: Compute and STORE all R_{tuv}^{(0)} for t+u+v <= L_total\n";
    std::cout << "Given pre-computed Boys function values F_N(T)\n\n";

    std::cout << "Number of integrals by L_total (tetrahedral number):\n";
    for (int L = 0; L <= 8; ++L) {
        std::cout << "  L=" << L << ": " << layered::tetrahedral(L) << " integrals\n";
    }
    std::cout << "\n";

    std::cout << "Implementations:\n";
    std::cout << "  TMP:            Original recursive template metaprogramming\n";
    std::cout << "  Layered:        Layer-by-layer with common subexpression elimination\n";
    std::cout << "  LayeredCodegen: Auto-generated optimized code with CSE (symbolic expansion)\n";
    std::cout << std::string(80, '=') << "\n\n";

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
