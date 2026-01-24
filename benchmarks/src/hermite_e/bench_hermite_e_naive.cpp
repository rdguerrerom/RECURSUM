/**
 * @file bench_hermite_e_naive.cpp
 * @brief Naive runtime-loop baseline for Hermite E coefficient benchmarks
 *
 * Implements a straightforward runtime loop-based evaluation of Hermite
 * expansion coefficients. This serves as a baseline to show how much
 * TMP and symbolic approaches improve over naive implementations.
 */

#include <benchmark/benchmark.h>
#include <vector>
#include <cmath>
#include "benchmark_common.hpp"

using namespace recursum::benchmark;

namespace {

/**
 * @brief Naive runtime evaluation of Hermite E coefficient
 *
 * Uses runtime recursion with memoization. This is the "textbook"
 * approach without compile-time optimization.
 *
 * Recurrence (Helgaker-Taylor 1992, Eq. 7):
 *   E^{i+1,j}_t = (1/2p) * E^{i,j}_{t-1} + PA * E^{i,j}_t + (t+1) * E^{i,j}_{t+1}
 *
 * For t=0:
 *   E^{i+1,j}_0 = PA * E^{i,j}_0 + E^{i,j}_1
 */
class NaiveHermiteE {
public:
    NaiveHermiteE(int max_nA, int max_nB)
        : max_nA_(max_nA), max_nB_(max_nB)
        , max_t_(max_nA + max_nB)
    {
        // Pre-allocate storage for memoization
        // Index: [nA][nB][t]
        cache_.resize((max_nA + 1) * (max_nB + 1) * (max_t_ + 1));
    }

    /**
     * @brief Compute E^{nA,nB}_t using runtime recursion
     */
    Vec8d compute(int nA, int nB, int t, Vec8d PA, Vec8d PB, Vec8d p) {
        // Clear cache for new computation
        std::fill(computed_.begin(), computed_.end(), false);
        computed_.resize((max_nA_ + 1) * (max_nB_ + 1) * (max_t_ + 1), false);

        PA_ = PA;
        PB_ = PB;
        one_over_2p_ = 1.0 / (2.0 * p);

        return computeRecursive(nA, nB, t);
    }

private:
    int max_nA_, max_nB_, max_t_;
    std::vector<Vec8d> cache_;
    std::vector<bool> computed_;
    Vec8d PA_, PB_, one_over_2p_;

    int index(int nA, int nB, int t) const {
        return nA * (max_nB_ + 1) * (max_t_ + 1) + nB * (max_t_ + 1) + t;
    }

    Vec8d computeRecursive(int nA, int nB, int t) {
        // Validity check
        if (nA < 0 || nB < 0 || t < 0 || t > nA + nB) {
            return Vec8d(0.0);
        }

        // Base case
        if (nA == 0 && nB == 0 && t == 0) {
            return Vec8d(1.0);
        }

        // Check memoization cache
        int idx = index(nA, nB, t);
        if (computed_[idx]) {
            return cache_[idx];
        }

        Vec8d result;

        if (nA > 0) {
            // A-side recurrence
            if (t == 0) {
                // E^{i+1,j}_0 = PA * E^{i,j}_0 + E^{i,j}_1
                result = PA_ * computeRecursive(nA - 1, nB, 0)
                       + computeRecursive(nA - 1, nB, 1);
            } else {
                // E^{i+1,j}_t = (1/2p) * E^{i,j}_{t-1} + PA * E^{i,j}_t + (t+1) * E^{i,j}_{t+1}
                result = one_over_2p_ * computeRecursive(nA - 1, nB, t - 1)
                       + PA_ * computeRecursive(nA - 1, nB, t)
                       + Vec8d(t + 1) * computeRecursive(nA - 1, nB, t + 1);
            }
        } else if (nB > 0) {
            // B-side recurrence
            if (t == 0) {
                result = PB_ * computeRecursive(0, nB - 1, 0)
                       + computeRecursive(0, nB - 1, 1);
            } else {
                result = one_over_2p_ * computeRecursive(0, nB - 1, t - 1)
                       + PB_ * computeRecursive(0, nB - 1, t)
                       + Vec8d(t + 1) * computeRecursive(0, nB - 1, t + 1);
            }
        } else {
            result = Vec8d(0.0);
        }

        // Memoize
        cache_[idx] = result;
        computed_[idx] = true;

        return result;
    }
};

/**
 * @brief Alternative: Bottom-up DP approach (no recursion overhead)
 */
class BottomUpHermiteE {
public:
    BottomUpHermiteE(int max_nA, int max_nB)
        : max_nA_(max_nA), max_nB_(max_nB)
        , stride_nB_(max_nB + 1)
        , stride_t_(max_nA + max_nB + 1)
    {
        table_.resize((max_nA + 1) * stride_nB_ * stride_t_);
    }

    /**
     * @brief Fill table bottom-up and return E^{nA,nB}_t
     */
    Vec8d compute(int nA, int nB, int t, Vec8d PA, Vec8d PB, Vec8d p) {
        Vec8d one_over_2p = 1.0 / (2.0 * p);

        // Initialize table to zero
        std::fill(table_.begin(), table_.end(), Vec8d(0.0));

        // Base case: E^{0,0}_0 = 1
        table_[index(0, 0, 0)] = Vec8d(1.0);

        // Fill A-side (nB = 0): E^{i,0}_t
        for (int i = 1; i <= max_nA_; ++i) {
            for (int ti = 0; ti <= i; ++ti) {
                if (ti == 0) {
                    table_[index(i, 0, 0)] = PA * table_[index(i-1, 0, 0)]
                                            + table_[index(i-1, 0, 1)];
                } else {
                    table_[index(i, 0, ti)] =
                        one_over_2p * table_[index(i-1, 0, ti-1)]
                        + PA * table_[index(i-1, 0, ti)]
                        + Vec8d(ti + 1) * table_[index(i-1, 0, ti+1)];
                }
            }
        }

        // Fill B-side (nA = 0): E^{0,j}_t
        for (int j = 1; j <= max_nB_; ++j) {
            for (int tj = 0; tj <= j; ++tj) {
                if (tj == 0) {
                    table_[index(0, j, 0)] = PB * table_[index(0, j-1, 0)]
                                            + table_[index(0, j-1, 1)];
                } else {
                    table_[index(0, j, tj)] =
                        one_over_2p * table_[index(0, j-1, tj-1)]
                        + PB * table_[index(0, j-1, tj)]
                        + Vec8d(tj + 1) * table_[index(0, j-1, tj+1)];
                }
            }
        }

        // Fill mixed (nA > 0, nB > 0): E^{i,j}_t
        for (int i = 1; i <= max_nA_; ++i) {
            for (int j = 1; j <= max_nB_; ++j) {
                for (int tij = 0; tij <= i + j; ++tij) {
                    if (tij == 0) {
                        table_[index(i, j, 0)] = PA * table_[index(i-1, j, 0)]
                                                + table_[index(i-1, j, 1)];
                    } else {
                        table_[index(i, j, tij)] =
                            one_over_2p * table_[index(i-1, j, tij-1)]
                            + PA * table_[index(i-1, j, tij)]
                            + Vec8d(tij + 1) * table_[index(i-1, j, tij+1)];
                    }
                }
            }
        }

        return table_[index(nA, nB, t)];
    }

private:
    int max_nA_, max_nB_;
    int stride_nB_, stride_t_;
    std::vector<Vec8d> table_;

    int index(int nA, int nB, int t) const {
        return nA * stride_nB_ * stride_t_ + nB * stride_t_ + t;
    }
};

} // anonymous namespace

// =============================================================================
// Naive Recursive Implementation Benchmarks
// =============================================================================

#define BENCH_HERMITE_E_NAIVE(NA, NB, T) \
static void BM_HermiteE_Naive_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    NaiveHermiteE solver(NA, NB); \
    Vec8d result; \
    int64_t total_geometries = 0; \
    for (auto _ : state) { \
        result = solver.compute(NA, NB, T, params.PA[0], params.PB[0], params.p); \
        benchmark::DoNotOptimize(result); \
        total_geometries += 8; \
    } \
    state.counters["Geometries/s"] = benchmark::Counter( \
        total_geometries, benchmark::Counter::kIsRate); \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
    state.counters["impl"] = 2; \
} \
BENCHMARK(BM_HermiteE_Naive_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0)

// Selected test cases (same as TMP benchmarks)
BENCH_HERMITE_E_NAIVE(0, 0, 0);
BENCH_HERMITE_E_NAIVE(1, 1, 2);
BENCH_HERMITE_E_NAIVE(2, 2, 4);
BENCH_HERMITE_E_NAIVE(3, 3, 6);

// =============================================================================
// Bottom-Up DP Implementation Benchmarks
// =============================================================================

#define BENCH_HERMITE_E_BOTTOMUP(NA, NB, T) \
static void BM_HermiteE_BottomUp_##NA##_##NB##_##T(benchmark::State& state) { \
    GeometryParams params(42); \
    BottomUpHermiteE solver(NA, NB); \
    Vec8d result; \
    int64_t total_geometries = 0; \
    for (auto _ : state) { \
        result = solver.compute(NA, NB, T, params.PA[0], params.PB[0], params.p); \
        benchmark::DoNotOptimize(result); \
        total_geometries += 8; \
    } \
    state.counters["Geometries/s"] = benchmark::Counter( \
        total_geometries, benchmark::Counter::kIsRate); \
    state.counters["nA"] = NA; \
    state.counters["nB"] = NB; \
    state.counters["t"] = T; \
    state.counters["impl"] = 3; \
} \
BENCHMARK(BM_HermiteE_BottomUp_##NA##_##NB##_##T)->Unit(benchmark::kNanosecond)->MinTime(1.0)

// Selected test cases
BENCH_HERMITE_E_BOTTOMUP(0, 0, 0);
BENCH_HERMITE_E_BOTTOMUP(1, 1, 2);
BENCH_HERMITE_E_BOTTOMUP(2, 2, 4);
BENCH_HERMITE_E_BOTTOMUP(3, 3, 6);

// =============================================================================
// Shell Pair Benchmarks (Naive)
// =============================================================================

static void BM_HermiteE_Naive_ShellPair_ff(benchmark::State& state) {
    GeometryParams params(42);
    NaiveHermiteE solver(3, 3);
    int64_t total_geometries = 0;

    for (auto _ : state) {
        // Compute all coefficients for ff shell pair
        Vec8d result(0.0);
        for (int nA = 0; nA <= 3; ++nA) {
            for (int nB = 0; nB <= 3; ++nB) {
                for (int t = 0; t <= nA + nB; ++t) {
                    result += solver.compute(nA, nB, t, params.PA[0], params.PB[0], params.p);
                }
            }
        }
        benchmark::DoNotOptimize(result);
        total_geometries += 8;
    }

    state.counters["Geometries/s"] = benchmark::Counter(total_geometries, benchmark::Counter::kIsRate);
    state.counters["LA"] = 3;
    state.counters["LB"] = 3;
}
BENCHMARK(BM_HermiteE_Naive_ShellPair_ff)->Unit(benchmark::kMicrosecond)->MinTime(1.0);

static void BM_HermiteE_BottomUp_ShellPair_ff(benchmark::State& state) {
    GeometryParams params(42);
    BottomUpHermiteE solver(3, 3);
    int64_t total_geometries = 0;

    for (auto _ : state) {
        // The bottom-up approach computes all coefficients in one pass
        // Just retrieve the highest one (others are computed along the way)
        auto result = solver.compute(3, 3, 6, params.PA[0], params.PB[0], params.p);
        benchmark::DoNotOptimize(result);
        total_geometries += 8;
    }

    state.counters["Geometries/s"] = benchmark::Counter(total_geometries, benchmark::Counter::kIsRate);
    state.counters["LA"] = 3;
    state.counters["LB"] = 3;
}
BENCHMARK(BM_HermiteE_BottomUp_ShellPair_ff)->Unit(benchmark::kMicrosecond)->MinTime(1.0);
