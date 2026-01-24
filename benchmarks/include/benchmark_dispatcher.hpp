/**
 * @file benchmark_dispatcher.hpp
 * @brief Comprehensive dispatchers for TMP and Symbolic Hermite E benchmarks
 *
 * AUTO-GENERATED - DO NOT EDIT MANUALLY
 *
 * Coverage: E^{0,0}_0 to E^{4,4}_{8} (125 coefficients)
 * Also includes gradient dispatchers for dE/dPA, dE/dPB
 */

#pragma once

#include <functional>
#include <array>
#include <stdexcept>

#ifndef RECURSUM_VEC_TYPE
#include <recursum/vectorclass.h>
#endif

#include <recursum/mcmd/hermite_e.hpp>
#include <recursum/mcmd/hermite_grad.hpp>
#include "hermite_e_symbolic.hpp"
#include "hermite_grad_symbolic.hpp"

namespace recursum {
namespace benchmark {

// Maximum angular momentum for benchmarks
constexpr int BENCH_MAX_L = 4;  // Up to g-type (L=4)
constexpr int BENCH_MAX_N = 2 * BENCH_MAX_L;  // Max t index = 8

/**
 * @brief Function signature for coefficient evaluation
 */
using CoeffFunction = std::function<Vec8d(const Vec8d&, const Vec8d&, const Vec8d&)>;

/**
 * @brief Function signature for gradient evaluation (5 params)
 */
using GradFunction = std::function<Vec8d(const Vec8d&, const Vec8d&, const Vec8d&, const Vec8d&, const Vec8d&)>;

// =============================================================================
// TMP DISPATCHER
// =============================================================================

/**
 * @class TMPDispatcher
 * @brief Dispatcher for TMP HermiteE template functions
 */
class TMPDispatcher {
public:
    static constexpr int MAX_NA = BENCH_MAX_L;
    static constexpr int MAX_NB = BENCH_MAX_L;
    static constexpr int MAX_N = BENCH_MAX_N;

    using NTable = std::array<CoeffFunction, MAX_N + 1>;
    using NBTable = std::array<NTable, MAX_NB + 1>;
    using DispatchTable = std::array<NBTable, MAX_NA + 1>;

    static const TMPDispatcher& instance() {
        static TMPDispatcher dispatcher;
        return dispatcher;
    }

    Vec8d compute(int nA, int nB, int t,
                  const Vec8d& PA, const Vec8d& PB, const Vec8d& p) const {
        if (nA < 0 || nA > MAX_NA || nB < 0 || nB > MAX_NB ||
            t < 0 || t > MAX_N || t > nA + nB) {
            return Vec8d(0.0);
        }
        const auto& func = table_[nA][nB][t];
        if (!func) return Vec8d(0.0);
        return func(PA, PB, p);
    }

    TMPDispatcher(const TMPDispatcher&) = delete;
    TMPDispatcher& operator=(const TMPDispatcher&) = delete;

private:
    TMPDispatcher() {
        for (auto& na_slice : table_) {
            for (auto& nb_slice : na_slice) {
                for (auto& func : nb_slice) {
                    func = nullptr;
                }
            }
        }
        registerAll();
    }

    template<int NA, int NB, int T>
    void registerCoeff() {
        table_[NA][NB][T] = [](const Vec8d& pa, const Vec8d& pb, const Vec8d& p) -> Vec8d {
            return mcmd::HermiteE<NA, NB, T>::compute(pa, pb, p);
        };
    }

    template<int NA, int NB, int T>
    struct RegisterT {
        static void apply(TMPDispatcher* d) {
            if constexpr (T <= NA + NB && T <= MAX_N) {
                d->registerCoeff<NA, NB, T>();
                if constexpr (T + 1 <= NA + NB && T + 1 <= MAX_N) {
                    RegisterT<NA, NB, T + 1>::apply(d);
                }
            }
        }
    };

    template<int NA, int NB>
    struct RegisterNB {
        static void apply(TMPDispatcher* d) {
            RegisterT<NA, NB, 0>::apply(d);
            if constexpr (NB + 1 <= MAX_NB) {
                RegisterNB<NA, NB + 1>::apply(d);
            }
        }
    };

    template<int NA>
    struct RegisterNA {
        static void apply(TMPDispatcher* d) {
            RegisterNB<NA, 0>::apply(d);
            if constexpr (NA + 1 <= MAX_NA) {
                RegisterNA<NA + 1>::apply(d);
            }
        }
    };

    void registerAll() {
        RegisterNA<0>::apply(this);
    }

    DispatchTable table_;
};

// =============================================================================
// SYMBOLIC DISPATCHER
// =============================================================================

/**
 * @class SymbolicDispatcher
 * @brief Dispatcher for Symbolic (SymPy-generated) functions
 */
class SymbolicDispatcher {
public:
    static constexpr int MAX_NA = BENCH_MAX_L;
    static constexpr int MAX_NB = BENCH_MAX_L;
    static constexpr int MAX_N = BENCH_MAX_N;

    using NTable = std::array<CoeffFunction, MAX_N + 1>;
    using NBTable = std::array<NTable, MAX_NB + 1>;
    using DispatchTable = std::array<NBTable, MAX_NA + 1>;

    static const SymbolicDispatcher& instance() {
        static SymbolicDispatcher dispatcher;
        return dispatcher;
    }

    Vec8d compute(int nA, int nB, int t,
                  const Vec8d& PA, const Vec8d& PB, const Vec8d& one_over_2p) const {
        if (nA < 0 || nA > MAX_NA || nB < 0 || nB > MAX_NB ||
            t < 0 || t > MAX_N || t > nA + nB) {
            return Vec8d(0.0);
        }
        const auto& func = table_[nA][nB][t];
        if (!func) return Vec8d(0.0);
        return func(PA, PB, one_over_2p);
    }

    SymbolicDispatcher(const SymbolicDispatcher&) = delete;
    SymbolicDispatcher& operator=(const SymbolicDispatcher&) = delete;

private:
    SymbolicDispatcher() {
        for (auto& na_slice : table_) {
            for (auto& nb_slice : na_slice) {
                for (auto& func : nb_slice) {
                    func = nullptr;
                }
            }
        }
        registerAll();
    }

    void registerAll() {
        table_[0][0][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_0_0(PA, PB, p);
        };
        table_[0][1][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_1_0(PA, PB, p);
        };
        table_[0][1][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_1_1(PA, PB, p);
        };
        table_[0][2][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_2_0(PA, PB, p);
        };
        table_[0][2][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_2_1(PA, PB, p);
        };
        table_[0][2][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_2_2(PA, PB, p);
        };
        table_[0][3][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_3_0(PA, PB, p);
        };
        table_[0][3][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_3_1(PA, PB, p);
        };
        table_[0][3][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_3_2(PA, PB, p);
        };
        table_[0][3][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_3_3(PA, PB, p);
        };
        table_[0][4][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_4_0(PA, PB, p);
        };
        table_[0][4][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_4_1(PA, PB, p);
        };
        table_[0][4][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_4_2(PA, PB, p);
        };
        table_[0][4][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_4_3(PA, PB, p);
        };
        table_[0][4][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_0_4_4(PA, PB, p);
        };
        table_[1][0][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_0_0(PA, PB, p);
        };
        table_[1][0][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_0_1(PA, PB, p);
        };
        table_[1][1][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_1_0(PA, PB, p);
        };
        table_[1][1][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_1_1(PA, PB, p);
        };
        table_[1][1][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_1_2(PA, PB, p);
        };
        table_[1][2][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_2_0(PA, PB, p);
        };
        table_[1][2][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_2_1(PA, PB, p);
        };
        table_[1][2][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_2_2(PA, PB, p);
        };
        table_[1][2][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_2_3(PA, PB, p);
        };
        table_[1][3][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_3_0(PA, PB, p);
        };
        table_[1][3][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_3_1(PA, PB, p);
        };
        table_[1][3][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_3_2(PA, PB, p);
        };
        table_[1][3][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_3_3(PA, PB, p);
        };
        table_[1][3][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_3_4(PA, PB, p);
        };
        table_[1][4][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_4_0(PA, PB, p);
        };
        table_[1][4][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_4_1(PA, PB, p);
        };
        table_[1][4][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_4_2(PA, PB, p);
        };
        table_[1][4][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_4_3(PA, PB, p);
        };
        table_[1][4][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_4_4(PA, PB, p);
        };
        table_[1][4][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_1_4_5(PA, PB, p);
        };
        table_[2][0][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_0_0(PA, PB, p);
        };
        table_[2][0][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_0_1(PA, PB, p);
        };
        table_[2][0][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_0_2(PA, PB, p);
        };
        table_[2][1][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_1_0(PA, PB, p);
        };
        table_[2][1][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_1_1(PA, PB, p);
        };
        table_[2][1][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_1_2(PA, PB, p);
        };
        table_[2][1][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_1_3(PA, PB, p);
        };
        table_[2][2][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_2_0(PA, PB, p);
        };
        table_[2][2][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_2_1(PA, PB, p);
        };
        table_[2][2][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_2_2(PA, PB, p);
        };
        table_[2][2][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_2_3(PA, PB, p);
        };
        table_[2][2][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_2_4(PA, PB, p);
        };
        table_[2][3][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_3_0(PA, PB, p);
        };
        table_[2][3][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_3_1(PA, PB, p);
        };
        table_[2][3][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_3_2(PA, PB, p);
        };
        table_[2][3][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_3_3(PA, PB, p);
        };
        table_[2][3][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_3_4(PA, PB, p);
        };
        table_[2][3][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_3_5(PA, PB, p);
        };
        table_[2][4][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_4_0(PA, PB, p);
        };
        table_[2][4][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_4_1(PA, PB, p);
        };
        table_[2][4][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_4_2(PA, PB, p);
        };
        table_[2][4][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_4_3(PA, PB, p);
        };
        table_[2][4][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_4_4(PA, PB, p);
        };
        table_[2][4][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_4_5(PA, PB, p);
        };
        table_[2][4][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_2_4_6(PA, PB, p);
        };
        table_[3][0][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_0_0(PA, PB, p);
        };
        table_[3][0][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_0_1(PA, PB, p);
        };
        table_[3][0][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_0_2(PA, PB, p);
        };
        table_[3][0][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_0_3(PA, PB, p);
        };
        table_[3][1][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_1_0(PA, PB, p);
        };
        table_[3][1][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_1_1(PA, PB, p);
        };
        table_[3][1][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_1_2(PA, PB, p);
        };
        table_[3][1][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_1_3(PA, PB, p);
        };
        table_[3][1][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_1_4(PA, PB, p);
        };
        table_[3][2][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_2_0(PA, PB, p);
        };
        table_[3][2][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_2_1(PA, PB, p);
        };
        table_[3][2][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_2_2(PA, PB, p);
        };
        table_[3][2][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_2_3(PA, PB, p);
        };
        table_[3][2][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_2_4(PA, PB, p);
        };
        table_[3][2][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_2_5(PA, PB, p);
        };
        table_[3][3][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_3_0(PA, PB, p);
        };
        table_[3][3][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_3_1(PA, PB, p);
        };
        table_[3][3][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_3_2(PA, PB, p);
        };
        table_[3][3][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_3_3(PA, PB, p);
        };
        table_[3][3][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_3_4(PA, PB, p);
        };
        table_[3][3][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_3_5(PA, PB, p);
        };
        table_[3][3][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_3_6(PA, PB, p);
        };
        table_[3][4][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_4_0(PA, PB, p);
        };
        table_[3][4][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_4_1(PA, PB, p);
        };
        table_[3][4][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_4_2(PA, PB, p);
        };
        table_[3][4][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_4_3(PA, PB, p);
        };
        table_[3][4][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_4_4(PA, PB, p);
        };
        table_[3][4][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_4_5(PA, PB, p);
        };
        table_[3][4][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_4_6(PA, PB, p);
        };
        table_[3][4][7] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_3_4_7(PA, PB, p);
        };
        table_[4][0][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_0_0(PA, PB, p);
        };
        table_[4][0][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_0_1(PA, PB, p);
        };
        table_[4][0][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_0_2(PA, PB, p);
        };
        table_[4][0][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_0_3(PA, PB, p);
        };
        table_[4][0][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_0_4(PA, PB, p);
        };
        table_[4][1][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_1_0(PA, PB, p);
        };
        table_[4][1][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_1_1(PA, PB, p);
        };
        table_[4][1][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_1_2(PA, PB, p);
        };
        table_[4][1][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_1_3(PA, PB, p);
        };
        table_[4][1][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_1_4(PA, PB, p);
        };
        table_[4][1][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_1_5(PA, PB, p);
        };
        table_[4][2][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_2_0(PA, PB, p);
        };
        table_[4][2][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_2_1(PA, PB, p);
        };
        table_[4][2][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_2_2(PA, PB, p);
        };
        table_[4][2][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_2_3(PA, PB, p);
        };
        table_[4][2][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_2_4(PA, PB, p);
        };
        table_[4][2][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_2_5(PA, PB, p);
        };
        table_[4][2][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_2_6(PA, PB, p);
        };
        table_[4][3][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_3_0(PA, PB, p);
        };
        table_[4][3][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_3_1(PA, PB, p);
        };
        table_[4][3][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_3_2(PA, PB, p);
        };
        table_[4][3][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_3_3(PA, PB, p);
        };
        table_[4][3][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_3_4(PA, PB, p);
        };
        table_[4][3][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_3_5(PA, PB, p);
        };
        table_[4][3][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_3_6(PA, PB, p);
        };
        table_[4][3][7] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_3_7(PA, PB, p);
        };
        table_[4][4][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_4_0(PA, PB, p);
        };
        table_[4][4][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_4_1(PA, PB, p);
        };
        table_[4][4][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_4_2(PA, PB, p);
        };
        table_[4][4][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_4_3(PA, PB, p);
        };
        table_[4][4][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_4_4(PA, PB, p);
        };
        table_[4][4][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_4_5(PA, PB, p);
        };
        table_[4][4][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_4_6(PA, PB, p);
        };
        table_[4][4][7] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_4_7(PA, PB, p);
        };
        table_[4][4][8] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_e_symbolic_4_4_8(PA, PB, p);
        };
    }

    DispatchTable table_;
};

// =============================================================================
// GRADIENT DISPATCHERS
// =============================================================================

/**
 * @class TMPGradPADispatcher
 * @brief Dispatcher for TMP dE/dPA gradients
 */
class TMPGradPADispatcher {
public:
    static constexpr int MAX_NA = BENCH_MAX_L;
    static constexpr int MAX_NB = BENCH_MAX_L;
    static constexpr int MAX_N = BENCH_MAX_N;

    using NTable = std::array<CoeffFunction, MAX_N + 1>;
    using NBTable = std::array<NTable, MAX_NB + 1>;
    using DispatchTable = std::array<NBTable, MAX_NA + 1>;

    static const TMPGradPADispatcher& instance() {
        static TMPGradPADispatcher dispatcher;
        return dispatcher;
    }

    Vec8d compute(int nA, int nB, int t,
                  const Vec8d& PA, const Vec8d& PB, const Vec8d& p) const {
        if (nA < 0 || nA > MAX_NA || nB < 0 || nB > MAX_NB ||
            t < 0 || t > MAX_N || t > nA + nB) {
            return Vec8d(0.0);
        }
        const auto& func = table_[nA][nB][t];
        if (!func) return Vec8d(0.0);
        return func(PA, PB, p);
    }

private:
    TMPGradPADispatcher() {
        for (auto& na_slice : table_) {
            for (auto& nb_slice : na_slice) {
                for (auto& func : nb_slice) {
                    func = nullptr;
                }
            }
        }
        registerAll();
    }

    template<int NA, int NB, int T>
    void registerGrad() {
        table_[NA][NB][T] = [](const Vec8d& pa, const Vec8d& pb, const Vec8d& p) -> Vec8d {
            return mcmd::HermiteDerivPA<NA, NB, T>::compute(pa, pb, p);
        };
    }

    template<int NA, int NB, int T>
    struct RegisterT {
        static void apply(TMPGradPADispatcher* d) {
            if constexpr (T <= NA + NB && T <= MAX_N) {
                d->registerGrad<NA, NB, T>();
                if constexpr (T + 1 <= NA + NB && T + 1 <= MAX_N) {
                    RegisterT<NA, NB, T + 1>::apply(d);
                }
            }
        }
    };

    template<int NA, int NB>
    struct RegisterNB {
        static void apply(TMPGradPADispatcher* d) {
            RegisterT<NA, NB, 0>::apply(d);
            if constexpr (NB + 1 <= MAX_NB) {
                RegisterNB<NA, NB + 1>::apply(d);
            }
        }
    };

    template<int NA>
    struct RegisterNA {
        static void apply(TMPGradPADispatcher* d) {
            RegisterNB<NA, 0>::apply(d);
            if constexpr (NA + 1 <= MAX_NA) {
                RegisterNA<NA + 1>::apply(d);
            }
        }
    };

    void registerAll() {
        RegisterNA<0>::apply(this);
    }

    DispatchTable table_;
};

/**
 * @class SymbolicGradPADispatcher
 * @brief Dispatcher for Symbolic dE/dPA gradients
 */
class SymbolicGradPADispatcher {
public:
    static constexpr int MAX_NA = BENCH_MAX_L;
    static constexpr int MAX_NB = BENCH_MAX_L;
    static constexpr int MAX_N = BENCH_MAX_N;

    using NTable = std::array<CoeffFunction, MAX_N + 1>;
    using NBTable = std::array<NTable, MAX_NB + 1>;
    using DispatchTable = std::array<NBTable, MAX_NA + 1>;

    static const SymbolicGradPADispatcher& instance() {
        static SymbolicGradPADispatcher dispatcher;
        return dispatcher;
    }

    Vec8d compute(int nA, int nB, int t,
                  const Vec8d& PA, const Vec8d& PB, const Vec8d& one_over_2p) const {
        if (nA < 0 || nA > MAX_NA || nB < 0 || nB > MAX_NB ||
            t < 0 || t > MAX_N || t > nA + nB) {
            return Vec8d(0.0);
        }
        const auto& func = table_[nA][nB][t];
        if (!func) return Vec8d(0.0);
        return func(PA, PB, one_over_2p);
    }

private:
    SymbolicGradPADispatcher() {
        for (auto& na_slice : table_) {
            for (auto& nb_slice : na_slice) {
                for (auto& func : nb_slice) {
                    func = nullptr;
                }
            }
        }
        registerAll();
    }

    void registerAll() {
        table_[0][0][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_0_0(PA, PB, p);
        };
        table_[0][1][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_1_0(PA, PB, p);
        };
        table_[0][1][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_1_1(PA, PB, p);
        };
        table_[0][2][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_2_0(PA, PB, p);
        };
        table_[0][2][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_2_1(PA, PB, p);
        };
        table_[0][2][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_2_2(PA, PB, p);
        };
        table_[0][3][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_3_0(PA, PB, p);
        };
        table_[0][3][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_3_1(PA, PB, p);
        };
        table_[0][3][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_3_2(PA, PB, p);
        };
        table_[0][3][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_3_3(PA, PB, p);
        };
        table_[0][4][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_4_0(PA, PB, p);
        };
        table_[0][4][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_4_1(PA, PB, p);
        };
        table_[0][4][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_4_2(PA, PB, p);
        };
        table_[0][4][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_4_3(PA, PB, p);
        };
        table_[0][4][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_0_4_4(PA, PB, p);
        };
        table_[1][0][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_0_0(PA, PB, p);
        };
        table_[1][0][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_0_1(PA, PB, p);
        };
        table_[1][1][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_1_0(PA, PB, p);
        };
        table_[1][1][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_1_1(PA, PB, p);
        };
        table_[1][1][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_1_2(PA, PB, p);
        };
        table_[1][2][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_2_0(PA, PB, p);
        };
        table_[1][2][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_2_1(PA, PB, p);
        };
        table_[1][2][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_2_2(PA, PB, p);
        };
        table_[1][2][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_2_3(PA, PB, p);
        };
        table_[1][3][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_3_0(PA, PB, p);
        };
        table_[1][3][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_3_1(PA, PB, p);
        };
        table_[1][3][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_3_2(PA, PB, p);
        };
        table_[1][3][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_3_3(PA, PB, p);
        };
        table_[1][3][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_3_4(PA, PB, p);
        };
        table_[1][4][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_4_0(PA, PB, p);
        };
        table_[1][4][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_4_1(PA, PB, p);
        };
        table_[1][4][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_4_2(PA, PB, p);
        };
        table_[1][4][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_4_3(PA, PB, p);
        };
        table_[1][4][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_4_4(PA, PB, p);
        };
        table_[1][4][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_1_4_5(PA, PB, p);
        };
        table_[2][0][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_0_0(PA, PB, p);
        };
        table_[2][0][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_0_1(PA, PB, p);
        };
        table_[2][0][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_0_2(PA, PB, p);
        };
        table_[2][1][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_1_0(PA, PB, p);
        };
        table_[2][1][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_1_1(PA, PB, p);
        };
        table_[2][1][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_1_2(PA, PB, p);
        };
        table_[2][1][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_1_3(PA, PB, p);
        };
        table_[2][2][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_2_0(PA, PB, p);
        };
        table_[2][2][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_2_1(PA, PB, p);
        };
        table_[2][2][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_2_2(PA, PB, p);
        };
        table_[2][2][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_2_3(PA, PB, p);
        };
        table_[2][2][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_2_4(PA, PB, p);
        };
        table_[2][3][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_3_0(PA, PB, p);
        };
        table_[2][3][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_3_1(PA, PB, p);
        };
        table_[2][3][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_3_2(PA, PB, p);
        };
        table_[2][3][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_3_3(PA, PB, p);
        };
        table_[2][3][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_3_4(PA, PB, p);
        };
        table_[2][3][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_3_5(PA, PB, p);
        };
        table_[2][4][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_4_0(PA, PB, p);
        };
        table_[2][4][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_4_1(PA, PB, p);
        };
        table_[2][4][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_4_2(PA, PB, p);
        };
        table_[2][4][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_4_3(PA, PB, p);
        };
        table_[2][4][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_4_4(PA, PB, p);
        };
        table_[2][4][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_4_5(PA, PB, p);
        };
        table_[2][4][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_2_4_6(PA, PB, p);
        };
        table_[3][0][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_0_0(PA, PB, p);
        };
        table_[3][0][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_0_1(PA, PB, p);
        };
        table_[3][0][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_0_2(PA, PB, p);
        };
        table_[3][0][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_0_3(PA, PB, p);
        };
        table_[3][1][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_1_0(PA, PB, p);
        };
        table_[3][1][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_1_1(PA, PB, p);
        };
        table_[3][1][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_1_2(PA, PB, p);
        };
        table_[3][1][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_1_3(PA, PB, p);
        };
        table_[3][1][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_1_4(PA, PB, p);
        };
        table_[3][2][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_2_0(PA, PB, p);
        };
        table_[3][2][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_2_1(PA, PB, p);
        };
        table_[3][2][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_2_2(PA, PB, p);
        };
        table_[3][2][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_2_3(PA, PB, p);
        };
        table_[3][2][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_2_4(PA, PB, p);
        };
        table_[3][2][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_2_5(PA, PB, p);
        };
        table_[3][3][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_3_0(PA, PB, p);
        };
        table_[3][3][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_3_1(PA, PB, p);
        };
        table_[3][3][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_3_2(PA, PB, p);
        };
        table_[3][3][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_3_3(PA, PB, p);
        };
        table_[3][3][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_3_4(PA, PB, p);
        };
        table_[3][3][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_3_5(PA, PB, p);
        };
        table_[3][3][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_3_6(PA, PB, p);
        };
        table_[3][4][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_4_0(PA, PB, p);
        };
        table_[3][4][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_4_1(PA, PB, p);
        };
        table_[3][4][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_4_2(PA, PB, p);
        };
        table_[3][4][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_4_3(PA, PB, p);
        };
        table_[3][4][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_4_4(PA, PB, p);
        };
        table_[3][4][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_4_5(PA, PB, p);
        };
        table_[3][4][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_4_6(PA, PB, p);
        };
        table_[3][4][7] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_3_4_7(PA, PB, p);
        };
        table_[4][0][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_0_0(PA, PB, p);
        };
        table_[4][0][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_0_1(PA, PB, p);
        };
        table_[4][0][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_0_2(PA, PB, p);
        };
        table_[4][0][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_0_3(PA, PB, p);
        };
        table_[4][0][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_0_4(PA, PB, p);
        };
        table_[4][1][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_1_0(PA, PB, p);
        };
        table_[4][1][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_1_1(PA, PB, p);
        };
        table_[4][1][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_1_2(PA, PB, p);
        };
        table_[4][1][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_1_3(PA, PB, p);
        };
        table_[4][1][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_1_4(PA, PB, p);
        };
        table_[4][1][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_1_5(PA, PB, p);
        };
        table_[4][2][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_2_0(PA, PB, p);
        };
        table_[4][2][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_2_1(PA, PB, p);
        };
        table_[4][2][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_2_2(PA, PB, p);
        };
        table_[4][2][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_2_3(PA, PB, p);
        };
        table_[4][2][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_2_4(PA, PB, p);
        };
        table_[4][2][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_2_5(PA, PB, p);
        };
        table_[4][2][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_2_6(PA, PB, p);
        };
        table_[4][3][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_3_0(PA, PB, p);
        };
        table_[4][3][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_3_1(PA, PB, p);
        };
        table_[4][3][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_3_2(PA, PB, p);
        };
        table_[4][3][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_3_3(PA, PB, p);
        };
        table_[4][3][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_3_4(PA, PB, p);
        };
        table_[4][3][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_3_5(PA, PB, p);
        };
        table_[4][3][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_3_6(PA, PB, p);
        };
        table_[4][3][7] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_3_7(PA, PB, p);
        };
        table_[4][4][0] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_4_0(PA, PB, p);
        };
        table_[4][4][1] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_4_1(PA, PB, p);
        };
        table_[4][4][2] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_4_2(PA, PB, p);
        };
        table_[4][4][3] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_4_3(PA, PB, p);
        };
        table_[4][4][4] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_4_4(PA, PB, p);
        };
        table_[4][4][5] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_4_5(PA, PB, p);
        };
        table_[4][4][6] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_4_6(PA, PB, p);
        };
        table_[4][4][7] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_4_7(PA, PB, p);
        };
        table_[4][4][8] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
            return symbolic::hermite_dE_dPA_4_4_8(PA, PB, p);
        };
    }

    DispatchTable table_;
};

// =============================================================================
// CONVENIENCE FUNCTIONS
// =============================================================================

inline Vec8d computeTMP(int nA, int nB, int t,
                        const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
    return TMPDispatcher::instance().compute(nA, nB, t, PA, PB, p);
}

inline Vec8d computeSymbolic(int nA, int nB, int t,
                             const Vec8d& PA, const Vec8d& PB, const Vec8d& one_over_2p) {
    return SymbolicDispatcher::instance().compute(nA, nB, t, PA, PB, one_over_2p);
}

inline Vec8d computeTMPGradPA(int nA, int nB, int t,
                              const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {
    return TMPGradPADispatcher::instance().compute(nA, nB, t, PA, PB, p);
}

inline Vec8d computeSymbolicGradPA(int nA, int nB, int t,
                                   const Vec8d& PA, const Vec8d& PB, const Vec8d& one_over_2p) {
    return SymbolicGradPADispatcher::instance().compute(nA, nB, t, PA, PB, one_over_2p);
}

} // namespace benchmark
} // namespace recursum
