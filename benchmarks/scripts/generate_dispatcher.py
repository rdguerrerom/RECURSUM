#!/usr/bin/env python3
"""
Generate comprehensive dispatcher header for all coefficients up to L=4 (gg).
"""

import os

def generate_dispatcher_header(max_L: int = 4):
    """Generate benchmark_dispatcher.hpp with all coefficients."""

    max_N = 2 * max_L  # Maximum t index

    # Collect all valid (nA, nB, t) combinations
    coeffs = []
    for nA in range(max_L + 1):
        for nB in range(max_L + 1):
            for t in range(nA + nB + 1):
                coeffs.append((nA, nB, t))

    header = f'''/**
 * @file benchmark_dispatcher.hpp
 * @brief Comprehensive dispatchers for TMP and Symbolic Hermite E benchmarks
 *
 * AUTO-GENERATED - DO NOT EDIT MANUALLY
 *
 * Coverage: E^{{0,0}}_0 to E^{{{max_L},{max_L}}}_{{{2*max_L}}} ({len(coeffs)} coefficients)
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

namespace recursum {{
namespace benchmark {{

// Maximum angular momentum for benchmarks
constexpr int BENCH_MAX_L = {max_L};  // Up to g-type (L=4)
constexpr int BENCH_MAX_N = 2 * BENCH_MAX_L;  // Max t index = {max_N}

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
class TMPDispatcher {{
public:
    static constexpr int MAX_NA = BENCH_MAX_L;
    static constexpr int MAX_NB = BENCH_MAX_L;
    static constexpr int MAX_N = BENCH_MAX_N;

    using NTable = std::array<CoeffFunction, MAX_N + 1>;
    using NBTable = std::array<NTable, MAX_NB + 1>;
    using DispatchTable = std::array<NBTable, MAX_NA + 1>;

    static const TMPDispatcher& instance() {{
        static TMPDispatcher dispatcher;
        return dispatcher;
    }}

    Vec8d compute(int nA, int nB, int t,
                  const Vec8d& PA, const Vec8d& PB, const Vec8d& p) const {{
        if (nA < 0 || nA > MAX_NA || nB < 0 || nB > MAX_NB ||
            t < 0 || t > MAX_N || t > nA + nB) {{
            return Vec8d(0.0);
        }}
        const auto& func = table_[nA][nB][t];
        if (!func) return Vec8d(0.0);
        return func(PA, PB, p);
    }}

    TMPDispatcher(const TMPDispatcher&) = delete;
    TMPDispatcher& operator=(const TMPDispatcher&) = delete;

private:
    TMPDispatcher() {{
        for (auto& na_slice : table_) {{
            for (auto& nb_slice : na_slice) {{
                for (auto& func : nb_slice) {{
                    func = nullptr;
                }}
            }}
        }}
        registerAll();
    }}

    template<int NA, int NB, int T>
    void registerCoeff() {{
        table_[NA][NB][T] = [](const Vec8d& pa, const Vec8d& pb, const Vec8d& p) -> Vec8d {{
            return mcmd::HermiteE<NA, NB, T>::compute(pa, pb, p);
        }};
    }}

    template<int NA, int NB, int T>
    struct RegisterT {{
        static void apply(TMPDispatcher* d) {{
            if constexpr (T <= NA + NB && T <= MAX_N) {{
                d->registerCoeff<NA, NB, T>();
                if constexpr (T + 1 <= NA + NB && T + 1 <= MAX_N) {{
                    RegisterT<NA, NB, T + 1>::apply(d);
                }}
            }}
        }}
    }};

    template<int NA, int NB>
    struct RegisterNB {{
        static void apply(TMPDispatcher* d) {{
            RegisterT<NA, NB, 0>::apply(d);
            if constexpr (NB + 1 <= MAX_NB) {{
                RegisterNB<NA, NB + 1>::apply(d);
            }}
        }}
    }};

    template<int NA>
    struct RegisterNA {{
        static void apply(TMPDispatcher* d) {{
            RegisterNB<NA, 0>::apply(d);
            if constexpr (NA + 1 <= MAX_NA) {{
                RegisterNA<NA + 1>::apply(d);
            }}
        }}
    }};

    void registerAll() {{
        RegisterNA<0>::apply(this);
    }}

    DispatchTable table_;
}};

// =============================================================================
// SYMBOLIC DISPATCHER
// =============================================================================

/**
 * @class SymbolicDispatcher
 * @brief Dispatcher for Symbolic (SymPy-generated) functions
 */
class SymbolicDispatcher {{
public:
    static constexpr int MAX_NA = BENCH_MAX_L;
    static constexpr int MAX_NB = BENCH_MAX_L;
    static constexpr int MAX_N = BENCH_MAX_N;

    using NTable = std::array<CoeffFunction, MAX_N + 1>;
    using NBTable = std::array<NTable, MAX_NB + 1>;
    using DispatchTable = std::array<NBTable, MAX_NA + 1>;

    static const SymbolicDispatcher& instance() {{
        static SymbolicDispatcher dispatcher;
        return dispatcher;
    }}

    Vec8d compute(int nA, int nB, int t,
                  const Vec8d& PA, const Vec8d& PB, const Vec8d& one_over_2p) const {{
        if (nA < 0 || nA > MAX_NA || nB < 0 || nB > MAX_NB ||
            t < 0 || t > MAX_N || t > nA + nB) {{
            return Vec8d(0.0);
        }}
        const auto& func = table_[nA][nB][t];
        if (!func) return Vec8d(0.0);
        return func(PA, PB, one_over_2p);
    }}

    SymbolicDispatcher(const SymbolicDispatcher&) = delete;
    SymbolicDispatcher& operator=(const SymbolicDispatcher&) = delete;

private:
    SymbolicDispatcher() {{
        for (auto& na_slice : table_) {{
            for (auto& nb_slice : na_slice) {{
                for (auto& func : nb_slice) {{
                    func = nullptr;
                }}
            }}
        }}
        registerAll();
    }}

    void registerAll() {{
'''

    # Generate registration for all coefficients
    for nA, nB, t in coeffs:
        header += f'''        table_[{nA}][{nB}][{t}] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {{
            return symbolic::hermite_e_symbolic_{nA}_{nB}_{t}(PA, PB, p);
        }};
'''

    header += '''    }

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
'''

    # Generate registration for all gradient coefficients
    for nA, nB, t in coeffs:
        header += f'''        table_[{nA}][{nB}][{t}] = [](const Vec8d& PA, const Vec8d& PB, const Vec8d& p) {{
            return symbolic::hermite_dE_dPA_{nA}_{nB}_{t}(PA, PB, p);
        }};
'''

    header += '''    }

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
'''

    return header


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, '../include/benchmark_dispatcher.hpp')

    print("Generating comprehensive dispatcher header...")
    header = generate_dispatcher_header(max_L=4)

    with open(output_path, 'w') as f:
        f.write(header)

    print(f"Generated: {output_path}")


if __name__ == '__main__':
    main()
