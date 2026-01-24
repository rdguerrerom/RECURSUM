#pragma once

#include <type_traits>
#include <recursum/vectorclass.h>

// Portable force-inline macro for performance-critical compute methods
// This ensures template instantiations are fully inlined at compile time
#ifndef RECURSUM_FORCEINLINE
  #ifdef _MSC_VER
    #define RECURSUM_FORCEINLINE __forceinline
  #elif defined(__GNUC__) || defined(__clang__)
    #define RECURSUM_FORCEINLINE inline __attribute__((always_inline))
  #else
    #define RECURSUM_FORCEINLINE inline
  #endif
#endif

namespace mcmd_coulomb {


template<int t, int u, int v, typename Enable = void>
struct CoulombRCoeffLayer {
    static constexpr int N_VALUES = 0;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d /*PCx*/, Vec8d /*PCy*/, Vec8d /*PCz*/, const Vec8d* /*Boys*/) {
        // Invalid indices: do nothing
    }
};

template<>
struct CoulombRCoeffLayer<0, 0, 0, void> {
    static constexpr int N_VALUES = 4;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d /*PCx*/, Vec8d /*PCy*/, Vec8d /*PCz*/, const Vec8d* /*Boys*/) {
        // Base case: R_{000}^{(N)} = Boys[N] for all N
        for (int N = 0; N < N_VALUES; ++N) {
            out[N] = Boys[N];
        }
    }
};

template<int t, int u, int v>
struct CoulombRCoeffLayer<
    t, u, v,
    typename std::enable_if<((t > 0)) && ((t >= 0) && (u >= 0) && (v >= 0))>::type
> {
    static constexpr int N_VALUES = (t + u + v) + 4;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
        // Compute previous layer 0: (t - 1, u, v)
        Vec8d prev_0[N_VALUES + 1] = {};  // Zero-init, sized for shifted access
        CoulombRCoeffLayer<t - 1, u, v>::compute(prev_0, PCx, PCy, PCz, Boys);
        // Compute previous layer 1: (t - 2, u, v)
        Vec8d prev_1[N_VALUES + 1] = {};  // Zero-init, sized for shifted access
        CoulombRCoeffLayer<t - 2, u, v>::compute(prev_1, PCx, PCy, PCz, Boys);

        // General case: t = 0 to N_VALUES-1
        for (int N = 0; N < N_VALUES; ++N) {
            out[N] = PCx * prev_0[N + 1] + Vec8d(t - 1) * prev_1[N + 1];
        }
    }
};

template<int t, int u, int v>
struct CoulombRCoeffLayer<
    t, u, v,
    typename std::enable_if<((t == 0) && (u > 0)) && ((t >= 0) && (u >= 0) && (v >= 0))>::type
> {
    static constexpr int N_VALUES = (t + u + v) + 4;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
        // Compute previous layer 0: (t, u - 1, v)
        Vec8d prev_0[N_VALUES + 1] = {};  // Zero-init, sized for shifted access
        CoulombRCoeffLayer<t, u - 1, v>::compute(prev_0, PCx, PCy, PCz, Boys);
        // Compute previous layer 1: (t, u - 2, v)
        Vec8d prev_1[N_VALUES + 1] = {};  // Zero-init, sized for shifted access
        CoulombRCoeffLayer<t, u - 2, v>::compute(prev_1, PCx, PCy, PCz, Boys);

        // General case: t = 0 to N_VALUES-1
        for (int N = 0; N < N_VALUES; ++N) {
            out[N] = PCy * prev_0[N + 1] + Vec8d(u - 1) * prev_1[N + 1];
        }
    }
};

template<int t, int u, int v>
struct CoulombRCoeffLayer<
    t, u, v,
    typename std::enable_if<((t == 0) && (u == 0) && (v > 0)) && ((t >= 0) && (u >= 0) && (v >= 0))>::type
> {
    static constexpr int N_VALUES = (t + u + v) + 4;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
        // Compute previous layer 0: (t, u, v - 1)
        Vec8d prev_0[N_VALUES + 1] = {};  // Zero-init, sized for shifted access
        CoulombRCoeffLayer<t, u, v - 1>::compute(prev_0, PCx, PCy, PCz, Boys);
        // Compute previous layer 1: (t, u, v - 2)
        Vec8d prev_1[N_VALUES + 1] = {};  // Zero-init, sized for shifted access
        CoulombRCoeffLayer<t, u, v - 2>::compute(prev_1, PCx, PCy, PCz, Boys);

        // General case: t = 0 to N_VALUES-1
        for (int N = 0; N < N_VALUES; ++N) {
            out[N] = PCz * prev_0[N + 1] + Vec8d(v - 1) * prev_1[N + 1];
        }
    }
};

// API compatibility: single-value accessor
template<int t, int u, int v, int N>
struct CoulombRCoeff {
    static RECURSUM_FORCEINLINE Vec8d compute(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
        Vec8d layer[(t + u + v) + 4];
        CoulombRCoeffLayer<t, u, v>::compute(layer, PCx, PCy, PCz, Boys);
        return layer[N];
    }
};

} // namespace mcmd_coulomb