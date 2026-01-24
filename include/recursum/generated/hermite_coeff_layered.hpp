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

namespace mcmd_hermite {


template<int nA, int nB, typename Enable = void>
struct HermiteECoeffLayer {
    static constexpr int N_VALUES = 0;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/) {
        // Invalid indices: do nothing
    }
};

template<>
struct HermiteECoeffLayer<0, 0, void> {
    static constexpr int N_VALUES = 1;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/) {
        out[0] = Vec8d(1.0);
    }
};

template<int nA, int nB>
struct HermiteECoeffLayer<
    nA, nB,
    typename std::enable_if<((nA > 0) && (nB == 0)) && ((nA >= 0) && (nB >= 0))>::type
> {
    static constexpr int N_VALUES = (nA + nB) + 1;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d PA, Vec8d PB, Vec8d p) {
        // Compute previous layer once (initialize to zero for out-of-bounds accesses)
        Vec8d prev[nA + nB + 2] = {};  // Zero-initialize
        HermiteECoeffLayer<nA - 1, nB>::compute(prev, PA, PB, p);

        // t = 0 special case
        out[0] = PA * prev[0] + Vec8d(1) * prev[1];

        // General case: t = 1 to N_VALUES-1
        for (int t = 1; t < N_VALUES; ++t) {
            out[t] = 0.5 / p * prev[t - 1] + PA * prev[t] + Vec8d(t + 1) * prev[t + 1];
        }
    }
};

template<int nA, int nB>
struct HermiteECoeffLayer<
    nA, nB,
    typename std::enable_if<((nA == 0) && (nB > 0)) && ((nA >= 0) && (nB >= 0))>::type
> {
    static constexpr int N_VALUES = (nA + nB) + 1;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d PA, Vec8d PB, Vec8d p) {
        // Compute previous layer once (initialize to zero for out-of-bounds accesses)
        Vec8d prev[nA + nB + 2] = {};  // Zero-initialize
        HermiteECoeffLayer<nA, nB - 1>::compute(prev, PA, PB, p);

        // t = 0 special case
        out[0] = PB * prev[0] + Vec8d(1) * prev[1];

        // General case: t = 1 to N_VALUES-1
        for (int t = 1; t < N_VALUES; ++t) {
            out[t] = 0.5 / p * prev[t - 1] + PB * prev[t] + Vec8d(t + 1) * prev[t + 1];
        }
    }
};

template<int nA, int nB>
struct HermiteECoeffLayer<
    nA, nB,
    typename std::enable_if<((nA > 0) && (nB > 0)) && ((nA >= 0) && (nB >= 0))>::type
> {
    static constexpr int N_VALUES = (nA + nB) + 1;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d PA, Vec8d PB, Vec8d p) {
        // Compute previous layer once (initialize to zero for out-of-bounds accesses)
        Vec8d prev[nA + nB + 2] = {};  // Zero-initialize
        HermiteECoeffLayer<nA - 1, nB>::compute(prev, PA, PB, p);

        // t = 0 special case
        out[0] = PA * prev[0] + Vec8d(1) * prev[1];

        // General case: t = 1 to N_VALUES-1
        for (int t = 1; t < N_VALUES; ++t) {
            out[t] = 0.5 / p * prev[t - 1] + PA * prev[t] + Vec8d(t + 1) * prev[t + 1];
        }
    }
};

// API compatibility: single-value accessor
template<int nA, int nB, int t>
struct HermiteECoeff {
    static RECURSUM_FORCEINLINE Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        Vec8d layer[(nA + nB) + 1];
        HermiteECoeffLayer<nA, nB>::compute(layer, PA, PB, p);
        return layer[t];
    }
};

} // namespace mcmd_hermite