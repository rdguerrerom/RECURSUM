/**
 * @file hermite_e_layered.hpp
 * @brief Layer-by-layer Hermite E computation with TRUE CSE
 *
 * This implementation computes ALL t-values for a given (nA, nB) simultaneously,
 * enabling genuine common subexpression elimination. Each "layer" depends only
 * on the previous layer, computed exactly once.
 *
 * Key optimization: E^{nA,nB} for all t values is computed from E^{nA-1,nB}
 * which is computed only ONCE and reused for all t values.
 *
 * Performance: O(nA * (nA+nB)) instead of O(3^(nA+nB)) for naive recursion.
 */

#pragma once

#include <array>
#include <type_traits>

#ifndef RECURSUM_VEC_TYPE
#include <recursum/vectorclass.h>
#endif

namespace recursum {
namespace mcmd {
namespace layered {

// Maximum supported angular momentum
constexpr int MAX_L = 4;
constexpr int MAX_T = 2 * MAX_L + 1;  // Max t index + 1

/**
 * @brief Layer holding all E^{nA,nB}_t for t = 0..nA+nB
 *
 * This struct holds a complete "row" of coefficients, enabling CSE
 * since the previous layer is computed only once.
 */
template<int nA, int nB>
struct HermiteELayer {
    static constexpr int N = nA + nB + 1;  // Number of valid t values

    /**
     * @brief Compute all E^{nA,nB}_t values at once
     * @return Array of E values indexed by t
     */
    static std::array<Vec8d, MAX_T> compute(Vec8d PA, Vec8d PB, Vec8d p) {
        std::array<Vec8d, MAX_T> result{};

        // Get previous layer (computed ONCE, reused for all t)
        auto prev = HermiteELayer<nA - 1, nB>::compute(PA, PB, p);

        Vec8d inv2p = Vec8d(0.5) / p;

        // Compute each t value using the recurrence:
        // E^{nA,nB}_t = inv2p * E^{nA-1,nB}_{t-1} + PA * E^{nA-1,nB}_t + (t+1) * E^{nA-1,nB}_{t+1}

        // t = 0 special case (no t-1 term)
        result[0] = PA * prev[0] + prev[1];

        // t = 1 to nA+nB-1
        for (int t = 1; t < nA + nB; ++t) {
            result[t] = inv2p * prev[t-1] + PA * prev[t] + Vec8d(t + 1) * prev[t + 1];
        }

        // t = nA+nB (no t+1 term)
        if (nA + nB > 0) {
            int t = nA + nB;
            result[t] = inv2p * prev[t-1] + PA * prev[t];
        }

        return result;
    }
};

// Base case: E^{0,0} layer
template<>
struct HermiteELayer<0, 0> {
    static std::array<Vec8d, MAX_T> compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/) {
        std::array<Vec8d, MAX_T> result{};
        result[0] = Vec8d(1.0);  // E^{0,0}_0 = 1
        return result;
    }
};

// B-side layers: E^{0,nB}
template<int nB>
struct HermiteELayer<0, nB> {
    static std::array<Vec8d, MAX_T> compute(Vec8d PA, Vec8d PB, Vec8d p) {
        static_assert(nB > 0, "nB must be positive");

        std::array<Vec8d, MAX_T> result{};
        auto prev = HermiteELayer<0, nB - 1>::compute(PA, PB, p);

        Vec8d inv2p = Vec8d(0.5) / p;

        // t = 0
        result[0] = PB * prev[0] + prev[1];

        // t = 1 to nB-1
        for (int t = 1; t < nB; ++t) {
            result[t] = inv2p * prev[t-1] + PB * prev[t] + Vec8d(t + 1) * prev[t + 1];
        }

        // t = nB
        result[nB] = inv2p * prev[nB-1] + PB * prev[nB];

        return result;
    }
};

/**
 * @brief Optimized single-coefficient accessor using layered computation
 *
 * This maintains API compatibility with HermiteE<nA, nB, t>::compute()
 * but uses the layered approach internally for efficiency.
 */
template<int nA, int nB, int t>
struct HermiteEOpt {
    static_assert(t >= 0 && t <= nA + nB, "Invalid t index");

    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        auto layer = HermiteELayer<nA, nB>::compute(PA, PB, p);
        return layer[t];
    }
};

// Convenience function
template<int nA, int nB, int t>
inline Vec8d hermite_e_opt(Vec8d PA, Vec8d PB, Vec8d p) {
    return HermiteEOpt<nA, nB, t>::compute(PA, PB, p);
}

} // namespace layered
} // namespace mcmd
} // namespace recursum
