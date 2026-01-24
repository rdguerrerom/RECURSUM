/**
 * @file hermite_e.hpp
 * @brief Unified Hermite E coefficient solver using Helgaker-Taylor 1992 recurrence
 *
 * RECURSUM alternative to McMD/coeff_solver.hpp with:
 * - Cleaner base case handling (E^{0,0}_0 = 1)
 * - Single-path recurrence (increment-i only) for consistency
 * - Proper handling of t=0 case including E_{t+1} term
 *
 * Mathematical Foundation (Helgaker-Taylor 1992, Eq. 7):
 * -------------------------------------------------------
 * E^{i+1,j}_t = (1/2p) × E^{i,j}_{t-1} + PA × E^{i,j}_t + (t+1) × E^{i,j}_{t+1}
 *
 * For t=0:
 * E^{i+1,j}_0 = PA × E^{i,j}_0 + E^{i,j}_1
 *
 * Base case: E^{0,0}_0 = 1
 *
 * This is the CORRECT formulation that includes the (t+1)×E_{t+1} term
 * even when t=0, which was a critical bug fix in the RECURSUM validation.
 *
 * Key differences from McMD/coeff_solver.hpp:
 * ------------------------------------------
 * 1. Uses `p` (combined exponent) instead of `aAB = 1/(2p)`
 * 2. Single-path recurrence (always decrement nA first) instead of two-branch averaging
 * 3. Explicit handling of t=0 vs t>0 cases
 *
 * Reference: Helgaker & Taylor, Theor. Chim. Acta 83 (1992) 177-183
 */

#pragma once

#include <type_traits>

#ifndef RECURSUM_VEC_TYPE
#include <recursum/vectorclass.h>
#endif

namespace recursum {
namespace mcmd {

/**
 * @brief Hermite expansion coefficient E^{nA,nB}_t
 *
 * Template parameters:
 * @tparam nA Angular momentum on center A (i)
 * @tparam nB Angular momentum on center B (j)
 * @tparam t  Hermite auxiliary index (0 <= t <= nA + nB)
 */
template<int nA, int nB, int t, typename Enable = void>
struct HermiteE {
    static Vec8d compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/) {
        return Vec8d(0.0);  // Invalid parameter combinations
    }
};

// Base case: E^{0,0}_0 = 1
template<>
struct HermiteE<0, 0, 0, void> {
    static Vec8d compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/) {
        return Vec8d(1.0);
    }
};

// A-side recurrence for t=0 (nA > 0, nB = 0)
// E^{i+1,0}_0 = PA × E^{i,0}_0 + E^{i,0}_1
template<int nA, int t>
struct HermiteE<
    nA, 0, t,
    typename std::enable_if<(nA > 0) && (t == 0)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        return PA * HermiteE<nA - 1, 0, t>::compute(PA, PB, p)
             + Vec8d(t + 1) * HermiteE<nA - 1, 0, t + 1>::compute(PA, PB, p);
    }
};

// A-side recurrence for t>0 (nA > 0, nB = 0)
// E^{i+1,0}_t = (1/2p) × E^{i,0}_{t-1} + PA × E^{i,0}_t + (t+1) × E^{i,0}_{t+1}
template<int nA, int t>
struct HermiteE<
    nA, 0, t,
    typename std::enable_if<(nA > 0) && (t > 0) && (t <= nA)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        return Vec8d(0.5) / p * HermiteE<nA - 1, 0, t - 1>::compute(PA, PB, p)
             + PA * HermiteE<nA - 1, 0, t>::compute(PA, PB, p)
             + Vec8d(t + 1) * HermiteE<nA - 1, 0, t + 1>::compute(PA, PB, p);
    }
};

// B-side recurrence for t=0 (nA = 0, nB > 0)
// E^{0,j+1}_0 = PB × E^{0,j}_0 + E^{0,j}_1
template<int nB, int t>
struct HermiteE<
    0, nB, t,
    typename std::enable_if<(nB > 0) && (t == 0)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        return PB * HermiteE<0, nB - 1, t>::compute(PA, PB, p)
             + Vec8d(t + 1) * HermiteE<0, nB - 1, t + 1>::compute(PA, PB, p);
    }
};

// B-side recurrence for t>0 (nA = 0, nB > 0)
// E^{0,j+1}_t = (1/2p) × E^{0,j}_{t-1} + PB × E^{0,j}_t + (t+1) × E^{0,j}_{t+1}
template<int nB, int t>
struct HermiteE<
    0, nB, t,
    typename std::enable_if<(nB > 0) && (t > 0) && (t <= nB)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        return Vec8d(0.5) / p * HermiteE<0, nB - 1, t - 1>::compute(PA, PB, p)
             + PB * HermiteE<0, nB - 1, t>::compute(PA, PB, p)
             + Vec8d(t + 1) * HermiteE<0, nB - 1, t + 1>::compute(PA, PB, p);
    }
};

// General case for t=0 (nA > 0, nB > 0)
// Use increment-i (decrement nA) recurrence only
// E^{i+1,j}_0 = PA × E^{i,j}_0 + E^{i,j}_1
template<int nA, int nB, int t>
struct HermiteE<
    nA, nB, t,
    typename std::enable_if<(nA > 0) && (nB > 0) && (t == 0)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        return PA * HermiteE<nA - 1, nB, t>::compute(PA, PB, p)
             + Vec8d(t + 1) * HermiteE<nA - 1, nB, t + 1>::compute(PA, PB, p);
    }
};

// General case for t>0 (nA > 0, nB > 0)
// Use increment-i (decrement nA) recurrence only
// E^{i+1,j}_t = (1/2p) × E^{i,j}_{t-1} + PA × E^{i,j}_t + (t+1) × E^{i,j}_{t+1}
template<int nA, int nB, int t>
struct HermiteE<
    nA, nB, t,
    typename std::enable_if<(nA > 0) && (nB > 0) && (t > 0) && (t <= nA + nB)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        return Vec8d(0.5) / p * HermiteE<nA - 1, nB, t - 1>::compute(PA, PB, p)
             + PA * HermiteE<nA - 1, nB, t>::compute(PA, PB, p)
             + Vec8d(t + 1) * HermiteE<nA - 1, nB, t + 1>::compute(PA, PB, p);
    }
};

// Convenience function wrapper
template<int nA, int nB, int t>
inline Vec8d hermite_e(Vec8d PA, Vec8d PB, Vec8d p) {
    return HermiteE<nA, nB, t>::compute(PA, PB, p);
}

} // namespace mcmd
} // namespace recursum
