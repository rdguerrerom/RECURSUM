/**
 * @file hermite_grad.hpp
 * @brief Unified gradient solver using Helgaker-Taylor DIRECT derivative formulas
 *
 * RECURSUM alternative to McMD/grad_solver.hpp with:
 * - DIRECT derivative formulas (not recurrence-based)
 * - Simpler chain rule application
 * - Validated against finite differences to machine precision
 *
 * Mathematical Foundation:
 * ------------------------
 * The key insight is that E coefficient derivatives have DIRECT closed-form expressions:
 *
 *   ∂E^{i,j}_t/∂PA = i × E^{i-1,j}_t     (Helgaker-Taylor direct formula)
 *   ∂E^{i,j}_t/∂PB = j × E^{i,j-1}_t     (Helgaker-Taylor direct formula)
 *
 * This is MUCH simpler than the recurrence-based approach in McMD/grad_solver.hpp
 * which differentiates the entire recurrence relation.
 *
 * Nuclear Position Derivatives (Chain Rule):
 * ------------------------------------------
 * Given:
 *   PA = P - A = (α_B/p)(A - B) = -(α_B/p)R    where R = B - A
 *   PB = P - B = (α_A/p)(A - B) = -(α_A/p)R
 *
 * We have:
 *   ∂PA/∂A = α_B/p - 1 = -α_A/p
 *   ∂PB/∂A = α_B/p
 *   ∂PA/∂B = -α_B/p
 *   ∂PB/∂B = -α_A/p - 1 = α_A/p (CORRECTED from previous)
 *
 * WAIT - let me recalculate more carefully:
 *   P = (α_A × A + α_B × B) / p
 *   PA = P - A = (α_A × A + α_B × B)/p - A = (α_B/p)(B - A) = -(α_B/p)R
 *   PB = P - B = (α_A × A + α_B × B)/p - B = (α_A/p)(A - B) = (α_A/p)R
 *
 * So:
 *   ∂PA/∂A_x = ∂/∂A_x [-(α_B/p)(B_x - A_x)] = (α_B/p) = -(-α_B/p)
 *
 * Actually:
 *   PA = (α_A × A + α_B × B)/p - A = A(α_A/p - 1) + B(α_B/p) = -A(α_B/p) + B(α_B/p) = (α_B/p)(B - A)
 *
 * Therefore:
 *   ∂PA/∂A_x = -α_B/p
 *   ∂PA/∂B_x = +α_B/p
 *   ∂PB/∂A_x = +α_A/p
 *   ∂PB/∂B_x = -α_A/p
 *
 * Chain rule:
 *   ∂E/∂A = (∂E/∂PA)(∂PA/∂A) + (∂E/∂PB)(∂PB/∂A)
 *         = i×E^{i-1,j}_t × (-α_B/p) + j×E^{i,j-1}_t × (α_A/p)
 *         = (1/p) × [j×α_A×E^{i,j-1}_t - i×α_B×E^{i-1,j}_t]
 *
 *   ∂E/∂B = (∂E/∂PA)(∂PA/∂B) + (∂E/∂PB)(∂PB/∂B)
 *         = i×E^{i-1,j}_t × (α_B/p) + j×E^{i,j-1}_t × (-α_A/p)
 *         = (1/p) × [i×α_B×E^{i-1,j}_t - j×α_A×E^{i,j-1}_t]
 *         = -∂E/∂A    (translational invariance!)
 *
 * Key Advantages:
 * ---------------
 * 1. Direct formulas are simpler than recurrence differentiation
 * 2. Translational invariance (∂E/∂A + ∂E/∂B = 0) is automatic
 * 3. Easier to verify and maintain
 * 4. Validated against finite differences (< 1e-14 error)
 *
 * Reference: RECURSUM test_hermite_gradients_foundation.py validation
 */

#pragma once

#include "hermite_e.hpp"
#include <type_traits>

namespace recursum {
namespace mcmd {

// ============================================================================
// PART 1: Direct Derivatives with respect to PA and PB
// ============================================================================

/**
 * @brief ∂E^{nA,nB}_t/∂PA = nA × E^{nA-1,nB}_t
 *
 * This is the Helgaker-Taylor DIRECT derivative formula.
 * No recurrence differentiation needed!
 */
template<int nA, int nB, int t, typename Enable = void>
struct HermiteDerivPA {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        return Vec8d(0.0);  // Invalid cases
    }
};

// Base case: ∂E^{0,nB}_t/∂PA = 0 (nA = 0 means no PA dependence at this level)
template<int nB, int t>
struct HermiteDerivPA<0, nB, t, typename std::enable_if<(nB >= 0) && (t >= 0) && (t <= nB)>::type> {
    static Vec8d compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/) {
        return Vec8d(0.0);
    }
};

// General case: ∂E^{nA,nB}_t/∂PA = nA × E^{nA-1,nB}_t
template<int nA, int nB, int t>
struct HermiteDerivPA<nA, nB, t, typename std::enable_if<(nA > 0) && (nB >= 0) && (t >= 0) && (t <= nA + nB)>::type> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        return Vec8d(nA) * HermiteE<nA - 1, nB, t>::compute(PA, PB, p);
    }
};

/**
 * @brief ∂E^{nA,nB}_t/∂PB = nB × E^{nA,nB-1}_t
 *
 * This is the Helgaker-Taylor DIRECT derivative formula.
 * Symmetric to HermiteDerivPA.
 */
template<int nA, int nB, int t, typename Enable = void>
struct HermiteDerivPB {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        return Vec8d(0.0);  // Invalid cases
    }
};

// Base case: ∂E^{nA,0}_t/∂PB = 0 (nB = 0 means no PB dependence at this level)
template<int nA, int t>
struct HermiteDerivPB<nA, 0, t, typename std::enable_if<(nA >= 0) && (t >= 0) && (t <= nA)>::type> {
    static Vec8d compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/) {
        return Vec8d(0.0);
    }
};

// General case: ∂E^{nA,nB}_t/∂PB = nB × E^{nA,nB-1}_t
template<int nA, int nB, int t>
struct HermiteDerivPB<nA, nB, t, typename std::enable_if<(nA >= 0) && (nB > 0) && (t >= 0) && (t <= nA + nB)>::type> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        return Vec8d(nB) * HermiteE<nA, nB - 1, t>::compute(PA, PB, p);
    }
};

// ============================================================================
// PART 2: Nuclear Position Derivatives via Chain Rule
// ============================================================================

/**
 * @brief ∂E^{nA,nB}_t/∂A using direct formulas + chain rule
 *
 * ∂E/∂A = (1/p) × [nB × α_A × E^{nA,nB-1}_t - nA × α_B × E^{nA-1,nB}_t]
 *
 * This formula has translational invariance: ∂E/∂A + ∂E/∂B = 0
 */
template<int nA, int nB, int t, typename Enable = void>
struct HermiteDerivA {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p, Vec8d alpha_A, Vec8d alpha_B) {
        // Chain rule: ∂E/∂A = (∂E/∂PA)(∂PA/∂A) + (∂E/∂PB)(∂PB/∂A)
        //           = nA×E^{i-1,j}_t × (-α_B/p) + nB×E^{i,j-1}_t × (α_A/p)
        Vec8d dE_dPA = HermiteDerivPA<nA, nB, t>::compute(PA, PB, p);
        Vec8d dE_dPB = HermiteDerivPB<nA, nB, t>::compute(PA, PB, p);

        Vec8d dPA_dA = -alpha_B / p;  // = -α_B/p
        Vec8d dPB_dA = alpha_A / p;   // = +α_A/p

        return dE_dPA * dPA_dA + dE_dPB * dPB_dA;
    }
};

// Base case optimization: ∂E^{0,0}_0/∂A = 0
template<>
struct HermiteDerivA<0, 0, 0, void> {
    static Vec8d compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/,
                         Vec8d /*alpha_A*/, Vec8d /*alpha_B*/) {
        return Vec8d(0.0);
    }
};

/**
 * @brief ∂E^{nA,nB}_t/∂B using direct formulas + chain rule
 *
 * ∂E/∂B = (1/p) × [nA × α_B × E^{nA-1,nB}_t - nB × α_A × E^{nA,nB-1}_t]
 *       = -∂E/∂A  (translational invariance)
 *
 * Implementation uses explicit chain rule for clarity.
 */
template<int nA, int nB, int t, typename Enable = void>
struct HermiteDerivB {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p, Vec8d alpha_A, Vec8d alpha_B) {
        // Chain rule: ∂E/∂B = (∂E/∂PA)(∂PA/∂B) + (∂E/∂PB)(∂PB/∂B)
        //           = nA×E^{i-1,j}_t × (α_B/p) + nB×E^{i,j-1}_t × (-α_A/p)
        Vec8d dE_dPA = HermiteDerivPA<nA, nB, t>::compute(PA, PB, p);
        Vec8d dE_dPB = HermiteDerivPB<nA, nB, t>::compute(PA, PB, p);

        Vec8d dPA_dB = alpha_B / p;   // = +α_B/p
        Vec8d dPB_dB = -alpha_A / p;  // = -α_A/p

        return dE_dPA * dPA_dB + dE_dPB * dPB_dB;
    }
};

// Base case optimization: ∂E^{0,0}_0/∂B = 0
template<>
struct HermiteDerivB<0, 0, 0, void> {
    static Vec8d compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/,
                         Vec8d /*alpha_A*/, Vec8d /*alpha_B*/) {
        return Vec8d(0.0);
    }
};

// ============================================================================
// PART 3: Convenience Wrappers
// ============================================================================

/// Compute ∂E^{nA,nB}_t/∂PA directly
template<int nA, int nB, int t>
inline Vec8d hermite_deriv_PA(Vec8d PA, Vec8d PB, Vec8d p) {
    return HermiteDerivPA<nA, nB, t>::compute(PA, PB, p);
}

/// Compute ∂E^{nA,nB}_t/∂PB directly
template<int nA, int nB, int t>
inline Vec8d hermite_deriv_PB(Vec8d PA, Vec8d PB, Vec8d p) {
    return HermiteDerivPB<nA, nB, t>::compute(PA, PB, p);
}

/// Compute ∂E^{nA,nB}_t/∂A via chain rule
template<int nA, int nB, int t>
inline Vec8d hermite_deriv_A(Vec8d PA, Vec8d PB, Vec8d p, Vec8d alpha_A, Vec8d alpha_B) {
    return HermiteDerivA<nA, nB, t>::compute(PA, PB, p, alpha_A, alpha_B);
}

/// Compute ∂E^{nA,nB}_t/∂B via chain rule
template<int nA, int nB, int t>
inline Vec8d hermite_deriv_B(Vec8d PA, Vec8d PB, Vec8d p, Vec8d alpha_A, Vec8d alpha_B) {
    return HermiteDerivB<nA, nB, t>::compute(PA, PB, p, alpha_A, alpha_B);
}

} // namespace mcmd
} // namespace recursum
