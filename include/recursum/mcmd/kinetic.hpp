/**
 * @file kinetic.hpp
 * @brief Kinetic energy integral using McMurchie-Davidson formula
 *
 * Implements the correct McMD kinetic energy formula from TeraChem SI (Eq. 13)
 * as an alternative to the Obara-Saika approach currently used in McMD.
 *
 * Mathematical Foundation (TeraChem SI Eq. S12-S13):
 * -------------------------------------------------
 * The kinetic energy integral T = ⟨μ|-½∇²|ν⟩ in the McMD scheme is:
 *
 * T_μν = C_μ × C_ν × K_AB × (π/p)^{3/2} × [
 *     kinetic_x × E_y × E_z +
 *     E_x × kinetic_y × E_z +
 *     E_x × E_y × kinetic_z
 * ]
 *
 * where each kinetic term applies the kinetic operator to one Cartesian direction:
 *
 *   kinetic(i,j) = -j(j-1)/2 × E^{i,j-2}_0 + (2j+1)×b × E^{i,j}_0 - 2b² × E^{i,j+2}_0
 *
 * Here:
 *   - b = α_B (exponent of the ket Gaussian)
 *   - E^{i,j}_0 is the Hermite coefficient with t=0
 *   - C_μ, C_ν are Cartesian normalization constants
 *   - K_AB = exp(-μ × R²_AB) is the Gaussian product factor
 *   - p = α_A + α_B is the combined exponent
 *
 * Key Insight:
 * -----------
 * The kinetic operator acts on the ket function (B-side), which is why the
 * formula involves shifting j (angular momentum on B) by ±2. The three terms
 * correspond to:
 *   1. -j(j-1)/2 × E^{i,j-2}_0: Lowering j by 2 (second derivative brings down j)
 *   2. (2j+1)×b × E^{i,j}_0: Diagonal contribution
 *   3. -2b² × E^{i,j+2}_0: Raising j by 2 (laplacian action on Gaussian)
 *
 * Comparison with Obara-Saika:
 * ---------------------------
 * The Obara-Saika formula (McMD/kinetic_obara_saika.cpp) uses:
 *   T = b×(2L+3)×S - 2b²×[S^{+2x} + S^{+2y} + S^{+2z}] - 0.5×[l(l-1)×S^{-2}]
 *
 * Both approaches are mathematically equivalent but:
 * - Obara-Saika computes overlaps with shifted angular momenta
 * - McMD uses Hermite coefficients directly (more unified with S, V integrals)
 *
 * Reference:
 *   - TeraChem SI (F_orbitals_article.pdf), Eq. S12-S13
 *   - python/gradient_reference/core/kinetic_integral_mcmd.py
 */

#pragma once

#include "hermite_e.hpp"
#include <cmath>

namespace recursum {
namespace mcmd {

/**
 * @brief Compute the 1D kinetic contribution for McMD formula
 *
 * kinetic(nA, nB) = -nB(nB-1)/2 × E^{nA,nB-2}_0 + (2nB+1)×b × E^{nA,nB}_0 - 2b² × E^{nA,nB+2}_0
 *
 * This computes the kinetic operator contribution in one Cartesian direction.
 *
 * @tparam nA Angular momentum on A in this direction
 * @tparam nB Angular momentum on B in this direction
 */
template<int nA, int nB>
struct KineticTerm1D {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p, Vec8d alpha_B) {
        Vec8d result(0.0);

        // Term 1: -nB(nB-1)/2 × E^{nA,nB-2}_0
        // Only contributes if nB >= 2
        if constexpr (nB >= 2) {
            result += Vec8d(-0.5 * nB * (nB - 1)) * HermiteE<nA, nB - 2, 0>::compute(PA, PB, p);
        }

        // Term 2: (2nB+1)×b × E^{nA,nB}_0
        result += Vec8d(2 * nB + 1) * alpha_B * HermiteE<nA, nB, 0>::compute(PA, PB, p);

        // Term 3: -2b² × E^{nA,nB+2}_0
        result += Vec8d(-2.0) * alpha_B * alpha_B * HermiteE<nA, nB + 2, 0>::compute(PA, PB, p);

        return result;
    }
};

/**
 * @brief Compute the 3D kinetic integral using McMD formula
 *
 * T = K_AB × (π/p)^{3/2} × [kinetic_x × E_y × E_z +
 *                           E_x × kinetic_y × E_z +
 *                           E_x × E_y × kinetic_z]
 *
 * Note: Normalization factors C_μ, C_ν should be applied externally.
 *
 * @tparam lA,mA,nA Cartesian angular momentum on A (x,y,z)
 * @tparam lB,mB,nB Cartesian angular momentum on B (x,y,z)
 */
template<int lA, int mA, int nA, int lB, int mB, int nB>
struct KineticIntegral3D {
    static Vec8d compute(Vec8d PAx, Vec8d PAy, Vec8d PAz,
                         Vec8d PBx, Vec8d PBy, Vec8d PBz,
                         Vec8d p, Vec8d alpha_B, Vec8d K_AB) {
        // E coefficients in each direction (t=0 for overlap integrals)
        Vec8d Ex = HermiteE<lA, lB, 0>::compute(PAx, PBx, p);
        Vec8d Ey = HermiteE<mA, mB, 0>::compute(PAy, PBy, p);
        Vec8d Ez = HermiteE<nA, nB, 0>::compute(PAz, PBz, p);

        // Kinetic contributions in each direction
        Vec8d Tx = KineticTerm1D<lA, lB>::compute(PAx, PBx, p, alpha_B);
        Vec8d Ty = KineticTerm1D<mA, mB>::compute(PAy, PBy, p, alpha_B);
        Vec8d Tz = KineticTerm1D<nA, nB>::compute(PAz, PBz, p, alpha_B);

        // Sum: Tx×Ey×Ez + Ex×Ty×Ez + Ex×Ey×Tz
        Vec8d T_value = Tx * Ey * Ez + Ex * Ty * Ez + Ex * Ey * Tz;

        // Apply Gaussian product prefactor: K_AB × (π/p)^{3/2}
        constexpr double PI = 3.14159265358979323846;
        Vec8d prefactor = K_AB * pow(Vec8d(PI) / p, Vec8d(1.5));

        return prefactor * T_value;
    }
};

/**
 * @brief Convenience function for kinetic integral computation
 *
 * Computes the primitive kinetic integral between two Cartesian Gaussians.
 * The normalization factors are NOT included - they should be applied externally
 * using GaussianNormalization::cartesianNorm3D().
 *
 * @param lA,mA,nA Cartesian angular momentum on A
 * @param lB,mB,nB Cartesian angular momentum on B
 * @param alpha_A  Exponent on center A
 * @param alpha_B  Exponent on center B
 * @param A_x,A_y,A_z Position of center A
 * @param B_x,B_y,B_z Position of center B
 */
template<int lA, int mA, int nA, int lB, int mB, int nB>
inline Vec8d kinetic_integral(
    Vec8d alpha_A, Vec8d alpha_B,
    Vec8d A_x, Vec8d A_y, Vec8d A_z,
    Vec8d B_x, Vec8d B_y, Vec8d B_z
) {
    // Combined exponent
    Vec8d p = alpha_A + alpha_B;

    // Gaussian product center P = (α_A × A + α_B × B) / p
    Vec8d P_x = (alpha_A * A_x + alpha_B * B_x) / p;
    Vec8d P_y = (alpha_A * A_y + alpha_B * B_y) / p;
    Vec8d P_z = (alpha_A * A_z + alpha_B * B_z) / p;

    // Displacement vectors PA and PB
    Vec8d PAx = P_x - A_x;
    Vec8d PAy = P_y - A_y;
    Vec8d PAz = P_z - A_z;
    Vec8d PBx = P_x - B_x;
    Vec8d PBy = P_y - B_y;
    Vec8d PBz = P_z - B_z;

    // Gaussian product factor K_AB = exp(-μ × R²_AB)
    Vec8d mu = alpha_A * alpha_B / p;
    Vec8d R_AB_sq = (A_x - B_x) * (A_x - B_x)
                  + (A_y - B_y) * (A_y - B_y)
                  + (A_z - B_z) * (A_z - B_z);
    Vec8d K_AB = exp(-mu * R_AB_sq);

    return KineticIntegral3D<lA, mA, nA, lB, mB, nB>::compute(
        PAx, PAy, PAz, PBx, PBy, PBz, p, alpha_B, K_AB
    );
}

} // namespace mcmd
} // namespace recursum
