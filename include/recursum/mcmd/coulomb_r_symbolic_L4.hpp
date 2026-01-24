#pragma once

#include <array>
#include <recursum/vectorclass.h>

namespace recursum {
namespace mcmd {
namespace coulomb_r_symbolic {

// Forward declaration of primary template
template<int L_total>
struct CoulombRSymbolic;

/**
 * Coulomb R auxiliary integrals for L_max = 4
 *
 * Computes all R_{tuv}^{(0)} for t+u+v <= 4
 *
 * Operations Analysis:
 *   Initial: 213 ops
 *   Final (with CSE): 213 ops
 *   Reduction: 0.0%
 */

template<>
struct CoulombRSymbolic<4> {
    static constexpr int N_VALUES = 35;

    static inline void compute(Vec8d* R, Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {

        // Common multipliers
        constexpr double m0 = 3;
        constexpr double m1 = 6;

        // Calculate R values
        R[0] = Boys[0];
        R[1] = Boys[1]*PCx;
        R[2] = Boys[1]*PCy;
        R[3] = Boys[1]*PCz;
        R[4] = Boys[1] + Boys[2]*PCx*PCx;
        R[5] = Boys[2]*PCx*PCy;
        R[6] = Boys[2]*PCx*PCz;
        R[7] = Boys[1] + Boys[2]*PCy*PCy;
        R[8] = Boys[2]*PCy*PCz;
        R[9] = Boys[1] + Boys[2]*PCz*PCz;
        R[10] = PCx*(m0*Boys[2] + Boys[3]*PCx*PCx);
        R[11] = PCy*(Boys[2] + Boys[3]*PCx*PCx);
        R[12] = PCz*(Boys[2] + Boys[3]*PCx*PCx);
        R[13] = PCx*(Boys[2] + Boys[3]*PCy*PCy);
        R[14] = Boys[3]*PCx*PCy*PCz;
        R[15] = PCx*(Boys[2] + Boys[3]*PCz*PCz);
        R[16] = PCy*(m0*Boys[2] + Boys[3]*PCy*PCy);
        R[17] = PCz*(Boys[2] + Boys[3]*PCy*PCy);
        R[18] = PCy*(Boys[2] + Boys[3]*PCz*PCz);
        R[19] = PCz*(m0*Boys[2] + Boys[3]*PCz*PCz);
        R[20] = m0*Boys[2] +m1*Boys[3]*PCx*PCx + Boys[4]*PCx*PCx*PCx*PCx;
        R[21] = PCx*PCy*(m0*Boys[3] + Boys[4]*PCx*PCx);
        R[22] = PCx*PCz*(m0*Boys[3] + Boys[4]*PCx*PCx);
        R[23] = Boys[2] + Boys[3]*PCx*PCx + Boys[3]*PCy*PCy + Boys[4]*PCx*PCx*PCy*PCy;
        R[24] = PCy*PCz*(Boys[3] + Boys[4]*PCx*PCx);
        R[25] = Boys[2] + Boys[3]*PCx*PCx + Boys[3]*PCz*PCz + Boys[4]*PCx*PCx*PCz*PCz;
        R[26] = PCx*PCy*(m0*Boys[3] + Boys[4]*PCy*PCy);
        R[27] = PCx*PCz*(Boys[3] + Boys[4]*PCy*PCy);
        R[28] = PCx*PCy*(Boys[3] + Boys[4]*PCz*PCz);
        R[29] = PCx*PCz*(m0*Boys[3] + Boys[4]*PCz*PCz);
        R[30] = m0*Boys[2] +m1*Boys[3]*PCy*PCy + Boys[4]*PCy*PCy*PCy*PCy;
        R[31] = PCy*PCz*(m0*Boys[3] + Boys[4]*PCy*PCy);
        R[32] = Boys[2] + Boys[3]*PCy*PCy + Boys[3]*PCz*PCz + Boys[4]*PCy*PCy*PCz*PCz;
        R[33] = PCy*PCz*(m0*Boys[3] + Boys[4]*PCz*PCz);
        R[34] = m0*Boys[2] +m1*Boys[3]*PCz*PCz + Boys[4]*PCz*PCz*PCz*PCz;
    }
};

} // namespace coulomb_r_symbolic
} // namespace mcmd
} // namespace recursum
