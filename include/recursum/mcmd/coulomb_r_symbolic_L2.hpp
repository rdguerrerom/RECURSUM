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
 * Coulomb R auxiliary integrals for L_max = 2
 *
 * Computes all R_{tuv}^{(0)} for t+u+v <= 2
 *
 * Operations Analysis:
 *   Initial: 21 ops
 *   Final (with CSE): 21 ops
 *   Reduction: 0.0%
 */

template<>
struct CoulombRSymbolic<2> {
    static constexpr int N_VALUES = 10;

    static inline void compute(Vec8d* R, Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {

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
    }
};

} // namespace coulomb_r_symbolic
} // namespace mcmd
} // namespace recursum
