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
 * Coulomb R auxiliary integrals for L_max = 1
 *
 * Computes all R_{tuv}^{(0)} for t+u+v <= 1
 *
 * Operations Analysis:
 *   Initial: 3 ops
 *   Final (with CSE): 3 ops
 *   Reduction: 0.0%
 */

template<>
struct CoulombRSymbolic<1> {
    static constexpr int N_VALUES = 4;

    static inline void compute(Vec8d* R, Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {

        // Calculate R values
        R[0] = Boys[0];
        R[1] = Boys[1]*PCx;
        R[2] = Boys[1]*PCy;
        R[3] = Boys[1]*PCz;
    }
};

} // namespace coulomb_r_symbolic
} // namespace mcmd
} // namespace recursum
