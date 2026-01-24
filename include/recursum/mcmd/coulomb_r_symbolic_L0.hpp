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
 * Coulomb R auxiliary integrals for L_max = 0
 *
 * Computes all R_{tuv}^{(0)} for t+u+v <= 0
 *
 * Operations Analysis:
 *   Initial: 0 ops
 *   Final (with CSE): 0 ops
 *   Reduction: 100.0%
 */

template<>
struct CoulombRSymbolic<0> {
    static constexpr int N_VALUES = 1;

    static inline void compute(Vec8d* R, Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {

        // Calculate R values
        R[0] = Boys[0];
    }
};

} // namespace coulomb_r_symbolic
} // namespace mcmd
} // namespace recursum
