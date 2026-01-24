/**
 * @file coulomb_r_layered.hpp
 * @brief Layer-by-layer Coulomb R integral computation with TRUE CSE
 *
 * This implementation computes ALL R_{tuv}^{(N)} values for a given maximum
 * angular momentum L_total simultaneously, enabling genuine common subexpression
 * elimination. Each N-layer depends only on the N+1 layer, computed exactly once.
 *
 * Key optimization: For a shell quartet with L_total = L_ab + L_cd, all R values
 * at auxiliary index N are computed from R values at N+1 which are computed
 * only ONCE and reused for all (t,u,v) combinations.
 *
 * Memory layout: Triangular indexing scheme for (t,u,v) with t+u+v <= L_max
 *   index(t,u,v) = tetrahedral_offset(t+u+v) + triangular_offset(u+v) + v
 *
 * Performance: O(L_total^3) instead of O(3^L_total) for naive recursion.
 *
 * Reference:
 *   - McMurchie & Davidson, J. Comput. Phys. 26 (1978) 218-231
 *   - Helgaker et al., "Molecular Electronic-Structure Theory" Ch. 9
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

// Maximum supported total angular momentum for R integrals
// For (gg|gg) integrals: L_total = 4 + 4 + 4 + 4 = 16
constexpr int MAX_R_L = 16;

// Tetrahedral number: count of (t,u,v) with t+u+v <= L
// T(L) = (L+1)(L+2)(L+3)/6
constexpr int tetrahedral(int L) {
    return (L + 1) * (L + 2) * (L + 3) / 6;
}

// Triangular number: count of (u,v) with u+v <= L
// Tri(L) = (L+1)(L+2)/2
constexpr int triangular(int L) {
    return (L + 1) * (L + 2) / 2;
}

// Maximum storage size for R array
constexpr int MAX_R_SIZE = tetrahedral(MAX_R_L);

/**
 * @brief Index into R array for (t,u,v) with t+u+v <= L_max
 *
 * The indexing scheme groups by total angular momentum L = t+u+v,
 * then by (u+v), then by v.
 */
constexpr int r_index(int t, int u, int v) {
    int L = t + u + v;
    // Offset to start of shell with total L
    int shell_offset = tetrahedral(L - 1);
    // Offset within shell for given (u+v)
    int uv = u + v;
    int uv_offset = triangular(uv - 1);
    return shell_offset + uv_offset + v;
}

/**
 * @brief Compute all R_{tuv}^{(N)} values for N=0 and t+u+v <= L_total
 *
 * @tparam L_total Maximum total angular momentum (t+u+v)
 *
 * The computation proceeds from N = L_total down to N = 0,
 * building up all required R values layer by layer.
 */
template<int L_total>
struct CoulombRLayer {
    static_assert(L_total >= 0 && L_total <= MAX_R_L, "L_total out of range");

    // Number of R values needed
    static constexpr int N_VALUES = tetrahedral(L_total);

    /**
     * @brief Compute all R_{tuv}^{(0)} values
     *
     * @param PCx, PCy, PCz  Distance from Gaussian product center P to point C
     * @param Boys           Pre-computed Boys function values F_0(T) to F_{L_total}(T)
     * @return Array of R values indexed by r_index(t,u,v)
     */
    static std::array<Vec8d, MAX_R_SIZE> compute(Vec8d PCx, Vec8d PCy, Vec8d PCz,
                                                  const Vec8d* Boys) {
        // We need two arrays to swap between N and N+1 layers
        std::array<Vec8d, MAX_R_SIZE> curr{};
        std::array<Vec8d, MAX_R_SIZE> prev{};

        // Initialize: R_{000}^{(L_total)} = F_{L_total}(T)
        prev[0] = Boys[L_total];

        // Build from N = L_total-1 down to N = 0
        for (int N = L_total - 1; N >= 0; --N) {
            // At level N, we need all (t,u,v) with t+u+v <= L_total - N
            int L_max_at_N = L_total - N;

            // Clear current layer
            for (int i = 0; i < tetrahedral(L_max_at_N); ++i) {
                curr[i] = Vec8d(0.0);
            }

            // Base case at this N: R_{000}^{(N)} = F_N(T)
            curr[0] = Boys[N];

            // Z-recurrence first: build R_{00v}^{(N)} for v > 0
            for (int v = 1; v <= L_max_at_N; ++v) {
                // R_{00v}^{(N)} = PCz * R_{00,v-1}^{(N+1)} + (v-1) * R_{00,v-2}^{(N+1)}
                Vec8d term1 = PCz * prev[r_index(0, 0, v - 1)];
                Vec8d term2 = (v > 1) ? Vec8d(v - 1) * prev[r_index(0, 0, v - 2)] : Vec8d(0.0);
                curr[r_index(0, 0, v)] = term1 + term2;
            }

            // Y-recurrence: build R_{0uv}^{(N)} for u > 0
            for (int u = 1; u <= L_max_at_N; ++u) {
                for (int v = 0; v <= L_max_at_N - u; ++v) {
                    // R_{0uv}^{(N)} = PCy * R_{0,u-1,v}^{(N+1)} + (u-1) * R_{0,u-2,v}^{(N+1)}
                    Vec8d term1 = PCy * prev[r_index(0, u - 1, v)];
                    Vec8d term2 = (u > 1) ? Vec8d(u - 1) * prev[r_index(0, u - 2, v)] : Vec8d(0.0);
                    curr[r_index(0, u, v)] = term1 + term2;
                }
            }

            // X-recurrence: build R_{tuv}^{(N)} for t > 0
            for (int t = 1; t <= L_max_at_N; ++t) {
                for (int u = 0; u <= L_max_at_N - t; ++u) {
                    for (int v = 0; v <= L_max_at_N - t - u; ++v) {
                        // R_{tuv}^{(N)} = PCx * R_{t-1,u,v}^{(N+1)} + (t-1) * R_{t-2,u,v}^{(N+1)}
                        Vec8d term1 = PCx * prev[r_index(t - 1, u, v)];
                        Vec8d term2 = (t > 1) ? Vec8d(t - 1) * prev[r_index(t - 2, u, v)] : Vec8d(0.0);
                        curr[r_index(t, u, v)] = term1 + term2;
                    }
                }
            }

            // Swap arrays for next iteration
            std::swap(curr, prev);
        }

        // After the loop, results are in prev (due to final swap)
        return prev;
    }
};

// Specialization for L_total = 0
template<>
struct CoulombRLayer<0> {
    static constexpr int N_VALUES = 1;

    static std::array<Vec8d, MAX_R_SIZE> compute(Vec8d /*PCx*/, Vec8d /*PCy*/, Vec8d /*PCz*/,
                                                  const Vec8d* Boys) {
        std::array<Vec8d, MAX_R_SIZE> result{};
        result[0] = Boys[0];  // R_{000}^{(0)} = F_0(T)
        return result;
    }
};

/**
 * @brief Convenience function to compute a single R_{tuv}^{(0)} using layered approach
 *
 * @note For best performance when multiple R values are needed, use CoulombRLayer
 *       directly to avoid recomputing the entire layer for each query.
 */
template<int t, int u, int v>
inline Vec8d coulomb_r_layered(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
    constexpr int L_total = t + u + v;
    auto R = CoulombRLayer<L_total>::compute(PCx, PCy, PCz, Boys);
    return R[r_index(t, u, v)];
}

/**
 * @brief Get R value from pre-computed layer array
 */
template<int t, int u, int v>
inline Vec8d get_r(const std::array<Vec8d, MAX_R_SIZE>& R) {
    return R[r_index(t, u, v)];
}

// Non-template version for runtime indexing
inline Vec8d get_r(const std::array<Vec8d, MAX_R_SIZE>& R, int t, int u, int v) {
    return R[r_index(t, u, v)];
}

} // namespace layered
} // namespace mcmd
} // namespace recursum
