/**
 * @file coulomb_aux.hpp
 * @brief Hermite Coulomb auxiliary integrals R_{tuv}^{(N)} and their derivatives
 *
 * RECURSUM alternative to McMD/coulomb_solver.hpp with:
 * - Unified namespace (recursum::mcmd)
 * - Cleaner documentation
 * - Consistent with hermite_e.hpp and hermite_grad.hpp
 *
 * Mathematical Foundation:
 * ------------------------
 * The Hermite Coulomb auxiliary integrals arise in:
 * - Nuclear attraction: ⟨μ|1/|r-C||ν⟩
 * - Electron repulsion: (μν|λσ)
 *
 * They are defined as:
 *   R_{tuv}^{(N)} = ∫ Λ_t(x) Λ_u(y) Λ_v(z) × (1/|r-C|) dr
 *
 * where Λ_t are Hermite Gaussians centered at the product center P.
 *
 * Recurrence Relations:
 * --------------------
 * Base case: R_{000}^{(N)} = F_N(T)  (Boys function)
 *   where T = p × |PC|² and PC = P - C
 *
 * X-recurrence (t > 0):
 *   R_{tuv}^{(N)} = PC_x × R_{t-1,u,v}^{(N+1)} + (t-1) × R_{t-2,u,v}^{(N+1)}
 *
 * Y-recurrence (t = 0, u > 0):
 *   R_{0uv}^{(N)} = PC_y × R_{0,u-1,v}^{(N+1)} + (u-1) × R_{0,u-2,v}^{(N+1)}
 *
 * Z-recurrence (t = u = 0, v > 0):
 *   R_{00v}^{(N)} = PC_z × R_{00,v-1}^{(N+1)} + (v-1) × R_{00,v-2}^{(N+1)}
 *
 * Derivative Recurrences:
 * ----------------------
 * The key insight is that the P-derivative of R integrals follows:
 *   ∂R_{tuv}/∂P_x = -R_{t+1,u,v}   (NEGATIVE sign!)
 *
 * This allows nuclear derivatives via chain rule:
 *   ∂V/∂A = (∂V/∂P)(∂P/∂A) + ...
 *
 * Reference:
 *   - McMurchie & Davidson, J. Comput. Phys. 26 (1978) 218-231
 *   - Helgaker et al., "Molecular Electronic-Structure Theory" Ch. 9
 */

#pragma once

#include <type_traits>

#ifndef RECURSUM_VEC_TYPE
#include <recursum/vectorclass.h>
#endif

namespace recursum {
namespace mcmd {

// ============================================================================
// PART 1: Auxiliary Integrals R_{tuv}^{(N)}
// ============================================================================

/**
 * @brief Hermite Coulomb auxiliary integral R_{t,u,v}^{(N)}
 *
 * @tparam t Hermite index in x-direction
 * @tparam u Hermite index in y-direction
 * @tparam v Hermite index in z-direction
 * @tparam N Boys function auxiliary index (starts at 0 for integrals)
 */
template<int t, int u, int v, int N, typename Enable = void>
struct CoulombR {
    static Vec8d compute(Vec8d /*PCx*/, Vec8d /*PCy*/, Vec8d /*PCz*/,
                         const Vec8d* /*Boys*/) {
        return Vec8d(0.0);  // Invalid cases
    }
};

// Base case: R_{000}^{(N)} = F_N(T)
template<int N>
struct CoulombR<0, 0, 0, N, typename std::enable_if<(N >= 0)>::type> {
    static Vec8d compute(Vec8d /*PCx*/, Vec8d /*PCy*/, Vec8d /*PCz*/,
                         const Vec8d* Boys) {
        return Boys[N];
    }
};

// X-recurrence: t > 0
// R_{tuv}^{(N)} = PC_x × R_{t-1,u,v}^{(N+1)} + (t-1) × R_{t-2,u,v}^{(N+1)}
template<int t, int u, int v, int N>
struct CoulombR<t, u, v, N,
    typename std::enable_if<(t > 0) && (u >= 0) && (v >= 0) && (N >= 0)>::type>
{
    static Vec8d compute(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
        Vec8d term1 = PCx * CoulombR<t - 1, u, v, N + 1>::compute(PCx, PCy, PCz, Boys);
        Vec8d term2 = Vec8d(t - 1) * CoulombR<t - 2, u, v, N + 1>::compute(PCx, PCy, PCz, Boys);
        return term1 + term2;
    }
};

// Y-recurrence: t = 0, u > 0
// R_{0uv}^{(N)} = PC_y × R_{0,u-1,v}^{(N+1)} + (u-1) × R_{0,u-2,v}^{(N+1)}
template<int u, int v, int N>
struct CoulombR<0, u, v, N,
    typename std::enable_if<(u > 0) && (v >= 0) && (N >= 0)>::type>
{
    static Vec8d compute(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
        Vec8d term1 = PCy * CoulombR<0, u - 1, v, N + 1>::compute(PCx, PCy, PCz, Boys);
        Vec8d term2 = Vec8d(u - 1) * CoulombR<0, u - 2, v, N + 1>::compute(PCx, PCy, PCz, Boys);
        return term1 + term2;
    }
};

// Z-recurrence: t = u = 0, v > 0
// R_{00v}^{(N)} = PC_z × R_{00,v-1}^{(N+1)} + (v-1) × R_{00,v-2}^{(N+1)}
template<int v, int N>
struct CoulombR<0, 0, v, N,
    typename std::enable_if<(v > 0) && (N >= 0)>::type>
{
    static Vec8d compute(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
        Vec8d term1 = PCz * CoulombR<0, 0, v - 1, N + 1>::compute(PCx, PCy, PCz, Boys);
        Vec8d term2 = Vec8d(v - 1) * CoulombR<0, 0, v - 2, N + 1>::compute(PCx, PCy, PCz, Boys);
        return term1 + term2;
    }
};

// ============================================================================
// PART 2: P-Derivatives Using Higher Angular Momentum
// ============================================================================

/**
 * @brief "Incremented angular momentum" form: R_{t+1,u,v}^{(N)}
 *
 * For gradient calculations, the derivative of V integral w.r.t. P_x
 * can be expressed using R integrals with incremented angular momentum:
 *
 *   ∂/∂P_x [Σ E_t R_t] = Σ E_t × ∂R_t/∂P_x = Σ E_t × (-R_{t+1})  (at fixed T)
 *
 * IMPORTANT: This formula assumes the Boys function argument T = p|PC|² is
 * held constant. In reality, T also depends on P, so the full derivative
 * requires the chain-rule approach in CoulombRDerivPCx.
 *
 * The "incremented angular momentum" form is useful for:
 * - Obara-Saika style gradient formulas
 * - Cases where T is treated as an independent variable
 */
template<int t, int u, int v, int N>
struct CoulombRIncX {
    static Vec8d compute(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
        // R_{t+1,u,v}^{(N)} - used in gradient formulas
        return CoulombR<t + 1, u, v, N>::compute(PCx, PCy, PCz, Boys);
    }
};

template<int t, int u, int v, int N>
struct CoulombRIncY {
    static Vec8d compute(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
        return CoulombR<t, u + 1, v, N>::compute(PCx, PCy, PCz, Boys);
    }
};

template<int t, int u, int v, int N>
struct CoulombRIncZ {
    static Vec8d compute(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
        return CoulombR<t, u, v + 1, N>::compute(PCx, PCy, PCz, Boys);
    }
};

// ============================================================================
// PART 3: Alternative Derivative Formulation (Chain Rule)
// ============================================================================

/**
 * @brief ∂R_{tuv}^{(N)}/∂(PC_x) using chain rule recurrence
 *
 * This is the derivative with respect to PC_x = P_x - C_x.
 * It uses the recurrence differentiation approach from McMD/coulomb_solver.hpp.
 *
 * For gradient calculations, the simpler form ∂R/∂P = -R_{t+1,...} is often
 * preferred, but this chain-rule form is kept for completeness.
 *
 * Base case: ∂R_{000}^{(N)}/∂(PC_x) = -2p × PC_x × F_{N+1}(T)
 *   (from chain rule on Boys function: dF_N/dT × dT/d(PC_x))
 *
 * X-recurrence (t > 0):
 *   ∂R_{tuv}/∂(PC_x) = R_{t-1,u,v}^{(N+1)}   <- EXTRA TERM from explicit PC_x!
 *                    + PC_x × ∂R_{t-1,u,v}^{(N+1)}/∂(PC_x)
 *                    + (t-1) × ∂R_{t-2,u,v}^{(N+1)}/∂(PC_x)
 */
template<int t, int u, int v, int N, typename Enable = void>
struct CoulombRDerivPCx {
    static Vec8d compute(Vec8d /*PCx*/, Vec8d /*PCy*/, Vec8d /*PCz*/,
                         const Vec8d* /*Boys*/, Vec8d /*p*/) {
        return Vec8d(0.0);
    }
};

// Base case: chain rule on Boys function
template<int N>
struct CoulombRDerivPCx<0, 0, 0, N, typename std::enable_if<(N >= 0)>::type> {
    static Vec8d compute(Vec8d PCx, Vec8d /*PCy*/, Vec8d /*PCz*/,
                         const Vec8d* Boys, Vec8d p) {
        // dF_N(T)/d(PC_x) = F'_N(T) × dT/d(PC_x)
        //                 = -F_{N+1}(T) × 2p × PC_x
        return Vec8d(-2.0) * p * PCx * Boys[N + 1];
    }
};

// X-recurrence with extra term
template<int t, int u, int v, int N>
struct CoulombRDerivPCx<t, u, v, N,
    typename std::enable_if<(t > 0) && (u >= 0) && (v >= 0) && (N >= 0)>::type>
{
    static Vec8d compute(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys, Vec8d p) {
        // Extra term from d(PC_x)/d(PC_x) = 1
        Vec8d extra = CoulombR<t - 1, u, v, N + 1>::compute(PCx, PCy, PCz, Boys);

        Vec8d term1 = PCx * CoulombRDerivPCx<t - 1, u, v, N + 1>::compute(PCx, PCy, PCz, Boys, p);
        Vec8d term2 = Vec8d(t - 1) * CoulombRDerivPCx<t - 2, u, v, N + 1>::compute(PCx, PCy, PCz, Boys, p);

        return extra + term1 + term2;
    }
};

// Y-recurrence (no extra term)
template<int u, int v, int N>
struct CoulombRDerivPCx<0, u, v, N,
    typename std::enable_if<(u > 0) && (v >= 0) && (N >= 0)>::type>
{
    static Vec8d compute(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys, Vec8d p) {
        Vec8d term1 = PCy * CoulombRDerivPCx<0, u - 1, v, N + 1>::compute(PCx, PCy, PCz, Boys, p);
        Vec8d term2 = Vec8d(u - 1) * CoulombRDerivPCx<0, u - 2, v, N + 1>::compute(PCx, PCy, PCz, Boys, p);
        return term1 + term2;
    }
};

// Z-recurrence (no extra term)
template<int v, int N>
struct CoulombRDerivPCx<0, 0, v, N,
    typename std::enable_if<(v > 0) && (N >= 0)>::type>
{
    static Vec8d compute(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys, Vec8d p) {
        Vec8d term1 = PCz * CoulombRDerivPCx<0, 0, v - 1, N + 1>::compute(PCx, PCy, PCz, Boys, p);
        Vec8d term2 = Vec8d(v - 1) * CoulombRDerivPCx<0, 0, v - 2, N + 1>::compute(PCx, PCy, PCz, Boys, p);
        return term1 + term2;
    }
};

// Similar structures for Y and Z derivatives (following coulomb_solver.hpp pattern)
// Omitted for brevity - the ∂R/∂P = -R_{t+1,...} form is recommended

// ============================================================================
// PART 4: Convenience Wrappers
// ============================================================================

/// Compute R_{tuv}^{(N)}
template<int t, int u, int v, int N = 0>
inline Vec8d coulomb_r(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
    return CoulombR<t, u, v, N>::compute(PCx, PCy, PCz, Boys);
}

/// Compute R_{t+1,u,v}^{(N)} (incremented x angular momentum)
template<int t, int u, int v, int N = 0>
inline Vec8d coulomb_r_inc_x(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
    return CoulombRIncX<t, u, v, N>::compute(PCx, PCy, PCz, Boys);
}

/// Compute R_{t,u+1,v}^{(N)} (incremented y angular momentum)
template<int t, int u, int v, int N = 0>
inline Vec8d coulomb_r_inc_y(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
    return CoulombRIncY<t, u, v, N>::compute(PCx, PCy, PCz, Boys);
}

/// Compute R_{t,u,v+1}^{(N)} (incremented z angular momentum)
template<int t, int u, int v, int N = 0>
inline Vec8d coulomb_r_inc_z(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys) {
    return CoulombRIncZ<t, u, v, N>::compute(PCx, PCy, PCz, Boys);
}

/// Compute ∂R_{tuv}^{(N)}/∂(PC_x) using chain rule (includes Boys function derivative)
template<int t, int u, int v, int N = 0>
inline Vec8d coulomb_r_deriv_PCx(Vec8d PCx, Vec8d PCy, Vec8d PCz, const Vec8d* Boys, Vec8d p) {
    return CoulombRDerivPCx<t, u, v, N>::compute(PCx, PCy, PCz, Boys, p);
}

} // namespace mcmd
} // namespace recursum
