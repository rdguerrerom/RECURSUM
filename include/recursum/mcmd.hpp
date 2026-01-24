/**
 * @file mcmd.hpp
 * @brief Unified McMurchie-Davidson integral framework for quantum chemistry
 *
 * This header provides the complete RECURSUM implementation of McMurchie-Davidson
 * recurrences for Gaussian integral evaluation. It serves as a cleaner alternative
 * to the individual solver headers in McMD:
 *
 *   - mcmd/hermite_e.hpp    -> Alternative to McMD/coeff_solver.hpp
 *   - mcmd/hermite_grad.hpp -> Alternative to McMD/grad_solver.hpp
 *   - mcmd/kinetic.hpp      -> McMD kinetic formula (TeraChem Eq. 13)
 *   - mcmd/coulomb_aux.hpp  -> Alternative to McMD/coulomb_solver.hpp
 *
 * Architecture Principles:
 * -----------------------
 * 1. **Unified base case handling**: All recurrences start from E^{0,0}_0 = 1
 * 2. **Direct derivative formulas**: Use ∂E/∂PA = i × E^{i-1,j}_t (Helgaker-Taylor)
 * 3. **Clear chain rules**: Nuclear derivatives via ∂E/∂A = (∂E/∂PA)(∂PA/∂A) + ...
 * 4. **Single-path recurrence**: Increment-i only (no two-branch averaging)
 * 5. **Translational invariance**: ∂E/∂A + ∂E/∂B = 0 is automatic
 *
 * Key Differences from McMD:
 * -------------------------
 * | Feature                  | McMD                    | RECURSUM               |
 * |--------------------------|-------------------------|------------------------|
 * | E coefficient parameter  | aAB = 1/(2p)           | p (combined exponent)  |
 * | General case recurrence  | Two-branch averaging    | Single-path (increment-i)|
 * | Derivative formulas      | Recurrence-based        | Direct (i × E^{i-1,j}) |
 * | Kinetic integral         | Obara-Saika             | McMD formula (Eq. 13)  |
 *
 * Validation:
 * ----------
 * All formulas have been validated against:
 * 1. Finite difference derivatives (< 1e-14 error)
 * 2. McMD recurrence-based approach (exact agreement)
 * 3. PySCF integral values (working McMD implementation)
 *
 * Usage:
 * -----
 * @code
 * #include <recursum/mcmd.hpp>
 *
 * using namespace recursum::mcmd;
 *
 * // Compute E coefficient
 * Vec8d E = hermite_e<1, 1, 0>(PA, PB, p);
 *
 * // Compute direct derivative
 * Vec8d dE_dPA = hermite_deriv_PA<1, 1, 0>(PA, PB, p);
 *
 * // Compute nuclear derivative
 * Vec8d dE_dA = hermite_deriv_A<1, 1, 0>(PA, PB, p, alpha_A, alpha_B);
 *
 * // Compute kinetic integral (normalized externally)
 * Vec8d T = kinetic_integral<0,0,0, 0,0,0>(alpha_A, alpha_B, ...);
 *
 * // Compute Coulomb auxiliary integral
 * Vec8d R = coulomb_r<1, 0, 0>(PCx, PCy, PCz, Boys);
 * @endcode
 *
 * Reference:
 *   - Helgaker & Taylor, Theor. Chim. Acta 83 (1992) 177-183
 *   - McMurchie & Davidson, J. Comput. Phys. 26 (1978) 218-231
 *   - TeraChem SI (F_orbitals_article.pdf), Eq. S12-S13
 */

#pragma once

// Core E coefficient solver
#include "mcmd/hermite_e.hpp"

// Gradient solver (direct formulas + chain rule)
#include "mcmd/hermite_grad.hpp"

// Kinetic energy integral (McMD formula)
#include "mcmd/kinetic.hpp"

// Coulomb auxiliary integrals
#include "mcmd/coulomb_aux.hpp"
