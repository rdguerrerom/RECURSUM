#pragma once

/**
 * @file coulomb_r_symbolic.hpp
 * @brief Optimized Coulomb R auxiliary integrals using symbolic expansion + CSE
 *
 * This is the symbolic/codegen approach (impl=3, equivalent to "LayeredCodegen") for
 * Coulomb R integrals. It uses SymPy symbolic expansion followed by Common
 * Subexpression Elimination to generate highly optimized code.
 *
 * Analogous to the Symbolic implementation (impl=2) for Hermite E coefficients,
 * but this is considered the "LayeredCodegen" approach for plotting consistency.
 */

#include <recursum/mcmd/coulomb_r_symbolic_L0.hpp>
#include <recursum/mcmd/coulomb_r_symbolic_L1.hpp>
#include <recursum/mcmd/coulomb_r_symbolic_L2.hpp>
#include <recursum/mcmd/coulomb_r_symbolic_L3.hpp>
#include <recursum/mcmd/coulomb_r_symbolic_L4.hpp>

namespace recursum {
namespace mcmd {

/**
 * Unified interface for Coulomb R symbolic computations.
 *
 * Usage:
 *   std::array<Vec8d, N> R_array;
 *   CoulombRSymbolic<L_total>::compute(R_array.data(), PCx, PCy, PCz, Boys);
 */
template<int L_total>
using CoulombRSymbolic = coulomb_r_symbolic::CoulombRSymbolic<L_total>;

} // namespace mcmd
} // namespace recursum
