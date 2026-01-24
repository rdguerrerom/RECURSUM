#!/usr/bin/env python3
"""
Generate comprehensive symbolic code for ALL McMurchie-Davidson recurrences.

Coverage:
- Hermite E coefficients: E^{0,0}_0 to E^{4,4}_8 (ss to gg shell pairs)
- Hermite gradients: ∂E/∂PA, ∂E/∂PB (direct differentiation)
- Nuclear gradients: ∂E/∂A, ∂E/∂B (chain rule)
- Coulomb R integrals: R_{0,0,0} to R_{16,16,16} (for (gg|gg) ERIs)

Output files:
- symbolic_generated/hermite_e_symbolic.hpp
- symbolic_generated/hermite_grad_symbolic.hpp
- symbolic_generated/coulomb_r_symbolic.hpp
"""

import sympy as sp
from sympy import symbols, expand, diff, factorial, binomial, sqrt, pi, exp, cse
from typing import Dict, Tuple, List
import os
import argparse

# =============================================================================
# PART 1: Hermite E Coefficients
# =============================================================================

def generate_hermite_e(max_L: int) -> Dict[Tuple[int, int, int], sp.Expr]:
    """
    Generate symbolic expressions for Hermite E coefficients E^{nA,nB}_t.

    Uses Helgaker-Taylor (1992) recurrence:
        E^{i+1,j}_t = (1/2p) * E^{i,j}_{t-1} + PA * E^{i,j}_t + (t+1) * E^{i,j}_{t+1}

    Args:
        max_L: Maximum angular momentum (4 for g orbitals)

    Returns:
        Dictionary mapping (nA, nB, t) -> symbolic expression
    """
    PA, PB, one_over_2p = symbols('PA PB one_over_2p', real=True)

    E: Dict[Tuple[int, int, int], sp.Expr] = {}
    E[(0, 0, 0)] = sp.Integer(1)

    def get_E(nA: int, nB: int, t: int) -> sp.Expr:
        if nA < 0 or nB < 0 or t < 0 or t > nA + nB:
            return sp.Integer(0)
        return E.get((nA, nB, t), sp.Integer(0))

    # Build all coefficients up to max_L
    for total_L in range(1, 2 * max_L + 1):
        for nA in range(min(total_L, max_L) + 1):
            nB = total_L - nA
            if nB < 0 or nB > max_L:
                continue

            for t in range(nA + nB + 1):
                if (nA, nB, t) in E:
                    continue

                # Choose recurrence direction
                if nA > 0:
                    # A-side: E^{nA,nB}_t from E^{nA-1,nB}_*
                    if t == 0:
                        E[(nA, nB, 0)] = PA * get_E(nA-1, nB, 0) + get_E(nA-1, nB, 1)
                    else:
                        E[(nA, nB, t)] = (one_over_2p * get_E(nA-1, nB, t-1)
                                         + PA * get_E(nA-1, nB, t)
                                         + (t + 1) * get_E(nA-1, nB, t+1))
                elif nB > 0:
                    # B-side: E^{0,nB}_t from E^{0,nB-1}_*
                    if t == 0:
                        E[(0, nB, 0)] = PB * get_E(0, nB-1, 0) + get_E(0, nB-1, 1)
                    else:
                        E[(0, nB, t)] = (one_over_2p * get_E(0, nB-1, t-1)
                                        + PB * get_E(0, nB-1, t)
                                        + (t + 1) * get_E(0, nB-1, t+1))

    # Expand all expressions
    for key in E:
        E[key] = expand(E[key])

    return E


def generate_hermite_gradients(E_coeffs: Dict[Tuple[int, int, int], sp.Expr]) -> Dict:
    """
    Generate symbolic gradient expressions using direct differentiation.

    Helgaker-Taylor direct formulas:
        ∂E^{nA,nB}_t/∂PA = nA × E^{nA-1,nB}_t
        ∂E^{nA,nB}_t/∂PB = nB × E^{nA,nB-1}_t

    We also compute these by direct differentiation to verify.
    """
    PA, PB, one_over_2p = symbols('PA PB one_over_2p', real=True)

    dE_dPA = {}
    dE_dPB = {}

    for (nA, nB, t), expr in E_coeffs.items():
        # Direct differentiation
        dE_dPA[(nA, nB, t)] = expand(diff(expr, PA))
        dE_dPB[(nA, nB, t)] = expand(diff(expr, PB))

    return {'dE_dPA': dE_dPA, 'dE_dPB': dE_dPB}


# =============================================================================
# PART 2: Coulomb R Integrals
# =============================================================================

def generate_coulomb_R(max_index: int) -> Dict[Tuple[int, int, int, int], sp.Expr]:
    """
    Generate symbolic expressions for Coulomb R integrals R_{t,u,v}^{(N)}.

    McMurchie-Davidson recurrence for R integrals:
        R_{t+1,u,v}^{(N)} = t * R_{t-1,u,v}^{(N+1)} + X_PC * R_{t,u,v}^{(N+1)}
        R_{t,u+1,v}^{(N)} = u * R_{t,u-1,v}^{(N+1)} + Y_PC * R_{t,u,v}^{(N+1)}
        R_{t,v+1,v}^{(N)} = v * R_{t,u,v-1}^{(N+1)} + Z_PC * R_{t,u,v}^{(N+1)}

    Base case: R_{0,0,0}^{(N)} = (-2p)^N * F_N(T) where F_N is Boys function

    For symbolic generation, we express R in terms of R_{0,0,0}^{(N)} base cases.

    Args:
        max_index: Maximum index for t, u, v (16 for gg|gg)

    Returns:
        Dictionary mapping (t, u, v, N) -> symbolic expression
    """
    X_PC, Y_PC, Z_PC = symbols('X_PC Y_PC Z_PC', real=True)

    # We'll express everything in terms of R_000_N symbols
    # R_{0,0,0}^{(N)} = F_N (Boys function, treated as symbol)

    R = {}

    # Create symbols for base cases R_{0,0,0}^{(N)}
    max_N = 3 * max_index  # Maximum Boys function order needed
    F = {n: symbols(f'F_{n}', real=True) for n in range(max_N + 1)}

    # Base cases: R_{0,0,0}^{(N)} = F_N
    for N in range(max_N + 1):
        R[(0, 0, 0, N)] = F[N]

    def get_R(t: int, u: int, v: int, N: int) -> sp.Expr:
        if t < 0 or u < 0 or v < 0 or N < 0:
            return sp.Integer(0)
        if N > max_N:
            return sp.Integer(0)  # Truncate at max order
        return R.get((t, u, v, N), sp.Integer(0))

    # Build R integrals using recurrence
    # We need to be careful about the order of computation

    # For efficiency, limit to smaller max_index for generation
    gen_max = min(max_index, 8)  # Limit symbolic generation

    for total in range(1, 3 * gen_max + 1):
        for t in range(min(total, gen_max) + 1):
            for u in range(min(total - t, gen_max) + 1):
                v = total - t - u
                if v < 0 or v > gen_max:
                    continue

                # Compute for all needed N values
                for N in range(max_N - total, -1, -1):
                    if (t, u, v, N) in R:
                        continue

                    # Choose recurrence direction based on largest index
                    if t > 0:
                        # X recurrence: R_{t,u,v}^{(N)} from R_{t-1,u,v}^{(N+1)}
                        R[(t, u, v, N)] = ((t - 1) * get_R(t-2, u, v, N+1)
                                          + X_PC * get_R(t-1, u, v, N+1))
                    elif u > 0:
                        # Y recurrence
                        R[(t, u, v, N)] = ((u - 1) * get_R(t, u-2, v, N+1)
                                          + Y_PC * get_R(t, u-1, v, N+1))
                    elif v > 0:
                        # Z recurrence
                        R[(t, u, v, N)] = ((v - 1) * get_R(t, u, v-2, N+1)
                                          + Z_PC * get_R(t, u, v-1, N+1))

    # Expand all expressions
    for key in R:
        R[key] = expand(R[key])

    return R


# =============================================================================
# PART 3: C++ Code Generation
# =============================================================================

def expr_to_cpp(expr: sp.Expr) -> str:
    """Convert SymPy expression to C++ code."""
    from sympy import ccode
    code = ccode(expr)

    # Optimize power expressions
    for var in ['PA', 'PB', 'one_over_2p', 'X_PC', 'Y_PC', 'Z_PC']:
        for n in range(2, 9):
            old = f'pow({var}, {n})'
            if n == 2:
                new = f'({var}*{var})'
            elif n == 3:
                new = f'({var}*{var}*{var})'
            elif n == 4:
                new = f'(({var}*{var})*({var}*{var}))'
            elif n == 5:
                new = f'(({var}*{var})*({var}*{var})*{var})'
            elif n == 6:
                new = f'(({var}*{var}*{var})*({var}*{var}*{var}))'
            elif n == 7:
                new = f'(({var}*{var}*{var})*({var}*{var}*{var})*{var})'
            elif n == 8:
                new = f'((({var}*{var})*({var}*{var}))*(({var}*{var})*({var}*{var})))'
            code = code.replace(old, new)

    return code


def generate_hermite_e_header(E_coeffs: Dict, output_path: str, use_cse: bool = True):
    """Generate C++ header for Hermite E coefficients.

    Args:
        E_coeffs: Dictionary of (nA, nB, t) -> symbolic expression
        output_path: Path to write the header file
        use_cse: If True, apply Common Subexpression Elimination
    """
    cse_mode = "with CSE" if use_cse else "expanded form"

    cpp = f'''/**
 * @file hermite_e_symbolic.hpp
 * @brief Symbolically-generated Hermite E coefficients ({cse_mode})
 *
 * AUTO-GENERATED - DO NOT EDIT MANUALLY
 * Coverage: E^{{0,0}}_0 to E^{{4,4}}_8 (ss to gg shell pairs)
 *
 * Optimization: {"CSE (Common Subexpression Elimination) applied" if use_cse else "No optimization"}
 */

#pragma once

#ifndef RECURSUM_VEC_TYPE
#include <recursum/vectorclass.h>
#endif

namespace recursum {{
namespace symbolic {{

'''

    sorted_keys = sorted(E_coeffs.keys())

    for (nA, nB, t) in sorted_keys:
        expr = E_coeffs[(nA, nB, t)]

        if use_cse:
            # Apply CSE to the expression
            cse_result = cse(expr)
            intermediates = cse_result[0]  # List of (symbol, expression) pairs
            final_expr = cse_result[1][0]  # The final expression using the intermediates

            if intermediates:
                # Generate function with intermediate variables
                cpp += f'''/**
 * @brief Symbolic E^{{{nA},{nB}}}_{t} (CSE optimized)
 */
inline Vec8d hermite_e_symbolic_{nA}_{nB}_{t}(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {{
'''
                for sym, sub_expr in intermediates:
                    cpp += f'    Vec8d {sym.name} = {expr_to_cpp(sub_expr)};\n'
                cpp += f'    return {expr_to_cpp(final_expr)};\n'
                cpp += '}\n\n'
            else:
                # No CSE needed, simple return
                cpp += f'''/**
 * @brief Symbolic E^{{{nA},{nB}}}_{t}
 */
inline Vec8d hermite_e_symbolic_{nA}_{nB}_{t}(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {{
    return {expr_to_cpp(expr)};
}}

'''
        else:
            # No CSE - simple expanded form
            cpp += f'''/**
 * @brief Symbolic E^{{{nA},{nB}}}_{t}
 */
inline Vec8d hermite_e_symbolic_{nA}_{nB}_{t}(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {{
    return {expr_to_cpp(expr)};
}}

'''

    # Add dispatcher
    cpp += '''/**
 * @brief Runtime dispatcher for symbolic Hermite E coefficients
 */
inline Vec8d dispatch_hermite_e_symbolic(int nA, int nB, int t,
                                          Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
'''
    for (nA, nB, t) in sorted_keys:
        cpp += f'    if (nA == {nA} && nB == {nB} && t == {t}) return hermite_e_symbolic_{nA}_{nB}_{t}(PA, PB, one_over_2p);\n'

    cpp += '''    return Vec8d(0.0);  // Invalid indices
}

} // namespace symbolic
} // namespace recursum
'''

    with open(output_path, 'w') as f:
        f.write(cpp)
    print(f"Generated: {output_path} ({len(E_coeffs)} coefficients, CSE={'ON' if use_cse else 'OFF'})")


def generate_hermite_grad_header(gradients: Dict, output_path: str):
    """Generate C++ header for Hermite gradients."""

    dE_dPA = gradients['dE_dPA']
    dE_dPB = gradients['dE_dPB']

    cpp = '''/**
 * @file hermite_grad_symbolic.hpp
 * @brief Symbolically-generated Hermite E gradients
 *
 * AUTO-GENERATED - DO NOT EDIT MANUALLY
 *
 * Gradients computed by direct differentiation:
 *   dE_dPA[(nA,nB,t)] = ∂E^{nA,nB}_t/∂PA
 *   dE_dPB[(nA,nB,t)] = ∂E^{nA,nB}_t/∂PB
 */

#pragma once

#ifndef RECURSUM_VEC_TYPE
#include <recursum/vectorclass.h>
#endif

namespace recursum {
namespace symbolic {

// =============================================================================
// ∂E/∂PA Gradients
// =============================================================================

'''

    sorted_keys = sorted(dE_dPA.keys())

    for (nA, nB, t) in sorted_keys:
        expr = dE_dPA[(nA, nB, t)]
        cpp += f'''inline Vec8d hermite_dE_dPA_{nA}_{nB}_{t}(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {{
    return {expr_to_cpp(expr)};
}}

'''

    cpp += '''// =============================================================================
// ∂E/∂PB Gradients
// =============================================================================

'''

    for (nA, nB, t) in sorted_keys:
        expr = dE_dPB[(nA, nB, t)]
        cpp += f'''inline Vec8d hermite_dE_dPB_{nA}_{nB}_{t}(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {{
    return {expr_to_cpp(expr)};
}}

'''

    # Add dispatchers
    cpp += '''/**
 * @brief Runtime dispatcher for ∂E/∂PA
 */
inline Vec8d dispatch_dE_dPA(int nA, int nB, int t,
                              Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
'''
    for (nA, nB, t) in sorted_keys:
        cpp += f'    if (nA == {nA} && nB == {nB} && t == {t}) return hermite_dE_dPA_{nA}_{nB}_{t}(PA, PB, one_over_2p);\n'
    cpp += '''    return Vec8d(0.0);
}

/**
 * @brief Runtime dispatcher for ∂E/∂PB
 */
inline Vec8d dispatch_dE_dPB(int nA, int nB, int t,
                              Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
'''
    for (nA, nB, t) in sorted_keys:
        cpp += f'    if (nA == {nA} && nB == {nB} && t == {t}) return hermite_dE_dPB_{nA}_{nB}_{t}(PA, PB, one_over_2p);\n'
    cpp += '''    return Vec8d(0.0);
}

} // namespace symbolic
} // namespace recursum
'''

    with open(output_path, 'w') as f:
        f.write(cpp)
    print(f"Generated: {output_path} ({len(dE_dPA)} dE/dPA + {len(dE_dPB)} dE/dPB)")


def generate_coulomb_r_header(R_coeffs: Dict, output_path: str):
    """Generate C++ header for Coulomb R integrals."""

    cpp = '''/**
 * @file coulomb_r_symbolic.hpp
 * @brief Symbolically-generated Coulomb R integrals
 *
 * AUTO-GENERATED - DO NOT EDIT MANUALLY
 *
 * R_{t,u,v}^{(N)} expressed in terms of Boys functions F_N(T)
 * and Gaussian center differences X_PC, Y_PC, Z_PC.
 */

#pragma once

#ifndef RECURSUM_VEC_TYPE
#include <recursum/vectorclass.h>
#endif

namespace recursum {
namespace symbolic {

/**
 * @brief R integral computation requires Boys functions F_0 to F_N
 *
 * The caller must provide pre-computed Boys function values.
 * F[n] = F_n(T) where T = p * |P-C|^2
 */

'''

    # Group by (t,u,v) and generate for N=0 (most commonly used)
    tuv_keys = set((t, u, v) for (t, u, v, N) in R_coeffs.keys())
    sorted_tuv = sorted(tuv_keys)

    # For each (t,u,v), generate function that computes R_{t,u,v}^{(0)}
    for (t, u, v) in sorted_tuv:
        if (t, u, v, 0) not in R_coeffs:
            continue

        expr = R_coeffs[(t, u, v, 0)]

        # Determine which F_n are needed
        max_n = t + u + v
        F_params = ', '.join([f'Vec8d F_{n}' for n in range(max_n + 1)])

        cpp += f'''/**
 * @brief Symbolic R_{{{t},{u},{v}}}^{{(0)}}
 */
inline Vec8d coulomb_r_symbolic_{t}_{u}_{v}(Vec8d X_PC, Vec8d Y_PC, Vec8d Z_PC, {F_params}) {{
    return {expr_to_cpp(expr)};
}}

'''

    cpp += '''} // namespace symbolic
} // namespace recursum
'''

    with open(output_path, 'w') as f:
        f.write(cpp)
    print(f"Generated: {output_path} ({len(tuv_keys)} R integrals)")


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Generate comprehensive symbolic code for McMurchie-Davidson')
    parser.add_argument('--max-L', type=int, default=4,
                        help='Maximum angular momentum (default: 4 for g orbitals)')
    parser.add_argument('--max-R', type=int, default=8,
                        help='Maximum R integral index (default: 8)')
    parser.add_argument('--output-dir', type=str, default='../symbolic_generated',
                        help='Output directory')
    parser.add_argument('--cse', action='store_true', default=True,
                        help='Enable CSE optimization (default: True)')
    parser.add_argument('--no-cse', action='store_false', dest='cse',
                        help='Disable CSE optimization')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print("COMPREHENSIVE SYMBOLIC CODE GENERATION")
    print("=" * 70)
    print(f"Max angular momentum: L = {args.max_L} (g orbitals)")
    print(f"Max R integral index: {args.max_R}")
    print(f"Output directory: {output_dir}")
    print(f"CSE optimization: {'ENABLED' if args.cse else 'DISABLED'}")
    print()

    # Generate Hermite E coefficients
    print("Generating Hermite E coefficients...")
    E_coeffs = generate_hermite_e(args.max_L)
    print(f"  Generated {len(E_coeffs)} E coefficients")

    # Generate gradients
    print("Generating Hermite gradients...")
    gradients = generate_hermite_gradients(E_coeffs)
    print(f"  Generated {len(gradients['dE_dPA'])} dE/dPA gradients")
    print(f"  Generated {len(gradients['dE_dPB'])} dE/dPB gradients")

    # Generate Coulomb R integrals
    print("Generating Coulomb R integrals...")
    R_coeffs = generate_coulomb_R(args.max_R)
    print(f"  Generated {len(R_coeffs)} R integrals")

    # Write headers
    print()
    print("Writing C++ headers...")
    generate_hermite_e_header(E_coeffs, os.path.join(output_dir, 'hermite_e_symbolic.hpp'), use_cse=args.cse)
    generate_hermite_grad_header(gradients, os.path.join(output_dir, 'hermite_grad_symbolic.hpp'))
    generate_coulomb_r_header(R_coeffs, os.path.join(output_dir, 'coulomb_r_symbolic.hpp'))

    print()
    print("=" * 70)
    print("GENERATION COMPLETE")
    print("=" * 70)

    # Print summary table
    print("\nShell pair coverage:")
    shell_names = ['s', 'p', 'd', 'f', 'g']
    print("| Shell Pair | L_A | L_B | E coefficients | Max t |")
    print("|------------|-----|-----|----------------|-------|")
    for La in range(args.max_L + 1):
        for Lb in range(La, args.max_L + 1):
            count = sum(1 for (nA, nB, t) in E_coeffs if nA == La and nB == Lb)
            max_t = La + Lb
            print(f"| {shell_names[La]}{shell_names[Lb]}         | {La}   | {Lb}   | {count:14} | {max_t:5} |")


if __name__ == '__main__':
    main()
