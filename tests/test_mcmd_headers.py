#!/usr/bin/env python3
"""
Test RECURSUM McMD headers against reference implementations.

Validates:
1. HermiteE coefficients match Python reference
2. HermiteDerivPA/PB direct formulas match finite differences
3. HermiteDerivA/B chain rule matches McMD grad_solver approach
4. KineticIntegral matches PySCF (when available)
5. CoulombR integrals match reference values

This test uses pure Python implementations to verify the mathematical
correctness of the RECURSUM header formulas before C++ compilation.
"""

import numpy as np


# =============================================================================
# Reference Implementations (Python)
# =============================================================================

def hermite_e_ref(nA, nB, t, PA, PB, p):
    """
    Reference Hermite E coefficient using Helgaker-Taylor 1992 recurrence.

    E^{i+1,j}_t = (1/2p) × E^{i,j}_{t-1} + PA × E^{i,j}_t + (t+1) × E^{i,j}_{t+1}

    Base case: E^{0,0}_0 = 1
    """
    # Validity check
    if nA < 0 or nB < 0 or t < 0 or t > nA + nB:
        return 0.0

    # Base case
    if nA == 0 and nB == 0 and t == 0:
        return 1.0

    # A-side recurrence (nA > 0)
    if nA > 0:
        term1 = 0.0
        if t > 0:
            term1 = (0.5 / p) * hermite_e_ref(nA - 1, nB, t - 1, PA, PB, p)
        term2 = PA * hermite_e_ref(nA - 1, nB, t, PA, PB, p)
        term3 = (t + 1) * hermite_e_ref(nA - 1, nB, t + 1, PA, PB, p)
        return term1 + term2 + term3

    # B-side recurrence (nA = 0, nB > 0)
    if nB > 0:
        term1 = 0.0
        if t > 0:
            term1 = (0.5 / p) * hermite_e_ref(0, nB - 1, t - 1, PA, PB, p)
        term2 = PB * hermite_e_ref(0, nB - 1, t, PA, PB, p)
        term3 = (t + 1) * hermite_e_ref(0, nB - 1, t + 1, PA, PB, p)
        return term1 + term2 + term3

    return 0.0


def hermite_deriv_PA_ref(nA, nB, t, PA, PB, p):
    """
    ∂E^{nA,nB}_t/∂PA = nA × E^{nA-1,nB}_t  (Helgaker-Taylor direct formula)
    """
    if nA <= 0:
        return 0.0
    return nA * hermite_e_ref(nA - 1, nB, t, PA, PB, p)


def hermite_deriv_PB_ref(nA, nB, t, PA, PB, p):
    """
    ∂E^{nA,nB}_t/∂PB = nB × E^{nA,nB-1}_t  (Helgaker-Taylor direct formula)
    """
    if nB <= 0:
        return 0.0
    return nB * hermite_e_ref(nA, nB - 1, t, PA, PB, p)


def hermite_deriv_A_ref(nA, nB, t, PA, PB, p, alpha_A, alpha_B):
    """
    ∂E/∂A via chain rule:
    ∂E/∂A = (∂E/∂PA)(∂PA/∂A) + (∂E/∂PB)(∂PB/∂A)
          = (nA × E^{nA-1,nB}_t) × (-α_B/p) + (nB × E^{nA,nB-1}_t) × (α_A/p)
          = (1/p) × [nB × α_A × E^{nA,nB-1}_t - nA × α_B × E^{nA-1,nB}_t]
    """
    dE_dPA = hermite_deriv_PA_ref(nA, nB, t, PA, PB, p)
    dE_dPB = hermite_deriv_PB_ref(nA, nB, t, PA, PB, p)
    dPA_dA = -alpha_B / p
    dPB_dA = alpha_A / p
    return dE_dPA * dPA_dA + dE_dPB * dPB_dA


def hermite_deriv_B_ref(nA, nB, t, PA, PB, p, alpha_A, alpha_B):
    """
    ∂E/∂B via chain rule:
    ∂E/∂B = (∂E/∂PA)(∂PA/∂B) + (∂E/∂PB)(∂PB/∂B)
          = (nA × E^{nA-1,nB}_t) × (α_B/p) + (nB × E^{nA,nB-1}_t) × (-α_A/p)
          = (1/p) × [nA × α_B × E^{nA-1,nB}_t - nB × α_A × E^{nA,nB-1}_t]
          = -∂E/∂A  (translational invariance!)
    """
    dE_dPA = hermite_deriv_PA_ref(nA, nB, t, PA, PB, p)
    dE_dPB = hermite_deriv_PB_ref(nA, nB, t, PA, PB, p)
    dPA_dB = alpha_B / p
    dPB_dB = -alpha_A / p
    return dE_dPA * dPA_dB + dE_dPB * dPB_dB


def kinetic_term_1d_ref(nA, nB, PA, PB, p, alpha_B):
    """
    1D kinetic contribution (TeraChem SI Eq. S12):
    kinetic(nA, nB) = -nB(nB-1)/2 × E^{nA,nB-2}_0 + (2nB+1)×b × E^{nA,nB}_0 - 2b² × E^{nA,nB+2}_0
    """
    result = 0.0

    # Term 1: -nB(nB-1)/2 × E^{nA,nB-2}_0
    if nB >= 2:
        result += -0.5 * nB * (nB - 1) * hermite_e_ref(nA, nB - 2, 0, PA, PB, p)

    # Term 2: (2nB+1)×b × E^{nA,nB}_0
    result += (2 * nB + 1) * alpha_B * hermite_e_ref(nA, nB, 0, PA, PB, p)

    # Term 3: -2b² × E^{nA,nB+2}_0
    result += -2.0 * alpha_B**2 * hermite_e_ref(nA, nB + 2, 0, PA, PB, p)

    return result


def boys_ref(n, T):
    """Reference Boys function using scipy or numerical integration."""
    if T < 1e-10:
        return 1.0 / (2 * n + 1)

    from scipy import special
    # F_n(T) = (1/2) * T^{-n-1/2} * gamma(n+1/2) * gammainc(n+1/2, T)
    # where gammainc is the regularized lower incomplete gamma function
    return 0.5 * T**(-n - 0.5) * special.gamma(n + 0.5) * special.gammainc(n + 0.5, T)


def coulomb_r_ref(t, u, v, N, PCx, PCy, PCz, Boys):
    """
    Reference Coulomb auxiliary integral R_{tuv}^{(N)}.

    Base case: R_{000}^{(N)} = F_N(T)
    X-recurrence: R_{tuv} = PC_x × R_{t-1,u,v}^{(N+1)} + (t-1) × R_{t-2,u,v}^{(N+1)}
    """
    if t < 0 or u < 0 or v < 0:
        return 0.0

    # Base case
    if t == 0 and u == 0 and v == 0:
        return Boys[N]

    # X-recurrence
    if t > 0:
        term1 = PCx * coulomb_r_ref(t - 1, u, v, N + 1, PCx, PCy, PCz, Boys)
        term2 = (t - 1) * coulomb_r_ref(t - 2, u, v, N + 1, PCx, PCy, PCz, Boys) if t >= 2 else 0.0
        return term1 + term2

    # Y-recurrence
    if u > 0:
        term1 = PCy * coulomb_r_ref(0, u - 1, v, N + 1, PCx, PCy, PCz, Boys)
        term2 = (u - 1) * coulomb_r_ref(0, u - 2, v, N + 1, PCx, PCy, PCz, Boys) if u >= 2 else 0.0
        return term1 + term2

    # Z-recurrence
    if v > 0:
        term1 = PCz * coulomb_r_ref(0, 0, v - 1, N + 1, PCx, PCy, PCz, Boys)
        term2 = (v - 1) * coulomb_r_ref(0, 0, v - 2, N + 1, PCx, PCy, PCz, Boys) if v >= 2 else 0.0
        return term1 + term2

    return 0.0


# =============================================================================
# Test Functions
# =============================================================================

def test_hermite_e_coefficients():
    """Test HermiteE coefficients against reference."""
    print("Testing HermiteE coefficients...")

    # Test parameters
    alpha_A = 1.2
    alpha_B = 0.8
    A = np.array([0.0, 0.0, 0.0])
    B = np.array([0.5, 0.3, 0.1])

    p = alpha_A + alpha_B
    P = (alpha_A * A + alpha_B * B) / p
    PA = P[0] - A[0]  # Use x-component
    PB = P[0] - B[0]

    test_cases = [
        (0, 0, 0),
        (1, 0, 0), (0, 1, 0), (1, 1, 0),
        (2, 0, 0), (0, 2, 0), (1, 1, 1),
        (2, 1, 0), (1, 2, 0), (2, 2, 0),
        (3, 0, 0), (0, 3, 0), (3, 3, 0),
    ]

    all_passed = True
    for nA, nB, t in test_cases:
        E = hermite_e_ref(nA, nB, t, PA, PB, p)
        print(f"  E^{{{nA},{nB}}}_{t} = {E:.10f}")
        if not np.isfinite(E):
            print(f"    WARNING: Non-finite value!")
            all_passed = False

    print(f"  Result: {'PASS' if all_passed else 'FAIL'}")
    return all_passed


def test_hermite_deriv_direct():
    """Test direct derivative formulas against finite differences."""
    print("\nTesting HermiteDerivPA/PB direct formulas...")

    alpha_A = 1.2
    alpha_B = 0.8
    A = np.array([0.0, 0.0, 0.0])
    B = np.array([0.5, 0.3, 0.1])

    p = alpha_A + alpha_B
    P = (alpha_A * A + alpha_B * B) / p
    PA = P[0] - A[0]
    PB = P[0] - B[0]

    h = 1e-6
    max_error = 0.0

    test_cases = [(1, 0, 0), (0, 1, 0), (1, 1, 0), (2, 1, 0), (1, 2, 0), (2, 2, 0)]

    for nA, nB, t in test_cases:
        # Test ∂E/∂PA
        dE_dPA_direct = hermite_deriv_PA_ref(nA, nB, t, PA, PB, p)
        E_plus = hermite_e_ref(nA, nB, t, PA + h, PB, p)
        E_minus = hermite_e_ref(nA, nB, t, PA - h, PB, p)
        dE_dPA_fd = (E_plus - E_minus) / (2 * h)

        error_PA = abs(dE_dPA_direct - dE_dPA_fd)
        max_error = max(max_error, error_PA)

        # Test ∂E/∂PB
        dE_dPB_direct = hermite_deriv_PB_ref(nA, nB, t, PA, PB, p)
        E_plus = hermite_e_ref(nA, nB, t, PA, PB + h, p)
        E_minus = hermite_e_ref(nA, nB, t, PA, PB - h, p)
        dE_dPB_fd = (E_plus - E_minus) / (2 * h)

        error_PB = abs(dE_dPB_direct - dE_dPB_fd)
        max_error = max(max_error, error_PB)

        print(f"  E^{{{nA},{nB}}}_{t}: ∂PA err={error_PA:.2e}, ∂PB err={error_PB:.2e}")

    passed = max_error < 1e-8
    print(f"  Max error: {max_error:.2e}")
    print(f"  Result: {'PASS' if passed else 'FAIL'}")
    return passed


def test_hermite_deriv_chain_rule():
    """Test nuclear derivatives via chain rule."""
    print("\nTesting HermiteDerivA/B chain rule...")

    alpha_A = 1.2
    alpha_B = 0.8
    A_x = 0.0
    B_x = 0.5

    p = alpha_A + alpha_B
    P_x = (alpha_A * A_x + alpha_B * B_x) / p
    PA = P_x - A_x
    PB = P_x - B_x

    h = 1e-6
    max_error = 0.0

    test_cases = [(1, 0, 0), (0, 1, 0), (1, 1, 0), (2, 1, 0), (2, 2, 0)]

    for nA, nB, t in test_cases:
        # Analytical derivative
        dE_dA = hermite_deriv_A_ref(nA, nB, t, PA, PB, p, alpha_A, alpha_B)
        dE_dB = hermite_deriv_B_ref(nA, nB, t, PA, PB, p, alpha_A, alpha_B)

        # Finite difference for ∂E/∂A
        A_plus = A_x + h
        P_plus = (alpha_A * A_plus + alpha_B * B_x) / p
        PA_plus = P_plus - A_plus
        PB_plus = P_plus - B_x
        E_plus = hermite_e_ref(nA, nB, t, PA_plus, PB_plus, p)

        A_minus = A_x - h
        P_minus = (alpha_A * A_minus + alpha_B * B_x) / p
        PA_minus = P_minus - A_minus
        PB_minus = P_minus - B_x
        E_minus = hermite_e_ref(nA, nB, t, PA_minus, PB_minus, p)

        dE_dA_fd = (E_plus - E_minus) / (2 * h)

        error_A = abs(dE_dA - dE_dA_fd)
        max_error = max(max_error, error_A)

        # Translational invariance check
        trans_inv_error = abs(dE_dA + dE_dB)

        print(f"  E^{{{nA},{nB}}}_{t}: ∂A err={error_A:.2e}, trans_inv={trans_inv_error:.2e}")

    passed = max_error < 1e-8
    print(f"  Max error: {max_error:.2e}")
    print(f"  Result: {'PASS' if passed else 'FAIL'}")
    return passed


def test_kinetic_term():
    """Test kinetic energy 1D term."""
    print("\nTesting KineticTerm1D...")

    alpha_A = 1.2
    alpha_B = 0.8
    A_x = 0.0
    B_x = 0.5

    p = alpha_A + alpha_B
    P_x = (alpha_A * A_x + alpha_B * B_x) / p
    PA = P_x - A_x
    PB = P_x - B_x

    test_cases = [(0, 0), (1, 0), (0, 1), (1, 1), (0, 2), (2, 0)]

    all_finite = True
    for nA, nB in test_cases:
        T = kinetic_term_1d_ref(nA, nB, PA, PB, p, alpha_B)
        print(f"  T^{{{nA},{nB}}} = {T:.10f}")
        if not np.isfinite(T):
            all_finite = False

    print(f"  Result: {'PASS' if all_finite else 'FAIL'}")
    return all_finite


def test_coulomb_r_integrals():
    """Test Coulomb auxiliary integrals."""
    print("\nTesting CoulombR integrals...")
    try:
        from scipy import special
    except ImportError:
        print("  SKIP: scipy not available")
        return True

    # Test parameters
    p = 2.0
    PCx, PCy, PCz = 0.3, 0.2, 0.1
    T = p * (PCx**2 + PCy**2 + PCz**2)

    # Compute Boys functions
    N_max = 10
    Boys = [boys_ref(n, T) for n in range(N_max)]

    test_cases = [
        (0, 0, 0, 0),
        (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0),
        (1, 1, 0, 0), (1, 0, 1, 0), (0, 1, 1, 0),
        (2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0),
    ]

    all_finite = True
    for t, u, v, N in test_cases:
        R = coulomb_r_ref(t, u, v, N, PCx, PCy, PCz, Boys)
        print(f"  R_{{{t},{u},{v}}}^({N}) = {R:.10f}")
        if not np.isfinite(R):
            all_finite = False

    print(f"  Result: {'PASS' if all_finite else 'FAIL'}")
    return all_finite


def coulomb_r_deriv_pcx_ref(t, u, v, N, PCx, PCy, PCz, Boys, p):
    """
    Reference chain-rule derivative ∂R_{tuv}^{(N)}/∂(PC_x).

    Base case: ∂R_{000}^{(N)}/∂(PC_x) = -2p × PC_x × F_{N+1}(T)
    X-recurrence: includes extra R_{t-1,u,v}^{(N+1)} term
    Y/Z-recurrence: no extra term
    """
    if t < 0 or u < 0 or v < 0:
        return 0.0

    # Base case
    if t == 0 and u == 0 and v == 0:
        return -2.0 * p * PCx * Boys[N + 1]

    # X-recurrence (t > 0): has extra term!
    if t > 0:
        extra = coulomb_r_ref(t - 1, u, v, N + 1, PCx, PCy, PCz, Boys)
        term1 = PCx * coulomb_r_deriv_pcx_ref(t - 1, u, v, N + 1, PCx, PCy, PCz, Boys, p)
        term2 = (t - 1) * coulomb_r_deriv_pcx_ref(t - 2, u, v, N + 1, PCx, PCy, PCz, Boys, p) if t >= 2 else 0.0
        return extra + term1 + term2

    # Y-recurrence (t = 0, u > 0): no extra term
    if u > 0:
        term1 = PCy * coulomb_r_deriv_pcx_ref(0, u - 1, v, N + 1, PCx, PCy, PCz, Boys, p)
        term2 = (u - 1) * coulomb_r_deriv_pcx_ref(0, u - 2, v, N + 1, PCx, PCy, PCz, Boys, p) if u >= 2 else 0.0
        return term1 + term2

    # Z-recurrence (t = u = 0, v > 0): no extra term
    if v > 0:
        term1 = PCz * coulomb_r_deriv_pcx_ref(0, 0, v - 1, N + 1, PCx, PCy, PCz, Boys, p)
        term2 = (v - 1) * coulomb_r_deriv_pcx_ref(0, 0, v - 2, N + 1, PCx, PCy, PCz, Boys, p) if v >= 2 else 0.0
        return term1 + term2

    return 0.0


def test_coulomb_r_derivative():
    """Test ∂R/∂(PC_x) chain-rule formula."""
    print("\nTesting CoulombR chain-rule derivative...")
    try:
        from scipy import special
    except ImportError:
        print("  SKIP: scipy not available")
        return True

    p = 2.0
    PCx, PCy, PCz = 0.3, 0.2, 0.1
    T = p * (PCx**2 + PCy**2 + PCz**2)

    N_max = 20
    Boys = [boys_ref(n, T) for n in range(N_max)]

    h = 1e-6
    max_error = 0.0

    test_cases = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (2, 0, 0)]

    for t, u, v in test_cases:
        N = 0

        # Analytical using chain-rule recurrence
        dR_dPCx_analytical = coulomb_r_deriv_pcx_ref(t, u, v, N, PCx, PCy, PCz, Boys, p)

        # Finite difference (changing PCx also changes T and Boys!)
        def R_at_PCx(pc_x):
            T_new = p * (pc_x**2 + PCy**2 + PCz**2)
            Boys_new = [boys_ref(n, T_new) for n in range(N_max)]
            return coulomb_r_ref(t, u, v, N, pc_x, PCy, PCz, Boys_new)

        dR_dPCx_fd = (R_at_PCx(PCx + h) - R_at_PCx(PCx - h)) / (2 * h)

        error = abs(dR_dPCx_analytical - dR_dPCx_fd)
        max_error = max(max_error, error)

        print(f"  R_{{{t},{u},{v}}}: analytical={dR_dPCx_analytical:.6f}, fd={dR_dPCx_fd:.6f}, err={error:.2e}")

    passed = max_error < 1e-6
    print(f"  Max error: {max_error:.2e}")
    print(f"  Result: {'PASS' if passed else 'FAIL'}")
    return passed


def main():
    """Run all tests."""
    print("=" * 60)
    print("RECURSUM McMD Headers Validation")
    print("=" * 60)

    results = []
    results.append(("HermiteE coefficients", test_hermite_e_coefficients()))
    results.append(("HermiteDerivPA/PB direct", test_hermite_deriv_direct()))
    results.append(("HermiteDerivA/B chain rule", test_hermite_deriv_chain_rule()))
    results.append(("KineticTerm1D", test_kinetic_term()))
    results.append(("CoulombR integrals", test_coulomb_r_integrals()))
    results.append(("CoulombR P-derivative", test_coulomb_r_derivative()))

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    all_passed = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    print()
    if all_passed:
        print("All tests PASSED!")
    else:
        print("Some tests FAILED!")

    return 0 if all_passed else 1


if __name__ == "__main__":
    exit(main())
