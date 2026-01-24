"""
McMD Gradient Recurrence Relations

E coefficient derivatives for McMurchie-Davidson gradient integrals.
Based on equations S4-S7 from F_orbital_revisionSI.pdf (TeraChem SI).

These recurrences implement the derivative of Gaussian overlap expansion
coefficients with respect to the separation coordinate Δ = A - B.
"""

from ..codegen.recurrence import Recurrence


def e_coefficient_derivative() -> Recurrence:
    """
    E^{i,j,1} derivative coefficients (Equations S4-S7).

    These are the derivatives of the McMurchie-Davidson E coefficients
    with respect to the separation coordinate Δτ = Aτ - Bτ.

    **Recurrence Relations:**

    Base case (Eq S6):
        E^{0,0,1}_{0,τ} = 2a(Pτ - Aτ)exp(-ab/(a+b)(Aτ-Bτ)²)

    Increment i (Eq S4):
        E^{i+1,j,1}_t = (1/2p)E^{i,j,1}_{t-1} - (b/(a+b))E^{i,j}_t
                        + (Pτ - Aτ)E^{i,j,1}_t + (t+1)E^{i,j,1}_{t+1}

    Increment j (Eq S5):
        E^{i,j+1,1}_t = (1/2p)E^{i,j,1}_{t-1} + (a/(a+b))E^{i,j}_t
                        + (Pτ - Bτ)E^{i,j,1}_t + (t+1)E^{i,j,1}_{t+1}

    Boundary (Eq S7):
        E^{i,j,1}_t = 0  if t < 0 or t > i + j

    **Parameters:**
    - i, j: Angular momentum indices (0 to 5 for f orbitals)
    - t: Hermite polynomial index (0 to i+j+1 for derivatives)
    - a, b: Gaussian exponents
    - p: Combined exponent = a + b
    - P_tau: Center-of-mass coordinate component
    - A_tau, B_tau: GTO center coordinate components
    - Delta_sq: (Aτ - Bτ)² for exponential term

    **Note on normalization:**
    The non-derivative E coefficients E^{i,j}_t are referenced in the
    recurrence. These must be computed separately with correct Cartesian
    normalization via GaussianNormalization::cartesianNorm3D(). RECURSUM
    handles only the recurrence logic; normalization is applied externally.
    """
    rec = Recurrence(
        name="E_deriv",
        indices=["i", "j", "t"],
        runtime_vars=["a", "b", "p", "P_tau", "A_tau", "B_tau", "Delta_sq"],
        namespace="mcmd",
        max_indices={"i": 3, "j": 3, "t": 7},  # Max L=3 (f-orbitals), t up to 6+1 for derivatives
        scipy_reference=None  # No SciPy equivalent for McMD E coefficients
    )

    # Base case (Eq S6): E^{0,0,1}_{0,τ} = 2a(Pτ - Aτ)exp(-ab/(a+b)(Aτ-Bτ)²)
    # Note: This can also be written as -2b(Pτ - Bτ)exp(...) (equivalent)
    rec.base(
        i=0, j=0, t=0,
        value="2.0 * a * (P_tau - A_tau) * exp(-(a * b) / (a + b) * Delta_sq)"
    )

    # Boundary conditions (Eq S7)
    # E^{i,j,1}_t = 0 if t < 0 or t > i + j
    # This is handled implicitly by RECURSUM - we simply don't define
    # out-of-bounds cases, and they default to zero.

    # Recurrence for incrementing i (Eq S4)
    # E^{i+1,j,1}_t = (1/2p)E^{i,j,1}_{t-1} - (b/(a+b))E^{i,j}_t
    #                 + (Pτ-Aτ)E^{i,j,1}_t + (t+1)E^{i,j,1}_{t+1}
    #
    # CRITICAL: The term E^{i,j}_t is the NON-derivative coefficient.
    # This must be provided as a separate recurrence or computed externally.
    rec.rule(
        "i > 0 && t >= 0 && t <= i + j",
        """(0.5 / p) * E_deriv[i-1, j, t-1]
            - (b / (a + b)) * E[i-1, j, t]
            + (P_tau - A_tau) * E_deriv[i-1, j, t]
            + (t + 1.0) * E_deriv[i-1, j, t+1]""",
        name="Increment i (Eq S4)"
    )

    # Recurrence for incrementing j (Eq S5)
    # E^{i,j+1,1}_t = (1/2p)E^{i,j,1}_{t-1} + (a/(a+b))E^{i,j}_t
    #                 + (Pτ-Bτ)E^{i,j,1}_t + (t+1)E^{i,j,1}_{t+1}
    rec.rule(
        "j > 0 && t >= 0 && t <= i + j",
        """(0.5 / p) * E_deriv[i, j-1, t-1]
            + (a / (a + b)) * E[i, j-1, t]
            + (P_tau - B_tau) * E_deriv[i, j-1, t]
            + (t + 1.0) * E_deriv[i, j-1, t+1]""",
        name="Increment j (Eq S5)"
    )

    return rec


def e_coefficient() -> Recurrence:
    """
    E^{i,j} non-derivative coefficients (standard McMurchie-Davidson).

    These are the standard Gaussian overlap expansion coefficients
    used in the McMurchie-Davidson algorithm. Required for computing
    E^{i,j,1} derivatives.

    **Recurrence Relations:**

    Base cases:
        E^{0,0}_0 = exp(-ab/(a+b)(A-B)²)
        E^{0,0}_t = 0 for t > 0

    Horizontal recursion:
        E^{i+1,j}_t = (1/2p)E^{i,j}_{t-1} + (Px - Ax)E^{i,j}_t + (t+1)E^{i,j}_{t+1}
        E^{i,j+1}_t = (1/2p)E^{i,j}_{t-1} + (Px - Bx)E^{i,j}_t + (t+1)E^{i,j}_{t+1}

    **Note:**
    This is a standard recurrence. If McMD already has E coefficients
    computed, this may not need to be implemented in RECURSUM. It's
    included here for completeness and to provide E[i,j,t] terms
    needed by E_deriv.
    """
    rec = Recurrence(
        name="E",
        indices=["i", "j", "t"],
        runtime_vars=["a", "b", "p", "P_tau", "A_tau", "B_tau", "Delta_sq"],
        namespace="mcmd",
        max_indices={"i": 3, "j": 3, "t": 6},  # Max L=3 (f-orbitals)
        scipy_reference=None
    )

    # Base case: E^{0,0}_0 = exp(-ab/(a+b)(A-B)²)
    rec.base(
        i=0, j=0, t=0,
        value="exp(-(a * b) / (a + b) * Delta_sq)"
    )

    # E^{0,0}_t = 0 for t > 0 (only t=0 survives for base case)
    # Handled implicitly by not defining other base cases

    # Increment i: E^{i+1,j}_t = (1/2p)E^{i,j}_{t-1} + (Px-Ax)E^{i,j}_t + (t+1)E^{i,j}_{t+1}
    rec.rule(
        "i > 0 && t >= 0 && t <= i + j",
        """(0.5 / p) * E[i-1, j, t-1]
            + (P_tau - A_tau) * E[i-1, j, t]
            + (t + 1.0) * E[i-1, j, t+1]""",
        name="Increment i"
    )

    # Increment j: E^{i,j+1}_t = (1/2p)E^{i,j}_{t-1} + (Px-Bx)E^{i,j}_t + (t+1)E^{i,j}_{t+1}
    rec.rule(
        "j > 0 && t >= 0 && t <= i + j",
        """(0.5 / p) * E[i, j-1, t-1]
            + (P_tau - B_tau) * E[i, j-1, t]
            + (t + 1.0) * E[i, j-1, t+1]""",
        name="Increment j"
    )

    return rec
