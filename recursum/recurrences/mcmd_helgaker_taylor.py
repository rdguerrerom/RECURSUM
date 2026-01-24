"""
Helgaker-Taylor derivative recurrences for E coefficients.

These are derivatives of E^{i,j}_t with respect to ATOMIC COORDINATES (A_τ, B_τ),
which include the derivative of the exponential Gaussian prefactor.

This is different from derivatives w.r.t. P_τ which do not include the exponential term.

References:
    - TeraChem SI Equations S4-S7 (F_orbital_revisionSI.pdf)
    - Helgaker et al., "Molecular Electronic-Structure Theory" Ch. 9
"""

from ..codegen import Recurrence


def hermite_e_deriv_A() -> Recurrence:
    """
    Hermite coefficient derivative ∂E^{i,j}_t / ∂A_τ (Helgaker-Taylor form).

    This is the FULL derivative including the exponential Gaussian factor:

    ∂/∂A_τ [E^{i,j}_t × exp(-μ(A-B)²)]

    where the exponential derivative contributes an additional term.

    Recurrence (TeraChem SI Eq S4):
    E^{i+1,j,1}_{t,τ} = (1/2p)E^{i,j,1}_{t-1,τ} - [b/(a+b)]E^{i,j}_{t,τ}
                         + (P_τ - A_τ)E^{i,j,1}_{t,τ} + (t+1)E^{i,j,1}_{t+1,τ}

    Base case (TeraChem SI Eq S6):
    E^{0,0,1}_{0,τ} = 2a(P_τ - A_τ) × exp(-ab/(a+b)(A_τ-B_τ)²)

    Note: The "1" superscript denotes derivative order (first derivative).
    The E^{i,j}_{t,τ} term (without derivative) comes from ∂/∂A_τ of the exponential.

    Parameters:
        i, j: Angular momentum quantum numbers for centers A, B
        t: Hermite index (0 ≤ t ≤ i+j for derivatives)
        p: a + b (sum of exponents)
        PA: P_τ - A_τ (P is Gaussian product center)
        PB: P_τ - B_τ
        aAB: a / (a+b) (reduced exponent ratio)
        bAB: b / (a+b)
    """
    rec = Recurrence(
        name="HermiteE_dA",
        indices=["i", "j", "t"],
        runtime_vars=["p", "PA", "PB", "aAB", "bAB"],
        namespace="mcmd_helgaker_taylor",
        max_indices={"i": 3, "j": 3, "t": 7}  # Max L=3 (f-orbitals), t can go up to i+j+1 for derivatives
    )

    rec.validity("i >= 0", "j >= 0", "t >= 0", "t <= i + j + 1")

    # Base case: E^{0,0,1}_0 = 2a(P - A) = 2a × PA
    # The exponential factor is absorbed into the coefficient
    # Note: We use PA directly since P - A = PA in our convention
    rec.base(i=0, j=0, t=0, value="2.0 * aAB * p * PA")

    # Boundary conditions: E^{i,j,1}_t = 0 for t < 0 or t > i+j
    rec.base("t < 0", value=0.0)
    rec.base("t > i + j", value=0.0)

    # Increment i recurrence (TeraChem SI Eq S4)
    # E^{i+1,j,1}_t = (1/2p)E^{i,j,1}_{t-1} - (b/p)E^{i,j}_t
    #                  + PA×E^{i,j,1}_t + (t+1)×E^{i,j,1}_{t+1}
    rec.rule(
        "i > 0",
        "0.5 / p * E_dA[i-1, j, t-1] - bAB * E[i-1, j, t] + "
        "PA * E_dA[i-1, j, t] + (t + 1) * E_dA[i-1, j, t+1]",
        name="Increment i (Eq S4)",
        requires=["E"]  # Requires standard E coefficient
    )

    # Increment j recurrence (TeraChem SI Eq S5)
    # E^{i,j+1,1}_t = (1/2p)E^{i,j,1}_{t-1} + (a/p)E^{i,j}_t
    #                  + PB×E^{i,j,1}_t + (t+1)×E^{i,j,1}_{t+1}
    rec.rule(
        "j > 0",
        "0.5 / p * E_dA[i, j-1, t-1] + aAB * E[i, j-1, t] + "
        "PB * E_dA[i, j-1, t] + (t + 1) * E_dA[i, j-1, t+1]",
        name="Increment j (Eq S5)",
        requires=["E"]
    )

    return rec


def hermite_e_deriv_B() -> Recurrence:
    """
    Hermite coefficient derivative ∂E^{i,j}_t / ∂B_τ (Helgaker-Taylor form).

    Similar to ∂E/∂A but with opposite sign from the exponential derivative.

    Recurrence (TeraChem SI Eq S5 adapted for B):
    The derivative w.r.t. B has opposite sign from A in the exponential term.

    Base case:
    E^{0,0,1}_{0,τ} w.r.t. B = -2b(P_τ - B_τ) × exp(...)
                              = -2b × PB

    Note: Since ∂/∂B_τ [exp(-μ(A-B)²)] = -∂/∂A_τ [exp(...)], the B derivative
    has opposite sign in the E^{i,j}_t term.
    """
    rec = Recurrence(
        name="HermiteE_dB",
        indices=["i", "j", "t"],
        runtime_vars=["p", "PA", "PB", "aAB", "bAB"],
        namespace="mcmd_helgaker_taylor",
        max_indices={"i": 3, "j": 3, "t": 7}  # Max L=3 (f-orbitals)
    )

    rec.validity("i >= 0", "j >= 0", "t >= 0", "t <= i + j + 1")

    # Base case: E^{0,0,1}_0 w.r.t. B = -2b × PB
    rec.base(i=0, j=0, t=0, value="-2.0 * bAB * p * PB")

    # Boundary conditions
    rec.base("t < 0", value=0.0)
    rec.base("t > i + j", value=0.0)

    # Increment i: same as dA but with opposite sign on E term
    rec.rule(
        "i > 0",
        "0.5 / p * E_dB[i-1, j, t-1] + bAB * E[i-1, j, t] + "
        "PA * E_dB[i-1, j, t] + (t + 1) * E_dB[i-1, j, t+1]",
        name="Increment i for dB",
        requires=["E"]
    )

    # Increment j: opposite sign on E term
    rec.rule(
        "j > 0",
        "0.5 / p * E_dB[i, j-1, t-1] - aAB * E[i, j-1, t] + "
        "PB * E_dB[i, j-1, t] + (t + 1) * E_dB[i, j-1, t+1]",
        name="Increment j for dB",
        requires=["E"]
    )

    return rec
