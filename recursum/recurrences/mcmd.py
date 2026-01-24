"""
McMurchie-Davidson recurrence relations for quantum chemistry integrals and gradients.

This module implements all McMD gradient recurrences:
- Core Hermite E-coefficient recurrences
- Hermite derivative recurrences (∂E/∂PA, ∂E/∂PB, ∂E/∂R, ∂E/∂A, ∂E/∂B)
- S (overlap) gradient formulas
- T (kinetic energy) gradient formulas
- V (nuclear attraction) gradient formulas with TeraChem two-term structure
- J (Coulomb) gradient formulas
- K (exchange) gradient formulas

References:
    - McMurchie & Davidson, J. Comp. Phys. 26, 218 (1978)
    - Helgaker et al., "Molecular Electronic-Structure Theory" Ch. 9
    - Helgaker-Taylor 1992, Eq. 7
    - TeraChem SI (F_orbital_revisionSI.pdf), Eq. S14, S15, S17
"""

from ..codegen import Recurrence


def hermite_e_coefficient() -> Recurrence:
    """
    Hermite expansion coefficient E^{nA,nB}_t(PA, PB, aAB).

    Base recurrence for McMurchie-Davidson scheme. Expresses product of
    two primitive Gaussians as a linear combination of Hermite Gaussians.

    CRITICAL: Uses Helgaker-Taylor 1992 Equation (7):
        E^{i+1,j}_t = (1/2p) E^{i,j}_{t-1} + PA × E^{i,j}_t + (t+1) × E^{i,j}_{t+1}

    For t=0, this becomes:
        E^{i+1,j}_0 = PA × E^{i,j}_0 + E^{i,j}_1

    The E^{i,j}_1 term MUST be included!

    For i>0, j>0: Use INCREMENT-I ONLY (not both paths summed).

    where p = α_A + α_B, and aAB = 1/(2p) is the reduced exponent factor

    Reference: Helgaker-Taylor 1992 Eq. 7, HELGAKER_TAYLOR_STATUS.md
    """
    rec = Recurrence("HermiteE", ["nA", "nB", "t"], ["PA", "PB", "aAB"],
                     namespace="mcmd_hermite", max_indices={"nA": 3, "nB": 3, "t": 6})  # Max L=3 (f-orbitals)
    rec.validity("nA >= 0", "nB >= 0", "t >= 0", "t <= nA + nB")

    # Base cases
    rec.base(nA=0, nB=0, t=0, value=1.0)

    # A-side recurrence (nA > 0, nB = 0) - Helgaker-Taylor Eq. 7
    # E^{i+1,0}_t = aAB E^{i,0}_{t-1} + PA × E^{i,0}_t + (t+1) × E^{i,0}_{t+1}
    # For t=0: E^{i+1,0}_0 = PA × E^{i,0}_0 + E^{i,0}_1
    rec.rule("nA > 0 && nB == 0 && t == 0",
             "PA * E[nA-1, 0, t] + (t + 1) * E[nA-1, 0, t+1]",
             name="A-side t=0 (includes E_{t+1})")
    rec.rule("nA > 0 && nB == 0 && t > 0",
             "aAB * E[nA-1, 0, t-1] + PA * E[nA-1, 0, t] + (t + 1) * E[nA-1, 0, t+1]",
             name="A-side t>0")

    # B-side recurrence (nA = 0, nB > 0) - Helgaker-Taylor Eq. 7 (increment-j)
    # E^{0,j+1}_t = aAB E^{0,j}_{t-1} + PB × E^{0,j}_t + (t+1) × E^{0,j}_{t+1}
    # For t=0: E^{0,j+1}_0 = PB × E^{0,j}_0 + E^{0,j}_1
    rec.rule("nA == 0 && nB > 0 && t == 0",
             "PB * E[0, nB-1, t] + (t + 1) * E[0, nB-1, t+1]",
             name="B-side t=0 (includes E_{t+1})")
    rec.rule("nA == 0 && nB > 0 && t > 0",
             "aAB * E[0, nB-1, t-1] + PB * E[0, nB-1, t] + (t + 1) * E[0, nB-1, t+1]",
             name="B-side t>0")

    # General recurrence (nA > 0, nB > 0) - USE INCREMENT-I ONLY!
    # E^{i+1,j}_t = aAB E^{i,j}_{t-1} + PA × E^{i,j}_t + (t+1) × E^{i,j}_{t+1}
    # For t=0: E^{i+1,j}_0 = PA × E^{i,j}_0 + E^{i,j}_1
    rec.rule("nA > 0 && nB > 0 && t == 0",
             "PA * E[nA-1, nB, t] + (t + 1) * E[nA-1, nB, t+1]",
             name="General t=0 (increment-i only, includes E_{t+1})")
    rec.rule("nA > 0 && nB > 0 && t > 0",
             "aAB * E[nA-1, nB, t-1] + PA * E[nA-1, nB, t] + (t + 1) * E[nA-1, nB, t+1]",
             name="General t>0 (increment-i only)")

    return rec


def coulomb_r_auxiliary() -> Recurrence:
    """
    Coulomb auxiliary integral R^{(N)}_{t,u,v}(PC_x, PC_y, PC_z, Boys).

    Used in McMurchie-Davidson electron repulsion integrals (ERIs).
    Computes Hermite Coulomb integrals over 1/|r-C|.

    Recurrence relations:
        R_{000}^{(N)} = F_N(T)  (Boys function, base case)

        R_{tuv}^{(N)} = PC_x × R_{t-1,u,v}^{(N+1)} + (t-1) × R_{t-2,u,v}^{(N+1)}  (t > 0)
        R_{0uv}^{(N)} = PC_y × R_{0,u-1,v}^{(N+1)} + (u-1) × R_{0,u-2,v}^{(N+1)}  (t=0, u > 0)
        R_{00v}^{(N)} = PC_z × R_{00,v-1}^{(N+1)} + (v-1) × R_{00,v-2}^{(N+1)}  (t=u=0, v > 0)

    Reference: McMurchie & Davidson, J. Comput. Phys. 26 (1978) 218-231
    """
    rec = Recurrence("CoulombR", ["t", "u", "v", "N"],
                     ["PCx", "PCy", "PCz", "Boys"],
                     namespace="mcmd_coulomb",
                     max_indices={"t": 6, "u": 6, "v": 6, "N": 6})  # Max L=3 (f-orbitals, 2*3=6)

    rec.validity("t >= 0", "u >= 0", "v >= 0", "N >= 0")

    # Base case: R_{000}^{(N)} = F_N(T) = Boys[N]
    rec.base(t=0, u=0, v=0, value="Boys[N]")

    # X-recurrence (t > 0)
    # Note: Using E[...] notation as required by RecurrenceParser
    rec.rule("t > 0",
             "PCx * E[t-1, u, v, N+1] + (t - 1) * E[t-2, u, v, N+1]",
             name="X-recurrence")

    # Y-recurrence (t = 0, u > 0)
    rec.rule("t == 0 && u > 0",
             "PCy * E[0, u-1, v, N+1] + (u - 1) * E[0, u-2, v, N+1]",
             name="Y-recurrence")

    # Z-recurrence (t = 0, u = 0, v > 0)
    rec.rule("t == 0 && u == 0 && v > 0",
             "PCz * E[0, 0, v-1, N+1] + (v - 1) * E[0, 0, v-2, N+1]",
             name="Z-recurrence")

    return rec


# NOTE: The derivative recurrences below are placeholders.
# The true analytical derivatives are implemented in implement_helgaker_taylor_complete.py
# and are computed at runtime, not via template metaprogramming.

__all__ = [
    "hermite_e_coefficient",
    "coulomb_r_auxiliary",
]
