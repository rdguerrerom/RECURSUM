"""Modified spherical Bessel function recurrences for STO integrals."""

from ..codegen import Recurrence


def modified_spherical_bessel_i() -> Recurrence:
    """Modified spherical Bessel function of the first kind i_n(x).

    DLMF section 10.47 - https://dlmf.nist.gov/10.47
    """
    rec = Recurrence("ModSphBesselI", ["n"], ["inv_x", "i0", "i1"],
                     namespace="bessel_sto", max_indices={"n": 15})
    rec.validity("n >= 0")
    rec.base(n=0, value="i0")
    rec.base(n=1, value="i1")
    rec.rule("n > 1", "E[n-2] + (-(2*n-1)) * inv_x * E[n-1]",
             name="Upward recurrence")
    return rec


def modified_spherical_bessel_k() -> Recurrence:
    """Modified spherical Bessel function of the second kind k_n(x).

    DLMF section 10.47 - https://dlmf.nist.gov/10.47
    """
    rec = Recurrence("ModSphBesselK", ["n"], ["inv_x", "k0", "k1"],
                     namespace="bessel_sto", max_indices={"n": 15})
    rec.validity("n >= 0")
    rec.base(n=0, value="k0")
    rec.base(n=1, value="k1")
    rec.rule("n > 1", "E[n-2] + (2*n-1) * inv_x * E[n-1]",
             name="Upward recurrence (stable)")
    return rec


def reduced_bessel_b() -> Recurrence:
    """Reduced modified spherical Bessel b_n(x) = e^{-x} i_n(x).

    DLMF section 10.47 - https://dlmf.nist.gov/10.47 (scaled form)
    """
    rec = Recurrence("ReducedBesselB", ["n"], ["inv_x", "b0", "b1"],
                     namespace="bessel_sto", max_indices={"n": 15})
    rec.validity("n >= 0")
    rec.base(n=0, value="b0")
    rec.base(n=1, value="b1")
    rec.rule("n > 1", "E[n-2] + (-(2*n-1)) * inv_x * E[n-1]",
             name="Upward recurrence")
    return rec


def reduced_bessel_a() -> Recurrence:
    """Reduced modified spherical Bessel a_n(x) = e^{x} k_n(x).

    DLMF section 10.47 - https://dlmf.nist.gov/10.47 (scaled form)
    """
    rec = Recurrence("ReducedBesselA", ["n"], ["inv_x", "a0", "a1"],
                     namespace="bessel_sto", max_indices={"n": 15})
    rec.validity("n >= 0")
    rec.base(n=0, value="a0")
    rec.base(n=1, value="a1")
    rec.rule("n > 1", "E[n-2] + (2*n-1) * inv_x * E[n-1]",
             name="Upward recurrence (stable)")
    return rec


__all__ = [
    "modified_spherical_bessel_i",
    "modified_spherical_bessel_k",
    "reduced_bessel_b",
    "reduced_bessel_a",
]
