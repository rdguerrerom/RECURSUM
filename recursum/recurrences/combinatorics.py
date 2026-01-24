"""Combinatorics recurrence relations."""

from ..codegen import Recurrence


def binomial_coefficients() -> Recurrence:
    """Binomial coefficients C(n,k)."""
    rec = Recurrence("Binomial", ["n", "k"], [], namespace="combinatorics",
                     max_indices={"n": 10, "k": 10})
    rec.validity("n >= 0", "k >= 0", "n >= k")
    rec.base(n=0, k=0, value=1.0)
    rec.rule("k == 0", "1 * E[n-1, k]", name="k=0 edge")
    rec.rule("n == k", "1 * E[n-1, k-1]", name="n=k edge")
    rec.rule("n > k && k > 0", "E[n-1, k-1] + E[n-1, k]", name="Pascal's rule")
    return rec


def fibonacci() -> Recurrence:
    """Fibonacci-like with parameter x."""
    rec = Recurrence("Fibonacci", ["n"], ["x"], namespace="sequences",
                     max_indices={"n": 20})
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.base(n=1, value="x")
    rec.rule("n > 1", "x * E[n-1] + E[n-2]", name="Fibonacci-like")
    return rec


__all__ = ["binomial_coefficients", "fibonacci"]
