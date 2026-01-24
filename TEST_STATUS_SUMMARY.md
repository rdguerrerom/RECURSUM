# RECURSUM Test Status Summary

**Date:** 2026-01-24
**Total Recurrences:** 35
**Total Tests:** 206 (166 passed, 40 failed, 0 skipped)
**Pass Rate:** 80.6%

## Status by Recurrence

### Fully Passing (22 recurrences - 63%)

These recurrences pass all generated tests:

1. **AssocLegendre** - Associated Legendre polynomials (2/2 tests)
2. **BesselY** - Bessel function Y_n (9/9 tests)
3. **Binomial** - Binomial coefficients (2/2 tests)
4. **Boys** - Boys function (2/2 tests)
5. **Fibonacci** - Fibonacci sequence (2/2 tests)
6. **Gaunt** - Gaunt coefficients (2/2 tests)
7. **Gegenbauer** - Gegenbauer polynomials (9/9 tests)
8. **Hermite** - Hermite coefficients (McMurchie-Davidson) (2/2 tests)
9. **HermiteH** - Hermite polynomials H_n (9/9 tests)
10. **HermiteHe** - Probabilist's Hermite polynomials (9/9 tests)
11. **Jacobi** - Jacobi polynomials (9/9 tests)
12. **Laguerre** - Laguerre polynomials (9/9 tests)
13. **ModifiedBesselK** - Modified Bessel K_n (9/9 tests)
14. **ModSphBesselK** - Modified spherical Bessel k_n (9/9 tests)
15. **ReducedBesselA** - Reduced Bessel A (2/2 tests)
16. **ReducedBesselB** - Reduced Bessel B (2/2 tests)
17. **Rys2D** - Rys 2D integrals (2/2 tests)
18. **RysHRR** - Rys horizontal recurrence (2/2 tests)
19. **RysPoly** - Rys polynomials (2/2 tests)
20. **RysVRRFull** - Rys VRR full (2/2 tests)
21. **SphericalBesselY** - Spherical Bessel y_n (9/9 tests)
22. **STOAuxB** - STO auxiliary B function (2/2 tests)

### Partially Failing (11 recurrences - 31%)

These recurrences have some failing tests, likely due to numerical precision issues at high orders or edge cases:

1. **AiryAi** - Airy Ai function (4/9 passed, 5/9 failed)
2. **AiryBi** - Airy Bi function (4/9 passed, 5/9 failed)
3. **Bernoulli** - Bernoulli polynomials (4/9 passed, 5/9 failed) **[INCORRECT RECURRENCE]**
4. **BesselJ** - Bessel function J_n (5/9 passed, 4/9 failed)
5. **ChebyshevT** - Chebyshev polynomials T_n (8/9 passed, 1/9 failed)
6. **ChebyshevU** - Chebyshev polynomials U_n (8/9 passed, 1/9 failed)
7. **Euler** - Euler polynomials (4/9 passed, 5/9 failed) **[INCORRECT RECURRENCE]**
8. **Legendre** - Legendre polynomials (8/9 passed, 1/9 failed)
9. **ModifiedBesselI** - Modified Bessel I_n (5/9 passed, 4/9 failed)
10. **ModSphBesselI** - Modified spherical Bessel i_n (5/9 passed, 4/9 failed)
11. **SphericalBesselJ** - Spherical Bessel j_n (4/9 passed, 5/9 failed)

### Skipped (2 recurrences - 6%)

1. **AssocLaguerre** - Associated Laguerre polynomials (0 tests run)
2. **HermiteE** - Hermite E coefficients (0 tests run)

## Critical Issues

### Euler and Bernoulli Polynomials

**Finding:** Euler and Bernoulli polynomials **do not have simple three-term recurrence relations** as implemented in `recursum/recurrences/special.py`.

**Evidence:**
- DLMF §24.5 provides only summation-based identities, not three-term recurrences
- DLMF §24.5.2 for Euler: ∑_{k=0}^n (n choose k) E_k(x) + E_n(x) = 2x^n
- DLMF §24.5.1 for Bernoulli: ∑_{k=0}^{n-1} (n choose k) B_k(x) = nx^{n-1}
- The implemented formula `E_n = (2x E_{n-1} - (n-1) E_{n-2}) / n` is **mathematically incorrect**
- Test comparisons against SymPy show errors of 10-1000x magnitude

**Recommendation:**
Euler and Bernoulli polynomials should be **removed from RECURSUM** or clearly marked as unsupported, since they do not fit the three-term recurrence paradigm that RECURSUM is designed for.

**Alternative Approaches:**
1. Use SymPy directly: `sympy.euler(n, x)` and `sympy.bernoulli(n, x)`
2. Implement the correct DLMF summation formulas (n-term recurrence)
3. Use closed-form expressions where available

### Numerical Precision Issues

Several Bessel-type functions show numerical precision failures at high orders or near zeros. This is a known issue with upward recurrence relations and may require:
- Backward recursion for stability
- Miller's algorithm
- Direct series evaluation at critical points

## References

- [DLMF §24.5 - Recurrence Relations for Bernoulli and Euler Polynomials](https://dlmf.nist.gov/24.5)
- [DLMF §24.4 - Derivatives and Differences](https://dlmf.nist.gov/24.4)

## Next Steps

1. Remove or document limitations of Euler and Bernoulli polynomial implementations
2. Investigate numerical stability issues in Bessel-type functions
3. Add stability warnings for high-order evaluations
4. Consider implementing backward recursion for stable evaluation
