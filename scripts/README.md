# RECURSUM: Recurrence Relation Code Generator

A fully automated Python DSL and C++ code generator for zero-overhead recurrence relation evaluation using template metaprogramming with seamless Python integration.

## Quick Start

```bash
# Clone and install
git clone <repository-url>
cd RECURSUM
pip install -e .[dev]

# Use in Python
python
>>> import recursum._orthogonal as ortho
>>> import numpy as np
>>> x = np.array([0.5, 1.0, 1.5])
>>> result = ortho.legendre(n=5, x=x)
```

## Features

- **Fully Automated Pipeline**: Code generation, compilation, and testing
- **Zero Runtime Overhead**: Compile-time template metaprogramming
- **SciPy Validation**: Automatic comparison tests against SciPy
- **Interactive Notebooks**: Auto-generated Jupyter notebooks for visualization
- **Hot-Reload**: Watch mode for rapid development
- **SIMD Optimized**: Vec8d vectorization for 8x performance

## Installation

### Requirements

- Python 3.8+
- CMake 3.18+
- C++ compiler with C++17 support
- pybind11 (installed automatically)

### Development Installation

```bash
pip install -e .[dev]
```

This will:
1. Generate C++ template headers from recurrence definitions
2. Generate runtime dispatchers
3. Generate pybind11 bindings
4. Compile C++ extensions
5. Generate pytest tests and Jupyter notebooks

---

# Legacy: Original DSL Documentation

A Python-based DSL and code generator that produces SFINAE-enabled C++ template implementations for solving arbitrary recurrence relations at compile time.

## Overview

This tool transforms declarative recurrence relation specifications into optimized C++ template metaprogramming code. The generated code uses SFINAE (`std::enable_if`) to select the appropriate specialization at compile time, resulting in zero runtime overhead for recurrence evaluation.

## Features

- **Einsum-inspired DSL**: Intuitive notation for expressing recurrence relations
- **Multi-index support**: Handle recurrences with arbitrary numbers of indices
- **Automatic SFINAE generation**: Constraints compiled to `enable_if` conditions
- **Branch averaging**: Built-in support for symmetric recurrences (e.g., McMurchie-Davidson)
- **Scaling factors**: Native support for divided recurrences like Legendre polynomials
- **SIMD-ready output**: Generated code uses `Vec8d` for vectorization compatibility

## Installation

No dependencies required beyond Python 3.7+. Simply copy `recurrence_codegen.py` to your project.

## Quick Start

```python
from recurrence_codegen import Recurrence

# Define a simple recurrence (Chebyshev T_n)
rec = Recurrence("ChebyshevT", ["n"], ["x"], namespace="chebyshev")
rec.validity("n >= 0")
rec.base(n=0, value=1.0)
rec.base(n=1, value="x")
rec.rule("n > 1", "(2) * x * E[n-1] + (-1) * E[n-2]")

# Generate C++ code
print(rec.generate())
```

## DSL Syntax

### Index Notation
```
E[i,j,t]           # Recurrence evaluated at indices i, j, t
E[i-1,j,t+1]       # With index shifts
E[n-2]             # Single-index recurrence
```

### Coefficient Types
```
x                  # Runtime variable
(n+1)              # Index-dependent coefficient
(2*n-1)            # Compound index expression
0.5                # Numeric constant
aAB                # Named parameter
```

### Expression Syntax
```python
# Sum of terms
"coeff * E[shifts] + coeff * E[shifts] + ..."

# Examples
"x * E[n-1] + E[n-2]"                                    # Fibonacci-like
"(2*n-1) * x * E[n-1] + (-(n-1)) * E[n-2]"              # Legendre (before scaling)
"aAB * E[nA-1, nB, N-1] + PA * E[nA-1, nB, N]"          # Hermite coefficient
```

### Constraint Syntax
```python
"n > 0"                    # Single constraint
"n >= 0 && m >= 0"         # Combined constraints
["nA > 0", "nB == 0"]      # List form
```

## API Reference

### `Recurrence(name, indices, runtime_vars, vec_type="Vec8d", namespace="")`

Create a new recurrence definition.

- `name`: Struct name prefix (e.g., "Hermite" → `HermiteCoeff`)
- `indices`: List of template parameter names (compile-time integers)
- `runtime_vars`: List of runtime parameter names (Vec8d values)
- `vec_type`: SIMD vector type (default: "Vec8d")
- `namespace`: C++ namespace for generated code

### `.validity(*constraints)`

Set validity constraints determining when results are non-zero.

```python
rec.validity("n >= 0", "l >= 0", "n >= l")
```

### `.base(value, **index_values)`

Define a base case with fixed index values.

```python
rec.base(n=0, value=1.0)
rec.base(n=1, value="x")
rec.base(l=0, m=0, value=1.0)
```

### `.rule(constraints, expression, scale=None, name="")`

Add a recurrence rule.

```python
# Simple rule
rec.rule("n > 1", "x * E[n-1] + E[n-2]")

# With scaling: result = expression / scale
rec.rule("n > 1", 
         "(2*n-1) * x * E[n-1] + (-(n-1)) * E[n-2]",
         scale="1/n",
         name="Three-term recurrence")
```

### `.branch_average(constraints, branches, name="")`

Average multiple equivalent branches (for numerical stability).

```python
rec.branch_average(
    "nA > 0 && nB > 0",
    ["aAB * E[nA, nB-1, N-1] + PB * E[nA, nB-1, N] + (N+1) * E[nA, nB-1, N+1]",
     "aAB * E[nA-1, nB, N-1] + PA * E[nA-1, nB, N] + (N+1) * E[nA-1, nB, N+1]"],
    name="Two-branch average"
)
```

### `.generate()`

Generate the C++ header file content as a string.

## Built-in Recurrences

### Orthogonal Polynomials

| Function | Description | Recurrence |
|----------|-------------|------------|
| `legendre_polynomials()` | Legendre $P_n(x)$ | $P_n = \frac{(2n-1)x P_{n-1} - (n-1)P_{n-2}}{n}$ |
| `chebyshev_T()` | Chebyshev 1st kind $T_n(x)$ | $T_n = 2x T_{n-1} - T_{n-2}$ |
| `chebyshev_U()` | Chebyshev 2nd kind $U_n(x)$ | $U_n = 2x U_{n-1} - U_{n-2}$ |
| `hermite_He()` | Probabilist's Hermite $He_n(x)$ | $He_n = x He_{n-1} - (n-1) He_{n-2}$ |
| `hermite_H()` | Physicist's Hermite $H_n(x)$ | $H_n = 2x H_{n-1} - 2(n-1) H_{n-2}$ |
| `laguerre()` | Laguerre $L_n(x)$ | $L_n = \frac{(2n-1-x)L_{n-1} - (n-1)L_{n-2}}{n}$ |
| `associated_legendre()` | Associated Legendre $P_l^m(x)$ | Multi-rule with diagonal, off-diagonal, general |

### Quantum Chemistry

| Function | Description | Use Case |
|----------|-------------|----------|
| `hermite_coefficients()` | McMurchie-Davidson $E^{i,j}_t$ | GTO molecular integrals |
| `boys_function()` | Boys function $F_n(T)$ | Gaussian integral auxiliary |
| `modified_spherical_bessel_i()` | $i_n(x) = \sqrt{\pi/(2x)} I_{n+1/2}(x)$ | STO integrals |
| `modified_spherical_bessel_k()` | $k_n(x) = \sqrt{\pi/(2x)} K_{n+1/2}(x) \cdot (2/\pi)$ | STO integrals |
| `reduced_bessel_b()` | $b_n(x) = e^{-x} i_n(x)$ | Overflow-safe STO integrals |
| `reduced_bessel_a()` | $a_n(x) = e^{x} k_n(x)$ | Polynomial form for STOs |
| `sto_auxiliary_B()` | Filter-Steinborn $B_{n,l}$ | STO translation coefficients |
| `gaunt_coefficients()` | $\int Y_{l_1}^{m_1} Y_{l_2}^{m_2} Y_{l_3}^{m_3} d\Omega$ | Angular momentum coupling |

### Rys Quadrature (ERIs)

| Function | Description | Use Case |
|----------|-------------|----------|
| `rys_2d_integral()` | $I_x(n, 0, m, 0)$ 2D integral | VRR for Rys quadrature |
| `rys_horizontal_transfer()` | HRR/transfer step | Angular momentum distribution |
| `rys_vrr_full()` | Full 4-index VRR | Modern integral codes |
| `rys_polynomial_recursion()` | $R_n(X, t)$ polynomials | Quadrature points |

### Combinatorics

| Function | Description | Recurrence |
|----------|-------------|------------|
| `binomial_coefficients()` | $\binom{n}{k}$ | Pascal's rule: $\binom{n}{k} = \binom{n-1}{k-1} + \binom{n-1}{k}$ |
| `fibonacci()` | Generalized Fibonacci | $F_n = x F_{n-1} + F_{n-2}$ |

## Modified Spherical Bessel Functions for STO Integrals

The modified spherical Bessel functions appear in Slater-type orbital (STO) multi-center integral evaluation. They are not part of the STO definition itself, but arise as auxiliary functions in analytical integral schemes (Guseinov, Harris-Michaels, Filter-Steinborn).

### Definitions

**Modified spherical Bessel of the first kind:**
$$i_n(x) = \sqrt{\frac{\pi}{2x}} I_{n+1/2}(x)$$

Base cases:
- $i_0(x) = \sinh(x)/x$
- $i_1(x) = \cosh(x)/x - \sinh(x)/x^2$

**Modified spherical Bessel of the second kind:**
$$k_n(x) = \sqrt{\frac{\pi}{2x}} K_{n+1/2}(x) \cdot \frac{2}{\pi}$$

Base cases:
- $k_0(x) = \frac{\pi}{2} e^{-x}/x$
- $k_1(x) = k_0(x)(1 + 1/x)$

### Recurrence Relations

```
i_n(x) = i_{n-2}(x) - (2n-1)/x · i_{n-1}(x)    # Upward (unstable for large n)
k_n(x) = k_{n-2}(x) + (2n-1)/x · k_{n-1}(x)    # Upward (stable)
```

### Reduced (Scaled) Forms

To avoid numerical overflow, the reduced Bessel functions factor out the exponential:

$$b_n(x) = e^{-x} i_n(x)$$
$$a_n(x) = e^{x} k_n(x)$$

The reduced forms satisfy the same recurrence relations as their parent functions (the exponential cancels). For $a_n(x)$, the result is purely polynomial in $1/x$:

$$a_n(x) = \frac{\pi}{2x} \sum_{k=0}^{n} \frac{(n+k)!}{k!(n-k)!} (2x)^{-k}$$

### Numerical Stability

- **$k_n$ and $a_n$**: Upward recurrence is stable
- **$i_n$ and $b_n$**: Upward recurrence is unstable for large $n$; use Miller's backward algorithm starting from an asymptotic estimate

## Rys Quadrature for Electron Repulsion Integrals

The Rys polynomial method (King & Dupuis 1976, Rys, Dupuis & King 1983) evaluates two-electron repulsion integrals via numerical quadrature using specialized orthogonal polynomials.

### The ERI Problem

The two-electron repulsion integral over Gaussian basis functions:

$$\langle ij | kl \rangle = \int\int \phi_i(1)\phi_j(1) \frac{1}{r_{12}} \phi_k(2)\phi_l(2) \, d\tau_1 d\tau_2$$

### Rys Quadrature Formula

The integral is computed as an N-point quadrature:

$$\langle ij | kl \rangle = 2\sqrt{\rho/\pi} \sum_{\alpha=0}^{N} I_x(t_\alpha) I_y(t_\alpha) I_z(t_\alpha) W_\alpha$$

where $t_\alpha$ and $W_\alpha$ are roots and weights of the N-th order Rys polynomial, and N is chosen as the smallest integer greater than $\lambda/2$ (where $\lambda$ is the total angular momentum).

### 2D Integral Vertical Recurrence (VRR)

The 2D integrals $I_x(n, 0, m, 0)$ are built recursively (Eqs. 15-16 in Augspurger et al.):

**Bra-side recurrence (increase n):**
```
I_x(n+1, 0, m, 0) = n·B₁₀·I_x(n-1, 0, m, 0)
                  + m·B₀₀·I_x(n, 0, m-1, 0)  
                  + C₀₀·I_x(n, 0, m, 0)
```

**Ket-side recurrence (increase m):**
```
I_x(n, 0, m+1, 0) = m·B'₁₀·I_x(n, 0, m-1, 0)
                  + n·B₀₀·I_x(n-1, 0, m, 0)
                  + C'₀₀·I_x(n, 0, m, 0)
```

### Recurrence Coefficients

The coefficients depend on the Rys root $t_\alpha$:

| Coefficient | Formula | Description |
|-------------|---------|-------------|
| $B_{00}$ | $\frac{t^2}{2(A+B)}$ | Bra-ket coupling |
| $B_{10}$ | $\frac{1}{2A} - \frac{B \cdot t^2}{2A(A+B)}$ | Bra vertical |
| $B_{01}$ | $\frac{1}{2B} - \frac{A \cdot t^2}{2B(A+B)}$ | Ket vertical |
| $C_{00}$ | $(P_x - x_i) + \frac{B(Q_x - P_x)}{A+B} t^2$ | Bra shift |
| $C'_{00}$ | $(Q_x - x_k) + \frac{A(P_x - Q_x)}{A+B} t^2$ | Ket shift |

where:
- $A = \alpha_i + \alpha_j$ (combined bra exponent)
- $B = \alpha_k + \alpha_l$ (combined ket exponent)
- $P = (\alpha_i r_i + \alpha_j r_j)/A$ (bra center)
- $Q = (\alpha_k r_k + \alpha_l r_l)/B$ (ket center)

### Horizontal Recurrence (HRR) / Transfer

After building $I_x(i_x+j_x, 0, k_x+l_x, 0)$, the angular momentum is "transferred" to individual indices (Eqs. 17-18):

**Bra transfer:**
```
I_x(i_x, j_x, m, 0) = I_x(i_x+1, j_x-1, m, 0) + (x_i - x_j)·I_x(i_x, j_x-1, m, 0)
```

**Ket transfer:**
```
I_x(i_x, j_x, k_x, l_x) = I_x(i_x, j_x, k_x+1, l_x-1) + (x_k - x_l)·I_x(i_x, j_x, k_x, l_x-1)
```

### Direct Summation Alternative (Eq. 19)

The transfer can be done more efficiently (2.5-3× faster) via direct summation:

$$I_x(i,j,m,0) = \sum_{n=0}^{j_x} q_n \cdot I_x(i+n, 0, m, 0)$$

where $q_n = \binom{j_x}{n} (x_j - x_i)^{j_x - n}$

### Rys Polynomial Recursion

The Rys polynomials $R_n(X, t)$ satisfy:

$$\beta_{n+1} R_{n+1} = (t^2 - \beta'_n) R_n - \beta_n R_{n-1}$$

where the coefficients $\beta$ are computed from weighted integrals over $[0,1]$. The roots of these polynomials provide the quadrature points.

## Generated Code Structure

The generator produces C++ headers with this structure:

```cpp
#pragma once
#include <type_traits>
#include <vectorclass.h>

namespace hermite {

// Primary template: returns 0 for invalid indices
template<int nA, int nB, int N, typename Enable = void>
struct HermiteCoeff {
    static Vec8d compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*aAB*/) {
        return Vec8d(0.0);
    }
};

// Base case: E^{0,0}_0 = 1
template<>
struct HermiteCoeff<0, 0, 0, void> {
    static Vec8d compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*aAB*/) {
        return Vec8d(1.0);
    }
};

// Rule with SFINAE constraints
template<int nA, int nB, int N>
struct HermiteCoeff<
    nA, nB, N,
    typename std::enable_if<(nA == 0) && (nB > 0) && ...>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d aAB) {
        return aAB * HermiteCoeff<nA, nB - 1, N - 1>::compute(PA, PB, aAB)
             + PB * HermiteCoeff<nA, nB - 1, N>::compute(PA, PB, aAB)
             + Vec8d(N+1) * HermiteCoeff<nA, nB - 1, N + 1>::compute(PA, PB, aAB);
    }
};

} // namespace hermite
```

## Command-Line Usage

```bash
# Print example recurrences to stdout
python recurrence_codegen.py

# Generate all built-in examples to a directory
python recurrence_codegen.py --generate ./output_headers/
```

## Example: Custom Recurrence

Define the Jacobi polynomials $P_n^{(\alpha,\beta)}(x)$:

```python
from recurrence_codegen import Recurrence

def jacobi_polynomials() -> Recurrence:
    """Jacobi polynomials P_n^{(α,β)}(x)."""
    rec = Recurrence("Jacobi", ["n"], ["x", "a1", "a2", "a3", "a4"],
                     namespace="jacobi")
    # a1..a4 are precomputed coefficients depending on α, β, n
    rec.validity("n >= 0")
    rec.base(n=0, value=1.0)
    rec.base(n=1, value="a1")  # (α-β)/2 + (α+β+2)x/2
    
    # General recurrence (coefficients passed at runtime)
    rec.rule("n > 1",
             "a2 * x * E[n-1] + a3 * E[n-1] + a4 * E[n-2]",
             name="Three-term recurrence")
    
    return rec
```

## Performance Characteristics

- **Compile-time**: $O(L^k)$ template instantiations for max index $L$ and $k$ indices
- **Runtime**: $O(1)$ constant-time evaluation per call
- **Memory**: Register-only computation, no heap allocation
- **SIMD**: Full vectorization preserved through `Vec8d` operations

## References

- McMurchie, L. E.; Davidson, E. R. J. Comput. Phys. 26 (1978) 218-231
- Helgaker, T.; Taylor, P. R. Theor. Chim. Acta 83 (1992) 177-183
- Filter, E.; Steinborn, E. O. Phys. Rev. A 18 (1978) 1-11
- Guseinov, I. I. J. Phys. B 3 (1970) 1399-1412
- Harris, F. E.; Michels, H. H. Adv. Chem. Phys. 13 (1967) 205-266
- King, H. F.; Dupuis, M. J. Comput. Phys. 21 (1976) 144-165
- Dupuis, M.; Rys, J.; King, H. F. J. Chem. Phys. 65 (1976) 111-116
- Rys, J.; Dupuis, M.; King, H. F. J. Comput. Chem. 4 (1983) 154-175
- Augspurger, J. D.; Bernholdt, D. E.; Dykstra, C. E. J. Comput. Chem. 11 (1990) 972-977
- Obara, S.; Saika, A. J. Chem. Phys. 84 (1986) 3963
- Head-Gordon, M.; Pople, J. A. J. Chem. Phys. 89 (1988) 5777-5786

## License

MIT License - See LICENSE file for details.
