#pragma once

#include <type_traits>
#include <vectorclass.h>

namespace sto_integrals {


template<int n, int l, typename Enable = void>
struct STOAuxBCoeff {
    static Vec8d compute(Vec8d /*x*/, Vec8d /*inv_x*/, Vec8d /*alpha*/, Vec8d /*beta*/) {
        return Vec8d(0.0);
    }
};

template<>
struct STOAuxBCoeff<0, 0, void> {
    static Vec8d compute(Vec8d /*x*/, Vec8d /*inv_x*/, Vec8d /*alpha*/, Vec8d /*beta*/) {
        return Vec8d(1.0);
    }
};

template<int n, int l>
struct STOAuxBCoeff<
    n, l,
    typename std::enable_if<(n == l) && (l > 0) && (n >= 0) && (l >= 0) && (n >= l)>::type
> {
    static Vec8d compute(Vec8d x, Vec8d inv_x, Vec8d alpha, Vec8d beta) {
        // Diagonal recursion
        return alpha * STOAuxBCoeff<n - 1, l - 1>::compute(x, inv_x, alpha, beta);
    }
};

template<int n, int l>
struct STOAuxBCoeff<
    n, l,
    typename std::enable_if<(n > l) && (l == 0) && (n >= 0) && (l >= 0) && (n >= l)>::type
> {
    static Vec8d compute(Vec8d x, Vec8d inv_x, Vec8d alpha, Vec8d beta) {
        // l=0 column
        return STOAuxBCoeff<n - 1, l>::compute(x, inv_x, alpha, beta)
                  + (Vec8d(2*n-1) * inv_x) * STOAuxBCoeff<n - 1, l>::compute(x, inv_x, alpha, beta);
    }
};

template<int n, int l>
struct STOAuxBCoeff<
    n, l,
    typename std::enable_if<(n > l) && (l > 0) && (n >= 0) && (l >= 0) && (n >= l)>::type
> {
    static Vec8d compute(Vec8d x, Vec8d inv_x, Vec8d alpha, Vec8d beta) {
        // General two-term
        return alpha * STOAuxBCoeff<n - 1, l - 1>::compute(x, inv_x, alpha, beta)
                  + beta * STOAuxBCoeff<n - 1, l + 1>::compute(x, inv_x, alpha, beta);
    }
};

} // namespace sto_integrals