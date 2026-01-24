#pragma once

#include <type_traits>
#include <vectorclass.h>

namespace bessel_sto {


template<int n, typename Enable = void>
struct ReducedBesselBCoeff {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d /*b0*/, Vec8d /*b1*/) {
        return Vec8d(0.0);
    }
};

template<>
struct ReducedBesselBCoeff<0, void> {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d b0, Vec8d /*b1*/) {
        return b0;
    }
};

template<>
struct ReducedBesselBCoeff<1, void> {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d /*b0*/, Vec8d b1) {
        return b1;
    }
};

template<int n>
struct ReducedBesselBCoeff<
    n,
    typename std::enable_if<(n > 1) && (n >= 0)>::type
> {
    static Vec8d compute(Vec8d inv_x, Vec8d b0, Vec8d b1) {
        // Upward recurrence
        return ReducedBesselBCoeff<n - 2>::compute(inv_x, b0, b1)
                  + (Vec8d(-(2*n-1)) * inv_x) * ReducedBesselBCoeff<n - 1>::compute(inv_x, b0, b1);
    }
};

} // namespace bessel_sto