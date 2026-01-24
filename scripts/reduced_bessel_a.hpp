#pragma once

#include <type_traits>
#include <vectorclass.h>

namespace bessel_sto {


template<int n, typename Enable = void>
struct ReducedBesselACoeff {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d /*a0*/, Vec8d /*a1*/) {
        return Vec8d(0.0);
    }
};

template<>
struct ReducedBesselACoeff<0, void> {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d a0, Vec8d /*a1*/) {
        return a0;
    }
};

template<>
struct ReducedBesselACoeff<1, void> {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d /*a0*/, Vec8d a1) {
        return a1;
    }
};

template<int n>
struct ReducedBesselACoeff<
    n,
    typename std::enable_if<(n > 1) && (n >= 0)>::type
> {
    static Vec8d compute(Vec8d inv_x, Vec8d a0, Vec8d a1) {
        // Upward recurrence (stable)
        return ReducedBesselACoeff<n - 2>::compute(inv_x, a0, a1)
                  + (Vec8d(2*n-1) * inv_x) * ReducedBesselACoeff<n - 1>::compute(inv_x, a0, a1);
    }
};

} // namespace bessel_sto