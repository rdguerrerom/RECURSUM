#pragma once

#include <type_traits>
#include <vectorclass.h>

namespace bessel_sto {


template<int n, typename Enable = void>
struct ModSphBesselICoeff {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d /*i0*/, Vec8d /*i1*/) {
        return Vec8d(0.0);
    }
};

template<>
struct ModSphBesselICoeff<0, void> {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d i0, Vec8d /*i1*/) {
        return i0;
    }
};

template<>
struct ModSphBesselICoeff<1, void> {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d /*i0*/, Vec8d i1) {
        return i1;
    }
};

template<int n>
struct ModSphBesselICoeff<
    n,
    typename std::enable_if<(n > 1) && (n >= 0)>::type
> {
    static Vec8d compute(Vec8d inv_x, Vec8d i0, Vec8d i1) {
        // Upward recurrence
        return ModSphBesselICoeff<n - 2>::compute(inv_x, i0, i1)
                  + (Vec8d(-(2*n-1)) * inv_x) * ModSphBesselICoeff<n - 1>::compute(inv_x, i0, i1);
    }
};

} // namespace bessel_sto