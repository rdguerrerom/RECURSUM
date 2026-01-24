#pragma once

#include <type_traits>
#include <vectorclass.h>

namespace bessel_sto {


template<int n, typename Enable = void>
struct ModSphBesselKCoeff {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d /*k0*/, Vec8d /*k1*/) {
        return Vec8d(0.0);
    }
};

template<>
struct ModSphBesselKCoeff<0, void> {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d k0, Vec8d /*k1*/) {
        return k0;
    }
};

template<>
struct ModSphBesselKCoeff<1, void> {
    static Vec8d compute(Vec8d /*inv_x*/, Vec8d /*k0*/, Vec8d k1) {
        return k1;
    }
};

template<int n>
struct ModSphBesselKCoeff<
    n,
    typename std::enable_if<(n > 1) && (n >= 0)>::type
> {
    static Vec8d compute(Vec8d inv_x, Vec8d k0, Vec8d k1) {
        // Upward recurrence (stable)
        return ModSphBesselKCoeff<n - 2>::compute(inv_x, k0, k1)
                  + (Vec8d(2*n-1) * inv_x) * ModSphBesselKCoeff<n - 1>::compute(inv_x, k0, k1);
    }
};

} // namespace bessel_sto