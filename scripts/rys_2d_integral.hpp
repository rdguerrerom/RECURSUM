#pragma once

#include <type_traits>
#include <vectorclass.h>

namespace rys_quadrature {


template<int n, int m, typename Enable = void>
struct Rys2DCoeff {
    static Vec8d compute(Vec8d /*B00*/, Vec8d /*B10*/, Vec8d /*B01*/, Vec8d /*C00*/, Vec8d /*C00p*/) {
        return Vec8d(0.0);
    }
};

template<>
struct Rys2DCoeff<0, 0, void> {
    static Vec8d compute(Vec8d /*B00*/, Vec8d /*B10*/, Vec8d /*B01*/, Vec8d /*C00*/, Vec8d /*C00p*/) {
        return Vec8d(1.0);
    }
};

template<int n, int m>
struct Rys2DCoeff<
    n, m,
    typename std::enable_if<(n > 0) && (m == 0) && (n >= 0) && (m >= 0)>::type
> {
    static Vec8d compute(Vec8d B00, Vec8d B10, Vec8d B01, Vec8d C00, Vec8d C00p) {
        // Bra VRR (m=0) [Eq. 15]
        return (Vec8d(n-1) * B10) * Rys2DCoeff<n - 2, m>::compute(B00, B10, B01, C00, C00p)
                  + C00 * Rys2DCoeff<n - 1, m>::compute(B00, B10, B01, C00, C00p);
    }
};

template<int n, int m>
struct Rys2DCoeff<
    n, m,
    typename std::enable_if<(m > 0) && (n == 0) && (n >= 0) && (m >= 0)>::type
> {
    static Vec8d compute(Vec8d B00, Vec8d B10, Vec8d B01, Vec8d C00, Vec8d C00p) {
        // Ket VRR (n=0) [Eq. 16]
        return (Vec8d(m-1) * B01) * Rys2DCoeff<n, m - 2>::compute(B00, B10, B01, C00, C00p)
                  + C00p * Rys2DCoeff<n, m - 1>::compute(B00, B10, B01, C00, C00p);
    }
};

template<int n, int m>
struct Rys2DCoeff<
    n, m,
    typename std::enable_if<(n > 0) && (m > 0) && (n >= 0) && (m >= 0)>::type
> {
    static Vec8d compute(Vec8d B00, Vec8d B10, Vec8d B01, Vec8d C00, Vec8d C00p) {
        // General VRR [Eq. 15]
        return (Vec8d(n-1) * B10) * Rys2DCoeff<n - 2, m>::compute(B00, B10, B01, C00, C00p)
                  + (Vec8d(m) * B00) * Rys2DCoeff<n - 1, m - 1>::compute(B00, B10, B01, C00, C00p)
                  + C00 * Rys2DCoeff<n - 1, m>::compute(B00, B10, B01, C00, C00p);
    }
};

} // namespace rys_quadrature