#pragma once

#include <type_traits>
#include <vectorclass.h>

namespace rys_quadrature {


template<int i, int j, int k, int l, typename Enable = void>
struct RysHRRCoeff {
    static Vec8d compute(Vec8d /*AB*/, Vec8d /*CD*/) {
        return Vec8d(0.0);
    }
};

template<>
struct RysHRRCoeff<0, 0, 0, 0, void> {
    static Vec8d compute(Vec8d /*AB*/, Vec8d /*CD*/) {
        return Vec8d(1.0);
    }
};

template<int i, int j, int k, int l>
struct RysHRRCoeff<
    i, j, k, l,
    typename std::enable_if<(j > 0) && (l == 0) && (i >= 0) && (j >= 0) && (k >= 0) && (l >= 0)>::type
> {
    static Vec8d compute(Vec8d AB, Vec8d CD) {
        // Bra HRR [Eq. 17]
        return RysHRRCoeff<i + 1, j - 1, k, l>::compute(AB, CD)
                  + AB * RysHRRCoeff<i, j - 1, k, l>::compute(AB, CD);
    }
};

template<int i, int j, int k, int l>
struct RysHRRCoeff<
    i, j, k, l,
    typename std::enable_if<(l > 0) && (i >= 0) && (j >= 0) && (k >= 0) && (l >= 0)>::type
> {
    static Vec8d compute(Vec8d AB, Vec8d CD) {
        // Ket HRR [Eq. 18]
        return RysHRRCoeff<i, j, k + 1, l - 1>::compute(AB, CD)
                  + CD * RysHRRCoeff<i, j, k, l - 1>::compute(AB, CD);
    }
};

} // namespace rys_quadrature