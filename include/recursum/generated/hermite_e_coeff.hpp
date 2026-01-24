#pragma once

#include <type_traits>
#include <recursum/vectorclass.h>

namespace mcmd_hermite {


template<int nA, int nB, int t, typename Enable = void>
struct HermiteECoeff {
    static Vec8d compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/) {
        return Vec8d(0.0);
    }
};

template<>
struct HermiteECoeff<0, 0, 0, void> {
    static Vec8d compute(Vec8d /*PA*/, Vec8d /*PB*/, Vec8d /*p*/) {
        return Vec8d(1.0);
    }
};

template<int nA, int nB, int t>
struct HermiteECoeff<
    nA, nB, t,
    typename std::enable_if<(nA > 0) && (nB == 0) && (t == 0) && (nA >= 0) && (nB >= 0) && (t >= 0) && (t <= nA + nB)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        // A-side t=0 (includes E_{t+1})
        return PA * HermiteECoeff<nA - 1, nB, t>::compute(PA, PB, p) + Vec8d(t + 1) * HermiteECoeff<nA - 1, nB, t + 1>::compute(PA, PB, p);
    }
};

template<int nA, int nB, int t>
struct HermiteECoeff<
    nA, nB, t,
    typename std::enable_if<(nA == 0) && (nB > 0) && (t == 0) && (nA >= 0) && (nB >= 0) && (t >= 0) && (t <= nA + nB)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        // B-side t=0 (includes E_{t+1})
        return PB * HermiteECoeff<nA, nB - 1, t>::compute(PA, PB, p) + Vec8d(t + 1) * HermiteECoeff<nA, nB - 1, t + 1>::compute(PA, PB, p);
    }
};

template<int nA, int nB, int t>
struct HermiteECoeff<
    nA, nB, t,
    typename std::enable_if<(nA > 0) && (nB == 0) && (t > 0) && (nA >= 0) && (nB >= 0) && (t >= 0) && (t <= nA + nB)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        // A-side t>0
        return 0.5 / p * HermiteECoeff<nA - 1, nB, t - 1>::compute(PA, PB, p) + PA * HermiteECoeff<nA - 1, nB, t>::compute(PA, PB, p) + Vec8d(t + 1) * HermiteECoeff<nA - 1, nB, t + 1>::compute(PA, PB, p);
    }
};

template<int nA, int nB, int t>
struct HermiteECoeff<
    nA, nB, t,
    typename std::enable_if<(nA == 0) && (nB > 0) && (t > 0) && (nA >= 0) && (nB >= 0) && (t >= 0) && (t <= nA + nB)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        // B-side t>0
        return 0.5 / p * HermiteECoeff<nA, nB - 1, t - 1>::compute(PA, PB, p) + PB * HermiteECoeff<nA, nB - 1, t>::compute(PA, PB, p) + Vec8d(t + 1) * HermiteECoeff<nA, nB - 1, t + 1>::compute(PA, PB, p);
    }
};

template<int nA, int nB, int t>
struct HermiteECoeff<
    nA, nB, t,
    typename std::enable_if<(nA > 0) && (nB > 0) && (t == 0) && (nA >= 0) && (nB >= 0) && (t >= 0) && (t <= nA + nB)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        // General t=0 (increment-i only, includes E_{t+1})
        return PA * HermiteECoeff<nA - 1, nB, t>::compute(PA, PB, p) + Vec8d(t + 1) * HermiteECoeff<nA - 1, nB, t + 1>::compute(PA, PB, p);
    }
};

template<int nA, int nB, int t>
struct HermiteECoeff<
    nA, nB, t,
    typename std::enable_if<(nA > 0) && (nB > 0) && (t > 0) && (nA >= 0) && (nB >= 0) && (t >= 0) && (t <= nA + nB)>::type
> {
    static Vec8d compute(Vec8d PA, Vec8d PB, Vec8d p) {
        // General t>0 (increment-i only)
        return 0.5 / p * HermiteECoeff<nA - 1, nB, t - 1>::compute(PA, PB, p) + PA * HermiteECoeff<nA - 1, nB, t>::compute(PA, PB, p) + Vec8d(t + 1) * HermiteECoeff<nA - 1, nB, t + 1>::compute(PA, PB, p);
    }
};

} // namespace mcmd_hermite