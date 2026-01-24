/**
 * @file hermite_grad_symbolic.hpp
 * @brief Symbolically-generated Hermite E gradients
 *
 * AUTO-GENERATED - DO NOT EDIT MANUALLY
 *
 * Gradients computed by direct differentiation:
 *   dE_dPA[(nA,nB,t)] = ∂E^{nA,nB}_t/∂PA
 *   dE_dPB[(nA,nB,t)] = ∂E^{nA,nB}_t/∂PB
 */

#pragma once

#ifndef RECURSUM_VEC_TYPE
#include <recursum/vectorclass.h>
#endif

namespace recursum {
namespace symbolic {

// =============================================================================
// ∂E/∂PA Gradients
// =============================================================================

inline Vec8d hermite_dE_dPA_0_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_0_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_1_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 1;
}

inline Vec8d hermite_dE_dPA_1_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_1_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PB;
}

inline Vec8d hermite_dE_dPA_1_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return one_over_2p;
}

inline Vec8d hermite_dE_dPA_1_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_1_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PB*PB) + one_over_2p;
}

inline Vec8d hermite_dE_dPA_1_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PB*one_over_2p;
}

inline Vec8d hermite_dE_dPA_1_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_1_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_1_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PB*PB*PB) + 3*PB*one_over_2p;
}

inline Vec8d hermite_dE_dPA_1_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PB*PB)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_1_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_1_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_1_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_1_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PB*PB)*(PB*PB)) + 6*(PB*PB)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_1_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PB*PB*PB)*one_over_2p + 12*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_1_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PB*PB)*(one_over_2p*one_over_2p) + 6*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_1_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_1_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_1_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_2_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA;
}

inline Vec8d hermite_dE_dPA_2_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*one_over_2p;
}

inline Vec8d hermite_dE_dPA_2_0_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_2_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*PB + 2*one_over_2p;
}

inline Vec8d hermite_dE_dPA_2_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*one_over_2p + 2*PB*one_over_2p;
}

inline Vec8d hermite_dE_dPA_2_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_1_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_2_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*(PB*PB) + 2*PA*one_over_2p + 4*PB*one_over_2p;
}

inline Vec8d hermite_dE_dPA_2_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*PB*one_over_2p + 2*(PB*PB)*one_over_2p + 6*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*(one_over_2p*one_over_2p) + 4*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_2_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_2_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*(PB*PB*PB) + 6*PA*PB*one_over_2p + 6*(PB*PB)*one_over_2p + 6*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*(PB*PB)*one_over_2p + 6*PA*(one_over_2p*one_over_2p) + 2*(PB*PB*PB)*one_over_2p + 18*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*PB*(one_over_2p*one_over_2p) + 6*(PB*PB)*(one_over_2p*one_over_2p) + 12*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*(one_over_2p*one_over_2p*one_over_2p) + 6*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_2_3_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_2_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*((PB*PB)*(PB*PB)) + 12*PA*(PB*PB)*one_over_2p + 6*PA*(one_over_2p*one_over_2p) + 8*(PB*PB*PB)*one_over_2p + 24*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 8*PA*(PB*PB*PB)*one_over_2p + 24*PA*PB*(one_over_2p*one_over_2p) + 2*((PB*PB)*(PB*PB))*one_over_2p + 36*(PB*PB)*(one_over_2p*one_over_2p) + 30*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PA*(PB*PB)*(one_over_2p*one_over_2p) + 12*PA*(one_over_2p*one_over_2p*one_over_2p) + 8*(PB*PB*PB)*(one_over_2p*one_over_2p) + 48*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 8*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 12*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 20*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_2_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 8*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_2_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPA_2_4_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_3_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA) + 3*one_over_2p;
}

inline Vec8d hermite_dE_dPA_3_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*one_over_2p;
}

inline Vec8d hermite_dE_dPA_3_0_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_0_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_3_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*PB + 6*PA*one_over_2p + 3*PB*one_over_2p;
}

inline Vec8d hermite_dE_dPA_3_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*one_over_2p + 6*PA*PB*one_over_2p + 9*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*(one_over_2p*one_over_2p) + 3*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_1_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_1_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_3_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*(PB*PB) + 3*(PA*PA)*one_over_2p + 12*PA*PB*one_over_2p + 3*(PB*PB)*one_over_2p + 9*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA)*PB*one_over_2p + 6*PA*(PB*PB)*one_over_2p + 18*PA*(one_over_2p*one_over_2p) + 18*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*(one_over_2p*one_over_2p) + 12*PA*PB*(one_over_2p*one_over_2p) + 3*(PB*PB)*(one_over_2p*one_over_2p) + 18*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*(one_over_2p*one_over_2p*one_over_2p) + 6*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_2_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_3_2_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_3_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*(PB*PB*PB) + 9*(PA*PA)*PB*one_over_2p + 18*PA*(PB*PB)*one_over_2p + 18*PA*(one_over_2p*one_over_2p) + 3*(PB*PB*PB)*one_over_2p + 27*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 9*(PA*PA)*(PB*PB)*one_over_2p + 9*(PA*PA)*(one_over_2p*one_over_2p) + 6*PA*(PB*PB*PB)*one_over_2p + 54*PA*PB*(one_over_2p*one_over_2p) + 27*(PB*PB)*(one_over_2p*one_over_2p) + 45*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 9*(PA*PA)*PB*(one_over_2p*one_over_2p) + 18*PA*(PB*PB)*(one_over_2p*one_over_2p) + 36*PA*(one_over_2p*one_over_2p*one_over_2p) + 3*(PB*PB*PB)*(one_over_2p*one_over_2p) + 54*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 18*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 9*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 30*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_3_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 9*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_3_3_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_3_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_3_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*((PB*PB)*(PB*PB)) + 18*(PA*PA)*(PB*PB)*one_over_2p + 9*(PA*PA)*(one_over_2p*one_over_2p) + 24*PA*(PB*PB*PB)*one_over_2p + 72*PA*PB*(one_over_2p*one_over_2p) + 3*((PB*PB)*(PB*PB))*one_over_2p + 54*(PB*PB)*(one_over_2p*one_over_2p) + 45*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*(PB*PB*PB)*one_over_2p + 36*(PA*PA)*PB*(one_over_2p*one_over_2p) + 6*PA*((PB*PB)*(PB*PB))*one_over_2p + 108*PA*(PB*PB)*(one_over_2p*one_over_2p) + 90*PA*(one_over_2p*one_over_2p*one_over_2p) + 36*(PB*PB*PB)*(one_over_2p*one_over_2p) + 180*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 18*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 18*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 24*PA*(PB*PB*PB)*(one_over_2p*one_over_2p) + 144*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 3*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p) + 108*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 135*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_3_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 36*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 60*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 12*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 120*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_3_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 24*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 18*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 45*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 12*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPA_3_4_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_3_4_7(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_4_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA) + 12*PA*one_over_2p;
}

inline Vec8d hermite_dE_dPA_4_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*one_over_2p + 12*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_0_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PA*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_0_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_0_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_4_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*PB + 12*(PA*PA)*one_over_2p + 12*PA*PB*one_over_2p + 12*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*one_over_2p + 12*(PA*PA)*PB*one_over_2p + 36*PA*(one_over_2p*one_over_2p) + 12*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*(one_over_2p*one_over_2p) + 12*PA*PB*(one_over_2p*one_over_2p) + 24*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_1_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PA*(one_over_2p*one_over_2p*one_over_2p) + 4*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_1_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_4_1_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_4_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*(PB*PB) + 4*(PA*PA*PA)*one_over_2p + 24*(PA*PA)*PB*one_over_2p + 12*PA*(PB*PB)*one_over_2p + 36*PA*(one_over_2p*one_over_2p) + 24*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 8*(PA*PA*PA)*PB*one_over_2p + 12*(PA*PA)*(PB*PB)*one_over_2p + 36*(PA*PA)*(one_over_2p*one_over_2p) + 72*PA*PB*(one_over_2p*one_over_2p) + 12*(PB*PB)*(one_over_2p*one_over_2p) + 60*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*(one_over_2p*one_over_2p) + 24*(PA*PA)*PB*(one_over_2p*one_over_2p) + 12*PA*(PB*PB)*(one_over_2p*one_over_2p) + 72*PA*(one_over_2p*one_over_2p*one_over_2p) + 48*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 24*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 4*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 40*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_4_2_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 8*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_4_2_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_2_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_4_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*(PB*PB*PB) + 12*(PA*PA*PA)*PB*one_over_2p + 36*(PA*PA)*(PB*PB)*one_over_2p + 36*(PA*PA)*(one_over_2p*one_over_2p) + 12*PA*(PB*PB*PB)*one_over_2p + 108*PA*PB*(one_over_2p*one_over_2p) + 36*(PB*PB)*(one_over_2p*one_over_2p) + 60*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA*PA)*(PB*PB)*one_over_2p + 12*(PA*PA*PA)*(one_over_2p*one_over_2p) + 12*(PA*PA)*(PB*PB*PB)*one_over_2p + 108*(PA*PA)*PB*(one_over_2p*one_over_2p) + 108*PA*(PB*PB)*(one_over_2p*one_over_2p) + 180*PA*(one_over_2p*one_over_2p*one_over_2p) + 12*(PB*PB*PB)*(one_over_2p*one_over_2p) + 180*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA*PA)*PB*(one_over_2p*one_over_2p) + 36*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 72*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 12*PA*(PB*PB*PB)*(one_over_2p*one_over_2p) + 216*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 72*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 180*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_4_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 36*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 36*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 120*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 4*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 120*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_4_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 36*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 12*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 60*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_3_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 12*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_3_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_4_3_7(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPA_4_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*((PB*PB)*(PB*PB)) + 24*(PA*PA*PA)*(PB*PB)*one_over_2p + 12*(PA*PA*PA)*(one_over_2p*one_over_2p) + 48*(PA*PA)*(PB*PB*PB)*one_over_2p + 144*(PA*PA)*PB*(one_over_2p*one_over_2p) + 12*PA*((PB*PB)*(PB*PB))*one_over_2p + 216*PA*(PB*PB)*(one_over_2p*one_over_2p) + 180*PA*(one_over_2p*one_over_2p*one_over_2p) + 48*(PB*PB*PB)*(one_over_2p*one_over_2p) + 240*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 16*(PA*PA*PA)*(PB*PB*PB)*one_over_2p + 48*(PA*PA*PA)*PB*(one_over_2p*one_over_2p) + 12*(PA*PA)*((PB*PB)*(PB*PB))*one_over_2p + 216*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 180*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 144*PA*(PB*PB*PB)*(one_over_2p*one_over_2p) + 720*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 12*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p) + 360*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 420*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_4_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 24*(PA*PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 24*(PA*PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 48*(PA*PA)*(PB*PB*PB)*(one_over_2p*one_over_2p) + 288*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 12*PA*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p) + 432*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 540*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 96*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 720*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_4_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 16*(PA*PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 72*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 120*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 48*PA*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 480*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 4*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p*one_over_2p) + 240*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 420*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 48*(PA*PA)*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 72*PA*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 180*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 16*(PB*PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 240*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 48*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 24*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 84*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_4_4_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PA*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)) + 16*PB*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPA_4_4_7(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPA_4_4_8(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

// =============================================================================
// ∂E/∂PB Gradients
// =============================================================================

inline Vec8d hermite_dE_dPB_0_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_0_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 1;
}

inline Vec8d hermite_dE_dPB_0_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_0_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PB;
}

inline Vec8d hermite_dE_dPB_0_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*one_over_2p;
}

inline Vec8d hermite_dE_dPB_0_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_0_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PB*PB) + 3*one_over_2p;
}

inline Vec8d hermite_dE_dPB_0_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PB*one_over_2p;
}

inline Vec8d hermite_dE_dPB_0_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_0_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_0_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PB*PB*PB) + 12*PB*one_over_2p;
}

inline Vec8d hermite_dE_dPB_0_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PB*PB)*one_over_2p + 12*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_0_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_0_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_0_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_1_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_1_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_1_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PA;
}

inline Vec8d hermite_dE_dPB_1_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return one_over_2p;
}

inline Vec8d hermite_dE_dPB_1_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_1_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*PB + 2*one_over_2p;
}

inline Vec8d hermite_dE_dPB_1_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*one_over_2p + 2*PB*one_over_2p;
}

inline Vec8d hermite_dE_dPB_1_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_1_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_1_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PA*(PB*PB) + 3*PA*one_over_2p + 6*PB*one_over_2p;
}

inline Vec8d hermite_dE_dPB_1_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*PB*one_over_2p + 3*(PB*PB)*one_over_2p + 9*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_1_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PA*(one_over_2p*one_over_2p) + 6*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_1_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_1_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_1_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*(PB*PB*PB) + 12*PA*PB*one_over_2p + 12*(PB*PB)*one_over_2p + 12*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_1_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PA*(PB*PB)*one_over_2p + 12*PA*(one_over_2p*one_over_2p) + 4*(PB*PB*PB)*one_over_2p + 36*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_1_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PA*PB*(one_over_2p*one_over_2p) + 12*(PB*PB)*(one_over_2p*one_over_2p) + 24*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_1_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*(one_over_2p*one_over_2p*one_over_2p) + 12*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_1_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_1_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_2_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_2_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_2_0_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_2_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA) + one_over_2p;
}

inline Vec8d hermite_dE_dPB_2_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*one_over_2p;
}

inline Vec8d hermite_dE_dPB_2_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_1_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_2_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*(PA*PA)*PB + 4*PA*one_over_2p + 2*PB*one_over_2p;
}

inline Vec8d hermite_dE_dPB_2_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*(PA*PA)*one_over_2p + 4*PA*PB*one_over_2p + 6*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*(one_over_2p*one_over_2p) + 2*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_2_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_2_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*(PB*PB) + 3*(PA*PA)*one_over_2p + 12*PA*PB*one_over_2p + 3*(PB*PB)*one_over_2p + 9*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA)*PB*one_over_2p + 6*PA*(PB*PB)*one_over_2p + 18*PA*(one_over_2p*one_over_2p) + 18*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*(one_over_2p*one_over_2p) + 12*PA*PB*(one_over_2p*one_over_2p) + 3*(PB*PB)*(one_over_2p*one_over_2p) + 18*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*(one_over_2p*one_over_2p*one_over_2p) + 6*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_2_3_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_2_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA)*(PB*PB*PB) + 12*(PA*PA)*PB*one_over_2p + 24*PA*(PB*PB)*one_over_2p + 24*PA*(one_over_2p*one_over_2p) + 4*(PB*PB*PB)*one_over_2p + 36*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*(PB*PB)*one_over_2p + 12*(PA*PA)*(one_over_2p*one_over_2p) + 8*PA*(PB*PB*PB)*one_over_2p + 72*PA*PB*(one_over_2p*one_over_2p) + 36*(PB*PB)*(one_over_2p*one_over_2p) + 60*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*PB*(one_over_2p*one_over_2p) + 24*PA*(PB*PB)*(one_over_2p*one_over_2p) + 48*PA*(one_over_2p*one_over_2p*one_over_2p) + 4*(PB*PB*PB)*(one_over_2p*one_over_2p) + 72*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 24*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 12*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 40*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_2_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 8*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 12*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_2_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPB_2_4_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_3_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_3_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_3_0_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_3_0_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_3_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA*PA) + 3*PA*one_over_2p;
}

inline Vec8d hermite_dE_dPB_3_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PA*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_1_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_1_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_3_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*(PA*PA*PA)*PB + 6*(PA*PA)*one_over_2p + 6*PA*PB*one_over_2p + 6*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*(PA*PA*PA)*one_over_2p + 6*(PA*PA)*PB*one_over_2p + 18*PA*(one_over_2p*one_over_2p) + 6*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA)*(one_over_2p*one_over_2p) + 6*PA*PB*(one_over_2p*one_over_2p) + 12*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*(one_over_2p*one_over_2p*one_over_2p) + 2*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_2_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_3_2_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_3_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA*PA)*(PB*PB) + 3*(PA*PA*PA)*one_over_2p + 18*(PA*PA)*PB*one_over_2p + 9*PA*(PB*PB)*one_over_2p + 27*PA*(one_over_2p*one_over_2p) + 18*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA*PA)*PB*one_over_2p + 9*(PA*PA)*(PB*PB)*one_over_2p + 27*(PA*PA)*(one_over_2p*one_over_2p) + 54*PA*PB*(one_over_2p*one_over_2p) + 9*(PB*PB)*(one_over_2p*one_over_2p) + 45*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA*PA)*(one_over_2p*one_over_2p) + 18*(PA*PA)*PB*(one_over_2p*one_over_2p) + 9*PA*(PB*PB)*(one_over_2p*one_over_2p) + 54*PA*(one_over_2p*one_over_2p*one_over_2p) + 36*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 9*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 18*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 3*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 30*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_3_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 9*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 6*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_3_3_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_3_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_3_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*(PB*PB*PB) + 12*(PA*PA*PA)*PB*one_over_2p + 36*(PA*PA)*(PB*PB)*one_over_2p + 36*(PA*PA)*(one_over_2p*one_over_2p) + 12*PA*(PB*PB*PB)*one_over_2p + 108*PA*PB*(one_over_2p*one_over_2p) + 36*(PB*PB)*(one_over_2p*one_over_2p) + 60*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA*PA)*(PB*PB)*one_over_2p + 12*(PA*PA*PA)*(one_over_2p*one_over_2p) + 12*(PA*PA)*(PB*PB*PB)*one_over_2p + 108*(PA*PA)*PB*(one_over_2p*one_over_2p) + 108*PA*(PB*PB)*(one_over_2p*one_over_2p) + 180*PA*(one_over_2p*one_over_2p*one_over_2p) + 12*(PB*PB*PB)*(one_over_2p*one_over_2p) + 180*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA*PA)*PB*(one_over_2p*one_over_2p) + 36*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 72*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 12*PA*(PB*PB*PB)*(one_over_2p*one_over_2p) + 216*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 72*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 180*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_3_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 36*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 36*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 120*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 4*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 120*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_3_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 36*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 12*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 60*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 12*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPB_3_4_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_3_4_7(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_4_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_4_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_4_0_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_4_0_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_4_0_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_4_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PA*PA)*(PA*PA)) + 6*(PA*PA)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*one_over_2p + 12*PA*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA)*(one_over_2p*one_over_2p) + 6*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_1_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_1_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_4_1_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_4_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*((PA*PA)*(PA*PA))*PB + 8*(PA*PA*PA)*one_over_2p + 12*(PA*PA)*PB*one_over_2p + 24*PA*(one_over_2p*one_over_2p) + 6*PB*(one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*((PA*PA)*(PA*PA))*one_over_2p + 8*(PA*PA*PA)*PB*one_over_2p + 36*(PA*PA)*(one_over_2p*one_over_2p) + 24*PA*PB*(one_over_2p*one_over_2p) + 30*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 8*(PA*PA*PA)*(one_over_2p*one_over_2p) + 12*(PA*PA)*PB*(one_over_2p*one_over_2p) + 48*PA*(one_over_2p*one_over_2p*one_over_2p) + 12*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 8*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 20*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_4_2_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 8*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 2*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_4_2_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_2_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_4_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*((PA*PA)*(PA*PA))*(PB*PB) + 3*((PA*PA)*(PA*PA))*one_over_2p + 24*(PA*PA*PA)*PB*one_over_2p + 18*(PA*PA)*(PB*PB)*one_over_2p + 54*(PA*PA)*(one_over_2p*one_over_2p) + 72*PA*PB*(one_over_2p*one_over_2p) + 9*(PB*PB)*(one_over_2p*one_over_2p) + 45*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*((PA*PA)*(PA*PA))*PB*one_over_2p + 12*(PA*PA*PA)*(PB*PB)*one_over_2p + 36*(PA*PA*PA)*(one_over_2p*one_over_2p) + 108*(PA*PA)*PB*(one_over_2p*one_over_2p) + 36*PA*(PB*PB)*(one_over_2p*one_over_2p) + 180*PA*(one_over_2p*one_over_2p*one_over_2p) + 90*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*((PA*PA)*(PA*PA))*(one_over_2p*one_over_2p) + 24*(PA*PA*PA)*PB*(one_over_2p*one_over_2p) + 18*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 108*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 144*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 18*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 135*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_4_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*(PA*PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 36*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 12*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 120*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 60*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_4_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 18*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 24*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 3*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 45*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_3_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 6*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_3_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_4_3_7(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

inline Vec8d hermite_dE_dPB_4_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((PA*PA)*(PA*PA))*(PB*PB*PB) + 12*((PA*PA)*(PA*PA))*PB*one_over_2p + 48*(PA*PA*PA)*(PB*PB)*one_over_2p + 48*(PA*PA*PA)*(one_over_2p*one_over_2p) + 24*(PA*PA)*(PB*PB*PB)*one_over_2p + 216*(PA*PA)*PB*(one_over_2p*one_over_2p) + 144*PA*(PB*PB)*(one_over_2p*one_over_2p) + 240*PA*(one_over_2p*one_over_2p*one_over_2p) + 12*(PB*PB*PB)*(one_over_2p*one_over_2p) + 180*PB*(one_over_2p*one_over_2p*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*((PA*PA)*(PA*PA))*(PB*PB)*one_over_2p + 12*((PA*PA)*(PA*PA))*(one_over_2p*one_over_2p) + 16*(PA*PA*PA)*(PB*PB*PB)*one_over_2p + 144*(PA*PA*PA)*PB*(one_over_2p*one_over_2p) + 216*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 360*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 48*PA*(PB*PB*PB)*(one_over_2p*one_over_2p) + 720*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 180*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 420*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_4_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 12*((PA*PA)*(PA*PA))*PB*(one_over_2p*one_over_2p) + 48*(PA*PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 96*(PA*PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 24*(PA*PA)*(PB*PB*PB)*(one_over_2p*one_over_2p) + 432*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 288*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 720*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 24*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 540*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_4_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((PA*PA)*(PA*PA))*(one_over_2p*one_over_2p*one_over_2p) + 48*(PA*PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 72*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 240*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 16*PA*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 480*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 120*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 420*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 16*(PA*PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 72*(PA*PA)*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 48*PA*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 240*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 4*(PB*PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 180*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 24*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 48*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 12*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 84*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_4_4_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 16*PA*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)) + 12*PB*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

inline Vec8d hermite_dE_dPB_4_4_7(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)*one_over_2p);
}

inline Vec8d hermite_dE_dPB_4_4_8(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 0;
}

/**
 * @brief Runtime dispatcher for ∂E/∂PA
 */
inline Vec8d dispatch_dE_dPA(int nA, int nB, int t,
                              Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    if (nA == 0 && nB == 0 && t == 0) return hermite_dE_dPA_0_0_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 1 && t == 0) return hermite_dE_dPA_0_1_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 1 && t == 1) return hermite_dE_dPA_0_1_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 2 && t == 0) return hermite_dE_dPA_0_2_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 2 && t == 1) return hermite_dE_dPA_0_2_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 2 && t == 2) return hermite_dE_dPA_0_2_2(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 0) return hermite_dE_dPA_0_3_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 1) return hermite_dE_dPA_0_3_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 2) return hermite_dE_dPA_0_3_2(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 3) return hermite_dE_dPA_0_3_3(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 0) return hermite_dE_dPA_0_4_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 1) return hermite_dE_dPA_0_4_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 2) return hermite_dE_dPA_0_4_2(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 3) return hermite_dE_dPA_0_4_3(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 4) return hermite_dE_dPA_0_4_4(PA, PB, one_over_2p);
    if (nA == 1 && nB == 0 && t == 0) return hermite_dE_dPA_1_0_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 0 && t == 1) return hermite_dE_dPA_1_0_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 1 && t == 0) return hermite_dE_dPA_1_1_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 1 && t == 1) return hermite_dE_dPA_1_1_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 1 && t == 2) return hermite_dE_dPA_1_1_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 0) return hermite_dE_dPA_1_2_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 1) return hermite_dE_dPA_1_2_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 2) return hermite_dE_dPA_1_2_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 3) return hermite_dE_dPA_1_2_3(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 0) return hermite_dE_dPA_1_3_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 1) return hermite_dE_dPA_1_3_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 2) return hermite_dE_dPA_1_3_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 3) return hermite_dE_dPA_1_3_3(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 4) return hermite_dE_dPA_1_3_4(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 0) return hermite_dE_dPA_1_4_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 1) return hermite_dE_dPA_1_4_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 2) return hermite_dE_dPA_1_4_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 3) return hermite_dE_dPA_1_4_3(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 4) return hermite_dE_dPA_1_4_4(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 5) return hermite_dE_dPA_1_4_5(PA, PB, one_over_2p);
    if (nA == 2 && nB == 0 && t == 0) return hermite_dE_dPA_2_0_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 0 && t == 1) return hermite_dE_dPA_2_0_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 0 && t == 2) return hermite_dE_dPA_2_0_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 0) return hermite_dE_dPA_2_1_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 1) return hermite_dE_dPA_2_1_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 2) return hermite_dE_dPA_2_1_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 3) return hermite_dE_dPA_2_1_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 0) return hermite_dE_dPA_2_2_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 1) return hermite_dE_dPA_2_2_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 2) return hermite_dE_dPA_2_2_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 3) return hermite_dE_dPA_2_2_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 4) return hermite_dE_dPA_2_2_4(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 0) return hermite_dE_dPA_2_3_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 1) return hermite_dE_dPA_2_3_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 2) return hermite_dE_dPA_2_3_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 3) return hermite_dE_dPA_2_3_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 4) return hermite_dE_dPA_2_3_4(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 5) return hermite_dE_dPA_2_3_5(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 0) return hermite_dE_dPA_2_4_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 1) return hermite_dE_dPA_2_4_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 2) return hermite_dE_dPA_2_4_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 3) return hermite_dE_dPA_2_4_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 4) return hermite_dE_dPA_2_4_4(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 5) return hermite_dE_dPA_2_4_5(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 6) return hermite_dE_dPA_2_4_6(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 0) return hermite_dE_dPA_3_0_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 1) return hermite_dE_dPA_3_0_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 2) return hermite_dE_dPA_3_0_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 3) return hermite_dE_dPA_3_0_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 0) return hermite_dE_dPA_3_1_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 1) return hermite_dE_dPA_3_1_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 2) return hermite_dE_dPA_3_1_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 3) return hermite_dE_dPA_3_1_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 4) return hermite_dE_dPA_3_1_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 0) return hermite_dE_dPA_3_2_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 1) return hermite_dE_dPA_3_2_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 2) return hermite_dE_dPA_3_2_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 3) return hermite_dE_dPA_3_2_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 4) return hermite_dE_dPA_3_2_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 5) return hermite_dE_dPA_3_2_5(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 0) return hermite_dE_dPA_3_3_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 1) return hermite_dE_dPA_3_3_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 2) return hermite_dE_dPA_3_3_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 3) return hermite_dE_dPA_3_3_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 4) return hermite_dE_dPA_3_3_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 5) return hermite_dE_dPA_3_3_5(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 6) return hermite_dE_dPA_3_3_6(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 0) return hermite_dE_dPA_3_4_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 1) return hermite_dE_dPA_3_4_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 2) return hermite_dE_dPA_3_4_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 3) return hermite_dE_dPA_3_4_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 4) return hermite_dE_dPA_3_4_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 5) return hermite_dE_dPA_3_4_5(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 6) return hermite_dE_dPA_3_4_6(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 7) return hermite_dE_dPA_3_4_7(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 0) return hermite_dE_dPA_4_0_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 1) return hermite_dE_dPA_4_0_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 2) return hermite_dE_dPA_4_0_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 3) return hermite_dE_dPA_4_0_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 4) return hermite_dE_dPA_4_0_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 0) return hermite_dE_dPA_4_1_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 1) return hermite_dE_dPA_4_1_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 2) return hermite_dE_dPA_4_1_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 3) return hermite_dE_dPA_4_1_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 4) return hermite_dE_dPA_4_1_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 5) return hermite_dE_dPA_4_1_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 0) return hermite_dE_dPA_4_2_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 1) return hermite_dE_dPA_4_2_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 2) return hermite_dE_dPA_4_2_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 3) return hermite_dE_dPA_4_2_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 4) return hermite_dE_dPA_4_2_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 5) return hermite_dE_dPA_4_2_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 6) return hermite_dE_dPA_4_2_6(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 0) return hermite_dE_dPA_4_3_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 1) return hermite_dE_dPA_4_3_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 2) return hermite_dE_dPA_4_3_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 3) return hermite_dE_dPA_4_3_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 4) return hermite_dE_dPA_4_3_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 5) return hermite_dE_dPA_4_3_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 6) return hermite_dE_dPA_4_3_6(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 7) return hermite_dE_dPA_4_3_7(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 0) return hermite_dE_dPA_4_4_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 1) return hermite_dE_dPA_4_4_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 2) return hermite_dE_dPA_4_4_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 3) return hermite_dE_dPA_4_4_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 4) return hermite_dE_dPA_4_4_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 5) return hermite_dE_dPA_4_4_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 6) return hermite_dE_dPA_4_4_6(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 7) return hermite_dE_dPA_4_4_7(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 8) return hermite_dE_dPA_4_4_8(PA, PB, one_over_2p);
    return Vec8d(0.0);
}

/**
 * @brief Runtime dispatcher for ∂E/∂PB
 */
inline Vec8d dispatch_dE_dPB(int nA, int nB, int t,
                              Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    if (nA == 0 && nB == 0 && t == 0) return hermite_dE_dPB_0_0_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 1 && t == 0) return hermite_dE_dPB_0_1_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 1 && t == 1) return hermite_dE_dPB_0_1_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 2 && t == 0) return hermite_dE_dPB_0_2_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 2 && t == 1) return hermite_dE_dPB_0_2_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 2 && t == 2) return hermite_dE_dPB_0_2_2(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 0) return hermite_dE_dPB_0_3_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 1) return hermite_dE_dPB_0_3_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 2) return hermite_dE_dPB_0_3_2(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 3) return hermite_dE_dPB_0_3_3(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 0) return hermite_dE_dPB_0_4_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 1) return hermite_dE_dPB_0_4_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 2) return hermite_dE_dPB_0_4_2(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 3) return hermite_dE_dPB_0_4_3(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 4) return hermite_dE_dPB_0_4_4(PA, PB, one_over_2p);
    if (nA == 1 && nB == 0 && t == 0) return hermite_dE_dPB_1_0_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 0 && t == 1) return hermite_dE_dPB_1_0_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 1 && t == 0) return hermite_dE_dPB_1_1_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 1 && t == 1) return hermite_dE_dPB_1_1_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 1 && t == 2) return hermite_dE_dPB_1_1_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 0) return hermite_dE_dPB_1_2_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 1) return hermite_dE_dPB_1_2_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 2) return hermite_dE_dPB_1_2_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 3) return hermite_dE_dPB_1_2_3(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 0) return hermite_dE_dPB_1_3_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 1) return hermite_dE_dPB_1_3_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 2) return hermite_dE_dPB_1_3_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 3) return hermite_dE_dPB_1_3_3(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 4) return hermite_dE_dPB_1_3_4(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 0) return hermite_dE_dPB_1_4_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 1) return hermite_dE_dPB_1_4_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 2) return hermite_dE_dPB_1_4_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 3) return hermite_dE_dPB_1_4_3(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 4) return hermite_dE_dPB_1_4_4(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 5) return hermite_dE_dPB_1_4_5(PA, PB, one_over_2p);
    if (nA == 2 && nB == 0 && t == 0) return hermite_dE_dPB_2_0_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 0 && t == 1) return hermite_dE_dPB_2_0_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 0 && t == 2) return hermite_dE_dPB_2_0_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 0) return hermite_dE_dPB_2_1_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 1) return hermite_dE_dPB_2_1_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 2) return hermite_dE_dPB_2_1_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 3) return hermite_dE_dPB_2_1_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 0) return hermite_dE_dPB_2_2_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 1) return hermite_dE_dPB_2_2_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 2) return hermite_dE_dPB_2_2_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 3) return hermite_dE_dPB_2_2_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 4) return hermite_dE_dPB_2_2_4(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 0) return hermite_dE_dPB_2_3_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 1) return hermite_dE_dPB_2_3_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 2) return hermite_dE_dPB_2_3_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 3) return hermite_dE_dPB_2_3_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 4) return hermite_dE_dPB_2_3_4(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 5) return hermite_dE_dPB_2_3_5(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 0) return hermite_dE_dPB_2_4_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 1) return hermite_dE_dPB_2_4_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 2) return hermite_dE_dPB_2_4_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 3) return hermite_dE_dPB_2_4_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 4) return hermite_dE_dPB_2_4_4(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 5) return hermite_dE_dPB_2_4_5(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 6) return hermite_dE_dPB_2_4_6(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 0) return hermite_dE_dPB_3_0_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 1) return hermite_dE_dPB_3_0_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 2) return hermite_dE_dPB_3_0_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 3) return hermite_dE_dPB_3_0_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 0) return hermite_dE_dPB_3_1_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 1) return hermite_dE_dPB_3_1_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 2) return hermite_dE_dPB_3_1_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 3) return hermite_dE_dPB_3_1_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 4) return hermite_dE_dPB_3_1_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 0) return hermite_dE_dPB_3_2_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 1) return hermite_dE_dPB_3_2_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 2) return hermite_dE_dPB_3_2_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 3) return hermite_dE_dPB_3_2_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 4) return hermite_dE_dPB_3_2_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 5) return hermite_dE_dPB_3_2_5(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 0) return hermite_dE_dPB_3_3_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 1) return hermite_dE_dPB_3_3_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 2) return hermite_dE_dPB_3_3_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 3) return hermite_dE_dPB_3_3_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 4) return hermite_dE_dPB_3_3_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 5) return hermite_dE_dPB_3_3_5(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 6) return hermite_dE_dPB_3_3_6(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 0) return hermite_dE_dPB_3_4_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 1) return hermite_dE_dPB_3_4_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 2) return hermite_dE_dPB_3_4_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 3) return hermite_dE_dPB_3_4_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 4) return hermite_dE_dPB_3_4_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 5) return hermite_dE_dPB_3_4_5(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 6) return hermite_dE_dPB_3_4_6(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 7) return hermite_dE_dPB_3_4_7(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 0) return hermite_dE_dPB_4_0_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 1) return hermite_dE_dPB_4_0_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 2) return hermite_dE_dPB_4_0_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 3) return hermite_dE_dPB_4_0_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 4) return hermite_dE_dPB_4_0_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 0) return hermite_dE_dPB_4_1_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 1) return hermite_dE_dPB_4_1_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 2) return hermite_dE_dPB_4_1_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 3) return hermite_dE_dPB_4_1_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 4) return hermite_dE_dPB_4_1_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 5) return hermite_dE_dPB_4_1_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 0) return hermite_dE_dPB_4_2_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 1) return hermite_dE_dPB_4_2_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 2) return hermite_dE_dPB_4_2_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 3) return hermite_dE_dPB_4_2_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 4) return hermite_dE_dPB_4_2_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 5) return hermite_dE_dPB_4_2_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 6) return hermite_dE_dPB_4_2_6(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 0) return hermite_dE_dPB_4_3_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 1) return hermite_dE_dPB_4_3_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 2) return hermite_dE_dPB_4_3_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 3) return hermite_dE_dPB_4_3_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 4) return hermite_dE_dPB_4_3_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 5) return hermite_dE_dPB_4_3_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 6) return hermite_dE_dPB_4_3_6(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 7) return hermite_dE_dPB_4_3_7(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 0) return hermite_dE_dPB_4_4_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 1) return hermite_dE_dPB_4_4_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 2) return hermite_dE_dPB_4_4_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 3) return hermite_dE_dPB_4_4_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 4) return hermite_dE_dPB_4_4_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 5) return hermite_dE_dPB_4_4_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 6) return hermite_dE_dPB_4_4_6(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 7) return hermite_dE_dPB_4_4_7(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 8) return hermite_dE_dPB_4_4_8(PA, PB, one_over_2p);
    return Vec8d(0.0);
}

} // namespace symbolic
} // namespace recursum
