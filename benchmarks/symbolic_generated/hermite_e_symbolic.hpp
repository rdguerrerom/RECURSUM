/**
 * @file hermite_e_symbolic.hpp
 * @brief Symbolically-generated Hermite E coefficients (expanded form)
 *
 * AUTO-GENERATED - DO NOT EDIT MANUALLY
 * Coverage: E^{0,0}_0 to E^{4,4}_8 (ss to gg shell pairs)
 */

#pragma once

#ifndef RECURSUM_VEC_TYPE
#include <recursum/vectorclass.h>
#endif

namespace recursum {
namespace symbolic {

/**
 * @brief Symbolic E^{0,0}_0
 */
inline Vec8d hermite_e_symbolic_0_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 1;
}

/**
 * @brief Symbolic E^{0,1}_0
 */
inline Vec8d hermite_e_symbolic_0_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PB;
}

/**
 * @brief Symbolic E^{0,1}_1
 */
inline Vec8d hermite_e_symbolic_0_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return one_over_2p;
}

/**
 * @brief Symbolic E^{0,2}_0
 */
inline Vec8d hermite_e_symbolic_0_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PB*PB) + one_over_2p;
}

/**
 * @brief Symbolic E^{0,2}_1
 */
inline Vec8d hermite_e_symbolic_0_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PB*one_over_2p;
}

/**
 * @brief Symbolic E^{0,2}_2
 */
inline Vec8d hermite_e_symbolic_0_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{0,3}_0
 */
inline Vec8d hermite_e_symbolic_0_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PB*PB*PB) + 3*PB*one_over_2p;
}

/**
 * @brief Symbolic E^{0,3}_1
 */
inline Vec8d hermite_e_symbolic_0_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PB*PB)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{0,3}_2
 */
inline Vec8d hermite_e_symbolic_0_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{0,3}_3
 */
inline Vec8d hermite_e_symbolic_0_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{0,4}_0
 */
inline Vec8d hermite_e_symbolic_0_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PB*PB)*(PB*PB)) + 6*(PB*PB)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{0,4}_1
 */
inline Vec8d hermite_e_symbolic_0_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PB*PB*PB)*one_over_2p + 12*PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{0,4}_2
 */
inline Vec8d hermite_e_symbolic_0_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PB*PB)*(one_over_2p*one_over_2p) + 6*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{0,4}_3
 */
inline Vec8d hermite_e_symbolic_0_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{0,4}_4
 */
inline Vec8d hermite_e_symbolic_0_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{1,0}_0
 */
inline Vec8d hermite_e_symbolic_1_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PA;
}

/**
 * @brief Symbolic E^{1,0}_1
 */
inline Vec8d hermite_e_symbolic_1_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return one_over_2p;
}

/**
 * @brief Symbolic E^{1,1}_0
 */
inline Vec8d hermite_e_symbolic_1_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PA*PB + one_over_2p;
}

/**
 * @brief Symbolic E^{1,1}_1
 */
inline Vec8d hermite_e_symbolic_1_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PA*one_over_2p + PB*one_over_2p;
}

/**
 * @brief Symbolic E^{1,1}_2
 */
inline Vec8d hermite_e_symbolic_1_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,2}_0
 */
inline Vec8d hermite_e_symbolic_1_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PA*(PB*PB) + PA*one_over_2p + 2*PB*one_over_2p;
}

/**
 * @brief Symbolic E^{1,2}_1
 */
inline Vec8d hermite_e_symbolic_1_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*PB*one_over_2p + (PB*PB)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,2}_2
 */
inline Vec8d hermite_e_symbolic_1_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PA*(one_over_2p*one_over_2p) + 2*PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,2}_3
 */
inline Vec8d hermite_e_symbolic_1_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,3}_0
 */
inline Vec8d hermite_e_symbolic_1_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PA*(PB*PB*PB) + 3*PA*PB*one_over_2p + 3*(PB*PB)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,3}_1
 */
inline Vec8d hermite_e_symbolic_1_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PA*(PB*PB)*one_over_2p + 3*PA*(one_over_2p*one_over_2p) + (PB*PB*PB)*one_over_2p + 9*PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,3}_2
 */
inline Vec8d hermite_e_symbolic_1_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PA*PB*(one_over_2p*one_over_2p) + 3*(PB*PB)*(one_over_2p*one_over_2p) + 6*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,3}_3
 */
inline Vec8d hermite_e_symbolic_1_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PA*(one_over_2p*one_over_2p*one_over_2p) + 3*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,3}_4
 */
inline Vec8d hermite_e_symbolic_1_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{1,4}_0
 */
inline Vec8d hermite_e_symbolic_1_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PA*((PB*PB)*(PB*PB)) + 6*PA*(PB*PB)*one_over_2p + 3*PA*(one_over_2p*one_over_2p) + 4*(PB*PB*PB)*one_over_2p + 12*PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,4}_1
 */
inline Vec8d hermite_e_symbolic_1_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*(PB*PB*PB)*one_over_2p + 12*PA*PB*(one_over_2p*one_over_2p) + ((PB*PB)*(PB*PB))*one_over_2p + 18*(PB*PB)*(one_over_2p*one_over_2p) + 15*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,4}_2
 */
inline Vec8d hermite_e_symbolic_1_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*PA*(PB*PB)*(one_over_2p*one_over_2p) + 6*PA*(one_over_2p*one_over_2p*one_over_2p) + 4*(PB*PB*PB)*(one_over_2p*one_over_2p) + 24*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{1,4}_3
 */
inline Vec8d hermite_e_symbolic_1_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 6*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 10*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{1,4}_4
 */
inline Vec8d hermite_e_symbolic_1_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 4*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{1,4}_5
 */
inline Vec8d hermite_e_symbolic_1_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{2,0}_0
 */
inline Vec8d hermite_e_symbolic_2_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA) + one_over_2p;
}

/**
 * @brief Symbolic E^{2,0}_1
 */
inline Vec8d hermite_e_symbolic_2_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*one_over_2p;
}

/**
 * @brief Symbolic E^{2,0}_2
 */
inline Vec8d hermite_e_symbolic_2_0_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,1}_0
 */
inline Vec8d hermite_e_symbolic_2_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA)*PB + 2*PA*one_over_2p + PB*one_over_2p;
}

/**
 * @brief Symbolic E^{2,1}_1
 */
inline Vec8d hermite_e_symbolic_2_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA)*one_over_2p + 2*PA*PB*one_over_2p + 3*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,1}_2
 */
inline Vec8d hermite_e_symbolic_2_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*(one_over_2p*one_over_2p) + PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,1}_3
 */
inline Vec8d hermite_e_symbolic_2_1_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,2}_0
 */
inline Vec8d hermite_e_symbolic_2_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA)*(PB*PB) + (PA*PA)*one_over_2p + 4*PA*PB*one_over_2p + (PB*PB)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,2}_1
 */
inline Vec8d hermite_e_symbolic_2_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*(PA*PA)*PB*one_over_2p + 2*PA*(PB*PB)*one_over_2p + 6*PA*(one_over_2p*one_over_2p) + 6*PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,2}_2
 */
inline Vec8d hermite_e_symbolic_2_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA)*(one_over_2p*one_over_2p) + 4*PA*PB*(one_over_2p*one_over_2p) + (PB*PB)*(one_over_2p*one_over_2p) + 6*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,2}_3
 */
inline Vec8d hermite_e_symbolic_2_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*(one_over_2p*one_over_2p*one_over_2p) + 2*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,2}_4
 */
inline Vec8d hermite_e_symbolic_2_2_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{2,3}_0
 */
inline Vec8d hermite_e_symbolic_2_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA)*(PB*PB*PB) + 3*(PA*PA)*PB*one_over_2p + 6*PA*(PB*PB)*one_over_2p + 6*PA*(one_over_2p*one_over_2p) + (PB*PB*PB)*one_over_2p + 9*PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,3}_1
 */
inline Vec8d hermite_e_symbolic_2_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*(PB*PB)*one_over_2p + 3*(PA*PA)*(one_over_2p*one_over_2p) + 2*PA*(PB*PB*PB)*one_over_2p + 18*PA*PB*(one_over_2p*one_over_2p) + 9*(PB*PB)*(one_over_2p*one_over_2p) + 15*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,3}_2
 */
inline Vec8d hermite_e_symbolic_2_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*PB*(one_over_2p*one_over_2p) + 6*PA*(PB*PB)*(one_over_2p*one_over_2p) + 12*PA*(one_over_2p*one_over_2p*one_over_2p) + (PB*PB*PB)*(one_over_2p*one_over_2p) + 18*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,3}_3
 */
inline Vec8d hermite_e_symbolic_2_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 6*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 3*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 10*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{2,3}_4
 */
inline Vec8d hermite_e_symbolic_2_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 3*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{2,3}_5
 */
inline Vec8d hermite_e_symbolic_2_3_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{2,4}_0
 */
inline Vec8d hermite_e_symbolic_2_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA)*((PB*PB)*(PB*PB)) + 6*(PA*PA)*(PB*PB)*one_over_2p + 3*(PA*PA)*(one_over_2p*one_over_2p) + 8*PA*(PB*PB*PB)*one_over_2p + 24*PA*PB*(one_over_2p*one_over_2p) + ((PB*PB)*(PB*PB))*one_over_2p + 18*(PB*PB)*(one_over_2p*one_over_2p) + 15*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,4}_1
 */
inline Vec8d hermite_e_symbolic_2_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA)*(PB*PB*PB)*one_over_2p + 12*(PA*PA)*PB*(one_over_2p*one_over_2p) + 2*PA*((PB*PB)*(PB*PB))*one_over_2p + 36*PA*(PB*PB)*(one_over_2p*one_over_2p) + 30*PA*(one_over_2p*one_over_2p*one_over_2p) + 12*(PB*PB*PB)*(one_over_2p*one_over_2p) + 60*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{2,4}_2
 */
inline Vec8d hermite_e_symbolic_2_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 6*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 8*PA*(PB*PB*PB)*(one_over_2p*one_over_2p) + 48*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + ((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p) + 36*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 45*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{2,4}_3
 */
inline Vec8d hermite_e_symbolic_2_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 12*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 20*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 4*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 40*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{2,4}_4
 */
inline Vec8d hermite_e_symbolic_2_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 8*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 6*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 15*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{2,4}_5
 */
inline Vec8d hermite_e_symbolic_2_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 4*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{2,4}_6
 */
inline Vec8d hermite_e_symbolic_2_4_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,0}_0
 */
inline Vec8d hermite_e_symbolic_3_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA*PA) + 3*PA*one_over_2p;
}

/**
 * @brief Symbolic E^{3,0}_1
 */
inline Vec8d hermite_e_symbolic_3_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,0}_2
 */
inline Vec8d hermite_e_symbolic_3_0_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PA*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,0}_3
 */
inline Vec8d hermite_e_symbolic_3_0_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,1}_0
 */
inline Vec8d hermite_e_symbolic_3_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA*PA)*PB + 3*(PA*PA)*one_over_2p + 3*PA*PB*one_over_2p + 3*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,1}_1
 */
inline Vec8d hermite_e_symbolic_3_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA*PA)*one_over_2p + 3*(PA*PA)*PB*one_over_2p + 9*PA*(one_over_2p*one_over_2p) + 3*PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,1}_2
 */
inline Vec8d hermite_e_symbolic_3_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*(one_over_2p*one_over_2p) + 3*PA*PB*(one_over_2p*one_over_2p) + 6*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,1}_3
 */
inline Vec8d hermite_e_symbolic_3_1_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PA*(one_over_2p*one_over_2p*one_over_2p) + PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,1}_4
 */
inline Vec8d hermite_e_symbolic_3_1_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,2}_0
 */
inline Vec8d hermite_e_symbolic_3_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA*PA)*(PB*PB) + (PA*PA*PA)*one_over_2p + 6*(PA*PA)*PB*one_over_2p + 3*PA*(PB*PB)*one_over_2p + 9*PA*(one_over_2p*one_over_2p) + 6*PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,2}_1
 */
inline Vec8d hermite_e_symbolic_3_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*(PA*PA*PA)*PB*one_over_2p + 3*(PA*PA)*(PB*PB)*one_over_2p + 9*(PA*PA)*(one_over_2p*one_over_2p) + 18*PA*PB*(one_over_2p*one_over_2p) + 3*(PB*PB)*(one_over_2p*one_over_2p) + 15*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,2}_2
 */
inline Vec8d hermite_e_symbolic_3_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA*PA)*(one_over_2p*one_over_2p) + 6*(PA*PA)*PB*(one_over_2p*one_over_2p) + 3*PA*(PB*PB)*(one_over_2p*one_over_2p) + 18*PA*(one_over_2p*one_over_2p*one_over_2p) + 12*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,2}_3
 */
inline Vec8d hermite_e_symbolic_3_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 6*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + (PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 10*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,2}_4
 */
inline Vec8d hermite_e_symbolic_3_2_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 2*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,2}_5
 */
inline Vec8d hermite_e_symbolic_3_2_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{3,3}_0
 */
inline Vec8d hermite_e_symbolic_3_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA*PA)*(PB*PB*PB) + 3*(PA*PA*PA)*PB*one_over_2p + 9*(PA*PA)*(PB*PB)*one_over_2p + 9*(PA*PA)*(one_over_2p*one_over_2p) + 3*PA*(PB*PB*PB)*one_over_2p + 27*PA*PB*(one_over_2p*one_over_2p) + 9*(PB*PB)*(one_over_2p*one_over_2p) + 15*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,3}_1
 */
inline Vec8d hermite_e_symbolic_3_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA*PA)*(PB*PB)*one_over_2p + 3*(PA*PA*PA)*(one_over_2p*one_over_2p) + 3*(PA*PA)*(PB*PB*PB)*one_over_2p + 27*(PA*PA)*PB*(one_over_2p*one_over_2p) + 27*PA*(PB*PB)*(one_over_2p*one_over_2p) + 45*PA*(one_over_2p*one_over_2p*one_over_2p) + 3*(PB*PB*PB)*(one_over_2p*one_over_2p) + 45*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,3}_2
 */
inline Vec8d hermite_e_symbolic_3_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA*PA)*PB*(one_over_2p*one_over_2p) + 9*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 18*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 3*PA*(PB*PB*PB)*(one_over_2p*one_over_2p) + 54*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 18*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 45*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,3}_3
 */
inline Vec8d hermite_e_symbolic_3_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 9*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 9*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 30*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + (PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 30*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,3}_4
 */
inline Vec8d hermite_e_symbolic_3_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 9*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 3*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 15*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{3,3}_5
 */
inline Vec8d hermite_e_symbolic_3_3_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 3*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{3,3}_6
 */
inline Vec8d hermite_e_symbolic_3_3_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,4}_0
 */
inline Vec8d hermite_e_symbolic_3_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA*PA)*((PB*PB)*(PB*PB)) + 6*(PA*PA*PA)*(PB*PB)*one_over_2p + 3*(PA*PA*PA)*(one_over_2p*one_over_2p) + 12*(PA*PA)*(PB*PB*PB)*one_over_2p + 36*(PA*PA)*PB*(one_over_2p*one_over_2p) + 3*PA*((PB*PB)*(PB*PB))*one_over_2p + 54*PA*(PB*PB)*(one_over_2p*one_over_2p) + 45*PA*(one_over_2p*one_over_2p*one_over_2p) + 12*(PB*PB*PB)*(one_over_2p*one_over_2p) + 60*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{3,4}_1
 */
inline Vec8d hermite_e_symbolic_3_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*(PB*PB*PB)*one_over_2p + 12*(PA*PA*PA)*PB*(one_over_2p*one_over_2p) + 3*(PA*PA)*((PB*PB)*(PB*PB))*one_over_2p + 54*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 45*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 36*PA*(PB*PB*PB)*(one_over_2p*one_over_2p) + 180*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 3*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p) + 90*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 105*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,4}_2
 */
inline Vec8d hermite_e_symbolic_3_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 6*(PA*PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 12*(PA*PA)*(PB*PB*PB)*(one_over_2p*one_over_2p) + 72*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 3*PA*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p) + 108*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 135*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 24*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 180*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,4}_3
 */
inline Vec8d hermite_e_symbolic_3_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 18*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 30*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 12*PA*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 120*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + ((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p*one_over_2p) + 60*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 105*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{3,4}_4
 */
inline Vec8d hermite_e_symbolic_3_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (PA*PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 12*(PA*PA)*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 18*PA*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 45*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 4*(PB*PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 60*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{3,4}_5
 */
inline Vec8d hermite_e_symbolic_3_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 12*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 6*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 21*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,4}_6
 */
inline Vec8d hermite_e_symbolic_3_4_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*PA*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)) + 4*PB*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{3,4}_7
 */
inline Vec8d hermite_e_symbolic_3_4_7(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,0}_0
 */
inline Vec8d hermite_e_symbolic_4_0_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PA*PA)*(PA*PA)) + 6*(PA*PA)*one_over_2p + 3*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{4,0}_1
 */
inline Vec8d hermite_e_symbolic_4_0_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*one_over_2p + 12*PA*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{4,0}_2
 */
inline Vec8d hermite_e_symbolic_4_0_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA)*(one_over_2p*one_over_2p) + 6*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{4,0}_3
 */
inline Vec8d hermite_e_symbolic_4_0_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{4,0}_4
 */
inline Vec8d hermite_e_symbolic_4_0_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,1}_0
 */
inline Vec8d hermite_e_symbolic_4_1_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PA*PA)*(PA*PA))*PB + 4*(PA*PA*PA)*one_over_2p + 6*(PA*PA)*PB*one_over_2p + 12*PA*(one_over_2p*one_over_2p) + 3*PB*(one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{4,1}_1
 */
inline Vec8d hermite_e_symbolic_4_1_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PA*PA)*(PA*PA))*one_over_2p + 4*(PA*PA*PA)*PB*one_over_2p + 18*(PA*PA)*(one_over_2p*one_over_2p) + 12*PA*PB*(one_over_2p*one_over_2p) + 15*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{4,1}_2
 */
inline Vec8d hermite_e_symbolic_4_1_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*(one_over_2p*one_over_2p) + 6*(PA*PA)*PB*(one_over_2p*one_over_2p) + 24*PA*(one_over_2p*one_over_2p*one_over_2p) + 6*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{4,1}_3
 */
inline Vec8d hermite_e_symbolic_4_1_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 4*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 10*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,1}_4
 */
inline Vec8d hermite_e_symbolic_4_1_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,1}_5
 */
inline Vec8d hermite_e_symbolic_4_1_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,2}_0
 */
inline Vec8d hermite_e_symbolic_4_2_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PA*PA)*(PA*PA))*(PB*PB) + ((PA*PA)*(PA*PA))*one_over_2p + 8*(PA*PA*PA)*PB*one_over_2p + 6*(PA*PA)*(PB*PB)*one_over_2p + 18*(PA*PA)*(one_over_2p*one_over_2p) + 24*PA*PB*(one_over_2p*one_over_2p) + 3*(PB*PB)*(one_over_2p*one_over_2p) + 15*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{4,2}_1
 */
inline Vec8d hermite_e_symbolic_4_2_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 2*((PA*PA)*(PA*PA))*PB*one_over_2p + 4*(PA*PA*PA)*(PB*PB)*one_over_2p + 12*(PA*PA*PA)*(one_over_2p*one_over_2p) + 36*(PA*PA)*PB*(one_over_2p*one_over_2p) + 12*PA*(PB*PB)*(one_over_2p*one_over_2p) + 60*PA*(one_over_2p*one_over_2p*one_over_2p) + 30*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{4,2}_2
 */
inline Vec8d hermite_e_symbolic_4_2_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PA*PA)*(PA*PA))*(one_over_2p*one_over_2p) + 8*(PA*PA*PA)*PB*(one_over_2p*one_over_2p) + 6*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 36*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 48*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 6*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 45*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,2}_3
 */
inline Vec8d hermite_e_symbolic_4_2_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 12*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 4*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 40*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 20*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,2}_4
 */
inline Vec8d hermite_e_symbolic_4_2_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 8*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + (PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 15*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,2}_5
 */
inline Vec8d hermite_e_symbolic_4_2_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 2*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,2}_6
 */
inline Vec8d hermite_e_symbolic_4_2_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,3}_0
 */
inline Vec8d hermite_e_symbolic_4_3_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PA*PA)*(PA*PA))*(PB*PB*PB) + 3*((PA*PA)*(PA*PA))*PB*one_over_2p + 12*(PA*PA*PA)*(PB*PB)*one_over_2p + 12*(PA*PA*PA)*(one_over_2p*one_over_2p) + 6*(PA*PA)*(PB*PB*PB)*one_over_2p + 54*(PA*PA)*PB*(one_over_2p*one_over_2p) + 36*PA*(PB*PB)*(one_over_2p*one_over_2p) + 60*PA*(one_over_2p*one_over_2p*one_over_2p) + 3*(PB*PB*PB)*(one_over_2p*one_over_2p) + 45*PB*(one_over_2p*one_over_2p*one_over_2p);
}

/**
 * @brief Symbolic E^{4,3}_1
 */
inline Vec8d hermite_e_symbolic_4_3_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*((PA*PA)*(PA*PA))*(PB*PB)*one_over_2p + 3*((PA*PA)*(PA*PA))*(one_over_2p*one_over_2p) + 4*(PA*PA*PA)*(PB*PB*PB)*one_over_2p + 36*(PA*PA*PA)*PB*(one_over_2p*one_over_2p) + 54*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 90*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 12*PA*(PB*PB*PB)*(one_over_2p*one_over_2p) + 180*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 45*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 105*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,3}_2
 */
inline Vec8d hermite_e_symbolic_4_3_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 3*((PA*PA)*(PA*PA))*PB*(one_over_2p*one_over_2p) + 12*(PA*PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 24*(PA*PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 6*(PA*PA)*(PB*PB*PB)*(one_over_2p*one_over_2p) + 108*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 72*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 180*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 6*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 135*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,3}_3
 */
inline Vec8d hermite_e_symbolic_4_3_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PA*PA)*(PA*PA))*(one_over_2p*one_over_2p*one_over_2p) + 12*(PA*PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 18*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 60*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 4*PA*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 120*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 30*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 105*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,3}_4
 */
inline Vec8d hermite_e_symbolic_4_3_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 18*(PA*PA)*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 12*PA*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 60*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + (PB*PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 45*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,3}_5
 */
inline Vec8d hermite_e_symbolic_4_3_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 12*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 3*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 21*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,3}_6
 */
inline Vec8d hermite_e_symbolic_4_3_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)) + 3*PB*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,3}_7
 */
inline Vec8d hermite_e_symbolic_4_3_7(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,4}_0
 */
inline Vec8d hermite_e_symbolic_4_4_0(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PA*PA)*(PA*PA))*((PB*PB)*(PB*PB)) + 6*((PA*PA)*(PA*PA))*(PB*PB)*one_over_2p + 3*((PA*PA)*(PA*PA))*(one_over_2p*one_over_2p) + 16*(PA*PA*PA)*(PB*PB*PB)*one_over_2p + 48*(PA*PA*PA)*PB*(one_over_2p*one_over_2p) + 6*(PA*PA)*((PB*PB)*(PB*PB))*one_over_2p + 108*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 90*(PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 48*PA*(PB*PB*PB)*(one_over_2p*one_over_2p) + 240*PA*PB*(one_over_2p*one_over_2p*one_over_2p) + 3*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p) + 90*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 105*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,4}_1
 */
inline Vec8d hermite_e_symbolic_4_4_1(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((PA*PA)*(PA*PA))*(PB*PB*PB)*one_over_2p + 12*((PA*PA)*(PA*PA))*PB*(one_over_2p*one_over_2p) + 4*(PA*PA*PA)*((PB*PB)*(PB*PB))*one_over_2p + 72*(PA*PA*PA)*(PB*PB)*(one_over_2p*one_over_2p) + 60*(PA*PA*PA)*(one_over_2p*one_over_2p*one_over_2p) + 72*(PA*PA)*(PB*PB*PB)*(one_over_2p*one_over_2p) + 360*(PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 12*PA*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p) + 360*PA*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 420*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 60*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 420*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,4}_2
 */
inline Vec8d hermite_e_symbolic_4_4_2(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*((PA*PA)*(PA*PA))*(PB*PB)*(one_over_2p*one_over_2p) + 6*((PA*PA)*(PA*PA))*(one_over_2p*one_over_2p*one_over_2p) + 16*(PA*PA*PA)*(PB*PB*PB)*(one_over_2p*one_over_2p) + 96*(PA*PA*PA)*PB*(one_over_2p*one_over_2p*one_over_2p) + 6*(PA*PA)*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p) + 216*(PA*PA)*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 270*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 96*PA*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 720*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 6*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p*one_over_2p) + 270*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 420*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,4}_3
 */
inline Vec8d hermite_e_symbolic_4_4_3(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*((PA*PA)*(PA*PA))*PB*(one_over_2p*one_over_2p*one_over_2p) + 24*(PA*PA*PA)*(PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 40*(PA*PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 24*(PA*PA)*(PB*PB*PB)*(one_over_2p*one_over_2p*one_over_2p) + 240*(PA*PA)*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 4*PA*((PB*PB)*(PB*PB))*(one_over_2p*one_over_2p*one_over_2p) + 240*PA*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 420*PA*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 40*(PB*PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 420*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,4}_4
 */
inline Vec8d hermite_e_symbolic_4_4_4(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return ((PA*PA)*(PA*PA))*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 16*(PA*PA*PA)*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 36*(PA*PA)*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 90*(PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 16*PA*(PB*PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 240*PA*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + ((PB*PB)*(PB*PB))*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)) + 90*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 210*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,4}_5
 */
inline Vec8d hermite_e_symbolic_4_4_5(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*(PA*PA*PA)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 24*(PA*PA)*PB*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 24*PA*(PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 84*PA*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)) + 4*(PB*PB*PB)*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)*one_over_2p) + 84*PB*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p));
}

/**
 * @brief Symbolic E^{4,4}_6
 */
inline Vec8d hermite_e_symbolic_4_4_6(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 6*(PA*PA)*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)) + 16*PA*PB*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)) + 6*(PB*PB)*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)) + 28*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,4}_7
 */
inline Vec8d hermite_e_symbolic_4_4_7(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return 4*PA*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)*one_over_2p) + 4*PB*((one_over_2p*one_over_2p*one_over_2p)*(one_over_2p*one_over_2p*one_over_2p)*one_over_2p);
}

/**
 * @brief Symbolic E^{4,4}_8
 */
inline Vec8d hermite_e_symbolic_4_4_8(Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    return (((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p))*((one_over_2p*one_over_2p)*(one_over_2p*one_over_2p)));
}

/**
 * @brief Runtime dispatcher for symbolic Hermite E coefficients
 */
inline Vec8d dispatch_hermite_e_symbolic(int nA, int nB, int t,
                                          Vec8d PA, Vec8d PB, Vec8d one_over_2p) {
    if (nA == 0 && nB == 0 && t == 0) return hermite_e_symbolic_0_0_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 1 && t == 0) return hermite_e_symbolic_0_1_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 1 && t == 1) return hermite_e_symbolic_0_1_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 2 && t == 0) return hermite_e_symbolic_0_2_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 2 && t == 1) return hermite_e_symbolic_0_2_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 2 && t == 2) return hermite_e_symbolic_0_2_2(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 0) return hermite_e_symbolic_0_3_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 1) return hermite_e_symbolic_0_3_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 2) return hermite_e_symbolic_0_3_2(PA, PB, one_over_2p);
    if (nA == 0 && nB == 3 && t == 3) return hermite_e_symbolic_0_3_3(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 0) return hermite_e_symbolic_0_4_0(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 1) return hermite_e_symbolic_0_4_1(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 2) return hermite_e_symbolic_0_4_2(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 3) return hermite_e_symbolic_0_4_3(PA, PB, one_over_2p);
    if (nA == 0 && nB == 4 && t == 4) return hermite_e_symbolic_0_4_4(PA, PB, one_over_2p);
    if (nA == 1 && nB == 0 && t == 0) return hermite_e_symbolic_1_0_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 0 && t == 1) return hermite_e_symbolic_1_0_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 1 && t == 0) return hermite_e_symbolic_1_1_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 1 && t == 1) return hermite_e_symbolic_1_1_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 1 && t == 2) return hermite_e_symbolic_1_1_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 0) return hermite_e_symbolic_1_2_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 1) return hermite_e_symbolic_1_2_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 2) return hermite_e_symbolic_1_2_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 2 && t == 3) return hermite_e_symbolic_1_2_3(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 0) return hermite_e_symbolic_1_3_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 1) return hermite_e_symbolic_1_3_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 2) return hermite_e_symbolic_1_3_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 3) return hermite_e_symbolic_1_3_3(PA, PB, one_over_2p);
    if (nA == 1 && nB == 3 && t == 4) return hermite_e_symbolic_1_3_4(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 0) return hermite_e_symbolic_1_4_0(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 1) return hermite_e_symbolic_1_4_1(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 2) return hermite_e_symbolic_1_4_2(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 3) return hermite_e_symbolic_1_4_3(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 4) return hermite_e_symbolic_1_4_4(PA, PB, one_over_2p);
    if (nA == 1 && nB == 4 && t == 5) return hermite_e_symbolic_1_4_5(PA, PB, one_over_2p);
    if (nA == 2 && nB == 0 && t == 0) return hermite_e_symbolic_2_0_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 0 && t == 1) return hermite_e_symbolic_2_0_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 0 && t == 2) return hermite_e_symbolic_2_0_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 0) return hermite_e_symbolic_2_1_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 1) return hermite_e_symbolic_2_1_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 2) return hermite_e_symbolic_2_1_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 1 && t == 3) return hermite_e_symbolic_2_1_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 0) return hermite_e_symbolic_2_2_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 1) return hermite_e_symbolic_2_2_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 2) return hermite_e_symbolic_2_2_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 3) return hermite_e_symbolic_2_2_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 2 && t == 4) return hermite_e_symbolic_2_2_4(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 0) return hermite_e_symbolic_2_3_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 1) return hermite_e_symbolic_2_3_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 2) return hermite_e_symbolic_2_3_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 3) return hermite_e_symbolic_2_3_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 4) return hermite_e_symbolic_2_3_4(PA, PB, one_over_2p);
    if (nA == 2 && nB == 3 && t == 5) return hermite_e_symbolic_2_3_5(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 0) return hermite_e_symbolic_2_4_0(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 1) return hermite_e_symbolic_2_4_1(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 2) return hermite_e_symbolic_2_4_2(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 3) return hermite_e_symbolic_2_4_3(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 4) return hermite_e_symbolic_2_4_4(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 5) return hermite_e_symbolic_2_4_5(PA, PB, one_over_2p);
    if (nA == 2 && nB == 4 && t == 6) return hermite_e_symbolic_2_4_6(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 0) return hermite_e_symbolic_3_0_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 1) return hermite_e_symbolic_3_0_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 2) return hermite_e_symbolic_3_0_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 0 && t == 3) return hermite_e_symbolic_3_0_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 0) return hermite_e_symbolic_3_1_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 1) return hermite_e_symbolic_3_1_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 2) return hermite_e_symbolic_3_1_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 3) return hermite_e_symbolic_3_1_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 1 && t == 4) return hermite_e_symbolic_3_1_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 0) return hermite_e_symbolic_3_2_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 1) return hermite_e_symbolic_3_2_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 2) return hermite_e_symbolic_3_2_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 3) return hermite_e_symbolic_3_2_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 4) return hermite_e_symbolic_3_2_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 2 && t == 5) return hermite_e_symbolic_3_2_5(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 0) return hermite_e_symbolic_3_3_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 1) return hermite_e_symbolic_3_3_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 2) return hermite_e_symbolic_3_3_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 3) return hermite_e_symbolic_3_3_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 4) return hermite_e_symbolic_3_3_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 5) return hermite_e_symbolic_3_3_5(PA, PB, one_over_2p);
    if (nA == 3 && nB == 3 && t == 6) return hermite_e_symbolic_3_3_6(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 0) return hermite_e_symbolic_3_4_0(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 1) return hermite_e_symbolic_3_4_1(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 2) return hermite_e_symbolic_3_4_2(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 3) return hermite_e_symbolic_3_4_3(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 4) return hermite_e_symbolic_3_4_4(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 5) return hermite_e_symbolic_3_4_5(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 6) return hermite_e_symbolic_3_4_6(PA, PB, one_over_2p);
    if (nA == 3 && nB == 4 && t == 7) return hermite_e_symbolic_3_4_7(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 0) return hermite_e_symbolic_4_0_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 1) return hermite_e_symbolic_4_0_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 2) return hermite_e_symbolic_4_0_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 3) return hermite_e_symbolic_4_0_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 0 && t == 4) return hermite_e_symbolic_4_0_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 0) return hermite_e_symbolic_4_1_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 1) return hermite_e_symbolic_4_1_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 2) return hermite_e_symbolic_4_1_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 3) return hermite_e_symbolic_4_1_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 4) return hermite_e_symbolic_4_1_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 1 && t == 5) return hermite_e_symbolic_4_1_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 0) return hermite_e_symbolic_4_2_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 1) return hermite_e_symbolic_4_2_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 2) return hermite_e_symbolic_4_2_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 3) return hermite_e_symbolic_4_2_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 4) return hermite_e_symbolic_4_2_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 5) return hermite_e_symbolic_4_2_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 2 && t == 6) return hermite_e_symbolic_4_2_6(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 0) return hermite_e_symbolic_4_3_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 1) return hermite_e_symbolic_4_3_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 2) return hermite_e_symbolic_4_3_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 3) return hermite_e_symbolic_4_3_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 4) return hermite_e_symbolic_4_3_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 5) return hermite_e_symbolic_4_3_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 6) return hermite_e_symbolic_4_3_6(PA, PB, one_over_2p);
    if (nA == 4 && nB == 3 && t == 7) return hermite_e_symbolic_4_3_7(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 0) return hermite_e_symbolic_4_4_0(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 1) return hermite_e_symbolic_4_4_1(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 2) return hermite_e_symbolic_4_4_2(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 3) return hermite_e_symbolic_4_4_3(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 4) return hermite_e_symbolic_4_4_4(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 5) return hermite_e_symbolic_4_4_5(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 6) return hermite_e_symbolic_4_4_6(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 7) return hermite_e_symbolic_4_4_7(PA, PB, one_over_2p);
    if (nA == 4 && nB == 4 && t == 8) return hermite_e_symbolic_4_4_8(PA, PB, one_over_2p);
    return Vec8d(0.0);  // Invalid indices
}

} // namespace symbolic
} // namespace recursum
