#include "recursum_eri_scalars.h"


// ERI class (ssss) [array]: 1 Cartesian integrals, 1 DAG nodes, peak liveness 0 slots
__attribute__((noinline)) void recursum_eri_ssss(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    out[0] = kf[0];
}

// ERI class (sssp) [array]: 3 Cartesian integrals, 8 DAG nodes, peak liveness 1 slots
__attribute__((noinline)) void recursum_eri_sssp(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[1];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    out[2] = sc[0] - s.CDz * kf[0];
    sc[0] = s.QCy * kf[0] + s.WQy * kf[1];
    out[1] = sc[0] - s.CDy * kf[0];
    sc[0] = s.QCx * kf[0] + s.WQx * kf[1];
    out[0] = sc[0] - s.CDx * kf[0];
}

// ERI class (ssps) [array]: 3 Cartesian integrals, 5 DAG nodes, peak liveness 0 slots
__attribute__((noinline)) void recursum_eri_ssps(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    out[2] = s.QCz * kf[0] + s.WQz * kf[1];
    out[1] = s.QCy * kf[0] + s.WQy * kf[1];
    out[0] = s.QCx * kf[0] + s.WQx * kf[1];
}

// ERI class (sspp) [array]: 9 Cartesian integrals, 24 DAG nodes, peak liveness 7 slots
__attribute__((noinline)) void recursum_eri_sspp(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[7];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    sc[1] = s.QCy * kf[0] + s.WQy * kf[1];
    sc[2] = s.QCx * kf[0] + s.WQx * kf[1];
    sc[3] = s.QCz * kf[1] + s.WQz * kf[2];
    sc[4] = s.QCy * kf[1] + s.WQy * kf[2];
    sc[5] = s.QCx * kf[1] + s.WQx * kf[2];
    sc[6] = s.QCz * sc[0] + s.WQz * sc[3] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    out[8] = sc[6] - s.CDz * sc[0];
    sc[6] = s.QCy * sc[0] + s.WQy * sc[3];
    out[7] = sc[6] - s.CDy * sc[0];
    out[5] = sc[6] - s.CDz * sc[1];
    sc[6] = s.QCy * sc[1] + s.WQy * sc[4] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    out[4] = sc[6] - s.CDy * sc[1];
    sc[6] = s.QCx * sc[0] + s.WQx * sc[3];
    out[6] = sc[6] - s.CDx * sc[0];
    out[2] = sc[6] - s.CDz * sc[2];
    sc[6] = s.QCx * sc[1] + s.WQx * sc[4];
    out[3] = sc[6] - s.CDx * sc[1];
    out[1] = sc[6] - s.CDy * sc[2];
    sc[6] = s.QCx * sc[2] + s.WQx * sc[5] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    out[0] = sc[6] - s.CDx * sc[2];
}

// ERI class (spss) [array]: 3 Cartesian integrals, 8 DAG nodes, peak liveness 1 slots
__attribute__((noinline)) void recursum_eri_spss(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[1];
    sc[0] = s.PAz * kf[0] + s.WPz * kf[1];
    out[2] = sc[0] - s.ABz * kf[0];
    sc[0] = s.PAy * kf[0] + s.WPy * kf[1];
    out[1] = sc[0] - s.ABy * kf[0];
    sc[0] = s.PAx * kf[0] + s.WPx * kf[1];
    out[0] = sc[0] - s.ABx * kf[0];
}

// ERI class (spsp) [array]: 9 Cartesian integrals, 42 DAG nodes, peak liveness 14 slots
__attribute__((noinline)) void recursum_eri_spsp(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[14];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    sc[1] = sc[0] - s.CDz * kf[0];
    sc[2] = s.QCy * kf[0] + s.WQy * kf[1];
    sc[3] = sc[2] - s.CDy * kf[0];
    sc[4] = s.QCx * kf[0] + s.WQx * kf[1];
    sc[5] = sc[4] - s.CDx * kf[0];
    sc[6] = s.PAz * kf[0] + s.WPz * kf[1];
    sc[7] = s.PAy * kf[0] + s.WPy * kf[1];
    sc[8] = s.PAx * kf[0] + s.WPx * kf[1];
    sc[9] = s.QCz * kf[1] + s.WQz * kf[2];
    sc[10] = s.QCy * kf[1] + s.WQy * kf[2];
    sc[11] = s.QCx * kf[1] + s.WQx * kf[2];
    sc[12] = s.PAz * sc[0] + s.WPz * sc[9] + 1.0 * s.inv_2zpq * kf[1];
    sc[13] = sc[12] - s.CDz * sc[6];
    out[8] = sc[13] - s.ABz * sc[1];
    sc[13] = s.PAz * sc[2] + s.WPz * sc[10];
    sc[12] = sc[13] - s.CDy * sc[6];
    out[7] = sc[12] - s.ABz * sc[3];
    sc[12] = s.PAz * sc[4] + s.WPz * sc[11];
    sc[13] = sc[12] - s.CDx * sc[6];
    out[6] = sc[13] - s.ABz * sc[5];
    sc[13] = s.PAy * sc[0] + s.WPy * sc[9];
    sc[12] = sc[13] - s.CDz * sc[7];
    out[5] = sc[12] - s.ABy * sc[1];
    sc[12] = s.PAy * sc[2] + s.WPy * sc[10] + 1.0 * s.inv_2zpq * kf[1];
    sc[13] = sc[12] - s.CDy * sc[7];
    out[4] = sc[13] - s.ABy * sc[3];
    sc[13] = s.PAy * sc[4] + s.WPy * sc[11];
    sc[12] = sc[13] - s.CDx * sc[7];
    out[3] = sc[12] - s.ABy * sc[5];
    sc[12] = s.PAx * sc[0] + s.WPx * sc[9];
    sc[9] = sc[12] - s.CDz * sc[8];
    out[2] = sc[9] - s.ABx * sc[1];
    sc[9] = s.PAx * sc[2] + s.WPx * sc[10];
    sc[10] = sc[9] - s.CDy * sc[8];
    out[1] = sc[10] - s.ABx * sc[3];
    sc[10] = s.PAx * sc[4] + s.WPx * sc[11] + 1.0 * s.inv_2zpq * kf[1];
    sc[11] = sc[10] - s.CDx * sc[8];
    out[0] = sc[11] - s.ABx * sc[5];
}

// ERI class (spps) [array]: 9 Cartesian integrals, 27 DAG nodes, peak liveness 7 slots
__attribute__((noinline)) void recursum_eri_spps(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[7];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    sc[1] = s.QCy * kf[0] + s.WQy * kf[1];
    sc[2] = s.QCx * kf[0] + s.WQx * kf[1];
    sc[3] = s.QCz * kf[1] + s.WQz * kf[2];
    sc[4] = s.QCy * kf[1] + s.WQy * kf[2];
    sc[5] = s.QCx * kf[1] + s.WQx * kf[2];
    sc[6] = s.PAz * sc[0] + s.WPz * sc[3] + 1.0 * s.inv_2zpq * kf[1];
    out[8] = sc[6] - s.ABz * sc[0];
    sc[6] = s.PAz * sc[1] + s.WPz * sc[4];
    out[7] = sc[6] - s.ABz * sc[1];
    sc[6] = s.PAz * sc[2] + s.WPz * sc[5];
    out[6] = sc[6] - s.ABz * sc[2];
    sc[6] = s.PAy * sc[0] + s.WPy * sc[3];
    out[5] = sc[6] - s.ABy * sc[0];
    sc[6] = s.PAy * sc[1] + s.WPy * sc[4] + 1.0 * s.inv_2zpq * kf[1];
    out[4] = sc[6] - s.ABy * sc[1];
    sc[6] = s.PAy * sc[2] + s.WPy * sc[5];
    out[3] = sc[6] - s.ABy * sc[2];
    sc[6] = s.PAx * sc[0] + s.WPx * sc[3];
    out[2] = sc[6] - s.ABx * sc[0];
    sc[6] = s.PAx * sc[1] + s.WPx * sc[4];
    out[1] = sc[6] - s.ABx * sc[1];
    sc[6] = s.PAx * sc[2] + s.WPx * sc[5] + 1.0 * s.inv_2zpq * kf[1];
    out[0] = sc[6] - s.ABx * sc[2];
}

// ERI class (sppp) [array]: 27 Cartesian integrals, 115 DAG nodes, peak liveness 35 slots
__attribute__((noinline)) void recursum_eri_sppp(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[35];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    sc[1] = s.QCy * kf[0] + s.WQy * kf[1];
    sc[2] = s.QCx * kf[0] + s.WQx * kf[1];
    sc[3] = s.QCz * kf[1] + s.WQz * kf[2];
    sc[4] = s.QCy * kf[1] + s.WQy * kf[2];
    sc[5] = s.QCx * kf[1] + s.WQx * kf[2];
    sc[6] = s.QCz * kf[2] + s.WQz * kf[3];
    sc[7] = s.QCy * kf[2] + s.WQy * kf[3];
    sc[8] = s.QCx * kf[2] + s.WQx * kf[3];
    sc[9] = s.QCz * sc[0] + s.WQz * sc[3] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    sc[10] = sc[9] - s.CDz * sc[0];
    sc[11] = s.QCy * sc[0] + s.WQy * sc[3];
    sc[12] = sc[11] - s.CDy * sc[0];
    sc[13] = sc[11] - s.CDz * sc[1];
    sc[14] = s.QCy * sc[1] + s.WQy * sc[4] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    sc[15] = sc[14] - s.CDy * sc[1];
    sc[16] = s.QCx * sc[0] + s.WQx * sc[3];
    sc[17] = sc[16] - s.CDx * sc[0];
    sc[18] = sc[16] - s.CDz * sc[2];
    sc[19] = s.QCx * sc[1] + s.WQx * sc[4];
    sc[20] = sc[19] - s.CDx * sc[1];
    sc[21] = sc[19] - s.CDy * sc[2];
    sc[22] = s.QCx * sc[2] + s.WQx * sc[5] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    sc[23] = sc[22] - s.CDx * sc[2];
    sc[24] = s.PAz * sc[0] + s.WPz * sc[3] + 1.0 * s.inv_2zpq * kf[1];
    sc[25] = s.PAz * sc[1] + s.WPz * sc[4];
    sc[26] = s.PAz * sc[2] + s.WPz * sc[5];
    sc[27] = s.PAy * sc[0] + s.WPy * sc[3];
    sc[28] = s.PAy * sc[1] + s.WPy * sc[4] + 1.0 * s.inv_2zpq * kf[1];
    sc[29] = s.PAy * sc[2] + s.WPy * sc[5];
    sc[30] = s.PAx * sc[0] + s.WPx * sc[3];
    sc[0] = s.PAx * sc[1] + s.WPx * sc[4];
    sc[1] = s.PAx * sc[2] + s.WPx * sc[5] + 1.0 * s.inv_2zpq * kf[1];
    sc[2] = s.QCz * sc[3] + s.WQz * sc[6] + 1.0 * s.inv_2zq * kf[1] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[2];
    sc[31] = s.QCy * sc[3] + s.WQy * sc[6];
    sc[32] = s.QCy * sc[4] + s.WQy * sc[7] + 1.0 * s.inv_2zq * kf[1] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[2];
    sc[33] = s.QCx * sc[3] + s.WQx * sc[6];
    sc[6] = s.QCx * sc[4] + s.WQx * sc[7];
    sc[7] = s.QCx * sc[5] + s.WQx * sc[8] + 1.0 * s.inv_2zq * kf[1] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[2];
    sc[8] = s.PAz * sc[9] + s.WPz * sc[2] + 2.0 * s.inv_2zpq * sc[3];
    sc[34] = sc[8] - s.CDz * sc[24];
    out[26] = sc[34] - s.ABz * sc[10];
    sc[34] = s.PAz * sc[11] + s.WPz * sc[31] + 1.0 * s.inv_2zpq * sc[4];
    sc[8] = sc[34] - s.CDy * sc[24];
    out[25] = sc[8] - s.ABz * sc[12];
    sc[8] = sc[34] - s.CDz * sc[25];
    out[23] = sc[8] - s.ABz * sc[13];
    sc[8] = s.PAz * sc[14] + s.WPz * sc[32];
    sc[34] = sc[8] - s.CDy * sc[25];
    out[22] = sc[34] - s.ABz * sc[15];
    sc[34] = s.PAz * sc[16] + s.WPz * sc[33] + 1.0 * s.inv_2zpq * sc[5];
    sc[8] = sc[34] - s.CDx * sc[24];
    out[24] = sc[8] - s.ABz * sc[17];
    sc[8] = sc[34] - s.CDz * sc[26];
    out[20] = sc[8] - s.ABz * sc[18];
    sc[8] = s.PAz * sc[19] + s.WPz * sc[6];
    sc[34] = sc[8] - s.CDx * sc[25];
    out[21] = sc[34] - s.ABz * sc[20];
    sc[34] = sc[8] - s.CDy * sc[26];
    out[19] = sc[34] - s.ABz * sc[21];
    sc[34] = s.PAz * sc[22] + s.WPz * sc[7];
    sc[8] = sc[34] - s.CDx * sc[26];
    out[18] = sc[8] - s.ABz * sc[23];
    sc[8] = s.PAy * sc[9] + s.WPy * sc[2];
    sc[34] = sc[8] - s.CDz * sc[27];
    out[17] = sc[34] - s.ABy * sc[10];
    sc[34] = s.PAy * sc[11] + s.WPy * sc[31] + 1.0 * s.inv_2zpq * sc[3];
    sc[8] = sc[34] - s.CDy * sc[27];
    out[16] = sc[8] - s.ABy * sc[12];
    sc[8] = sc[34] - s.CDz * sc[28];
    out[14] = sc[8] - s.ABy * sc[13];
    sc[8] = s.PAy * sc[14] + s.WPy * sc[32] + 2.0 * s.inv_2zpq * sc[4];
    sc[34] = sc[8] - s.CDy * sc[28];
    out[13] = sc[34] - s.ABy * sc[15];
    sc[34] = s.PAy * sc[16] + s.WPy * sc[33];
    sc[8] = sc[34] - s.CDx * sc[27];
    out[15] = sc[8] - s.ABy * sc[17];
    sc[8] = sc[34] - s.CDz * sc[29];
    out[11] = sc[8] - s.ABy * sc[18];
    sc[8] = s.PAy * sc[19] + s.WPy * sc[6] + 1.0 * s.inv_2zpq * sc[5];
    sc[34] = sc[8] - s.CDx * sc[28];
    out[12] = sc[34] - s.ABy * sc[20];
    sc[34] = sc[8] - s.CDy * sc[29];
    out[10] = sc[34] - s.ABy * sc[21];
    sc[34] = s.PAy * sc[22] + s.WPy * sc[7];
    sc[8] = sc[34] - s.CDx * sc[29];
    out[9] = sc[8] - s.ABy * sc[23];
    sc[8] = s.PAx * sc[9] + s.WPx * sc[2];
    sc[2] = sc[8] - s.CDz * sc[30];
    out[8] = sc[2] - s.ABx * sc[10];
    sc[2] = s.PAx * sc[11] + s.WPx * sc[31];
    sc[31] = sc[2] - s.CDy * sc[30];
    out[7] = sc[31] - s.ABx * sc[12];
    sc[31] = sc[2] - s.CDz * sc[0];
    out[5] = sc[31] - s.ABx * sc[13];
    sc[31] = s.PAx * sc[14] + s.WPx * sc[32];
    sc[32] = sc[31] - s.CDy * sc[0];
    out[4] = sc[32] - s.ABx * sc[15];
    sc[32] = s.PAx * sc[16] + s.WPx * sc[33] + 1.0 * s.inv_2zpq * sc[3];
    sc[33] = sc[32] - s.CDx * sc[30];
    out[6] = sc[33] - s.ABx * sc[17];
    sc[33] = sc[32] - s.CDz * sc[1];
    out[2] = sc[33] - s.ABx * sc[18];
    sc[33] = s.PAx * sc[19] + s.WPx * sc[6] + 1.0 * s.inv_2zpq * sc[4];
    sc[6] = sc[33] - s.CDx * sc[0];
    out[3] = sc[6] - s.ABx * sc[20];
    sc[6] = sc[33] - s.CDy * sc[1];
    out[1] = sc[6] - s.ABx * sc[21];
    sc[6] = s.PAx * sc[22] + s.WPx * sc[7] + 2.0 * s.inv_2zpq * sc[5];
    sc[7] = sc[6] - s.CDx * sc[1];
    out[0] = sc[7] - s.ABx * sc[23];
}

// ERI class (psss) [array]: 3 Cartesian integrals, 5 DAG nodes, peak liveness 0 slots
__attribute__((noinline)) void recursum_eri_psss(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    out[2] = s.PAz * kf[0] + s.WPz * kf[1];
    out[1] = s.PAy * kf[0] + s.WPy * kf[1];
    out[0] = s.PAx * kf[0] + s.WPx * kf[1];
}

// ERI class (pssp) [array]: 9 Cartesian integrals, 30 DAG nodes, peak liveness 10 slots
__attribute__((noinline)) void recursum_eri_pssp(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[10];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    sc[1] = s.QCy * kf[0] + s.WQy * kf[1];
    sc[2] = s.QCx * kf[0] + s.WQx * kf[1];
    sc[3] = s.PAz * kf[0] + s.WPz * kf[1];
    sc[4] = s.PAy * kf[0] + s.WPy * kf[1];
    sc[5] = s.PAx * kf[0] + s.WPx * kf[1];
    sc[6] = s.QCz * kf[1] + s.WQz * kf[2];
    sc[7] = s.QCy * kf[1] + s.WQy * kf[2];
    sc[8] = s.QCx * kf[1] + s.WQx * kf[2];
    sc[9] = s.PAz * sc[0] + s.WPz * sc[6] + 1.0 * s.inv_2zpq * kf[1];
    out[8] = sc[9] - s.CDz * sc[3];
    sc[9] = s.PAz * sc[1] + s.WPz * sc[7];
    out[7] = sc[9] - s.CDy * sc[3];
    sc[9] = s.PAz * sc[2] + s.WPz * sc[8];
    out[6] = sc[9] - s.CDx * sc[3];
    sc[9] = s.PAy * sc[0] + s.WPy * sc[6];
    out[5] = sc[9] - s.CDz * sc[4];
    sc[9] = s.PAy * sc[1] + s.WPy * sc[7] + 1.0 * s.inv_2zpq * kf[1];
    out[4] = sc[9] - s.CDy * sc[4];
    sc[9] = s.PAy * sc[2] + s.WPy * sc[8];
    out[3] = sc[9] - s.CDx * sc[4];
    sc[9] = s.PAx * sc[0] + s.WPx * sc[6];
    out[2] = sc[9] - s.CDz * sc[5];
    sc[9] = s.PAx * sc[1] + s.WPx * sc[7];
    out[1] = sc[9] - s.CDy * sc[5];
    sc[9] = s.PAx * sc[2] + s.WPx * sc[8] + 1.0 * s.inv_2zpq * kf[1];
    out[0] = sc[9] - s.CDx * sc[5];
}

// ERI class (psps) [array]: 9 Cartesian integrals, 18 DAG nodes, peak liveness 6 slots
__attribute__((noinline)) void recursum_eri_psps(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[6];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    sc[1] = s.QCy * kf[0] + s.WQy * kf[1];
    sc[2] = s.QCx * kf[0] + s.WQx * kf[1];
    sc[3] = s.QCz * kf[1] + s.WQz * kf[2];
    sc[4] = s.QCy * kf[1] + s.WQy * kf[2];
    sc[5] = s.QCx * kf[1] + s.WQx * kf[2];
    out[8] = s.PAz * sc[0] + s.WPz * sc[3] + 1.0 * s.inv_2zpq * kf[1];
    out[7] = s.PAz * sc[1] + s.WPz * sc[4];
    out[6] = s.PAz * sc[2] + s.WPz * sc[5];
    out[5] = s.PAy * sc[0] + s.WPy * sc[3];
    out[4] = s.PAy * sc[1] + s.WPy * sc[4] + 1.0 * s.inv_2zpq * kf[1];
    out[3] = s.PAy * sc[2] + s.WPy * sc[5];
    out[2] = s.PAx * sc[0] + s.WPx * sc[3];
    out[1] = s.PAx * sc[1] + s.WPx * sc[4];
    out[0] = s.PAx * sc[2] + s.WPx * sc[5] + 1.0 * s.inv_2zpq * kf[1];
}

// ERI class (pspp) [array]: 27 Cartesian integrals, 79 DAG nodes, peak liveness 25 slots
__attribute__((noinline)) void recursum_eri_pspp(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[25];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    sc[1] = s.QCy * kf[0] + s.WQy * kf[1];
    sc[2] = s.QCx * kf[0] + s.WQx * kf[1];
    sc[3] = s.QCz * kf[1] + s.WQz * kf[2];
    sc[4] = s.QCy * kf[1] + s.WQy * kf[2];
    sc[5] = s.QCx * kf[1] + s.WQx * kf[2];
    sc[6] = s.QCz * kf[2] + s.WQz * kf[3];
    sc[7] = s.QCy * kf[2] + s.WQy * kf[3];
    sc[8] = s.QCx * kf[2] + s.WQx * kf[3];
    sc[9] = s.QCz * sc[0] + s.WQz * sc[3] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    sc[10] = s.QCy * sc[0] + s.WQy * sc[3];
    sc[11] = s.QCy * sc[1] + s.WQy * sc[4] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    sc[12] = s.QCx * sc[0] + s.WQx * sc[3];
    sc[13] = s.QCx * sc[1] + s.WQx * sc[4];
    sc[14] = s.QCx * sc[2] + s.WQx * sc[5] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    sc[15] = s.PAz * sc[0] + s.WPz * sc[3] + 1.0 * s.inv_2zpq * kf[1];
    sc[16] = s.PAz * sc[1] + s.WPz * sc[4];
    sc[17] = s.PAz * sc[2] + s.WPz * sc[5];
    sc[18] = s.PAy * sc[0] + s.WPy * sc[3];
    sc[19] = s.PAy * sc[1] + s.WPy * sc[4] + 1.0 * s.inv_2zpq * kf[1];
    sc[20] = s.PAy * sc[2] + s.WPy * sc[5];
    sc[21] = s.PAx * sc[0] + s.WPx * sc[3];
    sc[0] = s.PAx * sc[1] + s.WPx * sc[4];
    sc[1] = s.PAx * sc[2] + s.WPx * sc[5] + 1.0 * s.inv_2zpq * kf[1];
    sc[2] = s.QCz * sc[3] + s.WQz * sc[6] + 1.0 * s.inv_2zq * kf[1] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[2];
    sc[22] = s.QCy * sc[3] + s.WQy * sc[6];
    sc[23] = s.QCy * sc[4] + s.WQy * sc[7] + 1.0 * s.inv_2zq * kf[1] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[2];
    sc[24] = s.QCx * sc[3] + s.WQx * sc[6];
    sc[6] = s.QCx * sc[4] + s.WQx * sc[7];
    sc[7] = s.QCx * sc[5] + s.WQx * sc[8] + 1.0 * s.inv_2zq * kf[1] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[2];
    sc[8] = s.PAz * sc[9] + s.WPz * sc[2] + 2.0 * s.inv_2zpq * sc[3];
    out[26] = sc[8] - s.CDz * sc[15];
    sc[8] = s.PAz * sc[10] + s.WPz * sc[22] + 1.0 * s.inv_2zpq * sc[4];
    out[25] = sc[8] - s.CDy * sc[15];
    out[23] = sc[8] - s.CDz * sc[16];
    sc[8] = s.PAz * sc[11] + s.WPz * sc[23];
    out[22] = sc[8] - s.CDy * sc[16];
    sc[8] = s.PAz * sc[12] + s.WPz * sc[24] + 1.0 * s.inv_2zpq * sc[5];
    out[24] = sc[8] - s.CDx * sc[15];
    out[20] = sc[8] - s.CDz * sc[17];
    sc[8] = s.PAz * sc[13] + s.WPz * sc[6];
    out[21] = sc[8] - s.CDx * sc[16];
    out[19] = sc[8] - s.CDy * sc[17];
    sc[8] = s.PAz * sc[14] + s.WPz * sc[7];
    out[18] = sc[8] - s.CDx * sc[17];
    sc[8] = s.PAy * sc[9] + s.WPy * sc[2];
    out[17] = sc[8] - s.CDz * sc[18];
    sc[8] = s.PAy * sc[10] + s.WPy * sc[22] + 1.0 * s.inv_2zpq * sc[3];
    out[16] = sc[8] - s.CDy * sc[18];
    out[14] = sc[8] - s.CDz * sc[19];
    sc[8] = s.PAy * sc[11] + s.WPy * sc[23] + 2.0 * s.inv_2zpq * sc[4];
    out[13] = sc[8] - s.CDy * sc[19];
    sc[8] = s.PAy * sc[12] + s.WPy * sc[24];
    out[15] = sc[8] - s.CDx * sc[18];
    out[11] = sc[8] - s.CDz * sc[20];
    sc[8] = s.PAy * sc[13] + s.WPy * sc[6] + 1.0 * s.inv_2zpq * sc[5];
    out[12] = sc[8] - s.CDx * sc[19];
    out[10] = sc[8] - s.CDy * sc[20];
    sc[8] = s.PAy * sc[14] + s.WPy * sc[7];
    out[9] = sc[8] - s.CDx * sc[20];
    sc[8] = s.PAx * sc[9] + s.WPx * sc[2];
    out[8] = sc[8] - s.CDz * sc[21];
    sc[8] = s.PAx * sc[10] + s.WPx * sc[22];
    out[7] = sc[8] - s.CDy * sc[21];
    out[5] = sc[8] - s.CDz * sc[0];
    sc[8] = s.PAx * sc[11] + s.WPx * sc[23];
    out[4] = sc[8] - s.CDy * sc[0];
    sc[8] = s.PAx * sc[12] + s.WPx * sc[24] + 1.0 * s.inv_2zpq * sc[3];
    out[6] = sc[8] - s.CDx * sc[21];
    out[2] = sc[8] - s.CDz * sc[1];
    sc[8] = s.PAx * sc[13] + s.WPx * sc[6] + 1.0 * s.inv_2zpq * sc[4];
    out[3] = sc[8] - s.CDx * sc[0];
    out[1] = sc[8] - s.CDy * sc[1];
    sc[8] = s.PAx * sc[14] + s.WPx * sc[7] + 2.0 * s.inv_2zpq * sc[5];
    out[0] = sc[8] - s.CDx * sc[1];
}

// ERI class (ppss) [array]: 9 Cartesian integrals, 24 DAG nodes, peak liveness 7 slots
__attribute__((noinline)) void recursum_eri_ppss(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[7];
    sc[0] = s.PAz * kf[0] + s.WPz * kf[1];
    sc[1] = s.PAy * kf[0] + s.WPy * kf[1];
    sc[2] = s.PAx * kf[0] + s.WPx * kf[1];
    sc[3] = s.PAz * kf[1] + s.WPz * kf[2];
    sc[4] = s.PAy * kf[1] + s.WPy * kf[2];
    sc[5] = s.PAx * kf[1] + s.WPx * kf[2];
    sc[6] = s.PAz * sc[0] + s.WPz * sc[3] + 1.0 * s.inv_2zp * kf[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * kf[1];
    out[8] = sc[6] - s.ABz * sc[0];
    sc[6] = s.PAy * sc[0] + s.WPy * sc[3];
    out[7] = sc[6] - s.ABy * sc[0];
    out[5] = sc[6] - s.ABz * sc[1];
    sc[6] = s.PAy * sc[1] + s.WPy * sc[4] + 1.0 * s.inv_2zp * kf[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * kf[1];
    out[4] = sc[6] - s.ABy * sc[1];
    sc[6] = s.PAx * sc[0] + s.WPx * sc[3];
    out[6] = sc[6] - s.ABx * sc[0];
    out[2] = sc[6] - s.ABz * sc[2];
    sc[6] = s.PAx * sc[1] + s.WPx * sc[4];
    out[3] = sc[6] - s.ABx * sc[1];
    out[1] = sc[6] - s.ABy * sc[2];
    sc[6] = s.PAx * sc[2] + s.WPx * sc[5] + 1.0 * s.inv_2zp * kf[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * kf[1];
    out[0] = sc[6] - s.ABx * sc[2];
}

// ERI class (ppsp) [array]: 27 Cartesian integrals, 115 DAG nodes, peak liveness 44 slots
__attribute__((noinline)) void recursum_eri_ppsp(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[44];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    sc[1] = s.QCy * kf[0] + s.WQy * kf[1];
    sc[2] = s.QCx * kf[0] + s.WQx * kf[1];
    sc[3] = s.PAz * kf[0] + s.WPz * kf[1];
    sc[4] = s.PAy * kf[0] + s.WPy * kf[1];
    sc[5] = s.PAx * kf[0] + s.WPx * kf[1];
    sc[6] = s.QCz * kf[1] + s.WQz * kf[2];
    sc[7] = s.QCy * kf[1] + s.WQy * kf[2];
    sc[8] = s.QCx * kf[1] + s.WQx * kf[2];
    sc[9] = s.PAz * kf[1] + s.WPz * kf[2];
    sc[10] = s.PAy * kf[1] + s.WPy * kf[2];
    sc[11] = s.PAx * kf[1] + s.WPx * kf[2];
    sc[12] = s.QCz * kf[2] + s.WQz * kf[3];
    sc[13] = s.QCy * kf[2] + s.WQy * kf[3];
    sc[14] = s.QCx * kf[2] + s.WQx * kf[3];
    sc[15] = s.PAz * sc[0] + s.WPz * sc[6] + 1.0 * s.inv_2zpq * kf[1];
    sc[16] = sc[15] - s.CDz * sc[3];
    sc[17] = s.PAz * sc[1] + s.WPz * sc[7];
    sc[18] = sc[17] - s.CDy * sc[3];
    sc[19] = s.PAz * sc[2] + s.WPz * sc[8];
    sc[20] = sc[19] - s.CDx * sc[3];
    sc[21] = s.PAz * sc[3] + s.WPz * sc[9] + 1.0 * s.inv_2zp * kf[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * kf[1];
    sc[22] = s.PAy * sc[0] + s.WPy * sc[6];
    sc[23] = sc[22] - s.CDz * sc[4];
    sc[24] = s.PAy * sc[1] + s.WPy * sc[7] + 1.0 * s.inv_2zpq * kf[1];
    sc[25] = sc[24] - s.CDy * sc[4];
    sc[26] = s.PAy * sc[2] + s.WPy * sc[8];
    sc[27] = sc[26] - s.CDx * sc[4];
    sc[28] = s.PAy * sc[3] + s.WPy * sc[9];
    sc[29] = s.PAy * sc[4] + s.WPy * sc[10] + 1.0 * s.inv_2zp * kf[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * kf[1];
    sc[30] = s.PAx * sc[0] + s.WPx * sc[6];
    sc[31] = sc[30] - s.CDz * sc[5];
    sc[32] = s.PAx * sc[1] + s.WPx * sc[7];
    sc[33] = sc[32] - s.CDy * sc[5];
    sc[34] = s.PAx * sc[2] + s.WPx * sc[8] + 1.0 * s.inv_2zpq * kf[1];
    sc[35] = sc[34] - s.CDx * sc[5];
    sc[36] = s.PAx * sc[3] + s.WPx * sc[9];
    sc[3] = s.PAx * sc[4] + s.WPx * sc[10];
    sc[4] = s.PAx * sc[5] + s.WPx * sc[11] + 1.0 * s.inv_2zp * kf[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * kf[1];
    sc[5] = s.PAz * sc[6] + s.WPz * sc[12] + 1.0 * s.inv_2zpq * kf[2];
    sc[37] = s.PAz * sc[7] + s.WPz * sc[13];
    sc[38] = s.PAz * sc[8] + s.WPz * sc[14];
    sc[39] = s.PAy * sc[6] + s.WPy * sc[12];
    sc[40] = s.PAy * sc[7] + s.WPy * sc[13] + 1.0 * s.inv_2zpq * kf[2];
    sc[41] = s.PAy * sc[8] + s.WPy * sc[14];
    sc[42] = s.PAx * sc[6] + s.WPx * sc[12];
    sc[12] = s.PAx * sc[7] + s.WPx * sc[13];
    sc[13] = s.PAx * sc[8] + s.WPx * sc[14] + 1.0 * s.inv_2zpq * kf[2];
    sc[14] = s.PAz * sc[15] + s.WPz * sc[5] + 1.0 * s.inv_2zp * sc[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[6] + 1.0 * s.inv_2zpq * sc[9];
    sc[43] = sc[14] - s.CDz * sc[21];
    out[26] = sc[43] - s.ABz * sc[16];
    sc[43] = s.PAz * sc[17] + s.WPz * sc[37] + 1.0 * s.inv_2zp * sc[1] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[7];
    sc[14] = sc[43] - s.CDy * sc[21];
    out[25] = sc[14] - s.ABz * sc[18];
    sc[14] = s.PAz * sc[19] + s.WPz * sc[38] + 1.0 * s.inv_2zp * sc[2] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[8];
    sc[43] = sc[14] - s.CDx * sc[21];
    out[24] = sc[43] - s.ABz * sc[20];
    sc[43] = s.PAy * sc[15] + s.WPy * sc[5];
    sc[14] = sc[43] - s.CDz * sc[28];
    out[23] = sc[14] - s.ABy * sc[16];
    out[17] = sc[14] - s.ABz * sc[23];
    sc[14] = s.PAy * sc[17] + s.WPy * sc[37] + 1.0 * s.inv_2zpq * sc[9];
    sc[43] = sc[14] - s.CDy * sc[28];
    out[22] = sc[43] - s.ABy * sc[18];
    out[16] = sc[43] - s.ABz * sc[25];
    sc[43] = s.PAy * sc[19] + s.WPy * sc[38];
    sc[14] = sc[43] - s.CDx * sc[28];
    out[21] = sc[14] - s.ABy * sc[20];
    out[15] = sc[14] - s.ABz * sc[27];
    sc[14] = s.PAy * sc[22] + s.WPy * sc[39] + 1.0 * s.inv_2zp * sc[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[6];
    sc[43] = sc[14] - s.CDz * sc[29];
    out[14] = sc[43] - s.ABy * sc[23];
    sc[43] = s.PAy * sc[24] + s.WPy * sc[40] + 1.0 * s.inv_2zp * sc[1] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[7] + 1.0 * s.inv_2zpq * sc[10];
    sc[14] = sc[43] - s.CDy * sc[29];
    out[13] = sc[14] - s.ABy * sc[25];
    sc[14] = s.PAy * sc[26] + s.WPy * sc[41] + 1.0 * s.inv_2zp * sc[2] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[8];
    sc[43] = sc[14] - s.CDx * sc[29];
    out[12] = sc[43] - s.ABy * sc[27];
    sc[43] = s.PAx * sc[15] + s.WPx * sc[5];
    sc[5] = sc[43] - s.CDz * sc[36];
    out[20] = sc[5] - s.ABx * sc[16];
    out[8] = sc[5] - s.ABz * sc[31];
    sc[5] = s.PAx * sc[17] + s.WPx * sc[37];
    sc[37] = sc[5] - s.CDy * sc[36];
    out[19] = sc[37] - s.ABx * sc[18];
    out[7] = sc[37] - s.ABz * sc[33];
    sc[37] = s.PAx * sc[19] + s.WPx * sc[38] + 1.0 * s.inv_2zpq * sc[9];
    sc[38] = sc[37] - s.CDx * sc[36];
    out[18] = sc[38] - s.ABx * sc[20];
    out[6] = sc[38] - s.ABz * sc[35];
    sc[38] = s.PAx * sc[22] + s.WPx * sc[39];
    sc[39] = sc[38] - s.CDz * sc[3];
    out[11] = sc[39] - s.ABx * sc[23];
    out[5] = sc[39] - s.ABy * sc[31];
    sc[39] = s.PAx * sc[24] + s.WPx * sc[40];
    sc[40] = sc[39] - s.CDy * sc[3];
    out[10] = sc[40] - s.ABx * sc[25];
    out[4] = sc[40] - s.ABy * sc[33];
    sc[40] = s.PAx * sc[26] + s.WPx * sc[41] + 1.0 * s.inv_2zpq * sc[10];
    sc[41] = sc[40] - s.CDx * sc[3];
    out[9] = sc[41] - s.ABx * sc[27];
    out[3] = sc[41] - s.ABy * sc[35];
    sc[41] = s.PAx * sc[30] + s.WPx * sc[42] + 1.0 * s.inv_2zp * sc[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[6];
    sc[42] = sc[41] - s.CDz * sc[4];
    out[2] = sc[42] - s.ABx * sc[31];
    sc[42] = s.PAx * sc[32] + s.WPx * sc[12] + 1.0 * s.inv_2zp * sc[1] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[7];
    sc[12] = sc[42] - s.CDy * sc[4];
    out[1] = sc[12] - s.ABx * sc[33];
    sc[12] = s.PAx * sc[34] + s.WPx * sc[13] + 1.0 * s.inv_2zp * sc[2] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[8] + 1.0 * s.inv_2zpq * sc[11];
    sc[13] = sc[12] - s.CDx * sc[4];
    out[0] = sc[13] - s.ABx * sc[35];
}

// ERI class (ppps) [array]: 27 Cartesian integrals, 79 DAG nodes, peak liveness 28 slots
__attribute__((noinline)) void recursum_eri_ppps(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[28];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    sc[1] = s.QCy * kf[0] + s.WQy * kf[1];
    sc[2] = s.QCx * kf[0] + s.WQx * kf[1];
    sc[3] = s.QCz * kf[1] + s.WQz * kf[2];
    sc[4] = s.QCy * kf[1] + s.WQy * kf[2];
    sc[5] = s.QCx * kf[1] + s.WQx * kf[2];
    sc[6] = s.PAz * kf[1] + s.WPz * kf[2];
    sc[7] = s.PAy * kf[1] + s.WPy * kf[2];
    sc[8] = s.PAx * kf[1] + s.WPx * kf[2];
    sc[9] = s.QCz * kf[2] + s.WQz * kf[3];
    sc[10] = s.QCy * kf[2] + s.WQy * kf[3];
    sc[11] = s.QCx * kf[2] + s.WQx * kf[3];
    sc[12] = s.PAz * sc[0] + s.WPz * sc[3] + 1.0 * s.inv_2zpq * kf[1];
    sc[13] = s.PAz * sc[1] + s.WPz * sc[4];
    sc[14] = s.PAz * sc[2] + s.WPz * sc[5];
    sc[15] = s.PAy * sc[0] + s.WPy * sc[3];
    sc[16] = s.PAy * sc[1] + s.WPy * sc[4] + 1.0 * s.inv_2zpq * kf[1];
    sc[17] = s.PAy * sc[2] + s.WPy * sc[5];
    sc[18] = s.PAx * sc[0] + s.WPx * sc[3];
    sc[19] = s.PAx * sc[1] + s.WPx * sc[4];
    sc[20] = s.PAx * sc[2] + s.WPx * sc[5] + 1.0 * s.inv_2zpq * kf[1];
    sc[21] = s.PAz * sc[3] + s.WPz * sc[9] + 1.0 * s.inv_2zpq * kf[2];
    sc[22] = s.PAz * sc[4] + s.WPz * sc[10];
    sc[23] = s.PAz * sc[5] + s.WPz * sc[11];
    sc[24] = s.PAy * sc[3] + s.WPy * sc[9];
    sc[25] = s.PAy * sc[4] + s.WPy * sc[10] + 1.0 * s.inv_2zpq * kf[2];
    sc[26] = s.PAy * sc[5] + s.WPy * sc[11];
    sc[27] = s.PAx * sc[3] + s.WPx * sc[9];
    sc[9] = s.PAx * sc[4] + s.WPx * sc[10];
    sc[10] = s.PAx * sc[5] + s.WPx * sc[11] + 1.0 * s.inv_2zpq * kf[2];
    sc[11] = s.PAz * sc[12] + s.WPz * sc[21] + 1.0 * s.inv_2zp * sc[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[3] + 1.0 * s.inv_2zpq * sc[6];
    out[26] = sc[11] - s.ABz * sc[12];
    sc[11] = s.PAz * sc[13] + s.WPz * sc[22] + 1.0 * s.inv_2zp * sc[1] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[4];
    out[25] = sc[11] - s.ABz * sc[13];
    sc[11] = s.PAz * sc[14] + s.WPz * sc[23] + 1.0 * s.inv_2zp * sc[2] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[5];
    out[24] = sc[11] - s.ABz * sc[14];
    sc[11] = s.PAy * sc[12] + s.WPy * sc[21];
    out[23] = sc[11] - s.ABy * sc[12];
    out[17] = sc[11] - s.ABz * sc[15];
    sc[11] = s.PAy * sc[13] + s.WPy * sc[22] + 1.0 * s.inv_2zpq * sc[6];
    out[22] = sc[11] - s.ABy * sc[13];
    out[16] = sc[11] - s.ABz * sc[16];
    sc[11] = s.PAy * sc[14] + s.WPy * sc[23];
    out[21] = sc[11] - s.ABy * sc[14];
    out[15] = sc[11] - s.ABz * sc[17];
    sc[11] = s.PAy * sc[15] + s.WPy * sc[24] + 1.0 * s.inv_2zp * sc[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[3];
    out[14] = sc[11] - s.ABy * sc[15];
    sc[11] = s.PAy * sc[16] + s.WPy * sc[25] + 1.0 * s.inv_2zp * sc[1] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[4] + 1.0 * s.inv_2zpq * sc[7];
    out[13] = sc[11] - s.ABy * sc[16];
    sc[11] = s.PAy * sc[17] + s.WPy * sc[26] + 1.0 * s.inv_2zp * sc[2] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[5];
    out[12] = sc[11] - s.ABy * sc[17];
    sc[11] = s.PAx * sc[12] + s.WPx * sc[21];
    out[20] = sc[11] - s.ABx * sc[12];
    out[8] = sc[11] - s.ABz * sc[18];
    sc[11] = s.PAx * sc[13] + s.WPx * sc[22];
    out[19] = sc[11] - s.ABx * sc[13];
    out[7] = sc[11] - s.ABz * sc[19];
    sc[11] = s.PAx * sc[14] + s.WPx * sc[23] + 1.0 * s.inv_2zpq * sc[6];
    out[18] = sc[11] - s.ABx * sc[14];
    out[6] = sc[11] - s.ABz * sc[20];
    sc[11] = s.PAx * sc[15] + s.WPx * sc[24];
    out[11] = sc[11] - s.ABx * sc[15];
    out[5] = sc[11] - s.ABy * sc[18];
    sc[11] = s.PAx * sc[16] + s.WPx * sc[25];
    out[10] = sc[11] - s.ABx * sc[16];
    out[4] = sc[11] - s.ABy * sc[19];
    sc[11] = s.PAx * sc[17] + s.WPx * sc[26] + 1.0 * s.inv_2zpq * sc[7];
    out[9] = sc[11] - s.ABx * sc[17];
    out[3] = sc[11] - s.ABy * sc[20];
    sc[11] = s.PAx * sc[18] + s.WPx * sc[27] + 1.0 * s.inv_2zp * sc[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[3];
    out[2] = sc[11] - s.ABx * sc[18];
    sc[11] = s.PAx * sc[19] + s.WPx * sc[9] + 1.0 * s.inv_2zp * sc[1] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[4];
    out[1] = sc[11] - s.ABx * sc[19];
    sc[11] = s.PAx * sc[20] + s.WPx * sc[10] + 1.0 * s.inv_2zp * sc[2] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[5] + 1.0 * s.inv_2zpq * sc[8];
    out[0] = sc[11] - s.ABx * sc[20];
}

// ERI class (pppp) [array]: 81 Cartesian integrals, 308 DAG nodes, peak liveness 106 slots
__attribute__((noinline)) void recursum_eri_pppp(
    const ScalarPack& __restrict__ s,
    const double* __restrict__ kf,   // kf[m] = K * F_m(T)
    double* __restrict__ out) {
    double sc[106];
    sc[0] = s.QCz * kf[0] + s.WQz * kf[1];
    sc[1] = s.QCy * kf[0] + s.WQy * kf[1];
    sc[2] = s.QCx * kf[0] + s.WQx * kf[1];
    sc[3] = s.QCz * kf[1] + s.WQz * kf[2];
    sc[4] = s.QCy * kf[1] + s.WQy * kf[2];
    sc[5] = s.QCx * kf[1] + s.WQx * kf[2];
    sc[6] = s.PAz * kf[1] + s.WPz * kf[2];
    sc[7] = s.PAy * kf[1] + s.WPy * kf[2];
    sc[8] = s.PAx * kf[1] + s.WPx * kf[2];
    sc[9] = s.QCz * kf[2] + s.WQz * kf[3];
    sc[10] = s.QCy * kf[2] + s.WQy * kf[3];
    sc[11] = s.QCx * kf[2] + s.WQx * kf[3];
    sc[12] = s.QCz * kf[3] + s.WQz * kf[4];
    sc[13] = s.QCy * kf[3] + s.WQy * kf[4];
    sc[14] = s.QCx * kf[3] + s.WQx * kf[4];
    sc[15] = s.QCz * sc[0] + s.WQz * sc[3] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    sc[16] = s.QCy * sc[0] + s.WQy * sc[3];
    sc[17] = s.QCy * sc[1] + s.WQy * sc[4] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    sc[18] = s.QCx * sc[0] + s.WQx * sc[3];
    sc[19] = s.QCx * sc[1] + s.WQx * sc[4];
    sc[20] = s.QCx * sc[2] + s.WQx * sc[5] + 1.0 * s.inv_2zq * kf[0] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[1];
    sc[21] = s.PAz * sc[0] + s.WPz * sc[3] + 1.0 * s.inv_2zpq * kf[1];
    sc[22] = s.PAz * sc[1] + s.WPz * sc[4];
    sc[23] = s.PAz * sc[2] + s.WPz * sc[5];
    sc[24] = s.PAy * sc[0] + s.WPy * sc[3];
    sc[25] = s.PAy * sc[1] + s.WPy * sc[4] + 1.0 * s.inv_2zpq * kf[1];
    sc[26] = s.PAy * sc[2] + s.WPy * sc[5];
    sc[27] = s.PAx * sc[0] + s.WPx * sc[3];
    sc[28] = s.PAx * sc[1] + s.WPx * sc[4];
    sc[29] = s.PAx * sc[2] + s.WPx * sc[5] + 1.0 * s.inv_2zpq * kf[1];
    sc[30] = s.QCz * sc[3] + s.WQz * sc[9] + 1.0 * s.inv_2zq * kf[1] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[2];
    sc[31] = s.QCy * sc[3] + s.WQy * sc[9];
    sc[32] = s.QCy * sc[4] + s.WQy * sc[10] + 1.0 * s.inv_2zq * kf[1] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[2];
    sc[33] = s.QCx * sc[3] + s.WQx * sc[9];
    sc[34] = s.QCx * sc[4] + s.WQx * sc[10];
    sc[35] = s.QCx * sc[5] + s.WQx * sc[11] + 1.0 * s.inv_2zq * kf[1] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[2];
    sc[36] = s.PAz * sc[3] + s.WPz * sc[9] + 1.0 * s.inv_2zpq * kf[2];
    sc[37] = s.PAz * sc[4] + s.WPz * sc[10];
    sc[38] = s.PAz * sc[5] + s.WPz * sc[11];
    sc[39] = s.PAy * sc[3] + s.WPy * sc[9];
    sc[40] = s.PAy * sc[4] + s.WPy * sc[10] + 1.0 * s.inv_2zpq * kf[2];
    sc[41] = s.PAy * sc[5] + s.WPy * sc[11];
    sc[42] = s.PAx * sc[3] + s.WPx * sc[9];
    sc[43] = s.PAx * sc[4] + s.WPx * sc[10];
    sc[44] = s.PAx * sc[5] + s.WPx * sc[11] + 1.0 * s.inv_2zpq * kf[2];
    sc[45] = s.QCz * sc[9] + s.WQz * sc[12] + 1.0 * s.inv_2zq * kf[2] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[3];
    sc[46] = s.QCy * sc[9] + s.WQy * sc[12];
    sc[47] = s.QCy * sc[10] + s.WQy * sc[13] + 1.0 * s.inv_2zq * kf[2] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[3];
    sc[48] = s.QCx * sc[9] + s.WQx * sc[12];
    sc[12] = s.QCx * sc[10] + s.WQx * sc[13];
    sc[13] = s.QCx * sc[11] + s.WQx * sc[14] + 1.0 * s.inv_2zq * kf[2] - (1.0 * s.inv_2zq) * s.frac_p_over_pq * kf[3];
    sc[14] = s.PAz * sc[15] + s.WPz * sc[30] + 2.0 * s.inv_2zpq * sc[3];
    sc[49] = sc[14] - s.CDz * sc[21];
    sc[50] = s.PAz * sc[16] + s.WPz * sc[31] + 1.0 * s.inv_2zpq * sc[4];
    sc[51] = sc[50] - s.CDy * sc[21];
    sc[52] = sc[50] - s.CDz * sc[22];
    sc[53] = s.PAz * sc[17] + s.WPz * sc[32];
    sc[54] = sc[53] - s.CDy * sc[22];
    sc[55] = s.PAz * sc[18] + s.WPz * sc[33] + 1.0 * s.inv_2zpq * sc[5];
    sc[56] = sc[55] - s.CDx * sc[21];
    sc[57] = sc[55] - s.CDz * sc[23];
    sc[58] = s.PAz * sc[19] + s.WPz * sc[34];
    sc[59] = sc[58] - s.CDx * sc[22];
    sc[60] = sc[58] - s.CDy * sc[23];
    sc[61] = s.PAz * sc[20] + s.WPz * sc[35];
    sc[62] = sc[61] - s.CDx * sc[23];
    sc[63] = s.PAz * sc[21] + s.WPz * sc[36] + 1.0 * s.inv_2zp * sc[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[3] + 1.0 * s.inv_2zpq * sc[6];
    sc[64] = s.PAz * sc[22] + s.WPz * sc[37] + 1.0 * s.inv_2zp * sc[1] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[4];
    sc[65] = s.PAz * sc[23] + s.WPz * sc[38] + 1.0 * s.inv_2zp * sc[2] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[5];
    sc[66] = s.PAy * sc[15] + s.WPy * sc[30];
    sc[67] = sc[66] - s.CDz * sc[24];
    sc[68] = s.PAy * sc[16] + s.WPy * sc[31] + 1.0 * s.inv_2zpq * sc[3];
    sc[69] = sc[68] - s.CDy * sc[24];
    sc[70] = sc[68] - s.CDz * sc[25];
    sc[71] = s.PAy * sc[17] + s.WPy * sc[32] + 2.0 * s.inv_2zpq * sc[4];
    sc[72] = sc[71] - s.CDy * sc[25];
    sc[73] = s.PAy * sc[18] + s.WPy * sc[33];
    sc[74] = sc[73] - s.CDx * sc[24];
    sc[75] = sc[73] - s.CDz * sc[26];
    sc[76] = s.PAy * sc[19] + s.WPy * sc[34] + 1.0 * s.inv_2zpq * sc[5];
    sc[77] = sc[76] - s.CDx * sc[25];
    sc[78] = sc[76] - s.CDy * sc[26];
    sc[79] = s.PAy * sc[20] + s.WPy * sc[35];
    sc[80] = sc[79] - s.CDx * sc[26];
    sc[81] = s.PAy * sc[21] + s.WPy * sc[36];
    sc[82] = s.PAy * sc[22] + s.WPy * sc[37] + 1.0 * s.inv_2zpq * sc[6];
    sc[83] = s.PAy * sc[23] + s.WPy * sc[38];
    sc[84] = s.PAy * sc[24] + s.WPy * sc[39] + 1.0 * s.inv_2zp * sc[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[3];
    sc[85] = s.PAy * sc[25] + s.WPy * sc[40] + 1.0 * s.inv_2zp * sc[1] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[4] + 1.0 * s.inv_2zpq * sc[7];
    sc[86] = s.PAy * sc[26] + s.WPy * sc[41] + 1.0 * s.inv_2zp * sc[2] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[5];
    sc[87] = s.PAx * sc[15] + s.WPx * sc[30];
    sc[88] = sc[87] - s.CDz * sc[27];
    sc[89] = s.PAx * sc[16] + s.WPx * sc[31];
    sc[90] = sc[89] - s.CDy * sc[27];
    sc[91] = sc[89] - s.CDz * sc[28];
    sc[92] = s.PAx * sc[17] + s.WPx * sc[32];
    sc[93] = sc[92] - s.CDy * sc[28];
    sc[94] = s.PAx * sc[18] + s.WPx * sc[33] + 1.0 * s.inv_2zpq * sc[3];
    sc[95] = sc[94] - s.CDx * sc[27];
    sc[96] = sc[94] - s.CDz * sc[29];
    sc[97] = s.PAx * sc[19] + s.WPx * sc[34] + 1.0 * s.inv_2zpq * sc[4];
    sc[98] = sc[97] - s.CDx * sc[28];
    sc[99] = sc[97] - s.CDy * sc[29];
    sc[100] = s.PAx * sc[20] + s.WPx * sc[35] + 2.0 * s.inv_2zpq * sc[5];
    sc[101] = sc[100] - s.CDx * sc[29];
    sc[102] = s.PAx * sc[21] + s.WPx * sc[36];
    sc[21] = s.PAx * sc[22] + s.WPx * sc[37];
    sc[22] = s.PAx * sc[23] + s.WPx * sc[38] + 1.0 * s.inv_2zpq * sc[6];
    sc[6] = s.PAx * sc[24] + s.WPx * sc[39];
    sc[24] = s.PAx * sc[25] + s.WPx * sc[40];
    sc[25] = s.PAx * sc[26] + s.WPx * sc[41] + 1.0 * s.inv_2zpq * sc[7];
    sc[7] = s.PAx * sc[27] + s.WPx * sc[42] + 1.0 * s.inv_2zp * sc[0] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[3];
    sc[27] = s.PAx * sc[28] + s.WPx * sc[43] + 1.0 * s.inv_2zp * sc[1] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[4];
    sc[28] = s.PAx * sc[29] + s.WPx * sc[44] + 1.0 * s.inv_2zp * sc[2] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[5] + 1.0 * s.inv_2zpq * sc[8];
    sc[8] = s.PAz * sc[30] + s.WPz * sc[45] + 2.0 * s.inv_2zpq * sc[9];
    sc[29] = s.PAz * sc[31] + s.WPz * sc[46] + 1.0 * s.inv_2zpq * sc[10];
    sc[5] = s.PAz * sc[32] + s.WPz * sc[47];
    sc[2] = s.PAz * sc[33] + s.WPz * sc[48] + 1.0 * s.inv_2zpq * sc[11];
    sc[4] = s.PAz * sc[34] + s.WPz * sc[12];
    sc[1] = s.PAz * sc[35] + s.WPz * sc[13];
    sc[3] = s.PAy * sc[30] + s.WPy * sc[45];
    sc[0] = s.PAy * sc[31] + s.WPy * sc[46] + 1.0 * s.inv_2zpq * sc[9];
    sc[26] = s.PAy * sc[32] + s.WPy * sc[47] + 2.0 * s.inv_2zpq * sc[10];
    sc[23] = s.PAy * sc[33] + s.WPy * sc[48];
    sc[103] = s.PAy * sc[34] + s.WPy * sc[12] + 1.0 * s.inv_2zpq * sc[11];
    sc[104] = s.PAy * sc[35] + s.WPy * sc[13];
    sc[105] = s.PAx * sc[30] + s.WPx * sc[45];
    sc[45] = s.PAx * sc[31] + s.WPx * sc[46];
    sc[46] = s.PAx * sc[32] + s.WPx * sc[47];
    sc[47] = s.PAx * sc[33] + s.WPx * sc[48] + 1.0 * s.inv_2zpq * sc[9];
    sc[48] = s.PAx * sc[34] + s.WPx * sc[12] + 1.0 * s.inv_2zpq * sc[10];
    sc[12] = s.PAx * sc[35] + s.WPx * sc[13] + 2.0 * s.inv_2zpq * sc[11];
    sc[13] = s.PAz * sc[14] + s.WPz * sc[8] + 1.0 * s.inv_2zp * sc[15] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[30] + 2.0 * s.inv_2zpq * sc[36];
    sc[11] = sc[13] - s.CDz * sc[63];
    out[80] = sc[11] - s.ABz * sc[49];
    sc[11] = s.PAz * sc[50] + s.WPz * sc[29] + 1.0 * s.inv_2zp * sc[16] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[31] + 1.0 * s.inv_2zpq * sc[37];
    sc[13] = sc[11] - s.CDy * sc[63];
    out[79] = sc[13] - s.ABz * sc[51];
    sc[13] = sc[11] - s.CDz * sc[64];
    out[77] = sc[13] - s.ABz * sc[52];
    sc[13] = s.PAz * sc[53] + s.WPz * sc[5] + 1.0 * s.inv_2zp * sc[17] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[32];
    sc[11] = sc[13] - s.CDy * sc[64];
    out[76] = sc[11] - s.ABz * sc[54];
    sc[11] = s.PAz * sc[55] + s.WPz * sc[2] + 1.0 * s.inv_2zp * sc[18] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[33] + 1.0 * s.inv_2zpq * sc[38];
    sc[13] = sc[11] - s.CDx * sc[63];
    out[78] = sc[13] - s.ABz * sc[56];
    sc[13] = sc[11] - s.CDz * sc[65];
    out[74] = sc[13] - s.ABz * sc[57];
    sc[13] = s.PAz * sc[58] + s.WPz * sc[4] + 1.0 * s.inv_2zp * sc[19] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[34];
    sc[11] = sc[13] - s.CDx * sc[64];
    out[75] = sc[11] - s.ABz * sc[59];
    sc[11] = sc[13] - s.CDy * sc[65];
    out[73] = sc[11] - s.ABz * sc[60];
    sc[11] = s.PAz * sc[61] + s.WPz * sc[1] + 1.0 * s.inv_2zp * sc[20] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[35];
    sc[13] = sc[11] - s.CDx * sc[65];
    out[72] = sc[13] - s.ABz * sc[62];
    sc[13] = s.PAy * sc[14] + s.WPy * sc[8];
    sc[11] = sc[13] - s.CDz * sc[81];
    out[71] = sc[11] - s.ABy * sc[49];
    out[53] = sc[11] - s.ABz * sc[67];
    sc[11] = s.PAy * sc[50] + s.WPy * sc[29] + 1.0 * s.inv_2zpq * sc[36];
    sc[13] = sc[11] - s.CDy * sc[81];
    out[70] = sc[13] - s.ABy * sc[51];
    out[52] = sc[13] - s.ABz * sc[69];
    sc[13] = sc[11] - s.CDz * sc[82];
    out[68] = sc[13] - s.ABy * sc[52];
    out[50] = sc[13] - s.ABz * sc[70];
    sc[13] = s.PAy * sc[53] + s.WPy * sc[5] + 2.0 * s.inv_2zpq * sc[37];
    sc[11] = sc[13] - s.CDy * sc[82];
    out[67] = sc[11] - s.ABy * sc[54];
    out[49] = sc[11] - s.ABz * sc[72];
    sc[11] = s.PAy * sc[55] + s.WPy * sc[2];
    sc[13] = sc[11] - s.CDx * sc[81];
    out[69] = sc[13] - s.ABy * sc[56];
    out[51] = sc[13] - s.ABz * sc[74];
    sc[13] = sc[11] - s.CDz * sc[83];
    out[65] = sc[13] - s.ABy * sc[57];
    out[47] = sc[13] - s.ABz * sc[75];
    sc[13] = s.PAy * sc[58] + s.WPy * sc[4] + 1.0 * s.inv_2zpq * sc[38];
    sc[11] = sc[13] - s.CDx * sc[82];
    out[66] = sc[11] - s.ABy * sc[59];
    out[48] = sc[11] - s.ABz * sc[77];
    sc[11] = sc[13] - s.CDy * sc[83];
    out[64] = sc[11] - s.ABy * sc[60];
    out[46] = sc[11] - s.ABz * sc[78];
    sc[11] = s.PAy * sc[61] + s.WPy * sc[1];
    sc[13] = sc[11] - s.CDx * sc[83];
    out[63] = sc[13] - s.ABy * sc[62];
    out[45] = sc[13] - s.ABz * sc[80];
    sc[13] = s.PAy * sc[66] + s.WPy * sc[3] + 1.0 * s.inv_2zp * sc[15] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[30];
    sc[11] = sc[13] - s.CDz * sc[84];
    out[44] = sc[11] - s.ABy * sc[67];
    sc[11] = s.PAy * sc[68] + s.WPy * sc[0] + 1.0 * s.inv_2zp * sc[16] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[31] + 1.0 * s.inv_2zpq * sc[39];
    sc[13] = sc[11] - s.CDy * sc[84];
    out[43] = sc[13] - s.ABy * sc[69];
    sc[13] = sc[11] - s.CDz * sc[85];
    out[41] = sc[13] - s.ABy * sc[70];
    sc[13] = s.PAy * sc[71] + s.WPy * sc[26] + 1.0 * s.inv_2zp * sc[17] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[32] + 2.0 * s.inv_2zpq * sc[40];
    sc[11] = sc[13] - s.CDy * sc[85];
    out[40] = sc[11] - s.ABy * sc[72];
    sc[11] = s.PAy * sc[73] + s.WPy * sc[23] + 1.0 * s.inv_2zp * sc[18] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[33];
    sc[13] = sc[11] - s.CDx * sc[84];
    out[42] = sc[13] - s.ABy * sc[74];
    sc[13] = sc[11] - s.CDz * sc[86];
    out[38] = sc[13] - s.ABy * sc[75];
    sc[13] = s.PAy * sc[76] + s.WPy * sc[103] + 1.0 * s.inv_2zp * sc[19] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[34] + 1.0 * s.inv_2zpq * sc[41];
    sc[11] = sc[13] - s.CDx * sc[85];
    out[39] = sc[11] - s.ABy * sc[77];
    sc[11] = sc[13] - s.CDy * sc[86];
    out[37] = sc[11] - s.ABy * sc[78];
    sc[11] = s.PAy * sc[79] + s.WPy * sc[104] + 1.0 * s.inv_2zp * sc[20] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[35];
    sc[13] = sc[11] - s.CDx * sc[86];
    out[36] = sc[13] - s.ABy * sc[80];
    sc[13] = s.PAx * sc[14] + s.WPx * sc[8];
    sc[8] = sc[13] - s.CDz * sc[102];
    out[62] = sc[8] - s.ABx * sc[49];
    out[26] = sc[8] - s.ABz * sc[88];
    sc[8] = s.PAx * sc[50] + s.WPx * sc[29];
    sc[29] = sc[8] - s.CDy * sc[102];
    out[61] = sc[29] - s.ABx * sc[51];
    out[25] = sc[29] - s.ABz * sc[90];
    sc[29] = sc[8] - s.CDz * sc[21];
    out[59] = sc[29] - s.ABx * sc[52];
    out[23] = sc[29] - s.ABz * sc[91];
    sc[29] = s.PAx * sc[53] + s.WPx * sc[5];
    sc[5] = sc[29] - s.CDy * sc[21];
    out[58] = sc[5] - s.ABx * sc[54];
    out[22] = sc[5] - s.ABz * sc[93];
    sc[5] = s.PAx * sc[55] + s.WPx * sc[2] + 1.0 * s.inv_2zpq * sc[36];
    sc[2] = sc[5] - s.CDx * sc[102];
    out[60] = sc[2] - s.ABx * sc[56];
    out[24] = sc[2] - s.ABz * sc[95];
    sc[2] = sc[5] - s.CDz * sc[22];
    out[56] = sc[2] - s.ABx * sc[57];
    out[20] = sc[2] - s.ABz * sc[96];
    sc[2] = s.PAx * sc[58] + s.WPx * sc[4] + 1.0 * s.inv_2zpq * sc[37];
    sc[4] = sc[2] - s.CDx * sc[21];
    out[57] = sc[4] - s.ABx * sc[59];
    out[21] = sc[4] - s.ABz * sc[98];
    sc[4] = sc[2] - s.CDy * sc[22];
    out[55] = sc[4] - s.ABx * sc[60];
    out[19] = sc[4] - s.ABz * sc[99];
    sc[4] = s.PAx * sc[61] + s.WPx * sc[1] + 2.0 * s.inv_2zpq * sc[38];
    sc[1] = sc[4] - s.CDx * sc[22];
    out[54] = sc[1] - s.ABx * sc[62];
    out[18] = sc[1] - s.ABz * sc[101];
    sc[1] = s.PAx * sc[66] + s.WPx * sc[3];
    sc[3] = sc[1] - s.CDz * sc[6];
    out[35] = sc[3] - s.ABx * sc[67];
    out[17] = sc[3] - s.ABy * sc[88];
    sc[3] = s.PAx * sc[68] + s.WPx * sc[0];
    sc[0] = sc[3] - s.CDy * sc[6];
    out[34] = sc[0] - s.ABx * sc[69];
    out[16] = sc[0] - s.ABy * sc[90];
    sc[0] = sc[3] - s.CDz * sc[24];
    out[32] = sc[0] - s.ABx * sc[70];
    out[14] = sc[0] - s.ABy * sc[91];
    sc[0] = s.PAx * sc[71] + s.WPx * sc[26];
    sc[26] = sc[0] - s.CDy * sc[24];
    out[31] = sc[26] - s.ABx * sc[72];
    out[13] = sc[26] - s.ABy * sc[93];
    sc[26] = s.PAx * sc[73] + s.WPx * sc[23] + 1.0 * s.inv_2zpq * sc[39];
    sc[23] = sc[26] - s.CDx * sc[6];
    out[33] = sc[23] - s.ABx * sc[74];
    out[15] = sc[23] - s.ABy * sc[95];
    sc[23] = sc[26] - s.CDz * sc[25];
    out[29] = sc[23] - s.ABx * sc[75];
    out[11] = sc[23] - s.ABy * sc[96];
    sc[23] = s.PAx * sc[76] + s.WPx * sc[103] + 1.0 * s.inv_2zpq * sc[40];
    sc[103] = sc[23] - s.CDx * sc[24];
    out[30] = sc[103] - s.ABx * sc[77];
    out[12] = sc[103] - s.ABy * sc[98];
    sc[103] = sc[23] - s.CDy * sc[25];
    out[28] = sc[103] - s.ABx * sc[78];
    out[10] = sc[103] - s.ABy * sc[99];
    sc[103] = s.PAx * sc[79] + s.WPx * sc[104] + 2.0 * s.inv_2zpq * sc[41];
    sc[104] = sc[103] - s.CDx * sc[25];
    out[27] = sc[104] - s.ABx * sc[80];
    out[9] = sc[104] - s.ABy * sc[101];
    sc[104] = s.PAx * sc[87] + s.WPx * sc[105] + 1.0 * s.inv_2zp * sc[15] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[30];
    sc[105] = sc[104] - s.CDz * sc[7];
    out[8] = sc[105] - s.ABx * sc[88];
    sc[105] = s.PAx * sc[89] + s.WPx * sc[45] + 1.0 * s.inv_2zp * sc[16] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[31];
    sc[45] = sc[105] - s.CDy * sc[7];
    out[7] = sc[45] - s.ABx * sc[90];
    sc[45] = sc[105] - s.CDz * sc[27];
    out[5] = sc[45] - s.ABx * sc[91];
    sc[45] = s.PAx * sc[92] + s.WPx * sc[46] + 1.0 * s.inv_2zp * sc[17] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[32];
    sc[46] = sc[45] - s.CDy * sc[27];
    out[4] = sc[46] - s.ABx * sc[93];
    sc[46] = s.PAx * sc[94] + s.WPx * sc[47] + 1.0 * s.inv_2zp * sc[18] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[33] + 1.0 * s.inv_2zpq * sc[42];
    sc[47] = sc[46] - s.CDx * sc[7];
    out[6] = sc[47] - s.ABx * sc[95];
    sc[47] = sc[46] - s.CDz * sc[28];
    out[2] = sc[47] - s.ABx * sc[96];
    sc[47] = s.PAx * sc[97] + s.WPx * sc[48] + 1.0 * s.inv_2zp * sc[19] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[34] + 1.0 * s.inv_2zpq * sc[43];
    sc[48] = sc[47] - s.CDx * sc[27];
    out[3] = sc[48] - s.ABx * sc[98];
    sc[48] = sc[47] - s.CDy * sc[28];
    out[1] = sc[48] - s.ABx * sc[99];
    sc[48] = s.PAx * sc[100] + s.WPx * sc[12] + 1.0 * s.inv_2zp * sc[20] - (1.0 * s.inv_2zp) * s.frac_q_over_pq * sc[35] + 2.0 * s.inv_2zpq * sc[44];
    sc[12] = sc[48] - s.CDx * sc[28];
    out[0] = sc[12] - s.ABx * sc[101];
}