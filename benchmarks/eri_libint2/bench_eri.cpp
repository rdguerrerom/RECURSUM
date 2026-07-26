// Apples-to-apples primitive Cartesian ERI benchmark: RECURSUM (DAG-topological
// emitter + peak-liveness slot reuse) vs libint2 2.8.1 (both built here with
// identical g++ -O3 -march=native -ffast-math, scalar double, VECLEN=1).
//
// Fairness design:
//   * Same fixed primitive quartet (exponents + centers) for both.
//   * Both time ONLY the VRR/HRR recurrence contraction. All per-quartet
//     prerequisites (PA/WP/oo2z scalars, prescaled (ss|ss)^m base integrals)
//     are computed ONCE outside the timed loop -- mirroring how each library is
//     used in a real contraction inner loop.
//   * RECURSUM:  recursum_eri_<cls>(scalars, kf, out)
//   * libint2 :  libint2_build_eri[la][lb][lc][ld](&inteval)
//   * Correctness of both is asserted against each other at startup before any
//     timing (they compute the SAME unnormalized primitive integrals: libint2
//     pfac == oracle K, so results match to ~1e-13).
#include <benchmark/benchmark.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

#include "recursum_eri_scalars.h"
#include "eri_classes.h"   // generated Cls table + recursum_call/dispatch

// ---- libint2 low-level C API + canonical populator -----------------------
#include <libint2.h>
#include <libint2/boys.h>
// prep_libint2.h declares these two extern; define them here.
libint2::FmEval_Chebyshev7<double> fmeval_chebyshev(28);
libint2::FmEval_Taylor<double, 6> fmeval_taylor(1e-15);
#include "prep_libint2.h"   // from libint-2.8.1/tests/eri (copied next to this file)

// ---------- fixed primitive quartet ---------------------------------------
struct Prim { double a, A[3]; };
static const double EXPA = 1.2, EXPB = 0.8, EXPC = 1.5, EXPD = 0.9;
static const double CA[3] = {0.0, 0.0, 0.0};
static const double CB[3] = {0.5, 0.1, -0.2};
static const double CC[3] = {0.1, 0.7, 0.3};
static const double CD[3] = {-0.3, 0.2, 0.9};

// ---------- RECURSUM prerequisite fill (mirror oracle.compute_scalars) -----
static void fill_scalars(ScalarPack& s, std::vector<double>& kf, int max_m) {
    double a = EXPA, b = EXPB, g = EXPC, d = EXPD;
    const double *A = CA, *B = CB, *C = CC, *D = CD;
    double zp = a + b, zq = g + d, zpq = zp + zq, rho = zp * zq / zpq;
    double prefAB = a * b / zp, prefCD = g * d / zq;
    double P[3], Q[3], W[3];
    for (int i = 0; i < 3; ++i) {
        P[i] = (a * A[i] + b * B[i]) / zp;
        Q[i] = (g * C[i] + d * D[i]) / zq;
        W[i] = (zp * P[i] + zq * Q[i]) / zpq;
    }
    s.PAx = P[0]-A[0]; s.PAy = P[1]-A[1]; s.PAz = P[2]-A[2];
    s.QCx = Q[0]-C[0]; s.QCy = Q[1]-C[1]; s.QCz = Q[2]-C[2];
    s.WPx = W[0]-P[0]; s.WPy = W[1]-P[1]; s.WPz = W[2]-P[2];
    s.WQx = W[0]-Q[0]; s.WQy = W[1]-Q[1]; s.WQz = W[2]-Q[2];
    s.ABx = B[0]-A[0]; s.ABy = B[1]-A[1]; s.ABz = B[2]-A[2];   // osx (B-A)
    s.CDx = D[0]-C[0]; s.CDy = D[1]-C[1]; s.CDz = D[2]-C[2];
    s.inv_2zp = 0.5/zp; s.inv_2zq = 0.5/zq; s.inv_2zpq = 0.5/zpq;
    s.frac_q_over_pq = zq/zpq; s.frac_p_over_pq = zp/zpq;

    double R2AB = 0, R2CD = 0, R2PQ = 0;
    for (int i = 0; i < 3; ++i) {
        R2AB += (A[i]-B[i])*(A[i]-B[i]);
        R2CD += (C[i]-D[i])*(C[i]-D[i]);
        R2PQ += (P[i]-Q[i])*(P[i]-Q[i]);
    }
    double geom = 34.986836655249725693 / (zp * zq * std::sqrt(zpq));
    double K = geom * std::exp(-prefAB*R2AB) * std::exp(-prefCD*R2CD);
    double T = rho * R2PQ;
    std::vector<double> F(max_m + 1);
    libint2::FmEval_Reference<double>::eval(F.data(), T, max_m);
    kf.resize(max_m + 1);
    for (int m = 0; m <= max_m; ++m) kf[m] = K * F[m];
}

// ---------- libint2 setup for one class -----------------------------------
static Libint_t* make_libint_eval(int la, int lb, int lc, int ld) {
    uint am[4] = {(uint)la, (uint)lb, (uint)lc, (uint)ld};
    RandomShellSet<4u> rs(am, 1, 1);           // veclen=1, contrdepth=1
    for (int i = 0; i < 3; ++i) {
        rs.R[0][i] = CA[i]; rs.R[1][i] = CB[i];
        rs.R[2][i] = CC[i]; rs.R[3][i] = CD[i];
    }
    rs.exp[0][0][0] = EXPA; rs.exp[1][0][0] = EXPB;
    rs.exp[2][0][0] = EXPC; rs.exp[3][0][0] = EXPD;
    rs.coef[0][0][0] = 1.0; rs.coef[1][0][0] = 1.0;
    rs.coef[2][0][0] = 1.0; rs.coef[3][0][0] = 1.0;

    int lmax = LIBINT2_MAX_AM_eri;
    Libint_t* ie = new Libint_t;
    LIBINT2_PREFIXED_NAME(libint2_init_eri)(ie, lmax, 0);
    prep_libint2(ie, rs, 0, 0);
    ie->contrdepth = 1;
    return ie;
}

// class table (ERI_CLASSES) and recursum_call() come from generated
// eri_classes.h -- libint2-CANONICAL orderings so both libraries compute the
// identical Cartesian integral set with no permutation.
#define CLASSES ERI_CLASSES
static inline void call_recursum(const char* c, const ScalarPack& s,
                                 const double* kf, double* out) {
    recursum_call(c, s, kf, out);
}

// ---------- correctness cross-check at startup ----------------------------
static bool verify_all() {
    bool ok = true;
    for (const auto& cl : CLASSES) {
        int max_m = cl.l[0]+cl.l[1]+cl.l[2]+cl.l[3];
        ScalarPack s; std::vector<double> kf;
        fill_scalars(s, kf, max_m);
        std::vector<double> rout(cl.nout);
        call_recursum(cl.name, s, kf.data(), rout.data());

        Libint_t* ie = make_libint_eval(cl.l[0],cl.l[1],cl.l[2],cl.l[3]);
        LIBINT2_PREFIXED_NAME(libint2_build_eri)[cl.l[0]][cl.l[1]][cl.l[2]][cl.l[3]](ie);
        const double* lout = ie->targets[0];

        double maxrel = 0, maxabs = 0;
        for (int i = 0; i < cl.nout; ++i) {
            double d = std::abs(rout[i] - lout[i]);
            maxabs = std::max(maxabs, d);
            maxrel = std::max(maxrel, d / (std::abs(lout[i]) + 1e-300));
        }
        printf("  %-5s n=%4d  max|Δ|=%.2e  maxrel=%.2e  %s\n",
               cl.name, cl.nout, maxabs, maxrel, maxrel < 1e-9 ? "OK" : "MISMATCH");
        ok &= (maxrel < 1e-9);
        LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(ie); delete ie;
    }
    return ok;
}

// ---------- benchmarks -----------------------------------------------------
static void BM_RECURSUM(benchmark::State& st, const Cls* cl) {
    int max_m = cl->l[0]+cl->l[1]+cl->l[2]+cl->l[3];
    ScalarPack s; std::vector<double> kf;
    fill_scalars(s, kf, max_m);
    std::vector<double> out(cl->nout);
    recursum_fn_t fn = recursum_dispatch(cl->name);  // resolve ONCE, like libint2's fn
    for (auto _ : st) {
        fn(s, kf.data(), out.data());
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }
    st.SetLabel("RECURSUM");
    st.counters["nout"] = cl->nout;
}

static void BM_libint2(benchmark::State& st, const Cls* cl) {
    Libint_t* ie = make_libint_eval(cl->l[0],cl->l[1],cl->l[2],cl->l[3]);
    auto fn = LIBINT2_PREFIXED_NAME(libint2_build_eri)[cl->l[0]][cl->l[1]][cl->l[2]][cl->l[3]];
    for (auto _ : st) {
        fn(ie);
        benchmark::DoNotOptimize(ie->targets[0]);
        benchmark::ClobberMemory();
    }
    st.SetLabel("libint2");
    st.counters["nout"] = cl->nout;
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(ie); delete ie;
}

int main(int argc, char** argv) {
    LIBINT2_PREFIXED_NAME(libint2_static_init)();
    printf("=== correctness cross-check (RECURSUM vs libint2) ===\n");
    if (!verify_all()) { printf("CORRECTNESS FAILED — aborting\n"); return 1; }
    printf("=== all classes agree; starting timing ===\n\n");

    for (const auto& cl : CLASSES) {
        benchmark::RegisterBenchmark((std::string("RECURSUM/")+cl.name).c_str(),
                                     BM_RECURSUM, &cl)->Unit(benchmark::kNanosecond)->MinTime(0.5);
        benchmark::RegisterBenchmark((std::string("libint2/")+cl.name).c_str(),
                                     BM_libint2, &cl)->Unit(benchmark::kNanosecond)->MinTime(0.5);
    }
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();
    return 0;
}
