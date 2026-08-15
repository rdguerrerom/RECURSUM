// FD-oracle correctness gate for the restored gradctr codegen (plan.md L76,
// SHRIKE_Research). Reads a class name + one primitive quartet's ScalarPack +
// kf[] + w2[]=2*alpha(A,B,C) from stdin, calls the generated
// recursum_gradctr_<cls> kernel, and prints every output section value so a
// Python driver can reconstruct d/dX_c = up(a+1_c) - a_c*dn(a-1_c) and compare
// against central finite difference of the ALREADY-VALIDATED value kernel
// (recursum_call, eri_classes.h) — the same methodology
// docs/RECURSUM_GRADIENT_RECOMMENDATIONS.md §4 used for the original
// (pre-regression) bit-exact-to-3.2e-11 validation.
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include "gradctr_dispatch.h"

int main() {
    char cls[8];
    int nout, mmax;
    if (scanf("%7s %d %d", cls, &nout, &mmax) != 3) { fprintf(stderr, "hdr read fail\n"); return 1; }
    ScalarPack s{};
    double* sp = reinterpret_cast<double*>(&s);
    for (int i = 0; i < 23; ++i) if (scanf("%lf", &sp[i]) != 1) { fprintf(stderr,"scalar read fail\n"); return 1; }
    double kf[64];
    for (int m = 0; m <= mmax; ++m) if (scanf("%lf", &kf[m]) != 1) { fprintf(stderr,"kf read fail\n"); return 1; }
    double w2[3];
    for (int i = 0; i < 3; ++i) if (scanf("%lf", &w2[i]) != 1) { fprintf(stderr,"w2 read fail\n"); return 1; }

    recursum_gradctr_t fn = recursum_gradctr_dispatch(cls);
    if (!fn) { printf("NOFN\n"); return 0; }
    int nscratch = recursum_grad_nscratch(cls);
    std::vector<double> out(nout, 0.0), sc(nscratch, 0.0);
    fn(&s, kf, w2, /*cd=*/1, /*mstride=*/mmax + 1, out.data(), sc.data());
    for (int i = 0; i < nout; ++i) printf("%.17e\n", out[i]);
    return 0;
}
