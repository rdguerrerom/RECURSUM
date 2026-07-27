// Numerically validate the emitted one-electron kernels against the
// PySCF-validated DAG blocks (onee_expected.h). Build:
//   g++ -O2 -std=c++17 onee_validate.cpp onee_kernels.cpp -o onee_validate
#include <cstdio>
#include <cmath>
#include "recursum_onee_scalars.h"
#include "onee_decls.h"
#include "onee_expected.h"

static int fails = 0;
static void chk(const char* cls, const char* op, const double* got,
                const double* exp, int n) {
    double e = 0, m = 0;
    for (int i = 0; i < n; i++) { m = std::fmax(m, std::fabs(exp[i]));
        e = std::fmax(e, std::fabs(got[i] - exp[i])); }
    double rel = e / (m > 1e-300 ? m : 1e-300);
    bool ok = rel < 1e-12;
    if (!ok) fails++;
    printf("  %-3s %s  relerr=%.1e  %s\n", cls, op, rel, ok ? "PASS" : "FAIL");
}

#define CHECK(nm, NO) do { \
    double o[NO], v[NO], k[NO]; \
    recursum_ovlp_##nm(TEST_S, TEST_S00, o); \
    recursum_nuc_##nm(TEST_S, TEST_KF, v); \
    recursum_kin_##nm(TEST_S, TEST_S00, k); \
    chk(#nm, "S", o, EXP_ovlp_##nm, NO); \
    chk(#nm, "V", v, EXP_nuc_##nm, NO); \
    chk(#nm, "T", k, EXP_kin_##nm, NO); \
} while (0)

int main() {
    printf("Emitted one-electron kernels vs PySCF-validated DAG:\n");
    CHECK(ss, 1);  CHECK(ps, 3);  CHECK(pp, 9);
    CHECK(ds, 6);  CHECK(dp, 18); CHECK(dd, 36);
    CHECK(fs, 10); CHECK(fp, 30); CHECK(fd, 60); CHECK(ff, 100);
    printf(fails ? "\n%d FAILURES\n" : "\nALL PASS\n", fails);
    return fails ? 1 : 0;
}
