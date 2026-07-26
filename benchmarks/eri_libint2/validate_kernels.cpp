// Validation harness: read class + scalars + kf[] from stdin, call the
// generated RECURSUM kernel, print outputs. Python compares to the oracle.
#include <cstdio>
#include <string>
#include "eri_classes.h"   // generated: recursum_call() dispatch

int main() {
    char cls[8];
    int max_m, n_out;
    if (scanf("%7s %d %d", cls, &max_m, &n_out) != 3) return 1;
    ScalarPack s;
    double* sp = reinterpret_cast<double*>(&s);
    for (int i = 0; i < 23; ++i) if(scanf("%lf", &sp[i])!=1) return 1;
    double kf[64];
    for (int m = 0; m <= max_m; ++m) if(scanf("%lf", &kf[m])!=1) return 1;
    double out[16384];
    recursum_call(cls, s, kf, out);
    for (int i = 0; i < n_out; ++i) printf("%.17e\n", out[i]);
    return 0;
}
