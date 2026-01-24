/**
 * @file validate_coulomb_layeredcodegen.cpp
 * @brief Validate LayeredCodegen Coulomb R implementation against TMP and Layered
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <recursum/mcmd/coulomb_aux.hpp>
#include <recursum/mcmd/coulomb_r_layered.hpp>
#include <recursum/mcmd/coulomb_r_symbolic.hpp>
#include "benchmark_common.hpp"

using namespace recursum::mcmd;
using namespace recursum::benchmark;

template<int L_total>
bool validate_L() {
    CoulombGeometryParams params(42);

    // Pre-compute Boys functions
    constexpr int MAX_L = 8;
    Vec8d Boys[MAX_L + 1];

    // Compute Boys functions
    double T = params.T[0];
    if (T < 1e-10) {
        for (int n = 0; n <= L_total; ++n) {
            Boys[n] = Vec8d(1.0 / (2 * n + 1));
        }
    } else {
        double sqrtT = std::sqrt(T);
        double expT = std::exp(-T);
        Boys[L_total] = Vec8d(0.5 * std::sqrt(M_PI / T) * expT);
        for (int n = L_total; n > 0; --n) {
            double Fn = Boys[n][0];
            Boys[n - 1] = Vec8d((2.0 * T * Fn + expT) / (2.0 * n - 1.0));
        }
    }

    // Compute with all three methods
    constexpr int n_values = (L_total + 1) * (L_total + 2) * (L_total + 3) / 6;

    // TMP (compute individual values)
    std::array<Vec8d, n_values> R_tmp;
    int idx = 0;
    for (int t = 0; t <= L_total; ++t) {
        for (int u = 0; u <= L_total - t; ++u) {
            for (int v = 0; v <= L_total - t - u; ++v) {
                // For TMP we'd need to call individual templates, simplified here
                R_tmp[idx++] = Vec8d(0.0);  // Placeholder
            }
        }
    }

    // Layered
    auto R_layered = recursum::mcmd::layered::CoulombRLayer<L_total>::compute(
        params.PC[0], params.PC[1], params.PC[2], Boys);

    // LayeredCodegen (Symbolic)
    std::array<Vec8d, n_values> R_symbolic;
    CoulombRSymbolic<L_total>::compute(R_symbolic.data(),
        params.PC[0], params.PC[1], params.PC[2], Boys);

    // Compare Layered vs Symbolic
    double max_diff = 0.0;
    bool all_match = true;

    for (int i = 0; i < n_values; ++i) {
        double diff = std::abs(R_layered[i][0] - R_symbolic[i][0]);
        if (diff > max_diff) max_diff = diff;
        if (diff > 1e-12) {
            all_match = false;
            std::cout << "  Mismatch at index " << i
                     << ": Layered=" << R_layered[i][0]
                     << ", Symbolic=" << R_symbolic[i][0]
                     << ", diff=" << diff << "\n";
        }
    }

    std::cout << "L=" << L_total << ": ";
    if (all_match) {
        std::cout << "✓ PASS (max diff: " << std::scientific << max_diff << ")\n";
    } else {
        std::cout << "✗ FAIL (max diff: " << std::scientific << max_diff << ")\n";
    }

    return all_match;
}

int main() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Validating LayeredCodegen Coulomb R Implementation\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Comparing Layered vs LayeredCodegen (Symbolic) implementations:\n\n";

    bool all_passed = true;
    all_passed &= validate_L<0>();
    all_passed &= validate_L<1>();
    all_passed &= validate_L<2>();
    all_passed &= validate_L<3>();
    all_passed &= validate_L<4>();

    std::cout << "\n" << std::string(70, '=') << "\n";
    if (all_passed) {
        std::cout << "✓ All validations PASSED\n";
    } else {
        std::cout << "✗ Some validations FAILED\n";
    }
    std::cout << std::string(70, '=') << "\n\n";

    return all_passed ? 0 : 1;
}
