/**
 * @file test_coulomb_r_layered.cpp
 * @brief Quick validation test for layered Coulomb R implementation
 */

#include <iostream>
#include <cmath>
#include <recursum/mcmd/coulomb_aux.hpp>
#include <recursum/mcmd/coulomb_r_layered.hpp>

using namespace recursum::mcmd;
using namespace recursum::mcmd::layered;

// Simple Boys function approximation for testing
// F_n(T) ≈ (2n-1)!! / (2T)^{n+1/2} for large T
// F_n(T) ≈ 1/(2n+1) for small T
void compute_boys(double T, int n_max, Vec8d* Boys) {
    if (T < 1e-10) {
        for (int n = 0; n <= n_max; ++n) {
            Boys[n] = Vec8d(1.0 / (2 * n + 1));
        }
        return;
    }

    // Use downward recursion from asymptotic approximation
    double sqrtT = std::sqrt(T);
    double expT = std::exp(-T);

    // Asymptotic for large n
    Boys[n_max] = Vec8d(0.5 * std::sqrt(M_PI / T) * expT);

    // Downward recursion: F_{n-1}(T) = (2T * F_n(T) + exp(-T)) / (2n - 1)
    for (int n = n_max; n > 0; --n) {
        double Fn = Boys[n][0];
        Boys[n - 1] = Vec8d((2.0 * T * Fn + expT) / (2.0 * n - 1.0));
    }
}

int main() {
    std::cout << "Testing Coulomb R Layered Implementation\n";
    std::cout << "=========================================\n\n";

    // Test parameters
    Vec8d PCx(0.5), PCy(0.3), PCz(0.2);
    double T = 1.0;  // Example argument

    // Test for different L_total values
    for (int L_total = 0; L_total <= 4; ++L_total) {
        std::cout << "L_total = " << L_total << ":\n";

        // Compute Boys functions
        Vec8d Boys[17];  // Up to L=16
        compute_boys(T, L_total, Boys);

        // Compute using original recursive method
        std::cout << "  Original vs Layered comparison:\n";

        // Only test a few representative cases
        if (L_total >= 0) {
            Vec8d orig_000 = CoulombR<0, 0, 0, 0>::compute(PCx, PCy, PCz, Boys);
            auto R = CoulombRLayer<0>::compute(PCx, PCy, PCz, Boys);
            Vec8d layered_000 = R[r_index(0, 0, 0)];
            double diff = std::abs(orig_000[0] - layered_000[0]);
            std::cout << "    R_{000}: orig=" << orig_000[0] << " layered=" << layered_000[0]
                      << " diff=" << diff << (diff < 1e-10 ? " OK" : " FAIL") << "\n";
        }

        if (L_total >= 1) {
            Vec8d orig_100 = CoulombR<1, 0, 0, 0>::compute(PCx, PCy, PCz, Boys);
            auto R = CoulombRLayer<1>::compute(PCx, PCy, PCz, Boys);
            Vec8d layered_100 = R[r_index(1, 0, 0)];
            double diff = std::abs(orig_100[0] - layered_100[0]);
            std::cout << "    R_{100}: orig=" << orig_100[0] << " layered=" << layered_100[0]
                      << " diff=" << diff << (diff < 1e-10 ? " OK" : " FAIL") << "\n";

            Vec8d orig_010 = CoulombR<0, 1, 0, 0>::compute(PCx, PCy, PCz, Boys);
            Vec8d layered_010 = R[r_index(0, 1, 0)];
            diff = std::abs(orig_010[0] - layered_010[0]);
            std::cout << "    R_{010}: orig=" << orig_010[0] << " layered=" << layered_010[0]
                      << " diff=" << diff << (diff < 1e-10 ? " OK" : " FAIL") << "\n";

            Vec8d orig_001 = CoulombR<0, 0, 1, 0>::compute(PCx, PCy, PCz, Boys);
            Vec8d layered_001 = R[r_index(0, 0, 1)];
            diff = std::abs(orig_001[0] - layered_001[0]);
            std::cout << "    R_{001}: orig=" << orig_001[0] << " layered=" << layered_001[0]
                      << " diff=" << diff << (diff < 1e-10 ? " OK" : " FAIL") << "\n";
        }

        if (L_total >= 2) {
            Vec8d orig_200 = CoulombR<2, 0, 0, 0>::compute(PCx, PCy, PCz, Boys);
            auto R = CoulombRLayer<2>::compute(PCx, PCy, PCz, Boys);
            Vec8d layered_200 = R[r_index(2, 0, 0)];
            double diff = std::abs(orig_200[0] - layered_200[0]);
            std::cout << "    R_{200}: orig=" << orig_200[0] << " layered=" << layered_200[0]
                      << " diff=" << diff << (diff < 1e-10 ? " OK" : " FAIL") << "\n";

            Vec8d orig_110 = CoulombR<1, 1, 0, 0>::compute(PCx, PCy, PCz, Boys);
            Vec8d layered_110 = R[r_index(1, 1, 0)];
            diff = std::abs(orig_110[0] - layered_110[0]);
            std::cout << "    R_{110}: orig=" << orig_110[0] << " layered=" << layered_110[0]
                      << " diff=" << diff << (diff < 1e-10 ? " OK" : " FAIL") << "\n";
        }

        if (L_total >= 3) {
            Vec8d orig_111 = CoulombR<1, 1, 1, 0>::compute(PCx, PCy, PCz, Boys);
            auto R = CoulombRLayer<3>::compute(PCx, PCy, PCz, Boys);
            Vec8d layered_111 = R[r_index(1, 1, 1)];
            double diff = std::abs(orig_111[0] - layered_111[0]);
            std::cout << "    R_{111}: orig=" << orig_111[0] << " layered=" << layered_111[0]
                      << " diff=" << diff << (diff < 1e-10 ? " OK" : " FAIL") << "\n";
        }

        if (L_total >= 4) {
            Vec8d orig_400 = CoulombR<4, 0, 0, 0>::compute(PCx, PCy, PCz, Boys);
            auto R = CoulombRLayer<4>::compute(PCx, PCy, PCz, Boys);
            Vec8d layered_400 = R[r_index(4, 0, 0)];
            double diff = std::abs(orig_400[0] - layered_400[0]);
            std::cout << "    R_{400}: orig=" << orig_400[0] << " layered=" << layered_400[0]
                      << " diff=" << diff << (diff < 1e-10 ? " OK" : " FAIL") << "\n";

            Vec8d orig_220 = CoulombR<2, 2, 0, 0>::compute(PCx, PCy, PCz, Boys);
            Vec8d layered_220 = R[r_index(2, 2, 0)];
            diff = std::abs(orig_220[0] - layered_220[0]);
            std::cout << "    R_{220}: orig=" << orig_220[0] << " layered=" << layered_220[0]
                      << " diff=" << diff << (diff < 1e-10 ? " OK" : " FAIL") << "\n";
        }

        std::cout << "\n";
    }

    // Test r_index function
    std::cout << "Testing r_index function:\n";
    std::cout << "  r_index(0,0,0) = " << r_index(0,0,0) << " (expected 0)\n";
    std::cout << "  r_index(1,0,0) = " << r_index(1,0,0) << " (expected 1)\n";
    std::cout << "  r_index(0,1,0) = " << r_index(0,1,0) << " (expected 2)\n";
    std::cout << "  r_index(0,0,1) = " << r_index(0,0,1) << " (expected 3)\n";
    std::cout << "  r_index(2,0,0) = " << r_index(2,0,0) << " (expected 4)\n";

    std::cout << "\nDone.\n";
    return 0;
}
