/**
 * @file benchmark_common.hpp
 * @brief Common utilities for RECURSUM benchmarks
 *
 * Provides RandomGenerator and GeometryParams for generating realistic
 * SIMD test data with 8 different geometries per Vec8d.
 */

#pragma once

#include <random>
#include <array>
#include <cmath>

#ifndef RECURSUM_VEC_TYPE
#include <recursum/vectorclass.h>
#endif

namespace recursum {
namespace benchmark {

/**
 * @brief Thread-safe random number generator for benchmark data
 */
class RandomGenerator {
public:
    static std::mt19937& engine() {
        thread_local std::mt19937 gen(std::random_device{}());
        return gen;
    }

    /**
     * @brief Generate Vec8d with 8 DIFFERENT random values
     *
     * Critical for proper SIMD benchmarking - each lane should have
     * unique data to avoid artificial cache effects.
     */
    static Vec8d randomVec8d(double min = -2.0, double max = 2.0) {
        std::uniform_real_distribution<double> dist(min, max);
        alignas(64) double data[8];
        for (int i = 0; i < 8; ++i) {
            data[i] = dist(engine());
        }
        Vec8d result;
        result.load(data);
        return result;
    }

    /**
     * @brief Generate random Gaussian exponents (realistic range)
     */
    static Vec8d randomExponent(double min = 0.1, double max = 10.0) {
        return randomVec8d(min, max);
    }

    /**
     * @brief Generate random coordinates (molecular scale)
     */
    static Vec8d randomCoord(double min = -3.0, double max = 3.0) {
        return randomVec8d(min, max);
    }
};

/**
 * @brief Geometry parameters for Hermite coefficient benchmarks
 *
 * Stores 8 different molecular geometries in SIMD vectors:
 * - PA, PB: Displacement vectors from product center to A and B
 * - p: Combined exponent (alpha_A + alpha_B)
 * - one_over_2p: 1/(2p) = aAB in McMD notation
 * - alpha_A, alpha_B: Individual Gaussian exponents
 */
struct GeometryParams {
    std::array<Vec8d, 3> PA;   // P - A displacement (x, y, z)
    std::array<Vec8d, 3> PB;   // P - B displacement (x, y, z)
    Vec8d p;                    // Combined exponent: alpha_A + alpha_B
    Vec8d one_over_2p;          // 1/(2p) = aAB
    Vec8d alpha_A, alpha_B;     // Individual exponents

    // Original centers (for gradient benchmarks)
    std::array<Vec8d, 3> A;
    std::array<Vec8d, 3> B;
    std::array<Vec8d, 3> P;    // Product center

    /**
     * @brief Generate random geometry parameters
     *
     * Creates 8 different geometries with realistic molecular parameters.
     */
    GeometryParams() {
        // Random exponents in realistic range
        alpha_A = RandomGenerator::randomExponent(0.2, 5.0);
        alpha_B = RandomGenerator::randomExponent(0.2, 5.0);
        p = alpha_A + alpha_B;
        one_over_2p = 1.0 / (2.0 * p);

        // Generate 8 different center A positions
        A[0] = RandomGenerator::randomCoord(-1.0, 1.0);
        A[1] = RandomGenerator::randomCoord(-1.0, 1.0);
        A[2] = RandomGenerator::randomCoord(-1.0, 1.0);

        // Generate 8 different center B positions (further spread)
        B[0] = RandomGenerator::randomCoord(-3.0, 3.0);
        B[1] = RandomGenerator::randomCoord(-3.0, 3.0);
        B[2] = RandomGenerator::randomCoord(-3.0, 3.0);

        // Compute product center P = (alpha_A * A + alpha_B * B) / p
        for (int dim = 0; dim < 3; ++dim) {
            P[dim] = (alpha_A * A[dim] + alpha_B * B[dim]) / p;
            PA[dim] = P[dim] - A[dim];
            PB[dim] = P[dim] - B[dim];
        }
    }

    /**
     * @brief Create with specific seed for reproducibility
     */
    explicit GeometryParams(unsigned int seed) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> exp_dist(0.2, 5.0);
        std::uniform_real_distribution<double> coord_dist(-3.0, 3.0);

        alignas(64) double alpha_A_arr[8], alpha_B_arr[8];
        alignas(64) double Ax[8], Ay[8], Az[8], Bx[8], By[8], Bz[8];

        for (int i = 0; i < 8; ++i) {
            alpha_A_arr[i] = exp_dist(gen);
            alpha_B_arr[i] = exp_dist(gen);
            Ax[i] = coord_dist(gen) * 0.5;  // Tighter for A
            Ay[i] = coord_dist(gen) * 0.5;
            Az[i] = coord_dist(gen) * 0.5;
            Bx[i] = coord_dist(gen);
            By[i] = coord_dist(gen);
            Bz[i] = coord_dist(gen);
        }

        alpha_A.load(alpha_A_arr);
        alpha_B.load(alpha_B_arr);
        A[0].load(Ax); A[1].load(Ay); A[2].load(Az);
        B[0].load(Bx); B[1].load(By); B[2].load(Bz);

        p = alpha_A + alpha_B;
        one_over_2p = 1.0 / (2.0 * p);

        for (int dim = 0; dim < 3; ++dim) {
            P[dim] = (alpha_A * A[dim] + alpha_B * B[dim]) / p;
            PA[dim] = P[dim] - A[dim];
            PB[dim] = P[dim] - B[dim];
        }
    }
};

/**
 * @brief Geometry parameters for Coulomb integral benchmarks
 *
 * Extends GeometryParams with PC vector (product center to nucleus C)
 * and Boys function argument T.
 */
struct CoulombGeometryParams : public GeometryParams {
    std::array<Vec8d, 3> C;    // Nucleus position
    std::array<Vec8d, 3> PC;   // P - C displacement
    Vec8d T;                    // Boys function argument: p * |PC|^2

    CoulombGeometryParams() : GeometryParams() {
        // Random nucleus position
        C[0] = RandomGenerator::randomCoord(-2.0, 2.0);
        C[1] = RandomGenerator::randomCoord(-2.0, 2.0);
        C[2] = RandomGenerator::randomCoord(-2.0, 2.0);

        // Compute PC and T
        Vec8d PC2(0.0);
        for (int dim = 0; dim < 3; ++dim) {
            PC[dim] = P[dim] - C[dim];
            PC2 += PC[dim] * PC[dim];
        }
        T = p * PC2;
    }

    explicit CoulombGeometryParams(unsigned int seed) : GeometryParams(seed) {
        std::mt19937 gen(seed + 1000);  // Different seed for C
        std::uniform_real_distribution<double> coord_dist(-2.0, 2.0);

        alignas(64) double Cx[8], Cy[8], Cz[8];
        for (int i = 0; i < 8; ++i) {
            Cx[i] = coord_dist(gen);
            Cy[i] = coord_dist(gen);
            Cz[i] = coord_dist(gen);
        }
        C[0].load(Cx); C[1].load(Cy); C[2].load(Cz);

        Vec8d PC2(0.0);
        for (int dim = 0; dim < 3; ++dim) {
            PC[dim] = P[dim] - C[dim];
            PC2 += PC[dim] * PC[dim];
        }
        T = p * PC2;
    }
};

/**
 * @brief Shell pair information for ERI benchmarks
 */
struct ShellPairInfo {
    int LA, LB;           // Angular momenta
    const char* name;      // Human-readable name (e.g., "ss", "pp", "ff")
    int num_E_coeffs;     // Number of E coefficients
    int num_R_integrals;  // Number of R integrals for this shell type

    constexpr ShellPairInfo(int la, int lb, const char* n, int e, int r)
        : LA(la), LB(lb), name(n), num_E_coeffs(e), num_R_integrals(r) {}
};

// Shell pair definitions for benchmarking
// Coverage: (ss|ss) through (ff|ff)
constexpr ShellPairInfo SHELL_PAIRS[] = {
    {0, 0, "ss", 1,   1},
    {0, 1, "sp", 2,   4},
    {1, 1, "pp", 4,  10},
    {0, 2, "sd", 3,  10},
    {1, 2, "pd", 6,  20},
    {2, 2, "dd", 9,  35},
    {0, 3, "sf", 4,  20},
    {1, 3, "pf", 12, 56},
    {2, 3, "df", 18, 84},
    {3, 3, "ff", 25, 120}
};
constexpr int NUM_SHELL_PAIRS = sizeof(SHELL_PAIRS) / sizeof(SHELL_PAIRS[0]);

/**
 * @brief Number of Cartesian functions for angular momentum L
 */
constexpr int numCartesian(int L) {
    return (L + 1) * (L + 2) / 2;
}

/**
 * @brief Total angular momentum range for shell pair
 */
constexpr int maxT(int LA, int LB) {
    return LA + LB;
}

} // namespace benchmark
} // namespace recursum
