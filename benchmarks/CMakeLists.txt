# benchmarks/CMakeLists.txt
# Google Benchmark integration for RECURSUM TMP vs Symbolic comparison

cmake_minimum_required(VERSION 3.18)
project(recursum_benchmarks VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(BENCH_USE_NATIVE_ARCH "Use -march=native for SIMD" ON)
option(BENCH_ENABLE_ADVISOR "Enable Intel Advisor instrumentation" OFF)

# ==========================================================================
# Google Benchmark (via FetchContent)
# ==========================================================================
include(FetchContent)
FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.3
)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(benchmark)

# ==========================================================================
# Include Directories
# ==========================================================================
set(RECURSUM_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
include_directories(
    ${RECURSUM_ROOT}/include
    ${RECURSUM_ROOT}/src/generated
    ${RECURSUM_ROOT}/external/vectorclass
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/symbolic_generated
)

# ==========================================================================
# Compiler Flags
# ==========================================================================
if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    # Intel oneAPI icpx with aggressive math optimizations
    set(BENCH_FLAGS
        -O3
        -xHost                  # Use best instructions for host CPU
        -fp-model=fast          # Aggressive floating-point model
        -ffast-math
        -funroll-loops
        -fno-omit-frame-pointer
        -g
    )
    message(STATUS "Using Intel oneAPI compiler with aggressive math optimizations")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(BENCH_FLAGS -O3 -ffast-math -funroll-loops -fno-omit-frame-pointer -g)
    if(BENCH_USE_NATIVE_ARCH)
        list(APPEND BENCH_FLAGS -march=native)
    endif()
elseif(MSVC)
    set(BENCH_FLAGS /O2 /arch:AVX2 /fp:fast)
endif()

# ==========================================================================
# Benchmark Executables
# ==========================================================================

# Hermite E Coefficient Benchmarks (TMP vs Symbolic vs Naive)
add_executable(bench_hermite_e
    src/hermite_e/bench_hermite_e_tmp.cpp
    src/hermite_e/bench_hermite_e_symbolic.cpp
    src/hermite_e/bench_hermite_e_naive.cpp
    src/hermite_e/bench_hermite_e_compare.cpp
)
target_link_libraries(bench_hermite_e PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_hermite_e PRIVATE ${BENCH_FLAGS})

# Output to bin directory
set_target_properties(
    bench_hermite_e
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# Comprehensive benchmark (all shell pairs ss to gg)
add_executable(bench_comprehensive
    src/hermite_e/bench_comprehensive.cpp
)
target_link_libraries(bench_comprehensive PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_comprehensive PRIVATE ${BENCH_FLAGS})
set_target_properties(
    bench_comprehensive
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# ALL coefficients benchmark (125 E coefficients + gradients)
add_executable(bench_all_coefficients
    src/hermite_e/bench_all_coefficients.cpp
)
target_link_libraries(bench_all_coefficients PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_all_coefficients PRIVATE ${BENCH_FLAGS})
set_target_properties(
    bench_all_coefficients
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# Direct comparison benchmark (NO dispatcher overhead, full compiler optimization)
add_executable(bench_direct_comparison
    src/hermite_e/bench_direct_comparison.cpp
)
target_link_libraries(bench_direct_comparison PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_direct_comparison PRIVATE ${BENCH_FLAGS})
set_target_properties(
    bench_direct_comparison
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# Layered CSE comparison benchmark (Original vs Layered vs Symbolic)
add_executable(bench_layered_comparison
    src/hermite_e/bench_layered_comparison.cpp
)
target_link_libraries(bench_layered_comparison PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_layered_comparison PRIVATE ${BENCH_FLAGS})
set_target_properties(
    bench_layered_comparison
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# Full layer benchmark (computing ALL t values at once)
add_executable(bench_full_layer
    src/hermite_e/bench_full_layer.cpp
)
target_link_libraries(bench_full_layer PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_full_layer PRIVATE ${BENCH_FLAGS})
set_target_properties(
    bench_full_layer
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# Realistic McMurchie-Davidson benchmark (full layers for all shell pairs)
add_executable(bench_mcmd_realistic
    src/hermite_e/bench_mcmd_realistic.cpp
)
target_link_libraries(bench_mcmd_realistic PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_mcmd_realistic PRIVATE ${BENCH_FLAGS})
set_target_properties(
    bench_mcmd_realistic
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# Test for Coulomb R layered implementation
add_executable(test_coulomb_r_layered
    src/coulomb_r/test_coulomb_r_layered.cpp
)
target_compile_options(test_coulomb_r_layered PRIVATE ${BENCH_FLAGS})
set_target_properties(
    test_coulomb_r_layered
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

add_executable(validate_coulomb_layeredcodegen
    src/coulomb_r/validate_coulomb_layeredcodegen.cpp
)
target_compile_options(validate_coulomb_layeredcodegen PRIVATE ${BENCH_FLAGS})
set_target_properties(
    validate_coulomb_layeredcodegen
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# ==========================================================================
# NEW: Proper ERI-component benchmarks
# ==========================================================================

# Hermite Expansion Coefficients E_N^{i,j} benchmark
add_executable(bench_hermite_coefficients
    src/hermite_e/bench_hermite_coefficients.cpp
)
target_link_libraries(bench_hermite_coefficients PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_hermite_coefficients PRIVATE ${BENCH_FLAGS})
set_target_properties(
    bench_hermite_coefficients
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# Coulomb Auxiliary Integrals R_{tuv} benchmark
add_executable(bench_coulomb_hermite
    src/coulomb_r/bench_coulomb_hermite.cpp
)
target_link_libraries(bench_coulomb_hermite PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_coulomb_hermite PRIVATE ${BENCH_FLAGS})
set_target_properties(
    bench_coulomb_hermite
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# ==========================================================================
# NEW: Layered Codegen Validation
# ==========================================================================

# Test LayeredCppGenerator output (TMP vs Layered Hand-written vs Layered Codegen)
add_executable(bench_layered_codegen_test
    src/hermite_e/bench_layered_codegen_test.cpp
)
target_link_libraries(bench_layered_codegen_test PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_layered_codegen_test PRIVATE ${BENCH_FLAGS})
set_target_properties(
    bench_layered_codegen_test
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# ==========================================================================
# NEW: J/K Matrix Benchmarks for Alkanes
# ==========================================================================

# J and K matrix benchmarks (alkane scaling CH4 through C4H10)
add_executable(bench_jk_alkanes
    src/jk_matrix/bench_jk_alkanes.cpp
)
target_link_libraries(bench_jk_alkanes PRIVATE benchmark::benchmark pthread)
target_compile_options(bench_jk_alkanes PRIVATE ${BENCH_FLAGS})
set_target_properties(
    bench_jk_alkanes
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
)

# ==========================================================================
# Print Configuration Summary
# ==========================================================================
message(STATUS "")
message(STATUS "RECURSUM Benchmark Configuration:")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Native Arch: ${BENCH_USE_NATIVE_ARCH}")
message(STATUS "  Advisor: ${BENCH_ENABLE_ADVISOR}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "")
