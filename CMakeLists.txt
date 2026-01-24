cmake_minimum_required(VERSION 3.18)
project(RECURSUM_Research VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(RECURSUM_BUILD_TESTS "Build C++ tests" OFF)
option(RECURSUM_USE_NATIVE_ARCH "Use -march=native for SIMD" ON)

# Find Python first
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Find pybind11 (prefer Python package over system)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE pybind11_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE pybind11_FIND_RESULT
)

if(pybind11_FIND_RESULT EQUAL 0)
    message(STATUS "Using pybind11 from Python: ${pybind11_DIR}")
    set(pybind11_DIR "${pybind11_DIR}" CACHE PATH "Path to pybind11 CMake files")
endif()

find_package(pybind11 CONFIG REQUIRED)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/generated
)

# Compiler flags for SIMD
if(RECURSUM_USE_NATIVE_ARCH)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(-march=native -ffast-math)
    elseif(MSVC)
        add_compile_options(/arch:AVX2)
    endif()
endif()

# Enable ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    message(STATUS "Using ccache: ${CCACHE_PROGRAM}")
endif()

# Add subdirectory for C++ modules
add_subdirectory(src)

# Benchmarks (optional)
option(RECURSUM_BUILD_BENCHMARKS "Build performance benchmarks" OFF)
if(RECURSUM_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
