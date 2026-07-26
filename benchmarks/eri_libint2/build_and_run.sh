#!/usr/bin/env bash
# Build + run the RECURSUM-vs-libint2 primitive ERI benchmark with IDENTICAL
# flags for both sides. Fails loudly if libint2.a is missing.
set -euo pipefail
cd "$(dirname "$0")"

LB=libint_build/libint-2.8.1
LIBINT_A=$(find "$LB/build" -name "libint2.a" | head -1)
BENCH_SRC=$(find ../benchmarks -path "*benchmark-src/include" -type d | head -1)
# prefer the gcc-built benchmark lib (build/), not the icpx one (build_intel/)
BENCH_A=$(find ../benchmarks/build/_deps -name "libbenchmark.a" 2>/dev/null | head -1)
[ -z "$BENCH_A" ] && BENCH_A=$(find ../benchmarks -name "libbenchmark.a" | head -1)

[ -z "$LIBINT_A" ] && { echo "ERROR: libint2.a not built yet"; exit 1; }
[ -z "$BENCH_A" ]  && { echo "ERROR: libbenchmark.a not found"; exit 1; }

FLAGS="-O3 -march=native -ffast-math -funroll-loops -std=c++17 -DNDEBUG"
echo "libint2.a : $LIBINT_A"
echo "benchmark : $BENCH_A"
echo "flags     : $FLAGS"
echo "compiler  : $(g++ --version | head -1)"

g++ $FLAGS \
    -I. -I"$LB/build/include" -I"$LB/include" -I"$LB/include/libint2" -I"$BENCH_SRC" \
    bench_eri.cpp "$LIBINT_A" "$BENCH_A" -lpthread -no-pie \
    -o bench_eri
echo "BUILD OK"

./bench_eri \
    --benchmark_repetitions=20 \
    --benchmark_report_aggregates_only=true \
    --benchmark_out=results_eri.json \
    --benchmark_out_format=json
echo "RESULTS -> results_eri.json"
