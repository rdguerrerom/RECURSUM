#!/bin/bash
# RECURSUM TMP vs Symbolic Benchmark Suite
# Runs benchmarks and collects real data in JSON format

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="${SCRIPT_DIR}/.."
BIN_DIR="${BENCH_DIR}/bin"
RESULTS_DIR="${BENCH_DIR}/results/raw"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p "${RESULTS_DIR}"

print_header() {
    echo -e "${BLUE}==================================================================${NC}"
    echo -e "${BLUE}  $1${NC}"
    echo -e "${BLUE}==================================================================${NC}"
}

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

# System info
print_header "RECURSUM Benchmark Suite: TMP vs Symbolic"
echo ""
print_info "Timestamp: ${TIMESTAMP}"
print_info "CPU: $(lscpu | grep 'Model name' | cut -d: -f2 | xargs)"
print_info "Cores: $(nproc)"
print_info "Memory: $(free -h | grep Mem | awk '{print $2}')"

# Check CPU governor
if [ -f /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor ]; then
    GOVERNOR=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor)
    print_info "CPU Governor: ${GOVERNOR}"
    if [ "${GOVERNOR}" != "performance" ]; then
        print_warn "Set governor to 'performance' for consistent results:"
        print_warn "  sudo cpupower frequency-set -g performance"
    fi
fi

# Save system info
cat > "${RESULTS_DIR}/system_info_${TIMESTAMP}.json" << EOF
{
    "timestamp": "${TIMESTAMP}",
    "cpu": "$(lscpu | grep 'Model name' | cut -d: -f2 | xargs)",
    "cores": $(nproc),
    "governor": "${GOVERNOR:-unknown}",
    "memory": "$(free -h | grep Mem | awk '{print $2}')"
}
EOF

echo ""

# Check if benchmarks exist
if [ ! -f "${BIN_DIR}/bench_hermite_e" ]; then
    echo -e "${RED}[ERROR]${NC} Benchmark binary not found: ${BIN_DIR}/bench_hermite_e"
    echo "Build benchmarks first with:"
    echo "  mkdir -p build && cd build"
    echo "  cmake .. -DRECURSUM_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release"
    echo "  cmake --build . -j\$(nproc)"
    exit 1
fi

run_benchmark() {
    local name=$1
    local exe=$2
    local filter=$3
    local output="${RESULTS_DIR}/${name}_${TIMESTAMP}"

    print_info "Running ${name}..."

    # Run benchmark with JSON output
    ${BIN_DIR}/${exe} \
        --benchmark_filter="${filter}" \
        --benchmark_out="${output}.json" \
        --benchmark_out_format=json \
        --benchmark_repetitions=10 \
        --benchmark_report_aggregates_only=false \
        --benchmark_min_time=1.0 \
        2>&1 | tee "${output}.txt"

    print_info "Results: ${output}.json"
    echo ""
}

# Run benchmark suites
print_header "1. Hermite E Coefficient Comparison (TMP vs Symbolic vs Naive)"
run_benchmark "hermite_e_all" "bench_hermite_e" "BM_"

print_header "2. Direct Comparison Benchmarks"
run_benchmark "hermite_e_compare" "bench_hermite_e" "BM_Compare"

print_header "3. Scaling Analysis"
run_benchmark "hermite_e_scaling" "bench_hermite_e" "BM_Scaling"

print_header "4. Memory Bandwidth Analysis"
run_benchmark "hermite_e_memory" "bench_hermite_e" "BM_.*Memory"

print_header "Benchmark Complete"
print_info "Raw results saved to: ${RESULTS_DIR}"
print_info ""
print_info "To analyze results, run:"
print_info "  python3 ${BENCH_DIR}/analysis/analyze_benchmarks.py ${RESULTS_DIR}/hermite_e_all_${TIMESTAMP}.json"
