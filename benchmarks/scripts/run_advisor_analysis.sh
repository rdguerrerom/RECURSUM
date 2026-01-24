#!/bin/bash
# Intel Advisor Analysis for RECURSUM TMP vs Symbolic Comparison
# Collects real Trip Counts, FLOPs, and Cache metrics

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="${SCRIPT_DIR}/.."
BIN_DIR="${BENCH_DIR}/bin"
OUTPUT_DIR="${BENCH_DIR}/results/advisor"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
ADVISOR_PROJECT="${OUTPUT_DIR}/recursum_advisor_${TIMESTAMP}"

mkdir -p "${OUTPUT_DIR}"

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

# Check Intel oneAPI environment
if ! command -v advixe-cl &> /dev/null; then
    print_warn "Intel Advisor not found. Attempting to source oneAPI..."
    if [ -f /opt/intel/oneapi/setvars.sh ]; then
        source /opt/intel/oneapi/setvars.sh
    else
        echo -e "${RED}[ERROR]${NC} Intel oneAPI not found."
        echo "Install from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/advisor.html"
        exit 1
    fi
fi

# Check benchmark binary
if [ ! -f "${BIN_DIR}/bench_hermite_e" ]; then
    echo -e "${RED}[ERROR]${NC} Benchmark binary not found. Build with -DRECURSUM_BUILD_BENCHMARKS=ON"
    exit 1
fi

print_header "Intel Advisor Analysis: RECURSUM TMP vs Symbolic"
echo ""
print_info "Timestamp: ${TIMESTAMP}"
print_info "Project: ${ADVISOR_PROJECT}"
print_info "Binary: ${BIN_DIR}/bench_hermite_e"
echo ""

# 1. Survey Analysis (Hotspots + FLOP counts)
print_header "1. Survey Analysis (Hotspots + Trip Counts + FLOPs)"
advixe-cl --collect=survey \
    --project-dir="${ADVISOR_PROJECT}" \
    --search-dir src:r="${BENCH_DIR}/.." \
    --stacks \
    --trip-counts \
    --flop \
    -- "${BIN_DIR}/bench_hermite_e" \
        --benchmark_filter="BM_Compare" \
        --benchmark_repetitions=5 \
        --benchmark_min_time=0.5

print_info "Survey analysis complete"
echo ""

# 2. Trip Counts with Cache Simulation
print_header "2. Trip Counts + Cache Simulation"
advixe-cl --collect=tripcounts \
    --project-dir="${ADVISOR_PROJECT}" \
    --flop \
    --stacks \
    --enable-cache-simulation \
    --cachesim-mode=all-levels \
    -- "${BIN_DIR}/bench_hermite_e" \
        --benchmark_filter="BM_Scaling" \
        --benchmark_repetitions=3 \
        --benchmark_min_time=0.5

print_info "Trip counts analysis complete"
echo ""

# 3. Memory Access Pattern Analysis
print_header "3. Memory Access Patterns (MAP)"
advixe-cl --collect=map \
    --project-dir="${ADVISOR_PROJECT}" \
    --stacks \
    -- "${BIN_DIR}/bench_hermite_e" \
        --benchmark_filter="BM_.*Memory" \
        --benchmark_repetitions=3 \
        --benchmark_min_time=0.5

print_info "MAP analysis complete"
echo ""

# 4. Generate Reports
print_header "4. Generating Reports"

# Roofline HTML
advixe-cl --report=roofline \
    --project-dir="${ADVISOR_PROJECT}" \
    --report-output="${OUTPUT_DIR}/roofline_${TIMESTAMP}.html" \
    --format=html 2>/dev/null || print_warn "Roofline report generation failed (may need GUI)"

# Survey CSV (for analysis script)
advixe-cl --report=survey \
    --project-dir="${ADVISOR_PROJECT}" \
    --report-output="${OUTPUT_DIR}/survey_${TIMESTAMP}.csv" \
    --format=csv

# Trip counts CSV
advixe-cl --report=tripcounts \
    --project-dir="${ADVISOR_PROJECT}" \
    --report-output="${OUTPUT_DIR}/tripcounts_${TIMESTAMP}.csv" \
    --format=csv \
    --show-all-columns 2>/dev/null || print_warn "Trip counts report limited"

# Summary
advixe-cl --report=summary \
    --project-dir="${ADVISOR_PROJECT}" \
    --report-output="${OUTPUT_DIR}/summary_${TIMESTAMP}.txt" \
    --format=text

print_header "Analysis Complete"
echo ""
print_info "Results saved to: ${OUTPUT_DIR}"
print_info ""
print_info "Files generated:"
print_info "  - survey_${TIMESTAMP}.csv"
print_info "  - tripcounts_${TIMESTAMP}.csv"
print_info "  - summary_${TIMESTAMP}.txt"
print_info "  - roofline_${TIMESTAMP}.html (if supported)"
print_info ""
print_info "Open GUI for interactive analysis:"
print_info "  advixe-gui ${ADVISOR_PROJECT}"
