#!/usr/bin/env python3
"""
RECURSUM Benchmark Analysis Script

IMPORTANT: All plots and statistics are generated EXCLUSIVELY from real
benchmark data (JSON output from Google Benchmark). No estimates or
theoretical values are used.

Usage:
    python analyze_benchmarks.py results/raw/hermite_e_all_*.json
    python analyze_benchmarks.py results/raw/*.json --output-dir results/figures
"""

import json
import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import statistics


@dataclass
class BenchmarkResult:
    """Single benchmark result from JSON."""
    name: str
    real_time: float  # nanoseconds
    cpu_time: float
    iterations: int
    impl: int  # 0=TMP, 1=Symbolic, 2=Naive, 3=BottomUp
    nA: int = 0
    nB: int = 0
    t: int = 0
    L_total: int = 0
    geometries_per_sec: float = 0.0


class BenchmarkAnalyzer:
    """Analyze benchmark results from JSON files."""

    IMPL_NAMES = {
        0: 'TMP',
        1: 'Symbolic',
        2: 'Naive',
        3: 'BottomUp'
    }

    def __init__(self, json_paths: List[str]):
        self.results: List[BenchmarkResult] = []
        self.context: Dict = {}

        for path in json_paths:
            self._load_json(path)

    def _load_json(self, path: str):
        """Load benchmark results from JSON file."""
        with open(path, 'r') as f:
            data = json.load(f)

        # Store context (system info)
        if 'context' in data:
            self.context.update(data['context'])

        # Parse benchmarks
        for bench in data.get('benchmarks', []):
            # Skip aggregate results (mean, median, stddev)
            if bench.get('aggregate_name'):
                continue

            result = BenchmarkResult(
                name=bench.get('name', ''),
                real_time=bench.get('real_time', 0),
                cpu_time=bench.get('cpu_time', 0),
                iterations=bench.get('iterations', 0),
                impl=int(bench.get('impl', -1)),
                nA=int(bench.get('nA', 0)),
                nB=int(bench.get('nB', 0)),
                t=int(bench.get('t', 0)),
                L_total=int(bench.get('L_total', 0)),
                geometries_per_sec=bench.get('Geometries/s', 0)
            )
            self.results.append(result)

    def get_results_by_impl(self, impl: int) -> List[BenchmarkResult]:
        """Get all results for a specific implementation."""
        return [r for r in self.results if r.impl == impl]

    def get_results_by_coefficient(self, nA: int, nB: int, t: int) -> Dict[int, List[BenchmarkResult]]:
        """Get results grouped by implementation for a specific coefficient."""
        results = {}
        for r in self.results:
            if r.nA == nA and r.nB == nB and r.t == t:
                if r.impl not in results:
                    results[r.impl] = []
                results[r.impl].append(r)
        return results

    def calculate_speedup(self, baseline_impl: int = 1, target_impl: int = 0) -> Dict[Tuple[int, int, int], float]:
        """
        Calculate speedup of target over baseline.

        Returns dict mapping (nA, nB, t) -> speedup factor.
        speedup > 1 means target is faster.
        """
        speedups = {}

        # Group by coefficient
        coeff_results = {}
        for r in self.results:
            key = (r.nA, r.nB, r.t)
            if key not in coeff_results:
                coeff_results[key] = {}
            if r.impl not in coeff_results[key]:
                coeff_results[key][r.impl] = []
            coeff_results[key][r.impl].append(r.real_time)

        # Calculate speedups
        for key, impls in coeff_results.items():
            if baseline_impl in impls and target_impl in impls:
                baseline_time = statistics.mean(impls[baseline_impl])
                target_time = statistics.mean(impls[target_impl])
                if target_time > 0:
                    speedups[key] = baseline_time / target_time

        return speedups

    def print_summary(self):
        """Print summary statistics from real data."""
        print("\n" + "=" * 80)
        print("BENCHMARK SUMMARY (from actual measurements)")
        print("=" * 80)

        if self.context:
            print(f"\nSystem: {self.context.get('host_name', 'Unknown')}")
            print(f"CPU: {self.context.get('num_cpus', '?')} cores")

        print(f"\nTotal benchmark results: {len(self.results)}")

        # Group by implementation
        for impl in sorted(set(r.impl for r in self.results if r.impl >= 0)):
            impl_results = self.get_results_by_impl(impl)
            if impl_results:
                times = [r.real_time for r in impl_results]
                print(f"\n{self.IMPL_NAMES.get(impl, f'Impl {impl}')}:")
                print(f"  Count: {len(impl_results)}")
                print(f"  Mean time: {statistics.mean(times):.2f} ns")
                print(f"  Std dev: {statistics.stdev(times) if len(times) > 1 else 0:.2f} ns")
                print(f"  Min: {min(times):.2f} ns")
                print(f"  Max: {max(times):.2f} ns")

        # TMP vs Symbolic speedup
        speedups = self.calculate_speedup(baseline_impl=1, target_impl=0)
        if speedups:
            values = list(speedups.values())
            print(f"\nTMP vs Symbolic Speedup (from {len(speedups)} coefficients):")
            print(f"  Mean: {statistics.mean(values):.2f}x")
            print(f"  Median: {statistics.median(values):.2f}x")
            print(f"  Range: {min(values):.2f}x - {max(values):.2f}x")

    def export_csv(self, output_path: str):
        """Export results to CSV for external analysis."""
        with open(output_path, 'w') as f:
            f.write("name,impl,nA,nB,t,L_total,real_time_ns,cpu_time_ns,iterations,geometries_per_sec\n")
            for r in self.results:
                f.write(f"{r.name},{r.impl},{r.nA},{r.nB},{r.t},{r.L_total},"
                        f"{r.real_time},{r.cpu_time},{r.iterations},{r.geometries_per_sec}\n")
        print(f"Exported to {output_path}")


def generate_plots(analyzer: BenchmarkAnalyzer, output_dir: str):
    """
    Generate publication-ready plots from REAL benchmark data.

    All plots are based exclusively on actual measurements from the JSON files.
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available. Install with: pip install matplotlib")
        return

    os.makedirs(output_dir, exist_ok=True)

    # Publication-quality settings
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
    })

    COLORS = {
        'TMP': '#0077BB',
        'Symbolic': '#EE7733',
        'Naive': '#33BBEE',
        'BottomUp': '#EE3377'
    }

    # Plot 1: Speedup by coefficient
    speedups = analyzer.calculate_speedup(baseline_impl=1, target_impl=0)
    if speedups:
        fig, ax = plt.subplots(figsize=(8, 5))

        # Sort by L_total = nA + nB
        sorted_keys = sorted(speedups.keys(), key=lambda x: (x[0] + x[1], x[2]))
        labels = [f"E({k[0]},{k[1]},{k[2]})" for k in sorted_keys]
        values = [speedups[k] for k in sorted_keys]

        x = np.arange(len(labels))
        bars = ax.bar(x, values, color=COLORS['TMP'], edgecolor='black', linewidth=0.5)

        ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, label='Equal performance')
        ax.set_xlabel('Hermite E Coefficient')
        ax.set_ylabel('Speedup (TMP / Symbolic)')
        ax.set_title('TMP Speedup over Symbolic (from actual measurements)')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, axis='y', alpha=0.3)

        # Annotate values > 2x or < 0.5x
        for i, v in enumerate(values):
            if v > 2.0 or v < 0.5:
                ax.annotate(f'{v:.1f}x', xy=(i, v), ha='center', va='bottom', fontsize=7)

        plt.tight_layout()
        output_path = os.path.join(output_dir, 'speedup_by_coefficient.png')
        plt.savefig(output_path)
        plt.close()
        print(f"Generated: {output_path}")

    # Plot 2: Scaling with angular momentum
    scaling_results = {}
    for r in analyzer.results:
        if 'Scaling' in r.name and r.L_total >= 0:
            key = (r.impl, r.L_total)
            if key not in scaling_results:
                scaling_results[key] = []
            scaling_results[key].append(r.real_time)

    if scaling_results:
        fig, ax = plt.subplots(figsize=(7, 5))

        for impl in [0, 1]:  # TMP and Symbolic
            L_values = []
            means = []
            stds = []

            for L in range(0, 8, 2):  # L=0,2,4,6
                key = (impl, L)
                if key in scaling_results:
                    times = scaling_results[key]
                    L_values.append(L)
                    means.append(statistics.mean(times))
                    stds.append(statistics.stdev(times) if len(times) > 1 else 0)

            if L_values:
                impl_name = analyzer.IMPL_NAMES.get(impl, f'Impl {impl}')
                ax.errorbar(L_values, means, yerr=stds, marker='o', capsize=5,
                           color=COLORS.get(impl_name, 'gray'), linewidth=2,
                           markersize=8, label=impl_name)

        ax.set_xlabel('Total Angular Momentum ($L_A + L_B$)')
        ax.set_ylabel('Execution Time (ns)')
        ax.set_title('Performance Scaling (from actual measurements)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')

        plt.tight_layout()
        output_path = os.path.join(output_dir, 'scaling_by_L.png')
        plt.savefig(output_path)
        plt.close()
        print(f"Generated: {output_path}")

    # Plot 3: Comparison bar chart
    compare_results = {}
    for r in analyzer.results:
        if 'Compare' in r.name:
            # Extract coefficient from name
            parts = r.name.split('_')
            if len(parts) >= 4:
                try:
                    coeff_key = f"{parts[-3]}_{parts[-2]}_{parts[-1]}"
                    if coeff_key not in compare_results:
                        compare_results[coeff_key] = {}
                    impl_name = analyzer.IMPL_NAMES.get(r.impl, f'Impl{r.impl}')
                    if impl_name not in compare_results[coeff_key]:
                        compare_results[coeff_key][impl_name] = []
                    compare_results[coeff_key][impl_name].append(r.real_time)
                except (ValueError, IndexError):
                    pass

    if compare_results:
        fig, ax = plt.subplots(figsize=(10, 6))

        coefficients = sorted(compare_results.keys())
        x = np.arange(len(coefficients))
        width = 0.35

        tmp_means = []
        sym_means = []
        tmp_stds = []
        sym_stds = []

        for coeff in coefficients:
            tmp_times = compare_results[coeff].get('TMP', [])
            sym_times = compare_results[coeff].get('Symbolic', [])
            tmp_means.append(statistics.mean(tmp_times) if tmp_times else 0)
            sym_means.append(statistics.mean(sym_times) if sym_times else 0)
            tmp_stds.append(statistics.stdev(tmp_times) if len(tmp_times) > 1 else 0)
            sym_stds.append(statistics.stdev(sym_times) if len(sym_times) > 1 else 0)

        ax.bar(x - width/2, tmp_means, width, yerr=tmp_stds, label='TMP',
               color=COLORS['TMP'], edgecolor='black', linewidth=0.5, capsize=3)
        ax.bar(x + width/2, sym_means, width, yerr=sym_stds, label='Symbolic',
               color=COLORS['Symbolic'], edgecolor='black', linewidth=0.5, capsize=3)

        # Annotate speedup
        for i, (tmp, sym) in enumerate(zip(tmp_means, sym_means)):
            if tmp > 0 and sym > 0:
                speedup = sym / tmp
                ax.annotate(f'{speedup:.1f}x', xy=(i, max(tmp, sym) * 1.1),
                           ha='center', fontsize=8, color='green' if speedup > 1 else 'red')

        ax.set_xlabel('Coefficient')
        ax.set_ylabel('Execution Time (ns)')
        ax.set_title('TMP vs Symbolic Comparison (from actual measurements)')
        ax.set_xticks(x)
        ax.set_xticklabels([f'E({c.replace("_", ",")})' for c in coefficients], rotation=45, ha='right')
        ax.legend()
        ax.grid(True, axis='y', alpha=0.3)

        plt.tight_layout()
        output_path = os.path.join(output_dir, 'comparison_bar.png')
        plt.savefig(output_path)
        plt.close()
        print(f"Generated: {output_path}")

    print(f"\nAll plots saved to: {output_dir}")
    print("NOTE: All plots generated from actual benchmark measurements only.")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze RECURSUM benchmark results (from actual JSON data only)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    %(prog)s results/raw/hermite_e_all_*.json
    %(prog)s results/raw/*.json --output-dir results/figures --plot
    %(prog)s results/raw/hermite_e_all_20240115.json --csv results.csv
"""
    )
    parser.add_argument('json_files', nargs='+', help='JSON benchmark result files')
    parser.add_argument('--output-dir', default='../results/figures',
                        help='Output directory for plots')
    parser.add_argument('--plot', action='store_true', help='Generate plots')
    parser.add_argument('--csv', type=str, help='Export results to CSV')
    parser.add_argument('--quiet', action='store_true', help='Suppress summary output')

    args = parser.parse_args()

    # Verify all input files exist
    for path in args.json_files:
        if not os.path.exists(path):
            print(f"Error: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    # Load and analyze
    analyzer = BenchmarkAnalyzer(args.json_files)

    if not args.quiet:
        analyzer.print_summary()

    if args.csv:
        analyzer.export_csv(args.csv)

    if args.plot:
        output_dir = os.path.join(os.path.dirname(__file__), args.output_dir)
        generate_plots(analyzer, output_dir)


if __name__ == '__main__':
    main()
