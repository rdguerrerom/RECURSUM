#!/usr/bin/env python3
"""
Regenerate all 5 benchmark plots with publication-quality formatting.

This script follows the scientific plotting guidelines and uses:
- Optimized font sizes (font_scale=1.04)
- Log10 scale for Y-axis on performance plots
- Okabe-Ito colorblind-safe palette
- Both PNG (300 DPI) and PDF outputs
- Error bars (±1σ standard deviation)
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from collections import defaultdict

# ============================================================================
# Publication Style Configuration
# ============================================================================

# Okabe-Ito colorblind-safe palette
OKABE_ITO_CYCLE = [
    '#0072B2',  # blue
    '#E69F00',  # orange
    '#009E73',  # green
    '#CC79A7',  # purple (reddish-purple)
    '#D55E00',  # vermillion
    '#56B4E9',  # sky blue
    '#F0E442',  # yellow
    '#000000',  # black
]

def configure_publication_style(font_scale=1.1):
    """Configure matplotlib for publication-quality figures with smaller fonts."""
    plt.rcParams.update({
        # Font configuration
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',

        # Font sizes - SMALLER than before
        'font.size': 9 * font_scale,
        'axes.labelsize': 10 * font_scale,
        'axes.titlesize': 11 * font_scale,
        'xtick.labelsize': 8 * font_scale,
        'ytick.labelsize': 8 * font_scale,
        'legend.fontsize': 8 * font_scale,

        # Axes configuration
        'axes.linewidth': 0.8,
        'axes.labelpad': 3,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.prop_cycle': plt.cycler(color=OKABE_ITO_CYCLE),

        # Tick configuration
        'xtick.major.width': 0.8,
        'xtick.minor.width': 0.6,
        'ytick.major.width': 0.8,
        'ytick.minor.width': 0.6,
        'xtick.major.size': 3.5,
        'xtick.minor.size': 2,
        'ytick.major.size': 3.5,
        'ytick.minor.size': 2,
        'xtick.direction': 'out',
        'ytick.direction': 'out',

        # Legend configuration
        'legend.frameon': False,
        'legend.borderpad': 0.3,
        'legend.handlelength': 1.5,
        'legend.handletextpad': 0.5,
        'legend.columnspacing': 1.0,

        # Line configuration
        'lines.linewidth': 1.5,
        'lines.markersize': 5,

        # Figure configuration
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
    })

def save_figure(fig, filename, output_dir):
    """Save figure in both PDF and PNG formats."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for fmt in ['pdf', 'png']:
        filepath = output_path / f"{filename}.{fmt}"
        fig.savefig(filepath, format=fmt, dpi=300, bbox_inches='tight', pad_inches=0.02)
        print(f"  Saved: {filepath}")

# ============================================================================
# Data Loading and Processing
# ============================================================================

def load_benchmark_data(filepath):
    """Load Google Benchmark JSON data."""
    with open(filepath, 'r') as f:
        data = json.load(f)
    return data['benchmarks']

def parse_hermite_data(benchmarks):
    """
    Parse hermite coefficients benchmark data.

    Returns:
        dict: Nested dict with structure [impl][shell_pair] = {'mean': float, 'std': float}
    """
    # Group data by implementation and shell pair
    data = defaultdict(lambda: defaultdict(list))

    for bench in benchmarks:
        if bench['run_type'] != 'iteration':
            continue

        impl = int(bench['impl'])
        nA = int(bench['nA'])
        nB = int(bench['nB'])
        L = int(bench['L'])
        time_ns = bench['real_time']

        # Create shell pair label
        shells = ['s', 'p', 'd', 'f', 'g', 'h', 'i']
        shell_pair = f"{shells[nA]}{shells[nB]}"

        data[impl][(L, shell_pair, nA, nB)].append(time_ns)

    # Calculate statistics
    stats = {}
    for impl in data:
        stats[impl] = {}
        for key, times in data[impl].items():
            L, shell_pair, nA, nB = key
            stats[impl][key] = {
                'L': L,
                'shell_pair': shell_pair,
                'nA': nA,
                'nB': nB,
                'mean': np.mean(times),
                'std': np.std(times, ddof=1),
                'times': times
            }

    return stats

def parse_coulomb_data(benchmarks):
    """
    Parse coulomb hermite benchmark data.

    Returns:
        dict: Nested dict with structure [impl][L_total] = {'mean': float, 'std': float, 'n_integrals': int}
    """
    data = defaultdict(lambda: defaultdict(list))
    n_integrals_map = {}

    for bench in benchmarks:
        if bench['run_type'] != 'iteration':
            continue

        impl = int(bench['impl'])
        L_total = int(bench['L_total'])
        time_ns = bench['real_time']
        n_integrals = int(bench['n_integrals'])

        data[impl][L_total].append(time_ns)
        n_integrals_map[(impl, L_total)] = n_integrals

    # Calculate statistics
    stats = {}
    for impl in data:
        stats[impl] = {}
        for L_total, times in data[impl].items():
            stats[impl][L_total] = {
                'mean': np.mean(times),
                'std': np.std(times, ddof=1),
                'n_integrals': n_integrals_map[(impl, L_total)],
                'times': times
            }

    return stats

# ============================================================================
# Plot 1: Hermite Coefficients Comparison (4-way comparison)
# ============================================================================

def plot_hermite_comparison(hermite_stats, output_dir):
    """
    Line plot comparing TMP, Layered, Symbolic, and LayeredCodegen implementations.
    """
    configure_publication_style(font_scale=1.04)

    # Create figure
    fig, ax = plt.subplots(figsize=(7.5, 4.5))

    # Define shell pairs to plot (ordered by complexity)
    shell_pairs_order = ['ss', 'sp', 'pp', 'sd', 'pd', 'dd', 'ff', 'gg']

    # Implementation labels and styles
    impl_info = {
        0: {'label': 'TMP (Original)', 'marker': 'o', 'linestyle': '-'},
        1: {'label': 'Layered (Hand-written)', 'marker': 's', 'linestyle': '-'},
        2: {'label': 'Symbolic (SymPy)', 'marker': '^', 'linestyle': '--'},
        3: {'label': 'LayeredCodegen (NEW)', 'marker': 'D', 'linestyle': '-'},
    }

    # Extract data for each implementation
    for impl_idx, impl in enumerate([0, 1, 2, 3]):
        if impl not in hermite_stats:
            continue

        means = []
        stds = []
        L_values = []

        for shell_pair in shell_pairs_order:
            # Find matching data
            found = False
            for key, stats in hermite_stats[impl].items():
                if stats['shell_pair'] == shell_pair:
                    means.append(stats['mean'])
                    stds.append(stats['std'])
                    L_values.append(stats['L'])
                    found = True
                    break

            if not found:
                means.append(np.nan)
                stds.append(0)
                L_values.append(0)

        # Create x-axis labels with L values
        x_labels = [f"{sp}\nL={L}" for sp, L in zip(shell_pairs_order, L_values)]

        # Plot line with error bars
        x = np.arange(len(shell_pairs_order))
        ax.errorbar(x, means, yerr=stds,
                   label=impl_info[impl]['label'],
                   color=OKABE_ITO_CYCLE[impl_idx],
                   marker=impl_info[impl]['marker'],
                   markersize=6,
                   linestyle=impl_info[impl]['linestyle'],
                   linewidth=1.5,
                   capsize=3,
                   capthick=0.8)

    # Log scale for y-axis
    ax.set_yscale('log')

    # Formatting
    ax.set_xlabel('Shell Pair and Total Angular Momentum L', fontsize=8)
    ax.set_ylabel('Time (nanoseconds)', fontsize=8)
    ax.set_title('Hermite Expansion Coefficients:\nTMP vs Layered vs LayeredCodegen vs Symbolic',
                 fontsize=9, pad=8)
    ax.set_xticks(np.arange(len(shell_pairs_order)))
    ax.set_xticklabels(x_labels, fontsize=7)
    ax.legend(loc='upper left', fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3, linewidth=0.5, which='both')
    ax.minorticks_on()

    # Save
    save_figure(fig, 'hermite_coefficients_comparison', output_dir)
    plt.close(fig)

# ============================================================================
# Plot 2: Hermite Coefficients vs L (Scaling with angular momentum)
# ============================================================================

def plot_hermite_vs_L(hermite_stats, output_dir):
    """
    Line plot showing execution time vs total angular momentum L.
    """
    configure_publication_style(font_scale=1.04)

    fig, ax = plt.subplots(figsize=(6, 4.5))

    impl_info = {
        0: {'label': 'TMP', 'marker': 'o'},
        1: {'label': 'Layered', 'marker': 's'},
        2: {'label': 'Symbolic', 'marker': '^'},
        3: {'label': 'LayeredCodegen', 'marker': 'D'},
    }

    for impl_idx, impl in enumerate([0, 1, 2, 3]):
        if impl not in hermite_stats:
            continue

        # Group by L
        L_data = defaultdict(list)
        for key, stats in hermite_stats[impl].items():
            L = stats['L']
            L_data[L].append(stats['mean'])

        # Calculate mean for each L
        L_values = sorted(L_data.keys())
        means = [np.mean(L_data[L]) for L in L_values]
        stds = [np.std(L_data[L]) if len(L_data[L]) > 1 else 0 for L in L_values]

        ax.errorbar(L_values, means, yerr=stds,
                   label=impl_info[impl]['label'],
                   color=OKABE_ITO_CYCLE[impl_idx],
                   marker=impl_info[impl]['marker'],
                   markersize=6,
                   linewidth=1.5,
                   capsize=3,
                   capthick=0.8)

    # Log scale for y-axis
    ax.set_yscale('log')

    # Formatting
    ax.set_xlabel('Total Angular Momentum L', fontsize=8)
    ax.set_ylabel('Time (nanoseconds)', fontsize=8)
    ax.set_title('Hermite Coefficients: Scaling with Angular Momentum', fontsize=9, pad=8)
    ax.legend(loc='upper left', fontsize=7)
    ax.grid(True, alpha=0.3, linewidth=0.5, which='both')
    ax.minorticks_on()

    save_figure(fig, 'hermite_coefficients_vs_L', output_dir)
    plt.close(fig)

# ============================================================================
# Plot 3: LayeredCodegen Speedup (Bar chart)
# ============================================================================

def plot_layered_codegen_speedup(hermite_stats, output_dir):
    """
    Bar chart showing LayeredCodegen speedup for ss shell.
    """
    configure_publication_style(font_scale=1.04)

    fig, ax = plt.subplots(figsize=(5, 4))

    # Get ss shell data for each implementation
    ss_times = {}
    for impl in [0, 1, 2, 3]:
        if impl not in hermite_stats:
            continue
        for key, stats in hermite_stats[impl].items():
            if stats['shell_pair'] == 'ss':
                ss_times[impl] = stats['mean']
                break

    # Calculate speedups (LayeredCodegen vs others)
    if 3 in ss_times:
        layered_codegen_time = ss_times[3]

        comparisons = []
        speedups = []
        colors = []

        if 1 in ss_times:
            comparisons.append('vs Layered')
            speedups.append(ss_times[1] / layered_codegen_time)
            colors.append(OKABE_ITO_CYCLE[1])

        if 0 in ss_times:
            comparisons.append('vs TMP')
            speedups.append(ss_times[0] / layered_codegen_time)
            colors.append(OKABE_ITO_CYCLE[0])

        if 2 in ss_times:
            comparisons.append('vs Symbolic')
            speedups.append(ss_times[2] / layered_codegen_time)
            colors.append(OKABE_ITO_CYCLE[2])

        # Plot bars
        x = np.arange(len(comparisons))
        bars = ax.bar(x, speedups, color=colors, width=0.6, edgecolor='black', linewidth=0.8)

        # Add value labels on bars
        for bar, speedup in zip(bars, speedups):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{speedup:.2f}×',
                   ha='center', va='bottom', fontsize=7)

        # Formatting
        ax.set_ylabel('Speedup Factor', fontsize=8)
        ax.set_xlabel('Comparison', fontsize=8)
        ax.set_title('LayeredCodegen Speedup (ss shell)', fontsize=9, pad=8)
        ax.set_xticks(x)
        ax.set_xticklabels(comparisons, fontsize=7)
        ax.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax.grid(axis='y', alpha=0.3, linewidth=0.5)

    save_figure(fig, 'hermite_layered_codegen_speedup', output_dir)
    plt.close(fig)

# ============================================================================
# Plot 4: Coulomb Hermite Comparison (Line plot)
# ============================================================================

def plot_coulomb_comparison(coulomb_stats, output_dir):
    """
    Line plot comparing TMP, Layered, and LayeredCodegen for Coulomb Hermite integrals.
    """
    configure_publication_style(font_scale=1.04)

    fig, ax = plt.subplots(figsize=(6, 4.5))

    impl_info = {
        0: {'label': 'TMP', 'marker': 'o'},
        1: {'label': 'Layered', 'marker': 's'},
        3: {'label': 'LayeredCodegen', 'marker': 'D'},
    }

    # Get L_total values
    all_L = set()
    for impl in [0, 1, 3]:
        if impl in coulomb_stats:
            all_L.update(coulomb_stats[impl].keys())
    L_values = sorted(all_L)

    for impl_idx, impl in enumerate([0, 1, 3]):
        if impl not in coulomb_stats:
            continue

        means = [coulomb_stats[impl].get(L, {'mean': np.nan})['mean'] for L in L_values]
        stds = [coulomb_stats[impl].get(L, {'std': 0})['std'] for L in L_values]

        ax.errorbar(L_values, means, yerr=stds,
                   label=impl_info[impl]['label'],
                   color=OKABE_ITO_CYCLE[impl_idx],
                   marker=impl_info[impl]['marker'],
                   markersize=6,
                   linewidth=1.5,
                   capsize=3,
                   capthick=0.8)

    # Log scale for y-axis
    ax.set_yscale('log')

    # Formatting
    ax.set_xlabel('Total Angular Momentum $L_{total}$', fontsize=8)
    ax.set_ylabel('Time (nanoseconds)', fontsize=8)
    ax.set_title('Coulomb Hermite Integrals: TMP vs Layered vs LayeredCodegen', fontsize=9, pad=8)
    ax.legend(loc='upper left', fontsize=7)
    ax.grid(True, alpha=0.3, linewidth=0.5, which='both')
    ax.minorticks_on()

    save_figure(fig, 'coulomb_hermite_comparison', output_dir)
    plt.close(fig)

# ============================================================================
# Plot 5: Coulomb Hermite Scaling (Log-log plot)
# ============================================================================

def plot_coulomb_scaling(coulomb_stats, output_dir):
    """
    Log-log plot showing scaling behavior with number of integrals.
    """
    configure_publication_style(font_scale=1.04)

    fig, ax = plt.subplots(figsize=(6, 4.5))

    impl_info = {
        0: {'label': 'TMP', 'marker': 'o'},
        1: {'label': 'Layered', 'marker': 's'},
        3: {'label': 'LayeredCodegen', 'marker': 'D'},
    }

    for impl_idx, impl in enumerate([0, 1, 3]):
        if impl not in coulomb_stats:
            continue

        L_values = sorted(coulomb_stats[impl].keys())

        n_integrals_list = []
        means = []
        stds = []

        for L in L_values:
            stats = coulomb_stats[impl][L]
            n_integrals_list.append(stats['n_integrals'])
            means.append(stats['mean'])
            stds.append(stats['std'])

        ax.errorbar(n_integrals_list, means, yerr=stds,
                   label=impl_info[impl]['label'],
                   color=OKABE_ITO_CYCLE[impl_idx],
                   marker=impl_info[impl]['marker'],
                   markersize=6,
                   linewidth=1.5,
                   capsize=3,
                   capthick=0.8)

    # Log-log scale
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Formatting
    ax.set_xlabel('Number of Integrals', fontsize=8)
    ax.set_ylabel('Time (nanoseconds)', fontsize=8)
    ax.set_title('Coulomb Hermite: Problem Size Scaling', fontsize=9, pad=8)
    ax.legend(loc='upper left', fontsize=7)
    ax.grid(True, alpha=0.3, linewidth=0.5, which='both')
    ax.minorticks_on()

    save_figure(fig, 'coulomb_hermite_scaling', output_dir)
    plt.close(fig)

# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Generate all 5 plots with publication-quality formatting."""
    # Paths
    data_dir = Path('/home/ruben/Research/Science/Projects/RECURSUM/benchmarks/results/raw')
    output_dir = Path('/home/ruben/Research/Science/Projects/RECURSUM/benchmarks/results/figures')

    hermite_file = data_dir / 'hermite_coefficients.json'
    coulomb_file = data_dir / 'coulomb_hermite.json'

    print("=" * 70)
    print("Regenerating All 5 Benchmark Plots with Larger Font Sizes")
    print("=" * 70)

    # Load and parse data
    print("\nLoading data...")
    hermite_benchmarks = load_benchmark_data(hermite_file)
    coulomb_benchmarks = load_benchmark_data(coulomb_file)

    print(f"  Hermite coefficients: {len(hermite_benchmarks)} benchmark entries")
    print(f"  Coulomb Hermite: {len(coulomb_benchmarks)} benchmark entries")

    print("\nParsing benchmark results...")
    hermite_stats = parse_hermite_data(hermite_benchmarks)
    coulomb_stats = parse_coulomb_data(coulomb_benchmarks)

    print(f"  Hermite implementations found: {sorted(hermite_stats.keys())}")
    print(f"  Coulomb implementations found: {sorted(coulomb_stats.keys())}")

    # Generate plots
    print("\nGenerating plots with smaller fonts...")

    print("\n[1/5] Hermite Coefficients Comparison (4-way)")
    plot_hermite_comparison(hermite_stats, output_dir)

    print("\n[2/5] Hermite Coefficients vs L")
    plot_hermite_vs_L(hermite_stats, output_dir)

    print("\n[3/5] LayeredCodegen Speedup")
    plot_layered_codegen_speedup(hermite_stats, output_dir)

    print("\n[4/5] Coulomb Hermite Comparison")
    plot_coulomb_comparison(coulomb_stats, output_dir)

    print("\n[5/5] Coulomb Hermite Scaling")
    plot_coulomb_scaling(coulomb_stats, output_dir)

    print("\n" + "=" * 70)
    print("All 5 plots regenerated successfully with larger font sizes (30% increase)!")
    print(f"Output directory: {output_dir}")
    print("=" * 70)

if __name__ == '__main__':
    main()
