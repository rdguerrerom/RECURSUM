#!/usr/bin/env python3
"""
Generate publication-ready plots for McMurchie-Davidson benchmark results.

This script follows the scientific plotting guidelines at:
~/Research/Writing/scientific_plotting_agent_guide.md
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from collections import defaultdict

# ============================================================================
# Publication Style Configuration (from plotting guide)
# ============================================================================

# Okabe-Ito colorblind-safe palette
OKABE_ITO_CYCLE = [
    '#0072B2',  # blue
    '#E69F00',  # orange
    '#009E73',  # green
    '#CC79A7',  # purple
    '#D55E00',  # vermillion
    '#56B4E9',  # sky blue
    '#F0E442',  # yellow
    '#000000',  # black
]

# Journal column widths (Nature)
COLUMN_WIDTH_SINGLE = 3.46  # inches
COLUMN_WIDTH_DOUBLE = 7.09  # inches

def configure_publication_style(font_scale=1.2):
    """Configure matplotlib for publication-quality figures."""
    plt.rcParams.update({
        # Font configuration
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',

        # Font sizes (scaled for potential LaTeX reduction)
        'font.size': 9 * font_scale,
        'axes.labelsize': 10 * font_scale,
        'axes.titlesize': 11 * font_scale,
        'xtick.labelsize': 9 * font_scale,
        'ytick.labelsize': 9 * font_scale,
        'legend.fontsize': 8 * font_scale,

        # Axes configuration
        'axes.linewidth': 0.8,
        'axes.labelpad': 4,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.prop_cycle': plt.cycler(color=OKABE_ITO_CYCLE),

        # Tick configuration
        'xtick.major.width': 0.8,
        'xtick.minor.width': 0.6,
        'ytick.major.width': 0.8,
        'ytick.minor.width': 0.6,
        'xtick.major.size': 4,
        'xtick.minor.size': 2,
        'ytick.major.size': 4,
        'ytick.minor.size': 2,
        'xtick.direction': 'out',
        'ytick.direction': 'out',

        # Legend configuration
        'legend.frameon': False,
        'legend.borderpad': 0.4,
        'legend.handlelength': 1.5,
        'legend.handletextpad': 0.5,
        'legend.columnspacing': 1.0,

        # Line configuration
        'lines.linewidth': 1.5,
        'lines.markersize': 6,

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
        print(f"Saved: {filepath}")

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
# Plot 1: Hermite Coefficients Comparison (Bar Chart)
# ============================================================================

def plot_hermite_comparison(hermite_stats, output_dir):
    """
    Bar chart comparing TMP, Layered CSE, and Symbolic implementations.
    """
    configure_publication_style(font_scale=1.2)

    # Create figure
    fig, ax = plt.subplots(figsize=(COLUMN_WIDTH_DOUBLE, COLUMN_WIDTH_DOUBLE * 0.4))

    # Define shell pairs to plot (ordered by complexity)
    shell_pairs_order = ['ss', 'sp', 'pp', 'sd', 'pd', 'dd', 'ff', 'gg']

    # Implementation labels
    impl_labels = {
        0: 'TMP',
        1: 'Layered CSE',
        2: 'Symbolic'
    }

    # Extract data for each implementation
    plot_data = {impl: [] for impl in [0, 1, 2]}

    for shell_pair in shell_pairs_order:
        for impl in [0, 1, 2]:
            # Find matching data
            found = False
            for key, stats in hermite_stats[impl].items():
                if stats['shell_pair'] == shell_pair:
                    plot_data[impl].append({
                        'mean': stats['mean'],
                        'std': stats['std'],
                        'shell_pair': shell_pair
                    })
                    found = True
                    break
            if not found:
                plot_data[impl].append({'mean': np.nan, 'std': np.nan, 'shell_pair': shell_pair})

    # Plot bars
    x = np.arange(len(shell_pairs_order))
    width = 0.25

    for i, impl in enumerate([0, 1, 2]):
        means = [d['mean'] for d in plot_data[impl]]
        stds = [d['std'] for d in plot_data[impl]]

        offset = (i - 1) * width
        ax.bar(x + offset, means, width, yerr=stds,
               label=impl_labels[impl], color=OKABE_ITO_CYCLE[i],
               capsize=3, error_kw={'linewidth': 0.8})

    # Log scale for y-axis
    ax.set_yscale('log')

    # Formatting
    ax.set_xlabel('Shell Pair')
    ax.set_ylabel('Execution Time (ns)')
    ax.set_xticks(x)
    ax.set_xticklabels(shell_pairs_order)
    ax.legend(loc='upper left', ncol=3)
    ax.grid(axis='y', alpha=0.3, linewidth=0.5)

    # Save
    save_figure(fig, 'hermite_coefficients_comparison', output_dir)
    plt.close(fig)

# ============================================================================
# Plot 2: Hermite Coefficients vs L (Line Plot)
# ============================================================================

def plot_hermite_vs_L(hermite_stats, output_dir):
    """
    Line plot showing execution time vs total angular momentum L.
    """
    configure_publication_style(font_scale=1.2)

    fig, ax = plt.subplots(figsize=(COLUMN_WIDTH_SINGLE, COLUMN_WIDTH_SINGLE * 0.75))

    impl_labels = {0: 'TMP', 1: 'Layered CSE', 2: 'Symbolic'}
    markers = ['o', 's', '^']

    for impl_idx, impl in enumerate([0, 1, 2]):
        # Group by L
        L_data = defaultdict(list)
        for key, stats in hermite_stats[impl].items():
            L = stats['L']
            L_data[L].append(stats['mean'])

        # Calculate mean for each L
        L_values = sorted(L_data.keys())
        means = [np.mean(L_data[L]) for L in L_values]
        stds = [np.std(L_data[L]) for L in L_values]

        ax.errorbar(L_values, means, yerr=stds,
                   label=impl_labels[impl],
                   color=OKABE_ITO_CYCLE[impl_idx],
                   marker=markers[impl_idx],
                   markersize=6,
                   linewidth=1.5,
                   capsize=3,
                   capthick=0.8)

    # Log scale for y-axis
    ax.set_yscale('log')

    # Formatting
    ax.set_xlabel('Total Angular Momentum $L$')
    ax.set_ylabel('Execution Time (ns)')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3, linewidth=0.5, which='both')
    ax.minorticks_on()

    save_figure(fig, 'hermite_coefficients_vs_L', output_dir)
    plt.close(fig)

# ============================================================================
# Plot 3: Coulomb Hermite Comparison (Bar Chart)
# ============================================================================

def plot_coulomb_comparison(coulomb_stats, output_dir):
    """
    Bar chart comparing TMP and Layered CSE for Coulomb Hermite integrals.
    """
    configure_publication_style(font_scale=1.2)

    fig, ax = plt.subplots(figsize=(COLUMN_WIDTH_DOUBLE, COLUMN_WIDTH_DOUBLE * 0.4))

    impl_labels = {0: 'TMP', 1: 'Layered CSE'}

    # Get L_total values
    L_values = sorted(set(coulomb_stats[0].keys()) | set(coulomb_stats[1].keys()))

    # Extract data
    x = np.arange(len(L_values))
    width = 0.35

    for impl_idx, impl in enumerate([0, 1]):
        means = [coulomb_stats[impl].get(L, {'mean': np.nan})['mean'] for L in L_values]
        stds = [coulomb_stats[impl].get(L, {'std': 0})['std'] for L in L_values]

        offset = (impl_idx - 0.5) * width
        ax.bar(x + offset, means, width, yerr=stds,
               label=impl_labels[impl], color=OKABE_ITO_CYCLE[impl_idx],
               capsize=3, error_kw={'linewidth': 0.8})

    # Log scale for y-axis
    ax.set_yscale('log')

    # Formatting
    ax.set_xlabel('Total Angular Momentum $L_{\\mathrm{total}}$')
    ax.set_ylabel('Execution Time (ns)')
    ax.set_xticks(x)
    ax.set_xticklabels([str(L) for L in L_values])
    ax.legend(loc='upper left')
    ax.grid(axis='y', alpha=0.3, linewidth=0.5)

    save_figure(fig, 'coulomb_hermite_comparison', output_dir)
    plt.close(fig)

# ============================================================================
# Plot 4: Coulomb Hermite Scaling (Line Plot)
# ============================================================================

def plot_coulomb_scaling(coulomb_stats, output_dir):
    """
    Line plot showing cost per integral vs L_total.
    """
    configure_publication_style(font_scale=1.2)

    fig, ax = plt.subplots(figsize=(COLUMN_WIDTH_SINGLE, COLUMN_WIDTH_SINGLE * 0.75))

    impl_labels = {0: 'TMP', 1: 'Layered CSE'}
    markers = ['o', 's']

    for impl_idx, impl in enumerate([0, 1]):
        L_values = sorted(coulomb_stats[impl].keys())

        # Calculate cost per integral
        cost_per_integral = []
        cost_per_integral_std = []

        for L in L_values:
            stats = coulomb_stats[impl][L]
            n_integrals = stats['n_integrals']
            mean_time = stats['mean']
            std_time = stats['std']

            cost_per_integral.append(mean_time / n_integrals)
            cost_per_integral_std.append(std_time / n_integrals)

        ax.errorbar(L_values, cost_per_integral, yerr=cost_per_integral_std,
                   label=impl_labels[impl],
                   color=OKABE_ITO_CYCLE[impl_idx],
                   marker=markers[impl_idx],
                   markersize=6,
                   linewidth=1.5,
                   capsize=3,
                   capthick=0.8)

    # Log scale for y-axis
    ax.set_yscale('log')

    # Formatting
    ax.set_xlabel('Total Angular Momentum $L_{\\mathrm{total}}$')
    ax.set_ylabel('Time/Integral (ns/integral)')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3, linewidth=0.5, which='both')
    ax.minorticks_on()

    save_figure(fig, 'coulomb_hermite_scaling', output_dir)
    plt.close(fig)

# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Generate all plots."""
    # Paths
    data_dir = Path('/home/ruben/Research/Science/Projects/RECURSUM/benchmarks/results/raw')
    output_dir = Path('/home/ruben/Research/Science/Projects/RECURSUM/benchmarks/results/figures')

    hermite_file = data_dir / 'hermite_coefficients.json'
    coulomb_file = data_dir / 'coulomb_hermite.json'

    print("=" * 70)
    print("Generating Publication-Ready McMurchie-Davidson Benchmark Plots")
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

    print(f"  Hermite implementations: {list(hermite_stats.keys())}")
    print(f"  Coulomb implementations: {list(coulomb_stats.keys())}")

    # Generate plots
    print("\nGenerating plots...")
    print("\n[1/4] Hermite Coefficients Comparison (Bar Chart)")
    plot_hermite_comparison(hermite_stats, output_dir)

    print("\n[2/4] Hermite Coefficients vs L (Line Plot)")
    plot_hermite_vs_L(hermite_stats, output_dir)

    print("\n[3/4] Coulomb Hermite Comparison (Bar Chart)")
    plot_coulomb_comparison(coulomb_stats, output_dir)

    print("\n[4/4] Coulomb Hermite Scaling (Line Plot)")
    plot_coulomb_scaling(coulomb_stats, output_dir)

    print("\n" + "=" * 70)
    print("All plots generated successfully!")
    print(f"Output directory: {output_dir}")
    print("=" * 70)

if __name__ == '__main__':
    main()
