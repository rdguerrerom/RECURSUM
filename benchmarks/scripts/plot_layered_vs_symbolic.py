#!/usr/bin/env python3
"""
Generate publication-ready plots comparing Layered CSE (impl=1) vs Symbolic (impl=2) implementations.
Excludes TMP (impl=0) from all visualizations.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# ============================================================================
# Configuration from Scientific Plotting Guide
# ============================================================================

# Okabe-Ito colorblind-safe palette
OKABE_ITO = {
    'black':     '#000000',
    'orange':    '#E69F00',
    'sky_blue':  '#56B4E9',
    'green':     '#009E73',
    'yellow':    '#F0E442',
    'blue':      '#0072B2',
    'vermillion':'#D55E00',
    'purple':    '#CC79A7',
}

# Use orange for Layered, green for Symbolic (colorblind-safe)
LAYERED_COLOR = OKABE_ITO['orange']
SYMBOLIC_COLOR = OKABE_ITO['green']

# Journal column widths
COLUMN_WIDTHS = {
    'nature': {'single': 3.46, 'onehalf': 4.76, 'double': 7.09},
    'science': {'single': 3.35, 'onehalf': 4.57, 'double': 7.01},
    'aps': {'single': 3.39, 'onehalf': 5.12, 'double': 7.01},
}

def configure_publication_style(font_scale=1.2):
    """Configure matplotlib for publication-quality figures."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',

        # Font sizes
        'axes.labelsize': 10 * font_scale,
        'axes.titlesize': 11 * font_scale,
        'xtick.labelsize': 9 * font_scale,
        'ytick.labelsize': 9 * font_scale,
        'legend.fontsize': 8 * font_scale,
        'font.size': 9 * font_scale,

        # Axes configuration
        'axes.linewidth': 0.8,
        'axes.labelpad': 4,
        'axes.spines.top': False,
        'axes.spines.right': False,

        # Tick configuration
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'xtick.direction': 'out',
        'ytick.direction': 'out',

        # Legend configuration
        'legend.frameon': False,
        'legend.borderpad': 0.4,
        'legend.handlelength': 1.5,

        # Figure configuration
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
    })

# ============================================================================
# Data Loading and Processing
# ============================================================================

def load_benchmark_data(json_path):
    """Load Google Benchmark JSON data."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data['benchmarks']

def extract_hermite_data(benchmarks):
    """
    Extract Hermite coefficient data for Layered (impl=1) and Symbolic (impl=2) ONLY.
    Returns dict with shell pairs as keys.
    """
    # Shell pair mapping
    shell_pairs = {
        0: 'ss', 1: 'sp', 2: 'pp',
        3: 'sd', 4: 'pd', 5: 'dd',
        6: 'ff', 7: 'gg'
    }

    data = {impl: {} for impl in [1, 2]}  # Only Layered and Symbolic

    for bench in benchmarks:
        impl_raw = bench.get('impl', -1)
        # Skip if NaN or not a valid number
        if np.isnan(impl_raw) if isinstance(impl_raw, float) else False:
            continue
        impl = int(impl_raw)
        if impl not in [1, 2]:  # Exclude TMP (impl=0)
            continue

        L = int(bench['L'])
        shell = shell_pairs.get(L, f'L{L}')

        if shell not in data[impl]:
            data[impl][shell] = {'times': [], 'L': L}

        data[impl][shell]['times'].append(bench['cpu_time'])

    # Calculate statistics
    for impl in [1, 2]:
        for shell in data[impl]:
            times = np.array(data[impl][shell]['times'])
            data[impl][shell]['mean'] = np.mean(times)
            data[impl][shell]['std'] = np.std(times)
            data[impl][shell]['n'] = len(times)

    return data

def extract_coulomb_data(benchmarks):
    """
    Extract Coulomb Hermite data for Layered implementation ONLY.
    Returns dict with L_total as keys.
    """
    data = {impl: {} for impl in [1]}  # Only Layered (TMP data incomplete)

    for bench in benchmarks:
        impl_raw = bench.get('impl', -1)
        # Skip if NaN or not a valid number
        if np.isnan(impl_raw) if isinstance(impl_raw, float) else False:
            continue
        impl = int(impl_raw)
        if impl != 1:  # Only Layered
            continue

        L = int(bench['L_total'])

        if L not in data[impl]:
            data[impl][L] = {'times': []}

        data[impl][L]['times'].append(bench['cpu_time'])

    # Calculate statistics
    for impl in [1]:
        for L in data[impl]:
            times = np.array(data[impl][L]['times'])
            data[impl][L]['mean'] = np.mean(times)
            data[impl][L]['std'] = np.std(times)
            data[impl][L]['n'] = len(times)

    return data

# ============================================================================
# Plot 1: Hermite Coefficients Bar Chart Comparison
# ============================================================================

def plot_hermite_bar_comparison(data, output_dir):
    """Bar chart comparing Layered CSE vs Symbolic for Hermite coefficients."""
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=(7.09, 4.5))

    # Extract shell pairs in order
    shell_order = ['ss', 'sp', 'pp', 'sd', 'pd', 'dd', 'ff', 'gg']

    # Prepare data arrays
    layered_means = []
    layered_stds = []
    symbolic_means = []
    symbolic_stds = []
    x_labels = []

    for shell in shell_order:
        if shell in data[1] and shell in data[2]:
            layered_means.append(data[1][shell]['mean'])
            layered_stds.append(data[1][shell]['std'])
            symbolic_means.append(data[2][shell]['mean'])
            symbolic_stds.append(data[2][shell]['std'])
            x_labels.append(shell)

    x = np.arange(len(x_labels))
    width = 0.35

    # Create bars
    bars1 = ax.bar(x - width/2, layered_means, width, yerr=layered_stds,
                   label='Layered CSE', color=LAYERED_COLOR, alpha=0.8,
                   capsize=3, error_kw={'linewidth': 1})
    bars2 = ax.bar(x + width/2, symbolic_means, width, yerr=symbolic_stds,
                   label='Symbolic', color=SYMBOLIC_COLOR, alpha=0.8,
                   capsize=3, error_kw={'linewidth': 1})

    # Formatting
    ax.set_xlabel('Shell Pair')
    ax.set_ylabel('Execution Time (ns)')
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels)
    ax.legend(loc='upper left')
    ax.grid(axis='y', alpha=0.3, linewidth=0.5)

    # Save
    fig.tight_layout()
    fig.savefig(output_dir / 'hermite_layered_vs_symbolic.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / 'hermite_layered_vs_symbolic.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_dir / 'hermite_layered_vs_symbolic.pdf'}")
    print(f"Saved: {output_dir / 'hermite_layered_vs_symbolic.png'}")

# ============================================================================
# Plot 2: Hermite Scaling Analysis (Log Scale)
# ============================================================================

def plot_hermite_scaling(data, output_dir):
    """Line plot showing scaling with angular momentum."""
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=(6.0, 4.5))

    # Extract L values and times
    L_layered = []
    mean_layered = []
    std_layered = []
    L_symbolic = []
    mean_symbolic = []
    std_symbolic = []

    for shell, shell_data in sorted(data[1].items(), key=lambda x: x[1]['L']):
        L_layered.append(shell_data['L'])
        mean_layered.append(shell_data['mean'])
        std_layered.append(shell_data['std'])

    for shell, shell_data in sorted(data[2].items(), key=lambda x: x[1]['L']):
        L_symbolic.append(shell_data['L'])
        mean_symbolic.append(shell_data['mean'])
        std_symbolic.append(shell_data['std'])

    # Plot with markers
    ax.errorbar(L_layered, mean_layered, yerr=std_layered,
                marker='o', markersize=6, linewidth=2, capsize=4,
                label='Layered CSE', color=LAYERED_COLOR)
    ax.errorbar(L_symbolic, mean_symbolic, yerr=std_symbolic,
                marker='s', markersize=6, linewidth=2, capsize=4,
                label='Symbolic', color=SYMBOLIC_COLOR)

    # Find and annotate crossover point
    L_common = sorted(set(L_layered) & set(L_symbolic))
    layered_interp = np.interp(L_common, L_layered, mean_layered)
    symbolic_interp = np.interp(L_common, L_symbolic, mean_symbolic)

    # Formatting
    ax.set_xlabel('Total Angular Momentum (L)')
    ax.set_ylabel('Execution Time (ns)')
    ax.set_yscale('log')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3, linewidth=0.5, which='both')
    ax.set_xticks(range(0, 9))

    # Save
    fig.tight_layout()
    fig.savefig(output_dir / 'hermite_layered_vs_symbolic_scaling.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / 'hermite_layered_vs_symbolic_scaling.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_dir / 'hermite_layered_vs_symbolic_scaling.pdf'}")
    print(f"Saved: {output_dir / 'hermite_layered_vs_symbolic_scaling.png'}")

# ============================================================================
# Plot 3: Coulomb Hermite Layered Only
# ============================================================================

def plot_coulomb_layered(data, output_dir):
    """Single-line plot showing Layered CSE performance for Coulomb Hermite."""
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=(6.0, 4.5))

    # Extract L_total values and times for Layered only
    L_values = []
    means = []
    stds = []

    for L in sorted(data[1].keys()):
        L_values.append(L)
        means.append(data[1][L]['mean'])
        stds.append(data[1][L]['std'])

    # Plot
    ax.errorbar(L_values, means, yerr=stds,
                marker='o', markersize=7, linewidth=2.5, capsize=4,
                label='Layered CSE', color=LAYERED_COLOR)

    # Formatting
    ax.set_xlabel('Total Angular Momentum ($L_{total}$)')
    ax.set_ylabel('Execution Time (ns)')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_xticks(range(0, max(L_values)+1))

    # Save
    fig.tight_layout()
    fig.savefig(output_dir / 'coulomb_layered_scaling.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / 'coulomb_layered_scaling.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_dir / 'coulomb_layered_scaling.pdf'}")
    print(f"Saved: {output_dir / 'coulomb_layered_scaling.png'}")

# ============================================================================
# Plot 4: Relative Performance (Speedup Ratio)
# ============================================================================

def plot_relative_performance(data, output_dir):
    """Bar chart showing Symbolic/Layered speedup ratio."""
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=(7.09, 4.5))

    # Extract shell pairs in order
    shell_order = ['ss', 'sp', 'pp', 'sd', 'pd', 'dd', 'ff', 'gg']

    speedup_ratios = []
    x_labels = []

    for shell in shell_order:
        if shell in data[1] and shell in data[2]:
            ratio = data[2][shell]['mean'] / data[1][shell]['mean']
            speedup_ratios.append(ratio)
            x_labels.append(shell)

    x = np.arange(len(x_labels))

    # Color bars based on speedup (> 1 means Layered is faster)
    colors = [LAYERED_COLOR if r > 1 else SYMBOLIC_COLOR for r in speedup_ratios]

    # Create bars
    bars = ax.bar(x, speedup_ratios, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)

    # Add horizontal line at y=1.0
    ax.axhline(y=1.0, color='black', linestyle='--', linewidth=1.5, label='Equal Performance')

    # Formatting
    ax.set_xlabel('Shell Pair')
    ax.set_ylabel('Speedup Factor (Symbolic / Layered)')
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels)
    ax.grid(axis='y', alpha=0.3, linewidth=0.5)
    ax.legend(loc='upper right')

    # Add text labels on bars
    for i, (bar, ratio) in enumerate(zip(bars, speedup_ratios)):
        height = bar.get_height()
        label = f'{ratio:.2f}×'
        ax.text(bar.get_x() + bar.get_width()/2., height,
                label, ha='center', va='bottom', fontsize=8)

    # Save
    fig.tight_layout()
    fig.savefig(output_dir / 'hermite_layered_symbolic_speedup.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / 'hermite_layered_symbolic_speedup.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_dir / 'hermite_layered_symbolic_speedup.pdf'}")
    print(f"Saved: {output_dir / 'hermite_layered_symbolic_speedup.png'}")

# ============================================================================
# Main Execution
# ============================================================================

def main():
    # Paths
    base_dir = Path('/home/ruben/Research/Science/Projects/RECURSUM/benchmarks')
    hermite_json = base_dir / 'results/raw/hermite_coefficients.json'
    coulomb_json = base_dir / 'results/raw/coulomb_hermite.json'
    output_dir = base_dir / 'results/figures'

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading benchmark data...")
    hermite_benchmarks = load_benchmark_data(hermite_json)
    coulomb_benchmarks = load_benchmark_data(coulomb_json)

    print("Extracting Hermite data (Layered and Symbolic only)...")
    hermite_data = extract_hermite_data(hermite_benchmarks)

    print("Extracting Coulomb data (Layered only)...")
    coulomb_data = extract_coulomb_data(coulomb_benchmarks)

    print("\nGenerating plots...")
    print("-" * 60)

    # Generate all 4 plots
    plot_hermite_bar_comparison(hermite_data, output_dir)
    plot_hermite_scaling(hermite_data, output_dir)
    plot_coulomb_layered(coulomb_data, output_dir)
    plot_relative_performance(hermite_data, output_dir)

    print("-" * 60)
    print(f"\nAll plots saved to: {output_dir}")
    print("\nGenerated files:")
    print("  1. hermite_layered_vs_symbolic.pdf/png")
    print("  2. hermite_layered_vs_symbolic_scaling.pdf/png")
    print("  3. coulomb_layered_scaling.pdf/png")
    print("  4. hermite_layered_symbolic_speedup.pdf/png")
    print("\nTotal: 8 files (4 PDF + 4 PNG)")

if __name__ == '__main__':
    main()
