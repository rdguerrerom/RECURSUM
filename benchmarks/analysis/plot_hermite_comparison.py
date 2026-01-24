#!/usr/bin/env python3
"""
Publication-Quality Plots for TMP vs Symbolic McMurchie-Davidson Recurrence Comparison

This script generates publication-ready figures from ACTUAL benchmark data.
All values are from real measurements - no theoretical or estimated values.

Following guidelines from: ~/Research/Writing/scientific_plotting_agent_guide.md
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from collections import defaultdict
import statistics


# =============================================================================
# PUBLICATION STYLE CONFIGURATION (from scientific_plotting_agent_guide.md)
# =============================================================================

# Journal column widths in inches
COLUMN_WIDTHS = {
    'nature': {'single': 3.46, 'onehalf': 4.76, 'double': 7.09},
    'science': {'single': 3.35, 'onehalf': 4.57, 'double': 7.01},
    'aps': {'single': 3.39, 'onehalf': 5.12, 'double': 7.01},
    'acs': {'single': 3.33, 'onehalf': 4.33, 'double': 7.00},
    'elsevier': {'single': 3.54, 'onehalf': 5.51, 'double': 7.48},
}

# Okabe-Ito colorblind-safe palette (MANDATORY per guide)
OKABE_ITO = {
    'blue':      '#0072B2',
    'orange':    '#E69F00',
    'green':     '#009E73',
    'purple':    '#CC79A7',
    'vermillion':'#D55E00',
    'sky_blue':  '#56B4E9',
    'yellow':    '#F0E442',
    'black':     '#000000',
}

# Colors for our implementations
IMPL_COLORS = {
    'TMP': OKABE_ITO['blue'],
    'Symbolic': OKABE_ITO['orange'],
}

# Markers for distinguishing features (supplement color per guide)
MARKERS = {
    'TMP': 'o',
    'Symbolic': 's',
}


def get_figsize(journal='aps', width='single', aspect=0.75):
    """Get figure dimensions for target journal."""
    w = COLUMN_WIDTHS[journal][width]
    h = w * aspect
    return (w, h)


def get_font_scale(fig_width_inches, target_min_font=6):
    """Calculate font scale factor based on figure width."""
    reference_width = 7.0  # double column
    scale = reference_width / fig_width_inches
    return max(1.0, min(scale, 2.5))


def configure_publication_style(journal='aps', font_scale=1.0, use_latex=False):
    """Configure matplotlib for publication-quality figures."""

    # Base font sizes (will be scaled)
    base_sizes = {
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 8,
        'figure.titlesize': 12,
        'font.size': 9,
    }

    # Apply scaling
    scaled_sizes = {k: v * font_scale for k, v in base_sizes.items()}

    plt.rcParams.update({
        # Font configuration - using serif for scientific publications
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'DejaVu Serif', 'Computer Modern Roman'],
        'mathtext.fontset': 'cm',  # Computer Modern for math

        # Font sizes
        **scaled_sizes,

        # Axes configuration
        'axes.linewidth': 0.8,
        'axes.labelpad': 4,
        'axes.spines.top': False,
        'axes.spines.right': False,

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

        # LaTeX rendering
        'text.usetex': use_latex,
    })


def add_panel_label(ax, label, fontsize=12, fontweight='bold', offset=(-0.12, 1.05)):
    """Add panel label (a), (b), etc. to axes."""
    ax.text(offset[0], offset[1], label,
            transform=ax.transAxes,
            fontsize=fontsize,
            fontweight=fontweight,
            va='top',
            ha='left')


# =============================================================================
# DATA LOADING
# =============================================================================

def load_benchmark_data(json_path):
    """
    Load benchmark results from JSON file.

    Returns only ACTUAL measured data, filtering out aggregate statistics.
    """
    with open(json_path, 'r') as f:
        data = json.load(f)

    context = data.get('context', {})
    benchmarks = []

    for bench in data.get('benchmarks', []):
        # Skip aggregate results (mean, median, stddev, cv)
        if bench.get('aggregate_name'):
            continue

        # Skip entries with NaN values
        impl = bench.get('impl', -1)
        if impl != impl:  # NaN check
            continue

        benchmarks.append({
            'name': bench.get('name', ''),
            'real_time': bench.get('real_time', 0),
            'cpu_time': bench.get('cpu_time', 0),
            'iterations': bench.get('iterations', 0),
            'impl': int(impl),
            'L_total': int(bench.get('L_total', 0)) if bench.get('L_total', 0) == bench.get('L_total', 0) else 0,
            'nA': int(bench.get('nA', 0)) if 'nA' in bench else 0,
            'nB': int(bench.get('nB', 0)) if 'nB' in bench else 0,
            't': int(bench.get('t', 0)) if 't' in bench else 0,
        })

    return context, benchmarks


def extract_coefficient_from_name(name):
    """Extract nA, nB, t from benchmark name like 'BM_Compare_E_0_0_0'."""
    parts = name.split('_')
    try:
        # Look for pattern with three consecutive integers at the end
        for i in range(len(parts) - 2):
            try:
                nA = int(parts[i])
                nB = int(parts[i + 1])
                t = int(parts[i + 2])
                return nA, nB, t
            except ValueError:
                continue
    except (ValueError, IndexError):
        pass
    return None, None, None


# =============================================================================
# PLOT GENERATION
# =============================================================================

def plot_performance_comparison(benchmarks, output_dir, journal='aps'):
    """
    Plot 1: TMP vs Symbolic Performance Comparison Bar Chart

    Compares execution times for each E coefficient using ACTUAL benchmark data.
    """
    fig_width = COLUMN_WIDTHS[journal]['onehalf']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal=journal, font_scale=font_scale)

    # Group data by coefficient and implementation
    coeff_data = defaultdict(lambda: {'TMP': [], 'Symbolic': []})

    for bench in benchmarks:
        if 'Compare' in bench['name'] or 'E_' in bench['name']:
            nA, nB, t = extract_coefficient_from_name(bench['name'])
            if nA is not None:
                impl_name = 'TMP' if bench['impl'] == 0 else 'Symbolic' if bench['impl'] == 1 else None
                if impl_name:
                    coeff_data[(nA, nB, t)][impl_name].append(bench['real_time'])

    if not coeff_data:
        print("Warning: No comparison data found for performance plot")
        return None

    # Sort coefficients by t value (to show crossover pattern)
    sorted_coeffs = sorted(coeff_data.keys(), key=lambda x: (x[2], x[0] + x[1]))

    # Prepare plot data
    labels = []
    tmp_means = []
    tmp_stds = []
    sym_means = []
    sym_stds = []

    for coeff in sorted_coeffs:
        tmp_times = coeff_data[coeff]['TMP']
        sym_times = coeff_data[coeff]['Symbolic']

        if tmp_times and sym_times:
            labels.append(f'$E^{{{coeff[0]},{coeff[1]}}}_{{{coeff[2]}}}$')
            tmp_means.append(statistics.mean(tmp_times))
            tmp_stds.append(statistics.stdev(tmp_times) if len(tmp_times) > 1 else 0)
            sym_means.append(statistics.mean(sym_times))
            sym_stds.append(statistics.stdev(sym_times) if len(sym_times) > 1 else 0)

    if not labels:
        print("Warning: No matching TMP/Symbolic pairs found")
        return None

    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width, fig_width * 0.7))

    x = np.arange(len(labels))
    width = 0.35

    # Bar plots with error bars
    bars_tmp = ax.bar(x - width/2, tmp_means, width,
                      yerr=tmp_stds,
                      label='TMP (Template Metaprogramming)',
                      color=IMPL_COLORS['TMP'],
                      edgecolor='black',
                      linewidth=0.5,
                      capsize=2,
                      error_kw={'linewidth': 0.8})

    bars_sym = ax.bar(x + width/2, sym_means, width,
                      yerr=sym_stds,
                      label='Symbolic (Code Generation)',
                      color=IMPL_COLORS['Symbolic'],
                      edgecolor='black',
                      linewidth=0.5,
                      capsize=2,
                      error_kw={'linewidth': 0.8})

    # Formatting
    ax.set_xlabel('Hermite E Coefficient')
    ax.set_ylabel('Execution Time (ns)')
    ax.set_title('TMP vs Symbolic Performance Comparison\n(Actual Measurements)', fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.legend(loc='upper left', fontsize=7)
    ax.grid(True, axis='y', alpha=0.3, linewidth=0.5)

    # Add minor gridlines
    ax.yaxis.set_minor_locator(plt.AutoMinorLocator())

    plt.tight_layout()

    # Save in multiple formats
    output_base = Path(output_dir) / 'fig1_performance_comparison'
    fig.savefig(f'{output_base}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_base}.pdf and .png")

    plt.close()
    return str(output_base)


def plot_speedup_analysis(benchmarks, output_dir, journal='aps'):
    """
    Plot 2: Speedup Analysis with Crossover Point

    Shows speedup ratio (Symbolic time / TMP time) highlighting crossover.
    Speedup > 1 means TMP is faster; < 1 means Symbolic is faster.
    """
    fig_width = COLUMN_WIDTHS[journal]['onehalf']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal=journal, font_scale=font_scale)

    # Group data by coefficient and implementation
    coeff_data = defaultdict(lambda: {'TMP': [], 'Symbolic': []})

    for bench in benchmarks:
        if 'Compare' in bench['name'] or 'E_' in bench['name']:
            nA, nB, t = extract_coefficient_from_name(bench['name'])
            if nA is not None:
                impl_name = 'TMP' if bench['impl'] == 0 else 'Symbolic' if bench['impl'] == 1 else None
                if impl_name:
                    coeff_data[(nA, nB, t)][impl_name].append(bench['real_time'])

    if not coeff_data:
        print("Warning: No comparison data found for speedup plot")
        return None

    # Calculate speedups (Symbolic/TMP - so >1 means TMP wins)
    speedup_data = {}
    for coeff, data in coeff_data.items():
        if data['TMP'] and data['Symbolic']:
            tmp_mean = statistics.mean(data['TMP'])
            sym_mean = statistics.mean(data['Symbolic'])
            if tmp_mean > 0:
                speedup_data[coeff] = sym_mean / tmp_mean

    # Sort by t value to show crossover
    sorted_coeffs = sorted(speedup_data.keys(), key=lambda x: (x[2], x[0] + x[1]))

    labels = [f'$E^{{{c[0]},{c[1]}}}_{{{c[2]}}}$' for c in sorted_coeffs]
    speedups = [speedup_data[c] for c in sorted_coeffs]
    t_values = [c[2] for c in sorted_coeffs]

    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width, fig_width * 0.65))

    x = np.arange(len(labels))

    # Color bars based on which method wins
    colors = [IMPL_COLORS['TMP'] if s >= 1 else IMPL_COLORS['Symbolic'] for s in speedups]

    bars = ax.bar(x, speedups, color=colors, edgecolor='black', linewidth=0.5)

    # Add reference line at speedup = 1
    ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=1.5,
               label='Equal Performance', zorder=0)

    # Find and highlight crossover region
    crossover_indices = []
    for i in range(len(speedups) - 1):
        if (speedups[i] < 1 and speedups[i+1] >= 1) or (speedups[i] >= 1 and speedups[i+1] < 1):
            crossover_indices.append(i)

    # Shade crossover region
    for idx in crossover_indices:
        ax.axvspan(idx + 0.5, idx + 1.5, alpha=0.2, color='gray',
                   label='Crossover Region' if idx == crossover_indices[0] else None)

    # Annotate speedup values
    for i, (bar, speedup) in enumerate(zip(bars, speedups)):
        height = bar.get_height()
        va = 'bottom' if height >= 1 else 'top'
        offset = 0.05 if height >= 1 else -0.05
        ax.annotate(f'{speedup:.2f}x',
                    xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3 if height >= 1 else -3),
                    textcoords='offset points',
                    ha='center', va=va,
                    fontsize=6, rotation=0)

    # Formatting
    ax.set_xlabel('Hermite E Coefficient (ordered by $t$)')
    ax.set_ylabel('Speedup Ratio (Symbolic Time / TMP Time)')
    ax.set_title('Speedup Analysis: TMP vs Symbolic\n(Actual Measurements)', fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha='right')
    ax.legend(loc='upper right', fontsize=7)
    ax.grid(True, axis='y', alpha=0.3, linewidth=0.5)

    # Add annotation explaining interpretation
    ax.text(0.02, 0.98, 'Speedup > 1: TMP faster\nSpeedup < 1: Symbolic faster',
            transform=ax.transAxes, fontsize=6,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8))

    plt.tight_layout()

    # Save
    output_base = Path(output_dir) / 'fig2_speedup_analysis'
    fig.savefig(f'{output_base}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_base}.pdf and .png")

    plt.close()
    return str(output_base)


def plot_scaling_with_L(benchmarks, output_dir, journal='aps'):
    """
    Plot 3: Performance Scaling with Angular Momentum

    Shows how execution time scales with L_total = nA + nB.
    Uses ACTUAL benchmark data from scaling tests.
    """
    fig_width = COLUMN_WIDTHS[journal]['onehalf']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal=journal, font_scale=font_scale)

    # Group data by implementation and L_total
    scaling_data = defaultdict(lambda: defaultdict(list))

    for bench in benchmarks:
        if 'Scaling' in bench['name']:
            impl = bench['impl']
            L = bench['L_total']
            if impl in [0, 1] and L >= 0:
                impl_name = 'TMP' if impl == 0 else 'Symbolic'
                scaling_data[impl_name][L].append(bench['real_time'])

    if not scaling_data:
        print("Warning: No scaling data found")
        return None

    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width, fig_width * 0.75))

    # Plot each implementation
    for impl_name in ['TMP', 'Symbolic']:
        if impl_name not in scaling_data:
            continue

        L_values = sorted(scaling_data[impl_name].keys())
        means = []
        stds = []

        for L in L_values:
            times = scaling_data[impl_name][L]
            means.append(statistics.mean(times))
            stds.append(statistics.stdev(times) if len(times) > 1 else 0)

        ax.errorbar(L_values, means, yerr=stds,
                    marker=MARKERS[impl_name],
                    color=IMPL_COLORS[impl_name],
                    linewidth=1.5,
                    markersize=7,
                    capsize=3,
                    capthick=1,
                    label=impl_name,
                    markeredgecolor='black',
                    markeredgewidth=0.5)

    # Formatting
    ax.set_xlabel(r'Total Angular Momentum $L_{total} = n_A + n_B$')
    ax.set_ylabel('Execution Time (ns)')
    ax.set_title('Performance Scaling with Angular Momentum\n(Actual Measurements)', fontsize=10)
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3, linewidth=0.5)

    # Use log scale if data spans multiple orders of magnitude
    all_means = []
    for impl_name in scaling_data:
        for L in scaling_data[impl_name]:
            all_means.extend(scaling_data[impl_name][L])

    if all_means and max(all_means) / min(all_means) > 10:
        ax.set_yscale('log')
        ax.set_ylabel('Execution Time (ns, log scale)')

    # Set integer ticks for L values
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    plt.tight_layout()

    # Save
    output_base = Path(output_dir) / 'fig3_scaling_with_L'
    fig.savefig(f'{output_base}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_base}.pdf and .png")

    plt.close()
    return str(output_base)


def plot_combined_figure(benchmarks, output_dir, journal='aps'):
    """
    Generate a combined multi-panel figure with all three plots.
    Following scientific_plotting_agent_guide.md for multi-panel layout.
    """
    fig_width = COLUMN_WIDTHS[journal]['double']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal=journal, font_scale=font_scale)

    # Create 1x3 subplot layout
    fig, axes = plt.subplots(1, 3, figsize=(fig_width, fig_width * 0.35))

    # =========================================================================
    # Panel (a): Performance Comparison
    # =========================================================================
    ax = axes[0]

    coeff_data = defaultdict(lambda: {'TMP': [], 'Symbolic': []})
    for bench in benchmarks:
        if 'Compare' in bench['name'] or 'E_' in bench['name']:
            nA, nB, t = extract_coefficient_from_name(bench['name'])
            if nA is not None:
                impl_name = 'TMP' if bench['impl'] == 0 else 'Symbolic' if bench['impl'] == 1 else None
                if impl_name:
                    coeff_data[(nA, nB, t)][impl_name].append(bench['real_time'])

    sorted_coeffs = sorted(coeff_data.keys(), key=lambda x: (x[2], x[0] + x[1]))

    labels = []
    tmp_means = []
    sym_means = []
    for coeff in sorted_coeffs:
        if coeff_data[coeff]['TMP'] and coeff_data[coeff]['Symbolic']:
            labels.append(f'$E^{{{coeff[0]},{coeff[1]}}}_{{{coeff[2]}}}$')
            tmp_means.append(statistics.mean(coeff_data[coeff]['TMP']))
            sym_means.append(statistics.mean(coeff_data[coeff]['Symbolic']))

    if labels:
        x = np.arange(len(labels))
        width = 0.35
        ax.bar(x - width/2, tmp_means, width, label='TMP', color=IMPL_COLORS['TMP'],
               edgecolor='black', linewidth=0.3)
        ax.bar(x + width/2, sym_means, width, label='Symbolic', color=IMPL_COLORS['Symbolic'],
               edgecolor='black', linewidth=0.3)
        ax.set_xlabel('Coefficient')
        ax.set_ylabel('Time (ns)')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=6, rotation=45, ha='right')
        ax.legend(loc='upper left', fontsize=6)
        ax.grid(True, axis='y', alpha=0.3)

    add_panel_label(ax, '(a)')

    # =========================================================================
    # Panel (b): Speedup Analysis
    # =========================================================================
    ax = axes[1]

    speedup_data = {}
    for coeff, data in coeff_data.items():
        if data['TMP'] and data['Symbolic']:
            tmp_mean = statistics.mean(data['TMP'])
            sym_mean = statistics.mean(data['Symbolic'])
            if tmp_mean > 0:
                speedup_data[coeff] = sym_mean / tmp_mean

    sorted_coeffs = sorted(speedup_data.keys(), key=lambda x: (x[2], x[0] + x[1]))
    labels = [f'$E^{{{c[0]},{c[1]}}}_{{{c[2]}}}$' for c in sorted_coeffs]
    speedups = [speedup_data[c] for c in sorted_coeffs]

    if speedups:
        x = np.arange(len(labels))
        colors = [IMPL_COLORS['TMP'] if s >= 1 else IMPL_COLORS['Symbolic'] for s in speedups]
        ax.bar(x, speedups, color=colors, edgecolor='black', linewidth=0.3)
        ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=1)
        ax.set_xlabel('Coefficient')
        ax.set_ylabel('Speedup (Sym/TMP)')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=6, rotation=45, ha='right')
        ax.grid(True, axis='y', alpha=0.3)

    add_panel_label(ax, '(b)')

    # =========================================================================
    # Panel (c): Scaling with L
    # =========================================================================
    ax = axes[2]

    scaling_data = defaultdict(lambda: defaultdict(list))
    for bench in benchmarks:
        if 'Scaling' in bench['name']:
            impl = bench['impl']
            L = bench['L_total']
            if impl in [0, 1] and L >= 0:
                impl_name = 'TMP' if impl == 0 else 'Symbolic'
                scaling_data[impl_name][L].append(bench['real_time'])

    for impl_name in ['TMP', 'Symbolic']:
        if impl_name not in scaling_data:
            continue
        L_values = sorted(scaling_data[impl_name].keys())
        means = [statistics.mean(scaling_data[impl_name][L]) for L in L_values]
        stds = [statistics.stdev(scaling_data[impl_name][L]) if len(scaling_data[impl_name][L]) > 1 else 0
                for L in L_values]
        ax.errorbar(L_values, means, yerr=stds,
                    marker=MARKERS[impl_name], color=IMPL_COLORS[impl_name],
                    linewidth=1.2, markersize=5, capsize=2, label=impl_name)

    ax.set_xlabel(r'$L_{total}$')
    ax.set_ylabel('Time (ns)')
    ax.legend(loc='upper left', fontsize=6)
    ax.grid(True, alpha=0.3)

    # Check if log scale needed
    all_means = []
    for impl_name in scaling_data:
        for L in scaling_data[impl_name]:
            all_means.extend(scaling_data[impl_name][L])
    if all_means and max(all_means) / min(all_means) > 10:
        ax.set_yscale('log')

    add_panel_label(ax, '(c)')

    # Finalize
    plt.tight_layout()

    output_base = Path(output_dir) / 'fig_combined_comparison'
    fig.savefig(f'{output_base}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_base}.pdf and .png")

    plt.close()
    return str(output_base)


def print_data_summary(benchmarks, context):
    """Print summary of loaded benchmark data."""
    print("\n" + "=" * 70)
    print("BENCHMARK DATA SUMMARY (Actual Measurements)")
    print("=" * 70)

    if context:
        print(f"\nSystem: {context.get('host_name', 'Unknown')}")
        print(f"Date: {context.get('date', 'Unknown')}")
        print(f"CPUs: {context.get('num_cpus', '?')}")
        print(f"Build: {context.get('library_build_type', 'Unknown')}")

    # Count by implementation
    impl_counts = defaultdict(int)
    impl_names = {0: 'TMP', 1: 'Symbolic', 2: 'Naive', 3: 'BottomUp'}
    for b in benchmarks:
        impl_counts[b['impl']] += 1

    print(f"\nTotal measurements: {len(benchmarks)}")
    for impl, count in sorted(impl_counts.items()):
        print(f"  {impl_names.get(impl, f'Impl {impl}')}: {count}")

    # Check for scaling and comparison data
    scaling_count = sum(1 for b in benchmarks if 'Scaling' in b['name'])
    compare_count = sum(1 for b in benchmarks if 'Compare' in b['name'] or 'E_' in b['name'])
    print(f"\nScaling benchmarks: {scaling_count}")
    print(f"Comparison benchmarks: {compare_count}")
    print("=" * 70 + "\n")


# =============================================================================
# MAIN
# =============================================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate publication-quality plots from ACTUAL benchmark data'
    )
    parser.add_argument('json_file', help='Path to benchmark JSON file')
    parser.add_argument('--output-dir', '-o', default=None,
                        help='Output directory for figures')
    parser.add_argument('--journal', '-j', default='aps',
                        choices=['nature', 'science', 'aps', 'acs', 'elsevier'],
                        help='Target journal for figure dimensions')

    args = parser.parse_args()

    # Set default output directory
    if args.output_dir is None:
        args.output_dir = str(Path(args.json_file).parent.parent / 'figures')

    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading benchmark data from: {args.json_file}")
    context, benchmarks = load_benchmark_data(args.json_file)

    # Print summary
    print_data_summary(benchmarks, context)

    # Generate plots
    print("Generating publication-quality figures...")
    print(f"Target journal format: {args.journal}")
    print(f"Output directory: {args.output_dir}\n")

    # Individual plots
    plot_performance_comparison(benchmarks, args.output_dir, args.journal)
    plot_speedup_analysis(benchmarks, args.output_dir, args.journal)
    plot_scaling_with_L(benchmarks, args.output_dir, args.journal)

    # Combined multi-panel figure
    plot_combined_figure(benchmarks, args.output_dir, args.journal)

    print("\nAll figures generated successfully!")
    print("NOTE: All plots use ACTUAL benchmark measurements only.")


if __name__ == '__main__':
    main()
