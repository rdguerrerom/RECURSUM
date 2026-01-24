#!/usr/bin/env python3
"""
Publication-quality benchmark comparison plots for McMurchie-Davidson algorithm components.

Generates figures following scientific plotting guide standards:
- Colorblind-safe Okabe-Ito palette
- Physical Review style formatting
- 300 DPI output in PDF and PNG formats

Author: Generated with Claude Code
Date: 2026-01-15
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from collections import defaultdict

# =============================================================================
# Configuration following scientific plotting guide
# =============================================================================

# Journal column widths in inches
COLUMN_WIDTHS = {
    'nature': {'single': 3.46, 'onehalf': 4.76, 'double': 7.09},
    'science': {'single': 3.35, 'onehalf': 4.57, 'double': 7.01},
    'aps': {'single': 3.39, 'onehalf': 5.12, 'double': 7.01},
    'acs': {'single': 3.33, 'onehalf': 4.33, 'double': 7.00},
    'elsevier': {'single': 3.54, 'onehalf': 5.51, 'double': 7.48},
}

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

# Line styles for additional differentiation
LINE_STYLES = ['-', '--', '-.', ':']
MARKERS = ['o', 's', '^', 'D', 'v', '<', '>', 'p']


def get_figsize(journal='aps', width='single', aspect=0.75):
    """Get figure dimensions for target journal."""
    w = COLUMN_WIDTHS[journal][width]
    h = w * aspect
    return (w, h)


def get_font_scale(fig_width_inches, target_min_font=6):
    """Calculate font scale factor based on figure width."""
    reference_width = 7.0
    scale = reference_width / fig_width_inches
    return max(1.0, min(scale, 2.5))


def configure_publication_style(journal='aps', font_scale=1.0):
    """Configure matplotlib for publication-quality figures."""
    base_sizes = {
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 8,
        'figure.titlesize': 12,
        'font.size': 9,
    }

    scaled_sizes = {k: v * font_scale for k, v in base_sizes.items()}

    font_families = {
        'nature': 'Helvetica',
        'science': 'Arial',
        'aps': 'Helvetica',
        'acs': 'Arial',
        'elsevier': 'Arial',
    }

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': [font_families.get(journal, 'Arial'), 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',
        **scaled_sizes,
        'axes.linewidth': 0.8,
        'axes.labelpad': 4,
        'axes.spines.top': False,
        'axes.spines.right': False,
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
        'legend.frameon': False,
        'legend.borderpad': 0.4,
        'legend.handlelength': 1.5,
        'legend.handletextpad': 0.5,
        'legend.columnspacing': 1.0,
        'lines.linewidth': 1.5,
        'lines.markersize': 5,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
    })

    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=OKABE_ITO_CYCLE)


def save_publication_figure(fig, filename, formats=('pdf', 'png'), dpi=300):
    """Save figure in multiple formats for publication."""
    for fmt in formats:
        output_path = f'{filename}.{fmt}'
        kwargs = {
            'bbox_inches': 'tight',
            'pad_inches': 0.02,
            'transparent': False,
        }
        if fmt in ('png', 'tiff'):
            kwargs['dpi'] = dpi
        elif fmt == 'pdf':
            kwargs['dpi'] = dpi

        fig.savefig(output_path, format=fmt, **kwargs)
        print(f'Saved: {output_path}')


# =============================================================================
# Data Loading and Parsing
# =============================================================================

def load_benchmark_data(json_path):
    """Load Google Benchmark JSON and extract aggregate statistics."""
    with open(json_path, 'r') as f:
        data = json.load(f)

    # Extract only aggregate entries (mean and stddev)
    aggregates = defaultdict(dict)

    for entry in data['benchmarks']:
        if entry.get('run_type') != 'aggregate':
            continue

        name = entry['run_name']
        agg_name = entry.get('aggregate_name', '')

        if agg_name in ('mean', 'stddev'):
            if name not in aggregates:
                aggregates[name] = {}
            aggregates[name][agg_name] = entry

    return aggregates


def parse_hermite_data(aggregates):
    """Parse Hermite coefficient benchmark data into structured format."""
    shell_pairs = ['ss', 'sp', 'pp', 'sd', 'pd', 'dd', 'ff', 'gg']
    impl_names = {0: 'TMP', 1: 'Layered', 2: 'Symbolic'}

    results = {impl: {'shells': [], 'L': [], 'mean': [], 'std': [], 'n_coeffs': []}
               for impl in impl_names.values()}

    for name, agg_data in aggregates.items():
        if 'mean' not in agg_data:
            continue

        mean_entry = agg_data['mean']
        std_entry = agg_data.get('stddev', {})

        # Extract implementation and shell pair from name
        # Format: HermiteE/TMP/ss/min_time:1.000
        parts = name.split('/')
        if len(parts) < 3:
            continue

        impl_str = parts[1]
        shell = parts[2]

        if shell not in shell_pairs:
            continue

        impl_id = int(mean_entry.get('impl', -1))
        if impl_id not in impl_names:
            continue

        impl_name = impl_names[impl_id]

        results[impl_name]['shells'].append(shell)
        results[impl_name]['L'].append(int(mean_entry.get('L', 0)))
        results[impl_name]['mean'].append(mean_entry['cpu_time'])
        results[impl_name]['std'].append(std_entry.get('cpu_time', 0))
        results[impl_name]['n_coeffs'].append(int(mean_entry.get('n_coeffs', 1)))

    # Sort by shell pair order
    for impl in results:
        if results[impl]['shells']:
            indices = [shell_pairs.index(s) for s in results[impl]['shells']]
            sorted_idx = np.argsort(indices)
            for key in results[impl]:
                results[impl][key] = [results[impl][key][i] for i in sorted_idx]

    return results, shell_pairs


def parse_coulomb_data(aggregates):
    """Parse Coulomb auxiliary integral benchmark data into structured format."""
    impl_names = {0: 'TMP', 1: 'Layered'}

    results = {impl: {'L_total': [], 'mean': [], 'std': [], 'n_integrals': []}
               for impl in impl_names.values()}

    for name, agg_data in aggregates.items():
        if 'mean' not in agg_data:
            continue

        mean_entry = agg_data['mean']
        std_entry = agg_data.get('stddev', {})

        impl_id = int(mean_entry.get('impl', -1))
        if impl_id not in impl_names:
            continue

        impl_name = impl_names[impl_id]

        results[impl_name]['L_total'].append(int(mean_entry.get('L_total', 0)))
        results[impl_name]['mean'].append(mean_entry['cpu_time'])
        results[impl_name]['std'].append(std_entry.get('cpu_time', 0))
        results[impl_name]['n_integrals'].append(int(mean_entry.get('n_integrals', 1)))

    # Sort by L_total
    for impl in results:
        if results[impl]['L_total']:
            sorted_idx = np.argsort(results[impl]['L_total'])
            for key in results[impl]:
                results[impl][key] = [results[impl][key][i] for i in sorted_idx]

    return results


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_hermite_bar_comparison(hermite_data, shell_pairs, output_dir):
    """
    Plot 1: Bar chart comparing TMP, Layered, and Symbolic implementations
    for Hermite expansion coefficients.
    """
    fig_width = COLUMN_WIDTHS['aps']['onehalf']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal='aps', font_scale=font_scale)

    fig, ax = plt.subplots(figsize=(fig_width, fig_width * 0.7))

    # Prepare data - ensure all implementations have same shells
    implementations = ['TMP', 'Layered', 'Symbolic']
    colors = OKABE_ITO_CYCLE[:3]

    # Get shells that exist in all implementations
    available_shells = []
    for shell in shell_pairs:
        if all(shell in hermite_data[impl]['shells'] for impl in implementations):
            available_shells.append(shell)

    n_shells = len(available_shells)
    n_impl = len(implementations)
    bar_width = 0.25
    x = np.arange(n_shells)

    for i, (impl, color) in enumerate(zip(implementations, colors)):
        means = []
        stds = []
        for shell in available_shells:
            idx = hermite_data[impl]['shells'].index(shell)
            means.append(hermite_data[impl]['mean'][idx])
            stds.append(hermite_data[impl]['std'][idx])

        offset = (i - n_impl/2 + 0.5) * bar_width
        bars = ax.bar(x + offset, means, bar_width, yerr=stds,
                      label=impl, color=color,
                      capsize=2, error_kw={'linewidth': 0.8})

    ax.set_xlabel('Shell pair')
    ax.set_ylabel('Execution time (ns)')
    ax.set_title(r'Hermite Expansion Coefficients $E_N^{i,j}$')
    ax.set_xticks(x)
    ax.set_xticklabels(available_shells)
    ax.legend(loc='upper left')

    # Use log scale for better visibility of differences
    ax.set_yscale('log')
    ax.set_ylim(bottom=0.1)

    fig.tight_layout()

    output_path = output_dir / 'hermite_coefficients_comparison'
    save_publication_figure(fig, str(output_path))
    plt.close(fig)

    return str(output_path)


def plot_hermite_scaling(hermite_data, output_dir):
    """
    Plot 2: Line plot showing timing vs total angular momentum L
    with log scale on y-axis.
    """
    fig_width = COLUMN_WIDTHS['aps']['onehalf']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal='aps', font_scale=font_scale)

    fig, ax = plt.subplots(figsize=(fig_width, fig_width * 0.7))

    implementations = ['TMP', 'Layered', 'Symbolic']
    colors = OKABE_ITO_CYCLE[:3]
    markers = MARKERS[:3]
    linestyles = ['-', '--', '-.']

    for impl, color, marker, ls in zip(implementations, colors, markers, linestyles):
        L_vals = np.array(hermite_data[impl]['L'])
        means = np.array(hermite_data[impl]['mean'])
        stds = np.array(hermite_data[impl]['std'])

        # Sort by L
        sorted_idx = np.argsort(L_vals)
        L_vals = L_vals[sorted_idx]
        means = means[sorted_idx]
        stds = stds[sorted_idx]

        ax.errorbar(L_vals, means, yerr=stds,
                    label=impl, color=color, marker=marker,
                    linestyle=ls, markersize=6, capsize=2,
                    markerfacecolor=color, markeredgecolor='white',
                    markeredgewidth=0.5, linewidth=1.5)

    ax.set_xlabel(r'Total angular momentum $L = n_A + n_B$')
    ax.set_ylabel('Execution time (ns)')
    ax.set_title('Scaling with Angular Momentum')
    ax.set_yscale('log')
    ax.legend(loc='upper left')

    # Add minor gridlines for log scale
    ax.grid(True, which='major', linestyle='-', alpha=0.3)
    ax.grid(True, which='minor', linestyle=':', alpha=0.2)

    # Set x-axis to integer ticks
    ax.set_xticks(np.arange(0, max(hermite_data['TMP']['L']) + 1, 1))

    fig.tight_layout()

    output_path = output_dir / 'hermite_coefficients_vs_L'
    save_publication_figure(fig, str(output_path))
    plt.close(fig)

    return str(output_path)


def plot_coulomb_bar_comparison(coulomb_data, output_dir):
    """
    Plot 3: Bar chart comparing TMP and Layered implementations
    for Coulomb auxiliary integrals.
    """
    fig_width = COLUMN_WIDTHS['aps']['onehalf']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal='aps', font_scale=font_scale)

    fig, ax = plt.subplots(figsize=(fig_width, fig_width * 0.7))

    implementations = ['TMP', 'Layered']
    colors = OKABE_ITO_CYCLE[:2]

    L_values = sorted(set(coulomb_data['TMP']['L_total']))
    n_L = len(L_values)
    n_impl = len(implementations)
    bar_width = 0.35
    x = np.arange(n_L)

    for i, (impl, color) in enumerate(zip(implementations, colors)):
        means = []
        stds = []
        for L in L_values:
            idx = coulomb_data[impl]['L_total'].index(L)
            means.append(coulomb_data[impl]['mean'][idx])
            stds.append(coulomb_data[impl]['std'][idx])

        offset = (i - n_impl/2 + 0.5) * bar_width
        bars = ax.bar(x + offset, means, bar_width, yerr=stds,
                      label=impl, color=color,
                      capsize=2, error_kw={'linewidth': 0.8})

    ax.set_xlabel(r'Total angular momentum $L_{total}$')
    ax.set_ylabel('Execution time (ns)')
    ax.set_title(r'Coulomb Auxiliary Integrals $R_{tuv}$')
    ax.set_xticks(x)
    ax.set_xticklabels([str(L) for L in L_values])
    ax.legend(loc='upper left')

    # Use log scale
    ax.set_yscale('log')

    fig.tight_layout()

    output_path = output_dir / 'coulomb_hermite_comparison'
    save_publication_figure(fig, str(output_path))
    plt.close(fig)

    return str(output_path)


def plot_coulomb_scaling(coulomb_data, output_dir):
    """
    Plot 4: Line plot showing cost per integral vs L_total.
    """
    fig_width = COLUMN_WIDTHS['aps']['onehalf']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal='aps', font_scale=font_scale)

    fig, ax = plt.subplots(figsize=(fig_width, fig_width * 0.7))

    implementations = ['TMP', 'Layered']
    colors = OKABE_ITO_CYCLE[:2]
    markers = MARKERS[:2]
    linestyles = ['-', '--']

    for impl, color, marker, ls in zip(implementations, colors, markers, linestyles):
        L_vals = np.array(coulomb_data[impl]['L_total'])
        means = np.array(coulomb_data[impl]['mean'])
        stds = np.array(coulomb_data[impl]['std'])
        n_integrals = np.array(coulomb_data[impl]['n_integrals'])

        # Calculate time per integral
        time_per_integral = means / n_integrals
        std_per_integral = stds / n_integrals

        # Sort by L
        sorted_idx = np.argsort(L_vals)
        L_vals = L_vals[sorted_idx]
        time_per_integral = time_per_integral[sorted_idx]
        std_per_integral = std_per_integral[sorted_idx]

        ax.errorbar(L_vals, time_per_integral, yerr=std_per_integral,
                    label=impl, color=color, marker=marker,
                    linestyle=ls, markersize=6, capsize=2,
                    markerfacecolor=color, markeredgecolor='white',
                    markeredgewidth=0.5, linewidth=1.5)

    ax.set_xlabel(r'Total angular momentum $L_{total}$')
    ax.set_ylabel('Time per integral (ns/integral)')
    ax.set_title('Cost per Integral vs Angular Momentum')

    # Use log scale to show both implementations clearly (orders of magnitude difference)
    ax.set_yscale('log')
    ax.legend(loc='upper right')

    # Add gridlines for log scale
    ax.grid(True, which='major', linestyle='-', alpha=0.3)
    ax.grid(True, which='minor', linestyle=':', alpha=0.2)

    # Set x-axis to integer ticks
    max_L = max(coulomb_data['TMP']['L_total'])
    ax.set_xticks(np.arange(0, max_L + 1, 1))

    fig.tight_layout()

    output_path = output_dir / 'coulomb_hermite_scaling'
    save_publication_figure(fig, str(output_path))
    plt.close(fig)

    return str(output_path)


# =============================================================================
# Main Execution
# =============================================================================

def main():
    """Generate all benchmark comparison plots."""
    # Paths
    base_dir = Path('/home/ruben/Research/Science/Projects/RECURSUM/benchmarks')
    hermite_json = base_dir / 'results/raw/hermite_coefficients.json'
    coulomb_json = base_dir / 'results/raw/coulomb_hermite.json'
    output_dir = base_dir / 'results/figures'

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("McMurchie-Davidson Benchmark Comparison Plots")
    print("=" * 60)

    # Load data
    print("\nLoading Hermite coefficient benchmarks...")
    hermite_aggregates = load_benchmark_data(hermite_json)
    hermite_data, shell_pairs = parse_hermite_data(hermite_aggregates)
    print(f"  Found {len(hermite_aggregates)} benchmark configurations")

    print("\nLoading Coulomb auxiliary integral benchmarks...")
    coulomb_aggregates = load_benchmark_data(coulomb_json)
    coulomb_data = parse_coulomb_data(coulomb_aggregates)
    print(f"  Found {len(coulomb_aggregates)} benchmark configurations")

    # Generate plots
    print("\n" + "-" * 60)
    print("Generating plots...")
    print("-" * 60)

    print("\nPlot 1: Hermite coefficients bar comparison")
    path1 = plot_hermite_bar_comparison(hermite_data, shell_pairs, output_dir)

    print("\nPlot 2: Hermite coefficients scaling with L")
    path2 = plot_hermite_scaling(hermite_data, output_dir)

    print("\nPlot 3: Coulomb integrals bar comparison")
    path3 = plot_coulomb_bar_comparison(coulomb_data, output_dir)

    print("\nPlot 4: Coulomb integrals scaling")
    path4 = plot_coulomb_scaling(coulomb_data, output_dir)

    print("\n" + "=" * 60)
    print("All plots generated successfully!")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir}")
    print("\nGenerated files:")
    print(f"  - {path1}.pdf/.png")
    print(f"  - {path2}.pdf/.png")
    print(f"  - {path3}.pdf/.png")
    print(f"  - {path4}.pdf/.png")


if __name__ == '__main__':
    main()
