#!/usr/bin/env python3
"""
Publication-quality plots for RECURSUM Hermite E coefficient benchmark comparison.
TMP (Template Metaprogramming) vs Symbolic (SymPy-generated) implementations.

Follows scientific plotting guidelines for Nature/Science/APS journals.
All data extracted from actual benchmark JSON - no estimated values.

Author: Scientific Plotting Agent
Date: 2026-01-14
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

# ==============================================================================
# Publication Style Configuration (from scientific_plotting_agent_guide.md)
# ==============================================================================

# Journal column widths in inches
COLUMN_WIDTHS = {
    'nature': {'single': 3.46, 'onehalf': 4.76, 'double': 7.09},
    'aps': {'single': 3.39, 'onehalf': 5.12, 'double': 7.01},
}

# Okabe-Ito colorblind-safe palette
OKABE_ITO = {
    'blue': '#0072B2',
    'orange': '#E69F00',
    'green': '#009E73',
    'purple': '#CC79A7',
    'vermillion': '#D55E00',
    'sky_blue': '#56B4E9',
    'yellow': '#F0E442',
    'black': '#000000',
}

# Color assignments for implementations
COLORS = {
    'TMP': OKABE_ITO['blue'],
    'Symbolic': OKABE_ITO['orange'],
}


def configure_publication_style(journal='nature', font_scale=1.0):
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

    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'DejaVu Serif', 'serif'],
        'mathtext.fontset': 'cm',  # Computer Modern for math

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

        'lines.linewidth': 1.5,
        'lines.markersize': 5,

        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
    })


def get_figsize(journal='nature', width='single', aspect=0.75):
    """Get figure dimensions for target journal."""
    w = COLUMN_WIDTHS[journal][width]
    h = w * aspect
    return (w, h)


def get_font_scale(fig_width_inches):
    """Calculate font scale factor based on figure width."""
    reference_width = 7.0
    scale = reference_width / fig_width_inches
    return max(1.0, min(scale, 2.5))


def add_panel_label(ax, label, fontsize=12, fontweight='bold', offset=(-0.12, 1.05)):
    """Add panel label (a), (b), etc. to axes."""
    ax.text(offset[0], offset[1], label, transform=ax.transAxes,
            fontsize=fontsize, fontweight=fontweight, va='top', ha='left')


# ==============================================================================
# Data Loading and Parsing
# ==============================================================================

def load_benchmark_data(json_path):
    """Load and parse benchmark JSON file.

    Handles NaN values which are invalid JSON by replacing them with null.
    """
    import re
    with open(json_path, 'r') as f:
        content = f.read()
    # Replace NaN values with null (valid JSON)
    content = re.sub(r':\s*-?nan\b', ': null', content, flags=re.IGNORECASE)
    data = json.loads(content)
    return data


def extract_comparison_benchmarks(data):
    """
    Extract comparison benchmark data (TMP vs Symbolic for specific coefficients).

    Returns dict with structure:
    {
        (nA, nB, t): {
            'TMP': {'mean': float, 'stddev': float},
            'Symbolic': {'mean': float, 'stddev': float}
        }
    }
    """
    results = defaultdict(lambda: defaultdict(dict))

    for bench in data['benchmarks']:
        name = bench.get('name', '')

        # Filter for comparison benchmarks with mean/stddev aggregates
        if 'BM_Compare' not in name:
            continue

        # Get aggregate type
        if '_mean' in name:
            agg_type = 'mean'
        elif '_stddev' in name:
            agg_type = 'stddev'
        else:
            continue

        # Extract implementation type and indices
        impl_label = bench.get('label', '')
        if impl_label not in ['TMP', 'Symbolic']:
            # Infer from name
            if 'TMP' in name:
                impl_label = 'TMP'
            elif 'Symbolic' in name:
                impl_label = 'Symbolic'
            else:
                continue

        nA = int(bench.get('nA', 0))
        nB = int(bench.get('nB', 0))
        t = int(bench.get('t', 0))

        real_time = bench.get('real_time', 0)

        results[(nA, nB, t)][impl_label][agg_type] = real_time

    return dict(results)


def extract_scaling_benchmarks(data):
    """
    Extract scaling benchmark data (performance vs L_total).

    Returns dict with structure:
    {
        L_total: {
            'TMP': {'mean': float, 'stddev': float},
            'Symbolic': {'mean': float, 'stddev': float}
        }
    }
    """
    results = defaultdict(lambda: defaultdict(dict))

    for bench in data['benchmarks']:
        name = bench.get('name', '')

        if 'BM_Scaling' not in name:
            continue

        # Get aggregate type
        if '_mean' in name:
            agg_type = 'mean'
        elif '_stddev' in name:
            agg_type = 'stddev'
        else:
            continue

        # Determine implementation
        if 'TMP' in name:
            impl_label = 'TMP'
        elif 'Symbolic' in name:
            impl_label = 'Symbolic'
        else:
            continue

        L_total = int(bench.get('L_total', 0))
        real_time = bench.get('real_time', 0)

        results[L_total][impl_label][agg_type] = real_time

    return dict(results)


# ==============================================================================
# Plotting Functions
# ==============================================================================

def plot_execution_time_comparison(comparison_data, output_dir, journal='nature'):
    """
    Plot 1: Bar chart comparing TMP vs Symbolic execution times.

    Shows key coefficients: E^{0,0}_0, E^{1,1}_0, E^{1,1}_2, E^{2,2}_2, E^{2,2}_4,
    and E^{3,3}_t for t=0,1,2,3,4,5
    """
    # Select representative coefficients for visualization
    # Focus on diagonal cases (nA=nB) which are most commonly used
    selected_coefficients = [
        (0, 0, 0),  # E^{0,0}_0 - simplest
        (1, 1, 0),  # E^{1,1}_0
        (1, 1, 2),  # E^{1,1}_2 - max t for (1,1)
        (2, 2, 0),  # E^{2,2}_0
        (2, 2, 2),  # E^{2,2}_2
        (2, 2, 4),  # E^{2,2}_4 - max t for (2,2)
        (3, 3, 0),  # E^{3,3}_0
        (3, 3, 1),  # E^{3,3}_1
        (3, 3, 2),  # E^{3,3}_2
        (3, 3, 3),  # E^{3,3}_3
        (3, 3, 4),  # E^{3,3}_4
        (3, 3, 5),  # E^{3,3}_5
    ]

    # Filter to only coefficients that exist in data
    available = [c for c in selected_coefficients if c in comparison_data]

    if not available:
        print("Warning: No comparison data found for selected coefficients")
        return None

    # Extract data
    labels = []
    tmp_means = []
    tmp_stds = []
    sym_means = []
    sym_stds = []

    for (nA, nB, t) in available:
        entry = comparison_data[(nA, nB, t)]

        if 'TMP' in entry and 'Symbolic' in entry:
            if 'mean' in entry['TMP'] and 'mean' in entry['Symbolic']:
                labels.append(f'$E^{{{nA},{nB}}}_{{{t}}}$')
                tmp_means.append(entry['TMP']['mean'])
                tmp_stds.append(entry['TMP'].get('stddev', 0))
                sym_means.append(entry['Symbolic']['mean'])
                sym_stds.append(entry['Symbolic'].get('stddev', 0))

    if not labels:
        print("Warning: No complete TMP/Symbolic pairs found")
        return None

    # Configure style
    fig_width = COLUMN_WIDTHS[journal]['double']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal=journal, font_scale=font_scale)

    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width, fig_width * 0.5))

    x = np.arange(len(labels))
    width = 0.35

    # Plot bars with error bars
    bars_tmp = ax.bar(x - width/2, tmp_means, width,
                      yerr=tmp_stds, capsize=2,
                      label='TMP', color=COLORS['TMP'],
                      edgecolor='none', alpha=0.85)
    bars_sym = ax.bar(x + width/2, sym_means, width,
                      yerr=sym_stds, capsize=2,
                      label='Symbolic', color=COLORS['Symbolic'],
                      edgecolor='none', alpha=0.85)

    # Formatting
    ax.set_xlabel('Hermite E Coefficient')
    ax.set_ylabel('Execution Time (ns)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.legend(loc='upper left')

    # Add minor gridlines for readability
    ax.yaxis.grid(True, linestyle='--', alpha=0.3, which='major')
    ax.set_axisbelow(True)

    plt.tight_layout()

    # Save
    output_path = output_dir / 'hermite_e_execution_times'
    fig.savefig(f'{output_path}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}.pdf and {output_path}.png")

    return fig


def plot_speedup_analysis(comparison_data, output_dir, journal='nature'):
    """
    Plot 2: Speedup chart showing Symbolic/TMP ratio.

    Identifies crossover point where TMP becomes faster than Symbolic.
    Values > 1 mean Symbolic is faster; < 1 mean TMP is faster.
    """
    # Focus on E^{3,3}_t series to show crossover clearly
    e33_coefficients = [(3, 3, t) for t in range(7)]
    available = [c for c in e33_coefficients if c in comparison_data]

    if not available:
        print("Warning: No E^{3,3}_t data found")
        return None

    t_values = []
    speedup_ratios = []  # TMP/Symbolic (>1 means TMP slower, Symbolic wins)

    for (nA, nB, t) in available:
        entry = comparison_data[(nA, nB, t)]
        if 'TMP' in entry and 'Symbolic' in entry:
            if 'mean' in entry['TMP'] and 'mean' in entry['Symbolic']:
                tmp_time = entry['TMP']['mean']
                sym_time = entry['Symbolic']['mean']
                if sym_time > 0:
                    t_values.append(t)
                    # Speedup = TMP_time / Symbolic_time
                    # >1 means Symbolic is faster
                    speedup_ratios.append(tmp_time / sym_time)

    if not t_values:
        print("Warning: Could not compute speedup ratios")
        return None

    # Configure style
    fig_width = COLUMN_WIDTHS[journal]['single']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal=journal, font_scale=font_scale)

    fig, ax = plt.subplots(figsize=(fig_width, fig_width * 0.85))

    # Plot speedup as bars with color indicating which is faster
    colors = [COLORS['Symbolic'] if s > 1 else COLORS['TMP'] for s in speedup_ratios]
    bars = ax.bar(t_values, speedup_ratios, color=colors, edgecolor='none', alpha=0.85)

    # Add reference line at speedup = 1
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1, label='Equal performance')

    # Add annotations for interpretation
    ax.text(0.02, 0.98, 'Symbolic faster', transform=ax.transAxes,
            fontsize=8, va='top', ha='left', color=COLORS['Symbolic'], fontweight='bold')
    ax.text(0.02, 0.02, 'TMP faster', transform=ax.transAxes,
            fontsize=8, va='bottom', ha='left', color=COLORS['TMP'], fontweight='bold')

    # Formatting
    ax.set_xlabel(r'Index $t$ in $E^{3,3}_t$')
    ax.set_ylabel('Speedup Ratio (TMP time / Symbolic time)')
    ax.set_xticks(t_values)

    # Add value labels on bars
    for bar, ratio in zip(bars, speedup_ratios):
        height = bar.get_height()
        ax.annotate(f'{ratio:.2f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=7)

    plt.tight_layout()

    # Save
    output_path = output_dir / 'hermite_e_speedup_e33'
    fig.savefig(f'{output_path}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}.pdf and {output_path}.png")

    return fig


def plot_scaling_analysis(scaling_data, output_dir, journal='nature'):
    """
    Plot 3: Scaling analysis - execution time vs L_total.

    Shows how both implementations scale with total angular momentum.
    """
    if not scaling_data:
        print("Warning: No scaling data available")
        return None

    L_values = sorted(scaling_data.keys())

    tmp_means = []
    tmp_stds = []
    sym_means = []
    sym_stds = []
    valid_L = []

    for L in L_values:
        entry = scaling_data[L]
        if 'TMP' in entry and 'Symbolic' in entry:
            if 'mean' in entry['TMP'] and 'mean' in entry['Symbolic']:
                valid_L.append(L)
                tmp_means.append(entry['TMP']['mean'])
                tmp_stds.append(entry['TMP'].get('stddev', 0))
                sym_means.append(entry['Symbolic']['mean'])
                sym_stds.append(entry['Symbolic'].get('stddev', 0))

    if not valid_L:
        print("Warning: No complete scaling data found")
        return None

    # Configure style
    fig_width = COLUMN_WIDTHS[journal]['single']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal=journal, font_scale=font_scale)

    fig, ax = plt.subplots(figsize=(fig_width, fig_width * 0.85))

    # Plot with error bars and markers
    ax.errorbar(valid_L, tmp_means, yerr=tmp_stds,
                color=COLORS['TMP'], marker='o', markersize=6,
                capsize=3, capthick=1, linewidth=1.5,
                label='TMP', linestyle='-')
    ax.errorbar(valid_L, sym_means, yerr=sym_stds,
                color=COLORS['Symbolic'], marker='s', markersize=6,
                capsize=3, capthick=1, linewidth=1.5,
                label='Symbolic', linestyle='--')

    # Formatting
    ax.set_xlabel(r'Total Angular Momentum $L_{\mathrm{total}}$')
    ax.set_ylabel('Execution Time (ns)')
    ax.set_xticks(valid_L)
    ax.legend(loc='upper left')

    # Add minor gridlines
    ax.yaxis.grid(True, linestyle='--', alpha=0.3, which='major')
    ax.set_axisbelow(True)

    plt.tight_layout()

    # Save
    output_path = output_dir / 'hermite_e_scaling'
    fig.savefig(f'{output_path}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}.pdf and {output_path}.png")

    return fig


def plot_comprehensive_comparison(comparison_data, scaling_data, output_dir, journal='nature'):
    """
    Combined multi-panel figure for publication.

    Panel (a): Bar chart of execution times for key coefficients
    Panel (b): Speedup analysis for E^{3,3}_t series
    Panel (c): Scaling with L_total
    """
    # Configure style
    fig_width = COLUMN_WIDTHS[journal]['double']
    font_scale = get_font_scale(fig_width)
    configure_publication_style(journal=journal, font_scale=font_scale)

    # Create 1x3 subplot layout
    fig, axes = plt.subplots(1, 3, figsize=(fig_width, fig_width * 0.35))

    # ===== Panel (a): Execution times for E^{3,3}_t =====
    ax = axes[0]
    e33_coefficients = [(3, 3, t) for t in range(6)]  # t=0 to 5

    labels = []
    tmp_means = []
    tmp_stds = []
    sym_means = []
    sym_stds = []

    for (nA, nB, t) in e33_coefficients:
        if (nA, nB, t) in comparison_data:
            entry = comparison_data[(nA, nB, t)]
            if 'TMP' in entry and 'Symbolic' in entry:
                if 'mean' in entry['TMP'] and 'mean' in entry['Symbolic']:
                    labels.append(f'$t={t}$')
                    tmp_means.append(entry['TMP']['mean'])
                    tmp_stds.append(entry['TMP'].get('stddev', 0))
                    sym_means.append(entry['Symbolic']['mean'])
                    sym_stds.append(entry['Symbolic'].get('stddev', 0))

    if labels:
        x = np.arange(len(labels))
        width = 0.35
        ax.bar(x - width/2, tmp_means, width, yerr=tmp_stds, capsize=2,
               label='TMP', color=COLORS['TMP'], edgecolor='none', alpha=0.85)
        ax.bar(x + width/2, sym_means, width, yerr=sym_stds, capsize=2,
               label='Symbolic', color=COLORS['Symbolic'], edgecolor='none', alpha=0.85)
        ax.set_xlabel(r'Index $t$ in $E^{3,3}_t$')
        ax.set_ylabel('Time (ns)')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(loc='upper left', fontsize=7)
        ax.yaxis.grid(True, linestyle='--', alpha=0.3)
        ax.set_axisbelow(True)

    add_panel_label(ax, '(a)', fontsize=10)

    # ===== Panel (b): Speedup ratios =====
    ax = axes[1]
    t_values = []
    speedup_ratios = []

    for t in range(6):
        if (3, 3, t) in comparison_data:
            entry = comparison_data[(3, 3, t)]
            if 'TMP' in entry and 'Symbolic' in entry:
                if 'mean' in entry['TMP'] and 'mean' in entry['Symbolic']:
                    tmp_time = entry['TMP']['mean']
                    sym_time = entry['Symbolic']['mean']
                    if sym_time > 0:
                        t_values.append(t)
                        speedup_ratios.append(tmp_time / sym_time)

    if t_values:
        colors = [COLORS['Symbolic'] if s > 1 else COLORS['TMP'] for s in speedup_ratios]
        bars = ax.bar(t_values, speedup_ratios, color=colors, edgecolor='none', alpha=0.85)
        ax.axhline(y=1, color='gray', linestyle='--', linewidth=1)
        ax.set_xlabel(r'Index $t$')
        ax.set_ylabel('Speedup (TMP/Symbolic)')
        ax.set_xticks(t_values)

        # Compact annotations
        ax.text(0.95, 0.95, 'Sym. faster', transform=ax.transAxes,
                fontsize=6, va='top', ha='right', color=COLORS['Symbolic'])
        ax.text(0.95, 0.05, 'TMP faster', transform=ax.transAxes,
                fontsize=6, va='bottom', ha='right', color=COLORS['TMP'])

    add_panel_label(ax, '(b)', fontsize=10)

    # ===== Panel (c): Scaling analysis =====
    ax = axes[2]

    if scaling_data:
        L_values = sorted(scaling_data.keys())
        tmp_means = []
        tmp_stds = []
        sym_means = []
        sym_stds = []
        valid_L = []

        for L in L_values:
            entry = scaling_data[L]
            if 'TMP' in entry and 'Symbolic' in entry:
                if 'mean' in entry['TMP'] and 'mean' in entry['Symbolic']:
                    valid_L.append(L)
                    tmp_means.append(entry['TMP']['mean'])
                    tmp_stds.append(entry['TMP'].get('stddev', 0))
                    sym_means.append(entry['Symbolic']['mean'])
                    sym_stds.append(entry['Symbolic'].get('stddev', 0))

        if valid_L:
            ax.errorbar(valid_L, tmp_means, yerr=tmp_stds,
                        color=COLORS['TMP'], marker='o', markersize=4,
                        capsize=2, linewidth=1.2, label='TMP')
            ax.errorbar(valid_L, sym_means, yerr=sym_stds,
                        color=COLORS['Symbolic'], marker='s', markersize=4,
                        capsize=2, linewidth=1.2, linestyle='--', label='Symbolic')
            ax.set_xlabel(r'$L_{\mathrm{total}}$')
            ax.set_ylabel('Time (ns)')
            ax.set_xticks(valid_L)
            ax.legend(loc='upper left', fontsize=7)
            ax.yaxis.grid(True, linestyle='--', alpha=0.3)
            ax.set_axisbelow(True)

    add_panel_label(ax, '(c)', fontsize=10)

    plt.tight_layout()

    # Save
    output_path = output_dir / 'hermite_e_comprehensive'
    fig.savefig(f'{output_path}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}.pdf and {output_path}.png")

    return fig


def print_data_summary(comparison_data, scaling_data):
    """Print summary of extracted benchmark data."""
    print("\n" + "="*70)
    print("BENCHMARK DATA SUMMARY")
    print("="*70)

    print("\n--- Comparison Benchmarks (TMP vs Symbolic) ---")
    print(f"Total coefficient combinations: {len(comparison_data)}")

    # Print E^{3,3}_t series details
    print("\nE^{3,3}_t series (key comparison):")
    print(f"{'Coeff':<12} {'TMP (ns)':<15} {'Symbolic (ns)':<15} {'Speedup':<10}")
    print("-" * 52)

    for t in range(7):
        key = (3, 3, t)
        if key in comparison_data:
            entry = comparison_data[key]
            if 'TMP' in entry and 'Symbolic' in entry:
                tmp = entry['TMP'].get('mean', 0)
                sym = entry['Symbolic'].get('mean', 0)
                speedup = tmp / sym if sym > 0 else 0
                winner = "Sym" if speedup > 1 else "TMP"
                print(f"E^{{3,3}}_{t:<5} {tmp:<15.3f} {sym:<15.3f} {speedup:<6.2f} ({winner})")

    print("\n--- Scaling Benchmarks ---")
    if scaling_data:
        print(f"L_total values: {sorted(scaling_data.keys())}")
        print("\nScaling data:")
        print(f"{'L_total':<10} {'TMP (ns)':<15} {'Symbolic (ns)':<15}")
        print("-" * 40)
        for L in sorted(scaling_data.keys()):
            entry = scaling_data[L]
            tmp = entry.get('TMP', {}).get('mean', 0)
            sym = entry.get('Symbolic', {}).get('mean', 0)
            print(f"{L:<10} {tmp:<15.3f} {sym:<15.3f}")
    else:
        print("No scaling data found")

    print("="*70 + "\n")


# ==============================================================================
# Main Execution
# ==============================================================================

def main():
    # Paths
    json_path = Path('/home/ruben/Research/Science/Projects/RECURSUM/benchmarks/results/raw/hermite_e_comparison.json')
    output_dir = Path('/home/ruben/Research/Science/Projects/RECURSUM/benchmarks/results/figures')

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading benchmark data...")
    data = load_benchmark_data(json_path)

    print("Extracting comparison benchmarks...")
    comparison_data = extract_comparison_benchmarks(data)

    print("Extracting scaling benchmarks...")
    scaling_data = extract_scaling_benchmarks(data)

    # Print data summary
    print_data_summary(comparison_data, scaling_data)

    # Generate plots
    print("Generating plots...")

    print("\n1. Execution time comparison bar chart...")
    plot_execution_time_comparison(comparison_data, output_dir)

    print("\n2. Speedup analysis for E^{3,3}_t...")
    plot_speedup_analysis(comparison_data, output_dir)

    print("\n3. Scaling analysis...")
    plot_scaling_analysis(scaling_data, output_dir)

    print("\n4. Comprehensive multi-panel figure...")
    plot_comprehensive_comparison(comparison_data, scaling_data, output_dir)

    print("\nAll plots generated successfully!")
    print(f"Output directory: {output_dir}")


if __name__ == '__main__':
    main()
