#!/usr/bin/env python3
"""
Generate publication-ready plots for McMurchie-Davidson benchmark results.
Compares TMP, Layered, LayeredCodegen, and Symbolic implementations.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from typing import Dict, List, Tuple

# ============================================================================
# Configuration from scientific_plotting_agent_guide.md
# ============================================================================

# Okabe-Ito colorblind-safe palette
OKABE_ITO_CYCLE = [
    '#0072B2',  # blue
    '#E69F00',  # orange
    '#009E73',  # green
    '#CC79A7',  # purple/pink
    '#D55E00',  # vermillion
    '#56B4E9',  # sky blue
    '#F0E442',  # yellow
    '#000000',  # black
]

# Line styles for additional differentiation
LINE_STYLES = ['-', '--', '-.', ':']
MARKERS = ['o', 's', '^', 'D', 'v', '<', '>', 'p']

# Journal column widths (Nature)
COLUMN_WIDTHS = {
    'nature': {'single': 3.46, 'onehalf': 4.76, 'double': 7.09},
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
        # Font configuration
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',

        # Font sizes
        **scaled_sizes,

        # Axes configuration
        'axes.linewidth': 0.8,
        'axes.labelpad': 4,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'grid.linewidth': 0.5,

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

        # Color cycle
        'axes.prop_cycle': plt.cycler(color=OKABE_ITO_CYCLE),
    })

def get_figsize(journal='nature', width='single', aspect=0.75):
    """Get figure dimensions for target journal."""
    w = COLUMN_WIDTHS[journal][width]
    h = w * aspect
    return (w, h)

# ============================================================================
# Data Loading and Processing
# ============================================================================

def load_benchmark_data(json_path: Path) -> Dict:
    """Load Google Benchmark JSON output."""
    with open(json_path, 'r') as f:
        return json.load(f)

def extract_hermite_statistics(data: Dict) -> Dict[str, Dict]:
    """
    Extract mean and stddev statistics for Hermite E coefficients.

    Returns dictionary keyed by (impl, shell_pair) with mean/stddev values.
    """
    results = {}

    for bench in data['benchmarks']:
        name = bench['name']

        # Only process mean statistics
        if '_mean' not in name:
            continue

        # Extract implementation and shell pair from name
        # Format: HermiteE/{impl}/{shell}/min_time:1.000_mean
        parts = name.split('/')
        impl_name = parts[1]  # TMP, Layered, LayeredCodegen, Symbolic
        shell = parts[2]  # ss, sp, pp, sd, pd, dd, ff, gg

        impl = int(bench['impl'])
        L = int(bench['L'])
        nA = int(bench['nA'])
        nB = int(bench['nB'])
        n_coeffs = int(bench['n_coeffs'])

        key = (impl, L, shell)

        # Find corresponding stddev
        stddev = 0.0
        stddev_name = name.replace('_mean', '_stddev')
        for b in data['benchmarks']:
            if b['name'] == stddev_name:
                stddev = b['real_time']
                break

        results[key] = {
            'impl': impl,
            'impl_name': impl_name,
            'shell': shell,
            'L': L,
            'nA': nA,
            'nB': nB,
            'n_coeffs': n_coeffs,
            'mean': bench['real_time'],
            'stddev': stddev,
            'time_unit': bench['time_unit'],
        }

    return results

def extract_coulomb_statistics(data: Dict) -> Dict[str, Dict]:
    """
    Extract mean and stddev statistics for Coulomb R integrals.

    Returns dictionary keyed by (impl, L_total) with mean/stddev values.
    """
    results = {}

    for bench in data['benchmarks']:
        name = bench['name']

        # Only process mean statistics
        if '_mean' not in name:
            continue

        # Extract implementation from name
        # Format: CoulombR/{impl}/L{L}/min_time:1.000_mean
        parts = name.split('/')
        impl_name = parts[1]  # TMP, Layered

        impl = int(bench['impl'])
        L_total = int(bench['L_total'])
        n_integrals = int(bench['n_integrals'])

        key = (impl, L_total)

        # Find corresponding stddev
        stddev = 0.0
        stddev_name = name.replace('_mean', '_stddev')
        for b in data['benchmarks']:
            if b['name'] == stddev_name:
                stddev = b['real_time']
                break

        results[key] = {
            'impl': impl,
            'impl_name': impl_name,
            'L_total': L_total,
            'n_integrals': n_integrals,
            'mean': bench['real_time'],
            'stddev': stddev,
            'time_unit': bench['time_unit'],
        }

    return results

# ============================================================================
# Plotting Functions
# ============================================================================

def plot_hermite_comparison(stats: Dict, output_dir: Path):
    """
    Plot 1: Hermite Coefficients Comparison (4 implementations).
    """
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=get_figsize('nature', 'onehalf', aspect=0.85))

    # Define shell pairs in order
    shell_order = ['ss', 'sp', 'pp', 'sd', 'pd', 'dd', 'ff', 'gg']

    # Implementation names and styling
    impl_info = {
        0: {'name': 'TMP (Original)', 'color': OKABE_ITO_CYCLE[0], 'marker': 'o', 'ls': '-'},
        1: {'name': 'Layered (Hand-written)', 'color': OKABE_ITO_CYCLE[1], 'marker': 's', 'ls': '--'},
        2: {'name': 'Symbolic (SymPy)', 'color': OKABE_ITO_CYCLE[2], 'marker': '^', 'ls': '-.'},
        3: {'name': 'LayeredCodegen (NEW)', 'color': OKABE_ITO_CYCLE[4], 'marker': 'D', 'ls': '-'},
    }

    # Organize data by implementation
    for impl in [0, 1, 2, 3]:
        x_positions = []
        y_means = []
        y_errors = []
        x_labels = []

        for idx, shell in enumerate(shell_order):
            # Find data for this impl and shell
            for key, data in stats.items():
                if data['impl'] == impl and data['shell'] == shell:
                    x_positions.append(idx)
                    y_means.append(data['mean'])
                    y_errors.append(data['stddev'])
                    if impl == 0:  # Only collect labels once
                        x_labels.append(f"{shell}\nL={data['L']}")
                    break

        if y_means:  # Only plot if data exists
            info = impl_info[impl]
            ax.errorbar(x_positions, y_means, yerr=y_errors,
                       label=info['name'], color=info['color'],
                       marker=info['marker'], linestyle=info['ls'],
                       markersize=6, capsize=3, linewidth=1.5,
                       markerfacecolor='white', markeredgewidth=1.5)

    # Formatting
    ax.set_yscale('log')
    ax.set_xticks(range(len(shell_order)))
    ax.set_xticklabels([f"{s}\nL={i//2}" for i, s in enumerate(shell_order)], fontsize=8)
    ax.set_xlabel('Shell Pair and Total Angular Momentum L', fontweight='bold')
    ax.set_ylabel('Time (nanoseconds)', fontweight='bold')
    ax.set_title('Hermite Expansion Coefficients:\nTMP vs Layered vs LayeredCodegen vs Symbolic',
                fontsize=11, fontweight='bold', pad=10)
    ax.legend(loc='upper left', fontsize=7, frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.3, which='both')

    # Save
    for fmt in ['png', 'pdf']:
        output_path = output_dir / f'hermite_coefficients_comparison.{fmt}'
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_path}")

    plt.close(fig)

def plot_hermite_scaling(stats: Dict, output_dir: Path):
    """
    Plot 2: Hermite Coefficients Scaling with L.
    """
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=get_figsize('nature', 'single', aspect=0.85))

    impl_info = {
        0: {'name': 'TMP', 'color': OKABE_ITO_CYCLE[0], 'marker': 'o', 'ls': '-'},
        1: {'name': 'Layered', 'color': OKABE_ITO_CYCLE[1], 'marker': 's', 'ls': '--'},
        2: {'name': 'Symbolic', 'color': OKABE_ITO_CYCLE[2], 'marker': '^', 'ls': '-.'},
        3: {'name': 'LayeredCodegen', 'color': OKABE_ITO_CYCLE[4], 'marker': 'D', 'ls': '-'},
    }

    # Organize data by implementation and L
    for impl in [0, 1, 2, 3]:
        L_values = []
        means = []
        errors = []

        # Collect data grouped by L
        L_dict = {}
        for key, data in stats.items():
            if data['impl'] == impl:
                L = data['L']
                if L not in L_dict:
                    L_dict[L] = []
                L_dict[L].append((data['mean'], data['stddev']))

        # Average over shells with same L if multiple exist
        for L in sorted(L_dict.keys()):
            L_values.append(L)
            values = L_dict[L]
            mean_of_means = np.mean([v[0] for v in values])
            # Propagate errors (quadrature sum divided by sqrt(n))
            avg_error = np.sqrt(np.sum([v[1]**2 for v in values])) / len(values)
            means.append(mean_of_means)
            errors.append(avg_error)

        if means:
            info = impl_info[impl]
            ax.errorbar(L_values, means, yerr=errors,
                       label=info['name'], color=info['color'],
                       marker=info['marker'], linestyle=info['ls'],
                       markersize=6, capsize=3, linewidth=1.5,
                       markerfacecolor='white', markeredgewidth=1.5)

    # Formatting
    ax.set_yscale('log')
    ax.set_xlabel('Total Angular Momentum L', fontweight='bold')
    ax.set_ylabel('Time (nanoseconds)', fontweight='bold')
    ax.set_title('Hermite Coefficient Performance Scaling', fontsize=11, fontweight='bold')
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3, which='both')

    # Save
    for fmt in ['png', 'pdf']:
        output_path = output_dir / f'hermite_coefficients_vs_L.{fmt}'
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_path}")

    plt.close(fig)

def plot_layered_codegen_speedup(stats: Dict, output_dir: Path):
    """
    Plot 3: LayeredCodegen Speedup over Hand-written Layered.
    """
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=get_figsize('nature', 'single', aspect=0.75))

    # Define shell pairs in order
    shell_order = ['ss', 'sp', 'pp', 'sd', 'pd', 'dd', 'ff', 'gg']

    x_positions = []
    speedups = []
    x_labels = []

    for idx, shell in enumerate(shell_order):
        layered_time = None
        codegen_time = None
        L_val = None

        # Find data for Layered (impl=1) and LayeredCodegen (impl=3)
        for key, data in stats.items():
            if data['shell'] == shell:
                if data['impl'] == 1:  # Layered
                    layered_time = data['mean']
                    L_val = data['L']
                elif data['impl'] == 3:  # LayeredCodegen
                    codegen_time = data['mean']

        if layered_time and codegen_time:
            speedup = layered_time / codegen_time
            x_positions.append(idx)
            speedups.append(speedup)
            x_labels.append(f"{shell}\nL={L_val}")

    # Plot bars
    bars = ax.bar(x_positions, speedups, color=OKABE_ITO_CYCLE[3],
                  edgecolor='black', linewidth=0.8, alpha=0.8)

    # Add reference line at 1.0 (no speedup)
    ax.axhline(y=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5, label='No speedup')

    # Add value labels on bars
    for i, (pos, speedup) in enumerate(zip(x_positions, speedups)):
        ax.text(pos, speedup + 0.2, f'{speedup:.1f}×',
               ha='center', va='bottom', fontsize=8, fontweight='bold')

    # Formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=8)
    ax.set_xlabel('Shell Pair and Total Angular Momentum L', fontweight='bold')
    ax.set_ylabel('Speedup Factor (Layered / LayeredCodegen)', fontweight='bold')
    ax.set_title('LayeredCodegen Speedup vs Hand-written Layered',
                fontsize=11, fontweight='bold', pad=10)
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(bottom=0)

    # Save
    for fmt in ['png', 'pdf']:
        output_path = output_dir / f'hermite_layered_codegen_speedup.{fmt}'
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_path}")

    plt.close(fig)

def plot_coulomb_comparison(stats: Dict, output_dir: Path):
    """
    Plot 4: Coulomb R Integrals Comparison.
    """
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=get_figsize('nature', 'single', aspect=0.85))

    impl_info = {
        0: {'name': 'TMP', 'color': OKABE_ITO_CYCLE[0], 'marker': 'o', 'ls': '-'},
        1: {'name': 'Layered', 'color': OKABE_ITO_CYCLE[1], 'marker': 's', 'ls': '--'},
    }

    for impl in [0, 1]:
        L_values = []
        means = []
        errors = []

        for key, data in stats.items():
            if data['impl'] == impl:
                L_values.append(data['L_total'])
                means.append(data['mean'])
                errors.append(data['stddev'])

        # Sort by L_total
        if L_values:
            sorted_data = sorted(zip(L_values, means, errors))
            L_values, means, errors = zip(*sorted_data)

            info = impl_info[impl]
            ax.errorbar(L_values, means, yerr=errors,
                       label=info['name'], color=info['color'],
                       marker=info['marker'], linestyle=info['ls'],
                       markersize=7, capsize=3, linewidth=1.5,
                       markerfacecolor='white', markeredgewidth=1.5)

    # Formatting
    ax.set_yscale('log')
    ax.set_xlabel('Total Angular Momentum L', fontweight='bold')
    ax.set_ylabel('Time (nanoseconds)', fontweight='bold')
    ax.set_title('Coulomb Auxiliary Integrals $R_{tuv}$:\nTMP vs Layered',
                fontsize=11, fontweight='bold', pad=10)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xticks(range(0, 9))

    # Save
    for fmt in ['png', 'pdf']:
        output_path = output_dir / f'coulomb_hermite_comparison.{fmt}'
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_path}")

    plt.close(fig)

def plot_coulomb_scaling(stats: Dict, output_dir: Path):
    """
    Plot 5: Coulomb R Scaling Analysis.
    """
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=get_figsize('nature', 'single', aspect=0.85))

    impl_info = {
        0: {'name': 'TMP', 'color': OKABE_ITO_CYCLE[0], 'marker': 'o', 'ls': '-'},
        1: {'name': 'Layered', 'color': OKABE_ITO_CYCLE[1], 'marker': 's', 'ls': '--'},
    }

    for impl in [0, 1]:
        n_integrals = []
        means = []
        errors = []

        for key, data in stats.items():
            if data['impl'] == impl:
                n_integrals.append(data['n_integrals'])
                means.append(data['mean'])
                errors.append(data['stddev'])

        # Sort by n_integrals
        if n_integrals:
            sorted_data = sorted(zip(n_integrals, means, errors))
            n_integrals, means, errors = zip(*sorted_data)

            info = impl_info[impl]
            ax.errorbar(n_integrals, means, yerr=errors,
                       label=info['name'], color=info['color'],
                       marker=info['marker'], linestyle=info['ls'],
                       markersize=7, capsize=3, linewidth=1.5,
                       markerfacecolor='white', markeredgewidth=1.5)

    # Formatting
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Number of Integrals', fontweight='bold')
    ax.set_ylabel('Time (nanoseconds)', fontweight='bold')
    ax.set_title('Coulomb Auxiliary Integrals:\nScaling with Problem Size',
                fontsize=11, fontweight='bold', pad=10)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

    # Save
    for fmt in ['png', 'pdf']:
        output_path = output_dir / f'coulomb_hermite_scaling.{fmt}'
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_path}")

    plt.close(fig)

# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Generate all benchmark plots."""
    # Setup paths
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / 'results' / 'raw'
    output_dir = base_dir / 'results' / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print("Loading benchmark data...")
    hermite_file = data_dir / 'hermite_coefficients.json'
    coulomb_file = data_dir / 'coulomb_hermite.json'

    hermite_data = load_benchmark_data(hermite_file)
    coulomb_data = load_benchmark_data(coulomb_file)

    # Extract statistics
    print("Processing Hermite coefficient benchmarks...")
    hermite_stats = extract_hermite_statistics(hermite_data)
    print(f"  Found {len(hermite_stats)} benchmark results")

    print("Processing Coulomb R integral benchmarks...")
    coulomb_stats = extract_coulomb_statistics(coulomb_data)
    print(f"  Found {len(coulomb_stats)} benchmark results")

    # Generate plots
    print("\nGenerating plots...")
    print("  [1/5] Hermite coefficients comparison...")
    plot_hermite_comparison(hermite_stats, output_dir)

    print("  [2/5] Hermite coefficients scaling with L...")
    plot_hermite_scaling(hermite_stats, output_dir)

    print("  [3/5] LayeredCodegen speedup analysis...")
    plot_layered_codegen_speedup(hermite_stats, output_dir)

    print("  [4/5] Coulomb R integrals comparison...")
    plot_coulomb_comparison(coulomb_stats, output_dir)

    print("  [5/5] Coulomb R scaling analysis...")
    plot_coulomb_scaling(coulomb_stats, output_dir)

    print(f"\nAll plots saved to: {output_dir}")
    print("\nKey findings:")

    # Calculate and report key statistics
    # Find LayeredCodegen vs Layered speedup
    layered_ss = hermite_stats.get((1, 0, 'ss'), {}).get('mean', None)
    codegen_ss = hermite_stats.get((3, 0, 'ss'), {}).get('mean', None)
    if layered_ss and codegen_ss:
        speedup = layered_ss / codegen_ss
        print(f"  - LayeredCodegen achieves {speedup:.1f}× speedup over hand-written Layered (ss shell)")

    # Compare TMP vs LayeredCodegen
    tmp_ss = hermite_stats.get((0, 0, 'ss'), {}).get('mean', None)
    if tmp_ss and codegen_ss:
        ratio = codegen_ss / tmp_ss
        print(f"  - LayeredCodegen is {ratio:.2f}× the time of TMP (ss shell)")
        if ratio < 1.1:
            print(f"    → LayeredCodegen MATCHES TMP performance!")

    print("\nDone!")

if __name__ == '__main__':
    main()
