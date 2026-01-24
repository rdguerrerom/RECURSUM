#!/usr/bin/env python3
"""
Publication-quality plots for McMurchie-Davidson benchmark results.

Compares three implementations for computing Hermite E coefficients:
- impl=0: Original TMP (separate recursive calls per coefficient)
- impl=1: Symbolic (SymPy-generated polynomial expressions)
- impl=2: Layered TMP with CSE (computes ALL t values at once, true CSE)

Following scientific plotting guide standards:
- Colorblind-safe Okabe-Ito palette
- Journal-appropriate figure dimensions (APS double column)
- Proper font scaling for publication
- Error bars showing uncertainty
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from collections import defaultdict

# =============================================================================
# Configuration following scientific_plotting_agent_guide.md
# =============================================================================

# Okabe-Ito colorblind-safe palette
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

# Implementation colors (3 distinct, colorblind-safe)
IMPL_COLORS = {
    0: OKABE_ITO['blue'],       # Original TMP
    1: OKABE_ITO['orange'],     # Symbolic
    2: OKABE_ITO['green'],      # Layered CSE
}

IMPL_NAMES = {
    0: 'Original TMP',
    1: 'Symbolic',
    2: 'Layered CSE',
}

# Journal column widths (APS Physical Review)
COLUMN_WIDTHS = {
    'single': 3.39,
    'onehalf': 5.12,
    'double': 7.01,
}

# Shell pair order by L value
SHELL_PAIR_ORDER = ['ss', 'sp', 'pp', 'sd', 'pd', 'sf', 'pf', 'dd', 'df', 'sg', 'pg', 'ff', 'dg', 'fg', 'gg']
SHELL_PAIR_L = {
    'ss': 0, 'sp': 1, 'pp': 2, 'sd': 2, 'pd': 3, 'sf': 3, 'pf': 4,
    'dd': 4, 'df': 5, 'sg': 4, 'pg': 5, 'ff': 6, 'dg': 6, 'fg': 7, 'gg': 8
}


def configure_publication_style(font_scale=1.2):
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
        'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',

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
        'lines.markersize': 5,

        # Figure configuration
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
    })


def load_benchmark_data(json_path):
    """Load and parse benchmark JSON data."""
    with open(json_path, 'r') as f:
        data = json.load(f)

    benchmarks = data['benchmarks']

    # Organize data by shell pair and implementation
    results = defaultdict(lambda: defaultdict(dict))

    for b in benchmarks:
        if b['run_type'] != 'aggregate':
            continue

        name = b['name']
        agg_name = b.get('aggregate_name', '')

        # Extract shell pair and impl from name
        shell_pair = None
        impl = None

        if 'Orig_' in name:
            shell_pair = name.split('Orig_')[1].split('/')[0]
            impl = 0
        elif 'CSE_' in name:
            shell_pair = name.split('CSE_')[1].split('/')[0]
            impl = 2
        elif 'Sym_' in name:
            shell_pair = name.split('Sym_')[1].split('/')[0]
            impl = 1

        if shell_pair is None or impl is None:
            continue

        # Handle aggregate name suffixes
        if '_mean' in name:
            results[shell_pair][impl]['mean'] = b['cpu_time']
            results[shell_pair][impl]['L'] = b.get('L', SHELL_PAIR_L.get(shell_pair, 0))
            results[shell_pair][impl]['n_coeffs'] = b.get('n_coeffs', 1)
        elif '_stddev' in name:
            results[shell_pair][impl]['stddev'] = b['cpu_time']
        elif '_median' in name:
            results[shell_pair][impl]['median'] = b['cpu_time']

    return results


def plot_shell_pair_comparison(results, output_path):
    """
    Plot 1: Bar chart comparison of implementations by shell pair.
    """
    configure_publication_style(font_scale=1.2)

    fig, ax = plt.subplots(figsize=(COLUMN_WIDTHS['double'], 4.0))

    # Filter to shell pairs that have all implementations
    shell_pairs = [sp for sp in SHELL_PAIR_ORDER if sp in results]
    n_pairs = len(shell_pairs)
    n_impls = 3

    x = np.arange(n_pairs)
    width = 0.25

    # Plot bars for each implementation
    for i, impl in enumerate([0, 2, 1]):  # Order: Original, CSE, Symbolic
        means = []
        stds = []
        for sp in shell_pairs:
            if impl in results[sp]:
                means.append(results[sp][impl].get('mean', 0))
                stds.append(results[sp][impl].get('stddev', 0))
            else:
                means.append(0)
                stds.append(0)

        offset = (i - 1) * width
        bars = ax.bar(x + offset, means, width, yerr=stds,
                     label=IMPL_NAMES[impl], color=IMPL_COLORS[impl],
                     capsize=2, error_kw={'linewidth': 0.8})

    ax.set_xlabel('Shell Pair')
    ax.set_ylabel('Time per Layer (ns)')
    ax.set_xticks(x)
    ax.set_xticklabels(shell_pairs, fontsize=8)
    ax.legend(loc='upper left')

    # Add secondary x-axis showing L values
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(x)
    L_labels = [str(SHELL_PAIR_L.get(sp, '')) for sp in shell_pairs]
    ax2.set_xticklabels(L_labels, fontsize=7)
    ax2.set_xlabel('Total Angular Momentum L', fontsize=9)
    ax2.spines['top'].set_visible(False)
    ax2.tick_params(length=0)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    # Also save PNG for preview
    png_path = str(output_path).replace('.pdf', '.png')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_path}")
    print(f"Saved: {png_path}")


def plot_speedup_vs_L(results, output_path):
    """
    Plot 2: Speedup analysis showing CSE speedup vs Original and Symbolic.
    """
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=(COLUMN_WIDTHS['single'], 3.0))

    # Calculate speedups
    shell_pairs = [sp for sp in SHELL_PAIR_ORDER if sp in results]

    L_values = []
    speedup_vs_orig = []
    speedup_vs_symb = []
    shell_labels = []

    for sp in shell_pairs:
        if 0 in results[sp] and 2 in results[sp]:
            L = SHELL_PAIR_L.get(sp, 0)
            L_values.append(L)
            shell_labels.append(sp)

            orig_time = results[sp][0].get('mean', 1)
            cse_time = results[sp][2].get('mean', 1)
            symb_time = results[sp].get(1, {}).get('mean', orig_time)

            speedup_vs_orig.append(orig_time / cse_time if cse_time > 0 else 1)
            speedup_vs_symb.append(symb_time / cse_time if cse_time > 0 else 1)

    # Plot speedup lines
    ax.plot(L_values, speedup_vs_orig, 'o-', color=OKABE_ITO['blue'],
            label='CSE vs Original', markersize=6, markerfacecolor='white', markeredgewidth=1.5)
    ax.plot(L_values, speedup_vs_symb, 's--', color=OKABE_ITO['orange'],
            label='CSE vs Symbolic', markersize=5, markerfacecolor='white', markeredgewidth=1.5)

    # Reference line at speedup = 1
    ax.axhline(y=1, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax.text(ax.get_xlim()[1] * 0.98, 1.05, 'No speedup', fontsize=7,
            ha='right', va='bottom', color='gray')

    ax.set_xlabel('Total Angular Momentum L')
    ax.set_ylabel('Speedup Factor')
    ax.legend(loc='upper left', fontsize=8)

    ax.set_xlim(-0.3, max(L_values) + 0.3)
    ax.set_xticks(sorted(set(L_values)))

    # Shade region where CSE wins
    ax.fill_between(ax.get_xlim(), 1, ax.get_ylim()[1], alpha=0.1, color=OKABE_ITO['green'])

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    # Also save PNG for preview
    png_path = str(output_path).replace('.pdf', '.png')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_path}")
    print(f"Saved: {png_path}")


def plot_timing_vs_L(results, output_path):
    """
    Plot 3: Scaling with L showing timing for all implementations on log scale.
    """
    configure_publication_style(font_scale=1.3)

    fig, ax = plt.subplots(figsize=(COLUMN_WIDTHS['single'], 3.5))

    shell_pairs = [sp for sp in SHELL_PAIR_ORDER if sp in results]

    # Collect data for each implementation
    impl_data = {0: {'L': [], 'time': [], 'sp': []},
                 1: {'L': [], 'time': [], 'sp': []},
                 2: {'L': [], 'time': [], 'sp': []}}

    for sp in shell_pairs:
        L = SHELL_PAIR_L.get(sp, 0)
        for impl in [0, 1, 2]:
            if impl in results[sp] and 'mean' in results[sp][impl]:
                impl_data[impl]['L'].append(L)
                impl_data[impl]['time'].append(results[sp][impl]['mean'])
                impl_data[impl]['sp'].append(sp)

    # Markers for each implementation
    markers = {0: 'o', 1: 's', 2: '^'}

    # Plot each implementation
    for impl in [0, 1, 2]:
        if impl_data[impl]['L']:
            ax.semilogy(impl_data[impl]['L'], impl_data[impl]['time'],
                       marker=markers[impl], linestyle='-' if impl != 1 else '--',
                       color=IMPL_COLORS[impl], label=IMPL_NAMES[impl],
                       markersize=5, markerfacecolor='white', markeredgewidth=1.2)

    # Add theoretical scaling reference lines
    L_ref = np.array([0, 2, 4, 6, 8])

    # O(3^L) scaling - exponential
    base_exp = impl_data[0]['time'][0] if impl_data[0]['time'] else 1
    exp_scaling = base_exp * (3 ** L_ref) / (3 ** L_ref[0])
    ax.semilogy(L_ref, exp_scaling, ':', color='gray', alpha=0.5, linewidth=1.5,
               label=r'$O(3^L)$ scaling')

    ax.set_xlabel('Total Angular Momentum L')
    ax.set_ylabel('Time per Layer (ns)')
    ax.legend(loc='upper left', fontsize=7)

    ax.set_xlim(-0.3, 8.3)
    ax.set_xticks(range(9))

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    # Also save PNG for preview
    png_path = str(output_path).replace('.pdf', '.png')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_path}")
    print(f"Saved: {png_path}")


def generate_summary_table(results, output_path):
    """
    Generate summary table showing winner for each shell pair.
    """
    shell_pairs = [sp for sp in SHELL_PAIR_ORDER if sp in results]

    lines = []
    lines.append("=" * 90)
    lines.append("McMurchie-Davidson Benchmark Summary: Winner by Shell Pair")
    lines.append("=" * 90)
    lines.append("")
    lines.append(f"{'Shell':<8} {'L':<4} {'Original (ns)':<16} {'Symbolic (ns)':<16} {'CSE (ns)':<16} {'Winner':<12} {'Improvement':<12}")
    lines.append("-" * 90)

    total_wins = {0: 0, 1: 0, 2: 0}

    for sp in shell_pairs:
        L = SHELL_PAIR_L.get(sp, 0)

        times = {}
        for impl in [0, 1, 2]:
            if impl in results[sp] and 'mean' in results[sp][impl]:
                times[impl] = results[sp][impl]['mean']

        if not times:
            continue

        # Find winner (minimum time)
        winner = min(times, key=times.get)
        total_wins[winner] += 1

        # Calculate improvement vs slowest
        slowest = max(times.values())
        fastest = times[winner]
        improvement = ((slowest - fastest) / slowest * 100) if slowest > 0 else 0

        orig_str = f"{times.get(0, 0):.3f}" if 0 in times else "N/A"
        symb_str = f"{times.get(1, 0):.3f}" if 1 in times else "N/A"
        cse_str = f"{times.get(2, 0):.3f}" if 2 in times else "N/A"

        lines.append(f"{sp:<8} {L:<4} {orig_str:<16} {symb_str:<16} {cse_str:<16} {IMPL_NAMES[winner]:<12} {improvement:>6.1f}%")

    lines.append("-" * 90)
    lines.append("")
    lines.append("Summary:")
    lines.append(f"  Original TMP wins: {total_wins[0]} shell pairs")
    lines.append(f"  Symbolic wins:     {total_wins[1]} shell pairs")
    lines.append(f"  Layered CSE wins:  {total_wins[2]} shell pairs")
    lines.append("")
    lines.append("=" * 90)

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))

    print(f"Saved: {output_path}")

    # Also print to console
    for line in lines:
        print(line)


def main():
    """Main entry point."""
    # Paths
    json_path = Path('/home/ruben/Research/Science/Projects/RECURSUM/benchmarks/results/raw/mcmd_realistic_100reps.json')
    output_dir = Path('/home/ruben/Research/Science/Projects/RECURSUM/benchmarks/results/figures')

    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print("Loading benchmark data...")
    results = load_benchmark_data(json_path)
    print(f"Loaded data for {len(results)} shell pairs")

    # Generate plots
    print("\nGenerating plots...")

    plot_shell_pair_comparison(
        results,
        output_dir / 'shell_pair_comparison.pdf'
    )

    plot_speedup_vs_L(
        results,
        output_dir / 'speedup_vs_L.pdf'
    )

    plot_timing_vs_L(
        results,
        output_dir / 'timing_vs_L.pdf'
    )

    generate_summary_table(
        results,
        output_dir / 'summary_table.txt'
    )

    print("\nAll plots generated successfully!")


if __name__ == '__main__':
    main()
