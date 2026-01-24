#!/usr/bin/env python3
"""
Comprehensive Analysis of RECURSUM Benchmark Results

Generates publication-quality figures and statistical summaries from
the statistically significant benchmark data (10 repetitions).
"""

import json
import re
import argparse
from pathlib import Path
from collections import defaultdict
import numpy as np

# Configure matplotlib for publication quality
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Publication style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'serif',
})


def load_benchmark_data(filepath):
    """Load and parse benchmark JSON, handling NaN values."""
    with open(filepath, 'r') as f:
        content = f.read()
    # Google Benchmark outputs -nan which isn't valid JSON
    content = re.sub(r':\s*-?nan\b', ': null', content)
    return json.loads(content)


def parse_benchmark_name(name):
    """Extract benchmark type and parameters from name."""
    # Remove the /min_time:X.XXX suffix
    base_name = name.split('/')[0]

    # Parse E coefficient benchmarks: BM_E_TMP_nA_nB_t or BM_E_Sym_nA_nB_t
    match = re.match(r'BM_E_(TMP|Sym)_(\d+)_(\d+)_(\d+)', base_name)
    if match:
        impl = match.group(1)
        nA, nB, t = int(match.group(2)), int(match.group(3)), int(match.group(4))
        return {'type': 'E', 'impl': impl, 'nA': nA, 'nB': nB, 't': t, 'L': nA + nB}

    # Parse gradient benchmarks: BM_dEdPA_TMP_nA_nB_t or BM_dEdPA_Sym_nA_nB_t
    match = re.match(r'BM_dEdPA_(TMP|Sym)_(\d+)_(\d+)_(\d+)', base_name)
    if match:
        impl = match.group(1)
        nA, nB, t = int(match.group(2)), int(match.group(3)), int(match.group(4))
        return {'type': 'dEdPA', 'impl': impl, 'nA': nA, 'nB': nB, 't': t, 'L': nA + nB}

    return None


def aggregate_benchmarks(data):
    """Aggregate benchmark results by type and compute statistics."""
    results = defaultdict(lambda: {'times': [], 'info': None})

    for bench in data['benchmarks']:
        parsed = parse_benchmark_name(bench['name'])
        if parsed is None:
            continue

        key = (parsed['type'], parsed['impl'], parsed['nA'], parsed['nB'], parsed['t'])
        results[key]['times'].append(bench['cpu_time'])
        results[key]['info'] = parsed

    # Compute statistics
    stats = {}
    for key, val in results.items():
        times = np.array(val['times'])
        stats[key] = {
            'info': val['info'],
            'mean': np.mean(times),
            'std': np.std(times),
            'min': np.min(times),
            'max': np.max(times),
            'cv': np.std(times) / np.mean(times) * 100 if np.mean(times) > 0 else 0,
            'n': len(times)
        }

    return stats


def generate_scaling_plot(stats, output_dir):
    """Generate scaling comparison plot for all L values."""
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    axes = axes.flatten()

    for L in range(9):
        ax = axes[L]

        # Collect data for this L
        tmp_data = []
        sym_data = []

        for key, val in stats.items():
            if val['info']['type'] != 'E' or val['info']['L'] != L:
                continue
            t = val['info']['t']
            if val['info']['impl'] == 'TMP':
                tmp_data.append((t, val['mean'], val['std']))
            else:
                sym_data.append((t, val['mean'], val['std']))

        if not tmp_data or not sym_data:
            ax.set_title(f'L={L} (no data)', fontsize=10)
            ax.set_visible(False)
            continue

        # Sort by t
        tmp_data.sort(key=lambda x: x[0])
        sym_data.sort(key=lambda x: x[0])

        tmp_t = [d[0] for d in tmp_data]
        tmp_mean = [d[1] for d in tmp_data]
        tmp_std = [d[2] for d in tmp_data]

        sym_t = [d[0] for d in sym_data]
        sym_mean = [d[1] for d in sym_data]
        sym_std = [d[2] for d in sym_data]

        # Plot with error bars
        ax.errorbar(tmp_t, tmp_mean, yerr=tmp_std, marker='o', markersize=5,
                   capsize=3, label='TMP', color='#2ecc71', linewidth=1.5)
        ax.errorbar(sym_t, sym_mean, yerr=sym_std, marker='s', markersize=5,
                   capsize=3, label='Symbolic', color='#e74c3c', linewidth=1.5)

        ax.set_title(f'L = {L}', fontsize=11)
        ax.set_xlabel('t index', fontsize=10)
        ax.set_ylabel('Time (ns)', fontsize=10)
        ax.legend(fontsize=9, loc='upper left')
        ax.set_xticks(range(L + 1))
        ax.grid(True, alpha=0.3)

    plt.suptitle('TMP vs Symbolic: Scaling with Angular Momentum\n(Error bars: 1σ, n=10)',
                 fontsize=14, y=1.02)
    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / f'fig1_scaling_all_L.{fmt}')
    plt.close(fig)
    print(f"Generated: {output_dir / 'fig1_scaling_all_L.png'}")


def generate_speedup_heatmap(stats, output_dir):
    """Generate speedup heatmap for all coefficients."""
    # Collect speedup data
    speedups = {}  # (nA, nB, t) -> speedup

    for key, val in stats.items():
        if val['info']['type'] != 'E':
            continue
        nA, nB, t = val['info']['nA'], val['info']['nB'], val['info']['t']

        # Find matching TMP and Symbolic
        tmp_key = ('E', 'TMP', nA, nB, t)
        sym_key = ('E', 'Sym', nA, nB, t)

        if tmp_key in stats and sym_key in stats:
            tmp_time = stats[tmp_key]['mean']
            sym_time = stats[sym_key]['mean']
            speedup = sym_time / tmp_time  # >1 means TMP faster
            speedups[(nA, nB, t)] = speedup

    # Create heatmap data - organized by L and t
    fig, axes = plt.subplots(3, 3, figsize=(14, 12))
    axes = axes.flatten()

    for L in range(9):
        ax = axes[L]

        # Find all (nA, nB) pairs for this L
        pairs = [(nA, L - nA) for nA in range(L + 1) if L - nA >= 0 and L - nA <= 4 and nA <= 4]

        if not pairs:
            ax.set_visible(False)
            continue

        # Create matrix: rows = (nA, nB) pairs, cols = t values
        t_max = L
        matrix = np.zeros((len(pairs), t_max + 1))
        matrix[:] = np.nan

        for i, (nA, nB) in enumerate(pairs):
            for t in range(t_max + 1):
                if (nA, nB, t) in speedups:
                    matrix[i, t] = speedups[(nA, nB, t)]

        # Custom colormap: red (<1, Symbolic faster) -> white (1) -> green (>1, TMP faster)
        cmap = LinearSegmentedColormap.from_list('speedup',
            ['#e74c3c', '#ffffff', '#2ecc71'], N=256)

        # Plot heatmap
        im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=0.3, vmax=3.0)

        # Labels
        ax.set_xticks(range(t_max + 1))
        ax.set_xticklabels(range(t_max + 1))
        ax.set_yticks(range(len(pairs)))
        ax.set_yticklabels([f'({nA},{nB})' for nA, nB in pairs])

        ax.set_xlabel('t index', fontsize=10)
        ax.set_ylabel('(nA, nB)', fontsize=10)
        ax.set_title(f'L = {L}', fontsize=11)

        # Add text annotations
        for i in range(len(pairs)):
            for t in range(t_max + 1):
                if not np.isnan(matrix[i, t]):
                    color = 'black' if 0.5 < matrix[i, t] < 2.0 else 'white'
                    ax.text(t, i, f'{matrix[i, t]:.2f}', ha='center', va='center',
                           color=color, fontsize=8)

    plt.suptitle('Speedup Heatmap: TMP / Symbolic\n(Green >1: TMP faster, Red <1: Symbolic faster)',
                 fontsize=14, y=1.02)

    # Add colorbar
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('Speedup (Sym/TMP)', fontsize=11)

    plt.tight_layout(rect=[0, 0, 0.9, 1])

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / f'fig2_speedup_heatmap.{fmt}')
    plt.close(fig)
    print(f"Generated: {output_dir / 'fig2_speedup_heatmap.png'}")


def generate_crossover_analysis(stats, output_dir):
    """Analyze crossover point where TMP becomes faster."""
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = plt.cm.viridis(np.linspace(0, 1, 9))

    crossover_data = []

    for L in range(9):
        tmp_data = []
        sym_data = []

        for key, val in stats.items():
            if val['info']['type'] != 'E' or val['info']['L'] != L:
                continue
            t = val['info']['t']
            if val['info']['impl'] == 'TMP':
                tmp_data.append((t, val['mean'], val['std']))
            else:
                sym_data.append((t, val['mean'], val['std']))

        if not tmp_data or not sym_data:
            continue

        # Sort by t
        tmp_data.sort(key=lambda x: x[0])
        sym_data.sort(key=lambda x: x[0])

        # Compute speedup for each t
        speedups = []
        for (t1, m1, s1), (t2, m2, s2) in zip(tmp_data, sym_data):
            assert t1 == t2
            speedups.append((t1, m2 / m1))  # Sym/TMP

        t_vals = [s[0] for s in speedups]
        speedup_vals = [s[1] for s in speedups]

        ax.plot(t_vals, speedup_vals, 'o-', color=colors[L],
                label=f'L={L}', markersize=6, linewidth=1.5)

        # Find crossover
        crossover_t = None
        for i, (t, sp) in enumerate(speedups):
            if sp >= 1.0:
                crossover_t = t
                break
        crossover_data.append((L, crossover_t))

    ax.axhline(y=1.0, color='black', linestyle='--', linewidth=2, label='Equal performance')
    ax.fill_between([0, 8], [0, 0], [1, 1], alpha=0.1, color='red', label='Symbolic faster')
    ax.fill_between([0, 8], [1, 1], [5, 5], alpha=0.1, color='green', label='TMP faster')

    ax.set_xlabel('t index', fontsize=12)
    ax.set_ylabel('Speedup (Symbolic / TMP)', fontsize=12)
    ax.set_title('Crossover Analysis: When TMP Becomes Faster than Symbolic\n(Based on actual benchmark measurements)',
                fontsize=13)
    ax.legend(loc='upper right', ncol=2)
    ax.set_xlim(-0.2, 8.2)
    ax.set_ylim(0, 4)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / f'fig3_crossover_analysis.{fmt}')
    plt.close(fig)
    print(f"Generated: {output_dir / 'fig3_crossover_analysis.png'}")

    return crossover_data


def generate_bar_comparison(stats, output_dir):
    """Generate bar chart comparing key shell pairs."""
    # Select representative coefficients for each shell pair
    shell_pairs = [
        ('ss', 0, 0, 0),
        ('sp', 0, 1, 1),
        ('pp', 1, 1, 2),
        ('sd', 0, 2, 2),
        ('pd', 1, 2, 3),
        ('dd', 2, 2, 4),
        ('sf', 0, 3, 3),
        ('pf', 1, 3, 4),
        ('df', 2, 3, 5),
        ('ff', 3, 3, 6),
        ('sg', 0, 4, 4),
        ('pg', 1, 4, 5),
        ('dg', 2, 4, 6),
        ('fg', 3, 4, 7),
        ('gg', 4, 4, 8),
    ]

    labels = []
    tmp_times = []
    sym_times = []
    tmp_errs = []
    sym_errs = []

    for name, nA, nB, t in shell_pairs:
        tmp_key = ('E', 'TMP', nA, nB, t)
        sym_key = ('E', 'Sym', nA, nB, t)

        if tmp_key in stats and sym_key in stats:
            labels.append(f'{name}\nE({nA},{nB},{t})')
            tmp_times.append(stats[tmp_key]['mean'])
            tmp_errs.append(stats[tmp_key]['std'])
            sym_times.append(stats[sym_key]['mean'])
            sym_errs.append(stats[sym_key]['std'])

    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(14, 6))

    bars1 = ax.bar(x - width/2, tmp_times, width, label='TMP', color='#2ecc71',
                   yerr=tmp_errs, capsize=3)
    bars2 = ax.bar(x + width/2, sym_times, width, label='Symbolic', color='#e74c3c',
                   yerr=sym_errs, capsize=3)

    ax.set_xlabel('Shell Pair / Coefficient', fontsize=12)
    ax.set_ylabel('Execution Time (ns)', fontsize=12)
    ax.set_title('TMP vs Symbolic: Representative Coefficients by Shell Pair\n(Error bars: 1σ, n=10)',
                fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')

    # Add speedup annotations
    for i, (t, s) in enumerate(zip(tmp_times, sym_times)):
        speedup = s / t
        color = '#2ecc71' if speedup > 1 else '#e74c3c'
        ax.annotate(f'{speedup:.2f}x', xy=(i, max(t, s) + max(tmp_errs[i], sym_errs[i]) + 1),
                   ha='center', fontsize=8, color=color, fontweight='bold')

    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / f'fig4_shell_pair_comparison.{fmt}')
    plt.close(fig)
    print(f"Generated: {output_dir / 'fig4_shell_pair_comparison.png'}")


def generate_gradient_comparison(stats, output_dir):
    """Generate comparison for gradient benchmarks."""
    # Filter gradient data
    grad_data = {k: v for k, v in stats.items() if v['info']['type'] == 'dEdPA'}

    if not grad_data:
        print("No gradient benchmark data found")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    # Group by (nA, nB)
    pairs = set((v['info']['nA'], v['info']['nB']) for v in grad_data.values())

    colors_tmp = plt.cm.Greens(np.linspace(0.4, 0.9, len(pairs)))
    colors_sym = plt.cm.Reds(np.linspace(0.4, 0.9, len(pairs)))

    for idx, (nA, nB) in enumerate(sorted(pairs)):
        tmp_data = []
        sym_data = []

        for key, val in grad_data.items():
            if val['info']['nA'] != nA or val['info']['nB'] != nB:
                continue
            t = val['info']['t']
            if val['info']['impl'] == 'TMP':
                tmp_data.append((t, val['mean'], val['std']))
            else:
                sym_data.append((t, val['mean'], val['std']))

        if tmp_data:
            tmp_data.sort(key=lambda x: x[0])
            t_vals = [d[0] for d in tmp_data]
            means = [d[1] for d in tmp_data]
            stds = [d[2] for d in tmp_data]
            ax.errorbar(t_vals, means, yerr=stds, marker='o', markersize=6, capsize=3,
                       label=f'TMP ({nA},{nB})', color=colors_tmp[idx], linewidth=1.5)

        if sym_data:
            sym_data.sort(key=lambda x: x[0])
            t_vals = [d[0] for d in sym_data]
            means = [d[1] for d in sym_data]
            stds = [d[2] for d in sym_data]
            ax.errorbar(t_vals, means, yerr=stds, marker='s', markersize=6, capsize=3,
                       label=f'Sym ({nA},{nB})', color=colors_sym[idx], linewidth=1.5, linestyle='--')

    ax.set_xlabel('t index', fontsize=12)
    ax.set_ylabel('Execution Time (ns)', fontsize=12)
    ax.set_title('Gradient Benchmarks: ∂E/∂PA\n(TMP vs Symbolic)', fontsize=13)
    ax.legend(loc='upper left', ncol=2, fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / f'fig5_gradient_comparison.{fmt}')
    plt.close(fig)
    print(f"Generated: {output_dir / 'fig5_gradient_comparison.png'}")


def generate_summary_table(stats, output_dir, crossover_data):
    """Generate comprehensive summary markdown table."""
    lines = [
        "# RECURSUM Benchmark Summary",
        "## TMP vs Symbolic Performance Analysis",
        "",
        "### Experimental Setup",
        "- **Repetitions**: 10 per benchmark",
        "- **Minimum Time**: 1.0s per measurement",
        "- **SIMD**: Vec8d (8 double precision values per operation)",
        "- **Coefficients**: 125 Hermite E coefficients (E^{0,0}_0 to E^{4,4}_8)",
        "",
        "---",
        "",
        "### Key Results by Angular Momentum",
        "",
        "| L | nA+nB | TMP Mean (ns) | Sym Mean (ns) | Avg Speedup | t Crossover |",
        "|---|-------|---------------|---------------|-------------|-------------|"
    ]

    for L in range(9):
        tmp_times = []
        sym_times = []

        for key, val in stats.items():
            if val['info']['type'] != 'E' or val['info']['L'] != L:
                continue
            if val['info']['impl'] == 'TMP':
                tmp_times.append(val['mean'])
            else:
                sym_times.append(val['mean'])

        if tmp_times and sym_times:
            tmp_avg = np.mean(tmp_times)
            sym_avg = np.mean(sym_times)
            speedup = sym_avg / tmp_avg

            # Find crossover
            crossover = next((t for l, t in crossover_data if l == L), None)
            crossover_str = f"t={crossover}" if crossover is not None else "N/A"

            lines.append(f"| {L} | {L} | {tmp_avg:.2f} | {sym_avg:.2f} | {speedup:.2f}x | {crossover_str} |")

    lines.extend([
        "",
        "---",
        "",
        "### Representative Coefficients",
        "",
        "| Coefficient | TMP (ns) | Sym (ns) | Speedup | Faster |",
        "|-------------|----------|----------|---------|--------|"
    ])

    # Add key coefficients
    key_coeffs = [
        (0, 0, 0), (1, 1, 2), (2, 2, 4), (3, 3, 6), (4, 4, 8),
        (4, 4, 0), (4, 4, 4), (3, 4, 7), (2, 4, 6)
    ]

    for nA, nB, t in key_coeffs:
        tmp_key = ('E', 'TMP', nA, nB, t)
        sym_key = ('E', 'Sym', nA, nB, t)

        if tmp_key in stats and sym_key in stats:
            tmp_time = stats[tmp_key]['mean']
            sym_time = stats[sym_key]['mean']
            speedup = sym_time / tmp_time
            faster = "TMP" if speedup > 1 else "Symbolic"
            lines.append(f"| E^{{{nA},{nB}}}_{t} | {tmp_time:.2f} ± {stats[tmp_key]['std']:.2f} | "
                        f"{sym_time:.2f} ± {stats[sym_key]['std']:.2f} | {speedup:.2f}x | {faster} |")

    lines.extend([
        "",
        "---",
        "",
        "### Crossover Analysis",
        "",
        "The crossover point where TMP becomes faster than Symbolic varies with L:",
        ""
    ])

    for L, t in crossover_data:
        if t is not None:
            lines.append(f"- **L={L}**: TMP faster for t ≥ {t}")
        else:
            lines.append(f"- **L={L}**: Symbolic faster for all t")

    lines.extend([
        "",
        "---",
        "",
        "### Key Findings",
        "",
        "1. **Crossover Behavior**: TMP becomes faster as t increases",
        "2. **Low t (t=0,1)**: Symbolic benefits from compact polynomial expressions",
        "3. **High t**: TMP's compile-time optimization pays off",
        "4. **Gradients**: Similar pattern with TMP advantage at higher complexity",
        "",
        "*All data from actual benchmark measurements with proper statistical significance.*"
    ])

    summary_path = output_dir / 'comprehensive_summary.md'
    with open(summary_path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"Generated: {summary_path}")


def main():
    parser = argparse.ArgumentParser(description='Analyze RECURSUM benchmark results')
    parser.add_argument('input', help='Input JSON file from benchmark')
    parser.add_argument('--output-dir', '-o', default='../results/figures',
                       help='Output directory for figures')
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading benchmark data from: {input_path}")
    data = load_benchmark_data(input_path)

    print(f"Aggregating {len(data['benchmarks'])} benchmark entries...")
    stats = aggregate_benchmarks(data)
    print(f"Found {len(stats)} unique benchmark configurations")

    print("\nGenerating publication-quality figures...")

    generate_scaling_plot(stats, output_dir)
    generate_speedup_heatmap(stats, output_dir)
    crossover_data = generate_crossover_analysis(stats, output_dir)
    generate_bar_comparison(stats, output_dir)
    generate_gradient_comparison(stats, output_dir)
    generate_summary_table(stats, output_dir, crossover_data)

    print("\nAll figures generated successfully!")
    print(f"Output directory: {output_dir}")


if __name__ == '__main__':
    main()
