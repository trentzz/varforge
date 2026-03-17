#!/usr/bin/env python3
"""Generate publication-quality benchmark figures for the VarForge paper."""

import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# Style: clean, publication-ready
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.linewidth': 0.8,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

BLUE = '#2E86AB'
RED = '#C0392B'
GREEN = '#27AE60'
ORANGE = '#E67E22'
GRAY = '#7F8C8D'
PURPLE = '#8E44AD'

with open('../../benchmark_output/results.json') as f:
    data = json.load(f)


def fig_coverage_scaling():
    """Figure 1: Wall time and memory vs coverage depth (dual y-axis)."""
    covs = [1, 5, 10, 30, 50, 100, 200]
    walls = []
    mems = []
    for c in covs:
        d = data['coverage_scaling'][f'cov_{c}x']
        walls.append(d['avg_wall_secs'])
        mems.append(d['avg_peak_mem_kb'] / 1024)  # MB

    fig, ax1 = plt.subplots(figsize=(5.5, 3.5))

    ax1.plot(covs, walls, 'o-', color=BLUE, linewidth=1.5, markersize=5,
             label='Wall time', zorder=3)
    ax1.set_xlabel('Coverage depth (×)')
    ax1.set_ylabel('Wall time (s)', color=BLUE)
    ax1.tick_params(axis='y', labelcolor=BLUE)
    ax1.set_xlim(-5, 210)

    ax2 = ax1.twinx()
    ax2.plot(covs, mems, 's--', color=RED, linewidth=1.5, markersize=5,
             label='Peak memory', zorder=3)
    ax2.set_ylabel('Peak memory (MB)', color=RED)
    ax2.tick_params(axis='y', labelcolor=RED)

    # Linear fit for wall time
    coeffs = np.polyfit(covs, walls, 1)
    fit_x = np.linspace(0, 210, 100)
    fit_y = np.polyval(coeffs, fit_x)
    ax1.plot(fit_x, fit_y, ':', color=BLUE, alpha=0.4, linewidth=1)

    # R² annotation
    ss_res = sum((w - np.polyval(coeffs, c))**2 for c, w in zip(covs, walls))
    ss_tot = sum((w - np.mean(walls))**2 for w in walls)
    r2 = 1 - ss_res / ss_tot
    ax1.text(0.05, 0.92, f'$R^2 = {r2:.4f}$', transform=ax1.transAxes,
             fontsize=9, color=GRAY)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left',
               framealpha=0.9)

    ax1.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig('fig_coverage_scaling.pdf')
    fig.savefig('fig_coverage_scaling.png')
    plt.close(fig)
    print('  fig_coverage_scaling.pdf')


def fig_feature_overhead():
    """Figure 2: Feature overhead bar chart."""
    features = ['Baseline', '+Variants\n(500)', '+UMI\nsimplex',
                '+FFPE/\noxoG', '+GC\nbias', 'All\ncombined']
    keys = ['feat_base', 'feat_variants', 'feat_umi',
            'feat_artifacts', 'feat_gc_bias', 'feat_all']
    times = [data['feature_overhead'][k]['avg_wall_secs'] for k in keys]
    base = times[0]
    overheads = [(t - base) / base * 100 for t in times]

    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    colors = [GRAY, BLUE, ORANGE, RED, GREEN, PURPLE]
    bars = ax.bar(range(len(features)), times, color=colors, width=0.65,
                  edgecolor='white', linewidth=0.5)

    # Overhead labels
    for i, (bar, ovh) in enumerate(zip(bars, overheads)):
        if i == 0:
            label = 'baseline'
        else:
            label = f'+{ovh:.0f}%'
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                label, ha='center', va='bottom', fontsize=8, color='#333')

    ax.set_xticks(range(len(features)))
    ax.set_xticklabels(features, fontsize=8)
    ax.set_ylabel('Wall time (s)')
    ax.set_ylim(0, max(times) * 1.2)
    ax.grid(axis='y', alpha=0.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig('fig_feature_overhead.pdf')
    fig.savefig('fig_feature_overhead.png')
    plt.close(fig)
    print('  fig_feature_overhead.pdf')


def fig_thread_scaling():
    """Figure 3: Thread scaling with I/O-bound annotation."""
    threads = [1, 2, 4, 6, 8, 12]
    walls = [data['thread_scaling'][f'threads_{t}']['avg_wall_secs']
             for t in threads]
    base = walls[0]
    speedups = [base / w for w in walls]

    fig, ax = plt.subplots(figsize=(4.5, 3.5))

    # Ideal line
    ideal = list(range(1, 13))
    ax.plot(ideal, ideal, '--', color=GRAY, linewidth=0.8, alpha=0.5,
            label='Ideal (linear)')

    # Actual
    ax.plot(threads, speedups, 'o-', color=BLUE, linewidth=1.5,
            markersize=6, label='VarForge (10 MB, 30×)', zorder=3)

    ax.set_xlabel('Thread count')
    ax.set_ylabel('Speedup (×)')
    ax.set_xlim(0.5, 12.5)
    ax.set_ylim(0.5, 3.5)
    ax.set_xticks(threads)

    # I/O bound annotation
    ax.annotate('I/O-bound\nplateau',
                xy=(6, speedups[3]), xytext=(8.5, 1.8),
                fontsize=8, color=GRAY,
                arrowprops=dict(arrowstyle='->', color=GRAY, lw=0.8))

    ax.legend(loc='upper left', framealpha=0.9)
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    fig.savefig('fig_thread_scaling.pdf')
    fig.savefig('fig_thread_scaling.png')
    plt.close(fig)
    print('  fig_thread_scaling.pdf')


def fig_memory_scaling():
    """Figure 4: Memory per read pair across coverage depths."""
    covs = [1, 5, 10, 30, 50, 100, 200]
    # Approximate read pairs for 1MB ref at each coverage
    # 1MB ref, 150bp reads, 300bp fragments → pairs = coverage * 1e6 / 300
    pairs = [int(c * 1e6 / (300/2)) for c in covs]
    mems = [data['coverage_scaling'][f'cov_{c}x']['avg_peak_mem_kb'] / 1024
            for c in covs]

    # KB per pair
    kb_per_pair = [data['coverage_scaling'][f'cov_{c}x']['avg_peak_mem_kb'] / p
                   for c, p in zip(covs, pairs)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 3.0))

    # Left: Memory vs read pairs
    ax1.plot([p/1000 for p in pairs], mems, 'o-', color=RED, linewidth=1.5,
             markersize=5)
    ax1.set_xlabel('Read pairs (×1000)')
    ax1.set_ylabel('Peak memory (MB)')
    ax1.grid(True, alpha=0.2)

    # Linear fit
    coeffs = np.polyfit(pairs, mems, 1)
    r2_mem = 1 - sum((m - np.polyval(coeffs, p))**2 for p, m in zip(pairs, mems)) / \
             sum((m - np.mean(mems))**2 for m in mems)
    ax1.text(0.05, 0.88, f'$R^2 = {r2_mem:.4f}$', transform=ax1.transAxes,
             fontsize=9, color=GRAY)
    ax1.text(0.05, 0.76, f'{coeffs[0]*1024:.1f} KB/pair',
             transform=ax1.transAxes, fontsize=9, color=GRAY)

    # Right: KB per pair (should be constant)
    ax2.bar(range(len(covs)), kb_per_pair, color=GREEN, width=0.6,
            edgecolor='white')
    ax2.set_xticks(range(len(covs)))
    ax2.set_xticklabels([f'{c}×' for c in covs], fontsize=8)
    ax2.set_xlabel('Coverage depth')
    ax2.set_ylabel('Memory per read pair (KB)')
    ax2.set_ylim(0, max(kb_per_pair) * 1.3)
    ax2.axhline(y=np.mean(kb_per_pair[2:]), color=GRAY, linestyle='--',
                linewidth=0.8, alpha=0.5)
    ax2.text(0.5, 0.88, f'Mean: {np.mean(kb_per_pair[2:]):.2f} KB/pair',
             transform=ax2.transAxes, fontsize=8, color=GRAY, ha='center')
    ax2.grid(axis='y', alpha=0.2)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig('fig_memory_scaling.pdf')
    fig.savefig('fig_memory_scaling.png')
    plt.close(fig)
    print('  fig_memory_scaling.pdf')


def fig_throughput_by_scenario():
    """Figure 5: Throughput across main benchmark scenarios."""
    scenarios = {
        '01 Baseline\n(10MB,30×)': ('01_baseline', 10e6, 30),
        '02 +Variants\n(10MB,30×)': ('02_with_variants', 10e6, 30),
        '04 High cov\n(1MB,200×)': ('04_very_high_coverage', 1e6, 200),
        '07 cfDNA\n(1MB,200×)': ('07_cfdna', 1e6, 200),
        '08 FFPE\n(10MB,30×)': ('08_ffpe_artifacts', 10e6, 30),
    }

    names = list(scenarios.keys())
    throughputs = []
    for name, (key, ref_size, cov) in scenarios.items():
        d = data['main_benchmarks'][key]
        pairs = cov * ref_size / 150  # approximate
        rps = pairs / d['avg_wall_secs']
        throughputs.append(rps)

    fig, ax = plt.subplots(figsize=(5.5, 3.0))
    colors = [BLUE, BLUE, ORANGE, GREEN, RED]
    bars = ax.barh(range(len(names)), [t/1000 for t in throughputs],
                   color=colors, height=0.6, edgecolor='white')

    for bar, t in zip(bars, throughputs):
        ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height()/2,
                f'{t/1000:.1f}K', va='center', fontsize=8, color='#333')

    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel('Throughput (×1000 read pairs/s)')
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig('fig_throughput_scenarios.pdf')
    fig.savefig('fig_throughput_scenarios.png')
    plt.close(fig)
    print('  fig_throughput_scenarios.pdf')


if __name__ == '__main__':
    print('Generating figures...')
    fig_coverage_scaling()
    fig_feature_overhead()
    fig_thread_scaling()
    fig_memory_scaling()
    fig_throughput_by_scenario()
    print('Done.')
