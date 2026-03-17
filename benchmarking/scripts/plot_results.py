#!/usr/bin/env python3
"""Generate publication-quality benchmark figures for VarForge.

Reads data from benchmarking/results/tables/.
Writes figures to benchmarking/results/graphs/.

Style follows ~/.claude/guides/graph-style.md:
  Okabe-Ito palette, white background, horizontal grid only,
  no top/right spines, serif font, PDF output.
"""

import csv
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# -- Okabe-Ito palette (colourblind-safe) --
PALETTE = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # green
    "#CC79A7",  # pink
    "#56B4E9",  # sky
    "#D55E00",  # vermillion
    "#F0E442",  # yellow
]

mpl.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "#FAFAFA",
    "axes.edgecolor": "#CCCCCC",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.color": "#E0E0E0",
    "grid.linewidth": 0.8,
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "lines.linewidth": 1.8,
    "axes.prop_cycle": mpl.cycler(color=PALETTE),
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
})

BLUE, ORANGE, GREEN, PINK, SKY, VERMILLION, YELLOW = PALETTE

TABLES = os.path.join(os.path.dirname(__file__), '..', 'results', 'tables')
GRAPHS = os.path.join(os.path.dirname(__file__), '..', 'results', 'graphs')
os.makedirs(GRAPHS, exist_ok=True)


def load_json():
    path = os.path.join(TABLES, 'results.json')
    with open(path) as f:
        return json.load(f)


def load_csv(name):
    path = os.path.join(TABLES, name)
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def save(fig, name):
    fig.savefig(os.path.join(GRAPHS, name + '.pdf'))
    fig.savefig(os.path.join(GRAPHS, name + '.png'))
    plt.close(fig)
    print(f'  {name}.pdf')


def fig_coverage_scaling(data):
    """Wall time and peak memory vs coverage depth."""
    covs = [1, 5, 10, 30, 50, 100, 200]
    walls, mems = [], []
    for c in covs:
        d = data['coverage_scaling'][f'cov_{c}x']
        walls.append(d['avg_wall_secs'])
        mems.append(d['avg_peak_mem_kb'] / 1024)

    fig, ax1 = plt.subplots(figsize=(84/25.4, 65/25.4))  # 84mm wide

    ax1.plot(covs, walls, 'o-', color=BLUE, markersize=4, label='Wall time')
    ax1.set_xlabel('Coverage depth (x)')
    ax1.set_ylabel('Wall time (s)', color=BLUE)
    ax1.tick_params(axis='y', labelcolor=BLUE)

    # Linear fit
    coeffs = np.polyfit(covs, walls, 1)
    fit_x = np.linspace(0, 210, 100)
    ax1.plot(fit_x, np.polyval(coeffs, fit_x), ':', color=BLUE, alpha=0.3)

    ss_res = sum((w - np.polyval(coeffs, c))**2 for c, w in zip(covs, walls))
    ss_tot = sum((w - np.mean(walls))**2 for w in walls)
    r2 = 1 - ss_res / ss_tot
    ax1.text(0.03, 0.93, f'$R^2 = {r2:.4f}$', transform=ax1.transAxes,
             fontsize=8, color='#666666')

    ax2 = ax1.twinx()
    ax2.plot(covs, mems, 's--', color=VERMILLION, markersize=4, label='Peak memory')
    ax2.set_ylabel('Peak memory (MB)', color=VERMILLION)
    ax2.tick_params(axis='y', labelcolor=VERMILLION)
    ax2.spines['right'].set_visible(True)
    ax2.spines['right'].set_color('#CCCCCC')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

    ax1.set_xlim(-5, 210)
    save(fig, 'coverage_scaling')


def fig_feature_overhead(data):
    """Bar chart of feature overhead."""
    labels = ['Baseline', '+Variants', '+UMI', '+FFPE/oxoG', '+GC bias', 'All']
    keys = ['feat_base', 'feat_variants', 'feat_umi',
            'feat_artifacts', 'feat_gc_bias', 'feat_all']
    times = [data['feature_overhead'][k]['avg_wall_secs'] for k in keys]
    base = times[0]

    fig, ax = plt.subplots(figsize=(84/25.4, 55/25.4))
    colors = [BLUE, ORANGE, GREEN, VERMILLION, SKY, PINK]
    bars = ax.bar(range(len(labels)), times, color=colors, width=0.65,
                  edgecolor='white', linewidth=0.5)

    for i, (bar, t) in enumerate(zip(bars, times)):
        if i == 0:
            label = 'baseline'
        else:
            label = f'+{(t - base) / base * 100:.0f}%'
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                label, ha='center', va='bottom', fontsize=7, color='#444444')

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7, rotation=15, ha='right')
    ax.set_ylabel('Wall time (s)')
    ax.set_ylim(0, max(times) * 1.2)
    save(fig, 'feature_overhead')


def fig_thread_scaling(data):
    """Thread scaling: shows I/O-bound slowdown with more threads."""
    threads = [1, 2, 4, 6, 8, 12]
    walls = [data['thread_scaling'][f'threads_{t}']['avg_wall_secs']
             for t in threads]
    base = walls[0]
    rel = [base / w for w in walls]  # values <1 mean slower than 1 thread

    fig, ax = plt.subplots(figsize=(70/25.4, 55/25.4))

    ax.axhline(y=1.0, color='#CCCCCC', linestyle='--', linewidth=0.8, label='1-thread baseline')
    ax.plot(threads, rel, 'o-', color=BLUE, markersize=5,
            label='VarForge (10 MB, 30x)')

    ax.set_xlabel('Thread count')
    ax.set_ylabel('Relative speed (1T = 1.0)')
    ax.set_xlim(0.5, 12.5)
    ax.set_ylim(0.6, 1.2)
    ax.set_xticks(threads)

    ax.annotate('I/O-bound:\nmore threads = slower', xy=(12, rel[-1]),
                xytext=(7, 0.72), fontsize=7, color='#888888',
                arrowprops=dict(arrowstyle='->', color='#AAAAAA', lw=0.8))

    ax.legend(loc='upper right', fontsize=8)
    save(fig, 'thread_scaling')


def fig_memory_per_pair(data):
    """Memory per read pair across coverage depths."""
    covs = [1, 5, 10, 30, 50, 100, 200]
    pairs = [int(c * 1e6 / 150) for c in covs]
    mems_mb = [data['coverage_scaling'][f'cov_{c}x']['avg_peak_mem_kb'] / 1024
               for c in covs]
    kb_per_pair = [data['coverage_scaling'][f'cov_{c}x']['avg_peak_mem_kb'] / p
                   for c, p in zip(covs, pairs)]

    fig, ax = plt.subplots(figsize=(84/25.4, 55/25.4))
    ax.bar(range(len(covs)), kb_per_pair, color=GREEN, width=0.6,
           edgecolor='white')
    ax.set_xticks(range(len(covs)))
    ax.set_xticklabels([f'{c}x' for c in covs], fontsize=8)
    ax.set_xlabel('Coverage depth')
    ax.set_ylabel('Memory per read pair (KB)')
    ax.set_ylim(0, max(kb_per_pair) * 1.3)

    mean_val = np.mean(kb_per_pair[2:])  # skip 1x and 5x (startup dominated)
    ax.axhline(y=mean_val, color='#AAAAAA', linestyle='--', linewidth=0.8)
    ax.text(len(covs) - 1, mean_val + 0.15,
            f'{mean_val:.1f} KB/pair', fontsize=8, color='#666666', ha='right')

    save(fig, 'memory_per_pair')


def fig_throughput_scenarios(data):
    """Throughput across main scenarios."""
    scenarios = [
        ('Baseline\n(10 MB, 30x)', '01_baseline', 10e6, 30),
        ('+Variants\n(10 MB, 30x)', '02_with_variants', 10e6, 30),
        ('High cov\n(1 MB, 200x)', '04_very_high_coverage', 1e6, 200),
        ('cfDNA\n(1 MB, 200x)', '07_cfdna', 1e6, 200),
        ('FFPE\n(10 MB, 30x)', '08_ffpe_artifacts', 10e6, 30),
        ('Panel+UMI\n(1 MB, 200x)', '11_panel_umi', 1e6, 200),
    ]

    names, throughputs = [], []
    for name, key, ref_size, cov in scenarios:
        d = data['main_benchmarks'][key]
        pairs = cov * ref_size / 150
        rps = pairs / d['avg_wall_secs']
        names.append(name)
        throughputs.append(rps)

    fig, ax = plt.subplots(figsize=(84/25.4, 60/25.4))
    colors = [BLUE, ORANGE, GREEN, PINK, VERMILLION, SKY]
    bars = ax.barh(range(len(names)), [t/1000 for t in throughputs],
                   color=colors, height=0.55, edgecolor='white')

    for bar, t in zip(bars, throughputs):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2,
                f'{t/1000:.1f}K', va='center', fontsize=7, color='#444444')

    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel('Throughput (x1000 read pairs/s)')
    ax.invert_yaxis()
    ax.grid(axis='x')
    save(fig, 'throughput_scenarios')


if __name__ == '__main__':
    print('Generating figures...')
    data = load_json()
    fig_coverage_scaling(data)
    fig_feature_overhead(data)
    fig_thread_scaling(data)
    fig_memory_per_pair(data)
    fig_throughput_scenarios(data)
    print('Done.')
