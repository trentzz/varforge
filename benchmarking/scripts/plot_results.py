#!/usr/bin/env python3
"""Generate publication-quality benchmark figures for VarForge.

Reads data from benchmarking/results/tables/.
Writes figures to benchmarking/results/graphs/.

Style: Okabe-Ito palette, white background, horizontal grid only,
no top/right spines, DejaVu Sans, PDF + PNG output at 600 DPI.

Design decisions:
- No dual y-axes. Use two-panel (gridspec) figures instead.
- Prefer Cleveland dot plots over bar charts when zero is not a
  meaningful baseline (e.g. runtime comparisons across scenarios).
- Direct annotation instead of legend where feasible.
- Export width = 84 mm (single column) for all figures.
- Thread scaling is moved to appendix; RAM and throughput are
  the primary performance messages.
"""

import csv
import json
import os

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

BLUE, ORANGE, GREEN, PINK, SKY, VERMILLION, YELLOW = PALETTE

COL_W = 84 / 25.4          # 84 mm single-column width in inches

mpl.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#BBBBBB",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.color": "#E8E8E8",
    "grid.linewidth": 0.6,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "lines.linewidth": 1.0,
    "axes.prop_cycle": mpl.cycler(color=PALETTE),
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

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
    with open(path) as f:
        return list(csv.DictReader(f))


def save(fig, name):
    fig.savefig(os.path.join(GRAPHS, name + '.pdf'))
    fig.savefig(os.path.join(GRAPHS, name + '.png'))
    plt.close(fig)
    print(f'  {name}.pdf + .png')


# ---------------------------------------------------------------------------
# Figure 1: Coverage scaling — two panels, no dual axis
# Takeaway: both wall time and RAM grow linearly with coverage.
# ---------------------------------------------------------------------------
def fig_coverage_scaling(data):
    covs = [1, 5, 10, 30, 50, 100, 200]
    walls, mems = [], []
    for c in covs:
        d = data['coverage_scaling'][f'cov_{c}x']
        walls.append(d['avg_wall_secs'])
        mems.append(d['avg_peak_mem_kb'] / 1024)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(COL_W, COL_W * 0.62),
                                    gridspec_kw={'wspace': 0.45})

    def _fit_annotate(ax, xs, ys):
        coeffs = np.polyfit(xs, ys, 1)
        fx = np.linspace(0, 210, 200)
        ax.plot(fx, np.polyval(coeffs, fx), '-', color='#BBBBBB',
                linewidth=0.8, zorder=0)
        ss_res = sum((y - np.polyval(coeffs, x))**2 for x, y in zip(xs, ys))
        ss_tot = sum((y - np.mean(ys))**2 for y in ys)
        r2 = 1 - ss_res / ss_tot
        ax.text(0.97, 0.06, f'$R^2={r2:.4f}$', transform=ax.transAxes,
                fontsize=6.5, color='#888888', ha='right')

    # Wall time
    ax1.plot(covs, walls, 'o', color=BLUE, markersize=4, zorder=3)
    _fit_annotate(ax1, covs, walls)
    ax1.set_xlabel('Coverage (x)')
    ax1.set_ylabel('Wall time / s')
    ax1.set_xlim(-5, 215)
    ax1.set_ylim(bottom=0)

    # Peak RAM
    ax2.plot(covs, mems, 'o', color=VERMILLION, markersize=4, zorder=3)
    _fit_annotate(ax2, covs, mems)
    ax2.set_xlabel('Coverage (x)')
    ax2.set_ylabel('Peak RAM / MB')
    ax2.set_xlim(-5, 215)
    ax2.set_ylim(bottom=0)
    ax2.yaxis.grid(True, color='#E8E8E8', linewidth=0.6)

    save(fig, 'coverage_scaling')


# ---------------------------------------------------------------------------
# Figure 2: Feature overhead — two-panel horizontal dot chart
# Takeaway: UMI dominates; FFPE adds 15%; variants and GC bias are negligible.
# ---------------------------------------------------------------------------
def fig_feature_overhead(data):
    labels = ['Baseline', '+Variants', '+GC bias', '+FFPE/oxoG', '+UMI', 'All features']
    keys   = ['feat_base', 'feat_variants', 'feat_gc_bias',
              'feat_artifacts', 'feat_umi', 'feat_all']

    walls = [data['feature_overhead'][k]['avg_wall_secs'] for k in keys]
    rams  = [data['feature_overhead'][k]['avg_peak_mem_kb'] / 1024 for k in keys]
    base_t, base_r = walls[0], rams[0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(COL_W, COL_W * 0.72),
                                    gridspec_kw={'wspace': 0.08})

    y = range(len(labels))
    colors = [BLUE, SKY, GREEN, ORANGE, VERMILLION, PINK]

    # Wall time panel
    ax1.set_xlim(0, max(walls) * 1.35)
    for i, (yy, t, c) in enumerate(zip(y, walls, colors)):
        ax1.plot([0, t], [yy, yy], '-', color='#DDDDDD', linewidth=1.5, zorder=1)
        ax1.plot(t, yy, 'o', color=c, markersize=6, zorder=3)
        if i == 0:
            label = f'{t:.1f} s'
        else:
            pct = (t - base_t) / base_t * 100
            label = f'+{pct:.0f}%'
        ax1.text(t + max(walls) * 0.04, yy, label,
                 va='center', fontsize=6.5, color='#444444')
    ax1.set_yticks(list(y))
    ax1.set_yticklabels(labels, fontsize=7)
    ax1.set_xlabel('Wall time / s')
    ax1.invert_yaxis()
    ax1.grid(axis='x', color='#E8E8E8', linewidth=0.6)
    ax1.grid(axis='y', visible=False)
    ax1.spines['left'].set_visible(False)
    ax1.tick_params(left=False)

    # RAM panel
    ax2.set_xlim(0, max(rams) * 1.35)
    for i, (yy, r, c) in enumerate(zip(y, rams, colors)):
        ax2.plot([0, r], [yy, yy], '-', color='#DDDDDD', linewidth=1.5, zorder=1)
        ax2.plot(r, yy, 's', color=c, markersize=6, zorder=3)
        label = f'{r:.0f} MB'
        ax2.text(r + max(rams) * 0.04, yy, label,
                 va='center', fontsize=6.5, color='#444444')
    ax2.set_yticks(list(y))
    ax2.set_yticklabels([])
    ax2.set_xlabel('Peak RAM / MB')
    ax2.invert_yaxis()
    ax2.grid(axis='x', color='#E8E8E8', linewidth=0.6)
    ax2.grid(axis='y', visible=False)
    ax2.spines['left'].set_visible(False)
    ax2.tick_params(left=False)

    save(fig, 'feature_overhead')


# ---------------------------------------------------------------------------
# Figure 3: hg38 chr22 benchmark — RAM and throughput per config
# Takeaway: RAM ranges from 0.56 GB (BED-filtered panel) to 4.2 GB (30x WGS).
# ---------------------------------------------------------------------------
def fig_hg38_benchmark(rows):
    labels = [
        'WGS baseline (5x)',
        'WGS + variants (5x)',
        'Panel UMI (5x)',
        'Twist duplex (5x)',
        'cfDNA (5x)',
        'FFPE tumour (5x)',
        'WGS baseline (30x)',
        'Panel BED 200x',
    ]
    rams   = [float(r['peak_mem_kb']) / 1024 / 1024 for r in rows]  # GB
    pairss = [float(r['read_pairs']) / float(r['wall_secs'])
              for r in rows]  # pairs/s

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(COL_W, COL_W * 0.90),
                                    gridspec_kw={'wspace': 0.08})

    y = np.arange(len(labels))
    # Colour: highlight 30x WGS and Twist duplex
    colors = [BLUE, BLUE, GREEN, ORANGE, SKY, VERMILLION, PINK, GREEN]

    # RAM panel — key message: memory limitations
    max_r = max(rams)
    ax1.set_xlim(0, max_r * 1.55)
    for yy, r, c in zip(y, rams, colors):
        ax1.plot([0, r], [yy, yy], '-', color='#DDDDDD', linewidth=1.5, zorder=1)
        ax1.plot(r, yy, 'o', color=c, markersize=5.5, zorder=3)
        ax1.text(r + max_r * 0.04, yy, f'{r:.2f} GB',
                 va='center', fontsize=6, color='#444444')
    # Machine limit annotation
    ax1.axvline(x=12, color='#CC0000', linewidth=0.6, linestyle='--', alpha=0.6)
    ax1.text(12.1, len(labels) - 0.5, '12 GB\nsystem', fontsize=5.5,
             color='#CC0000', va='top')
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=6.5)
    ax1.set_xlabel('Peak RAM / GB')
    ax1.invert_yaxis()
    ax1.grid(axis='x', color='#E8E8E8', linewidth=0.6)
    ax1.grid(axis='y', visible=False)
    ax1.spines['left'].set_visible(False)
    ax1.tick_params(left=False)

    # Throughput panel
    max_p = max(pairss)
    ax2.set_xlim(0, max_p * 1.40)
    for yy, p, c in zip(y, pairss, colors):
        ax2.plot([0, p], [yy, yy], '-', color='#DDDDDD', linewidth=1.5, zorder=1)
        ax2.plot(p, yy, 's', color=c, markersize=5.5, zorder=3)
        ax2.text(p + max_p * 0.04, yy, f'{p/1000:.1f}K',
                 va='center', fontsize=6, color='#444444')
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlabel('Throughput / pairs s\u207b\u00b9')
    ax2.invert_yaxis()
    ax2.grid(axis='x', color='#E8E8E8', linewidth=0.6)
    ax2.grid(axis='y', visible=False)
    ax2.spines['left'].set_visible(False)
    ax2.tick_params(left=False)

    save(fig, 'hg38_benchmark')


# ---------------------------------------------------------------------------
# Figure 4: Memory per read pair — dot plot showing bounded growth
# Takeaway: ~0.45 KB/pair is constant across coverage depths (1 MB reference).
# ---------------------------------------------------------------------------
def fig_memory_per_pair(data):
    covs  = [1, 5, 10, 30, 50, 100, 200]
    pairs = [c * 1_000_000 / 150 for c in covs]
    kb_pp = [data['coverage_scaling'][f'cov_{c}x']['avg_peak_mem_kb'] / p
             for c, p in zip(covs, pairs)]

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.62))

    ax.scatter([f'{c}x' for c in covs], kb_pp, color=BLUE, s=28, zorder=3)

    # Plateau mean (skip 1x startup noise)
    mean_val = np.mean(kb_pp[2:])
    ax.axhline(mean_val, color='#AAAAAA', linestyle='--', linewidth=0.8)
    ax.text(len(covs) - 1.1, mean_val + 0.06,
            f'{mean_val:.2f} KB/pair', fontsize=7, color='#666666', ha='right')

    ax.set_xlabel('Coverage depth')
    ax.set_ylabel('Memory per read pair / KB')
    ax.set_ylim(bottom=0, top=max(kb_pp) * 1.4)

    save(fig, 'memory_per_pair')


# ---------------------------------------------------------------------------
# Figure 5: Throughput across main scenarios — horizontal dot chart
# Takeaway: I/O-bound configs cluster at ~14 K pairs/s; UMI lowers this.
# ---------------------------------------------------------------------------
def fig_throughput_scenarios(data):
    scenarios = [
        ('Baseline (10 MB, 30x)',   '01_baseline',        10e6, 30),
        ('+Variants (10 MB, 30x)',  '02_with_variants',   10e6, 30),
        ('High cov (1 MB, 200x)',   '04_very_high_coverage', 1e6, 200),
        ('cfDNA (1 MB, 200x)',      '07_cfdna',           1e6, 200),
        ('FFPE (10 MB, 30x)',       '08_ffpe_artifacts',  10e6, 30),
        ('Panel+UMI (1 MB, 200x)', '11_panel_umi',        1e6, 200),
    ]

    names, throughputs = [], []
    for name, key, ref_size, cov in scenarios:
        d = data['main_benchmarks'][key]
        pairs = cov * ref_size / 150
        names.append(name)
        throughputs.append(pairs / d['avg_wall_secs'] / 1000)  # K pairs/s

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.70))
    colors = [BLUE, SKY, GREEN, PINK, ORANGE, VERMILLION]
    y = np.arange(len(names))

    max_t = max(throughputs)
    ax.set_xlim(0, max_t * 1.30)
    for yy, t, c in zip(y, throughputs, colors):
        ax.plot([0, t], [yy, yy], '-', color='#DDDDDD', linewidth=1.5, zorder=1)
        ax.plot(t, yy, 'o', color=c, markersize=6, zorder=3)
        ax.text(t + max_t * 0.03, yy, f'{t:.1f}K', va='center',
                fontsize=6.5, color='#444444')

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel('Throughput / 1000 pairs s\u207b\u00b9')
    ax.invert_yaxis()
    ax.grid(axis='x', color='#E8E8E8', linewidth=0.6)
    ax.grid(axis='y', visible=False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False)

    save(fig, 'throughput_scenarios')


# ---------------------------------------------------------------------------
# Figure 6: Thread scaling — kept for appendix, shown small
# Takeaway: I/O-bound workload degrades with more threads.
# ---------------------------------------------------------------------------
def fig_thread_scaling(data):
    threads = [1, 2, 4, 6, 8, 12]
    walls = [data['thread_scaling'][f'threads_{t}']['avg_wall_secs']
             for t in threads]
    base = walls[0]
    rel  = [base / w for w in walls]

    fig, ax = plt.subplots(figsize=(COL_W * 0.65, COL_W * 0.52))

    ax.axhline(y=1.0, color='#CCCCCC', linestyle='--', linewidth=0.6)
    ax.plot(threads, rel, 'o-', color=BLUE, markersize=3.5, linewidth=0.9)

    ax.set_xlabel('Thread count', fontsize=7.5)
    ax.set_ylabel('Relative speed', fontsize=7.5)
    ax.set_xlim(0.5, 12.5)
    ax.set_ylim(0.55, 1.15)
    ax.set_xticks(threads)
    ax.tick_params(labelsize=7)

    ax.text(0.5, 0.08, 'I/O-bound workload', fontsize=6.5, color='#888888',
            transform=ax.transAxes)

    save(fig, 'thread_scaling')


if __name__ == '__main__':
    print('Generating figures...')
    data = load_json()
    hg38 = load_csv('hg38_benchmarks.csv')

    fig_coverage_scaling(data)
    fig_feature_overhead(data)
    if hg38:
        fig_hg38_benchmark(hg38)
    fig_memory_per_pair(data)
    fig_throughput_scenarios(data)
    fig_thread_scaling(data)
    print('Done.')
