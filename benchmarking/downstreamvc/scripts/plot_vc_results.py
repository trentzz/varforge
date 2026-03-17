#!/usr/bin/env python3
"""Plot downstream variant caller benchmark results.

Reads results/run_*/summary.tsv and produces publication figures:
  - sensitivity_by_vaf.pdf: sensitivity per VAF tier and scenario
  - precision_recall.pdf: precision-recall by scenario
  - fp_by_scenario.pdf: false positive counts (artefact impact)

Run after run_benchmark.sh completes.

Usage:
  python3 plot_vc_results.py --results results/run_20260318/summary.tsv
"""

import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Match project graph style (see ~/.claude/guides/graph-style.md)
PALETTE = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # green
    "#CC79A7",  # pink
    "#56B4E9",  # sky
    "#D55E00",  # vermillion
]
BLUE, ORANGE, GREEN, PINK, SKY, VERMILLION = PALETTE
COL_W = 84 / 25.4  # 84 mm in inches

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
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

VAF_TIER_ORDER  = ["high", "medium", "low", "very_low"]
VAF_TIER_LABELS = [">20%", "10–20%", "5–10%", "1–5%"]

SCENARIO_ORDER = [
    "SNV baseline",
    "Low-VAF ctDNA",
    "FFPE artefacts",
    "Panel + UMI",
]


def load_results(path: str) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(fig, outdir: str, name: str):
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, name + ".pdf"))
    fig.savefig(os.path.join(outdir, name + ".png"))
    plt.close(fig)
    print(f"  {name}.pdf + .png")


def fig_sensitivity_by_vaf(rows: list[dict], outdir: str):
    """Grouped dot chart: sensitivity per VAF tier, one group per scenario.

    Takeaway: sensitivity degrades at lower VAF; how fast depends on the scenario.
    """
    # Only SNV rows, all VAF tiers
    snv = [r for r in rows if r.get("type") == "SNV" and r.get("tier") in VAF_TIER_ORDER]
    if not snv:
        snv = [r for r in rows if r.get("type") == "ALL" and r.get("tier") in VAF_TIER_ORDER]
    if not snv:
        print("  WARNING: no per-tier SNV data available. Skipping sensitivity_by_vaf.")
        return

    scenarios = [s for s in SCENARIO_ORDER if any(r["scenario"] == s for r in snv)]
    colors = [BLUE, ORANGE, GREEN, PINK][:len(scenarios)]

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.72))

    n_tiers = len(VAF_TIER_ORDER)
    n_scen  = len(scenarios)
    width   = 0.15
    x       = np.arange(n_tiers)

    for si, (scen, color) in enumerate(zip(scenarios, colors)):
        offsets = x + (si - n_scen / 2 + 0.5) * width
        sens = []
        for tier in VAF_TIER_ORDER:
            match = [r for r in snv if r["scenario"] == scen and r["tier"] == tier]
            sens.append(float(match[0]["sensitivity"]) if match else 0.0)
        ax.plot(offsets, sens, "o-", color=color, markersize=4,
                linewidth=0.8, label=scen)

    ax.set_xticks(x)
    ax.set_xticklabels(VAF_TIER_LABELS, fontsize=7.5)
    ax.set_xlabel("VAF tier")
    ax.set_ylabel("Sensitivity")
    ax.set_ylim(0, 1.08)
    ax.axhline(1.0, color="#DDDDDD", linewidth=0.6, linestyle="--")
    ax.legend(loc="lower left", fontsize=6.5)

    save(fig, outdir, "sensitivity_by_vaf")


def fig_precision_recall(rows: list[dict], outdir: str):
    """Precision-recall scatter by scenario (all-tier ALL variant type).

    Takeaway: where each scenario sits in precision-recall space.
    """
    summary = [r for r in rows if r.get("tier") == "all" and r.get("type") == "ALL"]
    if not summary:
        print("  WARNING: no all-tier summary data. Skipping precision_recall.")
        return

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.75))

    markers = ["o", "s", "^", "D"]
    for row, color, marker in zip(summary, PALETTE, markers):
        prec = float(row["precision"])
        sens = float(row["sensitivity"])
        ax.plot(sens, prec, marker=marker, color=color, markersize=7, zorder=3)
        # Put label above for low-precision points to avoid overlap
        if prec < 0.1:
            xytext = (5, 12)
        else:
            xytext = (5, 3)
        ax.annotate(row["scenario"],
                    xy=(sens, prec),
                    xytext=xytext,
                    textcoords="offset points",
                    fontsize=6.5,
                    color="#444444")

    ax.set_xlabel("Sensitivity (recall)")
    ax.set_ylabel("Precision")
    ax.set_xlim(-0.05, 1.10)
    ax.set_ylim(-0.05, 1.10)
    ax.axhline(1.0, color="#EEEEEE", linewidth=0.5)
    ax.axvline(1.0, color="#EEEEEE", linewidth=0.5)
    ax.grid(axis="both", color="#E8E8E8", linewidth=0.6)

    save(fig, outdir, "precision_recall")


def fig_fp_by_scenario(rows: list[dict], outdir: str):
    """Horizontal dot chart: false positive counts per scenario (log scale).

    Takeaway: FFPE artefacts overwhelm stock Mutect2 with >100k FPs; baseline and panel are clean.
    """
    summary = [r for r in rows if r.get("tier") == "all" and r.get("type") == "ALL"]
    if not summary:
        print("  WARNING: no all-tier summary data. Skipping fp_by_scenario.")
        return

    scenarios = [r["scenario"] for r in summary]
    fps       = [int(r["fp"]) for r in summary]
    # Replace 0 with 0.5 for log scale display
    fps_plot  = [max(fp, 0.5) for fp in fps]
    colors    = [PALETTE[i % len(PALETTE)] for i in range(len(scenarios))]

    fig, ax = plt.subplots(figsize=(COL_W, COL_W * 0.65))
    y = np.arange(len(scenarios))

    for yy, fp, fp_p, c in zip(y, fps, fps_plot, colors):
        ax.plot([0.5, fp_p], [yy, yy], "-", color="#DDDDDD", linewidth=1.5, zorder=1)
        ax.plot(fp_p, yy, "o", color=c, markersize=6, zorder=3)
        label = str(fp) if fp > 0 else "0"
        ax.text(fp_p * 1.8, yy, label, va="center", fontsize=6.5, color="#444444")

    ax.set_xscale("log")
    ax.set_xlim(0.4, max(fps_plot) * 12)
    ax.set_yticks(y)
    ax.set_yticklabels(scenarios, fontsize=7)
    ax.set_xlabel("False positives (log scale)")
    ax.invert_yaxis()
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.6)
    ax.grid(axis="y", visible=False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)

    save(fig, outdir, "fp_by_scenario")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True,
                        help="Path to summary.tsv from run_benchmark.sh")
    parser.add_argument("--outdir", default=None,
                        help="Output directory for figures (default: same dir as results)")
    args = parser.parse_args()

    if not os.path.exists(args.results):
        print(f"ERROR: results file not found: {args.results}")
        sys.exit(1)

    outdir = args.outdir or os.path.dirname(args.results)
    rows   = load_results(args.results)
    print(f"Loaded {len(rows)} rows from {args.results}")

    print("Generating figures...")
    fig_sensitivity_by_vaf(rows, outdir)
    fig_precision_recall(rows, outdir)
    fig_fp_by_scenario(rows, outdir)
    print("Done.")


if __name__ == "__main__":
    main()
