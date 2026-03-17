#!/usr/bin/env python3
"""Evaluate variant caller performance against VarForge truth VCF.

Computes sensitivity, precision, and F1 per variant type and VAF tier.
Reads rtg-tools vcfeval output (TP/FP/FN counts) or falls back to direct
VCF comparison if rtg output is not available.

Usage:
  # Single scenario evaluation
  python3 evaluate.py \
      --truth results/s1/generate/truth.vcf.gz \
      --calls results/s1/call/filtered.vcf.gz \
      --rtg-eval results/s1/eval/rtg \
      --scenario "SNV baseline" \
      --outdir results/s1/eval

  # Aggregate all scenarios
  python3 evaluate.py --aggregate results/run_20260318 --outdir results/run_20260318
"""

import argparse
import csv
import gzip
import json
import os
import sys
from collections import defaultdict
from pathlib import Path


# VAF tiers used throughout the benchmark.
VAF_TIERS = [
    ("high",     0.20, 1.00),
    ("medium",   0.10, 0.20),
    ("low",      0.05, 0.10),
    ("very_low", 0.00, 0.05),
]

VARIANT_TYPES = ["SNV", "INDEL", "MNV", "ALL"]


def parse_vcf(path: str) -> list[dict]:
    """Parse a VCF (plain or gzipped), returning a list of variant dicts."""
    variants = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, vid, ref, alt = fields[:5]
            info = fields[7] if len(fields) > 7 else ""
            fmt  = fields[8] if len(fields) > 8 else ""
            sample = fields[9] if len(fields) > 9 else ""

            # Variant type
            vtype = classify_variant(ref, alt)

            # Extract VAF from INFO (VarForge truth) or FORMAT (called VCF)
            vaf = extract_vaf(info, fmt, sample)

            # Filter flag (FILTER column)
            filt = fields[6] if len(fields) > 6 else "."

            variants.append({
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "type": vtype,
                "vaf": vaf,
                "filter": filt,
                "id": f"{chrom}:{pos}:{ref}:{alt}",
            })
    return variants


def classify_variant(ref: str, alt: str) -> str:
    alts = alt.split(",")[0]  # take first alt
    if len(ref) == 1 and len(alts) == 1:
        return "SNV"
    if len(ref) != len(alts):
        return "INDEL"
    return "MNV"


def extract_vaf(info: str, fmt: str, sample: str) -> float | None:
    """Extract VAF from INFO VAF= tag (truth) or FORMAT AF field (called)."""
    # Truth VCF: INFO contains VAF=0.25
    for token in info.split(";"):
        if token.startswith("VAF="):
            try:
                return float(token[4:])
            except ValueError:
                pass
    # Called VCF: FORMAT AF or AD field
    if fmt and sample:
        keys = fmt.split(":")
        vals = sample.split(":")
        fdict = dict(zip(keys, vals))
        if "AF" in fdict:
            try:
                return float(fdict["AF"].split(",")[0])
            except (ValueError, IndexError):
                pass
        if "AD" in fdict:
            try:
                ads = [int(x) for x in fdict["AD"].split(",")]
                total = sum(ads)
                if total > 0 and len(ads) >= 2:
                    return ads[1] / total
            except (ValueError, IndexError):
                pass
    return None


def vaf_tier(vaf: float | None) -> str:
    if vaf is None:
        return "unknown"
    for name, lo, hi in VAF_TIERS:
        if lo <= vaf < hi:
            return name
    return "unknown"


def compute_metrics(tp: int, fp: int, fn: int) -> dict:
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    precision   = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (2 * precision * sensitivity / (precision + sensitivity)
          if (precision + sensitivity) > 0 else 0.0)
    return {
        "tp": tp, "fp": fp, "fn": fn,
        "sensitivity": round(sensitivity, 4),
        "precision":   round(precision, 4),
        "f1":          round(f1, 4),
    }


def evaluate_vcfs(truth_path: str, calls_path: str) -> dict:
    """Direct VCF comparison (fallback when rtg-tools output is unavailable)."""
    truth  = {v["id"]: v for v in parse_vcf(truth_path)}
    called = {v["id"]: v for v in parse_vcf(calls_path)
              if v["filter"] in ("PASS", ".")}

    # Per tier, per type
    counts = defaultdict(lambda: {"tp": 0, "fp": 0, "fn": 0})

    for vid, tv in truth.items():
        tier = vaf_tier(tv["vaf"])
        vtype = tv["type"]
        if vid in called:
            counts[(tier, vtype)]["tp"] += 1
            counts[("all", vtype)]["tp"] += 1
            counts[(tier, "ALL")]["tp"] += 1
            counts[("all", "ALL")]["tp"] += 1
        else:
            counts[(tier, vtype)]["fn"] += 1
            counts[("all", vtype)]["fn"] += 1
            counts[(tier, "ALL")]["fn"] += 1
            counts[("all", "ALL")]["fn"] += 1

    for vid, cv in called.items():
        if vid not in truth:
            tier = vaf_tier(cv["vaf"])
            vtype = cv["type"]
            counts[(tier, vtype)]["fp"] += 1
            counts[("all", vtype)]["fp"] += 1
            counts[(tier, "ALL")]["fp"] += 1
            counts[("all", "ALL")]["fp"] += 1

    results = {}
    for (tier, vtype), c in counts.items():
        key = f"{tier}/{vtype}"
        results[key] = compute_metrics(c["tp"], c["fp"], c["fn"])
    return results


def load_rtg_summary(rtg_dir: str) -> dict | None:
    """Parse rtg-tools vcfeval summary.txt into TP/FP/FN counts."""
    summary = os.path.join(rtg_dir, "summary.txt")
    if not os.path.exists(summary):
        return None
    # rtg summary.txt has lines: Threshold  True-pos-baseline  ...
    with open(summary) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    # Find the "None" threshold row (unthresholded)
    for line in lines:
        if line.startswith("None") or line.lower().startswith("none"):
            parts = line.split()
            if len(parts) >= 6:
                try:
                    return {
                        "tp": int(parts[1]),
                        "fp": int(parts[3]),
                        "fn": int(parts[4]),
                    }
                except (ValueError, IndexError):
                    pass
    return None


def write_summary(rows: list[dict], outdir: str, filename: str = "summary.tsv"):
    path = os.path.join(outdir, filename)
    os.makedirs(outdir, exist_ok=True)
    if not rows:
        print("  WARNING: no results to write.")
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Written: {path}")


def write_json(data: dict, outdir: str, filename: str):
    path = os.path.join(outdir, filename)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def run_single(args):
    print(f"Evaluating scenario: {args.scenario}")

    # Try rtg-tools output first
    rtg_result = None
    if args.rtg_eval and os.path.isdir(args.rtg_eval):
        rtg_result = load_rtg_summary(args.rtg_eval)

    if rtg_result:
        print("  Using rtg-tools vcfeval summary (all variants combined).")
        m = compute_metrics(**rtg_result)
        rows = [{"scenario": args.scenario, "tier": "all", "type": "ALL", **m}]
    else:
        print("  rtg-tools output not available; using direct VCF comparison.")
        results = evaluate_vcfs(args.truth, args.calls)
        rows = []
        for key, m in sorted(results.items()):
            tier, vtype = key.split("/")
            rows.append({"scenario": args.scenario, "tier": tier, "type": vtype, **m})

    write_summary(rows, args.outdir)
    write_json({"scenario": args.scenario, "results": rows},
               args.outdir, "results.json")

    # Print a compact table to stdout
    print(f"\n  {'Tier':<12} {'Type':<8} {'Sens':>6} {'Prec':>6} {'F1':>6} "
          f"{'TP':>6} {'FP':>6} {'FN':>6}")
    print("  " + "-" * 60)
    for row in rows:
        if row.get("tier") == "all" and row.get("type") == "ALL":
            tier, vtype = row.get("tier", "-"), row.get("type", "-")
            print(f"  {tier:<12} {vtype:<8} "
                  f"{row['sensitivity']:>6.3f} {row['precision']:>6.3f} "
                  f"{row['f1']:>6.3f} "
                  f"{row['tp']:>6} {row['fp']:>6} {row['fn']:>6}")


def run_aggregate(args):
    """Collect all scenario summary.tsv files into one."""
    all_rows = []
    for sdir in sorted(Path(args.aggregate).glob("s[0-9]*")):
        tsv = sdir / "eval" / "summary.tsv"
        if not tsv.exists():
            continue
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            all_rows.extend(reader)

    if not all_rows:
        print("No scenario results found to aggregate.")
        return

    write_summary(all_rows, args.outdir, "summary.tsv")
    print(f"\nAggregated {len(all_rows)} rows from "
          f"{len(set(r.get('scenario','?') for r in all_rows))} scenarios.")


def main():
    parser = argparse.ArgumentParser(description="Evaluate VarForge downstream VC benchmark.")
    sub = parser.add_subparsers()

    # Single-scenario mode
    parser.add_argument("--truth",     help="Path to truth VCF (.vcf.gz)")
    parser.add_argument("--calls",     help="Path to caller VCF (.vcf.gz)")
    parser.add_argument("--rtg-eval",  dest="rtg_eval", default=None,
                        help="Path to rtg vcfeval output directory")
    parser.add_argument("--scenario",  default="unknown",
                        help="Scenario name for labelling")
    parser.add_argument("--outdir",    required=True, help="Output directory")

    # Aggregation mode
    parser.add_argument("--aggregate", default=None,
                        help="Aggregate: path to run directory containing s1/ s2/ ...")

    args = parser.parse_args()

    if args.aggregate:
        run_aggregate(args)
    elif args.truth and args.calls:
        run_single(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
