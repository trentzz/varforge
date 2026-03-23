# EPIC-TWIST-DUPLEX: Twist Biosciences UMI Duplex Sequencing

## Goal

Make VarForge generate publication-quality synthetic data that matches the characteristics
of Twist Biosciences hybrid-capture duplex UMI sequencing. The output must be realistic
enough to benchmark duplex-aware variant callers (fgbio, HUMID, Mutect2 with duplex mode)
on SNV, indel, and SV mutations at ultra-low VAFs.

## Motivation

Twist duplex sequencing is the dominant platform for ultra-sensitive ctDNA detection
(0.001%–0.1% VAF). Benchmarking callers on it requires:

1. Correct AB/BA strand family structure with independent sequencing errors per strand
2. Accurate duplex conversion rate (fraction of molecules with both strands observed)
3. Strand concordance for true positives (alt must appear in both AB and BA strands of the
   same molecule — the defining signal of a true variant in duplex calling)
4. Twist-specific capture characteristics: high uniformity (CV < 0.2), >95% on-target,
   panel-scale target regions
5. SVs represented as breakpoint-spanning reads in both strand orientations

## What Success Looks Like

A successful Twist duplex simulation produces output where:

| Metric | Target | Rationale |
|--------|--------|-----------|
| Duplex conversion rate | ≥ 85% | Fraction of unique molecules with both AB and BA families observed. Twist specs target 85–95%. |
| Strand concordance (true variants) | ≥ 99% | For each spiked-in variant, ≥ 99% of alt-carrying molecules have the alt in both AB and BA families. |
| VAF accuracy | Within 2× for VAF ≥ 0.001% | At 2000× coverage, a 0.01% VAF variant expects ~0.4 duplex molecules; Poisson variance is high but the simulation should be unbiased. |
| On-target rate | ≥ 95% | Fraction of read families overlapping the capture BED. |
| Coverage CV | ≤ 0.25 | Coefficient of variation across targets. Twist panels are among the most uniform available. |
| Family size distribution | LogNormal(μ=3.5, σ=1.5) | Typical Twist duplex family size distribution from published benchmarks. |
| AB/BA family size ratio | 0.8–1.2 | AB and BA families from the same molecule should have similar sizes. |
| SV duplex support | ≥ 80% of alt molecules have both strands at breakpoint | For deletions, duplications, and translocations. |

## Scope

- In scope: Twist-specific preset, duplex concordance tracking, duplex conversion rate,
  ultra-low VAF sampling correctness, SV duplex support, sim_report.json metrics.
- Out of scope: duplex consensus calling (that is the caller's job), methylation, base
  modification tags (MM/ML).

## Tasks

- T105: Add dedicated `twist` preset with Twist-specific parameters
- T106: Wire duplex_alt_count (track alt concordance across AB and BA strands per molecule)
- T107: Add strand concordance metric to sim_report.json
- T108: Add duplex conversion rate to sim_report.json
- T109: Validate SV reads appear in both AB and BA strand families at breakpoints
- T110: Verify ultra-low VAF Poisson sampling is unbiased at VAF < 0.01%
- T111: Add Twist capture uniformity parameters (CV target, on-target fraction) to preset
- T112: Write a benchmark scenario config and document the success criteria from this epic
