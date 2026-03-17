# Paper Intent: VarForge v0.1

This document describes the goals, scope, and editorial direction of the VarForge paper.
Update it when the target audience, claims, or scope change.

---

## What this paper is

A tool paper introducing VarForge: a single Rust binary that generates synthetic
cancer sequencing data driven by a YAML configuration. The paper is written for
bioinformatics engineers and researchers who benchmark variant callers, UMI
deduplication tools, and liquid biopsy algorithms.

---

## Primary claims (must be supported by evidence in the paper)

1. VarForge unifies read generation, variant spike-in, UMI/duplex simulation,
   cfDNA fragmentation, and library artefact injection in one tool. No other
   single tool covers this feature space.

2. Stochastic VAF (per-read Bernoulli sampling) produces more realistic coverage
   at variant sites than deterministic spike-in. This matters for sensitivity
   estimates at low VAF.

3. The streaming Rust architecture bounds memory by concurrent batch size, not
   dataset size. This makes full-chromosome and WGS simulations feasible on a
   standard workstation.

4. BED region filtering enables clinical-depth panel simulation (200x, 1 Mbp
   target on chr22) in under 90 seconds.

5. VarForge supports the Twist Biosciences duplex UMI protocol, enabling
   generation of realistic test datasets for liquid biopsy variant detection
   pipelines.

---

## What this paper is NOT

- A validation paper. Output correctness (VAF accuracy, cfDNA fragment
  distributions, UMI family size distributions) is not yet independently
  validated. This is stated openly as the main limitation and next step.
- A comparison paper. Direct throughput comparison with other tools is
  complicated by feature asymmetry. We note what comparable tools do, but
  do not claim "N times faster than X".
- A bioinformatics methods paper. We do not invent new algorithms. We engineer
  a practical tool that combines known models into a unified, usable system.

---

## Target audience

Primary: bioinformatics engineers who need ground-truth data to test tools
(variant callers, UMI deduplication, cfDNA analysis) and who have been stitching
together ART + BAMSurgeon + custom scripts.

Secondary: researchers developing new cancer sequencing methods who need
configurable reference datasets.

---

## Format

The canonical format is the two-column paper (`main_2col.tex`). Use this as the
primary working format. The single-column normal and big formats exist for
readability during review. The 2-pager is for conference abstract submissions.

---

## Current paper structure (main_2col.tex)

1. Abstract
2. Introduction: the gap (no single tool does UMI + cfDNA + variants), the
   contributions
3. Background: cfDNA biology, UMI/duplex sequencing, VAF statistics
4. Related Work: feature comparison table, two clear gaps
5. Method: architecture, read generation, VAF sampling, variant spike-in,
   tumour model, UMI/duplex, cfDNA, artefacts, implementation
6. Results: coverage scaling, feature overhead, throughput scenarios,
   hg38 benchmark, memory efficiency
7. Evaluation: claim-by-claim check against benchmark evidence
8. Discussion: UMI gap, cfDNA gap, stochastic VAF, performance for CI/CD
9. Future Work: VAF validation, cfDNA validation, downstream caller benchmarking
10. Conclusion
11. References
12. Appendix: extended feature table, coverage scaling data

---

## Figures

Priority figures (in the paper):
- `coverage_scaling`: two panels (wall time, RAM) vs coverage depth
- `feature_overhead`: two-panel dot chart (time overhead %, RAM per feature)
- `hg38_benchmark`: two-panel dot chart (peak RAM, throughput) per config
- `memory_per_pair`: dot plot, constant ~0.45 KB/pair confirms bounded memory
- `throughput_scenarios`: horizontal dot chart comparing scenario throughput

Appendix / lower priority:
- `thread_scaling`: I/O-bound degradation with threads (kept small, honest result)

The pipeline TikZ diagram stays in Method as the main architecture figure.

---

## Key editorial rules

- Thread scaling is an honest limitation, not a headline result. Keep the
  thread scaling figure small and in the Evaluation or appendix, not Results.
- RAM usage and memory limitations are a primary concern for real-world users.
  Show RAM prominently in the hg38 benchmark figure, with a reference line
  showing the machine limit.
- The Twist duplex section should be short (one paragraph + one figure or table).
  It is an example use case, not a separate contribution.
- Do not claim "first" without evidence. The feature comparison table is the
  evidence. Let it speak.
- Output correctness validation is the stated gap. Name it clearly and early.
  Hiding limitations damages credibility more than acknowledging them.

---

## Outstanding work (update this list as items are completed)

- [ ] VAF accuracy chi-squared test
- [ ] cfDNA fragment distribution extraction from BAM output
- [ ] Downstream variant caller benchmarking (Mutect2 or VarDict)
- [ ] Compute-bound thread scaling benchmark (UMI simulation at 1/2/4/6/12 threads)
- [ ] FragmentModel::Custom implementation (currently warns and falls back to Normal)
- [ ] umi.spacer_length config field for exact Twist protocol simulation
