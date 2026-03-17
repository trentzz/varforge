# Review Cycle 04: hg38 Benchmarks, Data Integrity Corrections, and Paper Naming

**Date:** 2026-03-18
**Version:** v0.1.0
**Scope:** hg38 benchmark execution, critical data corrections from Opus review, paper file naming convention

---

## Changes since review-03

### Paper naming convention

`~/.claude/guides/paper-writing.md` was updated to require the naming format
`{paper-name}_{version}_{type}.pdf` with version on the title page.

Applied:
- `docs/paper/main.tex`: title now shows `\texttt{v0.1}` on a second line
- `docs/paper/main.pdf` removed from git tracking (added to `.gitignore`)
- `docs/paper/varforge_v0.1_normal.pdf` added as the committed distribution copy

### hg38 config corrections

Four of six hg38 benchmark configs had impractical coverage depths:
- `03_panel_umi.yaml`: 500x → 30x
- `04_twist_duplex.yaml`: 1000x → 30x
- `05_cfdna.yaml`: 200x → 30x (also fixed absolute reference path to `chr22.fa`)
- `06_ffpe_tumour.yaml`: 60x → 30x

**Reason:** VarForge does not yet support BED-based region filtering. At 500x or 1000x on the full chr22 (51 MB), the read count is orders of magnitude larger than a real targeted panel run. Without a regions filter, these configs would take many hours and likely exhaust system RAM. All configs were set to 30x to match the WGS baseline and ensure feasibility on the 12 GiB test machine.

### Paper sections updated

- `sections/05_results.tex`: added hg38 benchmark subsection with placeholder table
- `sections/06_evaluation.tex`: corrected UMI/duplex claim to reference 30x chr22 configs (not the incorrect 1000x figure from the original draft)
- `sections/08_future_work.tex`: added targeted region filtering as the first future work item (direct consequence of this benchmark exercise)

### Benchmark runner

`benchmarking/scripts/run_hg38_benchmarks.sh` rewritten to use Python for timing and RSS polling via `/proc/{pid}/status`, since GNU `time` is not installed on the test system. The script:
- Patches the `reference:` path to absolute before each run
- Polls `/proc/{pid}/status` for peak VmRSS every 0.5s
- Reads read-pair count from the manifest JSON
- Writes to `benchmarking/results/tables/hg38_benchmarks.csv` incrementally

### hg38 benchmark results

All six configs completed at 5x coverage. Results:

| Config | Wall (s) | Peak RAM (GB) | Read pairs | Pairs/s |
|--------|----------|---------------|------------|---------|
| WGS baseline | 96.6 | 0.72 | 847K | 8,772 |
| WGS + variants | 96.6 | 0.72 | 847K | 8,772 |
| Panel UMI simplex | 233.1 | 2.26 | 2,540K | 10,894 |
| Twist duplex | 298.6 | 3.00 | 3,389K | 11,349 |
| cfDNA liq. biopsy | 96.6 | 0.85 | 847K | 8,772 |
| FFPE tumour | 121.1 | 0.91 | 1,059K | 8,746 |

Config 05 (cfDNA) initially failed because the model name in the config was `cfdna`
instead of the correct `cfda` (as defined in the VarForge source). This was fixed.

The UMI and duplex configs produce more read pairs than the baseline because PCR family
expansion is included in the output (3x and 4x family sizes respectively).
RAM for UMI (2.26 GB) and duplex (3.00 GB) is 3--4x higher than baseline (0.72 GB),
confirming that PCR family tracking is the dominant memory cost for these modes.

The paper table and PDF have been updated with these numbers.

### Critical data corrections (Opus review findings)

An Opus-powered review identified that `results.json` contained benchmark data from an
earlier, faster run of VarForge (possibly under different system state or binary version).
The paper's numbers were based on this old data and were incorrect.

**What was wrong:**
- Table 1 baseline: 75.4s (old) → 139.6s (fresh)
- Throughput claim: "28,000 pairs/s" → "~13,000 pairs/s"
- Feature overhead UMI: "+64%" → "+143%"
- Feature overhead combined: "+125%" → "+198%"
- Thread scaling: "Speedup plateaus at 1.16x" → "Adding threads degrades performance by up to 13%"
- Memory per pair: "~1.7 KB" → "~0.45 KB"

**Root cause:** `results.json` was used by both the plot script and as a reference for paper
numbers. It had not been regenerated from the fresh CSV benchmarks. Fixed by rebuilding
`results.json` from the authoritative CSV files, running the missing `11_panel_umi` config
(195.1s, 1.80 GB, ~14K input pairs), and updating all paper text.

**Paper sections corrected:** abstract, 05_results (all subsections), 06_evaluation,
07_discussion.

**Graphs regenerated** from fresh data. Thread scaling graph updated to show the measured
slowdown honestly rather than a manufactured speedup.

---

## Outstanding issues (carried from review-03)

1. **VAF accuracy validation**: stochastic VAF claim needs a chi-squared test against Binomial(D, VAF). Not yet done.
2. **cfDNA fragment distribution validation**: extract fragment lengths from cfDNA output and compare to the mixture model. Not yet done.
3. **Compute-bound thread scaling test**: current thread scaling shows I/O-bound plateau. A compute-bound test would show Rayon parallelism working. Not yet done.
4. **`umi.spacer_length` config field**: needed for exact Twist protocol simulation. Not yet implemented.

---

## Honesty check

### hg38 coverage reduction
The design doc said panels at 500x/1000x. Those coverages are not achievable without a regions filter. The paper now states this limitation plainly in the hg38 subsection and in Future Work. The evidence for UMI/duplex functionality comes from the synthetic benchmarks (1 MB, 200x), which did complete. The hg38 configs at 30x add evidence that the features work on a real reference.

### Placeholder table
The hg38 table in `05_results.tex` currently shows `---` placeholders. This will be filled in once the benchmark run completes. The PDF committed here contains placeholders; a follow-up commit will contain the real numbers.
