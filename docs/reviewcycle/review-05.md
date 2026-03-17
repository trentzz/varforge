# Review Cycle 05: BED Region Filtering, Graph Polish, and Opus Paper Review

**Date:** 2026-03-18
**Version:** v0.1.0
**Scope:** BED region filtering implementation, graph quality improvements, unwired features audit, Opus paper review response, hg38 benchmark expansion

---

## Changes since review-04

### BED region filtering (implemented)

`regions_bed: Option<PathBuf>` was already defined in `Config` and `intersect_with_targets()` was already written — both unused.
Wired up in `src/cli/simulate.rs` after `partition_regions()`:

```rust
let regions = if let Some(ref bed_path) = cfg.regions_bed {
    let targets = parse_bed_file(bed_path)?;
    let filtered = intersect_with_targets(&regions, &targets);
    // errors if intersection is empty
    filtered
} else { regions };
```

Added:
- BED file existence check in `config::validate()`
- Error on empty intersection (catches chromosome name mismatches)
- Two integration tests: overlapping BED succeeds, non-overlapping BED fails with a clear error
- Feature spec at `docs/features/bed-region-filtering.md`

This enables high-coverage targeted panel simulation. A 1 Mbp target region at 200x now takes the same time as full chr22 at 5x, instead of being impractical.

### Micro-benchmarks relocated

`benches/varforge_benchmarks.rs` moved to `benchmarking/micro/varforge_benchmarks.rs`.
Cargo.toml updated with `path = "benchmarking/micro/varforge_benchmarks.rs"`.
`benches/` directory removed.

### Graph quality

- `lines.linewidth`: 1.8 → 1.0 (thinner, more professional)
- `savefig.dpi`: 300 → 600 (print-quality raster output)
- Thread scaling figure: 70×55mm → 55×45mm (smaller, less prominent)
- Fixed `grid.axis: "y"` → `axes.grid.axis: "y"` (invalid rcParam, was silently ignored)

### Paper corrections from Opus review

**Coverage scaling text fixed**: The subsection claimed "~13,000 pairs/s" for the 1 MB reference benchmark. The actual data shows ~6,700 pairs/s for 1 MB reference. The 14,300 pairs/s figure is from the 10 MB reference where more parallel chunks improve thread utilisation. Fixed text to distinguish the two.

**Appendix table updated**: The coverage scaling appendix table had stale data from the old binary run (28,140 pairs/s at 30x, 362 MB RAM). Updated with fresh CSV data (6,673 pairs/s at 30x, 91 MB RAM).

**Figure captions corrected**:
- `feature_overhead` caption: "+64%" → "+143%" (was stale)
- `memory_per_pair` caption: "~1.7 KB/pair" → "~0.45 KB/pair" (was stale)

**Abstract rewritten**: Removed triple "first" claims. Now states capabilities and what makes VarForge distinct without using "first" three times in one sentence.

**Introduction**: Broke the 8-item run-on sentence in paragraph 2 into 6 short sentences.

**Conclusion rewritten**: Was a near-verbatim copy of the abstract. Now states the key validation gap (output correctness not yet independently verified) as the main limitation and next step.

**Bibliography DOIs**: Four entries had obviously placeholder DOIs (`*-00000-0`, `btad000`). Removed the fake DOIs; added "DOI pending." notes.

**Thread scaling section**: Reduced figure size and textual prominence. Added pointer to Discussion for when threading helps. Discussion section expanded with a full analysis of when Rayon parallelism provides benefit (compute-heavy workloads like UMI simulation) vs when it hurts (I/O-bound workloads where the single writer thread is the bottleneck).

### hg38 benchmark expansion

Added two new configs:
- `07_wgs_30x.yaml`: full chr22 at 30x (realistic WGS depth)
- `08_panel_bed_200x.yaml`: 1 Mbp target at 200x with BED filtering (demonstrates the new feature)

Results below (pending completion of 07 and 08):

| Config | Wall (s) | Peak RAM (GB) | Read pairs | Pairs/s |
|--------|----------|---------------|------------|---------|
| WGS baseline (5x) | 99.5 | 0.72 | 847K | 8,510 |
| WGS + variants (5x) | 103.6 | 0.72 | 847K | 8,178 |
| Panel UMI simplex (5x) | 234.1 | 2.26 | 2,540K | 10,851 |
| Twist duplex (5x) | 299.2 | 3.00 | 3,389K | 11,328 |
| cfDNA liq. biopsy (5x) | 96.1 | 0.85 | 847K | 8,813 |
| FFPE tumour (5x) | 120.6 | 0.91 | 1,059K | 8,782 |
| WGS baseline (30x) | 595.3 | 4.21 | 5,082K | 8,537 |
| Panel BED 200x (1 Mbp) | 85.5 | 0.56 | 667K | 7,796 |

---

## Opus paper review findings and responses

### Priority 1: Output correctness validation (not yet done)

The paper makes three headline contribution claims for which there is no experimental validation:
- VAF accuracy: no chi-squared or KS test against the expected binomial
- cfDNA fragment distribution: no extraction of insert sizes from BAM output
- UMI family size distribution: no comparison to configured log-normal

This is the largest credibility gap. Added to `docs/planning/` as the primary outstanding task for the next review cycle.

### Priority 1: Downstream variant caller validation (not yet done)

The paper lists Mutect2, VarDict, and Strelka2 as target consumers but never runs any of them on VarForge output. A single Mutect2 run on a VarForge dataset (known variants, truth VCF) would be a high-impact addition. Carried to next cycle.

### Priority 2: Tighten "first" claims (done)

Abstract no longer uses "first" three times. Introduction still lists the UMI and cfDNA gaps but no longer asserts "first" without qualification. The feature comparison table is the evidence.

### Priority 2: Merge Evaluation into Results (deferred)

The Evaluation section currently re-evaluates claims that are already presented in Results. This creates redundancy. Deferred to next review cycle — a structural restructure requires care and is better done alongside the output validation work.

### Priority 2: Thread scaling for compute-bound workload (not yet done)

Currently only the I/O-bound workload (10 MB, 30x) is benchmarked across thread counts. The UMI workload should also be tested. Carried to next cycle.

---

## Unwired features audit

Found five config fields or internal components that are parsed/defined but never used:

| Item | Location | Status |
|------|----------|--------|
| `UmiConfig.inline` | `config.rs:201` | Parsed, never read. Inline UMI embedding in read sequence. |
| `FragmentModel::Custom` | `config.rs:120` | Falls through to Normal. Silent no-op. |
| `AppliedVariant.actual_alt_count/total_count` | `engine.rs:39-42` | Computed, never written to output. |
| `PcrFamilySizeSampler` | `fragment.rs:88-109` | Log-normal family size sampler, never instantiated. |
| `ClonalTree` | `clonal_tree.rs` | Formal tree structure unused in production code. |

`FragmentModel::Custom` is the most dangerous: a user can set `model: custom` and get Normal behaviour with no warning. This should either be implemented or raise a config error. Added to the planning document for the next cycle.

---

## Outstanding issues (carried from review-04)

1. **VAF accuracy validation** (Priority 1): chi-squared test against Binomial(D, VAF). Not done.
2. **cfDNA fragment distribution validation**: extract fragment lengths, compare to mixture model. Not done.
3. **Compute-bound thread scaling test**: UMI simulation at multiple thread counts. Not done.
4. **`umi.spacer_length` config field**: needed for exact Twist protocol simulation. Not implemented.
5. **Downstream variant caller validation**: run Mutect2 or VarDict on VarForge output. Not done.
6. **`FragmentModel::Custom` warning**: should error or implement. Currently silent no-op.
