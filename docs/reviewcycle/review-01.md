# VarForge Code Review #1

**Date:** 2026-03-17
**Reviewer:** Senior Rust / Bioinformatics Engineer
**Scope:** All source files under `src/`
**Commit:** current working tree (pre-commit state)

---

## Executive Summary

VarForge is a well-structured Rust project with clean module boundaries, good test coverage for core algorithms, and correct domain modeling for cancer sequencing simulation. The codebase demonstrates solid engineering on the fundamentals: deterministic seeding, stochastic VAF sampling, proper coordinate handling, and extensible configuration.

The issues found below are predominantly in three areas: **(1)** performance-critical allocations in the simulation hot path, **(2)** correctness edge cases in genomic coordinate handling and variant spike-in, and **(3)** architecture improvements for long-term maintainability.

**Issue counts by severity:**
- Critical: 3
- High: 10
- Medium: 18
- Low: 8

---

## Critical Issues

### C-1. `apply_variant_to_seq` clones the entire fragment sequence per variant

**File:** `src/core/engine.rs`, line ~637
**Category:** Performance
**Severity:** Critical

The function creates a `Read` by cloning the full fragment sequence (`seq.clone()`) and a quality placeholder vec, applies the variant, then copies the modified sequence back. This happens inside the inner variant loop of the hot path (`generate_region`), meaning for a region with N read pairs and M variants, we do up to N*M full fragment-length allocations and copies.

```rust
let mut temp_read = Read::new(seq.clone(), qual_placeholder);
// ...
*seq = temp_read.seq;
```

**Suggested fix:** Refactor `spike_snv`, `spike_mnv`, and `spike_indel` to operate directly on `&mut [u8]` slices rather than requiring a `Read` wrapper. This eliminates two allocations and two copies per variant per read. A minimal change would be to add `_raw` variants of these functions that take `&mut Vec<u8>` directly.

### C-2. `ParametricQualityModel` is reconstructed inside the per-read loop

**File:** `src/core/engine.rs`, lines ~371-374
**Category:** Performance
**Severity:** Critical

Inside `generate_region`, when no empirical quality model is loaded, a new `ParametricQualityModel` is constructed for every read pair:

```rust
let parametric = ParametricQualityModel::new(
    self.config.quality.mean_quality,
    self.config.quality.tail_decay,
);
```

While the struct itself is small, the `Normal::new()` distribution inside is called once per construction. This should be hoisted above the read-generation loop and reused.

**Suggested fix:** Build the parametric model once before the `while pair_idx < n_pairs` loop. Store it alongside the empirical model in an enum to avoid the branch inside the loop as well.

### C-3. `ReferenceGenome::sequence()` called per fragment in hot loop with Mutex contention

**File:** `src/core/engine.rs`, line ~294; `src/io/reference.rs`, line ~144
**Category:** Performance
**Severity:** Critical

Each fragment generation calls `self.reference.sequence()`, which locks a `Mutex<IndexedReader>`. Although the `Clone` impl reopens the file to get a separate handle for parallel workers, the `generate_region` method takes `&mut self`, which holds a single `Arc<ReferenceGenome>`. The Mutex is locked and unlocked once per fragment, which involves a syscall overhead on contended paths.

More critically, each call does a fresh file seek and allocation for the sequence `Vec<u8>`. For a 1 Mbp region at 30x coverage with 150bp reads, that is ~100,000 individual file seeks.

**Suggested fix:** Pre-fetch the entire region's reference sequence once at the top of `generate_region` and slice into it for each fragment. This reduces file I/O from O(n_pairs) to O(1) per region and eliminates the Mutex entirely from the hot path:

```rust
let region_seq = self.reference.sequence(region)?;
// then for each fragment:
let frag_seq = region_seq[frag_offset..frag_offset + frag_len].to_vec();
```

---

## High Severity Issues

### H-1. `intersect_with_targets` is O(regions * targets) - quadratic for panel sequencing

**File:** `src/core/coverage.rs`, lines 36-50
**Category:** Performance
**Severity:** High

This function performs a nested loop over all regions and all targets. For a typical exome capture with ~200,000 target intervals and ~3,000 genome regions, this is 600 million iterations. The function should use an interval tree or at minimum sort targets by chromosome and binary-search.

**Suggested fix:** Use `HashMap<String, Vec<Region>>` keyed by chromosome, sort each chromosome's targets by start position, and use binary search to find overlaps.

### H-2. `CaptureModel::coverage_multiplier_at` is O(n_targets) per position

**File:** `src/core/capture.rs`, lines 118-133
**Category:** Performance
**Severity:** High

Linear scan through all target regions for every position lookup. For WES with 200k targets, this is called per-fragment in the hot loop.

**Suggested fix:** Build an interval index (sorted + binary search, or an interval tree) at construction time.

### H-3. cfDNA fragment sampler ignores config mean/sd parameters

**File:** `src/core/engine.rs`, lines 232-242
**Category:** Correctness
**Severity:** High

When `FragmentModel::Cfda` is selected, the fragment sampler is constructed with hardcoded parameters `(167.0, 334.0, 0.85, 0.0)` rather than using `self.config.fragment.mean` and `self.config.fragment.sd`:

```rust
FragmentModel::Cfda => Sampler::Cfdna(CfdnaFragmentSampler::new(
    167.0, 334.0, 0.85, 0.0,
)),
```

Users who set `fragment.mean: 155` in their YAML config (e.g., to model a different nucleosomal protection pattern) would silently get the default 167bp peak.

**Suggested fix:** Pass `self.config.fragment.mean`, compute di-peak as `2 * mean`, and derive ctdna_fraction from purity config.

### H-4. `find_cn_region` is O(n) linear scan per variant per fragment

**File:** `src/variants/cnv.rs`, lines 118-124; called from `src/core/engine.rs`, line 192
**Category:** Performance
**Severity:** High

Called once per region to adjust coverage, which is acceptable. However, the same pattern is used for per-variant lookup. For genomes with many CNV regions (e.g., highly rearranged cancers), this should use sorted lookup.

### H-5. `format!` allocation for every read name in hot loop

**File:** `src/core/engine.rs`, line ~383
**Category:** Performance
**Severity:** High

```rust
name: format!("{}_{}:{}", region.chrom, frag_start, pair_idx),
```

This allocates a new `String` for each of potentially millions of read pairs. The read name also clones `region.chrom` via the format macro.

**Suggested fix:** Pre-compute the chromosome prefix once, and use a reusable buffer with `write!` to avoid repeated allocations:

```rust
let mut name_buf = String::with_capacity(64);
// in loop:
name_buf.clear();
write!(&mut name_buf, "{}_{}:{}", chrom_prefix, frag_start, pair_idx).unwrap();
```

### H-6. `region.chrom.clone()` on every fragment and read pair

**File:** `src/core/engine.rs`, lines ~293, ~388
**Category:** Performance
**Severity:** High

Each `ReadPair` stores `chrom: String` and each fragment `Region::new` call clones the chromosome name. For millions of reads in a single region, these are all identical strings.

**Suggested fix:** Use `Arc<str>` or a reference-counted chromosome string. Alternatively, since all reads in a region share the same chrom, store it once in the output batch and remove `chrom` from `ReadPair`.

### H-7. `write_fastq_record` allocates a Vec for quality encoding on every record

**File:** `src/io/fastq.rs`, line 99
**Category:** Performance
**Severity:** High

```rust
let encoded: Vec<u8> = qual.iter().map(|&q| q + 33).collect();
```

For millions of reads, this creates a temporary Vec per record. Use a reusable buffer, or write the bytes one at a time (the write is already buffered by GzEncoder + BufWriter).

**Suggested fix:** Write the encoded quality bytes directly:
```rust
for &q in qual {
    writer.write_all(&[q + 33])?;
}
```

### H-8. Purity combination logic bug in edit subcommand config mode

**File:** `src/cli/edit.rs`, line 172
**Category:** Correctness
**Severity:** High

```rust
let purity = cfg.purity.min(purity_override);
```

This takes the *minimum* of the config purity and the override purity. The intent appears to be that the CLI override should take precedence, but `min` means the CLI default of `1.0` will never raise a config-specified purity. If the config says `0.3` and CLI says `1.0` (default), the result is `0.3` -- correct. But if config says `0.9` and user passes `--purity 0.5`, the result is `0.5` -- correct by accident. If config says `0.3` and user passes `--purity 0.7`, the result is `0.3` -- wrong, user intent ignored.

**Suggested fix:** Use a proper override pattern: if the CLI purity differs from its default value (`1.0`), use the CLI value; otherwise use the config value.

### H-9. `Read::new` panics with `assert_eq!` in non-test code

**File:** `src/core/types.rs`, line 47
**Category:** Correctness / Safety
**Severity:** High

```rust
assert_eq!(seq.len(), qual.len(), "sequence and quality length mismatch");
```

This is a debug assertion that panics in production. If any code path produces a mismatched seq/qual (e.g., an indel spike-in edge case), the entire simulation crashes rather than reporting an error.

**Suggested fix:** Return `Result<Read>` or use `debug_assert_eq!` and clamp quality length at runtime.

### H-10. `NormalFragmentSampler::new` panics with `expect` on invalid parameters

**File:** `src/core/fragment.rs`, line 18
**Category:** Correctness / Safety
**Severity:** High

```rust
dist: Normal::new(mean, sd).expect("invalid normal parameters"),
```

If a user provides `sd: 0` or `sd: -1` in their config, this panics. Same issue in `CfdnaFragmentSampler::new` (lines 50-51) and `PcrFamilySizeSampler::new` (line 101).

**Suggested fix:** Validate these parameters in `config::validate()` and return `Result` from constructors.

---

## Medium Severity Issues

### M-1. `serde_yaml` is deprecated, use `serde_yml`

**File:** `Cargo.toml`, line 48
**Category:** Rust Code Quality
**Severity:** Medium

The `serde_yaml` crate (v0.9) is unmaintained. The community successor is `serde_yml`. While functionally equivalent today, this should be migrated before it becomes a security concern.

### M-2. Duplicate `GcBiasConfig` type definitions

**File:** `src/io/config.rs` (lines 236-257) vs `src/core/gc_bias.rs` (lines 31-47)
**Category:** Architecture
**Severity:** Medium

There are two separate `GcBiasConfig` structs: one in `io::config` (with `model: String`) and one in `core::gc_bias` (with `model: GcBiasModelKind` enum). The engine's `build_gc_bias_model` manually translates between them. This duplication invites drift.

**Suggested fix:** Use a single `GcBiasConfig` with the enum type, implementing custom serde for the string-to-enum conversion.

### M-3. `Region::len()` can underflow when `start > end`

**File:** `src/core/types.rs`, line 18-20
**Category:** Correctness
**Severity:** Medium

```rust
pub fn len(&self) -> u64 {
    self.end - self.start
}
```

If `start > end` (e.g., due to a parsing bug), this wraps around to a huge value since both are `u64`. The `is_empty` check uses `>=`, which would catch `start == end` but not the underflow in `len()`.

**Suggested fix:** `self.end.saturating_sub(self.start)` or validate in the constructor.

### M-4. Missing validation: `coverage <= 0.0`

**File:** `src/io/config.rs`, `validate()` function
**Category:** User Features
**Severity:** Medium

Coverage is not validated to be positive. A user could specify `coverage: 0` or `coverage: -5` and get silently incorrect behavior (0 reads generated).

**Suggested fix:** Add `anyhow::ensure!(config.sample.coverage > 0.0, ...)` to `validate()`.

### M-5. Missing validation: `read_length == 0`

**File:** `src/io/config.rs`, `validate()` function
**Category:** User Features
**Severity:** Medium

Similarly, `read_length: 0` is not rejected but would produce empty reads.

### M-6. `UmiFamily.family_size` distribution is recomputed every call

**File:** `src/umi/families.rs`, lines 60-67
**Category:** Performance
**Severity:** Medium

`sample_family_size` recomputes the LogNormal parameters and creates a new distribution object every time. With UMI mode at 1000x coverage, this is called millions of times.

**Suggested fix:** Create the `LogNormal` distribution once in the engine setup and pass it in, or use the `PcrFamilySizeSampler` struct that already exists in `fragment.rs`.

### M-7. Indel spike-in uses hardcoded "NNNNNNNNNN" padding

**File:** `src/core/engine.rs`, line ~649
**Category:** Correctness
**Severity:** Medium

```rust
spike_indel(&mut temp_read, frag_start, &variant.mutation, b"NNNNNNNNNN");
```

For deletions longer than 10bp, the padding will be insufficient and the read will be shorter than expected or filled with fewer bases than needed. The padding should be sourced from the actual reference sequence downstream of the deletion.

**Suggested fix:** Pass the pre-fetched region sequence (from fix C-3) to provide proper reference padding.

### M-8. `Variant` has redundant `chrom` field alongside `MutationType::Sv`

**File:** `src/core/types.rs`, lines 84-90 and 69
**Category:** Architecture
**Severity:** Medium

`Variant.chrom` duplicates the chromosome stored inside `MutationType::Sv { chrom, ... }`. This opens the door to inconsistency.

**Suggested fix:** Either remove `chrom` from the Sv variant or use the outer `Variant.chrom` exclusively.

### M-9. No `Eq` derive on `Variant`, `MutationType` uses `f64`

**File:** `src/core/types.rs`
**Category:** Rust Code Quality
**Severity:** Medium

`MutationType` derives `PartialEq` but cannot derive `Eq` because `Variant` contains `f64` (expected_vaf). This is a known design tension, but it means Variants cannot be used in HashSets or as HashMap keys. Consider a newtype wrapper for VAF with explicit equality semantics.

### M-10. `config.umi.clone()` inside `generate_region`

**File:** `src/core/engine.rs`, line ~452
**Category:** Rust Code Quality
**Severity:** Medium

```rust
if let Some(umi_cfg) = &self.config.umi.clone() {
```

This clones the entire UMI config on every region invocation. Use `if let Some(ref umi_cfg) = self.config.umi` instead.

### M-11. `config.artifacts.clone()` inside `generate_region`

**File:** `src/core/engine.rs`, line ~483
**Category:** Rust Code Quality
**Severity:** Medium

Same unnecessary clone pattern as M-10.

### M-12. `inject_snv_into_md` is incorrect -- appends rather than inserting at position

**File:** `src/editor/read_modifier.rs`, lines 570-577
**Category:** Correctness
**Severity:** Medium

The function just appends the ref base to the MD string rather than splitting the numeric run at the correct offset. This produces incorrect MD tags on edited BAM records, which will confuse downstream tools (Mutect2, VarDict) that rely on MD for variant calling.

**Suggested fix:** Implement proper MD string parsing and reconstruction, or omit MD tags from edited records (clear them, as done for indels).

### M-13. `BamEditor::run()` loads all records into memory before writing

**File:** `src/editor/bam_editor.rs`
**Category:** Performance
**Severity:** Medium

The editor reads all BAM records into memory (as `RecordBuf` objects) and groups them by position before writing. For WGS BAM files this can consume 50-100+ GB of RAM.

**Suggested fix:** Implement a streaming approach: process records in coordinate-sorted order, buffering only reads at positions near variants.

### M-14. `end_enrichment` division by zero for `length == 1`

**File:** `src/artifacts/ffpe.rs`, line 59
**Category:** Correctness
**Severity:** Medium

```rust
let dist_from_end = pos.min(length - 1 - pos) as f64;
```

When `length == 1` and `pos == 0`: `length - 1 - pos = 0`, `dist_from_end = 0`, `norm_dist = 0 / 0.5 = 0.0` -- this works. But when `pos > 0` and `length == 1`, `length - 1 - pos` underflows (usize). In practice pos < length so `pos == 0` always when `length == 1`, but this is fragile.

### M-15. `flate2` with `rust_backend` is slower than system zlib

**File:** `Cargo.toml`, line 41
**Category:** Performance
**Severity:** Medium

The pure-Rust `miniz_oxide` backend is ~2-3x slower than system zlib for compression. For WGS simulations producing hundreds of GB of FASTQ, this is a significant cost. Consider using `features = ["zlib"]` or `features = ["zlib-ng"]` with a fallback feature flag.

### M-16. `sample_family_size` creates LogNormal distribution per-call

**File:** `src/umi/families.rs`, line 64
**Category:** Performance
**Severity:** Medium

Same issue as M-6 but located in a different module. The `LogNormal::new(mu, sigma).expect(...)` call is repeated for every UMI family.

### M-17. Multi-ALT VCF records: only first ALT allele is parsed

**File:** `src/io/vcf_input.rs`, line 148
**Category:** User Features
**Severity:** Medium

```rust
let alt_allele = alt_field.split(',').next().unwrap_or(alt_field);
```

Multi-allelic VCF records silently drop all but the first ALT allele. At minimum, log a warning when alleles are dropped. Ideally, generate one `Variant` per ALT allele.

### M-18. `--threads` flag accepted by `learn-profile` and `edit` but unused

**File:** `src/cli/edit.rs`, line 51; `src/cli/learn_profile.rs`, line 35
**Category:** User Features
**Severity:** Medium

Both functions accept `_threads: Option<usize>` but ignore it. This is misleading to users who expect parallelism from passing `--threads 8`.

---

## Low Severity Issues

### L-1. Excessive `#[allow(dead_code)]` annotations

**Files:** Multiple (types.rs, fragment.rs, structural.rs, barcode.rs, cnv.rs, etc.)
**Category:** Rust Code Quality
**Severity:** Low

Many public functions and types are annotated `#[allow(dead_code)]`. This suppresses useful compiler warnings and suggests these APIs may not yet be integrated. A periodic audit should remove unused code or integrate it.

### L-2. `noodles_base_to_u8` is a no-op identity function

**File:** `src/editor/read_modifier.rs`, line 603-605
**Category:** Rust Code Quality
**Severity:** Low

```rust
pub fn noodles_base_to_u8(b: &u8) -> u8 {
    *b
}
```

This function just dereferences the pointer. It was likely intended for a version of noodles that used a newtype, but in the current version it serves no purpose.

### L-3. `variant_position` duplicated across modules

**File:** `src/core/engine.rs` (line ~622) vs `src/io/truth_vcf.rs` (line ~146)
**Category:** Architecture
**Severity:** Low

Two identical implementations of position extraction from `Variant`. Use the `Variant::pos()` method already defined in `truth_vcf.rs`.

### L-4. Missing `#[must_use]` on pure functions

**Files:** `src/variants/vaf.rs::expected_vaf`, `src/core/coverage.rs::read_pairs_for_coverage`, etc.
**Category:** Rust Code Quality
**Severity:** Low

Several pure functions that compute values should have `#[must_use]` to catch accidentally discarded results.

### L-5. `current_timestamp_utc` hand-rolls date formatting

**File:** `src/io/manifest.rs`, lines 92-130
**Category:** Rust Code Quality
**Severity:** Low

The manual UTC timestamp formatting avoids a dependency but re-implements date arithmetic (leap years, month lengths). This is fragile. Consider using `std::time::SystemTime` with a minimal formatting crate, or at least add more thorough tests for edge cases (year boundaries, leap years).

### L-6. `ReadPair.name` is a `String` owned per pair

**File:** `src/core/types.rs`, line 30
**Category:** Performance
**Severity:** Low

Combined with H-5 and H-6, each ReadPair owns three Strings (name, chrom, and the fragment_start formatted into name). For very high-coverage simulations, this memory overhead adds up.

### L-7. Tests construct `Config` structs with full field lists

**Files:** Many test modules (cancer_presets, multi_sample, presets, manifest, etc.)
**Category:** Testing
**Severity:** Low

Config construction in tests is verbose and fragile -- adding a new field to `Config` requires updating every test helper. Consider providing a `Config::default_for_testing()` method or a builder.

### L-8. `simulate.rs` writer thread drops errors silently on channel closure

**File:** `src/cli/simulate.rs` (writer thread)
**Category:** Correctness
**Severity:** Low

If the writer thread encounters an error (e.g., disk full), it may not propagate back to the main thread cleanly depending on channel ordering.

---

## Testing Gaps

1. **No integration test for the full simulate pipeline.** Unit tests cover individual modules well, but there is no end-to-end test that runs `simulate` with a real (small) FASTA reference and verifies the output FASTQ/BAM/VCF contents are valid and contain the expected variants.

2. **No test for `generate_region` with empirical quality model.** The engine tests only exercise the parametric model path.

3. **No test for multi-threaded determinism.** The per-region seed derivation is tested, but there is no test that runs the full parallel pipeline twice with the same seed and asserts byte-identical output.

4. **No test for the streaming output pipeline** (`simulate.rs` with crossbeam channel). This is a complex concurrent system that should have at least a smoke test.

5. **Edge case missing:** Fragment length > region length. The engine clamps `max_start_offset` but does not test this boundary condition.

6. **Edge case missing:** Variant at position 0 (first base of chromosome). The 0-based to 1-based coordinate conversion in VCF writing should be tested at this boundary.

---

## Architecture Observations

1. **Good:** Clean separation of concerns between `core/` (simulation engine), `io/` (file I/O), `variants/` (mutation modeling), and `cli/` (user interface). Each module is independently testable.

2. **Good:** The `FragmentSampler` trait and `QualityModel` trait enable pluggable implementations. The empirical quality model cleanly extends the parametric one.

3. **Good:** Stochastic VAF sampling via binomial distribution is a clear improvement over deterministic spike-in tools.

4. **Improvement opportunity:** The `SimulationEngine` struct is becoming a "god object" -- `generate_region` is ~400 lines with deeply nested logic. Consider extracting a `FragmentGenerator`, `VariantSpiker`, and `ArtifactInjector` that each handle one phase of the pipeline.

5. **Improvement opportunity:** `ReadPair` carries too much per-pair data (owned Strings). A more memory-efficient design would separate read sequences from metadata, using indices or references into shared data.

6. **Improvement opportunity:** The `Variant` / `MutationType` type hierarchy conflates the variant specification (what mutation to spike) with the positional encoding. An SV variant stores its chromosome redundantly. Consider separating the coordinate from the mutation description.

---

## Summary of Top Priority Actions

1. **Pre-fetch region sequence** (C-3) -- single largest performance win, eliminates per-fragment file I/O.
2. **Eliminate `apply_variant_to_seq` clone** (C-1) -- remove per-variant-per-read allocation in hot path.
3. **Hoist quality model construction** (C-2) -- avoid per-read object creation.
4. **Fix cfDNA parameter passthrough** (H-3) -- silent correctness bug.
5. **Fix `Read::new` panic** (H-9) -- production crash risk.
6. **Fix edit purity override logic** (H-8) -- user-facing correctness bug.
7. **Add end-to-end integration test** -- validates the full pipeline.
