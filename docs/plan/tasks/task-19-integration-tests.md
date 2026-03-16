# Task 19: Integration Tests

## Phase
6 - Quality & Polish

## Dependencies
Task 06 (Simulate Command Orchestrator)

## Objective
Write comprehensive integration tests that validate end-to-end simulation correctness: the full pipeline from config to output files, with verification of output file formats, variant presence, and statistical properties.

## Context
Individual modules have unit tests, but we need integration tests that exercise the full pipeline and catch issues at module boundaries. These tests should cover every major use case.

## Key Files
- **Create**: `tests/integration/` directory with test files
- **Create**: `tests/fixtures/` with test reference, config files
- **Create**: `tests/common/mod.rs` with shared test utilities

## Requirements

### Test Fixtures
Create minimal but complete test data:
- Small reference FASTA (~10kb, 3 chromosomes) with known sequence
- Corresponding .fai index
- Config YAML files for each test scenario
- Optional: small VCF with known mutations

### Integration Test Scenarios

1. **`test_minimal_simulation`** - Smallest valid config (1 region, 1x coverage)
   - Verify: FASTQ files exist, are valid gzip, contain reads of correct length

2. **`test_snv_spike_in`** - 10 SNVs at 50% VAF, 100x coverage
   - Verify: Truth VCF has 10 variants, FASTQ reads at variant positions show ~50% alt alleles

3. **`test_indel_spike_in`** - Insertions and deletions at various sizes
   - Verify: Truth VCF correct, reads show indels

4. **`test_umi_mode`** - UMI enabled, simplex mode
   - Verify: Read headers contain RX:Z tag, PCR family copies exist

5. **`test_duplex_mode`** - UMI duplex mode
   - Verify: Duplex UMI pairs (A+B / B+A) present

6. **`test_cfdna_mode`** - cfDNA fragment model
   - Verify: Fragment size distribution peaks near 167bp

7. **`test_ffpe_artifacts`** - FFPE damage enabled
   - Verify: Elevated C>T transition rate

8. **`test_tumor_purity`** - 50% purity
   - Verify: VAF of clonal variants approximately halved

9. **`test_subclonal_variants`** - Clonal tree with subclones
   - Verify: Subclonal variants have lower VAF than clonal

10. **`test_targeted_panel`** - BED file with target regions
    - Verify: Reads only in targeted regions (plus off-target if modeled)

11. **`test_deterministic_seed`** - Same seed, same output
    - Verify: Byte-identical FASTQ output for same seed

12. **`test_bam_output`** - BAM output enabled
    - Verify: Valid BAM with correct header, flags, tags

13. **`test_high_coverage`** - 1000x coverage on small region
    - Verify: Approximately correct read count, reasonable memory usage

14. **`test_low_vaf_detection`** - Variants at 1% VAF, 1000x coverage
    - Verify: Alt reads present but rare (~10 per position)

15. **`test_multi_sample`** - Longitudinal simulation (if Task 15 complete)
    - Verify: Same mutations across samples, different VAFs

### Validation Utilities
```rust
// In tests/common/mod.rs
fn count_fastq_records(path: &Path) -> usize;
fn extract_variants_from_vcf(path: &Path) -> Vec<Variant>;
fn check_alt_allele_fraction(fastq_r1: &Path, pos: u64, alt: u8) -> f64;
fn verify_bam_header(path: &Path, expected_samples: &[&str]) -> bool;
fn fragment_size_distribution(fastq_r1: &Path, fastq_r2: &Path) -> Vec<usize>;
```

## Tests
This task IS the tests. All 15 scenarios above.

## Acceptance Criteria
- [ ] All 15 integration test scenarios pass
- [ ] Tests run in < 60 seconds total (small test data)
- [ ] Test fixtures are minimal and self-contained
- [ ] `cargo test --test integration` passes
- [ ] No flaky tests (deterministic with seed)
