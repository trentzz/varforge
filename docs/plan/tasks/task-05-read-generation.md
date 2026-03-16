# Task 05: Read Generation Engine

## Phase
2 - Core Pipeline

## Dependencies
Task 01 (Reference Genome Loading)

## Objective
Implement the core read generation pipeline that ties together fragment sampling, sequence extraction, variant spike-in, quality modeling, UMI attachment, and artifact injection to produce `ReadPair` structs for a given genomic region.

## Context
This is the central engine of VarForge. All the individual modules (fragment, quality, variants, UMI, artifacts) are already implemented with unit tests. This task wires them together into a coherent pipeline that processes one region at a time.

The pipeline for generating reads in a single region:
1. Calculate number of read pairs needed for target coverage
2. For each read pair:
   a. Sample fragment length from appropriate model (normal/cfDNA)
   b. Sample random position within region
   c. Extract reference sequence for fragment
   d. Check if any variants overlap this fragment; if so, apply spike-in (stochastic VAF sampling)
   e. Generate quality scores with position-dependent decay
   f. Inject sequencing errors based on quality scores
   g. If UMI mode: attach UMI barcode, generate PCR family copies
   h. If artifact mode: inject FFPE/oxoG damage
3. Collect all ReadPairs for the region

## Key Files
- **Create**: `src/core/engine.rs`
- **Modify**: `src/core/mod.rs` (add module export)
- **Existing context**: All files in `src/core/`, `src/variants/`, `src/umi/`, `src/artifacts/`

## Requirements

### Struct: `SimulationEngine`
```rust
pub struct SimulationEngine {
    config: SimulationConfig,
    reference: ReferenceGenome,
    rng: StdRng,
}
```

### Core Method
```rust
pub fn generate_region(
    &mut self,
    region: &Region,
    variants: &[Variant],
) -> Result<RegionOutput>
```

### Struct: `RegionOutput`
```rust
pub struct RegionOutput {
    pub read_pairs: Vec<ReadPair>,
    pub applied_variants: Vec<AppliedVariant>,  // variants that were actually spiked
}

pub struct AppliedVariant {
    pub variant: Variant,
    pub actual_alt_count: u32,
    pub actual_total_count: u32,
}
```

### Pipeline Steps (in order)
1. `coverage::read_pairs_for_coverage()` - How many pairs needed
2. For each pair:
   - `fragment::*Sampler::sample()` - Fragment length
   - Reference sequence extraction from `ReferenceGenome`
   - `variants::vaf::sample_alt_count()` + `variants::spike_in::spike_*()` - Variant injection
   - `core::quality::ParametricQualityModel::generate_qualities()` - Quality scores
   - `core::quality::inject_errors()` - Sequencing errors
3. Post-pair processing:
   - `umi::barcode::generate_umi()` - UMI attachment (if enabled)
   - `umi::families::generate_pcr_copies()` - PCR family expansion (if enabled)
   - `artifacts::ffpe::inject_ffpe_damage()` - FFPE artifacts (if enabled)
   - `artifacts::ffpe::inject_oxog_damage()` - OxoG artifacts (if enabled)
   - `artifacts::duplicates::select_duplicates()` - PCR duplicates (if enabled)

## Tests

### Unit Tests
1. `test_generate_region_basic` - Generates expected number of read pairs for a simple region
2. `test_read_lengths_correct` - All reads have configured read length
3. `test_variant_spiked` - A variant at 100% VAF appears in all reads overlapping its position
4. `test_variant_stochastic` - A variant at 50% VAF produces alt reads ~50% of the time (statistical test)
5. `test_no_variants` - Region with no variants produces reference-only reads
6. `test_umi_attached` - When UMI enabled, read pairs have UMI metadata
7. `test_pcr_families` - When UMI enabled, PCR copies produced for each original
8. `test_artifacts_applied` - When FFPE enabled, C>T transitions present at configured rate
9. `test_cfDNA_fragments` - cfDNA mode produces fragments with expected nucleosomal distribution
10. `test_deterministic_with_seed` - Same seed produces same output
11. `test_applied_variants_tracking` - AppliedVariant correctly tracks actual alt/total counts

### Test Fixtures
Use a small in-memory reference (e.g., 1000bp of known sequence) rather than a real genome file. Create a helper that builds a minimal `SimulationConfig`.

## Acceptance Criteria
- [ ] Pipeline generates ReadPairs for arbitrary regions
- [ ] Variants correctly spiked with stochastic VAF
- [ ] UMI, artifacts, PCR families all integrate correctly
- [ ] All 11 tests pass
- [ ] `cargo test --lib core::engine` passes
- [ ] `cargo clippy` clean
