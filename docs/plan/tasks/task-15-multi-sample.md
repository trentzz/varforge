# Task 15: Multi-Sample Longitudinal Simulation

## Phase
5 - Extended Simulation Modes

## Dependencies
Task 06 (Simulate Command Orchestrator)

## Objective
Support simulating multiple samples from the same patient with shared clonal architecture but different tumor fractions, coverage levels, and time points -- enabling longitudinal monitoring benchmarks.

## Context
Liquid biopsy monitoring tracks tumor burden over time. A patient might have samples at diagnosis (TF=5%), post-treatment (TF=0.1%), and relapse (TF=3%). All samples share the same clonal tree and mutations, but the tumor fraction and possibly clone proportions change.

No existing tool supports this natively.

## Key Files
- **Create**: `src/core/multi_sample.rs`
- **Modify**: `src/io/config.rs` (extend config for multi-sample)
- **Modify**: `src/cli/simulate.rs` (loop over samples)

## Requirements

### Config Extension
```yaml
samples:
  - name: "diagnosis"
    coverage: 200
    tumour_fraction: 0.05
    fragment_model: cfDNA
  - name: "post_treatment"
    coverage: 200
    tumour_fraction: 0.001
    fragment_model: cfDNA
  - name: "relapse"
    coverage: 200
    tumour_fraction: 0.03
    fragment_model: cfDNA
    clonal_shift:  # optional: adjust clone CCFs
      clone_A: 0.8  # increased
      clone_B: 0.1  # decreased
```

### Shared State
- Same reference genome
- Same mutation list
- Same clonal tree structure
- Per-sample: tumor fraction, coverage, fragment model, optional CCF adjustments

### Output
Each sample gets its own output directory:
```
output/
  diagnosis/
    diagnosis_R1.fastq.gz
    diagnosis_R2.fastq.gz
    truth.vcf
  post_treatment/
    ...
  relapse/
    ...
  manifest.json  # Combined manifest with all samples
```

## Tests

### Unit Tests
1. `test_multi_sample_config_parsing` - YAML with multiple samples parses correctly
2. `test_shared_mutations` - All samples have same mutations in truth VCF
3. `test_different_tumor_fractions` - VAFs differ between samples as expected
4. `test_clonal_shift` - CCF adjustments applied correctly
5. `test_per_sample_output` - Each sample has its own output directory
6. `test_combined_manifest` - Manifest includes all samples

## Acceptance Criteria
- [ ] Multiple samples simulated from shared clonal architecture
- [ ] Tumor fraction varies per sample
- [ ] Clonal CCF shifts supported
- [ ] All 6 tests pass
- [ ] `cargo clippy` clean
