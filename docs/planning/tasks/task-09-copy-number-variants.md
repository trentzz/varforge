# Task 09: Copy Number Variant Support

## Phase
3 - Extended Variant Support

## Dependencies
Task 05 (Read Generation Engine)

## Objective
Implement copy number variant simulation: amplifications and deletions that affect coverage depth across regions, integrated with the tumour purity/ploidy model.

## Context
CNVs are fundamental to cancer genomics. A region with copy number 4 (amplification) should have ~2x the expected coverage compared to diploid (CN=2). A region with copy number 1 (heterozygous deletion) should have ~0.5x. This interacts with tumour purity: coverage = depth × (purity × CN_tumor + (1-purity) × CN_normal) / 2.

Target tools: CNVkit, GATK CNV, ichorCNA, FACETS

## Key Files
- **Create**: `src/variants/cnv.rs`
- **Modify**: `src/variants/mod.rs`, `src/core/engine.rs` (coverage adjustment per region)
- **Existing context**: `src/tumour/clonal_tree.rs`, `src/variants/vaf.rs` (already has CN in VAF equation)

## Requirements

### Config Extension
```yaml
copy_number:
  - region: "chr7:55000000-55200000"
    tumor_cn: 4    # amplification (EGFR)
    normal_cn: 2
  - region: "chr17:7500000-7700000"
    tumor_cn: 1    # het deletion (TP53)
    normal_cn: 2
  - region: "chr13:32300000-32400000"
    tumor_cn: 0    # homozygous deletion (BRCA2)
    normal_cn: 2
```

### Coverage Adjustment
For a region with tumor CN and normal CN:
```
adjusted_coverage = base_coverage × (purity × tumor_cn + (1 - purity) × normal_cn) / ploidy
```

### Allele-Specific CN
Support major/minor allele copy numbers for LOH:
- CN=2 with major=2, minor=0 → LOH (all reads from one allele)
- CN=3 with major=2, minor=1 → unbalanced gain
- Affects VAF of heterozygous variants in the region

### Integration with VAF Model
The existing `expected_vaf()` function in `vaf.rs` already takes `cn_tumor` and `cn_normal` parameters. This task ensures the engine passes correct CN values per-region.

## Tests

### Unit Tests
1. `test_amplification_coverage` - CN=4 region has ~2x coverage vs CN=2
2. `test_deletion_coverage` - CN=1 region has ~0.5x coverage
3. `test_homodel_coverage` - CN=0 region has only normal contamination reads
4. `test_loh_allele_frequency` - LOH region: het SNP shows VAF ~1.0 for major allele
5. `test_cn_with_purity` - Coverage correctly adjusted for tumor purity
6. `test_vaf_in_amplified_region` - VAF shifts correctly for amplified variants
7. `test_cn_config_parsing` - YAML config with copy_number section parses correctly

## Acceptance Criteria
- [ ] Coverage depth adjusts based on copy number and purity
- [ ] Allele-specific CN affects variant allele frequencies
- [ ] LOH correctly modeled
- [ ] Homozygous deletions produce only normal-contamination reads
- [ ] All 7 tests pass
- [ ] `cargo clippy` clean
