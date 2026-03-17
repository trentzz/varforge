# Task 10: VCF Input Parsing

## Phase
3 - Extended Variant Support

## Dependencies
Task 05 (Read Generation Engine)

## Objective
Implement VCF file parsing so users can supply a list of mutations to spike in, rather than relying solely on random mutation generation.

## Context
Many users have a specific set of mutations they want to test (e.g., known cancer hotspot mutations, a list from COSMIC, or mutations detected in a real patient). The config already has `vcf_path` in `MutationConfig` but parsing is not implemented.

## Key Files
- **Create**: `src/io/vcf_input.rs`
- **Modify**: `src/io/mod.rs`, `src/core/engine.rs` (use parsed variants)
- **Existing context**: `src/io/config.rs` (MutationConfig has vcf_path), `src/core/types.rs` (Variant struct)

## Requirements

### Parsing
- Read VCF using `noodles-vcf`
- Extract: CHROM, POS, REF, ALT, and optionally INFO fields for VAF/clone assignment
- Support both VCF and VCF.gz (bgzipped) input
- Map VCF records to `Variant` structs

### VCF INFO Fields (optional, with defaults)
- `VAF` or `AF` → target VAF (default: use tumour model to calculate)
- `CLONE` → clone assignment (default: "founder")
- `CN` → copy number at this locus (default: 2)

### Variant Type Detection
- REF.len() == 1 && ALT.len() == 1 → SNV
- REF.len() > ALT.len() → Deletion
- REF.len() < ALT.len() → Insertion
- REF.len() > 1 && ALT.len() > 1 && REF.len() == ALT.len() → MNV
- ALT contains BND notation → SV (Task 08)

### Validation
- Verify REF allele matches reference genome at that position
- Warn (don't error) on REF mismatch (common with different reference builds)
- Skip variants on chromosomes not in reference
- Report count of loaded/skipped variants

## Tests

### Unit Tests
1. `test_parse_snv` - Single SNV from VCF
2. `test_parse_indel` - Insertion and deletion from VCF
3. `test_parse_mnv` - Multi-nucleotide variant
4. `test_parse_multiple` - Multiple variants from one VCF
5. `test_custom_vaf` - VAF from INFO field overrides default
6. `test_clone_assignment` - CLONE from INFO field
7. `test_missing_info_defaults` - Defaults applied when INFO fields absent
8. `test_ref_mismatch_warning` - Warning on REF mismatch (not error)
9. `test_skip_unknown_chrom` - Variants on unknown chromosomes skipped with warning
10. `test_bgzipped_vcf` - Parse .vcf.gz input

### Test Fixtures
Create small VCF files in test setup using string literals.

## Acceptance Criteria
- [ ] VCF and VCF.gz files parsed into Variant structs
- [ ] Custom VAF and clone assignment from INFO fields
- [ ] Ref allele validation with warning on mismatch
- [ ] All 10 tests pass
- [ ] `cargo clippy` clean
