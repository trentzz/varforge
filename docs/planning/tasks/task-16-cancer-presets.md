# Task 16: Cancer-Type Presets

## Phase
5 - Extended Simulation Modes

## Dependencies
Task 13 (CLI Overrides & Presets)

## Objective
Provide built-in cancer-type presets that configure realistic mutation spectra, variant counts, tumor purity ranges, and common driver mutations for specific cancer types.

## Context
MOV&RSim (2025) provides presets for 21 cancer types, which is a key differentiator. VarForge should offer similar capability with better performance and additional features (UMI, cfDNA).

Presets are based on COSMIC mutation signatures, typical mutation burdens, and common driver genes per cancer type.

## Key Files
- **Create**: `src/cli/cancer_presets.rs`
- **Modify**: `src/cli/presets.rs` (integrate cancer presets)

## Requirements

### Cancer Type Presets
```
--preset cancer:lung_adeno     # NSCLC: high TMB, EGFR/KRAS drivers, smoking signature (SBS4)
--preset cancer:colorectal     # CRC: MSI-H option, APC/KRAS/TP53, SBS1/SBS5
--preset cancer:breast_tnbc    # TNBC: high TMB, TP53/BRCA1, SBS3 (HRD)
--preset cancer:melanoma       # High TMB, BRAF V600E, UV signature (SBS7a/b)
--preset cancer:aml            # Low TMB, FLT3/NPM1/DNMT3A
--preset cancer:prostate       # Low TMB, TMPRSS2-ERG fusion, AR
--preset cancer:pancreatic     # Low purity, KRAS G12D, TP53, SMAD4
--preset cancer:glioblastoma   # IDH1/EGFR/PTEN, moderate TMB
```

### Preset Contents
Each preset defines:
- Typical mutation burden (mutations per Mb)
- Mutation type distribution (SNV/indel/SV/CNV ratios)
- Dominant COSMIC SBS signatures and their weights
- Common driver mutations (with configurable inclusion)
- Typical tumor purity range
- Typical ploidy
- Suggested coverage and sequencing mode

### Mutation Signature Application
- Weight the substitution type distribution by COSMIC signatures
- SBS4 (smoking): C>A dominant
- SBS7 (UV): C>T at dipyrimidines
- SBS3 (HRD): relatively flat spectrum
- SBS1/5 (aging): C>T at CpG

## Tests

### Unit Tests
1. `test_all_cancer_presets_valid` - Every preset produces a valid config
2. `test_lung_adeno_signatures` - Lung preset has SBS4 signature characteristics
3. `test_melanoma_high_tmb` - Melanoma has high mutation burden
4. `test_aml_low_tmb` - AML has low mutation burden
5. `test_driver_mutations_included` - Known drivers present in mutation list
6. `test_preset_override` - User can override preset values via CLI

## Acceptance Criteria
- [ ] At least 8 cancer type presets
- [ ] Mutation spectra reflect COSMIC signatures
- [ ] Common driver mutations included
- [ ] All 6 tests pass
- [ ] `cargo clippy` clean
