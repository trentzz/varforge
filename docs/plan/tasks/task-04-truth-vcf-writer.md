# Task 04: Truth VCF Writer

## Phase
1 - Foundation I/O

## Dependencies
None

## Objective
Implement a truth VCF writer that outputs all spiked-in variants with their expected VAFs, clone assignments, and variant types as a ground-truth file for benchmarking.

## Context
The truth VCF is essential for benchmarking: users compare their variant caller's output against VarForge's truth VCF to compute sensitivity, specificity, precision, and recall. Every variant spiked into the simulated reads must appear in the truth VCF.

The project uses `noodles-vcf` (pure Rust).

## Key Files
- **Create**: `src/io/truth_vcf.rs`
- **Modify**: `src/io/mod.rs` (add module export)
- **Existing context**: `src/core/types.rs` (Variant, MutationType), `src/variants/vaf.rs` (expected_vaf)

## Requirements

### Struct: `TruthVcfWriter`
```rust
pub struct TruthVcfWriter {
    // noodles VCF writer
}
```

### Methods
1. `TruthVcfWriter::new(path: &Path, sample_name: &str, contigs: &[(String, u64)]) -> Result<Self>` - Create truth VCF
2. `fn write_variant(&mut self, variant: &Variant, ref_allele: &[u8], alt_allele: &[u8]) -> Result<()>` - Write one variant
3. `fn finish(self) -> Result<()>` - Finalize

### VCF Header
- `##fileformat=VCFv4.3`
- `##source=VarForge`
- `##INFO=<ID=EXPECTED_VAF,Number=1,Type=Float,Description="Expected variant allele frequency">`
- `##INFO=<ID=CLONE,Number=1,Type=String,Description="Clone assignment">`
- `##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Variant type: SNV, INDEL, MNV, SV, CNV">`
- `##INFO=<ID=CCF,Number=1,Type=Float,Description="Cancer cell fraction of assigned clone">`
- `##contig=<ID={chrom},length={length}>` for each reference sequence
- `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`
- `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT {sample}`

### Record Format
- CHROM, POS from variant
- REF, ALT alleles
- QUAL = "." (not applicable for truth)
- FILTER = "PASS"
- INFO = `EXPECTED_VAF={vaf};CLONE={clone};VARTYPE={type};CCF={ccf}`
- FORMAT/SAMPLE = `GT` with `0/1` (heterozygous)

## Tests

### Unit Tests
1. `test_write_header` - Valid VCF header with all INFO fields
2. `test_write_snv` - SNV record with correct fields
3. `test_write_indel` - Insertion and deletion records
4. `test_write_mnv` - MNV record
5. `test_expected_vaf_in_info` - EXPECTED_VAF present and correct
6. `test_clone_assignment` - CLONE info field populated
7. `test_multiple_variants` - Multiple records in sorted order
8. `test_contig_headers` - Contig lines match reference

## Acceptance Criteria
- [ ] Valid VCF 4.3 format output
- [ ] All spiked variants represented with correct alleles
- [ ] EXPECTED_VAF, CLONE, VARTYPE, CCF in INFO field
- [ ] All 8 tests pass
- [ ] `cargo clippy` clean
