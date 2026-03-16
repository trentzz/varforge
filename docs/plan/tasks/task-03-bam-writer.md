# Task 03: BAM Writer

## Phase
1 - Foundation I/O

## Dependencies
None

## Objective
Implement a BAM writer that outputs aligned read pairs with proper headers, CIGAR strings, mate information, and optional UMI/duplex tags.

## Context
VarForge needs to output BAM files for users who want to skip the alignment step or need pre-aligned data for testing downstream tools. The project uses `noodles-bam` and `noodles-sam` (pure Rust).

BAM records need: proper header with @SQ lines from reference, @RG read groups, CIGAR from spike-in operations, mate pair linkage, and auxiliary tags (RX for UMI, MI for family ID).

## Key Files
- **Create**: `src/io/bam.rs`
- **Modify**: `src/io/mod.rs` (add module export)
- **Existing context**: `src/core/types.rs` (ReadPair), `src/io/config.rs` (SampleConfig for read groups)

## Requirements

### Struct: `BamWriter`
```rust
pub struct BamWriter {
    // noodles BAM writer
    // Reference sequence dictionary for header
}
```

### Methods
1. `BamWriter::new(path: &Path, ref_sequences: &[(String, u64)], sample_config: &SampleConfig) -> Result<Self>` - Create BAM with header
2. `fn write_pair(&mut self, pair: &ReadPair, ref_id: usize, pos: u64, cigar_r1: &str, cigar_r2: &str) -> Result<()>` - Write aligned pair
3. `fn finish(self) -> Result<()>` - Finalize BAM

### BAM Header
- `@HD VN:1.6 SO:coordinate`
- `@SQ SN:{chrom} LN:{length}` for each reference sequence
- `@RG ID:{sample} SM:{sample} PL:{platform} LB:{sample}`

### SAM Tags to Set
- `RG:Z:{read_group}` - Read group
- `RX:Z:{umi}` - UMI barcode (if UMI mode enabled)
- `MI:i:{family_id}` - Molecular identifier / family ID (if UMI mode)
- Proper mate flags (0x1 paired, 0x2 proper pair, 0x20 mate reverse, 0x40/0x80 first/second in pair)

## Tests

### Unit Tests
1. `test_write_header` - Valid BAM header with SQ and RG lines
2. `test_write_single_pair` - Write and read back one pair, verify fields
3. `test_mate_flags` - Proper paired-end flags set
4. `test_umi_tags` - RX and MI tags present when UMI data provided
5. `test_cigar_preserved` - CIGAR strings from spike-in preserved in output
6. `test_read_group` - RG tag matches header

### Notes
- Use `noodles-bam` reader in tests to verify written output
- Use `tempfile` for output paths

## Acceptance Criteria
- [ ] Valid BAM files produced with correct header
- [ ] Paired-end flags and mate information correct
- [ ] UMI tags (RX, MI) written when present
- [ ] CIGAR strings preserved from variant spike-in
- [ ] All 6 tests pass
- [ ] `cargo clippy` clean
