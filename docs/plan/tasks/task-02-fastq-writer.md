# Task 02: FASTQ Writer

## Phase
1 - Foundation I/O

## Dependencies
None

## Objective
Implement a streaming, gzip-compressed FASTQ writer that outputs paired-end reads (R1/R2) from `ReadPair` structs.

## Context
VarForge generates `ReadPair` structs (defined in `src/core/types.rs`) containing two `Read`s, each with a sequence (`Vec<u8>`) and quality scores (`Vec<u8>`). These need to be written as paired FASTQ files (sample_R1.fastq.gz, sample_R2.fastq.gz).

The project has `gzp` for multi-threaded gzip compression and `noodles-fastq` available.

## Key Files
- **Create**: `src/io/fastq.rs`
- **Modify**: `src/io/mod.rs` (add module export)
- **Existing context**: `src/core/types.rs` (ReadPair, Read structs)

## Requirements

### Struct: `FastqWriter`
```rust
pub struct FastqWriter {
    // Gzip-compressed writers for R1 and R2
}
```

### Methods
1. `FastqWriter::new(output_dir: &Path, sample_name: &str) -> Result<Self>` - Create R1/R2 gzipped files
2. `fn write_pair(&mut self, pair: &ReadPair, read_name: &str) -> Result<()>` - Write one pair
3. `fn write_pairs(&mut self, pairs: &[ReadPair], name_prefix: &str) -> Result<()>` - Write batch with auto-naming
4. `fn finish(self) -> Result<()>` - Flush and finalize gzip streams

### FASTQ Format
```
@{read_name}/1
{sequence}
+
{quality_as_ascii}
```

Quality scores are Phred+33 encoded (add 33 to raw quality value for ASCII).

### Read Naming Convention
`@{sample}:{region}:{index}/1` for R1, `/2` for R2. Include UMI in header if present:
`@{sample}:{region}:{index} RX:Z:{umi}/1`

## Tests

### Unit Tests
1. `test_write_single_pair` - Write one pair, decompress, verify FASTQ format
2. `test_write_multiple_pairs` - Write batch, verify count and ordering
3. `test_quality_encoding` - Phred+33 encoding correct (Q30 = '?', Q0 = '!')
4. `test_gzip_compressed` - Output files are valid gzip
5. `test_paired_consistency` - R1 and R2 have same number of records
6. `test_empty_writes` - Handle zero pairs gracefully
7. `test_output_paths` - Files named correctly (sample_R1.fastq.gz, sample_R2.fastq.gz)

### Test Helpers
Use `tempfile::TempDir` for output. Use `flate2::read::GzDecoder` to verify output.

## Acceptance Criteria
- [ ] Paired FASTQ files written with correct format
- [ ] Gzip compression working (multi-threaded via gzp)
- [ ] Quality scores properly Phred+33 encoded
- [ ] All 7 tests pass
- [ ] `cargo test --lib io::fastq` passes
- [ ] `cargo clippy` clean
