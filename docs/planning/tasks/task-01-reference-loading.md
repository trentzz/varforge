# Task 01: Reference Genome Loading

## Phase
1 - Foundation I/O

## Dependencies
None

## Objective
Implement a `ReferenceGenome` struct that loads an indexed FASTA file and provides efficient random-access sequence extraction for arbitrary genomic regions.

## Context
The simulation pipeline needs to read a reference genome (e.g., hg38.fa) to:
1. Extract sequences for fragment generation (reads are subsequences of reference + mutations)
2. Validate variant positions against reference alleles
3. Determine chromosome lengths for coverage calculations

The project already uses `noodles-fasta` (pure Rust). The FASTA must be indexed (.fai) for random access.

## Key Files
- **Create**: `src/io/reference.rs`
- **Modify**: `src/io/mod.rs` (add module export)
- **Existing context**: `src/core/types.rs` (Region struct), `src/core/coverage.rs` (uses regions)

## Requirements

### Struct: `ReferenceGenome`
```rust
pub struct ReferenceGenome {
    // Indexed FASTA reader
    // Chromosome name -> length mapping
}
```

### Methods
1. `ReferenceGenome::open(path: &Path) -> Result<Self>` - Open indexed FASTA, build chrom->length map
2. `fn sequence(&self, region: &Region) -> Result<Vec<u8>>` - Extract uppercase sequence for a region
3. `fn chromosome_lengths(&self) -> &HashMap<String, u64>` - Get all chromosome lengths
4. `fn contains_chromosome(&self, name: &str) -> bool` - Check if chromosome exists
5. `fn genome_size(&self) -> u64` - Total genome size

### Edge Cases
- Handle both "chr1" and "1" naming conventions
- Return error for out-of-bounds regions
- Uppercase all returned sequence (handle mixed-case references)
- Handle N-bases gracefully (they exist in reference)

## Tests

### Unit Tests (in `src/io/reference.rs`)
1. `test_open_missing_file` - Error on nonexistent file
2. `test_open_missing_index` - Error when .fai missing
3. `test_sequence_extraction` - Correct sequence for known region (use tempfile with small test FASTA)
4. `test_sequence_uppercase` - Mixed-case input returns uppercase
5. `test_chromosome_lengths` - Correct lengths from index
6. `test_out_of_bounds` - Error for region beyond chromosome length
7. `test_contains_chromosome` - True/false for known/unknown chroms
8. `test_genome_size` - Sum of all chromosome lengths

### Test Fixtures
Create a small test FASTA + .fai in a helper function:
```
>chr1
ACGTACGTACGTACGT
>chr2
NNNNACGTNNNN
```

## Acceptance Criteria
- [ ] `ReferenceGenome::open()` loads indexed FASTA
- [ ] `sequence()` returns correct bases for arbitrary regions
- [ ] All chromosome lengths accessible
- [ ] All 8 tests pass
- [ ] `cargo test --lib io::reference` passes
- [ ] `cargo clippy` clean
