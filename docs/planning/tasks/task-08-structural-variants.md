# Task 08: Structural Variant Support

## Phase
3 - Extended Variant Support

## Dependencies
Task 05 (Read Generation Engine)

## Objective
Add support for simulating structural variants: deletions (>50bp), insertions (>50bp), inversions, duplications, and translocations. Generate split reads and discordant pairs at SV breakpoints.

## Context
Structural variants are critical for cancer genomics (driver fusions like BCR-ABL, large deletions of tumor suppressors). Current spike-in supports only small indels (<50bp). SV simulation requires:
- Modifying the reference at breakpoints
- Generating split reads that span breakpoints
- Creating discordant read pairs (abnormal insert size or wrong orientation)
- Proper CIGAR strings with soft-clipping at breakpoints

Target tools: Manta, DELLY, LUMPY, SvABA, GRIDSS

## Key Files
- **Create**: `src/variants/structural.rs`
- **Modify**: `src/variants/mod.rs`, `src/core/types.rs` (extend MutationType)
- **Existing context**: `src/variants/spike_in.rs` (small variant model to follow)

## Requirements

### SV Types
```rust
pub enum StructuralVariant {
    Deletion { chrom: String, start: u64, end: u64 },
    Insertion { chrom: String, pos: u64, sequence: Vec<u8> },
    Inversion { chrom: String, start: u64, end: u64 },
    Duplication { chrom: String, start: u64, end: u64, copies: u32 },
    Translocation { chrom1: String, pos1: u64, chrom2: String, pos2: u64 },
}
```

### Read Generation at SVs
- **Deletion**: Reads spanning the deletion show a gap; split reads have soft-clipped segments
- **Insertion**: Reads at insertion point carry inserted sequence; excess bases soft-clipped
- **Inversion**: Reads in inverted region have reverse-complement sequence; discordant pairs at boundaries
- **Duplication**: Extra coverage in duplicated region; split reads at tandem dup boundaries
- **Translocation**: Chimeric reads spanning the breakpoint; discordant pairs mapping to different chromosomes

### Stochastic Sampling
SVs use the same VAF/CCF model as small variants -- binomial sampling determines which reads carry the SV allele.

## Tests

### Unit Tests
1. `test_deletion_split_reads` - Reads spanning deletion have correct soft-clipping
2. `test_deletion_coverage_drop` - Coverage drops in deleted region for affected allele
3. `test_insertion_reads` - Reads at insertion point carry inserted sequence
4. `test_inversion_orientation` - Reads in inverted region have reversed sequence
5. `test_duplication_coverage` - Duplicated region has elevated coverage
6. `test_translocation_chimeric` - Chimeric reads span breakpoint
7. `test_sv_vaf_stochastic` - SV at 50% VAF produces alt reads ~50% of time
8. `test_sv_truth_vcf` - SVs appear in truth VCF with correct breakpoint notation (BND format)

## Acceptance Criteria
- [ ] All 5 SV types generate correct read patterns
- [ ] Split reads and discordant pairs at breakpoints
- [ ] Stochastic VAF sampling for SVs
- [ ] SVs recorded in truth VCF
- [ ] All 8 tests pass
- [ ] `cargo clippy` clean
