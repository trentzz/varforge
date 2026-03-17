# Task 22: BAM Editor (Spike-in to Existing BAM)

## Phase
3 - Extended Variant Support (parallel to from-scratch simulation)

## Dependencies
Task 05 (Read Generation Engine - for variant spike-in logic)

## Objective
Implement a BAM editing mode that takes an existing BAM file and spikes in SNVs, indels, and SVs while preserving the original sequencing artifacts, error profiles, and read qualities. This is VarForge's answer to BAMSurgeon -- but faster (Rust), with stochastic VAF sampling, and supporting UMI-aware spike-in.

## Context
Many users have real sequencing data and want to add known variants for benchmarking. By editing a real BAM, all the natural sequencing characteristics are preserved:
- Real base quality distributions and error patterns
- Real fragment size distribution
- Real PCR duplicates and library artifacts
- Real GC bias and coverage variation
- Platform-specific systematic errors

BAMSurgeon does this but is very slow (~1 sec/variant), uses deterministic VAF (not stochastic), and has no UMI awareness. stochasticSim addresses VAF but only handles SNVs and is exome-only.

## Key Files
- **Create**: `src/editor/mod.rs`
- **Create**: `src/editor/bam_editor.rs` - Core BAM editing engine
- **Create**: `src/editor/read_modifier.rs` - Read-level modification logic
- **Create**: `src/cli/edit.rs` - `edit` subcommand
- **Modify**: `src/cli/mod.rs` (add `edit` subcommand)
- **Modify**: `src/main.rs` or `src/lib.rs` (add `editor` module)
- **Existing context**: `src/variants/spike_in.rs`, `src/variants/vaf.rs`, `src/io/truth_vcf.rs`

## Requirements

### CLI: `edit` subcommand
```
varforge edit \
  --bam input.bam \
  --reference hg38.fa \
  --variants mutations.vcf \
  --output output.bam \
  --seed 42 \
  --threads 4 \
  --purity 0.7 \
  --truth-vcf truth.vcf
```

Or with YAML config:
```
varforge edit --config edit_config.yaml
```

### Edit Config (YAML)
```yaml
mode: edit
input_bam: /path/to/input.bam
reference: /path/to/hg38.fa
output_bam: /path/to/output.bam
variants:
  vcf_path: /path/to/mutations.vcf
  # OR inline:
  # list:
  #   - { chrom: "chr7", pos: 55249071, ref: "T", alt: "A", vaf: 0.25 }  # EGFR T790M
tumour:
  purity: 0.7
  ploidy: 2
seed: 42
truth_vcf: /path/to/truth.vcf
```

### Core Algorithm

For each variant to spike in:

1. **Collect overlapping reads**: Query the input BAM for all reads overlapping the variant position
2. **Calculate target alt count**: Use binomial sampling based on expected VAF and total read depth at position
3. **Select reads to modify**: Randomly select `alt_count` reads to carry the variant, respecting:
   - Strand balance (proportional to existing F/R ratio)
   - Mapping quality threshold (skip MQ=0 reads)
   - Base quality threshold (skip BQ<20 at variant position)
4. **Apply modification**:
   - **SNV**: Change the base at the variant position, keep original quality score
   - **Indel**: Modify sequence and CIGAR string, adjust alignment
   - **SV**: For split reads at breakpoints, soft-clip and remap; for discordant pairs, adjust mate position/orientation
5. **Write output**: Copy all reads, with modified reads replacing originals
6. **Record truth**: Write applied variant with actual alt/total counts to truth VCF

### SNV Editing
- Replace base in SEQ field at the variant position
- Keep original QUAL score (preserves real quality profile)
- Update MD tag if present
- Update NM tag (edit distance)

### Indel Editing
- **Insertion**: Insert bases into SEQ and QUAL, update CIGAR (add I operation), adjust downstream alignment
- **Deletion**: Remove bases from SEQ and QUAL, update CIGAR (add D operation), adjust downstream positions
- Recalculate MD and NM tags
- For reads near the indel that extend past it: soft-clip if realignment would be complex

### SV Editing
- **Deletion (>50bp)**: For reads spanning the breakpoint, soft-clip at the breakpoint. Adjust mate pairs. Remove reads that fall entirely within the deletion (for tumor allele).
- **Insertion (>50bp)**: Insert sequence into spanning reads, soft-clip excess.
- **Inversion**: Reverse-complement the sequence in reads within the inverted region (for tumor allele reads).
- **Translocation**: Create chimeric reads at breakpoints with soft-clipping. Adjust supplementary alignments (SA tag).

### UMI-Aware Editing
When input BAM has UMI tags (RX):
- Group reads by UMI family before selecting which to modify
- Modify entire families together (a true variant should appear in all family members)
- This ensures downstream UMI consensus calling sees the variant correctly
- Artifactual variants (for testing error correction) can be spiked into individual reads within a family

### Stochastic VAF (key differentiator vs BAMSurgeon)
- Target VAF → expected alt count = VAF × total depth
- Actual alt count sampled from Binomial(total_depth, VAF)
- This produces realistic VAF scatter matching real sequencing

### Parallel Processing
- Process variants in parallel using rayon
- Sort variants by position to batch overlapping regions
- Use indexed BAM (BAI) for random access
- Thread-safe output via sorted merge of modified regions

## Tests

### Unit Tests
1. `test_snv_spike_in_bam` - Spike SNV, verify base changed at position, quality preserved
2. `test_snv_stochastic_vaf` - Spike at 50% VAF, verify alt fraction ~50% (statistical)
3. `test_indel_insertion_bam` - Spike insertion, verify CIGAR updated, sequence extended
4. `test_indel_deletion_bam` - Spike deletion, verify CIGAR updated, sequence shortened
5. `test_md_tag_updated` - MD tag correctly reflects the modification
6. `test_nm_tag_updated` - NM tag (edit distance) updated
7. `test_strand_balance` - Modified reads respect original strand distribution
8. `test_quality_preserved` - Base qualities not altered (except at modified positions for indels)
9. `test_umi_family_editing` - Entire UMI family modified together for true variants
10. `test_umi_artifact_editing` - Single read modified within family for artifact simulation
11. `test_truth_vcf_output` - Truth VCF contains all spiked variants with actual counts
12. `test_unmodified_reads_unchanged` - Reads not selected for modification are byte-identical to input
13. `test_deterministic_with_seed` - Same seed produces identical output

### Integration Tests
14. `test_edit_command_e2e` - Full `varforge edit` command produces valid output BAM
15. `test_edit_vs_simulate_consistency` - Variant at same VAF produces similar alt fractions in both modes

## Acceptance Criteria
- [ ] `varforge edit --bam input.bam --variants mutations.vcf` produces modified BAM
- [ ] SNV, indel spike-in with stochastic VAF
- [ ] Original read qualities and artifacts preserved
- [ ] CIGAR, MD, NM tags correctly updated
- [ ] UMI-aware family editing
- [ ] Truth VCF with actual applied counts
- [ ] All 15 tests pass
- [ ] `cargo clippy` clean
- [ ] Faster than BAMSurgeon for equivalent operations
