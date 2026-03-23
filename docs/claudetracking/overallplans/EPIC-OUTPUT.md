# EPIC-OUTPUT: Output Improvements

## Goal
Enrich output files so downstream tools can verify which reads carry mutations and how many alt molecules were observed per variant.

## Motivation
Truth VCF INFO fields and per-read variant tags enable direct benchmarking of variant callers without requiring a separate ground-truth pipeline. They make VarForge outputs self-documenting.

## Scope
- In scope: variant tags on read names, `variant_reads.tsv` sidecar, `N_ALT_MOL` and `N_DUPLEX_ALT` INFO fields in truth VCF.
- Out of scope: BAM tag injection (separate feature), CRAM output.

## Tasks
- T011: Aggregate N_ALT_MOL and N_DUPLEX_ALT per variant across regions
- T012: Write `N_ALT_MOL` and `N_DUPLEX_ALT` INFO fields to truth VCF
- T013: Add `variant_tags: Vec<VariantTag>` to `ReadPair`; populate during spike-in
- T014: Append variant tag string to FASTQ read name on write
- T015: Implement `variant_reads.tsv` sidecar writer
