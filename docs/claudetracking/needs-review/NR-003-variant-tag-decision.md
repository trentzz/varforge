# VariantTag: annotate FASTQ names or remove?

**Category**: feature
**Related epic**: EPIC-POLISH

## Context

The `VariantTag` struct is populated for every variant-carrying read pair in the
engine, but no output path consumes it. The struct and its population code add
complexity and allocations. Two paths forward:

## Options

1. **Annotate FASTQ read names**: append variant info to read names (e.g.
   `@readname:VAR:chr1:1000:SNV:0.15`). This lets users grep for variant reads
   without the truth VCF. Useful for debugging and lightweight validation.
   Tradeoff: makes read names longer, which increases FASTQ file size by ~5%.

2. **Remove VariantTag entirely**: delete the struct, the population code in
   `generate_region`, and the `variant_tags` field on `ReadPair`. Reduces
   complexity and per-read allocations. Tradeoff: loses the ability to
   annotate individual reads with their variant status.

## Recommendation

Option 1. Annotating read names is a low-cost way to add significant
debugging value. The 5% file size increase is negligible for a simulation tool
where output size is already configurable via coverage.
