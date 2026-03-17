# Feature: BED Region Filtering

**Status:** Implemented (v0.1.0)
**Author:** Trent Zeng
**Date:** 2026-03-18

---

## Motivation

Targeted sequencing panels (Twist, Agilent, IDT) cover a small fraction of the genome. Without region filtering, simulating a 500x panel on a full chromosome would require simulating billions of read pairs across non-target sequence, then discarding them. This is orders of magnitude more work than necessary and impractical on a workstation.

BED region filtering restricts simulation to intervals that overlap a supplied target BED file. Only regions that intersect at least one BED interval are simulated. This makes high-coverage panel simulation feasible and produces output that reflects real targeted sequencing data.

---

## Config Schema

Add the `regions_bed` field to the top-level config:

```yaml
regions_bed: /path/to/targets.bed
```

The path can be absolute or relative to the current working directory. The field is optional; omitting it simulates the full reference (previous behaviour, unchanged).

The BED file must have at least three tab-separated fields per line: `chrom`, `start` (0-based), `end`. Comment lines starting with `#`, `track`, or `browser` are skipped. Both tab and space delimiters are accepted as a fallback.

---

## Behaviour

1. VarForge partitions the reference into fixed-size chunks with `partition_regions()`.
2. If `regions_bed` is set, the chunks are intersected with the BED targets using `intersect_with_targets()`. Each retained chunk is clipped to the boundary of its overlapping target.
3. If the intersection is empty (e.g., chromosome names in the BED do not match the reference), VarForge exits with an error rather than silently producing empty output.
4. Only the retained, clipped chunks are dispatched to workers. Coverage depth applies over the simulated region length, not the full chromosome.

---

## Example

```yaml
reference: /data/hg38/chr22.fa
regions_bed: /data/panels/twist_exome_plus_chr22.bed

sample:
  name: panel_test
  read_length: 150
  coverage: 500.0

output:
  directory: output/panel_test
```

With a 51 Mbp chromosome and a 1 Mbp panel footprint, this simulates roughly 1/51 of the chunks and completes ~51x faster than a full-chromosome run at the same coverage.

---

## Implementation

All required infrastructure already existed before this feature was enabled:

| Component | Location | Notes |
|-----------|----------|-------|
| Config field | `src/io/config.rs:30` | `regions_bed: Option<PathBuf>` |
| BED parser | `src/cli/simulate.rs:958` | `parse_bed_file()` |
| Intersection | `src/core/coverage.rs:34` | `intersect_with_targets()` |
| Wiring | `src/cli/simulate.rs:195` | Inserted after `partition_regions()` |
| Validation | `src/io/config.rs:validate()` | Checks file exists; runtime check for empty intersection |

---

## Limitations

- No padding around targets. Reads that start before a target boundary but extend into it are not simulated. If overhanging reads are required, add padding to the BED file before running.
- Variant placement respects BED boundaries only indirectly: variants outside the retained chunks are not simulated. Variants inside retained chunks are simulated normally.
- Multi-sample mode applies BED filtering independently to each sample (each shares the same `regions_bed` from the top-level config).

---

## Future Work

- Per-sample BED overrides for multi-sample configs (different capture kits per sample).
- BED padding option (`regions_bed_padding: 50`) to include flanking sequence.
- Report the total targeted base count in the manifest.
