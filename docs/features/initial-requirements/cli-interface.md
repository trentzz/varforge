# CLI Interface Requirements

Command-line interface design.

## Subcommands

### `varforge simulate`

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-CLI-001 | Accept a YAML config file path as the primary argument | P0 |
| REQ-CLI-002 | Override individual config values via CLI flags (e.g., `--coverage 100`, `--seed 42`) | P1 |
| REQ-CLI-003 | `--dry-run` flag: validate config and report estimated output size without generating data | P1 |
| REQ-CLI-004 | `--random-mutations N` to auto-generate N random mutations without an input VCF | P0 |
| REQ-CLI-005 | `--vaf-range MIN-MAX` to control VAF range for random mutations (e.g., `0.0001-0.05`) | P0 |
| REQ-CLI-006 | `--preset {small,panel,wgs}` for quick common configurations | P1 |

### `varforge learn-profile`

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-CLI-010 | Accept a BAM file and output a learned error profile file | P1 |
| REQ-CLI-011 | `--sample-reads N` to limit the number of reads analyzed (default: 1 million) | P1 |
| REQ-CLI-012 | `--regions BED` to restrict profiling to specific genomic regions | P2 |

### `varforge validate`

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-CLI-020 | Accept a YAML config and report validation results (schema correctness, reference accessibility, region overlap, mutation conflicts) | P1 |
| REQ-CLI-021 | Exit with non-zero status on validation failure | P1 |

## Common Options

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-CLI-030 | `--threads N` to control parallelism (default: number of available cores) | P0 |
| REQ-CLI-031 | `--output-dir DIR` to override the output directory from config | P0 |
| REQ-CLI-032 | `--quiet` / `--verbose` flags to control log verbosity | P0 |

## Progress Reporting

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-CLI-040 | Display a progress bar or percentage when writing to a terminal (stderr) | P1 |
| REQ-CLI-041 | Machine-readable progress (JSON lines) when `--progress json` is specified | P2 |

## Logging

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-CLI-050 | Use structured logging (e.g., `tracing` crate) with levels: error, warn, info, debug, trace | P0 |
| REQ-CLI-051 | Default log level: info; configurable via `--log-level` or `VARFORGE_LOG` env var | P0 |
