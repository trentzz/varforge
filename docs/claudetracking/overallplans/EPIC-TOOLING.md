# EPIC-TOOLING: Developer Tooling

## Goal
Add developer and power-user tooling: a benchmarking suite CLI, amplicon sequencing support, tumour-in-normal contamination, config templating, and empirical error profile learning from BAM.

## Motivation
These features reduce friction for users running large-scale evaluations and for developers tuning simulation parameters from real data. They are lower priority than core biology features but increase VarForge's utility as a benchmarking platform.

## Scope
- In scope: `benchmark-suite` subcommand, amplicon fixed-length fragments, primer trimming, `normal.tumour_contamination`, `${variable}` substitution in config, `--set` CLI flag, `learn-profile` command, `.profile.tsv` serialisation.
- Out of scope: web UI, cloud execution, result aggregation dashboards.

## Tasks
- T043: Add `benchmark-suite` subcommand: VAF × coverage config grid + manifest
- T044: Add `capture.mode: amplicon`; implement fixed-length fragment sampling at amplicon coords
- T045: Implement `primer_trim` option for amplicon reads
- T046: Add `normal.tumour_contamination` field; apply somatic mutations at low level in normal simulation
- T047: Implement `${variable}` substitution in Config loading + `--set key=value` CLI flag
- T050: Complete `learn-profile` command: read quality-by-cycle from BAM
- T051: Serialise/deserialise quality profile to `.profile.tsv`
