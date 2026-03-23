# EPIC-BATCH: Multi-VAF Batch Mode and Dry-Run

## Goal
Allow users to sweep across multiple VAFs in a single command, and to preview expected molecule counts before running a full simulation.

## Motivation
Benchmarking tools requires datasets at many VAF levels. Manually cloning configs and running each is error-prone. A batch mode and a dry-run preview reduce friction and catch misconfigured runs before expensive I/O.

## Scope
- In scope: `vafs` list field in config, per-VAF output subdirectories, batch manifest TSV, `--dry-run` flag, expected count table.
- Out of scope: parallelism across batch runs (handled by existing rayon workers), GUI or interactive mode.

## Tasks
- T006: Add `vafs: Option<Vec<f64>>` field to `Config`
- T007: Implement per-VAF config forking in simulate command
- T008: Write batch manifest listing all VAF runs
- T009: Add `--dry-run` flag to simulate subcommand
- T010: Compute and print expected alt molecule count table
