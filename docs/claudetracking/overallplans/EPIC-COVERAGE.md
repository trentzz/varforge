# EPIC-COVERAGE: Coverage Models and QC Report

## Goal
Support per-target depth specified in the capture BED file, and produce a JSON QC report summarising coverage, family size distribution, and per-variant actual VAF after each run.

## Motivation
Real capture panels have non-uniform depth across targets. A per-target depth column lets users model this directly. The QC report closes the loop: users can verify the simulation matched the intended parameters.

## Scope
- In scope: optional depth column (col 4) in capture BED, `CaptureModel` per-target multipliers, per-chromosome coverage stats, family size distribution, `sim_report.json`.
- Out of scope: GC-bias models, insert-size-dependent coverage.

## Tasks
- T019: Parse optional depth column (col 4) from capture BED
- T020: Use stored per-target multipliers in `coverage_multiplier_at`
- T021: Track per-chromosome coverage statistics (mean, sd) during simulation
- T022: Compute family size distribution and duplicate rate
- T023: Write `sim_report.json` with coverage, VAF accuracy, family stats
