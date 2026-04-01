# EPIC-ERROR-MODEL: Sequencing Error and Noise Model (v0.2.0)

## Goal

Give VarForge a layered, composable sequencing error model that makes generated
reads indistinguishable from real instrument output in terms of error structure.
All new features live under `quality.sequencing_errors` in the YAML config and
are fully opt-in — absent config is identical to v0.1.x behaviour.

## Separation of concerns

`artifacts:` models sample preparation artifacts:
- FFPE deamination, oxoG damage, PCR duplicates, PCR substitution errors

`quality.sequencing_errors:` models instrument artifacts:
- Base-call errors (substitutions), sequencing indels, cycle-position effects,
  context-dependent error elevation, strand asymmetry, phasing bursts

## Architecture

All new error models are composed by an `ErrorOrchestrator` struct that runs
passes in a defined order per read pair:

1. Quality-score-driven substitutions (existing, unchanged)
2. Cycle-position additional substitutions (CycleErrorCurve)
3. Context-dependent multipliers (KmerErrorModel)
4. Sequencing indels (IndelErrorModel)
5. Strand bias quality offset to R2 (StrandBiasModel)
6. Correlated phasing burst errors (CorrelatedErrorModel)

The orchestrator is constructed once per region and shared across all read pairs
in the batch. All precomputed tables (cycle curve, k-mer hash) are built at init
time, keeping the per-read hot path allocation-free.

## What success looks like

- A simulated Illumina NovaSeq run shows: substitution error rate ~0.1%,
  indel rate ~0.005%, exponential error rise in the last 20% of cycles,
  elevated errors at poly-G and GGC contexts, R2 error rate 1.2–1.5× R1.
- All new features are off by default; existing tests pass unchanged.
- `varforge learn-profile` outputs a profile JSON that, when loaded, drives
  the full orchestrator automatically.
- Platform presets auto-populate sensible defaults for each chemistry.

## Tasks

- T151: IndelErrorModel — per-base sequencing indel injection
- T152: CycleErrorCurve — flat and exponential per-cycle error rate
- T153: KmerErrorModel — context-dependent error multipliers
- T154: StrandBiasModel + CorrelatedErrorModel — R2 asymmetry and phasing bursts
- T155: ErrorOrchestrator + SequencingErrorConfig + preset defaults
- T156: ProfileLearner extension — learn indel and cycle error rates from BAMs
