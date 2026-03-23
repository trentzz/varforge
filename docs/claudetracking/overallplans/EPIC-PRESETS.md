# EPIC-PRESETS: Chemistry Presets

## Goal
Allow users to select a named chemistry preset (e.g. `illumina-umi`, `pacbio-hifi`) that fills in sensible defaults without requiring every parameter to be specified in the YAML config.

## Motivation
New users should be able to get a working config with one line. Presets also act as canonical reference points for chemistry-specific parameter sets.

## Scope
- In scope: preset lookup table in `io/config.rs`, application of preset defaults after YAML deserialisation, integration test.
- Out of scope: user-defined custom presets, preset inheritance.

## Tasks
- T016: Define preset lookup table in `io/config.rs`
- T017: Apply preset defaults after YAML deserialisation
- T018: Add integration test: preset fills expected defaults
