# EPIC-TESTS: Test Coverage Gaps

## Goal

Add unit and integration tests for modules that currently have zero coverage.

## Motivation

bam_editor.rs, the edit/validate/benchmark-suite/learn-profile CLI subcommands, config variable
substitution edge cases, and ffpe boundary inputs are all untested. These are real execution
paths that carry bugs silently.

## Scope

- In scope: unit tests for each listed module, integration tests for each CLI subcommand.
- Out of scope: property-based testing, fuzzing (future work).

## Tasks

- T060: Tests for src/editor/bam_editor.rs
- T061: Integration tests for cli/edit.rs end-to-end
- T062: Tests for cli/validate.rs, cli/benchmark_suite.rs, cli/learn_profile.rs
- T063: Negative and edge-case tests for config variable substitution
- T064: Boundary tests for ffpe.rs end_enrichment
