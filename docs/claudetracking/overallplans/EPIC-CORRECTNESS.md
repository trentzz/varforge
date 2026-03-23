# EPIC-CORRECTNESS: Bug Fixes and Correctness Issues

## Goal

Fix known correctness bugs discovered during the March 2026 codebase review. All bugs produce
incorrect output rather than crashes, so they are silent data-quality problems.

## Motivation

Quality averaging overflow and operator precedence errors corrupt Phred scores in edited and
structural-variant reads. The MD tag stub produces invalid BAM records. The purity override
silently takes the minimum of two values rather than preferring the CLI value.

## Scope

- In scope: the four bugs listed below, regression tests for each.
- Out of scope: performance, new features, refactoring beyond the fix itself.

## Tasks

- T052: Fix quality averaging operator precedence in structural.rs and read_modifier.rs
- T053: Fix quality averaging u8 overflow in spike_in.rs
- T054: Implement inject_snv_into_md correctly
- T055: Fix purity merge in edit.rs to CLI-wins logic
