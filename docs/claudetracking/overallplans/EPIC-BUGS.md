# EPIC-BUGS: Bug Fixes

## Goal
Resolve known crashes and correctness failures in the simulation engine. All tasks in this epic are already complete.

## Motivation
SV boundary panics and incorrect duplex read generation produce invalid output and crash the simulator on real inputs. These must be fixed before other features are layered on top.

## Scope
- In scope: SV boundary panics, duplex strand logic, regression tests for both.
- Out of scope: new features, performance work.

## Tasks
- T001: Fix SV boundary panic: skip too-short post-deletion fragments
- T002: Fix SV boundary panic: clamp translocation partner-seq fetch to chromosome length
- T003: Add regression test for SV boundary panic
- T004: Fix duplex: implement BA-strand complement read pair
- T005: Add unit test: duplex simulation produces both AB and BA reads
