# EPIC-REFACTOR: Code Quality and Deduplication

## Goal

Remove duplicated logic, dead code, and structural problems identified in the March 2026 review.
No behaviour changes — this is purely internal quality work.

## Motivation

Duplicate `complement`/`reverse_complement` functions, three separate free functions extracting
variant position, a dead no-op wrapper, and near-identical BAM write methods all increase
maintenance cost and create divergence risk.

## Scope

- In scope: deduplication, dead code removal, consolidation of parallel implementations.
- Out of scope: new features, performance work, API changes visible to users.

## Tasks

- T056: Extract seq_utils module with complement and reverse_complement
- T057: Unify variant position and type helpers as methods on MutationType
- T058: Remove noodles_base_to_u8 dead wrapper
- T059: Unify write_pair and write_pair_with_umi in bam.rs
