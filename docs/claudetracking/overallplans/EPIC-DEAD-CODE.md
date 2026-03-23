# EPIC-DEAD-CODE: Dead Code Reduction

## Goal

Reduce the 67 `#[allow(dead_code)]` annotations to only those covering planned future work,
and remove genuinely unused items.

## Motivation

Dead code annotations hide real problems (clippy would catch regressions), increase binary
size, and create confusion about what is and is not part of the active API. Items tied to
planned tasks should be annotated with the task ID, not silently suppressed.

## Scope

- In scope: audit and remove unused items; annotate kept items with task references.
- Out of scope: refactoring the API surface (separate from dead code removal).

## Tasks

- T103: Audit all dead_code annotations; remove genuinely unused items; annotate kept items
  with task IDs
- T104: Wire PcrFamilySizeSampler into UMI family generation (log-normal family sizes)
