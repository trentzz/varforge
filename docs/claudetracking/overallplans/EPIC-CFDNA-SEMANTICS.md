# EPIC-CFDNA-SEMANTICS: cfDNA Model Correctness

## Goal

Fix the inverted ctDNA fraction parameter and add explicit cfDNA config fields so the
fragment model correctly represents liquid biopsy biology.

## Motivation

`ctdna_fraction = 1.0 - purity` is inverted: higher tumour burden means more ctDNA, not less.
This makes cfDNA simulations produce the opposite of intended VAF distributions for ctDNA
callers. The fix is a dedicated `cfdna.ctdna_fraction` field.

## Scope

- In scope: ctDNA fraction parameter, dedicated config field, SD override for fragment peaks.
- Out of scope: methylation patterns, fragment end motif enrichment (in EPIC-ADVANCED).

## Tasks

- T076: Fix inverted ctDNA fraction in cfDNA fragment sampler
- T077: Add explicit `cfdna.ctdna_fraction` config field
- T078: Allow user override of cfDNA fragment peak SDs; add source citation or reference
