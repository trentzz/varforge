# EPIC-MD-TAG: MD Tag Generation

## Goal

Generate correct MD tags on all simulated BAM records so that tools like VarDict,
samtools view, and Picard ValidateSamFile accept the output without error.

## Motivation

VarDict requires the MD tag for variant calling. Picard ValidateSamFile flags records missing
MD as MISSING_TAG_NM_MD. Currently all simulated BAM records lack MD entirely.

## Scope

- In scope: MD tag for reference-matching bases and spiked variants in the simulate path.
- Out of scope: the editor path inject_snv_into_md fix (covered by T054).

## Tasks

- T089: Generate MD tag for simulated BAM reads (reference-match runs + spiked alt bases)
