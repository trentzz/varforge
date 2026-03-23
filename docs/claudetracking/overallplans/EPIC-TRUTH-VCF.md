# EPIC-TRUTH-VCF: Truth VCF Completeness

## Goal

Make the truth VCF output complete and usable for benchmarking: correct genotypes,
accurate CCF, bgzip+tabix, and SV records.

## Motivation

Hardcoded `GT=0/1` is wrong for homozygous variants. Missing SV records mean the truth set
is incomplete for SV callers. Plain-text VCF is not accepted by bcftools or IGV without
manual compression.

## Scope

- In scope: genotype correctness, CCF accuracy, compression, SV records.
- Out of scope: GVCF format, multi-sample VCF merging.

## Tasks

- T079: Support homozygous genotypes (1/1) for variants with VAF ≥ 0.99
- T080: Compute CCF from ClonalTree once T067 is integrated (depends on T067)
- T081: Write bgzip-compressed truth VCF with tabix index
- T082: Include SV records in truth VCF (wire sv_vcf_info / sv_vcf_alt into writer)
