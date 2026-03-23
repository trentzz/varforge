# EPIC-GERMLINE: Germline Variants and Paired Tumour-Normal Mode

## Goal
Add germline SNP/indel simulation and a paired tumour-normal output mode where both samples share a germline background and the tumour additionally carries somatic mutations.

## Motivation
Most real variant calling workflows require a matched normal. Simulating paired samples enables end-to-end benchmarking of somatic callers that subtract germline signal.

## Scope
- In scope: `germline` config section, population-frequency SNP/indel generator, `germline_truth.vcf`, `mode: paired-tumour-normal`, `PairedSimulationMode`.
- Out of scope: structural germline variants, population stratification, trio simulation.

## Tasks
- T032: Add `germline` config section; implement population-frequency SNP/indel generator
- T033: Write `germline_truth.vcf` output
- T034: Add `mode: paired-tumour-normal` config; define paired config schema
- T035: Implement `PairedSimulationMode`: shared germline, split somatic spike-in
