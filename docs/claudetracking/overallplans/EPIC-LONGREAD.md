# EPIC-LONGREAD: Long-Read Sequencing Support

## Goal
Support PacBio HiFi and Oxford Nanopore R10 sequencing platforms, including appropriate fragment length distributions, quality models, and single-read BAM output.

## Motivation
Long-read sequencing is increasingly used for structural variant detection and phasing. VarForge must simulate these platforms to be useful for benchmarking long-read callers.

## Scope
- In scope: `LongReadFragmentSampler`, PacBio HiFi quality model, Nanopore R10 quality model, single-read BAM output, platform detection wiring.
- Out of scope: CCS consensus modelling, basecaller-specific error distributions, multi-pass reads.

## Tasks
- T024: Add `LongReadFragmentSampler` (log-normal, 5–50 kbp)
- T025: Add PacBio HiFi quality model (Q20–Q30, low indel rate)
- T026: Add Nanopore R10 quality model (Q15–Q25, homopolymer errors)
- T027: Single-read BAM output mode (no R2, template length = read length)
- T028: Wire platform detection into engine and fragment/quality selection
