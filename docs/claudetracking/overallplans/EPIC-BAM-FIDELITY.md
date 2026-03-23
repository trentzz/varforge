# EPIC-BAM-FIDELITY: BAM Output Correctness

## Goal

Fix BAM output issues that make generated files incompatible with downstream tools
(samtools, Picard, Mutect2, VarDict). All five bugs produce silently wrong BAM data.

## Motivation

Incorrect mate positions, missing NM tags, un-reversed R2 sequences, and hardcoded MAPQ
mean VarForge BAMs may fail validation or produce wrong variant calls in the very tools
VarForge is designed to benchmark.

## Scope

- In scope: FLAG correctness, mate position, NM tag, R2 orientation, MAPQ, TLEN overflow.
- Out of scope: full CIGAR generation for indel-containing reads (separate task).

## Tasks

- T071: Fix R2 mate position calculation in bam.rs
- T072: Add NM (edit distance) tag to generated BAM records
- T073: Reverse-complement R2 sequence and quality in BAM output
- T074: Guard template_length i32 overflow for long-read fragments
- T075: Wire configurable or model-based MAPQ (not hardcoded 60)
