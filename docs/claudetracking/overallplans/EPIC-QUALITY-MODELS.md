# EPIC-QUALITY-MODELS: Realistic Quality Score Profiles

## Goal

Replace placeholder uniform random quality models for PacBio HiFi and Nanopore R10 with
realistic position-dependent and context-dependent models.

## Motivation

The current models assign uniform random qualities within a range, which is unrealistic.
PacBio HiFi has a cycle-dependent profile; Nanopore has systematic errors in homopolymers.
Tools that use quality scores for filtering (e.g. fgbio, HUMID) will behave differently on
realistic vs uniform quality profiles.

## Scope

- In scope: position-dependent PacBio HiFi, homopolymer-aware Nanopore error model.
- Out of scope: CpG-specific error models, strand-specific asymmetries.

## Tasks

- T091: Implement position-dependent PacBio HiFi quality curve
- T092: Implement homopolymer-aware Nanopore R10 quality degradation
