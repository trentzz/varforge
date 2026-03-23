# EPIC-MULTI-SAMPLE: Multi-Sample Mode Improvements

## Goal

Make multi-sample (longitudinal) mode parallel and biologically consistent across samples
by sharing germline variants.

## Motivation

Sequential sample execution wastes CPU on multi-core machines. Each sample generating its
own independent germline set makes longitudinal tumour evolution analysis invalid — the same
patient should carry the same germline across all timepoints.

## Scope

- In scope: parallel sample execution, shared germline variant set.
- Out of scope: full tumour evolution modelling (clonal dynamics across timepoints).

## Tasks

- T098: Parallelise multi-sample runs across samples using rayon
- T099: Share germline variant set across all samples in a multi-sample run
