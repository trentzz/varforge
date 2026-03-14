# Artifact Simulation Requirements

Realistic library preparation and sequencing artifacts.

## FFPE Damage

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-ART-001 | Simulate formalin-induced cytosine deamination (C>T / G>A) at a configurable rate | P0 |
| REQ-ART-002 | FFPE damage is strand-specific: C>T on the affected strand, observed as G>A on the complement | P0 |
| REQ-ART-003 | Damage rate correlates with fragment ends (higher near break sites) | P1 |

## Oxidative Damage (oxoG)

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-ART-010 | Simulate 8-oxoguanine artifacts (G>T / C>A) at a configurable rate | P1 |
| REQ-ART-011 | oxoG artifacts appear predominantly on one strand (read 1 vs. read 2 asymmetry) | P1 |

## GC Bias

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-ART-020 | Model GC-dependent coverage variation using a configurable bias curve | P1 |
| REQ-ART-021 | Support a lookup table mapping GC% bins to depth multipliers | P1 |
| REQ-ART-022 | Default bias curve derived from typical Illumina whole-genome library prep | P2 |

## PCR Duplicates

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-ART-030 | Generate PCR duplicate reads sharing identical alignment coordinates | P0 |
| REQ-ART-031 | Duplicate rate configurable as a fraction of total reads (e.g., 0.15 for 15% duplication) | P0 |
| REQ-ART-032 | In UMI mode, PCR duplicates share the same UMI (see UMI requirements) | P0 |

## PCR Error Injection

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-ART-040 | Simulate PCR polymerase errors introduced during amplification | P1 |
| REQ-ART-041 | PCR errors accumulate per cycle: early-cycle errors propagate to more copies than late-cycle errors | P1 |
| REQ-ART-042 | Configurable per-cycle error rate (default ~1e-5 per bp per cycle) | P1 |
| REQ-ART-043 | PCR errors within a UMI family are shared across duplicates descended from the error-carrying copy | P1 |
