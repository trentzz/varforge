# Tumour Model Requirements

Defines tumour architecture: purity, clonal structure, ploidy, and copy number state.

## Purity

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-TUM-001 | Accept a tumour purity value (0.0–1.0) controlling the fraction of reads from tumour vs. normal | P0 |
| REQ-TUM-002 | Mix tumour and normal reads at the fragment level so that coverage statistics are realistic | P0 |

## Clonal Architecture

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-TUM-010 | Define a clonal tree as a set of clones with parent-child relationships | P0 |
| REQ-TUM-011 | Assign a cancer cell fraction (CCF) to each clone (root clone CCF = 1.0, subclones < parent) | P0 |
| REQ-TUM-012 | Validate that child clone CCFs do not exceed parent CCF | P0 |
| REQ-TUM-013 | Assign each mutation to exactly one clone; mutations are inherited by all descendant clones | P0 |
| REQ-TUM-014 | Support at least 10 clones per simulation | P1 |

## Ploidy and Copy Number

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-TUM-020 | Accept a genome-wide base ploidy for the tumour (default 2) | P0 |
| REQ-TUM-021 | Support per-locus copy number states defined as (total CN, minor CN) pairs | P0 |
| REQ-TUM-022 | Allow different CN profiles per clone | P1 |
| REQ-TUM-023 | Use CN state to modulate both read depth and expected VAF via the VAF equation | P0 |

## Normal Contamination

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-TUM-030 | Generate normal reads at the complement of tumour purity (normal fraction = 1 − purity) | P0 |
| REQ-TUM-031 | Normal reads carry reference alleles at all somatic variant sites | P0 |
| REQ-TUM-032 | Allow injection of germline variants into both tumour and normal reads | P1 |
