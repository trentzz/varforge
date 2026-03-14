# Liquid Biopsy Requirements

Simulation of cell-free DNA (cfDNA) and circulating tumour DNA (ctDNA) sequencing data.

## Nucleosomal Fragment Size Distribution

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-LBX-001 | Model cfDNA fragment sizes as a mixture of mononucleosomal (~167 bp) and dinucleosomal (~334 bp) peaks | P0 |
| REQ-LBX-002 | Include 10 bp periodicity sub-peaks reflecting the helical pitch of DNA wound around nucleosomes | P0 |
| REQ-LBX-003 | Configurable relative weights of mononucleosomal vs. dinucleosomal components | P1 |

## ctDNA Fragment Enrichment

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-LBX-010 | Tumour-derived (ctDNA) fragments use a shorter size distribution centered around ~143 bp | P0 |
| REQ-LBX-011 | Normal cfDNA fragments use the standard nucleosomal distribution (~167 bp peak) | P0 |
| REQ-LBX-012 | The size difference between ctDNA and normal cfDNA should be configurable | P1 |

## Low Tumour Fraction Support

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-LBX-020 | Support tumour fractions as low as 0.001 (0.1%) | P0 |
| REQ-LBX-021 | Correctly model expected VAFs at low tumour fractions (e.g., a clonal heterozygous SNV at 0.5% TF yields ~0.25% VAF) | P0 |
| REQ-LBX-022 | Generate sufficient total depth to make low-frequency variants detectable (e.g., 1000× for 0.5% VAF) | P0 |

## Multi-Sample Longitudinal Series

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-LBX-030 | Generate multiple samples sharing the same clonal tree but with different tumour fractions | P1 |
| REQ-LBX-031 | Accept an array of tumour fractions (e.g., [0.05, 0.02, 0.01, 0.005]) representing serial timepoints | P1 |
| REQ-LBX-032 | Each sample in a series uses the same mutation set and clonal structure, only TF varies | P1 |
| REQ-LBX-033 | Allow per-timepoint clone emergence or extinction (optional subclone CCF changes) | P2 |
