# EPIC-ADVANCED: Advanced Tumour Biology

## Goal
Simulate advanced biological phenomena: haplotype phasing, complex SV signatures (HRD, TDP, chromothripsis), cross-sample contamination, cfDNA fragment end motifs, and microsatellite instability.

## Motivation
These features cover the long tail of real tumour biology that basic SNV/indel simulation cannot reproduce. They are needed to benchmark specialised callers and to produce publication-quality synthetic datasets.

## Scope
- In scope: haplotype-aware spike-in, HRD/TDP/chromothripsis SV generators, SNP-profile and BAM-source contamination, 4-mer end motif rejection sampling, MSI indel enrichment, hg38 microsatellite loci BED.
- Out of scope: clonal evolution modelling, spatial heterogeneity, epigenetic marks.

## Tasks
- T036: Add `haplotype: Option<u8>` to `Variant`; implement haplotype-aware spike-in
- T037: Implement HRD deletion pattern generator (100 kbp–10 Mbp deletions with LOH)
- T038: Implement tandem duplicator phenotype (enriched short tandem dups 1–10 kbp)
- T039: Implement chromothripsis model (random cuts + re-ligation on one chromosome)
- T040: Add `contamination` config block; implement synthetic contamination from SNP profile VCF
- T041: Implement BAM-source contamination (randomly subsample reads from donor BAM)
- T042: Add 4-mer end motif frequency table; implement rejection sampling in fragment generation
- T048: Add `tumour.msi: true`; implement microsatellite-enriched indel generation
- T049: Bundle hg38 microsatellite loci BED as an optional static asset
