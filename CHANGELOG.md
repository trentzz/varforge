# Changelog

All notable changes to VarForge are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Version numbers follow semantic versioning: Z-bump for fixes and additions,
Y-bump for breaking changes.

---

## [Unreleased] — v0.1.1

### Added

- **Inline UMI FASTQ and BAM output** (`umi.inline: true`). The UMI sequence is
  prepended to the read sequence so tools such as fgbio `ExtractUmisFromBam` can
  strip it without a custom read-name parser.
- **Inline UMI spacer** (`umi.spacer`). An optional fixed nucleotide sequence
  appended after the inline UMI and before the template sequence. Matches the
  Twist Biosciences AT-spacer layout and any other chemistry that places a
  fixed adapter between the UMI and template.
- **Configurable duplex conversion rate** (`umi.duplex_conversion_rate`). Controls
  the fraction of duplex molecules for which both the AB and BA strands are
  recovered. The remainder produce only an AB-strand family, simulating real
  library preparation losses. Default: 1.0 (all molecules yield both strands).
- **UMI sequencing error injection** (`umi.error_rate`). Injects random base-call
  errors into UMI sequences at a configurable per-base rate. Produces near-miss
  UMI families to test the error-correction tolerance of deduplication tools
  (fgbio, HUMID, UMI-tools). Default: 0.0.
- **Updated Twist preset** (`--preset twist` / `examples/twist_duplex_benchmark.yaml`).
  Now uses 5 bp UMI, AT spacer, 90 % duplex conversion rate, and 0.1 % UMI
  error rate, matching the Twist Biosciences Comprehensive Exome Panel layout.

### Fixed

- BA strand `ref_seq` was written as the forward-strand sequence instead of
  the reverse complement. This caused incorrect base qualities and mismatches
  for BA-strand reads in BAM output. BA-strand reads now carry the correct
  reverse-complemented reference sequence.
- BAM R2 soft-clip position was calculated from the wrong read end, shifting
  the soft-clip CIGAR operation by the clip length. BAM R2 records now have
  the soft-clip at the correct position.

---

## [0.1.0] — 2025-03-01

Initial public release.

### Features

- Somatic SNV, indel, and MNV simulation with configurable VAF ranges.
- Structural variants: DEL, INS, INV, DUP, TRA with HRD, TDP, and
  CHROMOTHRIPSIS signatures.
- COSMIC SBS signature-weighted base selection.
- Tumour purity and multi-clone subclonal architecture.
- Paired tumour-normal simulation.
- Germline SNP/indel simulation with truth VCF output.
- cfDNA fragment model with mononucleosomal and dinucleosomal peaks and
  end-motif rejection sampling.
- Long-read fragment model (log-normal sampler).
- Duplex UMI barcodes with PCR amplification families.
- FFPE deamination and OxoG artefact simulation.
- Hybrid-capture and amplicon enrichment model.
- GC bias model.
- Copy number alteration with depth scaling.
- Microsatellite instability (MSI) mode.
- Multi-sample / longitudinal series mode.
- Empirical quality profile learning from real BAM files.
- YAML configuration with variable substitution and CLI overrides.
- Named presets for common scenarios and cancer types.
- Truth VCF, BAM, FASTQ, and manifest outputs.
- Pure Rust dependency stack; single binary, no C libraries required.
