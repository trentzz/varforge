# Changelog

All notable changes to VarForge are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Version numbers follow semantic versioning: Z-bump for fixes and additions,
Y-bump for breaking changes.

---

## [Unreleased] — v0.2.0

### Added

- **Sequencing error model (`ErrorOrchestrator`)**. A unified error-injection
  pipeline replaces the old flat quality model. It composes all error sources
  in a single pass: base substitutions, cycle-position decay, k-mer context
  modifiers, sequencing indels, strand bias, and phasing bursts.
- **Cycle-position error rates**. Three cycle error curve models: `flat`
  (uniform rate), `exponential` (tail rise), and `custom` (user-supplied TSV).
  Configurable via `quality.sequencing_errors.cycle_error_model`.
- **k-mer context errors**. Context-dependent error multipliers keyed on
  surrounding sequence context (k = 1–5). Inline rules or an external JSON
  profile via `quality.sequencing_errors.kmer_length` and `context_rules`.
- **Sequencing indels**. Sequencing-level insertion and deletion errors
  independent of somatic mutations, drawn from a geometric length distribution.
  Configurable via `quality.sequencing_errors.indel_rate`.
- **Strand bias model**. Per-read error rate asymmetry between R1 and R2 via
  `quality.sequencing_errors.r2_error_multiplier` and `r2_quality_offset`.
- **Correlated phasing burst errors**. Runs of correlated errors modelling
  phasing failures, controlled by `quality.sequencing_errors.burst_rate` and
  `burst_length_mean`.
- **ProfileLearner CIGAR-based indel extraction**. `varforge learn-profile`
  now counts CIGAR `I` and `D` operations per cycle (MAPQ ≥ 30 only) and
  exports `indel_error_profile` and `cycle_error_rates` fields to the profile
  JSON. Loading these fields auto-configures the `ErrorOrchestrator`.
- **Three platform presets**: `illumina_novaseq`, `pacbio_hifi`, and
  `nanopore_r10`. Each preset sets realistic error rates, indel rates, and
  cycle models for the target platform.

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
