# VarForge Extension Roadmap

Generated: 2026-03-23. Based on external feedback from trentzz/kam and independent
codebase analysis.

---

## Summary

The immediate work is one confirmed panic (bug 1) and one behavioural gap (duplex
tracking). Beyond those, there are nine feature requests from the kam benchmark run.
This document also adds fourteen items from independent analysis aimed at broadening
the range of supported use cases.

---

## Part A: Feedback from trentzz/kam

### A1. [BUG] Out-of-bounds panic on short post-deletion haplotype

**Priority**: P0 — blocks simulation at 0.25% VAF benchmarks.

**Symptom**: `range end index 150 out of range for slice of length 107`

**Likely cause**: After a large deletion near a chromosome boundary, the remaining
haplotype region is shorter than `read_length`. The fragment start is sampled within
this truncated region, then a 150 bp read is attempted against ~107 bp of remaining
sequence. The issue is in `engine.rs` `generate_region`: the fragment length sampled
from the distribution exceeds the actual region length after the deletion shrinks it.
Even though `frag_end` is clamped to `region.end`, SV spike-in (applied post-hoc to
already-constructed reads) may change the reads' effective haplotype length. The
translocation path also fetches partner sequence without checking chromosome bounds.

**Fix**:
- In the pre-SV fragment loop: after clamping `actual_frag_len`, skip fragments where
  `actual_frag_len < read_length / 2` (too short to produce useful reads) rather than
  padding with N.
- In `apply_sv_to_read` translocation branch: clamp `partner_end` to chromosome length
  before calling `reference.sequence`.
- Add a regression test: simulate with a deletion that leaves fewer than `read_length`
  bases before the chromosome end.

---

### A2. [FEATURE] Multi-VAF batch mode

**Priority**: P1

**Config addition**:
```yaml
vafs: [0.25, 0.5, 1.0, 2.0, 5.0]   # generates one output dir per VAF
```

When `vafs` is set, `vaf` in the mutations block is ignored. Each VAF level runs as an
independent simulation with output in `{output.directory}/{vaf_pct}pct/`. A top-level
manifest lists all runs.

**Implementation**:
- Add `vafs: Option<Vec<f64>>` to `Config`.
- In the `simulate` command, detect `vafs`, build per-VAF `Config` copies (override
  `mutations.random.vaf_min` and `vaf_max` to pin the VAF, and the output directory),
  then run each sequentially or in parallel.
- Single VCF input mode: scale all variant VAFs by the ratio `target_vaf / original_vaf`
  for each batch step.

---

### A3. [FEATURE] Dry-run with expected alt molecule counts

**Priority**: P1

**Behaviour**: `varforge simulate --dry-run config.yaml` prints a table without
generating any reads:

```
Variant       VAF     Coverage  Family size (mean)  Expected alt families
chr7:55191822  0.005   1000      3.0                 ~1.7
chr7:55191822  0.010   1000      3.0                 ~3.3
```

Formula: `n_alt_expected = coverage × vaf / family_size_mean`

**Implementation**:
- Add `--dry-run` flag to the CLI `simulate` subcommand.
- Before any I/O, compute the table for all variants × VAF levels and print to stdout.
- Exit 0 if all VAFs yield ≥1 expected alt family, exit 1 otherwise (so CI can gate
  on this).

---

### A4. [FEATURE] Variant-tagged read output

**Priority**: P1

**Behaviour**: Each read pair that carries a variant alt allele is tagged.

- FASTQ header: append ` VT:Z:chr7:55191822:SNV` to the read name.
- Sidecar TSV: `{output_dir}/variant_reads.tsv` with columns
  `read_name`, `chrom`, `pos`, `vartype`, `vaf`, `clone_id`.

**Implementation**:
- `ReadPair` needs an optional `Vec<VariantTag>` field.
- When a variant is spiked into a fragment in `generate_region`, push a tag.
- FASTQ writer appends the tag string to the read name.
- A new `VariantReadsWriter` streams tags to the TSV alongside existing writers.

---

### A5. [BUG/FEATURE] Duplex molecule tracking

**Priority**: P2

**Symptom**: `n_duplex_alt = 0` for all SV simulations even when `duplex: true`.

**Root cause (suspected)**: The UMI family expansion in `engine.rs` (~line 512)
generates PCR copies but does not pair AB and BA strand reads. `generate_pcr_copies`
in `umi/families.rs` creates copies from a single original, not from two complementary
originals. The duplex UMI format (`AAAA-BBBB`) is written to the read name, but no BA
strand is ever generated, so duplex consensus tools see only one strand and count zero
duplex families.

**Fix**:
- When `umi.duplex = true`, generate a second "partner" read pair for each original by
  reverse-complementing both reads and swapping R1/R2, then tagging with the swapped
  barcode. Both the AB and BA copies go into the family.
- Add `N_DUPLEX_ALT` INFO field to truth VCF (see A9).
- Add a unit test that verifies BA-strand reads appear when `duplex: true`.

---

### A6. [FEATURE] Per-target coverage model from BED depth column

**Priority**: P2

**Config addition**:
```yaml
capture:
  per_target_depth_bed: "targets_with_depth.bed"  # 4th column = depth multiplier
```

**Behaviour**: When the per-target BED is provided, use the 4th column as the absolute
depth multiplier for that target instead of the log-normal sampling model.

**Implementation**:
- Extend `CaptureModel` to optionally store a `HashMap<(chrom, start, end), f64>`
  per-target depth.
- In `coverage_multiplier_at`, look up the target and return the stored multiplier if
  present.

---

### A7. [FEATURE] Chemistry presets

**Priority**: P2

A named preset fills in UMI, fragment, and quality defaults in one field:

```yaml
preset: twist-umi-duplex
```

Built-in presets:

| Name | umi.length | umi.duplex | fragment.mean | fragment.sd | quality.mean |
|------|-----------|------------|---------------|-------------|--------------|
| twist-umi-duplex | 8 | true | 200 | 30 | 37 |
| illumina-wgs | — | false | 300 | 50 | 36 |
| illumina-wes | — | false | 200 | 40 | 35 |
| illumina-ctdna | 8 | false | 167 | 20 | 36 |
| pacbio-hifi | — | false | 15000 | 5000 | — |
| nanopore-r10 | — | false | 20000 | 10000 | — |

**Implementation**:
- Add `preset: Option<String>` to `Config`.
- After YAML deserialisation, apply preset defaults for any fields that are `None`.
- Presets are a static lookup table in `io/config.rs`; no external file needed.

---

### A8. [FEATURE] Empirical error profile from real BAM

**Priority**: P3

**Config addition**:
```yaml
quality:
  profile_bam: "real_run.bam"   # learn from this file
```

`varforge learn-profile` is already stubbed. This completes it:
- Read quality scores by cycle position from the BAM.
- Store a 2D array `[cycle][quality_bin]` of counts.
- Save as a `.profile.bin` (bincode-serialised) or `.profile.tsv`.
- `quality.profile_path` then loads the saved profile.

The BAM learning path just reads query quality strings; it does not need alignment.

---

### A9. [FEATURE] Simulated molecule counts in truth VCF INFO

**Priority**: P3

Add two INFO fields to the truth VCF:

```
N_ALT_MOL=<int>      # number of read pairs carrying the alt allele
N_DUPLEX_ALT=<int>   # number of AB+BA duplex families carrying the alt allele
```

**Implementation**:
- `AppliedVariant` already tracks `actual_alt_count` and `actual_total_count`.
- Aggregate these across regions in the simulate command.
- Pass the per-variant molecule counts to `TruthVcfWriter`.

---

## Part B: Independent Analysis

### B1. [FEATURE] Long-read sequencing support (PacBio HiFi, Nanopore)

**Priority**: P1 — unlocks a completely different set of target callers.

Long reads require:
- Fragment lengths of 5–50 kbp (sampled from log-normal or gamma).
- PacBio HiFi: high accuracy (Q20–Q30), low indel rate (~0.1%), near-zero systematic
  error. CCS consensus model.
- Nanopore R10: Q15–Q25, homopolymer-dependent error, strand bias (C>T on one strand).
- Single-read BAM output (no read2, template length = read length).
- No UMI (long reads use different deduplication).

**Implementation**:
- Add `platform: "pacbio-hifi" | "nanopore-r10"` to `SampleConfig`.
- Add `LongReadFragmentSampler` and `LongReadQualityModel`.
- BAM writer: detect single-read mode and omit R2 fields.
- The chemistry presets in A7 cover PacBio HiFi and Nanopore R10 as defaults.

---

### B2. [FEATURE] COSMIC mutational signatures

**Priority**: P1 — makes simulated somatic mutations biologically realistic.

**Config addition**:
```yaml
mutations:
  signatures:
    - name: SBS4    # tobacco smoking
      weight: 0.6
    - name: SBS2    # APOBEC
      weight: 0.4
```

VarForge uses the trinucleotide context probabilities from COSMIC v3.4 to sample the
alt base and position for each random SNV. The existing random mutation generator
currently picks positions uniformly; this replaces the base-selection step.

**Data**: Embed the 96-channel SBS, 78-channel DBS, and 83-channel ID signature
matrices as a compressed static asset (COSMIC licence allows this for research tools).

---

### B3. [FEATURE] Germline variant simulation

**Priority**: P1

Many somatic callers are tested against germline SNPs to evaluate false positive rates.
Currently VarForge has no germline variant model.

**Config addition**:
```yaml
germline:
  population_vcf: "gnomad_chr7.vcf.gz"   # or use built-in model
  het_snp_density: 1.0                    # per-kbp, used when no VCF
  hom_snp_density: 0.3
```

Germline variants are applied at VAF 0.5 (het) or 1.0 (hom), to both tumour and normal
samples simultaneously. They appear in a separate truth VCF (`germline_truth.vcf`).

---

### B4. [FEATURE] Paired tumour-normal simultaneous generation

**Priority**: P1

Running tumour and normal separately then merging is error-prone. Add:

```yaml
mode: paired-tumour-normal
tumour:
  purity: 0.7
  output_subdir: tumour
normal:
  output_subdir: normal
```

Both samples share the reference genome, germline variants, and read simulation
parameters (coverage may differ). Somatic mutations appear only in the tumour output.

**Implementation**:
- New `PairedSimulationMode` in `core/multi_sample.rs`.
- Reuse region iteration and reference loading; fork the per-read pipeline at the
  mutation-spike stage.

---

### B5. [FEATURE] Chromosome-level haplotype phasing

**Priority**: P2

Currently, mutations on heterozygous positions are assigned randomly to whichever read
covers the site. Phasing-aware callers (WhatsHap, SHAPEIT) and allele-specific
expression tools need reads that are correctly phased.

**Config addition**:
```yaml
mutations:
  phasing: random   # current behaviour
  # phasing: haplotype   # all mutations assigned to hap A or hap B deterministically
```

**Implementation**:
- Assign each variant a `haplotype: Option<u8>` field (0 = hap A, 1 = hap B).
- When spiking, only apply the variant if the fragment start was sampled from the
  correct haplotype (determined by RNG with fixed seed per-region).

---

### B6. [FEATURE] COSMIC SV signatures and complex SVs

**Priority**: P2

Single SVs are already supported. Add:

- **Tandem duplicator phenotype**: Enrichment of short tandem duplications (1–10 kbp).
- **HRD pattern**: Large deletions (100 kbp–10 Mbp) with LOH, simulating BRCA1/2
  deficiency.
- **Chromothripsis**: A configurable number of random cuts on one chromosome,
  with random re-ligation (subset of fragments inverted or deleted).

**Config addition**:
```yaml
mutations:
  sv_signature: HRD   # or TDP (tandem duplicator phenotype)
```

---

### B7. [FEATURE] Cross-sample contamination

**Priority**: P2

Useful for testing CONPAIR, VerifyBamID, Concordance callers.

```yaml
contamination:
  source_sample: "donor_B.bam"   # reads randomly drawn from this BAM
  fraction: 0.01                  # 1% contamination
```

Or synthetic contamination (no source BAM):
```yaml
contamination:
  fraction: 0.01
  snp_profile: "donor_b_snps.vcf"
```

---

### B8. [FEATURE] Fragment end motif simulation (cfDNA)

**Priority**: P2

cfDNA liquid biopsy tools (DELFI, Griffin, CancerSEEK) rely on fragment end motifs
(the 4-mer at each fragment end). Real cfDNA has a strong preference for certain end
contexts driven by DNase and MNase cleavage.

**Config addition**:
```yaml
fragment:
  end_motif_model: "plasma"   # uses empirical 4-mer frequencies from literature
```

**Implementation**:
- Store a 256-entry frequency table (4-mer counts for each possible end context).
- During fragment sampling, use rejection sampling to match the target end-motif
  distribution.

---

### B9. [FEATURE] Simulation QC report

**Priority**: P2

After each run, write `{output_dir}/sim_report.json` with:
- Achieved coverage per chromosome (mean, sd).
- Per-variant: target VAF, actual alt count, actual total count, actual VAF.
- Family size distribution (mean, median, 95th percentile).
- Duplicate rate achieved.
- GC bias curve (GC fraction vs normalised coverage).

This replaces the manual post-hoc analysis step that users currently have to do
themselves to verify the simulation ran as expected.

---

### B10. [FEATURE] Benchmarking mode

**Priority**: P3

```
varforge benchmark-suite \
  --config base.yaml \
  --vafs 0.1,0.5,1.0,2.0,5.0 \
  --coverages 100,500,1000 \
  --output benchmark_suite/
```

Generates a grid of configs and a manifest linking each output dir to its parameters.
Designed to feed directly into benchmarking pipelines like kam.

---

### B11. [FEATURE] Amplicon sequencing support

**Priority**: P3

Targeted amplicon panels (Ion Torrent, TruSight Oncology) have:
- Very tight fragment distribution (fragments are exactly amplicon length).
- Read pairs that always start at primer coordinates.
- Primer sequence contamination at read ends.
- No meaningful off-target reads.

```yaml
capture:
  mode: amplicon    # vs default "hybrid-capture"
  amplicons_bed: "amplicons.bed"
  primer_trim: true
```

---

### B12. [FEATURE] Tumour-in-normal contamination

**Priority**: P3

Simulates a normal sample with low-level somatic contamination (e.g., circulating
tumour cells in peripheral blood). Relevant for testing contamination-aware germline
callers and matched-normal validation.

```yaml
normal:
  tumour_contamination: 0.002   # 0.2% tumour DNA in normal sample
```

---

### B13. [FEATURE] Config templating

**Priority**: P3

Support `${variable}` substitution in YAML values:

```yaml
sample:
  name: "${sample_name}"
  coverage: ${coverage}
```

```
varforge simulate \
  --config template.yaml \
  --set sample_name=PATIENT_01 \
  --set coverage=500
```

Avoids maintaining N near-identical config files for N samples or N coverage levels.

---

### B14. [FEATURE] MSI simulation

**Priority**: P3

Microsatellite instability tumours have elevated indel rates at microsatellite loci.

**Config addition**:
```yaml
tumour:
  msi: true          # enables microsatellite-targeted indel enrichment
  msi_bed: null      # optional: restrict to specific loci
```

Requires a microsatellite loci annotation BED, either user-supplied or bundled for
common reference genomes (hg19, hg38, GRCh38).

---

## Implementation order

### Immediate (next sprint)
1. A1 — panic fix (SV deletion boundary)
2. A5 — duplex tracking fix
3. A2 — multi-VAF batch mode
4. A3 — dry-run expected counts
5. A9 — molecule counts in truth VCF

### Short term
6. A4 — variant-tagged reads
7. A7 — chemistry presets (prerequisite for B1)
8. B1 — long-read support
9. B3 — germline variants
10. B4 — paired tumour-normal mode

### Medium term
11. B2 — COSMIC signatures
12. A6 — per-target depth BED
13. B9 — QC report
14. A8 — empirical error profile from BAM
15. B5 — haplotype phasing

### Longer term
16. B6 — complex SVs / SV signatures
17. B7 — cross-sample contamination
18. B8 — cfDNA end motif model
19. B10 — benchmarking mode
20. B11 — amplicon support
21. B12 — tumour-in-normal contamination
22. B13 — config templating
23. B14 — MSI simulation
