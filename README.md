# VarForge

[![CI](https://github.com/your-org/varforge/actions/workflows/ci.yml/badge.svg)](https://github.com/your-org/varforge/actions)
[![Crates.io](https://img.shields.io/crates/v/varforge.svg)](https://crates.io/crates/varforge)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

VarForge is a fast, single-binary Rust tool for generating synthetic cancer sequencing test data with controlled ground truth. It produces realistic FASTQ/BAM files with known mutations, tumour parameters, UMI barcodes, and cfDNA fragment profiles for benchmarking bioinformatics pipelines.

---

## Features

| Feature | VarForge | BAMSurgeon | ART | NEAT |
|---------|----------|------------|-----|------|
| Single binary, no runtime deps | Yes | No | No | No |
| Controlled somatic mutations (SNV/indel/MNV) | Yes | Yes | No | Partial |
| Tumour purity / clonal architecture | Yes | No | No | No |
| cfDNA fragment model | Yes | No | No | No |
| Duplex UMI barcodes | Yes | No | No | No |
| FFPE / oxoG artefacts | Yes | No | No | No |
| Longitudinal / multi-sample series | Yes | No | No | No |
| Copy number alterations | Yes | No | No | No |
| GC bias model | Yes | No | Partial | No |
| Truth VCF output | Yes | Partial | No | Yes |
| YAML configuration | Yes | No | No | No |

---

## Installation

### From crates.io

```
cargo install varforge
```

### From source

```
git clone https://github.com/your-org/varforge
cd varforge
cargo build --release
./target/release/varforge --help
```

### Pre-built binaries

Pre-built binaries for Linux (x86_64, aarch64) and macOS (Apple Silicon) are available on the [Releases](https://github.com/your-org/varforge/releases) page.

---

## Quickstart

The following example generates a 30x WGS tumour sample with 5000 random somatic mutations in about 20–40 minutes on an 8-core machine.

**1. Create a config file (`quickstart.yaml`):**

```yaml
reference: /data/ref/hg38.fa

output:
  directory: out/quickstart
  fastq: true
  truth_vcf: true

sample:
  name: TUMOUR
  coverage: 30.0

tumour:
  purity: 0.70

mutations:
  random:
    count: 5000
    vaf_min: 0.05
    vaf_max: 0.60
    snv_fraction: 0.80
    indel_fraction: 0.15
    mnv_fraction: 0.05

seed: 42
```

**2. Validate the config:**

```
varforge validate --config quickstart.yaml
```

**3. Run the simulation:**

```
varforge simulate --config quickstart.yaml
```

**4. Inspect the output:**

```
out/quickstart/
  TUMOUR_R1.fastq.gz       # Read 1 FASTQ
  TUMOUR_R2.fastq.gz       # Read 2 FASTQ
  truth.vcf.gz             # Ground-truth VCF with all injected variants
  manifest.tsv             # Sample metadata (name, coverage, purity, paths)
```

**CLI overrides** — any config value can be overridden at the command line:

```
varforge simulate --config quickstart.yaml \
    --coverage 60 --purity 0.5 --seed 99
```

**Presets** — skip the config entirely for common scenarios:

```
varforge simulate --config quickstart.yaml --preset wgs
varforge simulate --config quickstart.yaml --preset cancer:melanoma
```

---

## CLI Reference

```
varforge [OPTIONS] <COMMAND>

Commands:
  simulate      Run a simulation from a YAML config
  validate      Validate a YAML config without running
  edit          Spike variants into an existing BAM file
  learn-profile Learn an error/quality profile from a real BAM file

Global options:
  -t, --threads <N>          Number of threads (default: all available cores)
      --log-level <LEVEL>    error | warn | info | debug | trace (default: info)
```

### `simulate`

```
varforge simulate --config <FILE> [OPTIONS]

Options:
  -c, --config <FILE>           Path to YAML configuration file (required)
  -o, --output-dir <DIR>        Override output directory
      --seed <N>                Override random seed
      --coverage <F>            Override coverage depth (x)
      --read-length <N>         Override read length (bp)
      --purity <F>              Override tumour purity (0.0–1.0)
      --fragment-mean <F>       Override fragment mean length (bp)
      --fragment-sd <F>         Override fragment length standard deviation (bp)
      --random-mutations <N>    Generate N random mutations (no VCF needed)
      --vaf-range <MIN-MAX>     VAF range for random mutations (e.g. 0.001-0.05)
      --preset <NAME>           Apply a named preset (see Presets section)
      --dry-run                 Validate config and estimate output size only
```

### `validate`

```
varforge validate --config <FILE>
```

Parses the YAML config and checks all fields for consistency. Exits with status 0 if valid, non-zero otherwise with a descriptive error message.

### `edit`

```
varforge edit --bam <IN.bam> --vcf <VARIANTS.vcf> --output <OUT.bam>
```

Spikes variants from a VCF directly into an existing BAM file without re-simulating reads. Useful for adding a handful of known mutations to a real or previously simulated dataset.

### `learn-profile`

```
varforge learn-profile --bam <BAM> --output <PROFILE.json>
```

Learns an empirical base-quality and error profile from a real BAM file. The resulting JSON can be referenced from the `quality.profile_path` config field to produce reads with a realistic, data-driven quality model instead of the parametric default.

---

## Configuration Reference

All simulation parameters are specified in a YAML file. Only `reference` and `output.directory` are required; everything else has a sensible default.

### Top-level fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `reference` | path | — | Path to FASTA reference genome (required) |
| `output` | OutputConfig | — | Output format and directory (required) |
| `sample` | SampleConfig | see below | Read generation parameters |
| `fragment` | FragmentConfig | see below | Insert size distribution |
| `quality` | QualityConfig | see below | Base quality model |
| `tumour` | TumourConfig | null | Tumour purity and clonal architecture |
| `mutations` | MutationConfig | null | Somatic mutation injection |
| `umi` | UmiConfig | null | UMI barcode configuration |
| `artifacts` | ArtifactConfig | null | Sequencing artefact simulation |
| `copy_number` | list of CopyNumberConfig | null | Copy number alterations |
| `gc_bias` | GcBiasConfig | null | GC content coverage bias |
| `capture` | CaptureConfig | null | Hybrid-capture enrichment model |
| `samples` | list of SampleEntry | null | Multi-sample / longitudinal series |
| `chromosomes` | list of strings | null (all) | Restrict simulation to named chromosomes |
| `regions_bed` | path | null | Restrict simulation to BED file regions |
| `seed` | integer | null (random) | Random seed for reproducibility |
| `threads` | integer | null (all cores) | Worker thread count |

---

### `output`

```yaml
output:
  directory: out/my_run   # required
  fastq: true             # write gzip-compressed FASTQ files (default: true)
  bam: false              # write coordinate-sorted BAM (default: false)
  truth_vcf: true         # write ground-truth VCF (default: true)
  manifest: true          # write manifest.tsv (default: true)
```

---

### `sample`

```yaml
sample:
  name: TUMOUR_01       # sample name used in file names and read headers (default: SAMPLE)
  read_length: 150      # read length in bp (default: 150)
  coverage: 30.0        # mean target coverage depth (default: 30.0)
  platform: illumina    # sequencing platform tag; written to BAM @RG header (optional)
```

---

### `fragment`

Controls the insert size (fragment length) distribution.

```yaml
fragment:
  model: normal    # normal | cfda | custom (default: normal)
  mean: 300.0      # mean fragment length in bp (default: 300.0)
  sd: 50.0         # standard deviation in bp (default: 50.0)
```

**Fragment models:**

- `normal` — Gaussian distribution. Suitable for standard library prep from fresh-frozen tissue or cell lines.
- `cfda` — Short, nucleosome-phased distribution reflecting cell-free DNA in plasma. Typical mean ~167 bp, sd ~20 bp.
- `custom` — Reserved for future user-supplied length histograms.

---

### `quality`

```yaml
quality:
  mean_quality: 36        # mean Phred quality score for the first cycle (default: 36)
  tail_decay: 0.003       # per-cycle quality decay rate (default: 0.003)
  profile_path: null      # optional path to empirical profile JSON from learn-profile
```

If `profile_path` is set, the empirical profile overrides `mean_quality` and `tail_decay`.

---

### `tumour`

```yaml
tumour:
  purity: 0.70    # fraction of cells that are tumour (0.0–1.0; default: 1.0)
  ploidy: 2       # tumour ploidy (default: 2)
  clones:         # optional list of clones for subclonal architecture
    - id: trunk
      ccf: 1.0           # cancer cell fraction (0.0–1.0)
    - id: subclone_a
      ccf: 0.40
      parent: trunk      # parent clone ID (optional; omit for founding clone)
```

When `clones` is empty, all mutations are assigned to a single clonal population at the specified `purity`. When clones are defined, mutations are distributed across the clone tree and their effective VAF is:

```
VAF = purity × CCF / ploidy
```

---

### `mutations`

```yaml
mutations:
  vcf: /path/to/variants.vcf.gz   # optional: inject specific variants from VCF
  random:                          # optional: add random somatic mutations
    count: 5000                    # number of mutations to generate
    vaf_min: 0.001                 # minimum VAF (default: 0.001)
    vaf_max: 0.50                  # maximum VAF (default: 0.5)
    snv_fraction: 0.80             # fraction that are SNVs (default: 0.80)
    indel_fraction: 0.15           # fraction that are indels (default: 0.15)
    mnv_fraction: 0.05             # fraction that are MNVs (default: 0.05)
```

`snv_fraction + indel_fraction + mnv_fraction` must sum to exactly 1.0.

Both `vcf` and `random` may be specified simultaneously; the VCF variants are injected first and then random mutations are added at non-overlapping positions.

---

### `umi`

```yaml
umi:
  length: 8              # UMI barcode length in bases (default: 8)
  duplex: false          # enable duplex (double-stranded) UMI mode (default: false)
  pcr_cycles: 10         # number of PCR amplification cycles (default: 10)
  family_size_mean: 3.0  # mean read family size (default: 3.0)
  family_size_sd: 1.5    # standard deviation of family size (default: 1.5)
  inline: true           # prepend UMI to read sequence (default: false)
```

When `inline: true`, the UMI is prepended to the read sequence (e.g. for fgbio `ExtractUmisFromBam`). When `inline: false`, the UMI is written into the read name (e.g. `@READ:ACGTACGT`).

When `duplex: true`, each molecule is tagged with a strand-specific UMI pair supporting duplex consensus calling tools such as fgbio `CallDuplexConsensusReads`.

---

### `artifacts`

```yaml
artifacts:
  ffpe_damage_rate: 0.02    # C>T deamination rate (0.0–1.0; null = disabled)
  oxog_rate: 0.01           # 8-oxoG C>A transversion rate (0.0–1.0; null = disabled)
  duplicate_rate: 0.15      # PCR duplicate fraction (0.0–1.0; null = disabled)
  pcr_error_rate: 0.001     # PCR substitution error rate per base (null = disabled)
```

All fields are optional. Omit the entire `artifacts` block (or set individual rates to `null`) to disable artefact simulation.

---

### `copy_number`

```yaml
copy_number:
  - region: "chr7:55000000-55200000"   # chrom:start-end (1-based, inclusive)
    tumor_cn: 4                        # tumour copy number (default: 2)
    normal_cn: 2                       # normal copy number (default: 2)
    major_cn: 3                        # major allele CN for LOH modeling (optional)
    minor_cn: 1                        # minor allele CN for LOH modeling (optional)
```

Multiple entries may be listed. Overlapping regions are applied in order (last wins). Read depth in each region is scaled proportionally to `tumor_cn / normal_cn`.

---

### `gc_bias`

```yaml
gc_bias:
  enabled: true      # apply GC bias model (default: true when block is present)
  model: default     # default | flat | custom (default: "default")
  severity: 1.0      # bias multiplier: 0 = none, 1 = realistic, 2 = extreme (default: 1.0)
```

The `default` model applies an empirical coverage reduction at GC extremes (< 30 % or > 70 % GC). Setting `severity: 0` disables the effect while keeping the block present.

---

### `capture`

```yaml
capture:
  enabled: true                             # activate capture model (default: true)
  targets_bed: /data/panels/panel.bed       # path to capture target BED (optional)
  off_target_fraction: 0.20                 # fraction of reads mapping off-target (default: 0.2)
  coverage_uniformity: 0.30                 # per-target LogNormal σ (0 = uniform; default: 0.3)
  edge_dropoff_bases: 50                    # exponential dropoff at target edges in bp (default: 50)
```

When `targets_bed` is omitted, the capture model distributes reads uniformly across whichever chromosomes or regions are active.

---

### `samples` (multi-sample / longitudinal)

When `samples` is present, VarForge generates one output sub-directory per entry and a combined `manifest.tsv`. Each entry shares the top-level `reference`, `mutations`, `tumour`, and `fragment` settings but can override coverage, tumour fraction, and fragment model independently.

```yaml
samples:
  - name: timepoint_1
    coverage: 1000.0
    tumour_fraction: 0.05        # ctDNA fraction for this sample (default: 1.0)
    fragment_model: cfda         # override fragment model (optional)
    clonal_shift:                # per-clone CCF adjustments at this timepoint (optional)
      subclone_a: 0.10
  - name: timepoint_2
    coverage: 1000.0
    tumour_fraction: 0.002
    fragment_model: cfda
```

---

## Presets Reference

Presets are named configuration bundles that set sensible defaults for common scenarios. A preset is applied *before* the YAML config values, so explicit YAML fields always win.

```
varforge simulate --config base.yaml --preset <NAME>
```

### Built-in presets

| Preset | Coverage | Fragment | Mutations | UMI | Notes |
|--------|----------|----------|-----------|-----|-------|
| `small` | 1x | normal | 100 random (chr22 only) | no | Smoke test; completes in ~30 s |
| `panel` | 500x | normal | 50 random | inline 8-mer | Targeted panel benchmarking |
| `wgs` | 30x | normal | 5 000 random | no | Whole-genome variant calling |
| `cfdna` | 200x | cfda (167 bp) | 200 random, VAF 0.1–5% | duplex | Liquid biopsy simulation |
| `ffpe` | 30x | normal | 500 random | no | FFPE artefacts enabled |
| `umi` | 1 000x | normal | 50 random | duplex 9-mer | High-depth duplex consensus |

### Cancer-type presets

Cancer presets are accessed with the `cancer:` namespace prefix. Each preset sets biologically realistic mutation counts, VAF ranges, tumour purity, and mutation-type fractions based on published COSMIC mutational signatures.

```
varforge simulate --config base.yaml --preset cancer:melanoma
```

| Preset | Cancer | Dominant Signature | Typical TMB | Purity |
|--------|--------|--------------------|-------------|--------|
| `cancer:lung_adeno` | Lung adenocarcinoma | SBS4 (smoking, C>A) | ~8 mut/Mb | 60% |
| `cancer:colorectal` | Colorectal (MSS) | SBS1/SBS5 (aging, C>T) | ~5 mut/Mb | 65% |
| `cancer:breast_tnbc` | Triple-negative breast | SBS3 (HRD, flat) | ~5 mut/Mb | 55% |
| `cancer:melanoma` | Cutaneous melanoma | SBS7a/b (UV, C>T) | ~30 mut/Mb | 70% |
| `cancer:aml` | Acute myeloid leukaemia | SBS1/SBS5 (aging) | ~1 mut/Mb | 80% |
| `cancer:prostate` | Prostate adenocarcinoma | SBS1/SBS5 (aging) | ~2 mut/Mb | 50% |
| `cancer:pancreatic` | Pancreatic ductal | SBS1/SBS5 (aging) | ~3 mut/Mb | 25% |
| `cancer:glioblastoma` | Glioblastoma (IDH-wt) | SBS1/SBS5 (aging) | ~4 mut/Mb | 65% |

---

## Example Configs

Ready-to-use example configs are in the `examples/` directory.

| File | Use case |
|------|----------|
| [`examples/minimal.yaml`](examples/minimal.yaml) | Simplest possible simulation (defaults only) |
| [`examples/wgs_30x.yaml`](examples/wgs_30x.yaml) | Standard 30x WGS tumour with random mutations |
| [`examples/panel_umi.yaml`](examples/panel_umi.yaml) | Targeted panel with inline 8-mer UMI |
| [`examples/cfdna_monitoring.yaml`](examples/cfdna_monitoring.yaml) | cfDNA longitudinal series (4 timepoints) |
| [`examples/ffpe_artifacts.yaml`](examples/ffpe_artifacts.yaml) | FFPE-damaged tumour sample |
| [`examples/tumor_normal.yaml`](examples/tumor_normal.yaml) | Matched tumour/normal pair |
| [`examples/subclonal.yaml`](examples/subclonal.yaml) | Four-clone tumour with copy number alterations |
| [`examples/high_depth.yaml`](examples/high_depth.yaml) | 1000x duplex UMI for low-VAF detection |
| [`examples/custom_mutations.yaml`](examples/custom_mutations.yaml) | Inject specific variants from a VCF file |

---

## Output Formats

### FASTQ headers

Read names follow the format:

```
@{SAMPLE}:{CHROM}:{POS}:{READ_NUM}[:UMI={BARCODE}]
```

Example: `@TUMOUR:chr7:55191822:1:UMI=ACGTACGT`

The UMI suffix is only present when a `umi` block is configured and `inline: false`.

### BAM tags

When BAM output is enabled (`output.bam: true`), the following non-standard tags are written:

| Tag | Type | Description |
|-----|------|-------------|
| `RX` | Z | Raw UMI sequence |
| `MI` | Z | Molecule ID (read family identifier) |
| `tp` | i | 1 if the read carries a simulated somatic variant, 0 otherwise |
| `cl` | Z | Clone ID the variant was assigned to |

### Truth VCF fields

The truth VCF written to `truth.vcf.gz` uses the following INFO fields:

| Field | Description |
|-------|-------------|
| `VAF` | Target variant allele frequency |
| `CLONE` | Clone ID the variant was assigned to |
| `CCF` | Cancer cell fraction of the assigned clone |
| `TYPE` | Variant type: `SNV`, `INDEL`, or `MNV` |

---

## Use Case Recipes

### Benchmarking a somatic variant caller

```yaml
# Generate matched tumour/normal with known variants.
reference: /data/ref/hg38.fa
output:
  directory: out/caller_bench
  bam: true
  truth_vcf: true
samples:
  - name: TUMOUR
    coverage: 60.0
    tumour_fraction: 1.0
  - name: NORMAL
    coverage: 30.0
    tumour_fraction: 0.0
tumour:
  purity: 0.65
mutations:
  random:
    count: 1000
    vaf_min: 0.05
    vaf_max: 0.60
    snv_fraction: 0.80
    indel_fraction: 0.15
    mnv_fraction: 0.05
seed: 42
```

Run your caller against `TUMOUR/` with `NORMAL/` as the matched normal. Evaluate with:

```
bcftools stats --apply-filters PASS \
    caller_output.vcf.gz truth.vcf.gz
```

### Benchmarking UMI deduplication (fgbio)

Use `examples/panel_umi.yaml` with `inline: true` and then pipe through fgbio:

```
fgbio ExtractUmisFromBam \
    --input out/panel_umi/PANEL_UMI.bam \
    --output extracted.bam \
    --read-structure 8M+T 8M+T

fgbio GroupReadsByUmi \
    --input extracted.bam \
    --output grouped.bam \
    --strategy paired

fgbio CallMolecularConsensusReads \
    --input grouped.bam \
    --output consensus.bam \
    --min-reads 1
```

### Liquid biopsy sensitivity curve

Generate cfDNA samples at multiple tumour fractions and measure detection rate at each:

```
for TF in 0.10 0.05 0.01 0.005 0.001; do
    varforge simulate --config examples/cfdna_monitoring.yaml \
        --output-dir out/tf_${TF} \
        --purity ${TF} \
        --seed 42
done
```

### FFPE artefact filter development

Use `examples/ffpe_artifacts.yaml` to generate data with realistic FFPE damage, then evaluate your artefact filter:

- True positives: variants present in `truth.vcf.gz`
- Artefacts: C>T calls not in the truth VCF with strand bias (use the `tp` BAM tag to distinguish)

---

## Performance Tuning

VarForge uses Rayon for data parallelism. Performance scales approximately linearly with thread count up to the number of chromosomes being simulated.

### Thread count

```
varforge simulate --config cfg.yaml --threads 16
```

Or set in the config:

```yaml
threads: 16
```

### Restricting scope

For development and testing, restrict to one or a few chromosomes:

```yaml
chromosomes:
  - chr22
```

Or to a BED file of target regions:

```yaml
regions_bed: /data/panels/hotspot_panel.bed
```

### Approximate runtimes (8-core laptop, hg38)

| Scenario | Coverage | Mutations | Time |
|----------|----------|-----------|------|
| `small` preset | 1x, chr22 | 100 | ~30 s |
| Panel (chr7, chr12, chr17) | 500x | 50 | ~2 min |
| WGS 30x | 30x, all | 5 000 | ~25 min |
| WGS 60x | 60x, all | 5 000 | ~50 min |
| Ultra-deep panel 1000x | 1000x, 3 chroms | 50 | ~8 min |

### Memory usage

Peak memory is approximately:

```
2 × read_length × threads × (coverage / 30) MB
```

For 30x WGS with 150 bp reads on 8 threads: ~600 MB. For 1000x panel on 8 threads: ~200 MB (limited region).

---

## Comparison with Other Tools

| Scenario | Recommended tool |
|----------|-----------------|
| Realistic Illumina base-quality profiles | VarForge (parametric or `learn-profile`) |
| Controlled somatic variant spike-in | VarForge or BAMSurgeon |
| Spike into a real patient BAM | BAMSurgeon (preserves real read background) |
| Simple read generation, no mutations | ART or NanoSim |
| Whole-genome de novo simulation | NEAT or VarForge |
| cfDNA / liquid biopsy data | VarForge (only tool with native cfDNA model) |
| UMI-tagged duplex sequencing | VarForge (only tool with native duplex model) |
| FFPE artefact simulation | VarForge (only tool with FFPE + oxoG model) |
| Multi-sample longitudinal series | VarForge (only tool with native time-series) |

VarForge is the right choice when you need a complete, reproducible, ground-truth dataset for pipeline benchmarking and do not have access to real patient sequencing data. BAMSurgeon is the better choice when you need to insert a small number of variants into an existing real-data BAM while preserving the authentic read background.

---

## License

[MIT](LICENSE)
