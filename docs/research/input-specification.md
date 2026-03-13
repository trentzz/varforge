# Input Specification Design

## Primary Format: YAML Configuration

YAML is the recommended primary input format for VarForge. It is expressive enough to describe complex simulation scenarios (clonal architecture, UMI parameters, cfDNA models) while being human-readable and easily validated with serde.

### Full Example Configuration

```yaml
# VarForge simulation configuration
reference: /path/to/hg38.fa
output_prefix: sim_liquid_biopsy
seed: 42

# Sample-level parameters
sample:
  purity: 0.02              # 2% tumour fraction (liquid biopsy)
  ploidy: 2.0
  coverage: 2000            # 2000x for liquid biopsy panel
  read_length: 150
  paired_end: true

# Fragment size model
fragments:
  model: cfDNA              # "standard" | "cfDNA" | "custom"
  # Standard model parameters (used when model = "standard"):
  # mean: 350
  # std: 50
  #
  # cfDNA model uses nucleosomal mixture automatically
  # Custom model allows explicit mixture components:
  # components:
  #   - { mean: 167, std: 15, weight: 0.85 }
  #   - { mean: 143, std: 10, weight: 0.08 }
  #   - { mean: 334, std: 25, weight: 0.05 }
  #   - { mean: 500, std: 30, weight: 0.02 }

# Clonal architecture
clones:
  - id: founder
    ccf: 1.0                # Present in all tumour cells
    mutations:
      - { chr: chr7,  pos: 140453136, ref: A, alt: T }     # BRAF V600E
      - { chr: chr17, pos: 7578406,   ref: C, alt: T }     # TP53
  - id: subclone_A
    ccf: 0.3                # Present in 30% of tumour cells
    parent: founder
    mutations:
      - { chr: chr12, pos: 25398284, ref: C, alt: A }      # KRAS
  - id: subclone_B
    ccf: 0.15
    parent: founder
    mutations:
      - { chr: chr3, pos: 178936091, ref: G, alt: A }      # PIK3CA

# Alternative: load mutations from a VCF file
# mutations_vcf: /path/to/mutations.vcf

# Regions to simulate (optional; default = whole genome)
# regions: /path/to/panel.bed

# Sequencing error profile
errors:
  model: parametric         # "parametric" | "learned"
  base_error_rate: 0.001    # Overall error rate
  # For learned model:
  # profile_bam: /path/to/real_data.bam

# Library preparation
library:
  pcr_cycles: 8
  duplicate_rate: 0.15
  gc_bias: true
  gc_bias_model: standard   # "standard" | "none" | learned from BAM

# Artifacts (all rates are per-eligible-base)
artifacts:
  ffpe_damage_rate: 0.0     # C>T deamination rate (0 = disabled)
  oxog_rate: 0.0            # G>T oxidative damage rate
  deamination_rate: 0.0     # General deamination

# UMI configuration
umi:
  enabled: true
  length: 8                 # UMI barcode length in bp
  type: duplex              # "simplex" | "duplex"
  position: inline_r1       # "inline_r1" | "inline_both" | "index_read"
  error_rate: 0.01          # Per-base UMI sequencing error rate
  family_size:
    distribution: lognormal  # "lognormal" | "negative_binomial" | "fixed"
    mean: 5
    std: 3
    min: 1

# Output formats
output:
  fastq: true
  bam: true
  truth_vcf: true
  truth_bed: false           # For CNV regions
  manifest: true             # JSON summary of simulation parameters

# Read group information
read_group:
  id: "SIM001"
  sample: "SIMULATED_TUMOR"
  library: "SIM_LIB_001"
  platform: "ILLUMINA"
  platform_unit: "SIM_FLOWCELL.1"
```

### Minimal Configuration

```yaml
reference: /path/to/hg38.fa
output_prefix: simple_test

sample:
  coverage: 100
  read_length: 150

clones:
  - id: tumour
    ccf: 1.0
    mutations:
      - { chr: chr7, pos: 140453136, ref: A, alt: T }
```

All unspecified fields use sensible defaults (purity=1.0, ploidy=2.0, standard fragment model, parametric errors, no UMIs, no artifacts, FASTQ output only).

---

## Mutation Specification via VCF

For users who prefer standard formats or have mutations from existing pipelines, VarForge accepts VCF files with custom INFO fields:

```vcf
##fileformat=VCFv4.3
##INFO=<ID=CCF,Number=1,Type=Float,Description="Cancer cell fraction">
##INFO=<ID=CLONE,Number=1,Type=String,Description="Clone identifier">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr7	140453136	.	A	T	.	.	CCF=1.0;CLONE=founder
chr17	7578406	.	C	T	.	.	CCF=1.0;CLONE=founder
chr12	25398284	.	C	A	.	.	CCF=0.3;CLONE=subclone_A
chr2	29416089	.	A	ACCCGT	.	.	CCF=1.0;CLONE=founder
chr5	112175770	.	ACTG	A	.	.	CCF=0.5;CLONE=subclone_B
chr9	5000000	.	N	<DEL>	.	.	SVTYPE=DEL;SVLEN=-50000;CCF=1.0;CLONE=founder
```

VAF is NOT specified in the VCF -- it is computed by VarForge from the CCF, sample purity, and copy number at each locus. This is more biologically correct than specifying VAF directly.

Referenced from YAML config:
```yaml
mutations_vcf: /path/to/mutations.vcf
```

---

## Key User-Controllable Parameters

| Category | Parameter | Type | Default |
|----------|-----------|------|---------|
| **Coverage** | `sample.coverage` | int | 30 |
| | `regions` (BED file) | path | whole genome |
| **Reads** | `sample.read_length` | int | 150 |
| | `sample.paired_end` | bool | true |
| **Fragments** | `fragments.model` | enum | standard |
| | `fragments.mean` / `fragments.std` | float | 350 / 50 |
| **Tumour** | `sample.purity` | float | 1.0 |
| | `sample.ploidy` | float | 2.0 |
| | `clones[].ccf` | float | required |
| **Mutations** | Inline YAML or VCF file | - | required |
| **Errors** | `errors.base_error_rate` | float | 0.001 |
| | `errors.model` | enum | parametric |
| **Library** | `library.pcr_cycles` | int | 8 |
| | `library.duplicate_rate` | float | 0.0 |
| | `library.gc_bias` | bool | false |
| **Artifacts** | `artifacts.ffpe_damage_rate` | float | 0.0 |
| | `artifacts.oxog_rate` | float | 0.0 |
| **UMI** | `umi.enabled` | bool | false |
| | `umi.length` | int | 8 |
| | `umi.type` | enum | simplex |
| | `umi.error_rate` | float | 0.01 |
| **Output** | `output.fastq` / `output.bam` | bool | true / false |
| | `output.truth_vcf` | bool | true |
| **Reproducibility** | `seed` | int | random |

---

## CLI Interface

```
varforge simulate --config simulation.yaml
varforge simulate --config simulation.yaml --threads 8
varforge learn-profile --bam real_data.bam --output profile.json
varforge validate --truth truth.vcf --calls calls.vcf
```

### Subcommands

| Command | Description |
|---------|-------------|
| `simulate` | Run simulation from YAML config |
| `learn-profile` | Learn error/quality profile from a real BAM file |
| `validate` | Compare variant calls against truth VCF (convenience wrapper) |
