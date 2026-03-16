# Comprehensive Tool Comparison (March 2025)

## Overview

This document catalogs every significant synthetic sequencing data generation tool, organized by category, with detailed analysis of capabilities, limitations, and maintenance status. This informs VarForge's design to fill the gaps no existing tool addresses.

---

## Mutation Spike-in Tools (BAM modification)

### BAMSurgeon
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/adamewing/bamsurgeon |
| **Language** | Python (85%), Shell (14%) |
| **Stars/Citations** | 247 stars; highly cited (DREAM Challenge papers in Nature Methods) |
| **Last Release** | v1.4.1 (August 2022); 622 commits |
| **Status** | Maintenance mode / low activity |

**Capabilities**: Adds SNVs, indels, SVs to existing BAMs via local contig assembly and read replacement. Gold standard for ICGC-TCGA DREAM Somatic Mutation Calling Challenge. Multiple aligner support (BWA, Novoalign, GSNAP, Bowtie2).

**Limitations**:
- **Very slow**: ~4-5 sec/variant for >20K SNVs; optimized to ~1 sec/variant but still prohibitive for WGS
- Requires >30x coverage BAM as input for local assembly
- **Deterministic spike-in bias**: 10% VAF at 30x = exactly 3 alt reads (no stochastic sampling)
- Actual breakpoints may differ from input due to contig assembly imprecision
- Cannot simulate from scratch -- needs existing BAM
- 39 open issues (header errors with samtools 1.10, WES region failures, quality score mismatches)
- No UMI/duplex, cfDNA, or artifact simulation

### SomatoSim
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/BieseckerLab/SomatoSim |
| **Language** | Python |
| **Stars** | 7 stars |
| **Last Activity** | June 2022 |
| **Status** | Inactive |

**Capabilities**: Precise control over VAF and variant positions. Preserves real sequencing environment. Models strand distribution. MQ/BQ thresholds. Extensive simulation reports.

**Limitations**: SNVs only (no indels, SVs, CNVs). Exome/targeted only. Requires existing BAM. Same deterministic spike-in as BAMSurgeon.

### stochasticSim
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/BrianOSullivanGit/stochasticSim |
| **Language** | C (HTSlib) |
| **Paper** | NAR Cancer 2023 |
| **Status** | Published 2023 |

**Capabilities**: Haplotype-aware, phased tumor simulations from 1000 Genomes. **Addresses stochastic sampling** (unlike BAMSurgeon). Includes FFPE deamination and 8-oxoG artifacts. Tracks source of all non-reference bases.

**Limitations**: No somatic indels, no SVs/CNVs. Exome-scale only. Small user base. No UMI or cfDNA.

### Xome-Blender
| Attribute | Detail |
|-----------|--------|
| **Paper** | PLOS ONE 2018 |

**Capabilities**: Mixes real tumour/normal BAMs at specified ratios for synthetic cancer genomes with user-defined subclone count and CNA count.

**Limitations**: Requires existing real tumour/normal BAMs. Cannot generate de novo.

---

## Read Simulators (FASTQ generation)

### ART
| Attribute | Detail |
|-----------|--------|
| **URL** | https://www.niehs.nih.gov/research/resources/software/biostatistics/art |
| **Language** | C++ |
| **Citations** | >3500 (one of the most cited read simulators) |
| **Last Release** | "GreatSmokyMountains" (April 2016) |
| **Status** | Unmaintained (original); art_modern actively maintained |

**Capabilities**: Empirical, technology-specific error models (position-dependent quality profiles). Supports Illumina, 454, SOLiD. Generates golden ground-truth aligned reads. Fast (chr17 at 10x in <12 min).

**Limitations**: **Strongly inflated error rates** in some datasets requiring recalibration. No variant injection. No cancer/tumor modeling, UMI, cfDNA, PCR duplicates, or GC bias. Single-threaded. Unmaintained since 2016.

### art_modern
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/YU-Zhejian/art_modern |
| **Language** | C++17 |
| **Stars** | 15 stars; 627 commits |
| **Status** | Actively maintained |

**Capabilities**: Community re-implementation of ART with MPI support, BAM/SAM output, RNA-Seq mode, multiple coverage modes, container support, modern build system.

### NEAT
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/ncsa/NEAT |
| **Language** | Python 3.10+ |
| **Stars** | 67 stars; v4.3.6 |
| **Status** | Actively maintained |

**Capabilities**: Learns error/mutation models from real datasets. Variant injection from VCF. SNPs, indels, inversions. Paired-end/single-end. Multi-threaded (v4.3.5+, 2-3x speedup).

**Limitations**: Python -- slower than compiled tools. No tumor/cancer modeling, UMI/duplex, cfDNA, or artifact simulation.

**Note**: **rusty-neat** (Rust reimplementation at github.com/ncsa/rusty-neat) exists at v1.2.0 (Sept 2025) with 2,506 commits, but incomplete -- no BAM output, no model learning, 20 open issues.

### DWGSIM
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/nh13/DWGSIM |
| **Language** | C (85%) |
| **Stars** | 103 stars; v0.1.16 (Nov 2025) |
| **Status** | Actively maintained |

**Capabilities**: Lightweight, fast WGS simulator. Configurable error rates (uniform/increasing/decreasing). Paired-end. Evaluation scripts for mapping/variant assessment.

**Limitations**: Simple error model (no empirical profiles). No cancer features, UMI, cfDNA, or artifacts.

### ReSeq
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/schmeing/ReSeq |
| **Language** | C++ |
| **Stars** | 51 stars |
| **Paper** | Genome Biology 2021 |
| **Status** | Maintained |

**Capabilities**: **Most realistic Illumina simulator.** Learns and reproduces systematic errors, fragment-based coverage, GC bias, PCR duplicates (negative binomial model), sampling matrices from real data.

**Limitations**: Slower than simpler simulators. No cancer/variant injection. No UMI or cfDNA.

### InSilicoSeq
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/HadrienG/InSilicoSeq |
| **Language** | Python |
| **Status** | Active |

**Capabilities**: Premade error models for MiSeq, HiSeq, NovaSeq, NextSeq. v2.0 adds amplicon simulation.

**Limitations**: Primarily metagenomics. No cancer features.

### wgsim
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/lh3/wgsim |
| **Language** | C |

Minimal, fast WGS read simulator with SNPs, indels, uniform errors. Very simple error model.

### NGSNGS
| Attribute | Detail |
|-----------|--------|
| **Paper** | Bioinformatics (Jan 2023) |
| **Status** | Published 2023 |

**Capabilities**: Multithreaded NGS read simulator. Faster than existing methods.

### Mason
| Attribute | Detail |
|-----------|--------|
| **Language** | C++ (SeqAn) |
| **Stars** | 8 stars |
| **Status** | Dormant (docs reference v2.0.0-beta1 from Feb 2014) |

Multi-tool suite (genome sim, VCF materializer, read sequencing, methylation, splicing). No cancer features, UMI, or cfDNA.

### Badread
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/rrwick/Badread |
| **Language** | Python |

Long-read simulator (ONT/PacBio) with chimeras and systematic errors. Not for Illumina.

---

## Cancer/Tumour Evolution Simulators

### VarSim
| Attribute | Detail |
|-----------|--------|
| **URL** | https://bioinform.github.io/varsim/ |
| **Language** | Java (70.6%), Python, Shell |
| **Stars** | 91 stars; v0.8.6 (Sept 2020) |
| **Status** | Maintenance mode |

**Capabilities**: End-to-end framework (variant generation + read simulation + validation). Leverages real variant databases (dbSNP, DGV). SNVs, indels, large SVs. Tumor/normal pairs. Built-in validation with interactive reports. Supports ART and DWGSIM backends.

**Limitations**: Complex setup (Java + Python + Maven + Ant). 36 open issues. No UMI, cfDNA, or artifacts. Heavy JVM dependency. Last release 2020.

### MOV&RSim (2025)
| Attribute | Detail |
|-----------|--------|
| **Paper** | BMC Bioinformatics (Nov 2025) |
| **Status** | New (2025) |

**Capabilities**: Cancer-specific simulator leveraging well-annotated variant databases. **Cancer-specific presets for 21 cancer types.** Full control over biological and technical parameters. Dockerized. Developed after analysis of 9 existing somatic simulators found none provided complete control.

**Limitations**: Very new, small community.

### GENOMICON-Seq (2024-2025)
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/Rounge-lab/GENOMICON-Seq |
| **Paper** | Scientific Reports (2025) |

**Capabilities**: **First framework to simulate full library prep pipeline** including PCR amplification errors, probe-capture enrichment, and Illumina sequencing biases. Models APOBEC3 edits and COSMIC SBS signatures. Tracks ground truth through entire pipeline. Docker-based.

### PSiTE
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/hchyang/PSiTE |
| **Language** | Python |

Phylogeny-guided simulator for tumour evolution. Modules for germline genome, somatic events, clone genomes, WGS/WES reads. Bulk and single-cell. Complex multi-module setup.

### HeteroGenesis + w-Wessim
WES-focused tumour evolution simulator with clonal architecture. Requires separate read simulator. Published Bioinformatics 2019.

### SVEngine
SV simulator with cancer clonal evolution. SV-focused only. Published GigaScience 2018.

### VISOR
Haplotype-aware SV simulator for short/linked/long reads. SV-focused, no small variants. Published 2019.

### CINner
R-based chromosomal instability modeler. Focal CNAs, missegregation, WGD. No raw reads.

### Synth4bench (2024)
Framework for benchmarking tumor-only somatic callers. Built on NEAT v3.3. Focuses on low-VAF variants (<=10%). Benchmarked Mutect2, FreeBayes, VarDict, VarScan2, LoFreq. R scripts for analysis.

---

## Mutation Signature Simulators

### SigProfilerSimulator
| Attribute | Detail |
|-----------|--------|
| **URL** | https://github.com/AlexandrovLab/SigProfilerSimulator |
| **Language** | Python |

Generates null-hypothesis mutational landscapes. Transforms real somatic catalogs maintaining burden and signature patterns. SBS, DBS, indel signatures. Operates at mutation catalog level, not read/BAM level.

### SomaticSiMu
Python-based. Generates substitutions and indels with biologically representative mutation signature probabilities. Genome-level, not read-level.

---

## UMI / Duplex / cfDNA Tools

### UMI-tools (simulation component)
Single-locus UMI amplification simulation. Not genome-wide.

### Kennedy Lab Duplex Pipeline
Analysis pipeline with basic test data generation. Not a dedicated simulator.

### shendurelab/cfDNA
Procedural cfDNA simulation using kmer frequencies and insert size distributions. Research-grade code. Fragmentomics focused, no variant spike-in.

### Fragmentstein
Converts tabulated fragment coordinates into paired-end BAM files. Does NOT generate fragments de novo. Published Bioinformatics 2024.

### FinaleToolkit
Analyzes real cfDNA fragmentation patterns. Does NOT simulate them. Published PMC 2024.

---

## Summary Comparison Matrix

| Tool | From Scratch | SNV | Indel | SV | CNV | Tumor Model | UMI | cfDNA | Artifacts | Stochastic VAF | Speed | Maintained |
|------|-------------|-----|-------|-----|-----|-------------|-----|-------|-----------|----------------|-------|------------|
| BAMSurgeon | No | ✓ | ✓ | ✓ | ✗ | Partial | ✗ | ✗ | ✗ | ✗ | Slow | Low |
| SomatoSim | No | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | Fast | No |
| stochasticSim | No | ✓ | ✗ | ✗ | ✗ | ✓ | ✗ | ✗ | ✓ | ✓ | Med | Low |
| ART | Yes | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | N/A | Fast | No |
| NEAT | Yes | ✓ | ✓ | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | Med | Yes |
| DWGSIM | Yes | ✓ | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | Fast | Yes |
| ReSeq | Yes | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✓* | N/A | Med | Yes |
| VarSim | Yes | ✓ | ✓ | ✓ | ✗ | Partial | ✗ | ✗ | ✗ | ✗ | Med | No |
| MOV&RSim | Yes | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | ✗ | ? | Med | New |
| GENOMICON-Seq | Yes | ✓ | ✗ | ✗ | ✗ | Partial | ✗ | ✗ | ✓ | ✗ | Med | New |
| PSiTE | Yes | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | ✗ | ✗ | Slow | Low |
| **VarForge** | **Yes** | **✓** | **✓** | **✓** | **✓** | **✓** | **✓** | **✓** | **✓** | **✓** | **Fast** | **Active** |

*ReSeq models PCR duplicates and GC bias but not FFPE/oxoG artifacts specifically.

---

## How Target Consumers Currently Benchmark

### Somatic Variant Callers (Mutect2, VarDict, Strelka2, VarScan2, FreeBayes, LoFreq)
- BAMSurgeon spike-ins into real BAMs (slow, deterministic bias)
- Genome in a Bottle (GIAB) NA12878/NA24385 mixtures
- SEQC2 HCC1395 truth set (~37K SNVs, ~58K indels)
- Synth4bench framework using NEAT-generated data
- PrecisionFDA Truth Challenge data
- Performance varies dramatically with mutation frequency: VarDict > Mutect2/Strelka2 at low VAFs

### UMI/Duplex Tools (fgbio, UMI-tools, HUMID)
- **Circular benchmarking**: fgbio itself annotates BAMs with synthetic UMI sequences
- Seraseq ctDNA reference materials (wet-lab controls)
- Ad hoc synthetic data from manual construction
- No standard simulation framework exists

### cfDNA Analysis (ichorCNA, DELFI)
- Real patient samples with known tumor fractions
- Low-coverage WGS at various depths
- ichorCNA: excels at high purity (>=50%), limited at low tumor fractions
- DELFI: ML on fragmentation patterns, validated with clinical data
- No in silico simulation tool exists

### Read Mappers (BWA, minimap2, samtools)
- ART and DWGSIM-generated reads (simple but fast)
- GIAB truth sets for accuracy assessment

---

## Key Takeaways for VarForge Design

1. **No existing tool occupies VarForge's planned niche**: integrated from-scratch simulation with cancer modeling + UMI + cfDNA + artifacts
2. **UMI/duplex simulation is completely unserved** -- VarForge would be first
3. **cfDNA simulation is completely unserved** for in silico generation
4. **Stochastic VAF sampling** is rare (only stochasticSim), VarForge already implements it
5. **Rust performance** would be a differentiator for WGS-scale workloads
6. **YAML config + modern CLI** would address the "poor documentation and integration" complaint
7. **MOV&RSim and GENOMICON-Seq** are the closest competitors but lack UMI, cfDNA, and Rust performance
8. **ReSeq** sets the bar for realism -- VarForge should aim for comparable or better fidelity
