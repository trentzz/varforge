# Existing Tools for Synthetic Cancer Sequencing Data

## Mutation Spike-in Tools (BAM modification)

### BAMSurgeon
- **URL**: https://github.com/adamewing/bamsurgeon
- **Language**: Python (wraps samtools/picard/BWA)
- **What it does**: Adds SNVs, indels, and SVs to existing BAM files at controlled VAFs. Gold standard for the ICGC-TCGA DREAM Somatic Mutation Calling Challenge.
- **I/O**: BAM in -> BAM + truth VCF out
- **Maintained**: Low activity; mature
- **Limitations**:
  - **Slow**: ~56 minutes for runs SomatoSim does in ~44 seconds
  - **Spike-in bias**: A 10% VAF at 30x depth always gets exactly 3 alt reads instead of stochastic sampling -- unrealistic determinism
  - No cfDNA fragment simulation
  - No UMI support

### SomatoSim
- **URL**: https://github.com/BieseckerLab/SomatoSim
- **Language**: Python
- **What it does**: Simulates somatic SNVs in BAM files with control over position, count, and VAF. Models forward/reverse strand distribution from original data.
- **I/O**: BAM in -> BAM + truth VCF out
- **Maintained**: Published 2021; moderate activity
- **Limitations**: SNVs only (no indels, SVs, CNVs). Same deterministic spike-in as BAMSurgeon.

### stochasticSim
- **URL**: https://github.com/BrianOSullivanGit/stochasticSim
- **Paper**: NAR Cancer 2023
- **What it does**: Comprehensive tumour simulation based on phased, personalized genomes. Unlike BAMSurgeon, accounts for stochastic sampling errors (avoids spike-in bias). Tracks the source of all non-reference bases.
- **I/O**: Phased genome + variant specs in -> simulated data + truth out
- **Maintained**: Published 2023
- **Limitations**: Newer tool, smaller community.

### Xome-Blender
- **Paper**: PLOS ONE 2018
- **What it does**: Mixes real tumour and normal BAMs at specified ratios to create synthetic cancer genomes with user-defined features (subclone count, CNA count).
- **I/O**: BAM files in -> mixed BAM + truth VCF out
- **Limitations**: Requires existing real tumour/normal BAMs; doesn't generate de novo.

---

## Read Simulators (FASTQ generation)

### ART
- **URL**: https://www.niehs.nih.gov/research/resources/software/biostatistics/art
- **Modernized fork**: https://github.com/YU-Zhejian/art_modern
- **Language**: C++
- **What it does**: Generates synthetic Illumina, 454, and SOLiD reads using empirical error models from large real datasets. Outputs ground-truth alignments.
- **I/O**: FASTA in -> FASTQ + ALN/SAM out
- **Maintained**: Original is mature/static. `art_modern` actively maintained.
- **Limitations**: No variant injection. No UMI or cfDNA support.

### DWGSIM
- **URL**: https://github.com/nh13/DWGSIM
- **Language**: C
- **What it does**: WGS simulator with SNP/indel polymorphisms and uniform substitution errors.
- **I/O**: FASTA in -> FASTQ + truth out
- **Limitations**: Simpler error model than ART. No SV support.

### NEAT (NExt-generation sequencing Analysis Toolkit)
- **URL**: https://github.com/ncsa/NEAT (v4.x)
- **Language**: Python
- **What it does**: Fine-grained read simulator using models learned from real datasets. Simulates substitution errors, GC-bias coverage, mutation models.
- **I/O**: FASTA + optional VCF in -> FASTQ + BAM + VCF out
- **Maintained**: Yes (v4.3.5+)
- **Limitations**: Historically hard-coded error models requiring calibration.

### InSilicoSeq
- **URL**: https://github.com/HadrienG/InSilicoSeq
- **Language**: Python
- **What it does**: Produces realistic Illumina reads with premade error models for MiSeq, HiSeq, NovaSeq, NextSeq. v2.0 adds amplicon simulation.
- **I/O**: FASTA in -> FASTQ out
- **Limitations**: Primarily designed for metagenomics. No cancer features.

### wgsim
- **URL**: https://github.com/lh3/wgsim
- **Language**: C
- **What it does**: Minimal, fast WGS read simulator with SNPs, indels, and uniform errors.
- **Limitations**: Very simple error model.

### ReSeq
- **URL**: https://github.com/schmeing/ReSeq
- **Paper**: Genome Biology 2021
- **Language**: C++
- **What it does**: High-fidelity simulator that learns and reproduces systematic errors and fragment-based coverage from real data.
- **I/O**: BAM (training) + FASTA in -> FASTQ out
- **Limitations**: Slower than simpler simulators.

### Mason
- **Language**: C++ (SeqAn library)
- **What it does**: Read simulator for second-gen sequencing with ground-truth alignments.
- **Limitations**: Tied to SeqAn ecosystem.

### Badread
- **URL**: https://github.com/rrwick/Badread
- **Language**: Python
- **What it does**: Long-read simulator (ONT/PacBio) with chimeras, low-quality regions, systematic errors.
- **Limitations**: Long reads only; not for Illumina.

---

## Tumour Evolution / Heterogeneity Simulators

### HeteroGenesis + w-Wessim
- **Paper**: Bioinformatics 2019
- **What it does**: Generates realistically evolved tumour genomes with clonal evolution. Paired with w-Wessim for WES simulation. Supports multi-region and longitudinal data.
- **Limitations**: WES-focused; requires separate read simulator.

### PSiTE (Phylogeny-guided Simulator for Tumor Evolution)
- **URL**: https://github.com/hchyang/PSiTE
- **Language**: Python
- **What it does**: Simulates NGS data from tumour evolution guided by phylogenetic trees. Modules for germline genome, somatic events, clone genomes, WGS/WES reads. Supports bulk and single-cell.
- **I/O**: Phylogeny + params in -> FASTA + FASTQ + VCF + BAM out
- **Limitations**: Complex multi-module setup.

### SVEngine
- **Paper**: GigaScience 2018
- **What it does**: SV simulator with cancer clonal evolution features. Simulates locus-specific VAF to mimic phylogeny.
- **Limitations**: SVs only; not for small variants.

### VISOR
- **URL**: https://github.com/davidebolo1993/VISOR
- **Language**: Python
- **What it does**: Haplotype-aware SV simulator for short/linked/long reads. Supports complex SVs, diploid/polyploid, bulk and single-cell.
- **Limitations**: SV-focused; no small variants or clonal evolution.

### VarSim
- **URL**: https://bioinform.github.io/varsim/
- **Language**: Java + Python
- **What it does**: High-fidelity framework synthesizing diploid genomes with germline and somatic mutations from realistic models. Automates simulation + validation.
- **Limitations**: Complex setup. Depends on external read simulators.

### CINner
- **URL**: https://github.com/dinhngockhanh/CINner
- **Language**: R
- **What it does**: Models chromosomal instability at single-cell resolution. Focal CNAs, missegregation, WGD, drivers.
- **Limitations**: Focused on CIN/CNA; no raw reads.

---

## Mutation Signature / Catalog Simulators

### SigProfilerSimulator
- **URL**: https://github.com/AlexandrovLab/SigProfilerSimulator
- **Language**: Python
- **What it does**: Generates null-hypothesis mutational landscapes. Transforms real somatic catalogs into simulated ones maintaining burden and signature patterns. SBS, DBS, indel signatures.
- **I/O**: VCF/MAF in -> simulated catalogs out
- **Limitations**: Operates at mutation catalog level, not read/BAM level.

### SomaticSiMu
- **URL**: https://github.com/HillLab/SomaticSiMu
- **Language**: Python
- **What it does**: Generates substitutions and indels with biologically representative mutation signature probabilities.
- **Limitations**: Genome-level mutation, not read-level.

---

## UMI / Duplex / cfDNA Simulators

### UMI-tools (simulation component)
- **URL**: https://umi-tools.readthedocs.io/
- **Language**: Python
- **What it does**: Primarily UMI processing, but includes limited simulation of UMI amplification for benchmarking at a single locus.
- **Limitations**: Single-locus only; not genome-wide.

### Duplex Sequencing Pipeline (Kennedy Lab)
- **URL**: https://github.com/Kennedy-Lab-UW/Duplex-Seq-Pipeline
- **Language**: Python
- **What it does**: Analysis pipeline for duplex sequencing. Includes test data generation with simulated PCR errors.
- **Limitations**: Analysis pipeline, not a dedicated simulator.

### shendurelab/cfDNA
- **URL**: https://github.com/shendurelab/cfDNA
- **Language**: Python
- **What it does**: Procedural simulation of aligned cfDNA data. Models fragment initiation using kmer frequencies, samples length from insert size distribution (~160-170bp).
- **Limitations**: Research-grade code. Focused on fragmentomics, not variant spike-in.

---

## Benchmarking Frameworks and Truth Sets

### ICGC-TCGA DREAM Challenge
- **Paper**: Nature Methods 2015
- **What it does**: Crowdsourced benchmark using BAMSurgeon-generated synthetic tumours. 248 SNV submissions, 204 SV submissions evaluated.
- **Data**: Available via EGA (EGAS00001002092)

### Genome in a Bottle (GIAB)
- Provides high-confidence truth sets for reference samples (HG001-HG007)
- Somatic benchmarking via NA12878/NA24385 mixtures at known ratios
- Maintained by NIST

### SEQC2 Consortium
- **Reference**: HCC1395 (TNBC cell line) + HCC1395BL (matched normal)
- **Truth set**: ~37,000 SNVs + ~58,000 indels via multi-platform consensus
- Mixtures at controlled tumour fractions available

### nf-core/readsimulator
- **URL**: https://github.com/nf-core/readsimulator
- Nextflow pipeline wrapping multiple read simulators for amplicon, target capture, metagenome, and WGS data.
