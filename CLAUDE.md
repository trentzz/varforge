# VarForge

Rust tool for generating synthetic cancer sequencing test data.

## Project Overview

VarForge generates realistic synthetic sequencing data (FASTQ/BAM) with controlled mutations, tumour parameters, UMI tags, and cfDNA fragment profiles for benchmarking bioinformatics tools.

## Research Documentation

All research and design documents are in `docs/research/`:

- `existing-tools.md` - Catalog of existing simulation tools (BAMSurgeon, ART, NEAT, etc.)
- `gap-analysis.md` - Gaps in current tooling that VarForge addresses
- `technical-requirements.md` - Implementation details: error profiles, VAF modeling, PCR duplicates, etc.
- `cancer-genomics-details.md` - Domain knowledge: variant types, cfDNA, UMI/duplex sequencing
- `rust-ecosystem.md` - Rust crates for BAM/FASTQ/VCF handling
- `input-specification.md` - Proposed YAML configuration format

## Key Design Decisions

- **Rust** for performance (high-coverage WGS simulation is compute-intensive) and single-binary distribution
- **YAML primary config** with optional VCF for mutation lists
- **rust-htslib** for BAM/VCF I/O (mature, wraps htslib C library)
- **noodles-fasta** for reference genome access (pure Rust)
- **seq_io** for fast FASTQ writing

## Target consumers of generated data

km, samtools, HUMID, fgbio, UMI-tools, Mutect2, VarDict, Strelka2, CNVkit, and similar tools.
