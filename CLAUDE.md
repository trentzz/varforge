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

## Feature Requirements

Initial requirements documentation is in `docs/features/initial-requirements/`:

- `core-simulation.md` - Read generation, fragment models, coverage, quality scores
- `variant-spike-in.md` - SNV, indel, MNV, SV, CNV injection with stochastic VAF
- `tumour-model.md` - Purity, clonal architecture, ploidy, copy number
- `umi-duplex.md` - UMI barcodes, simplex/duplex modes, PCR families
- `liquid-biopsy.md` - cfDNA fragmentation, ctDNA enrichment, low TF support
- `artifacts.md` - FFPE damage, oxoG, GC bias, PCR duplicates/errors
- `io-formats.md` - YAML config, VCF/BED input, FASTQ/BAM/truth VCF output
- `cli-interface.md` - Subcommands, threading, progress, logging
- `performance.md` - Parallelization, streaming, memory, throughput targets

## Key Design Decisions

- **Rust** for performance (high-coverage WGS simulation is compute-intensive) and single-binary distribution
- **YAML primary config** with optional VCF for mutation lists
- **rust-htslib** for BAM/VCF I/O (mature, wraps htslib C library)
- **noodles-fasta** for reference genome access (pure Rust)
- **seq_io** for fast FASTQ writing

## Target consumers of generated data

km, samtools, HUMID, fgbio, UMI-tools, Mutect2, VarDict, Strelka2, CNVkit, and similar tools.
