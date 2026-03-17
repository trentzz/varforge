# VarForge

Rust tool for generating synthetic cancer sequencing test data.

## Overview

VarForge generates realistic synthetic sequencing data (FASTQ/BAM) with controlled mutations, tumour parameters, UMI tags, and cfDNA fragment profiles for benchmarking bioinformatics tools.

## Project Structure

```
varforge/
├── src/                    # Implementation
├── tests/                  # Integration tests
├── examples/               # Runnable YAML configs
├── benchmarking/           # Performance evaluation
│   ├── design/             # What is measured and why
│   ├── micro/              # Criterion micro-benchmarks (cargo bench)
│   ├── results/
│   │   ├── graphs/         # Publication figures (PDF/PNG)
│   │   └── tables/         # Raw data (JSON/CSV)
│   └── scripts/            # run_benchmarks.sh, plot_results.py
└── docs/
    ├── research/           # Literature context and gap analysis
    ├── features/           # Feature specifications
    ├── planning/           # Milestones, tasks, decision log
    ├── reviewcycle/        # Review notes and responses
    ├── paper/              # LaTeX source (imports graphs from benchmarking/)
    │   ├── main.tex
    │   ├── references.bib
    │   └── sections/       # One .tex per section
    └── archive/            # Superseded versions
```

## Key Design Decisions

- **Pure Rust** dependency stack (noodles for BAM/VCF/FASTA, seq_io for FASTQ, flate2 for compression). No C libraries.
- **YAML primary config** with optional VCF for mutation lists.
- **Streaming architecture**: bounded crossbeam-channel between rayon workers and a dedicated writer thread. Memory scales with channel depth, not dataset size.
- **Stochastic VAF**: per-read Bernoulli sampling, not deterministic spike-in.

## Paper PDFs

All final PDFs live in `docs/paper/pdfs/`. Naming scheme: `varforge_v{version}_{type}.pdf`.
Types: `normal` (single-col), `normal_2col` (two-col, canonical format), `big` (single-col with ToC), `2p` (two-page summary).
Delete old versioned PDFs when rebuilding. Do not accumulate stale versions.

The two-column format (`main_2col.tex`) is the canonical working format.

## Pre-commit checks

Before committing, run and fix all failures:

```
cargo fmt -- --check
cargo clippy -- -D warnings
cargo test
```

## Target consumers

km, samtools, HUMID, fgbio, UMI-tools, Mutect2, VarDict, Strelka2, CNVkit.
