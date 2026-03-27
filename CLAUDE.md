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
cargo clippy --all-targets -- -D warnings
cargo test
```

## Target consumers

km, samtools, HUMID, fgbio, UMI-tools, Mutect2, VarDict, Strelka2, CNVkit.

---

## Task Tracking

All development tasks are tracked in `docs/claudetracking/KANBAN.md` as a kanban board
with four columns: TODO, IN_PROGRESS, DONE, BLOCKED.

### Rules

- Before starting any task, move it from TODO to **IN_PROGRESS** in the kanban.
- Mark **DONE** only after `cargo fmt -- --check`, `cargo clippy --all-targets -- -D warnings`, and
  `cargo test` all pass cleanly.
- If a task is blocked, move it to **BLOCKED** with a one-line note explaining why.
- Check the `Depends on` column before starting: all listed dependencies must be DONE first.
- When running parallel agents, each agent is responsible for moving its own tasks to
  the correct column when it finishes.

### Format

Each task row: `ID | Title | Priority | Depends on | Notes`

Priority levels: P0 (critical bug), P1 (high), P2 (medium), P3 (low).

Task IDs are stable: never reuse an ID, even after a task is deleted.
