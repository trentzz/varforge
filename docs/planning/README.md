# VarForge Implementation Plan

## Architecture Overview

VarForge is a single-binary Rust tool that generates synthetic cancer sequencing data from a YAML configuration. The simulation pipeline flows:

```
YAML Config → Reference Loading → Region Partitioning → Per-Region Read Generation → Output Writing
                                                              ↓
                                                    Fragment Sampling
                                                    Quality Modeling
                                                    Variant Spike-in
                                                    UMI Attachment
                                                    Artifact Injection
                                                    PCR Duplication
```

## Design Principles

1. **Config-driven**: Every simulation parameter is set in YAML, with CLI overrides
2. **Stochastic-first**: All sampling uses proper distributions (binomial VAF, log-normal families, etc.)
3. **Streaming I/O**: Never hold full genome in memory; process region-by-region
4. **Parallel by default**: rayon work-stealing across genomic regions
5. **Composable modules**: Each feature (UMI, cfDNA, artifacts) is independently toggleable
6. **Test-driven**: Every module has unit tests; integration tests validate end-to-end correctness

## Target Users & Use Cases

### Primary Users
1. **Bioinformatics pipeline developers**: Need ground-truth data to validate variant callers (Mutect2, VarDict, Strelka2), UMI tools (fgbio, HUMID), cfDNA analyzers (ichorCNA, DELFI)
2. **Clinical genomics labs**: Need synthetic data for assay validation without patient data (privacy, availability)
3. **Tool benchmarking researchers**: Need controlled experiments comparing caller sensitivity/specificity across VAF ranges, tumor purities, and coverage depths
4. **CI/CD pipelines**: Need fast, deterministic (seeded) small-dataset generation for automated testing

### Key Use Cases
- **Variant caller benchmarking**: Generate tumor/normal pairs with known mutations at controlled VAFs, run callers, compare to truth VCF
- **UMI tool validation**: Generate reads with UMI barcodes, known PCR families, and known true/artifactual variants to test consensus calling
- **Liquid biopsy assay development**: Generate cfDNA-like fragments with low tumor fractions (0.1-5%) to test ctDNA detection limits
- **FFPE artifact testing**: Generate reads with controlled FFPE damage rates to test artifact filtering
- **Panel design validation**: Generate targeted panel data with known variants to validate capture and calling
- **Sensitivity analysis**: Sweep across VAFs (0.1% to 50%), coverages (50x to 5000x), tumor purities to map caller performance landscapes
- **Longitudinal monitoring**: Generate time-series samples with shared clonal architecture and varying tumor fractions
- **Education/training**: Provide realistic but controlled data for teaching cancer genomics concepts

## Phase Structure

### Phase 1: Foundation I/O (Tasks 01-04)
Wire up the basic I/O: read a reference genome, write FASTQ, write BAM, write truth VCF. These are prerequisites for everything else.

### Phase 2: Core Pipeline (Tasks 05-07)
Build the read generation engine that ties all existing modules together, implement the `simulate` command orchestrator, and add parallel processing.

### Phase 3: Extended Variant Support (Tasks 08-10)
Add structural variants, copy number variants, and VCF input parsing for user-supplied mutation lists.

### Phase 4: Advanced Features (Tasks 11-14)
GC bias modeling, learned error profiles from real data, CLI flag overrides/presets, and simulation manifest output.

### Phase 5: Extended Simulation Modes (Tasks 15-18)
Multi-sample/longitudinal simulation, cancer-type presets with COSMIC signatures, the learn-profile subcommand, and capture/panel efficiency modeling.

### Phase 6: Quality & Polish (Tasks 19-21)
Integration tests, performance benchmarking/optimization, and documentation with example configs.

## Task Index

| Task | Name | Phase | Dependencies | File |
|------|------|-------|-------------|------|
| 01 | Reference Genome Loading | 1 | None | [task-01-reference-loading.md](tasks/task-01-reference-loading.md) |
| 02 | FASTQ Writer | 1 | None | [task-02-fastq-writer.md](tasks/task-02-fastq-writer.md) |
| 03 | BAM Writer | 1 | None | [task-03-bam-writer.md](tasks/task-03-bam-writer.md) |
| 04 | Truth VCF Writer | 1 | None | [task-04-truth-vcf-writer.md](tasks/task-04-truth-vcf-writer.md) |
| 05 | Read Generation Engine | 2 | 01 | [task-05-read-generation.md](tasks/task-05-read-generation.md) |
| 06 | Simulate Command Orchestrator | 2 | 01-05 | [task-06-simulate-command.md](tasks/task-06-simulate-command.md) |
| 07 | Parallel Region Processing | 2 | 05, 06 | [task-07-parallel-processing.md](tasks/task-07-parallel-processing.md) |
| 08 | Structural Variant Support | 3 | 05 | [task-08-structural-variants.md](tasks/task-08-structural-variants.md) |
| 09 | Copy Number Variant Support | 3 | 05 | [task-09-copy-number-variants.md](tasks/task-09-copy-number-variants.md) |
| 10 | VCF Input Parsing | 3 | 05 | [task-10-vcf-input-parsing.md](tasks/task-10-vcf-input-parsing.md) |
| 11 | GC Bias Modeling | 4 | 05 | [task-11-gc-bias.md](tasks/task-11-gc-bias.md) |
| 12 | Learned Error Profiles | 4 | 05 | [task-12-error-profiles.md](tasks/task-12-error-profiles.md) |
| 13 | CLI Overrides & Presets | 4 | 06 | [task-13-cli-overrides-presets.md](tasks/task-13-cli-overrides-presets.md) |
| 14 | Simulation Manifest Output | 4 | 06 | [task-14-manifest-output.md](tasks/task-14-manifest-output.md) |
| 15 | Multi-Sample Longitudinal | 5 | 06 | [task-15-multi-sample.md](tasks/task-15-multi-sample.md) |
| 16 | Cancer-Type Presets | 5 | 13 | [task-16-cancer-presets.md](tasks/task-16-cancer-presets.md) |
| 17 | Learn-Profile Subcommand | 5 | 12 | [task-17-learn-profile.md](tasks/task-17-learn-profile.md) |
| 18 | Capture Efficiency Model | 5 | 05 | [task-18-capture-efficiency.md](tasks/task-18-capture-efficiency.md) |
| 19 | Integration Tests | 6 | 06 | [task-19-integration-tests.md](tasks/task-19-integration-tests.md) |
| 20 | Performance Benchmarking | 6 | 07 | [task-20-performance.md](tasks/task-20-performance.md) |
| 21 | Documentation & Examples | 6 | All | [task-21-documentation.md](tasks/task-21-documentation.md) |
| 22 | BAM Editor (Spike-in to Existing BAM) | 3 | 05 | [task-22-bam-editor.md](tasks/task-22-bam-editor.md) |

## Loop Command

To execute tasks iteratively with subagents:

```bash
# From project root, run the task loop
./scripts/task-loop.sh
```

Each task is self-contained with its own context file, acceptance criteria, and test requirements. A subagent can pick up any task whose dependencies are satisfied and implement it independently.
