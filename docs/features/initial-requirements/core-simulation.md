# Core Simulation Requirements

Read generation engine that produces realistic sequencing data from a reference genome.

## FASTQ Generation

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-SIM-001 | Generate paired-end FASTQ reads from an indexed reference FASTA | P0 |
| REQ-SIM-002 | Support configurable read lengths (50–300 bp, default 150 bp) | P0 |
| REQ-SIM-003 | Assign read names encoding tile, lane, and pair information compatible with Illumina naming conventions | P1 |

## Fragment Size Models

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-SIM-010 | Standard WGS fragment size model: Gaussian distribution with configurable mean and standard deviation (default mean=300, sd=50) | P0 |
| REQ-SIM-011 | cfDNA fragment size model: mixture distribution with mononucleosomal peak (~167 bp), dinucleosomal peak (~334 bp), and 10 bp periodicity | P0 |
| REQ-SIM-012 | Custom mixture model: user-specified weighted Gaussian components | P1 |
| REQ-SIM-013 | ctDNA shorter fragment enrichment mode: shift tumour-derived fragments toward ~143 bp | P0 |

## Base Quality Scores

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-SIM-020 | Parametric quality model: position-dependent quality decay with configurable mean quality, tail slope, and read-position curve | P0 |
| REQ-SIM-021 | Learned quality model: load an empirical error profile derived from a real BAM via `varforge learn-profile` | P1 |
| REQ-SIM-022 | Inject base-call errors at rates consistent with the quality scores (Phred-scaled probability) | P0 |
| REQ-SIM-023 | Support platform-specific error profiles: Illumina (substitution-dominated), with hooks for future Nanopore/PacBio profiles | P2 |

## Coverage Control

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-SIM-030 | Uniform coverage mode: target a specified mean depth genome-wide | P0 |
| REQ-SIM-031 | GC bias mode: modulate local coverage as a function of regional GC content using a configurable bias curve | P1 |
| REQ-SIM-032 | Targeted panel mode: restrict read generation to regions defined in a BED file, with configurable on-target/off-target depth ratio | P0 |
| REQ-SIM-033 | Per-region depth override: allow a BED file with per-interval depth multipliers | P2 |

## Small Dataset / Quick Mode

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-SIM-035 | Chromosome subset mode: restrict simulation to a list of chromosomes (e.g., `chromosomes: [chr22]`) for fast iteration | P0 |
| REQ-SIM-036 | Region-only mode: simulate reads only within BED-defined regions with no off-target reads, for minimal output size | P0 |
| REQ-SIM-037 | Low-coverage quick mode: allow fractional coverage (e.g., `coverage: 0.1`) for smoke-testing pipelines | P0 |
| REQ-SIM-038 | Synthetic mini-reference mode: generate a small synthetic reference (e.g., 1 Mb) with spiked variants, avoiding the need for a full human genome download | P1 |
| REQ-SIM-039 | Preset profiles: `--preset small` (chr22, 30×), `--preset panel` (BED, 500×), `--preset wgs` (full genome, 30×) for common use cases | P1 |

## Read Groups

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-SIM-040 | Assign read group tags (ID, SM, LB, PL, PU) to all generated reads | P0 |
| REQ-SIM-041 | Support multiple read groups within a single simulation (e.g., multi-lane) | P1 |

## Reproducibility

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-SIM-050 | Accept a user-specified random seed for deterministic output | P0 |
| REQ-SIM-051 | Record the effective seed in the simulation manifest so runs can be reproduced | P0 |
| REQ-SIM-052 | Guarantee identical output for identical config + seed + VarForge version, regardless of thread count | P0 |
