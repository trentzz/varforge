# I/O Format Requirements

Input configuration and output file formats.

## Input: YAML Configuration

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-IOF-001 | Accept a YAML file as the primary simulation configuration | P0 |
| REQ-IOF-002 | Validate YAML against a schema on load; emit clear error messages for invalid fields | P0 |
| REQ-IOF-003 | Support a minimal config (reference + output path + coverage) with sensible defaults for all other parameters | P0 |

## Input: VCF Mutation List

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-IOF-010 | Accept an optional VCF file specifying mutations to spike in | P0 |
| REQ-IOF-011 | Read standard VCF fields (CHROM, POS, REF, ALT) for variant definition | P0 |
| REQ-IOF-012 | Support custom INFO fields for clone assignment (`CLONE`), target VAF override (`VAF`), and copy number context (`CN`) | P1 |

## Input: BED Regions

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-IOF-020 | Accept a BED file defining targeted capture regions (parsed via `noodles-bed` or `rust-bio`) | P0 |
| REQ-IOF-021 | Support an optional fourth column for per-region depth multipliers | P2 |

## Output: FASTQ

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-IOF-030 | Write gzipped paired-end FASTQ files (R1 + R2) using `seq_io` for writing and `gzp` for multi-threaded compression | P0 |
| REQ-IOF-031 | Configurable gzip compression level (default 6) | P1 |
| REQ-IOF-032 | Multi-threaded gzip compression via `gzp` with libdeflater backend (pigz-style block parallelism) | P1 |

## Output: BAM

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-IOF-040 | Optionally write a coordinate-sorted, indexed BAM file | P0 |
| REQ-IOF-041 | Include read group, mate information, and proper CIGAR strings | P0 |
| REQ-IOF-042 | Multi-threaded BAM compression via rust-htslib thread pool | P1 |
| REQ-IOF-043 | Write BAM index (.bai) alongside the BAM file | P0 |

## Output: Truth VCF

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-IOF-050 | Emit a truth VCF containing all spiked variants | P0 |
| REQ-IOF-051 | Include INFO fields: expected VAF, achieved VAF, clone, CCF, CN state, REF/ALT read counts | P0 |
| REQ-IOF-052 | Truth VCF must be coordinate-sorted and bgzipped with tabix index | P1 |

## Output: Simulation Manifest

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-IOF-060 | Write a JSON manifest summarizing the simulation run | P0 |
| REQ-IOF-061 | Manifest includes: VarForge version, effective random seed, input config hash, output file paths, total reads generated, runtime, and achieved mean coverage | P0 |

## Learned Error Profiles

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-IOF-070 | `learn-profile` writes a binary or JSON profile file capturing position-dependent error rates, quality distributions, and context-dependent substitution frequencies | P1 |
| REQ-IOF-071 | The `simulate` command accepts the profile file as an input parameter | P1 |
