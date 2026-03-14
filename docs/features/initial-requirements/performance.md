# Performance Requirements

Throughput, memory, and scalability targets.

## Parallelization

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-PRF-001 | Partition the genome into independent regions and simulate in parallel using rayon | P0 |
| REQ-PRF-002 | Region boundaries must not split variant sites or read pairs across workers | P0 |
| REQ-PRF-003 | Scale near-linearly up to at least 16 threads | P0 |
| REQ-PRF-004 | Thread count configurable at runtime (default: all available cores) | P0 |

## Streaming Output

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-PRF-010 | Stream reads to output files as they are generated; do not accumulate the full dataset in memory | P0 |
| REQ-PRF-011 | Use bounded channel or buffer between generation and I/O threads to decouple compute from disk | P0 |
| REQ-PRF-012 | BAM sorting via external merge sort if BAM output is requested (avoid holding all records in memory) | P1 |

## Reference Genome Access

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-PRF-020 | Use memory-mapped or indexed access to the reference FASTA (do not load the entire genome into RAM) | P0 |
| REQ-PRF-021 | Support .fai indexed FASTA files | P0 |

## Throughput Targets

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-PRF-030 | Generate 30× WGS paired-end data (human genome) in under 60 minutes on a 16-core machine | P0 |
| REQ-PRF-031 | Generate a targeted panel simulation (1000 genes, 500× depth) in under 5 minutes | P1 |

## Memory

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-PRF-040 | Peak memory usage for 30× WGS under 8 GB | P0 |
| REQ-PRF-041 | Targeted panel simulations should use under 2 GB | P1 |

## Compression

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-PRF-050 | Multi-threaded BAM compression using the htslib thread pool (via `rust-htslib`) | P1 |
| REQ-PRF-051 | Multi-threaded FASTQ gzip compression using `gzp` (pigz-style block compression with libdeflater backend) | P1 |

## Crate Leverage Policy

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-PRF-060 | Use existing Rust crates for all I/O formats, compression, distributions, parallelism, CLI, and logging — do not reimplement functionality available in the ecosystem | P0 |
| REQ-PRF-061 | Only build custom implementations for domain-specific logic with no existing crate: GC bias models, error profiles, cfDNA fragmentation, mixture distributions, variant spike-in, tumour models, UMI simulation, and artifact injection | P0 |
