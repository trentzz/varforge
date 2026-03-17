# hg38 Benchmark Design

## Purpose

The synthetic reference benchmarks (1 MB, 10 MB) measure scaling behaviour in a controlled environment.
The hg38 benchmarks measure correctness and feature coverage on a real genome with realistic GC content, repeat structure, and chromosome length.
This closes the gap between synthetic benchmarks and what users will actually encounter.

## Reference

Chromosome 22 from hg38 (~51 Mbp).
Chosen because it is the smallest autosome, large enough to be realistic, and small enough to complete benchmarks in reasonable time on a 12 GiB system.

## Coverage tiers

**Tier 1 (feature validation, 5x)**: Confirms each feature runs without error on a real reference. Low coverage completes quickly and uses minimal RAM.

**Tier 2 (realistic WGS, 30x)**: Full chr22 at clinical WGS depth. ~5M read pairs, ~4 GB peak RAM, ~10 min. Validates throughput at real scale and confirms bounded memory holds at realistic pair counts.

**Tier 3 (targeted panel with BED, 200x)**: Restricted to a synthetic 1 Mbp target region using BED filtering. High coverage on a small footprint is now feasible since BED filtering was implemented in v0.1.0. Simulates a realistic panel run without simulating non-target regions.

The synthetic benchmarks (1 MB, 10 MB) remain the primary throughput and scaling evidence. The hg38 benchmarks confirm correct behaviour on real sequence and add a realistic RAM data point at WGS scale.

## Memory guidance for constrained systems

VarForge's peak RAM scales with the number of in-flight read pairs in the streaming channel, not with total pairs generated.
For systems with limited RAM, two settings reduce memory:
- `performance.output_buffer_regions`: reduce from the default of 64 to 8 or 16.
- `regions_bed`: restrict simulation to a target region, reducing the total work and the steady-state buffer occupancy.
- Fewer threads: each rayon worker holds a batch of reads; fewer workers means fewer concurrent batches.

## Experiments

### 1. WGS baseline (5x)

Standard whole-genome simulation on chr22.
No optional features.
Measures baseline behaviour on a real reference.

### 2. WGS with variants (5x, 500 mutations)

Same as baseline but with 500 random somatic mutations.
Confirms variant injection works on real sequence with realistic GC content.

### 3. Targeted panel (5x, 50 mutations, UMI simplex)

Simulates a targeted panel sequencing run.
UMI simplex, 8-mer, inline, 10 PCR cycles.
Confirms UMI family generation and tagging work on chr22.

### 4. Twist-style UMI duplex (5x, 50 mutations)

Simulates Twist Biosciences UMI duplex protocol.
5-mer UMI (matching Twist read structure), duplex mode.
Confirms dual-strand tagging works on a real reference.

### 5. cfDNA liquid biopsy (5x, 200 mutations, 2% purity)

cfDNA fragment model on a real reference.
Low tumour fraction to test the sensitivity use case.
Confirms nucleosomal fragment length distribution is applied to real sequence.

### 6. FFPE tumour (5x, 500 mutations, FFPE + oxoG)

Simulates a formalin-fixed tumour sample.
Confirms FFPE deamination and oxidative damage are applied to real sequence.

## Metrics

For each experiment:
- Wall time (seconds)
- Peak RSS (MB)
- Read pairs generated
- Output file size (FASTQ.gz, bytes)
- Throughput (read pairs per second)

## Hardware

Same as synthetic benchmarks: AMD Ryzen 5 7640HS, 12 GiB DDR5, NVMe SSD.
