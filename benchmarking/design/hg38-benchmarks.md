# hg38 Benchmark Design

## Purpose

The synthetic reference benchmarks (1 MB, 10 MB) measure scaling behaviour in a controlled environment.
The hg38 benchmarks measure correctness and feature coverage on a real genome with realistic GC content, repeat structure, and chromosome length.
This closes the gap between synthetic benchmarks and what users will actually encounter.

## Reference

Chromosome 22 from hg38 (~51 Mbp).
Chosen because it is the smallest autosome, large enough to be realistic, and small enough to complete benchmarks in reasonable time on a 12 GiB system.

## Coverage note

All configs run at 5x coverage.
The purpose of these benchmarks is to confirm that each feature (UMI, duplex, cfDNA, FFPE) runs without error on a real reference and produces plausible output.
Throughput and scaling data are derived from the synthetic benchmarks.
Higher clinical depths (500x for panels, 200x for cfDNA) require BED-based region filtering, which is not yet implemented.

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
