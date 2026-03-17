# hg38 Benchmark Design

## Purpose

The synthetic reference benchmarks (1 MB, 10 MB) measure scaling behaviour in a controlled environment.
The hg38 benchmarks measure performance on a real genome with realistic GC content, repeat structure, and chromosome length.
This closes the gap between synthetic benchmarks and what users will actually encounter.

## Reference

Chromosome 22 from hg38 (~51 Mbp).
Chosen because it is the smallest autosome, large enough to be realistic, and small enough to complete benchmarks in reasonable time on a 12 GiB system.

## Experiments

### 1. WGS baseline (30x)

Standard whole-genome simulation on chr22.
No optional features.
Measures baseline throughput on a real reference.

### 2. WGS with variants (30x, 500 mutations)

Same as baseline but with 500 random somatic mutations.
Measures variant injection overhead on a real reference.

### 3. Targeted panel (500x, 50 mutations, UMI simplex)

Simulates a targeted panel sequencing run.
UMI simplex, 8-mer, inline, 10 PCR cycles.
Restricted to a subset of chr22 via regions.

### 4. Twist-style UMI duplex (1000x, 50 mutations)

Simulates Twist Biosciences UMI duplex protocol.
5-mer UMI (matching Twist read structure), duplex mode.
High depth to test UMI family generation at scale.

### 5. cfDNA liquid biopsy (200x, 200 mutations, 2% purity)

cfDNA fragment model on a real reference.
Low tumour fraction to test sensitivity use case.

### 6. FFPE tumour (60x, 500 mutations, FFPE + oxoG)

Simulates a formalin-fixed tumour sample.
Typical clinical coverage depth.
Tests artefact injection on real sequence.

### 7. Tumour-normal pair (60x tumour, 30x normal, 1000 mutations)

Multi-sample mode producing matched pair.
Tests the longitudinal/multi-sample pipeline.

## Metrics

For each experiment:
- Wall time (seconds)
- Peak RSS (MB)
- Read pairs generated
- Output file size (FASTQ.gz, bytes)
- Throughput (read pairs per second)

## Hardware

Same as synthetic benchmarks: AMD Ryzen 5 7640HS, 12 GiB DDR5, NVMe SSD.
