# VarForge Benchmark Design

## What matters for a sequencing data simulator

VarForge generates synthetic cancer sequencing data. Users care about three things: how fast it runs, how much memory it needs, and whether the output is correct. The benchmarks test all three.

## Experiments

### 1. Coverage scaling (throughput and memory)

**Hypothesis**: Wall time and peak memory scale linearly with coverage depth. The streaming architecture keeps per-pair memory constant regardless of total dataset size.

**What it measures**: Read pairs per second, wall time, and peak RSS across seven coverage levels (1x to 200x) on a fixed 1 MB reference with 100 mutations.

**What it does not measure**: I/O throughput at WGS scale. The 1 MB reference keeps the workload compute-bound, isolating the simulation engine from storage bottlenecks.

**Baseline**: Each coverage level is its own baseline. The linearity of the curve is the result.

### 2. Feature overhead

**Hypothesis**: Each feature adds a bounded, predictable cost. Variant injection and GC bias are cheap. UMI simulation is the most expensive feature because it tracks PCR family genealogies.

**What it measures**: Wall time for each feature in isolation and all features combined, relative to a no-feature baseline. Same reference and coverage across all configs.

**What it does not measure**: Interaction effects beyond the all-combined case. Two-way feature interactions are not tested individually.

**Baseline**: The no-feature configuration (1 MB, 30x, no mutations, no UMI, no artifacts, no GC bias).

### 3. Thread scaling

**Hypothesis**: Scaling depends on the compute-to-I/O ratio. I/O-bound workloads (small reference, low coverage) plateau early. Compute-bound workloads (high coverage, many features) scale better.

**What it measures**: Wall time from 1 to 12 threads on a fixed config. Speedup and parallel efficiency.

**What it does not measure**: Scaling beyond 12 threads. The test system has 6 cores / 12 threads.

**Baseline**: Single-thread execution of the same config.

### 4. VAF accuracy

**Hypothesis**: VarForge's stochastic VAF sampling produces alt-read counts that follow the expected Binomial distribution. This is the core methodological contribution over deterministic tools like BAMSurgeon.

**What it measures**: For variants at known VAF and known coverage, the observed alt-read count distribution. Compare the empirical distribution to Binomial(D, VAF) using a chi-squared goodness-of-fit test.

**What it does not measure**: Whether the VAF values in the truth VCF match the configured target. That is a configuration validation test, not a statistical accuracy test.

**Baseline**: The theoretical Binomial(D, VAF) distribution.

### 5. cfDNA fragment distribution

**Hypothesis**: The cfDNA fragment model produces fragment lengths matching the four-component Gaussian mixture with 10 bp periodicity below 167 bp.

**What it measures**: Fragment length histogram from cfDNA simulation output. Visual comparison against the expected mixture density. Peak positions at ~143, ~167, and ~334 bp.

**What it does not measure**: Whether the ctDNA shift produces the correct differential enrichment below 167 bp. That requires comparing tumour-fraction-specific fragment pools.

**Baseline**: The parametric mixture model defined in the paper (Equation 3).

### 6. Output size

**Hypothesis**: Compressed FASTQ output size is predictable from coverage and read length.

**What it measures**: FASTQ.gz file sizes across configurations. Compression ratio relative to uncompressed size.

**What it does not measure**: BAM file sizes (BAM output is optional and not the primary output format).

## Hardware

AMD Ryzen 5 7640HS (6 cores / 12 threads), 12 GiB DDR5, NVMe SSD, Debian 13, Rust 1.94.0.

## Reproducibility

All runs use `--seed N` where N is the iteration number. Scripts in `benchmarking/scripts/` produce all results. Raw timing data is stored in `benchmarking/results/tables/`.
