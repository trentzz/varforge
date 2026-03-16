# Benchmark Results

Hardware, software, and methodology details are in
[hardware.md](hardware.md) and [methodology.md](methodology.md).

All raw numbers are also available in machine-readable form in
`benchmark_output/results.json`.

---

## Main Scenario Benchmarks

All runs use all available cores (12) unless noted.

| Config | Ref | Coverage | Variants | Features | Iters | Avg Wall (s) | Min Wall (s) | Avg Peak RAM (MB) |
|--------|-----|----------|----------|----------|-------|-------------|-------------|-------------------|
| 01_baseline | 10 MB | 30x | none | none | 3 | 75.4 | 74.7 | 3,421 |
| 02_with_variants | 10 MB | 30x | 500 | none | 3 | 75.9 | 75.2 | 3,423 |
| 03_high_coverage | 10 MB | 100x | 500 | none | 1 | 251.1 | 251.1 | 11,352 |
| 04_very_high_coverage | 1 MB | 200x | 100 | none | 3 | 47.8 | 47.5 | 2,271 |
| 05_umi_simplex | 10 MB | 100x | none | UMI simplex | 1 | 376.1 | 376.1 | 35,552 ⚠️ |
| 06_umi_duplex | 10 MB | 100x | none | UMI duplex | 1 | 373.6 | 373.6 | 35,553 ⚠️ |
| 07_cfdna | 1 MB | 200x | 200 | cfDNA, 2% purity | 3 | 47.5 | 47.1 | 2,382 |
| 08_ffpe_artifacts | 10 MB | 30x | 500 | FFPE+oxoG | 3 | 87.5 | 86.1 | 4,547 |
| 09_all_features | 10 MB | 100x | 500 | all | 1 | 471.0 | 471.0 | 39,521 ⚠️ |
| 10_large_genome | 50 MB | 30x | 1,000 | none | 1 | 374.5 | 374.5 | 16,846 ⚠️ |
| 11_panel_simulation | 1 MB | 200x | 50 | UMI simplex | 3 | 76.8 | 74.7 | 7,161 |
| 12_ultra_deep | 1 MB | 500x | 20 | UMI duplex | 1 | 190.5 | 190.5 | 17,901 ⚠️ |

⚠️ Peak RAM exceeds 12 GiB physical memory; system used swap. Times are degraded.

### Key Observations

- **Baseline throughput**: ~13,400 read pairs/sec at 30x on a 10 MB reference using 12
  cores.  This corresponds to ~1 GB compressed FASTQ (R1+R2) in 75 seconds.
- **Variant spike-in overhead**: Adding 500 random mutations adds < 1 second (+0.6%)
  to a 30x 10 MB run — negligible.
- **FFPE artifact overhead**: FFPE + oxoG damage adds ~11.6 seconds (+15%) at 30x on
  10 MB, producing 1.26M read pairs with damage patterns.
- **UMI memory regression**: UMI mode at 10 MB × 100x (10.1M read pairs) requires
  ~35 GB RAM — far exceeding the test system's 12 GB.  This indicates that VarForge
  currently buffers all UMI families in memory before writing output.  This is a
  known scalability issue to be addressed in a future release.
- **cfDNA and panel targets are efficient**: Using a 1 MB reference at 200x produces
  ~673K read pairs in ~48 seconds with only 2.4 GB RAM.
- **50 MB genome**: 5M read pairs at 30x in 374 seconds (1.36 hours estimated for a
  3 Gb human genome at 30x with proportional scaling).

---

## Thread Scaling

Configuration: `02_with_variants` (10 MB reference, 30x, 500 variants).

| Threads | Avg Wall (s) | Speedup vs 1T | Efficiency |
|---------|-------------|---------------|------------|
| 1       | 86.0        | 1.00×         | 100%       |
| 2       | 75.8        | 1.13×         | 57%        |
| 4       | 74.0        | 1.16×         | 29%        |
| 6       | 74.3        | 1.16×         | 19%        |
| 8       | 75.1        | 1.14×         | 14%        |
| 12      | 75.7        | 1.14×         | 10%        |

### Interpretation

Thread scaling is nearly flat beyond 2 threads for the 10 MB / 30x configuration.
The likely explanation is that the workload at this scale is **I/O-bound rather than
compute-bound**: writing ~280 MB of compressed FASTQ (R1+R2) saturates the storage
pipeline regardless of how many cores are computing sequences.

The single-thread run is ~14% slower than multi-thread, suggesting a modest parallelism
benefit for the compute portion (fragment generation, quality scoring, variant application).

For larger workloads (high-coverage, larger genomes) where the compute-to-I/O ratio is
higher, thread scaling is expected to improve.  The UMI/all-features runs at 10 MB × 100x
produce ~10× more data per run; their single-pass throughput would benefit more from
additional cores.

---

## Coverage Scaling

Configuration: 1 MB reference, 100 random mutations, varying coverage.
Three iterations each.

| Coverage | Avg Read Pairs | Avg Wall (s) | Avg Peak RAM (MB) | Reads/sec |
|----------|---------------|-------------|-------------------|-----------|
| 1x       | ~6,735        | 0.3         | 37                | ~22,450   |
| 5x       | ~33,675       | 1.3         | 81                | ~25,900   |
| 10x      | ~67,350       | 2.5         | 137               | ~26,940   |
| 30x      | ~202,050      | 7.2         | 362               | ~28,063   |
| 50x      | ~336,750      | 11.9        | 586               | ~28,382   |
| 100x     | ~673,335      | 23.8        | 1,146             | ~28,291   |
| 200x     | ~1,346,670    | 47.7        | 2,271             | ~28,231   |

### Interpretation

- **Linear scaling**: Wall time and memory both scale linearly with coverage depth
  (R² > 0.999 for both), confirming O(N) complexity in the number of read pairs.
- **Throughput stabilises** at ~28,000 read pairs/sec above 10x (startup/warmup
  overhead dominates at very low coverage).
- **Memory per read pair**: approximately 1.69 KB/read-pair, consistent across depths.
  This will need to decrease for the system to handle high-coverage large-genome runs
  within typical HPC memory limits.

---

## Feature Overhead

Configuration: 1 MB reference, 30x (~100K read pairs).
Three iterations each.  Values show mean wall time and overhead relative to the no-feature
baseline.

| Feature config | Avg Wall (s) | Overhead | Notes |
|----------------|-------------|----------|-------|
| Base (no features) | 7.2 | — | |
| +Variants (500) | 7.5 | +0.2s (+3%) | Variant lookup + application |
| +UMI simplex | 11.9 | +4.6s (+64%) | PCR family simulation |
| +FFPE/oxoG artifacts | 8.4 | +1.2s (+16%) | Per-base damage pass |
| +GC bias | 7.4 | +0.2s (+2%) | Lightweight per-window sampling |
| All combined | 16.3 | +9.1s (+125%) | Slightly super-additive |

### Interpretation

- **Variants** are essentially free at 500 mutations in a 1 MB region.
- **GC bias** is negligible (window-based lookup table).
- **FFPE/oxoG artifacts** add ~16% — a per-base probabilistic pass over all sequence data.
- **UMI simulation** is the most expensive individual feature at +64%.  PCR family
  simulation requires tracking molecule genealogies and assigning barcodes, adding both
  CPU and memory overhead.
- **Combined overhead is super-additive** (+125% vs sum-of-parts ~85%) due to interaction
  effects: UMI families must carry variant and artifact state through each PCR cycle.

---

## Throughput Summary

| Scenario | Read pairs / sec | Context |
|----------|-----------------|---------|
| Baseline, no features, 10 MB, 30x | 13,400 | I/O-bound (FASTQ write) |
| With 500 variants, 10 MB, 30x | 13,300 | +0% vs baseline |
| cfDNA mode, 1 MB, 200x | 14,100 | Slightly faster (shorter fragments) |
| FFPE artifacts, 10 MB, 30x | 14,400 | More reads per fragment (shorter inserts) |
| Coverage scaling plateau (1 MB) | ~28,000 | Compute-bound on small reference |
| 50 MB genome, 30x | 13,400 | Consistent with 10 MB at same depth |
| UMI simplex, 10 MB, 100x | 26,900 | Reported by app (excludes swap overhead) |

Note: reads/sec reported by VarForge internally (excluding FASTQ write time) differs from
wall-clock throughput when I/O is the bottleneck.

---

## Projected Scaling to Full Human Genome (3 Gb)

Extrapolating from 10 MB × 30x at ~13,400 read pairs/sec:

| Depth | Estimated read pairs | Estimated wall time |
|-------|---------------------|---------------------|
| 30x   | 600M                | ~12.4 hours         |
| 100x  | 2B                  | ~41.4 hours (swap risk) |

These projections assume no memory pressure.  With sufficient RAM (≥64 GB for 30x whole
genome), wall time is expected to scale linearly.  I/O throughput to fast SSD is the
primary bottleneck.

A future streaming output mode (writing reads as they are generated rather than buffering
by region) would dramatically reduce peak memory and likely improve throughput by
eliminating swap.
