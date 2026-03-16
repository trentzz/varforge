# Feature: Streaming Output Pipeline

## Status
**Proposed** — Addresses the #1 scalability bottleneck identified in benchmarking.

## Problem Statement

VarForge currently buffers all generated reads in memory before writing to disk. The pipeline is:

```
Region 1 → generate all reads → Vec<ReadPair>  ─┐
Region 2 → generate all reads → Vec<ReadPair>  ─┤  collect ALL  →  write ALL to FASTQ/BAM
Region 3 → generate all reads → Vec<ReadPair>  ─┤
...                                              ─┘
```

This architecture causes **memory usage proportional to the total number of read pairs**:

| Scenario | Read Pairs | Memory (measured) | Memory per pair |
|----------|-----------|-------------------|-----------------|
| 10MB, 30x, baseline | 1.0M | 3.4 GiB | 1.69 KB |
| 10MB, 100x, UMI | 10.1M | 35.6 GiB | 3.5 KB |
| 10MB, 100x, all features | 10.1M | 39.5 GiB | 3.9 KB |

UMI mode is particularly expensive because each original molecule spawns a PCR family of 3-10 copies, multiplying the read count 3-10x.

**For a full human genome at 30x**: ~600M read pairs × 1.69 KB = **~1 TB of buffer memory**, which is obviously infeasible.

The root cause is in two places:

### 1. `SimulationEngine::generate_region()` (src/core/engine.rs)

Returns `RegionOutput { read_pairs: Vec<ReadPair>, ... }` — all reads for a region are materialized in a vector.

### 2. `run_single_sample()` (src/cli/simulate.rs)

Collects all `RegionOutput` results from parallel region processing into a `Vec<Result<(Region, RegionOutput)>>`, then iterates to write sequentially:

```rust
let region_results: Vec<Result<(Region, RegionOutput)>> = regions
    .par_iter()
    .enumerate()
    .map(|(idx, region)| { /* generate reads */ })
    .collect();  // <-- ALL reads in memory at once

for result in region_results {
    // write to FASTQ/BAM
}
```

## Solution: Channel-Based Streaming Pipeline

Replace the collect-then-write pattern with a **bounded channel** where worker threads send reads directly to a writer thread. Reads are written to disk as they are generated, keeping only a small buffer in memory.

### Architecture

```
                                     ┌─────────────┐
Region 1 ──→ Engine Worker ──→──┐    │             │
Region 2 ──→ Engine Worker ──→──┤    │  Bounded    │    ┌──────────────┐
Region 3 ──→ Engine Worker ──→──┼──→ │  Channel    │──→ │ Writer Thread│──→ FASTQ.gz
Region 4 ──→ Engine Worker ──→──┤    │  (N slots)  │    │              │──→ BAM
   ...                    ──→──┘    │             │    │              │──→ Truth VCF
                                     └─────────────┘    └──────────────┘
        rayon thread pool                                  single thread
```

### Key Design Decisions

#### Bounded Channel (crossbeam or std::sync::mpsc)

```rust
use crossbeam_channel::{bounded, Sender, Receiver};

// Buffer at most 64 regions worth of reads at a time
let (tx, rx): (Sender<RegionBatch>, Receiver<RegionBatch>) = bounded(64);
```

The bounded channel provides **backpressure**: if the writer thread falls behind, worker threads block on send, preventing unbounded memory growth. The buffer size (e.g., 64 region batches) controls the memory ceiling.

#### RegionBatch — the unit of transfer

```rust
pub struct RegionBatch {
    pub region: Region,
    pub read_pairs: Vec<ReadPair>,
    pub applied_variants: Vec<AppliedVariant>,
}
```

Each region's reads are still generated as a batch (the per-region Vec is bounded by per-region coverage, typically a few thousand pairs for a 1 Mbp chunk at 30x). The key change is that batches are **sent to the writer immediately** rather than accumulated.

#### Per-region memory bound

For a 1 Mbp chunk at 30x with 150bp reads:
- read_pairs = 30 × 1,000,000 / (2 × 150) = ~100,000 pairs
- Memory per batch ≈ 100K × 1.69 KB ≈ 169 MB

With 64 batches buffered: 64 × 169 MB ≈ **10.5 GiB max** — acceptable for most systems. And this is for 30x WGS; targeted panels would use far less.

The buffer size can be tuned based on available memory:
```yaml
performance:
  output_buffer_regions: 64  # default, tune down for low-memory systems
```

#### Writer Thread

A dedicated writer thread consumes `RegionBatch` items from the channel and writes them to FASTQ/BAM/VCF:

```rust
fn writer_thread(
    rx: Receiver<RegionBatch>,
    mut fastq: FastqWriter,
    mut bam: Option<BamWriter>,
    mut truth_vcf: Option<TruthVcfWriter>,
) -> Result<WriterStats> {
    let mut stats = WriterStats::default();
    for batch in rx {
        // Write each read pair to FASTQ
        for (i, pair) in batch.read_pairs.iter().enumerate() {
            fastq.write_pair(pair, &pair.name)?;
        }
        // Write to BAM if enabled
        if let Some(ref mut bam_w) = bam {
            for pair in &batch.read_pairs {
                bam_w.write_pair(/* ... */)?;
            }
        }
        // Record applied variants for truth VCF
        for av in &batch.applied_variants {
            if let Some(ref mut vcf_w) = truth_vcf {
                vcf_w.write_variant(/* ... */)?;
            }
        }
        stats.total_read_pairs += batch.read_pairs.len();
        stats.regions_completed += 1;
    }
    // Finalize all writers
    fastq.finish()?;
    if let Some(bam_w) = bam { bam_w.finish()?; }
    if let Some(vcf_w) = truth_vcf { vcf_w.finish()?; }
    Ok(stats)
}
```

#### Modified Parallel Region Loop

```rust
let (tx, rx) = bounded(buffer_size);

// Spawn writer thread
let writer_handle = std::thread::spawn(move || {
    writer_thread(rx, fastq_writer, bam_writer, truth_vcf_writer)
});

// Worker threads send batches as they complete
regions.par_iter().enumerate().try_for_each(|(idx, region)| {
    let mut engine = SimulationEngine::new_with_shared_reference(/* ... */);
    let output = engine.generate_region(region, &variants_for_region)?;
    tx.send(RegionBatch {
        region: region.clone(),
        read_pairs: output.read_pairs,
        applied_variants: output.applied_variants,
    }).map_err(|e| anyhow::anyhow!("writer channel closed: {}", e))?;
    pb.inc(1);
    Ok::<(), anyhow::Error>(())
})?;

// Signal completion by dropping the sender
drop(tx);

// Wait for writer to finish
let stats = writer_handle.join()
    .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;
```

### UMI-Specific Optimization

UMI mode is the worst offender because PCR family expansion multiplies reads 3-10x. Two additional optimizations:

#### 1. Deferred PCR Expansion

Currently, `generate_region()` generates PCR copies for every UMI family and stores them all in the `Vec<ReadPair>`. Instead, defer expansion to the writer thread:

```rust
pub struct DeferredUmiFamily {
    pub original: ReadPair,
    pub umi: String,
    pub family_size: usize,
    pub pcr_error_rate: f64,
}

pub struct RegionBatch {
    pub region: Region,
    pub read_pairs: Vec<ReadPair>,           // non-UMI reads
    pub umi_families: Vec<DeferredUmiFamily>, // UMI: expand during write
    pub applied_variants: Vec<AppliedVariant>,
}
```

The writer thread expands PCR families one at a time, writes each copy immediately, and discards it. This reduces peak memory by the PCR amplification factor (typically 3-10x).

**Expected memory savings for UMI mode:**
- Current: 10.1M read pairs × 3.5 KB = 35.6 GiB
- With deferred expansion: ~1.0M originals × 3.5 KB + ~170 MB channel buffer = **~3.7 GiB**
- **~10x memory reduction**

#### 2. Intra-Region Streaming for Very High Coverage

For ultra-deep sequencing (>500x), even a single region's reads may be large. Split generation within a region into sub-batches:

```rust
pub fn generate_region_streaming(
    &mut self,
    region: &Region,
    variants: &[Variant],
    batch_size: usize,  // e.g., 10,000 pairs per sub-batch
    tx: &Sender<RegionBatch>,
) -> Result<Vec<AppliedVariant>> {
    let n_pairs = self.calculate_pairs(region);
    let mut applied = Vec::new();

    for chunk_start in (0..n_pairs).step_by(batch_size) {
        let chunk_end = (chunk_start + batch_size).min(n_pairs);
        let pairs = self.generate_pairs(region, variants, chunk_start, chunk_end)?;
        // ... track applied variants ...
        tx.send(RegionBatch { read_pairs: pairs, /* ... */ })?;
    }
    Ok(applied)
}
```

### Memory Budget Analysis

After streaming + deferred UMI expansion:

| Scenario | Current Memory | Streaming Memory | Reduction |
|----------|---------------|-----------------|-----------|
| 10MB, 30x, baseline | 3.4 GiB | ~0.5 GiB | 7x |
| 10MB, 100x, UMI | 35.6 GiB | ~3.7 GiB | 10x |
| 10MB, 100x, all | 39.5 GiB | ~4.2 GiB | 9x |
| 3Gb, 30x, baseline | ~1 TB (est.) | ~10.5 GiB | ~100x |
| 3Gb, 30x, UMI | infeasible | ~15 GiB | ∞ |

### BAM Sorting Consideration

Streaming output means reads arrive in region order (which is coordinate-sorted by region start). Within a region, read positions are random. Options:

1. **Unsorted BAM** (default for FASTQ workflows — users will re-align anyway)
2. **Region-sorted BAM** — sort within each region batch before writing. Gives approximately-sorted output, suitable for many tools.
3. **External sort** — write unsorted BAM, then use a post-processing sort step. Most robust but adds I/O overhead.

Recommendation: Default to region-sorted (option 2), which is sufficient for most tools and has no extra I/O cost.

### Implementation Plan

#### Phase 1: Core Streaming (highest impact)
1. Add `crossbeam-channel` dependency
2. Refactor `run_single_sample()` to use channel-based pipeline
3. Spawn dedicated writer thread
4. Keep `generate_region()` API unchanged (returns Vec)
5. Add `performance.output_buffer_regions` config option
6. **Expected outcome**: 7-100x memory reduction for non-UMI workloads

#### Phase 2: Deferred UMI Expansion
1. Add `DeferredUmiFamily` struct
2. Split UMI family generation from PCR expansion in engine
3. Writer thread performs PCR expansion inline
4. **Expected outcome**: 10x additional memory reduction for UMI workloads

#### Phase 3: Intra-Region Streaming
1. Add `generate_region_streaming()` method to engine
2. Sub-batch generation within very large regions
3. **Expected outcome**: Enables ultra-deep (>2000x) on targeted panels without memory issues

#### Phase 4: Sorted BAM Streaming
1. Implement per-batch coordinate sorting
2. Write sorted batches sequentially for approximately-sorted BAM output
3. Optional external sort post-processing step

### Dependencies

- `crossbeam-channel` (mature, widely used, lock-free bounded channels)
- No other new dependencies

### Testing Strategy

1. **Correctness**: Same seed must produce byte-identical FASTQ output with and without streaming (just memory/ordering changes, not content)
2. **Memory regression test**: Monitor peak RSS during simulation, assert it stays below threshold
3. **Throughput test**: Streaming should not reduce throughput (may improve it due to overlapping compute+I/O)
4. **Backpressure test**: Slow writer (throttled) should not cause OOM
5. **UMI deferred expansion test**: PCR copies match direct expansion (same seed, same output)

### Metrics to Track

After implementation, re-run the benchmark suite and compare:
- Peak RSS at each scenario
- Wall time (should be similar or slightly better)
- Throughput (reads/sec)
- Thread scaling (should improve — writer thread decouples I/O from compute)

### Expected Impact on Benchmark Results

| Metric | Before | After (projected) |
|--------|--------|-------------------|
| 10MB 30x peak RAM | 3.4 GiB | ~0.5 GiB |
| 10MB 100x UMI peak RAM | 35.6 GiB | ~3.7 GiB |
| 3Gb 30x peak RAM | ~1 TB | ~10.5 GiB |
| 3Gb 30x feasibility | Infeasible (12 GiB system) | **Feasible** |
| Thread scaling | Flat (I/O-bound) | **Improved** (I/O overlaps compute) |
| 10MB 30x wall time | 75s | ~70s (slight improvement from overlap) |

The most important outcome: **full human genome simulation becomes feasible on a standard workstation** (16-64 GiB RAM), which is currently impossible.
