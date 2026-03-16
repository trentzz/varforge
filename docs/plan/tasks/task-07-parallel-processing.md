# Task 07: Parallel Region Processing

## Phase
2 - Core Pipeline

## Dependencies
Tasks 05, 06 (Read Generation Engine, Simulate Command)

## Objective
Add rayon-based parallel processing so regions are simulated concurrently across available CPU cores, with thread-safe output writing and proper progress tracking.

## Context
The project already has `rayon` as a dependency and `coverage::partition_regions()` to split the genome into chunks. This task wires up parallel execution in the simulate command.

Key challenge: output writers (FASTQ, BAM) need thread-safe access. Strategy: generate reads in parallel, collect results, write sequentially (I/O is the bottleneck anyway and sequential writing avoids interleaving).

## Key Files
- **Modify**: `src/cli/simulate.rs` (add parallel region processing)
- **Modify**: `src/core/engine.rs` (ensure per-region RNG is deterministic regardless of thread scheduling)
- **Existing context**: `src/core/coverage.rs` (partition_regions)

## Requirements

### Parallel Strategy
```rust
use rayon::prelude::*;

let results: Vec<RegionOutput> = regions
    .par_iter()
    .map(|region| {
        // Each region gets its own RNG seeded deterministically from master seed + region index
        let mut engine = SimulationEngine::new_for_region(config, reference, seed, region_idx);
        engine.generate_region(region, &variants_for_region)
    })
    .collect::<Result<Vec<_>>>()?;
```

### Determinism
- Master seed from config → per-region seed = hash(master_seed, region_index)
- Output order is deterministic regardless of thread scheduling (sort by region before writing)
- Same seed + same thread count = identical output

### Thread Control
- Respect `--threads` CLI option (set rayon global thread pool size)
- Default to available CPU cores

### Progress
- `indicatif::ProgressBar` with region count, updated atomically from parallel workers
- Multi-progress bar: one per thread showing current region (optional, P1)

### Memory Management
- Process regions in batches if total read count would exceed memory threshold
- Configurable batch size based on available memory (default: process all if <8GB estimated)
- Flush completed batches to disk before starting next

## Tests

### Unit Tests
1. `test_deterministic_parallel` - Same seed produces same output with 1 thread and 4 threads
2. `test_region_seeds_unique` - Each region gets a unique but deterministic RNG seed
3. `test_thread_count_respected` - Thread pool size matches config

### Integration Tests
4. `test_parallel_speedup` - 4 threads is measurably faster than 1 thread for a medium workload
5. `test_output_order_stable` - Output read ordering is identical across runs regardless of thread scheduling

## Acceptance Criteria
- [ ] Simulation uses rayon for parallel region processing
- [ ] `--threads N` controls parallelism
- [ ] Output is deterministic for same seed regardless of thread count
- [ ] Progress bar updates from parallel workers
- [ ] All 5 tests pass
- [ ] `cargo clippy` clean
