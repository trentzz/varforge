# Task 20: Performance Benchmarking & Optimization

## Phase
6 - Quality & Polish

## Dependencies
Task 07 (Parallel Region Processing)

## Objective
Benchmark VarForge's performance against targets and existing tools, identify bottlenecks, and optimize for the throughput and memory targets defined in requirements.

## Context
Performance targets from `docs/features/initial-requirements/performance.md`:
- 30x WGS in <60 minutes (16 threads)
- Targeted panel in <5 minutes
- Memory: <8GB for WGS, <2GB for panels
- Near-linear scaling up to 16 threads

## Key Files
- **Create**: `benches/` directory with criterion benchmarks
- **Create**: `scripts/benchmark.sh` for end-to-end timing
- **Modify**: Any hot-path code identified as bottleneck

## Requirements

### Micro-benchmarks (criterion)
1. Fragment sampling throughput (fragments/sec)
2. Quality score generation throughput
3. Variant spike-in throughput (variants/sec)
4. FASTQ writing throughput (reads/sec, MB/sec)
5. BAM writing throughput
6. UMI generation + PCR family expansion throughput

### Macro-benchmarks
1. Small panel: 50 targets, 500x, UMI → time and memory
2. WES: ~60Mb targets, 100x → time and memory
3. WGS estimate: chr22 at 30x → extrapolate to full genome

### Thread Scaling
- Benchmark with 1, 2, 4, 8, 16 threads
- Plot speedup curve
- Identify serialization bottlenecks

### Optimization Targets
- **FASTQ gzip**: Ensure gzp multi-threaded compression is engaged
- **Memory**: Streaming output (don't hold all reads in memory)
- **Allocation**: Minimize Vec allocations in hot loops (pre-allocate)
- **SIMD**: Consider for quality score generation if beneficial

### Comparison
Run ART and DWGSIM on same region/coverage for baseline comparison. VarForge should be within 2x of ART for pure read generation (with more features).

## Tests
Benchmarks are not traditional pass/fail tests but should have performance regression guards:
1. `bench_fragment_sampling` - >1M fragments/sec
2. `bench_quality_generation` - >500K reads/sec
3. `bench_fastq_writing` - >100MB/sec compressed
4. `bench_full_panel` - <5 minutes for standard panel

## Acceptance Criteria
- [ ] Criterion benchmarks for all hot paths
- [ ] End-to-end benchmark script
- [ ] Thread scaling measured and documented
- [ ] Performance within 2x of pure C read simulators
- [ ] Memory usage within targets
