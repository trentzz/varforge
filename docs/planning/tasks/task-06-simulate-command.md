# Task 06: Simulate Command Orchestrator

## Phase
2 - Core Pipeline

## Dependencies
Tasks 01-05 (all Foundation I/O + Read Generation Engine)

## Objective
Implement the `simulate` subcommand that orchestrates the full simulation: load config, load reference, partition regions, generate reads per region, and write all output files.

## Context
Currently `src/cli/simulate.rs` is a stub that prints "not yet implemented". This task fills it in to orchestrate the complete pipeline. It coordinates between config loading, the simulation engine, and all output writers.

## Key Files
- **Modify**: `src/cli/simulate.rs` (replace stub with full implementation)
- **Modify**: `src/cli/mod.rs` (if needed for new options)
- **Existing context**: `src/io/config.rs`, `src/core/engine.rs` (Task 05), `src/io/reference.rs` (Task 01), `src/io/fastq.rs` (Task 02), `src/io/bam.rs` (Task 03), `src/io/truth_vcf.rs` (Task 04)

## Requirements

### Orchestration Flow
```
1. Load and validate config (existing)
2. Open reference genome
3. Load variants (from config random gen or VCF)
4. Determine regions to simulate:
   - If BED file provided: use those regions
   - If whole genome: all chromosomes from reference
5. Partition regions for parallel processing
6. Initialize output writers (FASTQ, optionally BAM, truth VCF)
7. For each region partition:
   a. Run SimulationEngine::generate_region()
   b. Write ReadPairs to FASTQ
   c. If BAM enabled: write to BAM
   d. Collect applied variants
8. Write truth VCF with all applied variants
9. Write simulation manifest (summary)
10. Report statistics (reads generated, variants applied, time elapsed)
```

### Progress Reporting
- Use `indicatif` progress bar showing regions processed
- Log key milestones with `tracing` (config loaded, reference opened, simulation started, complete)
- Report final statistics: total reads, total variants, wall time, throughput (reads/sec)

### Error Handling
- Validate reference genome has all chromosomes referenced in config/BED
- Validate variants fall within reference bounds
- Report clear errors for missing files, invalid configs

### Dry Run Mode
When `--dry-run` is set:
- Load and validate everything
- Calculate and report expected read count, estimated file sizes
- Do not generate any reads or output files

## Tests

### Integration Tests (in `tests/` directory)
1. `test_simulate_minimal` - Smallest possible simulation (tiny reference, 1 region, 1x coverage) produces valid FASTQ
2. `test_simulate_with_variants` - Simulation with variants produces truth VCF
3. `test_simulate_dry_run` - Dry run produces no output files but reports stats
4. `test_simulate_missing_reference` - Clear error for missing reference
5. `test_simulate_bam_output` - BAM output flag produces valid BAM
6. `test_simulate_umi_mode` - UMI config produces reads with UMI in headers
7. `test_simulate_deterministic` - Same seed produces identical output files

### Test Fixtures
Create minimal test FASTA files in `tests/fixtures/` or generate them in test setup.

## Acceptance Criteria
- [ ] `varforge simulate --config test.yaml` produces FASTQ output
- [ ] Progress bar shows during simulation
- [ ] Truth VCF generated with all spiked variants
- [ ] Dry run mode works correctly
- [ ] All 7 tests pass
- [ ] `cargo clippy` clean
