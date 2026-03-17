# Task 14: Simulation Manifest Output

## Phase
4 - Advanced Features

## Dependencies
Task 06 (Simulate Command Orchestrator)

## Objective
Generate a simulation manifest (JSON) that records all parameters, generated file paths, variant counts, and statistics for reproducibility and provenance tracking.

## Context
Users need to know exactly what parameters produced a given dataset, especially when running many simulations for benchmarking sweeps. The manifest is a machine-readable record of everything.

## Key Files
- **Create**: `src/io/manifest.rs`
- **Modify**: `src/io/mod.rs`, `src/cli/simulate.rs` (write manifest after simulation)

## Requirements

### Manifest Format (JSON)
```json
{
  "varforge_version": "0.1.0",
  "timestamp": "2025-03-16T14:30:00Z",
  "seed": 42,
  "config": { /* full resolved config */ },
  "reference": {
    "path": "/path/to/hg38.fa",
    "genome_size": 3088286401,
    "chromosomes": 25
  },
  "output_files": {
    "fastq_r1": "output/sample_R1.fastq.gz",
    "fastq_r2": "output/sample_R2.fastq.gz",
    "bam": "output/sample.bam",
    "truth_vcf": "output/truth.vcf.gz"
  },
  "statistics": {
    "total_read_pairs": 450000000,
    "total_bases": 135000000000,
    "variants_spiked": 5000,
    "variants_by_type": { "SNV": 4000, "Indel": 750, "MNV": 250 },
    "mean_coverage_achieved": 30.2,
    "regions_simulated": 1500,
    "wall_time_seconds": 2400,
    "reads_per_second": 187500
  }
}
```

### Methods
1. `Manifest::new(config, reference_info) -> Self`
2. `fn add_output_file(&mut self, key: &str, path: &Path)`
3. `fn set_statistics(&mut self, stats: SimulationStatistics)`
4. `fn write(&self, path: &Path) -> Result<()>`

## Tests

### Unit Tests
1. `test_manifest_creation` - Valid manifest with all fields
2. `test_manifest_json_valid` - Output is parseable JSON
3. `test_statistics_populated` - All statistics fields present
4. `test_config_included` - Full config serialized in manifest
5. `test_file_paths_relative` - Output paths relative to output directory

## Acceptance Criteria
- [ ] JSON manifest written after every simulation
- [ ] Contains full config, file paths, and statistics
- [ ] All 5 tests pass
- [ ] `cargo clippy` clean
