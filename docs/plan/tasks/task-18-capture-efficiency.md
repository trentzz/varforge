# Task 18: Capture Efficiency Model

## Phase
5 - Extended Simulation Modes

## Dependencies
Task 05 (Read Generation Engine)

## Objective
Model target capture efficiency for panel/WES simulation: uneven coverage across targets, off-target reads, and coverage dropoff at target boundaries.

## Context
Targeted sequencing (panels, WES) uses hybridization capture or amplicon-based enrichment. Coverage is not uniform: some targets capture better than others, coverage drops at target edges, and some reads map off-target. GENOMICON-Seq is the only tool that models probe-capture enrichment, but it's Docker-only and limited.

## Key Files
- **Create**: `src/core/capture.rs`
- **Modify**: `src/core/mod.rs`, `src/core/engine.rs` (apply capture model)

## Requirements

### Capture Model
```rust
pub struct CaptureModel {
    target_regions: Vec<Region>,
    off_target_fraction: f64,       // default: 0.2 (20% off-target)
    coverage_uniformity: f64,        // 0.0 = perfectly uniform, 1.0 = highly variable
    edge_dropoff_bases: u32,         // bases of dropoff at target edges (default: 50)
}
```

### Coverage Variation
- Each target gets a coverage multiplier sampled from LogNormal(0, uniformity)
- Target edges have exponential coverage dropoff
- Off-target regions get `off_target_fraction × mean_coverage`

### Config
```yaml
capture:
  enabled: true
  targets_bed: "panel.bed"
  off_target_fraction: 0.2
  coverage_uniformity: 0.3
  edge_dropoff_bases: 50
```

## Tests

### Unit Tests
1. `test_uniform_targets` - uniformity=0 gives equal coverage across targets
2. `test_variable_targets` - uniformity>0 gives variable coverage
3. `test_edge_dropoff` - Coverage decreases at target boundaries
4. `test_off_target_reads` - Off-target fraction approximately correct
5. `test_no_capture` - Disabled capture gives uniform coverage

## Acceptance Criteria
- [ ] Per-target coverage variation
- [ ] Edge coverage dropoff
- [ ] Off-target reads at configured fraction
- [ ] All 5 tests pass
- [ ] `cargo clippy` clean
