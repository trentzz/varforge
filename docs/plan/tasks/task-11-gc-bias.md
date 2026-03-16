# Task 11: GC Bias Modeling

## Phase
4 - Advanced Features

## Dependencies
Task 05 (Read Generation Engine)

## Objective
Implement GC-content-based coverage bias so that regions with extreme GC content have reduced coverage, matching real Illumina sequencing behavior.

## Context
Real sequencing has systematic GC bias: AT-rich and GC-rich regions have lower coverage than regions near 50% GC. This affects variant calling sensitivity in those regions. ReSeq is the gold standard for realistic GC bias modeling.

## Key Files
- **Create**: `src/core/gc_bias.rs`
- **Modify**: `src/core/mod.rs`, `src/core/engine.rs` (apply GC bias to coverage)

## Requirements

### GC Bias Model
Default model: bell-shaped curve centered at ~50% GC content.

```rust
pub struct GcBiasModel {
    // Parameters for the bias curve
}

impl GcBiasModel {
    /// Returns coverage multiplier (0.0 to 1.0) for a given GC fraction
    pub fn coverage_multiplier(&self, gc_fraction: f64) -> f64;

    /// Calculate GC fraction of a sequence
    pub fn gc_fraction(sequence: &[u8]) -> f64;
}
```

### Default Bias Curve
Modeled as a quadratic or Gaussian centered at 0.5:
- GC = 0.5 → multiplier = 1.0
- GC = 0.2 or 0.8 → multiplier ~0.7
- GC = 0.1 or 0.9 → multiplier ~0.3

### Application
For each fragment:
1. Calculate GC fraction of the fragment sequence
2. Get coverage multiplier from model
3. Use multiplier as probability of keeping the fragment (rejection sampling)

### Configurable
```yaml
gc_bias:
  enabled: true
  model: default  # or "flat" (no bias) or "custom"
  severity: 1.0   # multiplier on bias effect (0 = no bias, 2 = extreme bias)
```

## Tests

### Unit Tests
1. `test_gc_fraction` - Correct GC% calculation for known sequences
2. `test_multiplier_at_50_percent` - ~1.0 at 50% GC
3. `test_multiplier_at_extremes` - Reduced at 10% and 90% GC
4. `test_flat_model` - Flat model returns 1.0 everywhere
5. `test_severity_scaling` - Higher severity amplifies bias
6. `test_coverage_distribution` - After GC bias, AT-rich regions have fewer reads (statistical)
7. `test_gc_fraction_n_bases` - N bases excluded from GC calculation

## Acceptance Criteria
- [ ] GC bias model reduces coverage at extreme GC content
- [ ] Configurable severity (can be disabled)
- [ ] All 7 tests pass
- [ ] `cargo clippy` clean
