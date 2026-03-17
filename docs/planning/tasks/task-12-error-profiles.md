# Task 12: Learned Error Profiles

## Phase
4 - Advanced Features

## Dependencies
Task 05 (Read Generation Engine)

## Objective
Support loading empirical error profiles (quality score distributions, substitution matrices) from a JSON/YAML profile file, enabling more realistic simulation that matches a specific sequencer or library prep.

## Context
The current `ParametricQualityModel` uses a simple decay function. Real sequencers have position-specific quality distributions and context-dependent error patterns (e.g., after GGG runs on NovaSeq). This task adds support for loading learned profiles.

## Key Files
- **Create**: `src/core/error_profile.rs`
- **Modify**: `src/core/mod.rs`, `src/core/quality.rs` (use loaded profiles), `src/core/engine.rs`
- **Existing context**: `src/core/quality.rs` (ParametricQualityModel)

## Requirements

### Profile Format (JSON)
```json
{
  "platform": "NovaSeq 6000",
  "read_length": 150,
  "quality_distribution": {
    "read1": [[37, 0.8], [36, 0.1], [35, 0.05], ...],  // per-position quality distributions
    "read2": [[36, 0.7], [35, 0.15], ...]
  },
  "substitution_matrix": {
    "A>C": 0.001, "A>G": 0.003, "A>T": 0.0005,
    "C>A": 0.001, "C>G": 0.0005, "C>T": 0.004,
    ...
  },
  "context_effects": {
    "GGG": { "quality_penalty": 5 }  // NovaSeq poly-G quality drop
  }
}
```

### EmpiricalQualityModel
```rust
pub struct EmpiricalQualityModel {
    quality_distributions: Vec<Vec<(u8, f64)>>,  // per-position weighted quality values
    substitution_matrix: HashMap<(u8, u8), f64>,  // (from, to) -> probability
    context_effects: HashMap<Vec<u8>, ContextEffect>,
}
```

### Integration
- If `quality.profile_path` is set in config, load empirical model instead of parametric
- Falls back to parametric model if no profile provided
- Both models implement the same trait

### Trait Abstraction
```rust
pub trait QualityModel {
    fn generate_qualities(&self, length: usize, rng: &mut impl Rng) -> Vec<u8>;
    fn inject_errors(&self, sequence: &mut [u8], qualities: &[u8], rng: &mut impl Rng);
}
```

## Tests

### Unit Tests
1. `test_load_profile` - Parse a valid profile JSON
2. `test_invalid_profile` - Error on malformed profile
3. `test_empirical_quality_distribution` - Generated qualities match loaded distribution (chi-square)
4. `test_substitution_matrix` - Error types match loaded matrix proportions
5. `test_context_effects` - Quality drops in poly-G context
6. `test_trait_compatibility` - Both parametric and empirical implement QualityModel
7. `test_fallback_to_parametric` - No profile path → parametric model used

## Acceptance Criteria
- [ ] Empirical profiles loaded from JSON
- [ ] Quality and error models follow loaded distributions
- [ ] Context-dependent effects applied
- [ ] QualityModel trait abstracts both models
- [ ] All 7 tests pass
- [ ] `cargo clippy` clean
