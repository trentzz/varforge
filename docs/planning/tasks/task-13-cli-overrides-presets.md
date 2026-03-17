# Task 13: CLI Overrides & Presets

## Phase
4 - Advanced Features

## Dependencies
Task 06 (Simulate Command Orchestrator)

## Objective
Implement CLI flag overrides for common config parameters and built-in presets for quick setup.

## Context
Users need to quickly override config values without editing YAML (e.g., `--coverage 100` to test higher coverage). Presets provide one-command setups for common scenarios.

The CLI structure in `src/cli/mod.rs` already has `--coverage`, `--seed`, `--preset` flags defined but not wired up.

## Key Files
- **Modify**: `src/cli/mod.rs` (ensure all override flags defined)
- **Modify**: `src/cli/simulate.rs` (apply overrides to loaded config)
- **Create**: `src/cli/presets.rs` (preset configurations)

## Requirements

### CLI Overrides
These flags override the corresponding YAML config value:
- `--coverage <N>` → `sample.coverage`
- `--read-length <N>` → `sample.read_length`
- `--seed <N>` → reproducibility seed
- `--output-dir <PATH>` → `output.directory`
- `--threads <N>` → thread count
- `--purity <FLOAT>` → `tumour.purity`
- `--vaf-range <MIN-MAX>` → `mutations.random.vaf_range`
- `--random-mutations <N>` → `mutations.random.count`
- `--fragment-mean <N>` → `fragment.mean`
- `--fragment-sd <N>` → `fragment.sd`

### Override Application
```rust
fn apply_overrides(config: &mut SimulationConfig, args: &SimulateArgs) {
    if let Some(cov) = args.coverage { config.sample.coverage = cov; }
    // ... etc
}
```

### Presets
```
--preset small    → 1x coverage, chr22 only, 100 random mutations, ~30 sec runtime
--preset panel    → 500x coverage, BED targets, UMI enabled, 50 mutations
--preset wgs      → 30x coverage, whole genome, 5000 mutations
--preset cfdna    → 200x coverage, cfDNA fragment model, low TF (2%), UMI enabled
--preset ffpe     → 30x coverage, FFPE damage 0.02, oxoG 0.01
--preset umi      → 1000x coverage, UMI duplex mode, panel targets
```

Each preset is a partial config that gets merged with defaults before user overrides.

### Override Precedence
`defaults < preset < YAML config < CLI flags`

## Tests

### Unit Tests
1. `test_coverage_override` - CLI flag overrides config coverage
2. `test_seed_override` - CLI seed overrides config seed
3. `test_preset_small` - Small preset has expected values
4. `test_preset_cfdna` - cfDNA preset enables cfDNA fragment model
5. `test_precedence` - CLI flag > YAML > preset > default
6. `test_no_overrides` - Config unchanged when no flags provided
7. `test_all_presets_valid` - Every preset produces a valid config

## Acceptance Criteria
- [ ] CLI flags override YAML config values
- [ ] All 6 presets defined and functional
- [ ] Correct precedence order
- [ ] All 7 tests pass
- [ ] `cargo clippy` clean
