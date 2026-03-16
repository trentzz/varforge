# Task 21: Documentation & Examples

## Phase
6 - Quality & Polish

## Dependencies
All previous tasks

## Objective
Write user-facing documentation with example configs for every major use case, a quickstart guide, and annotated output explanations.

## Key Files
- **Modify**: `README.md` (comprehensive user guide)
- **Create**: `examples/` directory with example configs
- **Create**: `docs/user-guide/` with detailed docs

## Requirements

### Example Configs
```
examples/
  minimal.yaml           # Simplest possible simulation
  wgs_30x.yaml           # Standard WGS with random mutations
  panel_umi.yaml         # Targeted panel with UMI
  cfdna_monitoring.yaml  # cfDNA longitudinal series
  ffpe_artifacts.yaml    # FFPE-damaged sample
  tumor_normal.yaml      # Matched tumor/normal pair
  subclonal.yaml         # Complex clonal architecture
  high_depth.yaml        # 1000x for low-VAF detection
  custom_mutations.yaml  # VCF input with specific mutations
```

### User Guide Sections
1. Installation
2. Quickstart (3-minute example)
3. Configuration reference (every YAML field documented)
4. Use case recipes (one page per use case)
5. Output format reference (FASTQ headers, BAM tags, VCF fields)
6. Presets reference
7. Performance tuning guide
8. Comparison with other tools (when to use VarForge vs alternatives)

### README
- Badge for CI status, version, license
- One-paragraph description
- Feature highlights with comparison table
- Installation (cargo install, binary download)
- Quickstart example
- Links to user guide

## Acceptance Criteria
- [ ] 9 example configs, all valid and tested
- [ ] User guide covers all major features
- [ ] README is clear and compelling
- [ ] All examples produce valid output (`cargo test` verifies)
