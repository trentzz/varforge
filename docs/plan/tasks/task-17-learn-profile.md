# Task 17: Learn-Profile Subcommand

## Phase
5 - Extended Simulation Modes

## Dependencies
Task 12 (Learned Error Profiles)

## Objective
Implement the `learn-profile` subcommand that analyzes a real BAM file and extracts quality score distributions, substitution error patterns, fragment size distributions, and GC bias curves into a reusable profile file.

## Context
The CLI already defines `learn-profile` as a planned subcommand. This enables users to make VarForge's output match their specific sequencer/library prep. ReSeq pioneered this approach and showed it dramatically improves simulation realism.

## Key Files
- **Create**: `src/cli/learn_profile.rs`
- **Create**: `src/core/profile_learner.rs`
- **Modify**: `src/cli/mod.rs` (add subcommand)

## Requirements

### CLI
```
varforge learn-profile \
  --bam input.bam \
  --output profile.json \
  --sample-size 1000000    # reads to sample (default)
  --threads 4
```

### What to Learn
1. **Quality scores**: Per-position quality distribution for R1 and R2 (histograms)
2. **Substitution errors**: Base-by-base error rates using MD tags or reference comparison
3. **Fragment sizes**: Insert size distribution from proper pairs (mean, sd, full histogram)
4. **GC bias**: GC content vs coverage depth curve (bin by GC%, compute mean depth per bin)
5. **Context effects**: Tri-nucleotide context error rates (if sufficient data)

### Output Format
Same JSON format as Task 12's error profile, so it can be directly used with `quality.profile_path` in config.

### Sampling Strategy
- Randomly sample N reads from BAM (default 1M) for efficiency
- Use only properly paired, primary alignments
- Skip duplicate-flagged reads
- Report percentage of reads sampled and confidence

## Tests

### Unit Tests
1. `test_quality_extraction` - Correct quality distributions from synthetic BAM
2. `test_fragment_size_extraction` - Insert sizes match expected distribution
3. `test_gc_bias_extraction` - GC bias curve reasonable for synthetic data
4. `test_output_format` - Profile JSON matches expected schema
5. `test_sampling` - Sampling respects requested count
6. `test_round_trip` - Learned profile used for simulation produces similar quality distribution

## Acceptance Criteria
- [ ] `varforge learn-profile --bam input.bam` produces profile JSON
- [ ] Profile captures quality, errors, fragment sizes, GC bias
- [ ] Profile usable with `quality.profile_path` config
- [ ] All 6 tests pass
- [ ] `cargo clippy` clean
