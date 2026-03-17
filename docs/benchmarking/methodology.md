# Benchmark Methodology

## Experiments

Four complementary experiments characterise VarForge performance:

1. **Main scenario benchmarks** — seven configurations within 12 GiB memory
2. **Coverage scaling** — wall time and memory vs. depth (1×–200×)
3. **Feature overhead** — incremental cost of each feature
4. **Thread scaling** — speedup vs. thread count

## Synthetic References

Generated with ~40% GC content via `benchmark_output/gen_refs.py`:

| File | Size | Chromosomes | Purpose |
|------|------|-------------|---------|
| `ref_1mb.fa` | ~1 MB | 6 | Coverage/feature scaling |
| `ref_10mb.fa` | ~10 MB | 8 | Main scenarios |
| `ref_50mb.fa` | ~50 MB | 12 | Large genome projection |

## Main Scenarios

| Config | Reference | Coverage | Variants | Features | Iterations |
|--------|-----------|----------|----------|----------|------------|
| Baseline | 10 MB | 30× | none | none | 3 |
| +Variants | 10 MB | 30× | 500 | none | 3 |
| High coverage | 10 MB | 100× | 500 | none | 1 |
| Very high cov | 1 MB | 200× | 100 | none | 3 |
| cfDNA | 1 MB | 200× | 200 | cfDNA, 2% purity | 3 |
| FFPE artifacts | 10 MB | 30× | 500 | FFPE + oxoG | 3 |
| Panel + UMI | 1 MB | 200× | 50 | UMI simplex | 3 |

## Measurement

- Wall time and peak RSS via `/usr/bin/time -v`
- 3 iterations for configs under ~2 min; 1 for larger
- Seeds differ between iterations (seed = iteration number)
- Output files deleted after timing capture

## Figure Generation

Benchmark figures are generated from `benchmark_output/results.json` using `docs/paper/generate_figures.py` (matplotlib). Figures are saved as PDF for LaTeX inclusion and PNG for documentation.
