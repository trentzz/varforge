# Downstream Variant Caller Benchmarking

This directory contains the design, configuration, and scripts for evaluating
how well standard variant callers (Mutect2, VarDict) perform on VarForge output.

The benchmark addresses the primary outstanding validation gap: does VarForge
produce reads that behave like real tumour sequencing data when presented to
a downstream caller?

---

## Quick start

```bash
# 1. Install tools (see setup.md)
# 2. Build the VarForge binary
cargo build --release

# 3. Run the full pipeline
cd benchmarking/downstreamvc
./scripts/run_benchmark.sh \
    --varforge ../../target/release/varforge \
    --ref ../../chr22.fa \
    --outdir results/run_$(date +%Y%m%d) \
    --threads 12
```

Results are written to `results/run_<date>/summary.tsv` and a per-scenario
breakdown in `results/run_<date>/eval/`.

---

## What is measured

For each scenario, VarForge generates reads with a known set of injected
variants (the truth VCF). The pipeline:

1. Aligns reads to the reference with bwa-mem2.
2. Calls somatic variants with Mutect2.
3. Compares called variants against the truth VCF with rtg-tools vcfeval.
4. Reports sensitivity, precision, and F1 per VAF tier.

This answers: "Given that VarForge injected a variant at VAF X, does Mutect2
detect it?"

---

## Scenarios

| # | Name | Coverage | VAF range | Features |
|---|------|----------|-----------|----------|
| 1 | SNV baseline | 30x | 0.05--0.50 | none |
| 2 | Low-VAF ctDNA | 30x | 0.01--0.10 | cfDNA frags |
| 3 | FFPE artefacts | 30x | 0.05--0.50 | FFPE/oxoG |
| 4 | Panel high-depth | 200x, 1 Mbp BED | 0.01--0.20 | UMI simplex |

See `design.md` for the full rationale.

---

## Directory layout

```
downstreamvc/
├── README.md          this file
├── design.md          experimental design, metrics, expected outcomes
├── setup.md           tool installation instructions
├── configs/           VarForge YAML configs (one per scenario)
├── scripts/
│   ├── run_benchmark.sh   end-to-end pipeline
│   └── evaluate.py        VCF comparison and metric computation
└── results/           output directory (gitignored except .gitkeep)
```
