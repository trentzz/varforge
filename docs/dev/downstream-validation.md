# Downstream Validation Guide

This guide covers how to set up and run the downstream variant caller validation
for VarForge. The validation answers the key question: does VarForge produce reads
that behave like real tumour sequencing data when presented to a somatic variant caller?

See `benchmarking/downstreamvc/design.md` for the full experimental design and
expected outcomes.

---

## What it validates

The benchmark injects known variants into synthetic reads, aligns with bwa-mem2,
calls somatic variants with Mutect2, then compares calls against the truth VCF
using rtg-tools vcfeval. Sensitivity, precision, and F1 are reported per VAF tier
and variant type.

Four scenarios are defined:

| # | Name | Coverage | VAF range | Features |
|---|------|----------|-----------|----------|
| 1 | SNV baseline | 30x chr22 | 0.05--0.50 | none |
| 2 | Low-VAF ctDNA | 30x chr22 | 0.01--0.10 | cfDNA fragmentation |
| 3 | FFPE artefacts | 30x chr22 | 0.05--0.50 | FFPE/oxoG damage |
| 4 | Panel + UMI | 200x, 1 Mbp BED | 0.01--0.20 | UMI simplex |

Scenario 1 is the primary sanity check. Scenarios 2--4 test specific capabilities.

A separate fgbio UMI validation (T142) checks that VarForge duplex UMI output
produces correct family size distributions after `GroupReadsByUmi`.

---

## Prerequisites

**Software**: conda or mamba (miniforge recommended). All other tools are installed
into the conda environment by the setup script.

**Disk**: ~15 GB for all four scenarios. ~4 GB for scenario 1 alone.

**Reference genome**: hg38 chr22. The setup script downloads this automatically.

---

## Installation

Run the setup script from the repository root. It creates the `varforge-vc` conda
environment and downloads the chr22 reference.

```bash
bash scripts/setup_downstream.sh
```

The script installs: bwa-mem2, samtools (>=1.17), gatk4 (>=4.4), rtg-tools, fgbio,
and Python (>=3.10). The reference is placed at `benchmarking/refs/chr22.fa`.

After the script completes, activate the environment and index the reference:

```bash
conda activate varforge-vc

samtools faidx benchmarking/refs/chr22.fa
bwa-mem2.avx2 index benchmarking/refs/chr22.fa
gatk CreateSequenceDictionary -R benchmarking/refs/chr22.fa
rtg format -o benchmarking/refs/chr22.sdf benchmarking/refs/chr22.fa
```

Indexing takes roughly 10--20 minutes (bwa-mem2 index dominates). Run it once
and the indices are reused for all subsequent benchmark runs.

---

## Running the benchmark

Build VarForge first:

```bash
cargo build --release
```

Then run the benchmark script. The script handles generation, alignment, variant
calling, and evaluation for all scenarios.

```bash
cd benchmarking/downstreamvc

./scripts/run_benchmark.sh \
    --varforge ../../target/release/varforge \
    --ref ../../benchmarking/refs/chr22.fa \
    --sdf ../../benchmarking/refs/chr22.sdf \
    --outdir results/run_$(date +%Y%m%d) \
    --threads 12
```

To run a single scenario:

```bash
./scripts/run_benchmark.sh \
    --varforge ../../target/release/varforge \
    --ref ../../benchmarking/refs/chr22.fa \
    --sdf ../../benchmarking/refs/chr22.sdf \
    --scenarios 1 \
    --threads 12
```

To skip re-generating reads (reuse existing FASTQs from a previous run):

```bash
./scripts/run_benchmark.sh ... --skip-generate
```

---

## Running the fgbio UMI validation (T142)

This validation checks that UMI family sizes from VarForge duplex output are
consistent with the configured `family_size_mean`.

Generate a duplex UMI BAM using the twist-umi-duplex preset:

```bash
./target/release/varforge simulate \
    --config benchmarking/downstreamvc/configs/04_panel_umi.yaml \
    --threads 8
```

Align the reads:

```bash
bwa-mem2.avx2 mem -t 8 \
    -R "@RG\tID:1\tSM:SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    benchmarking/refs/chr22.fa R1.fastq.gz R2.fastq.gz \
| samtools sort -n \
| samtools fixmate -m - - \
| samtools sort \
| samtools markdup - aligned.bam

samtools index aligned.bam
```

Run fgbio family grouping:

```bash
fgbio GroupReadsByUmi \
    -i aligned.bam \
    -o grouped.bam \
    -f family-sizes.txt \
    --strategy adjacency
```

The histogram in `family-sizes.txt` reports read counts per family. Compare the
mean family size against the `family_size_mean` set in the VarForge config. A
good match (within ~20%) confirms UMI tag generation is working correctly.

---

## Outputs

Each run writes to `benchmarking/downstreamvc/results/run_<timestamp>/`:

```
run_20260328/
├── versions.txt           tool versions and run parameters
├── summary.tsv            aggregated results across all scenarios
└── s1/                    per-scenario directory (s1 through s4)
    ├── config.yaml        patched VarForge config
    ├── generate/          FASTQ output and truth VCF from VarForge
    ├── align/             sorted, markdup BAM
    ├── call/              Mutect2 raw and filtered VCF
    └── eval/
        ├── rtg/           rtg vcfeval output
        └── summary.tsv    per-scenario metrics
```

The `summary.tsv` columns are:
`scenario, tier, variant_type, tp, fp, fn, sensitivity, precision, f1`

VAF tiers: `high` (0.20--0.50), `medium` (0.10--0.20), `low` (0.05--0.10),
`very_low` (0.01--0.05).

---

## Interpreting results

Key checks for scenario 1 (SNV baseline):

- Sensitivity >90% for VAF >0.20 confirms high-VAF variants are detectable.
- Sensitivity >75% for VAF 0.05--0.10 is the minimum acceptable floor.
- Precision >95% confirms no systematic false positives.

For scenario 2 (ctDNA): expect sensitivity to drop below 0.05 VAF at 30x. This
is expected and quantifies the coverage floor, not a VarForge defect.

For scenario 3 (FFPE): expect more false positives than scenario 1. Mutect2
artefact filtering should reduce but not eliminate them.

Log files for each step are in the per-scenario directories. Start with
`generate.log`, `bwa.log`, `mutect2.log`, and `rtg.log` when debugging failures.

---

## Blocked status

T135 and T142 are currently blocked because this machine lacks the required tools
and reference. Run `scripts/setup_downstream.sh` to unblock both tasks.
