# Twist Biosciences Duplex Benchmarking

This document defines the reference scenario and success criteria for benchmarking
duplex consensus pipelines and low-VAF variant callers against VarForge-generated
Twist Biosciences hybrid-capture data.

## Sequencing model

Twist Biosciences Comprehensive Exome Panel duplex sequencing uses 9 bp inline
dual UMIs. Each original molecule produces two read families: AB (forward strand,
UMI `A-B`) and BA (reverse strand, UMI `B-A`). After duplex consensus calling,
each molecule contributes one consensus read. At 2 000× raw coverage with a mean
family size of 3.5×, the expected duplex consensus depth is approximately 500–600×.

Key parameters in the benchmark config:

| Parameter | Value |
|---|---|
| Raw coverage | 2 000× |
| Read length | 150 bp |
| Fragment mean | 170 bp |
| Fragment SD | 30 bp |
| UMI length | 9 bp |
| UMI mode | Inline dual, duplex |
| PCR cycles | 10 |
| Mean family size | 3.5× |
| On-target fraction | >= 95 % |
| Coverage CV target | <= 0.25 |

## Variant load

| Type | Count | VAF range |
|---|---|---|
| SNV | ~50 | 0.001–0.05 |
| Indel | ~10 | 0.001–0.05 |
| SV (HRD deletions) | 5 | 0.01–0.1 |

The SNV/indel VAF range covers the ultra-low ctDNA detection window where duplex
sequencing provides a meaningful sensitivity advantage over standard sequencing.

## Success metrics

These criteria define a passing run for benchmarking purposes. A tool that meets
all criteria is considered validated against the Twist duplex model.

| Metric | Target | Source |
|---|---|---|
| Duplex conversion rate | >= 0.90 | `sim_report.json` → `umi.duplex_conversion_rate` |
| AB/BA strand concordance (per variant) | >= 0.95 | `sim_report.json` → `strand_concordance.per_variant[*].strand_concordance` |
| Overall strand concordance | >= 0.95 | `sim_report.json` → `strand_concordance.overall` |
| VAF accuracy (all variants) | |obs - exp| < 0.20 × exp | truth VCF `N_ALT_MOL` vs `N_DUPLEX_ALT` |
| On-target fraction | >= 0.95 | `sim_report.json` → `capture.achieved_on_target_fraction` |
| Coverage CV | <= 0.25 | `sim_report.json` → `capture.achieved_coverage_cv` |

## How to run

1. Install VarForge and obtain a reference genome (hg38 recommended):

```
cargo install --path .
```

2. Run the benchmark scenario:

```bash
varforge simulate --config examples/twist_duplex_benchmark.yaml \
  --set reference=/path/to/hg38.fa \
  --set output_dir=/tmp/twist_benchmark \
  --set targets_bed=examples/twist_duplex_benchmark_targets.bed
```

3. Inspect the simulation report:

```bash
cat /tmp/twist_benchmark/sim_report.json | python3 -m json.tool
```

## Interpreting sim_report.json

The simulation report is written to `<output_dir>/sim_report.json`. For a duplex
run, it contains the following fields relevant to the benchmark:

```json
{
  "umi": {
    "enabled": true,
    "family_size_mean": 3.5,
    "duplex_total_molecules": 1234567,
    "duplex_conversion_rate": 1.0
  },
  "strand_concordance": {
    "overall": 1.0,
    "per_variant": [
      {
        "chrom": "chr1",
        "pos": 12345678,
        "n_alt_mol": 42,
        "n_duplex_alt": 42,
        "strand_concordance": 1.0
      }
    ]
  }
}
```

`duplex_conversion_rate` is the fraction of simulated molecules for which both
AB and BA families were generated. In VarForge, this is always 1.0 because both
strands are simulated deterministically. In real data, conversion rates of 0.85–0.95
are typical due to library preparation losses.

`strand_concordance` per variant measures the fraction of alt molecules where both
strands carry the variant. A value below 0.95 indicates that a variant may be
filtered by duplex callers as a strand artefact.

## Duplex consensus calling

### fgbio

```bash
# Step 1: sort by queryname for GroupReadsByUmi
samtools sort -n -o qname_sorted.bam output/TWIST_DUPLEX.bam

# Step 2: group reads into UMI families
fgbio GroupReadsByUmi \
  --input=qname_sorted.bam \
  --strategy=adjacency --edits=1 \
  --output=grouped.bam

# Step 3: call duplex consensus reads
fgbio CallDuplexConsensusReads \
  --input=grouped.bam \
  --output=consensus.bam \
  --min-reads=1 1 1 \
  --min-base-quality=30

# Step 4: align consensus reads
bwa mem hg38.fa consensus.bam | samtools sort -o consensus_aligned.bam
samtools index consensus_aligned.bam
```

### HUMID

```bash
# HUMID handles duplex consensus natively from FASTQ
humid --umi-length 9 \
  output/TWIST_DUPLEX_R1.fastq.gz \
  output/TWIST_DUPLEX_R2.fastq.gz \
  --out-prefix humid_consensus
```

## Variant calling on consensus output

After duplex consensus calling, run a standard somatic variant caller:

```bash
# Mutect2 (GATK4)
gatk Mutect2 \
  -R hg38.fa \
  -I consensus_aligned.bam \
  -O raw_calls.vcf.gz \
  --min-base-quality-score 30

# Evaluate against truth VCF
bcftools stats \
  --regions-file examples/twist_duplex_benchmark_targets.bed \
  output/truth.vcf raw_calls.vcf.gz
```

## Notes on the synthetic BED file

`examples/twist_duplex_benchmark_targets.bed` contains 100 synthetic capture
targets on chr1, each 150 bp wide, spaced every 1 Mbp starting at position
10 000 000. These coordinates are present in hg38. Replace this file with the
actual Twist panel BED when benchmarking against real data.
