# Downstream Variant Caller Benchmark Design

## Purpose

VarForge claims to produce synthetic reads with ground-truth variants suitable
for benchmarking somatic variant callers. This benchmark validates that claim
directly: inject known variants, call with Mutect2, measure how many are found.

This is the highest-value validation experiment because it mimics exactly how
a real user would use VarForge.

---

## Metrics

### Primary metrics

**Sensitivity** (recall): fraction of truth variants detected.

    sensitivity = TP / (TP + FN)

A variant is a true positive (TP) if the caller reports a variant at the same
position, with the same alt allele, and the genotype (VAF estimate) is within
a reasonable tolerance.

**Precision**: fraction of calls that correspond to a truth variant.

    precision = TP / (TP + FP)

**F1**: harmonic mean of sensitivity and precision.

    F1 = 2 * (precision * sensitivity) / (precision + sensitivity)

### Stratification

All three metrics are computed per VAF tier:

| Tier | VAF range | Clinical analogue |
|------|-----------|-------------------|
| High | 0.20--0.50 | Primary tumour, high purity |
| Medium | 0.10--0.20 | Subclonal or lower purity |
| Low | 0.05--0.10 | Early detection, low purity |
| Very low | 0.01--0.05 | Liquid biopsy, ctDNA |

Stratification matters because callers have different sensitivities at different
VAFs. A tool that detects everything above 20% but misses sub-5% variants is
not useful for liquid biopsy workflows.

### Variant type breakdown

Report separately for:
- SNVs (the most common variant type, most callers handle well)
- Indels (harder to align and call accurately)
- MNVs (often split into separate SNVs by callers)

---

## Scenarios

### Scenario 1: SNV baseline (30x WGS)

**Goal**: establish the baseline sensitivity and precision of Mutect2 on clean
VarForge output with no artefacts.

**Config**: `configs/01_snv_baseline.yaml`

- Reference: chr22 (hg38, ~51 Mbp)
- Coverage: 30x
- Mutations: 500 random SNVs/indels, VAF uniform in [0.05, 0.50]
- Purity: 0.80 (high, to avoid purity masking low-VAF variants)
- Features: none (baseline)
- Seed: 42

**Expected outcome**: sensitivity >90% for VAF >0.20, >75% for VAF 0.05--0.10.
Precision >95% (no artefacts to produce false positives).

---

### Scenario 2: Low-VAF ctDNA (30x, cfDNA fragmentation)

**Goal**: test detection at very low VAF using cfDNA-like reads, simulating
circulating tumour DNA in liquid biopsy.

**Config**: `configs/02_low_vaf_ctdna.yaml`

- Coverage: 30x
- Mutations: 200 SNVs, VAF uniform in [0.01, 0.10]
- Purity: 0.05 (5% tumour fraction, realistic for ctDNA)
- Features: cfDNA fragmentation (4-component Gaussian mixture)
- Seed: 43

**Expected outcome**: sensitivity will drop significantly below VAF 0.05 at 30x;
this quantifies the coverage floor for ctDNA detection with standard callers.
The point of this scenario is to show what coverage is needed, not to prove
callers are perfect.

---

### Scenario 3: FFPE artefacts (30x)

**Goal**: test whether FFPE artefacts injected by VarForge cause false positives
in Mutect2, and whether the caller's built-in artefact filters handle them.

**Config**: `configs/03_ffpe_artefacts.yaml`

- Coverage: 30x
- Mutations: 300 SNVs, VAF in [0.10, 0.50]
- Features: FFPE damage rate 2%, oxoG rate 1%
- Seed: 44

**Expected outcome**: false positive rate increases relative to scenario 1.
Mutect2 artefact filtering (--germline-resource, orientation bias filter) should
reduce but not eliminate FFPE FPs. This characterises the artefact-rejection
capability of the caller on VarForge data.

---

### Scenario 4: Panel + UMI, high depth (200x, 1 Mbp BED target)

**Goal**: test high-depth targeted panel calling where UMI-based error correction
would be applied upstream of the caller.

**Config**: `configs/04_panel_umi.yaml`

- Coverage: 200x on a 1 Mbp BED target (chr22:30000000-31000000)
- Mutations: 50 SNVs within the target, VAF in [0.01, 0.20]
- Features: UMI simplex
- Seed: 45

**Note**: This scenario generates UMI-tagged reads. In a real pipeline, reads
would be grouped and deduplicated (fgbio or HUMID) before alignment. This
benchmark runs Mutect2 directly on the raw aligned reads (no UMI deduplication)
to establish a baseline. A follow-up should add the deduplication step and
compare.

---

## Pipeline

```
VarForge → FASTQ (R1, R2) + truth VCF
    ↓
bwa-mem2 mem -t {threads} {ref} R1.fq.gz R2.fq.gz
    ↓
samtools sort | samtools addreplacerg -r ID:1 -r SM:SAMPLE -r PL:ILLUMINA
    ↓
samtools markdup (flag duplicates, do not remove for Mutect2)
    ↓
samtools index
    ↓
gatk Mutect2 -R {ref} -I {bam} -tumor SAMPLE --output {vcf}
    ↓
gatk FilterMutectCalls -R {ref} -V {vcf} --output {filtered_vcf}
    ↓
rtg-tools vcfeval -b {truth_vcf} -c {filtered_vcf} -t {ref.sdf} -o {eval_dir}
    ↓
evaluate.py → summary.tsv
```

### Reference genome preparation

```bash
# Index for bwa-mem2
bwa-mem2 index chr22.fa

# Create sequence dictionary for GATK
gatk CreateSequenceDictionary -R chr22.fa

# Create RTG-tools SDF (sequence data format)
rtg format -o chr22.sdf chr22.fa
```

---

## Tool versions

All tools should be recorded in `results/run_<date>/versions.txt`. The pipeline
script does this automatically. Expected versions:

- bwa-mem2: ≥2.2.1
- samtools: ≥1.17
- gatk: ≥4.4.0
- rtg-tools: ≥3.12

---

## Interpreting results

The summary table (`summary.tsv`) has columns:
`scenario, tier, variant_type, tp, fp, fn, sensitivity, precision, f1`

Key questions to answer:
1. Does VarForge output pass through the full pipeline without errors? (correctness check)
2. Is sensitivity at VAF >0.20 above 90%? (sanity check: high-VAF variants should be trivially detectable)
3. Does sensitivity degrade gracefully at lower VAF, or does it cliff? (indicates alignment/calling issues)
4. Does the FFPE scenario show more FPs than baseline? By how much? (validates artefact simulation)

---

## Limitations

- This benchmark uses the same reference for generation and alignment. A caller
  optimised for a reference will have an advantage over callers tested on
  different references.
- Mutect2 is tested without a matched normal. In clinical practice, a matched
  normal eliminates many germline variants. Single-sample mode is more
  challenging and more conservative.
- UMI deduplication is not applied in scenario 4. Results for that scenario
  should be interpreted as a lower bound.
- The truth VCF from VarForge reflects injected variants. Some injected variants
  may fall in low-complexity regions where alignment itself fails; these will
  appear as false negatives regardless of caller quality.
