# EPIC-KAM-COMPAT: kam Downstream Compatibility

## Goal

Make VarForge's Twist duplex output directly consumable by
[kam](https://github.com/trentzz/kam) without post-processing.

kam is an alignment-free variant caller that processes Twist duplex UMI
sequencing data. It expects inline UMIs at the start of both R1 and R2 reads,
a realistic duplex conversion rate, and UMI-level sequencing errors. VarForge
currently lacks all three.

## Motivation

The downstream user of VarForge synthetic data is kam. For that pipeline to
work end-to-end:

1. R1 and R2 must carry the UMI (and optional spacer) prepended to the template
   sequence, matching Twist's `5M2S+T` FASTQ layout.
2. Only ~85–95% of molecules should produce both AB and BA strand families.
   Simulating 100% duplex conversion inflates duplex depth and produces
   unrealistically clean data.
3. UMI bases have sequencing errors (~0.1%), which affects UMI clustering.
   Without this noise, UMI error-correction algorithms are never exercised.

## What Success Looks Like

| Check | Target |
|-------|--------|
| R1 and R2 sequences begin with UMI bytes when `umi.inline: true` | exact |
| Optional spacer (e.g. `AT`) follows UMI in FASTQ | exact |
| `umi.duplex_conversion_rate` drives BA strand dropout | Bernoulli sampling |
| Reported `duplex_conversion_rate` in sim_report.json matches config | within 2% at n=10000 |
| `umi.error_rate` injects substitution errors into UMI sequences | rate within 2× |
| Twist preset updated: 5 bp UMI, `AT` spacer, 90% conversion, 0.1% UMI error | exact |

## Scope

- In scope: inline UMI in FASTQ and BAM, duplex conversion rate, UMI error
  injection, Twist preset update, integration tests.
- Out of scope: UMI allow-lists, per-position UMI quality models, methylation.

## Tasks

- T146: Implement inline UMI FASTQ output with optional spacer
- T147: Add configurable duplex conversion rate (stochastic BA strand dropout)
- T148: Wire UMI sequencing error injection into production path
- T149: Update Twist preset for kam compatibility
- T150: Integration tests for inline UMI, duplex conversion rate, UMI errors
