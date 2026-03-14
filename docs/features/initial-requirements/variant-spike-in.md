# Variant Spike-in Requirements

Controlled injection of mutations into synthetic reads with stochastic allele-frequency modeling.

## VAF Model

The expected VAF for a somatic variant is:

```
VAF = CCF × multiplicity × purity / (purity × CN_tumour + (1 − purity) × CN_normal)
```

Where:
- **CCF** — cancer cell fraction of the clone carrying the mutation
- **multiplicity** — number of copies of the mutant allele in tumour cells
- **purity** — tumour content of the sample
- **CN_tumour** — total copy number at that locus in tumour
- **CN_normal** — total copy number at that locus in normal (typically 2)

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-VAR-001 | Compute expected VAF from the equation above using per-variant tumour model parameters | P0 |
| REQ-VAR-002 | Sample observed variant read counts from a **binomial distribution** around the expected VAF (not deterministic rounding) | P0 |
| REQ-VAR-003 | Allow per-variant VAF override in the input VCF (bypass equation) | P1 |

## Random Mutation Generation

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-VAR-004 | Generate random mutations without requiring an input VCF — user specifies count, type distribution, and target VAF range | P0 |
| REQ-VAR-005 | Accept parameters: `--random-mutations N` with `--vaf-range MIN-MAX` (e.g., `--random-mutations 100 --vaf-range 0.0001-0.05`) | P0 |
| REQ-VAR-006 | Distribute random mutations across the genome (or target regions if BED provided), avoiding centromeres, gaps, and low-complexity regions | P0 |
| REQ-VAR-007 | Configurable mutation type mix: fraction of SNVs, indels, MNVs (default: 80% SNV, 15% indel, 5% MNV) | P1 |
| REQ-VAR-008 | Random indel size distribution: geometric or user-specified, with configurable max length | P1 |
| REQ-VAR-009 | All randomly generated mutations are recorded in the truth VCF identically to VCF-specified mutations | P0 |

## SNVs

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-VAR-010 | Spike single-nucleotide variants into reads overlapping the target position | P0 |
| REQ-VAR-011 | Modify base quality at spiked positions to match surrounding context (avoid artificial quality peaks) | P0 |

## Indels

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-VAR-020 | Spike insertions and deletions (1–1000 bp) into overlapping reads | P0 |
| REQ-VAR-021 | Recalculate CIGAR strings for reads carrying indels | P0 |
| REQ-VAR-022 | Adjust mate alignment coordinates when indels shift read boundaries | P0 |

## MNVs

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-VAR-030 | Spike multi-nucleotide variants (2+ adjacent substitutions) phased on the same reads | P0 |
| REQ-VAR-031 | Ensure MNV constituent bases are always co-present on a read (never partially spiked) | P0 |

## Structural Variants

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-VAR-040 | Generate split reads at SV breakpoints | P1 |
| REQ-VAR-041 | Generate discordant read pairs spanning SV breakpoints | P1 |
| REQ-VAR-042 | Produce supplementary alignments (SA tag) for split-read evidence | P1 |
| REQ-VAR-043 | Support SV types: deletions, duplications, inversions, translocations | P1 |
| REQ-VAR-044 | Support insertions of novel sequence at breakpoints | P2 |

## Copy Number Variants

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-VAR-050 | Simulate CNV gains by increasing local read depth proportional to the copy number ratio | P0 |
| REQ-VAR-051 | Simulate CNV losses (including LOH) by decreasing local read depth | P0 |
| REQ-VAR-052 | Apply CNV depth changes per-clone according to the clonal tree | P1 |

## Truth Output

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-VAR-060 | Emit a truth VCF listing every spiked variant with the **expected** and **achieved** VAF | P0 |
| REQ-VAR-061 | Include per-variant read support counts (REF, ALT) in truth VCF INFO/FORMAT fields | P0 |
| REQ-VAR-062 | Annotate truth VCF with clone assignment, CCF, and CN context | P1 |
