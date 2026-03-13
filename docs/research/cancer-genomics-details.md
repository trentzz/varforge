# Cancer Genomics Domain Knowledge

## Variant Types in Sequencing Data

### SNVs (Single Nucleotide Variants)
Single-base mismatches vs. reference. In BAM, appear as quality-weighted mismatches at a genomic coordinate. In FASTQ, indistinguishable from sequencing errors until statistical analysis. A fraction of reads carry the alt allele (= VAF).

### Small Indels
Represented via CIGAR operators `I` (insertion) and `D` (deletion). Complexities:
- **Left-alignment normalization**: Same indel can be represented at multiple positions in repetitive regions
- Reads spanning indels may be soft-clipped rather than correctly aligned
- Local realignment needed to resolve ambiguous placements
- Length distribution: mostly 1-2bp, up to ~50bp for "small" indels

### MNVs (Multi-Nucleotide Variants)
Two or more adjacent SNVs on the same haplotype. Critical for codon-level interpretation: two adjacent SNVs reported separately may predict different amino acid changes than when considered together. The key simulation requirement is that alt alleles must appear on the **same reads** (phased in cis), not independently distributed.

### Structural Variants

| SV Type | Discordant Pairs | Split Reads | Read Depth |
|---------|-----------------|-------------|------------|
| Deletion | Insert size larger than expected | Breakpoint-spanning reads split | Decreased in deleted region |
| Tandem Duplication | Unexpected orientation/insert size | Split at duplication boundaries | Increased in duplicated region |
| Inversion | Same-strand orientation | Primary + supplementary in opposite orientations | No change |
| Translocation | Mates on different chromosomes | Soft-clipped at breakpoints; supplementary on different chr | No local change |

Key BAM features: supplementary alignments (SA tag), soft-clipping (S in CIGAR), discordant read pairs.

### Copy Number Variants
Detected through read depth analysis:
- Log2 ratio: `log2(tumor_depth / normal_depth)` -- ~0.58 for 1-copy gain, large negative for deletions
- B-allele frequency (BAF) at het SNP sites distinguishes allele-specific events (LOH vs balanced gains)
- Segmentation algorithms (CBS) identify contiguous regions of consistent copy number

---

## The VAF Equation

```
Expected_VAF = (CCF * multiplicity * purity) / (purity * CN_tumor + (1 - purity) * CN_normal)
```

| Parameter | Description | Typical values |
|-----------|-------------|----------------|
| CCF | Cancer cell fraction | 0-1 (1 = clonal) |
| multiplicity | Copies of mutant allele per carrier cell | 1 (het), 2 (hom) |
| purity | Tumour cell fraction in sample | 0.3-0.9 (tissue), 0.001-0.05 (liquid biopsy) |
| CN_tumor | Total copy number at locus in tumour | 2 (diploid), 1-6+ |
| CN_normal | Total copy number in normal | 2 |

**Examples**:
- Clonal het SNV, diploid, 60% purity: VAF = 1 * 1 * 0.6 / (0.6 * 2 + 0.4 * 2) = 0.30
- Clonal het SNV, diploid, 1% tumour fraction (liquid biopsy): VAF = 0.005
- Subclonal (CCF=0.3) het SNV, diploid, 60% purity: VAF = 0.3 * 0.6 / 2 = 0.09

---

## Liquid Biopsy / cfDNA

### Fragment Size Distribution
Follows a nucleosomal protection pattern:

- **Primary peak: ~167bp** -- DNA wrapped around one nucleosome (~147bp) + ~20bp linker DNA
- **Sub-nucleosomal shoulder: ~143-147bp** -- nucleosome core without linker
- **Di-nucleosomal: ~334bp**
- **Tri-nucleosomal: ~500bp**
- **10bp periodicity** below 167bp reflecting DNA helical pitch on nucleosome surface

### ctDNA vs Normal cfDNA
- ctDNA fragments are **shorter** (enriched below 167bp, peak ~143bp)
- Distinct end-motif preferences
- Aberrant nucleosome positioning from cancer cells
- Similar pattern to fetal cfDNA vs maternal cfDNA

### Typical Tumour Fractions
| Context | Typical TF |
|---------|-----------|
| Advanced/metastatic | 5-30% |
| Stage II-IV detection | 1-10% |
| Stage I detection | 0.1-1% (often below detection) |
| MRD monitoring | 0.001-0.01% |
| Clinical "high TF" threshold | >1% |

Below 1% TF, standard sequencing error rates (~0.1-1%) mask true variants. UMI/duplex sequencing is required.

---

## UMI and Duplex Sequencing

### UMI Encoding
UMIs are short random barcodes (typically 4-12bp) ligated before PCR:
- **Inline**: First N bases of the read are the UMI
- **Index read**: UMI in a separate Illumina index read
- **Read name**: Often encoded as `@READNAME:UMISEQUENCE`

Read structure notation: `M` = molecular barcode, `S` = skip/spacer, `T` = template
- Example: `8M+T,+T` = 8bp UMI + template on R1, template only on R2

### Duplex Sequencing
Both strands of the original molecule are tagged with coordinated barcodes:

1. Ligate adapters with random UMIs to both ends of each fragment
2. Top strand gets UMI pair `A-B`, bottom strand gets `B-A` (swapped)
3. After PCR and sequencing, group reads by UMI + position
4. Build single-strand consensus (SSCS) for each UMI family
5. Compare alpha (A-B) and beta (B-A) SSCS to create duplex consensus (DCS)
6. True mutations appear on **both** strands; errors appear on one

### Error Suppression Performance

| Method | Error Rate |
|--------|-----------|
| Raw Illumina | ~10^-3 |
| Single-strand consensus (SSCS) | ~3.4 x 10^-5 |
| Duplex consensus (DCS) | ~10^-7 to 10^-10 |

DCS eliminates first-cycle PCR artifacts and DNA damage artifacts that SSCS cannot catch.

### Family Size Distribution
- Without error correction: dominated by singletons (unique barcodes)
- With UMI error correction (1-3 mismatch tolerance): singletons collapse into families
- Typical peak: 2-10 reads per family, right-skewed distribution
- Depends on: library complexity, sequencing depth, PCR cycles

### Common UMI Lengths

| Length | Unique Sequences | Typical Use |
|--------|-----------------|-------------|
| 4bp per end | 256 per end | Original duplex (Schmitt et al.) |
| 6-8bp | 4K-65K | Standard simplex dedup |
| 8bp | 65,536 | Most commercial panels (IDT xGen) |
| 12bp | ~16.7M | High-complexity applications |

---

## Tools That Consume Synthetic Data

### km (iric-soft/km)
- K-mer-based targeted variant detection
- Input: Jellyfish k-mer databases from FASTQ + target FASTA sequences
- Use cases: NPM1 insertions, FLT3-ITD in leukemia
- Test edge cases: variable-length insertions at same locus, low-VAF variants

### samtools
- General BAM manipulation, depth, mpileup, stats
- Test edge cases: supplementary alignments, unmapped mates, duplicate flags, UMI tags

### HUMID
- Reference-free FASTQ UMI deduplication
- Input: FASTQ with inline UMIs
- Test edge cases: UMI collisions, UMI sequencing errors, family size variation

### fgbio
- Full UMI pipeline: extract, group, consensus call, filter
- Input: uBAM or aligned BAM with RX tags
- Test edge cases: mixed family sizes, UMI errors causing fragmentation, supplementary alignments

### UMI-tools
- UMI extraction (FASTQ) and network-based deduplication (BAM)
- Test edge cases: high collision rates, PCR chimeras

### Variant Callers

| Caller | Strengths | Weaknesses | Key Test Scenarios |
|--------|-----------|------------|-------------------|
| Mutect2 | Best at VAF <10%, Bayesian model | Slower, needs matched normal or PoN | Low-VAF (0.1-1%), clustered mutations, multiallelic sites |
| VarDict | Good low-VAF, native MNV support | Higher false-positive rate | Very low VAF, complex indels, MNVs in same codon |
| Strelka2 | Fast (20x faster than Mutect2), good at VAF >=20% | Degrades below 10% VAF, needs matched normal | Performance across VAF spectrum |

---

## Benchmarking Standards

### GIAB
- High-confidence truth sets for HG001-HG007
- Somatic benchmarking: NA12878/NA24385 mixtures at known ratios
- Maintained by NIST

### SEQC2
- HCC1395 (TNBC) + HCC1395BL (matched normal)
- ~37K SNVs + ~58K indels truth set
- Controlled tumour fraction mixtures available

### ICGC-TCGA DREAM Challenge
- BAMSurgeon-generated synthetic tumours
- 248 SNV submissions, 204 SV submissions evaluated
- Three synthetic tumour-normal pairs of increasing difficulty
- Data: EGA EGAS00001002092
