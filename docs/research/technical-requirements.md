# Technical Requirements

## Core Simulation Pipeline

### 1. Reference Genome Access
- Indexed FASTA (.fa + .fai) random access
- Extract sequences at mutation sites and read generation regions
- Crate: `noodles-fasta` (pure Rust, indexed access) or `rust-htslib`

### 2. Read Generation with Controlled Error Profiles

**Parametric model (default)**:
- Position-dependent base quality scores (quality degrades toward 3' end)
- Read 2 has higher error rates than Read 1
- Substitution errors are context-dependent (influenced by preceding/following bases)
- Quality score sampled from a distribution per position, then error applied probabilistically based on quality

**Empirical model (learned from real data)**:
- Parse a real BAM file to extract per-position quality score distributions
- Build substitution probability matrices per (position, context, quality) tuple
- Sample from these learned distributions during simulation

**Implementation**: Pre-compute error probability tables. Use `rand_distr` for sampling quality scores. Apply errors by sampling a uniform random against the position-specific error probability.

### 3. Variant Spike-in at Controlled VAFs

**The VAF equation**:
```
Expected_VAF = (CCF * multiplicity * purity) / (purity * CN_tumor + (1 - purity) * CN_normal)
```

Where:
- CCF = cancer cell fraction (0-1)
- multiplicity = copies of the mutant allele per cell carrying it
- purity = tumour cell fraction in sample
- CN_tumor = total copy number at locus in tumour
- CN_normal = total copy number in normal (usually 2)

**Simple case** (het SNV, diploid, clonal): `VAF = purity / 2`

**Stochastic sampling**: For each read overlapping a variant site, draw from Bernoulli(VAF) to decide if this read carries the alt allele. This naturally produces binomial-distributed alt read counts, unlike deterministic spike-in tools.

**For SNVs**: Replace the reference base with the alt base in selected reads. Adjust quality score if desired (matching the quality distribution of the surrounding bases).

**For indels**: Modify the read sequence and CIGAR string. Insertions add bases (CIGAR: `I`), deletions remove bases (CIGAR: `D`). After modification, realign locally using Smith-Waterman (rust-bio) to produce a realistic CIGAR string.

**For SVs**: Generate split reads (soft-clipped at breakpoints with supplementary alignments) and discordant read pairs (unexpected insert size or mate chromosome). Set SA tags for supplementary alignments.

### 4. Tumour Purity / Clonal Architecture

**Approach**: Define a clonal tree where each node (clone) has:
- A cellular prevalence (CCF)
- A set of mutations
- A parent clone (mutations are cumulative down the tree)

For a sample with purity `p`:
- Generate `p * coverage` reads from the tumour genome
- Generate `(1-p) * coverage` reads from the normal genome
- Mutations in clone `i` with CCF `c_i` appear at VAF = `c_i * p * multiplicity / (p * CN_tumor + (1-p) * 2)`

### 5. Fragment Size Distribution

**Standard WGS**: Normal distribution with configurable mean (typically ~350bp) and std dev (~50bp).

**cfDNA (liquid biopsy)**: Mixture of Gaussians modeling nucleosomal protection:
- Mono-nucleosomal peak: mean ~167bp, std ~15bp (dominant component, weight ~0.85)
- Sub-nucleosomal shoulder: mean ~143bp, std ~10bp (weight ~0.08)
- Di-nucleosomal peak: mean ~334bp, std ~25bp (weight ~0.05)
- Tri-nucleosomal: mean ~500bp, std ~30bp (weight ~0.02)

For ctDNA (tumour-derived fragments): shift peaks ~10-20bp shorter, increase weight of sub-nucleosomal component.

The 10bp periodicity (nucleosome winding) can be modeled as a sinusoidal modulation on the fragment length probability.

### 6. PCR Duplicate Simulation

- Each original template molecule gets an amplification factor from a log-normal or negative binomial distribution
- Parameters derived from the target duplicate rate
- All copies of a template share: same insert coordinates (start/end), same UMI (if applicable)
- Each copy gets independent sequencing errors
- Duplicate flag (`0x400`) is NOT set in output (that's the deduplicator's job)

### 7. UMI Tag Generation

**Simplex UMIs**:
- Generate random N-mer UMIs (configurable length, typically 8bp)
- Prepend to Read 1 (inline) or store in separate index read
- Each original molecule gets one UMI
- PCR copies share the same UMI
- Apply UMI sequencing errors at configurable rate (typically 1-2% per base)

**Duplex UMIs**:
- Each molecule gets a UMI pair: (A, B)
- Top strand reads tagged with `A+B` (or `A-B` with separator)
- Bottom strand reads tagged with `B+A` (swapped)
- Downstream tools (fgbio) use the swapped pattern to identify complementary strand families
- True mutations appear on both strands; artifacts appear on one

**BAM tags**: Write `RX` tag (raw UMI sequence), optionally `MI` tag (molecule index after grouping).

### 8. Library Prep Artifacts

**FFPE damage**: C>T / G>A transitions from cytosine deamination. Apply at configurable rate (e.g., 0-5% of cytosines). Strand-specific: deamination on one strand only.

**Oxidative damage (oxoG)**: G>T / C>A transversions at 8-oxoguanine sites. Configurable rate.

**GC bias**: Coverage multiplier as a function of fragment GC content. Model with a polynomial or lookup table. Can be learned from real data.

### 9. Output Formats

- **FASTQ**: Paired-end R1/R2 files with optional UMI in read name or inline
- **BAM**: Aligned reads with proper headers, read groups, CIGAR strings, mate information
- **Truth VCF**: Ground truth variants with actual achieved VAFs
- **Truth BED**: For CNVs, the regions with altered copy number
- **Manifest**: JSON/YAML summary of simulation parameters and statistics

---

## Performance Considerations

### Why Rust Matters

Simulating 200x WGS coverage of a 3Gb genome:
- ~2 billion read pairs at 150bp
- At 1M reads/sec (optimistic Python), that's ~33 minutes just for read generation
- Rust with rayon parallelism can realistically achieve 10-50M reads/sec
- BAM writing is I/O bound; rust-htslib supports multi-threaded BAM compression

### Parallelization Strategy

- Divide genome into non-overlapping regions (e.g., per chromosome or per megabase)
- Generate reads for each region independently using rayon
- Merge sorted output streams
- Random seed per region for reproducibility: `region_seed = hash(global_seed, region_id)`

### Memory Management

- Stream reads to output; don't hold entire dataset in memory
- Reference genome: memory-map the FASTA file
- Error profile tables: pre-compute and share across threads (read-only)
- Fragment buffer per thread: bounded queue of generated fragments

---

## Recommended Crate Stack

| Crate | Purpose |
|-------|---------|
| `rust-htslib` | BAM/CRAM read/write, VCF/BCF read/write |
| `noodles-fasta` | Reference genome indexed access (pure Rust) |
| `seq_io` | Fast FASTQ writing |
| `rust-bio` | Pairwise alignment (Smith-Waterman) for CIGAR recalculation |
| `rand` + `rand_distr` | All random sampling (fragment sizes, errors, VAF, UMIs) |
| `rayon` | Parallel read generation across genomic regions |
| `clap` | CLI argument parsing |
| `serde` + `serde_yaml` | Configuration file parsing |
| `log` + `env_logger` | Logging |
