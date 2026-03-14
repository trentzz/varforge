# Rust Bioinformatics Ecosystem

## BAM/SAM/CRAM

### rust-htslib
- **URL**: https://github.com/rust-bio/rust-htslib
- **Type**: FFI bindings to the C htslib library (same engine as samtools/bcftools)
- **Features**:
  - BAM, SAM, CRAM read/write
  - Multi-threaded BAM compression/decompression
  - Indexed random access (fetch by region)
  - Serde serialization for `bam::Record`
  - Full CRAM support (bzip2 + lzma)
  - VCF/BCF read/write via `bcf` module
- **Maturity**: High. Battle-tested, wraps the reference C implementation.
- **Downside**: C dependency (requires htslib, zlib, bzip2, lzma to compile)
- **Recommended for VarForge**: Yes -- primary choice for BAM/VCF I/O

### noodles
- **URL**: https://github.com/zaeleus/noodles | https://docs.rs/noodles
- **Type**: Pure Rust implementation (no C dependencies)
- **Supported formats**: BAM 1.6, SAM 1.6, CRAM 3.0/3.1, BGZF, VCF 4.3/4.4, BCF 2.2, FASTA, FASTQ, BED, GFF3, GTF, CSI, tabix, htsget, refget
- **Features**: Async I/O with Tokio, spec-compliant
- **Maturity**: Experimental (API explicitly described as experimental, frequent breaking changes)
- **Recommended for VarForge**: Use selectively -- `noodles-fasta` for reference access (pure Rust, no C dep), but prefer rust-htslib for BAM/VCF

---

## FASTQ

| Crate | Speed | Notes |
|-------|-------|-------|
| **seq_io** | Fastest | Zero-allocation iteration. Also writes. Best for raw throughput. |
| **needletail** | Fast | MIT-licensed, minimal-copy parser with k-mer support. |
| **bio::io::fastq** | Moderate | Part of rust-bio. Simple API. Slower due to allocations and UTF-8 checks. |
| **noodles-fastq** | Good | Part of noodles ecosystem. Async support. |
| **fastq** (fastq-rs) | Fast | Built-in threaded decompression. |

**Recommended for VarForge**: `seq_io` for FASTQ writing performance.

---

## VCF/BCF

### rust-htslib::bcf
- Full VCF/BCF read/write via htslib bindings
- Indexed random access, multi-threaded I/O
- BCF has zero parsing overhead
- **Recommended**: Yes, for reading mutation specification VCFs and writing truth VCFs

### noodles-vcf / noodles-bcf
- Pure Rust VCF 4.3/4.4 and BCF 2.2
- Async I/O
- **Maturity**: Experimental

---

## FASTA (Reference Genome)

| Crate | Notes |
|-------|-------|
| **noodles-fasta** | Pure Rust, indexed FASTA (.fai) support, async. **Recommended.** |
| **bio::io::fasta** | Simple reader/writer with index support. Part of rust-bio. |
| **seq_io** | Fast FASTA parsing. |

All support `.fai` indexed random access, essential for extracting reference sequences at mutation sites.

---

## Algorithms and Utilities

### rust-bio
- **URL**: https://rust-bio.github.io/ | https://github.com/rust-bio
- Pairwise alignment (Smith-Waterman, Needleman-Wunsch, banded)
- Sequence utilities, pattern matching
- FASTA/FASTQ/BED I/O
- **Recommended**: Yes, for local realignment after mutation spiking (CIGAR recalculation)

### rand + rand_distr
- **URL**: https://crates.io/crates/rand, https://crates.io/crates/rand_distr
- Random number generation with distributions: Normal, LogNormal, Beta, Gamma, Poisson, Binomial, Bernoulli, etc.
- **Critical for**: Fragment size sampling, error injection, VAF stochastic sampling, UMI generation, PCR amplification modeling
- **Recommended**: Yes, core dependency

### statrs
- Statistical distributions and functions (complementary to rand_distr)
- Probability density functions, CDF, quantile functions

### rayon
- Data parallelism (parallel iterators)
- **Critical for**: Parallel read generation across genomic regions
- **Recommended**: Yes

---

## Compression

| Crate | Purpose | Notes |
|-------|---------|-------|
| **flate2** | DEFLATE/gzip/zlib | Standard Rust compression. Supports miniz_oxide (pure Rust) and zlib-ng backends. Single-threaded. |
| **gzp** | Multi-threaded gzip/BGZF | Drop-in replacement for `Write` that compresses blocks in parallel (pigz-style). Supports BGZF block format (64 KB blocks) used in BAM. Uses libdeflater backend by default for best performance. **Recommended for FASTQ gzip and BGZF.** |
| **libdeflater** | High-performance DEFLATE | Faster than flate2 for known-size inputs. Used as gzp backend. |
| **bgzip** | BGZF format | Pure Rust BGZF reader/writer. Alternative to gzp for BAM-style block gzip. |

---

## BED File Parsing

| Crate | Purpose | Notes |
|-------|---------|-------|
| **noodles-bed** | BED reader/writer | Part of noodles ecosystem. Pure Rust, spec-compliant. **Recommended.** |
| **rust-bio** | BED I/O | Basic BED reader in `bio::io::bed`. Simpler API but less feature-complete. |

---

## CLI and Configuration

| Crate | Purpose |
|-------|---------|
| **clap** | CLI argument parsing with derive macros |
| **serde** | Serialization/deserialization framework |
| **serde_yaml** | YAML config parsing |
| **toml** | TOML config parsing (alternative) |

---

## Progress Reporting and Logging

| Crate | Purpose | Notes |
|-------|---------|-------|
| **indicatif** | Progress bars and spinners | Configurable styles, multi-bar support, ETA estimation. Widely used. **Recommended.** |
| **tracing** | Structured logging/diagnostics | Span-based instrumentation, multiple subscribers, async-aware. Superior to log+env_logger for complex apps. **Recommended over log.** |
| **tracing-subscriber** | tracing output formatters | JSON output, filtering, layered subscribers. |
| **env_logger** | Simple env-based log filtering | Simpler alternative if tracing is overkill. |

---

## Existing Rust Bioinformatics Tools (Reference)

### Varlociraptor
- **URL**: https://github.com/varlociraptor/varlociraptor
- Unified variant caller for SNVs, MNVs, indels, SVs in Rust
- Grammar-based scenario configuration for tumour/normal, pedigrees
- Most significant Rust tool in cancer variant calling
- Demonstrates that complex cancer bioinformatics tools can be built successfully in Rust

### biotest
- **URL**: https://github.com/natir/biotest
- Generates random test data for bioinformatics (FASTA, FASTQ, VCF, sequences, quality strings)
- Very basic; not cancer-specific
- Part of the rust-bio ecosystem

### rustybam
- **URL**: https://github.com/vollgerlab/rustybam
- Bioinformatics toolkit for manipulation of BAM, BED, FASTA, FASTQ files
- Useful reference for Rust BAM/BED manipulation patterns

### Notable absence
**No Rust-based cancer sequencing simulator exists.** The ecosystem has strong infrastructure (I/O, calling) but no simulation tools. VarForge would be the first.

---

## What VarForge Must Build (No Existing Crate)

These capabilities have no existing Rust crate and must be implemented:

- **GC bias coverage modeling** — no Rust crate models GC-dependent coverage curves; must implement from literature (e.g., ReSeq/NEAT approach)
- **Illumina error profiles** — no Rust crate provides position-dependent quality/error models; must learn from BAM or use parametric model
- **cfDNA fragment size distribution** — nucleosomal periodicity model must be built from scratch
- **Mixture distributions** — rand_distr provides individual distributions (Normal, LogNormal, Binomial, NegativeBinomial) but no mixture/composite; must implement weighted mixture sampling
- **Variant spike-in engine** — core novel contribution; no existing crate
- **Tumour clonal architecture** — no existing Rust crate
- **UMI/duplex simulation** — no existing crate for UMI family generation with PCR amplification modeling
- **FFPE/oxoG artifact injection** — no existing crate

---

## Recommended Stack for VarForge

```toml
[dependencies]
rust-htslib = "0.47"        # BAM/VCF I/O, CIGAR manipulation, aux tags
noodles-fasta = "0.42"      # Reference genome indexed access (pure Rust)
noodles-bed = "0.15"        # BED file parsing (pure Rust)
seq_io = "0.3"              # Fast FASTQ writing
rust-bio = "1.6"            # Pairwise alignment, BED I/O
rand = "0.8"                # RNG with seedable generators
rand_distr = "0.4"          # Normal, LogNormal, Binomial, NegBinomial, etc.
statrs = "0.17"             # PDF/CDF/quantile functions
rayon = "1.10"              # Data parallelism (per-region work stealing)
gzp = "0.11"                # Multi-threaded gzip/BGZF compression
flate2 = "1.0"              # Single-threaded gzip (fallback)
clap = { version = "4", features = ["derive"] }
serde = { version = "1", features = ["derive"] }
serde_yaml = "0.9"
serde_json = "1"            # Simulation manifest output
tracing = "0.1"             # Structured logging
tracing-subscriber = "0.3"  # Log output formatting
indicatif = "0.17"          # Progress bars
```

Note: Version numbers are approximate; check crates.io for latest versions at implementation time.
