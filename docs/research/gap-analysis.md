# Gap Analysis: What's Missing from Current Tools

## Critical Gaps

### 1. No tool combines cfDNA fragmentation + variant spike-in

No existing tool generates liquid biopsy data with both realistic cfDNA fragment size distributions (nucleosomal peaks at ~167bp, ~334bp) AND controlled somatic variant spike-in at tumour fractions typical of clinical cfDNA (0.1-5%). The closest is shendurelab/cfDNA for fragment simulation, but it has no variant spike-in capability.

**Impact**: Developers of ctDNA variant callers and liquid biopsy analysis pipelines cannot generate realistic end-to-end test data.

### 2. No dedicated UMI/duplex sequencing read simulator

There is no tool that generates realistic FASTQ data with:
- Inline UMI barcodes of configurable length
- Duplex strand pairing (alpha/beta UMI swapping)
- Realistic PCR family size distributions
- UMI sequencing errors at configurable rates
- Variant spike-in that correctly appears on both strands (true mutations) or one strand (artifacts)

UMI-tools has a single-locus simulation component. The Kennedy Lab pipeline has basic test data generation. Neither is a general-purpose simulator.

**Impact**: Tools like HUMID, fgbio, UMI-tools, and duplex consensus callers lack comprehensive test data.

### 3. Spike-in bias in existing tools

BAMSurgeon and SomatoSim use deterministic spike-in: a variant at 10% VAF with 30x depth always produces exactly 3 alt reads. Real sequencing involves stochastic sampling -- the actual number follows a binomial distribution. Only stochasticSim (2023) addresses this, but it has a small user base.

**Impact**: Benchmarks may overestimate variant caller performance because the test data is "too clean."

### 4. No Rust-based simulator

The entire simulation tool landscape is Python/R/C++. There is no Rust-based tool despite:
- Rust's performance advantages (critical for high-coverage WGS simulation -- e.g., 200x WGS generates billions of reads)
- Single-binary distribution (no Python environment management)
- Memory safety guarantees
- Excellent parallelism via rayon

**Impact**: Simulation remains a bottleneck in CI/CD pipelines for bioinformatics tools.

### 5. No single tool handles the full pipeline

Current workflow requires chaining multiple tools:
1. Generate mutated reference genome (SigProfilerSimulator, SomaticSiMu)
2. Simulate reads (ART, NEAT, ReSeq)
3. Spike mutations into BAM (BAMSurgeon, SomatoSim)
4. Add UMIs (manual scripting)
5. Simulate cfDNA fragmentation (manual scripting)

Each tool has its own input format, dependencies, and quirks. Integration is fragile.

---

## Secondary Gaps

### 6. Library prep artifact simulation is ad hoc

FFPE damage (C>T deamination), oxidative damage (G>T at 8-oxoG sites), and other library prep artifacts are important for testing error correction tools, but no simulator provides configurable, biologically-motivated artifact injection.

### 7. No integrated clonal architecture + read generation

Tools like PSiTE and HeteroGenesis model clonal evolution, but their read generation is delegated to external simulators (ART, wgsim). The handoff between evolution modeling and read generation is manual and error-prone.

### 8. Poor support for targeted panels

Most simulators assume WGS. Targeted panel simulation (with capture efficiency variation, off-target reads, target boundary effects) requires additional manual setup.

### 9. No support for multi-sample / longitudinal simulation

Liquid biopsy monitoring involves serial samples from the same patient. Simulating a time series with consistent clonal architecture but varying tumour fractions is not supported by any single tool.

---

## What VarForge Should Address

| Gap | Priority | Approach |
|-----|----------|----------|
| cfDNA + variant spike-in | **P0** | Nucleosomal fragment model + VAF-aware read generation |
| UMI/duplex simulation | **P0** | Inline UMI generation, family size modeling, strand pairing |
| Stochastic VAF sampling | **P0** | Binomial sampling for alt read count given target VAF and depth |
| Single Rust binary | **P0** | Rust implementation with rust-htslib, noodles, seq_io |
| Full pipeline in one tool | **P1** | Reference -> mutated genome -> reads -> BAM + truth VCF |
| Library artifacts | **P1** | Configurable FFPE, oxoG, deamination rates |
| Clonal architecture | **P1** | Tree-based clone definition with per-clone mutations and CCFs |
| Targeted panels | **P2** | BED-file-driven region selection with capture efficiency model |
| Multi-sample series | **P2** | Shared clonal tree, per-sample tumour fraction |
| Error profile learning | **P2** | Learn quality/error profiles from user-provided real BAM |
