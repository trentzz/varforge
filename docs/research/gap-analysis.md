# Gap Analysis: What's Missing from Current Tools

*Updated March 2025 with comprehensive tool landscape research*

## Critical Gaps

### 1. No tool combines cfDNA fragmentation + variant spike-in

**No dedicated cfDNA read simulator exists** that models the characteristic ~167bp nucleosomal fragmentation pattern, jagged ends, tumor fraction mixing, or fragment-end motifs. Existing benchmarking relies on wet-lab reference materials (Seraseq ctDNA controls, synthetic spike-ins into DNA-free plasma) rather than in silico simulation.

- Fragmentstein converts tabulated fragment coordinates into BAM but does not generate fragments de novo
- FinaleToolkit analyzes real cfDNA patterns but does not simulate them
- shendurelab/cfDNA extracts distributions from real data, no variant spike-in
- Sonicated DNA has "end properties different from naturally occurring cfDNA" -- naive fragmentation is insufficient
- ichorCNA needs low-coverage WGS with realistic coverage variation and tumor fraction
- DELFI needs whole-genome fragmentation profiles with realistic nucleosome positioning

**Impact**: ctDNA variant callers, liquid biopsy tools (ichorCNA, DELFI), and cfDNA fragmentomics tools have no in silico simulation capability for comprehensive benchmarking.

### 2. No dedicated UMI/duplex sequencing read simulator

**Zero tools** provide genome-wide UMI/duplex sequencing simulation:
- Inline UMI barcodes of configurable length
- Duplex strand pairing (alpha/beta UMI swapping)
- Realistic PCR family size distributions (log-normal)
- UMI sequencing errors at configurable rates
- Variant spike-in that correctly appears on both strands (true mutations) vs one strand (artifacts)

- UMI-tools has single-locus simulation only
- Kennedy Lab pipeline has basic test data, not a general simulator
- scReadSim handles UMI for scRNA-seq only, not targeted DNA
- BMC Genomics 2024 UMI benchmarking paper used manual ad hoc construction

**Impact**: fgbio, UMI-tools, HUMID, UMIErrorCorrect, UMI-VarCal, and duplex consensus callers use circular benchmarking or ad hoc data.

### 3. Spike-in bias (deterministic VAF) in most tools

BAMSurgeon and SomatoSim use deterministic spike-in: 10% VAF at 30x = exactly 3 alt reads. Real sequencing is stochastic (binomial distribution). Only stochasticSim (2023) addresses this but is limited to SNVs only and exome-scale.

MOV&RSim (2025) analyzed 9 somatic simulators and found **none provided complete control over both biological and technical parameters**.

**Impact**: Benchmarks overestimate variant caller performance with unrealistically clean data.

### 4. No Rust-based simulator exists

The entire landscape is Python/R/C++/Java:
- Python tools (BAMSurgeon, NEAT, SomatoSim): slow for WGS-scale
- C++ tools (ART, ReSeq): fast but unmaintained or limited scope
- Java tools (VarSim): complex dependencies, JVM overhead
- No single-binary distribution tool exists

rusty-neat exists but is incomplete (no BAM output, no model learning, 20 open issues).

**Impact**: Simulation remains a bottleneck in CI/CD. Users must manage Python environments or C build dependencies.

### 5. No single tool handles the full pipeline

Current workflow chains 3-5 separate tools:
1. Generate/mutate reference genome (SigProfilerSimulator, SomaticSiMu)
2. Simulate reads (ART, NEAT, ReSeq)
3. Spike mutations into BAM (BAMSurgeon, SomatoSim)
4. Add UMIs (manual scripting)
5. Simulate cfDNA fragmentation (manual scripting)

Each tool has different formats, dependencies, and documentation quality. Integration is fragile.

**Impact**: Hours of manual pipeline setup for each benchmarking experiment.

---

## Secondary Gaps

### 6. Library prep artifact simulation is fragmented

- FFPE damage (C>T deamination): only stochasticSim models this
- Oxidative damage (8-oxoG, G>T): only stochasticSim
- GC bias: only ReSeq models this realistically
- PCR duplicates: only ReSeq uses negative binomial model
- PCR errors during amplification: only GENOMICON-Seq models the full amplification pipeline
- No tool combines all artifact types in one configurable system

### 7. No integrated clonal architecture + from-scratch read generation

PSiTE and HeteroGenesis model clonal evolution but delegate read generation to external simulators (ART, wgsim). The handoff is manual and error-prone. VarSim integrates somewhat but is Java-based and unmaintained.

### 8. Poor support for targeted panels and capture simulation

Most simulators assume WGS. Targeted panel simulation needs:
- BED-defined target regions with capture efficiency variation
- Off-target read fractions
- Target boundary coverage dropoff
- Probe-specific biases (only GENOMICON-Seq models this for amplicon/WES)

### 9. No multi-sample / longitudinal simulation

Liquid biopsy monitoring involves serial samples. Simulating a time series with consistent clonal architecture but varying tumour fractions is not supported by any single tool.

### 10. Most tools are unmaintained

| Tool | Last Release | Status |
|------|-------------|--------|
| ART | 2016 | Abandoned |
| pIRS | 2016 | Abandoned |
| Mason | 2014 | Dormant |
| VarSim | 2020 | Stale |
| BAMSurgeon | 2022 | Low activity |
| SomatoSim | 2022 | Inactive |

The field needs actively maintained tools with modern interfaces.

---

## What VarForge Addresses

| Gap | Priority | VarForge Approach | Status |
|-----|----------|-------------------|--------|
| cfDNA + variant spike-in | **P0** | Nucleosomal mixture fragment model + VAF-aware generation | Core model implemented |
| UMI/duplex simulation | **P0** | Inline UMI, family size modeling, strand pairing | Core model implemented |
| Stochastic VAF sampling | **P0** | Binomial sampling for alt read count | Implemented |
| Single Rust binary | **P0** | Pure Rust with noodles, rayon, gzp | Scaffolded |
| Full pipeline in one tool | **P0** | Config → reads → FASTQ/BAM + truth VCF | Pipeline not yet wired |
| Library artifacts | **P1** | Configurable FFPE, oxoG, PCR errors/duplicates | Models implemented |
| Clonal architecture | **P1** | Tree-based clone definitions with CCF validation | Implemented |
| Targeted panels | **P1** | BED-file region selection with intersection | Implemented |
| Multi-sample series | **P2** | Shared clonal tree, per-sample tumour fraction | Not yet implemented |
| Error profile learning | **P2** | Learn quality/error profiles from real BAM | Not yet implemented |
| GC bias modeling | **P2** | Fragment-level GC content bias | Not yet implemented |
| SV/CNV support | **P2** | Structural variants, copy number variants | Not yet implemented |
| Cancer-type presets | **P2** | Predefined configs for common cancer types | Not yet implemented |

## Competitive Positioning

VarForge's unique value proposition vs every existing tool:

1. **vs BAMSurgeon**: From-scratch generation (no input BAM needed), stochastic VAF, 100x faster (Rust), UMI + cfDNA support
2. **vs ART/DWGSIM**: Variant injection, cancer modeling, UMI, cfDNA, artifacts -- not just raw read simulation
3. **vs NEAT**: Compiled Rust (10-100x faster), cancer modeling, UMI, cfDNA, artifacts
4. **vs ReSeq**: Variant injection, cancer modeling, UMI, cfDNA -- ReSeq only learns and replays profiles
5. **vs stochasticSim**: Indels + SVs + CNVs (not just SNVs), UMI, cfDNA, from-scratch generation, much faster
6. **vs VarSim**: Single binary (no Java/Python), UMI, cfDNA, artifacts, actively maintained
7. **vs MOV&RSim**: UMI/duplex support, cfDNA fragmentation, Rust performance, comprehensive artifact modeling
8. **vs GENOMICON-Seq**: UMI/duplex, cfDNA, clonal architecture, broader variant types, Rust performance
9. **vs PSiTE**: Integrated read generation (no external simulator needed), UMI, cfDNA, single binary
