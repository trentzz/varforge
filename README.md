# VarForge

A Rust tool for generating synthetic cancer sequencing test data for bioinformatics tool validation.

## Motivation

Testing bioinformatics tools (variant callers, UMI deduplicators, error correction pipelines, etc.) requires well-characterized test data with known ground truth. Existing simulators are fragmented across many tools, mostly written in Python/R, and none adequately handle the full spectrum of modern assay types -- particularly UMI-tagged duplex sequencing and liquid biopsy (cfDNA) data.

VarForge aims to be a single, fast, comprehensive tool that generates realistic synthetic sequencing data with:

- Controlled mutation spike-in (SNVs, indels, MNVs, SVs, CNVs) at specified VAFs
- Tumour purity/fraction and clonal architecture modeling
- Realistic Illumina error profiles (learned from real data or parametric)
- UMI-tagged reads (simplex and duplex sequencing)
- cfDNA fragment size distributions for liquid biopsy simulation
- Library prep artifact simulation (FFPE damage, oxidative damage, GC bias)
- PCR duplicate modeling with realistic family size distributions
- Multiple output formats: FASTQ, BAM, truth VCF

## Target Use Cases

Testing and benchmarking tools such as:

- **Variant callers**: Mutect2, VarDict, Strelka2, Varlociraptor
- **K-mer analysis**: km
- **BAM manipulation**: samtools
- **UMI deduplication**: HUMID, fgbio, UMI-tools
- **Error correction**: Duplex consensus pipelines
- **Copy number**: CNVkit, FACETS

## Project Status

**Research phase.** See [docs/research/](docs/research/) for detailed analysis of existing tools, identified gaps, and technical requirements.

## Research Documentation

| Document | Description |
|----------|-------------|
| [Existing Tools](docs/research/existing-tools.md) | Comprehensive catalog of current simulation/benchmarking tools |
| [Gap Analysis](docs/research/gap-analysis.md) | What's missing from the current landscape |
| [Technical Requirements](docs/research/technical-requirements.md) | What's needed to build VarForge |
| [Cancer Genomics Details](docs/research/cancer-genomics-details.md) | Domain knowledge: variants, cfDNA, UMIs, duplex sequencing |
| [Rust Ecosystem](docs/research/rust-ecosystem.md) | Available Rust crates for bioinformatics |
| [Input Specification Design](docs/research/input-specification.md) | Proposed configuration format |

## License

[MIT](LICENSE)
