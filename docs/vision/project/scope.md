# VarForge: Scope

## In scope

- Variant simulation: SNVs, indels, MNVs, structural variants.
- Tumour purity and subclonal architecture.
- UMI simulation: simplex and duplex barcoding.
- cfDNA fragment profiles: configurable fragment length distributions.
- Sequencing artefacts: FFPE, oxidative damage (oxoG).
- Copy number alterations.
- GC bias simulation.
- Capture and panel simulation.
- Multi-sample and longitudinal series.
- BAM editing for spike-in (T066).
- Profile learning from real sequencing data.
- Cancer-type presets with driver mutations (T068).
- Clonal tree integration (T067).

## Out of scope

- Platform-specific error modelling beyond Illumina short-read basics. ART handles this.
- Population-scale simulation.
- Complex structural rearrangements: chromothripsis, templated insertions, and similar.
- RNA-seq simulation.
- Methylation simulation.
- Native long-read error models. PacBio and Nanopore have basic read-length support only. Deep platform modelling is deferred.
- Contamination simulation. The config field exists but the feature is not implemented. Deferred to a future version.
