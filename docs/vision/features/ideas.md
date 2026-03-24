# VarForge: Post-release Ideas

These are not planned for v1. Collect them here so they are not lost.

- **Contamination simulation**: cross-sample contamination at a configurable rate. The config stub already exists.
- **Custom fragment length models**: user-supplied fragment length distributions as a lookup table or empirical file.
- **Chromosome name aliasing**: automatic translation between `chr1` and `1` naming styles.
- **Automatic BAM and tabix indexing**: run `samtools index` and `tabix` automatically after writing output files.
- **Cloud-native output**: write BAM and FASTQ directly to S3 or GCS buckets without local staging.
- **Web UI for config generation**: a browser-based form that produces a valid YAML config. Lowers the barrier for non-developer users.
- **VCF-to-VCF benchmarking mode**: skip alignment entirely and compare a truth VCF against a calls VCF directly. Useful for pure variant calling benchmarks.
- **Clonal evolution models**: simulate tumour evolution over time with branching clonal trees and selection pressure.
- **Native long-read error models**: proper PacBio kinetics and Nanopore squiggle-level artefacts.
