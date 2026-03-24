# VarForge: Target Audience

## Primary: bioinformatics developers

Teams building or benchmarking variant callers, UMI tools, and cfDNA analysis pipelines. They need ground truth data with known mutations at known VAFs to verify that their tools work correctly. Target tools include Mutect2, VarDict, Strelka2, fgbio, UMI-tools, CNVkit, and samtools.

Their pain: generating realistic cancer sequencing test data currently requires chaining multiple tools (ART + BAMSurgeon + custom scripts). VarForge replaces that pipeline with a single command.

## Secondary: clinical laboratory teams

Teams validating sequencing pipelines for regulatory or accreditation purposes. They need reproducible synthetic data that exercises the full pipeline without touching real patient samples. Seed-based reproducibility and a stable truth VCF are critical for this use case.

## What both groups need

- Output compatibility with standard bioinformatics tools (BAM, FASTQ, VCF formats).
- Seed-based reproducibility: same config, same seed, same output.
- Realistic error profiles so tools behave as they would on real data.
- Comprehensive ground truth: a VCF listing every simulated variant with VAF, depth, and strand information.
