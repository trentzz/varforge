# VarForge: Project Goals

VarForge generates realistic synthetic sequencing data (FASTQ and BAM) with controlled ground truth for benchmarking variant calling pipelines. Success means bioinformatics teams can validate their tools against known truth without needing real patient samples.

## What it is

A single-binary, pure-Rust simulation tool. It produces FASTQ and BAM with configurable mutations, tumour parameters, UMI barcodes, cfDNA profiles, and sequencing artefacts. Users declare their simulation in a YAML config file and get reproducible output.

## What it is not

- Not a platform-level error modeller. ART and PBSIM own that space.
- Not a population genetics simulator.
- Not a clinical data generator for patient records or clinical workflows.

## Current state

Feature-complete core simulation. Paper drafted. 106 tasks completed. Three features are in-progress (BAM editor, ClonalTree integration, cancer presets). Several critical bugs need fixing. The next phase is stability, thorough testing against downstream tools, and release alongside the paper.

## Measure of success

Bioinformatics developers and laboratory teams can run Mutect2, fgbio, UMI-tools, hap.py, and other standard tools against VarForge output and get expected results. The ground truth VCF matches what the tools recover.
