# Downstream validation scope for release

**Category**: scope
**Related epic**: EPIC-RELEASE-GATE

## Context

The vision docs state: "No release until [Mutect2, fgbio, hap.py] experiments pass."
Scripts exist in `benchmarking/downstreamvc/` but have never been run. No results
exist. This is the primary quality gate.

Running all three tools requires:
- hg38 reference genome (~3 GB)
- GATK/Mutect2 (~5 GB with resources)
- fgbio (~200 MB)
- hap.py (~1 GB with dependencies)

## Options

1. **Full validation (all three tools)**: highest confidence but requires significant setup and compute time.
2. **Mutect2 only**: validates the core use case (somatic variant calling). fgbio and hap.py are secondary.
3. **Unit-level validation only**: verify VarForge output is well-formed (valid BAM headers, sorted coordinates, correct tags) without running downstream tools. Fastest but weakest.

## Recommendation

Option 2 first, then option 1 if time permits. Mutect2 is the most important
validation because it tests the full pipeline: BAM correctness, variant
positions, and truth VCF accuracy. fgbio validation (T142) can follow as a
separate task.
