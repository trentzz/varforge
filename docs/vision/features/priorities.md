# VarForge: Feature Priorities

## Must finish before release

These block the release. Nothing ships until all of these are done.

**In-progress features**
- T066: Streaming BAM editor. Required for spike-in use case.
- T067: ClonalTree integration into the simulation engine. Required for subclonal architecture.
- T068: Cancer preset driver mutations. Required for preset-based workflows.
- T004: Wire duplex family count into AppliedVariant. Required for correct truth VCF output.

**Critical bug fixes**
- UMI RX:Z tag not written to BAM in production. Breaks fgbio compatibility.
- Truth VCF not sorted correctly. Breaks tabix indexing and downstream tools.
- `Normal::new().unwrap()` panics when standard deviation is zero. Crashes on edge-case configs.
- Platform name case not normalised. GATK rejects BAMs with incorrect `@RG PL` values.
- NM tag inaccurate for variant reads. Causes incorrect edit distance reporting.
- Remove ContaminationConfig or replace with a stub that errors clearly. Currently it silently does nothing.

## Important but not blocking

Needed before the release is considered stable and paper-ready.

- Integration tests for all major features: SNV, indel, SV, UMI, cfDNA, artefacts, copy number.
- Module-level documentation for the 22 undocumented modules.
- Remove or justify the 69 `#[allow(dead_code)]` annotations.
- Add `@PG` program record to BAM header for provenance.
- Automatic BAM index and tabix index generation after output is written.
