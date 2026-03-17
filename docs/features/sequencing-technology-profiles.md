# Sequencing Technology Profiles

VarForge must simulate data that matches real library preparation protocols.
Different protocols produce different UMI structures, fragment size distributions, and read structures.
This document specifies the technology profiles VarForge should support.

## Current support

VarForge currently supports:
- Configurable UMI length (1 to 16 bp)
- Simplex and duplex UMI modes
- Inline UMI (prepended to read) or read-name UMI
- Configurable PCR cycle count and family size distribution
- cfDNA nucleosomal mixture model
- FFPE and oxoG artefact injection

## Technology-specific gaps

### Twist Biosciences UMI Duplex

Twist uses a 5 bp random UMI followed by a 2 bp spacer on each read.
The read structure is `5M2S+T` (5 bases molecular barcode, 2 bases skip, then template).
Duplex identification uses canonical pairing of the two UMI sequences.

What VarForge needs:
- **Spacer support**: configurable spacer length between UMI and template (currently not modelled)
- **Read structure string**: accept a read structure specification like `5M2S+T` to describe UMI, spacer, and template layout
- The existing 5 bp UMI + duplex mode covers the core functionality

### IDT xGen UMI Adapters

IDT uses a 3 bp fixed UMI in the i5 index read, not inline.
Some versions use 3+3 bp dual-indexed UMI.

What VarForge needs:
- **Index-read UMI**: UMI stored in a separate index read rather than inline (currently only inline or read-name)
- This is a lower priority since most analysis pipelines extract UMIs to the RX tag before processing

### Illumina TruSight Oncology (TSO 500)

TSO 500 uses dual UMI: 7 bp UMI in each adapter.
Read structure includes a fixed 4 bp anchor sequence.

What VarForge needs:
- **Dual UMI**: independent UMI on R1 and R2 adapters (partially supported via duplex mode)
- **Anchor sequences**: fixed sequences adjacent to UMI (not modelled)

### Single-strand library prep

Used for ancient DNA and some cfDNA protocols.
All molecules are single-stranded before adapter ligation.
No duplex information is available.

VarForge's simplex mode already handles this case.
No changes needed.

## Priority order

1. **Spacer support** (low effort, enables Twist protocol simulation)
2. **Read structure specification** (medium effort, generalises UMI layout)
3. **Index-read UMI** (medium effort, enables IDT xGen)
4. **Anchor sequences** (low effort, enables TSO 500)

## Current workarounds

For Twist duplex simulation today, users can:
- Set `umi.length: 5` and `umi.duplex: true`
- The 2 bp spacer is not simulated but does not affect downstream analysis if the pipeline is configured to skip the spacer bases

For IDT xGen simulation today:
- Set `umi.inline: false` (UMI stored in read name)
- The analysis pipeline extracts UMIs from read names, which matches the xGen workflow
