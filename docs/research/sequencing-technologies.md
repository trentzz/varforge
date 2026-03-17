# Sequencing Technology UMI Profiles

Research on UMI and duplex sequencing protocols relevant to VarForge simulation.

## Twist Bioscience UMI Adapter System

Twist adapters contain a 5 bp random UMI and a 10 bp unique dual index per adapter.
The read structure is `5M2S+T`: 5 bases molecular barcode, 2 bases spacer (skipped), then template.
Duplex identification uses canonical pairing of the UMI sequences from both adapter ends.

Recommended workflow (from Twist tech note DOC-001337):
1. fgbio `FastqToBam` to extract UMIs from read sequences
2. BWA-MEM alignment
3. fgbio `GroupReadsByUmi` with `--strategy paired`
4. fgbio `CallMolecularConsensusReads` for simplex consensus
5. fgbio `CallDuplexConsensusReads` for duplex consensus
6. Variant calling on consensus BAM

VarForge simulation parameters for Twist:
- `umi.length: 5`
- `umi.duplex: true`
- `umi.inline: true`
- Missing: 2 bp spacer between UMI and template (needs `umi.spacer_length` field)

Sources:
- https://www.twistbioscience.com/products/ngs/library-preparation/twist-umi-adapter-system
- https://www.twistbioscience.com/resources/guideguideline/processing-sequencing-data-utilizing-twist-unique-molecular-identifier-umi

## IDT xGen UDI-UMI Adapters

IDT uses a 9 bp degenerate UMI positioned adjacent to the i7 index.
The UMI is sequenced as part of the i7 index read (requires 9 extra index cycles).
Dual 8 bp unique dual indexes (UDI) handle index hopping.

The UMI is not inline in the read sequence. It is in the index read.

VarForge simulation parameters for IDT xGen:
- `umi.length: 9`
- `umi.duplex: false` (simplex UMI only)
- `umi.inline: false` (UMI in read name, simulating extraction from index read)

Sources:
- https://www.idtdna.com/pages/support/faqs/how-do-i-sequence-the-umi-in-the-xgen-udi-umi-adapters
- https://www.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-library-preparation/ngs-adapters-indexing-primers/adapters-indexing-primers-for-illumina

## Illumina TruSight Oncology 500 (TSO 500)

TSO 500 uses dual UMIs, one from each adapter end.
UMI sequences are extracted during FASTQ generation (BCL Convert, not bcl2fastq).
The combined UMI appears in the read name as `UMI1+UMI2` (e.g., `CGAACGT+GTACACG`).
Each UMI is 7 bp.

VarForge simulation parameters for TSO 500:
- `umi.length: 7`
- `umi.duplex: true` (dual UMI from both adapter ends)
- `umi.inline: false` (UMI extracted to read name)

Sources:
- https://www.illumina.com/products/by-type/clinical-research-products/trusight-oncology-umi.html
- https://support.illumina.com/sequencing/sequencing_kits/trusight-oncology-500.html

## Agilent SureSelect XT HS2

Uses molecular barcodes (MBCs) of 6 bp on each adapter end (6+6 = 12 bp total).
MBCs are inline at the 5' end of each read.
Supports both simplex and duplex consensus calling.

VarForge simulation:
- `umi.length: 6`
- `umi.duplex: true`
- `umi.inline: true`

## NEB NEBNext UltraFS with UMI

Uses 4 bp inline UMI on each read.
Combined with unique dual indexing.
Typical applications: low-input cfDNA, FFPE.

VarForge simulation:
- `umi.length: 4`
- `umi.duplex: true`
- `umi.inline: true`

## Safe-SeqS

Uses 12-14 bp random UMI (endogenous or exogenous).
Simplex only. Single-end or paired-end.
Original method by Kinde et al. (2011).

VarForge simulation:
- `umi.length: 14`
- `umi.duplex: false`
- `umi.inline: true`

## Summary of UMI configurations

| Protocol | UMI length | Position | Duplex | Spacer |
|----------|-----------|----------|--------|--------|
| Twist UMI | 5 bp | inline | yes | 2 bp |
| IDT xGen | 9 bp | index read | no | none |
| TSO 500 | 7 bp | read name | yes | none |
| Agilent HS2 | 6 bp | inline | yes | none |
| NEB UltraFS | 4 bp | inline | yes | none |
| Safe-SeqS | 12-14 bp | inline | no | none |

## What VarForge needs to add

1. `umi.spacer_length` config field (integer, default 0). Inserts N random bases between UMI and template. Enables Twist protocol.
2. All other protocols are already supported by the existing `umi.length`, `umi.duplex`, and `umi.inline` configuration fields.
