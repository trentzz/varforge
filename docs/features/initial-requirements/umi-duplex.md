# UMI & Duplex Sequencing Requirements

Simulation of unique molecular identifiers and duplex sequencing protocols.

## UMI Barcode Generation

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-UMI-001 | Generate random inline UMI barcodes with configurable length (default 8 bp) | P0 |
| REQ-UMI-002 | Prepend UMIs to the read sequence (inline format) or store in BAM RX tag, configurable per run | P0 |
| REQ-UMI-003 | Support fixed UMI barcode lists (allow-list mode) in addition to random generation | P1 |

## Simplex UMI Mode

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-UMI-010 | Tag each original molecule with a single UMI | P0 |
| REQ-UMI-011 | All PCR copies of a molecule share the same UMI | P0 |

## Duplex UMI Mode

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-UMI-020 | Tag each original molecule with a paired UMI (AB + BA strand convention) | P0 |
| REQ-UMI-021 | Read 1 of alpha strand carries UMI-A+UMI-B; read 1 of beta strand carries UMI-B+UMI-A | P0 |
| REQ-UMI-022 | Both strands of a duplex molecule originate from the same genomic fragment | P0 |
| REQ-UMI-023 | Configurable alpha/beta strand recovery rate (not all molecules yield both strands) | P1 |

## PCR Amplification Families

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-UMI-030 | Simulate PCR duplication by generating multiple copies of each source molecule | P0 |
| REQ-UMI-031 | Family size drawn from a **lognormal** distribution (default) with configurable mean and sigma | P0 |
| REQ-UMI-032 | Alternative family size distribution: **negative binomial** with configurable r and p | P1 |
| REQ-UMI-033 | PCR copies share identical alignment positions and UMI, with optional small per-copy coordinate jitter | P0 |

## UMI Sequencing Errors

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-UMI-040 | Inject base-call errors into UMI sequences at a configurable rate (default matching read error rate) | P1 |
| REQ-UMI-041 | UMI errors should produce near-miss UMI families that test deduplication tool error correction | P1 |

## BAM Tag Conventions

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-UMI-050 | Write UMI sequences to the `RX` BAM tag (SAM spec standard) | P0 |
| REQ-UMI-051 | Optionally assign molecular identifiers to the `MI` BAM tag (integer grouping ID) | P1 |
| REQ-UMI-052 | Mark PCR duplicates via the SAM flag (0x400) when requested | P1 |

## Tool Compatibility

| ID | Requirement | Priority |
|----|-------------|----------|
| REQ-UMI-060 | Output must be consumable by **fgbio** GroupReadsByUmi and CallMolecularConsensusReads | P0 |
| REQ-UMI-061 | Output must be consumable by **UMI-tools** group and dedup | P0 |
| REQ-UMI-062 | Output must be consumable by **HUMID** for UMI clustering | P0 |
