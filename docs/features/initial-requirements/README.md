# VarForge Initial Requirements

Feature requirements derived from the [research documentation](../../research/).

## Documents

| Document | Scope |
|----------|-------|
| [Core Simulation](core-simulation.md) | Read generation, fragment models, coverage, quality scores, small dataset modes |
| [Variant Spike-in](variant-spike-in.md) | SNV, indel, MNV, SV, CNV injection with stochastic VAF, random mutation generation |
| [Tumour Model](tumour-model.md) | Purity, clonal architecture, ploidy, copy number |
| [UMI & Duplex](umi-duplex.md) | UMI barcodes, simplex/duplex modes, PCR families |
| [Liquid Biopsy](liquid-biopsy.md) | cfDNA fragmentation, ctDNA enrichment, low TF support |
| [Artifacts](artifacts.md) | FFPE damage, oxoG, GC bias, PCR duplicates/errors |
| [I/O Formats](io-formats.md) | YAML config, VCF/BED input, FASTQ/BAM/truth VCF output |
| [CLI Interface](cli-interface.md) | Subcommands, threading, progress, logging |
| [Performance](performance.md) | Parallelization, streaming, memory, throughput targets, crate leverage policy |

## Requirement ID Convention

Each requirement uses the prefix `REQ-<AREA>-<NNN>`:

- `SIM` — core simulation
- `VAR` — variant spike-in
- `TUM` — tumour model
- `UMI` — UMI / duplex
- `LBX` — liquid biopsy
- `ART` — artifacts
- `IOF` — I/O formats
- `CLI` — CLI interface
- `PRF` — performance

## Priority Levels

- **P0** — Must-have for initial release
- **P1** — Important, implement if feasible in first release
- **P2** — Planned for future release
