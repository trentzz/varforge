# VarForge: Comparison Strategy

## BAMSurgeon

The closest competitor for tumour simulation. BAMSurgeon injects variants into existing BAMs and supports SNVs, indels, and SVs. It requires a real input BAM, does not handle UMI barcodes, cfDNA profiles, or artefact injection, and requires Python dependencies.

Claim: VarForge matches BAMSurgeon's variant injection quality while adding UMI, cfDNA, artefacts, and a simpler single-command workflow.

## ART

The standard Illumina read simulator. ART produces realistic read error profiles but has no cancer-specific features: no VAF control, no tumour purity, no UMI, no cfDNA.

Claim: VarForge handles cancer-specific simulation that ART cannot. ART remains better for platform error modelling when that is the primary concern.

## NEAT

Another synthetic data tool with some variant injection capability. Less widely used than BAMSurgeon. Include for completeness.

## Comparison axes

| Axis | VarForge | BAMSurgeon | ART | NEAT |
|------|----------|------------|-----|------|
| SNV/indel/SV simulation | | | | |
| Tumour purity and subclones | | | | |
| UMI barcoding | | | | |
| cfDNA fragment profiles | | | | |
| FFPE/oxoG artefacts | | | | |
| Cancer presets | | | | |
| Ground truth VCF output | | | | |
| Single-command workflow | | | | |
| Seed reproducibility | | | | |
| Pure-binary install | | | | |

Fill this table with actual values when running comparisons.
