# VarForge Paper: Focus

## The story

Existing synthetic data tools either lack cancer-specific features (ART, NEAT) or require complex multi-tool pipelines (BAMSurgeon + ART + custom scripts). VarForge provides an integrated, single-command solution covering the full spectrum of cancer sequencing simulation with proper ground truth.

## Key claims

1. VarForge is the first single-tool solution covering tumour simulation, UMI barcoding, cfDNA profiles, and artefact injection together.
2. Output is compatible with standard downstream tools. Demonstrate with Mutect2, fgbio, and hap.py.
3. The streaming architecture scales to whole-genome depths without proportional memory growth.
4. Per-read Bernoulli VAF sampling is more realistic than deterministic spike-in. Observed VAF follows expected binomial variance.

## What to emphasise

- Breadth of features in a single tool.
- Downstream compatibility with real pipelines.
- Reproducibility via seed-based control.
- Ease of use: YAML config and cancer-type presets.

## What to downplay

- Error model sophistication. VarForge is not competing with ART on platform-level read modelling.
- Long-read support. Basic only. Do not invite comparison with NanoSim or PBSIM3.

## Tone

State what VarForge does and show it works. Let the feature table and validation results carry the argument. Do not oversell.
