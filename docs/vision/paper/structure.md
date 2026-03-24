# VarForge Paper: Structure

## Narrative arc

Problem: no single tool handles the full cancer sequencing simulation stack.
Gap: existing tools require multi-step pipelines and lack UMI, cfDNA, and artefact support.
Solution: VarForge integrates everything in one tool with a single config file.
Proof: downstream tools (Mutect2, fgbio, hap.py) produce expected results on VarForge output.
Limitations: honest account of what VarForge does not do.

## Sections

1. **Introduction**: the problem space. Variant callers, UMI tools, and cfDNA pipelines need ground truth data. Real patient data cannot be shared. Existing synthetic tools fall short.

2. **Background**: survey of existing tools. BAMSurgeon (tumour spike-in, requires real BAM input), ART (read simulation, no cancer features), NEAT (synthetic reads, limited cancer modelling). What each lacks.

3. **Method**: VarForge architecture.
   - Streaming pipeline: bounded channel between rayon workers and a writer thread.
   - Config system: YAML with optional VCF for variant lists.
   - Variant injection: per-read Bernoulli sampling for stochastic VAF.
   - UMI model: simplex and duplex families.
   - cfDNA fragment model: configurable length distributions.
   - Artefact model: FFPE and oxoG.
   - Cancer presets: driver mutations by cancer type.

4. **Results**:
   - Mutect2 variant recovery across VAF range.
   - fgbio UMI grouping accuracy on duplex data.
   - hap.py concordance against truth VCF.
   - Throughput and memory scaling benchmarks.

5. **Evaluation**: feature comparison table. VarForge vs BAMSurgeon vs ART vs NEAT.

6. **Discussion**: design trade-offs (pure Rust, streaming, Bernoulli VAF). Limitations. When to use VarForge and when not to.

7. **Future Work**: clonal evolution models, native long-read support, cloud output.

8. **Conclusion**: single paragraph. What was done and why it matters.
