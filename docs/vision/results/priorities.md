# VarForge: Results Priorities

## Must have

**1. Downstream tool compatibility**

This is the primary quality gate. Show that:
- Mutect2 calls variants from VarForge BAMs at expected sensitivity and precision.
- fgbio groups UMI families correctly on duplex VarForge output.
- hap.py concordance against the truth VCF is high.

If these experiments fail, VarForge is not ready to release.

**2. Performance benchmarks**

- Throughput in reads per second at 30x, 100x, and 300x WGS depths.
- Memory usage at each depth.
- Thread scaling curve.

These support the streaming architecture claim in the paper.

**3. Feature comparison table**

VarForge vs BAMSurgeon vs ART vs NEAT across all major feature dimensions. This is the paper's main evaluation figure.

## Lower priority

- Detailed error model accuracy. VarForge is not competing on this axis. Basic validation is enough.
- Absolute runtime speed comparisons. We claim to be fast enough, not fastest.
- Comparison with tools outside the main three competitors.
