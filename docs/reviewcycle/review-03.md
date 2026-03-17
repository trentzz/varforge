# Review Cycle 03: Structure, Benchmarks, and Sequencing Technology Audit

**Date:** 2026-03-17
**Version:** v0.1.0 (post-restructure)
**Scope:** Project structure compliance, benchmark coverage, sequencing technology support

---

## Structure audit

The project was restructured to follow `~/.claude/guides/research-project.md`.

**Compliant:**
- `benchmarking/` is top-level with `design/`, `results/graphs/`, `results/tables/`, `scripts/`
- `docs/planning/` replaces `docs/plan/`
- Paper uses `sections/` with one `.tex` per section
- Paper imports graphs from `../../benchmarking/results/graphs/` via relative paths
- Graph style follows Okabe-Ito palette, white background, no top/right spines
- Vancouver numbered citations via biblatex/biber
- Git committer set to `Trent Zeng <dev@trentz.me>`

**Still needed:**
- `benchmarking/archive/` exists but is empty. First benchmark set should be archived when superseded.
- `docs/dev/` and `docs/dev/scripts/` are empty. Build/test/run instructions should go here.
- No `docs/paper/figures/` directory for non-graph figures (pipeline diagram is TikZ, so this is acceptable for now).

---

## Benchmark coverage

### What is tested

1. **Coverage scaling** (1 MB synthetic ref, 1x to 200x): linear scaling confirmed, R^2 > 0.999
2. **Feature overhead** (1 MB synthetic ref, 30x): variant +3%, GC +2%, FFPE +16%, UMI +64%, all +125%
3. **Thread scaling** (10 MB synthetic ref, 30x): plateaus at 1.16x (I/O bound)
4. **hg38 benchmarks** (chr22, ~51 MB): WGS baseline, variants, panel UMI, Twist duplex, cfDNA, FFPE (running)

### What is not tested and should be

1. **VAF accuracy**: the stochastic VAF claim is a core contribution but has no quantitative validation in the benchmarks. Need to generate data at known VAF, count alt reads, and show the distribution matches Binomial(D, VAF).
2. **cfDNA fragment distribution**: the mixture model claim is central but not validated. Need to extract fragment lengths from cfDNA output and show they match the expected distribution.
3. **Truth VCF correctness**: no test verifies that all injected variants appear in the truth VCF with correct positions and alleles.
4. **Determinism**: no test verifies that the same seed produces identical output across runs.
5. **Multi-sample mode**: not benchmarked. Tumour-normal pair generation is listed in the benchmark design but not yet run.

### Benchmark methodology gap

The thread scaling benchmark shows poor results because the test workload is I/O bound.
This is an honest result for the tested configuration, but it should be supplemented with a compute-bound thread scaling test (e.g., 1 MB with all features enabled) to show that Rayon parallelism works when the bottleneck is compute.

---

## Sequencing technology audit

### Currently supported protocols

| Protocol | UMI len | Mode | Position | Fragment | Status |
|----------|---------|------|----------|----------|--------|
| Generic simplex | 1-16 bp | simplex | inline or read-name | configurable | supported |
| Generic duplex | 1-16 bp | duplex | inline or read-name | configurable | supported |
| cfDNA | N/A | N/A | N/A | nucleosomal mixture | supported |
| FFPE | N/A | N/A | N/A | normal (shortened) | supported |

### Protocols that need additional support

| Protocol | What is missing | Priority |
|----------|----------------|----------|
| Twist UMI Duplex | 2 bp spacer between UMI and template | medium |
| IDT xGen | UMI in index read (not inline) | low |
| Illumina TSO 500 | 4 bp anchor sequence adjacent to UMI | low |
| Agilent SureSelect XT HS2 | dual molecular barcode (6+6 bp) | low |

### Assessment

VarForge's configurable UMI length and duplex mode cover the core functionality of all major protocols.
The main gap is spacer/anchor sequences between UMI and template.
This is a small implementation change (add a `spacer_length` config field and prepend N bases after the UMI).

For Twist duplex specifically: setting `umi.length: 5` and `umi.duplex: true` produces data that is functionally correct.
The missing 2 bp spacer does not affect downstream analysis when the pipeline is configured to skip those bases.

---

## Correctness

- All 57 unit tests pass (`cargo test --all-targets`)
- Integration tests pass
- Clippy is clean
- No unsafe code
- Build succeeds on CI (self-hosted runner)

---

## Honesty check

### Claims that are well-supported
- "First genome-wide UMI and duplex simulator": no other tool does this. Supported.
- "First in-silico cfDNA generator with variant support": no other tool does this. Supported.
- "Streaming architecture bounds memory": coverage scaling data confirms 1.7 KB/pair constant. Supported.
- "Linear scaling with coverage": R^2 > 0.999. Supported.

### Claims that need stronger evidence
- "Stochastic VAF matching Binomial distribution": claimed but not yet validated quantitatively. Need a VAF accuracy benchmark.
- "Throughput substantially higher than Python-based alternatives": no direct comparison run. The claim is plausible (28K pairs/s in Rust vs typical Python throughput) but unquantified.

### Known limitations stated plainly
- Thread scaling is poor for I/O-bound workloads. The paper reports this honestly.
- UMI mode uses significantly more memory than baseline. The paper reports this.
- No long-read support. Listed in Future Work.
- No methylation support. Listed in Future Work.

---

## Planned changes

1. Add a VAF accuracy validation to the benchmark suite (generate data, count alt reads, chi-squared test against Binomial)
2. Add a cfDNA fragment distribution validation (extract lengths, compare to mixture model)
3. Add a compute-bound thread scaling test
4. Add `umi.spacer_length` config field (enables Twist protocol simulation)
5. Fill in `docs/dev/` with build/test/run instructions
