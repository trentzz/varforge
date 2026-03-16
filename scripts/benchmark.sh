#!/usr/bin/env bash
# End-to-end benchmark script for VarForge.
#
# Usage:
#   ./scripts/benchmark.sh [--threads N] [--output-dir DIR]
#
# What it measures:
#   1. Criterion micro-benchmarks (fragment sampling, quality, variants, FASTQ, UMI)
#   2. Thread-scaling: runs the micro-benchmarks at 1/2/4/8/16 threads
#   3. Macro-benchmarks: simulates small-panel, WES, and chr22-WGS scenarios using
#      `varforge simulate` (if a reference FASTA is available) and records
#      wall-clock time + peak RSS.
#
# Requirements:
#   - cargo (Rust toolchain)
#   - GNU time (for memory measurement, optional)
#   - A reference FASTA at ${REFERENCE_FASTA} (optional; skipped if absent)

set -euo pipefail

###############################################################################
# Configuration
###############################################################################

THREADS="${VARFORGE_THREADS:-16}"
OUTPUT_DIR="${VARFORGE_BENCH_OUT:-/tmp/varforge_bench_$(date +%Y%m%d_%H%M%S)}"
REFERENCE_FASTA="${VARFORGE_REFERENCE:-}"  # e.g. /data/hg38/hg38.fa
BINARY="${CARGO_TARGET_DIR:-target}/release/varforge"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

###############################################################################
# Parse arguments
###############################################################################

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)   THREADS="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --reference) REFERENCE_FASTA="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "${OUTPUT_DIR}"
LOG="${OUTPUT_DIR}/benchmark.log"
echo "VarForge benchmark run — $(date)" | tee -a "${LOG}"
echo "Output directory: ${OUTPUT_DIR}" | tee -a "${LOG}"
echo "" | tee -a "${LOG}"

###############################################################################
# Helper: portable time + peak RSS
###############################################################################

# Check whether GNU time is available (macOS ships /usr/bin/time which lacks -f)
if command -v gtime &>/dev/null; then
    TIME_CMD=(gtime -f "  elapsed=%e s  peak_rss=%M KB")
elif /usr/bin/time --version &>/dev/null 2>&1; then
    TIME_CMD=(/usr/bin/time -f "  elapsed=%e s  peak_rss=%M KB")
else
    # Fallback: no memory measurement
    TIME_CMD=(time)
fi

###############################################################################
# 1. Build release binary
###############################################################################

echo "=== Building release binary ===" | tee -a "${LOG}"
cd "${REPO_ROOT}"
cargo build --release 2>&1 | tail -5 | tee -a "${LOG}"
echo "" | tee -a "${LOG}"

###############################################################################
# 2. Criterion micro-benchmarks
###############################################################################

echo "=== Criterion micro-benchmarks ===" | tee -a "${LOG}"
echo "Running: cargo bench" | tee -a "${LOG}"
cargo bench 2>&1 | tee "${OUTPUT_DIR}/criterion.log"
echo "Criterion output saved to ${OUTPUT_DIR}/criterion.log" | tee -a "${LOG}"
echo "" | tee -a "${LOG}"

###############################################################################
# 3. Thread-scaling micro-benchmark
#    VarForge uses rayon; RAYON_NUM_THREADS controls the thread pool size.
#    We run the fragment_sampling and quality_generation groups across thread
#    counts so the scheduling overhead is visible even in micro-benchmarks.
###############################################################################

echo "=== Thread-scaling (fragment_sampling + quality_generation) ===" | tee -a "${LOG}"
SCALING_LOG="${OUTPUT_DIR}/thread_scaling.log"
echo "threads,group,median_ns" > "${SCALING_LOG}"

for t in 1 2 4 8 ${THREADS}; do
    echo "  threads=$t" | tee -a "${LOG}"
    RAYON_NUM_THREADS="${t}" cargo bench \
        -- fragment_sampling quality_generation \
        2>&1 \
        | grep -E 'time:\s+\[' \
        | awk -v t="${t}" '
            /fragment_sampling/ { group="fragment_sampling" }
            /quality_generation/ { group="quality_generation" }
            /time:/ {
                # Extract median from "[ lower  median  upper ]"
                match($0, /\[([0-9.]+) ([a-z]+)  ([0-9.]+) ([a-z]+)  ([0-9.]+)/, m)
                print t "," group "," m[3] m[4]
            }
        ' >> "${SCALING_LOG}" 2>/dev/null || true
done
echo "Thread-scaling results saved to ${SCALING_LOG}" | tee -a "${LOG}"
echo "" | tee -a "${LOG}"

###############################################################################
# 4. Macro-benchmarks (requires reference FASTA and varforge binary)
###############################################################################

run_macro() {
    local name="$1"
    local config="$2"
    local out_dir="${OUTPUT_DIR}/macro_${name}"
    mkdir -p "${out_dir}"
    echo "  Running macro-benchmark: ${name}" | tee -a "${LOG}"
    "${TIME_CMD[@]}" "${BINARY}" simulate \
        --config "${config}" \
        --output-dir "${out_dir}" \
        --threads "${THREADS}" \
        2>&1 | tee "${out_dir}/run.log" | tail -5 | tee -a "${LOG}"
    echo "" | tee -a "${LOG}"
}

generate_panel_config() {
    local ref="$1"
    local out="$2"
    cat > "${out}" <<YAML
reference: "${ref}"
output:
  prefix: panel_bench
sample:
  name: panel_bench
  purity: 0.5
  coverage: 500
  read_length: 150
fragment:
  model: normal
  mean: 300
  sd: 50
targets:
  - { chrom: chr17, start: 7_669_000, end: 7_676_000 }   # TP53 region (~7 kb)
  - { chrom: chr7,  start: 55_086_000, end: 55_088_000 } # EGFR region (~2 kb)
  - { chrom: chr12, start: 25_378_000, end: 25_380_000 } # KRAS region (~2 kb)
variants:
  - { chrom: chr17, pos: 7674220, ref: C, alt: T, vaf: 0.25 }
  - { chrom: chr7,  pos: 55086971, ref: T, alt: A, vaf: 0.10 }
umi:
  enabled: true
  length: 8
  duplex: false
YAML
}

generate_wes_config() {
    local ref="$1"
    local out="$2"
    cat > "${out}" <<YAML
reference: "${ref}"
output:
  prefix: wes_bench
sample:
  name: wes_bench
  purity: 0.3
  coverage: 100
  read_length: 150
fragment:
  model: normal
  mean: 300
  sd: 50
targets:
  # Approximate WES: chr17 coding exons (~1 Mb as a proxy for a fast run)
  - { chrom: chr17, start: 1_000_000, end: 2_000_000 }
variants:
  - { chrom: chr17, pos: 1_500_000, ref: A, alt: G, vaf: 0.15 }
YAML
}

generate_chr22_config() {
    local ref="$1"
    local out="$2"
    cat > "${out}" <<YAML
reference: "${ref}"
output:
  prefix: chr22_30x
sample:
  name: chr22_30x
  purity: 1.0
  coverage: 30
  read_length: 150
fragment:
  model: normal
  mean: 300
  sd: 50
YAML
}

echo "=== Macro-benchmarks ===" | tee -a "${LOG}"

if [[ -z "${REFERENCE_FASTA}" ]]; then
    echo "  SKIPPED: set VARFORGE_REFERENCE=/path/to/hg38.fa to enable macro-benchmarks" | tee -a "${LOG}"
else
    if [[ ! -f "${REFERENCE_FASTA}" ]]; then
        echo "  SKIPPED: reference not found at ${REFERENCE_FASTA}" | tee -a "${LOG}"
    else
        PANEL_CFG="${OUTPUT_DIR}/panel_config.yaml"
        WES_CFG="${OUTPUT_DIR}/wes_config.yaml"
        CHR22_CFG="${OUTPUT_DIR}/chr22_config.yaml"

        generate_panel_config "${REFERENCE_FASTA}" "${PANEL_CFG}"
        generate_wes_config   "${REFERENCE_FASTA}" "${WES_CFG}"
        generate_chr22_config "${REFERENCE_FASTA}" "${CHR22_CFG}"

        run_macro "small_panel"  "${PANEL_CFG}"
        run_macro "wes_100x"     "${WES_CFG}"
        run_macro "chr22_30x"    "${CHR22_CFG}"
    fi
fi
echo "" | tee -a "${LOG}"

###############################################################################
# 5. Summary
###############################################################################

echo "=== Summary ===" | tee -a "${LOG}"
echo "All outputs written to: ${OUTPUT_DIR}" | tee -a "${LOG}"
echo "" | tee -a "${LOG}"

# Extract and print the key throughput lines from criterion output
if [[ -f "${OUTPUT_DIR}/criterion.log" ]]; then
    echo "Key throughput results (from Criterion):" | tee -a "${LOG}"
    grep -E 'thrpt:' "${OUTPUT_DIR}/criterion.log" \
        | sed 's/^/  /' \
        | tee -a "${LOG}" \
        || echo "  (no throughput measurements found)" | tee -a "${LOG}"
fi

echo "" | tee -a "${LOG}"
echo "Benchmark complete — $(date)" | tee -a "${LOG}"
