#!/usr/bin/env bash
# Downstream variant caller benchmark pipeline.
# Generates reads with VarForge, aligns with bwa-mem2, calls with Mutect2,
# and evaluates sensitivity/precision against the truth VCF.
#
# Usage:
#   ./run_benchmark.sh [options]
#
# Options:
#   --varforge PATH   path to varforge binary (default: ../../target/release/varforge)
#   --ref PATH        path to chr22.fa (default: ../../benchmarking/chr22.fa)
#   --sdf PATH        path to chr22.sdf for rtg-tools (default: ../../benchmarking/chr22.sdf)
#   --outdir PATH     output directory (default: ../results/run_YYYYMMDD)
#   --threads N       thread count (default: 12)
#   --scenarios LIST  comma-separated list of scenarios to run (default: 1,2,3,4)
#   --skip-generate   skip VarForge generation (reuse existing FASTQs)
#
# Requires: bwa-mem2, samtools, gatk, rtg-tools, python3
# See setup.md for installation instructions.

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
VARFORGE="${VARFORGE:-$ROOT/target/release/varforge}"
REF="${REF:-$ROOT/benchmarking/chr22.fa}"
SDF="${SDF:-$ROOT/benchmarking/chr22.sdf}"
OUTDIR="${OUTDIR:-$SCRIPT_DIR/../results/run_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-12}"
SCENARIOS="${SCENARIOS:-1,2,3,4}"
SKIP_GENERATE="${SKIP_GENERATE:-false}"
CONFIG_DIR="$SCRIPT_DIR/../configs"

# ── Parse args ─────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --varforge)  VARFORGE="$2";      shift 2 ;;
        --ref)       REF="$2";           shift 2 ;;
        --sdf)       SDF="$2";           shift 2 ;;
        --outdir)    OUTDIR="$2";        shift 2 ;;
        --threads)   THREADS="$2";       shift 2 ;;
        --scenarios) SCENARIOS="$2";     shift 2 ;;
        --skip-generate) SKIP_GENERATE=true; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Checks ─────────────────────────────────────────────────────────────────────
check_tool() { command -v "$1" &>/dev/null || { echo "ERROR: $1 not found. See setup.md."; exit 1; }; }
check_tool bwa-mem2
check_tool samtools
check_tool gatk
check_tool rtg
check_tool python3

[[ -f "$VARFORGE" ]] || { echo "ERROR: varforge binary not found at $VARFORGE"; exit 1; }
[[ -f "$REF" ]]      || { echo "ERROR: reference not found at $REF"; exit 1; }
[[ -f "${REF}.bwt.2bit.64" ]] || { echo "ERROR: bwa-mem2 index missing. Run: bwa-mem2 index $REF"; exit 1; }
[[ -d "$SDF" ]] || { echo "ERROR: RTG SDF not found at $SDF. Run: rtg format -o $SDF $REF"; exit 1; }

mkdir -p "$OUTDIR"

# Record tool versions
{
    echo "varforge: $($VARFORGE --version 2>&1 | head -1)"
    echo "bwa-mem2: $(bwa-mem2 version 2>&1 | head -1)"
    echo "samtools: $(samtools --version | head -1)"
    echo "gatk:     $(gatk --version 2>&1 | head -1)"
    echo "rtg:      $(rtg version 2>&1 | grep Version | head -1)"
    echo "date:     $(date -Iseconds)"
    echo "ref:      $REF"
    echo "threads:  $THREADS"
} > "$OUTDIR/versions.txt"
echo "Tool versions recorded in $OUTDIR/versions.txt"

# ── Scenario definitions ────────────────────────────────────────────────────────
declare -A SCENARIO_CONFIG=(
    [1]="01_snv_baseline.yaml"
    [2]="02_low_vaf_ctdna.yaml"
    [3]="03_ffpe_artefacts.yaml"
    [4]="04_panel_umi.yaml"
)
declare -A SCENARIO_NAME=(
    [1]="SNV baseline"
    [2]="Low-VAF ctDNA"
    [3]="FFPE artefacts"
    [4]="Panel + UMI"
)

# ── Main loop ──────────────────────────────────────────────────────────────────
IFS=',' read -ra SELECTED <<< "$SCENARIOS"

for SNUM in "${SELECTED[@]}"; do
    SNUM="${SNUM// /}"
    CONFIG="${SCENARIO_CONFIG[$SNUM]:-}"
    NAME="${SCENARIO_NAME[$SNUM]:-}"
    [[ -z "$CONFIG" ]] && { echo "Unknown scenario: $SNUM (valid: 1,2,3,4)"; continue; }

    echo ""
    echo "══════════════════════════════════════════════"
    echo "  Scenario $SNUM: $NAME"
    echo "══════════════════════════════════════════════"

    SDIR="$OUTDIR/s${SNUM}"
    mkdir -p "$SDIR"/{generate,align,call,eval}

    CFGFILE="$SDIR/config.yaml"
    # Patch OUTPUT_DIR placeholder in config
    sed "s|OUTPUT_DIR|$SDIR/generate|g" "$CONFIG_DIR/$CONFIG" > "$CFGFILE"

    # ── Step 1: Generate reads ────────────────────────────────────────────────
    if [[ "$SKIP_GENERATE" == "true" ]] && [[ -f "$SDIR/generate/R1.fastq.gz" ]]; then
        echo "  [1/4] Skipping generation (reusing existing FASTQs)"
    else
        echo "  [1/4] Generating reads..."
        "$VARFORGE" simulate \
            --config "$CFGFILE" \
            --threads "$THREADS" \
            2>&1 | tee "$SDIR/generate.log"
        echo "        Done."
    fi

    # Locate output files (VarForge names them based on sample name)
    R1=$(ls "$SDIR"/generate/*.R1.fastq.gz 2>/dev/null | head -1)
    R2=$(ls "$SDIR"/generate/*.R2.fastq.gz 2>/dev/null | head -1)
    TRUTH=$(ls "$SDIR"/generate/*.truth.vcf.gz 2>/dev/null | head -1)
    [[ -z "$R1" ]]    && { echo "ERROR: R1 FASTQ not found in $SDIR/generate/"; exit 1; }
    [[ -z "$TRUTH" ]] && { echo "ERROR: truth VCF not found in $SDIR/generate/"; exit 1; }

    # ── Step 2: Align ─────────────────────────────────────────────────────────
    echo "  [2/4] Aligning with bwa-mem2..."
    BAM="$SDIR/align/aligned.bam"
    bwa-mem2 mem -t "$THREADS" \
        -R "@RG\tID:1\tSM:SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
        "$REF" "$R1" "$R2" 2>"$SDIR/align/bwa.log" \
    | samtools sort -@ "$THREADS" -o "$SDIR/align/sorted.bam"

    samtools markdup -@ "$THREADS" \
        "$SDIR/align/sorted.bam" "$BAM" \
        2>"$SDIR/align/markdup.log"

    samtools index -@ "$THREADS" "$BAM"
    echo "        Done. BAM: $BAM"

    # ── Step 3: Call variants with Mutect2 ────────────────────────────────────
    echo "  [3/4] Calling variants with Mutect2..."
    RAW_VCF="$SDIR/call/raw.vcf.gz"
    FILT_VCF="$SDIR/call/filtered.vcf.gz"

    gatk Mutect2 \
        -R "$REF" \
        -I "$BAM" \
        -tumor SAMPLE \
        --output "$RAW_VCF" \
        --native-pair-hmm-threads "$THREADS" \
        2>"$SDIR/call/mutect2.log"

    gatk FilterMutectCalls \
        -R "$REF" \
        -V "$RAW_VCF" \
        --output "$FILT_VCF" \
        2>"$SDIR/call/filter.log"

    echo "        Done. Filtered VCF: $FILT_VCF"

    # ── Step 4: Evaluate ──────────────────────────────────────────────────────
    echo "  [4/4] Evaluating against truth VCF..."
    rtg vcfeval \
        -b "$TRUTH" \
        -c "$FILT_VCF" \
        -t "$SDF" \
        -o "$SDIR/eval/rtg" \
        --output-mode annotate \
        2>"$SDIR/eval/rtg.log" || true  # rtg exits non-zero even on partial success

    python3 "$SCRIPT_DIR/evaluate.py" \
        --truth "$TRUTH" \
        --calls "$FILT_VCF" \
        --rtg-eval "$SDIR/eval/rtg" \
        --scenario "$NAME" \
        --outdir "$SDIR/eval"

    echo "        Done. Results: $SDIR/eval/summary.tsv"
done

# ── Aggregate summary ─────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════"
echo "  Aggregating results"
echo "══════════════════════════════════════════════"
python3 "$SCRIPT_DIR/evaluate.py" \
    --aggregate "$OUTDIR" \
    --outdir "$OUTDIR"

echo ""
echo "All done. Full summary: $OUTDIR/summary.tsv"
