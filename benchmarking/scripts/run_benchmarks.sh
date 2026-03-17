#!/usr/bin/env bash
# VarForge Comprehensive Benchmark Runner
# Run from: 
set -euo pipefail

VARFORGE="./target/release/varforge"
BENCH_DIR="./benchmark_output"
CONFIGS_DIR="${BENCH_DIR}/configs"
TIMING_DIR="${BENCH_DIR}/timing"
RUNS_DIR="${BENCH_DIR}/runs"

# Use busybox time if GNU time is not available
if command -v /usr/bin/time &>/dev/null; then
    TIME_BIN="/usr/bin/time"
elif [ -x "/tmp/time" ]; then
    TIME_BIN="/tmp/time"
elif [ -x "/mnt/tzeng-local/tzeng-thesis/fastp_0.23.4_dir/usr/bin/time" ]; then
    # Copy to /tmp and make symlink named "time" (busybox needs applet name)
    cp /mnt/tzeng-local/tzeng-thesis/fastp_0.23.4_dir/usr/bin/time /tmp/busybox_time_bin
    ln -sf /tmp/busybox_time_bin /tmp/time
    TIME_BIN="/tmp/time"
else
    echo "ERROR: no suitable 'time' binary found" >&2
    exit 1
fi
echo "Using time binary: ${TIME_BIN}"

mkdir -p "${TIMING_DIR}" "${RUNS_DIR}"

# Helper: run one simulation with timing via -o flag
# Usage: run_sim CONFIG_PATH OUTPUT_DIR TIMING_FILE [EXTRA_ARGS...]
run_sim() {
    local config="$1"
    local out_dir="$2"
    local timing_file="$3"
    shift 3
    local extra_args=("$@")

    rm -rf "${out_dir}"
    mkdir -p "$(dirname "${out_dir}")"

    if "${TIME_BIN}" -v -o "${timing_file}" "${VARFORGE}" simulate \
        --config "${config}" \
        --output-dir "${out_dir}" \
        "${extra_args[@]}" ; then
        return 0
    else
        echo "  ERROR: simulation failed for ${config}" >&2
        return 1
    fi
}

# Helper: extract wall clock seconds from time -v output
wall_secs_from_timing() {
    local f="$1"
    grep "Elapsed (wall clock)" "${f}" | sed 's/.*: //' | awk -F'[ m:s]' '
    {
        # Format is like "1m 14.85s" or "0m 37.73s"
        for(i=1; i<=NF; i++) {
            if ($i ~ /m$/) { mins = substr($i, 1, length($i)-1) }
            if ($i ~ /s$/) { secs = substr($i, 1, length($i)-1) }
        }
        print mins*60 + secs
    }' 2>/dev/null || \
    grep "Elapsed (wall clock)" "${f}" | sed 's/.*: //' | awk -F'[: ]' '{
        if (NF==3) print $1*60 + $2
        else if (NF==4) print $1*3600 + $2*60 + $3
        else print 0
    }'
}

peak_mem_kb_from_timing() {
    local f="$1"
    grep "Maximum resident" "${f}" | awk '{print $NF}'
}

echo "============================================================"
echo "VarForge Benchmark Suite"
echo "Date: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "CPU: $(grep 'model name' /proc/cpuinfo | head -1 | cut -d: -f2 | xargs)"
echo "Cores: $(nproc)"
echo "RAM: $(free -h | awk '/Mem:/{print $2}')"
echo "VarForge: $(${VARFORGE} --version 2>&1)"
echo "============================================================"
echo ""

# ============================================================
# STEP 3: Main benchmark runs (all 12 configs)
# Configs that take >2 min get 1 iteration; others get 3.
# ============================================================
echo "=== STEP 3: Main Benchmark Runs ==="

# Format: "config_name:iterations"
# 03,05,06,09 are 10MB 100x = ~4min each -> 1 iter
# 10 is 50MB 30x = ~6min each -> 1 iter
# 12 is 1MB 500x = ~3min each -> 1 iter
declare -A CONFIG_ITERS=(
    ["01_baseline"]=3
    ["02_with_variants"]=3
    ["03_high_coverage"]=1
    ["04_very_high_coverage"]=3
    ["05_umi_simplex"]=1
    ["06_umi_duplex"]=1
    ["07_cfdna"]=3
    ["08_ffpe_artifacts"]=3
    ["09_all_features"]=1
    ["10_large_genome"]=1
    ["11_panel_simulation"]=3
    ["12_ultra_deep"]=1
)

for config in "${CONFIGS_DIR}"/*.yaml; do
    cfgname=$(basename "${config}" .yaml)
    iters="${CONFIG_ITERS[${cfgname}]:-3}"
    echo ""
    echo "--- Config: ${cfgname} (${iters} iterations) ---"

    for iter in $(seq 1 "${iters}"); do
        out_dir="${RUNS_DIR}/${cfgname}_iter${iter}"
        timing_file="${TIMING_DIR}/${cfgname}_iter${iter}.txt"

        echo -n "  Iteration ${iter}... "

        if run_sim "${config}" "${out_dir}" "${timing_file}" --seed "${iter}"; then
            wall=$(wall_secs_from_timing "${timing_file}")
            mem=$(peak_mem_kb_from_timing "${timing_file}")

            # Count read pairs from manifest
            read_pairs=$(grep -o '"total_read_pairs":[0-9]*' "${out_dir}"/*.manifest.json 2>/dev/null | head -1 | cut -d: -f2 || echo "N/A")
            r1_size=$(du -b "${out_dir}"/*_R1.fastq.gz 2>/dev/null | awk '{print $1}' | head -1 || echo "0")

            echo "wall=${wall}s mem=${mem}KB reads=${read_pairs} r1=${r1_size}B"
            rm -rf "${out_dir}"
        else
            echo "FAILED"
        fi
    done
done

# ============================================================
# STEP 4: Thread Scaling
# Using 02_with_variants config (10MB 30x = ~75s per run)
# 6 thread counts x 3 iterations = 18 runs x ~75s = ~23min
# ============================================================
echo ""
echo "=== STEP 4: Thread Scaling ==="
THREAD_COUNTS=(1 2 4 6 8 12)
THREAD_CONFIG="${CONFIGS_DIR}/02_with_variants.yaml"

for threads in "${THREAD_COUNTS[@]}"; do
    echo ""
    echo "--- Threads: ${threads} ---"
    for iter in 1 2 3; do
        out_dir="${RUNS_DIR}/thread_scale_t${threads}_iter${iter}"
        timing_file="${TIMING_DIR}/thread_scale_t${threads}_iter${iter}.txt"

        echo -n "  Iter ${iter}... "

        if "${TIME_BIN}" -v -o "${timing_file}" "${VARFORGE}" simulate \
            --config "${THREAD_CONFIG}" \
            --output-dir "${out_dir}" \
            --seed "${iter}" \
            -t "${threads}" ; then

            wall=$(wall_secs_from_timing "${timing_file}")
            mem=$(peak_mem_kb_from_timing "${timing_file}")
            echo "wall=${wall}s mem=${mem}KB"
            rm -rf "${out_dir}"
        else
            echo "FAILED"
        fi
    done
done

# ============================================================
# STEP 5: Coverage Scaling
# Use 1MB ref to keep memory manageable: ~20s per run max at 200x
# 7 coverages x 3 iterations = 21 runs x ~20s avg = ~7min
# ============================================================
echo ""
echo "=== STEP 5: Coverage Scaling ==="
COVERAGES=(1 5 10 30 50 100 200)

# Create a coverage-scaling config based on 1MB ref with variants
COV_CONFIG_DIR="${BENCH_DIR}/configs_cov"
mkdir -p "${COV_CONFIG_DIR}"
cat > "${COV_CONFIG_DIR}/cov_scale_base.yaml" << COVYAML
reference: benchmark_output/ref_1mb.fa
output:
  directory: /tmp/cov_scale_placeholder
  fastq: true
  bam: false
  truth_vcf: false
  manifest: false
sample:
  name: COV_SCALE
  read_length: 150
  coverage: 30.0
fragment:
  model: normal
  mean: 300.0
  sd: 50.0
mutations:
  random:
    count: 100
    vaf_min: 0.05
    vaf_max: 0.50
    snv_fraction: 0.80
    indel_fraction: 0.15
    mnv_fraction: 0.05
seed: 1
COVYAML

for cov in "${COVERAGES[@]}"; do
    echo ""
    echo "--- Coverage: ${cov}x ---"
    for iter in 1 2 3; do
        out_dir="${RUNS_DIR}/cov_scale_${cov}x_iter${iter}"
        timing_file="${TIMING_DIR}/cov_scale_${cov}x_iter${iter}.txt"

        echo -n "  Iter ${iter}... "

        if "${TIME_BIN}" -v -o "${timing_file}" "${VARFORGE}" simulate \
            --config "${COV_CONFIG_DIR}/cov_scale_base.yaml" \
            --output-dir "${out_dir}" \
            --seed "${iter}" \
            --coverage "${cov}" ; then

            wall=$(wall_secs_from_timing "${timing_file}")
            mem=$(peak_mem_kb_from_timing "${timing_file}")
            r1_size=$(du -b "${out_dir}"/*_R1.fastq.gz 2>/dev/null | awk '{print $1}' | head -1 || echo "0")
            echo "wall=${wall}s mem=${mem}KB r1=${r1_size}B"
            rm -rf "${out_dir}"
        else
            echo "FAILED"
        fi
    done
done

# ============================================================
# STEP 6: Feature Overhead
# All use 1MB ref at 30x for fast comparisons (~5s each)
# 6 features x 3 iterations = 18 runs x ~5s = ~90s total
# ============================================================
echo ""
echo "=== STEP 6: Feature Overhead ==="

FEAT_CONFIGS_DIR="${BENCH_DIR}/configs_features"
mkdir -p "${FEAT_CONFIGS_DIR}"
REF1MB="benchmark_output/ref_1mb.fa"

# Base (no features)
cat > "${FEAT_CONFIGS_DIR}/feat_base.yaml" << YAML
reference: ${REF1MB}
output:
  directory: /tmp/feat_placeholder
  fastq: true
  bam: false
  truth_vcf: false
  manifest: false
sample:
  name: FEAT_BASE
  read_length: 150
  coverage: 30.0
fragment:
  model: normal
  mean: 300.0
  sd: 50.0
seed: 42
YAML

# +variants
cat > "${FEAT_CONFIGS_DIR}/feat_variants.yaml" << YAML
reference: ${REF1MB}
output:
  directory: /tmp/feat_placeholder
  fastq: true
  bam: false
  truth_vcf: false
  manifest: false
sample:
  name: FEAT_VARIANTS
  read_length: 150
  coverage: 30.0
fragment:
  model: normal
  mean: 300.0
  sd: 50.0
mutations:
  random:
    count: 500
    vaf_min: 0.05
    vaf_max: 0.50
    snv_fraction: 0.80
    indel_fraction: 0.15
    mnv_fraction: 0.05
seed: 42
YAML

# +UMI
cat > "${FEAT_CONFIGS_DIR}/feat_umi.yaml" << YAML
reference: ${REF1MB}
output:
  directory: /tmp/feat_placeholder
  fastq: true
  bam: false
  truth_vcf: false
  manifest: false
sample:
  name: FEAT_UMI
  read_length: 150
  coverage: 30.0
fragment:
  model: normal
  mean: 300.0
  sd: 50.0
umi:
  length: 8
  duplex: false
  pcr_cycles: 10
  family_size_mean: 3.0
  family_size_sd: 1.5
  inline: true
seed: 42
YAML

# +artifacts
cat > "${FEAT_CONFIGS_DIR}/feat_artifacts.yaml" << YAML
reference: ${REF1MB}
output:
  directory: /tmp/feat_placeholder
  fastq: true
  bam: false
  truth_vcf: false
  manifest: false
sample:
  name: FEAT_ARTIFACTS
  read_length: 150
  coverage: 30.0
fragment:
  model: normal
  mean: 300.0
  sd: 50.0
artifacts:
  ffpe_damage_rate: 0.02
  oxog_rate: 0.01
  duplicate_rate: 0.15
  pcr_error_rate: 0.001
seed: 42
YAML

# +GC bias
cat > "${FEAT_CONFIGS_DIR}/feat_gc_bias.yaml" << YAML
reference: ${REF1MB}
output:
  directory: /tmp/feat_placeholder
  fastq: true
  bam: false
  truth_vcf: false
  manifest: false
sample:
  name: FEAT_GC_BIAS
  read_length: 150
  coverage: 30.0
fragment:
  model: normal
  mean: 300.0
  sd: 50.0
gc_bias:
  enabled: true
  model: default
  severity: 1.0
seed: 42
YAML

# all combined
cat > "${FEAT_CONFIGS_DIR}/feat_all.yaml" << YAML
reference: ${REF1MB}
output:
  directory: /tmp/feat_placeholder
  fastq: true
  bam: false
  truth_vcf: false
  manifest: false
sample:
  name: FEAT_ALL
  read_length: 150
  coverage: 30.0
fragment:
  model: normal
  mean: 300.0
  sd: 50.0
mutations:
  random:
    count: 500
    vaf_min: 0.05
    vaf_max: 0.50
    snv_fraction: 0.80
    indel_fraction: 0.15
    mnv_fraction: 0.05
umi:
  length: 8
  duplex: false
  pcr_cycles: 10
  family_size_mean: 3.0
  family_size_sd: 1.5
  inline: true
artifacts:
  ffpe_damage_rate: 0.02
  oxog_rate: 0.01
  duplicate_rate: 0.15
  pcr_error_rate: 0.001
gc_bias:
  enabled: true
  model: default
  severity: 1.0
seed: 42
YAML

for feat_cfg in feat_base feat_variants feat_umi feat_artifacts feat_gc_bias feat_all; do
    echo ""
    echo "--- Feature: ${feat_cfg} ---"
    for iter in 1 2 3; do
        out_dir="${RUNS_DIR}/feat_${feat_cfg}_iter${iter}"
        timing_file="${TIMING_DIR}/feat_${feat_cfg}_iter${iter}.txt"

        echo -n "  Iter ${iter}... "

        if "${TIME_BIN}" -v -o "${timing_file}" "${VARFORGE}" simulate \
            --config "${FEAT_CONFIGS_DIR}/${feat_cfg}.yaml" \
            --output-dir "${out_dir}" \
            --seed "${iter}" ; then

            wall=$(wall_secs_from_timing "${timing_file}")
            mem=$(peak_mem_kb_from_timing "${timing_file}")
            echo "wall=${wall}s mem=${mem}KB"
            rm -rf "${out_dir}"
        else
            echo "FAILED"
        fi
    done
done

# ============================================================
# STEP 7: Compile results to JSON
# ============================================================
echo ""
echo "=== STEP 7: Compiling Results ==="

python3 - <<'PYEOF'
import json, os, re, glob

bench_dir = "./benchmark_output"
timing_dir = os.path.join(bench_dir, "timing")

def parse_timing_file(path):
    """Parse time -v output file."""
    result = {}
    if not os.path.exists(path):
        return result
    with open(path) as f:
        content = f.read()

    # Wall clock time - formats: "1m 14.85s" or "0m 37.73s" or "1:14.85" or "0:37.73"
    m = re.search(r'Elapsed \(wall clock\) time[^:]*: (\d+)m\s+([\d.]+)s', content)
    if m:
        result['wall_secs'] = int(m.group(1)) * 60 + float(m.group(2))
    else:
        m = re.search(r'Elapsed \(wall clock\) time[^:]*: (\d+):(\d+):(\d+)\.(\d+)', content)
        if m:
            result['wall_secs'] = int(m.group(1))*3600 + int(m.group(2))*60 + int(m.group(3)) + int(m.group(4))/100
        else:
            m = re.search(r'Elapsed \(wall clock\) time[^:]*: (\d+):([\d.]+)', content)
            if m:
                result['wall_secs'] = int(m.group(1)) * 60 + float(m.group(2))

    # Peak memory
    m = re.search(r'Maximum resident set size.*?: (\d+)', content)
    if m:
        result['peak_mem_kb'] = int(m.group(1))

    # CPU usage
    m = re.search(r'Percent of CPU this job got: (\d+)%', content)
    if m:
        result['cpu_pct'] = int(m.group(1))

    return result

def avg(vals):
    vals = [v for v in vals if v is not None]
    return round(sum(vals) / len(vals), 4) if vals else None

def parse_group(prefix, max_iterations=3):
    runs = []
    for i in range(1, max_iterations + 1):
        fname = os.path.join(timing_dir, f"{prefix}_iter{i}.txt")
        if os.path.exists(fname):
            t = parse_timing_file(fname)
            if t:
                runs.append(t)

    walls = [r.get('wall_secs') for r in runs]
    mems = [r.get('peak_mem_kb') for r in runs]
    return {
        "iterations": runs,
        "n_iterations": len(runs),
        "avg_wall_secs": avg(walls),
        "avg_peak_mem_kb": avg(mems),
        "min_wall_secs": min((v for v in walls if v is not None), default=None),
        "max_wall_secs": max((v for v in walls if v is not None), default=None),
    }

import subprocess
try:
    ver = subprocess.check_output(["./target/release/varforge", "--version"],
                                   stderr=subprocess.STDOUT).decode().strip()
except:
    ver = "unknown"

results = {
    "metadata": {
        "date": "2026-03-16",
        "varforge_version": ver,
        "cpu": open("/proc/cpuinfo").read().split("model name")[1].split("\n")[0].strip(": ") if os.path.exists("/proc/cpuinfo") else "unknown",
        "cores": os.cpu_count(),
        "notes": "High-coverage configs (03,05,06,09,10,12) ran 1 iteration due to memory/time constraints on 12GB RAM system",
    },
    "main_benchmarks": {},
    "thread_scaling": {},
    "coverage_scaling": {},
    "feature_overhead": {},
}

# Main benchmarks
configs = sorted(glob.glob(os.path.join(bench_dir, "configs", "*.yaml")))
for cfg in configs:
    name = os.path.basename(cfg).replace(".yaml", "")
    results["main_benchmarks"][name] = parse_group(name)

# Thread scaling
for t in [1, 2, 4, 6, 8, 12]:
    key = f"threads_{t}"
    results["thread_scaling"][key] = parse_group(f"thread_scale_t{t}")

# Coverage scaling
for cov in [1, 5, 10, 30, 50, 100, 200]:
    key = f"cov_{cov}x"
    results["coverage_scaling"][key] = parse_group(f"cov_scale_{cov}x")

# Feature overhead
for feat in ["feat_base", "feat_variants", "feat_umi", "feat_artifacts", "feat_gc_bias", "feat_all"]:
    results["feature_overhead"][feat] = parse_group(feat)

output_path = os.path.join(bench_dir, "results.json")
with open(output_path, "w") as f:
    json.dump(results, f, indent=2)
print(f"Results written to {output_path}")
PYEOF

echo ""
echo "============================================================"
echo "Benchmark suite complete!"
echo "Results: ${BENCH_DIR}/results.json"
echo "Timing data: ${TIMING_DIR}/"
echo "============================================================"
