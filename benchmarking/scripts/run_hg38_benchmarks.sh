#!/usr/bin/env bash
# VarForge hg38 (chr22) Benchmark Runner
# Runs from the repo root: bash benchmarking/scripts/run_hg38_benchmarks.sh
#
# Requires:
#   - ./target/release/varforge (release build)
#   - ./benchmarking/chr22.fa   (not committed to git)
#
# Outputs:
#   - benchmarking/results/tables/hg38_benchmarks.csv
set -euo pipefail

VARFORGE="./target/release/varforge"
CONFIGS_DIR="./benchmarking/scripts/configs_hg38"
TIMING_DIR="/tmp/hg38_timing"
OUT_DIR="benchmarking/results/tables"

mkdir -p "${TIMING_DIR}" "${OUT_DIR}"

# Verify binary exists
if [[ ! -x "${VARFORGE}" ]]; then
    echo "Building varforge..."
    cargo build --release --quiet
fi

# Verify reference exists
if [[ ! -f "./benchmarking/chr22.fa" ]]; then
    echo "ERROR: benchmarking/chr22.fa not found" >&2
    echo "Download with:" >&2
    echo "  wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz | gunzip > benchmarking/chr22.fa" >&2
    exit 1
fi

echo "============================================================"
echo "VarForge hg38 (chr22) Benchmark"
echo "Date: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "CPU:  $(grep 'model name' /proc/cpuinfo | head -1 | cut -d: -f2 | xargs)"
echo "RAM:  $(free -h | awk '/Mem:/{print $2}')"
echo "Ref:  $(du -h ./benchmarking/chr22.fa | cut -f1) (chr22.fa)"
echo "============================================================"
echo ""

# Write CSV header
CSV="${OUT_DIR}/hg38_benchmarks.csv"
echo "config,wall_secs,peak_mem_kb,read_pairs,r1_gz_bytes" > "${CSV}"

# Python runner: times the simulation and polls /proc for peak RSS
python3 - "${VARFORGE}" "${CONFIGS_DIR}" "${TIMING_DIR}" "${OUT_DIR}" "${CSV}" "$(pwd)" <<'PYEOF'
import subprocess, sys, time, os, re, json, glob, shutil

varforge, configs_dir, timing_dir, out_dir, csv_path, cwd = sys.argv[1:]

configs = [
    "01_wgs_baseline",
    "02_wgs_variants",
    "03_panel_umi",
    "04_twist_duplex",
    "05_cfdna",
    "06_ffpe_tumour",
]

def peak_rss_kb(pid):
    """Read current RSS from /proc/{pid}/status."""
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except (FileNotFoundError, ValueError):
        pass
    return 0

def run_config(name, seed=42):
    print(f"  Running {name}...", flush=True)
    orig = os.path.join(configs_dir, f"{name}.yaml")
    tmp_cfg = os.path.join(timing_dir, f"{name}.yaml")
    tmp_out = f"/tmp/hg38_run_{name}"

    # Patch reference path to absolute
    with open(orig) as f:
        content = f.read()
    content = content.replace("reference: chr22.fa",
                              f"reference: {cwd}/benchmarking/chr22.fa")
    with open(tmp_cfg, "w") as f:
        f.write(content)

    shutil.rmtree(tmp_out, ignore_errors=True)
    os.makedirs(tmp_out, exist_ok=True)

    cmd = [varforge, "simulate",
           "--config", tmp_cfg,
           "--output-dir", tmp_out,
           "--seed", str(seed)]

    t0 = time.monotonic()
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    peak_kb = 0
    while proc.poll() is None:
        rss = peak_rss_kb(proc.pid)
        if rss > peak_kb:
            peak_kb = rss
        time.sleep(0.5)

    elapsed = time.monotonic() - t0
    rc = proc.returncode

    if rc != 0:
        print(f"    FAILED (exit {rc})", flush=True)
        return f"{name},ERROR,ERROR,ERROR,ERROR"

    wall = round(elapsed, 2)
    mem  = peak_kb

    # Read pairs from manifest
    pairs = 0
    for mf in glob.glob(os.path.join(tmp_out, "*.manifest.json")):
        try:
            with open(mf) as f:
                d = json.load(f)
            pairs = d.get("total_read_pairs", 0)
            break
        except Exception:
            pass

    # R1 file size
    r1gz = 0
    for gz in glob.glob(os.path.join(tmp_out, "*_R1.fastq.gz")):
        r1gz = os.path.getsize(gz)
        break

    print(f"    wall={wall}s  mem={mem}KB  pairs={pairs}  r1={r1gz}B", flush=True)
    shutil.rmtree(tmp_out, ignore_errors=True)
    os.remove(tmp_cfg)
    return f"{name},{wall},{mem},{pairs},{r1gz}"

with open(csv_path, "a") as out:
    for cfg in configs:
        row = run_config(cfg)
        out.write(row + "\n")
        out.flush()

print(f"\nResults written to {csv_path}")
PYEOF

echo ""
echo "============================================================"
cat "${CSV}"
echo "============================================================"
