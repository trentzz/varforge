#!/usr/bin/env python3
"""
Compile VarForge benchmark timing results into results.json.
Run from 
"""
import json, os, re, glob, subprocess

bench_dir = "./benchmark_output"
timing_dir = os.path.join(bench_dir, "timing")


def parse_timing_file(path):
    """Parse busybox/GNU time -v output file."""
    result = {}
    if not os.path.exists(path):
        return result
    with open(path) as f:
        content = f.read()

    # Wall clock time: "1m 15.53s" or "0m 37.73s"
    m = re.search(r'(\d+)m\s+([\d.]+)s', content)
    if m:
        result['wall_secs'] = round(int(m.group(1)) * 60 + float(m.group(2)), 2)
    else:
        # Try h:mm:ss.cs format
        m = re.search(r'(\d+):(\d+):([\d.]+)', content)
        if m:
            result['wall_secs'] = round(int(m.group(1))*3600 + int(m.group(2))*60 + float(m.group(3)), 2)
        else:
            m = re.search(r'(\d+):([\d.]+)', content)
            if m:
                result['wall_secs'] = round(int(m.group(1)) * 60 + float(m.group(2)), 2)

    # Peak memory
    m = re.search(r'Maximum resident set size.*?: (\d+)', content)
    if m:
        result['peak_mem_kb'] = int(m.group(1))

    # CPU usage
    m = re.search(r'Percent of CPU this job got: (\d+)%', content)
    if m:
        result['cpu_pct'] = int(m.group(1))

    # User + system time
    m = re.search(r'User time \(seconds\): ([\d.]+)', content)
    if m:
        result['user_time_secs'] = float(m.group(1))
    m = re.search(r'System time \(seconds\): ([\d.]+)', content)
    if m:
        result['system_time_secs'] = float(m.group(1))

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


try:
    ver = subprocess.check_output(["./target/release/varforge", "--version"],
                                   stderr=subprocess.STDOUT).decode().strip()
except Exception:
    ver = "unknown"

cpu_info = "unknown"
if os.path.exists("/proc/cpuinfo"):
    with open("/proc/cpuinfo") as f:
        for line in f:
            if "model name" in line:
                cpu_info = line.split(":")[1].strip()
                break

results = {
    "metadata": {
        "date": "2026-03-16",
        "varforge_version": ver,
        "cpu": cpu_info,
        "cores": os.cpu_count(),
        "platform": "linux",
        "notes": [
            "High-coverage configs (03,05,06,09,10,12) ran 1 iteration due to memory/time constraints (12GB RAM system)",
            "04_very_high_coverage: reduced from 10MB/500x to 1MB/200x",
            "07_cfdna: reference reduced from 10MB to 1MB",
            "12_ultra_deep: reduced from 2000x to 500x",
            "Coverage scaling experiment uses 1MB reference for memory feasibility",
            "Feature overhead experiment uses 1MB reference at 30x",
        ],
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

# Feature overhead (benchmark script writes timing files as feat_feat_<name>_iterN.txt)
for feat in ["feat_base", "feat_variants", "feat_umi", "feat_artifacts", "feat_gc_bias", "feat_all"]:
    # Try both naming conventions
    data = parse_group(f"feat_{feat}")
    if data.get("n_iterations", 0) == 0:
        data = parse_group(feat)
    results["feature_overhead"][feat] = data

output_path = os.path.join(bench_dir, "results.json")
with open(output_path, "w") as f:
    json.dump(results, f, indent=2)
print(f"Results written to {output_path}")

# Print summary table
print("\n=== SUMMARY ===")
print(f"\nMain Benchmarks:")
print(f"{'Config':<35} {'Iters':>5} {'Avg Wall (s)':>14} {'Min Wall (s)':>14} {'Avg Mem (MB)':>14}")
print("-" * 86)
for name, data in results["main_benchmarks"].items():
    n = data.get("n_iterations", 0)
    avg_w = data.get("avg_wall_secs")
    min_w = data.get("min_wall_secs")
    avg_m = data.get("avg_peak_mem_kb")
    avg_w_str = f"{avg_w:.1f}" if avg_w else "N/A"
    min_w_str = f"{min_w:.1f}" if min_w else "N/A"
    avg_m_str = f"{avg_m/1024:.0f}" if avg_m else "N/A"
    print(f"{name:<35} {n:>5} {avg_w_str:>14} {min_w_str:>14} {avg_m_str:>14}")

print(f"\nThread Scaling (02_with_variants, 10MB 30x):")
print(f"{'Threads':>8} {'Iters':>5} {'Avg Wall (s)':>14} {'Speedup':>10}")
print("-" * 42)
base_wall = None
for t in [1, 2, 4, 6, 8, 12]:
    key = f"threads_{t}"
    data = results["thread_scaling"][key]
    avg_w = data.get("avg_wall_secs")
    n = data.get("n_iterations", 0)
    if t == 1 and avg_w:
        base_wall = avg_w
    speedup = f"{base_wall/avg_w:.2f}x" if (base_wall and avg_w) else "N/A"
    avg_w_str = f"{avg_w:.1f}" if avg_w else "N/A"
    print(f"{t:>8} {n:>5} {avg_w_str:>14} {speedup:>10}")

print(f"\nCoverage Scaling (1MB ref):")
print(f"{'Coverage':>10} {'Iters':>5} {'Avg Wall (s)':>14} {'Avg Mem (MB)':>14}")
print("-" * 48)
for cov in [1, 5, 10, 30, 50, 100, 200]:
    key = f"cov_{cov}x"
    data = results["coverage_scaling"][key]
    avg_w = data.get("avg_wall_secs")
    avg_m = data.get("avg_peak_mem_kb")
    n = data.get("n_iterations", 0)
    avg_w_str = f"{avg_w:.1f}" if avg_w else "N/A"
    avg_m_str = f"{avg_m/1024:.0f}" if avg_m else "N/A"
    print(f"{cov:>9}x {n:>5} {avg_w_str:>14} {avg_m_str:>14}")

print(f"\nFeature Overhead (1MB ref, 30x):")
print(f"{'Feature':<20} {'Iters':>5} {'Avg Wall (s)':>14} {'Overhead':>12}")
print("-" * 56)
base_feat_wall = None
for feat in ["feat_base", "feat_variants", "feat_umi", "feat_artifacts", "feat_gc_bias", "feat_all"]:
    data = results["feature_overhead"].get(feat, {})
    avg_w = data.get("avg_wall_secs")
    n = data.get("n_iterations", 0)
    if feat == "feat_base" and avg_w:
        base_feat_wall = avg_w
    overhead = f"+{(avg_w - base_feat_wall):.1f}s ({(avg_w/base_feat_wall - 1)*100:.0f}%)" if (base_feat_wall and avg_w and feat != "feat_base") else "baseline"
    avg_w_str = f"{avg_w:.1f}" if avg_w else "N/A"
    print(f"{feat:<20} {n:>5} {avg_w_str:>14} {overhead:>12}")
