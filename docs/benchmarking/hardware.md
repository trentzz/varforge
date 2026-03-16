# Benchmark Hardware

All VarForge benchmarks were executed on the following hardware configuration.

## System Specifications

| Component     | Details                                                     |
|---------------|-------------------------------------------------------------|
| CPU           | AMD Ryzen 5 7640HS w/ Radeon 760M Graphics                 |
| Cores / Threads | 6 cores / 12 threads                                      |
| Max Clock     | 5.0 GHz (boost)                                             |
| RAM           | 12 GiB DDR5                                                 |
| Storage       | 454 GB SSD (NVMe)                                           |
| OS            | Debian GNU/Linux 13 (trixie), kernel 6.12.74+deb13+1-amd64 |
| Rust toolchain | 1.94.0 (stable)                                            |
| VarForge      | 0.1.0 (`--release` build, no PGO)                          |

## Date

Benchmarks collected: 2026-03-16

## Notes

- The system has 12 GiB physical RAM. Several high-coverage scenarios (configs 05, 06, 09)
  exceeded physical memory and required swap, resulting in degraded throughput and inflated
  peak RSS figures reported by `time -v`. These results are annotated in the results file.
- The CPU is a mobile APU running on AC power. Boost behaviour is consistent across runs
  because the workload is sustained-compute rather than bursty.
- No other significant processes were running during benchmarks.
