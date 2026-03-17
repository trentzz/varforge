# Benchmark Hardware

**Date:** 2026-03-17

| Component | Specification |
|-----------|---------------|
| CPU | AMD Ryzen 5 7640HS (6 cores / 12 threads, 5.0 GHz boost) |
| RAM | 12 GiB DDR5 |
| Storage | 454 GB NVMe SSD |
| OS | Debian GNU/Linux 13, kernel 6.12.74+deb13+1-amd64 |
| Rust | 1.94.0 stable (2026-03-02) |
| VarForge | 0.1.0, `--release` mode (LTO thin, symbols stripped) |

## Notes

- High-coverage UMI configurations exceed 12 GiB physical memory and require swap. Those results are reported separately as memory projections.
- The NVMe SSD bandwidth limits I/O-bound workloads (~280 MB compressed FASTQ at 10 MB reference / 30× coverage saturates the write pipeline).
- All benchmarks use the pure-Rust `flate2` backend for gzip compression (no system zlib).
- No other significant processes were running during benchmarks.
