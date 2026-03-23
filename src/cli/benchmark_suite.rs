//! `benchmark-suite` subcommand: generate a VAF × coverage grid of simulations.
//!
//! Each combination of VAF and coverage is run as an independent single-sample
//! simulation.  Output goes into subdirectories named `vaf{VAF}_cov{COV}x/`
//! under the base output directory.  A TSV manifest summarising all runs is
//! written to `{base}/benchmark_manifest.tsv`.

use anyhow::{Context, Result};

use super::BenchmarkSuiteOpts;
use crate::io::config;
use crate::io::reference::ReferenceGenome;

/// Run the VAF × coverage benchmark grid.
pub fn run(opts: BenchmarkSuiteOpts, cli_threads: Option<usize>) -> Result<()> {
    let mut base_cfg = config::load(&opts.config)
        .with_context(|| format!("failed to load config: {}", opts.config.display()))?;

    if let Some(ref out) = opts.output {
        base_cfg.output.directory = out.clone();
    }

    anyhow::ensure!(
        !opts.vafs.is_empty(),
        "at least one --vafs value is required"
    );
    anyhow::ensure!(
        !opts.coverages.is_empty(),
        "at least one --coverages value is required"
    );

    let base_output = base_cfg.output.directory.clone();
    std::fs::create_dir_all(&base_output).with_context(|| {
        format!(
            "failed to create output directory: {}",
            base_output.display()
        )
    })?;

    // Build the (vaf, coverage) grid.
    let grid: Vec<(f64, f64)> = opts
        .vafs
        .iter()
        .flat_map(|&v| opts.coverages.iter().map(move |&c| (v, c)))
        .collect();

    eprintln!(
        "Running {}×{} benchmark grid ({} total runs)",
        opts.vafs.len(),
        opts.coverages.len(),
        grid.len()
    );

    let mut manifest_rows: Vec<String> = Vec::new();
    manifest_rows.push("vaf\tcoverage\toutput_dir\tseed\tstatus".to_string());

    for (vaf, coverage) in &grid {
        let mut run_cfg = base_cfg.clone();

        // Override VAF in random mutations config.
        if let Some(ref mut mutations) = run_cfg.mutations {
            if let Some(ref mut random) = mutations.random {
                random.vaf_min = *vaf;
                // Add a tiny epsilon so vaf_min < vaf_max, as the validator requires.
                random.vaf_max = vaf + 1e-9;
            }
        }

        // Override coverage.
        run_cfg.sample.coverage = *coverage;

        // Set output subdirectory: {base}/vaf{VAF}_cov{COV}x/
        let vaf_label = format!("{:.4}", vaf)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();
        let cov_label = format!("{:.0}", coverage);
        let run_dir = base_output.join(format!("vaf{}_cov{}x", vaf_label, cov_label));
        run_cfg.output.directory = run_dir.clone();

        // Apply CLI thread pool config.
        if let Some(t) = cli_threads {
            run_cfg.threads = Some(t);
        }

        eprintln!(
            "  VAF={:.4} coverage={}x -> {}",
            vaf,
            cov_label,
            run_dir.display()
        );

        let seed = run_cfg.seed;

        let reference = ReferenceGenome::open(&run_cfg.reference).with_context(|| {
            format!("failed to open reference: {}", run_cfg.reference.display())
        })?;

        let status = match crate::cli::simulate::run_single_sample(
            run_cfg,
            reference,
            std::time::Instant::now(),
            false,
        ) {
            Ok(()) => "ok".to_string(),
            Err(e) => format!("err:{}", e),
        };

        manifest_rows.push(format!(
            "{:.6}\t{:.1}\t{}\t{}\t{}",
            vaf,
            coverage,
            run_dir.display(),
            seed.map(|s| s.to_string())
                .unwrap_or_else(|| "-".to_string()),
            status
        ));
    }

    // Write TSV manifest summarising all runs.
    let manifest_path = base_output.join("benchmark_manifest.tsv");
    std::fs::write(&manifest_path, manifest_rows.join("\n") + "\n")
        .with_context(|| format!("failed to write manifest: {}", manifest_path.display()))?;

    eprintln!("Benchmark manifest written to {}", manifest_path.display());
    Ok(())
}
