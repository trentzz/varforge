//! `varforge learn-profile` subcommand.
//!
//! Analyzes a real BAM file and emits an error/quality profile that can
//! be fed back into VarForge via `quality.profile_path`.
//! Output format is either JSON (default) or TSV.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use tracing::info;

use crate::core::profile_learner::{
    read_quality_tsv, write_quality_tsv, LearnerConfig, ProfileLearner,
};

/// Output format for the learned quality profile.
#[derive(Debug, Clone, PartialEq, clap::ValueEnum)]
pub enum ProfileFormat {
    /// JSON format (default). Compatible with `quality.profile_path`.
    Json,
    /// TSV format with one row per cycle and one column per quality value.
    Tsv,
}

/// CLI options for the `learn-profile` subcommand.
#[derive(Parser, Debug)]
pub struct LearnProfileOpts {
    /// Input BAM file to learn from
    #[arg(long)]
    pub bam: PathBuf,

    /// Output profile path (JSON or TSV depending on --format)
    #[arg(short, long)]
    pub output: PathBuf,

    /// Number of reads to sample (default: 1,000,000)
    #[arg(long, default_value_t = 1_000_000)]
    pub sample_size: usize,

    /// Minimum mapping quality for a read to be included
    #[arg(long, default_value_t = 20)]
    pub min_mapq: u8,

    /// Output format: `json` (default) or `tsv`
    #[arg(long, value_enum, default_value_t = ProfileFormat::Json)]
    pub format: ProfileFormat,
}

/// Entry point called from `main`.
pub fn run(opts: LearnProfileOpts, _threads: Option<usize>) -> Result<()> {
    info!(
        bam = %opts.bam.display(),
        output = %opts.output.display(),
        sample_size = opts.sample_size,
        "starting learn-profile"
    );

    let config = LearnerConfig {
        sample_size: opts.sample_size,
        min_mapq: opts.min_mapq,
    };

    let learner = ProfileLearner::new(config);
    let profile = learner
        .learn_from_bam(&opts.bam)
        .with_context(|| format!("failed to learn profile from {}", opts.bam.display()))?;

    match opts.format {
        ProfileFormat::Json => {
            let json = serde_json::to_string_pretty(&profile)
                .context("failed to serialise profile to JSON")?;
            std::fs::write(&opts.output, json.as_bytes())
                .with_context(|| format!("failed to write profile to {}", opts.output.display()))?;
        }
        ProfileFormat::Tsv => {
            write_quality_tsv(&profile, &opts.output).with_context(|| {
                format!("failed to write TSV profile to {}", opts.output.display())
            })?;
        }
    }

    info!(
        output = %opts.output.display(),
        read_length = profile.read_length,
        format = ?opts.format,
        "profile written"
    );

    Ok(())
}

/// Load a quality profile from a TSV file and return the per-cycle distribution.
///
/// This is a convenience wrapper used by the engine when `quality.profile_path`
/// points to a `.tsv` file.
// Not yet called from production code; reserved for a future engine hook.
#[allow(dead_code)]
pub fn load_quality_tsv(path: &std::path::Path) -> Result<Vec<Vec<[f64; 2]>>> {
    read_quality_tsv(path)
}
