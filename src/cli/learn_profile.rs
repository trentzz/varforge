//! `varforge learn-profile` subcommand.
//!
//! Analyzes a real BAM file and emits an error/quality profile JSON that can
//! be fed back into VarForge via `quality.profile_path`.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use tracing::info;

use crate::core::profile_learner::{LearnerConfig, ProfileLearner};

/// CLI options for the `learn-profile` subcommand.
#[derive(Parser, Debug)]
pub struct LearnProfileOpts {
    /// Input BAM file to learn from
    #[arg(long)]
    pub bam: PathBuf,

    /// Output profile JSON path
    #[arg(short, long)]
    pub output: PathBuf,

    /// Number of reads to sample (default: 1,000,000)
    #[arg(long, default_value_t = 1_000_000)]
    pub sample_size: usize,

    /// Minimum mapping quality for a read to be included
    #[arg(long, default_value_t = 20)]
    pub min_mapq: u8,
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

    let json =
        serde_json::to_string_pretty(&profile).context("failed to serialize profile to JSON")?;

    std::fs::write(&opts.output, json.as_bytes())
        .with_context(|| format!("failed to write profile to {}", opts.output.display()))?;

    info!(
        output = %opts.output.display(),
        read_length = profile.read_length,
        "profile written"
    );

    Ok(())
}
