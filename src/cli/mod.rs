pub mod simulate;
pub mod validate;

use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "varforge", version, about = "Synthetic cancer sequencing data generator")]
pub struct Args {
    #[command(subcommand)]
    pub command: Command,

    /// Log level: error, warn, info, debug, trace
    #[arg(long, global = true, default_value = "info")]
    pub log_level: String,

    /// Number of threads (default: all available cores)
    #[arg(short = 't', long, global = true)]
    pub threads: Option<usize>,
}

impl Args {
    pub fn log_level(&self) -> &str {
        &self.log_level
    }
}

#[derive(Subcommand)]
pub enum Command {
    /// Run a simulation from a YAML config
    Simulate(SimulateOpts),
    /// Validate a YAML config without running
    Validate(ValidateOpts),
}

#[derive(Parser)]
pub struct SimulateOpts {
    /// Path to YAML configuration file
    #[arg(short, long)]
    pub config: PathBuf,

    /// Override output directory
    #[arg(short, long)]
    pub output_dir: Option<PathBuf>,

    /// Override random seed
    #[arg(long)]
    pub seed: Option<u64>,

    /// Override coverage depth
    #[arg(long)]
    pub coverage: Option<f64>,

    /// Generate N random mutations (no VCF needed)
    #[arg(long)]
    pub random_mutations: Option<usize>,

    /// VAF range for random mutations (e.g., 0.001-0.05)
    #[arg(long)]
    pub vaf_range: Option<String>,

    /// Simulation preset: small, panel, wgs
    #[arg(long)]
    pub preset: Option<String>,

    /// Validate config only, report estimated output size
    #[arg(long)]
    pub dry_run: bool,
}

#[derive(Parser)]
pub struct ValidateOpts {
    /// Path to YAML configuration file
    #[arg(short, long)]
    pub config: PathBuf,
}
