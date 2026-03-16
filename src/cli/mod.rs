pub mod cancer_presets;
pub mod edit;
pub mod learn_profile;
pub mod presets;
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
    /// Spike variants into an existing BAM file
    Edit(edit::EditOpts),
    /// Learn an error/quality profile from a real BAM file
    LearnProfile(learn_profile::LearnProfileOpts),
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

    /// Override read length (bp)
    #[arg(long)]
    pub read_length: Option<usize>,

    /// Override tumour purity (0.0–1.0)
    #[arg(long)]
    pub purity: Option<f64>,

    /// Override fragment mean length (bp)
    #[arg(long)]
    pub fragment_mean: Option<f64>,

    /// Override fragment length standard deviation (bp)
    #[arg(long)]
    pub fragment_sd: Option<f64>,

    /// Generate N random mutations (no VCF needed)
    #[arg(long)]
    pub random_mutations: Option<usize>,

    /// VAF range for random mutations (e.g., 0.001-0.05)
    #[arg(long)]
    pub vaf_range: Option<String>,

    /// Simulation preset: small, panel, wgs, cfdna, ffpe, umi
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
