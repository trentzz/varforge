mod artifacts;
mod cli;
mod core;
mod editor;
mod io;
mod seq_utils;
mod tumour;
mod umi;
mod variants;

use anyhow::Result;
use clap::Parser;

fn main() -> Result<()> {
    let args = cli::Args::parse();

    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_env("VARFORGE_LOG")
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new(args.log_level())),
        )
        .with_writer(std::io::stderr)
        .init();

    let threads = args.threads;
    match args.command {
        cli::Command::Simulate(opts) => cli::simulate::run(opts, threads),
        cli::Command::Validate(opts) => cli::validate::run(opts),
        cli::Command::Edit(opts) => cli::edit::run(opts, threads),
        cli::Command::LearnProfile(opts) => cli::learn_profile::run(opts, threads),
        cli::Command::BenchmarkSuite(opts) => cli::benchmark_suite::run(opts, threads),
    }
}
