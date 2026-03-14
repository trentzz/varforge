mod cli;
mod core;
mod variants;
mod tumour;
mod umi;
mod artifacts;
mod io;

use anyhow::Result;
use clap::Parser;

fn main() -> Result<()> {
    let args = cli::Args::parse();

    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_env("VARFORGE_LOG")
                .unwrap_or_else(|_| {
                    tracing_subscriber::EnvFilter::new(args.log_level())
                }),
        )
        .with_writer(std::io::stderr)
        .init();

    match args.command {
        cli::Command::Simulate(opts) => cli::simulate::run(opts),
        cli::Command::Validate(opts) => cli::validate::run(opts),
    }
}
