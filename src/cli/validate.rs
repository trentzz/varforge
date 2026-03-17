use super::ValidateOpts;
use crate::io::config;
use anyhow::Result;

pub fn run(opts: ValidateOpts) -> Result<()> {
    let cfg = config::load(&opts.config)?;
    config::validate(&cfg)?;
    tracing::info!("configuration is valid");
    Ok(())
}
