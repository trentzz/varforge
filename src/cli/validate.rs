use anyhow::Result;
use super::ValidateOpts;
use crate::io::config;

pub fn run(opts: ValidateOpts) -> Result<()> {
    let cfg = config::load(&opts.config)?;
    config::validate(&cfg)?;
    tracing::info!("configuration is valid");
    Ok(())
}
