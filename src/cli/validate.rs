use super::ValidateOpts;
use crate::io::config;
use anyhow::Result;

pub fn run(opts: ValidateOpts) -> Result<()> {
    let cfg = config::load(&opts.config)?;
    config::validate(&cfg)?;

    // Pre-flight checks: attempt to open key resources and report diagnostics.

    // Reference genome: attempt to open and read chromosome list.
    let ref_status = match crate::io::reference::ReferenceGenome::open(&cfg.reference) {
        Ok(genome) => {
            let n_chroms = genome.chromosome_lengths().len();
            format!("readable ({n_chroms} chromosomes)")
        }
        Err(e) => format!("ERROR: {e}"),
    };

    // BED targets: parse and count regions.
    let bed_status = if let Some(ref bed_path) = cfg.regions_bed {
        match parse_bed_count(bed_path) {
            Ok(n) => format!("{n} targets"),
            Err(e) => format!("ERROR: {e}"),
        }
    } else {
        "not set".to_string()
    };

    // Mutations VCF: parse and count variants.
    let vcf_status = if let Some(ref mutations) = cfg.mutations {
        if let Some(ref vcf_path) = mutations.vcf {
            match crate::io::vcf_input::parse_vcf(vcf_path, None, None) {
                Ok(result) => format!("{} variants", result.variants.len()),
                Err(e) => format!("ERROR: {e}"),
            }
        } else {
            "not set".to_string()
        }
    } else {
        "not set".to_string()
    };

    println!(
        "Config OK. Reference: {}. BED: {}. VCF: {}.",
        ref_status, bed_status, vcf_status
    );

    Ok(())
}

/// Parse a BED file and return the count of valid target regions.
fn parse_bed_count(path: &std::path::Path) -> Result<usize> {
    use std::io::BufRead;

    let file =
        std::fs::File::open(path).map_err(|e| anyhow::anyhow!("failed to open BED file: {e}"))?;
    let reader = std::io::BufReader::new(file);
    let mut count = 0usize;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| anyhow::anyhow!("BED read error at line {line_no}: {e}"))?;
        let trimmed = line.trim();
        // Skip blank lines and track/browser/comment headers.
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track ")
            || trimmed.starts_with("browser ")
        {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 3 {
            anyhow::bail!("BED line {} has fewer than 3 fields", line_no + 1);
        }
        count += 1;
    }

    Ok(count)
}
