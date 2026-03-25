//! Implementation of the `validate` subcommand.
//!
//! Loads and structurally validates a YAML config, then performs pre-flight
//! checks on the reference genome, BED targets, and mutations VCF.
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_yaml(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    /// run() returns Ok for a config that passes structural validation.
    ///
    /// The pre-flight reference/BED/VCF checks are best-effort (they print errors
    /// rather than returning Err), so a config with a non-existent reference still
    /// returns Ok from run().
    #[test]
    fn test_run_valid_config_returns_ok() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/varforge_test_validate_ok
"#;
        let f = write_yaml(yaml);
        let opts = ValidateOpts {
            config: f.path().to_path_buf(),
        };
        let result = run(opts);
        assert!(
            result.is_ok(),
            "valid config should return Ok: {:?}",
            result
        );
    }

    /// run() returns Err when the config file does not exist.
    #[test]
    fn test_run_missing_config_errors() {
        let opts = ValidateOpts {
            config: std::path::PathBuf::from("/nonexistent/path/config.yaml"),
        };
        assert!(run(opts).is_err(), "missing config file should return Err");
    }

    /// run() returns Err for a malformed YAML config (missing required fields).
    #[test]
    fn test_run_invalid_yaml_errors() {
        let yaml = "not: valid: yaml: [[[";
        let f = write_yaml(yaml);
        let opts = ValidateOpts {
            config: f.path().to_path_buf(),
        };
        assert!(run(opts).is_err(), "invalid YAML should return Err");
    }

    // T062: Tests for parse_bed_count.

    #[test]
    fn test_parse_bed_count_basic() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t100\t200").unwrap();
        writeln!(f, "chr1\t300\t400").unwrap();
        let count = parse_bed_count(f.path()).unwrap();
        assert_eq!(count, 2);
    }

    #[test]
    fn test_parse_bed_count_skips_comments() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "# comment").unwrap();
        writeln!(f, "track name=foo").unwrap();
        writeln!(f, "browser position chr1:1-100").unwrap();
        writeln!(f, "chr1\t100\t200").unwrap();
        let count = parse_bed_count(f.path()).unwrap();
        assert_eq!(count, 1);
    }

    #[test]
    fn test_parse_bed_count_missing_file_errors() {
        let result = parse_bed_count(std::path::Path::new("/nonexistent/bed.bed"));
        assert!(result.is_err());
    }
}
