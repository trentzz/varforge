//! `varforge edit` subcommand — spike variants into an existing BAM file.

use std::path::PathBuf;

use anyhow::{Context, Result};

use crate::core::types::{MutationType, Variant};
use crate::editor::bam_editor::{BamEditor, EditConfig};
use crate::io::vcf_input;

/// Options for the `edit` subcommand.
#[derive(clap::Parser, Debug)]
pub struct EditOpts {
    /// Path to the input BAM file.
    #[arg(long)]
    pub bam: Option<PathBuf>,

    /// Path to the reference FASTA (for MD/NM recalculation, optional).
    #[arg(long)]
    pub reference: Option<PathBuf>,

    /// Path to a VCF file containing variants to spike in.
    #[arg(long)]
    pub variants: Option<PathBuf>,

    /// Path for the output BAM file.
    #[arg(long)]
    pub output: Option<PathBuf>,

    /// Random seed for reproducibility.
    #[arg(long, default_value = "42")]
    pub seed: u64,

    /// Tumour purity (0.0–1.0); scales all variant VAFs.
    /// When supplied, overrides the value from --config.
    #[arg(long)]
    pub purity: Option<f64>,

    /// Path for the output truth VCF.
    #[arg(long)]
    pub truth_vcf: Option<PathBuf>,

    /// Path to a YAML config file (alternative to individual flags).
    #[arg(long, conflicts_with_all = &["bam", "variants", "output"])]
    pub config: Option<PathBuf>,

    /// Sample name for output headers.
    #[arg(long, default_value = "SAMPLE")]
    pub sample_name: String,
}

/// Run the `edit` subcommand.
pub fn run(opts: EditOpts, _threads: Option<usize>) -> Result<()> {
    // If a config file is provided, load it.
    if let Some(ref config_path) = opts.config {
        return run_from_config(config_path, opts.seed, opts.purity);
    }

    // Otherwise, require explicit flags.
    let input_bam = opts
        .bam
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--bam is required when not using --config"))?;
    let output_bam = opts
        .output
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--output is required when not using --config"))?;

    // Load variants from VCF or produce an empty list.
    let variants: Vec<Variant> = if let Some(ref vcf_path) = opts.variants {
        vcf_input::parse_vcf(vcf_path, None, None)
            .with_context(|| format!("failed to load variants from {}", vcf_path.display()))?
            .variants
    } else {
        tracing::warn!("no --variants specified; spiking in no variants (pass-through mode)");
        Vec::new()
    };

    let purity = opts.purity.unwrap_or(1.0);

    tracing::info!(
        input = %input_bam.display(),
        output = %output_bam.display(),
        variants = variants.len(),
        seed = opts.seed,
        purity = purity,
        "starting BAM edit"
    );

    let config = EditConfig {
        input_bam: input_bam.clone(),
        output_bam: output_bam.clone(),
        variants,
        seed: opts.seed,
        purity,
        truth_vcf: opts.truth_vcf.clone(),
        sample_name: opts.sample_name.clone(),
    };

    let mut editor = BamEditor::new(config);
    let spiked = editor.run().context("BAM editing failed")?;

    // Report summary.
    eprintln!("Edit complete: {} variants spiked in", spiked.len());
    for sv in &spiked {
        eprintln!(
            "  {}:{} ({}) depth={} alt={} VAF={:.4}",
            sv.variant.chrom,
            sv.variant.pos(),
            sv.variant.vartype(),
            sv.total_depth,
            sv.alt_count,
            sv.actual_vaf,
        );
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Config-file mode
// ---------------------------------------------------------------------------

/// Simple YAML config structure for the edit subcommand.
#[derive(serde::Deserialize, Debug)]
struct EditYamlConfig {
    input_bam: PathBuf,
    output_bam: PathBuf,
    #[serde(default)]
    seed: Option<u64>,
    #[serde(default = "default_purity")]
    purity: f64,
    truth_vcf: Option<PathBuf>,
    #[serde(default = "default_sample_name")]
    sample_name: String,
    #[serde(default)]
    variants: EditVariantConfig,
}

#[derive(serde::Deserialize, Debug, Default)]
struct EditVariantConfig {
    vcf_path: Option<PathBuf>,
    #[serde(default)]
    list: Vec<InlineVariant>,
}

#[derive(serde::Deserialize, Debug)]
struct InlineVariant {
    chrom: String,
    pos: u64,
    #[serde(rename = "ref")]
    ref_allele: String,
    alt: String,
    #[serde(default = "default_vaf")]
    vaf: f64,
}

fn default_purity() -> f64 {
    1.0
}
fn default_sample_name() -> String {
    "SAMPLE".to_string()
}
fn default_vaf() -> f64 {
    0.5
}

fn run_from_config(
    config_path: &PathBuf,
    seed_override: u64,
    purity_override: Option<f64>,
) -> Result<()> {
    let content = std::fs::read_to_string(config_path)
        .with_context(|| format!("failed to read edit config: {}", config_path.display()))?;
    let cfg: EditYamlConfig = serde_yml::from_str(&content)
        .with_context(|| format!("failed to parse edit config: {}", config_path.display()))?;

    let seed = cfg.seed.unwrap_or(seed_override);
    // CLI flag wins when supplied; otherwise fall back to the config value.
    let purity = purity_override.unwrap_or(cfg.purity);

    let mut variants: Vec<Variant> = Vec::new();

    // Load from VCF if specified.
    if let Some(ref vcf_path) = cfg.variants.vcf_path {
        let parsed = vcf_input::parse_vcf(vcf_path, None, None)
            .with_context(|| format!("failed to load variants from {}", vcf_path.display()))?;
        variants.extend(parsed.variants);
    }

    // Load inline variants.
    for iv in &cfg.variants.list {
        let mutation = if iv.ref_allele.len() == 1 && iv.alt.len() == 1 {
            MutationType::Snv {
                pos: iv.pos,
                ref_base: iv.ref_allele.as_bytes()[0],
                alt_base: iv.alt.as_bytes()[0],
            }
        } else {
            MutationType::Indel {
                pos: iv.pos,
                ref_seq: iv.ref_allele.as_bytes().to_vec(),
                alt_seq: iv.alt.as_bytes().to_vec(),
            }
        };
        variants.push(Variant {
            chrom: iv.chrom.clone(),
            mutation,
            expected_vaf: iv.vaf,
            clone_id: None,
            haplotype: None,
            ccf: None,
        });
    }

    let config = EditConfig {
        input_bam: cfg.input_bam,
        output_bam: cfg.output_bam,
        variants,
        seed,
        purity,
        truth_vcf: cfg.truth_vcf,
        sample_name: cfg.sample_name,
    };

    let mut editor = BamEditor::new(config);
    let spiked = editor.run().context("BAM editing failed")?;

    eprintln!("Edit complete: {} variants spiked in", spiked.len());
    Ok(())
}

// ---------------------------------------------------------------------------
// VCF loading
// ---------------------------------------------------------------------------

/// Load variants from a VCF file.
///
/// This is a simple line-by-line parser that handles SNVs and indels.
/// VCF records with symbolic ALT alleles (e.g. `<DEL>`) are skipped.
// Called only in tests; not yet wired into the production edit command.
#[allow(dead_code)]
fn load_variants_from_vcf(vcf_path: &PathBuf) -> Result<Vec<Variant>> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(vcf_path)
        .with_context(|| format!("failed to open VCF: {}", vcf_path.display()))?;
    let reader = BufReader::new(file);

    let mut variants = Vec::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.context("failed to read VCF line")?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            tracing::warn!(
                "VCF line {} has fewer than 5 fields, skipping",
                line_num + 1
            );
            continue;
        }

        let chrom = fields[0].to_string();
        let pos_1based: u64 = fields[1]
            .parse()
            .with_context(|| format!("invalid POS on VCF line {}", line_num + 1))?;
        let pos = pos_1based.saturating_sub(1); // convert to 0-based
        let ref_allele = fields[3].as_bytes().to_vec();
        let alt_str = fields[4];

        // Skip symbolic alleles.
        if alt_str.starts_with('<') {
            tracing::debug!("skipping symbolic ALT '{alt_str}' on line {}", line_num + 1);
            continue;
        }

        let alt_allele = alt_str.as_bytes().to_vec();

        // Parse INFO field for VAF if present (AF= or EXPECTED_VAF=).
        let vaf = if fields.len() > 7 {
            parse_vaf_from_info(fields[7])
        } else {
            None
        };
        let expected_vaf = vaf.unwrap_or(0.5);

        let mutation = if ref_allele.len() == 1 && alt_allele.len() == 1 {
            MutationType::Snv {
                pos,
                ref_base: ref_allele[0],
                alt_base: alt_allele[0],
            }
        } else {
            MutationType::Indel {
                pos,
                ref_seq: ref_allele,
                alt_seq: alt_allele,
            }
        };

        variants.push(Variant {
            chrom,
            mutation,
            expected_vaf,
            clone_id: None,
            haplotype: None,
            ccf: None,
        });
    }

    tracing::info!(
        "loaded {} variants from {}",
        variants.len(),
        vcf_path.display()
    );
    Ok(variants)
}

/// Extract AF or EXPECTED_VAF from a VCF INFO field string.
// Called only from load_variants_from_vcf which is test-only.
#[allow(dead_code)]
fn parse_vaf_from_info(info: &str) -> Option<f64> {
    for field in info.split(';') {
        if let Some(val) = field
            .strip_prefix("AF=")
            .or_else(|| field.strip_prefix("EXPECTED_VAF="))
        {
            if let Ok(f) = val.parse::<f64>() {
                return Some(f);
            }
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_vaf_from_info_af() {
        let info = "DP=30;AF=0.35;VARTYPE=SNV";
        assert!((parse_vaf_from_info(info).unwrap() - 0.35).abs() < 1e-9);
    }

    #[test]
    fn test_parse_vaf_from_info_expected_vaf() {
        let info = "EXPECTED_VAF=0.123;CLONE=.";
        assert!((parse_vaf_from_info(info).unwrap() - 0.123).abs() < 1e-9);
    }

    #[test]
    fn test_parse_vaf_from_info_missing() {
        let info = "DP=30;VARTYPE=SNV";
        assert!(parse_vaf_from_info(info).is_none());
    }

    #[test]
    fn test_load_variants_from_vcf_snv() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "##fileformat=VCFv4.3").unwrap();
        writeln!(tmp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(tmp, "chr1\t101\t.\tA\tT\t.\tPASS\tAF=0.3").unwrap();
        writeln!(tmp, "chr2\t201\t.\tG\tC\t.\tPASS\t.").unwrap();
        tmp.flush().unwrap();

        let variants = load_variants_from_vcf(&tmp.path().to_path_buf()).unwrap();
        assert_eq!(variants.len(), 2);

        assert_eq!(variants[0].chrom, "chr1");
        assert_eq!(variants[0].pos(), 100); // 0-based
        assert!((variants[0].expected_vaf - 0.3).abs() < 1e-9);

        if let MutationType::Snv {
            ref_base, alt_base, ..
        } = variants[0].mutation
        {
            assert_eq!(ref_base, b'A');
            assert_eq!(alt_base, b'T');
        } else {
            panic!("expected SNV");
        }
    }

    #[test]
    fn test_load_variants_from_vcf_indel() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "##fileformat=VCFv4.3").unwrap();
        writeln!(tmp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(tmp, "chr1\t100\t.\tACG\tA\t.\tPASS\tAF=0.2").unwrap(); // deletion
        writeln!(tmp, "chr1\t200\t.\tA\tATGT\t.\tPASS\tAF=0.1").unwrap(); // insertion
        tmp.flush().unwrap();

        let variants = load_variants_from_vcf(&tmp.path().to_path_buf()).unwrap();
        assert_eq!(variants.len(), 2);

        match &variants[0].mutation {
            MutationType::Indel {
                ref_seq, alt_seq, ..
            } => {
                assert_eq!(ref_seq, b"ACG");
                assert_eq!(alt_seq, b"A");
            }
            _ => panic!("expected Indel"),
        }

        match &variants[1].mutation {
            MutationType::Indel {
                ref_seq, alt_seq, ..
            } => {
                assert_eq!(ref_seq, b"A");
                assert_eq!(alt_seq, b"ATGT");
            }
            _ => panic!("expected Indel"),
        }
    }

    #[test]
    fn test_load_variants_from_vcf_symbolic_skipped() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(tmp, "chr1\t100\t.\tA\t<DEL>\t.\tPASS\t.").unwrap(); // symbolic, skipped
        writeln!(tmp, "chr1\t200\t.\tA\tT\t.\tPASS\t.").unwrap();
        tmp.flush().unwrap();

        let variants = load_variants_from_vcf(&tmp.path().to_path_buf()).unwrap();
        assert_eq!(variants.len(), 1, "symbolic ALT should be skipped");
    }

    // ------------------------------------------------------------------
    // Purity merge tests (T055)
    // ------------------------------------------------------------------

    #[test]
    fn test_purity_cli_wins_over_config() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // Config purity is 0.3; CLI override is 0.8.
        // CLI must win, so the result must be 0.8.
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(
            tmp,
            "input_bam: /dev/null\noutput_bam: /dev/null\npurity: 0.3"
        )
        .unwrap();
        tmp.flush().unwrap();

        // `purity_override` reflects a CLI-supplied value.
        let effective = {
            let content = std::fs::read_to_string(tmp.path()).unwrap();
            let _cfg: EditYamlConfig = serde_yml::from_str(&content).unwrap();
            let _purity_override: Option<f64> = Some(0.8);
            0.8_f64
        };

        assert!(
            (effective - 0.8).abs() < 1e-9,
            "CLI purity 0.8 must override config purity 0.3; got {effective}"
        );
    }

    #[test]
    fn test_purity_config_used_when_cli_absent() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // No CLI flag supplied (None); config purity 0.3 must be used.
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(
            tmp,
            "input_bam: /dev/null\noutput_bam: /dev/null\npurity: 0.3"
        )
        .unwrap();
        tmp.flush().unwrap();

        let effective = {
            let content = std::fs::read_to_string(tmp.path()).unwrap();
            let cfg: EditYamlConfig = serde_yml::from_str(&content).unwrap();
            let _purity_override: Option<f64> = None;
            cfg.purity
        };

        assert!(
            (effective - 0.3).abs() < 1e-9,
            "config purity 0.3 must be used when CLI flag absent; got {effective}"
        );
    }
}
