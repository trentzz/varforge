use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub reference: PathBuf,
    pub output: OutputConfig,
    #[serde(default)]
    pub sample: SampleConfig,
    #[serde(default)]
    pub fragment: FragmentConfig,
    #[serde(default)]
    pub quality: QualityConfig,
    #[serde(default)]
    pub tumour: Option<TumourConfig>,
    #[serde(default)]
    pub mutations: Option<MutationConfig>,
    #[serde(default)]
    pub umi: Option<UmiConfig>,
    #[serde(default)]
    pub artifacts: Option<ArtifactConfig>,
    #[serde(default)]
    pub seed: Option<u64>,
    #[serde(default)]
    pub threads: Option<usize>,
    #[serde(default)]
    pub chromosomes: Option<Vec<String>>,
    #[serde(default)]
    pub regions_bed: Option<PathBuf>,
    #[serde(default)]
    pub copy_number: Option<Vec<CopyNumberConfig>>,
    #[serde(default)]
    pub gc_bias: Option<GcBiasConfig>,
    /// Optional multi-sample configuration for longitudinal simulation.
    ///
    /// When present, VarForge will simulate each entry independently from a
    /// shared clonal architecture, writing per-sample output sub-directories
    /// plus a combined manifest.
    #[serde(default)]
    pub samples: Option<Vec<SampleEntry>>,
    #[serde(default)]
    pub capture: Option<CaptureConfig>,
    #[serde(default)]
    pub performance: PerformanceConfig,
}

/// Performance tuning parameters for the streaming output pipeline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceConfig {
    /// Maximum number of region batches buffered in the streaming channel.
    ///
    /// Higher values use more memory but provide more overlap between compute
    /// and I/O. Lower values reduce peak memory. Default: 64.
    #[serde(default = "default_output_buffer_regions")]
    pub output_buffer_regions: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputConfig {
    pub directory: PathBuf,
    #[serde(default = "default_true")]
    pub fastq: bool,
    #[serde(default)]
    pub bam: bool,
    #[serde(default = "default_true")]
    pub truth_vcf: bool,
    #[serde(default = "default_true")]
    pub manifest: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleConfig {
    #[serde(default = "default_sample_name")]
    pub name: String,
    #[serde(default = "default_read_length")]
    pub read_length: usize,
    #[serde(default = "default_coverage")]
    pub coverage: f64,
    #[serde(default)]
    pub platform: Option<String>,
}

impl Default for SampleConfig {
    fn default() -> Self {
        Self {
            name: default_sample_name(),
            read_length: default_read_length(),
            coverage: default_coverage(),
            platform: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentConfig {
    #[serde(default = "default_fragment_model")]
    pub model: FragmentModel,
    #[serde(default = "default_fragment_mean")]
    pub mean: f64,
    #[serde(default = "default_fragment_sd")]
    pub sd: f64,
}

impl Default for FragmentConfig {
    fn default() -> Self {
        Self {
            model: FragmentModel::Normal,
            mean: default_fragment_mean(),
            sd: default_fragment_sd(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FragmentModel {
    Normal,
    Cfda,
    Custom,
}

fn default_fragment_model() -> FragmentModel {
    FragmentModel::Normal
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityConfig {
    #[serde(default = "default_mean_quality")]
    pub mean_quality: u8,
    #[serde(default = "default_quality_decay")]
    pub tail_decay: f64,
    #[serde(default)]
    pub profile_path: Option<PathBuf>,
}

impl Default for QualityConfig {
    fn default() -> Self {
        Self {
            mean_quality: default_mean_quality(),
            tail_decay: default_quality_decay(),
            profile_path: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TumourConfig {
    #[serde(default = "default_purity")]
    pub purity: f64,
    #[serde(default = "default_ploidy")]
    pub ploidy: u32,
    #[serde(default)]
    pub clones: Vec<CloneConfig>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CloneConfig {
    pub id: String,
    pub ccf: f64,
    #[serde(default)]
    pub parent: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MutationConfig {
    #[serde(default)]
    pub vcf: Option<PathBuf>,
    #[serde(default)]
    pub random: Option<RandomMutationConfig>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RandomMutationConfig {
    pub count: usize,
    #[serde(default = "default_vaf_min")]
    pub vaf_min: f64,
    #[serde(default = "default_vaf_max")]
    pub vaf_max: f64,
    #[serde(default = "default_snv_fraction")]
    pub snv_fraction: f64,
    #[serde(default = "default_indel_fraction")]
    pub indel_fraction: f64,
    #[serde(default = "default_mnv_fraction")]
    pub mnv_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UmiConfig {
    #[serde(default = "default_umi_length")]
    pub length: usize,
    #[serde(default)]
    pub duplex: bool,
    #[serde(default = "default_pcr_cycles")]
    pub pcr_cycles: u32,
    #[serde(default = "default_family_size_mean")]
    pub family_size_mean: f64,
    #[serde(default = "default_family_size_sd")]
    pub family_size_sd: f64,
    #[serde(default)]
    pub inline: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArtifactConfig {
    #[serde(default)]
    pub ffpe_damage_rate: Option<f64>,
    #[serde(default)]
    pub oxog_rate: Option<f64>,
    #[serde(default)]
    pub duplicate_rate: Option<f64>,
    #[serde(default)]
    pub pcr_error_rate: Option<f64>,
}

/// A single copy number alteration entry from the YAML config.
///
/// The `region` field uses the compact `chr:start-end` notation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CopyNumberConfig {
    /// Genomic region in `chrom:start-end` format, e.g. `"chr7:55000000-55200000"`.
    pub region: String,
    #[serde(default = "default_normal_cn")]
    pub tumor_cn: u32,
    #[serde(default = "default_normal_cn")]
    pub normal_cn: u32,
    /// Major allele copy number (optional; for allele-specific / LOH modeling).
    #[serde(default)]
    pub major_cn: Option<u32>,
    /// Minor allele copy number (optional; for allele-specific / LOH modeling).
    #[serde(default)]
    pub minor_cn: Option<u32>,
}

/// Configuration for GC-content coverage bias (mirrors `gc_bias:` YAML section).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GcBiasConfig {
    /// Whether GC bias is applied.
    #[serde(default = "default_true")]
    pub enabled: bool,
    /// Bias curve: `"default"`, `"flat"`, or `"custom"`.
    #[serde(default = "default_gc_bias_model")]
    pub model: String,
    /// Severity multiplier (0 = no bias, 1 = realistic, 2 = extreme).
    #[serde(default = "default_gc_bias_severity")]
    pub severity: f64,
}

impl Default for GcBiasConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            model: default_gc_bias_model(),
            severity: default_gc_bias_severity(),
        }
    }
}

/// A single sample entry in the `samples:` list for multi-sample / longitudinal
/// simulation.
///
/// Each entry overrides the per-sample coverage, tumour fraction, and fragment
/// model while sharing the same reference, mutation list, and clonal tree
/// defined at the top level.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleEntry {
    /// Output sub-directory name and prefix for FASTQ / VCF files.
    pub name: String,
    /// Target sequencing coverage for this sample.
    #[serde(default = "default_coverage")]
    pub coverage: f64,
    /// Tumour fraction (ctDNA fraction) for this time point.
    #[serde(default = "default_tumour_fraction")]
    pub tumour_fraction: f64,
    /// Fragment model override (`normal` or `cfDNA`).  Falls back to the
    /// top-level `fragment.model` when absent.
    #[serde(default)]
    pub fragment_model: Option<FragmentModel>,
    /// Optional per-clone CCF adjustments for this time point.
    ///
    /// Keys are clone IDs; values are the new CCF (0.0–1.0).
    #[serde(default)]
    pub clonal_shift: std::collections::HashMap<String, f64>,
}

fn default_tumour_fraction() -> f64 {
    1.0
}

/// Configuration for target-capture efficiency (mirrors `capture:` YAML section).
///
/// ```yaml
/// capture:
///   enabled: true
///   targets_bed: "panel.bed"
///   off_target_fraction: 0.2
///   coverage_uniformity: 0.3
///   edge_dropoff_bases: 50
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CaptureConfig {
    /// Whether the capture model is active.
    #[serde(default = "default_true")]
    pub enabled: bool,
    /// Path to BED file listing capture target regions.
    #[serde(default)]
    pub targets_bed: Option<PathBuf>,
    /// Fraction of reads that map off-target (0.0–1.0; default 0.2).
    #[serde(default = "default_off_target_fraction")]
    pub off_target_fraction: f64,
    /// Per-target coverage variation via LogNormal σ (0 = uniform; default 0.3).
    #[serde(default = "default_coverage_uniformity")]
    pub coverage_uniformity: f64,
    /// Bases of exponential coverage drop-off at target edges (default 50).
    #[serde(default = "default_edge_dropoff_bases")]
    pub edge_dropoff_bases: u32,
}

impl Default for CaptureConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            targets_bed: None,
            off_target_fraction: default_off_target_fraction(),
            coverage_uniformity: default_coverage_uniformity(),
            edge_dropoff_bases: default_edge_dropoff_bases(),
        }
    }
}

// Default value functions
fn default_true() -> bool {
    true
}
fn default_normal_cn() -> u32 {
    2
}
fn default_sample_name() -> String {
    "SAMPLE".to_string()
}
fn default_read_length() -> usize {
    150
}
pub fn default_coverage() -> f64 {
    30.0
}
pub fn default_fragment_mean() -> f64 {
    300.0
}
pub fn default_fragment_sd() -> f64 {
    50.0
}
fn default_mean_quality() -> u8 {
    36
}
fn default_quality_decay() -> f64 {
    0.003
}
fn default_purity() -> f64 {
    1.0
}
fn default_ploidy() -> u32 {
    2
}
fn default_vaf_min() -> f64 {
    0.001
}
fn default_vaf_max() -> f64 {
    0.5
}
fn default_snv_fraction() -> f64 {
    0.80
}
fn default_indel_fraction() -> f64 {
    0.15
}
fn default_mnv_fraction() -> f64 {
    0.05
}
fn default_umi_length() -> usize {
    8
}
fn default_pcr_cycles() -> u32 {
    10
}
fn default_family_size_mean() -> f64 {
    3.0
}
fn default_family_size_sd() -> f64 {
    1.5
}
fn default_gc_bias_model() -> String {
    "default".to_string()
}
fn default_gc_bias_severity() -> f64 {
    1.0
}
impl Default for PerformanceConfig {
    fn default() -> Self {
        Self {
            output_buffer_regions: default_output_buffer_regions(),
        }
    }
}

fn default_output_buffer_regions() -> usize {
    64
}
fn default_off_target_fraction() -> f64 {
    0.2
}
fn default_coverage_uniformity() -> f64 {
    0.3
}
fn default_edge_dropoff_bases() -> u32 {
    50
}

pub fn load(path: &Path) -> Result<Config> {
    let contents = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read config file: {}", path.display()))?;
    let config: Config = serde_yaml::from_str(&contents)
        .with_context(|| format!("failed to parse config file: {}", path.display()))?;
    Ok(config)
}

pub fn validate(config: &Config) -> Result<()> {
    anyhow::ensure!(
        config.reference.exists(),
        "reference file not found: {}",
        config.reference.display()
    );

    if let Some(ref bed) = config.regions_bed {
        anyhow::ensure!(bed.exists(), "regions_bed file not found: {}", bed.display());
    }

    anyhow::ensure!(
        config.sample.coverage > 0.0,
        "coverage must be positive, got {}",
        config.sample.coverage
    );
    anyhow::ensure!(
        config.sample.read_length > 0,
        "read_length must be positive, got {}",
        config.sample.read_length
    );
    anyhow::ensure!(
        config.fragment.mean > 0.0,
        "fragment mean must be positive, got {}",
        config.fragment.mean
    );
    anyhow::ensure!(
        config.fragment.sd >= 0.0,
        "fragment sd must be non-negative, got {}",
        config.fragment.sd
    );

    if let Some(tumour) = &config.tumour {
        anyhow::ensure!(
            (0.0..=1.0).contains(&tumour.purity),
            "tumour purity must be between 0.0 and 1.0, got {}",
            tumour.purity
        );
        for clone in &tumour.clones {
            anyhow::ensure!(
                (0.0..=1.0).contains(&clone.ccf),
                "clone {} CCF must be between 0.0 and 1.0, got {}",
                clone.id,
                clone.ccf
            );
        }
    }

    if let Some(mutations) = &config.mutations {
        if let Some(random) = &mutations.random {
            anyhow::ensure!(random.count > 0, "random mutation count must be > 0");
            anyhow::ensure!(
                random.vaf_min < random.vaf_max,
                "vaf_min ({}) must be less than vaf_max ({})",
                random.vaf_min,
                random.vaf_max
            );
            let total = random.snv_fraction + random.indel_fraction + random.mnv_fraction;
            anyhow::ensure!(
                (total - 1.0).abs() < 1e-6,
                "mutation type fractions must sum to 1.0, got {}",
                total
            );
        }
    }

    if let Some(umi) = &config.umi {
        anyhow::ensure!(umi.length > 0, "UMI length must be > 0");
    }

    Ok(())
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

    #[test]
    fn test_minimal_config_parse() {
        let yaml = r#"
reference: /tmp/ref.fa
output:
  directory: /tmp/out
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert_eq!(cfg.sample.read_length, 150);
        assert_eq!(cfg.sample.coverage, 30.0);
        assert_eq!(cfg.fragment.mean, 300.0);
    }

    #[test]
    fn test_full_config_parse() {
        let yaml = r#"
reference: /tmp/ref.fa
output:
  directory: /tmp/out
  fastq: true
  bam: true
  truth_vcf: true
sample:
  name: TUMOUR_01
  read_length: 100
  coverage: 60.0
fragment:
  model: normal
  mean: 250.0
  sd: 40.0
tumour:
  purity: 0.7
  ploidy: 2
  clones:
    - id: clone_a
      ccf: 1.0
    - id: clone_b
      ccf: 0.3
      parent: clone_a
mutations:
  random:
    count: 50
    vaf_min: 0.01
    vaf_max: 0.5
umi:
  length: 12
  duplex: true
  pcr_cycles: 8
artifacts:
  ffpe_damage_rate: 0.01
  duplicate_rate: 0.15
seed: 42
chromosomes:
  - chr22
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert_eq!(cfg.sample.name, "TUMOUR_01");
        assert_eq!(cfg.sample.read_length, 100);
        assert_eq!(cfg.tumour.as_ref().unwrap().purity, 0.7);
        assert_eq!(cfg.tumour.as_ref().unwrap().clones.len(), 2);
        assert_eq!(
            cfg.mutations
                .as_ref()
                .unwrap()
                .random
                .as_ref()
                .unwrap()
                .count,
            50
        );
        assert!(cfg.umi.as_ref().unwrap().duplex);
        assert_eq!(cfg.seed, Some(42));
    }

    #[test]
    fn test_validate_bad_purity() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/out
tumour:
  purity: 1.5
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert!(validate(&cfg).is_err());
    }

    #[test]
    fn test_validate_bad_vaf_range() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/out
mutations:
  random:
    count: 10
    vaf_min: 0.5
    vaf_max: 0.1
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert!(validate(&cfg).is_err());
    }

    #[test]
    fn test_validate_bad_type_fractions() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/out
mutations:
  random:
    count: 10
    snv_fraction: 0.5
    indel_fraction: 0.5
    mnv_fraction: 0.5
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert!(validate(&cfg).is_err());
    }

    #[test]
    fn test_cfda_fragment_model() {
        let yaml = r#"
reference: /tmp/ref.fa
output:
  directory: /tmp/out
fragment:
  model: cfda
  mean: 167.0
  sd: 20.0
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert!(matches!(cfg.fragment.model, FragmentModel::Cfda));
    }
}
