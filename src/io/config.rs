//! YAML configuration structs for all simulation parameters, with serde deserialisation and validation.

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
    /// Optional chemistry preset name. When set, fills any fields that are
    /// absent or at their zero default with values appropriate for that
    /// sequencing chemistry. Explicit YAML values always take precedence.
    #[serde(default)]
    pub preset: Option<String>,
    /// When set, drives batch mode: one simulation run per VAF value.
    /// Overrides `mutations.random.vaf_min`/`vaf_max`.
    #[serde(default)]
    pub vafs: Option<Vec<f64>>,
    /// Optional germline variant simulation configuration.
    ///
    /// When present, germline SNPs and indels at population-frequency densities
    /// are added to every simulated sample alongside any somatic variants.
    #[serde(default)]
    pub germline: Option<GermlineConfig>,
    /// Optional paired tumour-normal simulation configuration.
    ///
    /// When present, VarForge runs two simulations: one tumour sample (with all
    /// configured somatic and germline variants) and one normal sample (germline
    /// variants only). Outputs are written to `tumour/` and `normal/`
    /// sub-directories under `output.directory`.
    #[serde(default)]
    pub paired: Option<PairedConfig>,
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
    /// Write BAM in single-read mode (one record per read, not paired).
    ///
    /// Use this for long-read platforms where each molecule is a single read.
    #[serde(default)]
    pub single_read_bam: bool,
    /// Write `germline_truth.vcf` alongside the somatic truth VCF when
    /// germline variants were simulated.
    #[serde(default = "default_true")]
    pub germline_vcf: bool,
    /// Mapping quality written to every BAM record.
    ///
    /// Default 60 is appropriate for short-read simulated data where
    /// alignment is perfect. Lower values suit error-prone platforms.
    #[serde(default = "default_mapq")]
    pub mapq: u8,
    /// Annotate FASTQ read names with variant information for reads that carry
    /// a spiked-in variant.
    ///
    /// When enabled, each read name gains one or more space-separated tags of
    /// the form `VT:Z:<chrom>:<pos>:<type>` (e.g. `VT:Z:chr1:1000:SNV`).
    /// Disabled by default to keep read names clean for production use.
    /// Enable for debugging or truth-labelled benchmarking datasets.
    #[serde(default)]
    pub annotate_reads: bool,
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
    /// Long-read fragment length model. When set, the log-normal sampler is
    /// used instead of the standard normal or cfDNA sampler.
    #[serde(default)]
    pub long_read: Option<LongReadFragmentConfig>,
    /// End motif model for fragment selection.
    ///
    /// When set to `"plasma"`, fragment start positions are accepted using
    /// rejection sampling against an empirical plasma cfDNA 4-mer end motif
    /// frequency table. Fragments whose 5' end motif has low frequency in
    /// plasma are rejected and resampled.
    #[serde(default)]
    pub end_motif_model: Option<String>,
    /// Explicit circulating tumour DNA fraction for the cfDNA fragment model.
    ///
    /// When set, this overrides the purity-derived ctDNA fraction. Must be in
    /// [0.0, 1.0]. A value of 0.05 means 5% of fragments are tumour-derived
    /// and will be drawn from the shorter ctDNA distribution (~143 bp) rather
    /// than the nucleosomal peaks.
    #[serde(default)]
    pub ctdna_fraction: Option<f64>,
    /// Standard deviation of the mononucleosomal fragment peak (bp).
    ///
    /// Defaults to 20 bp. Source: Cristiano et al. 2019 Science (DELFI study)
    /// fragment size distributions in healthy controls and cancer patients.
    #[serde(default)]
    pub mono_sd: Option<f64>,
    /// Standard deviation of the dinucleosomal fragment peak (bp).
    ///
    /// Defaults to 30 bp. Source: Cristiano et al. 2019 Science (DELFI study)
    /// fragment size distributions in healthy controls and cancer patients.
    #[serde(default)]
    pub di_sd: Option<f64>,
}

impl Default for FragmentConfig {
    fn default() -> Self {
        Self {
            model: FragmentModel::Normal,
            mean: default_fragment_mean(),
            sd: default_fragment_sd(),
            long_read: None,
            end_motif_model: None,
            ctdna_fraction: None,
            mono_sd: None,
            di_sd: None,
        }
    }
}

/// Fragment length configuration for long-read platforms (PacBio, Nanopore).
///
/// Fragment lengths are sampled from a log-normal distribution parameterised
/// by the linear-space mean and standard deviation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LongReadFragmentConfig {
    /// Mean fragment length in base pairs.
    #[serde(default = "default_lr_mean")]
    pub mean: usize,
    /// Standard deviation of fragment length in base pairs.
    #[serde(default = "default_lr_sd")]
    pub sd: usize,
    /// Minimum fragment length in base pairs.
    #[serde(default = "default_lr_min")]
    pub min_len: usize,
    /// Maximum fragment length in base pairs.
    #[serde(default = "default_lr_max")]
    pub max_len: usize,
}

fn default_lr_mean() -> usize {
    15000
}
fn default_lr_sd() -> usize {
    5000
}
fn default_lr_min() -> usize {
    1000
}
fn default_lr_max() -> usize {
    100000
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FragmentModel {
    Normal,
    Cfda,
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
    /// Enable microsatellite instability (MSI) mode.
    ///
    /// When true, indel rates at homopolymer and dinucleotide repeat loci are
    /// elevated to simulate MSI-high tumours.
    #[serde(default)]
    pub msi: bool,
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
    /// SV signature to apply: `"HRD"`, `"TDP"`, or `"CHROMOTHRIPSIS"`.
    ///
    /// When set, generates SVs with the pattern characteristic of the named
    /// biological phenotype.
    #[serde(default)]
    pub sv_signature: Option<String>,
    /// Number of SVs to generate for the configured signature (default: 10).
    #[serde(default = "default_sv_count")]
    pub sv_count: usize,
    /// When true, inject canonical driver mutations from the active cancer
    /// preset into the variant list.  Only mutations with fully specified
    /// genomic coordinates are injected; structural or fusion events with
    /// missing allele fields are silently skipped.
    ///
    /// Set automatically to `true` when a `cancer:` preset is applied.
    /// Has no effect when `preset` is not a cancer preset.
    #[serde(default)]
    pub include_driver_mutations: bool,
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
    /// COSMIC SBS signature name for weighted alt base selection (e.g. "SBS1").
    /// When set, replaces uniform base selection with signature-weighted probabilities.
    #[serde(default)]
    pub signature: Option<String>,
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
    /// Sequencing mode: `"panel"` (default) or `"amplicon"`.
    ///
    /// In amplicon mode, fragments exactly span each target region; no
    /// off-target reads are generated.
    #[serde(default = "default_capture_mode")]
    pub mode: String,
    /// Number of bases to soft-clip from each read end when in amplicon mode.
    ///
    /// Simulates primer trimming: the first and last `primer_trim` bases of each
    /// read are removed from the sequence and quality strings before output.
    /// Only applied when `mode == "amplicon"` and `primer_trim > 0`.
    #[serde(default)]
    pub primer_trim: usize,
    /// Target coefficient of variation for per-target coverage.
    ///
    /// When set, the simulation emits a `tracing::warn!` if the achieved CV
    /// exceeds this value. Used to validate uniformity requirements (e.g.
    /// Twist recommends CV ≤ 0.25).
    #[serde(default)]
    pub coverage_cv_target: Option<f64>,
    /// Target on-target fraction (0.0–1.0).
    ///
    /// When set, the simulation emits a `tracing::warn!` if the achieved
    /// on-target fraction falls below this value. Twist panels expect > 0.95.
    #[serde(default)]
    pub on_target_fraction_target: Option<f64>,
}

impl Default for CaptureConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            targets_bed: None,
            off_target_fraction: default_off_target_fraction(),
            coverage_uniformity: default_coverage_uniformity(),
            edge_dropoff_bases: default_edge_dropoff_bases(),
            mode: default_capture_mode(),
            primer_trim: 0,
            coverage_cv_target: None,
            on_target_fraction_target: None,
        }
    }
}

/// Configuration for germline variant simulation.
///
/// Germline variants are distributed across the simulated regions at the
/// given density per kilobase pair. They are assigned VAF 0.5 (heterozygous)
/// or 1.0 (homozygous). When `vcf` is set, the random generator is bypassed
/// and the listed variants are used instead.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GermlineConfig {
    /// Heterozygous SNP density per kbp (default: 0.6).
    #[serde(default = "default_het_snp_density")]
    pub het_snp_density: f64,
    /// Homozygous SNP density per kbp (default: 0.3).
    #[serde(default = "default_hom_snp_density")]
    pub hom_snp_density: f64,
    /// Heterozygous indel density per kbp (default: 0.05).
    #[serde(default = "default_het_indel_density")]
    pub het_indel_density: f64,
    /// Optional VCF of germline variants to use instead of random generation.
    #[serde(default)]
    pub vcf: Option<std::path::PathBuf>,
}

impl Default for GermlineConfig {
    fn default() -> Self {
        Self {
            het_snp_density: default_het_snp_density(),
            hom_snp_density: default_hom_snp_density(),
            het_indel_density: default_het_indel_density(),
            vcf: None,
        }
    }
}

fn default_het_snp_density() -> f64 {
    0.6
}
fn default_hom_snp_density() -> f64 {
    0.3
}
fn default_het_indel_density() -> f64 {
    0.05
}

/// Configuration for paired tumour-normal simulation mode.
///
/// When present in a config file, VarForge runs a tumour simulation
/// (with somatic and germline variants) and a normal simulation (with
/// germline variants only). Outputs land in `tumour/` and `normal/`
/// sub-directories.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairedConfig {
    /// Coverage for the normal sample (default: 30.0).
    #[serde(default = "default_coverage")]
    pub normal_coverage: f64,
    /// Sample name for the normal output (default: "NORMAL").
    #[serde(default = "default_normal_sample_name")]
    pub normal_sample_name: String,
    /// Tumour contamination fraction in the normal sample (0.0–1.0, default: 0.0).
    #[serde(default)]
    pub tumour_contamination_in_normal: f64,
}

impl Default for PairedConfig {
    fn default() -> Self {
        Self {
            normal_coverage: default_coverage(),
            normal_sample_name: default_normal_sample_name(),
            tumour_contamination_in_normal: 0.0,
        }
    }
}

fn default_normal_sample_name() -> String {
    "NORMAL".to_string()
}

// Default value functions
fn default_true() -> bool {
    true
}
fn default_mapq() -> u8 {
    60
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
fn default_sv_count() -> usize {
    10
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

fn default_capture_mode() -> String {
    "panel".to_string()
}

/// Known sequencing chemistry presets.
///
/// Each variant carries fixed defaults for UMI, fragment, and quality
/// parameters appropriate for that platform. Fields left as `None` in
/// the preset mean the platform does not use UMIs.
#[derive(Debug, Clone, PartialEq)]
pub enum ChemistryPreset {
    TwistUmiDuplex,
    IlluminaWgs,
    IlluminaWes,
    IlluminaCtdna,
    PacbioHifi,
    NanoporeR10,
}

impl ChemistryPreset {
    /// Parse a preset name from a string, returning `None` for unknown names.
    ///
    /// # Examples
    ///
    /// ```
    /// use varforge::io::config::ChemistryPreset;
    ///
    /// assert!(ChemistryPreset::from_name("illumina-wgs").is_some());
    /// assert!(ChemistryPreset::from_name("twist-umi-duplex").is_some());
    /// assert!(ChemistryPreset::from_name("unknown-preset").is_none());
    /// ```
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "twist-umi-duplex" => Some(Self::TwistUmiDuplex),
            "illumina-wgs" => Some(Self::IlluminaWgs),
            "illumina-wes" => Some(Self::IlluminaWes),
            "illumina-ctdna" => Some(Self::IlluminaCtdna),
            "pacbio-hifi" => Some(Self::PacbioHifi),
            "nanopore-r10" => Some(Self::NanoporeR10),
            _ => None,
        }
    }
}

/// Apply preset defaults to a config after deserialisation.
///
/// Only fills fields that are `None` or at their zero/default value.
/// Explicit values from the YAML always take precedence.
pub fn apply_preset(config: &mut Config, preset: &ChemistryPreset) {
    // Default values used to detect "not explicitly set".
    let default_fragment = FragmentConfig::default();
    let default_quality = QualityConfig::default();

    match preset {
        ChemistryPreset::TwistUmiDuplex => {
            fill_fragment(
                config,
                &default_fragment,
                FragmentModel::Normal,
                200.0,
                30.0,
            );
            fill_quality(config, &default_quality, 37);
            fill_umi(config, 8, true, false);
        }
        ChemistryPreset::IlluminaWgs => {
            fill_fragment(
                config,
                &default_fragment,
                FragmentModel::Normal,
                300.0,
                50.0,
            );
            fill_quality(config, &default_quality, 36);
            // No UMI for this chemistry.
        }
        ChemistryPreset::IlluminaWes => {
            fill_fragment(
                config,
                &default_fragment,
                FragmentModel::Normal,
                200.0,
                40.0,
            );
            fill_quality(config, &default_quality, 35);
            // No UMI for this chemistry.
        }
        ChemistryPreset::IlluminaCtdna => {
            fill_fragment(config, &default_fragment, FragmentModel::Cfda, 167.0, 20.0);
            fill_quality(config, &default_quality, 36);
            fill_umi(config, 8, false, false);
        }
        ChemistryPreset::PacbioHifi => {
            fill_fragment(
                config,
                &default_fragment,
                FragmentModel::Normal,
                15000.0,
                5000.0,
            );
            fill_quality(config, &default_quality, 25);
            // No UMI for this chemistry.
        }
        ChemistryPreset::NanoporeR10 => {
            fill_fragment(
                config,
                &default_fragment,
                FragmentModel::Normal,
                20000.0,
                10000.0,
            );
            fill_quality(config, &default_quality, 20);
            // No UMI for this chemistry.
        }
    }
}

/// Fill fragment config fields that are still at their defaults.
///
/// The model is always applied from the preset. The serialised default
/// (`normal`) is identical to an explicit `model: normal` in YAML, so
/// there is no way to distinguish them after deserialisation.
///
/// Mean and SD are only overwritten when they still equal the compiled
/// defaults, indicating the user did not set them.
fn fill_fragment(
    config: &mut Config,
    default: &FragmentConfig,
    model: FragmentModel,
    mean: f64,
    sd: f64,
) {
    config.fragment.model = model;
    if (config.fragment.mean - default.mean).abs() < f64::EPSILON {
        config.fragment.mean = mean;
    }
    if (config.fragment.sd - default.sd).abs() < f64::EPSILON {
        config.fragment.sd = sd;
    }
}

/// Fill quality config fields that are still at their defaults.
fn fill_quality(config: &mut Config, default: &QualityConfig, mean_quality: u8) {
    if config.quality.mean_quality == default.mean_quality {
        config.quality.mean_quality = mean_quality;
    }
}

/// Fill UMI config fields, creating the struct if it does not exist.
///
/// If `config.umi` is `None`, creates a new `UmiConfig` with the given
/// values. If it already exists, only fills fields that are still at
/// their zero/default values.
fn fill_umi(config: &mut Config, length: usize, duplex: bool, inline: bool) {
    if config.umi.is_none() {
        config.umi = Some(UmiConfig {
            length,
            duplex,
            inline,
            pcr_cycles: default_pcr_cycles(),
            family_size_mean: default_family_size_mean(),
            family_size_sd: default_family_size_sd(),
        });
    } else if let Some(umi) = config.umi.as_mut() {
        // `default_umi_length()` is 8, so we cannot distinguish "not set"
        // from "set to 8". We only overwrite zero-valued length.
        if umi.length == 0 {
            umi.length = length;
        }
        // Boolean fields: only fill if still at the false zero-value.
        if !umi.duplex {
            umi.duplex = duplex;
        }
        if !umi.inline {
            umi.inline = inline;
        }
    }
}

pub fn load(path: &Path) -> Result<Config> {
    let contents = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read config file: {}", path.display()))?;
    let mut config: Config = serde_yaml::from_str(&contents)
        .with_context(|| format!("failed to parse config file: {}", path.display()))?;
    if let Some(preset_name) = config.preset.clone() {
        if let Some(preset) = ChemistryPreset::from_name(&preset_name) {
            apply_preset(&mut config, &preset);
        } else {
            anyhow::bail!("unknown chemistry preset: {}", preset_name);
        }
    }
    Ok(config)
}

/// Load a YAML config file, substituting `${key}` placeholders with values from `vars`.
///
/// Unresolved placeholders cause an error. Call this instead of `load` when
/// `--set` values are provided.
pub fn load_with_vars(
    path: &Path,
    vars: &std::collections::HashMap<String, String>,
) -> Result<Config> {
    let raw = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read config file: {}", path.display()))?;
    let substituted = substitute_vars(&raw, vars)?;
    let mut config: Config = serde_yaml::from_str(&substituted)
        .with_context(|| format!("failed to parse config file: {}", path.display()))?;
    if let Some(preset_name) = config.preset.clone() {
        if let Some(preset) = ChemistryPreset::from_name(&preset_name) {
            apply_preset(&mut config, &preset);
        } else {
            anyhow::bail!("unknown chemistry preset: {}", preset_name);
        }
    }
    Ok(config)
}

/// Replace `${key}` placeholders in `text` with values from `vars`.
///
/// Returns an error if a placeholder has no corresponding value.
fn substitute_vars(text: &str, vars: &std::collections::HashMap<String, String>) -> Result<String> {
    let mut result = text.to_string();
    let mut i = 0;
    while let Some(start) = result[i..].find("${") {
        let abs_start = i + start;
        if let Some(end_offset) = result[abs_start + 2..].find('}') {
            let abs_end = abs_start + 2 + end_offset;
            let key = result[abs_start + 2..abs_end].to_string();
            if let Some(val) = vars.get(&key) {
                let placeholder = format!("${{{key}}}");
                result = result.replacen(&placeholder, val, 1);
                // Don't advance i; re-scan from abs_start in case substitution
                // introduced another placeholder.
            } else {
                anyhow::bail!(
                    "config placeholder '${{{key}}}' has no --set value; \
                     supply --set {key}=<value>"
                );
            }
        } else {
            // Unterminated placeholder: skip past the "${".
            i = abs_start + 2;
        }
    }
    Ok(result)
}

/// Parse a region string in `chrom:start-end` format.
///
/// Returns `(chrom, start, end)` on success. Returns an error if the string
/// is missing a colon, missing a dash, or if the coordinates are not valid
/// integers.
///
/// # Examples
///
/// ```
/// use varforge::io::config::parse_region;
///
/// let (chrom, start, end) = parse_region("chr7:55000000-55200000").unwrap();
/// assert_eq!(chrom, "chr7");
/// assert_eq!(start, 55_000_000);
/// assert_eq!(end, 55_200_000);
///
/// // Missing colon is an error.
/// assert!(parse_region("chr1_1000_2000").is_err());
///
/// // Non-integer coordinates are an error.
/// assert!(parse_region("chr1:abc-2000").is_err());
/// ```
pub fn parse_region(s: &str) -> Result<(String, u64, u64)> {
    let colon = s
        .find(':')
        .ok_or_else(|| anyhow::anyhow!("expected 'chrom:start-end' but found no ':' in '{}'", s))?;
    let chrom = s[..colon].to_string();
    let coords = &s[colon + 1..];
    let dash = coords.find('-').ok_or_else(|| {
        anyhow::anyhow!(
            "expected 'chrom:start-end' but found no '-' after ':' in '{}'",
            s
        )
    })?;
    let start: u64 = coords[..dash].parse().map_err(|_| {
        anyhow::anyhow!(
            "start coordinate '{}' in '{}' is not a valid integer",
            &coords[..dash],
            s
        )
    })?;
    let end: u64 = coords[dash + 1..].parse().map_err(|_| {
        anyhow::anyhow!(
            "end coordinate '{}' in '{}' is not a valid integer",
            &coords[dash + 1..],
            s
        )
    })?;
    Ok((chrom, start, end))
}

pub fn validate(config: &Config) -> Result<()> {
    anyhow::ensure!(
        config.reference.exists(),
        "reference file not found: {}",
        config.reference.display()
    );

    if let Some(ref bed) = config.regions_bed {
        anyhow::ensure!(
            bed.exists(),
            "regions_bed file not found: {}",
            bed.display()
        );
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
        config.fragment.sd > 0.0,
        "fragment.sd must be greater than zero, got {}",
        config.fragment.sd
    );

    // T115: mono_sd, di_sd, and long_read.sd must be positive so that
    // Normal::new() and LogNormal::new() succeed.
    if let Some(mono_sd) = config.fragment.mono_sd {
        anyhow::ensure!(
            mono_sd > 0.0,
            "fragment.mono_sd must be greater than zero, got {}",
            mono_sd
        );
    }
    if let Some(di_sd) = config.fragment.di_sd {
        anyhow::ensure!(
            di_sd > 0.0,
            "fragment.di_sd must be greater than zero, got {}",
            di_sd
        );
    }
    if let Some(ref lr) = config.fragment.long_read {
        anyhow::ensure!(
            lr.sd > 0,
            "fragment.long_read.sd must be greater than zero, got {}",
            lr.sd
        );
    }

    if let Some(ctdna_frac) = config.fragment.ctdna_fraction {
        anyhow::ensure!(
            (0.0..=1.0).contains(&ctdna_frac),
            "fragment.ctdna_fraction must be in [0.0, 1.0], got {}",
            ctdna_frac
        );
    }

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

    anyhow::ensure!(
        config.output.mapq <= 254,
        "mapq must be 0-254 (255 is reserved in SAM spec)"
    );

    // T083: warn when reads are longer than fragments (produces artifactual overlapping pairs).
    if config.sample.read_length > config.fragment.mean as usize {
        tracing::warn!(
            read_length = config.sample.read_length,
            fragment_mean = config.fragment.mean,
            "read_length exceeds fragment mean: reads will overlap, \
             producing artifactual paired-end overlap"
        );
    }

    // T084: capture BED path existence check.
    if let Some(ref capture) = config.capture {
        if let Some(ref bed) = capture.targets_bed {
            anyhow::ensure!(
                bed.exists(),
                "capture targets_bed file not found: {}",
                bed.display()
            );
        }
    }

    // T085: UMI length must be less than read length.
    if let Some(umi) = &config.umi {
        anyhow::ensure!(umi.length > 0, "UMI length must be > 0");
        anyhow::ensure!(
            umi.length < config.sample.read_length,
            "UMI length ({}) must be less than read_length ({})",
            umi.length,
            config.sample.read_length
        );

        // T134: inline UMI trimming is not implemented in v1. Reject early
        // so users get a clear message rather than silently incorrect output.
        anyhow::ensure!(
            !umi.inline,
            "inline UMI mode is not supported in v1; set umi.inline to false or omit it"
        );
    }

    // T086: validate copy number region strings.
    if let Some(ref cn_regions) = config.copy_number {
        for cn in cn_regions {
            parse_region(&cn.region).with_context(|| {
                format!(
                    "invalid copy_number region '{}': must be chrom:start-end",
                    cn.region
                )
            })?;
        }
    }

    // T088: mutations VCF path existence check.
    if let Some(ref mutations) = config.mutations {
        if let Some(ref vcf_path) = mutations.vcf {
            anyhow::ensure!(
                vcf_path.exists(),
                "mutations vcf file not found: {}",
                vcf_path.display()
            );
        }
    }

    if let Some(vafs) = &config.vafs {
        for &vaf in vafs {
            anyhow::ensure!(
                vaf > 0.0 && vaf <= 1.0,
                "each VAF in `vafs` must be in (0.0, 1.0], got {}",
                vaf
            );
        }
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

    /// Preset `twist-umi-duplex` fills UMI length, duplex flag, and quality.
    #[test]
    fn test_preset_twist_umi_duplex_fills_defaults() {
        let yaml = r#"
reference: /tmp/ref.fa
output:
  directory: /tmp/out
preset: twist-umi-duplex
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        let umi = cfg.umi.expect("umi should be set by preset");
        assert_eq!(umi.length, 8);
        assert!(umi.duplex);
        assert!(!umi.inline);
        assert_eq!(cfg.fragment.mean, 200.0);
        assert_eq!(cfg.fragment.sd, 30.0);
        assert_eq!(cfg.quality.mean_quality, 37);
    }

    /// Preset `illumina-wgs` applies fragment defaults and leaves UMI as `None`.
    #[test]
    fn test_preset_illumina_wgs_no_umi() {
        let yaml = r#"
reference: /tmp/ref.fa
output:
  directory: /tmp/out
preset: illumina-wgs
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert!(
            cfg.umi.is_none(),
            "illumina-wgs should not create UMI config"
        );
        assert_eq!(cfg.fragment.mean, 300.0);
        assert_eq!(cfg.fragment.sd, 50.0);
        assert_eq!(cfg.quality.mean_quality, 36);
    }

    /// An explicit `umi.length` in the YAML beats the preset default.
    #[test]
    fn test_explicit_umi_length_beats_preset() {
        let yaml = r#"
reference: /tmp/ref.fa
output:
  directory: /tmp/out
preset: twist-umi-duplex
umi:
  length: 12
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        let umi = cfg.umi.expect("umi should be set");
        // Explicit length 12 overrides the preset default of 8.
        assert_eq!(umi.length, 12);
    }

    /// `umi.inline: true` is rejected by validation with a clear error.
    #[test]
    fn test_inline_umi_rejected() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/out
umi:
  length: 8
  inline: true
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        let err = validate(&cfg).unwrap_err();
        assert!(
            err.to_string().contains("inline UMI mode is not supported"),
            "unexpected error message: {err}"
        );
    }

    /// Explicit `fragment.ctdna_fraction` is parsed and validated correctly.
    #[test]
    fn test_ctdna_fraction_field_parses() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/out
fragment:
  model: cfda
  mean: 167.0
  sd: 20.0
  ctdna_fraction: 0.03
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert_eq!(cfg.fragment.ctdna_fraction, Some(0.03));
    }

    /// `fragment.ctdna_fraction` outside [0.0, 1.0] fails validation.
    #[test]
    fn test_ctdna_fraction_out_of_range_fails() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/out
fragment:
  ctdna_fraction: 1.5
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert!(validate(&cfg).is_err());
    }

    /// `fragment.mono_sd` and `fragment.di_sd` are parsed correctly.
    #[test]
    fn test_mono_di_sd_fields_parse() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/out
fragment:
  model: cfda
  mean: 167.0
  sd: 20.0
  mono_sd: 15.0
  di_sd: 25.0
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert_eq!(cfg.fragment.mono_sd, Some(15.0));
        assert_eq!(cfg.fragment.di_sd, Some(25.0));
    }

    // T063: Edge-case tests for substitute_vars.

    #[test]
    fn test_substitute_vars_unknown_key_errors() {
        let vars = std::collections::HashMap::new();
        let result = substitute_vars("prefix_${unknown}", &vars);
        assert!(result.is_err(), "unknown key should return an error");
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("unknown"), "error should name the missing key");
    }

    #[test]
    fn test_substitute_vars_empty_value() {
        let mut vars = std::collections::HashMap::new();
        vars.insert("key".to_string(), String::new());
        let result = substitute_vars("prefix_${key}_suffix", &vars).unwrap();
        assert_eq!(result, "prefix__suffix");
    }

    #[test]
    fn test_substitute_vars_adjacent_placeholders() {
        let mut vars = std::collections::HashMap::new();
        vars.insert("a".to_string(), "hello".to_string());
        vars.insert("b".to_string(), "world".to_string());
        let result = substitute_vars("${a}${b}", &vars).unwrap();
        assert_eq!(result, "helloworld");
    }

    #[test]
    fn test_substitute_vars_unterminated_passes_through() {
        // An unterminated placeholder ("${key" with no '}') is skipped over.
        let vars = std::collections::HashMap::new();
        let result = substitute_vars("prefix_${key", &vars).unwrap();
        assert_eq!(
            result, "prefix_${key",
            "unterminated placeholder should pass through"
        );
    }

    #[test]
    fn test_substitute_vars_no_placeholders() {
        let vars = std::collections::HashMap::new();
        let result = substitute_vars("plain string", &vars).unwrap();
        assert_eq!(result, "plain string");
    }

    // T115: zero standard deviation must fail validation.

    /// `fragment.sd: 0.0` fails validation with a clear message.
    #[test]
    fn test_fragment_sd_zero_fails_validation() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/out
fragment:
  sd: 0.0
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        let err = validate(&cfg).unwrap_err();
        assert!(
            err.to_string().contains("greater than zero"),
            "error should mention 'greater than zero', got: {}",
            err
        );
    }

    /// `fragment.mono_sd: 0.0` fails validation.
    #[test]
    fn test_fragment_mono_sd_zero_fails_validation() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/out
fragment:
  sd: 20.0
  mono_sd: 0.0
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert!(
            validate(&cfg).is_err(),
            "mono_sd: 0.0 should fail validation"
        );
    }

    /// `fragment.di_sd: 0.0` fails validation.
    #[test]
    fn test_fragment_di_sd_zero_fails_validation() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/out
fragment:
  sd: 20.0
  di_sd: 0.0
"#;
        let f = write_yaml(yaml);
        let cfg = load(f.path()).unwrap();
        assert!(validate(&cfg).is_err(), "di_sd: 0.0 should fail validation");
    }
}
