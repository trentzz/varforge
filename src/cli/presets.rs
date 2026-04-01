//! Built-in simulation presets for common use-cases.
//!
//! Each preset is a *partial* configuration that is applied on top of
//! the YAML defaults and before any CLI flag overrides.  The precedence
//! order is therefore:
//!
//! ```text
//! defaults < preset < YAML config < CLI flags
//! ```
//!
//! A preset is represented as a plain struct so that callers can inspect
//! individual fields without deserialising YAML.

use crate::io::config::{
    ArtifactConfig, CaptureConfig, ContextRuleConfig, FragmentConfig, FragmentModel,
    MutationConfig, RandomMutationConfig, SampleConfig, SequencingErrorConfig, UmiConfig,
};

/// A partial configuration overlay produced by a preset.
///
/// Only the fields that a preset wishes to set are `Some`; all other
/// fields are `None` and leave the base config untouched.
#[derive(Debug, Clone, Default)]
pub struct PresetOverlay {
    pub coverage: Option<f64>,
    pub read_length: Option<usize>,
    pub chromosomes: Option<Vec<String>>,
    pub fragment: Option<FragmentConfig>,
    pub mutations: Option<MutationConfig>,
    pub umi: Option<UmiConfig>,
    pub artifacts: Option<ArtifactConfig>,
    pub purity: Option<f64>,
    pub capture: Option<CaptureConfig>,
    /// Optional sequencing error model configuration to apply as a preset default.
    pub sequencing_errors: Option<SequencingErrorConfig>,
}

/// Return the overlay for a named preset, or an error if the name is unknown.
///
/// Preset names in the `cancer:` namespace (e.g. `cancer:lung_adeno`) are
/// routed to [`crate::cli::cancer_presets::get`].
pub fn get(name: &str) -> anyhow::Result<PresetOverlay> {
    if let Some(cancer_name) = name.strip_prefix("cancer:") {
        let mut overlay = crate::cli::cancer_presets::get(cancer_name)?;
        // Enable driver mutation injection for all cancer presets by default.
        if let Some(ref mut muts) = overlay.mutations {
            muts.include_driver_mutations = true;
        }
        return Ok(overlay);
    }
    match name {
        "small" => Ok(preset_small()),
        "panel" => Ok(preset_panel()),
        "wgs" => Ok(preset_wgs()),
        "cfdna" => Ok(preset_cfdna()),
        "ffpe" => Ok(preset_ffpe()),
        "umi" => Ok(preset_umi()),
        "twist" => Ok(preset_twist()),
        "illumina_novaseq" => Ok(preset_illumina_novaseq()),
        "pacbio_hifi" => Ok(preset_pacbio_hifi()),
        "nanopore_r10" => Ok(preset_nanopore_r10()),
        other => anyhow::bail!(
            "unknown preset '{}'; valid choices: small, panel, wgs, cfdna, ffpe, umi, twist, \
             illumina_novaseq, pacbio_hifi, nanopore_r10, \
             or cancer:<type> (e.g. cancer:lung_adeno)",
            other
        ),
    }
}

/// Return the names of all built-in presets (base presets only; cancer presets
/// are accessible via the `cancer:` prefix — see
/// [`crate::cli::cancer_presets::all_names`]).
pub fn all_names() -> &'static [&'static str] {
    &[
        "small",
        "panel",
        "wgs",
        "cfdna",
        "ffpe",
        "umi",
        "twist",
        "illumina_novaseq",
        "pacbio_hifi",
        "nanopore_r10",
    ]
}

// ---------------------------------------------------------------------------
// Individual preset definitions
// ---------------------------------------------------------------------------

/// `small` – fast smoke-test preset.
///
/// 1× coverage, chr22 only, 100 random mutations.  Designed to finish in
/// ~30 seconds on a laptop.
fn preset_small() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(1.0),
        chromosomes: Some(vec!["chr22".to_string()]),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                count: 100,
                vaf_min: 0.05,
                vaf_max: 0.5,
                snv_fraction: 0.80,
                indel_fraction: 0.15,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_count: 0,
            sv_signature: None,
            include_driver_mutations: false,
        }),
        ..Default::default()
    }
}

/// `panel` – targeted-sequencing preset.
///
/// 500× coverage, UMI enabled, 50 mutations.
fn preset_panel() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(500.0),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                count: 50,
                vaf_min: 0.001,
                vaf_max: 0.5,
                snv_fraction: 0.80,
                indel_fraction: 0.15,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_count: 0,
            sv_signature: None,
            include_driver_mutations: false,
        }),
        umi: Some(UmiConfig {
            length: 8,
            duplex: false,
            pcr_cycles: 10,
            family_size_mean: 3.0,
            family_size_sd: 1.5,
            inline: true,
            spacer: None,
            duplex_conversion_rate: None,
            error_rate: None,
        }),
        ..Default::default()
    }
}

/// `wgs` – whole-genome sequencing preset.
///
/// 30× coverage, 5 000 random mutations across the whole genome.
fn preset_wgs() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                count: 5000,
                vaf_min: 0.001,
                vaf_max: 0.5,
                snv_fraction: 0.80,
                indel_fraction: 0.15,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_count: 0,
            sv_signature: None,
            include_driver_mutations: false,
        }),
        ..Default::default()
    }
}

/// `cfdna` – cell-free DNA / liquid biopsy preset.
///
/// 200× coverage, cfDNA fragment model, low tumour fraction (2 %), UMI enabled.
fn preset_cfdna() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(200.0),
        fragment: Some(FragmentConfig {
            model: FragmentModel::Cfda,
            mean: 167.0,
            sd: 20.0,
            long_read: None,
            end_motif_model: None,
            ctdna_fraction: None,
            mono_sd: None,
            di_sd: None,
        }),
        purity: Some(0.02),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                count: 200,
                vaf_min: 0.001,
                vaf_max: 0.05,
                snv_fraction: 0.80,
                indel_fraction: 0.15,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_count: 0,
            sv_signature: None,
            include_driver_mutations: false,
        }),
        umi: Some(UmiConfig {
            length: 8,
            duplex: true,
            pcr_cycles: 10,
            family_size_mean: 3.0,
            family_size_sd: 1.5,
            inline: false,
            spacer: None,
            duplex_conversion_rate: None,
            error_rate: None,
        }),
        ..Default::default()
    }
}

/// `ffpe` – formalin-fixed paraffin-embedded sample preset.
///
/// 30× coverage with FFPE deamination and oxoG artefacts.
fn preset_ffpe() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        artifacts: Some(ArtifactConfig {
            ffpe_damage_rate: Some(0.02),
            oxog_rate: Some(0.01),
            duplicate_rate: None,
            pcr_error_rate: None,
        }),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                count: 500,
                vaf_min: 0.001,
                vaf_max: 0.5,
                snv_fraction: 0.80,
                indel_fraction: 0.15,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_count: 0,
            sv_signature: None,
            include_driver_mutations: false,
        }),
        ..Default::default()
    }
}

/// `umi` – high-depth duplex UMI sequencing preset.
///
/// 1 000× coverage, duplex UMI mode, panel-scale mutation set.
fn preset_umi() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(1000.0),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                count: 50,
                vaf_min: 0.001,
                vaf_max: 0.5,
                snv_fraction: 0.80,
                indel_fraction: 0.15,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_count: 0,
            sv_signature: None,
            include_driver_mutations: false,
        }),
        umi: Some(UmiConfig {
            length: 9,
            duplex: true,
            pcr_cycles: 12,
            family_size_mean: 4.0,
            family_size_sd: 1.5,
            inline: false,
            spacer: None,
            duplex_conversion_rate: None,
            error_rate: None,
        }),
        ..Default::default()
    }
}

/// `twist` – Twist Biosciences hybrid-capture duplex UMI panel preset.
///
/// 2 000× coverage, 150 bp reads, cfDNA-range fragment sizes. Uses a 5 bp
/// inline dual UMI with a 2 bp `AT` spacer (Twist layout), 90 % duplex
/// conversion rate, and 0.1 % UMI sequencing error rate. VAF range tuned for
/// ultra-low ctDNA detection.
fn preset_twist() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(2000.0),
        read_length: Some(150),
        fragment: Some(FragmentConfig {
            model: FragmentModel::Normal,
            mean: 170.0,
            sd: 30.0,
            long_read: None,
            end_motif_model: None,
            ctdna_fraction: None,
            mono_sd: None,
            di_sd: None,
        }),
        umi: Some(UmiConfig {
            length: 5,
            duplex: true,
            pcr_cycles: 10,
            family_size_mean: 3.5,
            family_size_sd: 1.5,
            inline: true,
            spacer: Some("AT".to_string()),
            duplex_conversion_rate: Some(0.90),
            error_rate: Some(0.001),
        }),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                count: 50,
                vaf_min: 0.0001,
                vaf_max: 0.1,
                snv_fraction: 0.80,
                indel_fraction: 0.15,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_count: 0,
            sv_signature: None,
            include_driver_mutations: false,
        }),
        capture: Some(CaptureConfig {
            enabled: true,
            targets_bed: None,
            off_target_fraction: 0.03,
            coverage_uniformity: 0.15,
            edge_dropoff_bases: 50,
            mode: "panel".to_string(),
            primer_trim: 0,
            coverage_cv_target: Some(0.25),
            on_target_fraction_target: Some(0.95),
        }),
        ..Default::default()
    }
}

/// `illumina_novaseq` – Illumina NovaSeq 6000/X platform preset.
///
/// 30× coverage with realistic NovaSeq sequencing error parameters: 0.1%
/// base error rate, exponential tail rise over the last 20 % of cycles,
/// low sequencing indel rate, GGG and GGC context elevation, and R2
/// asymmetry (1.3× error rate, Q2 quality offset).
fn preset_illumina_novaseq() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        sequencing_errors: Some(SequencingErrorConfig {
            base_error_rate: Some(0.001),
            cycle_error_model: Some("exponential".to_string()),
            tail_start_fraction: Some(0.8),
            tail_rate_multiplier: Some(5.0),
            indel_rate: Some(0.00005),
            indel_insertion_fraction: Some(0.5),
            max_indel_length: Some(1),
            r2_error_multiplier: Some(1.3),
            r2_quality_offset: Some(2),
            context_rules: vec![
                ContextRuleConfig {
                    context: "GGG".to_string(),
                    sub_multiplier: 3.0,
                    indel_multiplier: 2.0,
                },
                ContextRuleConfig {
                    context: "GGC".to_string(),
                    sub_multiplier: 2.0,
                    indel_multiplier: 1.5,
                },
            ],
            kmer_length: Some(3),
            ..Default::default()
        }),
        ..Default::default()
    }
}

/// `pacbio_hifi` – PacBio HiFi (CCS) long-read platform preset.
///
/// 15× coverage with HiFi-typical sequencing parameters: flat 0.01%
/// base error rate, low indel rate biased towards insertions (consistent
/// with PacBio's characteristic error mode).
fn preset_pacbio_hifi() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(15.0),
        sequencing_errors: Some(SequencingErrorConfig {
            base_error_rate: Some(0.0001),
            cycle_error_model: Some("flat".to_string()),
            indel_rate: Some(0.001),
            indel_insertion_fraction: Some(0.6),
            max_indel_length: Some(2),
            ..Default::default()
        }),
        ..Default::default()
    }
}

/// `nanopore_r10` – Oxford Nanopore R10.4 pore long-read platform preset.
///
/// 30× coverage with R10-typical error parameters: higher indel rate
/// biased towards insertions, reflecting R10.4's characteristic deletion
/// and insertion error modes at homopolymers.
fn preset_nanopore_r10() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        sequencing_errors: Some(SequencingErrorConfig {
            base_error_rate: Some(0.005),
            cycle_error_model: Some("flat".to_string()),
            indel_rate: Some(0.03),
            indel_insertion_fraction: Some(0.55),
            max_indel_length: Some(5),
            ..Default::default()
        }),
        ..Default::default()
    }
}

// ---------------------------------------------------------------------------
// Apply a preset overlay onto a Config
// ---------------------------------------------------------------------------

use crate::io::config::Config;

/// Apply `overlay` to `config`, respecting the rule that the overlay takes
/// lower precedence than explicit YAML values.
///
/// Because YAML values are already loaded into `config` at this point we use
/// a *coarse merge*: only fields that the preset touches **and** that still
/// hold their default values are overwritten.  For fields without clear
/// sentinel defaults (e.g. `mutations`) we do not overwrite if the YAML has
/// already provided a value.
///
/// The intent is:  **defaults < preset < YAML config < CLI flags**.
/// Call this function *before* applying CLI flag overrides.
pub fn apply_preset_to_config(config: &mut Config, overlay: &PresetOverlay) {
    use crate::io::config::{default_coverage, default_fragment_mean, default_fragment_sd};

    if let Some(cov) = overlay.coverage {
        // Only apply if coverage is still the compiled default.
        if (config.sample.coverage - default_coverage()).abs() < 1e-9 {
            config.sample.coverage = cov;
        }
    }

    if let Some(rl) = overlay.read_length {
        if config.sample.read_length == SampleConfig::default().read_length {
            config.sample.read_length = rl;
        }
    }

    if let Some(ref chroms) = overlay.chromosomes {
        if config.chromosomes.is_none() {
            config.chromosomes = Some(chroms.clone());
        }
    }

    if let Some(ref frag) = overlay.fragment {
        // Apply if both mean and sd are still at their defaults.
        if (config.fragment.mean - default_fragment_mean()).abs() < 1e-9
            && (config.fragment.sd - default_fragment_sd()).abs() < 1e-9
        {
            config.fragment = frag.clone();
        }
    }

    if let Some(ref muts) = overlay.mutations {
        if config.mutations.is_none() {
            config.mutations = Some(muts.clone());
        }
    }

    if let Some(ref umi) = overlay.umi {
        if config.umi.is_none() {
            config.umi = Some(umi.clone());
        }
    }

    if let Some(ref arts) = overlay.artifacts {
        if config.artifacts.is_none() {
            config.artifacts = Some(arts.clone());
        }
    }

    if let Some(purity) = overlay.purity {
        // Apply preset purity only when no tumour block is defined in YAML.
        if config.tumour.is_none() {
            config.tumour = Some(crate::io::config::TumourConfig {
                purity,
                ploidy: 2,
                clones: Vec::new(),
                msi: false,
            });
        }
    }

    if let Some(ref cap) = overlay.capture {
        if config.capture.is_none() {
            config.capture = Some(cap.clone());
        }
    }

    if let Some(ref se) = overlay.sequencing_errors {
        // Only apply preset defaults when the YAML has not configured sequencing_errors.
        if config.quality.sequencing_errors.is_none() {
            config.quality.sequencing_errors = Some(se.clone());
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_preset_small() {
        let overlay = get("small").unwrap();
        assert_eq!(overlay.coverage, Some(1.0));
        assert_eq!(
            overlay.chromosomes.as_deref(),
            Some(["chr22".to_string()].as_slice())
        );
        let muts = overlay.mutations.unwrap();
        let rand = muts.random.unwrap();
        assert_eq!(rand.count, 100);
    }

    #[test]
    fn test_preset_cfdna() {
        let overlay = get("cfdna").unwrap();
        assert_eq!(overlay.coverage, Some(200.0));
        let frag = overlay.fragment.unwrap();
        assert!(matches!(frag.model, FragmentModel::Cfda));
        assert!((frag.mean - 167.0).abs() < 1e-9);
        // purity should be the low-TF value
        assert!((overlay.purity.unwrap() - 0.02).abs() < 1e-9);
        // UMI duplex should be enabled
        assert!(overlay.umi.unwrap().duplex);
    }

    #[test]
    fn test_all_presets_valid() {
        use crate::io::config::{
            Config, FragmentConfig, OutputConfig, QualityConfig, SampleConfig,
        };
        use std::path::PathBuf;

        let base_config = || Config {
            reference: PathBuf::from("/dev/null"),
            output: OutputConfig {
                directory: PathBuf::from("/tmp"),
                fastq: true,
                bam: false,
                truth_vcf: false,
                manifest: false,
                germline_vcf: false,
                single_read_bam: false,
                mapq: 60,
                annotate_reads: false,
            },
            sample: SampleConfig::default(),
            fragment: FragmentConfig::default(),
            quality: QualityConfig::default(),
            tumour: None,
            mutations: None,
            umi: None,
            artifacts: None,
            seed: None,
            threads: None,
            chromosomes: None,
            regions_bed: None,
            copy_number: None,
            gc_bias: None,
            samples: None,
            capture: None,
            performance: Default::default(),
            preset: None,
            vafs: None,
            germline: None,
            paired: None,
        };

        for name in all_names() {
            let overlay = get(name).expect(name);
            let mut cfg = base_config();
            apply_preset_to_config(&mut cfg, &overlay);
            // Basic sanity: coverage must be positive after applying preset.
            assert!(
                cfg.sample.coverage > 0.0,
                "preset '{name}' gave non-positive coverage"
            );
        }
    }

    #[test]
    fn test_unknown_preset_errors() {
        assert!(get("nonexistent").is_err());
    }

    #[test]
    fn test_preset_precedence_yaml_wins() {
        // When YAML already set a value, the preset must NOT overwrite it.
        use crate::io::config::{
            Config, FragmentConfig, OutputConfig, QualityConfig, SampleConfig,
        };
        use std::path::PathBuf;

        let mut cfg = Config {
            reference: PathBuf::from("/dev/null"),
            output: OutputConfig {
                directory: PathBuf::from("/tmp"),
                fastq: true,
                bam: false,
                truth_vcf: false,
                manifest: false,
                germline_vcf: false,
                single_read_bam: false,
                mapq: 60,
                annotate_reads: false,
            },
            // User set coverage to 60x in YAML.
            sample: SampleConfig {
                coverage: 60.0,
                ..SampleConfig::default()
            },
            fragment: FragmentConfig::default(),
            quality: QualityConfig::default(),
            tumour: None,
            mutations: None,
            umi: None,
            artifacts: None,
            seed: None,
            threads: None,
            chromosomes: None,
            regions_bed: None,
            copy_number: None,
            gc_bias: None,
            samples: None,
            capture: None,
            performance: Default::default(),
            preset: None,
            vafs: None,
            germline: None,
            paired: None,
        };

        let overlay = get("small").unwrap(); // small preset wants 1x
        apply_preset_to_config(&mut cfg, &overlay);

        // YAML value of 60x must survive.
        assert!(
            (cfg.sample.coverage - 60.0).abs() < 1e-9,
            "YAML coverage should not be overwritten by preset"
        );
    }
}
