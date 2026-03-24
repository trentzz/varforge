//! Multi-sample longitudinal simulation support.
//!
//! This module provides types and logic for simulating several samples from a
//! shared clonal architecture while allowing each sample to have its own:
//!
//! - Tumour fraction (ctDNA fraction / purity)
//! - Coverage level
//! - Fragment model
//! - Optional per-clone CCF adjustments (`clonal_shift`)
//!
//! # Usage
//!
//! Build a [`MultiSamplePlan`] from the top-level `Config`, then call
//! [`MultiSamplePlan::per_sample_configs`] to obtain one `Config` per sample.
//! The caller is responsible for running the simulation engine for each config
//! and writing per-sample output directories.

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;

use crate::io::config::{
    CloneConfig, Config, FragmentConfig, FragmentModel, SampleConfig, SampleEntry, TumourConfig,
};

// ---------------------------------------------------------------------------
// Plan types
// ---------------------------------------------------------------------------

/// A fully resolved simulation configuration for a single sample in a
/// multi-sample run.
#[derive(Debug, Clone)]
pub struct ResolvedSample {
    /// Sample name (used as sub-directory and file prefix).
    pub name: String,
    /// The full `Config` to pass to the simulation engine for this sample.
    pub config: Config,
    /// Output directory for this sample's files.
    pub output_dir: PathBuf,
}

/// Describes how a single sample differs from the shared base configuration.
#[derive(Debug, Clone)]
pub struct SampleDelta {
    /// Sample name.
    pub name: String,
    /// Per-sample coverage override.
    pub coverage: f64,
    /// Effective tumour fraction (replaces `tumour.purity`).
    pub tumour_fraction: f64,
    /// Optional fragment model override.
    pub fragment_model: Option<FragmentModel>,
    /// Per-clone CCF adjustments.
    pub clonal_shift: HashMap<String, f64>,
}

impl SampleDelta {
    /// Build a [`SampleDelta`] from a config [`SampleEntry`].
    pub fn from_entry(entry: &SampleEntry) -> Self {
        Self {
            name: entry.name.clone(),
            coverage: entry.coverage,
            tumour_fraction: entry.tumour_fraction,
            fragment_model: entry.fragment_model.clone(),
            clonal_shift: entry.clonal_shift.clone(),
        }
    }
}

// ---------------------------------------------------------------------------
// Plan
// ---------------------------------------------------------------------------

/// Orchestrates multi-sample simulation from a shared base configuration.
pub struct MultiSamplePlan {
    /// Base configuration shared by all samples.
    pub base: Config,
    /// Per-sample deltas derived from `base.samples`.
    pub samples: Vec<SampleDelta>,
}

impl MultiSamplePlan {
    /// Build a plan from a configuration that contains `samples:`.
    ///
    /// Returns `None` when `config.samples` is absent or empty so the caller
    /// can fall back to single-sample mode.
    pub fn from_config(config: Config) -> Option<Self> {
        let entries = config.samples.as_ref()?;
        if entries.is_empty() {
            return None;
        }
        let samples = entries.iter().map(SampleDelta::from_entry).collect();
        Some(Self {
            base: config,
            samples,
        })
    }

    /// Resolve the full `Config` for every sample.
    ///
    /// Each resolved config:
    /// 1. Copies all shared fields from `base`.
    /// 2. Overrides `sample.name`, `sample.coverage`, and `tumour.purity` with
    ///    per-sample values.
    /// 3. Applies `clonal_shift` to clone CCFs.
    /// 4. Overrides `fragment.model` if the sample specifies one.
    /// 5. Points `output.directory` to a sub-directory named after the sample.
    pub fn per_sample_configs(&self) -> Result<Vec<ResolvedSample>> {
        let root_dir = self.base.output.directory.clone();
        let mut out = Vec::with_capacity(self.samples.len());

        for delta in &self.samples {
            let mut cfg = self.base.clone();

            // Clear the samples list so the engine runs in single-sample mode.
            cfg.samples = None;

            // Per-sample output sub-directory.
            let sample_dir = root_dir.join(&delta.name);
            cfg.output.directory = sample_dir.clone();

            // Override sample metadata.
            cfg.sample = SampleConfig {
                name: delta.name.clone(),
                coverage: delta.coverage,
                ..cfg.sample.clone()
            };

            // Override tumour fraction (purity).
            apply_tumour_fraction(&mut cfg, delta.tumour_fraction, &delta.clonal_shift)?;

            // Override fragment model when specified.
            if let Some(ref model) = delta.fragment_model {
                cfg.fragment = FragmentConfig {
                    model: model.clone(),
                    ..cfg.fragment
                };
            }

            out.push(ResolvedSample {
                name: delta.name.clone(),
                config: cfg,
                output_dir: sample_dir,
            });
        }

        Ok(out)
    }

    /// Return the number of samples in this plan.
    // Not yet called from production code; retained as a standard collection API method.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Return `true` if the plan has no samples.
    // Not yet called from production code; retained as a standard collection API method.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Apply the per-sample tumour fraction and optional clonal CCF shifts.
///
/// - Sets `tumour.purity` to `tumour_fraction`.
/// - For each clone ID in `clonal_shift`, updates the matching clone's CCF.
/// - If `cfg.tumour` is absent, creates a minimal `TumourConfig`.
fn apply_tumour_fraction(
    cfg: &mut Config,
    tumour_fraction: f64,
    clonal_shift: &HashMap<String, f64>,
) -> Result<()> {
    match &mut cfg.tumour {
        Some(tumour) => {
            tumour.purity = tumour_fraction;
            for clone in &mut tumour.clones {
                if let Some(&new_ccf) = clonal_shift.get(&clone.id) {
                    anyhow::ensure!(
                        (0.0..=1.0).contains(&new_ccf),
                        "clonal_shift CCF for clone '{}' must be in [0, 1], got {}",
                        clone.id,
                        new_ccf
                    );
                    clone.ccf = new_ccf;
                }
            }
        }
        None => {
            // Build a minimal TumourConfig from the clonal_shift map.
            let clones: Vec<CloneConfig> = clonal_shift
                .iter()
                .map(|(id, &ccf)| CloneConfig {
                    id: id.clone(),
                    ccf,
                    parent: None,
                })
                .collect();
            cfg.tumour = Some(TumourConfig {
                purity: tumour_fraction,
                ploidy: 2,
                clones,
                msi: false,
            });
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Combined manifest
// ---------------------------------------------------------------------------

/// Metadata for a single simulated sample, written into the combined manifest.
#[derive(Debug, Clone, serde::Serialize)]
pub struct SampleManifestEntry {
    pub name: String,
    pub output_dir: String,
    pub coverage: f64,
    pub tumour_fraction: f64,
    pub total_read_pairs: u64,
    pub variants_applied: usize,
}

/// Write a combined manifest JSON file at `<root_dir>/manifest.json`.
pub fn write_combined_manifest(
    root_dir: &std::path::Path,
    entries: &[SampleManifestEntry],
    varforge_version: &str,
) -> Result<()> {
    let manifest_path = root_dir.join("manifest.json");

    let manifest = serde_json::json!({
        "varforge_version": varforge_version,
        "multi_sample": true,
        "sample_count": entries.len(),
        "samples": entries,
    });

    std::fs::write(
        &manifest_path,
        serde_json::to_string_pretty(&manifest)
            .map_err(|e| anyhow::anyhow!("failed to serialize combined manifest: {}", e))?,
    )
    .map_err(|e| {
        anyhow::anyhow!(
            "failed to write combined manifest {}: {}",
            manifest_path.display(),
            e
        )
    })?;

    tracing::info!("combined manifest written to {}", manifest_path.display());
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    use crate::io::config::{
        Config, FragmentConfig, FragmentModel, OutputConfig, QualityConfig, SampleConfig,
    };

    // -----------------------------------------------------------------------
    // Test helpers
    // -----------------------------------------------------------------------

    fn base_config(out_dir: &str) -> Config {
        Config {
            reference: std::path::PathBuf::from("/dev/null"),
            output: OutputConfig {
                directory: std::path::PathBuf::from(out_dir),
                fastq: true,
                bam: false,
                truth_vcf: true,
                manifest: true,
                germline_vcf: false,
                single_read_bam: false,
                mapq: 60,
            },
            sample: SampleConfig {
                name: "BASE".to_string(),
                read_length: 150,
                coverage: 30.0,
                platform: None,
            },
            fragment: FragmentConfig {
                model: FragmentModel::Normal,
                mean: 300.0,
                sd: 50.0,
                long_read: None,
                end_motif_model: None,
                ctdna_fraction: None,
                mono_sd: None,
                di_sd: None,
            },
            quality: QualityConfig {
                mean_quality: 36,
                tail_decay: 0.003,
                profile_path: None,
            },
            tumour: Some(TumourConfig {
                purity: 0.5,
                ploidy: 2,
                clones: vec![
                    CloneConfig {
                        id: "clone_A".to_string(),
                        ccf: 1.0,
                        parent: None,
                    },
                    CloneConfig {
                        id: "clone_B".to_string(),
                        ccf: 0.4,
                        parent: Some("clone_A".to_string()),
                    },
                ],
                msi: false,
            }),
            mutations: None,
            umi: None,
            artifacts: None,
            seed: Some(42),
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
        }
    }

    fn write_yaml(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    // -----------------------------------------------------------------------
    // 1. test_multi_sample_config_parsing
    // -----------------------------------------------------------------------
    #[test]
    fn test_multi_sample_config_parsing() {
        let yaml = r#"
reference: /dev/null
output:
  directory: /tmp/multi
samples:
  - name: "diagnosis"
    coverage: 200
    tumour_fraction: 0.05
    fragment_model: cfda
  - name: "post_treatment"
    coverage: 200
    tumour_fraction: 0.001
    fragment_model: cfda
  - name: "relapse"
    coverage: 200
    tumour_fraction: 0.03
    fragment_model: cfda
    clonal_shift:
      clone_A: 0.8
      clone_B: 0.1
"#;
        let f = write_yaml(yaml);
        let cfg = crate::io::config::load(f.path()).unwrap();

        let entries = cfg.samples.as_ref().expect("samples should be present");
        assert_eq!(entries.len(), 3, "should have 3 sample entries");

        assert_eq!(entries[0].name, "diagnosis");
        assert!((entries[0].coverage - 200.0).abs() < 1e-9);
        assert!((entries[0].tumour_fraction - 0.05).abs() < 1e-9);
        assert!(matches!(
            entries[0].fragment_model,
            Some(FragmentModel::Cfda)
        ));

        assert_eq!(entries[2].name, "relapse");
        assert_eq!(entries[2].clonal_shift.len(), 2);
        assert!((entries[2].clonal_shift["clone_A"] - 0.8).abs() < 1e-9);
        assert!((entries[2].clonal_shift["clone_B"] - 0.1).abs() < 1e-9);
    }

    // -----------------------------------------------------------------------
    // 2. test_shared_mutations
    // -----------------------------------------------------------------------
    /// All resolved configs reference the same mutations config (shared clonal
    /// architecture).  Since we clone from the base config, the mutations field
    /// is identical across all samples.
    #[test]
    fn test_shared_mutations() {
        use crate::io::config::MutationConfig;

        let mut cfg = base_config("/tmp/shared");
        cfg.mutations = Some(MutationConfig {
            vcf: None,
            random: None,
            sv_count: 0,
            sv_signature: None,
            include_driver_mutations: false,
        });
        cfg.samples = Some(vec![
            SampleEntry {
                name: "s1".to_string(),
                coverage: 100.0,
                tumour_fraction: 0.1,
                fragment_model: None,
                clonal_shift: HashMap::new(),
            },
            SampleEntry {
                name: "s2".to_string(),
                coverage: 50.0,
                tumour_fraction: 0.01,
                fragment_model: None,
                clonal_shift: HashMap::new(),
            },
        ]);

        let plan = MultiSamplePlan::from_config(cfg).expect("plan should be created");
        let resolved = plan.per_sample_configs().unwrap();
        assert_eq!(resolved.len(), 2);

        // Both samples share the same mutations config (or lack thereof).
        let m0 = &resolved[0].config.mutations;
        let m1 = &resolved[1].config.mutations;
        assert!(
            m0.is_some() && m1.is_some(),
            "both should have mutations config"
        );
    }

    // -----------------------------------------------------------------------
    // 3. test_different_tumor_fractions
    // -----------------------------------------------------------------------
    /// Each sample's resolved config carries its own purity (= tumour fraction).
    #[test]
    fn test_different_tumor_fractions() {
        let mut cfg = base_config("/tmp/tf_test");
        cfg.samples = Some(vec![
            SampleEntry {
                name: "diagnosis".to_string(),
                coverage: 200.0,
                tumour_fraction: 0.05,
                fragment_model: None,
                clonal_shift: HashMap::new(),
            },
            SampleEntry {
                name: "post_treatment".to_string(),
                coverage: 200.0,
                tumour_fraction: 0.001,
                fragment_model: None,
                clonal_shift: HashMap::new(),
            },
            SampleEntry {
                name: "relapse".to_string(),
                coverage: 200.0,
                tumour_fraction: 0.03,
                fragment_model: None,
                clonal_shift: HashMap::new(),
            },
        ]);

        let plan = MultiSamplePlan::from_config(cfg).expect("plan should be created");
        let resolved = plan.per_sample_configs().unwrap();

        let purity = |idx: usize| {
            resolved[idx]
                .config
                .tumour
                .as_ref()
                .expect("tumour config must be present")
                .purity
        };

        assert!(
            (purity(0) - 0.05).abs() < 1e-9,
            "diagnosis purity should be 0.05"
        );
        assert!(
            (purity(1) - 0.001).abs() < 1e-9,
            "post_treatment purity should be 0.001"
        );
        assert!(
            (purity(2) - 0.03).abs() < 1e-9,
            "relapse purity should be 0.03"
        );

        // All different from each other.
        assert_ne!(purity(0).to_bits(), purity(1).to_bits());
        assert_ne!(purity(1).to_bits(), purity(2).to_bits());
    }

    // -----------------------------------------------------------------------
    // 4. test_clonal_shift
    // -----------------------------------------------------------------------
    /// Clonal CCF adjustments are correctly applied to the matching clones.
    #[test]
    fn test_clonal_shift() {
        let mut cfg = base_config("/tmp/clonal_shift");
        let mut shift = HashMap::new();
        shift.insert("clone_A".to_string(), 0.8);
        shift.insert("clone_B".to_string(), 0.1);

        cfg.samples = Some(vec![SampleEntry {
            name: "relapse".to_string(),
            coverage: 200.0,
            tumour_fraction: 0.03,
            fragment_model: None,
            clonal_shift: shift,
        }]);

        let plan = MultiSamplePlan::from_config(cfg).expect("plan should be created");
        let resolved = plan.per_sample_configs().unwrap();
        assert_eq!(resolved.len(), 1);

        let tumour = resolved[0]
            .config
            .tumour
            .as_ref()
            .expect("tumour must be present");

        let clone_a = tumour
            .clones
            .iter()
            .find(|c| c.id == "clone_A")
            .expect("clone_A must exist");
        let clone_b = tumour
            .clones
            .iter()
            .find(|c| c.id == "clone_B")
            .expect("clone_B must exist");

        assert!(
            (clone_a.ccf - 0.8).abs() < 1e-9,
            "clone_A CCF should be 0.8 after shift"
        );
        assert!(
            (clone_b.ccf - 0.1).abs() < 1e-9,
            "clone_B CCF should be 0.1 after shift"
        );
    }

    // -----------------------------------------------------------------------
    // 5. test_per_sample_output
    // -----------------------------------------------------------------------
    /// Each resolved sample has its own output sub-directory.
    #[test]
    fn test_per_sample_output() {
        let root = "/tmp/per_sample_test";
        let mut cfg = base_config(root);
        cfg.samples = Some(vec![
            SampleEntry {
                name: "diagnosis".to_string(),
                coverage: 200.0,
                tumour_fraction: 0.05,
                fragment_model: None,
                clonal_shift: HashMap::new(),
            },
            SampleEntry {
                name: "relapse".to_string(),
                coverage: 200.0,
                tumour_fraction: 0.03,
                fragment_model: None,
                clonal_shift: HashMap::new(),
            },
        ]);

        let plan = MultiSamplePlan::from_config(cfg).expect("plan should be created");
        let resolved = plan.per_sample_configs().unwrap();

        let dir0 = &resolved[0].output_dir;
        let dir1 = &resolved[1].output_dir;

        // Each sample directory is a sub-directory of root named after the sample.
        assert_eq!(dir0, &std::path::PathBuf::from(root).join("diagnosis"));
        assert_eq!(dir1, &std::path::PathBuf::from(root).join("relapse"));

        // The resolved config's output.directory matches the output_dir field.
        assert_eq!(resolved[0].config.output.directory, *dir0);
        assert_eq!(resolved[1].config.output.directory, *dir1);

        // Directories are distinct.
        assert_ne!(dir0, dir1);
    }

    // -----------------------------------------------------------------------
    // 6. test_combined_manifest
    // -----------------------------------------------------------------------
    /// The combined manifest JSON is written and includes all sample entries.
    #[test]
    fn test_combined_manifest() {
        let tmp = tempfile::TempDir::new().unwrap();
        let root = tmp.path();

        let entries = vec![
            SampleManifestEntry {
                name: "diagnosis".to_string(),
                output_dir: root.join("diagnosis").to_string_lossy().into_owned(),
                coverage: 200.0,
                tumour_fraction: 0.05,
                total_read_pairs: 1_000_000,
                variants_applied: 42,
            },
            SampleManifestEntry {
                name: "relapse".to_string(),
                output_dir: root.join("relapse").to_string_lossy().into_owned(),
                coverage: 200.0,
                tumour_fraction: 0.03,
                total_read_pairs: 1_000_000,
                variants_applied: 42,
            },
        ];

        write_combined_manifest(root, &entries, "0.1.0").unwrap();

        let manifest_path = root.join("manifest.json");
        assert!(manifest_path.exists(), "manifest.json should be written");

        let content = std::fs::read_to_string(&manifest_path).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();

        assert_eq!(parsed["multi_sample"], true);
        assert_eq!(parsed["sample_count"], 2);

        let samples = parsed["samples"]
            .as_array()
            .expect("samples should be an array");
        assert_eq!(samples.len(), 2);
        assert_eq!(samples[0]["name"], "diagnosis");
        assert_eq!(samples[1]["name"], "relapse");
        assert!((samples[0]["tumour_fraction"].as_f64().unwrap() - 0.05).abs() < 1e-9);
    }
}
