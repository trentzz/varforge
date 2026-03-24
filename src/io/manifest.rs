//! Simulation manifest: records config, reference metadata, output file paths, and run statistics
//! as a machine-readable JSON file for reproducibility and auditing.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::io::config::Config;

/// Information about the reference genome used in the simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceInfo {
    pub path: PathBuf,
    pub genome_size: u64,
    pub chromosomes: usize,
}

/// Statistics collected during simulation.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SimulationStatistics {
    pub total_read_pairs: u64,
    pub total_bases: u64,
    pub variants_spiked: u64,
    pub variants_by_type: HashMap<String, u64>,
    pub mean_coverage_achieved: f64,
    pub regions_simulated: usize,
    pub wall_time_seconds: f64,
    pub reads_per_second: f64,
}

/// Machine-readable record of everything produced by a simulation run.
///
/// The manifest is serialised as pretty-printed JSON and written to the
/// output directory after every successful simulation so that users can
/// reproduce or audit any dataset.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Manifest {
    pub varforge_version: String,
    pub timestamp: String,
    pub seed: Option<u64>,
    pub config: Config,
    pub reference: ReferenceInfo,
    pub output_files: HashMap<String, PathBuf>,
    pub statistics: SimulationStatistics,
}

impl Manifest {
    /// Create a new manifest from the resolved config and reference metadata.
    ///
    /// The timestamp is set to the current UTC time in ISO 8601 format.
    /// `output_files` and `statistics` start empty and are populated via
    /// [`add_output_file`] and [`set_statistics`] respectively.
    pub fn new(config: Config, reference_info: ReferenceInfo) -> Self {
        let seed = config.seed;
        let timestamp = current_timestamp_utc();
        Self {
            varforge_version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp,
            seed,
            config,
            reference: reference_info,
            output_files: HashMap::new(),
            statistics: SimulationStatistics::default(),
        }
    }

    /// Record an output file path under the given key (e.g. `"fastq_r1"`).
    pub fn add_output_file(&mut self, key: &str, path: &Path) {
        self.output_files
            .insert(key.to_string(), path.to_path_buf());
    }

    /// Replace the statistics section with the provided values.
    pub fn set_statistics(&mut self, stats: SimulationStatistics) {
        self.statistics = stats;
    }

    /// Serialise the manifest to pretty-printed JSON and write it to `path`.
    // Called only in tests; production code serialises via serde_json directly.
    #[allow(dead_code)]
    pub fn write(&self, path: &Path) -> Result<()> {
        let json =
            serde_json::to_string_pretty(self).context("failed to serialise manifest to JSON")?;
        std::fs::write(path, json)
            .with_context(|| format!("failed to write manifest to {}", path.display()))?;
        Ok(())
    }
}

/// Return the current wall-clock time as an ISO 8601 UTC string.
///
/// We avoid pulling in the `chrono` or `time` crates just for this; instead
/// we format the UNIX timestamp ourselves to produce a valid RFC 3339 string
/// (e.g. `"2025-03-16T14:30:00Z"`).
fn current_timestamp_utc() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);

    // Simple manual UTC conversion (no leap-second awareness needed here).
    let mut days = secs / 86400;
    let time_of_day = secs % 86400;
    let hh = time_of_day / 3600;
    let mm = (time_of_day % 3600) / 60;
    let ss = time_of_day % 60;

    // Determine year, month, day from `days` (days since 1970-01-01).
    let mut year = 1970u32;
    loop {
        let days_in_year = if is_leap(year) { 366 } else { 365 };
        if days < days_in_year {
            break;
        }
        days -= days_in_year;
        year += 1;
    }
    let months = [
        31u64,
        if is_leap(year) { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    let mut month = 1u32;
    for &m in &months {
        if days < m {
            break;
        }
        days -= m;
        month += 1;
    }
    let day = days + 1;

    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        year, month, day, hh, mm, ss
    )
}

fn is_leap(year: u32) -> bool {
    (year % 4 == 0 && year % 100 != 0) || year % 400 == 0
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use tempfile::TempDir;

    use crate::io::config::{
        Config, FragmentConfig, FragmentModel, OutputConfig, QualityConfig, SampleConfig,
    };

    /// Build a minimal Config that does not require files to exist on disk.
    fn make_config() -> Config {
        Config {
            reference: PathBuf::from("/tmp/ref.fa"),
            output: OutputConfig {
                directory: PathBuf::from("/tmp/out"),
                fastq: true,
                bam: false,
                truth_vcf: true,
                manifest: true,
                germline_vcf: false,
                single_read_bam: false,
                mapq: 60,
            },
            sample: SampleConfig {
                name: "SAMPLE".to_string(),
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
            tumour: None,
            mutations: None,
            umi: None,
            artifacts: None,
            copy_number: None,
            gc_bias: None,
            samples: None,
            capture: None,
            seed: Some(42),
            threads: None,
            chromosomes: None,
            regions_bed: None,
            performance: Default::default(),
            preset: None,
            vafs: None,
            germline: None,
            paired: None,
        }
    }

    fn make_reference_info() -> ReferenceInfo {
        ReferenceInfo {
            path: PathBuf::from("/tmp/ref.fa"),
            genome_size: 3_088_286_401,
            chromosomes: 25,
        }
    }

    // ------------------------------------------------------------------
    // Test 1: Valid manifest with all fields populated
    // ------------------------------------------------------------------
    #[test]
    fn test_manifest_creation() {
        let cfg = make_config();
        let ref_info = make_reference_info();

        let mut manifest = Manifest::new(cfg, ref_info.clone());

        // version and seed
        assert_eq!(manifest.varforge_version, env!("CARGO_PKG_VERSION"));
        assert_eq!(manifest.seed, Some(42));

        // reference round-trips correctly
        assert_eq!(manifest.reference.genome_size, 3_088_286_401);
        assert_eq!(manifest.reference.chromosomes, 25);

        // add an output file
        manifest.add_output_file("fastq_r1", Path::new("output/sample_R1.fastq.gz"));
        assert!(manifest.output_files.contains_key("fastq_r1"));

        // set statistics
        let stats = SimulationStatistics {
            total_read_pairs: 450_000_000,
            variants_spiked: 5_000,
            wall_time_seconds: 2_400.0,
            ..Default::default()
        };
        manifest.set_statistics(stats);

        assert_eq!(manifest.statistics.total_read_pairs, 450_000_000);
        assert_eq!(manifest.statistics.variants_spiked, 5_000);
    }

    // ------------------------------------------------------------------
    // Test 2: Output is parseable JSON
    // ------------------------------------------------------------------
    #[test]
    fn test_manifest_json_valid() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("manifest.json");

        let cfg = make_config();
        let ref_info = make_reference_info();
        let manifest = Manifest::new(cfg, ref_info);

        manifest.write(&path).expect("write should succeed");

        let raw = std::fs::read_to_string(&path).expect("file should exist");
        let parsed: serde_json::Value =
            serde_json::from_str(&raw).expect("output must be valid JSON");

        // Top-level keys must be present
        assert!(parsed.get("varforge_version").is_some());
        assert!(parsed.get("timestamp").is_some());
        assert!(parsed.get("config").is_some());
        assert!(parsed.get("reference").is_some());
        assert!(parsed.get("output_files").is_some());
        assert!(parsed.get("statistics").is_some());
    }

    // ------------------------------------------------------------------
    // Test 3: All statistics fields are present after set_statistics
    // ------------------------------------------------------------------
    #[test]
    fn test_statistics_populated() {
        let cfg = make_config();
        let ref_info = make_reference_info();
        let mut manifest = Manifest::new(cfg, ref_info);

        let mut by_type = HashMap::new();
        by_type.insert("SNV".to_string(), 4_000u64);
        by_type.insert("Indel".to_string(), 750u64);
        by_type.insert("MNV".to_string(), 250u64);

        let stats = SimulationStatistics {
            total_read_pairs: 450_000_000,
            total_bases: 135_000_000_000,
            variants_spiked: 5_000,
            variants_by_type: by_type,
            mean_coverage_achieved: 30.2,
            regions_simulated: 1_500,
            wall_time_seconds: 2_400.0,
            reads_per_second: 187_500.0,
        };

        manifest.set_statistics(stats);

        let s = &manifest.statistics;
        assert_eq!(s.total_read_pairs, 450_000_000);
        assert_eq!(s.total_bases, 135_000_000_000);
        assert_eq!(s.variants_spiked, 5_000);
        assert_eq!(s.variants_by_type["SNV"], 4_000);
        assert_eq!(s.variants_by_type["Indel"], 750);
        assert_eq!(s.variants_by_type["MNV"], 250);
        assert!((s.mean_coverage_achieved - 30.2).abs() < 1e-9);
        assert_eq!(s.regions_simulated, 1_500);
        assert!((s.wall_time_seconds - 2_400.0).abs() < 1e-9);
        assert!((s.reads_per_second - 187_500.0).abs() < 1e-9);
    }

    // ------------------------------------------------------------------
    // Test 4: Full config is serialised inside the manifest
    // ------------------------------------------------------------------
    #[test]
    fn test_config_included() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("manifest.json");

        let cfg = make_config();
        let ref_info = make_reference_info();
        let manifest = Manifest::new(cfg, ref_info);

        manifest.write(&path).unwrap();

        let raw = std::fs::read_to_string(&path).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&raw).unwrap();

        let config_section = parsed.get("config").expect("config key must exist");

        // Verify a few representative fields round-tripped correctly.
        assert_eq!(config_section["sample"]["name"].as_str().unwrap(), "SAMPLE");
        assert_eq!(
            config_section["sample"]["read_length"].as_u64().unwrap(),
            150
        );
        assert!((config_section["sample"]["coverage"].as_f64().unwrap() - 30.0).abs() < 1e-9);
        assert_eq!(config_section["seed"].as_u64().unwrap(), 42);
    }

    // ------------------------------------------------------------------
    // Test 5: Output paths relative to output directory
    // ------------------------------------------------------------------
    #[test]
    fn test_file_paths_relative() {
        let cfg = make_config();
        let ref_info = make_reference_info();
        let mut manifest = Manifest::new(cfg, ref_info);

        // Add relative paths (as a real caller would produce).
        manifest.add_output_file("fastq_r1", Path::new("output/sample_R1.fastq.gz"));
        manifest.add_output_file("fastq_r2", Path::new("output/sample_R2.fastq.gz"));
        manifest.add_output_file("bam", Path::new("output/sample.bam"));
        manifest.add_output_file("truth_vcf", Path::new("output/truth.vcf.gz"));

        assert_eq!(manifest.output_files.len(), 4);

        // Every stored path must be relative (not start with '/').
        for (key, p) in &manifest.output_files {
            assert!(
                p.is_relative(),
                "output file '{}' has an absolute path: {}",
                key,
                p.display()
            );
        }

        // Round-trip through JSON to confirm paths survive serialisation.
        let json = serde_json::to_string(&manifest).unwrap();
        let back: Manifest = serde_json::from_str(&json).unwrap();
        assert_eq!(
            back.output_files["fastq_r1"],
            PathBuf::from("output/sample_R1.fastq.gz")
        );
    }
}
