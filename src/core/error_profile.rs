//! Empirical error profiles loaded from JSON.
//!
//! Profiles capture position-specific quality distributions, substitution
//! matrices, and context-dependent effects (e.g., poly-G quality drop on
//! NovaSeq) from real sequencing runs.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::core::quality::{ParametricQualityModel, QualityModel};

// ---------------------------------------------------------------------------
// On-disk profile format
// ---------------------------------------------------------------------------

/// Raw JSON structure for a profile file.
///
/// `quality_distribution.read1` is a list of per-position entries. Each
/// position entry is itself a list of `[quality, weight]` pairs, e.g.:
///
/// ```json
/// "read1": [
///   [[37, 0.8], [36, 0.1], [35, 0.1]],   // position 0
///   [[36, 0.7], [35, 0.2], [34, 0.1]],   // position 1
/// ]
/// ```
#[derive(Debug, Deserialize, Serialize)]
pub struct ProfileJson {
    pub platform: Option<String>,
    pub read_length: usize,
    pub quality_distribution: QualityDistributionJson,
    /// Substitution probabilities keyed by `"A>C"` style strings.
    #[serde(default)]
    pub substitution_matrix: HashMap<String, f64>,
    #[serde(default)]
    pub context_effects: HashMap<String, ContextEffectJson>,
}

/// Per-position quality distributions for each read in a pair.
///
/// Each element of `read1` / `read2` corresponds to one read position and
/// contains a list of `[quality_value, weight]` pairs.
#[derive(Debug, Deserialize, Serialize)]
pub struct QualityDistributionJson {
    /// `read1[pos]` is a list of `[quality_value, weight]` pairs.
    pub read1: Vec<Vec<[f64; 2]>>,
    /// Optional per-position list for read 2. Falls back to `read1` when absent.
    #[serde(default)]
    pub read2: Option<Vec<Vec<[f64; 2]>>>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct ContextEffectJson {
    /// Number of Phred points to subtract from quality when context matches.
    pub quality_penalty: u8,
}

// ---------------------------------------------------------------------------
// Internal representation
// ---------------------------------------------------------------------------

/// A context-dependent quality penalty.
#[derive(Debug, Clone)]
pub struct ContextEffect {
    /// Quality score penalty (subtracted from Phred score, floored at 2).
    pub quality_penalty: u8,
}

// ---------------------------------------------------------------------------
// EmpiricalQualityModel
// ---------------------------------------------------------------------------

/// Quality and error model derived from a real sequencing run.
pub struct EmpiricalQualityModel {
    /// Per-position CDF tables for read 1.
    ///
    /// `quality_distributions[pos]` is a list of `(quality, cumulative_weight)`
    /// pairs sorted ascending by cumulative weight.
    quality_distributions: Vec<Vec<(u8, f64)>>,
    /// `(from_base, to_base)` → substitution probability.
    substitution_matrix: HashMap<(u8, u8), f64>,
    /// Context sequences (as byte vecs) → quality penalty.
    context_effects: HashMap<Vec<u8>, ContextEffect>,
    /// Original platform string (for diagnostics).
    // Stored for future diagnostic output; not yet read by any caller.
    #[allow(dead_code)]
    pub platform: Option<String>,
}

impl EmpiricalQualityModel {
    /// Load an empirical model from a JSON profile file.
    pub fn from_file(path: &Path) -> Result<Self> {
        let contents = std::fs::read_to_string(path)
            .with_context(|| format!("cannot read profile: {}", path.display()))?;
        Self::from_json_str(&contents)
    }

    /// Parse a model from a JSON string (useful in tests).
    pub fn from_json_str(json: &str) -> Result<Self> {
        let profile: ProfileJson =
            serde_json::from_str(json).context("failed to parse error profile JSON")?;
        Self::from_profile_json(profile)
    }

    fn from_profile_json(profile: ProfileJson) -> Result<Self> {
        let raw_read1 = &profile.quality_distribution.read1;
        anyhow::ensure!(!raw_read1.is_empty(), "quality_distribution.read1 is empty");

        let quality_distributions = Self::build_distributions(raw_read1)?;

        // Parse "A>C" style keys into (from_byte, to_byte) pairs.
        let mut substitution_matrix: HashMap<(u8, u8), f64> =
            HashMap::with_capacity(profile.substitution_matrix.len());
        for (key, prob) in &profile.substitution_matrix {
            let (from, to) = Self::parse_substitution_key(key)?;
            substitution_matrix.insert((from, to), *prob);
        }

        let mut context_effects: HashMap<Vec<u8>, ContextEffect> =
            HashMap::with_capacity(profile.context_effects.len());
        for (ctx_str, effect_json) in profile.context_effects {
            context_effects.insert(
                ctx_str.into_bytes(),
                ContextEffect {
                    quality_penalty: effect_json.quality_penalty,
                },
            );
        }

        Ok(Self {
            quality_distributions,
            substitution_matrix,
            context_effects,
            platform: profile.platform,
        })
    }

    /// Convert raw per-position `[[quality, weight], ...]` lists into CDF tables.
    fn build_distributions(raw: &[Vec<[f64; 2]>]) -> Result<Vec<Vec<(u8, f64)>>> {
        raw.iter()
            .enumerate()
            .map(|(pos, pairs)| {
                anyhow::ensure!(
                    !pairs.is_empty(),
                    "quality_distribution entry at position {} is empty",
                    pos
                );
                let mut entries: Vec<(u8, f64)> =
                    pairs.iter().map(|pair| (pair[0] as u8, pair[1])).collect();

                let total: f64 = entries.iter().map(|(_, w)| w).sum();
                anyhow::ensure!(
                    total > 0.0,
                    "quality distribution weights must sum to > 0 at position {}",
                    pos
                );
                let mut cumulative = 0.0_f64;
                for (_, w) in &mut entries {
                    cumulative += *w / total;
                    *w = cumulative;
                }
                // Guarantee the last entry is exactly 1.0 (guard f64 rounding).
                if let Some(last) = entries.last_mut() {
                    last.1 = 1.0;
                }
                Ok(entries)
            })
            .collect()
    }

    /// Parse `"A>C"` into `(b'A', b'C')`.
    fn parse_substitution_key(key: &str) -> Result<(u8, u8)> {
        let bytes = key.as_bytes();
        anyhow::ensure!(
            bytes.len() == 3 && bytes[1] == b'>',
            "substitution key must be in 'X>Y' format, got: {}",
            key
        );
        Ok((bytes[0], bytes[2]))
    }

    /// Sample a quality score from the CDF table for `pos`.
    ///
    /// Wraps to the last position if `pos` exceeds the table length.
    fn sample_quality_at<R: Rng>(&self, pos: usize, rng: &mut R) -> u8 {
        let table_pos = pos.min(self.quality_distributions.len() - 1);
        let cdf = &self.quality_distributions[table_pos];
        let r: f64 = rng.gen();
        for &(q, cum) in cdf {
            if r <= cum {
                return q;
            }
        }
        cdf.last().map(|&(q, _)| q).unwrap_or(30)
    }

    /// Sum of all substitution probabilities from a given base.
    fn total_substitution_prob(&self, from: u8) -> f64 {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        BASES
            .iter()
            .filter(|&&b| b != from)
            .map(|&to| {
                self.substitution_matrix
                    .get(&(from, to))
                    .copied()
                    .unwrap_or(0.0)
            })
            .sum()
    }

    /// Choose a substitution target for `from` according to the loaded matrix.
    fn sample_substitution<R: Rng>(&self, from: u8, rng: &mut R) -> u8 {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let total = self.total_substitution_prob(from);
        if total <= 0.0 {
            return from;
        }
        let r: f64 = rng.gen::<f64>() * total;
        let mut cumulative = 0.0;
        for &to in BASES.iter().filter(|&&b| b != from) {
            cumulative += self
                .substitution_matrix
                .get(&(from, to))
                .copied()
                .unwrap_or(0.0);
            if r <= cumulative {
                return to;
            }
        }
        *BASES.iter().find(|&&b| b != from).unwrap_or(&b'A')
    }

    /// Apply context-effect penalties to a quality array in-place.
    ///
    /// Scans `seq` for every registered context and subtracts the penalty
    /// from the quality of the base immediately following the context.
    pub fn apply_context_effects(&self, seq: &[u8], qualities: &mut [u8]) {
        for (ctx_bytes, effect) in &self.context_effects {
            let ctx_len = ctx_bytes.len();
            for end in ctx_len..=seq.len() {
                let window = &seq[end - ctx_len..end];
                if window == ctx_bytes.as_slice() && end < qualities.len() {
                    qualities[end] = qualities[end].saturating_sub(effect.quality_penalty).max(2);
                }
            }
        }
    }
}

impl QualityModel for EmpiricalQualityModel {
    fn generate_qualities<R: Rng>(&self, read_length: usize, rng: &mut R) -> Vec<u8> {
        (0..read_length)
            .map(|pos| self.sample_quality_at(pos, rng))
            .collect()
    }

    fn inject_errors<R: Rng>(&self, sequence: &mut [u8], qualities: &[u8], rng: &mut R) {
        let mut effective = qualities.to_vec();
        self.apply_context_effects(sequence, &mut effective);

        for i in 0..sequence.len() {
            let from = sequence[i];
            let error_prob = ParametricQualityModel::error_probability(effective[i]);

            if rng.gen::<f64>() < error_prob {
                let total_sub = self.total_substitution_prob(from);
                if total_sub > 0.0 {
                    sequence[i] = self.sample_substitution(from, rng);
                } else {
                    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
                    loop {
                        let new_base = BASES[rng.gen_range(0..4)];
                        if new_base != from {
                            sequence[i] = new_base;
                            break;
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Factory: choose model from config
// ---------------------------------------------------------------------------

/// Return an `EmpiricalQualityModel` if `profile_path` is `Some`, else `None`.
///
/// `None` signals the caller to fall back to `ParametricQualityModel`.
// Called only in tests; production code uses `build_empirical_quality` in engine.rs.
#[cfg(test)]
pub fn load_from_config(profile_path: Option<&Path>) -> Result<Option<EmpiricalQualityModel>> {
    match profile_path {
        Some(path) => Ok(Some(EmpiricalQualityModel::from_file(path)?)),
        None => Ok(None),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    use crate::core::quality::{ParametricQualityModel, QualityModel};

    fn minimal_profile_json() -> &'static str {
        r#"{
            "platform": "NovaSeq 6000",
            "read_length": 150,
            "quality_distribution": {
                "read1": [
                    [[37, 0.8], [36, 0.1], [35, 0.1]],
                    [[36, 0.7], [35, 0.2], [34, 0.1]],
                    [[35, 0.6], [34, 0.3], [33, 0.1]]
                ]
            },
            "substitution_matrix": {
                "A>C": 0.001,
                "A>G": 0.003,
                "A>T": 0.0005,
                "C>A": 0.001,
                "C>G": 0.0005,
                "C>T": 0.004,
                "G>A": 0.003,
                "G>C": 0.0005,
                "G>T": 0.001,
                "T>A": 0.0005,
                "T>C": 0.004,
                "T>G": 0.001
            },
            "context_effects": {
                "GGG": { "quality_penalty": 5 }
            }
        }"#
    }

    // Test 1 – load a valid profile
    #[test]
    fn test_load_profile() {
        let model = EmpiricalQualityModel::from_json_str(minimal_profile_json())
            .expect("should parse valid profile");
        assert_eq!(model.platform.as_deref(), Some("NovaSeq 6000"));
        assert_eq!(model.quality_distributions.len(), 3);
        assert!(!model.substitution_matrix.is_empty());
        assert!(model.context_effects.contains_key(b"GGG".as_ref()));
    }

    // Test 2 – reject malformed profiles
    #[test]
    fn test_invalid_profile() {
        // Missing required `read_length` field
        let bad_json = r#"{ "platform": "test", "quality_distribution": { "read1": [] } }"#;
        assert!(
            EmpiricalQualityModel::from_json_str(bad_json).is_err(),
            "should fail on missing read_length"
        );

        // Empty read1 array
        let empty_read1 = r#"{
            "platform": "test",
            "read_length": 100,
            "quality_distribution": { "read1": [] }
        }"#;
        assert!(
            EmpiricalQualityModel::from_json_str(empty_read1).is_err(),
            "should fail on empty read1"
        );

        // Completely malformed JSON
        assert!(EmpiricalQualityModel::from_json_str("not json").is_err());
    }

    // Test 3 – generated qualities match loaded distribution
    #[test]
    fn test_empirical_quality_distribution() {
        // Single position: Q30 50 %, Q20 30 %, Q10 20 %
        let json = r#"{
            "platform": "test",
            "read_length": 1,
            "quality_distribution": {
                "read1": [
                    [[30, 0.5], [20, 0.3], [10, 0.2]]
                ]
            }
        }"#;
        let model = EmpiricalQualityModel::from_json_str(json).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        let n = 10_000usize;
        let mut counts: HashMap<u8, usize> = HashMap::new();
        for _ in 0..n {
            let quals = model.generate_qualities(1, &mut rng);
            *counts.entry(quals[0]).or_insert(0) += 1;
        }

        let freq_30 = *counts.get(&30).unwrap_or(&0) as f64 / n as f64;
        let freq_20 = *counts.get(&20).unwrap_or(&0) as f64 / n as f64;
        let freq_10 = *counts.get(&10).unwrap_or(&0) as f64 / n as f64;

        assert!(
            (freq_30 - 0.50).abs() < 0.04,
            "Q30 frequency {freq_30} not close to 0.50"
        );
        assert!(
            (freq_20 - 0.30).abs() < 0.04,
            "Q20 frequency {freq_20} not close to 0.30"
        );
        assert!(
            (freq_10 - 0.20).abs() < 0.04,
            "Q10 frequency {freq_10} not close to 0.20"
        );
    }

    // Test 4 – error types match loaded substitution matrix proportions
    #[test]
    fn test_substitution_matrix() {
        // Q10 = 10 % error rate; C>T is 4x more likely than C>A.
        let json = r#"{
            "platform": "test",
            "read_length": 1,
            "quality_distribution": {
                "read1": [
                    [[10, 1.0]]
                ]
            },
            "substitution_matrix": {
                "C>T": 0.08,
                "C>A": 0.02
            }
        }"#;
        let model = EmpiricalQualityModel::from_json_str(json).unwrap();
        let mut rng = StdRng::seed_from_u64(99);

        let n = 100_000usize;
        let mut ct_count = 0usize;
        let mut ca_count = 0usize;

        for _ in 0..n {
            let mut seq = vec![b'C'];
            let quals = vec![10u8];
            model.inject_errors(&mut seq, &quals, &mut rng);
            match seq[0] {
                b'T' => ct_count += 1,
                b'A' => ca_count += 1,
                _ => {}
            }
        }

        let total_errors = ct_count + ca_count;
        assert!(total_errors > 0, "expected some errors at Q10");

        let ct_fraction = ct_count as f64 / total_errors as f64;
        assert!(
            (ct_fraction - 0.80).abs() < 0.05,
            "C>T fraction {ct_fraction} not close to 0.80"
        );
    }

    // Test 5 – quality drops in poly-G context
    #[test]
    fn test_context_effects() {
        let json = r#"{
            "platform": "test",
            "read_length": 4,
            "quality_distribution": {
                "read1": [
                    [[30, 1.0]],
                    [[30, 1.0]],
                    [[30, 1.0]],
                    [[30, 1.0]]
                ]
            },
            "context_effects": {
                "GGG": { "quality_penalty": 5 }
            }
        }"#;
        let model = EmpiricalQualityModel::from_json_str(json).unwrap();
        let mut rng = StdRng::seed_from_u64(7);

        let seq = b"GGGA";
        let mut qualities = model.generate_qualities(4, &mut rng);
        assert!(
            qualities.iter().all(|&q| q == 30),
            "all should start at Q30"
        );

        model.apply_context_effects(seq, &mut qualities);

        // Position 3 (A, after GGG) should drop to 25.
        assert_eq!(
            qualities[3], 25,
            "expected Q25 after GGG context, got {}",
            qualities[3]
        );
        assert_eq!(qualities[0], 30);
        assert_eq!(qualities[1], 30);
        assert_eq!(qualities[2], 30);
    }

    // Test 6 – both models implement QualityModel trait
    #[test]
    fn test_trait_compatibility() {
        fn use_model<M: QualityModel>(model: &M, rng: &mut StdRng) {
            let quals = model.generate_qualities(50, rng);
            assert_eq!(quals.len(), 50);
            let mut seq = vec![b'A'; 50];
            model.inject_errors(&mut seq, &quals, rng);
            for &b in &seq {
                assert!(matches!(b, b'A' | b'C' | b'G' | b'T'));
            }
        }

        let mut rng = StdRng::seed_from_u64(1);

        let parametric = ParametricQualityModel::new(36, 0.003);
        use_model(&parametric, &mut rng);

        let empirical = EmpiricalQualityModel::from_json_str(minimal_profile_json()).unwrap();
        use_model(&empirical, &mut rng);
    }

    // Test 7 – no profile_path → parametric model used
    #[test]
    fn test_fallback_to_parametric() {
        let result = load_from_config(None).expect("should not error with no path");
        assert!(
            result.is_none(),
            "expected None when no profile path is set"
        );

        let model = ParametricQualityModel::new(36, 0.003);
        let mut rng = StdRng::seed_from_u64(42);
        let quals = model.generate_qualities(100, &mut rng);
        assert_eq!(quals.len(), 100);
        let mut seq = vec![b'A'; 100];
        model.inject_errors(&mut seq, &quals, &mut rng);
        for &b in &seq {
            assert!(matches!(b, b'A' | b'C' | b'G' | b'T'));
        }
    }
}
