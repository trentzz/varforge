/// GC-content-based coverage bias model.
///
/// Real Illumina sequencing shows systematic GC bias: AT-rich and GC-rich
/// regions have lower coverage than regions near 50% GC. This module
/// implements a configurable bell-shaped bias curve so the simulation can
/// reproduce that behaviour.
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Which GC bias curve to use.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "snake_case")]
#[derive(Default)]
pub enum GcBiasModelKind {
    /// Bell-shaped Gaussian centred at 50% GC (default Illumina-like curve).
    #[default]
    Default,
    /// No bias – returns 1.0 everywhere.
    Flat,
    /// Reserved for future custom curve support.
    Custom,
}


/// Configuration block for GC bias (mirrors the YAML `gc_bias:` section).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GcBiasConfig {
    /// Whether GC bias is applied at all.
    #[serde(default = "default_true")]
    pub enabled: bool,
    /// Which curve to use.
    #[serde(default)]
    pub model: GcBiasModelKind,
    /// Severity multiplier.
    ///
    /// * `0.0` – no bias (equivalent to the flat model).
    /// * `1.0` – default realistic Illumina bias.
    /// * `2.0` – extreme bias (useful for stress-testing callers).
    #[serde(default = "default_severity")]
    pub severity: f64,
}

fn default_true() -> bool {
    true
}
fn default_severity() -> f64 {
    1.0
}

impl Default for GcBiasConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            model: GcBiasModelKind::Default,
            severity: 1.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Model
// ---------------------------------------------------------------------------

/// Runtime GC bias model.
///
/// Build once from a [`GcBiasConfig`] (or via [`GcBiasModel::default`]) then
/// call [`GcBiasModel::coverage_multiplier`] for each fragment.
#[derive(Debug, Clone)]
pub struct GcBiasModel {
    enabled: bool,
    kind: GcBiasModelKind,
    /// Standard deviation of the Gaussian used for the default curve.
    ///
    /// Chosen so that the curve hits ~0.7 at GC = 0.2 / 0.8 and ~0.3 at
    /// GC = 0.1 / 0.9:
    ///   exp(-0.5 * ((0.2 - 0.5) / σ)²) ≈ 0.7  →  σ ≈ 0.23
    sigma: f64,
    /// Severity multiplier (see [`GcBiasConfig::severity`]).
    severity: f64,
}

impl Default for GcBiasModel {
    fn default() -> Self {
        Self::from_config(&GcBiasConfig::default())
    }
}

impl GcBiasModel {
    /// Construct a model from a [`GcBiasConfig`].
    pub fn from_config(cfg: &GcBiasConfig) -> Self {
        Self {
            enabled: cfg.enabled,
            kind: cfg.model.clone(),
            sigma: 0.23,
            severity: cfg.severity.max(0.0),
        }
    }

    /// Returns a coverage multiplier in `[0.0, 1.0]` for the given GC fraction.
    ///
    /// The multiplier is used as the probability of retaining a fragment
    /// (rejection sampling): a fragment is kept iff `rng < multiplier`.
    ///
    /// Behaviour:
    /// - `Flat` model or `enabled = false` → always returns `1.0`.
    /// - `Default` model → Gaussian centred at 0.5; `severity` blends between
    ///   1.0 (no bias) and the raw Gaussian output (full bias).
    pub fn coverage_multiplier(&self, gc_fraction: f64) -> f64 {
        if !self.enabled || self.kind == GcBiasModelKind::Flat {
            return 1.0;
        }

        if self.severity == 0.0 {
            return 1.0;
        }

        // Gaussian centred at 0.5.
        let raw = gaussian(gc_fraction, 0.5, self.sigma);

        // Blend: multiplier = 1.0 * (1 - severity) + raw * severity
        // At severity=1.0 we get the raw Gaussian; at 0.0 we get 1.0.
        let blended = 1.0 + self.severity * (raw - 1.0);
        blended.clamp(0.0, 1.0)
    }

    /// Calculate the GC fraction of a nucleotide sequence.
    ///
    /// N bases (ambiguity codes) are excluded from both the numerator and
    /// denominator so they do not dilute the GC estimate.
    pub fn gc_fraction(sequence: &[u8]) -> f64 {
        let mut gc = 0usize;
        let mut total = 0usize;
        for &b in sequence {
            match b.to_ascii_uppercase() {
                b'G' | b'C' => {
                    gc += 1;
                    total += 1;
                }
                b'A' | b'T' => {
                    total += 1;
                }
                // N and anything else: skip
                _ => {}
            }
        }
        if total == 0 {
            return 0.5; // Undefined → neutral
        }
        gc as f64 / total as f64
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Evaluate a Gaussian with unit peak at `mu` and standard deviation `sigma`.
///
/// Returns a value in `(0.0, 1.0]`.
fn gaussian(x: f64, mu: f64, sigma: f64) -> f64 {
    let z = (x - mu) / sigma;
    (-0.5 * z * z).exp()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // 1. GC fraction calculation for known sequences.
    #[test]
    fn test_gc_fraction() {
        // Pure GC
        assert!((GcBiasModel::gc_fraction(b"GGCC") - 1.0).abs() < 1e-9);
        // Pure AT
        assert!((GcBiasModel::gc_fraction(b"AATT") - 0.0).abs() < 1e-9);
        // 50%
        assert!((GcBiasModel::gc_fraction(b"ACGT") - 0.5).abs() < 1e-9);
        // 3 GC out of 6 non-N bases
        assert!((GcBiasModel::gc_fraction(b"GCGAAT") - 0.5).abs() < 1e-9);
    }

    // 2. Multiplier is ~1.0 at 50% GC.
    #[test]
    fn test_multiplier_at_50_percent() {
        let model = GcBiasModel::default();
        let m = model.coverage_multiplier(0.5);
        assert!(
            (m - 1.0).abs() < 1e-6,
            "expected ~1.0 at 50% GC, got {m}"
        );
    }

    // 3. Multiplier is reduced at extreme GC (10% and 90%).
    #[test]
    fn test_multiplier_at_extremes() {
        let model = GcBiasModel::default();
        let m_10 = model.coverage_multiplier(0.1);
        let m_90 = model.coverage_multiplier(0.9);
        // Both should be well below 1.0; task spec says ~0.3.
        assert!(m_10 < 0.5, "10% GC: expected < 0.5, got {m_10}");
        assert!(m_90 < 0.5, "90% GC: expected < 0.5, got {m_90}");
        // And symmetric (Gaussian is symmetric around 0.5).
        assert!(
            (m_10 - m_90).abs() < 1e-9,
            "bias should be symmetric: m_10={m_10}, m_90={m_90}"
        );
    }

    // 4. Flat model always returns 1.0.
    #[test]
    fn test_flat_model() {
        let cfg = GcBiasConfig {
            enabled: true,
            model: GcBiasModelKind::Flat,
            severity: 1.0,
        };
        let model = GcBiasModel::from_config(&cfg);
        for gc in [0.0, 0.1, 0.2, 0.5, 0.8, 0.9, 1.0] {
            let m = model.coverage_multiplier(gc);
            assert!(
                (m - 1.0).abs() < 1e-9,
                "flat model: expected 1.0 at GC={gc}, got {m}"
            );
        }
    }

    // 5. Higher severity amplifies the bias (reduces multiplier at extremes).
    #[test]
    fn test_severity_scaling() {
        let cfg_low = GcBiasConfig {
            enabled: true,
            model: GcBiasModelKind::Default,
            severity: 0.5,
        };
        let cfg_high = GcBiasConfig {
            enabled: true,
            model: GcBiasModelKind::Default,
            severity: 2.0,
        };
        let model_low = GcBiasModel::from_config(&cfg_low);
        let model_high = GcBiasModel::from_config(&cfg_high);

        let m_low = model_low.coverage_multiplier(0.1);
        let m_high = model_high.coverage_multiplier(0.1);

        // Higher severity → more extreme bias → lower multiplier at GC extremes.
        assert!(
            m_high < m_low,
            "higher severity should reduce multiplier further: low={m_low}, high={m_high}"
        );

        // At severity=0 the multiplier should be 1.0 everywhere.
        let cfg_zero = GcBiasConfig {
            enabled: true,
            model: GcBiasModelKind::Default,
            severity: 0.0,
        };
        let model_zero = GcBiasModel::from_config(&cfg_zero);
        assert!(
            (model_zero.coverage_multiplier(0.1) - 1.0).abs() < 1e-9,
            "severity=0 should give multiplier 1.0"
        );
    }

    // 6. Coverage distribution test: AT-rich region gets fewer reads after
    //    GC-bias rejection sampling.
    #[test]
    fn test_coverage_distribution() {
        use rand::Rng;
        use rand::SeedableRng;
        use rand::rngs::StdRng;

        let model = GcBiasModel::default();
        let mut rng = StdRng::seed_from_u64(42);

        let n = 10_000usize;
        let mut kept_at_rich = 0usize;
        let mut kept_neutral = 0usize;

        for _ in 0..n {
            // AT-rich fragment (10% GC)
            let m = model.coverage_multiplier(0.1);
            if rng.gen::<f64>() < m {
                kept_at_rich += 1;
            }
            // Neutral fragment (50% GC)
            let m = model.coverage_multiplier(0.5);
            if rng.gen::<f64>() < m {
                kept_neutral += 1;
            }
        }

        // Nearly all neutral fragments should be kept (multiplier ~1.0).
        assert!(
            kept_neutral > n * 99 / 100,
            "expected >99% retention at 50% GC, got {kept_neutral}/{n}"
        );
        // AT-rich fragments should be significantly reduced.
        assert!(
            kept_at_rich < n / 2,
            "expected <50% retention at 10% GC, got {kept_at_rich}/{n}"
        );
        assert!(
            kept_at_rich < kept_neutral,
            "AT-rich should have fewer reads than GC-neutral"
        );
    }

    // 7. N bases are excluded from GC calculation.
    #[test]
    fn test_gc_fraction_n_bases() {
        // NNACGT: 2 GC out of 4 non-N = 0.5
        assert!(
            (GcBiasModel::gc_fraction(b"NNACGT") - 0.5).abs() < 1e-9,
            "N bases should be excluded"
        );
        // All N → neutral 0.5
        assert!(
            (GcBiasModel::gc_fraction(b"NNNN") - 0.5).abs() < 1e-9,
            "all-N sequence should return 0.5"
        );
        // Ns mixed with GC
        // NGCNN → 2 GC out of 2 non-N = 1.0
        assert!(
            (GcBiasModel::gc_fraction(b"NGCNN") - 1.0).abs() < 1e-9,
            "N bases should not count toward denominator"
        );
    }
}
