//! Capture efficiency model for targeted sequencing (panel / WES).
//!
//! Models three phenomena common to hybridisation-capture enrichment:
//! 1. **Per-target coverage variation** – each target gets a coverage
//!    multiplier sampled from LogNormal(0, uniformity).
//! 2. **Edge drop-off** – exponential decay of coverage within
//!    `edge_dropoff_bases` of each target boundary.
//! 3. **Off-target fraction** – some fraction of reads lands outside all
//!    target regions and receives `off_target_fraction × mean_coverage`.

use rand::Rng;
use rand_distr::{Distribution, LogNormal};

use super::types::Region;

/// Model for target-capture efficiency applied during coverage calculation.
pub struct CaptureModel {
    /// The target regions for this panel / WES experiment.
    pub target_regions: Vec<Region>,
    /// Optional fixed depth multiplier per target, parallel to `target_regions`.
    ///
    /// When `Some(d)`, the target uses `d` directly instead of sampling from
    /// LogNormal.  When `None`, the normal LogNormal sampling path is used.
    pub target_depths: Vec<Option<f64>>,
    /// Fraction of reads that map off-target (default: 0.2).
    pub off_target_fraction: f64,
    /// Controls per-target coverage variation via LogNormal σ.
    /// 0.0 = perfectly uniform, larger values = more variable.
    pub coverage_uniformity: f64,
    /// Number of bases at each target edge over which coverage decays (default: 50).
    pub edge_dropoff_bases: u32,
    /// Sequencing mode: `"panel"` or `"amplicon"`.
    pub mode: String,
}

impl CaptureModel {
    /// Create a new CaptureModel.
    ///
    /// `target_depths` is parallel to `target_regions`.  Pass an empty `Vec`
    /// or a `Vec` of `None`s to fall back to LogNormal sampling for all targets.
    pub fn new(
        target_regions: Vec<Region>,
        target_depths: Vec<Option<f64>>,
        off_target_fraction: f64,
        coverage_uniformity: f64,
        edge_dropoff_bases: u32,
        mode: String,
    ) -> Self {
        Self {
            target_regions,
            target_depths,
            off_target_fraction,
            coverage_uniformity,
            edge_dropoff_bases,
            mode,
        }
    }

    /// Returns true when this model is in amplicon mode.
    ///
    /// In amplicon mode, fragments exactly span each target region rather than
    /// being sampled from a random position within the broader capture window.
    pub fn is_amplicon(&self) -> bool {
        self.mode == "amplicon"
    }

    /// Sample a per-target coverage multiplier from LogNormal(0, uniformity).
    ///
    /// When `uniformity == 0.0` the distribution degenerates to a point mass
    /// at 1.0, so every target gets an identical multiplier.
    pub fn sample_target_multiplier<R: Rng>(&self, rng: &mut R) -> f64 {
        if self.coverage_uniformity == 0.0 {
            return 1.0;
        }
        // LogNormal(μ=0, σ=uniformity) has median 1.0.
        let dist =
            LogNormal::new(0.0, self.coverage_uniformity).expect("invalid LogNormal parameters");
        dist.sample(rng)
    }

    /// Sample a per-target coverage multiplier for the target at `index`.
    ///
    /// If a fixed depth was recorded for this target in `target_depths`, that
    /// value is returned directly.  Otherwise the LogNormal sampling path is
    /// used.
    fn sample_target_multiplier_at<R: Rng>(&self, index: usize, rng: &mut R) -> f64 {
        if let Some(Some(d)) = self.target_depths.get(index) {
            return *d;
        }
        self.sample_target_multiplier(rng)
    }

    /// Compute the edge drop-off multiplier at `position` within a target.
    ///
    /// `target_start` and `target_end` are 0-based half-open coordinates.
    /// Returns a value in (0, 1]: 1.0 in the interior, decaying exponentially
    /// within `edge_dropoff_bases` of either boundary.
    ///
    /// The decay formula is: `e^(-distance_from_edge / decay_constant)` where
    /// `decay_constant = edge_dropoff_bases / 3.0` (so that coverage reaches
    /// ~5 % at the very edge).
    pub fn edge_multiplier(&self, position: u64, target_start: u64, target_end: u64) -> f64 {
        if self.edge_dropoff_bases == 0 {
            return 1.0;
        }
        let dropoff = self.edge_dropoff_bases as u64;
        let decay_const = (self.edge_dropoff_bases as f64) / 3.0;

        // Distance from the nearest boundary.
        let dist_from_start = position.saturating_sub(target_start);
        let dist_from_end = target_end.saturating_sub(position + 1);
        let dist_from_edge = dist_from_start.min(dist_from_end);

        if dist_from_edge >= dropoff {
            1.0
        } else {
            // Saturating exponential: 0 at the edge, approaches 1 in the interior.
            // multiplier = 1 - exp(-dist_from_edge / decay_const)
            let raw = 1.0 - (-(dist_from_edge as f64) / decay_const).exp();
            // Clamp to (0, 1] — the above is always in [0, 1), so add epsilon floor.
            raw.max(f64::EPSILON)
        }
    }

    /// Effective coverage for a position that is **on-target**.
    ///
    /// `mean_coverage` is the experiment-wide mean coverage.
    /// `target_multiplier` is the per-target LogNormal sample for this target.
    // Called only in tests; production code calls coverage_multiplier directly.
    #[cfg(test)]
    pub fn on_target_coverage(
        &self,
        mean_coverage: f64,
        target_multiplier: f64,
        edge_mult: f64,
    ) -> f64 {
        mean_coverage * target_multiplier * edge_mult
    }

    /// Effective coverage for a position that is **off-target**.
    // Called only in tests; production code calls coverage_multiplier directly.
    #[cfg(test)]
    pub fn off_target_coverage(&self, mean_coverage: f64) -> f64 {
        mean_coverage * self.off_target_fraction
    }

    /// Return the effective coverage multiplier for `position` on `chrom`.
    ///
    /// Looks up which (if any) target region contains the position, applies
    /// the per-target multiplier and edge drop-off.  If the position is not
    /// covered by any target, returns the off-target fraction.
    ///
    /// `target_multipliers` must be parallel to `self.target_regions`.
    pub fn coverage_multiplier_at(
        &self,
        chrom: &str,
        position: u64,
        target_multipliers: &[f64],
    ) -> f64 {
        for (i, target) in self.target_regions.iter().enumerate() {
            if target.chrom == chrom && position >= target.start && position < target.end {
                let edge_mult = self.edge_multiplier(position, target.start, target.end);
                return target_multipliers[i] * edge_mult;
            }
        }
        // Off-target position: return the off-target fraction directly as the
        // multiplier (caller multiplies this by mean_coverage).
        self.off_target_fraction
    }

    /// Sample one multiplier per target region, returned in the same order as
    /// `self.target_regions`.
    ///
    /// Targets that have a fixed depth in `target_depths` use that value
    /// directly.  All others are sampled from LogNormal.
    pub fn sample_all_target_multipliers<R: Rng>(&self, rng: &mut R) -> Vec<f64> {
        (0..self.target_regions.len())
            .map(|i| self.sample_target_multiplier_at(i, rng))
            .collect()
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

    fn make_regions() -> Vec<Region> {
        vec![
            Region::new("chr1", 0, 500),
            Region::new("chr1", 1000, 1500),
            Region::new("chr2", 0, 300),
        ]
    }

    // -----------------------------------------------------------------------
    // 1. test_uniform_targets
    //    uniformity=0 must give multiplier 1.0 for every target.
    // -----------------------------------------------------------------------
    #[test]
    fn test_uniform_targets() {
        let model = CaptureModel::new(make_regions(), vec![], 0.2, 0.0, 0, "panel".to_string());
        let mut rng = StdRng::seed_from_u64(1);

        // Draw many multipliers – they must all be exactly 1.0.
        for _ in 0..200 {
            let m = model.sample_target_multiplier(&mut rng);
            assert_eq!(m, 1.0, "uniformity=0 should give multiplier 1.0, got {m}");
        }
    }

    // -----------------------------------------------------------------------
    // 2. test_variable_targets
    //    uniformity>0 must produce a spread of multipliers across targets.
    // -----------------------------------------------------------------------
    #[test]
    fn test_variable_targets() {
        let model = CaptureModel::new(make_regions(), vec![], 0.2, 0.5, 0, "panel".to_string());
        let mut rng = StdRng::seed_from_u64(2);

        let multipliers: Vec<f64> = (0..500)
            .map(|_| model.sample_target_multiplier(&mut rng))
            .collect();

        let min = multipliers.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = multipliers
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);

        // With σ=0.5, the distribution has a coefficient of variation > 50%.
        // The range across 500 samples should be substantial (at least 2×).
        assert!(
            max / min > 2.0,
            "expected variable coverage (max/min > 2.0), got min={min:.3} max={max:.3}"
        );

        // The geometric mean should be close to 1.0 (since μ=0 in log-space).
        let log_mean: f64 =
            multipliers.iter().map(|x| x.ln()).sum::<f64>() / multipliers.len() as f64;
        assert!(
            log_mean.abs() < 0.15,
            "geometric mean should be near 1.0 (log mean near 0), got {log_mean:.4}"
        );
    }

    // -----------------------------------------------------------------------
    // 3. test_edge_dropoff
    //    Coverage decreases at target boundaries.
    // -----------------------------------------------------------------------
    #[test]
    fn test_edge_dropoff() {
        let dropoff_bases = 50u32;
        let model = CaptureModel::new(
            make_regions(),
            vec![],
            0.2,
            0.0,
            dropoff_bases,
            "panel".to_string(),
        );

        let target_start = 0u64;
        let target_end = 500u64;

        // Interior position: should be 1.0.
        let interior = model.edge_multiplier(250, target_start, target_end);
        assert_eq!(interior, 1.0, "interior multiplier should be 1.0");

        // Position right at the start boundary: distance from edge = 0.
        let at_edge = model.edge_multiplier(0, target_start, target_end);
        assert!(
            at_edge < 0.5,
            "multiplier at edge should be < 0.5, got {at_edge:.4}"
        );

        // Coverage should increase monotonically as we move away from the edge.
        let m0 = model.edge_multiplier(0, target_start, target_end);
        let m10 = model.edge_multiplier(10, target_start, target_end);
        let m30 = model.edge_multiplier(30, target_start, target_end);
        let m60 = model.edge_multiplier(60, target_start, target_end);

        assert!(
            m0 < m10,
            "coverage should increase moving away from edge: m0={m0:.4} m10={m10:.4}"
        );
        assert!(
            m10 < m30,
            "coverage should increase moving away from edge: m10={m10:.4} m30={m30:.4}"
        );
        assert!(
            m30 < m60,
            "coverage should continue increasing past edge zone: m30={m30:.4} m60={m60:.4}"
        );
        assert_eq!(
            m60, 1.0,
            "position beyond edge_dropoff_bases should have multiplier 1.0"
        );
    }

    // -----------------------------------------------------------------------
    // 4. test_off_target_reads
    //    Off-target fraction approximately correct.
    // -----------------------------------------------------------------------
    #[test]
    fn test_off_target_reads() {
        let off_target_fraction = 0.2_f64;
        let model = CaptureModel::new(
            make_regions(),
            vec![],
            off_target_fraction,
            0.0,
            0,
            "panel".to_string(),
        );
        let mean_coverage = 100.0_f64;

        // Simulate a mix of on-target and off-target positions.
        // We probe positions on chr3 (not in targets) – all off-target.
        let mut rng = StdRng::seed_from_u64(4);
        let multipliers = model.sample_all_target_multipliers(&mut rng);

        let off_target_cov = model.off_target_coverage(mean_coverage);
        let on_target_cov = model.on_target_coverage(mean_coverage, 1.0, 1.0);

        // Off-target coverage must equal mean * off_target_fraction.
        let expected_off = mean_coverage * off_target_fraction;
        assert!(
            (off_target_cov - expected_off).abs() < 1e-9,
            "off-target coverage should be {expected_off:.2}, got {off_target_cov:.2}"
        );

        // On-target coverage (with uniform multiplier=1, no edge) should be mean.
        assert!(
            (on_target_cov - mean_coverage).abs() < 1e-9,
            "on-target coverage should equal mean {mean_coverage:.2}, got {on_target_cov:.2}"
        );

        // coverage_multiplier_at for an off-target position (chr3) should return off_target_fraction.
        let off_mult = model.coverage_multiplier_at("chr3", 0, &multipliers);
        assert!(
            (off_mult - off_target_fraction).abs() < 1e-9,
            "off-target multiplier should be {off_target_fraction:.3}, got {off_mult:.6}"
        );
    }

    // -----------------------------------------------------------------------
    // 5. test_no_capture
    //    With uniformity=0, no edge dropoff, and off_target_fraction=1.0,
    //    every position (on- or off-target) gets the same multiplier (1.0),
    //    simulating uniform coverage – equivalent to capture being disabled.
    // -----------------------------------------------------------------------
    #[test]
    fn test_no_capture() {
        // off_target_fraction=1.0 means off-target positions get full coverage.
        // uniformity=0 means all targets get multiplier 1.0.
        // edge_dropoff_bases=0 means no edge effect.
        let regions = make_regions();
        let n = regions.len();
        let model = CaptureModel::new(regions, vec![], 1.0, 0.0, 0, "panel".to_string());
        let mut rng = StdRng::seed_from_u64(5);

        let multipliers = model.sample_all_target_multipliers(&mut rng);
        assert_eq!(multipliers.len(), n);

        // All on-target multipliers should be 1.0.
        for &m in &multipliers {
            assert_eq!(m, 1.0, "no-capture mode: expected multiplier 1.0, got {m}");
        }

        // Off-target positions should also get multiplier 1.0.
        let off = model.coverage_multiplier_at("chr99", 0, &multipliers);
        assert_eq!(
            off, 1.0,
            "no-capture off-target should also be 1.0, got {off}"
        );

        // All on-target positions (interior) should be 1.0.
        let on = model.coverage_multiplier_at("chr1", 250, &multipliers);
        assert_eq!(on, 1.0, "no-capture on-target should be 1.0, got {on}");
    }
}
