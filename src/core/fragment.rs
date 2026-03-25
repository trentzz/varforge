//! Fragment size sampling: Gaussian WGS, cfDNA nucleosomal, log-normal long-read, and PCR family size models.

use anyhow::{Context, Result};
use rand::Rng;
use rand_distr::{Distribution, LogNormal, Normal};

use crate::io::config::LongReadFragmentConfig;

/// Samples fragment sizes according to a configured model.
pub trait FragmentSampler: Send + Sync {
    fn sample<R: Rng>(&self, rng: &mut R) -> usize;
}

/// Standard WGS fragment size: Gaussian distribution.
pub struct NormalFragmentSampler {
    dist: Normal<f64>,
    min_size: usize,
}

impl NormalFragmentSampler {
    /// Create a Gaussian fragment size sampler.
    ///
    /// Returns an error if `sd` is NaN or infinite. The minimum sampled fragment
    /// size is clamped to 50 bp.
    ///
    /// # Examples
    ///
    /// ```
    /// use varforge::core::fragment::{NormalFragmentSampler, FragmentSampler};
    /// use rand::SeedableRng;
    /// use rand::rngs::StdRng;
    ///
    /// let sampler = NormalFragmentSampler::new(300.0, 50.0).unwrap();
    /// let mut rng = StdRng::seed_from_u64(0);
    /// let size = sampler.sample(&mut rng);
    /// // Fragment size is always at least the minimum (50 bp).
    /// assert!(size >= 50);
    ///
    /// // A non-finite sd produces an error.
    /// assert!(NormalFragmentSampler::new(300.0, f64::INFINITY).is_err());
    /// ```
    pub fn new(mean: f64, sd: f64) -> Result<Self> {
        Ok(Self {
            dist: Normal::new(mean, sd)
                .context("invalid normal distribution parameters (mean or sd is NaN/Inf)")?,
            min_size: 50,
        })
    }
}

impl FragmentSampler for NormalFragmentSampler {
    fn sample<R: Rng>(&self, rng: &mut R) -> usize {
        let size = self.dist.sample(rng).round() as i64;
        size.max(self.min_size as i64) as usize
    }
}

/// cfDNA fragment size: mixture of nucleosomal peaks with 10bp periodicity.
pub struct CfdnaFragmentSampler {
    mono_dist: Normal<f64>,
    di_dist: Normal<f64>,
    mono_weight: f64,
    ctdna_dist: Normal<f64>,
    ctdna_fraction: f64,
    min_size: usize,
}

impl CfdnaFragmentSampler {
    /// Create a cfDNA sampler.
    ///
    /// - `mono_peak`: mononucleosomal peak (default ~167)
    /// - `di_peak`: dinucleosomal peak (default ~334)
    /// - `mono_weight`: weight of mono vs di peak (default 0.85)
    /// - `ctdna_fraction`: fraction of fragments that are tumour-derived (shorter, ~143 bp)
    /// - `mono_sd`: SD of the mononucleosomal peak in bp (default 20.0)
    /// - `di_sd`: SD of the dinucleosomal peak in bp (default 30.0)
    ///
    /// Default SDs are from Cristiano et al. 2019 Science (DELFI study), which
    /// characterised cfDNA fragment size distributions across cancer types.
    pub fn new(
        mono_peak: f64,
        di_peak: f64,
        mono_weight: f64,
        ctdna_fraction: f64,
        mono_sd: f64,
        di_sd: f64,
    ) -> Result<Self> {
        Ok(Self {
            mono_dist: Normal::new(mono_peak, mono_sd)
                .context("invalid mononucleosomal distribution parameters")?,
            di_dist: Normal::new(di_peak, di_sd)
                .context("invalid dinucleosomal distribution parameters")?,
            mono_weight,
            ctdna_dist: Normal::new(143.0, 15.0)
                .context("invalid ctDNA distribution parameters")?,
            ctdna_fraction,
            min_size: 50,
        })
    }

    /// Add 10bp periodicity sub-peaks to a sampled size.
    fn apply_periodicity<R: Rng>(&self, size: f64, rng: &mut R) -> f64 {
        // DNA wraps around nucleosomes with ~10bp pitch
        // Add a small sinusoidal modulation
        let period_phase = (size / 10.0).fract();
        let adjustment = (period_phase * 2.0 * std::f64::consts::PI).sin() * 2.0;
        let noise: f64 = rng.random_range(-1.0..1.0);
        size + adjustment + noise
    }
}

impl FragmentSampler for CfdnaFragmentSampler {
    fn sample<R: Rng>(&self, rng: &mut R) -> usize {
        let is_ctdna = rng.random::<f64>() < self.ctdna_fraction;

        let raw_size = if is_ctdna {
            self.ctdna_dist.sample(rng)
        } else if rng.random::<f64>() < self.mono_weight {
            self.mono_dist.sample(rng)
        } else {
            self.di_dist.sample(rng)
        };

        let size = self.apply_periodicity(raw_size, rng);
        (size.round() as i64).max(self.min_size as i64) as usize
    }
}

/// Sample a fragment length for a long-read platform using a log-normal distribution.
///
/// The distribution is parameterised by the linear-space mean and standard
/// deviation. These are converted to the underlying normal parameters (mu, sigma)
/// before sampling. The result is clamped to [min_len, max_len].
pub fn sample_long_read_length<R: Rng>(cfg: &LongReadFragmentConfig, rng: &mut R) -> Result<usize> {
    let mean = cfg.mean as f64;
    let sd = cfg.sd as f64;
    let variance = sd * sd;
    // Convert linear-space mean/sd to log-normal mu/sigma.
    let mu = (mean * mean / (mean * mean + variance).sqrt()).ln();
    let sigma = (1.0_f64 + variance / (mean * mean)).ln().sqrt();
    let dist = LogNormal::new(mu, sigma)
        .context("invalid log-normal parameters for long-read fragment size")?;
    let sample = dist.sample(rng) as usize;
    Ok(sample.clamp(cfg.min_len, cfg.max_len))
}

/// PCR family size sampler using a log-normal distribution.
pub struct PcrFamilySizeSampler {
    dist: LogNormal<f64>,
}

impl PcrFamilySizeSampler {
    pub fn new(mean: f64, sd: f64) -> Result<Self> {
        // Convert mean/sd of the actual family size to log-space parameters
        let variance = sd * sd;
        let mu = (mean * mean / (mean * mean + variance).sqrt()).ln();
        let sigma = (1.0 + variance / (mean * mean)).ln().sqrt();
        Ok(Self {
            dist: LogNormal::new(mu, sigma)
                .context("invalid log-normal parameters for PCR family size")?,
        })
    }

    pub fn sample<R: Rng>(&self, rng: &mut R) -> usize {
        let size = self.dist.sample(rng).round() as usize;
        size.max(1) // at least 1 copy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_normal_fragment_sampler() {
        let sampler = NormalFragmentSampler::new(300.0, 50.0).unwrap();
        let mut rng = StdRng::seed_from_u64(42);
        let sizes: Vec<usize> = (0..1000).map(|_| sampler.sample(&mut rng)).collect();

        let mean = sizes.iter().sum::<usize>() as f64 / sizes.len() as f64;
        assert!(
            (mean - 300.0).abs() < 10.0,
            "mean {} too far from 300",
            mean
        );
        assert!(sizes.iter().all(|&s| s >= 50), "no fragment below minimum");
    }

    #[test]
    fn test_cfdna_fragment_sampler() {
        let sampler = CfdnaFragmentSampler::new(167.0, 334.0, 0.85, 0.0, 20.0, 30.0).unwrap();
        let mut rng = StdRng::seed_from_u64(42);
        let sizes: Vec<usize> = (0..10000).map(|_| sampler.sample(&mut rng)).collect();

        let mean = sizes.iter().sum::<usize>() as f64 / sizes.len() as f64;
        // With 85% mono (167) and 15% di (334), expected mean ~192
        assert!(
            mean > 150.0 && mean < 230.0,
            "cfDNA mean {} unexpected",
            mean
        );
    }

    #[test]
    fn test_cfdna_ctdna_shorter() {
        let sampler_normal =
            CfdnaFragmentSampler::new(167.0, 334.0, 0.85, 0.0, 20.0, 30.0).unwrap();
        let sampler_ctdna = CfdnaFragmentSampler::new(167.0, 334.0, 0.85, 1.0, 20.0, 30.0).unwrap();
        let mut rng = StdRng::seed_from_u64(42);

        let normal_mean: f64 = (0..5000)
            .map(|_| sampler_normal.sample(&mut rng) as f64)
            .sum::<f64>()
            / 5000.0;
        let ctdna_mean: f64 = (0..5000)
            .map(|_| sampler_ctdna.sample(&mut rng) as f64)
            .sum::<f64>()
            / 5000.0;

        assert!(
            ctdna_mean < normal_mean,
            "ctDNA mean ({}) should be shorter than normal cfDNA ({})",
            ctdna_mean,
            normal_mean
        );
    }

    #[test]
    fn test_pcr_family_size_sampler() {
        let sampler = PcrFamilySizeSampler::new(3.0, 1.5).unwrap();
        let mut rng = StdRng::seed_from_u64(42);
        let sizes: Vec<usize> = (0..10000).map(|_| sampler.sample(&mut rng)).collect();

        assert!(sizes.iter().all(|&s| s >= 1), "min family size is 1");
        let mean = sizes.iter().sum::<usize>() as f64 / sizes.len() as f64;
        assert!((mean - 3.0).abs() < 1.0, "mean {} too far from 3.0", mean);
    }

    #[test]
    fn test_long_read_length_range() {
        use crate::io::config::LongReadFragmentConfig;
        let cfg = LongReadFragmentConfig {
            mean: 15000,
            sd: 5000,
            min_len: 1000,
            max_len: 100000,
        };
        let mut rng = StdRng::seed_from_u64(42);
        let sizes: Vec<usize> = (0..1000)
            .map(|_| sample_long_read_length(&cfg, &mut rng).unwrap())
            .collect();
        assert!(
            sizes.iter().all(|&s| (1000..=100000).contains(&s)),
            "all samples must be within [min_len, max_len]"
        );
        let mean = sizes.iter().sum::<usize>() as f64 / sizes.len() as f64;
        // Log-normal with these parameters centres around 15000; allow wide tolerance.
        assert!(
            mean > 5000.0 && mean < 30000.0,
            "mean {} outside expected range",
            mean
        );
    }

    /// Regression test for T076: ctdna_fraction=0.01 should yield mostly normal-length
    /// fragments, not short ctDNA fragments. The fraction of samples near the ctDNA
    /// peak (~143 bp) should be approximately 1%, not ~99%.
    #[test]
    fn test_cfdna_low_ctdna_fraction_mostly_normal_lengths() {
        // 1% ctDNA: almost all fragments should be at nucleosomal peaks (>150 bp)
        let sampler = CfdnaFragmentSampler::new(167.0, 334.0, 0.85, 0.01, 20.0, 30.0).unwrap();
        let mut rng = StdRng::seed_from_u64(99);
        let n = 10_000usize;
        // Count fragments that are clearly in the short ctDNA range (<150 bp, well below
        // the mononucleosomal peak of 167 bp).
        // Use 130 bp as the threshold: this is well below the mononucleosomal
        // peak (167 bp ± 20 sd), so the mononucleosomal distribution contributes
        // only ~3% false positives, keeping the total well under the 10% limit.
        let short_count = (0..n).filter(|_| sampler.sample(&mut rng) < 130).count();
        let short_fraction = short_count as f64 / n as f64;
        // With ctdna_fraction=0.01, short-fragment fraction should be near 0.01, not ~0.99.
        assert!(
            short_fraction < 0.10,
            "expected <10% fragments below 130 bp at ctdna_fraction=0.01, got {:.1}%",
            short_fraction * 100.0
        );
    }

    /// Explicit ctdna_fraction=1.0 should yield entirely short fragments.
    #[test]
    fn test_cfdna_full_ctdna_fraction_all_short() {
        let sampler = CfdnaFragmentSampler::new(167.0, 334.0, 0.85, 1.0, 20.0, 30.0).unwrap();
        let mut rng = StdRng::seed_from_u64(77);
        let mean: f64 = (0..5_000)
            .map(|_| sampler.sample(&mut rng) as f64)
            .sum::<f64>()
            / 5_000.0;
        // ctDNA peak is ~143 bp; mean should be well below the nucleosomal peak of 167 bp.
        assert!(
            mean < 160.0,
            "expected mean <160 bp at ctdna_fraction=1.0, got {:.1}",
            mean
        );
    }

    #[test]
    fn test_deterministic_with_seed() {
        let sampler = NormalFragmentSampler::new(300.0, 50.0).unwrap();
        let mut rng1 = StdRng::seed_from_u64(123);
        let mut rng2 = StdRng::seed_from_u64(123);

        let sizes1: Vec<usize> = (0..100).map(|_| sampler.sample(&mut rng1)).collect();
        let sizes2: Vec<usize> = (0..100).map(|_| sampler.sample(&mut rng2)).collect();
        assert_eq!(sizes1, sizes2, "same seed must produce same results");
    }
}
