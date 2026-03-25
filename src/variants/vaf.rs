//! Variant allele frequency (VAF) calculation and stochastic alt-read sampling via binomial draws.

use rand::Rng;
use rand_distr::{Binomial, Distribution};

/// Compute expected VAF from tumour model parameters.
///
/// VAF = CCF * multiplicity * purity / (purity * CN_tumour + (1 - purity) * CN_normal)
// Called from tests and from vaf_in_amplified_region in cnv.rs (also test-only).
#[cfg(test)]
#[must_use]
pub fn expected_vaf(
    ccf: f64,
    multiplicity: u32,
    purity: f64,
    cn_tumour: u32,
    cn_normal: u32,
) -> f64 {
    let numerator = ccf * (multiplicity as f64) * purity;
    let denominator = purity * (cn_tumour as f64) + (1.0 - purity) * (cn_normal as f64);
    if denominator == 0.0 {
        return 0.0;
    }
    numerator / denominator
}

/// Sample the number of variant reads from a binomial distribution.
///
/// This avoids the deterministic spike-in bias present in tools like BAMSurgeon,
/// where a 10% VAF at 30x always gives exactly 3 alt reads.
///
/// # Examples
///
/// ```
/// use rand::SeedableRng;
/// use rand::rngs::StdRng;
/// use varforge::variants::vaf::sample_alt_count;
///
/// let mut rng = StdRng::seed_from_u64(42);
///
/// // At VAF=0.0, no alt reads are produced.
/// assert_eq!(sample_alt_count(100, 0.0, &mut rng), 0);
///
/// // At VAF=1.0, all reads are alt.
/// assert_eq!(sample_alt_count(50, 1.0, &mut rng), 50);
///
/// // At an intermediate VAF, the result is stochastic but bounded.
/// let count = sample_alt_count(100, 0.1, &mut rng);
/// assert!(count <= 100);
/// ```
pub fn sample_alt_count<R: Rng>(total_depth: u32, expected_vaf: f64, rng: &mut R) -> u32 {
    if expected_vaf <= 0.0 {
        return 0;
    }
    if expected_vaf >= 1.0 {
        return total_depth;
    }
    let dist = Binomial::new(total_depth as u64, expected_vaf).expect("invalid binomial params");
    dist.sample(rng) as u32
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_expected_vaf_simple() {
        // Clonal, heterozygous, pure tumour: CCF=1, mult=1, purity=1, CN_t=2, CN_n=2
        let vaf = expected_vaf(1.0, 1, 1.0, 2, 2);
        assert!((vaf - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_expected_vaf_with_purity() {
        // CCF=1, mult=1, purity=0.5, CN_t=2, CN_n=2
        // VAF = 1 * 1 * 0.5 / (0.5 * 2 + 0.5 * 2) = 0.5 / 2.0 = 0.25
        let vaf = expected_vaf(1.0, 1, 0.5, 2, 2);
        assert!((vaf - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_expected_vaf_subclone() {
        // CCF=0.3, mult=1, purity=1.0, CN_t=2, CN_n=2
        // VAF = 0.3 * 1 * 1.0 / (1.0 * 2) = 0.15
        let vaf = expected_vaf(0.3, 1, 1.0, 2, 2);
        assert!((vaf - 0.15).abs() < 1e-10);
    }

    #[test]
    fn test_expected_vaf_cnv_gain() {
        // CCF=1, mult=1, purity=1, CN_t=4, CN_n=2
        // VAF = 1 * 1 * 1 / (1 * 4) = 0.25
        let vaf = expected_vaf(1.0, 1, 1.0, 4, 2);
        assert!((vaf - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_expected_vaf_amplified_allele() {
        // CCF=1, mult=3, purity=1, CN_t=4, CN_n=2
        // VAF = 1 * 3 * 1 / (1 * 4) = 0.75
        let vaf = expected_vaf(1.0, 3, 1.0, 4, 2);
        assert!((vaf - 0.75).abs() < 1e-10);
    }

    #[test]
    fn test_expected_vaf_low_tf() {
        // Liquid biopsy: CCF=1, mult=1, purity=0.005 (0.5% TF), CN_t=2, CN_n=2
        // VAF = 1 * 1 * 0.005 / (0.005 * 2 + 0.995 * 2) = 0.005 / 2.0 = 0.0025
        let vaf = expected_vaf(1.0, 1, 0.005, 2, 2);
        assert!((vaf - 0.0025).abs() < 1e-10);
    }

    #[test]
    fn test_binomial_sampling_mean() {
        let mut rng = StdRng::seed_from_u64(42);
        let n = 100_000;
        let depth = 100u32;
        let vaf = 0.1;
        let total: u32 = (0..n).map(|_| sample_alt_count(depth, vaf, &mut rng)).sum();
        let mean = total as f64 / n as f64;
        assert!(
            (mean - 10.0).abs() < 0.5,
            "binomial mean {} should be close to 10",
            mean
        );
    }

    #[test]
    fn test_binomial_sampling_variance() {
        let mut rng = StdRng::seed_from_u64(42);
        let n = 100_000;
        let depth = 30u32;
        let vaf = 0.1;
        let counts: Vec<u32> = (0..n)
            .map(|_| sample_alt_count(depth, vaf, &mut rng))
            .collect();
        let mean = counts.iter().sum::<u32>() as f64 / n as f64;
        let variance = counts
            .iter()
            .map(|&c| (c as f64 - mean).powi(2))
            .sum::<f64>()
            / n as f64;
        // Expected variance = n*p*(1-p) = 30 * 0.1 * 0.9 = 2.7
        assert!(
            (variance - 2.7).abs() < 0.5,
            "variance {} should be close to 2.7",
            variance
        );
    }

    #[test]
    fn test_binomial_not_deterministic() {
        // The whole point: same VAF/depth should give DIFFERENT counts
        let mut rng = StdRng::seed_from_u64(42);
        let counts: Vec<u32> = (0..100)
            .map(|_| sample_alt_count(30, 0.1, &mut rng))
            .collect();
        let unique: std::collections::HashSet<_> = counts.iter().collect();
        assert!(
            unique.len() > 1,
            "binomial sampling should produce variable counts"
        );
    }

    #[test]
    fn test_ultra_low_vaf_unbiased_2000x() {
        // At depth=2000, VAF=0.0001, the expected mean alt count is 0.2 per read.
        // Verify the binomial sampler is unbiased at ultra-low VAF over 100,000 trials.
        let mut rng = StdRng::seed_from_u64(42);
        let n = 100_000u64;
        let total: u64 = (0..n)
            .map(|_| sample_alt_count(2000, 0.0001, &mut rng) as u64)
            .sum();
        let mean = total as f64 / n as f64;
        assert!(
            (0.15..=0.25).contains(&mean),
            "mean {} should be 0.2 ± 0.05 for ultra-low VAF at 2000x",
            mean
        );
    }

    #[test]
    fn test_ultra_low_vaf_unbiased_5000x() {
        // At depth=5000, VAF=0.00001, the expected mean alt count is 0.05 per read.
        // Verify the binomial sampler is unbiased at extremely low VAF over 100,000 trials.
        let mut rng = StdRng::seed_from_u64(43);
        let n = 100_000u64;
        let total: u64 = (0..n)
            .map(|_| sample_alt_count(5000, 0.00001, &mut rng) as u64)
            .sum();
        let mean = total as f64 / n as f64;
        assert!(
            (0.03..=0.07).contains(&mean),
            "mean {} should be 0.05 ± 0.02 for ultra-low VAF at 5000x",
            mean
        );
    }

    #[test]
    fn test_vaf_1_all_alt() {
        // At VAF=1.0 every read must be an alt. The sampler should short-circuit and
        // return total_depth without invoking the binomial distribution.
        let mut rng = StdRng::seed_from_u64(44);
        for _ in 0..100 {
            let count = sample_alt_count(100, 1.0, &mut rng);
            assert_eq!(count, 100, "VAF=1.0 must always yield all alt reads");
        }
    }
}
