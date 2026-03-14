use rand::Rng;
use rand_distr::{Distribution, Normal};

/// Generates base quality scores for a read.
pub trait QualityModel: Send + Sync {
    /// Generate quality scores for a read of the given length.
    fn generate_qualities<R: Rng>(&self, read_length: usize, rng: &mut R) -> Vec<u8>;

    /// Given a quality score, return the probability of an error.
    fn error_probability(quality: u8) -> f64 {
        10.0_f64.powf(-(quality as f64) / 10.0)
    }
}

/// Parametric quality model with position-dependent decay.
///
/// Quality starts high and decays toward the 3' end of the read,
/// mimicking Illumina's characteristic quality profile.
pub struct ParametricQualityModel {
    mean_quality: f64,
    tail_decay: f64,
    noise_sd: f64,
}

impl ParametricQualityModel {
    pub fn new(mean_quality: u8, tail_decay: f64) -> Self {
        Self {
            mean_quality: mean_quality as f64,
            tail_decay,
            noise_sd: 3.0,
        }
    }
}

impl QualityModel for ParametricQualityModel {
    fn generate_qualities<R: Rng>(&self, read_length: usize, rng: &mut R) -> Vec<u8> {
        let noise = Normal::new(0.0, self.noise_sd).unwrap();
        (0..read_length)
            .map(|pos| {
                let decay = self.tail_decay * (pos as f64) * (pos as f64);
                let q = self.mean_quality - decay + noise.sample(rng);
                q.round().clamp(2.0, 41.0) as u8
            })
            .collect()
    }
}

/// Inject base-call errors consistent with quality scores.
///
/// For each base, rolls a random number against the Phred error probability.
/// If an error occurs, substitutes with a random different base.
pub fn inject_errors<R: Rng>(seq: &mut [u8], qual: &[u8], rng: &mut R) {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

    for i in 0..seq.len() {
        let error_prob = ParametricQualityModel::error_probability(qual[i]);
        if rng.gen::<f64>() < error_prob {
            let original = seq[i];
            loop {
                let new_base = BASES[rng.gen_range(0..4)];
                if new_base != original {
                    seq[i] = new_base;
                    break;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_parametric_quality_length() {
        let model = ParametricQualityModel::new(36, 0.003);
        let mut rng = StdRng::seed_from_u64(42);
        let quals = model.generate_qualities(150, &mut rng);
        assert_eq!(quals.len(), 150);
    }

    #[test]
    fn test_quality_decay_toward_end() {
        let model = ParametricQualityModel::new(36, 0.003);
        let mut rng = StdRng::seed_from_u64(42);
        // Average over many runs to smooth out noise
        let n = 1000;
        let mut start_sum = 0u64;
        let mut end_sum = 0u64;
        for _ in 0..n {
            let quals = model.generate_qualities(150, &mut rng);
            start_sum += quals[..10].iter().map(|&q| q as u64).sum::<u64>();
            end_sum += quals[140..].iter().map(|&q| q as u64).sum::<u64>();
        }
        let start_mean = start_sum as f64 / (n * 10) as f64;
        let end_mean = end_sum as f64 / (n * 10) as f64;
        assert!(
            start_mean > end_mean,
            "start quality ({}) should be higher than end ({})",
            start_mean,
            end_mean
        );
    }

    #[test]
    fn test_quality_bounds() {
        let model = ParametricQualityModel::new(36, 0.003);
        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..100 {
            let quals = model.generate_qualities(150, &mut rng);
            for &q in &quals {
                assert!(q >= 2 && q <= 41, "quality {} out of bounds", q);
            }
        }
    }

    #[test]
    fn test_error_probability() {
        let p30 = ParametricQualityModel::error_probability(30);
        assert!((p30 - 0.001).abs() < 1e-6, "Q30 should be 0.001 error rate");

        let p20 = ParametricQualityModel::error_probability(20);
        assert!((p20 - 0.01).abs() < 1e-6, "Q20 should be 0.01 error rate");

        let p10 = ParametricQualityModel::error_probability(10);
        assert!((p10 - 0.1).abs() < 1e-6, "Q10 should be 0.1 error rate");
    }

    #[test]
    fn test_inject_errors_rate() {
        let mut rng = StdRng::seed_from_u64(42);
        let n = 100_000;
        let mut total_errors = 0usize;

        // Use Q20 (1% error rate) for all positions
        let qual = vec![20u8; 100];
        for _ in 0..n {
            let original = vec![b'A'; 100];
            let mut seq = original.clone();
            inject_errors(&mut seq, &qual, &mut rng);
            total_errors += seq.iter().zip(original.iter()).filter(|(a, b)| a != b).count();
        }

        let observed_rate = total_errors as f64 / (n * 100) as f64;
        assert!(
            (observed_rate - 0.01).abs() < 0.002,
            "error rate {} should be close to 0.01",
            observed_rate
        );
    }

    #[test]
    fn test_inject_errors_different_base() {
        let mut rng = StdRng::seed_from_u64(42);
        let qual = vec![10u8; 1000]; // high error rate
        let mut seq = vec![b'A'; 1000];
        inject_errors(&mut seq, &qual, &mut rng);

        for &base in &seq {
            assert!(
                matches!(base, b'A' | b'C' | b'G' | b'T'),
                "invalid base {}",
                base as char
            );
        }
    }
}
