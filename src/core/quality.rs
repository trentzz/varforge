//! Quality score models and error injection for Illumina, PacBio HiFi, and Nanopore R10 reads.

use rand::Rng;
use rand_distr::{Distribution, Normal};

/// Generates base quality scores and injects sequencing errors into reads.
pub trait QualityModel: Send + Sync {
    /// Generate quality scores for a read of the given length.
    fn generate_qualities<R: Rng>(&self, read_length: usize, rng: &mut R) -> Vec<u8>;

    /// Inject base-call errors into `sequence` consistent with `qualities`.
    ///
    /// Implementations should mutate bases according to their model's error
    /// distribution. The default implementation performs a uniform random
    /// substitution weighted only by the Phred error probability.
    fn inject_errors<R: Rng>(&self, sequence: &mut [u8], qualities: &[u8], rng: &mut R) {
        inject_errors(sequence, qualities, rng);
    }

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

/// Sample per-base qualities for a PacBio HiFi read with a position-dependent curve.
///
/// Real HiFi reads show highest quality in the interior (~Q35) and somewhat lower
/// quality near both ends (~Q25) due to subread stitching at the start and end.
/// Quality transitions linearly over ~50 bp at each end, then holds flat in the middle.
///
/// Reference: PacBio SMRT sequencing quality specification (Q20 minimum, median Q35).
pub fn sample_pacbio_hifi_qualities<R: Rng>(length: usize, rng: &mut R) -> Vec<u8> {
    const PEAK_Q: f64 = 35.0;
    const END_Q: f64 = 25.0;
    const TRANSITION_BP: usize = 50;
    const NOISE_SD: f64 = 1.5;

    let noise = Normal::new(0.0, NOISE_SD).unwrap();

    (0..length)
        .map(|pos| {
            // Distance from the nearer end.
            let dist_from_end = pos.min(length.saturating_sub(1 + pos));
            // Linear ramp from END_Q to PEAK_Q over TRANSITION_BP bases.
            let base_q = if dist_from_end >= TRANSITION_BP {
                PEAK_Q
            } else {
                END_Q + (PEAK_Q - END_Q) * (dist_from_end as f64 / TRANSITION_BP as f64)
            };
            let q = base_q + noise.sample(rng);
            q.round().clamp(2.0, 41.0) as u8
        })
        .collect()
}

/// Returns true if position `pos` falls inside a homopolymer run of length >= `min_len`.
///
/// A homopolymer run is a maximal stretch of the same base. The position is considered
/// inside such a run if the run that covers `pos` is at least `min_len` bases long.
pub fn is_homopolymer_run(seq: &[u8], pos: usize, min_len: usize) -> bool {
    if pos >= seq.len() {
        return false;
    }
    let base = seq[pos];
    // Find the start of the run.
    let mut start = pos;
    while start > 0 && seq[start - 1] == base {
        start -= 1;
    }
    // Find the end of the run.
    let mut end = pos;
    while end + 1 < seq.len() && seq[end + 1] == base {
        end += 1;
    }
    (end - start + 1) >= min_len
}

/// Sample per-base qualities for a Nanopore R10 read.
///
/// Produces a profile centred around Q20 with higher variance than HiFi.
/// Nanopore R10 achieves approximately Q20–Q25 accuracy.
///
/// Positions inside homopolymer runs of length >= 4 receive a penalty of 5–10
/// Phred units, floored at Q10. This reflects the well-known Nanopore systematic
/// error mode in homopolymeric contexts.
pub fn sample_nanopore_r10_qualities<R: Rng>(length: usize, seq: &[u8], rng: &mut R) -> Vec<u8> {
    const MIN_HOMOPOLYMER_LEN: usize = 4;
    const HOMOPOLYMER_PENALTY_MIN: f64 = 5.0;
    const HOMOPOLYMER_PENALTY_MAX: f64 = 10.0;
    const HOMOPOLYMER_FLOOR: f64 = 10.0;

    (0..length)
        .map(|pos| {
            // Q15–Q25 range; mode around Q20.
            let q: f64 = rng.gen_range(15.0..=25.0_f64);
            let q = if pos < seq.len() && is_homopolymer_run(seq, pos, MIN_HOMOPOLYMER_LEN) {
                let penalty: f64 = rng.gen_range(HOMOPOLYMER_PENALTY_MIN..=HOMOPOLYMER_PENALTY_MAX);
                (q - penalty).max(HOMOPOLYMER_FLOOR)
            } else {
                q
            };
            q.round() as u8
        })
        .collect()
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
    use rand::rngs::StdRng;
    use rand::SeedableRng;

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
                assert!((2..=41).contains(&q), "quality {} out of bounds", q);
            }
        }
    }

    #[test]
    fn test_pacbio_hifi_quality_length() {
        let mut rng = StdRng::seed_from_u64(42);
        let quals = sample_pacbio_hifi_qualities(1000, &mut rng);
        assert_eq!(quals.len(), 1000);
    }

    #[test]
    fn test_pacbio_hifi_quality_bounds() {
        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..100 {
            let quals = sample_pacbio_hifi_qualities(200, &mut rng);
            for &q in &quals {
                assert!(
                    (2..=41).contains(&q),
                    "PacBio HiFi quality {} out of bounds",
                    q
                );
            }
        }
    }

    #[test]
    fn test_pacbio_hifi_position_dependent() {
        // Mean quality at position 0 (near end) should be lower than at the midpoint.
        let length = 200usize;
        let n = 10_000;
        let mut q_at_0: f64 = 0.0;
        let mut q_at_mid: f64 = 0.0;
        let mid = length / 2;

        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..n {
            let quals = sample_pacbio_hifi_qualities(length, &mut rng);
            q_at_0 += quals[0] as f64;
            q_at_mid += quals[mid] as f64;
        }
        let mean_at_0 = q_at_0 / n as f64;
        let mean_at_mid = q_at_mid / n as f64;
        assert!(
            mean_at_0 < mean_at_mid,
            "mean quality at position 0 ({:.2}) should be less than at mid ({:.2})",
            mean_at_0,
            mean_at_mid
        );
    }

    #[test]
    fn test_is_homopolymer_run() {
        // "AAACGTTTTTAG": 3-bp poly-A at 0-2, 5-bp poly-T at 5-9.
        let seq = b"AAACGTTTTTAG";
        // Positions 0,1,2 are in a 3-bp poly-A — not >= 4.
        assert!(!is_homopolymer_run(seq, 0, 4));
        assert!(!is_homopolymer_run(seq, 2, 4));
        // Positions 6,7,8,9 are inside the 5-bp poly-T run.
        assert!(is_homopolymer_run(seq, 6, 4));
        assert!(is_homopolymer_run(seq, 9, 4));
        // G at position 4 is a run of 1.
        assert!(!is_homopolymer_run(seq, 4, 4));
        // Out-of-bounds position returns false.
        assert!(!is_homopolymer_run(seq, 100, 4));
    }

    #[test]
    fn test_nanopore_r10_homopolymer_quality_degradation() {
        // Build a 200-bp sequence: ACGT repeating as background (no homopolymer
        // runs), with positions 50-59 replaced by poly-A to create a 10-bp run.
        let mut seq: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        for b in &mut seq[50..60] {
            *b = b'A';
        }

        let n = 10_000;
        let mut q_homo: f64 = 0.0;
        let mut q_normal: f64 = 0.0;

        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..n {
            let quals = sample_nanopore_r10_qualities(200, &seq, &mut rng);
            // Mean over the 10-bp homopolymer region.
            for &q in &quals[50..60] {
                q_homo += q as f64;
            }
            // Mean over a plain region at positions 100-109.
            for &q in &quals[100..110] {
                q_normal += q as f64;
            }
        }
        let mean_homo = q_homo / (n * 10) as f64;
        let mean_normal = q_normal / (n * 10) as f64;
        assert!(
            mean_homo < mean_normal,
            "homopolymer quality ({:.2}) should be lower than normal region ({:.2})",
            mean_homo,
            mean_normal
        );
    }

    #[test]
    fn test_nanopore_r10_quality_range() {
        // Use a repeating ACGT sequence so no position is inside a homopolymer
        // run, isolating the base quality distribution from homopolymer penalties.
        let seq: Vec<u8> = (0..1000).map(|i| b"ACGT"[i % 4]).collect();
        let mut rng = StdRng::seed_from_u64(42);
        let quals = sample_nanopore_r10_qualities(1000, &seq, &mut rng);
        assert_eq!(quals.len(), 1000);
        // All qualities should be within valid Phred bounds.
        for &q in &quals {
            assert!(q >= 10, "Nanopore R10 quality {} below floor Q10", q);
            assert!(q <= 25, "Nanopore R10 quality {} above Q25", q);
        }
        let mean = quals.iter().map(|&q| q as f64).sum::<f64>() / quals.len() as f64;
        // Uniform over [15, 25] gives mean 20; allow ±1.5 with 1000 samples.
        assert!(
            (mean - 20.0).abs() < 1.5,
            "Nanopore R10 mean quality {} not near 20",
            mean
        );
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
            total_errors += seq
                .iter()
                .zip(original.iter())
                .filter(|(a, b)| a != b)
                .count();
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
