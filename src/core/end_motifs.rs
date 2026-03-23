//! Empirical 4-mer end motif frequencies for plasma cfDNA fragments.
//!
//! Source: Chandrananda et al. (2015) and Chan et al. (2016) plasma cfDNA
//! studies. Values are approximate relative frequencies. Motifs with elevated
//! plasma frequency are listed explicitly; all others fall back to a uniform
//! baseline.

use std::collections::HashMap;
use std::sync::OnceLock;

use rand::Rng;

static PLASMA_TABLE: OnceLock<HashMap<[u8; 4], f64>> = OnceLock::new();

fn plasma_table() -> &'static HashMap<[u8; 4], f64> {
    PLASMA_TABLE.get_or_init(build_plasma_table)
}

/// Build the plasma cfDNA 4-mer end motif frequency table.
///
/// Values are approximate relative frequencies from Chandrananda et al. 2015
/// and Chan et al. 2016. Motifs with high plasma prevalence are listed
/// explicitly. All other 4-mers fall back to a 1/256 baseline.
fn build_plasma_table() -> HashMap<[u8; 4], f64> {
    // Most frequent plasma cfDNA end motifs. Frequencies are relative and
    // approximate; they are not required to sum to 1.
    let entries: &[([u8; 4], f64)] = &[
        (*b"CCCA", 0.040),
        (*b"CCCG", 0.035),
        (*b"CCCT", 0.038),
        (*b"CCCC", 0.032),
        (*b"TCCC", 0.030),
        (*b"ACCC", 0.028),
        (*b"GCCC", 0.027),
        (*b"TTTC", 0.025),
        (*b"TTCA", 0.024),
        (*b"TTCG", 0.022),
        (*b"CTTT", 0.020),
        (*b"ATTT", 0.019),
        (*b"GTTT", 0.018),
        (*b"TTTT", 0.017),
        (*b"CAAA", 0.016),
        (*b"TAAA", 0.015),
        (*b"GAAA", 0.015),
        (*b"AAAA", 0.014),
        (*b"CGGG", 0.013),
        (*b"TGGG", 0.012),
        (*b"AGGG", 0.012),
        (*b"GGGG", 0.011),
    ];
    let mut map: HashMap<[u8; 4], f64> = HashMap::new();
    for &(motif, freq) in entries {
        map.insert(motif, freq);
    }
    map
}

/// Approximate maximum single-motif frequency in the plasma table.
///
/// Used as the upper bound for rejection sampling.
const MAX_PLASMA_FREQ: f64 = 0.040;

/// Relative frequency of a plasma cfDNA 4-mer end motif.
///
/// Returns the tabulated frequency for a known motif, or a baseline of
/// 1/256 for motifs not in the table.
pub fn plasma_end_motif_freq(motif: &[u8; 4]) -> f64 {
    plasma_table().get(motif).copied().unwrap_or(1.0 / 256.0)
}

/// Accept or reject a fragment based on its 5' end motif via rejection sampling.
///
/// Returns `true` if the fragment should be kept. Callers should discard the
/// fragment and resample when this returns `false`. When `motif_5p` is `None`
/// (e.g. near a chromosome boundary where the motif cannot be extracted), the
/// function always accepts.
pub fn accept_fragment_by_end_motif<R: Rng>(motif_5p: Option<[u8; 4]>, rng: &mut R) -> bool {
    match motif_5p {
        None => true,
        Some(ref m) => {
            let freq = plasma_end_motif_freq(m);
            rng.gen::<f64>() < freq / MAX_PLASMA_FREQ
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_known_motif_frequency() {
        // CCCA is the highest-frequency motif in the table.
        let freq = plasma_end_motif_freq(b"CCCA");
        assert!(
            (freq - 0.040).abs() < f64::EPSILON,
            "expected 0.040 for CCCA, got {freq}"
        );
    }

    #[test]
    fn test_unknown_motif_baseline() {
        // A motif not in the table should return the baseline 1/256.
        let freq = plasma_end_motif_freq(b"ACGT");
        let baseline = 1.0 / 256.0;
        assert!(
            (freq - baseline).abs() < 1e-10,
            "expected baseline {baseline} for unknown motif, got {freq}"
        );
    }

    #[test]
    fn test_none_motif_always_accepts() {
        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..100 {
            assert!(
                accept_fragment_by_end_motif(None, &mut rng),
                "None motif must always accept"
            );
        }
    }

    #[test]
    fn test_high_freq_motif_accepts_often() {
        let mut rng = StdRng::seed_from_u64(42);
        let motif = Some(*b"CCCA"); // freq = MAX_PLASMA_FREQ → accept rate = 1.0
        let accepted = (0..100)
            .filter(|_| accept_fragment_by_end_motif(motif, &mut rng))
            .count();
        // All should be accepted since freq / MAX_PLASMA_FREQ = 1.0.
        assert_eq!(
            accepted, 100,
            "CCCA at max freq should accept all fragments"
        );
    }

    #[test]
    fn test_low_freq_motif_rejects_often() {
        let mut rng = StdRng::seed_from_u64(1);
        let motif = Some(*b"ACGT"); // baseline freq ≈ 1/256 → very low accept rate
        let accepted = (0..1000)
            .filter(|_| accept_fragment_by_end_motif(motif, &mut rng))
            .count();
        // Baseline freq is 1/256 ≈ 0.0039. Accept rate = freq / MAX_PLASMA_FREQ
        // ≈ 0.097, so about 97/1000 expected. Verify it is well below 50% (500/1000).
        assert!(
            accepted < 200,
            "low-freq motif should rarely be accepted, got {accepted}/1000"
        );
    }
}
