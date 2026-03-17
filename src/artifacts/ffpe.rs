use rand::Rng;

/// Inject FFPE (formalin-fixed, paraffin-embedded) deamination artifacts.
///
/// FFPE damage causes C>T transitions on the affected strand (observed as G>A on complement).
/// Damage rate is typically higher near fragment ends.
pub fn inject_ffpe_damage<R: Rng>(
    seq: &mut [u8],
    damage_rate: f64,
    is_reverse_strand: bool,
    rng: &mut R,
) {
    let seq_len = seq.len();
    for (pos, base) in seq.iter_mut().enumerate() {
        // Damage is slightly higher near fragment ends
        let position_factor = end_enrichment(pos, seq_len);
        let effective_rate = damage_rate * position_factor;

        if is_reverse_strand {
            // On reverse strand, deamination of C on the original strand
            // appears as G>A
            if *base == b'G' && rng.gen::<f64>() < effective_rate {
                *base = b'A';
            }
        } else {
            // On forward strand, C>T deamination
            if *base == b'C' && rng.gen::<f64>() < effective_rate {
                *base = b'T';
            }
        }
    }
}

/// Inject oxidative damage (8-oxoguanine / oxoG artifacts).
///
/// oxoG causes G>T transversions, predominantly on one strand.
pub fn inject_oxog_damage<R: Rng>(seq: &mut [u8], damage_rate: f64, is_read1: bool, rng: &mut R) {
    // oxoG artifacts show strand asymmetry: primarily affect read 1
    let effective_rate = if is_read1 {
        damage_rate
    } else {
        damage_rate * 0.1
    };

    for base in seq.iter_mut() {
        if *base == b'G' && rng.gen::<f64>() < effective_rate {
            *base = b'T';
        }
    }
}

/// Enrichment factor for damage near fragment ends.
/// Returns 1.0-3.0, with higher values near ends.
fn end_enrichment(pos: usize, length: usize) -> f64 {
    if length == 0 {
        return 1.0;
    }
    let dist_from_end = pos.min(length - 1 - pos) as f64;
    let norm_dist = dist_from_end / (length as f64 / 2.0);
    // Exponential enrichment toward ends
    1.0 + 2.0 * (-3.0 * norm_dist).exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_ffpe_forward_strand_c_to_t() {
        let mut rng = StdRng::seed_from_u64(42);
        let original = vec![b'C'; 10000];
        let mut seq = original.clone();
        inject_ffpe_damage(&mut seq, 0.05, false, &mut rng);

        let mutations: usize = seq
            .iter()
            .zip(original.iter())
            .filter(|(a, b)| a != b)
            .count();
        assert!(mutations > 0, "should have some C>T damage");

        // All mutations should be C>T
        for (new, old) in seq.iter().zip(original.iter()) {
            if new != old {
                assert_eq!(*old, b'C');
                assert_eq!(*new, b'T');
            }
        }
    }

    #[test]
    fn test_ffpe_reverse_strand_g_to_a() {
        let mut rng = StdRng::seed_from_u64(42);
        let original = vec![b'G'; 10000];
        let mut seq = original.clone();
        inject_ffpe_damage(&mut seq, 0.05, true, &mut rng);

        let mutations: usize = seq
            .iter()
            .zip(original.iter())
            .filter(|(a, b)| a != b)
            .count();
        assert!(mutations > 0, "should have some G>A damage");

        for (new, old) in seq.iter().zip(original.iter()) {
            if new != old {
                assert_eq!(*old, b'G');
                assert_eq!(*new, b'A');
            }
        }
    }

    #[test]
    fn test_ffpe_no_damage_at_zero_rate() {
        let mut rng = StdRng::seed_from_u64(42);
        let original = vec![b'C'; 1000];
        let mut seq = original.clone();
        inject_ffpe_damage(&mut seq, 0.0, false, &mut rng);
        assert_eq!(seq, original);
    }

    #[test]
    fn test_ffpe_only_affects_c_forward() {
        let mut rng = StdRng::seed_from_u64(42);
        let original = vec![b'A'; 1000];
        let mut seq = original.clone();
        inject_ffpe_damage(&mut seq, 1.0, false, &mut rng);
        assert_eq!(
            seq, original,
            "A bases should not be affected by FFPE on forward strand"
        );
    }

    #[test]
    fn test_oxog_g_to_t() {
        let mut rng = StdRng::seed_from_u64(42);
        let original = vec![b'G'; 10000];
        let mut seq = original.clone();
        inject_oxog_damage(&mut seq, 0.05, true, &mut rng);

        for (new, old) in seq.iter().zip(original.iter()) {
            if new != old {
                assert_eq!(*old, b'G');
                assert_eq!(*new, b'T');
            }
        }
    }

    #[test]
    fn test_oxog_strand_asymmetry() {
        let mut rng1 = StdRng::seed_from_u64(42);
        let mut rng2 = StdRng::seed_from_u64(42);

        let mut seq_r1 = vec![b'G'; 100000];
        let mut seq_r2 = vec![b'G'; 100000];
        inject_oxog_damage(&mut seq_r1, 0.05, true, &mut rng1);
        inject_oxog_damage(&mut seq_r2, 0.05, false, &mut rng2);

        let r1_mutations = seq_r1.iter().filter(|&&b| b == b'T').count();
        let r2_mutations = seq_r2.iter().filter(|&&b| b == b'T').count();

        assert!(
            r1_mutations > r2_mutations * 5,
            "read1 ({}) should have much more oxoG than read2 ({})",
            r1_mutations,
            r2_mutations
        );
    }

    #[test]
    fn test_end_enrichment() {
        let start = end_enrichment(0, 100);
        let middle = end_enrichment(50, 100);
        assert!(
            start > middle,
            "ends ({}) should have more damage than middle ({})",
            start,
            middle
        );
    }
}
