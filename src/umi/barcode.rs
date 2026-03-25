//! UMI barcode generation: random sequences, duplex UMI pairs, error injection, and RX tag formatting.

use rand::Rng;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Generate a random UMI barcode of the specified length.
pub fn generate_umi<R: Rng>(length: usize, rng: &mut R) -> Vec<u8> {
    (0..length).map(|_| BASES[rng.random_range(0..4)]).collect()
}

/// Generate a duplex UMI pair (AB, BA).
///
/// In duplex sequencing, each molecule gets two UMIs: A and B.
/// - Alpha strand read 1 carries A+B
/// - Beta strand read 1 carries B+A
pub fn generate_duplex_umi_pair<R: Rng>(length: usize, rng: &mut R) -> (Vec<u8>, Vec<u8>) {
    let umi_a = generate_umi(length, rng);
    let umi_b = generate_umi(length, rng);

    let mut alpha = Vec::with_capacity(length * 2);
    alpha.extend_from_slice(&umi_a);
    alpha.extend_from_slice(&umi_b);

    let mut beta = Vec::with_capacity(length * 2);
    beta.extend_from_slice(&umi_b);
    beta.extend_from_slice(&umi_a);

    (alpha, beta)
}

/// Inject sequencing errors into a UMI sequence.
// Called only in tests; not yet wired into the simulation write path.
#[cfg(test)]
pub fn inject_umi_errors<R: Rng>(umi: &mut [u8], error_rate: f64, rng: &mut R) {
    for base in umi.iter_mut() {
        if rng.random::<f64>() < error_rate {
            let original = *base;
            loop {
                let new_base = BASES[rng.random_range(0..4)];
                if new_base != original {
                    *base = new_base;
                    break;
                }
            }
        }
    }
}

/// Format a UMI for the RX BAM tag (SAM spec).
/// For duplex: "AAAA-BBBB" format.
// Called only in tests; production UMI tagging uses the RX tag inline.
#[cfg(test)]
pub fn format_rx_tag(umi: &[u8], umi_length: usize) -> String {
    if umi.len() == umi_length {
        String::from_utf8_lossy(umi).to_string()
    } else if umi.len() == umi_length * 2 {
        format!(
            "{}-{}",
            String::from_utf8_lossy(&umi[..umi_length]),
            String::from_utf8_lossy(&umi[umi_length..])
        )
    } else {
        String::from_utf8_lossy(umi).to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_generate_umi_length() {
        let mut rng = StdRng::seed_from_u64(42);
        let umi = generate_umi(8, &mut rng);
        assert_eq!(umi.len(), 8);
    }

    #[test]
    fn test_umi_valid_bases() {
        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..100 {
            let umi = generate_umi(12, &mut rng);
            for &b in &umi {
                assert!(matches!(b, b'A' | b'C' | b'G' | b'T'));
            }
        }
    }

    #[test]
    fn test_duplex_umi_pair() {
        let mut rng = StdRng::seed_from_u64(42);
        let (alpha, beta) = generate_duplex_umi_pair(8, &mut rng);

        assert_eq!(alpha.len(), 16);
        assert_eq!(beta.len(), 16);

        // Alpha = A+B, Beta = B+A
        assert_eq!(&alpha[..8], &beta[8..]);
        assert_eq!(&alpha[8..], &beta[..8]);
    }

    #[test]
    fn test_inject_umi_errors() {
        let mut rng = StdRng::seed_from_u64(42);
        let original = generate_umi(1000, &mut rng);
        let mut errored = original.clone();
        inject_umi_errors(&mut errored, 0.01, &mut rng);

        let errors = original
            .iter()
            .zip(errored.iter())
            .filter(|(a, b)| a != b)
            .count();
        let rate = errors as f64 / 1000.0;
        assert!(
            (rate - 0.01).abs() < 0.01,
            "error rate {} should be close to 0.01",
            rate
        );
    }

    #[test]
    fn test_format_rx_simplex() {
        let umi = b"ACGTACGT";
        let tag = format_rx_tag(umi, 8);
        assert_eq!(tag, "ACGTACGT");
    }

    #[test]
    fn test_format_rx_duplex() {
        let umi = b"ACGTACGTTTTTAAAA";
        let tag = format_rx_tag(umi, 8);
        assert_eq!(tag, "ACGTACGT-TTTTAAAA");
    }

    #[test]
    fn test_deterministic_umi() {
        let mut rng1 = StdRng::seed_from_u64(42);
        let mut rng2 = StdRng::seed_from_u64(42);
        assert_eq!(generate_umi(8, &mut rng1), generate_umi(8, &mut rng2));
    }
}
