use rand::Rng;
use rand_distr::{Distribution, LogNormal};

use crate::core::types::ReadPair;

/// A UMI family: one original molecule and its PCR duplicates.
#[derive(Debug, Clone)]
pub struct UmiFamily {
    #[allow(dead_code)]
    pub umi: Vec<u8>,
    pub original: ReadPair,
    pub family_size: usize,
}

/// Generate PCR duplicate copies of a read pair within a UMI family.
///
/// All copies share the same alignment position and UMI.
/// PCR errors can be injected per-cycle to model error accumulation.
pub fn generate_pcr_copies(
    family: &UmiFamily,
    pcr_error_rate: f64,
    _pcr_cycles: u32,
    rng: &mut impl Rng,
) -> Vec<ReadPair> {
    let mut copies = Vec::with_capacity(family.family_size);

    for i in 0..family.family_size {
        let mut copy = family.original.clone();
        copy.name = format!("{}:pcr:{}", family.original.name, i);

        // Inject PCR errors (simplified: per-base rate across all copies)
        if pcr_error_rate > 0.0 {
            inject_pcr_errors(&mut copy.read1.seq, pcr_error_rate, rng);
            inject_pcr_errors(&mut copy.read2.seq, pcr_error_rate, rng);
        }

        copies.push(copy);
    }

    copies
}

fn inject_pcr_errors(seq: &mut [u8], rate: f64, rng: &mut impl Rng) {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    for base in seq.iter_mut() {
        if rng.gen::<f64>() < rate {
            let original = *base;
            loop {
                let new_base = BASES[rng.gen_range(0..4)];
                if new_base != original {
                    *base = new_base;
                    break;
                }
            }
        }
    }
}

/// Sample a family size from a log-normal distribution.
pub fn sample_family_size(mean: f64, sd: f64, rng: &mut impl Rng) -> usize {
    let variance = sd * sd;
    let mu = (mean * mean / (mean * mean + variance).sqrt()).ln();
    let sigma = (1.0 + variance / (mean * mean)).ln().sqrt();
    let dist = LogNormal::new(mu, sigma).expect("invalid lognormal params");
    let size = dist.sample(rng).round() as usize;
    size.max(1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::types::Read;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn make_family(family_size: usize) -> UmiFamily {
        UmiFamily {
            umi: b"ACGTACGT".to_vec(),
            original: ReadPair {
                name: "read_001".to_string(),
                read1: Read::new(vec![b'A'; 100], vec![30; 100]),
                read2: Read::new(vec![b'T'; 100], vec![30; 100]),
                fragment_start: 1000,
                fragment_length: 300,
                chrom: "chr1".to_string(),
            },
            family_size,
        }
    }

    #[test]
    fn test_generate_pcr_copies_count() {
        let mut rng = StdRng::seed_from_u64(42);
        let family = make_family(5);
        let copies = generate_pcr_copies(&family, 0.0, 10, &mut rng);
        assert_eq!(copies.len(), 5);
    }

    #[test]
    fn test_pcr_copies_same_position() {
        let mut rng = StdRng::seed_from_u64(42);
        let family = make_family(3);
        let copies = generate_pcr_copies(&family, 0.0, 10, &mut rng);
        for copy in &copies {
            assert_eq!(copy.fragment_start, family.original.fragment_start);
            assert_eq!(copy.chrom, family.original.chrom);
        }
    }

    #[test]
    fn test_pcr_errors_injected() {
        let mut rng = StdRng::seed_from_u64(42);
        let family = make_family(100);
        let copies = generate_pcr_copies(&family, 0.01, 10, &mut rng);

        let total_errors: usize = copies.iter().map(|c| {
            c.read1.seq.iter().zip(family.original.read1.seq.iter())
                .filter(|(a, b)| a != b).count()
        }).sum();

        // With 100 copies * 100 bases * 0.01 rate ≈ 100 errors
        assert!(total_errors > 0, "should have some PCR errors");
    }

    #[test]
    fn test_no_pcr_errors_when_rate_zero() {
        let mut rng = StdRng::seed_from_u64(42);
        let family = make_family(10);
        let copies = generate_pcr_copies(&family, 0.0, 10, &mut rng);
        for copy in &copies {
            assert_eq!(copy.read1.seq, family.original.read1.seq);
            assert_eq!(copy.read2.seq, family.original.read2.seq);
        }
    }

    #[test]
    fn test_sample_family_size() {
        let mut rng = StdRng::seed_from_u64(42);
        let sizes: Vec<usize> = (0..10000).map(|_| sample_family_size(3.0, 1.5, &mut rng)).collect();

        assert!(sizes.iter().all(|&s| s >= 1));
        let mean = sizes.iter().sum::<usize>() as f64 / sizes.len() as f64;
        assert!((mean - 3.0).abs() < 1.0, "mean {} too far from 3.0", mean);
    }

    #[test]
    fn test_unique_copy_names() {
        let mut rng = StdRng::seed_from_u64(42);
        let family = make_family(5);
        let copies = generate_pcr_copies(&family, 0.0, 10, &mut rng);
        let names: std::collections::HashSet<_> = copies.iter().map(|c| &c.name).collect();
        assert_eq!(names.len(), 5, "each copy should have a unique name");
    }
}
