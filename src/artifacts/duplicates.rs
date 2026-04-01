//! PCR duplicate simulation: selects reads for duplication and creates copies with identical alignment coordinates.

use crate::core::types::ReadPair;
use rand::Rng;

/// Select reads to duplicate based on a target duplication rate.
///
/// Returns indices of reads that should be duplicated.
pub fn select_duplicates<R: Rng>(
    read_count: usize,
    duplicate_rate: f64,
    rng: &mut R,
) -> Vec<usize> {
    if read_count == 0 {
        return Vec::new();
    }
    let n_dups = (read_count as f64 * duplicate_rate).round() as usize;
    let mut indices = Vec::with_capacity(n_dups);
    for _ in 0..n_dups {
        indices.push(rng.random_range(0..read_count));
    }
    indices
}

/// Create a duplicate of a read pair with the same alignment coordinates.
pub fn duplicate_read_pair(original: &ReadPair, dup_index: usize) -> ReadPair {
    let mut dup = original.clone();
    dup.name = format!("{}:dup:{}", original.name, dup_index);
    dup
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::types::Read;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_select_duplicates_count() {
        let mut rng = StdRng::seed_from_u64(42);
        let indices = select_duplicates(1000, 0.15, &mut rng);
        assert_eq!(indices.len(), 150);
    }

    #[test]
    fn test_select_duplicates_in_range() {
        let mut rng = StdRng::seed_from_u64(42);
        let indices = select_duplicates(100, 0.1, &mut rng);
        for &idx in &indices {
            assert!(idx < 100);
        }
    }

    #[test]
    fn test_duplicate_read_pair() {
        let original = ReadPair {
            name: "read_001".to_string(),
            read1: Read::new(vec![b'A'; 100], vec![30; 100]),
            read2: Read::new(vec![b'T'; 100], vec![30; 100]),
            fragment_start: 1000,
            fragment_length: 300,
            chrom: "chr1".to_string(),
            variant_tags: Vec::new(),
            ref_seq_r1: Vec::new(),
            ref_seq_r2: Vec::new(),
            inline_prefix_r1: None,
            inline_prefix_r2: None,
        };

        let dup = duplicate_read_pair(&original, 0);
        assert_eq!(dup.fragment_start, original.fragment_start);
        assert_eq!(dup.read1.seq, original.read1.seq);
        assert_ne!(dup.name, original.name);
    }

    #[test]
    fn test_zero_duplicate_rate() {
        let mut rng = StdRng::seed_from_u64(42);
        let indices = select_duplicates(1000, 0.0, &mut rng);
        assert!(indices.is_empty());
    }
}
