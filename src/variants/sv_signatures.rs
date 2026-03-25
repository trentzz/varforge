//! SV signature generators for known genomic instability phenotypes.
//!
//! Each generator produces a set of structural variant `Variant` records
//! with the size and type distribution characteristic of the phenotype.

use rand::Rng;
use rand_distr::{Distribution, LogNormal};

use crate::core::types::{MutationType, SvType, Variant};

/// Generate HRD-typical large deletions (100 kbp – 10 Mbp).
///
/// HRD deletions are distributed across all chromosomes; each is assigned
/// a random genomic position. LOH is indicated by setting the clone_id to
/// `"hrd_loh"`.
pub fn generate_hrd_deletions<R: Rng>(
    chrom_lengths: &[(String, u64)],
    count: usize,
    rng: &mut R,
) -> Vec<Variant> {
    let mut variants = Vec::with_capacity(count);

    if chrom_lengths.is_empty() || count == 0 {
        return variants;
    }

    // HRD deletions: log-normal size distribution centred around ~1 Mbp.
    // LogNormal(ln(500_000), 1.2) gives a range roughly from 100 kbp to 10 Mbp.
    let size_dist = LogNormal::new(13.1, 1.2).expect("valid log-normal params");

    for _ in 0..count {
        // Pick a random chromosome weighted by length.
        let (chrom, chrom_len) = match pick_chrom(chrom_lengths, rng) {
            Some(v) => v,
            None => continue,
        };

        // Sample size clamped to [100_000, 10_000_000].
        let size = (size_dist.sample(rng) as u64).clamp(100_000, 10_000_000);

        if chrom_len <= size + 1 {
            continue; // Chromosome too small for this deletion.
        }

        let start = rng.gen_range(0..chrom_len - size);
        let end = start + size;

        variants.push(Variant {
            chrom: chrom.clone(),
            mutation: MutationType::Sv {
                sv_type: SvType::Deletion,
                chrom,
                start,
                end,
            },
            expected_vaf: 0.5,
            clone_id: Some("hrd_loh".to_string()),
            haplotype: None,
            ccf: None,
        });
    }

    variants
}

/// Generate TDP-typical short tandem duplications (1 kbp – 10 kbp).
///
/// The tandem duplicator phenotype produces an elevated rate of tandem
/// duplications in the 1–10 kbp range.
pub fn generate_tdp_duplications<R: Rng>(
    chrom_lengths: &[(String, u64)],
    count: usize,
    rng: &mut R,
) -> Vec<Variant> {
    let mut variants = Vec::with_capacity(count);

    if chrom_lengths.is_empty() || count == 0 {
        return variants;
    }

    // TDP duplications: log-normal size distribution centred around ~3 kbp.
    let size_dist = LogNormal::new(8.0, 0.8).expect("valid log-normal params");

    for _ in 0..count {
        let (chrom, chrom_len) = match pick_chrom(chrom_lengths, rng) {
            Some(v) => v,
            None => continue,
        };

        let size = (size_dist.sample(rng) as u64).clamp(1_000, 10_000);

        if chrom_len <= size + 1 {
            continue;
        }

        let start = rng.gen_range(0..chrom_len - size);
        let end = start + size;

        variants.push(Variant {
            chrom: chrom.clone(),
            mutation: MutationType::Sv {
                sv_type: SvType::Duplication,
                chrom,
                start,
                end,
            },
            expected_vaf: 0.5,
            clone_id: Some("tdp".to_string()),
            haplotype: None,
            ccf: None,
        });
    }

    variants
}

/// Generate chromothripsis-pattern rearrangements.
///
/// Chromothripsis produces multiple rearrangements of a single chromosomal
/// region: a mix of deletions, inversions, and duplications from a localised
/// "shattering" event.
pub fn generate_chromothripsis<R: Rng>(
    chrom_lengths: &[(String, u64)],
    count: usize,
    rng: &mut R,
) -> Vec<Variant> {
    let mut variants = Vec::with_capacity(count);

    if chrom_lengths.is_empty() || count == 0 {
        return variants;
    }

    // Pick one chromosome to shatter.
    let (shattered_chrom, chrom_len) = match pick_chrom(chrom_lengths, rng) {
        Some(v) => v,
        None => return variants,
    };

    // The shattering window is a 10–50 Mbp region.
    let window_size = rng.gen_range(10_000_000u64..=50_000_000).min(chrom_len);
    let window_start = if chrom_len > window_size {
        rng.gen_range(0..chrom_len - window_size)
    } else {
        0
    };
    let window_end = window_start + window_size;

    // The window must be large enough to place at least one rearrangement.
    // We need window_start < window_end - 1000 for the range to be non-empty.
    if window_size < 2_000 {
        return variants;
    }

    // Generate `count` rearrangements within the shattering window.
    for _ in 0..count {
        let start = rng.gen_range(window_start..window_end.saturating_sub(1_000));
        let size = rng
            .gen_range(100_000u64..=5_000_000)
            .min(window_end - start);
        let end = start + size;

        // Randomly pick an SV type.
        let sv_type = match rng.gen_range(0..3u8) {
            0 => SvType::Deletion,
            1 => SvType::Inversion,
            _ => SvType::Duplication,
        };

        variants.push(Variant {
            chrom: shattered_chrom.clone(),
            mutation: MutationType::Sv {
                sv_type,
                chrom: shattered_chrom.clone(),
                start,
                end,
            },
            expected_vaf: 0.5,
            clone_id: Some("chromothripsis".to_string()),
            haplotype: None,
            ccf: None,
        });
    }

    variants
}

/// Pick a chromosome weighted by length.
///
/// Returns `None` if `chrom_lengths` is empty or all lengths are zero.
fn pick_chrom<R: Rng>(chrom_lengths: &[(String, u64)], rng: &mut R) -> Option<(String, u64)> {
    let total_len: u64 = chrom_lengths.iter().map(|(_, l)| l).sum();
    if total_len == 0 {
        return None;
    }
    let pick = rng.gen_range(0..total_len);
    let mut cumulative = 0u64;
    for (chrom, len) in chrom_lengths {
        cumulative += len;
        if pick < cumulative {
            return Some((chrom.clone(), *len));
        }
    }
    // Fallback: last chromosome.
    chrom_lengths.last().map(|(c, l)| (c.clone(), *l))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn test_chroms() -> Vec<(String, u64)> {
        vec![
            ("chr1".to_string(), 248_956_422),
            ("chr2".to_string(), 242_193_529),
            ("chr3".to_string(), 198_295_559),
        ]
    }

    #[test]
    fn test_hrd_deletions_count() {
        let mut rng = StdRng::seed_from_u64(42);
        let chroms = test_chroms();
        let variants = generate_hrd_deletions(&chroms, 10, &mut rng);
        // May be fewer than 10 if chromosomes are too small, but for realistic
        // sizes all 10 should be generated.
        assert!(!variants.is_empty());
        assert!(variants.len() <= 10);
    }

    #[test]
    fn test_hrd_deletions_size_range() {
        let mut rng = StdRng::seed_from_u64(1);
        let chroms = test_chroms();
        let variants = generate_hrd_deletions(&chroms, 100, &mut rng);
        for v in &variants {
            if let MutationType::Sv { start, end, .. } = &v.mutation {
                let size = end - start;
                assert!(
                    (100_000..=10_000_000).contains(&size),
                    "HRD deletion size {size} out of expected range"
                );
            }
        }
    }

    #[test]
    fn test_tdp_duplications_size_range() {
        let mut rng = StdRng::seed_from_u64(2);
        let chroms = test_chroms();
        let variants = generate_tdp_duplications(&chroms, 100, &mut rng);
        for v in &variants {
            if let MutationType::Sv {
                sv_type,
                start,
                end,
                ..
            } = &v.mutation
            {
                assert_eq!(*sv_type, SvType::Duplication);
                let size = end - start;
                assert!(
                    (1_000..=10_000).contains(&size),
                    "TDP duplication size {size} out of expected range"
                );
            }
        }
    }

    #[test]
    fn test_chromothripsis_single_chrom() {
        let mut rng = StdRng::seed_from_u64(3);
        let chroms = test_chroms();
        let variants = generate_chromothripsis(&chroms, 20, &mut rng);
        assert!(!variants.is_empty());
        // All variants must be on the same chromosome (the shattered one).
        let chrom_name = &variants[0].chrom;
        for v in &variants {
            assert_eq!(
                &v.chrom, chrom_name,
                "chromothripsis variants must all be on the same chromosome"
            );
        }
    }

    #[test]
    fn test_empty_inputs() {
        let mut rng = StdRng::seed_from_u64(0);
        assert!(generate_hrd_deletions(&[], 10, &mut rng).is_empty());
        assert!(generate_tdp_duplications(&[], 10, &mut rng).is_empty());
        assert!(generate_chromothripsis(&[], 10, &mut rng).is_empty());

        let chroms = test_chroms();
        assert!(generate_hrd_deletions(&chroms, 0, &mut rng).is_empty());
        assert!(generate_tdp_duplications(&chroms, 0, &mut rng).is_empty());
        assert!(generate_chromothripsis(&chroms, 0, &mut rng).is_empty());
    }

    #[test]
    fn test_clone_ids() {
        let mut rng = StdRng::seed_from_u64(99);
        let chroms = test_chroms();

        let hrd = generate_hrd_deletions(&chroms, 5, &mut rng);
        for v in &hrd {
            assert_eq!(v.clone_id.as_deref(), Some("hrd_loh"));
        }

        let tdp = generate_tdp_duplications(&chroms, 5, &mut rng);
        for v in &tdp {
            assert_eq!(v.clone_id.as_deref(), Some("tdp"));
        }

        let ct = generate_chromothripsis(&chroms, 5, &mut rng);
        for v in &ct {
            assert_eq!(v.clone_id.as_deref(), Some("chromothripsis"));
        }
    }
}
