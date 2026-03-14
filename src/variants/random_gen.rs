use rand::Rng;
use crate::core::types::{MutationType, Region, Variant};

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Generate random mutations distributed across the given regions.
pub fn generate_random_mutations<R: Rng>(
    regions: &[Region],
    count: usize,
    vaf_min: f64,
    vaf_max: f64,
    snv_fraction: f64,
    indel_fraction: f64,
    _mnv_fraction: f64,
    reference_lookup: &dyn Fn(&str, u64) -> Option<u8>,
    rng: &mut R,
) -> Vec<Variant> {
    let total_length: u64 = regions.iter().map(|r| r.len()).sum();
    if total_length == 0 || count == 0 {
        return Vec::new();
    }

    let mut variants = Vec::with_capacity(count);

    for _ in 0..count {
        // Pick a random position within the regions
        let rand_offset: u64 = rng.gen_range(0..total_length);
        let (chrom, pos) = offset_to_position(regions, rand_offset);

        // Pick mutation type
        let type_roll: f64 = rng.gen();
        let mutation = if type_roll < snv_fraction {
            generate_random_snv(&chrom, pos, reference_lookup, rng)
        } else if type_roll < snv_fraction + indel_fraction {
            generate_random_indel(&chrom, pos, reference_lookup, rng)
        } else {
            generate_random_mnv(&chrom, pos, reference_lookup, rng)
        };

        if let Some(mutation) = mutation {
            let vaf = rng.gen_range(vaf_min..vaf_max);
            variants.push(Variant {
                chrom,
                mutation,
                expected_vaf: vaf,
                clone_id: None,
            });
        }
    }

    variants
}

fn offset_to_position(regions: &[Region], offset: u64) -> (String, u64) {
    let mut remaining = offset;
    for region in regions {
        let len = region.len();
        if remaining < len {
            return (region.chrom.clone(), region.start + remaining);
        }
        remaining -= len;
    }
    // Fallback to last region
    let last = regions.last().unwrap();
    (last.chrom.clone(), last.end - 1)
}

fn generate_random_snv<R: Rng>(
    _chrom: &str,
    pos: u64,
    reference_lookup: &dyn Fn(&str, u64) -> Option<u8>,
    rng: &mut R,
) -> Option<MutationType> {
    let ref_base = reference_lookup(_chrom, pos).unwrap_or(b'N');
    if ref_base == b'N' {
        return None;
    }
    let alt_base = loop {
        let b = BASES[rng.gen_range(0..4)];
        if b != ref_base {
            break b;
        }
    };
    Some(MutationType::Snv { pos, ref_base, alt_base })
}

fn generate_random_indel<R: Rng>(
    _chrom: &str,
    pos: u64,
    reference_lookup: &dyn Fn(&str, u64) -> Option<u8>,
    rng: &mut R,
) -> Option<MutationType> {
    let ref_base = reference_lookup(_chrom, pos).unwrap_or(b'N');
    if ref_base == b'N' {
        return None;
    }

    let is_insertion = rng.gen_bool(0.5);
    let indel_len: usize = rng.gen_range(1..=10);

    if is_insertion {
        let mut alt_seq = vec![ref_base];
        for _ in 0..indel_len {
            alt_seq.push(BASES[rng.gen_range(0..4)]);
        }
        Some(MutationType::Indel {
            pos,
            ref_seq: vec![ref_base],
            alt_seq,
        })
    } else {
        // Deletion: need reference bases
        let mut ref_seq = vec![ref_base];
        for i in 1..=indel_len {
            let b = reference_lookup(_chrom, pos + i as u64).unwrap_or(b'N');
            if b == b'N' {
                return None;
            }
            ref_seq.push(b);
        }
        Some(MutationType::Indel {
            pos,
            ref_seq,
            alt_seq: vec![ref_base],
        })
    }
}

fn generate_random_mnv<R: Rng>(
    _chrom: &str,
    pos: u64,
    reference_lookup: &dyn Fn(&str, u64) -> Option<u8>,
    rng: &mut R,
) -> Option<MutationType> {
    let mnv_len: usize = rng.gen_range(2..=3);
    let mut ref_seq = Vec::with_capacity(mnv_len);
    let mut alt_seq = Vec::with_capacity(mnv_len);

    for i in 0..mnv_len {
        let ref_base = reference_lookup(_chrom, pos + i as u64).unwrap_or(b'N');
        if ref_base == b'N' {
            return None;
        }
        ref_seq.push(ref_base);
        let alt_base = loop {
            let b = BASES[rng.gen_range(0..4)];
            if b != ref_base {
                break b;
            }
        };
        alt_seq.push(alt_base);
    }

    Some(MutationType::Mnv { pos, ref_seq, alt_seq })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn simple_ref_lookup(chrom: &str, pos: u64) -> Option<u8> {
        // Return a deterministic base based on position
        let _ = chrom;
        Some(BASES[(pos % 4) as usize])
    }

    #[test]
    fn test_generate_random_mutations_count() {
        let mut rng = StdRng::seed_from_u64(42);
        let regions = vec![Region::new("chr1", 0, 10000)];
        let variants = generate_random_mutations(
            &regions, 100, 0.01, 0.5, 0.8, 0.15, 0.05,
            &simple_ref_lookup, &mut rng,
        );
        assert_eq!(variants.len(), 100);
    }

    #[test]
    fn test_vaf_range() {
        let mut rng = StdRng::seed_from_u64(42);
        let regions = vec![Region::new("chr1", 0, 10000)];
        let variants = generate_random_mutations(
            &regions, 1000, 0.001, 0.05, 0.8, 0.15, 0.05,
            &simple_ref_lookup, &mut rng,
        );
        for v in &variants {
            assert!(v.expected_vaf >= 0.001, "VAF {} below min", v.expected_vaf);
            assert!(v.expected_vaf < 0.05, "VAF {} above max", v.expected_vaf);
        }
    }

    #[test]
    fn test_mutation_type_distribution() {
        let mut rng = StdRng::seed_from_u64(42);
        let regions = vec![Region::new("chr1", 0, 100000)];
        let variants = generate_random_mutations(
            &regions, 10000, 0.01, 0.5, 0.8, 0.15, 0.05,
            &simple_ref_lookup, &mut rng,
        );

        let snv_count = variants.iter().filter(|v| matches!(v.mutation, MutationType::Snv { .. })).count();
        let indel_count = variants.iter().filter(|v| matches!(v.mutation, MutationType::Indel { .. })).count();
        let mnv_count = variants.iter().filter(|v| matches!(v.mutation, MutationType::Mnv { .. })).count();

        let snv_frac = snv_count as f64 / variants.len() as f64;
        let indel_frac = indel_count as f64 / variants.len() as f64;
        let mnv_frac = mnv_count as f64 / variants.len() as f64;

        assert!((snv_frac - 0.8).abs() < 0.05, "SNV fraction {} far from 0.8", snv_frac);
        assert!((indel_frac - 0.15).abs() < 0.05, "indel fraction {} far from 0.15", indel_frac);
        assert!((mnv_frac - 0.05).abs() < 0.03, "MNV fraction {} far from 0.05", mnv_frac);
    }

    #[test]
    fn test_snv_different_from_ref() {
        let mut rng = StdRng::seed_from_u64(42);
        let regions = vec![Region::new("chr1", 0, 10000)];
        let variants = generate_random_mutations(
            &regions, 1000, 0.01, 0.5, 1.0, 0.0, 0.0,
            &simple_ref_lookup, &mut rng,
        );
        for v in &variants {
            if let MutationType::Snv { ref_base, alt_base, .. } = v.mutation {
                assert_ne!(ref_base, alt_base, "alt must differ from ref");
            }
        }
    }

    #[test]
    fn test_empty_regions() {
        let mut rng = StdRng::seed_from_u64(42);
        let variants = generate_random_mutations(
            &[], 100, 0.01, 0.5, 0.8, 0.15, 0.05,
            &simple_ref_lookup, &mut rng,
        );
        assert!(variants.is_empty());
    }

    #[test]
    fn test_deterministic_with_seed() {
        let regions = vec![Region::new("chr1", 0, 10000)];
        let mut rng1 = StdRng::seed_from_u64(99);
        let mut rng2 = StdRng::seed_from_u64(99);

        let v1 = generate_random_mutations(
            &regions, 50, 0.01, 0.5, 0.8, 0.15, 0.05,
            &simple_ref_lookup, &mut rng1,
        );
        let v2 = generate_random_mutations(
            &regions, 50, 0.01, 0.5, 0.8, 0.15, 0.05,
            &simple_ref_lookup, &mut rng2,
        );

        for (a, b) in v1.iter().zip(v2.iter()) {
            assert_eq!(a.expected_vaf, b.expected_vaf);
            assert_eq!(a.chrom, b.chrom);
        }
    }
}
