use crate::core::types::{MutationType, Region, Variant};
use crate::variants::signatures::{builtin_signature, sbs96_index, SignatureVec};
use rand::Rng;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Generate random mutations distributed across the given regions.
///
/// When `signature` is `Some`, SNV alt base selection is weighted by the SBS96
/// probability vector for the trinucleotide context at each position.
/// For indels and MNVs, alt bases are chosen uniformly regardless of signature.
#[allow(clippy::too_many_arguments)]
pub fn generate_random_mutations<R: Rng>(
    regions: &[Region],
    count: usize,
    vaf_min: f64,
    vaf_max: f64,
    snv_fraction: f64,
    indel_fraction: f64,
    _mnv_fraction: f64,
    reference_lookup: &dyn Fn(&str, u64) -> Option<u8>,
    signature: Option<&SignatureVec>,
    rng: &mut R,
) -> Vec<Variant> {
    let total_length: u64 = regions.iter().map(|r| r.len()).sum();
    if total_length == 0 || count == 0 {
        return Vec::new();
    }

    let mut variants = Vec::with_capacity(count);

    for _ in 0..count {
        // Pick a random position within the regions.
        let rand_offset: u64 = rng.gen_range(0..total_length);
        let (chrom, pos) = offset_to_position(regions, rand_offset);

        // Pick mutation type.
        let type_roll: f64 = rng.gen();
        let mutation = if type_roll < snv_fraction {
            generate_random_snv(&chrom, pos, reference_lookup, signature, rng)
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
                haplotype: None,
                ccf: None,
            });
        }
    }

    variants
}

/// Resolve a named COSMIC SBS signature to its probability vector.
///
/// Returns `None` if the name is not recognised.
pub fn resolve_signature(name: &str) -> Option<&'static SignatureVec> {
    builtin_signature(name)
}

/// Choose an alt base using the SBS96 signature probabilities for the given context.
///
/// Returns the weighted alt base, or falls back to uniform selection if the
/// signature or context is unavailable.
pub fn signature_weighted_alt_base<R: Rng>(
    ref_base: u8,
    ctx_5p: u8,
    ctx_3p: u8,
    signature: &SignatureVec,
    rng: &mut R,
) -> u8 {
    let alts: Vec<u8> = BASES
        .iter()
        .copied()
        .filter(|&b| b != ref_base.to_ascii_uppercase())
        .collect();

    // For each alt, get the SBS96 weight.
    let weights: Vec<f64> = alts
        .iter()
        .map(|&alt| {
            if let Some(idx) = sbs96_index(ctx_5p, ref_base, alt, ctx_3p) {
                signature[idx]
            } else {
                1.0 / 3.0 // uniform fallback
            }
        })
        .collect();

    let total: f64 = weights.iter().sum();
    if total <= 0.0 {
        // All weights zero: fall back to uniform.
        return alts[rng.gen_range(0..alts.len())];
    }

    let r: f64 = rng.gen_range(0.0..total);
    let mut cumulative = 0.0;
    for (i, &w) in weights.iter().enumerate() {
        cumulative += w;
        if r < cumulative {
            return alts[i];
        }
    }
    *alts.last().unwrap()
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
    signature: Option<&SignatureVec>,
    rng: &mut R,
) -> Option<MutationType> {
    let ref_base = reference_lookup(_chrom, pos).unwrap_or(b'N');
    if ref_base == b'N' {
        return None;
    }

    let alt_base = if let Some(sig) = signature {
        // Fetch flanking bases for trinucleotide context.
        let ctx_5p = reference_lookup(_chrom, pos.saturating_sub(1)).unwrap_or(b'N');
        let ctx_3p = reference_lookup(_chrom, pos + 1).unwrap_or(b'N');
        if ctx_5p != b'N' && ctx_3p != b'N' && pos > 0 {
            signature_weighted_alt_base(ref_base, ctx_5p, ctx_3p, sig, rng)
        } else {
            // Context unavailable: fall back to uniform.
            loop {
                let b = BASES[rng.gen_range(0..4)];
                if b != ref_base {
                    break b;
                }
            }
        }
    } else {
        loop {
            let b = BASES[rng.gen_range(0..4)];
            if b != ref_base {
                break b;
            }
        }
    };

    Some(MutationType::Snv {
        pos,
        ref_base,
        alt_base,
    })
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

    Some(MutationType::Mnv {
        pos,
        ref_seq,
        alt_seq,
    })
}

/// Returns true if position `pos` in `seq` falls within a homopolymer run
/// (four or more identical bases) or a dinucleotide repeat (three or more
/// consecutive copies of a 2-mer).
///
/// Used to identify microsatellite-unstable loci for MSI indel generation.
#[allow(dead_code)]
pub fn is_repeat_locus(seq: &[u8], pos: usize) -> bool {
    if seq.is_empty() || pos >= seq.len() {
        return false;
    }

    // Check for a homopolymer run of length ≥ 4 that contains pos.
    let base = seq[pos];
    let run_start = seq[..pos]
        .iter()
        .rposition(|&b| b != base)
        .map(|i| i + 1)
        .unwrap_or(0);
    let run_end = seq[pos..]
        .iter()
        .position(|&b| b != base)
        .map(|i| i + pos)
        .unwrap_or(seq.len());
    if run_end - run_start >= 4 {
        return true;
    }

    // Check for a dinucleotide repeat: the 2-mer starting at pos repeated ≥ 3 times.
    if pos + 1 < seq.len() {
        let di = [seq[pos], seq[pos + 1]];
        let mut repeat_count = 1usize;
        let mut check_pos = pos + 2;
        while check_pos + 1 < seq.len() && seq[check_pos..check_pos + 2] == di {
            repeat_count += 1;
            check_pos += 2;
        }
        if repeat_count >= 3 {
            return true;
        }
    }

    false
}

/// Generate MSI-pattern indels: 1-bp insertions at homopolymer and
/// dinucleotide repeat loci across the given regions.
///
/// Accepts up to `target_count` repeat-locus positions found within
/// `max_attempts` random draws. Each variant is a 1-bp insertion that
/// extends the repeat by one copy of the anchor base, at VAF 0.3.
#[allow(dead_code)]
pub fn generate_msi_indels<R: Rng>(
    regions: &[Region],
    reference_lookup: &dyn Fn(&str, u64) -> Option<u8>,
    rng: &mut R,
    target_count: usize,
    max_attempts: usize,
) -> Vec<Variant> {
    let total_length: u64 = regions.iter().map(|r| r.len()).sum();
    if total_length == 0 || target_count == 0 {
        return Vec::new();
    }

    let mut variants = Vec::new();
    let mut attempts = 0usize;

    while variants.len() < target_count && attempts < max_attempts {
        attempts += 1;

        let rand_offset: u64 = rng.gen_range(0..total_length);
        let (chrom, pos) = offset_to_position(regions, rand_offset);

        // Fetch a small window around the position for repeat context detection.
        let window_start = pos.saturating_sub(5);
        let window_seq: Vec<u8> = (window_start..window_start + 12)
            .filter_map(|p| reference_lookup(&chrom, p))
            .collect();

        let local_pos = (pos - window_start) as usize;
        if local_pos >= window_seq.len() {
            continue;
        }

        if !is_repeat_locus(&window_seq, local_pos) {
            continue;
        }

        let ref_base = match reference_lookup(&chrom, pos) {
            Some(b) if b != b'N' => b,
            _ => continue,
        };

        // 1-bp insertion extending the repeat by one copy of the anchor base.
        let alt_seq = vec![ref_base, ref_base];
        variants.push(Variant {
            chrom,
            mutation: MutationType::Indel {
                pos,
                ref_seq: vec![ref_base],
                alt_seq,
            },
            expected_vaf: 0.3,
            clone_id: Some("msi".to_string()),
            haplotype: None,
            ccf: None,
        });
    }

    variants
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

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
            &regions,
            100,
            0.01,
            0.5,
            0.8,
            0.15,
            0.05,
            &simple_ref_lookup,
            None,
            &mut rng,
        );
        assert_eq!(variants.len(), 100);
    }

    #[test]
    fn test_vaf_range() {
        let mut rng = StdRng::seed_from_u64(42);
        let regions = vec![Region::new("chr1", 0, 10000)];
        let variants = generate_random_mutations(
            &regions,
            1000,
            0.001,
            0.05,
            0.8,
            0.15,
            0.05,
            &simple_ref_lookup,
            None,
            &mut rng,
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
            &regions,
            10000,
            0.01,
            0.5,
            0.8,
            0.15,
            0.05,
            &simple_ref_lookup,
            None,
            &mut rng,
        );

        let snv_count = variants
            .iter()
            .filter(|v| matches!(v.mutation, MutationType::Snv { .. }))
            .count();
        let indel_count = variants
            .iter()
            .filter(|v| matches!(v.mutation, MutationType::Indel { .. }))
            .count();
        let mnv_count = variants
            .iter()
            .filter(|v| matches!(v.mutation, MutationType::Mnv { .. }))
            .count();

        let snv_frac = snv_count as f64 / variants.len() as f64;
        let indel_frac = indel_count as f64 / variants.len() as f64;
        let mnv_frac = mnv_count as f64 / variants.len() as f64;

        assert!(
            (snv_frac - 0.8).abs() < 0.05,
            "SNV fraction {} far from 0.8",
            snv_frac
        );
        assert!(
            (indel_frac - 0.15).abs() < 0.05,
            "indel fraction {} far from 0.15",
            indel_frac
        );
        assert!(
            (mnv_frac - 0.05).abs() < 0.03,
            "MNV fraction {} far from 0.05",
            mnv_frac
        );
    }

    #[test]
    fn test_snv_different_from_ref() {
        let mut rng = StdRng::seed_from_u64(42);
        let regions = vec![Region::new("chr1", 0, 10000)];
        let variants = generate_random_mutations(
            &regions,
            1000,
            0.01,
            0.5,
            1.0,
            0.0,
            0.0,
            &simple_ref_lookup,
            None,
            &mut rng,
        );
        for v in &variants {
            if let MutationType::Snv {
                ref_base, alt_base, ..
            } = v.mutation
            {
                assert_ne!(ref_base, alt_base, "alt must differ from ref");
            }
        }
    }

    #[test]
    fn test_empty_regions() {
        let mut rng = StdRng::seed_from_u64(42);
        let variants = generate_random_mutations(
            &[],
            100,
            0.01,
            0.5,
            0.8,
            0.15,
            0.05,
            &simple_ref_lookup,
            None,
            &mut rng,
        );
        assert!(variants.is_empty());
    }

    #[test]
    fn test_deterministic_with_seed() {
        let regions = vec![Region::new("chr1", 0, 10000)];
        let mut rng1 = StdRng::seed_from_u64(99);
        let mut rng2 = StdRng::seed_from_u64(99);

        let v1 = generate_random_mutations(
            &regions,
            50,
            0.01,
            0.5,
            0.8,
            0.15,
            0.05,
            &simple_ref_lookup,
            None,
            &mut rng1,
        );
        let v2 = generate_random_mutations(
            &regions,
            50,
            0.01,
            0.5,
            0.8,
            0.15,
            0.05,
            &simple_ref_lookup,
            None,
            &mut rng2,
        );

        for (a, b) in v1.iter().zip(v2.iter()) {
            assert_eq!(a.expected_vaf, b.expected_vaf);
            assert_eq!(a.chrom, b.chrom);
        }
    }
}
