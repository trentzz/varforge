//! Germline SNP and indel generation at population-frequency densities.
//!
//! Variants are distributed across supplied regions at configurable densities
//! (SNPs per kbp, indels per kbp). Heterozygous variants receive VAF 0.5;
//! homozygous variants receive VAF 1.0. Clone IDs are set to `"germline_het"`
//! or `"germline_hom"` so downstream code can filter them from somatic variants.

use rand::Rng;

use crate::core::types::{MutationType, Region, Variant};
use crate::io::config::GermlineConfig;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Generate germline variants distributed across the given regions.
///
/// Returns a `Vec<Variant>` with VAF 0.5 (heterozygous) or 1.0 (homozygous).
/// Positions falling on an `N` base are skipped.
pub fn generate_germline_variants<R: Rng>(
    regions: &[Region],
    cfg: &GermlineConfig,
    reference_lookup: &dyn Fn(&str, u64) -> Option<u8>,
    rng: &mut R,
) -> Vec<Variant> {
    let mut variants = Vec::new();

    for region in regions {
        let region_kbp = region.len() as f64 / 1000.0;

        // Heterozygous SNPs.
        let het_count = (region_kbp * cfg.het_snp_density).round() as usize;
        for _ in 0..het_count {
            let pos = region.start + rng.gen_range(0..region.len());
            if let Some(ref_base) = reference_lookup(&region.chrom, pos) {
                if ref_base == b'N' {
                    continue;
                }
                let alts: Vec<u8> = BASES.iter().copied().filter(|&b| b != ref_base).collect();
                let alt = alts[rng.gen_range(0..alts.len())];
                variants.push(Variant {
                    chrom: region.chrom.clone(),
                    mutation: MutationType::Snv {
                        pos,
                        ref_base,
                        alt_base: alt,
                    },
                    expected_vaf: 0.5,
                    clone_id: Some("germline_het".to_string()),
                    haplotype: None,
                });
            }
        }

        // Homozygous SNPs.
        let hom_count = (region_kbp * cfg.hom_snp_density).round() as usize;
        for _ in 0..hom_count {
            let pos = region.start + rng.gen_range(0..region.len());
            if let Some(ref_base) = reference_lookup(&region.chrom, pos) {
                if ref_base == b'N' {
                    continue;
                }
                let alts: Vec<u8> = BASES.iter().copied().filter(|&b| b != ref_base).collect();
                let alt = alts[rng.gen_range(0..alts.len())];
                variants.push(Variant {
                    chrom: region.chrom.clone(),
                    mutation: MutationType::Snv {
                        pos,
                        ref_base,
                        alt_base: alt,
                    },
                    expected_vaf: 1.0,
                    clone_id: Some("germline_hom".to_string()),
                    haplotype: None,
                });
            }
        }

        // Heterozygous indels.
        let indel_count = (region_kbp * cfg.het_indel_density).round() as usize;
        for _ in 0..indel_count {
            let pos = region.start + rng.gen_range(0..region.len().saturating_sub(5));
            if rng.gen_bool(0.5) {
                // Insertion: anchor base + 1–3 random bases.
                let ins_len = rng.gen_range(1usize..=3);
                let inserted: Vec<u8> = (0..ins_len).map(|_| BASES[rng.gen_range(0..4)]).collect();
                if let Some(ref_base) = reference_lookup(&region.chrom, pos) {
                    let mut alt = vec![ref_base];
                    alt.extend_from_slice(&inserted);
                    variants.push(Variant {
                        chrom: region.chrom.clone(),
                        mutation: MutationType::Indel {
                            pos,
                            ref_seq: vec![ref_base],
                            alt_seq: alt,
                        },
                        expected_vaf: 0.5,
                        clone_id: Some("germline_het".to_string()),
                        haplotype: None,
                    });
                }
            } else {
                // Deletion: 1–3 bp.
                let del_len = rng.gen_range(1usize..=3);
                let mut ref_seq = Vec::new();
                for i in 0..=del_len {
                    if let Some(b) = reference_lookup(&region.chrom, pos + i as u64) {
                        ref_seq.push(b);
                    }
                }
                if ref_seq.len() == del_len + 1 {
                    let alt_seq = vec![ref_seq[0]];
                    variants.push(Variant {
                        chrom: region.chrom.clone(),
                        mutation: MutationType::Indel {
                            pos,
                            ref_seq,
                            alt_seq,
                        },
                        expected_vaf: 0.5,
                        clone_id: Some("germline_het".to_string()),
                        haplotype: None,
                    });
                }
            }
        }
    }

    variants
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn simple_lookup(_chrom: &str, pos: u64) -> Option<u8> {
        Some(BASES[(pos % 4) as usize])
    }

    #[test]
    fn test_germline_variant_counts() {
        let mut rng = StdRng::seed_from_u64(42);
        // 10 kbp region → expect ~6 het SNPs, ~3 hom SNPs, ~0–1 indels.
        let regions = vec![Region::new("chr1", 0, 10_000)];
        let cfg = GermlineConfig::default();
        let variants = generate_germline_variants(&regions, &cfg, &simple_lookup, &mut rng);
        assert!(
            !variants.is_empty(),
            "expected germline variants for a 10 kbp region"
        );
    }

    #[test]
    fn test_het_snps_have_vaf_half() {
        let mut rng = StdRng::seed_from_u64(1);
        let regions = vec![Region::new("chr1", 0, 100_000)];
        let cfg = GermlineConfig::default();
        let variants = generate_germline_variants(&regions, &cfg, &simple_lookup, &mut rng);
        for v in variants
            .iter()
            .filter(|v| v.clone_id.as_deref() == Some("germline_het"))
        {
            assert!(
                (v.expected_vaf - 0.5).abs() < f64::EPSILON,
                "het germline variant must have VAF 0.5, got {}",
                v.expected_vaf
            );
        }
    }

    #[test]
    fn test_hom_snps_have_vaf_one() {
        let mut rng = StdRng::seed_from_u64(2);
        let regions = vec![Region::new("chr1", 0, 100_000)];
        let cfg = GermlineConfig::default();
        let variants = generate_germline_variants(&regions, &cfg, &simple_lookup, &mut rng);
        for v in variants
            .iter()
            .filter(|v| v.clone_id.as_deref() == Some("germline_hom"))
        {
            assert!(
                (v.expected_vaf - 1.0).abs() < f64::EPSILON,
                "hom germline variant must have VAF 1.0, got {}",
                v.expected_vaf
            );
        }
    }

    #[test]
    fn test_snv_alt_differs_from_ref() {
        let mut rng = StdRng::seed_from_u64(7);
        let regions = vec![Region::new("chr1", 0, 100_000)];
        let cfg = GermlineConfig::default();
        let variants = generate_germline_variants(&regions, &cfg, &simple_lookup, &mut rng);
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
    fn test_empty_regions_produces_no_variants() {
        let mut rng = StdRng::seed_from_u64(99);
        let cfg = GermlineConfig::default();
        let variants = generate_germline_variants(&[], &cfg, &simple_lookup, &mut rng);
        assert!(variants.is_empty());
    }

    #[test]
    fn test_deterministic_with_seed() {
        let regions = vec![Region::new("chr1", 0, 50_000)];
        let cfg = GermlineConfig::default();
        let mut rng1 = StdRng::seed_from_u64(55);
        let mut rng2 = StdRng::seed_from_u64(55);
        let v1 = generate_germline_variants(&regions, &cfg, &simple_lookup, &mut rng1);
        let v2 = generate_germline_variants(&regions, &cfg, &simple_lookup, &mut rng2);
        assert_eq!(
            v1.len(),
            v2.len(),
            "same seed must produce the same number of variants"
        );
        for (a, b) in v1.iter().zip(v2.iter()) {
            assert_eq!(a.chrom, b.chrom);
            assert_eq!(a.expected_vaf, b.expected_vaf);
        }
    }
}
