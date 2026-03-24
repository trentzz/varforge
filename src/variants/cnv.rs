//! Copy number variant (CNV) modeling: amplifications, deletions, and LOH.
//!
//! Coverage adjustment formula:
//!   adjusted_coverage = base_coverage × (purity × tumor_cn + (1 - purity) × normal_cn) / ploidy
//!
//! For allele-specific copy numbers, major/minor alleles affect the expected
//! VAF of heterozygous variants in the CNV region.

/// A genomic region with assigned copy number state.
#[derive(Debug, Clone, PartialEq)]
pub struct CopyNumberRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub tumor_cn: u32,
    pub normal_cn: u32,
    /// Major allele copy number (for LOH / allele-specific modeling).
    pub major_cn: Option<u32>,
    /// Minor allele copy number (for LOH / allele-specific modeling).
    pub minor_cn: Option<u32>,
}

impl CopyNumberRegion {
    /// Create a new CNV region without allele-specific copy numbers.
    pub fn new(
        chrom: impl Into<String>,
        start: u64,
        end: u64,
        tumor_cn: u32,
        normal_cn: u32,
    ) -> Self {
        Self {
            chrom: chrom.into(),
            start,
            end,
            tumor_cn,
            normal_cn,
            major_cn: None,
            minor_cn: None,
        }
    }

    /// Create a new CNV region with allele-specific copy numbers.
    ///
    /// `major_cn + minor_cn` must equal `tumor_cn`.
    pub fn with_allele_specific(
        chrom: impl Into<String>,
        start: u64,
        end: u64,
        tumor_cn: u32,
        normal_cn: u32,
        major_cn: u32,
        minor_cn: u32,
    ) -> Self {
        debug_assert_eq!(
            major_cn + minor_cn,
            tumor_cn,
            "major_cn + minor_cn must equal tumor_cn"
        );
        Self {
            chrom: chrom.into(),
            start,
            end,
            tumor_cn,
            normal_cn,
            major_cn: Some(major_cn),
            minor_cn: Some(minor_cn),
        }
    }

    /// Return true if the given (chrom, pos) falls within this CNV region.
    pub fn contains(&self, chrom: &str, pos: u64) -> bool {
        self.chrom == chrom && pos >= self.start && pos < self.end
    }
}

/// Adjust base coverage for copy number state and tumour purity.
///
/// adjusted_coverage = base_coverage × (purity × tumor_cn + (1 - purity) × normal_cn) / ploidy
///
/// Returns 0.0 for homozygous tumour deletions with pure tumours (tumor_cn=0, purity=1).
pub fn adjusted_coverage(
    base_coverage: f64,
    purity: f64,
    tumor_cn: u32,
    normal_cn: u32,
    ploidy: u32,
) -> f64 {
    if ploidy == 0 {
        return 0.0;
    }
    let effective = purity * (tumor_cn as f64) + (1.0 - purity) * (normal_cn as f64);
    base_coverage * effective / (ploidy as f64)
}

/// Compute the expected VAF of a heterozygous SNP in a LOH region.
///
/// For a region with allele-specific CN (major_cn, minor_cn):
/// - The fraction of reads from the major allele from tumour cells:
///   `major_fraction = purity × major_cn / (purity × total_cn + (1-purity) × normal_cn)`
/// - For a het SNP on the major allele this is the expected VAF.
///
/// Returns `None` if the region has no allele-specific CN information.
// Called only in tests.
#[allow(dead_code)]
pub fn loh_vaf(purity: f64, major_cn: u32, minor_cn: u32, normal_cn: u32) -> f64 {
    let total_tumor_cn = major_cn + minor_cn;
    let denominator = purity * (total_tumor_cn as f64) + (1.0 - purity) * (normal_cn as f64);
    if denominator == 0.0 {
        return 0.0;
    }
    // A het SNP on the major allele in a pure tumour contributes major_cn copies
    // from tumour cells plus ~0.5 × normal_cn from normal cells (heterozygous normal).
    let numerator = purity * (major_cn as f64) + (1.0 - purity) * 0.5 * (normal_cn as f64);
    numerator / denominator
}

/// Compute the expected VAF of a variant in an amplified region given a
/// specific multiplicity (number of copies of the mutant allele).
///
/// This delegates to the formula from vaf.rs:
///   VAF = CCF × multiplicity × purity / (purity × cn_tumor + (1-purity) × cn_normal)
///
/// With CCF=1.0 (clonal) and provided multiplicity.
// Called only in tests.
#[allow(dead_code)]
pub fn vaf_in_amplified_region(
    purity: f64,
    multiplicity: u32,
    tumor_cn: u32,
    normal_cn: u32,
) -> f64 {
    crate::variants::vaf::expected_vaf(1.0, multiplicity, purity, tumor_cn, normal_cn)
}

/// Find the CNV region that a genomic position falls in, if any.
pub fn find_cn_region<'a>(
    regions: &'a [CopyNumberRegion],
    chrom: &str,
    pos: u64,
) -> Option<&'a CopyNumberRegion> {
    regions.iter().find(|r| r.contains(chrom, pos))
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // 1. Amplification coverage: CN=4 should produce ~2x vs CN=2
    // -----------------------------------------------------------------------
    #[test]
    fn test_amplification_coverage() {
        let base = 30.0;
        let purity = 1.0;
        let ploidy = 2;

        let diploid = adjusted_coverage(base, purity, 2, 2, ploidy);
        let amplified = adjusted_coverage(base, purity, 4, 2, ploidy);

        // diploid: 30 × (1.0 × 2) / 2 = 30
        assert!(
            (diploid - 30.0).abs() < 1e-10,
            "diploid coverage should be 30, got {diploid}"
        );
        // amplified: 30 × (1.0 × 4) / 2 = 60 → ~2x diploid
        assert!(
            (amplified - 60.0).abs() < 1e-10,
            "CN=4 coverage should be 60, got {amplified}"
        );
        assert!(
            (amplified / diploid - 2.0).abs() < 1e-10,
            "CN=4 should be ~2x CN=2, ratio was {}",
            amplified / diploid
        );
    }

    // -----------------------------------------------------------------------
    // 2. Deletion coverage: CN=1 should produce ~0.5x vs CN=2
    // -----------------------------------------------------------------------
    #[test]
    fn test_deletion_coverage() {
        let base = 30.0;
        let purity = 1.0;
        let ploidy = 2;

        let diploid = adjusted_coverage(base, purity, 2, 2, ploidy);
        let deleted = adjusted_coverage(base, purity, 1, 2, ploidy);

        // deleted: 30 × (1.0 × 1) / 2 = 15
        assert!(
            (deleted - 15.0).abs() < 1e-10,
            "CN=1 coverage should be 15, got {deleted}"
        );
        assert!(
            (deleted / diploid - 0.5).abs() < 1e-10,
            "CN=1 should be ~0.5x CN=2, ratio was {}",
            deleted / diploid
        );
    }

    // -----------------------------------------------------------------------
    // 3. Homozygous deletion: CN=0 with purity=1 → only normal contamination reads
    // -----------------------------------------------------------------------
    #[test]
    fn test_homodel_coverage() {
        let base = 30.0;
        let ploidy = 2;

        // Pure tumour: purity=1.0, tumor_cn=0 → no reads from tumour or normal
        let pure_homdel = adjusted_coverage(base, 1.0, 0, 2, ploidy);
        assert!(
            pure_homdel < 1e-10,
            "pure homdel should produce ~0 coverage, got {pure_homdel}"
        );

        // Impure tumour: purity=0.7, tumor_cn=0, normal_cn=2
        // coverage = 30 × (0.7 × 0 + 0.3 × 2) / 2 = 30 × 0.6 / 2 = 9.0
        let impure_homdel = adjusted_coverage(base, 0.7, 0, 2, ploidy);
        assert!(
            (impure_homdel - 9.0).abs() < 1e-10,
            "impure homdel should have only normal-contamination coverage (~9.0), got {impure_homdel}"
        );
    }

    // -----------------------------------------------------------------------
    // 4. LOH allele frequency: het SNP shows VAF ~1.0 for major allele (pure LOH)
    // -----------------------------------------------------------------------
    #[test]
    fn test_loh_allele_frequency() {
        // LOH: CN=2, major=2, minor=0 → all tumour reads from one allele
        // purity=1.0, normal_cn=2
        // VAF = (1.0 × 2 + 0.0 × 0.5 × 2) / (1.0 × 2) = 2/2 = 1.0
        let vaf = loh_vaf(1.0, 2, 0, 2);
        assert!(
            (vaf - 1.0).abs() < 1e-10,
            "pure LOH major allele VAF should be ~1.0, got {vaf}"
        );

        // With purity=0.5: het SNP on minor allele (minor_cn=0) should be ~0.0
        // The normal cells contribute 0.5 × 0.5 × 2 = 0.5 reads to the alt allele
        // from their het copy, out of total = 0.5 × 2 + 0.5 × 2 = 2.0
        // But we model the major allele here, which is cn=2 in tumour.
        // VAF(major) = (0.5×2 + 0.5×0.5×2) / (0.5×2 + 0.5×2)
        //            = (1.0 + 0.5) / 2.0 = 0.75
        let vaf_impure = loh_vaf(0.5, 2, 0, 2);
        assert!(
            (vaf_impure - 0.75).abs() < 1e-10,
            "impure LOH major allele VAF should be 0.75, got {vaf_impure}"
        );
    }

    // -----------------------------------------------------------------------
    // 5. Coverage correctly adjusted for tumour purity
    // -----------------------------------------------------------------------
    #[test]
    fn test_cn_with_purity() {
        let base = 100.0;
        let ploidy = 2;

        // purity=0.7, tumor_cn=4, normal_cn=2
        // coverage = 100 × (0.7×4 + 0.3×2) / 2 = 100 × (2.8 + 0.6) / 2 = 100 × 1.7 = 170
        let cov = adjusted_coverage(base, 0.7, 4, 2, ploidy);
        assert!(
            (cov - 170.0).abs() < 1e-10,
            "coverage with purity=0.7, CN=4 should be 170.0, got {cov}"
        );

        // purity=0.5, tumor_cn=2, normal_cn=2 → same as normal: 100 × (0.5×2+0.5×2)/2 = 100
        let cov_balanced = adjusted_coverage(base, 0.5, 2, 2, ploidy);
        assert!(
            (cov_balanced - 100.0).abs() < 1e-10,
            "balanced CN with purity=0.5 should give base coverage, got {cov_balanced}"
        );
    }

    // -----------------------------------------------------------------------
    // 6. VAF shifts correctly for amplified variants
    // -----------------------------------------------------------------------
    #[test]
    fn test_vaf_in_amplified_region() {
        // CN=4, purity=1.0, multiplicity=1 (one copy of alt allele)
        // VAF = 1 × 1 × 1 / (1 × 4) = 0.25
        let vaf_1 = vaf_in_amplified_region(1.0, 1, 4, 2);
        assert!(
            (vaf_1 - 0.25).abs() < 1e-10,
            "mult=1 in CN=4 should give VAF=0.25, got {vaf_1}"
        );

        // CN=4, purity=1.0, multiplicity=3 (three copies of alt)
        // VAF = 1 × 3 × 1 / (1 × 4) = 0.75
        let vaf_3 = vaf_in_amplified_region(1.0, 3, 4, 2);
        assert!(
            (vaf_3 - 0.75).abs() < 1e-10,
            "mult=3 in CN=4 should give VAF=0.75, got {vaf_3}"
        );

        // Unbalanced gain CN=3, major=2, minor=1, purity=1.0, mult=2
        // VAF = 1 × 2 × 1 / (1 × 3) = 2/3
        let vaf_unbalanced = vaf_in_amplified_region(1.0, 2, 3, 2);
        assert!(
            (vaf_unbalanced - 2.0 / 3.0).abs() < 1e-10,
            "mult=2 in CN=3 should give VAF=2/3, got {vaf_unbalanced}"
        );
    }

    // -----------------------------------------------------------------------
    // 7. YAML config with copy_number section parses correctly
    // -----------------------------------------------------------------------
    #[test]
    fn test_cn_config_parsing() {
        use crate::io::config::load;
        use std::io::Write;
        use tempfile::NamedTempFile;

        let yaml = r#"
reference: /tmp/ref.fa
output:
  directory: /tmp/out
copy_number:
  - region: "chr7:55000000-55200000"
    tumor_cn: 4
    normal_cn: 2
  - region: "chr17:7500000-7700000"
    tumor_cn: 1
    normal_cn: 2
  - region: "chr13:32300000-32400000"
    tumor_cn: 0
    normal_cn: 2
"#;

        let mut f = NamedTempFile::new().unwrap();
        f.write_all(yaml.as_bytes()).unwrap();
        let cfg = load(f.path()).unwrap();

        let cn = cfg
            .copy_number
            .as_ref()
            .expect("copy_number should be present");
        assert_eq!(cn.len(), 3);

        assert_eq!(cn[0].region, "chr7:55000000-55200000");
        assert_eq!(cn[0].tumor_cn, 4);
        assert_eq!(cn[0].normal_cn, 2);

        assert_eq!(cn[1].region, "chr17:7500000-7700000");
        assert_eq!(cn[1].tumor_cn, 1);
        assert_eq!(cn[1].normal_cn, 2);

        assert_eq!(cn[2].region, "chr13:32300000-32400000");
        assert_eq!(cn[2].tumor_cn, 0);
        assert_eq!(cn[2].normal_cn, 2);
    }
}
