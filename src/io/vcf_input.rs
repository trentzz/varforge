//! VCF input parsing: load user-supplied mutation lists from plain or bgzipped VCF.
//!
//! Supports `VAF`/`AF`, `CLONE`, and `CN` INFO fields.  Unknown chromosomes and
//! variants whose REF allele mismatches the reference are handled gracefully
//! (warn + skip / warn + keep respectively).

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use noodles_bgzf as bgzf;
use tracing::warn;

use crate::core::types::{MutationType, Variant};

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Result of parsing a VCF file.
#[derive(Debug, Default)]
pub struct VcfParseResult {
    /// Variants successfully parsed and mapped to `Variant` structs.
    pub variants: Vec<Variant>,
    /// Number of input records that were skipped (unknown chrom, missing ALT, etc.).
    pub skipped: usize,
}

/// Optional per-chromosome reference sequences used for REF allele validation.
///
/// Keys are chromosome names; values are the full sequence bytes.  Pass an
/// empty map to skip REF validation.
pub type RefSequences = HashMap<String, Vec<u8>>;

/// Parse variants from a VCF or bgzipped VCF file.
///
/// * `path` – path to the `.vcf` or `.vcf.gz` file.
/// * `known_chroms` – set of chromosome names present in the reference; variants
///   on chromosomes not in this set are skipped with a warning.  Pass `None` to
///   accept all chromosomes.
/// * `ref_seqs` – optional reference sequences used to validate REF alleles.
pub fn parse_vcf(
    path: &Path,
    known_chroms: Option<&[String]>,
    ref_seqs: Option<&RefSequences>,
) -> Result<VcfParseResult> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");

    let is_bgzipped = matches!(ext, "gz" | "bgz");

    if is_bgzipped {
        let file = std::fs::File::open(path)
            .with_context(|| format!("failed to open VCF: {}", path.display()))?;
        let bgzf_reader = bgzf::Reader::new(file);
        let buf_reader = BufReader::new(bgzf_reader);
        parse_from_reader(buf_reader, known_chroms, ref_seqs)
    } else {
        let file = std::fs::File::open(path)
            .with_context(|| format!("failed to open VCF: {}", path.display()))?;
        let buf_reader = BufReader::new(file);
        parse_from_reader(buf_reader, known_chroms, ref_seqs)
    }
}

/// Parse variants from any `BufRead` source (plain VCF text).
pub fn parse_from_reader<R: BufRead>(
    reader: R,
    known_chroms: Option<&[String]>,
    ref_seqs: Option<&RefSequences>,
) -> Result<VcfParseResult> {
    let mut result = VcfParseResult::default();

    for line in reader.lines() {
        let line = line.context("I/O error reading VCF")?;

        // Skip header lines.
        if line.starts_with('#') {
            continue;
        }
        // Skip blank lines.
        if line.trim().is_empty() {
            continue;
        }

        match parse_record(&line, known_chroms, ref_seqs) {
            Ok(Some(variant)) => result.variants.push(variant),
            Ok(None) => result.skipped += 1,
            Err(e) => {
                warn!("skipping malformed VCF record: {e}");
                result.skipped += 1;
            }
        }
    }

    tracing::debug!(
        loaded = result.variants.len(),
        skipped = result.skipped,
        "VCF parse complete"
    );

    Ok(result)
}

/// Parse variants from a string (convenience wrapper for tests).
// Called only in tests.
#[cfg(test)]
pub fn parse_vcf_str(
    input: &str,
    known_chroms: Option<&[String]>,
    ref_seqs: Option<&RefSequences>,
) -> Result<VcfParseResult> {
    parse_from_reader(input.as_bytes(), known_chroms, ref_seqs)
}

// ---------------------------------------------------------------------------
// Internal record parsing
// ---------------------------------------------------------------------------

/// Parse a single non-header VCF line into a `Variant`.
///
/// Returns `Ok(None)` when the record should be silently skipped (e.g., unknown
/// chromosome).  Logs warnings as appropriate before returning `None`.
fn parse_record(
    line: &str,
    known_chroms: Option<&[String]>,
    ref_seqs: Option<&RefSequences>,
) -> Result<Option<Variant>> {
    let fields: Vec<&str> = line.splitn(9, '\t').collect();
    anyhow::ensure!(fields.len() >= 5, "VCF record has fewer than 5 fields");

    let chrom = fields[0];
    let pos_1based: u64 = fields[1]
        .parse()
        .with_context(|| format!("invalid POS '{}'", fields[1]))?;
    let ref_allele = fields[3];
    let alt_field = fields[4];

    // Skip records with no ALT (monomorphic reference, `.`).
    if alt_field == "." {
        warn!(
            chrom = chrom,
            pos = pos_1based,
            "skipping record with no ALT allele"
        );
        return Ok(None);
    }

    // Take only the first ALT allele; warn if multiple are present.
    let mut alt_iter = alt_field.splitn(2, ',');
    let alt_allele = alt_iter.next().unwrap_or(alt_field);
    if alt_iter.next().is_some() {
        warn!(
            chrom = chrom,
            pos = pos_1based,
            alts = alt_field,
            "multi-ALT record detected; only the first ALT allele will be used"
        );
    }

    // Skip BND / SV notation (task-08 territory).
    if alt_allele.contains('[')
        || alt_allele.contains(']')
        || alt_allele.starts_with('.')
        || alt_allele.ends_with('.')
    {
        warn!(
            chrom = chrom,
            pos = pos_1based,
            alt = alt_allele,
            "skipping SV/BND record"
        );
        return Ok(None);
    }

    // Unknown chromosome check.
    if let Some(chroms) = known_chroms {
        if !chroms.iter().any(|c| c == chrom) {
            warn!(
                chrom = chrom,
                "skipping variant on chromosome not in reference"
            );
            return Ok(None);
        }
    }

    // VCF positions are 1-based; convert to 0-based internally.
    let pos_0based: u64 = pos_1based
        .checked_sub(1)
        .context("VCF POS is 0, which is invalid")?;

    // REF allele validation.
    if let Some(seqs) = ref_seqs {
        if let Some(seq) = seqs.get(chrom) {
            let start = pos_0based as usize;
            let end = start + ref_allele.len();
            if end <= seq.len() {
                let genome_ref = std::str::from_utf8(&seq[start..end])
                    .unwrap_or("")
                    .to_uppercase();
                if genome_ref != ref_allele.to_uppercase() {
                    warn!(
                        chrom = chrom,
                        pos = pos_1based,
                        vcf_ref = ref_allele,
                        genome_ref = genome_ref.as_str(),
                        "REF allele mismatch (possible reference build difference); keeping variant"
                    );
                }
            }
        }
    }

    // Parse INFO field (field index 7 if present).
    let info_str = if fields.len() >= 8 { fields[7] } else { "." };
    let (vaf, clone_id) = parse_info(info_str);

    // Determine mutation type.
    let mutation = classify_mutation(pos_0based, ref_allele, alt_allele)?;

    Ok(Some(Variant {
        chrom: chrom.to_string(),
        mutation,
        expected_vaf: vaf,
        clone_id,
        haplotype: None,
        ccf: None,
    }))
}

/// Extract VAF and CLONE from the INFO field string.
///
/// Recognised keys:
/// - `VAF=<float>` or `AF=<float>` → target VAF (default 0.5)
/// - `CLONE=<str>` → clone assignment (default `"founder"`)
fn parse_info(info: &str) -> (f64, Option<String>) {
    const DEFAULT_VAF: f64 = 0.5;
    const DEFAULT_CLONE: &str = "founder";

    if info == "." || info.is_empty() {
        return (DEFAULT_VAF, Some(DEFAULT_CLONE.to_string()));
    }

    let mut vaf: Option<f64> = None;
    let mut clone: Option<String> = None;

    for field in info.split(';') {
        if let Some(val) = field.strip_prefix("VAF=") {
            vaf = val.parse().ok();
        } else if vaf.is_none() {
            if let Some(val) = field.strip_prefix("AF=") {
                // AF may be comma-separated (multi-allelic); take first.
                let first = val.split(',').next().unwrap_or(val);
                vaf = first.parse().ok();
            }
        }
        if let Some(val) = field.strip_prefix("CLONE=") {
            clone = Some(val.to_string());
        }
    }

    let vaf = vaf.unwrap_or(DEFAULT_VAF);
    let clone = Some(clone.unwrap_or_else(|| DEFAULT_CLONE.to_string()));

    (vaf, clone)
}

/// Classify a variant by comparing REF and ALT lengths.
fn classify_mutation(pos: u64, ref_allele: &str, alt_allele: &str) -> Result<MutationType> {
    let ref_bytes = ref_allele.as_bytes().to_vec();
    let alt_bytes = alt_allele.as_bytes().to_vec();

    let ref_len = ref_bytes.len();
    let alt_len = alt_bytes.len();

    anyhow::ensure!(ref_len > 0, "REF allele is empty");
    anyhow::ensure!(alt_len > 0, "ALT allele is empty");

    let mutation = if ref_len == 1 && alt_len == 1 {
        // SNV
        MutationType::Snv {
            pos,
            ref_base: ref_bytes[0],
            alt_base: alt_bytes[0],
        }
    } else if ref_len == alt_len && ref_len > 1 {
        // MNV
        MutationType::Mnv {
            pos,
            ref_seq: ref_bytes,
            alt_seq: alt_bytes,
        }
    } else {
        // Insertion or deletion (both use Indel)
        MutationType::Indel {
            pos,
            ref_seq: ref_bytes,
            alt_seq: alt_bytes,
        }
    };

    Ok(mutation)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    // -----------------------------------------------------------------------
    // Helper: minimal VCF header
    // -----------------------------------------------------------------------
    const HEADER: &str = "\
##fileformat=VCFv4.3
##INFO=<ID=VAF,Number=1,Type=Float,Description=\"VAF\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">
##INFO=<ID=CLONE,Number=1,Type=String,Description=\"Clone\">
##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

    fn known_chroms() -> Vec<String> {
        vec!["chr1".to_string(), "chr2".to_string(), "chrX".to_string()]
    }

    // -----------------------------------------------------------------------
    // 1. test_parse_snv
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_snv() {
        let vcf = format!("{HEADER}chr1\t1000\t.\tA\tT\t.\tPASS\t.\n");
        let result = parse_vcf_str(&vcf, Some(&known_chroms()), None).unwrap();

        assert_eq!(result.variants.len(), 1);
        assert_eq!(result.skipped, 0);

        let v = &result.variants[0];
        assert_eq!(v.chrom, "chr1");
        assert!(
            matches!(
                &v.mutation,
                MutationType::Snv {
                    pos: 999,
                    ref_base: b'A',
                    alt_base: b'T'
                }
            ),
            "unexpected mutation: {:?}",
            v.mutation
        );
    }

    // -----------------------------------------------------------------------
    // 2. test_parse_indel
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_indel() {
        let vcf = format!(
            "{HEADER}\
chr1\t100\t.\tA\tATG\t.\tPASS\t.\n\
chr1\t200\t.\tATG\tA\t.\tPASS\t.\n"
        );
        let result = parse_vcf_str(&vcf, Some(&known_chroms()), None).unwrap();
        assert_eq!(result.variants.len(), 2);

        // Insertion
        let ins = &result.variants[0];
        assert!(
            matches!(&ins.mutation, MutationType::Indel { pos: 99, .. }),
            "expected insertion indel at pos 99, got {:?}",
            ins.mutation
        );
        if let MutationType::Indel {
            ref_seq, alt_seq, ..
        } = &ins.mutation
        {
            assert_eq!(ref_seq, b"A");
            assert_eq!(alt_seq, b"ATG");
        }

        // Deletion
        let del = &result.variants[1];
        assert!(
            matches!(&del.mutation, MutationType::Indel { pos: 199, .. }),
            "expected deletion indel, got {:?}",
            del.mutation
        );
        if let MutationType::Indel {
            ref_seq, alt_seq, ..
        } = &del.mutation
        {
            assert_eq!(ref_seq, b"ATG");
            assert_eq!(alt_seq, b"A");
        }
    }

    // -----------------------------------------------------------------------
    // 3. test_parse_mnv
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_mnv() {
        let vcf = format!("{HEADER}chr1\t500\t.\tAC\tTG\t.\tPASS\t.\n");
        let result = parse_vcf_str(&vcf, Some(&known_chroms()), None).unwrap();
        assert_eq!(result.variants.len(), 1);

        let v = &result.variants[0];
        assert!(
            matches!(&v.mutation, MutationType::Mnv { pos: 499, .. }),
            "expected MNV, got {:?}",
            v.mutation
        );
        if let MutationType::Mnv {
            ref_seq, alt_seq, ..
        } = &v.mutation
        {
            assert_eq!(ref_seq, b"AC");
            assert_eq!(alt_seq, b"TG");
        }
    }

    // -----------------------------------------------------------------------
    // 4. test_parse_multiple
    // -----------------------------------------------------------------------
    #[test]
    fn test_parse_multiple() {
        let vcf = format!(
            "{HEADER}\
chr1\t100\t.\tA\tC\t.\tPASS\t.\n\
chr1\t200\t.\tG\tT\t.\tPASS\t.\n\
chr2\t50\t.\tA\tATTT\t.\tPASS\t.\n\
chr2\t300\t.\tGCTA\tACTG\t.\tPASS\t.\n"
        );
        let result = parse_vcf_str(&vcf, Some(&known_chroms()), None).unwrap();
        assert_eq!(result.variants.len(), 4);
        assert_eq!(result.skipped, 0);

        // Order preserved.
        assert_eq!(result.variants[0].chrom, "chr1");
        assert_eq!(result.variants[2].chrom, "chr2");
    }

    // -----------------------------------------------------------------------
    // 5. test_custom_vaf
    // -----------------------------------------------------------------------
    #[test]
    fn test_custom_vaf() {
        // VAF from VAF= field.
        let vcf_vaf = format!("{HEADER}chr1\t1000\t.\tA\tT\t.\tPASS\tVAF=0.12\n");
        let result = parse_vcf_str(&vcf_vaf, Some(&known_chroms()), None).unwrap();
        assert!((result.variants[0].expected_vaf - 0.12).abs() < 1e-9);

        // AF= as fallback.
        let vcf_af = format!("{HEADER}chr1\t1000\t.\tA\tT\t.\tPASS\tAF=0.33\n");
        let result = parse_vcf_str(&vcf_af, Some(&known_chroms()), None).unwrap();
        assert!((result.variants[0].expected_vaf - 0.33).abs() < 1e-9);

        // VAF= takes precedence over AF=.
        let vcf_both = format!("{HEADER}chr1\t1000\t.\tA\tT\t.\tPASS\tVAF=0.25;AF=0.99\n");
        let result = parse_vcf_str(&vcf_both, Some(&known_chroms()), None).unwrap();
        assert!((result.variants[0].expected_vaf - 0.25).abs() < 1e-9);
    }

    // -----------------------------------------------------------------------
    // 6. test_clone_assignment
    // -----------------------------------------------------------------------
    #[test]
    fn test_clone_assignment() {
        let vcf = format!("{HEADER}chr1\t1000\t.\tA\tT\t.\tPASS\tCLONE=subclone_B\n");
        let result = parse_vcf_str(&vcf, Some(&known_chroms()), None).unwrap();
        assert_eq!(result.variants[0].clone_id.as_deref(), Some("subclone_B"));
    }

    // -----------------------------------------------------------------------
    // 7. test_missing_info_defaults
    // -----------------------------------------------------------------------
    #[test]
    fn test_missing_info_defaults() {
        // INFO is `.` – defaults should kick in.
        let vcf = format!("{HEADER}chr1\t1000\t.\tA\tT\t.\tPASS\t.\n");
        let result = parse_vcf_str(&vcf, Some(&known_chroms()), None).unwrap();
        let v = &result.variants[0];

        // Default VAF is 0.5.
        assert!(
            (v.expected_vaf - 0.5).abs() < 1e-9,
            "expected default VAF 0.5, got {}",
            v.expected_vaf
        );
        // Default clone_id is "founder".
        assert_eq!(v.clone_id.as_deref(), Some("founder"));
    }

    // -----------------------------------------------------------------------
    // 8. test_ref_mismatch_warning
    // -----------------------------------------------------------------------
    #[test]
    fn test_ref_mismatch_warning() {
        // Reference says position 999 is 'C', but VCF says REF=A.
        // The variant should still be returned (warn, not error).
        let mut ref_seqs: RefSequences = HashMap::new();
        // Build a short sequence: "GCATCGATCG..." with 'C' at index 999.
        let mut seq = vec![b'G'; 2000];
        seq[999] = b'C'; // position 1000 in 1-based = index 999 in 0-based
        ref_seqs.insert("chr1".to_string(), seq);

        let vcf = format!("{HEADER}chr1\t1000\t.\tA\tT\t.\tPASS\t.\n");
        let result = parse_vcf_str(&vcf, Some(&known_chroms()), Some(&ref_seqs)).unwrap();

        // Should still have 1 variant despite the mismatch.
        assert_eq!(
            result.variants.len(),
            1,
            "REF mismatch should warn, not skip"
        );
    }

    // -----------------------------------------------------------------------
    // 9. test_skip_unknown_chrom
    // -----------------------------------------------------------------------
    #[test]
    fn test_skip_unknown_chrom() {
        let vcf = format!(
            "{HEADER}\
chr1\t1000\t.\tA\tT\t.\tPASS\t.\n\
chrUNKNOWN\t500\t.\tG\tC\t.\tPASS\t.\n\
chr2\t200\t.\tA\tG\t.\tPASS\t.\n"
        );
        let result = parse_vcf_str(&vcf, Some(&known_chroms()), None).unwrap();
        assert_eq!(
            result.variants.len(),
            2,
            "only chr1 and chr2 should be kept"
        );
        assert_eq!(result.skipped, 1, "chrUNKNOWN record should be skipped");
    }

    // -----------------------------------------------------------------------
    // 10. test_bgzipped_vcf
    // -----------------------------------------------------------------------
    #[test]
    fn test_bgzipped_vcf() {
        let vcf_content = format!("{HEADER}chr1\t1000\t.\tA\tT\t.\tPASS\tVAF=0.42\n");

        // Write a proper BGZF file using noodles_bgzf::Writer.
        let mut tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        {
            let mut writer = bgzf::Writer::new(&mut tmp);
            writer.write_all(vcf_content.as_bytes()).unwrap();
            writer.try_finish().unwrap();
        }

        let result = parse_vcf(tmp.path(), Some(&known_chroms()), None).unwrap();
        assert_eq!(
            result.variants.len(),
            1,
            "bgzipped VCF should parse correctly"
        );
        let v = &result.variants[0];
        assert!(
            (v.expected_vaf - 0.42).abs() < 1e-9,
            "VAF mismatch: {}",
            v.expected_vaf
        );
        assert!(matches!(&v.mutation, MutationType::Snv { pos: 999, .. }));
    }
}
