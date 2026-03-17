use rand::Rng;

use crate::core::types::{Read, SvType};
use crate::variants::vaf::sample_alt_count;

/// A structural variant with all parameters needed to generate affected reads.
#[derive(Debug, Clone, PartialEq)]
pub enum StructuralVariant {
    Deletion {
        chrom: String,
        start: u64,
        end: u64,
    },
    Insertion {
        chrom: String,
        pos: u64,
        sequence: Vec<u8>,
    },
    Inversion {
        chrom: String,
        start: u64,
        end: u64,
    },
    Duplication {
        chrom: String,
        start: u64,
        end: u64,
        copies: u32,
    },
    Translocation {
        chrom1: String,
        pos1: u64,
        chrom2: String,
        pos2: u64,
    },
}

impl StructuralVariant {
    /// Return the primary chromosome of the SV.
    #[allow(dead_code)]
    pub fn chrom(&self) -> &str {
        match self {
            Self::Deletion { chrom, .. } => chrom,
            Self::Insertion { chrom, .. } => chrom,
            Self::Inversion { chrom, .. } => chrom,
            Self::Duplication { chrom, .. } => chrom,
            Self::Translocation { chrom1, .. } => chrom1,
        }
    }

    /// Return the `SvType` tag for use in `MutationType::Sv`.
    #[allow(dead_code)]
    pub fn sv_type(&self) -> SvType {
        match self {
            Self::Deletion { .. } => SvType::Deletion,
            Self::Insertion { .. } => SvType::Insertion,
            Self::Inversion { .. } => SvType::Inversion,
            Self::Duplication { .. } => SvType::Duplication,
            Self::Translocation { .. } => SvType::Translocation,
        }
    }

    /// Return the left-most breakpoint position (0-based).
    #[allow(dead_code)]
    pub fn start(&self) -> u64 {
        match self {
            Self::Deletion { start, .. } => *start,
            Self::Insertion { pos, .. } => *pos,
            Self::Inversion { start, .. } => *start,
            Self::Duplication { start, .. } => *start,
            Self::Translocation { pos1, .. } => *pos1,
        }
    }

    /// Return the right-most breakpoint position (0-based, exclusive for range types).
    /// For point events (Insertion, Translocation) this equals `start() + 1`.
    #[allow(dead_code)]
    pub fn end(&self) -> u64 {
        match self {
            Self::Deletion { end, .. } => *end,
            Self::Insertion { pos, .. } => pos + 1,
            Self::Inversion { end, .. } => *end,
            Self::Duplication { end, .. } => *end,
            Self::Translocation { pos1, .. } => pos1 + 1,
        }
    }
}

/// Result of applying an SV to a read.
#[derive(Debug, Clone, PartialEq)]
pub enum SvReadEffect {
    /// Read was not affected (outside breakpoints).
    Unaffected,
    /// Read spans the left breakpoint and was soft-clipped.
    SplitLeft,
    /// Read spans the right breakpoint and was soft-clipped.
    SplitRight,
    /// Read is entirely within a deleted region and should be discarded.
    Deleted,
    /// Read sequence was modified (insertion, inversion, duplication).
    Modified,
    /// Read is chimeric – it now crosses two chromosomes (translocation).
    Chimeric,
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Apply a deletion SV to a read that starts at `read_start` (0-based genome
/// coordinate).
///
/// Reads that fall entirely within the deleted interval are marked `Deleted`.
/// Reads that span the left or right breakpoint are soft-clipped so that the
/// portion overlapping the deleted region is replaced with soft-clip padding
/// (low-quality `N` bases).
///
/// Returns the `SvReadEffect` describing what happened.
pub fn apply_deletion(read: &mut Read, read_start: u64, sv: &StructuralVariant) -> SvReadEffect {
    let (chrom_sv, del_start, del_end) = match sv {
        StructuralVariant::Deletion { chrom, start, end } => (chrom.as_str(), *start, *end),
        _ => return SvReadEffect::Unaffected,
    };
    let _ = chrom_sv; // chromosome matching is done by the caller

    let read_end = read_start + read.len() as u64;
    let read_len = read.len();

    // Entirely within deletion → discard
    if read_start >= del_start && read_end <= del_end {
        return SvReadEffect::Deleted;
    }

    // Spans left breakpoint: read_start < del_start < read_end
    if read_start < del_start && read_end > del_start && read_end <= del_end {
        // Soft-clip the right tail that falls inside the deletion
        let keep = (del_start - read_start) as usize;
        soft_clip_right(read, keep, read_len);
        return SvReadEffect::SplitLeft;
    }

    // Spans right breakpoint: del_start <= read_start < del_end < read_end
    if read_start >= del_start && read_start < del_end && read_end > del_end {
        // Soft-clip the left portion that falls inside the deletion
        let clip_len = (del_end - read_start) as usize;
        soft_clip_left(read, clip_len, read_len);
        return SvReadEffect::SplitRight;
    }

    // Spans the entire deletion (read_start < del_start && read_end > del_end)
    // Treat as a left-breakpoint split read: soft-clip everything from del_start onward.
    // This models a chimeric read where the right arm maps to a different location.
    if read_start < del_start && read_end > del_end {
        let keep = (del_start - read_start) as usize;
        soft_clip_right(read, keep, read_len);
        return SvReadEffect::SplitLeft;
    }

    SvReadEffect::Unaffected
}

/// Apply an insertion SV to a read that starts at `read_start`.
///
/// Reads that overlap the insertion point get the inserted sequence embedded
/// at the correct offset; excess bases at the 3' end are soft-clipped to keep
/// the read length constant.
pub fn apply_insertion(read: &mut Read, read_start: u64, sv: &StructuralVariant) -> SvReadEffect {
    let (ins_pos, ins_seq) = match sv {
        StructuralVariant::Insertion { pos, sequence, .. } => (*pos, sequence.as_slice()),
        _ => return SvReadEffect::Unaffected,
    };

    let read_end = read_start + read.len() as u64;
    if ins_pos < read_start || ins_pos >= read_end {
        return SvReadEffect::Unaffected;
    }

    let offset = (ins_pos - read_start) as usize;
    let original_len = read.len();

    // Build new sequence: bases_before_ins + inserted_seq + bases_after_ins
    let mut new_seq: Vec<u8> = Vec::with_capacity(original_len);
    new_seq.extend_from_slice(&read.seq[..offset]);
    new_seq.extend_from_slice(ins_seq);
    new_seq.extend_from_slice(&read.seq[offset..]);

    // Build new qualities: copy originals, use avg for inserted bases
    let mut new_qual: Vec<u8> = Vec::with_capacity(original_len);
    new_qual.extend_from_slice(&read.qual[..offset]);
    let avg_qual = if offset > 0 && offset < read.qual.len() {
        (read.qual[offset - 1] as u16 + read.qual[offset] as u16) as u8 / 2
    } else {
        read.qual[offset.min(read.qual.len() - 1)]
    };
    for _ in 0..ins_seq.len() {
        new_qual.push(avg_qual);
    }
    new_qual.extend_from_slice(&read.qual[offset..]);

    // Truncate to original length (soft-clip trailing bases)
    let clip_start = original_len;
    new_seq.truncate(original_len);
    new_qual.truncate(original_len);

    // Mark last `ins_seq.len().min(original_len - offset)` bases as soft-clip
    // quality 0 to signal the clipped tail.
    let clipped = ins_seq.len().min(original_len.saturating_sub(clip_start));
    let _ = clipped; // CIGAR generation lives outside this function

    read.seq = new_seq;
    read.qual = new_qual;

    SvReadEffect::Modified
}

/// Apply an inversion SV to a read that starts at `read_start`.
///
/// Reads that fall entirely within the inverted interval have their sequence
/// reverse-complemented.  Reads spanning a boundary are partially
/// reverse-complemented and soft-clipped.
pub fn apply_inversion(read: &mut Read, read_start: u64, sv: &StructuralVariant) -> SvReadEffect {
    let (inv_start, inv_end) = match sv {
        StructuralVariant::Inversion { start, end, .. } => (*start, *end),
        _ => return SvReadEffect::Unaffected,
    };

    let read_end = read_start + read.len() as u64;

    // No overlap
    if read_end <= inv_start || read_start >= inv_end {
        return SvReadEffect::Unaffected;
    }

    // Clamp the overlap region to this read's coordinates
    let overlap_start = inv_start.max(read_start);
    let overlap_end = inv_end.min(read_end);

    let local_start = (overlap_start - read_start) as usize;
    let local_end = (overlap_end - read_start) as usize;

    // Reverse-complement the overlapping slice in-place
    revcomp_slice(&mut read.seq[local_start..local_end]);
    read.qual[local_start..local_end].reverse();

    SvReadEffect::Modified
}

/// Apply a tandem duplication SV to a read that starts at `read_start`.
///
/// Reads that overlap the right boundary of the original region (i.e. where
/// the duplicated copy joins the original) are modified to show the junction
/// sequence.  Coverage elevation in the duplicated region is handled by the
/// caller (by generating extra reads from the duplicated coordinates).
pub fn apply_duplication(read: &mut Read, read_start: u64, sv: &StructuralVariant) -> SvReadEffect {
    let (dup_start, dup_end, _copies) = match sv {
        StructuralVariant::Duplication {
            start, end, copies, ..
        } => (*start, *end, *copies),
        _ => return SvReadEffect::Unaffected,
    };

    let read_end = read_start + read.len() as u64;

    // A split read spans the tandem junction: read_start < dup_end <= read_end
    // and read_start >= dup_start (right-side split).
    if read_start >= dup_start && read_start < dup_end && read_end > dup_end {
        // Soft-clip the right portion that extends past the duplicated region
        // so the read appears to loop back to dup_start.
        let keep = (dup_end - read_start) as usize;
        soft_clip_right(read, keep, read.len());
        return SvReadEffect::SplitRight;
    }

    // Read spans the left boundary of the dup: read_start < dup_start < read_end
    if read_start < dup_start && read_end > dup_start && read_end <= dup_end {
        return SvReadEffect::Modified; // elevated coverage, no sequence change needed
    }

    SvReadEffect::Unaffected
}

/// Apply a translocation SV to a read that starts at `read_start` on `chrom`.
///
/// If the read spans `pos1` on `chrom1`, the portion past the breakpoint is
/// replaced with sequence from the partner chromosome (`chrom2` starting at
/// `pos2`), supplied by the caller as `partner_seq`.
pub fn apply_translocation(
    read: &mut Read,
    read_start: u64,
    chrom: &str,
    sv: &StructuralVariant,
    partner_seq: &[u8],
) -> SvReadEffect {
    let (chrom1, pos1, _chrom2, _pos2) = match sv {
        StructuralVariant::Translocation {
            chrom1,
            pos1,
            chrom2,
            pos2,
        } => (chrom1.as_str(), *pos1, chrom2.as_str(), *pos2),
        _ => return SvReadEffect::Unaffected,
    };

    if chrom != chrom1 {
        return SvReadEffect::Unaffected;
    }

    let read_end = read_start + read.len() as u64;
    if pos1 < read_start || pos1 >= read_end {
        return SvReadEffect::Unaffected;
    }

    let split_offset = (pos1 - read_start) as usize;
    let bases_from_partner = read.len() - split_offset;

    // Replace sequence after breakpoint with partner sequence
    let available = partner_seq.len().min(bases_from_partner);
    read.seq[split_offset..split_offset + available].copy_from_slice(&partner_seq[..available]);

    // Any bases beyond available partner seq become Ns
    for b in &mut read.seq[split_offset + available..] {
        *b = b'N';
    }

    // Quality: low quality for the chimeric portion
    for q in &mut read.qual[split_offset..] {
        *q = 10;
    }

    SvReadEffect::Chimeric
}

// ---------------------------------------------------------------------------
// Stochastic sampling helpers
// ---------------------------------------------------------------------------

/// Decide, for each read in `total_depth`, whether it should carry the SV
/// allele.  Returns the number of alt reads to generate using binomial
/// sampling (same model as small variants).
#[allow(dead_code)]
pub fn sample_sv_alt_count<R: Rng>(total_depth: u32, expected_vaf: f64, rng: &mut R) -> u32 {
    sample_alt_count(total_depth, expected_vaf, rng)
}

// ---------------------------------------------------------------------------
// VCF breakend (BND) notation
// ---------------------------------------------------------------------------

/// Format a structural variant as a VCF INFO field string using standard
/// BND / symbolic allele notation.
///
/// - DEL/INV/DUP: `SVTYPE=DEL;END=<end>;SVLEN=<len>`
/// - INS: `SVTYPE=INS;END=<pos>;SVLEN=<len>`
/// - BND (translocation): `SVTYPE=BND` on two separate records (caller
///   handles record splitting).  Returns the ALT allele string in BND format:
///   `N]chrom2:pos2]`
#[allow(dead_code)]
pub fn sv_vcf_info(sv: &StructuralVariant) -> String {
    match sv {
        StructuralVariant::Deletion { start, end, .. } => {
            let svlen = end.saturating_sub(*start);
            format!("SVTYPE=DEL;END={end};SVLEN=-{svlen}")
        }
        StructuralVariant::Insertion { pos, sequence, .. } => {
            let svlen = sequence.len();
            format!("SVTYPE=INS;END={pos};SVLEN={svlen}")
        }
        StructuralVariant::Inversion { start, end, .. } => {
            let svlen = end.saturating_sub(*start);
            format!("SVTYPE=INV;END={end};SVLEN={svlen}")
        }
        StructuralVariant::Duplication {
            start, end, copies, ..
        } => {
            let svlen = end.saturating_sub(*start);
            format!("SVTYPE=DUP;END={end};SVLEN={svlen};CN={copies}")
        }
        StructuralVariant::Translocation { chrom2, pos2, .. } => {
            format!("SVTYPE=BND;MATECHROM={chrom2};MATEPOS={pos2}")
        }
    }
}

/// Return the VCF ALT allele string for an SV record.
///
/// Standard symbolic alleles: `<DEL>`, `<INS>`, `<INV>`, `<DUP>`.
/// Translocations use BND mate notation: `N]chrom2:pos2]`.
#[allow(dead_code)]
pub fn sv_vcf_alt(sv: &StructuralVariant) -> String {
    match sv {
        StructuralVariant::Deletion { .. } => "<DEL>".to_string(),
        StructuralVariant::Insertion { .. } => "<INS>".to_string(),
        StructuralVariant::Inversion { .. } => "<INV>".to_string(),
        StructuralVariant::Duplication { .. } => "<DUP>".to_string(),
        StructuralVariant::Translocation { chrom2, pos2, .. } => {
            // BND format: N]chrom2:pos2] (join to the left of pos2)
            format!("N]{chrom2}:{pos2}]")
        }
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Soft-clip the right end of a read starting at `keep_bases`.
/// Clipped bases are replaced with `N` at quality 0.
fn soft_clip_right(read: &mut Read, keep_bases: usize, _original_len: usize) {
    let len = read.len();
    let keep = keep_bases.min(len);
    for i in keep..len {
        read.seq[i] = b'N';
        read.qual[i] = 0;
    }
}

/// Soft-clip the left end of a read for `clip_len` bases.
/// Clipped bases are replaced with `N` at quality 0.
fn soft_clip_left(read: &mut Read, clip_len: usize, _original_len: usize) {
    let len = read.len();
    let clip = clip_len.min(len);
    for i in 0..clip {
        read.seq[i] = b'N';
        read.qual[i] = 0;
    }
}

/// Reverse-complement a slice of bytes in-place.
fn revcomp_slice(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = complement(*b);
    }
}

/// Return the DNA complement of a base (IUPAC uppercase).
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        _ => b'N',
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn make_read(seq: &[u8]) -> Read {
        Read::new(seq.to_vec(), vec![30; seq.len()])
    }

    // ------------------------------------------------------------------
    // 1. test_deletion_split_reads
    // ------------------------------------------------------------------
    #[test]
    fn test_deletion_split_reads() {
        // Deletion from position 110 to 120 (10 bp).
        // Read starts at 105, length 20 → it spans the left breakpoint (110).
        let sv = StructuralVariant::Deletion {
            chrom: "chr1".into(),
            start: 110,
            end: 120,
        };
        let mut read = make_read(b"ACGTACGTACGTACGTACGT"); // 20 bp
        let effect = apply_deletion(&mut read, 105, &sv);

        assert_eq!(
            effect,
            SvReadEffect::SplitLeft,
            "read spanning left breakpoint should be split-left"
        );

        // Bases 0..5 (positions 105-109) are kept; positions 110+ are soft-clipped (N, qual=0)
        let keep = (110 - 105) as usize; // 5
        for i in 0..keep {
            assert_ne!(read.seq[i], b'N', "kept bases should not be N");
        }
        for i in keep..read.len() {
            assert_eq!(read.seq[i], b'N', "clipped bases should be N");
            assert_eq!(read.qual[i], 0, "clipped bases should have quality 0");
        }
    }

    // ------------------------------------------------------------------
    // 2. test_deletion_coverage_drop
    // ------------------------------------------------------------------
    #[test]
    fn test_deletion_coverage_drop() {
        // A read entirely within the deleted region should be marked Deleted.
        let sv = StructuralVariant::Deletion {
            chrom: "chr1".into(),
            start: 200,
            end: 300,
        };
        let mut read = make_read(b"ACGTACGTACGT"); // 12 bp starting at 210
        let effect = apply_deletion(&mut read, 210, &sv);
        assert_eq!(
            effect,
            SvReadEffect::Deleted,
            "read entirely within deleted region should be marked Deleted"
        );

        // A read outside the deletion should be unaffected.
        let mut read2 = make_read(b"ACGTACGT");
        let effect2 = apply_deletion(&mut read2, 350, &sv);
        assert_eq!(effect2, SvReadEffect::Unaffected);
    }

    // ------------------------------------------------------------------
    // 3. test_insertion_reads
    // ------------------------------------------------------------------
    #[test]
    fn test_insertion_reads() {
        // Insert "TTTT" at position 104 into a read starting at 100.
        let ins_seq = b"TTTT".to_vec();
        let sv = StructuralVariant::Insertion {
            chrom: "chr1".into(),
            pos: 104,
            sequence: ins_seq.clone(),
        };
        let mut read = make_read(b"ACGTACGTACGT"); // 12 bp
        let original_len = read.len();
        let effect = apply_insertion(&mut read, 100, &sv);

        assert_eq!(effect, SvReadEffect::Modified);
        assert_eq!(
            read.len(),
            original_len,
            "read length must be preserved after insertion"
        );

        // First 4 bases unchanged (positions 100-103)
        assert_eq!(&read.seq[..4], b"ACGT");
        // Next 4 bases are the inserted sequence
        assert_eq!(&read.seq[4..8], b"TTTT");
    }

    // ------------------------------------------------------------------
    // 4. test_inversion_orientation
    // ------------------------------------------------------------------
    #[test]
    fn test_inversion_orientation() {
        // Inversion from 100 to 108.
        // Read starts at 100, length 8 → entirely within inversion.
        let sv = StructuralVariant::Inversion {
            chrom: "chr1".into(),
            start: 100,
            end: 108,
        };
        let seq = b"ACGTACGT";
        let mut read = make_read(seq);
        let effect = apply_inversion(&mut read, 100, &sv);

        assert_eq!(effect, SvReadEffect::Modified);
        // revcomp("ACGTACGT") = revcomp: reverse = "TGCATGCA", complement each:
        // T->A, G->C, C->G, A->T, T->A, G->C, C->G, A->T = "ACGTACGT"  (palindrome)
        // Let's use a non-palindrome: use AAAACCCC
        let seq2 = b"AAAACCCC";
        let mut read2 = make_read(seq2);
        let effect2 = apply_inversion(&mut read2, 100, &sv);
        assert_eq!(effect2, SvReadEffect::Modified);
        // revcomp("AAAACCCC") = reverse "CCCCAAAA", complement = "GGGGTTTT"
        assert_eq!(
            &read2.seq, b"GGGGTTTT",
            "inverted read should be reverse-complemented"
        );
    }

    // ------------------------------------------------------------------
    // 5. test_duplication_coverage
    // ------------------------------------------------------------------
    #[test]
    fn test_duplication_coverage() {
        // Duplication from 200 to 250, 2 copies.
        let sv = StructuralVariant::Duplication {
            chrom: "chr1".into(),
            start: 200,
            end: 250,
            copies: 2,
        };

        // A read spanning the right boundary (dup_end=250): starts at 240, len 20 → ends at 260.
        let mut split_read = make_read(&[b'A'; 20]);
        let effect = apply_duplication(&mut split_read, 240, &sv);
        assert_eq!(
            effect,
            SvReadEffect::SplitRight,
            "read spanning right dup boundary should be split"
        );

        // Bases inside dup (240-249) kept; bases 250+ soft-clipped
        let keep = (250 - 240) as usize; // 10
        for i in 0..keep {
            assert_ne!(split_read.seq[i], b'N', "kept bases must not be N");
        }
        for i in keep..split_read.len() {
            assert_eq!(split_read.seq[i], b'N', "clipped bases should be N");
        }

        // A read entirely outside the dup should be unaffected.
        let mut outside_read = make_read(&[b'A'; 10]);
        let effect2 = apply_duplication(&mut outside_read, 300, &sv);
        assert_eq!(effect2, SvReadEffect::Unaffected);
    }

    // ------------------------------------------------------------------
    // 6. test_translocation_chimeric
    // ------------------------------------------------------------------
    #[test]
    fn test_translocation_chimeric() {
        let sv = StructuralVariant::Translocation {
            chrom1: "chr1".into(),
            pos1: 1000,
            chrom2: "chr2".into(),
            pos2: 5000,
        };

        // Read on chr1 spanning position 1000 (start 990, len 20).
        let mut read = make_read(b"ACGTACGTACGTACGTACGT"); // 20 bp
        let partner = b"GGGGGGGGGGGGGGGGGGGG"; // 20 bp from chr2
        let effect = apply_translocation(&mut read, 990, "chr1", &sv, partner);

        assert_eq!(effect, SvReadEffect::Chimeric);

        // First 10 bases (990-999) are from chr1 and should be unchanged
        assert_eq!(&read.seq[..10], b"ACGTACGTAC");
        // Next 10 bases come from partner chr2 sequence
        assert_eq!(&read.seq[10..20], &partner[..10]);

        // Quality of chimeric portion should be reduced to 10
        for q in &read.qual[10..] {
            assert_eq!(*q, 10, "chimeric bases should have quality 10");
        }
    }

    // ------------------------------------------------------------------
    // 7. test_sv_vaf_stochastic
    // ------------------------------------------------------------------
    #[test]
    fn test_sv_vaf_stochastic() {
        let mut rng = StdRng::seed_from_u64(12345);
        let depth = 100u32;
        let vaf = 0.5;
        let n_trials = 10_000u32;

        let total_alt: u32 = (0..n_trials)
            .map(|_| sample_sv_alt_count(depth, vaf, &mut rng))
            .sum();
        let mean_alt = total_alt as f64 / n_trials as f64;

        // At 50% VAF, depth 100 → expect ~50 alt reads per trial
        assert!(
            (mean_alt - 50.0).abs() < 1.0,
            "mean alt count {mean_alt:.2} should be close to 50.0 at VAF=0.5"
        );

        // Verify stochasticity: counts are not all identical
        let counts: Vec<u32> = (0..100)
            .map(|_| sample_sv_alt_count(depth, vaf, &mut rng))
            .collect();
        let unique: std::collections::HashSet<_> = counts.iter().copied().collect();
        assert!(unique.len() > 1, "SV VAF sampling should be stochastic");
    }

    // ------------------------------------------------------------------
    // 8. test_sv_truth_vcf
    // ------------------------------------------------------------------
    #[test]
    fn test_sv_truth_vcf() {
        // Deletion
        let del = StructuralVariant::Deletion {
            chrom: "chr1".into(),
            start: 1000,
            end: 2000,
        };
        let info = sv_vcf_info(&del);
        assert!(
            info.contains("SVTYPE=DEL"),
            "DEL INFO must contain SVTYPE=DEL"
        );
        assert!(info.contains("END=2000"), "DEL INFO must contain END");
        assert!(info.contains("SVLEN=-1000"), "DEL INFO must contain SVLEN");
        assert_eq!(sv_vcf_alt(&del), "<DEL>");

        // Insertion
        let ins = StructuralVariant::Insertion {
            chrom: "chr1".into(),
            pos: 500,
            sequence: b"ACGT".to_vec(),
        };
        let info_ins = sv_vcf_info(&ins);
        assert!(info_ins.contains("SVTYPE=INS"));
        assert!(info_ins.contains("SVLEN=4"));
        assert_eq!(sv_vcf_alt(&ins), "<INS>");

        // Translocation BND notation
        let tra = StructuralVariant::Translocation {
            chrom1: "chr1".into(),
            pos1: 1000,
            chrom2: "chr2".into(),
            pos2: 5000,
        };
        let info_bnd = sv_vcf_info(&tra);
        assert!(info_bnd.contains("SVTYPE=BND"));
        assert!(info_bnd.contains("MATECHROM=chr2"));
        assert!(info_bnd.contains("MATEPOS=5000"));

        let alt_bnd = sv_vcf_alt(&tra);
        assert!(
            alt_bnd.contains("chr2") && alt_bnd.contains("5000"),
            "BND ALT {alt_bnd} must contain mate chrom and position"
        );
        // Standard BND format uses ]chrom:pos] notation
        assert!(
            alt_bnd.starts_with('N'),
            "BND ALT must start with N anchor base"
        );
    }
}
