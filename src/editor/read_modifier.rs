//! Read-level modification logic for BAM editing.
//!
//! This module provides functions to apply variants (SNVs, indels) to individual
//! BAM records represented as [`RecordBuf`]. It handles updating the sequence,
//! quality scores, CIGAR string, MD tag, and NM tag.

use noodles_sam::alignment::{
    record::{
        cigar::{op::Kind, Op},
        data::field::Tag,
    },
    record_buf::{data::field::Value, Cigar, Data, QualityScores, Sequence},
    RecordBuf,
};

use crate::core::types::MutationType;

/// Result of applying a modification to a read.
#[derive(Debug, Clone, PartialEq)]
pub enum ModifyResult {
    /// Read was successfully modified.
    Modified,
    /// Read was not modified (no overlap, wrong ref base, etc.).
    Unchanged,
    /// Modification was attempted but skipped (e.g., low MQ, low BQ).
    Skipped,
}

/// Apply an SNV to a [`RecordBuf`].
///
/// Replaces the base at the variant position, keeping the original quality score.
/// Updates MD and NM tags if present.
///
/// Returns [`ModifyResult::Modified`] if the read was changed.
pub fn apply_snv(record: &mut RecordBuf, mutation: &MutationType) -> ModifyResult {
    let (pos, ref_base, alt_base) = match mutation {
        MutationType::Snv {
            pos,
            ref_base,
            alt_base,
        } => (*pos, *ref_base, *alt_base),
        _ => return ModifyResult::Unchanged,
    };

    // Get the 0-based alignment start position.
    let align_start = match record.alignment_start() {
        Some(p) => usize::from(p) as u64 - 1,   // noodles is 1-based
        None => return ModifyResult::Unchanged, // unmapped
    };

    // Compute the read length from the sequence.
    let seq_len = record.sequence().len();
    if seq_len == 0 {
        return ModifyResult::Unchanged;
    }

    // Find which read offset corresponds to the variant position, accounting for CIGAR.
    let offset = match cigar_ref_to_read_offset(record.cigar().as_ref(), align_start, pos) {
        Some(o) => o,
        None => return ModifyResult::Unchanged,
    };

    if offset >= seq_len {
        return ModifyResult::Unchanged;
    }

    // Check base quality threshold (BQ < 20 → skip).
    let bq = record
        .quality_scores()
        .as_ref()
        .get(offset)
        .copied()
        .unwrap_or(0);
    if bq < 20 {
        return ModifyResult::Skipped;
    }

    // Check the reference base matches.
    let current_base = sequence_base(record.sequence(), offset);
    if current_base != ref_base {
        return ModifyResult::Unchanged;
    }

    // Apply the substitution.
    let mut new_seq: Vec<u8> = record.sequence().as_ref().to_vec();
    new_seq[offset] = alt_base;

    let new_quals: Vec<u8> = record.quality_scores().as_ref().to_vec();

    *record.sequence_mut() = Sequence::from(new_seq.as_slice());
    *record.quality_scores_mut() = QualityScores::from(new_quals);

    // Update NM tag: increment by 1 for the introduced mismatch.
    update_nm_tag(record, 1);

    // Update MD tag to reflect the substitution.
    update_md_for_snv(record, offset, ref_base, align_start);

    ModifyResult::Modified
}

/// Apply an insertion to a [`RecordBuf`].
///
/// Inserts bases into the sequence at the variant position, updates CIGAR, MD, NM.
/// The read is soft-clipped at the end to maintain a similar reference span.
pub fn apply_insertion(record: &mut RecordBuf, mutation: &MutationType) -> ModifyResult {
    let (pos, ref_seq, alt_seq) = match mutation {
        MutationType::Indel {
            pos,
            ref_seq,
            alt_seq,
        } if alt_seq.len() > ref_seq.len() => (*pos, ref_seq.clone(), alt_seq.clone()),
        _ => return ModifyResult::Unchanged,
    };

    let align_start = match record.alignment_start() {
        Some(p) => usize::from(p) as u64 - 1,
        None => return ModifyResult::Unchanged,
    };

    let seq_len = record.sequence().len();
    if seq_len == 0 {
        return ModifyResult::Unchanged;
    }

    let offset = match cigar_ref_to_read_offset(record.cigar().as_ref(), align_start, pos) {
        Some(o) => o,
        None => return ModifyResult::Unchanged,
    };

    if offset >= seq_len {
        return ModifyResult::Unchanged;
    }

    // Check reference base at offset.
    if offset < ref_seq.len() || sequence_base(record.sequence(), offset) != ref_seq[0] {
        // Only do a soft check; if the anchor base doesn't match, skip.
        if !ref_seq.is_empty() && sequence_base(record.sequence(), offset) != ref_seq[0] {
            return ModifyResult::Unchanged;
        }
    }

    let insert_len = alt_seq.len() - ref_seq.len();

    let orig_seq: Vec<u8> = record.sequence().as_ref().to_vec();
    let orig_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

    // Build new sequence with insertion.
    let mut new_seq = Vec::with_capacity(seq_len + insert_len);
    new_seq.extend_from_slice(&orig_seq[..offset]);
    new_seq.extend_from_slice(&alt_seq);
    if offset + ref_seq.len() < orig_seq.len() {
        new_seq.extend_from_slice(&orig_seq[offset + ref_seq.len()..]);
    }

    // Build new quality scores (use average neighbor quality for inserted bases).
    let avg_qual = if offset > 0 && offset < orig_qual.len() {
        ((orig_qual[offset - 1] as u16 + orig_qual[offset] as u16) / 2) as u8
    } else {
        orig_qual.get(offset).copied().unwrap_or(30)
    };

    let mut new_qual = Vec::with_capacity(seq_len + insert_len);
    new_qual.extend_from_slice(&orig_qual[..offset]);
    for _ in 0..alt_seq.len() {
        new_qual.push(avg_qual);
    }
    if offset + ref_seq.len() < orig_qual.len() {
        new_qual.extend_from_slice(&orig_qual[offset + ref_seq.len()..]);
    }

    // Truncate to original length (soft-clip excess at end).
    new_seq.truncate(seq_len);
    new_qual.truncate(seq_len);

    *record.sequence_mut() = Sequence::from(new_seq.as_slice());
    *record.quality_scores_mut() = QualityScores::from(new_qual);

    // Update CIGAR to include the insertion.
    update_cigar_for_insertion(record, offset, insert_len, seq_len);

    // NM increases by insert_len.
    update_nm_tag(record, insert_len as i32);

    // Clear MD tag (recalculation is complex after indels; just remove it).
    clear_md_tag(record);

    ModifyResult::Modified
}

/// Apply a deletion to a [`RecordBuf`].
///
/// Removes bases from the sequence at the variant position, updates CIGAR, MD, NM.
pub fn apply_deletion(record: &mut RecordBuf, mutation: &MutationType) -> ModifyResult {
    let (pos, ref_seq, alt_seq) = match mutation {
        MutationType::Indel {
            pos,
            ref_seq,
            alt_seq,
        } if ref_seq.len() > alt_seq.len() => (*pos, ref_seq.clone(), alt_seq.clone()),
        _ => return ModifyResult::Unchanged,
    };

    let align_start = match record.alignment_start() {
        Some(p) => usize::from(p) as u64 - 1,
        None => return ModifyResult::Unchanged,
    };

    let seq_len = record.sequence().len();
    if seq_len == 0 {
        return ModifyResult::Unchanged;
    }

    let offset = match cigar_ref_to_read_offset(record.cigar().as_ref(), align_start, pos) {
        Some(o) => o,
        None => return ModifyResult::Unchanged,
    };

    if offset >= seq_len {
        return ModifyResult::Unchanged;
    }

    let del_len = ref_seq.len() - alt_seq.len();

    let orig_seq: Vec<u8> = record.sequence().as_ref().to_vec();
    let orig_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

    // Build new sequence with deletion (remove del_len bases after the anchor).
    let anchor_end = offset + alt_seq.len();
    let deleted_end = offset + ref_seq.len();

    let mut new_seq = Vec::with_capacity(seq_len);
    new_seq.extend_from_slice(&orig_seq[..anchor_end]);
    if deleted_end < orig_seq.len() {
        new_seq.extend_from_slice(&orig_seq[deleted_end..]);
    }

    let mut new_qual = Vec::with_capacity(seq_len);
    new_qual.extend_from_slice(&orig_qual[..anchor_end]);
    if deleted_end < orig_qual.len() {
        new_qual.extend_from_slice(&orig_qual[deleted_end..]);
    }

    // Pad with low-quality Ns to maintain original length.
    let fill_qual = orig_qual.last().copied().unwrap_or(30);
    while new_seq.len() < seq_len {
        new_seq.push(b'N');
        new_qual.push(fill_qual);
    }
    new_seq.truncate(seq_len);
    new_qual.truncate(seq_len);

    *record.sequence_mut() = Sequence::from(new_seq.as_slice());
    *record.quality_scores_mut() = QualityScores::from(new_qual);

    // Update CIGAR to include the deletion.
    update_cigar_for_deletion(record, offset, del_len, seq_len);

    // NM increases by del_len.
    update_nm_tag(record, del_len as i32);

    // Clear MD tag.
    clear_md_tag(record);

    ModifyResult::Modified
}

// ---------------------------------------------------------------------------
// CIGAR helpers
// ---------------------------------------------------------------------------

/// Given a CIGAR and alignment start, compute the read offset for a given
/// reference position. Returns `None` if the position is not covered.
pub fn cigar_ref_to_read_offset(cigar: &[Op], align_start: u64, ref_pos: u64) -> Option<usize> {
    if ref_pos < align_start {
        return None;
    }

    let mut ref_cursor = align_start;
    let mut read_cursor: usize = 0;

    for op in cigar {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // Advances both reference and read.
                if ref_pos >= ref_cursor && ref_pos < ref_cursor + len as u64 {
                    let delta = (ref_pos - ref_cursor) as usize;
                    return Some(read_cursor + delta);
                }
                ref_cursor += len as u64;
                read_cursor += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                // Only advances read.
                read_cursor += len;
            }
            Kind::Deletion | Kind::Skip => {
                // Only advances reference.
                if ref_pos >= ref_cursor && ref_pos < ref_cursor + len as u64 {
                    return None; // position is in a deletion
                }
                ref_cursor += len as u64;
            }
            Kind::HardClip | Kind::Pad => {
                // Advances neither.
            }
        }
    }

    None
}

/// Compute the reference span consumed by a list of CIGAR ops.
#[allow(dead_code)]
pub fn cigar_ref_span(cigar: &[Op]) -> u64 {
    cigar
        .iter()
        .filter_map(|op| op.kind().consumes_reference().then_some(op.len() as u64))
        .sum()
}

/// Compute the read (query) span consumed by a list of CIGAR ops.
#[allow(dead_code)]
pub fn cigar_read_span(cigar: &[Op]) -> usize {
    cigar
        .iter()
        .filter_map(|op| op.kind().consumes_read().then_some(op.len()))
        .sum()
}

/// Format CIGAR ops into a string like "5M2I143M".
#[allow(dead_code)]
pub fn cigar_to_string(ops: &[Op]) -> String {
    ops.iter()
        .map(|op| {
            let ch = match op.kind() {
                Kind::Match => 'M',
                Kind::Insertion => 'I',
                Kind::Deletion => 'D',
                Kind::Skip => 'N',
                Kind::SoftClip => 'S',
                Kind::HardClip => 'H',
                Kind::Pad => 'P',
                Kind::SequenceMatch => '=',
                Kind::SequenceMismatch => 'X',
            };
            format!("{}{}", op.len(), ch)
        })
        .collect::<Vec<_>>()
        .join("")
}

// ---------------------------------------------------------------------------
// CIGAR update helpers
// ---------------------------------------------------------------------------

/// Update the CIGAR string to include an insertion at `read_offset`.
fn update_cigar_for_insertion(
    record: &mut RecordBuf,
    read_offset: usize,
    insert_len: usize,
    read_len: usize,
) {
    let old_ops: Vec<Op> = record.cigar().as_ref().to_vec();
    let mut new_ops: Vec<Op> = Vec::new();

    let mut read_cursor: usize = 0;
    let mut inserted = false;

    for op in &old_ops {
        if inserted {
            new_ops.push(*op);
            continue;
        }

        let op_read_len = if op.kind().consumes_read() {
            op.len()
        } else {
            0
        };

        if read_cursor + op_read_len > read_offset && !inserted {
            // Split this op at the insertion point.
            let before = read_offset - read_cursor;
            let after = op_read_len - before;

            if before > 0 {
                new_ops.push(Op::new(op.kind(), before));
            }
            new_ops.push(Op::new(Kind::Insertion, insert_len));
            inserted = true;

            // The "after" part might need to be truncated because we inserted
            // bases and the read was truncated to original length.
            let remaining_read = read_len.saturating_sub(read_offset + insert_len);
            if op.kind().consumes_read() {
                let after_clamped = after.min(remaining_read);
                if after_clamped > 0 {
                    new_ops.push(Op::new(op.kind(), after_clamped));
                }
            } else if op.kind().consumes_reference() {
                // Deletion/skip: keep as-is after insertion.
                new_ops.push(*op);
            }

            read_cursor += op_read_len;
        } else {
            new_ops.push(*op);
            read_cursor += op_read_len;
        }
    }

    if !inserted {
        // Append insertion at the end.
        new_ops.push(Op::new(Kind::Insertion, insert_len));
    }

    // Merge consecutive ops of the same kind.
    let merged = merge_cigar_ops(new_ops);
    *record.cigar_mut() = merged.into_iter().collect::<Cigar>();
}

/// Update the CIGAR string to include a deletion at `read_offset`.
fn update_cigar_for_deletion(
    record: &mut RecordBuf,
    read_offset: usize,
    del_len: usize,
    _read_len: usize,
) {
    let old_ops: Vec<Op> = record.cigar().as_ref().to_vec();
    let mut new_ops: Vec<Op> = Vec::new();

    let mut read_cursor: usize = 0;
    let mut deleted = false;

    for op in &old_ops {
        if deleted {
            new_ops.push(*op);
            continue;
        }

        let op_read_len = if op.kind().consumes_read() {
            op.len()
        } else {
            0
        };

        if op.kind().consumes_read() && read_cursor + op_read_len > read_offset && !deleted {
            let before = read_offset - read_cursor;
            let after = op_read_len - before;

            if before > 0 {
                new_ops.push(Op::new(op.kind(), before));
            }
            new_ops.push(Op::new(Kind::Deletion, del_len));
            deleted = true;

            if after > 0 {
                new_ops.push(Op::new(op.kind(), after));
            }

            read_cursor += op_read_len;
        } else {
            new_ops.push(*op);
            read_cursor += op_read_len;
        }
    }

    if !deleted {
        new_ops.push(Op::new(Kind::Deletion, del_len));
    }

    let merged = merge_cigar_ops(new_ops);
    *record.cigar_mut() = merged.into_iter().collect::<Cigar>();
}

/// Merge consecutive CIGAR ops of the same kind.
fn merge_cigar_ops(ops: Vec<Op>) -> Vec<Op> {
    let mut result: Vec<Op> = Vec::new();
    for op in ops {
        match result.last_mut() {
            Some(last) if last.kind() == op.kind() => {
                *last = Op::new(last.kind(), last.len() + op.len());
            }
            _ => result.push(op),
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Tag update helpers
// ---------------------------------------------------------------------------

/// Update (increment) the NM tag by `delta`. If missing, set to `delta`.
fn update_nm_tag(record: &mut RecordBuf, delta: i32) {
    let current = record
        .data()
        .get(&Tag::EDIT_DISTANCE)
        .and_then(|v| v.as_int())
        .unwrap_or(0);
    let new_nm = (current + delta as i64).max(0) as i32;
    // Rebuild data with updated NM.
    let mut data: Data = record
        .data()
        .iter()
        .filter(|(t, _)| *t != Tag::EDIT_DISTANCE)
        .map(|(t, v)| (t, v.clone()))
        .collect();
    data.insert(Tag::EDIT_DISTANCE, Value::from(new_nm));
    *record.data_mut() = data;
}

/// Update the MD tag for an SNV substitution.
///
/// This is a simplified updater: if an existing MD tag is present and parseable,
/// it replaces the ref base at the position with the substitution marker.
/// If absent or complex, clears the MD tag.
fn update_md_for_snv(record: &mut RecordBuf, read_offset: usize, ref_base: u8, _align_start: u64) {
    // Get current MD string.
    let md_str = match record.data().get(&Tag::MISMATCHED_POSITIONS) {
        Some(Value::String(s)) => String::from_utf8(s.iter().copied().collect()).ok(),
        _ => None,
    };

    if let Some(md) = md_str {
        // Attempt to inject the ref base substitution into the MD string.
        let new_md = inject_snv_into_md(&md, read_offset, ref_base);
        let mut data: Data = record
            .data()
            .iter()
            .filter(|(t, _)| *t != Tag::MISMATCHED_POSITIONS)
            .map(|(t, v)| (t, v.clone()))
            .collect();
        data.insert(
            Tag::MISMATCHED_POSITIONS,
            Value::String(new_md.into_bytes().into()),
        );
        *record.data_mut() = data;
    }
    // If no MD tag exists, leave absent (no-op).
}

/// A single token from a parsed MD string.
#[derive(Debug, PartialEq)]
enum MdToken {
    /// A run of matching bases of the given length.
    Run(usize),
    /// A single-base mismatch showing the reference base.
    Mismatch(u8),
    /// Deleted reference bases (the sequence between `^` and the next digit).
    Deletion(Vec<u8>),
}

/// Parse an MD string into a list of tokens.
///
/// The MD format is `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`.
/// A numeric token represents that many consecutive matching bases.
/// An uppercase letter is a single-base mismatch (the reference base).
/// `^` followed by uppercase letters is a deleted reference sequence.
fn parse_md(md: &str) -> Vec<MdToken> {
    let bytes = md.as_bytes();
    let mut tokens = Vec::new();
    let mut i = 0;

    while i < bytes.len() {
        if bytes[i].is_ascii_digit() {
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            let n: usize = md[start..i].parse().unwrap_or(0);
            tokens.push(MdToken::Run(n));
        } else if bytes[i] == b'^' {
            i += 1;
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_uppercase() {
                i += 1;
            }
            tokens.push(MdToken::Deletion(bytes[start..i].to_vec()));
        } else if bytes[i].is_ascii_uppercase() {
            tokens.push(MdToken::Mismatch(bytes[i]));
            i += 1;
        } else {
            i += 1;
        }
    }

    tokens
}

/// Serialise a list of MD tokens back into a string.
///
/// Inserts a zero-length run between adjacent non-run tokens to keep the
/// output valid per the SAM spec.
fn serialise_md(tokens: &[MdToken]) -> String {
    let mut out = String::new();
    let mut last_was_nonrun = false;

    for token in tokens {
        match token {
            MdToken::Run(n) => {
                out.push_str(&n.to_string());
                last_was_nonrun = false;
            }
            MdToken::Mismatch(b) => {
                if last_was_nonrun {
                    out.push('0');
                }
                out.push(*b as char);
                last_was_nonrun = true;
            }
            MdToken::Deletion(bases) => {
                if last_was_nonrun {
                    out.push('0');
                }
                out.push('^');
                for &b in bases {
                    out.push(b as char);
                }
                last_was_nonrun = true;
            }
        }
    }

    out
}

/// Inject a ref-base SNV substitution into an MD string at `read_offset`.
///
/// Walks the token list counting consumed read bases. When the run that
/// covers `read_offset` is found, it is split and the ref base is inserted.
///
/// Example: MD="150", offset=50, ref_base=A → "50A99".
fn inject_snv_into_md(md: &str, read_offset: usize, ref_base: u8) -> String {
    let tokens = parse_md(md);
    let mut new_tokens: Vec<MdToken> = Vec::new();
    // `cursor` counts read bases consumed so far. Deletions do not consume
    // read bases, so they are passed through without advancing the cursor.
    let mut cursor: usize = 0;
    let mut injected = false;

    for token in tokens {
        if injected {
            new_tokens.push(token);
            continue;
        }

        match token {
            MdToken::Run(n) => {
                if cursor + n > read_offset {
                    // The target offset falls inside this run.
                    let before = read_offset - cursor;
                    let after = n - before - 1; // -1 for the substituted base
                                                // Always emit the run before the mismatch, even when
                                                // before==0, so the MD string starts with a digit per spec.
                    new_tokens.push(MdToken::Run(before));
                    new_tokens.push(MdToken::Mismatch(ref_base));
                    new_tokens.push(MdToken::Run(after));
                    injected = true;
                    cursor += n;
                } else {
                    cursor += n;
                    new_tokens.push(MdToken::Run(n));
                }
            }
            MdToken::Mismatch(b) => {
                // A mismatch consumes one read base.
                if cursor == read_offset {
                    // The SNV lands on an already-recorded mismatch; update
                    // the ref base.
                    new_tokens.push(MdToken::Mismatch(ref_base));
                    injected = true;
                } else {
                    new_tokens.push(MdToken::Mismatch(b));
                }
                cursor += 1;
            }
            MdToken::Deletion(bases) => {
                // Deletions do not consume read bases.
                new_tokens.push(MdToken::Deletion(bases));
            }
        }
    }

    serialise_md(&new_tokens)
}

/// Remove the MD tag from the record.
fn clear_md_tag(record: &mut RecordBuf) {
    let data: Data = record
        .data()
        .iter()
        .filter(|(t, _)| *t != Tag::MISMATCHED_POSITIONS)
        .map(|(t, v)| (t, v.clone()))
        .collect();
    *record.data_mut() = data;
}

// ---------------------------------------------------------------------------
// Sequence helpers
// ---------------------------------------------------------------------------

/// Get a base at `offset` from a noodles [`Sequence`] as a `u8`.
fn sequence_base(seq: &Sequence, offset: usize) -> u8 {
    seq.as_ref().get(offset).copied().unwrap_or(b'N')
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::{
        record::Flags,
        record_buf::{Cigar, Sequence},
        RecordBuf,
    };

    fn make_record(seq: &[u8], align_start: u64, cigar_str: &str) -> RecordBuf {
        let cigar_ops = parse_cigar_str(cigar_str);
        RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_alignment_start(Position::new(align_start as usize + 1).unwrap())
            .set_cigar(cigar_ops.into_iter().collect::<Cigar>())
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(vec![30u8; seq.len()]))
            .build()
    }

    fn make_record_with_nm(seq: &[u8], align_start: u64, cigar_str: &str, nm: i32) -> RecordBuf {
        let cigar_ops = parse_cigar_str(cigar_str);
        let mut data = Data::default();
        data.insert(Tag::EDIT_DISTANCE, Value::from(nm));
        RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_alignment_start(Position::new(align_start as usize + 1).unwrap())
            .set_cigar(cigar_ops.into_iter().collect::<Cigar>())
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(vec![30u8; seq.len()]))
            .set_data(data)
            .build()
    }

    fn parse_cigar_str(s: &str) -> Vec<Op> {
        crate::io::bam::parse_cigar(s).expect("valid cigar")
    }

    fn get_seq(record: &RecordBuf) -> Vec<u8> {
        record.sequence().as_ref().to_vec()
    }

    fn get_cigar_str(record: &RecordBuf) -> String {
        let ops: Vec<Op> = record.cigar().as_ref().to_vec();
        cigar_to_string(&ops)
    }

    // ------------------------------------------------------------------
    // SNV tests
    // ------------------------------------------------------------------

    #[test]
    fn test_apply_snv_basic() {
        // Read: ACGTACGT at position 100
        let mut record = make_record(b"ACGTACGT", 100, "8M");
        let mutation = MutationType::Snv {
            pos: 102,
            ref_base: b'G',
            alt_base: b'T',
        };
        let result = apply_snv(&mut record, &mutation);
        assert_eq!(result, ModifyResult::Modified);
        let seq = get_seq(&record);
        assert_eq!(&seq, b"ACTTACGT", "base at offset 2 should change G->T");
    }

    #[test]
    fn test_apply_snv_preserves_quality() {
        let mut record = make_record(b"ACGTACGT", 100, "8M");
        let orig_quals: Vec<u8> = record.quality_scores().as_ref().to_vec();
        let mutation = MutationType::Snv {
            pos: 102,
            ref_base: b'G',
            alt_base: b'T',
        };
        apply_snv(&mut record, &mutation);
        let new_quals: Vec<u8> = record.quality_scores().as_ref().to_vec();
        assert_eq!(
            orig_quals, new_quals,
            "quality scores must be preserved for SNV"
        );
    }

    #[test]
    fn test_apply_snv_wrong_ref_base_unchanged() {
        let mut record = make_record(b"ACGTACGT", 100, "8M");
        let mutation = MutationType::Snv {
            pos: 102,
            ref_base: b'T',
            alt_base: b'A',
        };
        let result = apply_snv(&mut record, &mutation);
        assert_eq!(result, ModifyResult::Unchanged);
        let seq = get_seq(&record);
        assert_eq!(&seq, b"ACGTACGT", "unchanged when ref base doesn't match");
    }

    #[test]
    fn test_apply_snv_out_of_range() {
        let mut record = make_record(b"ACGTACGT", 100, "8M");
        let mutation = MutationType::Snv {
            pos: 200,
            ref_base: b'A',
            alt_base: b'T',
        };
        let result = apply_snv(&mut record, &mutation);
        assert_eq!(result, ModifyResult::Unchanged);
    }

    #[test]
    fn test_apply_snv_updates_nm_tag() {
        let mut record = make_record_with_nm(b"ACGTACGT", 100, "8M", 0);
        let mutation = MutationType::Snv {
            pos: 100,
            ref_base: b'A',
            alt_base: b'C',
        };
        apply_snv(&mut record, &mutation);
        let nm = record
            .data()
            .get(&Tag::EDIT_DISTANCE)
            .and_then(|v| v.as_int())
            .unwrap_or(-1);
        assert_eq!(nm, 1, "NM should be incremented to 1 after SNV");
    }

    #[test]
    fn test_apply_snv_nm_increments_from_existing() {
        let mut record = make_record_with_nm(b"ACGTACGT", 100, "8M", 2);
        let mutation = MutationType::Snv {
            pos: 100,
            ref_base: b'A',
            alt_base: b'C',
        };
        apply_snv(&mut record, &mutation);
        let nm = record
            .data()
            .get(&Tag::EDIT_DISTANCE)
            .and_then(|v| v.as_int())
            .unwrap_or(-1);
        assert_eq!(nm, 3, "NM should increment from existing value");
    }

    // ------------------------------------------------------------------
    // Indel tests
    // ------------------------------------------------------------------

    #[test]
    fn test_apply_insertion_basic() {
        // Read: ACGTACGT at position 100, insert "TT" after offset 2 (pos 102)
        let mut record = make_record(b"ACGTACGT", 100, "8M");
        let mutation = MutationType::Indel {
            pos: 102,
            ref_seq: vec![b'G'],
            alt_seq: vec![b'G', b'T', b'T'],
        };
        let result = apply_insertion(&mut record, &mutation);
        assert_eq!(result, ModifyResult::Modified);
        let seq = get_seq(&record);
        // Original: ACGTACGT
        // Insert GTT at offset 2: AC + GTT + TACGT → truncated to 8: ACGTTTAC
        assert_eq!(seq.len(), 8, "read length preserved after insertion");
        // Check the insertion is in the sequence.
        assert_eq!(&seq[..5], b"ACGTT", "insertion visible at offset 2");
    }

    #[test]
    fn test_apply_insertion_updates_cigar() {
        let mut record = make_record(b"ACGTACGT", 100, "8M");
        let mutation = MutationType::Indel {
            pos: 102,
            ref_seq: vec![b'G'],
            alt_seq: vec![b'G', b'T', b'T'],
        };
        apply_insertion(&mut record, &mutation);
        let cigar = get_cigar_str(&record);
        assert!(
            cigar.contains('I'),
            "CIGAR should contain an insertion op: {cigar}"
        );
    }

    #[test]
    fn test_apply_deletion_basic() {
        // Read: ACGTACGT at position 100, delete GT at position 102
        let mut record = make_record(b"ACGTACGT", 100, "8M");
        let mutation = MutationType::Indel {
            pos: 102,
            ref_seq: vec![b'G', b'T'],
            alt_seq: vec![b'G'],
        };
        let result = apply_deletion(&mut record, &mutation);
        assert_eq!(result, ModifyResult::Modified);
        let seq = get_seq(&record);
        assert_eq!(seq.len(), 8, "read length preserved after deletion");
        // After deletion: ACG + ACGT + N (padding) = ACGACGTN
        assert_eq!(&seq[..3], b"ACG", "bases before deletion unchanged");
    }

    #[test]
    fn test_apply_deletion_updates_cigar() {
        let mut record = make_record(b"ACGTACGT", 100, "8M");
        let mutation = MutationType::Indel {
            pos: 102,
            ref_seq: vec![b'G', b'T'],
            alt_seq: vec![b'G'],
        };
        apply_deletion(&mut record, &mutation);
        let cigar = get_cigar_str(&record);
        assert!(
            cigar.contains('D'),
            "CIGAR should contain a deletion op: {cigar}"
        );
    }

    // ------------------------------------------------------------------
    // CIGAR offset calculation tests
    // ------------------------------------------------------------------

    #[test]
    fn test_cigar_ref_to_read_offset_simple_match() {
        let ops = parse_cigar_str("10M");
        // Position 105 in a read aligned at 100 → offset 5.
        let offset = cigar_ref_to_read_offset(&ops, 100, 105);
        assert_eq!(offset, Some(5));
    }

    #[test]
    fn test_cigar_ref_to_read_offset_with_soft_clip() {
        // 2S8M: 2 soft-clipped bases then 8M
        let ops = parse_cigar_str("2S8M");
        // Read starts at ref 100 (soft clip doesn't count for alignment start).
        // Offset for pos 102 = 2(softclip) + 2(match) = 4.
        let offset = cigar_ref_to_read_offset(&ops, 100, 102);
        assert_eq!(offset, Some(4));
    }

    #[test]
    fn test_cigar_ref_to_read_offset_deletion_returns_none() {
        // 5M2D5M: 5 matches, 2 deleted, 5 matches
        let ops = parse_cigar_str("5M2D5M");
        // Position in the deleted region should return None.
        let offset = cigar_ref_to_read_offset(&ops, 100, 106); // pos 106 is in 2D
        assert_eq!(offset, None);
    }

    #[test]
    fn test_cigar_ref_to_read_offset_before_start() {
        let ops = parse_cigar_str("10M");
        let offset = cigar_ref_to_read_offset(&ops, 100, 50);
        assert_eq!(offset, None);
    }

    // ------------------------------------------------------------------
    // Utility tests
    // ------------------------------------------------------------------

    #[test]
    fn test_cigar_to_string() {
        let ops = parse_cigar_str("5M2I143M");
        let s = cigar_to_string(&ops);
        assert_eq!(s, "5M2I143M");
    }

    #[test]
    fn test_merge_cigar_ops() {
        let ops = vec![
            Op::new(Kind::Match, 5),
            Op::new(Kind::Match, 3),
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::Match, 5),
        ];
        let merged = merge_cigar_ops(ops);
        assert_eq!(merged.len(), 3);
        assert_eq!(merged[0], Op::new(Kind::Match, 8));
        assert_eq!(merged[1], Op::new(Kind::Insertion, 2));
        assert_eq!(merged[2], Op::new(Kind::Match, 5));
    }

    // ------------------------------------------------------------------
    // Quality averaging overflow (T052)
    // ------------------------------------------------------------------

    #[test]
    fn test_insertion_quality_avg_no_overflow() {
        // Qualities of 200 on both sides of the insertion point must average to
        // 200. The old expression `(a as u16 + b as u16) as u8 / 2` would cast
        // 400 to u8 (wrapping to 144) then divide, giving 72.
        let cigar_ops = parse_cigar_str("8M");
        let mut record = RecordBuf::builder()
            .set_flags(noodles_sam::alignment::record::Flags::empty())
            .set_alignment_start(Position::new(101).unwrap())
            .set_cigar(cigar_ops.into_iter().collect::<Cigar>())
            .set_sequence(Sequence::from(b"ACGTACGT".as_ref()))
            .set_quality_scores(QualityScores::from(vec![200u8; 8]))
            .build();
        let mutation = MutationType::Indel {
            pos: 102,
            ref_seq: vec![b'G'],
            alt_seq: vec![b'G', b'T'],
        };
        let result = apply_insertion(&mut record, &mutation);
        assert_eq!(result, ModifyResult::Modified);
        // All quality scores must remain 200 after padding inserted bases.
        let quals: Vec<u8> = record.quality_scores().as_ref().to_vec();
        assert!(
            quals.iter().all(|&q| q == 200),
            "all qualities must be 200 after insertion; got {quals:?}"
        );
    }

    // ------------------------------------------------------------------
    // inject_snv_into_md tests (T054)
    // ------------------------------------------------------------------

    #[test]
    fn test_inject_snv_at_position_zero() {
        // MD="150", offset=0, ref=T → "0T149"
        let result = inject_snv_into_md("150", 0, b'T');
        assert_eq!(result, "0T149");
    }

    #[test]
    fn test_inject_snv_in_middle_of_all_match() {
        // MD="150", offset=50, ref=A → "50A99"
        let result = inject_snv_into_md("150", 50, b'A');
        assert_eq!(result, "50A99");
    }

    #[test]
    fn test_inject_snv_at_last_base() {
        // MD="150", offset=149, ref=C → "149C0"
        let result = inject_snv_into_md("150", 149, b'C');
        assert_eq!(result, "149C0");
    }

    #[test]
    fn test_inject_snv_into_read_with_existing_mismatch() {
        // MD="75A74", read offsets: 0..74 match, offset 75 is mismatch A,
        // offsets 76..149 match. Inject at offset=100, ref=G → "75A24G49"
        let result = inject_snv_into_md("75A74", 100, b'G');
        assert_eq!(result, "75A24G49");
    }

    #[test]
    fn test_inject_snv_into_read_with_deletion() {
        // MD="10^AC20", offset=25 (past the deletion).
        // Token breakdown: Run(10), Deletion(AC), Run(20).
        // Read offsets: 0..9 → Run(10), deletion consumes no read bases,
        // then 10..29 → Run(20). Offset 25 is at read position 25 - 10 = 15
        // inside the second run.
        // Expected: "10^AC15G4"
        let result = inject_snv_into_md("10^AC20", 25, b'G');
        assert_eq!(result, "10^AC15G4");
    }

    #[test]
    fn test_inject_snv_into_read_with_deletion_before_offset() {
        // MD="5^GT10", inject at offset=12, ref=C.
        // Tokens: Run(5), Deletion(GT), Run(10).
        // Offsets: 0..4 match (Run 5), deletion skipped, 5..14 match (Run 10).
        // Offset 12 → 7 into the second run → before=7, after=2.
        // Expected: "5^GT7C2"
        let result = inject_snv_into_md("5^GT10", 12, b'C');
        assert_eq!(result, "5^GT7C2");
    }
}
