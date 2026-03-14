use crate::core::types::{MutationType, Read};

/// Apply a mutation to a read if it overlaps the mutation position.
///
/// Returns true if the read was modified.
pub fn spike_snv(read: &mut Read, read_start: u64, mutation: &MutationType) -> bool {
    match mutation {
        MutationType::Snv { pos, ref_base, alt_base } => {
            let read_end = read_start + read.len() as u64;
            if *pos >= read_start && *pos < read_end {
                let offset = (*pos - read_start) as usize;
                if read.seq[offset] == *ref_base {
                    read.seq[offset] = *alt_base;
                    return true;
                }
            }
            false
        }
        _ => false,
    }
}

/// Apply an MNV (multi-nucleotide variant) to a read.
/// All bases must be on the same read — never partially spike.
pub fn spike_mnv(read: &mut Read, read_start: u64, mutation: &MutationType) -> bool {
    match mutation {
        MutationType::Mnv { pos, ref_seq, alt_seq } => {
            let read_end = read_start + read.len() as u64;
            let mnv_end = *pos + ref_seq.len() as u64;

            // MNV must be fully contained within the read
            if *pos >= read_start && mnv_end <= read_end {
                let offset = (*pos - read_start) as usize;
                // Verify reference bases match
                if read.seq[offset..offset + ref_seq.len()] == ref_seq[..] {
                    read.seq[offset..offset + alt_seq.len()].copy_from_slice(alt_seq);
                    return true;
                }
            }
            false
        }
        _ => false,
    }
}

/// Apply an indel to a read.
///
/// For insertions: insert bases and truncate to maintain read length.
/// For deletions: remove bases and pad with reference sequence.
pub fn spike_indel(
    read: &mut Read,
    read_start: u64,
    mutation: &MutationType,
    ref_seq_after: &[u8],
) -> bool {
    match mutation {
        MutationType::Indel { pos, ref_seq, alt_seq } => {
            let read_end = read_start + read.len() as u64;
            if *pos < read_start || *pos >= read_end {
                return false;
            }

            let offset = (*pos - read_start) as usize;
            let original_len = read.len();

            if alt_seq.len() > ref_seq.len() {
                // Insertion
                let insert_len = alt_seq.len() - ref_seq.len();
                let mut new_seq = Vec::with_capacity(original_len);
                new_seq.extend_from_slice(&read.seq[..offset]);
                new_seq.extend_from_slice(alt_seq);
                if offset + ref_seq.len() < read.seq.len() {
                    new_seq.extend_from_slice(&read.seq[offset + ref_seq.len()..]);
                }
                new_seq.truncate(original_len);
                // Pad qualities for inserted bases (use average of neighbors)
                let mut new_qual = Vec::with_capacity(original_len);
                new_qual.extend_from_slice(&read.qual[..offset]);
                let avg_qual = if offset > 0 && offset < read.qual.len() {
                    (read.qual[offset - 1] + read.qual[offset]) / 2
                } else {
                    read.qual[offset.min(read.qual.len() - 1)]
                };
                for _ in 0..alt_seq.len() {
                    new_qual.push(avg_qual);
                }
                if offset + ref_seq.len() < read.qual.len() {
                    new_qual.extend_from_slice(&read.qual[offset + ref_seq.len()..]);
                }
                new_qual.truncate(original_len);

                read.seq = new_seq;
                read.qual = new_qual;
            } else {
                // Deletion
                let del_len = ref_seq.len() - alt_seq.len();
                let mut new_seq = Vec::with_capacity(original_len);
                new_seq.extend_from_slice(&read.seq[..offset]);
                new_seq.extend_from_slice(alt_seq);
                if offset + ref_seq.len() < read.seq.len() {
                    new_seq.extend_from_slice(&read.seq[offset + ref_seq.len()..]);
                }
                // Pad with reference sequence to maintain read length
                let pad_needed = original_len.saturating_sub(new_seq.len());
                new_seq.extend_from_slice(&ref_seq_after[..pad_needed.min(ref_seq_after.len())]);
                new_seq.truncate(original_len);

                let mut new_qual = Vec::with_capacity(original_len);
                new_qual.extend_from_slice(&read.qual[..offset]);
                let avg_qual = read.qual[offset.min(read.qual.len() - 1)];
                for _ in 0..alt_seq.len() {
                    new_qual.push(avg_qual);
                }
                if offset + ref_seq.len() < read.qual.len() {
                    new_qual.extend_from_slice(&read.qual[offset + ref_seq.len()..]);
                }
                while new_qual.len() < original_len {
                    new_qual.push(avg_qual);
                }
                new_qual.truncate(original_len);

                read.seq = new_seq;
                read.qual = new_qual;
            }

            true
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_read(seq: &[u8]) -> Read {
        Read::new(seq.to_vec(), vec![30; seq.len()])
    }

    #[test]
    fn test_spike_snv_hit() {
        let mut read = make_read(b"ACGTACGT");
        let mutation = MutationType::Snv { pos: 102, ref_base: b'G', alt_base: b'T' };
        assert!(spike_snv(&mut read, 100, &mutation));
        assert_eq!(&read.seq, b"ACTTACGT");
    }

    #[test]
    fn test_spike_snv_miss() {
        let mut read = make_read(b"ACGTACGT");
        let mutation = MutationType::Snv { pos: 200, ref_base: b'G', alt_base: b'T' };
        assert!(!spike_snv(&mut read, 100, &mutation));
        assert_eq!(&read.seq, b"ACGTACGT"); // unchanged
    }

    #[test]
    fn test_spike_snv_wrong_ref() {
        let mut read = make_read(b"ACGTACGT");
        let mutation = MutationType::Snv { pos: 100, ref_base: b'T', alt_base: b'G' };
        assert!(!spike_snv(&mut read, 100, &mutation));
        assert_eq!(&read.seq, b"ACGTACGT"); // unchanged, ref didn't match
    }

    #[test]
    fn test_spike_mnv() {
        let mut read = make_read(b"ACGTACGT");
        let mutation = MutationType::Mnv {
            pos: 101,
            ref_seq: vec![b'C', b'G'],
            alt_seq: vec![b'T', b'A'],
        };
        assert!(spike_mnv(&mut read, 100, &mutation));
        assert_eq!(&read.seq, b"ATATACGT"); // pos 101 offset 1, CG -> TA
    }

    #[test]
    fn test_spike_mnv_partial_overlap_rejected() {
        let mut read = make_read(b"ACGT");
        // MNV extends beyond read end
        let mutation = MutationType::Mnv {
            pos: 103,
            ref_seq: vec![b'T', b'X'],
            alt_seq: vec![b'A', b'A'],
        };
        assert!(!spike_mnv(&mut read, 100, &mutation));
        assert_eq!(&read.seq, b"ACGT"); // unchanged
    }

    #[test]
    fn test_spike_insertion() {
        let mut read = make_read(b"ACGTACGT");
        let mutation = MutationType::Indel {
            pos: 102,
            ref_seq: vec![b'G'],
            alt_seq: vec![b'G', b'T', b'T'],
        };
        assert!(spike_indel(&mut read, 100, &mutation, b""));
        // After insertion: AC + GTT + TACGT -> ACGTTTAC (truncated to 8)
        assert_eq!(read.seq.len(), 8);
        assert_eq!(&read.seq[..5], b"ACGTT"); // the insertion is visible
    }

    #[test]
    fn test_spike_deletion() {
        let mut read = make_read(b"ACGTACGT");
        let mutation = MutationType::Indel {
            pos: 102,
            ref_seq: vec![b'G', b'T'],
            alt_seq: vec![b'G'],
        };
        let ref_after = b"NNNNNNNN"; // padding
        assert!(spike_indel(&mut read, 100, &mutation, ref_after));
        assert_eq!(read.seq.len(), 8);
    }

    #[test]
    fn test_spike_indel_no_overlap() {
        let mut read = make_read(b"ACGT");
        let mutation = MutationType::Indel {
            pos: 200,
            ref_seq: vec![b'A'],
            alt_seq: vec![b'A', b'T'],
        };
        assert!(!spike_indel(&mut read, 100, &mutation, b""));
    }

    #[test]
    fn test_read_length_preserved_after_indel() {
        let mut read = make_read(b"ACGTACGTACGT");
        let original_len = read.len();
        let mutation = MutationType::Indel {
            pos: 103,
            ref_seq: vec![b'T'],
            alt_seq: vec![b'T', b'A', b'A', b'A'],
        };
        spike_indel(&mut read, 100, &mutation, b"NNNN");
        assert_eq!(read.len(), original_len, "read length must be preserved");
        assert_eq!(read.seq.len(), read.qual.len(), "seq and qual must match");
    }
}
