//! Sequencing error injection: indels, cycle-position errors, context-dependent
//! multipliers, strand bias, and phasing bursts.

use rand::Rng;

/// Parameters controlling sequencing indel error injection.
///
/// Indel errors are distinct from somatic or germline variants: they are
/// instrument artefacts that cause spurious insertions or deletions in the
/// read sequence at a very low per-base rate.
// Used by downstream tasks in EPIC-ERROR-MODEL (T155 ErrorOrchestrator).
#[allow(dead_code)]
pub struct IndelErrorModel {
    /// Per-base probability that an indel event occurs at each position.
    pub indel_rate: f64,
    /// Fraction of indel events that are insertions; the remainder are deletions.
    pub insertion_fraction: f64,
    /// Maximum indel length. Lengths are drawn from Geometric(0.7), capped here.
    pub max_length: usize,
}

/// Inject sequencing indel errors into a read.
///
/// Iterates positions in reverse order (avoiding index-shift bugs), rolls
/// each position against `model.indel_rate`, and records events. Events are
/// then applied in reverse position order.
///
/// After all events the sequence and quality vectors are truncated or padded
/// to exactly `read_length`, preserving the fixed-length Illumina contract.
/// Insertions trim the 3' end; deletions pad with `b'N'` / quality 0.
///
/// The caller must ensure `seq.len() == qual.len() == read_length` on entry.
// Used by downstream tasks in EPIC-ERROR-MODEL (T155 ErrorOrchestrator).
#[allow(dead_code)]
pub fn inject_indel_errors(
    seq: &mut Vec<u8>,
    qual: &mut Vec<u8>,
    read_length: usize,
    model: &IndelErrorModel,
    rng: &mut impl Rng,
) {
    // Collect events. Iterate in reverse so that applying them back in reverse
    // position order is straightforward and avoids shifting already-processed
    // indices.
    let mut events: Vec<(usize, bool, usize)> = Vec::new(); // (pos, is_insertion, length)

    for pos in (0..seq.len()).rev() {
        if rng.random::<f64>() < model.indel_rate {
            let is_insertion = rng.random::<f64>() < model.insertion_fraction;
            // Draw length from Geometric(0.7) capped at max_length.
            let mut len = 1usize;
            while len < model.max_length && rng.random::<f64>() < 0.3 {
                len += 1;
            }
            events.push((pos, is_insertion, len));
        }
    }

    // Events are already in reverse position order (we iterated in reverse).
    // Apply each event so that earlier positions are not disturbed by later ones.
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

    for (pos, is_insertion, len) in events {
        // Earlier deletions may have shortened the sequence so that `pos` is
        // now out of bounds. Skip stale events rather than underflowing.
        if pos >= seq.len() {
            continue;
        }

        if is_insertion {
            // Insert `len` random bases into seq and copies of qual[pos] into
            // qual immediately after position `pos`.
            let insert_pos = (pos + 1).min(seq.len());
            let q = qual[pos];
            for k in 0..len {
                let base = BASES[rng.random_range(0..4)];
                seq.insert(insert_pos + k, base);
                qual.insert(insert_pos + k, q);
            }
        } else {
            // Delete min(len, remaining) bases starting at `pos`.
            let del_count = len.min(seq.len() - pos);
            seq.drain(pos..pos + del_count);
            qual.drain(pos..pos + del_count);
        }
    }

    // Enforce fixed-length contract.
    if seq.len() > read_length {
        seq.truncate(read_length);
        qual.truncate(read_length);
    } else {
        while seq.len() < read_length {
            seq.push(b'N');
            qual.push(0);
        }
    }

    debug_assert_eq!(
        seq.len(),
        read_length,
        "seq length mismatch after indel injection"
    );
    debug_assert_eq!(
        qual.len(),
        read_length,
        "qual length mismatch after indel injection"
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn make_read(len: usize) -> (Vec<u8>, Vec<u8>) {
        let seq = vec![b'A'; len];
        let qual = vec![30u8; len];
        (seq, qual)
    }

    #[test]
    fn test_indel_rate_zero_leaves_read_unchanged() {
        // At rate 0.0 no events fire; the read must be identical to the input.
        let model = IndelErrorModel {
            indel_rate: 0.0,
            insertion_fraction: 0.5,
            max_length: 3,
        };
        let read_length = 50;
        let mut rng = StdRng::seed_from_u64(42);
        let original_seq = vec![b'A'; read_length];
        let original_qual = vec![30u8; read_length];
        let mut seq = original_seq.clone();
        let mut qual = original_qual.clone();
        inject_indel_errors(&mut seq, &mut qual, read_length, &model, &mut rng);
        assert_eq!(seq, original_seq, "seq should be unchanged at rate 0.0");
        assert_eq!(qual, original_qual, "qual should be unchanged at rate 0.0");
    }

    #[test]
    fn test_high_indel_rate_modifies_reads() {
        // At rate 0.5 on 20-bp reads almost every read will be modified.
        let model = IndelErrorModel {
            indel_rate: 0.5,
            insertion_fraction: 0.5,
            max_length: 2,
        };
        let read_length = 20;
        let mut rng = StdRng::seed_from_u64(99);
        let mut changed = 0usize;
        let n_reads = 1_000;
        for _ in 0..n_reads {
            let original = vec![b'A'; read_length];
            let mut seq = original.clone();
            let mut qual = vec![30u8; read_length];
            inject_indel_errors(&mut seq, &mut qual, read_length, &model, &mut rng);
            if seq != original {
                changed += 1;
            }
        }
        assert!(
            changed > 900,
            "expected >90% reads modified at rate 0.5, got {}/{}",
            changed,
            n_reads
        );
    }

    #[test]
    fn test_fixed_length_contract() {
        let model = IndelErrorModel {
            indel_rate: 0.1,
            insertion_fraction: 0.5,
            max_length: 3,
        };
        let read_length = 100;
        let mut rng = StdRng::seed_from_u64(7);

        for _ in 0..10_000 {
            let (mut seq, mut qual) = make_read(read_length);
            inject_indel_errors(&mut seq, &mut qual, read_length, &model, &mut rng);
            assert_eq!(
                seq.len(),
                read_length,
                "seq length {} != {}",
                seq.len(),
                read_length
            );
            assert_eq!(
                qual.len(),
                read_length,
                "qual length {} != {}",
                qual.len(),
                read_length
            );
        }
    }

    #[test]
    fn test_only_insertions_when_fraction_one() {
        // insertion_fraction: 1.0 means all events are insertions. With random
        // bases inserted, at least some output reads should contain non-A bases.
        let model = IndelErrorModel {
            indel_rate: 0.3,
            insertion_fraction: 1.0,
            max_length: 1,
        };
        let read_length = 20;
        let mut rng = StdRng::seed_from_u64(7);
        let mut has_non_a = false;
        for _ in 0..100 {
            let mut seq = vec![b'A'; read_length];
            let mut qual = vec![30u8; read_length];
            inject_indel_errors(&mut seq, &mut qual, read_length, &model, &mut rng);
            assert_eq!(seq.len(), read_length, "fixed-length contract violated");
            if seq.iter().any(|&b| b != b'A') {
                has_non_a = true;
            }
        }
        assert!(
            has_non_a,
            "insertions should produce non-A bases in at least one read"
        );
    }

    #[test]
    fn test_length_distribution() {
        // With indel_rate 1.0 and insertion_fraction 1.0, every position triggers
        // an insertion. Insertions grow the sequence before truncation. Because the
        // Geometric(0.7) distribution has P(len=1) = 0.7, reads with max_length 5
        // will almost always be extended by exactly 1 base at each event. After
        // truncation back to read_length, the output must still differ from all-A
        // input in over 50% of runs (inserted random bases replace trailing A's).
        let model = IndelErrorModel {
            indel_rate: 1.0,
            insertion_fraction: 1.0,
            max_length: 5,
        };
        let read_length = 20;
        let n_reads = 1_000;
        let mut rng = StdRng::seed_from_u64(13);
        let mut modified = 0usize;

        for _ in 0..n_reads {
            let original = vec![b'A'; read_length];
            let mut seq = original.clone();
            let mut qual = vec![30u8; read_length];
            inject_indel_errors(&mut seq, &mut qual, read_length, &model, &mut rng);
            assert_eq!(seq.len(), read_length, "fixed-length contract violated");
            if seq != original {
                modified += 1;
            }
        }

        // At rate 1.0 with insertions only, virtually every read gets modified.
        assert!(
            modified > n_reads / 2,
            "expected >50% reads modified at rate 1.0, got {}/{}",
            modified,
            n_reads
        );
    }
}
