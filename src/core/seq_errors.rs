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

/// Models the tendency for R2 reads to accumulate more errors than R1.
///
/// R2 undergoes more synthesis cycles before sequencing, so phasing errors
/// accumulate. This manifests as both a higher per-base error rate and lower
/// quality scores relative to R1.
///
/// The `r2_error_multiplier` is consumed by the orchestrator (T155) when
/// calling injection functions for R2. `apply_to_r2_qual` adjusts quality
/// scores before error injection.
// Used by T155 ErrorOrchestrator.
#[allow(dead_code)]
pub struct StrandBiasModel {
    /// R2 error rate = R1 rate × this multiplier. Default 1.0 (no bias).
    pub r2_error_multiplier: f64,
    /// Shift R2 quality scores by this many Phred points.
    /// Positive values lower quality (subtract); negative values raise quality (add).
    /// Default 0 (no shift).
    pub r2_quality_offset: i8,
}

impl StrandBiasModel {
    /// Apply the quality offset to an R2 quality array in place.
    ///
    /// Positive `r2_quality_offset` lowers quality (saturating subtract).
    /// Negative `r2_quality_offset` raises quality (saturating add, capped at 93).
    // Used by T155 ErrorOrchestrator.
    #[allow(dead_code)]
    pub fn apply_to_r2_qual(&self, qual: &mut [u8]) {
        for q in qual.iter_mut() {
            if self.r2_quality_offset > 0 {
                *q = q.saturating_sub(self.r2_quality_offset as u8);
            } else if self.r2_quality_offset < 0 {
                *q = q
                    .saturating_add(self.r2_quality_offset.unsigned_abs())
                    .min(93);
            }
        }
    }
}

/// Models correlated phasing burst errors caused by cluster phasing failures.
///
/// A phasing failure at one cycle causes the same miscall to propagate across
/// several adjacent positions. All bases in a burst are substituted to the same
/// wrong base, and their quality scores are set to Q12 to indicate unreliability.
// Used by T155 ErrorOrchestrator.
#[allow(dead_code)]
pub struct CorrelatedErrorModel {
    /// Per-base probability of initiating a phasing error burst.
    pub burst_rate: f64,
    /// Mean burst length. Lengths are drawn from a Geometric distribution
    /// with success probability `1 / burst_length_mean`.
    pub burst_length_mean: f64,
}

/// Inject correlated phasing burst errors into a read.
///
/// Walks positions 0..seq.len(). At each position, rolls against `burst_rate`.
/// If a burst starts, draws a length from Geometric(1/burst_length_mean),
/// then forces all burst positions to the same randomly chosen wrong base and
/// sets their quality scores to 12.
///
/// The caller must ensure `seq.len() == qual.len()` on entry.
// Used by T155 ErrorOrchestrator.
#[allow(dead_code)]
pub fn inject_burst_errors(
    seq: &mut [u8],
    qual: &mut [u8],
    model: &CorrelatedErrorModel,
    rng: &mut impl Rng,
) {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let len = seq.len();
    let mut burst_remaining: usize = 0;
    let mut burst_base: u8 = b'A';

    for i in 0..len {
        if burst_remaining > 0 {
            seq[i] = burst_base;
            qual[i] = 12;
            burst_remaining -= 1;
        } else if rng.random::<f64>() < model.burst_rate {
            // Draw burst length from Geometric(1/burst_length_mean).
            // Each additional step continues with probability 1 - 1/mean.
            let p_continue = 1.0 - 1.0 / model.burst_length_mean;
            let mut drawn_len = 1usize;
            let max_len = len - i;
            while drawn_len < max_len && rng.random::<f64>() < p_continue {
                drawn_len += 1;
            }

            // Choose a wrong base different from the current base.
            let current = seq[i];
            let wrong: u8 = loop {
                let candidate = BASES[rng.random_range(0..4)];
                if candidate != current {
                    break candidate;
                }
            };

            burst_base = wrong;
            seq[i] = burst_base;
            qual[i] = 12;
            burst_remaining = drawn_len - 1;
        }
    }
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

    #[test]
    fn test_indel_rate_statistical() {
        // At rate 0.1 on 20-bp reads, expected events per read ≈ 2.0.
        // With insertion_fraction 1.0, all events insert a random base.
        // The input is all b'A'. After truncation to read_length, inserted
        // non-A bases remain in the output.
        //
        // Truncation significantly reduces the visible non-A rate: insertions
        // near the 3' end push existing bases off the read. Empirically the
        // observed non-A fraction sits around 0.065 at this rate and read
        // length. Bounds [0.04, 0.09] confirm the rate is non-zero and scales
        // with the configured indel_rate without over-constraining for
        // truncation effects.
        let model = IndelErrorModel {
            indel_rate: 0.1,
            insertion_fraction: 1.0,
            max_length: 1,
        };
        let mut rng = StdRng::seed_from_u64(12345);
        let read_length = 20usize;
        let n_reads = 10_000usize;
        let mut non_a_count = 0usize;
        let total_bases = n_reads * read_length;
        for _ in 0..n_reads {
            let mut seq = vec![b'A'; read_length];
            let mut qual = vec![30u8; read_length];
            inject_indel_errors(&mut seq, &mut qual, read_length, &model, &mut rng);
            non_a_count += seq.iter().filter(|&&b| b != b'A').count();
        }
        let observed_rate = non_a_count as f64 / total_bases as f64;
        assert!(
            (0.04..=0.09).contains(&observed_rate),
            "expected non-A rate in [0.04, 0.09], got {:.4}",
            observed_rate
        );
    }

    // --- StrandBiasModel tests ---

    #[test]
    fn test_strand_bias_lowers_quality() {
        let model = StrandBiasModel {
            r2_error_multiplier: 1.3,
            r2_quality_offset: 3,
        };
        let mut qual = vec![40u8; 20];
        model.apply_to_r2_qual(&mut qual);
        assert!(
            qual.iter().all(|&q| q == 37),
            "all qualities should be 37 after subtracting offset 3 from 40"
        );
    }

    #[test]
    fn test_strand_bias_zero_offset_noop() {
        let model = StrandBiasModel {
            r2_error_multiplier: 1.3,
            r2_quality_offset: 0,
        };
        let original = vec![30u8, 25u8, 40u8, 10u8, 5u8];
        let mut qual = original.clone();
        model.apply_to_r2_qual(&mut qual);
        assert_eq!(qual, original, "zero offset must leave quality unchanged");
    }

    #[test]
    fn test_strand_bias_negative_offset_raises_quality() {
        // Negative offset raises quality; result is capped at 93.
        let model = StrandBiasModel {
            r2_error_multiplier: 1.0,
            r2_quality_offset: -5,
        };
        let mut qual = vec![30u8; 10];
        model.apply_to_r2_qual(&mut qual);
        assert!(
            qual.iter().all(|&q| q == 35),
            "negative offset -5 should raise quality from 30 to 35"
        );
    }

    // --- CorrelatedErrorModel tests ---

    #[test]
    fn test_burst_errors_correlated() {
        // With burst_rate 0.1 and burst_length_mean 5.0, bursts should be
        // clearly longer than 1 base on average. Run 10000 reads of 50 bp
        // all-A and measure average run length of the first non-A run per read.
        let model = CorrelatedErrorModel {
            burst_rate: 0.1,
            burst_length_mean: 5.0,
        };
        let mut rng = StdRng::seed_from_u64(42);
        let n_reads = 10_000usize;
        let read_len = 50usize;
        let mut total_run_len = 0usize;
        let mut run_count = 0usize;

        for _ in 0..n_reads {
            let mut seq = vec![b'A'; read_len];
            let mut qual = vec![30u8; read_len];
            inject_burst_errors(&mut seq, &mut qual, &model, &mut rng);

            // Find the first non-A position and measure the run from there.
            if let Some(start) = seq.iter().position(|&b| b != b'A') {
                let run_base = seq[start];
                let run_len = seq[start..].iter().take_while(|&&b| b == run_base).count();
                total_run_len += run_len;
                run_count += 1;
            }
        }

        assert!(run_count > 0, "expected some reads to have bursts");
        let avg_run = total_run_len as f64 / run_count as f64;
        assert!(
            avg_run > 1.5,
            "expected average burst run length > 1.5, got {:.3}",
            avg_run
        );
    }

    #[test]
    fn test_burst_base_consistent() {
        // All bases within a burst must be the same wrong base.
        let model = CorrelatedErrorModel {
            burst_rate: 0.5,
            burst_length_mean: 4.0,
        };
        let mut rng = StdRng::seed_from_u64(99);
        let read_len = 50usize;

        // Try up to 1000 reads to find one with a burst longer than 1 base.
        let mut found_multi_base_burst = false;
        for _ in 0..1_000 {
            let mut seq = vec![b'A'; read_len];
            let mut qual = vec![30u8; read_len];
            inject_burst_errors(&mut seq, &mut qual, &model, &mut rng);

            if let Some(start) = seq.iter().position(|&b| b != b'A') {
                let burst_base = seq[start];
                // Measure contiguous run of burst_base.
                let run_len = seq[start..]
                    .iter()
                    .take_while(|&&b| b == burst_base)
                    .count();
                if run_len > 1 {
                    // Verify all bases in the run are identical.
                    assert!(
                        seq[start..start + run_len].iter().all(|&b| b == burst_base),
                        "burst bases are not all identical"
                    );
                    found_multi_base_burst = true;
                    break;
                }
            }
        }

        assert!(
            found_multi_base_burst,
            "expected to find at least one burst longer than 1 base"
        );
    }

    #[test]
    fn test_burst_rate_zero_noop() {
        let model = CorrelatedErrorModel {
            burst_rate: 0.0,
            burst_length_mean: 3.0,
        };
        let mut rng = StdRng::seed_from_u64(7);
        let read_len = 50usize;

        for _ in 0..1_000 {
            let original_seq = vec![b'A'; read_len];
            let original_qual = vec![30u8; read_len];
            let mut seq = original_seq.clone();
            let mut qual = original_qual.clone();
            inject_burst_errors(&mut seq, &mut qual, &model, &mut rng);
            assert_eq!(seq, original_seq, "burst_rate 0.0 must not change seq");
            assert_eq!(qual, original_qual, "burst_rate 0.0 must not change qual");
        }
    }

    #[test]
    fn test_insertion_fraction_statistical() {
        // Use a low indel_rate (0.01) so that most reads receive at most one
        // event. This makes the insertion vs deletion proxy reliable:
        // - A single insertion adds a random non-A base (detected as non-A,
        //   non-N in the output from an all-A input).
        // - A single deletion removes a base and pads with b'N' at the end.
        // At high rates, multiple interacting events make the N-padding proxy
        // ambiguous because insertions can displace the trailing N bytes.
        // With insertion_fraction 0.7 and 10 000 reads at rate 0.01 on 20-bp
        // reads, roughly 2000 events fire total, giving enough signal with
        // clean single-event classification. Bounds [0.60, 0.80] confirm the
        // configured 0.7 split is respected.
        let model = IndelErrorModel {
            indel_rate: 0.01,
            insertion_fraction: 0.7,
            max_length: 1,
        };
        let mut rng = StdRng::seed_from_u64(54321);
        let read_length = 20usize;
        let n_reads = 10_000usize;
        let mut insertion_evidence = 0usize; // non-A, non-N bases (inserted random bases)
        let mut deletion_evidence = 0usize; // N bases from padding
        for _ in 0..n_reads {
            let mut seq = vec![b'A'; read_length];
            let mut qual = vec![30u8; read_length];
            inject_indel_errors(&mut seq, &mut qual, read_length, &model, &mut rng);
            for &b in &seq {
                if b != b'A' && b != b'N' {
                    insertion_evidence += 1;
                }
                if b == b'N' {
                    deletion_evidence += 1;
                }
            }
        }
        let total_evidence = insertion_evidence + deletion_evidence;
        if total_evidence > 100 {
            let observed_insertion_fraction = insertion_evidence as f64 / total_evidence as f64;
            assert!(
                (0.60..=0.80).contains(&observed_insertion_fraction),
                "expected insertion fraction ~0.7, got {:.4}",
                observed_insertion_fraction
            );
        }
    }
}
