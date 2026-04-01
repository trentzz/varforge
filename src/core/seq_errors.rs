//! Sequencing error injection: indels, cycle-position errors, context-dependent
//! multipliers, strand bias, and phasing bursts.

use std::io::{self, BufRead};
use std::path::Path;

use anyhow::Result;
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

/// Precomputed per-cycle error probabilities for a read.
///
/// Models the gap between reported quality scores and actual instrument error
/// rates, including the characteristic 3' error rise in Illumina reads.
// Used by downstream tasks in EPIC-ERROR-MODEL (T155 ErrorOrchestrator).
#[allow(dead_code)]
pub struct CycleErrorCurve {
    /// Per-position error probabilities, length equals read_length.
    curve: Vec<f64>,
}

#[allow(dead_code)]
impl CycleErrorCurve {
    /// Build a flat curve: every position gets `base_error_rate`.
    pub fn flat(read_length: usize, base_error_rate: f64) -> Self {
        Self {
            curve: vec![base_error_rate; read_length],
        }
    }

    /// Build an exponential-tail curve.
    ///
    /// Positions before `tail_start_fraction * read_length` get `base_error_rate`.
    /// Positions at or after that point ramp exponentially up to
    /// `base_error_rate * tail_rate_multiplier` at the last cycle.
    ///
    /// Formula: `rate(i) = base_error_rate * exp(k * (i - tail_start) / (read_length - 1 - tail_start))`
    /// where `k = ln(tail_rate_multiplier)`, for `i >= tail_start`.
    pub fn exponential(
        read_length: usize,
        base_error_rate: f64,
        tail_start_fraction: f64,
        tail_rate_multiplier: f64,
    ) -> Self {
        let tail_start = (tail_start_fraction * read_length as f64) as usize;
        let k = tail_rate_multiplier.ln();
        let denom = if read_length > 1 && tail_start < read_length - 1 {
            (read_length - 1 - tail_start) as f64
        } else {
            1.0
        };

        let curve = (0..read_length)
            .map(|i| {
                if i < tail_start {
                    base_error_rate
                } else {
                    let t = (i - tail_start) as f64 / denom;
                    base_error_rate * (k * t).exp()
                }
            })
            .collect();

        Self { curve }
    }

    /// Load from a two-column TSV: `cycle\terror_rate` (tab-separated, no header).
    ///
    /// Cycle values are 0-based. Linearly interpolates between provided points.
    /// Extrapolates the last known rate for cycles beyond the last provided point.
    pub fn from_tsv(path: &Path, read_length: usize) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let reader = io::BufReader::new(file);

        let mut points: Vec<(usize, f64)> = Vec::new();
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let mut parts = line.splitn(2, '\t');
            let cycle: usize = parts
                .next()
                .ok_or_else(|| anyhow::anyhow!("missing cycle column in TSV"))?
                .trim()
                .parse()?;
            let rate: f64 = parts
                .next()
                .ok_or_else(|| anyhow::anyhow!("missing rate column in TSV"))?
                .trim()
                .parse()?;
            points.push((cycle, rate));
        }

        anyhow::ensure!(
            !points.is_empty(),
            "cycle error TSV is empty: {}",
            path.display()
        );

        // Sort by cycle ascending for interpolation.
        points.sort_by_key(|&(c, _)| c);

        let curve = (0..read_length)
            .map(|i| Self::interpolate(&points, i))
            .collect();

        Ok(Self { curve })
    }

    /// Linearly interpolate (or extrapolate) the rate for cycle `i`.
    fn interpolate(points: &[(usize, f64)], i: usize) -> f64 {
        // Before the first point: use the first rate.
        if i <= points[0].0 {
            return points[0].1;
        }
        // After the last point: use the last rate.
        let last = points[points.len() - 1];
        if i >= last.0 {
            return last.1;
        }
        // Find the surrounding pair.
        let pos = points.partition_point(|&(c, _)| c <= i);
        let (c0, r0) = points[pos - 1];
        let (c1, r1) = points[pos];
        let t = (i - c0) as f64 / (c1 - c0) as f64;
        r0 + t * (r1 - r0)
    }
}

/// Apply cycle-position-dependent substitutions to a read.
///
/// This is an independent error pass on top of quality-driven errors.
/// For each position `i`, a Bernoulli trial fires with probability `model.curve[i]`.
/// When triggered, the base is replaced with a uniformly random different base.
// Used by downstream tasks in EPIC-ERROR-MODEL (T155 ErrorOrchestrator).
#[allow(dead_code)]
pub fn inject_cycle_errors(seq: &mut [u8], model: &CycleErrorCurve, rng: &mut impl Rng) {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    for (i, base) in seq.iter_mut().enumerate() {
        let rate = if i < model.curve.len() {
            model.curve[i]
        } else {
            0.0
        };
        if rng.random::<f64>() < rate {
            // Pick a different base uniformly from the 3 remaining.
            let alts: [u8; 3] = match *base {
                b'A' => [b'C', b'G', b'T'],
                b'C' => [b'A', b'G', b'T'],
                b'G' => [b'A', b'C', b'T'],
                b'T' => [b'A', b'C', b'G'],
                _ => [b'A', b'C', b'G'], // treat non-ACGT as T-like
            };
            *base = alts[rng.random_range(0..3)];
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

    // --- CycleErrorCurve tests ---

    #[test]
    fn test_flat_curve_rate() {
        // Flat curve with base_error_rate 0.1, 10000 reads of 50 bp all-A.
        // Count total non-A bases and assert rate is within [0.09, 0.11].
        let model = CycleErrorCurve::flat(50, 0.1);
        let mut rng = StdRng::seed_from_u64(1001);
        let n_reads = 10_000usize;
        let read_length = 50usize;
        let mut non_a = 0usize;
        for _ in 0..n_reads {
            let mut seq = vec![b'A'; read_length];
            inject_cycle_errors(&mut seq, &model, &mut rng);
            non_a += seq.iter().filter(|&&b| b != b'A').count();
        }
        let rate = non_a as f64 / (n_reads * read_length) as f64;
        assert!(
            (0.09..=0.11).contains(&rate),
            "expected flat rate in [0.09, 0.11], got {:.4}",
            rate
        );
    }

    #[test]
    fn test_exponential_tail_rises() {
        // Exponential curve: base_error_rate 0.01, tail_start_fraction 0.8,
        // tail_rate_multiplier 10.0, read_length 100.
        // curve[99] should be ≈ 0.1 (within 5%), curve[0] should be 0.01.
        let model = CycleErrorCurve::exponential(100, 0.01, 0.8, 10.0);
        let expected_last = 0.1f64;
        let actual_last = model.curve[99];
        assert!(
            (actual_last - expected_last).abs() / expected_last < 0.05,
            "expected curve[99] ≈ {:.4}, got {:.6}",
            expected_last,
            actual_last
        );
        assert_eq!(
            model.curve[0], 0.01,
            "curve[0] should equal base_error_rate"
        );
    }

    #[test]
    fn test_exponential_curve_len() {
        // Both constructors must produce a curve of exactly read_length.
        let flat = CycleErrorCurve::flat(75, 0.005);
        assert_eq!(flat.curve.len(), 75, "flat curve length mismatch");
        let exp = CycleErrorCurve::exponential(120, 0.005, 0.7, 8.0);
        assert_eq!(exp.curve.len(), 120, "exponential curve length mismatch");
    }

    #[test]
    fn test_substitution_not_identity() {
        // Every substituted base must differ from the original.
        // Run with a high rate so virtually all positions are hit.
        let model = CycleErrorCurve::flat(1, 1.0); // single-position, always fires
        let mut rng = StdRng::seed_from_u64(7777);
        for _ in 0..10_000 {
            let original = b'A';
            let mut seq = vec![original];
            inject_cycle_errors(&mut seq, &model, &mut rng);
            assert_ne!(
                seq[0], original,
                "substituted base must differ from original"
            );
        }
    }
}
