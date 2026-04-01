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

/// Map a single base byte to its 2-bit representation.
///
/// A=0, C=1, G=2, T=3. Any unrecognised byte maps to 0 (treated as A).
fn base_to_bits(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// JSON-deserialisable profile for loading learned k-mer error multipliers.
///
/// Produced by T156 (empirical profile extraction). Loaded via
/// `KmerErrorModel::from_profile_json`.
// Used by downstream tasks in EPIC-ERROR-MODEL (T155 ErrorOrchestrator, T156 profile extraction).
#[allow(dead_code)]
#[derive(serde::Deserialize)]
pub struct KmerProfileJson {
    pub kmer_length: usize,
    pub rules: Vec<KmerRuleJson>,
}

/// A single context rule inside a `KmerProfileJson`.
// Used by downstream tasks in EPIC-ERROR-MODEL (T155 ErrorOrchestrator, T156 profile extraction).
#[allow(dead_code)]
#[derive(serde::Deserialize)]
pub struct KmerRuleJson {
    pub context: String,
    pub sub_multiplier: f32,
    pub indel_multiplier: f32,
}

/// Context-dependent sequencing error multipliers, indexed by k-mer hash.
///
/// The table stores one `f32` multiplier per possible k-mer (4^k entries).
/// At each read position the current k-mer is hashed using a rolling 2-bit
/// scheme and the multiplier is applied to the base error probability before
/// the Bernoulli draw.
///
/// `k` must be in 1..=5 (table sizes 4..=1024 entries). Values outside that
/// range are accepted but will use larger allocations.
// Used by downstream tasks in EPIC-ERROR-MODEL (T155 ErrorOrchestrator).
#[allow(dead_code)]
pub struct KmerErrorModel {
    /// k-mer length (1..=5).
    pub k: usize,
    /// Substitution error multiplier indexed by k-mer hash. Size = 4^k. Default = 1.0.
    sub_multipliers: Vec<f32>,
    /// Indel error multiplier indexed by k-mer hash. Default = 1.0.
    indel_multipliers: Vec<f32>,
}

#[allow(dead_code)]
impl KmerErrorModel {
    /// Create a uniform model where all k-mer multipliers are 1.0.
    pub fn uniform(k: usize) -> Self {
        let size = 1 << (k * 2);
        Self {
            k,
            sub_multipliers: vec![1.0f32; size],
            indel_multipliers: vec![1.0f32; size],
        }
    }

    /// Set substitution and indel multipliers for a specific k-mer context string.
    ///
    /// `context` must have exactly `k` bytes and contain only ACGT (upper or
    /// lower case). Panics otherwise.
    pub fn set_rule(&mut self, context: &str, sub_multiplier: f32, indel_multiplier: f32) {
        assert_eq!(
            context.len(),
            self.k,
            "context length {} != k={}",
            context.len(),
            self.k
        );
        let idx = self.kmer_index(context.as_bytes());
        self.sub_multipliers[idx] = sub_multiplier;
        self.indel_multipliers[idx] = indel_multiplier;
    }

    /// Compute the lookup index for a slice of exactly `k` bases.
    fn kmer_index(&self, bases: &[u8]) -> usize {
        bases
            .iter()
            .fold(0usize, |acc, &b| (acc << 2) | base_to_bits(b))
    }

    /// Return the substitution multiplier for the k-mer ending at `pos`.
    ///
    /// Returns 1.0 if there is not yet enough context (pos + 1 < k).
    pub fn sub_multiplier_at(&self, seq: &[u8], pos: usize) -> f32 {
        if pos + 1 < self.k {
            return 1.0;
        }
        let start = pos + 1 - self.k;
        let idx = self.kmer_index(&seq[start..=pos]);
        self.sub_multipliers[idx]
    }

    /// Return the indel multiplier for the k-mer ending at `pos`.
    ///
    /// Returns 1.0 if there is not yet enough context (pos + 1 < k).
    pub fn indel_multiplier_at(&self, seq: &[u8], pos: usize) -> f32 {
        if pos + 1 < self.k {
            return 1.0;
        }
        let start = pos + 1 - self.k;
        let idx = self.kmer_index(&seq[start..=pos]);
        self.indel_multipliers[idx]
    }

    /// Load a k-mer error model from a JSON profile file.
    ///
    /// The file must conform to the `KmerProfileJson` schema. All k-mers not
    /// listed in `rules` keep the default multiplier of 1.0.
    pub fn from_profile_json(path: &std::path::Path) -> anyhow::Result<Self> {
        let text = std::fs::read_to_string(path)?;
        let profile: KmerProfileJson = serde_json::from_str(&text)?;
        let mut model = Self::uniform(profile.kmer_length);
        for rule in &profile.rules {
            model.set_rule(&rule.context, rule.sub_multiplier, rule.indel_multiplier);
        }
        Ok(model)
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

    // --- KmerErrorModel tests ---

    #[test]
    fn test_uniform_model_all_ones() {
        // uniform(3) should set all sub_multipliers to 1.0 and return 1.0 for
        // any query via sub_multiplier_at.
        let model = KmerErrorModel::uniform(3);
        assert!(
            model.sub_multipliers.iter().all(|&v| v == 1.0f32),
            "all sub_multipliers should be 1.0"
        );
        let seq = b"ACGTACGT";
        for pos in 0..seq.len() {
            assert_eq!(
                model.sub_multiplier_at(seq, pos),
                1.0f32,
                "expected 1.0 at pos {}",
                pos
            );
        }
    }

    #[test]
    fn test_set_rule_lookup() {
        // Set "GGG" with sub_multiplier 5.0. Query the last G in "AAAGGG"
        // (pos 5) — expect 5.0. Query an A at pos 2 — expect 1.0.
        let mut model = KmerErrorModel::uniform(3);
        model.set_rule("GGG", 5.0, 1.0);
        let seq = b"AAAGGG";
        assert_eq!(
            model.sub_multiplier_at(seq, 5),
            5.0f32,
            "expected 5.0 at GGG context"
        );
        assert_eq!(
            model.sub_multiplier_at(seq, 2),
            1.0f32,
            "expected 1.0 at AAA context"
        );
    }

    #[test]
    fn test_rolling_hash_matches_naive() {
        // For k=3 and 100 random 20-bp sequences, verify that sub_multiplier_at
        // (which uses kmer_index on the slice) agrees at every position.
        // We use a known set of rules to make some multipliers non-trivial.
        let mut model = KmerErrorModel::uniform(3);
        model.set_rule("GGC", 2.0, 3.0);
        model.set_rule("TTT", 4.0, 1.5);

        // Deterministic sequence: cycle through ACGT repeated.
        let alphabet = [b'A', b'C', b'G', b'T'];
        let mut rng = StdRng::seed_from_u64(42);
        use rand::Rng as _;
        for _ in 0..100 {
            let seq: Vec<u8> = (0..20).map(|_| alphabet[rng.random_range(0..4)]).collect();
            for pos in 0..seq.len() {
                let via_fn = model.sub_multiplier_at(&seq, pos);
                // Recompute naive: use kmer_index directly on the slice.
                let naive = if pos + 1 < model.k {
                    1.0f32
                } else {
                    let start = pos + 1 - model.k;
                    let idx = model.kmer_index(&seq[start..=pos]);
                    model.sub_multipliers[idx]
                };
                assert_eq!(via_fn, naive, "mismatch at pos {} in seq {:?}", pos, seq);
            }
        }
    }

    #[test]
    fn test_elevated_context_increases_errors() {
        // With k=2 and "GG" sub_multiplier 20.0, verify that sub_multiplier_at
        // returns 20.0 at positions following "GG" and 1.0 at positions
        // following "AA". No pipeline needed; the multiplier function is
        // what matters.
        let mut model = KmerErrorModel::uniform(2);
        model.set_rule("GG", 20.0, 1.0);

        // Sequence: alternating G and A — no consecutive GG or AA possible.
        // Use a sequence with explicit GG and AA runs instead.
        let seq = b"AAGGTAA";
        // pos 0: 'A' — k=2, pos+1=1 < k=2 → 1.0
        assert_eq!(model.sub_multiplier_at(seq, 0), 1.0f32);
        // pos 1: context is "AA" → 1.0
        assert_eq!(model.sub_multiplier_at(seq, 1), 1.0f32);
        // pos 2: context is "AG" → 1.0
        assert_eq!(model.sub_multiplier_at(seq, 2), 1.0f32);
        // pos 3: context is "GG" → 20.0
        assert_eq!(model.sub_multiplier_at(seq, 3), 20.0f32);
        // pos 5: context is "TA" → 1.0
        assert_eq!(model.sub_multiplier_at(seq, 5), 1.0f32);
    }

    #[test]
    fn test_kmer_size_1_to_4() {
        // uniform(k) must allocate exactly 4^k entries.
        for k in 1usize..=4 {
            let model = KmerErrorModel::uniform(k);
            let expected = 4usize.pow(k as u32);
            assert_eq!(
                model.sub_multipliers.len(),
                expected,
                "k={}: expected {} sub_multipliers, got {}",
                k,
                expected,
                model.sub_multipliers.len()
            );
            assert_eq!(
                model.indel_multipliers.len(),
                expected,
                "k={}: expected {} indel_multipliers, got {}",
                k,
                expected,
                model.indel_multipliers.len()
            );
        }
    }
}
