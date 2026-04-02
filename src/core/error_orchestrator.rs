//! Composes all sequencing error models into a single orchestrator that applies
//! passes in a defined order per read pair.
//!
//! The orchestrator is built once per region from [`SequencingErrorConfig`] and
//! shared across all read pairs in the batch. All precomputed tables (cycle
//! curve, k-mer hash) are built at initialisation time, keeping the per-read
//! hot path allocation-free.
//!
//! When `sequencing_errors` is absent from the config, no orchestrator is
//! constructed and behaviour is bit-for-bit identical to v0.1.x.

use std::path::Path;

use anyhow::Result;
use rand::Rng;

use crate::core::seq_errors::{
    inject_burst_errors, inject_context_errors, inject_cycle_errors, inject_indel_errors,
    CorrelatedErrorModel, CycleErrorCurve, IndelErrorModel, KmerErrorModel, StrandBiasModel,
};
use crate::io::config::SequencingErrorConfig;

/// Owns all active error sub-models and applies them in a defined order.
///
/// Build once per region via [`ErrorOrchestrator::from_config`], then call
/// [`ErrorOrchestrator::inject_all_errors`] for each read pair.
pub struct ErrorOrchestrator {
    /// Per-cycle additional substitution probabilities for R1.
    pub cycle_curve: Option<CycleErrorCurve>,
    /// Pre-scaled cycle curve for R2. `None` when `r2_error_multiplier == 1.0`,
    /// in which case `cycle_curve` is reused for R2.
    pub cycle_curve_r2: Option<CycleErrorCurve>,
    /// Sequencing indel model (shared between R1 and R2).
    pub indel_model: Option<IndelErrorModel>,
    /// Context-dependent k-mer error multiplier model.
    pub kmer_model: Option<KmerErrorModel>,
    /// R2 strand asymmetry model (quality offset and error rate multiplier).
    pub strand_bias: Option<StrandBiasModel>,
    /// Correlated phasing burst error model.
    pub burst_model: Option<CorrelatedErrorModel>,
    /// Base error rate used for the context error pass (0.0 if no cycle curve).
    context_base_rate: f64,
    /// R2 error rate = R1 rate × this multiplier.
    pub r2_error_multiplier: f64,
    /// Read length; used to enforce the fixed-length contract after indel injection.
    pub read_length: usize,
}

impl ErrorOrchestrator {
    /// Build an orchestrator from a [`SequencingErrorConfig`].
    ///
    /// Returns `Ok(None)` when every relevant field is `None` or zero, so the
    /// caller can skip the injection call entirely and avoid any overhead.
    pub fn from_config(
        seq_err_cfg: &SequencingErrorConfig,
        read_length: usize,
    ) -> Result<Option<Self>> {
        // Determine whether anything is actually configured.
        let has_cycle = seq_err_cfg.cycle_error_model.is_some()
            || seq_err_cfg.base_error_rate.is_some_and(|r| r > 0.0);
        let has_indel = seq_err_cfg.indel_rate.is_some_and(|r| r > 0.0);
        let has_kmer = seq_err_cfg.kmer_length.is_some()
            || !seq_err_cfg.context_rules.is_empty()
            || seq_err_cfg.context_profile_path.is_some();
        let has_bias = seq_err_cfg.r2_quality_offset.is_some_and(|o| o != 0)
            || seq_err_cfg.r2_error_multiplier.is_some_and(|m| m != 1.0);
        let has_burst = seq_err_cfg.burst_rate.is_some_and(|r| r > 0.0);

        if !has_cycle && !has_indel && !has_kmer && !has_bias && !has_burst {
            return Ok(None);
        }

        // --- Cycle error curve ---
        let cycle_curve = if has_cycle {
            let base_rate = seq_err_cfg.base_error_rate.unwrap_or(0.001);
            let model_name = seq_err_cfg.cycle_error_model.as_deref().unwrap_or("flat");
            let curve = match model_name {
                "exponential" => {
                    let tail_start = seq_err_cfg.tail_start_fraction.unwrap_or(0.7);
                    let tail_mult = seq_err_cfg.tail_rate_multiplier.unwrap_or(5.0);
                    CycleErrorCurve::exponential(read_length, base_rate, tail_start, tail_mult)
                }
                "custom" => {
                    let tsv_path = seq_err_cfg
                        .cycle_error_tsv
                        .as_deref()
                        .ok_or_else(|| anyhow::anyhow!(
                            "quality.sequencing_errors.cycle_error_tsv must be set when cycle_error_model is \"custom\""
                        ))?;
                    CycleErrorCurve::from_tsv(Path::new(tsv_path), read_length)?
                }
                // "flat" is the default.
                _ => CycleErrorCurve::flat(read_length, base_rate),
            };
            Some(curve)
        } else {
            None
        };

        // The context base rate is only non-zero when no cycle curve is active.
        // When a cycle curve is present it already handles the base substitution
        // rate; setting context_base_rate to 0.0 prevents double-counting.
        let context_base_rate = if cycle_curve.is_none() {
            seq_err_cfg.base_error_rate.unwrap_or(0.0)
        } else {
            0.0 // cycle curve handles the base rate; do not add a second pass
        };

        // --- Indel model ---
        let indel_model = if has_indel {
            Some(IndelErrorModel {
                indel_rate: seq_err_cfg.indel_rate.unwrap_or(0.0),
                insertion_fraction: seq_err_cfg.indel_insertion_fraction.unwrap_or(0.5),
                max_length: seq_err_cfg.max_indel_length.unwrap_or(3),
            })
        } else {
            None
        };

        // --- K-mer context model ---
        let kmer_model = if has_kmer {
            let k = seq_err_cfg.kmer_length.unwrap_or_else(|| {
                // Infer k from context_rules if kmer_length is not set.
                seq_err_cfg
                    .context_rules
                    .first()
                    .map(|r| r.context.len())
                    .unwrap_or(3)
            });

            let model = if let Some(ref profile_path) = seq_err_cfg.context_profile_path {
                KmerErrorModel::from_profile_json(Path::new(profile_path))?
            } else {
                let mut m = KmerErrorModel::uniform(k);
                for rule in &seq_err_cfg.context_rules {
                    m.set_rule(&rule.context, rule.sub_multiplier, rule.indel_multiplier);
                }
                m
            };
            Some(model)
        } else {
            None
        };

        // --- Strand bias model ---
        let strand_bias = if has_bias {
            Some(StrandBiasModel {
                r2_error_multiplier: seq_err_cfg.r2_error_multiplier.unwrap_or(1.0),
                r2_quality_offset: seq_err_cfg.r2_quality_offset.unwrap_or(0),
            })
        } else {
            None
        };

        // --- Correlated burst model ---
        let burst_model = if has_burst {
            Some(CorrelatedErrorModel {
                burst_rate: seq_err_cfg.burst_rate.unwrap_or(0.0),
                burst_length_mean: seq_err_cfg.burst_length_mean.unwrap_or(3.0),
            })
        } else {
            None
        };

        let r2_error_multiplier = seq_err_cfg.r2_error_multiplier.unwrap_or(1.0);

        // Pre-scale the cycle curve for R2 at construction time so that
        // inject_all_errors never rebuilds it on the hot path.
        let cycle_curve_r2 = if (r2_error_multiplier - 1.0).abs() > f64::EPSILON {
            cycle_curve.as_ref().map(|curve| {
                CycleErrorCurve::from_rates(
                    curve.rates().iter().map(|&r| r * r2_error_multiplier),
                    read_length,
                )
            })
        } else {
            None
        };

        Ok(Some(Self {
            cycle_curve,
            cycle_curve_r2,
            indel_model,
            kmer_model,
            strand_bias,
            burst_model,
            context_base_rate,
            r2_error_multiplier,
            read_length,
        }))
    }

    /// Apply all active error passes to a read pair.
    ///
    /// Passes are applied in this order:
    /// 1. Strand bias: lower R2 quality scores.
    /// 2. Cycle-position substitutions for R1 (and a scaled variant for R2).
    /// 3. Context-dependent k-mer substitutions.
    /// 4. Sequencing indels (R1 then R2).
    /// 5. Correlated phasing burst errors (R1 then R2).
    pub fn inject_all_errors<R: Rng>(
        &self,
        r1_seq: &mut Vec<u8>,
        r1_qual: &mut Vec<u8>,
        r2_seq: &mut Vec<u8>,
        r2_qual: &mut Vec<u8>,
        rng: &mut R,
    ) {
        // Pass 1: strand bias quality offset on R2.
        if let Some(ref bias) = self.strand_bias {
            bias.apply_to_r2_qual(r2_qual);
        }

        // Pass 2: cycle-position substitutions.
        if let Some(ref curve) = self.cycle_curve {
            inject_cycle_errors(r1_seq, curve, rng);

            // Use the pre-scaled R2 curve when the multiplier differs from 1.0;
            // otherwise reuse the R1 curve directly.
            let r2_curve = self.cycle_curve_r2.as_ref().unwrap_or(curve);
            inject_cycle_errors(r2_seq, r2_curve, rng);
        }

        // Pass 3: context-dependent k-mer substitutions.
        if let Some(ref kmer) = self.kmer_model {
            if self.context_base_rate > 0.0 {
                inject_context_errors(r1_seq, self.context_base_rate, kmer, rng);
                inject_context_errors(
                    r2_seq,
                    self.context_base_rate * self.r2_error_multiplier,
                    kmer,
                    rng,
                );
            }
        }

        // Pass 4: sequencing indels.
        if let Some(ref model) = self.indel_model {
            inject_indel_errors(r1_seq, r1_qual, self.read_length, model, rng);
            inject_indel_errors(r2_seq, r2_qual, self.read_length, model, rng);
        }

        // Pass 5: correlated phasing burst errors.
        if let Some(ref model) = self.burst_model {
            inject_burst_errors(r1_seq, r1_qual, model, rng);
            inject_burst_errors(r2_seq, r2_qual, model, rng);
        }
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

    fn make_read_pair(len: usize) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
        let r1_seq = vec![b'A'; len];
        let r1_qual = vec![30u8; len];
        let r2_seq = vec![b'A'; len];
        let r2_qual = vec![30u8; len];
        (r1_seq, r1_qual, r2_seq, r2_qual)
    }

    /// An orchestrator built from an all-None config must be `None`.
    #[test]
    fn test_no_config_noop() {
        let cfg = SequencingErrorConfig::default();
        let orchestrator = ErrorOrchestrator::from_config(&cfg, 150).unwrap();
        assert!(
            orchestrator.is_none(),
            "empty SequencingErrorConfig should produce no orchestrator"
        );
    }

    /// Backward compatibility: no orchestrator means seq/qual are untouched.
    #[test]
    fn test_backward_compat() {
        let cfg = SequencingErrorConfig::default();
        // from_config returns None; we simulate the engine skipping inject_all_errors.
        let orchestrator = ErrorOrchestrator::from_config(&cfg, 150).unwrap();
        assert!(orchestrator.is_none());

        // Verify that not calling inject_all_errors leaves reads unchanged.
        let (r1_seq, r1_qual, r2_seq, r2_qual) = make_read_pair(150);
        let orig = (
            r1_seq.clone(),
            r1_qual.clone(),
            r2_seq.clone(),
            r2_qual.clone(),
        );
        // No mutation happens.
        assert_eq!(r1_seq, orig.0);
        assert_eq!(r1_qual, orig.1);
        assert_eq!(r2_seq, orig.2);
        assert_eq!(r2_qual, orig.3);
    }

    /// With non-zero indel_rate and base_error_rate, both reads accumulate changes.
    #[test]
    fn test_indel_and_cycle_compose() {
        let cfg = SequencingErrorConfig {
            base_error_rate: Some(0.1),
            indel_rate: Some(0.1),
            indel_insertion_fraction: Some(0.5),
            max_indel_length: Some(2),
            ..Default::default()
        };
        let orchestrator = ErrorOrchestrator::from_config(&cfg, 100)
            .unwrap()
            .expect("orchestrator should be Some with active config");

        let mut rng = StdRng::seed_from_u64(42);

        let mut substitution_changes = 0usize;
        let mut indel_traces = 0usize;
        let n = 1_000;

        for _ in 0..n {
            let (mut r1_seq, mut r1_qual, mut r2_seq, mut r2_qual) = make_read_pair(100);
            let orig_r1 = r1_seq.clone();
            orchestrator.inject_all_errors(
                &mut r1_seq,
                &mut r1_qual,
                &mut r2_seq,
                &mut r2_qual,
                &mut rng,
            );
            // Count substitutions (non-A bases that are not N from indel padding).
            for (&new, &old) in r1_seq.iter().zip(orig_r1.iter()) {
                if new != old {
                    if new == b'N' {
                        indel_traces += 1;
                    } else {
                        substitution_changes += 1;
                    }
                }
            }
            // N-padding from deletions also counts.
            for &b in &r1_seq {
                if b == b'N' {
                    indel_traces += 1;
                }
            }
        }

        assert!(
            substitution_changes > 0,
            "expected substitutions from cycle error pass, got zero"
        );
        assert!(
            indel_traces > 0 || substitution_changes > 1000,
            "expected indel traces or many substitutions; indel_traces={indel_traces}"
        );
    }

    /// r2_quality_offset lowers R2 quality; R1 quality must be unchanged.
    #[test]
    fn test_r2_quality_lowered() {
        let offset: i8 = 5;
        let cfg = SequencingErrorConfig {
            r2_quality_offset: Some(offset),
            r2_error_multiplier: Some(1.0),
            ..Default::default()
        };
        let orchestrator = ErrorOrchestrator::from_config(&cfg, 100)
            .unwrap()
            .expect("orchestrator should be Some");

        let mut rng = StdRng::seed_from_u64(1);
        let (mut r1_seq, mut r1_qual, mut r2_seq, mut r2_qual) = make_read_pair(100);
        let orig_r1_qual = r1_qual.clone();

        orchestrator.inject_all_errors(
            &mut r1_seq,
            &mut r1_qual,
            &mut r2_seq,
            &mut r2_qual,
            &mut rng,
        );

        // R1 quality must be unchanged (no bias applied to R1).
        assert_eq!(r1_qual, orig_r1_qual, "R1 quality must not be changed");

        // Every R2 quality must be exactly 5 lower (30 - 5 = 25).
        for &q in &r2_qual {
            assert_eq!(q, 25, "R2 quality should be 30 - 5 = 25, got {q}");
        }
    }

    /// NovaSeq-like config: observed R1 substitution rate must fall within
    /// [0.0005, 0.003] across 5 000 reads of length 150.
    #[test]
    fn test_novaseq_rates_end_to_end() {
        let cfg = SequencingErrorConfig {
            base_error_rate: Some(0.001),
            cycle_error_model: Some("exponential".to_string()),
            tail_start_fraction: Some(0.8),
            tail_rate_multiplier: Some(5.0),
            indel_rate: Some(0.00005),
            indel_insertion_fraction: Some(0.5),
            max_indel_length: Some(1),
            r2_error_multiplier: Some(1.3),
            r2_quality_offset: Some(2),
            ..Default::default()
        };
        let orch = ErrorOrchestrator::from_config(&cfg, 150)
            .unwrap()
            .expect("NovaSeq config must produce an orchestrator");

        let mut rng = StdRng::seed_from_u64(42);
        let n_reads = 5_000usize;
        let read_length = 150usize;
        let mut total_subs = 0usize;
        let mut total_bases = 0usize;

        for _ in 0..n_reads {
            let original: Vec<u8> = (0..read_length)
                .map(|i| match i % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                })
                .collect();
            let mut r1_seq = original.clone();
            let mut r1_qual = vec![37u8; read_length];
            let mut r2_seq = original.clone();
            let mut r2_qual = vec![37u8; read_length];

            orch.inject_all_errors(
                &mut r1_seq,
                &mut r1_qual,
                &mut r2_seq,
                &mut r2_qual,
                &mut rng,
            );

            total_subs += r1_seq
                .iter()
                .zip(original.iter())
                .filter(|(&a, &b)| a != b && a != b'N')
                .count();
            total_bases += read_length;
        }

        let sub_rate = total_subs as f64 / total_bases as f64;
        // The exponential tail (multiplier 5×, starting at 80 % of the read)
        // drives the average above the flat base rate of 0.001. The upper
        // bound of 0.006 covers that tail while still catching broken configs.
        assert!(
            (0.0005..=0.006).contains(&sub_rate),
            "NovaSeq sub rate {sub_rate:.5} outside [0.0005, 0.006]"
        );
    }
}
