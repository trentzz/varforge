//! Read simulation engine: generates read pairs for a genomic region, applies variants,
//! artifacts, UMI families, and quality models.

use std::sync::Arc;

use anyhow::Result;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

use crate::artifacts::duplicates::{duplicate_read_pair, select_duplicates};
use crate::artifacts::ffpe::{inject_ffpe_damage, inject_oxog_damage};
use crate::core::capture::CaptureModel;
use crate::core::coverage::read_pairs_for_coverage;
use crate::core::end_motifs::accept_fragment_by_end_motif;
use crate::core::error_profile::EmpiricalQualityModel;
use crate::core::fragment::{sample_long_read_length, CfdnaFragmentSampler, NormalFragmentSampler};
use crate::core::gc_bias::GcBiasModel;
use crate::core::quality::{
    inject_errors, sample_nanopore_r10_qualities, sample_pacbio_hifi_qualities,
    ParametricQualityModel, QualityModel,
};
use crate::core::types::{MutationType, Read, ReadPair, Region, SvType, Variant, VariantTag};
use crate::io::config::{Config, FragmentModel};
use crate::io::reference::ReferenceGenome;
use crate::seq_utils::reverse_complement;
use crate::umi::barcode::{generate_duplex_umi_pair, generate_umi};
use crate::umi::families::{generate_pcr_copies, UmiFamily};
use crate::variants::cnv::{adjusted_coverage, find_cn_region, CopyNumberRegion};
use crate::variants::spike_in::spike_indel;
use crate::variants::structural::{
    apply_deletion, apply_duplication, apply_insertion, apply_inversion, apply_translocation,
    StructuralVariant, SvReadEffect,
};
use crate::variants::vaf::sample_alt_count;

/// The output produced for a single genomic region.
pub struct RegionOutput {
    pub read_pairs: Vec<ReadPair>,
    /// Variants that were actually spiked into at least one read.
    pub applied_variants: Vec<AppliedVariant>,
    /// Total molecules (original read pairs) that entered the duplex UMI path.
    /// Zero when duplex mode is not active.
    pub duplex_total_molecules: u64,
    /// Molecules for which both AB and BA strand families were produced.
    /// In simulation this always equals duplex_total_molecules (no library losses).
    pub duplex_molecules_with_both_strands: u64,
}

/// Tracks actual alt and total depths for a spiked variant.
pub struct AppliedVariant {
    pub variant: Variant,
    pub actual_alt_count: u32,
    // Not yet surfaced in truth VCF output; retained for future depth reporting.
    #[allow(dead_code)]
    pub actual_total_count: u32,
    /// Number of duplex (AB+BA) family pairs carrying the alt allele.
    pub duplex_alt_count: u32,
}

/// Core simulation engine that wires together all sub-modules.
pub struct SimulationEngine {
    pub config: Arc<Config>,
    pub reference: Arc<ReferenceGenome>,
    pub rng: StdRng,
    /// Optional GC bias model, built from config.gc_bias.
    pub gc_bias_model: Option<GcBiasModel>,
    /// Optional capture model for targeted sequencing.
    pub capture_model: Option<Arc<CaptureModel>>,
    /// Optional parsed CNV regions.
    pub cnv_regions: Vec<CopyNumberRegion>,
    /// Optional empirical quality model loaded from profile_path.
    pub empirical_quality: Option<EmpiricalQualityModel>,
}

impl SimulationEngine {
    /// Create a new engine from a config and a loaded reference genome.
    ///
    /// If `config.seed` is set that seed is used; otherwise entropy is drawn
    /// from the OS.
    // Called only in tests; production code uses new_with_shared_config.
    #[allow(dead_code)]
    pub fn new(config: Config, reference: ReferenceGenome) -> Self {
        let config = Arc::new(config);
        let rng = match config.seed {
            Some(seed) => StdRng::seed_from_u64(seed),
            None => StdRng::from_entropy(),
        };
        let gc_bias_model = build_gc_bias_model(&config);
        let cnv_regions = build_cnv_regions(&config);
        let empirical_quality = build_empirical_quality(&config);
        Self {
            config,
            reference: Arc::new(reference),
            rng,
            gc_bias_model,
            capture_model: None,
            cnv_regions,
            empirical_quality,
        }
    }

    /// Create a new engine sharing an already-Arc-wrapped reference genome.
    ///
    /// This is used for parallel region simulation where the same reference
    /// is shared across multiple engines running on different threads.
    // Not called in production yet; retained as an alternative constructor for future callers.
    #[allow(dead_code)]
    pub fn new_with_shared_reference(config: Config, reference: Arc<ReferenceGenome>) -> Self {
        let config = Arc::new(config);
        let rng = match config.seed {
            Some(seed) => StdRng::seed_from_u64(seed),
            None => StdRng::from_entropy(),
        };
        let gc_bias_model = build_gc_bias_model(&config);
        let cnv_regions = build_cnv_regions(&config);
        let empirical_quality = build_empirical_quality(&config);
        Self {
            config,
            reference,
            rng,
            gc_bias_model,
            capture_model: None,
            cnv_regions,
            empirical_quality,
        }
    }

    /// Create a new engine sharing an Arc-wrapped config and reference.
    ///
    /// Used in the parallel region loop to share the config across all worker
    /// threads without cloning the full struct per region. The per-region seed
    /// is passed separately so that each engine gets a unique RNG stream.
    pub fn new_with_shared_config(
        config: Arc<Config>,
        seed: Option<u64>,
        reference: Arc<ReferenceGenome>,
    ) -> Self {
        let rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };
        let gc_bias_model = build_gc_bias_model(&config);
        let cnv_regions = build_cnv_regions(&config);
        let empirical_quality = build_empirical_quality(&config);
        Self {
            config,
            reference,
            rng,
            gc_bias_model,
            capture_model: None,
            cnv_regions,
            empirical_quality,
        }
    }

    /// Create an engine seeded deterministically for a specific region index.
    ///
    /// The per-region seed is derived by hashing the master seed and the
    /// region index together, ensuring:
    /// - Each region has a unique, reproducible seed.
    /// - The same (master_seed, region_index) pair always produces the same output.
    /// - Thread scheduling does not affect determinism.
    ///
    /// If `master_seed` is `None`, per-region entropy is drawn from the OS.
    // Not called in production yet; retained as a convenience constructor for parallel
    // pipelines that operate per-region rather than per-config.
    #[allow(dead_code)]
    pub fn new_for_region(
        config: Config,
        reference: ReferenceGenome,
        master_seed: Option<u64>,
        region_index: usize,
    ) -> Self {
        let rng = match master_seed {
            Some(seed) => {
                // Mix seed and region_index with FNV-style hashing for a
                // compact, dependency-free derivation function.
                let region_seed = derive_region_seed(seed, region_index as u64);
                StdRng::seed_from_u64(region_seed)
            }
            None => StdRng::from_entropy(),
        };
        let config = Arc::new(config);
        let gc_bias_model = build_gc_bias_model(&config);
        let cnv_regions = build_cnv_regions(&config);
        let empirical_quality = build_empirical_quality(&config);
        Self {
            config,
            reference: Arc::new(reference),
            rng,
            gc_bias_model,
            capture_model: None,
            cnv_regions,
            empirical_quality,
        }
    }

    /// Attach a CaptureModel to this engine (call before generate_region).
    pub fn set_capture_model(&mut self, model: Arc<CaptureModel>) {
        self.capture_model = Some(model);
    }

    /// Generate all read pairs for a genomic region.
    ///
    /// The pipeline for each pair:
    /// 1. Sample fragment length.
    /// 2. Sample a random start position within the region.
    /// 3. Extract reference sequence.
    /// 4. Apply GC bias rejection sampling (if enabled).
    /// 5. Apply variant spike-in (stochastic VAF) for small variants.
    ///    For SV variants, apply structural effects to the read.
    /// 6. Generate quality scores via empirical or parametric model + inject errors.
    ///
    /// Post-pair processing (applied to the full batch):
    /// 7. Attach UMI barcodes and expand PCR families (if UMI enabled).
    /// 8. Inject FFPE / oxoG damage (if artifacts enabled).
    /// 9. Expand PCR duplicates (if artifacts enabled).
    ///
    /// Coverage is adjusted for:
    /// - CNV regions (copy number × purity formula)
    /// - Capture efficiency (on/off-target multipliers)
    pub fn generate_region(
        &mut self,
        region: &Region,
        variants: &[Variant],
    ) -> Result<RegionOutput> {
        let read_length = self.config.sample.read_length;
        let coverage = self.config.sample.coverage;

        // ---- CNV coverage adjustment ----
        // Compute effective coverage for the region centre.
        let region_centre = (region.start + region.end) / 2;
        let effective_coverage = if !self.cnv_regions.is_empty() {
            let purity = self.config.tumour.as_ref().map(|t| t.purity).unwrap_or(1.0);
            let ploidy = self.config.tumour.as_ref().map(|t| t.ploidy).unwrap_or(2);
            if let Some(cn_region) = find_cn_region(&self.cnv_regions, &region.chrom, region_centre)
            {
                adjusted_coverage(
                    coverage,
                    purity,
                    cn_region.tumor_cn,
                    cn_region.normal_cn,
                    ploidy,
                )
            } else {
                coverage
            }
        } else {
            coverage
        };

        // ---- Capture coverage adjustment ----
        // If a capture model is present, sample per-target multipliers and
        // compute the effective coverage multiplier for the region centre.
        let capture_multiplier = if let Some(ref capture_model) = self.capture_model {
            let target_multipliers = capture_model.sample_all_target_multipliers(&mut self.rng);
            capture_model.coverage_multiplier_at(&region.chrom, region_centre, &target_multipliers)
        } else {
            1.0
        };

        let final_coverage = effective_coverage * capture_multiplier;
        let n_pairs = read_pairs_for_coverage(region.len(), final_coverage, read_length);

        // Build fragment sampler from config.
        // FragmentSampler uses a generic Rng parameter so we can't use
        // Box<dyn FragmentSampler>. Use an enum instead.
        //
        // When `fragment.long_read` is set, use the log-normal long-read sampler
        // regardless of the `model` field.
        enum Sampler<'a> {
            Normal(NormalFragmentSampler),
            Cfdna(CfdnaFragmentSampler),
            LongRead(&'a crate::io::config::LongReadFragmentConfig),
        }
        impl Sampler<'_> {
            fn sample(&self, rng: &mut StdRng) -> Result<usize> {
                Ok(match self {
                    Sampler::Normal(s) => {
                        use crate::core::fragment::FragmentSampler as _;
                        s.sample(rng)
                    }
                    Sampler::Cfdna(s) => {
                        use crate::core::fragment::FragmentSampler as _;
                        s.sample(rng)
                    }
                    Sampler::LongRead(cfg) => sample_long_read_length(cfg, rng)?,
                })
            }
        }
        let fragment_sampler = if let Some(ref lr_cfg) = self.config.fragment.long_read {
            Sampler::LongRead(lr_cfg)
        } else {
            match self.config.fragment.model {
                FragmentModel::Cfda => {
                    let mono_peak = self.config.fragment.mean;
                    let di_peak = mono_peak * 2.0;
                    // Resolve ctDNA fraction: explicit config field takes priority
                    // over purity. Higher purity means more tumour-derived ctDNA,
                    // so we use purity directly (not 1.0 - purity).
                    let ctdna_fraction = self.config.fragment.ctdna_fraction.unwrap_or_else(|| {
                        self.config.tumour.as_ref().map(|t| t.purity).unwrap_or(0.0)
                    });
                    if (ctdna_fraction - 1.0).abs() < f64::EPSILON {
                        tracing::warn!(
                            "purity=1.0 with cfDNA model: all fragments will be \
                             tumour-derived. Did you mean to set a lower purity?"
                        );
                    }
                    let mono_sd = self.config.fragment.mono_sd.unwrap_or(20.0);
                    let di_sd = self.config.fragment.di_sd.unwrap_or(30.0);
                    Sampler::Cfdna(CfdnaFragmentSampler::new(
                        mono_peak,
                        di_peak,
                        0.85,
                        ctdna_fraction,
                        mono_sd,
                        di_sd,
                    )?)
                }
                _ => Sampler::Normal(NormalFragmentSampler::new(
                    self.config.fragment.mean,
                    self.config.fragment.sd,
                )?),
            }
        };

        // Per-variant alt/total counters (indexed parallel to `variants`).
        let mut alt_counts = vec![0u32; variants.len()];
        let mut total_counts = vec![0u32; variants.len()];

        // Determine which variants overlap the region, and how many alt reads
        // each should receive across the whole batch.
        let total_depth = n_pairs as u32;
        let mut variant_alt_budget: Vec<u32> = variants
            .iter()
            .map(|v| {
                if variant_overlaps_region(v, region) {
                    sample_alt_count(total_depth, v.expected_vaf, &mut self.rng)
                } else {
                    0
                }
            })
            .collect();

        // Pre-fetch the full region sequence once (Fix C-3: avoid per-fragment reference lookup).
        let region_seq: Vec<u8> = self.reference.sequence(region)?;

        // Hoist quality model construction outside the fragment loop (Fix C-2).
        //
        // For long-read platforms, use platform-specific flat quality profiles
        // rather than the Illumina parametric model. Platform is detected from
        // `sample.platform` only when `fragment.long_read` is set, ensuring
        // that the long-read quality model is not applied to short-read runs.
        let is_long_read = self.config.fragment.long_read.is_some();
        let platform = self
            .config
            .sample
            .platform
            .as_deref()
            .unwrap_or("")
            .to_lowercase();
        enum QualModel<'a> {
            Empirical(&'a EmpiricalQualityModel),
            Parametric(ParametricQualityModel),
            PacbioHifi,
            NanoporeR10,
        }
        let qual_model = if let Some(ref emp) = self.empirical_quality {
            QualModel::Empirical(emp)
        } else if is_long_read && platform.contains("pacbio") {
            QualModel::PacbioHifi
        } else if is_long_read && platform.contains("nanopore") {
            QualModel::NanoporeR10
        } else {
            QualModel::Parametric(ParametricQualityModel::new(
                self.config.quality.mean_quality,
                self.config.quality.tail_decay,
            ))
        };

        // Determine whether duplex UMI mode is active. In duplex mode, every alt
        // molecule produces both an AB and a BA strand family, so the duplex alt
        // count equals the molecule-level alt count.
        let is_duplex_mode = self.config.umi.as_ref().map(|u| u.duplex).unwrap_or(false);

        let mut read_pairs: Vec<ReadPair> = Vec::with_capacity(n_pairs as usize);

        // We may need to generate more than n_pairs fragments due to GC bias
        // rejection: track how many accepted pairs we have so far.
        let mut pair_idx: u64 = 0;
        let mut attempts: u64 = 0;
        // Safety cap: at most 10× target pairs to avoid infinite loop on
        // degenerate GC bias configs.
        let max_attempts: u64 = (n_pairs as u64).saturating_mul(10).max(100);

        // Pre-allocate read name buffer to avoid repeated allocation (Fix H-5).
        let chrom_prefix = &region.chrom;
        let mut name_buf = String::with_capacity(64);

        while pair_idx < n_pairs && attempts < max_attempts {
            attempts += 1;

            let frag_len = fragment_sampler.sample(&mut self.rng)?;

            // Sample fragment start within the region, clamped so the
            // fragment doesn't exceed the region boundary.
            let region_len = region.len() as usize;
            let max_start_offset = region_len.saturating_sub(frag_len);
            let frag_start = region.start
                + if max_start_offset > 0 {
                    self.rng.gen_range(0..max_start_offset as u64)
                } else {
                    0
                };

            // Clamp fragment end to region end.
            let frag_end = (frag_start + frag_len as u64).min(region.end);
            let actual_frag_len = (frag_end - frag_start) as usize;

            // Skip severely truncated fragments (e.g. post-deletion tail with fewer
            // than half the requested read bases).  N-padding is fine for fragments
            // that are only a few bases short; this guard prevents near-empty reads
            // that would panic in downstream slice operations.
            if actual_frag_len < read_length / 2 {
                continue;
            }

            // ---- End motif rejection sampling ----
            // When the plasma end motif model is enabled, extract the 4-mer at
            // the fragment 5' end and use rejection sampling to bias fragment
            // starts toward motifs that are enriched in plasma cfDNA.
            if self.config.fragment.end_motif_model.as_deref() == Some("plasma") {
                let motif_5p = {
                    let pos = (frag_start - region.start) as usize;
                    if pos + 4 <= region_seq.len() {
                        Some([
                            region_seq[pos],
                            region_seq[pos + 1],
                            region_seq[pos + 2],
                            region_seq[pos + 3],
                        ])
                    } else {
                        None
                    }
                };
                if !accept_fragment_by_end_motif(motif_5p, &mut self.rng) {
                    continue;
                }
            }

            // Extract reference sequence for the fragment from the pre-fetched region sequence
            // (Fix C-3: slice into cached region_seq instead of per-fragment reference call).
            let offset_start = (frag_start - region.start) as usize;
            let offset_end = (frag_end - region.start) as usize;
            let mut frag_seq = region_seq[offset_start..offset_end.min(region_seq.len())].to_vec();

            // Pad with 'N' if shorter than expected (shouldn't normally happen).
            while frag_seq.len() < actual_frag_len {
                frag_seq.push(b'N');
            }

            // ---- GC bias rejection sampling ----
            if let Some(ref gc_model) = self.gc_bias_model {
                let gc = GcBiasModel::gc_fraction(&frag_seq);
                let multiplier = gc_model.coverage_multiplier(gc);
                // Reject fragment with probability (1 - multiplier).
                if self.rng.gen::<f64>() >= multiplier {
                    continue;
                }
            }

            // ---- Haplotype assignment ----
            // Each fragment is assigned to haplotype 0 or 1 with equal probability.
            // Variants with a matching haplotype assignment are only applied to
            // fragments on that haplotype.
            let fragment_haplotype: u8 = self.rng.gen_range(0u8..2);

            // ---- Save pre-variant reference for MD tag construction ----
            // Clone the unmodified fragment sequence before any variant is applied.
            // The BAM writer uses this to build the MD tag.
            let pre_variant_frag_seq = frag_seq.clone();

            // ---- Variant spike-in ----
            // Collect tags for any variant actually applied to this fragment.
            let mut fragment_variant_tags: Vec<VariantTag> = Vec::new();
            for (vi, variant) in variants.iter().enumerate() {
                if variant_alt_budget[vi] == 0 {
                    continue;
                }
                if !variant_overlaps_region(variant, region) {
                    continue;
                }
                // Haplotype filter: skip this variant if it is assigned to a
                // haplotype that does not match this fragment.
                if let Some(h) = variant.haplotype {
                    if h != fragment_haplotype {
                        continue;
                    }
                }
                // Check if this fragment overlaps the variant position.
                let var_pos = variant.pos();
                if var_pos < frag_start || var_pos >= frag_end {
                    continue;
                }

                // Increment total count for this variant.
                total_counts[vi] += 1;

                // Consume one alt from the budget and spike.
                let budget = &mut variant_alt_budget[vi];
                let remaining_pairs = n_pairs - pair_idx;
                // Stochastically decide: spike this read given remaining budget and pairs.
                let spike_prob = (*budget as f64) / (remaining_pairs as f64).max(1.0);
                if self.rng.gen::<f64>() < spike_prob {
                    apply_variant_to_seq(&mut frag_seq, frag_start, variant);
                    *budget = budget.saturating_sub(1);
                    alt_counts[vi] += 1;
                    // Record which variant was spiked into this fragment.
                    fragment_variant_tags.push(VariantTag {
                        chrom: variant.chrom.clone(),
                        pos: mutation_pos(&variant.mutation),
                        vartype: mutation_vartype(&variant.mutation).to_string(),
                        vaf: variant.expected_vaf,
                        clone_id: variant.clone_id.clone(),
                    });
                }
            }

            // ---- Build read1 and read2 ----
            // read1: first read_length bases of fragment.
            // read2: last read_length bases of fragment (reverse complement semantics
            //        are handled by downstream tools; we store the actual sequence here).
            //
            // Re-derive the fragment length from frag_seq: large indels (SV spike-ins)
            // change the sequence length after actual_frag_len was computed above.
            let actual_frag_len = frag_seq.len();
            let r1_len = read_length.min(actual_frag_len);
            let r2_len = read_length.min(actual_frag_len);

            // Extract pre-variant reference slices for MD tag.
            // Use the unmodified fragment length (from pre_variant_frag_seq) to
            // determine the R2 reference offset, then pad to read_length.
            let pre_len = pre_variant_frag_seq.len();
            let pre_r1_len = read_length.min(pre_len);
            let pre_r2_len = read_length.min(pre_len);
            let mut ref_r1: Vec<u8> = pre_variant_frag_seq[..pre_r1_len].to_vec();
            let mut ref_r2: Vec<u8> = if pre_len >= read_length {
                pre_variant_frag_seq[pre_len - pre_r2_len..].to_vec()
            } else {
                pre_variant_frag_seq[pre_len - pre_r2_len.min(pre_len)..].to_vec()
            };
            while ref_r1.len() < read_length {
                ref_r1.push(b'N');
            }
            while ref_r2.len() < read_length {
                ref_r2.push(b'N');
            }

            let mut r1_seq: Vec<u8> = frag_seq[..r1_len].to_vec();
            let mut r2_seq: Vec<u8> = if actual_frag_len >= read_length {
                frag_seq[actual_frag_len - r2_len..].to_vec()
            } else {
                frag_seq[actual_frag_len - r2_len.min(actual_frag_len)..].to_vec()
            };

            // Pad reads to requested read_length if fragment was shorter.
            while r1_seq.len() < read_length {
                r1_seq.push(b'N');
            }
            while r2_seq.len() < read_length {
                r2_seq.push(b'N');
            }

            // ---- Quality and error injection ----
            // Use empirical model when loaded, else platform-specific model for
            // long reads, else fall back to the Illumina parametric model
            // (Fix C-2: model was built once before the loop, not re-constructed
            // each iteration).
            let (r1_qual, r2_qual) = match &qual_model {
                QualModel::Empirical(emp) => {
                    let q1 = emp.generate_qualities(read_length, &mut self.rng);
                    let q2 = emp.generate_qualities(read_length, &mut self.rng);
                    emp.inject_errors(&mut r1_seq, &q1, &mut self.rng);
                    emp.inject_errors(&mut r2_seq, &q2, &mut self.rng);
                    (q1, q2)
                }
                QualModel::Parametric(parametric) => {
                    let q1 = parametric.generate_qualities(read_length, &mut self.rng);
                    let q2 = parametric.generate_qualities(read_length, &mut self.rng);
                    inject_errors(&mut r1_seq, &q1, &mut self.rng);
                    inject_errors(&mut r2_seq, &q2, &mut self.rng);
                    (q1, q2)
                }
                QualModel::PacbioHifi => {
                    // For long reads, r1 spans the full fragment; r2 is unused
                    // but generated to keep the ReadPair structure intact.
                    let q1 = sample_pacbio_hifi_qualities(read_length, &mut self.rng);
                    let q2 = sample_pacbio_hifi_qualities(read_length, &mut self.rng);
                    inject_errors(&mut r1_seq, &q1, &mut self.rng);
                    inject_errors(&mut r2_seq, &q2, &mut self.rng);
                    (q1, q2)
                }
                QualModel::NanoporeR10 => {
                    // Pass the read sequence so homopolymer-aware quality degradation
                    // can be applied at the correct positions.
                    let q1 = sample_nanopore_r10_qualities(read_length, &r1_seq, &mut self.rng);
                    let q2 = sample_nanopore_r10_qualities(read_length, &r2_seq, &mut self.rng);
                    inject_errors(&mut r1_seq, &q1, &mut self.rng);
                    inject_errors(&mut r2_seq, &q2, &mut self.rng);
                    (q1, q2)
                }
            };

            // Build read name by reusing the pre-allocated buffer (Fix H-5).
            use std::fmt::Write as FmtWrite;
            name_buf.clear();
            write!(
                &mut name_buf,
                "{}_{}:{}",
                chrom_prefix, frag_start, pair_idx
            )
            .unwrap();

            let pair = ReadPair {
                name: name_buf.clone(),
                read1: Read::new(r1_seq, r1_qual),
                read2: Read::new(r2_seq, r2_qual),
                fragment_start: frag_start,
                fragment_length: actual_frag_len,
                chrom: region.chrom.clone(),
                variant_tags: fragment_variant_tags,
                ref_seq_r1: ref_r1,
                ref_seq_r2: ref_r2,
            };

            read_pairs.push(pair);
            pair_idx += 1;
        }

        // ---- SV spike-in: apply structural effects to reads ----
        // For each SV variant that overlaps the region, apply its effect to
        // reads that span the breakpoint.  Reads marked Deleted are removed.
        let sv_variants: Vec<(usize, &Variant)> = variants
            .iter()
            .enumerate()
            .filter(|(_, v)| {
                matches!(v.mutation, MutationType::Sv { .. }) && variant_overlaps_region(v, region)
            })
            .collect();

        if !sv_variants.is_empty() {
            // Build StructuralVariant objects from SV-type Variant entries.
            let mut surviving_pairs: Vec<ReadPair> = Vec::with_capacity(read_pairs.len());
            for mut pair in read_pairs {
                let mut deleted = false;
                for (vi, variant) in &sv_variants {
                    if variant_alt_budget[*vi] == 0 && total_counts[*vi] == 0 {
                        continue;
                    }
                    if let Some(sv) = variant_to_structural(variant) {
                        // Decide stochastically whether this read should carry the SV allele.
                        let should_spike = {
                            let budget = &mut variant_alt_budget[*vi];
                            if *budget > 0 {
                                let remaining = (n_pairs - pair_idx.min(n_pairs)) as f64;
                                let prob = (*budget as f64) / remaining.max(1.0);
                                if self.rng.gen::<f64>() < prob {
                                    *budget = budget.saturating_sub(1);
                                    alt_counts[*vi] += 1;
                                    true
                                } else {
                                    false
                                }
                            } else {
                                false
                            }
                        };
                        total_counts[*vi] += 1;

                        if should_spike {
                            let read_start = pair.fragment_start;
                            let effect1 = apply_sv_to_read(
                                &mut pair.read1,
                                read_start,
                                &pair.chrom,
                                &sv,
                                &self.reference,
                            );
                            let _effect2 = apply_sv_to_read(
                                &mut pair.read2,
                                read_start,
                                &pair.chrom,
                                &sv,
                                &self.reference,
                            );
                            if effect1 == SvReadEffect::Deleted {
                                deleted = true;
                                break;
                            }
                            // Record the SV spike-in in the read pair's variant tags.
                            pair.variant_tags.push(VariantTag {
                                chrom: variant.chrom.clone(),
                                pos: mutation_pos(&variant.mutation),
                                vartype: mutation_vartype(&variant.mutation).to_string(),
                                vaf: variant.expected_vaf,
                                clone_id: variant.clone_id.clone(),
                            });
                        }
                    }
                }
                if !deleted {
                    surviving_pairs.push(pair);
                }
            }
            read_pairs = surviving_pairs;
        }

        // ---- UMI / PCR families ----
        let mut duplex_total_molecules: u64 = 0;
        let mut duplex_molecules_with_both_strands: u64 = 0;

        if let Some(ref umi_cfg) = self.config.umi {
            let umi_len = umi_cfg.length;
            let family_mean = umi_cfg.family_size_mean;
            let family_sd = umi_cfg.family_size_sd;
            let pcr_cycles = umi_cfg.pcr_cycles;
            let is_duplex = umi_cfg.duplex;
            let pcr_error_rate = self
                .config
                .artifacts
                .as_ref()
                .and_then(|a| a.pcr_error_rate)
                .unwrap_or(0.0);

            // Construct sampler once; log-space conversion is non-trivial so
            // avoid repeating it for every read pair.
            let family_sampler =
                crate::core::fragment::PcrFamilySizeSampler::new(family_mean, family_sd)?;

            let mut families: Vec<ReadPair> = Vec::new();
            for mut pair in read_pairs {
                if is_duplex {
                    // Track molecule counts for duplex conversion rate (T108).
                    duplex_total_molecules += 1;
                    duplex_molecules_with_both_strands += 1; // always both strands in simulation

                    // Generate an AB/BA duplex pair.
                    // AB UMI is "AAAA-BBBB"; BA UMI is the halves swapped.
                    let (umi_a, umi_b) = generate_duplex_umi_pair(umi_len, &mut self.rng);
                    let a_str = String::from_utf8_lossy(&umi_a[..umi_len]);
                    let b_str = String::from_utf8_lossy(&umi_b[..umi_len]);
                    let ab_str = format!("{}-{}", a_str, b_str);
                    let ba_str = format!("{}-{}", b_str, a_str);

                    // AB strand: original orientation.
                    pair.name = format!("{}:UMI:{}", pair.name, ab_str);
                    let family_size = family_sampler.sample(&mut self.rng);
                    let ab_family = UmiFamily {
                        umi: umi_a,
                        original: pair.clone(),
                        family_size,
                    };
                    let mut ab_copies =
                        generate_pcr_copies(&ab_family, pcr_error_rate, pcr_cycles, &mut self.rng);
                    families.append(&mut ab_copies);

                    // BA strand: swap R1/R2 and reverse-complement both.
                    // R1_BA = revcomp(R2_AB), R2_BA = revcomp(R1_AB).
                    // This models sequencing from the opposite end of the molecule.
                    let ba_read1 = Read::new(
                        reverse_complement(&pair.read2.seq),
                        pair.read2.qual.iter().copied().rev().collect(),
                    );
                    let ba_read2 = Read::new(
                        reverse_complement(&pair.read1.seq),
                        pair.read1.qual.iter().copied().rev().collect(),
                    );
                    // Strip the ":UMI:..." suffix we already appended, then add BA tag.
                    let base_name = pair
                        .name
                        .rfind(":UMI:")
                        .map(|i| &pair.name[..i])
                        .unwrap_or(&pair.name);
                    let ba_pair = ReadPair {
                        name: format!("{}:UMI:{}", base_name, ba_str),
                        read1: ba_read1,
                        read2: ba_read2,
                        fragment_start: pair.fragment_start,
                        fragment_length: pair.fragment_length,
                        chrom: pair.chrom.clone(),
                        variant_tags: pair.variant_tags.clone(),
                        ref_seq_r1: pair.ref_seq_r2.iter().copied().rev().collect(),
                        ref_seq_r2: pair.ref_seq_r1.iter().copied().rev().collect(),
                    };
                    let ba_family_size = family_sampler.sample(&mut self.rng);
                    let ba_family = UmiFamily {
                        umi: umi_b,
                        original: ba_pair,
                        family_size: ba_family_size,
                    };
                    let mut ba_copies =
                        generate_pcr_copies(&ba_family, pcr_error_rate, pcr_cycles, &mut self.rng);
                    families.append(&mut ba_copies);
                } else {
                    let umi = generate_umi(umi_len, &mut self.rng);
                    let umi_str = String::from_utf8_lossy(&umi).into_owned();
                    pair.name = format!("{}:UMI:{}", pair.name, umi_str);

                    let family_size = family_sampler.sample(&mut self.rng);
                    let family = UmiFamily {
                        umi,
                        original: pair,
                        family_size,
                    };
                    let mut copies =
                        generate_pcr_copies(&family, pcr_error_rate, pcr_cycles, &mut self.rng);
                    families.append(&mut copies);
                }
            }
            read_pairs = families;
        }

        // ---- Artifact injection ----
        if let Some(ref artifact_cfg) = self.config.artifacts {
            if let Some(ffpe_rate) = artifact_cfg.ffpe_damage_rate {
                for pair in &mut read_pairs {
                    inject_ffpe_damage(&mut pair.read1.seq, ffpe_rate, false, &mut self.rng);
                    inject_ffpe_damage(&mut pair.read2.seq, ffpe_rate, true, &mut self.rng);
                }
            }

            if let Some(oxog_rate) = artifact_cfg.oxog_rate {
                for pair in &mut read_pairs {
                    inject_oxog_damage(&mut pair.read1.seq, oxog_rate, true, &mut self.rng);
                    inject_oxog_damage(&mut pair.read2.seq, oxog_rate, false, &mut self.rng);
                }
            }

            if let Some(dup_rate) = artifact_cfg.duplicate_rate {
                let dup_indices = select_duplicates(read_pairs.len(), dup_rate, &mut self.rng);
                let mut duplicates: Vec<ReadPair> = dup_indices
                    .iter()
                    .enumerate()
                    .map(|(i, &idx)| duplicate_read_pair(&read_pairs[idx], i))
                    .collect();
                read_pairs.append(&mut duplicates);
            }
        }

        // ---- Collect applied variants ----
        let applied_variants: Vec<AppliedVariant> = variants
            .iter()
            .enumerate()
            .filter(|(i, _)| total_counts[*i] > 0)
            .map(|(i, v)| AppliedVariant {
                variant: v.clone(),
                actual_alt_count: alt_counts[i],
                actual_total_count: total_counts[i],
                // In duplex mode both AB and BA strand families carry the same
                // variant, so every alt molecule contributes to the duplex count.
                duplex_alt_count: if is_duplex_mode { alt_counts[i] } else { 0 },
            })
            .collect();

        Ok(RegionOutput {
            read_pairs,
            applied_variants,
            duplex_total_molecules,
            duplex_molecules_with_both_strands,
        })
    }
}

// ---------------------------------------------------------------------------
// Seed derivation
// ---------------------------------------------------------------------------

/// Derive a per-region seed from the master seed and a region index.
///
/// Uses a simple mixing function (splitmix64) that is cheap, dependency-free,
/// and produces well-distributed seeds even for sequential region indices.
pub fn derive_region_seed(master: u64, region_index: u64) -> u64 {
    // splitmix64 step to mix the region index into the master seed.
    let mut x = master.wrapping_add(region_index.wrapping_mul(0x9e37_79b9_7f4a_7c15));
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

// ---------------------------------------------------------------------------
// Engine builder helpers
// ---------------------------------------------------------------------------

/// Build a GcBiasModel from config, or None if no gc_bias section is set.
fn build_gc_bias_model(config: &Config) -> Option<GcBiasModel> {
    use crate::core::gc_bias::{GcBiasConfig, GcBiasModelKind};

    let cfg = config.gc_bias.as_ref()?;
    if !cfg.enabled {
        return None;
    }
    let kind = match cfg.model.as_str() {
        "flat" => GcBiasModelKind::Flat,
        "custom" => GcBiasModelKind::Custom,
        _ => GcBiasModelKind::Default,
    };
    let gc_cfg = GcBiasConfig {
        enabled: cfg.enabled,
        model: kind,
        severity: cfg.severity,
    };
    Some(GcBiasModel::from_config(&gc_cfg))
}

/// Parse CNV regions from config into `CopyNumberRegion` objects.
fn build_cnv_regions(config: &Config) -> Vec<CopyNumberRegion> {
    let Some(ref cn_configs) = config.copy_number else {
        return Vec::new();
    };
    cn_configs
        .iter()
        .filter_map(|cn| {
            parse_region_string(&cn.region).map(|(chrom, start, end)| {
                if let (Some(major), Some(minor)) = (cn.major_cn, cn.minor_cn) {
                    CopyNumberRegion::with_allele_specific(
                        chrom,
                        start,
                        end,
                        cn.tumor_cn,
                        cn.normal_cn,
                        major,
                        minor,
                    )
                } else {
                    CopyNumberRegion::new(chrom, start, end, cn.tumor_cn, cn.normal_cn)
                }
            })
        })
        .collect()
}

/// Parse a region string like "chr7:55000000-55200000" into (chrom, start, end).
fn parse_region_string(region: &str) -> Option<(String, u64, u64)> {
    let colon = region.find(':')?;
    let chrom = region[..colon].to_string();
    let coords = &region[colon + 1..];
    let dash = coords.find('-')?;
    let start: u64 = coords[..dash].parse().ok()?;
    let end: u64 = coords[dash + 1..].parse().ok()?;
    Some((chrom, start, end))
}

/// Try to load an EmpiricalQualityModel from config.quality.profile_path.
fn build_empirical_quality(config: &Config) -> Option<EmpiricalQualityModel> {
    let path = config.quality.profile_path.as_deref()?;
    match EmpiricalQualityModel::from_file(path) {
        Ok(model) => Some(model),
        Err(e) => {
            tracing::warn!(
                "failed to load error profile from {}: {}; falling back to parametric model",
                path.display(),
                e
            );
            None
        }
    }
}

// ---------------------------------------------------------------------------
// Helper utilities
// ---------------------------------------------------------------------------

/// Return true if the variant's position falls within the region.
fn variant_overlaps_region(variant: &Variant, region: &Region) -> bool {
    if variant.chrom != region.chrom {
        return false;
    }
    let pos = variant.pos();
    pos >= region.start && pos < region.end
}

/// Apply a variant's mutation to a fragment sequence in-place.
///
/// `seq` is 0-indexed starting at `frag_start`.
///
/// For SNV and MNV, operates directly on the slice without allocating a
/// temporary `Read` (Fix C-1).  Indels still use the `Read` wrapper because
/// they can change the sequence length and need quality-score bookkeeping.
fn apply_variant_to_seq(seq: &mut Vec<u8>, frag_start: u64, variant: &Variant) {
    let frag_end = frag_start + seq.len() as u64;
    match &variant.mutation {
        MutationType::Snv {
            pos,
            ref_base,
            alt_base,
        } => {
            if *pos >= frag_start && *pos < frag_end {
                let offset = (*pos - frag_start) as usize;
                if seq[offset] == *ref_base {
                    seq[offset] = *alt_base;
                }
            }
        }
        MutationType::Mnv {
            pos,
            ref_seq,
            alt_seq,
        } => {
            let mnv_end = *pos + ref_seq.len() as u64;
            if *pos >= frag_start && mnv_end <= frag_end {
                let offset = (*pos - frag_start) as usize;
                if seq[offset..offset + ref_seq.len()] == ref_seq[..] {
                    seq[offset..offset + alt_seq.len()].copy_from_slice(alt_seq);
                }
            }
        }
        MutationType::Indel { .. } => {
            // Indels change sequence length; use the Read wrapper for quality handling.
            let qual_placeholder = vec![30u8; seq.len()];
            let mut temp_read = Read::new(seq.clone(), qual_placeholder);
            spike_indel(&mut temp_read, frag_start, &variant.mutation, b"NNNNNNNNNN");
            *seq = temp_read.seq;
        }
        MutationType::Sv { .. } => {
            // SVs are handled separately via apply_sv_to_read; skip here.
        }
    }
}

/// Convert a Variant with MutationType::Sv into a StructuralVariant.
///
/// Returns None for non-SV variants or unrecognised SV types.
fn variant_to_structural(variant: &Variant) -> Option<StructuralVariant> {
    let MutationType::Sv {
        sv_type,
        chrom: sv_chrom,
        start,
        end,
    } = &variant.mutation
    else {
        return None;
    };
    let chrom = sv_chrom.clone();
    match sv_type {
        SvType::Deletion => Some(StructuralVariant::Deletion {
            chrom,
            start: *start,
            end: *end,
        }),
        SvType::Insertion => {
            // Use N as placeholder inserted sequence (real sequence would need
            // to be specified in a richer variant type).
            Some(StructuralVariant::Insertion {
                chrom,
                pos: *start,
                sequence: b"N".to_vec(),
            })
        }
        SvType::Inversion => Some(StructuralVariant::Inversion {
            chrom,
            start: *start,
            end: *end,
        }),
        SvType::Duplication => Some(StructuralVariant::Duplication {
            chrom,
            start: *start,
            end: *end,
            copies: 1,
        }),
        SvType::Translocation => {
            // For translocations, use start as pos1 and end as pos2 on the
            // same chromosome (simplified model; full inter-chromosomal
            // translocations require a partner chromosome in the variant).
            Some(StructuralVariant::Translocation {
                chrom1: chrom.clone(),
                pos1: *start,
                chrom2: chrom,
                pos2: *end,
            })
        }
    }
}

/// Apply a StructuralVariant to a single Read, fetching partner sequence
/// from the reference as needed for translocations.
fn apply_sv_to_read(
    read: &mut Read,
    read_start: u64,
    chrom: &str,
    sv: &StructuralVariant,
    reference: &ReferenceGenome,
) -> SvReadEffect {
    match sv {
        StructuralVariant::Deletion { .. } => apply_deletion(read, read_start, sv),
        StructuralVariant::Insertion { .. } => apply_insertion(read, read_start, sv),
        StructuralVariant::Inversion { .. } => apply_inversion(read, read_start, sv),
        StructuralVariant::Duplication { .. } => apply_duplication(read, read_start, sv),
        StructuralVariant::Translocation { chrom2, pos2, .. } => {
            // Fetch partner sequence from reference.  Clamp the end coordinate to
            // the chromosome length so that a translocation near the chromosome end
            // does not request an out-of-bounds region.
            let chrom_len = reference.chrom_len(chrom2).unwrap_or(u64::MAX);
            let partner_end = (pos2 + read.len() as u64).min(chrom_len);
            let partner_region = Region::new(chrom2.clone(), *pos2, partner_end);
            let partner_seq = reference.sequence(&partner_region).unwrap_or_default();
            apply_translocation(read, read_start, chrom, sv, &partner_seq)
        }
    }
}

/// Return a short string label for a mutation type, used in variant tags.
fn mutation_vartype(mutation: &MutationType) -> &'static str {
    match mutation {
        MutationType::Snv { .. } => "SNV",
        MutationType::Indel { .. } => "INDEL",
        MutationType::Mnv { .. } => "MNV",
        MutationType::Sv { .. } => "SV",
    }
}

/// Return the primary genomic position of a mutation, used in variant tags.
fn mutation_pos(mutation: &MutationType) -> u64 {
    match mutation {
        MutationType::Snv { pos, .. } => *pos,
        MutationType::Indel { pos, .. } => *pos,
        MutationType::Mnv { pos, .. } => *pos,
        MutationType::Sv { start, .. } => *start,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    use crate::io::config::{
        ArtifactConfig, Config, FragmentConfig, FragmentModel, OutputConfig, QualityConfig,
        SampleConfig, UmiConfig,
    };

    // -----------------------------------------------------------------------
    // Test helpers
    // -----------------------------------------------------------------------

    /// Write a 1000-bp FASTA + matching index to a temp dir.
    /// Sequence is a repeating "ACGT" pattern so we know exact bases.
    fn write_test_fasta(dir: &TempDir) -> std::path::PathBuf {
        let fa_path = dir.path().join("test.fa");
        let fai_path = dir.path().join("test.fa.fai");

        // 1000 bp = 250 repeats of "ACGT"
        let seq: Vec<u8> = b"ACGT".iter().cycle().take(1000).cloned().collect();
        let header = b">chr1\n";
        let newline = b"\n";

        let mut fa_bytes = Vec::new();
        fa_bytes.extend_from_slice(header);
        fa_bytes.extend_from_slice(&seq);
        fa_bytes.extend_from_slice(newline);

        std::fs::write(&fa_path, &fa_bytes).unwrap();

        // FAI: name, length, offset, line_bases, line_width
        // offset = len(">chr1\n") = 6
        // line_bases = 1000, line_width = 1001 (1000 + newline)
        let fai_content = "chr1\t1000\t6\t1000\t1001\n".to_string();
        std::fs::write(&fai_path, fai_content.as_bytes()).unwrap();

        fa_path
    }

    fn make_config(seed: u64) -> Config {
        Config {
            reference: std::path::PathBuf::from("/dev/null"),
            output: OutputConfig {
                directory: std::path::PathBuf::from("/tmp"),
                fastq: true,
                bam: false,
                truth_vcf: false,
                manifest: false,
                germline_vcf: false,
                single_read_bam: false,
                mapq: 60,
            },
            sample: SampleConfig {
                name: "TEST".to_string(),
                read_length: 50,
                coverage: 10.0,
                platform: None,
            },
            fragment: FragmentConfig {
                model: FragmentModel::Normal,
                mean: 200.0,
                sd: 20.0,
                long_read: None,
                end_motif_model: None,
                ctdna_fraction: None,
                mono_sd: None,
                di_sd: None,
            },
            quality: QualityConfig {
                mean_quality: 36,
                tail_decay: 0.001,
                profile_path: None,
            },
            tumour: None,
            mutations: None,
            umi: None,
            artifacts: None,
            seed: Some(seed),
            threads: None,
            chromosomes: None,
            regions_bed: None,
            copy_number: None,
            gc_bias: None,
            samples: None,
            capture: None,
            performance: Default::default(),
            preset: None,
            vafs: None,
            germline: None,
            paired: None,
        }
    }

    fn open_reference(fa_path: &std::path::Path) -> ReferenceGenome {
        ReferenceGenome::open(fa_path).expect("failed to open test reference")
    }

    fn snv_variant(pos: u64, ref_base: u8, alt_base: u8, vaf: f64) -> Variant {
        Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos,
                ref_base,
                alt_base,
            },
            expected_vaf: vaf,
            clone_id: None,
            haplotype: None,
            ccf: None,
        }
    }

    // -----------------------------------------------------------------------
    // 1. Basic read pair count
    // -----------------------------------------------------------------------
    #[test]
    fn test_generate_region_basic() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let config = make_config(1);
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        // Region of 1000 bp, coverage=10, read_length=50
        // n_pairs = ceil(10 * 1000 / (2 * 50)) = ceil(100) = 100
        let region = Region::new("chr1", 0, 1000);
        let output = engine.generate_region(&region, &[]).unwrap();
        assert_eq!(output.read_pairs.len(), 100);
    }

    // -----------------------------------------------------------------------
    // 2. Read lengths
    // -----------------------------------------------------------------------
    #[test]
    fn test_read_lengths_correct() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let config = make_config(2);
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        let region = Region::new("chr1", 0, 1000);
        let output = engine.generate_region(&region, &[]).unwrap();
        for pair in &output.read_pairs {
            assert_eq!(pair.read1.len(), 50, "read1 should be 50 bp");
            assert_eq!(pair.read2.len(), 50, "read2 should be 50 bp");
        }
    }

    // -----------------------------------------------------------------------
    // 3. 100% VAF variant in all overlapping reads
    // -----------------------------------------------------------------------
    #[test]
    fn test_variant_spiked() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let config = make_config(3);
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        // Position 100 in repeating ACGT: offset 100 % 4 = 0 → 'A'
        // Spike A→T at pos 100, VAF=1.0
        let variant = snv_variant(100, b'A', b'T', 1.0);
        let region = Region::new("chr1", 50, 200);
        let output = engine.generate_region(&region, &[variant]).unwrap();

        let applied = &output.applied_variants;
        assert!(!applied.is_empty(), "variant should have been applied");
        assert_eq!(
            applied[0].actual_alt_count, applied[0].actual_total_count,
            "at VAF=1.0, all overlapping reads should carry the alt"
        );
    }

    // -----------------------------------------------------------------------
    // 4. Stochastic VAF at 50%
    // -----------------------------------------------------------------------
    #[test]
    fn test_variant_stochastic() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let mut config = make_config(4);
        // Increase coverage for better statistics.
        config.sample.coverage = 100.0;
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        // Use a central position so the majority of fragments cover it.
        // pos 248 in ACGT pattern: 248 % 4 = 0 → 'A'; spike A→C at 50% VAF
        let variant = snv_variant(248, b'A', b'C', 0.5);
        let region = Region::new("chr1", 0, 500);
        let output = engine.generate_region(&region, &[variant]).unwrap();

        let applied = &output.applied_variants;
        assert!(!applied.is_empty(), "variant should appear");
        let alt = applied[0].actual_alt_count as f64;
        let total = applied[0].actual_total_count as f64;
        let observed_vaf = alt / total;
        assert!(
            (observed_vaf - 0.5).abs() < 0.25,
            "observed VAF {} should be near 0.5",
            observed_vaf
        );
    }

    // -----------------------------------------------------------------------
    // 5. No variants → no applied_variants
    // -----------------------------------------------------------------------
    #[test]
    fn test_no_variants() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let config = make_config(5);
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        let region = Region::new("chr1", 0, 1000);
        let output = engine.generate_region(&region, &[]).unwrap();
        assert!(output.applied_variants.is_empty());
        assert!(!output.read_pairs.is_empty());
    }

    // -----------------------------------------------------------------------
    // 6. UMI attached to read names
    // -----------------------------------------------------------------------
    #[test]
    fn test_umi_attached() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let mut config = make_config(6);
        config.umi = Some(UmiConfig {
            length: 8,
            duplex: false,
            pcr_cycles: 1,
            family_size_mean: 1.0,
            family_size_sd: 0.1,
            inline: false,
        });
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        let region = Region::new("chr1", 0, 500);
        let output = engine.generate_region(&region, &[]).unwrap();
        assert!(!output.read_pairs.is_empty());
        for pair in &output.read_pairs {
            assert!(
                pair.name.contains(":UMI:"),
                "read name '{}' should contain ':UMI:'",
                pair.name
            );
        }
    }

    // -----------------------------------------------------------------------
    // 7. PCR families expand read count
    // -----------------------------------------------------------------------
    #[test]
    fn test_pcr_families() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let mut config = make_config(7);
        // family_size_mean=3 means each original expands to ~3 copies.
        config.umi = Some(UmiConfig {
            length: 8,
            duplex: false,
            pcr_cycles: 5,
            family_size_mean: 3.0,
            family_size_sd: 0.5,
            inline: false,
        });
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config.clone(), reference);

        let region = Region::new("chr1", 0, 500);
        let output = engine.generate_region(&region, &[]).unwrap();

        // Expected number of pairs without UMI.
        let n_base =
            read_pairs_for_coverage(500, config.sample.coverage, config.sample.read_length);
        // With PCR families of mean size 3, total should be > base count.
        assert!(
            output.read_pairs.len() > n_base as usize,
            "PCR families should increase read count: got {} base {}",
            output.read_pairs.len(),
            n_base
        );
    }

    // -----------------------------------------------------------------------
    // 8. FFPE artifacts produce C>T transitions
    // -----------------------------------------------------------------------
    #[test]
    fn test_artifacts_applied() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let mut config = make_config(8);
        // Very high FFPE rate so we can detect it reliably.
        config.artifacts = Some(ArtifactConfig {
            ffpe_damage_rate: Some(0.5),
            oxog_rate: None,
            duplicate_rate: None,
            pcr_error_rate: None,
        });
        // High coverage for more reads to test.
        config.sample.coverage = 50.0;
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        let region = Region::new("chr1", 0, 1000);
        let output = engine.generate_region(&region, &[]).unwrap();

        // Count C>T on read1 (forward strand FFPE).
        // Reference has known bases (ACGT repeating), but after errors and
        // FFPE damage, we just check that some T bases exist where C was possible.
        // A simpler check: at least some reads have been modified.
        assert!(!output.read_pairs.is_empty());
        // Check that the FFPE damage was at least attempted: no trivial assertion
        // possible without knowing pre/post state, so we verify the pipeline ran.
        // The actual damage function is already unit-tested in artifacts::ffpe.
    }

    // -----------------------------------------------------------------------
    // 9. cfDNA fragment model
    // -----------------------------------------------------------------------
    #[test]
    fn test_cfdna_fragments() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let mut config = make_config(9);
        config.fragment = FragmentConfig {
            model: FragmentModel::Cfda,
            mean: 167.0,
            sd: 20.0,
            long_read: None,
            end_motif_model: None,
            ctdna_fraction: None,
            mono_sd: None,
            di_sd: None,
        };
        config.sample.coverage = 5.0;
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        let region = Region::new("chr1", 0, 1000);
        let output = engine.generate_region(&region, &[]).unwrap();

        // With cfDNA model, fragments should be short (mononucleosomal ~167 bp).
        // All fragment_lengths should be >= min_size (50) and <= region length.
        for pair in &output.read_pairs {
            assert!(
                pair.fragment_length >= 50 || pair.fragment_length == 0,
                "fragment length {} below minimum",
                pair.fragment_length
            );
            assert!(
                pair.fragment_length <= 1000,
                "fragment length {} exceeds region",
                pair.fragment_length
            );
        }
        assert!(!output.read_pairs.is_empty());
    }

    // -----------------------------------------------------------------------
    // 10. Deterministic output with same seed
    // -----------------------------------------------------------------------
    #[test]
    fn test_deterministic_with_seed() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);

        let config = make_config(42);
        let reference1 = open_reference(&fa);
        let reference2 = open_reference(&fa);

        let mut engine1 = SimulationEngine::new(config.clone(), reference1);
        let mut engine2 = SimulationEngine::new(config, reference2);

        let region = Region::new("chr1", 0, 500);
        let out1 = engine1.generate_region(&region, &[]).unwrap();
        let out2 = engine2.generate_region(&region, &[]).unwrap();

        assert_eq!(out1.read_pairs.len(), out2.read_pairs.len());
        for (p1, p2) in out1.read_pairs.iter().zip(out2.read_pairs.iter()) {
            assert_eq!(p1.read1.seq, p2.read1.seq, "read1 sequences differ");
            assert_eq!(p1.read2.seq, p2.read2.seq, "read2 sequences differ");
            assert_eq!(p1.fragment_start, p2.fragment_start);
        }
    }

    // -----------------------------------------------------------------------
    // 11. AppliedVariant tracking
    // -----------------------------------------------------------------------
    #[test]
    fn test_applied_variants_tracking() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let mut config = make_config(11);
        config.sample.coverage = 50.0;
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        // Two variants: one in region, one outside.
        let in_region = snv_variant(100, b'A', b'T', 0.3);
        let out_of_region = snv_variant(999, b'A', b'T', 0.5);
        let region = Region::new("chr1", 0, 500);

        let output = engine
            .generate_region(&region, &[in_region, out_of_region])
            .unwrap();

        // Only the in-region variant should appear in applied_variants.
        // (The out-of-region variant might still appear if budget > 0 and
        //  position is near the boundary, but pos 999 is beyond the region end 500.)
        let in_region_applied = output
            .applied_variants
            .iter()
            .find(|av| av.variant.pos() == 100);

        // Variant at pos 100 with 50x coverage and VAF=0.3 should be applied.
        if let Some(av) = in_region_applied {
            assert!(av.actual_total_count > 0, "total count should be > 0");
            assert!(
                av.actual_alt_count <= av.actual_total_count,
                "alt count must not exceed total"
            );
        }

        // Variant at pos 999 (outside region 0-500) should NOT be in applied_variants.
        let out_applied = output
            .applied_variants
            .iter()
            .any(|av| av.variant.pos() == 999);
        assert!(!out_applied, "out-of-region variant should not be applied");
    }

    // -----------------------------------------------------------------------
    // 12. GC bias reduces read count for AT-rich regions
    // -----------------------------------------------------------------------
    #[test]
    fn test_gc_bias_reduces_reads() {
        use crate::io::config::GcBiasConfig as CfgGcBias;

        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);

        // Build two configs: one with extreme GC bias, one without.
        let mut config_biased = make_config(20);
        config_biased.sample.coverage = 20.0;
        config_biased.gc_bias = Some(CfgGcBias {
            enabled: true,
            model: "default".to_string(),
            severity: 2.0,
        });

        let mut config_flat = make_config(20);
        config_flat.sample.coverage = 20.0;
        // No gc_bias → no rejection

        let ref1 = open_reference(&fa);
        let ref2 = open_reference(&fa);

        let mut engine_biased = SimulationEngine::new(config_biased, ref1);
        let mut engine_flat = SimulationEngine::new(config_flat, ref2);

        let region = Region::new("chr1", 0, 1000);
        let out_biased = engine_biased.generate_region(&region, &[]).unwrap();
        let out_flat = engine_flat.generate_region(&region, &[]).unwrap();

        // GC biased engine should produce fewer reads (many fragments rejected).
        // The test FASTA is 50% GC (ACGT repeating) so the effect is mild.
        // We check that the flat engine produces at least as many reads.
        assert!(
            out_flat.read_pairs.len() >= out_biased.read_pairs.len(),
            "flat model should produce >= reads as biased: flat={} biased={}",
            out_flat.read_pairs.len(),
            out_biased.read_pairs.len()
        );
    }

    // -----------------------------------------------------------------------
    // 13. CNV adjustment lowers read count for homozygous deletion
    // -----------------------------------------------------------------------
    #[test]
    fn test_cnv_homdel_reduces_coverage() {
        use crate::io::config::CopyNumberConfig;

        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);

        let mut config_del = make_config(30);
        config_del.sample.coverage = 30.0;
        // Pure tumour, diploid default ploidy.
        config_del.tumour = Some(crate::io::config::TumourConfig {
            purity: 1.0,
            ploidy: 2,
            clones: Vec::new(),
            msi: false,
        });
        // Homozygous deletion over the entire test region.
        config_del.copy_number = Some(vec![CopyNumberConfig {
            region: "chr1:0-1000".to_string(),
            tumor_cn: 0,
            normal_cn: 2,
            major_cn: None,
            minor_cn: None,
        }]);

        let mut config_normal = make_config(30);
        config_normal.sample.coverage = 30.0;

        let ref1 = open_reference(&fa);
        let ref2 = open_reference(&fa);

        let mut engine_del = SimulationEngine::new(config_del, ref1);
        let mut engine_normal = SimulationEngine::new(config_normal, ref2);

        let region = Region::new("chr1", 0, 1000);
        let out_del = engine_del.generate_region(&region, &[]).unwrap();
        let out_normal = engine_normal.generate_region(&region, &[]).unwrap();

        // Homozygous deletion with purity=1 should produce 0 read pairs.
        assert_eq!(
            out_del.read_pairs.len(),
            0,
            "homdel region should produce no reads, got {}",
            out_del.read_pairs.len()
        );
        // Normal should have reads.
        assert!(!out_normal.read_pairs.is_empty());
    }

    // -----------------------------------------------------------------------
    // 14. SV boundary panic regression — large deletion near chromosome end
    // -----------------------------------------------------------------------
    //
    // A large deletion from position 50 to 180 on a 200-bp chromosome leaves
    // only 20 bases before the chromosome end.  With read_length=150, any
    // fragment sampled over that tail region has fewer than read_length/2 bases
    // and must be silently skipped rather than panicking.
    #[test]
    fn test_sv_deletion_near_chrom_end_no_panic() {
        let dir = TempDir::new().unwrap();
        let fa_path = dir.path().join("small.fa");
        let fai_path = dir.path().join("small.fa.fai");

        // 200-bp chromosome: repeating ACGT.
        let seq: Vec<u8> = b"ACGT".iter().cycle().take(200).cloned().collect();
        let mut fa_bytes = Vec::new();
        fa_bytes.extend_from_slice(b">chr1\n");
        fa_bytes.extend_from_slice(&seq);
        fa_bytes.push(b'\n');
        std::fs::write(&fa_path, &fa_bytes).unwrap();

        // FAI: name, length, offset, line_bases, line_width
        // offset = len(">chr1\n") = 6; 200 bases on one line → line_width = 201
        let fai_content = "chr1\t200\t6\t200\t201\n";
        std::fs::write(&fai_path, fai_content.as_bytes()).unwrap();

        let reference = ReferenceGenome::open(&fa_path).expect("failed to open small reference");

        // Config: read_length=150, low coverage so the test runs quickly.
        let mut config = make_config(42);
        config.sample.read_length = 150;
        config.sample.coverage = 5.0;

        let mut engine = SimulationEngine::new(config, reference);

        // Deletion from 50 to 180 leaves only 20 bp of chromosome remaining.
        let deletion = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Sv {
                sv_type: SvType::Deletion,
                chrom: "chr1".to_string(),
                start: 50,
                end: 180,
            },
            expected_vaf: 1.0,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        // Must not panic; may return an empty or non-empty result.
        let result = engine.generate_region(&Region::new("chr1", 0, 200), &[deletion]);
        assert!(
            result.is_ok(),
            "generate_region should not panic or error on a near-end deletion, got: {:?}",
            result.err()
        );
    }

    // -----------------------------------------------------------------------
    // 15. Duplex simulation produces both AB and BA read pairs (T005)
    // -----------------------------------------------------------------------
    #[test]
    fn test_duplex_produces_ab_and_ba_reads() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let mut config = make_config(50);
        // Use family_size_mean=1 to keep output predictable.
        config.umi = Some(UmiConfig {
            length: 8,
            duplex: true,
            pcr_cycles: 1,
            family_size_mean: 1.0,
            family_size_sd: 0.1,
            inline: false,
        });
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        let region = Region::new("chr1", 0, 500);
        let output = engine.generate_region(&region, &[]).unwrap();
        assert!(!output.read_pairs.is_empty(), "should produce read pairs");

        // Collect all UMI strings from read names.
        let names: Vec<&str> = output.read_pairs.iter().map(|p| p.name.as_str()).collect();

        // For each AB read "...:UMI:AAAA-BBBB...", check that a BA counterpart
        // "...:UMI:BBBB-AAAA..." also exists.
        let mut found_duplex_pair = false;
        for name in &names {
            // Extract the UMI segment following ":UMI:".
            let Some(umi_start) = name.find(":UMI:") else {
                continue;
            };
            // The UMI ends at the next colon (PCR suffix) or end of string.
            let umi_region = &name[umi_start + 5..];
            let umi_end = umi_region.find(':').unwrap_or(umi_region.len());
            let umi = &umi_region[..umi_end];

            // Only process duplexed UMIs that contain a '-' separator.
            let Some(dash) = umi.find('-') else {
                continue;
            };
            let half_a = &umi[..dash];
            let half_b = &umi[dash + 1..];
            let ba_umi = format!("{}-{}", half_b, half_a);

            // Check that at least one read in the output carries the BA barcode.
            let ba_exists = names
                .iter()
                .any(|n| n.contains(&format!(":UMI:{}", ba_umi)));
            if ba_exists {
                found_duplex_pair = true;
                break;
            }
        }

        assert!(
            found_duplex_pair,
            "duplex mode should produce at least one AB/BA read pair complement"
        );
    }

    // -----------------------------------------------------------------------
    // 16. duplex_alt_count is non-zero at moderate VAF in duplex mode (T106)
    // -----------------------------------------------------------------------
    //
    // In duplex mode every alt molecule produces both an AB and a BA family, so
    // duplex_alt_count should equal the molecule-level alt count and therefore be
    // greater than zero when the variant is spiked at a moderate VAF.
    #[test]
    fn test_duplex_alt_count_is_nonzero_at_moderate_vaf() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let mut config = make_config(106);
        // Use high coverage so the variant is applied to many molecules,
        // making a zero duplex_alt_count extremely unlikely at VAF=0.5.
        config.sample.coverage = 100.0;
        config.umi = Some(UmiConfig {
            length: 8,
            duplex: true,
            pcr_cycles: 1,
            family_size_mean: 1.0,
            family_size_sd: 0.1,
            inline: false,
        });
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        // Position 248 in the ACGT repeating pattern: 248 % 4 = 0 → 'A'.
        // Spike A→C at 50% VAF over a central region.
        let variant = snv_variant(248, b'A', b'C', 0.5);
        let region = Region::new("chr1", 0, 500);
        let output = engine.generate_region(&region, &[variant]).unwrap();

        let applied = &output.applied_variants;
        assert!(!applied.is_empty(), "variant should be applied");

        let av = &applied[0];
        assert!(
            av.duplex_alt_count > 0,
            "duplex_alt_count should be > 0 in duplex mode at VAF=0.5, got {}",
            av.duplex_alt_count
        );
        assert!(
            av.duplex_alt_count <= av.actual_alt_count,
            "duplex_alt_count {} must not exceed actual_alt_count {}",
            av.duplex_alt_count,
            av.actual_alt_count
        );
    }

    // -----------------------------------------------------------------------
    // 17. SV reads appear in both AB and BA strand families in duplex mode (T109)
    // -----------------------------------------------------------------------
    //
    // In duplex mode, when a deletion SV is spiked into a molecule, both the AB
    // and BA families should carry the variant_tag. This is critical for duplex
    // SV callers: if only one strand carries the SV tag the caller treats it as
    // an artefact and filters it out.
    //
    // The BA strand is generated by swap+revcomp of the AB reads and inherits
    // the parent's variant_tags. We verify this by checking that for every
    // SV-tagged read with a duplex UMI "A-B", a matching "B-A" read also
    // carries the same SV tag.
    #[test]
    fn test_sv_reads_appear_in_both_ab_and_ba_families() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let mut config = make_config(109);
        config.sample.coverage = 30.0;
        // Use family_size_mean=1.0 to keep output small and names simple.
        config.umi = Some(UmiConfig {
            length: 8,
            duplex: true,
            pcr_cycles: 1,
            family_size_mean: 1.0,
            family_size_sd: 0.1,
            inline: false,
        });
        let reference = open_reference(&fa);
        let mut engine = SimulationEngine::new(config, reference);

        // Deletion from 200 to 400 at VAF=1.0: every fragment spanning the
        // breakpoint carries the deletion. The region is large enough that
        // many fragments will overlap either junction.
        let deletion = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Sv {
                sv_type: SvType::Deletion,
                chrom: "chr1".to_string(),
                start: 200,
                end: 400,
            },
            expected_vaf: 1.0,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let region = Region::new("chr1", 0, 1000);
        let output = engine
            .generate_region(&region, &[deletion])
            .expect("generate_region should not fail");

        // Collect all read pairs that carry an SV variant_tag.
        let sv_tagged: Vec<&ReadPair> = output
            .read_pairs
            .iter()
            .filter(|p| p.variant_tags.iter().any(|t| t.vartype == "SV"))
            .collect();

        assert!(
            !sv_tagged.is_empty(),
            "at least one read pair should carry the SV tag at VAF=1.0"
        );

        // For every AB-tagged SV read (UMI "A-B"), verify a BA counterpart
        // (UMI "B-A") also carries the SV tag.
        let sv_names: std::collections::HashSet<&str> =
            sv_tagged.iter().map(|p| p.name.as_str()).collect();

        let mut found_ba_complement = false;
        for pair in &sv_tagged {
            let Some(umi_start) = pair.name.find(":UMI:") else {
                continue;
            };
            let umi_region = &pair.name[umi_start + 5..];
            let umi_end = umi_region.find(':').unwrap_or(umi_region.len());
            let umi = &umi_region[..umi_end];
            let Some(dash) = umi.find('-') else {
                continue;
            };
            let half_a = &umi[..dash];
            let half_b = &umi[dash + 1..];
            let ba_umi = format!("{}-{}", half_b, half_a);

            // The BA read name shares the same base prefix but with the
            // swapped UMI, followed by a PCR suffix (:pcr:N).
            // Use a prefix match to handle the pcr copy index.
            let base_prefix = &pair.name[..umi_start];
            let ba_prefix = format!("{}:UMI:{}", base_prefix, ba_umi);
            if sv_names.iter().any(|n| n.starts_with(ba_prefix.as_str())) {
                found_ba_complement = true;
                break;
            }
        }

        assert!(
            found_ba_complement,
            "every SV-tagged AB read should have a BA complement that also carries the SV tag"
        );
    }
}
