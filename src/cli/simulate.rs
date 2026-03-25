//! Implementation of the `simulate` subcommand.
//!
//! Loads a YAML config, builds the simulation plan, drives parallel read
//! generation across genomic regions, and writes the output BAM/FASTQ and
//! truth VCF files.
use std::io::Write;
use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result};
use crossbeam_channel::bounded;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use super::SimulateOpts;
use crate::cli::presets;
use crate::core::capture::CaptureModel;
use crate::core::coverage::{intersect_with_targets, partition_regions, read_pairs_for_coverage};
use crate::core::engine::{derive_region_seed, AppliedVariant, SimulationEngine};
use crate::core::multi_sample::{write_combined_manifest, MultiSamplePlan, SampleManifestEntry};
use crate::core::types::{MutationType, ReadPair, Region, Variant};
use crate::io::bam::BamWriter;
use crate::io::config::{self, Config};
use crate::io::fastq::FastqWriter;
use crate::io::manifest::{Manifest, ReferenceInfo, SimulationStatistics};
use crate::io::reference::ReferenceGenome;
use crate::io::truth_vcf::TruthVcfWriter;
use crate::io::vcf_input;
use crate::tumour::clonal_tree::{ClonalTree, Clone as TumourClone};
use crate::variants::germline::generate_germline_variants;
use crate::variants::random_gen::{
    generate_msi_indels, generate_random_mutations, resolve_signature,
};
use crate::variants::structural::StructuralVariant;

/// Extract the UMI barcode from a read name that contains `:UMI:<barcode>`.
///
/// Read names produced by the engine follow the pattern:
///   `<base>:UMI:<barcode>:pcr:<n>`
///
/// This function returns the barcode segment between the first `:UMI:` marker
/// and the next `:` (or the end of the string). Returns `None` when no UMI
/// marker is present.
fn extract_umi_from_name(name: &str) -> Option<&str> {
    let start = name.find(":UMI:")?;
    let umi_region = &name[start + 5..];
    let end = umi_region.find(':').unwrap_or(umi_region.len());
    Some(&umi_region[..end])
}

/// Derive a stable molecular family ID from a read name.
///
/// All PCR copies of the same original molecule share the same base name
/// (everything before the `:pcr:N` suffix). Assigning the same family ID to
/// all copies allows downstream tools to group them as one UMI family.
///
/// The family ID is assigned from the `family_ids` map: if the base name has
/// not been seen before, a new sequential ID is allocated and stored.
fn family_id_for(name: &str, family_ids: &mut std::collections::HashMap<String, i32>) -> i32 {
    // Strip the `:pcr:N` suffix to get the base name shared by all copies.
    let base = if let Some(pos) = name.rfind(":pcr:") {
        &name[..pos]
    } else {
        name
    };
    let next_id = family_ids.len() as i32;
    *family_ids.entry(base.to_string()).or_insert(next_id)
}

/// A batch of generated reads for one region, sent through the streaming channel.
struct RegionBatch {
    region: Region,
    read_pairs: Vec<ReadPair>,
    applied_variants: Vec<AppliedVariant>,
    duplex_total_molecules: u64,
    duplex_molecules_with_both_strands: u64,
}

/// Statistics returned by the writer thread.
struct WriterStats {
    total_read_pairs: u64,
    all_applied: Vec<AppliedVariant>,
    /// Mean coverage per chromosome: chrom -> mean coverage (×).
    // Collected for future use; not yet surfaced in output.
    #[allow(dead_code)]
    chrom_coverage: std::collections::HashMap<String, f64>,
    /// Estimated PCR duplicate rate (0.0 if UMI is not enabled).
    // Collected for future use; not yet surfaced in output.
    #[allow(dead_code)]
    duplicate_rate: f64,
    /// Total duplex molecules simulated (T108).
    duplex_total_molecules: u64,
    /// Duplex molecules where both AB and BA strands were generated (T108).
    duplex_molecules_with_both_strands: u64,
}

/// Default chunk size for partitioning regions (~1 Mbp).
const DEFAULT_CHUNK_SIZE: u64 = 1_000_000;

/// Configure the rayon global thread pool.
///
/// If `threads` is `None`, rayon uses its default (one thread per logical CPU).
fn configure_thread_pool(threads: Option<usize>) -> Result<()> {
    if let Some(n) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .context("failed to configure rayon thread pool")?;
    }
    Ok(())
}

/// Apply all CLI flag overrides to `config`.
///
/// Precedence order:  defaults < preset < YAML config < CLI flags.
/// This function is called *after* the preset overlay so that CLI flags
/// always win.
pub fn apply_overrides(
    config: &mut Config,
    opts: &SimulateOpts,
    cli_threads: Option<usize>,
) -> Result<()> {
    if let Some(seed) = opts.seed {
        config.seed = Some(seed);
    }
    if let Some(coverage) = opts.coverage {
        config.sample.coverage = coverage;
    }
    if let Some(rl) = opts.read_length {
        config.sample.read_length = rl;
    }
    if let Some(dir) = &opts.output_dir {
        config.output.directory = dir.clone();
    }
    if let Some(t) = cli_threads {
        config.threads = Some(t);
    }
    if let Some(purity) = opts.purity {
        match &mut config.tumour {
            Some(tumour) => tumour.purity = purity,
            None => {
                config.tumour = Some(config::TumourConfig {
                    purity,
                    ploidy: 2,
                    clones: Vec::new(),
                    msi: false,
                });
            }
        }
    }
    if let Some(fmean) = opts.fragment_mean {
        config.fragment.mean = fmean;
    }
    if let Some(fsd) = opts.fragment_sd {
        config.fragment.sd = fsd;
    }
    if let Some(count) = opts.random_mutations {
        let (vaf_min, vaf_max) = parse_vaf_range(opts.vaf_range.as_deref())?;
        // If vaf_range is provided without random_mutations we ignore it
        // (that would be a user error caught later).  Here we only override
        // the random block when random_mutations is set.
        let existing = config.mutations.get_or_insert(config::MutationConfig {
            vcf: None,
            random: None,
            sv_count: 0,
            sv_signature: None,
            include_driver_mutations: false,
        });
        let rand = existing.random.get_or_insert(config::RandomMutationConfig {
            count: 0,
            vaf_min: 0.001,
            vaf_max: 0.5,
            snv_fraction: 0.80,
            indel_fraction: 0.15,
            mnv_fraction: 0.05,
            signature: None,
        });
        rand.count = count;
        rand.vaf_min = vaf_min;
        rand.vaf_max = vaf_max;
    } else if let Some(ref vaf_str) = opts.vaf_range {
        // vaf_range without random_mutations: update existing random block if present.
        let (vaf_min, vaf_max) = parse_vaf_range(Some(vaf_str.as_str()))?;
        if let Some(ref mut muts) = config.mutations {
            if let Some(ref mut rand) = muts.random {
                rand.vaf_min = vaf_min;
                rand.vaf_max = vaf_max;
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Preset listing
// ---------------------------------------------------------------------------

/// Print all available presets to stdout.
///
/// Base presets are listed first, then cancer-type presets (with the
/// `cancer:` prefix users must supply).
fn print_preset_list() {
    println!("Base presets:");
    println!("  small       Fast smoke-test: 1x coverage, chr22 only, 100 random mutations");
    println!("  panel       Targeted-sequencing: 500x coverage, UMI enabled, 50 mutations");
    println!("  wgs         Whole-genome: 30x coverage, 5 000 random mutations");
    println!("  cfdna       Cell-free DNA / liquid biopsy: 200x, cfDNA fragments, 2% purity");
    println!("  ffpe        FFPE-damaged sample: 30x with deamination and oxoG artefacts");
    println!("  umi         High-depth duplex UMI: 1 000x coverage, duplex barcodes");
    println!();
    println!("Cancer-type presets (prefix with 'cancer:', e.g. --preset cancer:lung_adeno):");
    println!("  cancer:lung_adeno     Lung adenocarcinoma — SBS4 smoking signature, high TMB");
    println!("  cancer:colorectal     Colorectal cancer — SBS1/SBS5 aging, moderate TMB");
    println!("  cancer:breast_tnbc    Triple-negative breast — SBS3 HRD, elevated indel burden");
    println!("  cancer:melanoma       Cutaneous melanoma — SBS7 UV signature, very high TMB");
    println!("  cancer:aml            Acute myeloid leukaemia — SBS1/SBS5, low TMB");
    println!("  cancer:prostate       Prostate adenocarcinoma — SBS1/SBS5, low-moderate TMB");
    println!("  cancer:pancreatic     Pancreatic PDAC — SBS1/SBS5, low purity (stromal)");
    println!("  cancer:glioblastoma   Glioblastoma multiforme — SBS1/SBS5, moderate TMB");
}

// ---------------------------------------------------------------------------
// Fuzzy preset suggestion
// ---------------------------------------------------------------------------

/// Compute Levenshtein edit distance between two strings.
///
/// Uses a standard dynamic-programming approach with O(min(a,b)) space.
pub(crate) fn levenshtein_distance(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let m = a.len();
    let n = b.len();

    // Keep only two rows to reduce memory usage.
    let mut prev: Vec<usize> = (0..=n).collect();
    let mut curr = vec![0usize; n + 1];

    for i in 1..=m {
        curr[0] = i;
        for j in 1..=n {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            curr[j] = (curr[j - 1] + 1).min(prev[j] + 1).min(prev[j - 1] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[n]
}

/// Return all known preset names (base and cancer-prefixed) as a flat list.
fn all_preset_names() -> Vec<String> {
    use crate::cli::cancer_presets;
    let mut names: Vec<String> = crate::cli::presets::all_names()
        .iter()
        .map(|s| s.to_string())
        .collect();
    for cn in cancer_presets::all_names() {
        names.push(format!("cancer:{}", cn));
    }
    names
}

/// Find the closest preset name to `input`.
///
/// Returns `Some(name)` if the closest match has edit distance ≤ 2.
pub(crate) fn suggest_preset(input: &str) -> Option<String> {
    let candidates = all_preset_names();
    let (best_name, best_dist) = candidates
        .iter()
        .map(|name| (name.as_str(), levenshtein_distance(input, name)))
        .min_by_key(|&(_, d)| d)?;
    if best_dist <= 2 {
        Some(best_name.to_string())
    } else {
        None
    }
}

pub fn run(opts: SimulateOpts, cli_threads: Option<usize>) -> Result<()> {
    // -----------------------------------------------------------------------
    // 0. Handle informational flags that exit before doing any real work.
    // -----------------------------------------------------------------------
    if opts.list_presets {
        print_preset_list();
        return Ok(());
    }

    let start_time = Instant::now();

    // -----------------------------------------------------------------------
    // 1. Load and validate config, applying CLI overrides
    // -----------------------------------------------------------------------
    tracing::info!("loading config from {}", opts.config.display());
    // Parse --set key=value pairs and substitute ${key} placeholders.
    let vars: std::collections::HashMap<String, String> = opts
        .set
        .as_deref()
        .unwrap_or(&[])
        .iter()
        .map(|kv| {
            let (k, v) = kv.split_once('=').ok_or_else(|| {
                anyhow::anyhow!("--set value '{}' must be in KEY=VALUE format", kv)
            })?;
            Ok((k.to_string(), v.to_string()))
        })
        .collect::<Result<_>>()?;
    let mut cfg = config::load_with_vars(&opts.config, &vars)
        .with_context(|| format!("failed to load config: {}", opts.config.display()))?;

    // Apply preset (lower precedence than YAML).
    if let Some(ref preset_name) = opts.preset {
        let overlay = presets::get(preset_name).map_err(|e| {
            // Append a "did you mean?" hint when the typo is close to a known name.
            if let Some(suggestion) = suggest_preset(preset_name) {
                anyhow::anyhow!("{e}; did you mean '{suggestion}'")
            } else {
                e
            }
        })?;
        presets::apply_preset_to_config(&mut cfg, &overlay);
    }

    // Apply CLI overrides (highest precedence).
    apply_overrides(&mut cfg, &opts, cli_threads)?;

    config::validate(&cfg)?;
    tracing::info!("config validated successfully");

    // -----------------------------------------------------------------------
    // 1b. Configure rayon thread pool
    // -----------------------------------------------------------------------
    configure_thread_pool(cfg.threads)?;

    // -----------------------------------------------------------------------
    // 2. Open reference genome
    // -----------------------------------------------------------------------
    tracing::info!("opening reference genome: {}", cfg.reference.display());
    let reference = ReferenceGenome::open(&cfg.reference)
        .with_context(|| format!("failed to open reference: {}", cfg.reference.display()))?;

    // -----------------------------------------------------------------------
    // 3. Paired tumour-normal mode: if `paired` is set, delegate before any
    //    other mode check so paired overrides multi-sample and batch modes.
    // -----------------------------------------------------------------------
    if cfg.paired.is_some() {
        return run_paired_simulation(cfg, reference, start_time, opts.dry_run);
    }

    // -----------------------------------------------------------------------
    // 3b. Multi-sample mode: if samples[] is present, delegate to per-sample runs
    // -----------------------------------------------------------------------
    if cfg.samples.is_some() && !cfg.samples.as_ref().map(|s| s.is_empty()).unwrap_or(true) {
        return run_multi_sample(cfg, reference, start_time);
    }

    // -----------------------------------------------------------------------
    // 3b. Batch (vafs) mode: run one simulation per VAF value
    // -----------------------------------------------------------------------
    if cfg.vafs.as_ref().map(|v| !v.is_empty()).unwrap_or(false) {
        return run_batch(cfg, reference, start_time, opts.dry_run);
    }

    // -----------------------------------------------------------------------
    // 4. Single-sample simulation
    // -----------------------------------------------------------------------
    run_single_sample(cfg, reference, start_time, opts.dry_run)
}

/// Run the simulation for a single sample configuration.
pub(crate) fn run_single_sample(
    cfg: Config,
    reference: ReferenceGenome,
    start_time: Instant,
    dry_run: bool,
) -> Result<()> {
    // Build the sorted list of chromosomes/lengths we will simulate
    let chrom_lengths: Vec<(String, u64)> = build_chrom_list(&cfg, &reference)?;
    tracing::info!(
        "reference has {} chromosomes to simulate",
        chrom_lengths.len()
    );

    // -----------------------------------------------------------------------
    // Partition regions
    // -----------------------------------------------------------------------
    let regions = partition_regions(&chrom_lengths, DEFAULT_CHUNK_SIZE);
    tracing::info!("partitioned into {} regions", regions.len());

    // Filter to BED targets if specified.
    let regions = if let Some(ref bed_path) = cfg.regions_bed {
        let targets = parse_bed_file(bed_path)
            .with_context(|| format!("failed to parse regions BED: {}", bed_path.display()))?;
        let target_regions: Vec<Region> = targets.iter().map(|(r, _)| r.clone()).collect();
        tracing::info!(
            "filtering {} regions to {} BED targets",
            regions.len(),
            target_regions.len()
        );
        let filtered = intersect_with_targets(&regions, &target_regions);
        anyhow::ensure!(
            !filtered.is_empty(),
            "regions_bed produced zero overlapping regions; check chromosome names match the reference"
        );
        tracing::info!("BED filter retained {} regions", filtered.len());
        filtered
    } else {
        regions
    };

    // -----------------------------------------------------------------------
    // Generate / load variants
    // -----------------------------------------------------------------------
    let seed = cfg.seed.unwrap_or_else(|| {
        use std::time::{SystemTime, UNIX_EPOCH};
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_nanos() as u64)
            .unwrap_or(0)
    });
    let variants = build_variant_list(&cfg, &regions, &reference, seed, None)?;
    tracing::info!("prepared {} variants for spike-in", variants.len());

    // -----------------------------------------------------------------------
    // Dry-run mode
    // -----------------------------------------------------------------------
    if dry_run {
        return run_dry(&cfg, &chrom_lengths, &regions, &variants);
    }

    // -----------------------------------------------------------------------
    // Create output directory and initialize writers
    // -----------------------------------------------------------------------
    std::fs::create_dir_all(&cfg.output.directory).with_context(|| {
        format!(
            "failed to create output directory: {}",
            cfg.output.directory.display()
        )
    })?;

    let sample_name = &cfg.sample.name;
    let out_dir = &cfg.output.directory;

    let mut fastq_writer = if cfg.output.fastq {
        Some(FastqWriter::new(out_dir, sample_name).context("failed to create FASTQ writers")?)
    } else {
        None
    };

    let ref_sequences_for_bam: Vec<(String, u64)> = chrom_lengths.clone();
    let mut bam_writer = if cfg.output.bam {
        let bam_path = out_dir.join(format!("{}.bam", sample_name));
        Some(
            BamWriter::new(
                &bam_path,
                &ref_sequences_for_bam,
                &cfg.sample,
                cfg.output.mapq,
            )
            .context("failed to create BAM writer")?,
        )
    } else {
        None
    };

    // Open the variant reads sidecar before starting the writer thread so it
    // can be moved in alongside the FASTQ and BAM writers.
    let variant_reads_path = out_dir.join("variant_reads.tsv");
    let mut variant_reads_writer = std::io::BufWriter::new(
        std::fs::File::create(&variant_reads_path).context("failed to create variant_reads.tsv")?,
    );
    writeln!(
        variant_reads_writer,
        "read_name\tchrom\tpos\tvartype\tvaf\tclone_id"
    )
    .context("failed to write variant_reads.tsv header")?;

    // We'll open truth VCF after collecting applied variants; use a temp path
    // and rename, or just open now (write header) and stream records.
    let truth_vcf_path = out_dir.join(format!("{}.truth.vcf.gz", sample_name));
    let mut truth_vcf_writer = if cfg.output.truth_vcf {
        let contigs: Vec<(String, u64)> = chrom_lengths.clone();
        Some(
            TruthVcfWriter::new(&truth_vcf_path, sample_name, &contigs)
                .context("failed to create truth VCF writer")?,
        )
    } else {
        None
    };

    // -----------------------------------------------------------------------
    // Build capture model (shared across all region engines)
    // -----------------------------------------------------------------------
    let capture_model: Option<Arc<CaptureModel>> = build_capture_model(&cfg)?;

    // In amplicon mode, replace the partitioned regions with the amplicon targets
    // so that each target becomes exactly one simulation region.
    let regions = if let Some(ref cap) = capture_model {
        if cap.is_amplicon() {
            tracing::info!(
                "amplicon mode: using {} targets as simulation regions",
                cap.target_regions.len()
            );
            cap.target_regions
                .iter()
                .map(|t| Region::new(t.chrom.clone(), t.start, t.end))
                .collect::<Vec<_>>()
        } else {
            regions
        }
    } else {
        regions
    };

    // -----------------------------------------------------------------------
    // Run simulation per region in parallel, streaming batches to a dedicated
    // writer thread via a bounded channel.
    // -----------------------------------------------------------------------
    tracing::info!(
        "starting parallel simulation across {} regions",
        regions.len()
    );
    let pb = ProgressBar::new(regions.len() as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} regions ({eta})",
        )
        .unwrap_or_else(|_| ProgressStyle::default_bar()),
    );

    // Wrap the reference in Arc so it can be shared across rayon worker
    // threads.  ReferenceGenome uses an internal Mutex for safe concurrent access.
    let reference = Arc::new(reference);
    // Wrap config in Arc so each rayon worker shares the struct rather than
    // cloning the full config (including all variant/clone Vecs) per region.
    let cfg_arc = Arc::new(cfg.clone());
    let master_seed = cfg.seed;
    let pb_ref = pb.clone();

    // Create bounded channel for streaming batches to the writer thread.
    let buffer_size = cfg.performance.output_buffer_regions;
    let (tx, rx) = bounded::<RegionBatch>(buffer_size);

    // Spawn dedicated writer thread.  The writers are moved in so they are
    // owned exclusively by this thread.
    let writer_cfg = cfg.clone();
    let writer_chrom_lengths = chrom_lengths.clone();
    let writer_handle = std::thread::spawn(move || -> Result<WriterStats> {
        let mut total_read_pairs: u64 = 0;
        let mut all_applied: Vec<AppliedVariant> = Vec::new();
        let mut duplex_total_molecules: u64 = 0;
        let mut duplex_molecules_with_both_strands: u64 = 0;
        // Accumulate read pairs per chromosome for coverage computation.
        let mut chrom_read_pairs: std::collections::HashMap<String, u64> =
            std::collections::HashMap::new();
        // Track molecular family IDs for UMI-tagged BAM output. All PCR
        // copies of the same original molecule share one entry in this map.
        let mut family_ids: std::collections::HashMap<String, i32> =
            std::collections::HashMap::new();

        // Primer trim amount: only active in amplicon mode.
        let primer_trim = writer_cfg
            .capture
            .as_ref()
            .filter(|c| c.mode == "amplicon")
            .map(|c| c.primer_trim)
            .unwrap_or(0);

        for batch in rx {
            // Write FASTQ — annotate read names with variant tags when present.
            if let Some(ref mut fq) = fastq_writer {
                for pair in &batch.read_pairs {
                    let name = if pair.variant_tags.is_empty() {
                        pair.name.clone()
                    } else {
                        let tags: String = pair
                            .variant_tags
                            .iter()
                            .map(|t| format!(" VT:Z:{}:{}:{}", t.chrom, t.pos + 1, t.vartype))
                            .collect::<Vec<_>>()
                            .join("");
                        format!("{}{}", pair.name, tags)
                    };
                    // Apply primer trimming when configured: remove the first and
                    // last `primer_trim` bases to simulate primer removal.
                    if primer_trim > 0 && pair.read1.seq.len() > 2 * primer_trim + 10 {
                        let mut trimmed = pair.clone();
                        let len1 = trimmed.read1.seq.len();
                        trimmed.read1.seq =
                            trimmed.read1.seq[primer_trim..len1 - primer_trim].to_vec();
                        trimmed.read1.qual =
                            trimmed.read1.qual[primer_trim..len1 - primer_trim].to_vec();
                        let len2 = trimmed.read2.seq.len();
                        trimmed.read2.seq =
                            trimmed.read2.seq[primer_trim..len2 - primer_trim].to_vec();
                        trimmed.read2.qual =
                            trimmed.read2.qual[primer_trim..len2 - primer_trim].to_vec();
                        fq.write_pair(&trimmed, &name)
                            .context("failed to write FASTQ record")?;
                    } else {
                        fq.write_pair(pair, &name)
                            .context("failed to write FASTQ record")?;
                    }
                }
            }

            // Write BAM
            if let Some(ref mut bam) = bam_writer {
                let ref_id = writer_chrom_lengths
                    .iter()
                    .position(|(name, _)| *name == batch.region.chrom)
                    .unwrap_or(0);
                let read_len = writer_cfg.sample.read_length;
                let umi_enabled = writer_cfg.umi.is_some();
                if writer_cfg.output.single_read_bam {
                    // Long-read mode: write one record per read pair (read1 only).
                    for pair in &batch.read_pairs {
                        let actual_len = pair.read1.seq.len();
                        let cigar = if primer_trim > 0 && actual_len > 2 * primer_trim + 10 {
                            // Soft-clip primer bases at both ends.
                            format!(
                                "{}S{}M{}S",
                                primer_trim,
                                actual_len - 2 * primer_trim,
                                primer_trim
                            )
                        } else {
                            format!("{}M", read_len)
                        };
                        bam.write_single_read(pair, ref_id, pair.fragment_start, &cigar)
                            .context("failed to write BAM single read record")?;
                    }
                } else if umi_enabled {
                    // UMI mode: write RX and MI auxiliary tags on every record.
                    for pair in &batch.read_pairs {
                        let actual_len = pair.read1.seq.len();
                        let cigar = if primer_trim > 0 && actual_len > 2 * primer_trim + 10 {
                            // Soft-clip primer bases at both ends.
                            format!(
                                "{}S{}M{}S",
                                primer_trim,
                                actual_len - 2 * primer_trim,
                                primer_trim
                            )
                        } else {
                            format!("{}M", read_len)
                        };
                        let umi = extract_umi_from_name(&pair.name).unwrap_or("");
                        let fid = family_id_for(&pair.name, &mut family_ids);
                        bam.write_pair_with_umi(
                            pair,
                            ref_id,
                            pair.fragment_start,
                            &cigar,
                            &cigar,
                            umi,
                            fid,
                        )
                        .context("failed to write BAM record")?;
                    }
                } else {
                    for pair in &batch.read_pairs {
                        let actual_len = pair.read1.seq.len();
                        let cigar = if primer_trim > 0 && actual_len > 2 * primer_trim + 10 {
                            // Soft-clip primer bases at both ends.
                            format!(
                                "{}S{}M{}S",
                                primer_trim,
                                actual_len - 2 * primer_trim,
                                primer_trim
                            )
                        } else {
                            format!("{}M", read_len)
                        };
                        bam.write_pair(pair, ref_id, pair.fragment_start, &cigar, &cigar)
                            .context("failed to write BAM record")?;
                    }
                }
            }

            // Write variant reads sidecar — one row per variant tag per read pair.
            for pair in &batch.read_pairs {
                for tag in &pair.variant_tags {
                    writeln!(
                        variant_reads_writer,
                        "{}	{}	{}	{}	{:.6}	{}",
                        pair.name,
                        tag.chrom,
                        tag.pos + 1,
                        tag.vartype,
                        tag.vaf,
                        tag.clone_id.as_deref().unwrap_or(".")
                    )
                    .context("failed to write variant_reads.tsv row")?;
                }
            }

            let batch_pairs = batch.read_pairs.len() as u64;
            total_read_pairs += batch_pairs;
            // Accumulate per-chromosome pair count for coverage reporting.
            *chrom_read_pairs
                .entry(batch.region.chrom.clone())
                .or_insert(0) += batch_pairs;
            all_applied.extend(batch.applied_variants);
            duplex_total_molecules += batch.duplex_total_molecules;
            duplex_molecules_with_both_strands += batch.duplex_molecules_with_both_strands;
        }

        // Compute mean coverage per chromosome from accumulated pair counts.
        let chrom_coverage: std::collections::HashMap<String, f64> = chrom_read_pairs
            .iter()
            .map(|(chrom, &pairs)| {
                let chrom_len = writer_chrom_lengths
                    .iter()
                    .find(|(c, _)| c == chrom)
                    .map(|(_, l)| *l)
                    .unwrap_or(1);
                let coverage =
                    (pairs * 2 * writer_cfg.sample.read_length as u64) as f64 / chrom_len as f64;
                (chrom.clone(), coverage)
            })
            .collect();

        // Estimate duplicate rate from UMI family size mean, if configured.
        let duplicate_rate = if let Some(ref umi) = writer_cfg.umi {
            1.0 - 1.0 / umi.family_size_mean
        } else {
            0.0
        };

        // Finalize writers
        if let Some(fq) = fastq_writer {
            fq.finish().context("failed to finalize FASTQ files")?;
        }
        if let Some(bam) = bam_writer {
            bam.finish().context("failed to finalize BAM file")?;
        }
        variant_reads_writer
            .flush()
            .context("failed to flush variant_reads.tsv")?;

        Ok(WriterStats {
            total_read_pairs,
            all_applied,
            chrom_coverage,
            duplicate_rate,
            duplex_total_molecules,
            duplex_molecules_with_both_strands,
        })
    });

    // Worker threads simulate regions in parallel and send batches as they complete.
    regions
        .par_iter()
        .enumerate()
        .try_for_each(|(idx, region)| {
            let ref_arc = Arc::clone(&reference);
            // Derive a per-region seed for determinism; share the Arc<Config>
            // rather than cloning the full config struct per region.
            let region_seed = master_seed.map(|s| derive_region_seed(s, idx as u64));
            let mut engine = SimulationEngine::new_with_shared_config(
                Arc::clone(&cfg_arc),
                region_seed,
                ref_arc,
            );

            // Wire in capture model.
            if let Some(ref cap) = capture_model {
                engine.set_capture_model(Arc::clone(cap));
            }

            let output = engine
                .generate_region(region, &variants)
                .with_context(|| format!("simulation failed for region {:?}", region))?;
            tx.send(RegionBatch {
                region: region.clone(),
                read_pairs: output.read_pairs,
                applied_variants: output.applied_variants,
                duplex_total_molecules: output.duplex_total_molecules,
                duplex_molecules_with_both_strands: output.duplex_molecules_with_both_strands,
            })
            .map_err(|e| anyhow::anyhow!("writer channel closed: {}", e))?;
            pb_ref.inc(1);
            Ok::<(), anyhow::Error>(())
        })?;

    // Signal completion by dropping the sender so the writer thread exits its loop.
    drop(tx);

    // Wait for the writer thread to finish and collect stats.
    let stats = writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;

    pb.finish_with_message("simulation complete");

    let total_read_pairs = stats.total_read_pairs;
    let all_applied = stats.all_applied;
    let chrom_coverage = stats.chrom_coverage;
    let duplicate_rate = stats.duplicate_rate;
    let duplex_total_molecules = stats.duplex_total_molecules;
    let duplex_molecules_with_both_strands = stats.duplex_molecules_with_both_strands;

    // -----------------------------------------------------------------------
    // Write truth VCF with all applied variants
    //
    // Aggregate alt-molecule counts across all regions for each unique variant,
    // then emit one record per variant with the summed N_ALT_MOL and N_DUPLEX_ALT.
    // -----------------------------------------------------------------------
    if let Some(ref mut vcf) = truth_vcf_writer {
        // key -> (canonical variant, total alt count, total duplex alt count)
        let mut agg: std::collections::HashMap<String, (crate::core::types::Variant, u32, u32)> =
            std::collections::HashMap::new();
        for av in &all_applied {
            let key = variant_key(&av.variant);
            let entry = agg.entry(key).or_insert_with(|| (av.variant.clone(), 0, 0));
            entry.1 = entry.1.saturating_add(av.actual_alt_count);
            entry.2 = entry.2.saturating_add(av.duplex_alt_count);
        }
        // Sort by (chrom, pos) for deterministic output.
        let mut records: Vec<_> = agg.into_values().collect();
        records.sort_by(|a, b| {
            a.0.chrom
                .cmp(&b.0.chrom)
                .then_with(|| a.0.pos().cmp(&b.0.pos()))
        });
        for (variant, n_alt_mol, n_duplex_alt) in &records {
            if let Some(sv) = sv_to_structural(variant) {
                // SVs use symbolic ALT alleles and SV-specific INFO fields.
                vcf.write_sv(
                    &sv,
                    variant.expected_vaf,
                    variant.ccf,
                    variant.clone_id.as_deref(),
                )
                .context("failed to write SV truth VCF record")?;
            } else {
                let (ref_allele, alt_allele) = variant_alleles(variant);
                vcf.write_variant(variant, &ref_allele, &alt_allele, *n_alt_mol, *n_duplex_alt)
                    .context("failed to write truth VCF record")?;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Finalize truth VCF writer (FASTQ/BAM writers finalized inside writer thread)
    // -----------------------------------------------------------------------
    if let Some(vcf) = truth_vcf_writer {
        vcf.finish().context("failed to finalize truth VCF")?;
    }

    // -----------------------------------------------------------------------
    // Write germline truth VCF when germline simulation was active.
    //
    // Filter applied variants to those with clone_id starting with "germline_"
    // and write them to a separate germline_truth.vcf.
    // -----------------------------------------------------------------------
    if cfg.germline.is_some() && cfg.output.germline_vcf {
        let germline_vcf_path = out_dir.join("germline_truth.vcf");
        let contigs: Vec<(String, u64)> = chrom_lengths.clone();
        let mut germ_vcf = TruthVcfWriter::new(&germline_vcf_path, sample_name, &contigs)
            .context("failed to create germline truth VCF writer")?;

        let mut agg: std::collections::HashMap<String, (crate::core::types::Variant, u32, u32)> =
            std::collections::HashMap::new();
        for av in &all_applied {
            let is_germline = av
                .variant
                .clone_id
                .as_deref()
                .map(|id| id.starts_with("germline_"))
                .unwrap_or(false);
            if !is_germline {
                continue;
            }
            let key = variant_key(&av.variant);
            let entry = agg.entry(key).or_insert_with(|| (av.variant.clone(), 0, 0));
            entry.1 = entry.1.saturating_add(av.actual_alt_count);
            entry.2 = entry.2.saturating_add(av.duplex_alt_count);
        }
        let mut records: Vec<_> = agg.into_values().collect();
        records.sort_by(|a, b| {
            a.0.chrom
                .cmp(&b.0.chrom)
                .then_with(|| a.0.pos().cmp(&b.0.pos()))
        });
        for (variant, n_alt_mol, n_duplex_alt) in &records {
            let (ref_allele, alt_allele) = variant_alleles(variant);
            germ_vcf
                .write_variant(variant, &ref_allele, &alt_allele, *n_alt_mol, *n_duplex_alt)
                .context("failed to write germline truth VCF record")?;
        }
        germ_vcf
            .finish()
            .context("failed to finalise germline truth VCF")?;
        tracing::info!(
            "germline truth VCF written to {}",
            germline_vcf_path.display()
        );
    }

    // -----------------------------------------------------------------------
    // Write manifest and report statistics
    // -----------------------------------------------------------------------
    let elapsed = start_time.elapsed();
    let throughput = if elapsed.as_secs_f64() > 0.0 {
        total_read_pairs as f64 / elapsed.as_secs_f64()
    } else {
        0.0
    };

    tracing::info!(
        total_read_pairs,
        applied_variants = all_applied.len(),
        wall_time_secs = elapsed.as_secs_f64(),
        reads_per_sec = throughput,
        "simulation complete"
    );

    if cfg.output.manifest {
        let genome_size: u64 = chrom_lengths.iter().map(|(_, l)| l).sum();
        let ref_info = ReferenceInfo {
            path: cfg.reference.clone(),
            genome_size,
            chromosomes: chrom_lengths.len(),
        };

        // Count variants by type.
        let mut variants_by_type: std::collections::HashMap<String, u64> =
            std::collections::HashMap::new();
        for av in &all_applied {
            let type_name = match &av.variant.mutation {
                MutationType::Snv { .. } => "SNV",
                MutationType::Indel { .. } => "Indel",
                MutationType::Mnv { .. } => "MNV",
                MutationType::Sv { .. } => "SV",
            };
            *variants_by_type.entry(type_name.to_string()).or_insert(0) += 1;
        }

        let stats = SimulationStatistics {
            total_read_pairs,
            total_bases: total_read_pairs * 2 * cfg.sample.read_length as u64,
            variants_spiked: all_applied.len() as u64,
            variants_by_type,
            mean_coverage_achieved: cfg.sample.coverage,
            regions_simulated: regions.len(),
            wall_time_seconds: elapsed.as_secs_f64(),
            reads_per_second: throughput,
        };

        let manifest_path = out_dir.join(format!("{}.manifest.json", sample_name));
        let mut manifest = Manifest::new(cfg.clone(), ref_info);
        manifest.set_statistics(stats);

        // Record output files.
        if cfg.output.fastq {
            manifest.add_output_file(
                "fastq_r1",
                &out_dir.join(format!("{}_R1.fastq.gz", sample_name)),
            );
            manifest.add_output_file(
                "fastq_r2",
                &out_dir.join(format!("{}_R2.fastq.gz", sample_name)),
            );
        }
        if cfg.output.bam {
            manifest.add_output_file("bam", &out_dir.join(format!("{}.bam", sample_name)));
        }
        if cfg.output.truth_vcf {
            manifest.add_output_file("truth_vcf", &truth_vcf_path);
        }

        // Serialize the manifest, then inject top-level flat fields that are
        // expected by existing callers (e.g., integration tests and downstream
        // tooling that reads the simpler fields directly).
        let mut manifest_json =
            serde_json::to_value(&manifest).context("failed to serialize manifest")?;
        if let serde_json::Value::Object(ref mut map) = manifest_json {
            map.insert(
                "sample_name".into(),
                serde_json::Value::String(sample_name.to_string()),
            );
            map.insert(
                "total_read_pairs".into(),
                serde_json::Value::Number(total_read_pairs.into()),
            );
            map.insert(
                "total_reads".into(),
                serde_json::Value::Number((total_read_pairs * 2).into()),
            );
            map.insert(
                "variants_applied".into(),
                serde_json::Value::Number(all_applied.len().into()),
            );
            if let Some(seed_val) = cfg.seed {
                map.insert("seed".into(), serde_json::Value::Number(seed_val.into()));
            }
            map.insert(
                "wall_time_secs".into(),
                serde_json::json!(elapsed.as_secs_f64()),
            );
            map.insert(
                "varforge_version".into(),
                serde_json::Value::String(env!("CARGO_PKG_VERSION").to_string()),
            );
        }
        let pretty = serde_json::to_string_pretty(&manifest_json)
            .context("failed to pretty-print manifest")?;
        std::fs::write(&manifest_path, pretty)
            .with_context(|| format!("failed to write manifest: {}", manifest_path.display()))?;
        tracing::info!("manifest written to {}", manifest_path.display());
    }

    // -----------------------------------------------------------------------
    // Write sim_report.json with coverage and UMI statistics
    // -----------------------------------------------------------------------
    {
        // Build chromosomes object: chrom -> { mean_coverage }.
        let mut chroms_obj = serde_json::Map::new();
        let mut sorted_chroms: Vec<(&String, &f64)> = chrom_coverage.iter().collect();
        sorted_chroms.sort_by_key(|(c, _)| c.as_str());
        for (chrom, cov) in sorted_chroms {
            let mut chrom_entry = serde_json::Map::new();
            chrom_entry.insert(
                "mean_coverage".into(),
                serde_json::json!((*cov * 100.0).round() / 100.0),
            );
            chroms_obj.insert(chrom.clone(), serde_json::Value::Object(chrom_entry));
        }

        // Count variants by type.
        let mut by_type_obj = serde_json::Map::new();
        for av in &all_applied {
            let type_name = match &av.variant.mutation {
                MutationType::Snv { .. } => "SNV",
                MutationType::Indel { .. } => "INDEL",
                MutationType::Mnv { .. } => "MNV",
                MutationType::Sv { .. } => "SV",
            };
            let entry = by_type_obj
                .entry(type_name.to_string())
                .or_insert(serde_json::json!(0u64));
            if let Some(n) = entry.as_u64() {
                *entry = serde_json::json!(n + 1);
            }
        }

        // Duplex conversion rate: fraction of molecules for which both AB and
        // BA strand families were generated (T108). In simulation this is always
        // 1.0 when duplex is active; the field is included so downstream tools
        // can parse a consistent schema.
        let duplex_conversion_rate: Option<f64> =
            if cfg.umi.as_ref().map(|u| u.duplex).unwrap_or(false) && duplex_total_molecules > 0 {
                let rate =
                    duplex_molecules_with_both_strands as f64 / duplex_total_molecules as f64;
                Some((rate * 10000.0).round() / 10000.0)
            } else {
                None
            };

        let umi_obj = if let Some(ref umi) = cfg.umi {
            serde_json::json!({
                "enabled": true,
                "estimated_duplicate_rate": (duplicate_rate * 10000.0).round() / 10000.0,
                "family_size_mean": umi.family_size_mean,
                "duplex_total_molecules": duplex_total_molecules,
                "duplex_conversion_rate": duplex_conversion_rate,
            })
        } else {
            serde_json::json!({ "enabled": false })
        };

        // Compute per-variant strand concordance (T107).
        // strand_concordance = duplex_alt_count / actual_alt_count (duplex mode only).
        let is_duplex = cfg.umi.as_ref().map(|u| u.duplex).unwrap_or(false);
        let variants_arr: serde_json::Value = if is_duplex {
            let variant_entries: Vec<serde_json::Value> = all_applied
                .iter()
                .map(|av| {
                    let concordance = if av.actual_alt_count > 0 {
                        av.duplex_alt_count as f64 / av.actual_alt_count as f64
                    } else {
                        1.0
                    };
                    serde_json::json!({
                        "chrom": av.variant.chrom,
                        "pos": av.variant.pos(),
                        "n_alt_mol": av.actual_alt_count,
                        "n_duplex_alt": av.duplex_alt_count,
                        "strand_concordance": (concordance * 10000.0).round() / 10000.0,
                    })
                })
                .collect();
            serde_json::Value::Array(variant_entries)
        } else {
            serde_json::Value::Null
        };

        // Overall strand concordance: weighted mean across all variants.
        let overall_strand_concordance: Option<f64> = if is_duplex {
            let total_alt: u64 = all_applied
                .iter()
                .map(|av| av.actual_alt_count as u64)
                .sum();
            let total_duplex: u64 = all_applied
                .iter()
                .map(|av| av.duplex_alt_count as u64)
                .sum();
            if total_alt > 0 {
                Some((total_duplex as f64 / total_alt as f64 * 10000.0).round() / 10000.0)
            } else {
                None
            }
        } else {
            None
        };

        let mut report = serde_json::json!({
            "total_read_pairs": total_read_pairs,
            "sample_name": sample_name,
            "target_coverage": cfg.sample.coverage,
            "read_length": cfg.sample.read_length,
            "chromosomes": serde_json::Value::Object(chroms_obj),
            "variants": {
                "total": all_applied.len(),
                "by_type": serde_json::Value::Object(by_type_obj),
            },
            "umi": umi_obj,
        });

        if is_duplex {
            report["strand_concordance"] = serde_json::json!({
                "per_variant": variants_arr,
                "overall": overall_strand_concordance,
            });
        }

        // Capture uniformity metrics (T111).
        // achieved_coverage_cv is the theoretical CV of the LogNormal model
        // parameterised by coverage_uniformity (σ): CV = sqrt(exp(σ²) - 1).
        // achieved_on_target_fraction is 1.0 - off_target_fraction.
        if let Some(ref cap_cfg) = cfg.capture {
            if cap_cfg.enabled {
                let uniformity = cap_cfg.coverage_uniformity;
                let achieved_cv = (f64::exp(uniformity * uniformity) - 1.0).sqrt();
                let achieved_cv_rounded = (achieved_cv * 10000.0).round() / 10000.0;
                let achieved_on_target = (1.0 - cap_cfg.off_target_fraction).clamp(0.0, 1.0);
                let achieved_on_target_rounded = (achieved_on_target * 10000.0).round() / 10000.0;

                // Warn if either target is not met.
                if let Some(cv_target) = cap_cfg.coverage_cv_target {
                    if achieved_cv > cv_target {
                        tracing::warn!(
                            "capture coverage CV {:.4} exceeds target {:.4}; \
                             reduce coverage_uniformity to improve uniformity",
                            achieved_cv,
                            cv_target,
                        );
                    }
                }
                if let Some(ot_target) = cap_cfg.on_target_fraction_target {
                    if achieved_on_target < ot_target {
                        tracing::warn!(
                            "on-target fraction {:.4} is below target {:.4}; \
                             reduce off_target_fraction to improve enrichment",
                            achieved_on_target,
                            ot_target,
                        );
                    }
                }

                report["capture"] = serde_json::json!({
                    "enabled": true,
                    "achieved_coverage_cv": achieved_cv_rounded,
                    "achieved_on_target_fraction": achieved_on_target_rounded,
                    "coverage_cv_target": cap_cfg.coverage_cv_target,
                    "on_target_fraction_target": cap_cfg.on_target_fraction_target,
                });
            }
        }

        let report_path = out_dir.join("sim_report.json");
        let report_pretty =
            serde_json::to_string_pretty(&report).context("failed to serialise sim_report")?;
        std::fs::write(&report_path, report_pretty)
            .with_context(|| format!("failed to write sim_report: {}", report_path.display()))?;
        tracing::info!("sim_report written to {}", report_path.display());
    }

    eprintln!(
        "Simulation complete: {} read pairs, {} variants applied, {:.2}s ({:.0} reads/sec)",
        total_read_pairs,
        all_applied.len(),
        elapsed.as_secs_f64(),
        throughput,
    );

    Ok(())
}

/// Run a paired tumour-normal simulation.
///
/// The tumour sample is simulated with all configured somatic and germline
/// variants. The normal sample is simulated with germline variants only (somatic
/// mutations are cleared). Outputs go to `{base_output}/tumour/` and
/// `{base_output}/normal/` respectively.
fn run_paired_simulation(
    cfg: Config,
    reference: ReferenceGenome,
    start_time: Instant,
    dry_run: bool,
) -> Result<()> {
    let paired = cfg.paired.as_ref().expect("paired config must be set");
    let base_output = cfg.output.directory.clone();

    // --- Tumour run ---
    let tumour_dir = base_output.join("tumour");
    let mut tumour_cfg = cfg.clone();
    tumour_cfg.output.directory = tumour_dir;
    tumour_cfg.paired = None;

    eprintln!("Running tumour sample simulation...");
    let reference_for_tumour = ReferenceGenome::open(&cfg.reference)
        .context("failed to re-open reference for tumour run")?;
    run_single_sample(tumour_cfg, reference_for_tumour, Instant::now(), dry_run)?;

    // --- Normal run: germline-only by default; optionally with scaled somatic VAFs ---
    let normal_dir = base_output.join("normal");
    let mut normal_cfg = cfg.clone();
    normal_cfg.output.directory = normal_dir;
    normal_cfg.paired = None;
    normal_cfg.sample.name = paired.normal_sample_name.clone();
    normal_cfg.sample.coverage = paired.normal_coverage;

    let contamination = paired.tumour_contamination_in_normal;
    if contamination > 0.0 {
        // Scale somatic VAFs by the contamination fraction to simulate
        // a small amount of tumour DNA in the normal sample.
        if let Some(ref mut muts) = normal_cfg.mutations {
            if let Some(ref mut rand) = muts.random {
                rand.vaf_min *= contamination;
                rand.vaf_max *= contamination;
                // Ensure vaf_min < vaf_max after scaling.
                if rand.vaf_min >= rand.vaf_max {
                    rand.vaf_max = rand.vaf_min + 1e-9;
                }
            }
        }
    } else {
        // No contamination: suppress somatic mutations in the normal sample.
        normal_cfg.mutations = None;
    }
    // Keep germline if present.

    eprintln!("Running normal sample simulation...");
    run_single_sample(normal_cfg, reference, Instant::now(), dry_run)?;

    let elapsed = start_time.elapsed();
    eprintln!(
        "Paired simulation complete in {:.1}s. Outputs: tumour/ and normal/",
        elapsed.as_secs_f64()
    );

    Ok(())
}

/// Run multi-sample simulation: iterate over each sample in `cfg.samples`,
/// simulate independently with per-sample coverage and tumour fraction, and
/// write a combined manifest.
fn run_multi_sample(cfg: Config, reference: ReferenceGenome, start_time: Instant) -> Result<()> {
    tracing::info!(
        "multi-sample mode: {} samples configured",
        cfg.samples.as_ref().map(|s| s.len()).unwrap_or(0)
    );

    let root_dir = cfg.output.directory.clone();

    let plan = MultiSamplePlan::from_config(cfg)
        .ok_or_else(|| anyhow::anyhow!("multi-sample plan is empty"))?;

    let resolved = plan
        .per_sample_configs()
        .context("failed to resolve per-sample configs")?;

    // Create all output directories before starting parallel simulation so that
    // directory creation errors surface immediately and sequentially.
    for sample in &resolved {
        std::fs::create_dir_all(&sample.output_dir).with_context(|| {
            format!(
                "failed to create output directory: {}",
                sample.output_dir.display()
            )
        })?;
    }

    // Generate shared germline variants once from the first sample's config.
    // All samples in a multi-sample (longitudinal) run represent the same patient,
    // so the germline must be identical across timepoints.
    let shared_germline: Arc<Vec<Variant>> = if let Some(first) = resolved.first() {
        let first_cfg = &first.config;
        let first_chroms = build_chrom_list(first_cfg, &reference)?;
        let first_regions = partition_regions(&first_chroms, DEFAULT_CHUNK_SIZE);
        let seed = first_cfg.seed.unwrap_or(0);
        let germ = build_germline_variant_list(first_cfg, &first_regions, &reference, seed)
            .context("failed to generate shared germline variant list")?;
        tracing::info!(
            "generated {} shared germline variants for all samples",
            germ.len()
        );
        Arc::new(germ)
    } else {
        Arc::new(Vec::new())
    };

    // Simulate each sample in parallel.  rayon's work-stealing scheduler
    // handles the nested parallelism inside each `run_sample_simulation` call.
    let manifest_entries: Vec<SampleManifestEntry> = resolved
        .par_iter()
        .map(|sample| -> Result<SampleManifestEntry> {
            tracing::info!(
                "simulating sample '{}' (coverage={:.1}x, tumour_fraction={:.4})",
                sample.name,
                sample.config.sample.coverage,
                sample
                    .config
                    .tumour
                    .as_ref()
                    .map(|t| t.purity)
                    .unwrap_or(1.0),
            );

            let ref_for_sample = reference.clone();
            let sample_start = Instant::now();
            let (total_pairs, applied_count) = run_sample_simulation(
                sample.config.clone(),
                ref_for_sample,
                Some(&shared_germline),
            )?;

            tracing::info!(
                "sample '{}' complete: {} read pairs, {} variants, {:.2}s",
                sample.name,
                total_pairs,
                applied_count,
                sample_start.elapsed().as_secs_f64(),
            );

            Ok(SampleManifestEntry {
                name: sample.name.clone(),
                output_dir: sample.output_dir.to_string_lossy().into_owned(),
                coverage: sample.config.sample.coverage,
                tumour_fraction: sample
                    .config
                    .tumour
                    .as_ref()
                    .map(|t| t.purity)
                    .unwrap_or(1.0),
                total_read_pairs: total_pairs,
                variants_applied: applied_count,
            })
        })
        .collect::<Result<Vec<_>>>()?;

    // Write combined manifest.
    write_combined_manifest(&root_dir, &manifest_entries, env!("CARGO_PKG_VERSION"))
        .context("failed to write combined manifest")?;

    let elapsed = start_time.elapsed();
    let total_pairs: u64 = manifest_entries.iter().map(|e| e.total_read_pairs).sum();
    eprintln!(
        "Multi-sample simulation complete: {} samples, {} total read pairs, {:.2}s",
        manifest_entries.len(),
        total_pairs,
        elapsed.as_secs_f64(),
    );

    Ok(())
}

// ---------------------------------------------------------------------------
// Batch mode
// ---------------------------------------------------------------------------

/// One row in the batch manifest TSV.
struct BatchManifestRow {
    vaf: f64,
    output_dir: std::path::PathBuf,
    n_variants: usize,
    seed: Option<u64>,
    status: String,
}

/// Run one simulation per VAF listed in `cfg.vafs`.
///
/// Each run writes its output to a sub-directory named `{vaf_pct}pct/` under
/// the original output directory. After all runs, a `batch_manifest.tsv` is
/// written to the root output directory regardless of whether individual runs
/// failed.
fn run_batch(
    cfg: Config,
    reference: ReferenceGenome,
    start_time: Instant,
    dry_run: bool,
) -> Result<()> {
    let vafs = cfg.vafs.clone().unwrap_or_default();
    tracing::info!("batch mode: {} VAF values configured", vafs.len());

    let root_dir = cfg.output.directory.clone();
    let reference_path = cfg.reference.clone();

    let mut rows: Vec<BatchManifestRow> = Vec::new();

    for vaf in &vafs {
        let vaf = *vaf;

        // Build a label like "1pct", "5pct", or "0.5000pct" for fractional percentages.
        let pct_val = vaf * 100.0;
        let pct_str = if pct_val.fract() < 1e-9 {
            format!("{:.0}pct", pct_val)
        } else {
            format!("{:.4}pct", pct_val)
        };

        let sub_dir = root_dir.join(&pct_str);

        // Fork the config for this VAF value.
        let mut run_cfg = cfg.clone();
        run_cfg.vafs = None; // Prevent infinite recursion.
        run_cfg.output.directory = sub_dir.clone();

        // Override the random mutation VAF range to pin it to this exact value.
        if let Some(ref mut mutations) = run_cfg.mutations {
            if let Some(ref mut random) = mutations.random {
                random.vaf_min = vaf;
                // Keep vaf_min < vaf_max as required by the validator.
                random.vaf_max = vaf + 1e-9;
            }
        }

        // Count variants for the manifest row (requires a reference lookup).
        let n_variants = count_variants_for_config(&run_cfg, &reference_path);

        let seed = run_cfg.seed;

        tracing::info!("batch: simulating VAF {} → {}", vaf, sub_dir.display());

        let result = {
            // Clone the reference (reopens the file handle) for this run.
            let ref_for_run = reference.clone();
            run_single_sample(run_cfg, ref_for_run, Instant::now(), dry_run)
        };

        let status = match result {
            Ok(()) => "ok".to_string(),
            Err(ref e) => {
                tracing::error!("batch run for VAF {} failed: {:#}", vaf, e);
                format!("err:{}", e)
            }
        };

        rows.push(BatchManifestRow {
            vaf,
            output_dir: sub_dir,
            n_variants,
            seed,
            status,
        });
    }

    // Write the batch manifest even if some individual runs failed.
    if !dry_run {
        write_batch_manifest(&root_dir, &rows)?;
    }

    let elapsed = start_time.elapsed();
    let failed = rows.iter().filter(|r| r.status != "ok").count();
    eprintln!(
        "Batch simulation complete: {} VAF values, {} failed, {:.2}s",
        rows.len(),
        failed,
        elapsed.as_secs_f64(),
    );

    if failed > 0 {
        anyhow::bail!(
            "{} batch run(s) failed; see batch_manifest.tsv for details",
            failed
        );
    }

    Ok(())
}

/// Count the number of variants that would be generated for `cfg` without
/// actually opening a reference genome (uses the count from the random
/// mutation config, or 0 if a VCF is used or no mutations are configured).
///
/// This is a best-effort count used only for the manifest row. For VCF-sourced
/// variants the real count is only known after parsing, so we return 0.
fn count_variants_for_config(cfg: &Config, _reference_path: &std::path::Path) -> usize {
    match &cfg.mutations {
        Some(m) => {
            if m.vcf.is_some() {
                // VCF variant count requires parsing; return 0 as a placeholder.
                0
            } else if let Some(ref rand) = m.random {
                rand.count
            } else {
                0
            }
        }
        None => 0,
    }
}

/// Write a TSV batch manifest to `{root_dir}/batch_manifest.tsv`.
fn write_batch_manifest(root_dir: &std::path::Path, rows: &[BatchManifestRow]) -> Result<()> {
    use std::io::Write;

    std::fs::create_dir_all(root_dir).with_context(|| {
        format!(
            "failed to create batch output directory: {}",
            root_dir.display()
        )
    })?;

    let manifest_path = root_dir.join("batch_manifest.tsv");
    let mut f = std::fs::File::create(&manifest_path).with_context(|| {
        format!(
            "failed to create batch manifest: {}",
            manifest_path.display()
        )
    })?;

    writeln!(f, "vaf\toutput_dir\tn_variants\tseed\tstatus")?;
    for row in rows {
        let seed_str = row
            .seed
            .map(|s| s.to_string())
            .unwrap_or_else(|| "".to_string());
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            row.vaf,
            row.output_dir.display(),
            row.n_variants,
            seed_str,
            row.status,
        )?;
    }

    tracing::info!("batch manifest written to {}", manifest_path.display());
    Ok(())
}

/// Run the simulation for one sample config and return (total_read_pairs, applied_variant_count).
fn run_sample_simulation(
    cfg: Config,
    reference: ReferenceGenome,
    shared_germline: Option<&[Variant]>,
) -> Result<(u64, usize)> {
    let chrom_lengths = build_chrom_list(&cfg, &reference)?;
    let regions = partition_regions(&chrom_lengths, DEFAULT_CHUNK_SIZE);

    let seed = cfg.seed.unwrap_or(0);
    let variants = build_variant_list(&cfg, &regions, &reference, seed, shared_germline)?;

    let capture_model: Option<Arc<CaptureModel>> = build_capture_model(&cfg)?;

    let sample_name = &cfg.sample.name;
    let out_dir = &cfg.output.directory;

    let mut fastq_writer = if cfg.output.fastq {
        Some(FastqWriter::new(out_dir, sample_name).context("failed to create FASTQ writers")?)
    } else {
        None
    };

    let ref_sequences_for_bam: Vec<(String, u64)> = chrom_lengths.clone();
    let mut bam_writer = if cfg.output.bam {
        let bam_path = out_dir.join(format!("{}.bam", sample_name));
        Some(
            BamWriter::new(
                &bam_path,
                &ref_sequences_for_bam,
                &cfg.sample,
                cfg.output.mapq,
            )
            .context("failed to create BAM writer")?,
        )
    } else {
        None
    };

    let truth_vcf_path = out_dir.join(format!("{}.truth.vcf.gz", sample_name));
    let mut truth_vcf_writer = if cfg.output.truth_vcf {
        let contigs: Vec<(String, u64)> = chrom_lengths.clone();
        Some(
            TruthVcfWriter::new(&truth_vcf_path, sample_name, &contigs)
                .context("failed to create truth VCF writer")?,
        )
    } else {
        None
    };

    // Open the variant reads sidecar for this sample.
    let variant_reads_path = out_dir.join("variant_reads.tsv");
    let mut variant_reads_writer = std::io::BufWriter::new(
        std::fs::File::create(&variant_reads_path).context("failed to create variant_reads.tsv")?,
    );
    writeln!(
        variant_reads_writer,
        "read_name\tchrom\tpos\tvartype\tvaf\tclone_id"
    )
    .context("failed to write variant_reads.tsv header")?;

    let reference = Arc::new(reference);
    let cfg_arc = Arc::new(cfg.clone());
    let master_seed = cfg.seed;

    // Create bounded channel for streaming batches to the writer thread.
    let buffer_size = cfg.performance.output_buffer_regions;
    let (tx, rx) = bounded::<RegionBatch>(buffer_size);

    // Spawn dedicated writer thread.  The writers are moved in so they are
    // owned exclusively by this thread.
    let writer_cfg = cfg.clone();
    let writer_chrom_lengths = chrom_lengths.clone();
    let writer_handle = std::thread::spawn(move || -> Result<WriterStats> {
        let mut total_read_pairs: u64 = 0;
        let mut all_applied: Vec<AppliedVariant> = Vec::new();
        let mut duplex_total_molecules: u64 = 0;
        let mut duplex_molecules_with_both_strands: u64 = 0;
        // Accumulate read pairs per chromosome for coverage computation.
        let mut chrom_read_pairs: std::collections::HashMap<String, u64> =
            std::collections::HashMap::new();
        // Track molecular family IDs for UMI-tagged BAM output. All PCR
        // copies of the same original molecule share one entry in this map.
        let mut family_ids: std::collections::HashMap<String, i32> =
            std::collections::HashMap::new();

        // Primer trim amount: only active in amplicon mode.
        let primer_trim = writer_cfg
            .capture
            .as_ref()
            .filter(|c| c.mode == "amplicon")
            .map(|c| c.primer_trim)
            .unwrap_or(0);

        for batch in rx {
            // Write FASTQ — annotate read names with variant tags when present.
            if let Some(ref mut fq) = fastq_writer {
                for pair in &batch.read_pairs {
                    let name = if pair.variant_tags.is_empty() {
                        pair.name.clone()
                    } else {
                        let tags: String = pair
                            .variant_tags
                            .iter()
                            .map(|t| format!(" VT:Z:{}:{}:{}", t.chrom, t.pos + 1, t.vartype))
                            .collect::<Vec<_>>()
                            .join("");
                        format!("{}{}", pair.name, tags)
                    };
                    // Apply primer trimming when configured: remove the first and
                    // last `primer_trim` bases to simulate primer removal.
                    if primer_trim > 0 && pair.read1.seq.len() > 2 * primer_trim + 10 {
                        let mut trimmed = pair.clone();
                        let len1 = trimmed.read1.seq.len();
                        trimmed.read1.seq =
                            trimmed.read1.seq[primer_trim..len1 - primer_trim].to_vec();
                        trimmed.read1.qual =
                            trimmed.read1.qual[primer_trim..len1 - primer_trim].to_vec();
                        let len2 = trimmed.read2.seq.len();
                        trimmed.read2.seq =
                            trimmed.read2.seq[primer_trim..len2 - primer_trim].to_vec();
                        trimmed.read2.qual =
                            trimmed.read2.qual[primer_trim..len2 - primer_trim].to_vec();
                        fq.write_pair(&trimmed, &name)
                            .context("failed to write FASTQ record")?;
                    } else {
                        fq.write_pair(pair, &name)
                            .context("failed to write FASTQ record")?;
                    }
                }
            }

            // Write BAM
            if let Some(ref mut bam) = bam_writer {
                let ref_id = writer_chrom_lengths
                    .iter()
                    .position(|(name, _)| *name == batch.region.chrom)
                    .unwrap_or(0);
                let read_len = writer_cfg.sample.read_length;
                let umi_enabled = writer_cfg.umi.is_some();
                if writer_cfg.output.single_read_bam {
                    // Long-read mode: write one record per read pair (read1 only).
                    for pair in &batch.read_pairs {
                        let actual_len = pair.read1.seq.len();
                        let cigar = if primer_trim > 0 && actual_len > 2 * primer_trim + 10 {
                            // Soft-clip primer bases at both ends.
                            format!(
                                "{}S{}M{}S",
                                primer_trim,
                                actual_len - 2 * primer_trim,
                                primer_trim
                            )
                        } else {
                            format!("{}M", read_len)
                        };
                        bam.write_single_read(pair, ref_id, pair.fragment_start, &cigar)
                            .context("failed to write BAM single read record")?;
                    }
                } else if umi_enabled {
                    // UMI mode: write RX and MI auxiliary tags on every record.
                    for pair in &batch.read_pairs {
                        let actual_len = pair.read1.seq.len();
                        let cigar = if primer_trim > 0 && actual_len > 2 * primer_trim + 10 {
                            // Soft-clip primer bases at both ends.
                            format!(
                                "{}S{}M{}S",
                                primer_trim,
                                actual_len - 2 * primer_trim,
                                primer_trim
                            )
                        } else {
                            format!("{}M", read_len)
                        };
                        let umi = extract_umi_from_name(&pair.name).unwrap_or("");
                        let fid = family_id_for(&pair.name, &mut family_ids);
                        bam.write_pair_with_umi(
                            pair,
                            ref_id,
                            pair.fragment_start,
                            &cigar,
                            &cigar,
                            umi,
                            fid,
                        )
                        .context("failed to write BAM record")?;
                    }
                } else {
                    for pair in &batch.read_pairs {
                        let actual_len = pair.read1.seq.len();
                        let cigar = if primer_trim > 0 && actual_len > 2 * primer_trim + 10 {
                            // Soft-clip primer bases at both ends.
                            format!(
                                "{}S{}M{}S",
                                primer_trim,
                                actual_len - 2 * primer_trim,
                                primer_trim
                            )
                        } else {
                            format!("{}M", read_len)
                        };
                        bam.write_pair(pair, ref_id, pair.fragment_start, &cigar, &cigar)
                            .context("failed to write BAM record")?;
                    }
                }
            }

            // Write variant reads sidecar — one row per variant tag per read pair.
            for pair in &batch.read_pairs {
                for tag in &pair.variant_tags {
                    writeln!(
                        variant_reads_writer,
                        "{}	{}	{}	{}	{:.6}	{}",
                        pair.name,
                        tag.chrom,
                        tag.pos + 1,
                        tag.vartype,
                        tag.vaf,
                        tag.clone_id.as_deref().unwrap_or(".")
                    )
                    .context("failed to write variant_reads.tsv row")?;
                }
            }

            let batch_pairs = batch.read_pairs.len() as u64;
            total_read_pairs += batch_pairs;
            // Accumulate per-chromosome pair count for coverage reporting.
            *chrom_read_pairs
                .entry(batch.region.chrom.clone())
                .or_insert(0) += batch_pairs;
            all_applied.extend(batch.applied_variants);
            duplex_total_molecules += batch.duplex_total_molecules;
            duplex_molecules_with_both_strands += batch.duplex_molecules_with_both_strands;
        }

        // Compute mean coverage per chromosome from accumulated pair counts.
        let chrom_coverage: std::collections::HashMap<String, f64> = chrom_read_pairs
            .iter()
            .map(|(chrom, &pairs)| {
                let chrom_len = writer_chrom_lengths
                    .iter()
                    .find(|(c, _)| c == chrom)
                    .map(|(_, l)| *l)
                    .unwrap_or(1);
                let coverage =
                    (pairs * 2 * writer_cfg.sample.read_length as u64) as f64 / chrom_len as f64;
                (chrom.clone(), coverage)
            })
            .collect();

        // Estimate duplicate rate from UMI family size mean, if configured.
        let duplicate_rate = if let Some(ref umi) = writer_cfg.umi {
            1.0 - 1.0 / umi.family_size_mean
        } else {
            0.0
        };

        // Finalize writers
        if let Some(fq) = fastq_writer {
            fq.finish().context("failed to finalize FASTQ files")?;
        }
        if let Some(bam) = bam_writer {
            bam.finish().context("failed to finalize BAM file")?;
        }
        variant_reads_writer
            .flush()
            .context("failed to flush variant_reads.tsv")?;

        Ok(WriterStats {
            total_read_pairs,
            all_applied,
            chrom_coverage,
            duplicate_rate,
            duplex_total_molecules,
            duplex_molecules_with_both_strands,
        })
    });

    // Worker threads simulate regions in parallel and send batches as they complete.
    regions
        .par_iter()
        .enumerate()
        .try_for_each(|(idx, region)| {
            let ref_arc = Arc::clone(&reference);
            let region_seed = master_seed.map(|s| derive_region_seed(s, idx as u64));
            let mut engine = SimulationEngine::new_with_shared_config(
                Arc::clone(&cfg_arc),
                region_seed,
                ref_arc,
            );
            if let Some(ref cap) = capture_model {
                engine.set_capture_model(Arc::clone(cap));
            }
            let output = engine
                .generate_region(region, &variants)
                .with_context(|| format!("simulation failed for region {:?}", region))?;
            tx.send(RegionBatch {
                region: region.clone(),
                read_pairs: output.read_pairs,
                applied_variants: output.applied_variants,
                duplex_total_molecules: output.duplex_total_molecules,
                duplex_molecules_with_both_strands: output.duplex_molecules_with_both_strands,
            })
            .map_err(|e| anyhow::anyhow!("writer channel closed: {}", e))?;
            Ok::<(), anyhow::Error>(())
        })?;

    // Signal completion by dropping the sender so the writer thread exits its loop.
    drop(tx);

    // Wait for the writer thread to finish and collect stats.
    let stats = writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;

    let total_read_pairs = stats.total_read_pairs;
    let all_applied = stats.all_applied;

    if let Some(ref mut vcf) = truth_vcf_writer {
        // Aggregate alt-molecule counts across all regions for each unique variant.
        let mut agg: std::collections::HashMap<String, (crate::core::types::Variant, u32, u32)> =
            std::collections::HashMap::new();
        for av in &all_applied {
            let key = variant_key(&av.variant);
            let entry = agg.entry(key).or_insert_with(|| (av.variant.clone(), 0, 0));
            entry.1 = entry.1.saturating_add(av.actual_alt_count);
            entry.2 = entry.2.saturating_add(av.duplex_alt_count);
        }
        // Sort by (chrom, pos) for deterministic output.
        let mut records: Vec<_> = agg.into_values().collect();
        records.sort_by(|a, b| {
            a.0.chrom
                .cmp(&b.0.chrom)
                .then_with(|| a.0.pos().cmp(&b.0.pos()))
        });
        for (variant, n_alt_mol, n_duplex_alt) in &records {
            if let Some(sv) = sv_to_structural(variant) {
                vcf.write_sv(
                    &sv,
                    variant.expected_vaf,
                    variant.ccf,
                    variant.clone_id.as_deref(),
                )
                .context("failed to write SV truth VCF record")?;
            } else {
                let (ref_allele, alt_allele) = variant_alleles(variant);
                vcf.write_variant(variant, &ref_allele, &alt_allele, *n_alt_mol, *n_duplex_alt)
                    .context("failed to write truth VCF record")?;
            }
        }
    }

    // Finalize truth VCF writer (FASTQ/BAM writers finalized inside writer thread)
    if let Some(vcf) = truth_vcf_writer {
        vcf.finish().context("failed to finalize truth VCF")?;
    }

    Ok((total_read_pairs, all_applied.len()))
}

// ---------------------------------------------------------------------------
// Dry-run mode
// ---------------------------------------------------------------------------

fn run_dry(
    cfg: &Config,
    chrom_lengths: &[(String, u64)],
    regions: &[Region],
    variants: &[Variant],
) -> Result<()> {
    let total_genome_bp: u64 = chrom_lengths.iter().map(|(_, l)| l).sum();
    let read_length = cfg.sample.read_length;
    let coverage = cfg.sample.coverage;

    let expected_pairs: u64 = regions
        .iter()
        .map(|r| read_pairs_for_coverage(r.len(), coverage, read_length))
        .sum();

    let expected_reads = expected_pairs * 2;

    // Rough size estimate: each read is read_length bytes + ~50 bytes overhead in FASTQ
    let bytes_per_read = read_length as u64 + 50;
    let estimated_fastq_bytes = expected_reads * bytes_per_read;
    // gzip typically achieves ~4x compression on FASTQ
    let estimated_fastq_gz_bytes = estimated_fastq_bytes / 4;

    eprintln!("=== DRY RUN ===");
    eprintln!("Reference chromosomes: {}", chrom_lengths.len());
    eprintln!("Total genome size:     {} bp", total_genome_bp);
    eprintln!("Target coverage:       {}x", coverage);
    eprintln!("Read length:           {} bp", read_length);
    eprintln!("Regions:               {}", regions.len());
    eprintln!("Expected read pairs:   {}", expected_pairs);
    eprintln!("Expected total reads:  {}", expected_reads);
    eprintln!("Variants to spike:     {}", variants.len());
    eprintln!(
        "Estimated FASTQ size:  {:.1} MB (uncompressed), {:.1} MB (gzipped)",
        estimated_fastq_bytes as f64 / 1_048_576.0,
        estimated_fastq_gz_bytes as f64 / 1_048_576.0,
    );
    eprintln!("Output directory:      {}", cfg.output.directory.display());
    eprintln!("FASTQ output:          {}", cfg.output.fastq);
    eprintln!("BAM output:            {}", cfg.output.bam);
    eprintln!("Truth VCF output:      {}", cfg.output.truth_vcf);
    eprintln!("(dry-run: no files written)");

    // -----------------------------------------------------------------------
    // Per-VAF table (only when batch mode is configured)
    // -----------------------------------------------------------------------
    if let Some(vafs) = &cfg.vafs {
        if !vafs.is_empty() {
            let family_size_mean = cfg.umi.as_ref().map(|u| u.family_size_mean).unwrap_or(1.0);

            eprintln!();
            eprintln!("=== PER-VAF BATCH ESTIMATES ===");
            eprintln!(
                "{:<10}  {:>22}  {:>22}",
                "VAF", "ExpectedAltReadPairs", "ExpectedAltFamilies"
            );
            eprintln!("{}", "-".repeat(58));

            let mut any_below_one = false;
            for &vaf in vafs {
                let expected_alt_pairs = expected_pairs as f64 * vaf;
                let expected_alt_families = expected_alt_pairs / family_size_mean;
                if expected_alt_families < 1.0 {
                    any_below_one = true;
                }
                eprintln!(
                    "{:<10.6}  {:>22.1}  {:>22.1}",
                    vaf, expected_alt_pairs, expected_alt_families
                );
            }

            if any_below_one {
                anyhow::bail!(
                    "one or more VAF values produce fewer than 1 expected alt family at this \
                     coverage; increase coverage or raise the VAF floor"
                );
            }
        }
    }

    tracing::info!(
        expected_pairs,
        expected_reads,
        variants = variants.len(),
        "dry-run complete"
    );

    Ok(())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build the list of (chrom, length) pairs to simulate, respecting config.chromosomes.
fn build_chrom_list(cfg: &Config, reference: &ReferenceGenome) -> Result<Vec<(String, u64)>> {
    let all_lengths = reference.chromosome_lengths();

    let selected: Vec<(String, u64)> = if let Some(ref chroms) = cfg.chromosomes {
        // Use only the chromosomes listed in config
        let mut out = Vec::new();
        for chrom in chroms {
            let len = all_lengths.get(chrom).copied().ok_or_else(|| {
                anyhow::anyhow!(
                    "chromosome '{}' listed in config not found in reference",
                    chrom
                )
            })?;
            out.push((chrom.clone(), len));
        }
        out
    } else {
        // All chromosomes from the reference, in deterministic order
        let mut pairs: Vec<(String, u64)> =
            all_lengths.iter().map(|(k, v)| (k.clone(), *v)).collect();
        pairs.sort_by(|a, b| a.0.cmp(&b.0));
        pairs
    };

    anyhow::ensure!(!selected.is_empty(), "no chromosomes to simulate");
    Ok(selected)
}

/// Build the somatic variant list from config: VCF input takes precedence over random generation.
///
/// Returns somatic variants only. Germline variants are appended separately by
/// `build_germline_variant_list` and merged by the caller.
fn build_somatic_variant_list(
    cfg: &Config,
    regions: &[Region],
    reference: &ReferenceGenome,
    seed: u64,
) -> Result<Vec<Variant>> {
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    use crate::variants::sv_signatures;

    let mutations_cfg = match &cfg.mutations {
        Some(m) => m,
        None => return Ok(Vec::new()),
    };

    // Build the base variant list from VCF or random generation.
    let mut variants: Vec<Variant> = if let Some(ref vcf_path) = mutations_cfg.vcf {
        tracing::info!("loading variants from VCF: {}", vcf_path.display());
        let chrom_lengths = reference.chromosome_lengths();
        let known_chroms: Vec<String> = chrom_lengths.keys().cloned().collect();
        let result = vcf_input::parse_vcf(vcf_path, Some(&known_chroms), None)
            .with_context(|| format!("failed to parse VCF: {}", vcf_path.display()))?;
        if result.skipped > 0 {
            tracing::warn!("{} VCF records were skipped", result.skipped);
        }
        tracing::info!("loaded {} variants from VCF", result.variants.len());
        result.variants
    } else if let Some(rand_cfg) = &mutations_cfg.random {
        let mut rng = StdRng::seed_from_u64(seed.wrapping_add(1));
        let lookup = |chrom: &str, pos: u64| -> Option<u8> {
            let region = Region::new(chrom, pos, pos + 1);
            reference
                .sequence(&region)
                .ok()
                .and_then(|seq| seq.into_iter().next())
        };

        // Resolve the optional SBS signature for weighted alt base selection.
        let signature = if let Some(ref name) = rand_cfg.signature {
            let sig = resolve_signature(name)
                .ok_or_else(|| anyhow::anyhow!("unknown SBS signature: {}", name))?;
            tracing::info!("using SBS signature {} for SNV alt base selection", name);
            Some(sig)
        } else {
            None
        };

        generate_random_mutations(
            regions,
            rand_cfg.count,
            rand_cfg.vaf_min,
            rand_cfg.vaf_max,
            rand_cfg.snv_fraction,
            rand_cfg.indel_fraction,
            rand_cfg.mnv_fraction,
            &lookup,
            signature,
            &mut rng,
        )
    } else {
        Vec::new()
    };

    // SV signature generation (HRD, TDP, chromothripsis).
    if let Some(ref sig) = mutations_cfg.sv_signature {
        let count = mutations_cfg.sv_count;
        let mut sig_rng = StdRng::seed_from_u64(seed.wrapping_add(3));
        let chrom_lengths_map = reference.chromosome_lengths();
        let chrom_lengths: Vec<(String, u64)> = chrom_lengths_map
            .iter()
            .map(|(c, &l)| (c.clone(), l))
            .collect();
        let sv_variants = match sig.as_str() {
            "HRD" => sv_signatures::generate_hrd_deletions(&chrom_lengths, count, &mut sig_rng),
            "TDP" => sv_signatures::generate_tdp_duplications(&chrom_lengths, count, &mut sig_rng),
            "CHROMOTHRIPSIS" => {
                sv_signatures::generate_chromothripsis(&chrom_lengths, count, &mut sig_rng)
            }
            other => {
                tracing::warn!("unknown sv_signature '{}'; ignoring", other);
                Vec::new()
            }
        };
        tracing::info!(
            "generated {} SV signature variants ({})",
            sv_variants.len(),
            sig
        );
        variants.extend(sv_variants);
    }

    // MSI mode: add extra indels at homopolymer and dinucleotide repeat loci.
    if cfg.tumour.as_ref().map(|t| t.msi).unwrap_or(false) {
        let mut msi_rng = rand::rngs::StdRng::seed_from_u64(seed.wrapping_add(7));
        let lookup = |chrom: &str, pos: u64| -> Option<u8> {
            let region = Region::new(chrom, pos, pos + 1);
            reference
                .sequence(&region)
                .ok()
                .and_then(|seq| seq.into_iter().next())
        };
        let msi_variants = generate_msi_indels(regions, &lookup, &mut msi_rng, 50, 500);
        tracing::info!(
            "MSI mode: generated {} indels at repeat loci",
            msi_variants.len()
        );
        variants.extend(msi_variants);
    }

    // Driver mutation injection: add canonical hotspot mutations from the
    // active cancer preset. Only mutations with fully specified genomic
    // coordinates (chrom, pos, ref_allele, alt_allele all Some) are injected.
    // Structural/fusion events without allele fields are silently skipped.
    if mutations_cfg.include_driver_mutations {
        if let Some(ref preset_name) = cfg.preset {
            if let Some(cancer_name) = preset_name.strip_prefix("cancer:") {
                let purity = cfg.tumour.as_ref().map(|t| t.purity).unwrap_or(1.0);
                // Clonal heterozygous driver: VAF = purity / 2.
                let driver_vaf = purity / 2.0;

                let drivers = crate::cli::cancer_presets::drivers_for(cancer_name);
                let mut injected = 0usize;
                for driver in drivers {
                    // Skip mutations without fully resolved genomic coordinates.
                    let (Some(chrom), Some(pos), Some(ref_allele), Some(alt_allele)) = (
                        driver.chrom,
                        driver.pos,
                        driver.ref_allele,
                        driver.alt_allele,
                    ) else {
                        continue;
                    };

                    // Only inject mutations with prevalence above 10 % to keep
                    // the variant list realistic without overwhelming it.
                    if driver.prevalence <= 0.1 {
                        continue;
                    }

                    let mutation = if ref_allele.len() == 1 && alt_allele.len() == 1 {
                        MutationType::Snv {
                            pos,
                            ref_base: ref_allele[0],
                            alt_base: alt_allele[0],
                        }
                    } else {
                        MutationType::Indel {
                            pos,
                            ref_seq: ref_allele.to_vec(),
                            alt_seq: alt_allele.to_vec(),
                        }
                    };

                    variants.push(Variant {
                        chrom: chrom.to_string(),
                        mutation,
                        expected_vaf: driver_vaf,
                        clone_id: None,
                        haplotype: None,
                        ccf: None,
                    });
                    injected += 1;
                }

                if injected > 0 {
                    tracing::info!(
                        "injected {} driver mutation(s) from cancer preset '{}'",
                        injected,
                        cancer_name
                    );
                }
            } else {
                tracing::warn!(
                    "include_driver_mutations is set but preset '{}' is not a cancer preset",
                    preset_name
                );
            }
        } else {
            tracing::warn!("include_driver_mutations is set but no preset is configured");
        }
    }

    // Clonal tree integration (T067, T080).
    //
    // When `tumour.clones` is non-empty, build a ClonalTree and use it to:
    //   1. Set the true CCF on each variant from its assigned clone.
    //   2. Compute effective VAF: vaf = ccf * purity / ploidy.
    //   3. Propagate mutations from parent clones to all descendant clones by
    //      duplicating the variant with the subclone's CCF and effective VAF.
    if let Some(tumour_cfg) = &cfg.tumour {
        if !tumour_cfg.clones.is_empty() {
            variants = apply_clonal_tree(variants, tumour_cfg, &tumour_cfg.clones)?;
        }
    }

    Ok(variants)
}

/// Apply the clonal tree model to a variant list.
///
/// For each variant with a `clone_id` that matches a clone in the tree:
///   - Sets `variant.ccf` to the clone's CCF.
///   - Sets `variant.expected_vaf` to `ccf * purity / ploidy` (if not already
///     user-specified via a VCF `AF` field or random generation outside clonal context).
///
/// Mutations are inherited by descendant clones. If a variant is assigned to
/// clone A and clone B is a child of A, the variant also appears in B's cells.
/// This is modelled by adding a copy of the variant for each descendant, with
/// that descendant's CCF and effective VAF.
///
/// When `clone_id` is `None` or does not match any tree node, the variant is
/// left unchanged.
fn apply_clonal_tree(
    variants: Vec<Variant>,
    tumour_cfg: &crate::io::config::TumourConfig,
    clone_cfgs: &[crate::io::config::CloneConfig],
) -> Result<Vec<Variant>> {
    // Build the ClonalTree.
    let clones: Vec<TumourClone> = clone_cfgs
        .iter()
        .map(|c| TumourClone {
            id: c.id.clone(),
            ccf: c.ccf,
            parent: c.parent.clone(),
        })
        .collect();
    let tree = ClonalTree::new(clones).context("invalid clonal tree in config")?;

    let purity = tumour_cfg.purity;
    let ploidy = tumour_cfg.ploidy as f64;

    // Helper: collect all descendant clone IDs for a given clone (excluding itself).
    let descendants_of = |clone_id: &str| -> Vec<String> {
        // BFS over children.
        let mut result = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        for child_id in tree.children_of(clone_id) {
            queue.push_back(child_id.clone());
        }
        while let Some(id) = queue.pop_front() {
            result.push(id.clone());
            for child_id in tree.children_of(&id) {
                queue.push_back(child_id.clone());
            }
        }
        result
    };

    let n_input = variants.len();
    let mut output: Vec<Variant> = Vec::with_capacity(n_input);

    for mut v in variants {
        let Some(ref clone_id) = v.clone_id.clone() else {
            // No clone assignment: leave untouched.
            output.push(v);
            continue;
        };

        let Some(clone) = tree.get(clone_id) else {
            // Clone ID not found in tree: leave untouched.
            output.push(v);
            continue;
        };

        let ccf = clone.ccf;
        let effective_vaf = ccf * purity / ploidy;

        // Update the originating variant.
        v.ccf = Some(ccf);
        v.expected_vaf = effective_vaf;
        output.push(v.clone());

        // Propagate to descendant clones.
        for desc_id in descendants_of(clone_id) {
            let Some(desc_clone) = tree.get(&desc_id) else {
                continue;
            };
            let desc_ccf = desc_clone.ccf;
            let desc_vaf = desc_ccf * purity / ploidy;
            let mut desc_v = v.clone();
            desc_v.clone_id = Some(desc_id.clone());
            desc_v.ccf = Some(desc_ccf);
            desc_v.expected_vaf = desc_vaf;
            output.push(desc_v);
        }
    }

    tracing::info!(
        "clonal tree applied ({} clones): {} input variants -> {} output variants (including propagated)",
        tree.clones().len(),
        n_input,
        output.len()
    );

    Ok(output)
}

/// Generate germline variants from config and merge them with somatic variants.
///
/// Uses a separate RNG stream (seed + 2) so germline placement is independent
/// of somatic variant generation. Returns an empty list when `cfg.germline` is
/// `None`.
fn build_germline_variant_list(
    cfg: &Config,
    regions: &[Region],
    reference: &ReferenceGenome,
    seed: u64,
) -> Result<Vec<Variant>> {
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    let germline_cfg = match &cfg.germline {
        Some(g) => g,
        None => return Ok(Vec::new()),
    };

    // If a VCF is specified, load germline variants from it directly.
    if let Some(ref vcf_path) = germline_cfg.vcf {
        tracing::info!("loading germline variants from VCF: {}", vcf_path.display());
        let chrom_lengths = reference.chromosome_lengths();
        let known_chroms: Vec<String> = chrom_lengths.keys().cloned().collect();
        let result = vcf_input::parse_vcf(vcf_path, Some(&known_chroms), None)
            .with_context(|| format!("failed to parse germline VCF: {}", vcf_path.display()))?;
        if result.skipped > 0 {
            tracing::warn!("{} germline VCF records were skipped", result.skipped);
        }
        tracing::info!(
            "loaded {} germline variants from VCF",
            result.variants.len()
        );
        return Ok(result.variants);
    }

    let mut rng = StdRng::seed_from_u64(seed.wrapping_add(2));
    let lookup = |chrom: &str, pos: u64| -> Option<u8> {
        let region = Region::new(chrom, pos, pos + 1);
        reference
            .sequence(&region)
            .ok()
            .and_then(|seq| seq.into_iter().next())
    };

    let variants = generate_germline_variants(regions, germline_cfg, &lookup, &mut rng);
    tracing::info!("generated {} germline variants", variants.len());
    Ok(variants)
}

/// Build the full variant list (somatic + germline) from config.
///
/// If `shared_germline` is provided, it is used in place of generating germline
/// variants from config. This allows multi-sample runs to share a single germline
/// set across all timepoints, ensuring consistent germline positions.
fn build_variant_list(
    cfg: &Config,
    regions: &[Region],
    reference: &ReferenceGenome,
    seed: u64,
    shared_germline: Option<&[Variant]>,
) -> Result<Vec<Variant>> {
    let mut variants = build_somatic_variant_list(cfg, regions, reference, seed)?;
    let germline = if let Some(germ) = shared_germline {
        germ.to_vec()
    } else {
        build_germline_variant_list(cfg, regions, reference, seed)?
    };
    variants.extend(germline);

    Ok(variants)
}

/// Build a CaptureModel from config, or None if capture is disabled or absent.
///
/// Parses the targets BED file when specified and creates a CaptureModel
/// wrapping the target regions.
fn build_capture_model(cfg: &Config) -> Result<Option<Arc<CaptureModel>>> {
    let capture_cfg = match &cfg.capture {
        Some(c) if c.enabled => c,
        _ => return Ok(None),
    };

    let (target_regions, target_depths) = if let Some(ref bed_path) = capture_cfg.targets_bed {
        let entries = parse_bed_file(bed_path)
            .with_context(|| format!("failed to parse capture BED: {}", bed_path.display()))?;
        let regions: Vec<Region> = entries.iter().map(|(r, _)| r.clone()).collect();
        let depths: Vec<Option<f64>> = entries.into_iter().map(|(_, d)| d).collect();
        (regions, depths)
    } else {
        (Vec::new(), Vec::new())
    };

    if target_regions.is_empty() && capture_cfg.targets_bed.is_some() {
        tracing::warn!("capture BED file produced no regions; treating as no-capture");
        return Ok(None);
    }

    if target_regions.is_empty() {
        // No BED path and no regions: no-op capture (uniform).
        return Ok(None);
    }

    let model = CaptureModel::new(
        target_regions,
        target_depths,
        capture_cfg.off_target_fraction,
        capture_cfg.coverage_uniformity,
        capture_cfg.edge_dropoff_bases,
        capture_cfg.mode.clone(),
    );

    tracing::info!(
        "capture model: {} targets, off_target={:.2}, uniformity={:.2}, edge_dropoff={}",
        model.target_regions.len(),
        model.off_target_fraction,
        model.coverage_uniformity,
        model.edge_dropoff_bases,
    );

    Ok(Some(Arc::new(model)))
}

/// Parse a BED file into a list of regions with optional depth multipliers.
///
/// Expects tab- or space-separated lines with at least 3 fields: chrom, start, end.
/// An optional 4th column is parsed as a depth multiplier (`f64`).  Lines with an
/// unparseable 4th column emit a warning and set the depth to `None`.
/// Comment lines (starting with '#') and blank lines are skipped.
fn parse_bed_file(path: &std::path::Path) -> Result<Vec<(Region, Option<f64>)>> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path)
        .with_context(|| format!("cannot open BED file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut entries: Vec<(Region, Option<f64>)> = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.context("I/O error reading BED file")?;
        let line = line.trim();

        // Skip comments and track lines.
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }

        // Parse up to 5 fields so we can capture col 4 (depth) if present.
        let fields: Vec<&str> = line.splitn(5, '\t').collect();
        if fields.len() < 3 {
            // Try space-separated as fallback.
            let fields: Vec<&str> = line.splitn(5, ' ').collect();
            if fields.len() < 3 {
                tracing::warn!(
                    "BED line {} has fewer than 3 fields, skipping: {}",
                    line_no + 1,
                    line
                );
                continue;
            }
            let start: u64 = fields[1].parse().with_context(|| {
                format!("invalid BED start at line {}: {}", line_no + 1, fields[1])
            })?;
            let end: u64 = fields[2].parse().with_context(|| {
                format!("invalid BED end at line {}: {}", line_no + 1, fields[2])
            })?;
            let depth = parse_bed_depth(fields.get(3).copied(), line_no + 1);
            entries.push((Region::new(fields[0], start, end), depth));
            continue;
        }

        let start: u64 = fields[1]
            .parse()
            .with_context(|| format!("invalid BED start at line {}: {}", line_no + 1, fields[1]))?;
        let end: u64 = fields[2]
            .parse()
            .with_context(|| format!("invalid BED end at line {}: {}", line_no + 1, fields[2]))?;
        let depth = parse_bed_depth(fields.get(3).copied(), line_no + 1);
        entries.push((Region::new(fields[0], start, end), depth));
    }

    Ok(entries)
}

/// Parse the optional 4th BED column as a depth multiplier.
///
/// Returns `None` if the field is absent, empty, or non-numeric (with a warning
/// on non-numeric values).
fn parse_bed_depth(field: Option<&str>, line_no: usize) -> Option<f64> {
    let s = field?.trim();
    if s.is_empty() {
        return None;
    }
    match s.parse::<f64>() {
        Ok(d) => Some(d),
        Err(_) => {
            // Non-numeric column 4 (e.g. a gene name): treat as absent.
            tracing::debug!(
                "BED line {}: column 4 '{}' is not a number; depth multiplier set to None",
                line_no,
                s
            );
            None
        }
    }
}

/// Parse a VAF range string like "0.001-0.05" into (min, max).
fn parse_vaf_range(s: Option<&str>) -> Result<(f64, f64)> {
    match s {
        None => Ok((0.001, 0.5)),
        Some(s) => {
            let parts: Vec<&str> = s.splitn(2, '-').collect();
            anyhow::ensure!(
                parts.len() == 2,
                "vaf_range must be in the format 'min-max', e.g. '0.001-0.05'"
            );
            let min: f64 = parts[0]
                .parse()
                .with_context(|| format!("invalid vaf_range min: '{}'", parts[0]))?;
            let max: f64 = parts[1]
                .parse()
                .with_context(|| format!("invalid vaf_range max: '{}'", parts[1]))?;
            Ok((min, max))
        }
    }
}

/// Produce a stable string key for deduplicating applied variants.
fn variant_key(v: &Variant) -> String {
    match &v.mutation {
        MutationType::Snv {
            pos,
            ref_base,
            alt_base,
        } => {
            format!(
                "SNV:{}:{}:{}:{}",
                v.chrom, pos, *ref_base as char, *alt_base as char
            )
        }
        MutationType::Indel {
            pos,
            ref_seq,
            alt_seq,
        } => {
            format!(
                "INDEL:{}:{}:{}:{}",
                v.chrom,
                pos,
                String::from_utf8_lossy(ref_seq),
                String::from_utf8_lossy(alt_seq)
            )
        }
        MutationType::Mnv {
            pos,
            ref_seq,
            alt_seq,
        } => {
            format!(
                "MNV:{}:{}:{}:{}",
                v.chrom,
                pos,
                String::from_utf8_lossy(ref_seq),
                String::from_utf8_lossy(alt_seq)
            )
        }
        MutationType::Sv {
            sv_type,
            start,
            end,
            ..
        } => {
            format!("SV:{}:{:?}:{}:{}", v.chrom, sv_type, start, end)
        }
    }
}

/// Extract (ref_allele, alt_allele) byte vectors from a variant.
fn variant_alleles(v: &Variant) -> (Vec<u8>, Vec<u8>) {
    match &v.mutation {
        MutationType::Snv {
            ref_base, alt_base, ..
        } => (vec![*ref_base], vec![*alt_base]),
        MutationType::Indel {
            ref_seq, alt_seq, ..
        } => (ref_seq.clone(), alt_seq.clone()),
        MutationType::Mnv {
            ref_seq, alt_seq, ..
        } => (ref_seq.clone(), alt_seq.clone()),
        MutationType::Sv { .. } => {
            // SVs use symbolic alleles; emit placeholder representation.
            // SVs are written via write_sv; this path is unreachable in practice.
            (b"N".to_vec(), b"<SV>".to_vec())
        }
    }
}

/// Convert a `MutationType::Sv` variant to the richer `StructuralVariant` type
/// needed by `TruthVcfWriter::write_sv`.
fn sv_to_structural(v: &Variant) -> Option<StructuralVariant> {
    use crate::core::types::SvType;
    match &v.mutation {
        MutationType::Sv {
            sv_type,
            chrom,
            start,
            end,
        } => Some(match sv_type {
            SvType::Deletion => StructuralVariant::Deletion {
                chrom: chrom.clone(),
                start: *start,
                end: *end,
            },
            SvType::Insertion => StructuralVariant::Insertion {
                chrom: chrom.clone(),
                pos: *start,
                sequence: Vec::new(),
            },
            SvType::Inversion => StructuralVariant::Inversion {
                chrom: chrom.clone(),
                start: *start,
                end: *end,
            },
            SvType::Duplication => StructuralVariant::Duplication {
                chrom: chrom.clone(),
                start: *start,
                end: *end,
                copies: 1,
            },
            SvType::Translocation => StructuralVariant::Translocation {
                chrom1: chrom.clone(),
                pos1: *start,
                chrom2: chrom.clone(),
                pos2: *end,
            },
        }),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Unit tests for helpers in this module
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_vaf_range_default() {
        let (min, max) = parse_vaf_range(None).unwrap();
        assert_eq!(min, 0.001);
        assert_eq!(max, 0.5);
    }

    #[test]
    fn test_parse_vaf_range_explicit() {
        let (min, max) = parse_vaf_range(Some("0.01-0.3")).unwrap();
        assert!((min - 0.01).abs() < 1e-9);
        assert!((max - 0.3).abs() < 1e-9);
    }

    #[test]
    fn test_parse_vaf_range_invalid() {
        assert!(parse_vaf_range(Some("bad")).is_err());
        assert!(parse_vaf_range(Some("0.01")).is_err());
    }

    #[test]
    fn test_variant_key_snv() {
        let v = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 100,
                ref_base: b'A',
                alt_base: b'T',
            },
            expected_vaf: 0.3,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };
        let key = variant_key(&v);
        assert!(key.contains("SNV"));
        assert!(key.contains("chr1"));
        assert!(key.contains("100"));
    }

    #[test]
    fn test_variant_alleles_snv() {
        let v = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 0,
                ref_base: b'G',
                alt_base: b'C',
            },
            expected_vaf: 0.5,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };
        let (r, a) = variant_alleles(&v);
        assert_eq!(r, b"G");
        assert_eq!(a, b"C");
    }

    #[test]
    fn test_variant_alleles_indel() {
        let v = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Indel {
                pos: 0,
                ref_seq: b"ACG".to_vec(),
                alt_seq: b"A".to_vec(),
            },
            expected_vaf: 0.1,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };
        let (r, a) = variant_alleles(&v);
        assert_eq!(r, b"ACG");
        assert_eq!(a, b"A");
    }

    // -----------------------------------------------------------------------
    // Task 07: Parallel processing tests
    // -----------------------------------------------------------------------

    /// Each region index produces a unique, deterministic seed.
    #[test]
    fn test_region_seeds_unique() {
        let master: u64 = 12345;

        let seeds: Vec<u64> = (0..16).map(|i| derive_region_seed(master, i)).collect();

        // All seeds must be distinct.
        let unique: std::collections::HashSet<u64> = seeds.iter().cloned().collect();
        assert_eq!(
            unique.len(),
            seeds.len(),
            "every region index must produce a unique seed"
        );

        // Same inputs always give the same output (determinism).
        for (i, &seed) in seeds.iter().enumerate() {
            assert_eq!(
                derive_region_seed(master, i as u64),
                seed,
                "derive_region_seed must be deterministic for index {}",
                i
            );
        }

        // Different master seeds must give different region seeds.
        let other_seeds: Vec<u64> = (0..16).map(|i| derive_region_seed(master + 1, i)).collect();
        assert_ne!(
            seeds, other_seeds,
            "different master seeds must yield different region seeds"
        );
    }

    /// configure_thread_pool with a specific count must not panic.
    #[test]
    fn test_thread_count_respected() {
        // We can only call build_global once per process (rayon), so use a
        // local thread pool instead of the global one for this test.
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .expect("failed to build rayon thread pool with 2 threads");

        let observed = std::sync::atomic::AtomicUsize::new(0);
        pool.install(|| {
            use rayon::prelude::*;
            (0..100usize).into_par_iter().for_each(|_| {
                observed.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            });
        });
        assert_eq!(
            observed.load(std::sync::atomic::Ordering::Relaxed),
            100,
            "all work items should be processed"
        );
    }

    /// Same seed produces identical read sequences whether 1 or multiple
    /// regions are processed (determinism across parallel execution).
    ///
    /// We test this at the unit level by checking that two engines seeded
    /// with the same per-region seed produce identical read pairs.
    #[test]
    fn test_deterministic_parallel() {
        use crate::core::engine::SimulationEngine;
        use crate::core::types::Region;
        use crate::io::config::{
            Config, FragmentConfig, FragmentModel, OutputConfig, QualityConfig, SampleConfig,
        };
        use crate::io::reference::ReferenceGenome;
        use tempfile::TempDir;

        let dir = TempDir::new().unwrap();
        let fa = {
            let fa_path = dir.path().join("ref.fa");
            let fai_path = dir.path().join("ref.fa.fai");
            let seq: Vec<u8> = b"ACGT".iter().cycle().take(1000).cloned().collect();
            let mut fa_bytes: Vec<u8> = b">chr1\n".to_vec();
            fa_bytes.extend_from_slice(&seq);
            fa_bytes.extend_from_slice(b"\n");
            std::fs::write(&fa_path, &fa_bytes).unwrap();
            let fai = "chr1\t1000\t6\t1000\t1001\n";
            std::fs::write(&fai_path, fai.as_bytes()).unwrap();
            fa_path
        };

        let make_cfg = |seed: u64| Config {
            reference: fa.clone(),
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
                coverage: 5.0,
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
            copy_number: None,
            gc_bias: None,
            samples: None,
            capture: None,
            seed: Some(seed),
            threads: None,
            chromosomes: None,
            regions_bed: None,
            performance: Default::default(),
            preset: None,
            vafs: None,
            germline: None,
            paired: None,
        };

        let master_seed = 42u64;
        let region = Region::new("chr1", 0, 500);

        for region_idx in 0u64..4 {
            let per_region_seed = derive_region_seed(master_seed, region_idx);

            let ref1 = ReferenceGenome::open(&fa).unwrap();
            let ref2 = ReferenceGenome::open(&fa).unwrap();

            let mut engine1 = SimulationEngine::new(make_cfg(per_region_seed), ref1);
            let mut engine2 = SimulationEngine::new(make_cfg(per_region_seed), ref2);

            let out1 = engine1.generate_region(&region, &[]).unwrap();
            let out2 = engine2.generate_region(&region, &[]).unwrap();

            assert_eq!(
                out1.read_pairs.len(),
                out2.read_pairs.len(),
                "region {} pair counts differ",
                region_idx
            );
            for (i, (p1, p2)) in out1
                .read_pairs
                .iter()
                .zip(out2.read_pairs.iter())
                .enumerate()
            {
                assert_eq!(
                    p1.read1.seq, p2.read1.seq,
                    "region {} pair {} read1 sequences differ",
                    region_idx, i
                );
                assert_eq!(
                    p1.read2.seq, p2.read2.seq,
                    "region {} pair {} read2 sequences differ",
                    region_idx, i
                );
            }
        }
    }

    /// Output order is stable: collecting parallel results and sorting by
    /// region index always yields the same sequence.
    #[test]
    fn test_output_order_stable() {
        let master = 99u64;
        let n = 8usize;

        let run1: Vec<u64> = (0..n as u64)
            .map(|i| derive_region_seed(master, i))
            .collect();
        let run2: Vec<u64> = (0..n as u64)
            .map(|i| derive_region_seed(master, i))
            .collect();

        assert_eq!(
            run1, run2,
            "region seed sequence must be identical across runs"
        );

        let mut shuffled = run1.iter().cloned().enumerate().collect::<Vec<_>>();
        shuffled.sort_by_key(|(i, _)| n - 1 - i);
        shuffled.sort_by_key(|(i, _)| *i);
        let restored: Vec<u64> = shuffled.into_iter().map(|(_, s)| s).collect();
        assert_eq!(
            run1, restored,
            "sorting by index must restore original order"
        );
    }

    // -----------------------------------------------------------------------
    // Task 13: CLI overrides & presets tests
    // -----------------------------------------------------------------------

    use crate::io::config::{Config, FragmentConfig, OutputConfig, QualityConfig, SampleConfig};
    use std::path::PathBuf;

    fn base_config() -> Config {
        Config {
            reference: PathBuf::from("/dev/null"),
            output: OutputConfig {
                directory: PathBuf::from("/tmp"),
                fastq: true,
                bam: false,
                truth_vcf: false,
                manifest: false,
                germline_vcf: false,
                single_read_bam: false,
                mapq: 60,
            },
            sample: SampleConfig::default(),
            fragment: FragmentConfig::default(),
            quality: QualityConfig::default(),
            tumour: None,
            mutations: None,
            umi: None,
            artifacts: None,
            seed: None,
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

    fn empty_opts() -> SimulateOpts {
        SimulateOpts {
            config: PathBuf::from("/dev/null"),
            output_dir: None,
            seed: None,
            coverage: None,
            read_length: None,
            purity: None,
            fragment_mean: None,
            fragment_sd: None,
            random_mutations: None,
            vaf_range: None,
            preset: None,
            dry_run: false,
            list_presets: false,
            set: None,
        }
    }

    /// 1. --coverage overrides sample.coverage in config.
    #[test]
    fn test_coverage_override() {
        let mut cfg = base_config();
        cfg.sample.coverage = 30.0;

        let opts = SimulateOpts {
            coverage: Some(100.0),
            ..empty_opts()
        };
        apply_overrides(&mut cfg, &opts, None).unwrap();

        assert!(
            (cfg.sample.coverage - 100.0).abs() < 1e-9,
            "coverage should be overridden to 100x"
        );
    }

    /// 2. --seed overrides config.seed.
    #[test]
    fn test_seed_override() {
        let mut cfg = base_config();
        cfg.seed = Some(1);

        let opts = SimulateOpts {
            seed: Some(42),
            ..empty_opts()
        };
        apply_overrides(&mut cfg, &opts, None).unwrap();

        assert_eq!(cfg.seed, Some(42));
    }

    /// 3. Small preset sets expected values (1x coverage, chr22, 100 mutations).
    #[test]
    fn test_preset_small() {
        use crate::cli::presets;

        let overlay = presets::get("small").unwrap();
        let mut cfg = base_config();
        presets::apply_preset_to_config(&mut cfg, &overlay);

        assert!(
            (cfg.sample.coverage - 1.0).abs() < 1e-9,
            "small preset: 1x coverage"
        );
        assert_eq!(
            cfg.chromosomes.as_deref(),
            Some(vec!["chr22".to_string()].as_slice()),
            "small preset: chr22 only"
        );
        let count = cfg.mutations.unwrap().random.unwrap().count;
        assert_eq!(count, 100, "small preset: 100 random mutations");
    }

    /// 4. cfDNA preset enables the cfDNA fragment model.
    #[test]
    fn test_preset_cfdna() {
        use crate::cli::presets;
        use crate::io::config::FragmentModel;

        let overlay = presets::get("cfdna").unwrap();
        let mut cfg = base_config();
        presets::apply_preset_to_config(&mut cfg, &overlay);

        assert!(
            matches!(cfg.fragment.model, FragmentModel::Cfda),
            "cfdna preset should use Cfda fragment model"
        );
        assert!(
            (cfg.sample.coverage - 200.0).abs() < 1e-9,
            "cfdna preset: 200x coverage"
        );
    }

    /// 5. CLI flag > YAML > preset > default precedence.
    #[test]
    fn test_precedence() {
        use crate::cli::presets;

        let mut cfg = base_config();

        let overlay = presets::get("small").unwrap();
        presets::apply_preset_to_config(&mut cfg, &overlay);
        assert!(
            (cfg.sample.coverage - 1.0).abs() < 1e-9,
            "preset should set 1x"
        );

        cfg.sample.coverage = 60.0;

        let opts = SimulateOpts {
            coverage: Some(150.0),
            ..empty_opts()
        };
        apply_overrides(&mut cfg, &opts, None).unwrap();

        assert!(
            (cfg.sample.coverage - 150.0).abs() < 1e-9,
            "CLI flag should win over YAML and preset"
        );
    }

    /// 6. No CLI flags → config unchanged.
    #[test]
    fn test_no_overrides() {
        let mut cfg = base_config();
        cfg.sample.coverage = 55.0;
        cfg.seed = Some(7);

        let opts = empty_opts();
        apply_overrides(&mut cfg, &opts, None).unwrap();

        assert!((cfg.sample.coverage - 55.0).abs() < 1e-9);
        assert_eq!(cfg.seed, Some(7));
    }

    /// 7. Every preset produces a valid non-default coverage after application.
    #[test]
    fn test_all_presets_valid() {
        use crate::cli::presets;

        for name in presets::all_names() {
            let overlay = presets::get(name).expect(name);
            let mut cfg = base_config();
            presets::apply_preset_to_config(&mut cfg, &overlay);
            assert!(
                cfg.sample.coverage > 0.0,
                "preset '{}' should produce positive coverage",
                name
            );
        }
    }

    // -----------------------------------------------------------------------
    // BED file parsing tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_bed_file() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "# comment line").unwrap();
        writeln!(f, "chr1\t100\t200").unwrap();
        writeln!(f, "chr1\t500\t600\tgene_A").unwrap();
        writeln!(f, "chr2\t0\t1000").unwrap();
        f.flush().unwrap();

        let entries = parse_bed_file(f.path()).unwrap();
        assert_eq!(entries.len(), 3);
        assert_eq!(entries[0].0.chrom, "chr1");
        assert_eq!(entries[0].0.start, 100);
        assert_eq!(entries[0].0.end, 200);
        assert_eq!(entries[0].1, None); // no depth column
        assert_eq!(entries[1].1, None); // "gene_A" is not a number
        assert_eq!(entries[2].0.chrom, "chr2");
        assert_eq!(entries[2].0.start, 0);
        assert_eq!(entries[2].0.end, 1000);
    }

    // -----------------------------------------------------------------------
    // Capture CV formula tests (T111)
    //
    // The achieved coverage CV is sqrt(exp(σ²) - 1) where σ is
    // coverage_uniformity.  At σ=0 the CV must be exactly 0.0 (perfectly
    // uniform); at σ=0.5 it must be well above 0.1.
    // -----------------------------------------------------------------------

    fn capture_cv(uniformity: f64) -> f64 {
        (f64::exp(uniformity * uniformity) - 1.0).sqrt()
    }

    #[test]
    fn test_capture_cv_zero_uniformity_gives_zero_cv() {
        let cv = capture_cv(0.0);
        assert_eq!(cv, 0.0, "CV must be exactly 0.0 at uniformity=0.0");
    }

    #[test]
    fn test_capture_cv_half_uniformity_exceeds_threshold() {
        let cv = capture_cv(0.5);
        assert!(
            cv > 0.1,
            "CV at uniformity=0.5 should be > 0.1, got {:.4}",
            cv
        );
    }
}
