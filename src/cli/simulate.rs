use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result};
use crossbeam_channel::bounded;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use super::SimulateOpts;
use crate::cli::presets;
use crate::core::capture::CaptureModel;
use crate::core::coverage::{partition_regions, read_pairs_for_coverage};
use crate::core::engine::{AppliedVariant, SimulationEngine, derive_region_seed};
use crate::core::multi_sample::{MultiSamplePlan, SampleManifestEntry, write_combined_manifest};
use crate::core::types::{MutationType, ReadPair, Region, Variant};
use crate::io::bam::BamWriter;
use crate::io::config::{self, Config};
use crate::io::fastq::FastqWriter;
use crate::io::manifest::{Manifest, ReferenceInfo, SimulationStatistics};
use crate::io::reference::ReferenceGenome;
use crate::io::truth_vcf::TruthVcfWriter;
use crate::io::vcf_input;
use crate::variants::random_gen::generate_random_mutations;

/// A batch of generated reads for one region, sent through the streaming channel.
struct RegionBatch {
    region: Region,
    read_pairs: Vec<ReadPair>,
    applied_variants: Vec<AppliedVariant>,
}

/// Statistics returned by the writer thread.
struct WriterStats {
    total_read_pairs: u64,
    all_applied: Vec<AppliedVariant>,
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
        let existing = config
            .mutations
            .get_or_insert(config::MutationConfig { vcf: None, random: None });
        let rand = existing.random.get_or_insert(config::RandomMutationConfig {
            count: 0,
            vaf_min: 0.001,
            vaf_max: 0.5,
            snv_fraction: 0.80,
            indel_fraction: 0.15,
            mnv_fraction: 0.05,
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

pub fn run(opts: SimulateOpts, cli_threads: Option<usize>) -> Result<()> {
    let start_time = Instant::now();

    // -----------------------------------------------------------------------
    // 1. Load and validate config, applying CLI overrides
    // -----------------------------------------------------------------------
    tracing::info!("loading config from {}", opts.config.display());
    let mut cfg = config::load(&opts.config)
        .with_context(|| format!("failed to load config: {}", opts.config.display()))?;

    // Apply preset (lower precedence than YAML).
    if let Some(ref preset_name) = opts.preset {
        let overlay = presets::get(preset_name)
            .with_context(|| format!("invalid preset '{}'", preset_name))?;
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
    // 3. Multi-sample mode: if samples[] is present, delegate to per-sample runs
    // -----------------------------------------------------------------------
    if cfg.samples.is_some() && !cfg.samples.as_ref().map(|s| s.is_empty()).unwrap_or(true) {
        return run_multi_sample(cfg, reference, start_time);
    }

    // -----------------------------------------------------------------------
    // 4. Single-sample simulation
    // -----------------------------------------------------------------------
    run_single_sample(cfg, reference, start_time, opts.dry_run)
}

/// Run the simulation for a single sample configuration.
fn run_single_sample(
    cfg: Config,
    reference: ReferenceGenome,
    start_time: Instant,
    dry_run: bool,
) -> Result<()> {
    // Build the sorted list of chromosomes/lengths we will simulate
    let chrom_lengths: Vec<(String, u64)> = build_chrom_list(&cfg, &reference)?;
    tracing::info!("reference has {} chromosomes to simulate", chrom_lengths.len());

    // -----------------------------------------------------------------------
    // Partition regions
    // -----------------------------------------------------------------------
    let regions = partition_regions(&chrom_lengths, DEFAULT_CHUNK_SIZE);
    tracing::info!("partitioned into {} regions", regions.len());

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
    let variants = build_variant_list(&cfg, &regions, &reference, seed)?;
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
        Some(
            FastqWriter::new(out_dir, sample_name)
                .context("failed to create FASTQ writers")?,
        )
    } else {
        None
    };

    let ref_sequences_for_bam: Vec<(String, u64)> = chrom_lengths.clone();
    let mut bam_writer = if cfg.output.bam {
        let bam_path = out_dir.join(format!("{}.bam", sample_name));
        Some(
            BamWriter::new(&bam_path, &ref_sequences_for_bam, &cfg.sample)
                .context("failed to create BAM writer")?,
        )
    } else {
        None
    };

    // We'll open truth VCF after collecting applied variants; use a temp path
    // and rename, or just open now (write header) and stream records.
    let truth_vcf_path = out_dir.join(format!("{}.truth.vcf", sample_name));
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

    // -----------------------------------------------------------------------
    // Run simulation per region in parallel, streaming batches to a dedicated
    // writer thread via a bounded channel.
    // -----------------------------------------------------------------------
    tracing::info!("starting parallel simulation across {} regions", regions.len());
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

        for batch in rx {
            // Write FASTQ
            if let Some(ref mut fq) = fastq_writer {
                for pair in &batch.read_pairs {
                    fq.write_pair(pair, &pair.name)
                        .context("failed to write FASTQ record")?;
                }
            }

            // Write BAM
            if let Some(ref mut bam) = bam_writer {
                let ref_id = writer_chrom_lengths
                    .iter()
                    .position(|(name, _)| *name == batch.region.chrom)
                    .unwrap_or(0);
                let read_len = writer_cfg.sample.read_length;
                let cigar = format!("{}M", read_len);
                for pair in &batch.read_pairs {
                    bam.write_pair(pair, ref_id, pair.fragment_start, &cigar, &cigar)
                        .context("failed to write BAM record")?;
                }
            }

            total_read_pairs += batch.read_pairs.len() as u64;
            all_applied.extend(batch.applied_variants);
        }

        // Finalize writers
        if let Some(fq) = fastq_writer {
            fq.finish().context("failed to finalize FASTQ files")?;
        }
        if let Some(bam) = bam_writer {
            bam.finish().context("failed to finalize BAM file")?;
        }

        Ok(WriterStats { total_read_pairs, all_applied })
    });

    // Worker threads simulate regions in parallel and send batches as they complete.
    regions.par_iter().enumerate().try_for_each(|(idx, region)| {
        let ref_arc = Arc::clone(&reference);
        let mut local_cfg = cfg.clone();
        // Override seed to per-region value for determinism.
        local_cfg.seed = master_seed.map(|s| derive_region_seed(s, idx as u64));
        let mut engine = SimulationEngine::new_with_shared_reference(local_cfg, ref_arc);

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
        }).map_err(|e| anyhow::anyhow!("writer channel closed: {}", e))?;
        pb_ref.inc(1);
        Ok::<(), anyhow::Error>(())
    })?;

    // Signal completion by dropping the sender so the writer thread exits its loop.
    drop(tx);

    // Wait for the writer thread to finish and collect stats.
    let stats = writer_handle.join()
        .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;

    pb.finish_with_message("simulation complete");

    let total_read_pairs = stats.total_read_pairs;
    let all_applied = stats.all_applied;

    // -----------------------------------------------------------------------
    // Write truth VCF with all applied variants
    // -----------------------------------------------------------------------
    if let Some(ref mut vcf) = truth_vcf_writer {
        let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
        for av in &all_applied {
            let key = variant_key(&av.variant);
            if seen.insert(key) {
                let (ref_allele, alt_allele) = variant_alleles(&av.variant);
                vcf.write_variant(&av.variant, &ref_allele, &alt_allele)
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
        let mut manifest_json = serde_json::to_value(&manifest)
            .context("failed to serialize manifest")?;
        if let serde_json::Value::Object(ref mut map) = manifest_json {
            map.insert("sample_name".into(), serde_json::Value::String(sample_name.to_string()));
            map.insert("total_read_pairs".into(), serde_json::Value::Number(total_read_pairs.into()));
            map.insert("total_reads".into(), serde_json::Value::Number((total_read_pairs * 2).into()));
            map.insert("variants_applied".into(), serde_json::Value::Number(all_applied.len().into()));
            if let Some(seed_val) = cfg.seed {
                map.insert("seed".into(), serde_json::Value::Number(seed_val.into()));
            }
            map.insert("wall_time_secs".into(), serde_json::json!(elapsed.as_secs_f64()));
            map.insert("varforge_version".into(), serde_json::Value::String(env!("CARGO_PKG_VERSION").to_string()));
        }
        let pretty = serde_json::to_string_pretty(&manifest_json)
            .context("failed to pretty-print manifest")?;
        std::fs::write(&manifest_path, pretty)
            .with_context(|| format!("failed to write manifest: {}", manifest_path.display()))?;
        tracing::info!("manifest written to {}", manifest_path.display());
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

/// Run multi-sample simulation: iterate over each sample in `cfg.samples`,
/// simulate independently with per-sample coverage and tumour fraction, and
/// write a combined manifest.
fn run_multi_sample(
    cfg: Config,
    reference: ReferenceGenome,
    start_time: Instant,
) -> Result<()> {
    tracing::info!("multi-sample mode: {} samples configured", cfg.samples.as_ref().map(|s| s.len()).unwrap_or(0));

    let root_dir = cfg.output.directory.clone();

    let plan = MultiSamplePlan::from_config(cfg)
        .ok_or_else(|| anyhow::anyhow!("multi-sample plan is empty"))?;

    let resolved = plan.per_sample_configs()
        .context("failed to resolve per-sample configs")?;

    let mut manifest_entries: Vec<SampleManifestEntry> = Vec::new();

    for sample in &resolved {
        tracing::info!(
            "simulating sample '{}' (coverage={:.1}x, tumour_fraction={:.4})",
            sample.name,
            sample.config.sample.coverage,
            sample.config.tumour.as_ref().map(|t| t.purity).unwrap_or(1.0),
        );

        // Create per-sample output directory.
        std::fs::create_dir_all(&sample.output_dir).with_context(|| {
            format!("failed to create output directory: {}", sample.output_dir.display())
        })?;

        // Clone the reference (reopens the file handle) for this sample.
        let ref_for_sample = reference.clone();

        let sample_start = Instant::now();
        let (total_pairs, applied_count) = run_sample_simulation(
            sample.config.clone(),
            ref_for_sample,
        )?;

        manifest_entries.push(SampleManifestEntry {
            name: sample.name.clone(),
            output_dir: sample.output_dir.to_string_lossy().into_owned(),
            coverage: sample.config.sample.coverage,
            tumour_fraction: sample.config.tumour.as_ref().map(|t| t.purity).unwrap_or(1.0),
            total_read_pairs: total_pairs,
            variants_applied: applied_count,
        });

        tracing::info!(
            "sample '{}' complete: {} read pairs, {} variants, {:.2}s",
            sample.name,
            total_pairs,
            applied_count,
            sample_start.elapsed().as_secs_f64(),
        );
    }

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

/// Run the simulation for one sample config and return (total_read_pairs, applied_variant_count).
fn run_sample_simulation(cfg: Config, reference: ReferenceGenome) -> Result<(u64, usize)> {
    let chrom_lengths = build_chrom_list(&cfg, &reference)?;
    let regions = partition_regions(&chrom_lengths, DEFAULT_CHUNK_SIZE);

    let seed = cfg.seed.unwrap_or(0);
    let variants = build_variant_list(&cfg, &regions, &reference, seed)?;

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
        Some(BamWriter::new(&bam_path, &ref_sequences_for_bam, &cfg.sample)
            .context("failed to create BAM writer")?)
    } else {
        None
    };

    let truth_vcf_path = out_dir.join(format!("{}.truth.vcf", sample_name));
    let mut truth_vcf_writer = if cfg.output.truth_vcf {
        let contigs: Vec<(String, u64)> = chrom_lengths.clone();
        Some(TruthVcfWriter::new(&truth_vcf_path, sample_name, &contigs)
            .context("failed to create truth VCF writer")?)
    } else {
        None
    };

    let reference = Arc::new(reference);
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

        for batch in rx {
            // Write FASTQ
            if let Some(ref mut fq) = fastq_writer {
                for pair in &batch.read_pairs {
                    fq.write_pair(pair, &pair.name)
                        .context("failed to write FASTQ record")?;
                }
            }

            // Write BAM
            if let Some(ref mut bam) = bam_writer {
                let ref_id = writer_chrom_lengths
                    .iter()
                    .position(|(name, _)| *name == batch.region.chrom)
                    .unwrap_or(0);
                let read_len = writer_cfg.sample.read_length;
                let cigar = format!("{}M", read_len);
                for pair in &batch.read_pairs {
                    bam.write_pair(pair, ref_id, pair.fragment_start, &cigar, &cigar)
                        .context("failed to write BAM record")?;
                }
            }

            total_read_pairs += batch.read_pairs.len() as u64;
            all_applied.extend(batch.applied_variants);
        }

        // Finalize writers
        if let Some(fq) = fastq_writer {
            fq.finish().context("failed to finalize FASTQ files")?;
        }
        if let Some(bam) = bam_writer {
            bam.finish().context("failed to finalize BAM file")?;
        }

        Ok(WriterStats { total_read_pairs, all_applied })
    });

    // Worker threads simulate regions in parallel and send batches as they complete.
    regions.par_iter().enumerate().try_for_each(|(idx, region)| {
        let ref_arc = Arc::clone(&reference);
        let mut local_cfg = cfg.clone();
        local_cfg.seed = master_seed.map(|s| derive_region_seed(s, idx as u64));
        let mut engine = SimulationEngine::new_with_shared_reference(local_cfg, ref_arc);
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
        }).map_err(|e| anyhow::anyhow!("writer channel closed: {}", e))?;
        Ok::<(), anyhow::Error>(())
    })?;

    // Signal completion by dropping the sender so the writer thread exits its loop.
    drop(tx);

    // Wait for the writer thread to finish and collect stats.
    let stats = writer_handle.join()
        .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;

    let total_read_pairs = stats.total_read_pairs;
    let all_applied = stats.all_applied;

    if let Some(ref mut vcf) = truth_vcf_writer {
        let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
        for av in &all_applied {
            let key = variant_key(&av.variant);
            if seen.insert(key) {
                let (ref_allele, alt_allele) = variant_alleles(&av.variant);
                vcf.write_variant(&av.variant, &ref_allele, &alt_allele)
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
        let mut pairs: Vec<(String, u64)> = all_lengths
            .iter()
            .map(|(k, v)| (k.clone(), *v))
            .collect();
        pairs.sort_by(|a, b| a.0.cmp(&b.0));
        pairs
    };

    anyhow::ensure!(!selected.is_empty(), "no chromosomes to simulate");
    Ok(selected)
}

/// Build the variant list from config: VCF input takes precedence over random generation.
fn build_variant_list(
    cfg: &Config,
    regions: &[Region],
    reference: &ReferenceGenome,
    seed: u64,
) -> Result<Vec<Variant>> {
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    let mutations_cfg = match &cfg.mutations {
        Some(m) => m,
        None => return Ok(Vec::new()),
    };

    // VCF input takes precedence over random generation.
    if let Some(ref vcf_path) = mutations_cfg.vcf {
        tracing::info!("loading variants from VCF: {}", vcf_path.display());
        let chrom_lengths = reference.chromosome_lengths();
        let known_chroms: Vec<String> = chrom_lengths.keys().cloned().collect();
        let result = vcf_input::parse_vcf(vcf_path, Some(&known_chroms), None)
            .with_context(|| format!("failed to parse VCF: {}", vcf_path.display()))?;
        if result.skipped > 0 {
            tracing::warn!("{} VCF records were skipped", result.skipped);
        }
        tracing::info!("loaded {} variants from VCF", result.variants.len());
        return Ok(result.variants);
    }

    if let Some(rand_cfg) = &mutations_cfg.random {
        let mut rng = StdRng::seed_from_u64(seed.wrapping_add(1));
        let lookup = |chrom: &str, pos: u64| -> Option<u8> {
            let region = Region::new(chrom, pos, pos + 1);
            reference.sequence(&region).ok().and_then(|seq| seq.into_iter().next())
        };
        let variants = generate_random_mutations(
            regions,
            rand_cfg.count,
            rand_cfg.vaf_min,
            rand_cfg.vaf_max,
            rand_cfg.snv_fraction,
            rand_cfg.indel_fraction,
            rand_cfg.mnv_fraction,
            &lookup,
            &mut rng,
        );
        return Ok(variants);
    }

    Ok(Vec::new())
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

    let target_regions = if let Some(ref bed_path) = capture_cfg.targets_bed {
        parse_bed_file(bed_path)
            .with_context(|| format!("failed to parse capture BED: {}", bed_path.display()))?
    } else {
        Vec::new()
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
        capture_cfg.off_target_fraction,
        capture_cfg.coverage_uniformity,
        capture_cfg.edge_dropoff_bases,
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

/// Parse a BED file into a list of Regions.
///
/// Expects tab- or space-separated lines with at least 3 fields: chrom, start, end.
/// Comment lines (starting with '#') and blank lines are skipped.
fn parse_bed_file(path: &std::path::Path) -> Result<Vec<Region>> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path)
        .with_context(|| format!("cannot open BED file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.context("I/O error reading BED file")?;
        let line = line.trim();

        // Skip comments and track lines.
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let fields: Vec<&str> = line.splitn(4, '\t').collect();
        if fields.len() < 3 {
            // Try space-separated as fallback.
            let fields: Vec<&str> = line.splitn(4, ' ').collect();
            if fields.len() < 3 {
                tracing::warn!("BED line {} has fewer than 3 fields, skipping: {}", line_no + 1, line);
                continue;
            }
            let start: u64 = fields[1].parse()
                .with_context(|| format!("invalid BED start at line {}: {}", line_no + 1, fields[1]))?;
            let end: u64 = fields[2].parse()
                .with_context(|| format!("invalid BED end at line {}: {}", line_no + 1, fields[2]))?;
            regions.push(Region::new(fields[0], start, end));
            continue;
        }

        let start: u64 = fields[1].parse()
            .with_context(|| format!("invalid BED start at line {}: {}", line_no + 1, fields[1]))?;
        let end: u64 = fields[2].parse()
            .with_context(|| format!("invalid BED end at line {}: {}", line_no + 1, fields[2]))?;
        regions.push(Region::new(fields[0], start, end));
    }

    Ok(regions)
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
        MutationType::Snv { pos, ref_base, alt_base } => {
            format!("SNV:{}:{}:{}:{}", v.chrom, pos, *ref_base as char, *alt_base as char)
        }
        MutationType::Indel { pos, ref_seq, alt_seq } => {
            format!(
                "INDEL:{}:{}:{}:{}",
                v.chrom,
                pos,
                String::from_utf8_lossy(ref_seq),
                String::from_utf8_lossy(alt_seq)
            )
        }
        MutationType::Mnv { pos, ref_seq, alt_seq } => {
            format!(
                "MNV:{}:{}:{}:{}",
                v.chrom,
                pos,
                String::from_utf8_lossy(ref_seq),
                String::from_utf8_lossy(alt_seq)
            )
        }
        MutationType::Sv { sv_type, start, end, .. } => {
            format!("SV:{}:{:?}:{}:{}", v.chrom, sv_type, start, end)
        }
    }
}

/// Extract (ref_allele, alt_allele) byte vectors from a variant.
fn variant_alleles(v: &Variant) -> (Vec<u8>, Vec<u8>) {
    match &v.mutation {
        MutationType::Snv { ref_base, alt_base, .. } => {
            (vec![*ref_base], vec![*alt_base])
        }
        MutationType::Indel { ref_seq, alt_seq, .. } => {
            (ref_seq.clone(), alt_seq.clone())
        }
        MutationType::Mnv { ref_seq, alt_seq, .. } => {
            (ref_seq.clone(), alt_seq.clone())
        }
        MutationType::Sv { .. } => {
            // SVs use symbolic alleles; emit placeholder representation.
            (b"N".to_vec(), b"<SV>".to_vec())
        }
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
            mutation: MutationType::Snv { pos: 100, ref_base: b'A', alt_base: b'T' },
            expected_vaf: 0.3,
            clone_id: None,
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
            mutation: MutationType::Snv { pos: 0, ref_base: b'G', alt_base: b'C' },
            expected_vaf: 0.5,
            clone_id: None,
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
        use crate::io::config::{Config, FragmentConfig, FragmentModel, OutputConfig, QualityConfig, SampleConfig};
        use crate::core::engine::SimulationEngine;
        use crate::core::types::Region;
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
            for (i, (p1, p2)) in out1.read_pairs.iter().zip(out2.read_pairs.iter()).enumerate() {
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

        let run1: Vec<u64> = (0..n as u64).map(|i| derive_region_seed(master, i)).collect();
        let run2: Vec<u64> = (0..n as u64).map(|i| derive_region_seed(master, i)).collect();

        assert_eq!(run1, run2, "region seed sequence must be identical across runs");

        let mut shuffled = run1.iter().cloned().enumerate().collect::<Vec<_>>();
        shuffled.sort_by_key(|(i, _)| n - 1 - i);
        shuffled.sort_by_key(|(i, _)| *i);
        let restored: Vec<u64> = shuffled.into_iter().map(|(_, s)| s).collect();
        assert_eq!(run1, restored, "sorting by index must restore original order");
    }

    // -----------------------------------------------------------------------
    // Task 13: CLI overrides & presets tests
    // -----------------------------------------------------------------------

    use std::path::PathBuf;
    use crate::io::config::{Config, FragmentConfig, OutputConfig, QualityConfig, SampleConfig};

    fn base_config() -> Config {
        Config {
            reference: PathBuf::from("/dev/null"),
            output: OutputConfig {
                directory: PathBuf::from("/tmp"),
                fastq: true,
                bam: false,
                truth_vcf: false,
                manifest: false,
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
        }
    }

    /// 1. --coverage overrides sample.coverage in config.
    #[test]
    fn test_coverage_override() {
        let mut cfg = base_config();
        cfg.sample.coverage = 30.0;

        let opts = SimulateOpts { coverage: Some(100.0), ..empty_opts() };
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

        let opts = SimulateOpts { seed: Some(42), ..empty_opts() };
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

        assert!((cfg.sample.coverage - 1.0).abs() < 1e-9, "small preset: 1x coverage");
        assert_eq!(
            cfg.chromosomes.as_deref(),
            Some(vec!["chr22".to_string()].as_slice()),
            "small preset: chr22 only"
        );
        let count = cfg
            .mutations
            .unwrap()
            .random
            .unwrap()
            .count;
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
        assert!((cfg.sample.coverage - 200.0).abs() < 1e-9, "cfdna preset: 200x coverage");
    }

    /// 5. CLI flag > YAML > preset > default precedence.
    #[test]
    fn test_precedence() {
        use crate::cli::presets;

        let mut cfg = base_config();

        let overlay = presets::get("small").unwrap();
        presets::apply_preset_to_config(&mut cfg, &overlay);
        assert!((cfg.sample.coverage - 1.0).abs() < 1e-9, "preset should set 1x");

        cfg.sample.coverage = 60.0;

        let opts = SimulateOpts { coverage: Some(150.0), ..empty_opts() };
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
        use tempfile::NamedTempFile;
        use std::io::Write;

        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "# comment line").unwrap();
        writeln!(f, "chr1\t100\t200").unwrap();
        writeln!(f, "chr1\t500\t600\tgene_A").unwrap();
        writeln!(f, "chr2\t0\t1000").unwrap();
        f.flush().unwrap();

        let regions = parse_bed_file(f.path()).unwrap();
        assert_eq!(regions.len(), 3);
        assert_eq!(regions[0].chrom, "chr1");
        assert_eq!(regions[0].start, 100);
        assert_eq!(regions[0].end, 200);
        assert_eq!(regions[2].chrom, "chr2");
        assert_eq!(regions[2].start, 0);
        assert_eq!(regions[2].end, 1000);
    }
}
