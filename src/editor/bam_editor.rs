//! Core BAM editing engine.
//!
//! Takes an existing BAM file, spikes in variants at specified positions with
//! stochastic VAF sampling, and writes a modified BAM file.
//!
//! The engine reads all records into memory, sorts them by
//! (reference_sequence_id, alignment_start), then processes them with a
//! single-pass sliding-window approach. Sorting up front makes the editor
//! robust to unsorted or partially sorted input, which the simulator does
//! not guarantee. This is acceptable because the BAM editor targets small,
//! targeted BAMs rather than whole-genome files.

use std::collections::HashMap;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::{record::data::field::Tag, RecordBuf},
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::core::types::{MutationType, Variant};
use crate::editor::read_modifier::{apply_deletion, apply_insertion, apply_snv, ModifyResult};
use crate::variants::vaf::sample_alt_count;

/// Configuration for the BAM editor.
#[derive(Debug, Clone)]
pub struct EditConfig {
    /// Path to the input BAM file.
    pub input_bam: std::path::PathBuf,
    /// Path to the output BAM file.
    pub output_bam: std::path::PathBuf,
    /// Variants to spike in.
    pub variants: Vec<Variant>,
    /// Random seed for reproducibility.
    pub seed: u64,
    /// Tumour purity (0.0-1.0). Scales all VAFs.
    pub purity: f64,
    /// Optional path to write a truth VCF.
    pub truth_vcf: Option<std::path::PathBuf>,
    /// Sample name for output headers.
    pub sample_name: String,
}

/// Result of spiking a single variant into the BAM.
#[derive(Debug, Clone)]
pub struct SpikedVariant {
    pub variant: Variant,
    pub total_depth: u32,
    pub alt_count: u32,
    pub actual_vaf: f64,
}

/// Core BAM editing engine.
pub struct BamEditor {
    config: EditConfig,
    rng: StdRng,
}

impl BamEditor {
    /// Create a new `BamEditor` from the given config.
    pub fn new(config: EditConfig) -> Self {
        let rng = StdRng::seed_from_u64(config.seed);
        Self { config, rng }
    }

    /// Run the editing pipeline:
    ///
    /// 1. Read all records from the input BAM into memory.
    /// 2. Sort records by (reference_sequence_id, alignment_start).
    /// 3. Sort variants by (chrom_id, pos).
    /// 4. Process records in order with a sliding window for active variants.
    /// 5. Once all reads for a variant region are buffered, spike that variant.
    /// 6. Flush writes as soon as records are no longer needed by any variant.
    /// 7. Optionally write truth VCF.
    pub fn run(&mut self) -> Result<Vec<SpikedVariant>> {
        let input_path = self.config.input_bam.clone();
        let output_path = self.config.output_bam.clone();
        let variants = self.config.variants.clone();
        let purity = self.config.purity;

        // Open input BAM and read header.
        let file = std::fs::File::open(&input_path)
            .with_context(|| format!("failed to open BAM file: {}", input_path.display()))?;
        let mut reader = bam::io::Reader::new(file);
        let header = reader
            .read_header()
            .with_context(|| format!("failed to read BAM header: {}", input_path.display()))?;

        // Open output BAM writer.
        let out_file = std::fs::File::create(&output_path)
            .with_context(|| format!("failed to create BAM file: {}", output_path.display()))?;
        let mut writer = bam::io::Writer::new(out_file);
        writer
            .write_header(&header)
            .context("failed to write BAM header")?;

        // Sort variants by (chrom_id, pos). Variants on unknown chromosomes are
        // dropped with a warning.
        let mut sorted_variants: Vec<(usize, u64, Variant)> = Vec::new();
        for v in &variants {
            match header
                .reference_sequences()
                .get_index_of(v.chrom.as_bytes())
            {
                Some(rid) => {
                    sorted_variants.push((rid, v.pos(), v.clone()));
                }
                None => {
                    tracing::warn!(
                        "chromosome '{}' not found in BAM header, skipping variant at pos {}",
                        v.chrom,
                        v.pos()
                    );
                }
            }
        }
        sorted_variants.sort_by_key(|&(rid, pos, _)| (rid, pos));

        // Read all records into memory, then sort by (rid, alignment_start).
        // This makes the editor robust to unsorted input from the simulator.
        let mut all_records: Vec<RecordBuf> = Vec::new();
        for result in reader.record_bufs(&header) {
            let record = result.context("failed to read BAM record")?;
            all_records.push(record);
        }
        all_records.sort_by(|a, b| {
            let rid_a = a.reference_sequence_id();
            let rid_b = b.reference_sequence_id();
            let start_a = a.alignment_start().map(|p| usize::from(p) as u64);
            let start_b = b.alignment_start().map(|p| usize::from(p) as u64);
            rid_a.cmp(&rid_b).then(start_a.cmp(&start_b))
        });

        // Process sorted records with a sliding-window approach.
        let mut spiked: Vec<SpikedVariant> = Vec::new();
        let mut variant_idx = 0usize;

        // Buffer holds records until all variants that could overlap them have
        // been processed. Modifications happen in place.
        let mut buffer: Vec<RecordBuf> = Vec::new();

        // `min_active` returns the (rid, pos) of the earliest unprocessed variant.
        // Records that cannot overlap any remaining variant can be flushed.
        let min_active = |vi: usize, sorted: &[(usize, u64, Variant)]| -> Option<(usize, u64)> {
            sorted.get(vi).map(|(rid, pos, _)| (*rid, *pos))
        };

        use noodles_sam::alignment::io::Write as _;

        for record in all_records {
            // Determine this record's alignment span.
            let rec_rid = record.reference_sequence_id();
            let rec_start = record.alignment_start().map(|p| usize::from(p) as u64 - 1);

            // Push into buffer unconditionally. We will flush eligible records
            // after processing any variants that have become fully covered.
            buffer.push(record);

            // Process variants whose region is fully covered by the buffer. A
            // variant at (vrid, vpos) is fully covered when the most recently
            // read record has alignment start > vpos on the same or a later
            // chromosome. At that point, no future record can overlap vpos.
            loop {
                let Some((vrid, vpos, _)) = sorted_variants.get(variant_idx) else {
                    break;
                };
                let vrid = *vrid;
                let vpos = *vpos;

                // Check whether the current record is past this variant.
                let past = match (rec_rid, rec_start) {
                    (Some(rr), Some(rs)) => rr > vrid || (rr == vrid && rs > vpos),
                    _ => false,
                };

                if !past {
                    // The current record may still overlap this variant; stop
                    // processing variants for now.
                    break;
                }

                // All reads for this variant have been buffered. Spike it.
                let variant = sorted_variants[variant_idx].2.clone();
                let effective_vaf = variant.expected_vaf * purity;
                if let Some(sv) =
                    self.spike_variant_in_buffer(&mut buffer, &header, &variant, effective_vaf)?
                {
                    spiked.push(sv);
                }

                variant_idx += 1;
            }

            // Flush records that cannot overlap any remaining variant.
            let min = min_active(variant_idx, &sorted_variants);
            let flush_up_to = find_flush_boundary(&buffer, min);

            for record in buffer.drain(..flush_up_to) {
                writer
                    .write_alignment_record(&header, &record)
                    .context("failed to write BAM record")?;
            }
        }

        // End of input: process all remaining variants against the buffer.
        while variant_idx < sorted_variants.len() {
            let variant = sorted_variants[variant_idx].2.clone();
            let effective_vaf = variant.expected_vaf * purity;
            if let Some(sv) =
                self.spike_variant_in_buffer(&mut buffer, &header, &variant, effective_vaf)?
            {
                spiked.push(sv);
            }
            variant_idx += 1;
        }

        // Flush all remaining buffered records.
        for record in &buffer {
            writer
                .write_alignment_record(&header, record)
                .context("failed to write BAM record")?;
        }

        writer.try_finish().context("failed to finalise BAM file")?;

        // Write truth VCF if requested.
        if let Some(ref vcf_path) = self.config.truth_vcf.clone() {
            let contigs: Vec<(String, u64)> = header
                .reference_sequences()
                .iter()
                .map(|(name, map)| {
                    let name = String::from_utf8_lossy(name).into_owned();
                    let len = usize::from(map.length()) as u64;
                    (name, len)
                })
                .collect();

            write_truth_vcf(vcf_path, &self.config.sample_name, &contigs, &spiked)?;
        }

        Ok(spiked)
    }

    /// Spike a single variant into records held in the buffer.
    ///
    /// Only records that overlap the variant position are eligible. The method
    /// modifies records in place.
    fn spike_variant_in_buffer(
        &mut self,
        buffer: &mut [RecordBuf],
        header: &sam::Header,
        variant: &Variant,
        effective_vaf: f64,
    ) -> Result<Option<SpikedVariant>> {
        let var_pos = variant.pos();
        let chrom = &variant.chrom;

        let ref_id = match header.reference_sequences().get_index_of(chrom.as_bytes()) {
            Some(id) => id,
            None => {
                tracing::warn!("chromosome '{}' not found in BAM header", chrom);
                return Ok(None);
            }
        };

        // Collect indices of overlapping records in the buffer.
        let overlapping: Vec<usize> = buffer
            .iter()
            .enumerate()
            .filter_map(|(i, record)| {
                let rid = record.reference_sequence_id()?;
                if rid != ref_id {
                    return None;
                }
                let align_start = usize::from(record.alignment_start()?) as u64 - 1;
                let read_len = record.sequence().len() as u64;
                let align_end = align_start + read_len;

                // Skip MQ=0 reads.
                let mq = record.mapping_quality().map(u8::from).unwrap_or(0);
                if mq == 0 {
                    return None;
                }

                if var_pos >= align_start && var_pos < align_end {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();

        let total_depth = overlapping.len() as u32;
        if total_depth == 0 {
            tracing::debug!("variant at {}:{} has zero depth, skipping", chrom, var_pos);
            return Ok(None);
        }

        // Sample alt count using binomial distribution (stochastic VAF).
        let alt_count = sample_alt_count(total_depth, effective_vaf, &mut self.rng);
        if alt_count == 0 {
            return Ok(Some(SpikedVariant {
                variant: variant.clone(),
                total_depth,
                alt_count: 0,
                actual_vaf: 0.0,
            }));
        }

        // Select which reads to modify, respecting strand balance.
        let to_modify =
            select_reads_for_modification(buffer, &overlapping, alt_count as usize, &mut self.rng);

        // Apply the modification to selected reads.
        let mut actual_modified: u32 = 0;
        for idx in to_modify {
            let record = &mut buffer[idx];
            let result = match &variant.mutation {
                MutationType::Snv { .. } => apply_snv(record, &variant.mutation),
                MutationType::Indel {
                    ref_seq, alt_seq, ..
                } => {
                    if alt_seq.len() > ref_seq.len() {
                        apply_insertion(record, &variant.mutation)
                    } else {
                        apply_deletion(record, &variant.mutation)
                    }
                }
                _ => ModifyResult::Unchanged,
            };
            if result == ModifyResult::Modified {
                actual_modified += 1;
            }
        }

        let actual_vaf = if total_depth > 0 {
            actual_modified as f64 / total_depth as f64
        } else {
            0.0
        };

        Ok(Some(SpikedVariant {
            variant: variant.clone(),
            total_depth,
            alt_count: actual_modified,
            actual_vaf,
        }))
    }
}

// ---------------------------------------------------------------------------
// Flush boundary helper
// ---------------------------------------------------------------------------

/// Return the number of leading buffer entries that can be safely flushed.
///
/// A record can be flushed if it cannot overlap any remaining variant. Given
/// the earliest remaining variant at `min_active` (rid, pos), a record is safe
/// to flush if:
///   - its reference id is less than min_active.rid, or
///   - its reference id equals min_active.rid and its alignment end <= min_active.pos.
///
/// When there are no remaining variants, all records are safe to flush.
fn find_flush_boundary(buffer: &[RecordBuf], min_active: Option<(usize, u64)>) -> usize {
    let Some((min_rid, min_pos)) = min_active else {
        // No variants left; flush everything.
        return buffer.len();
    };

    let mut count = 0;
    for record in buffer {
        let rid = match record.reference_sequence_id() {
            Some(r) => r,
            None => {
                // Unmapped record; safe to flush.
                count += 1;
                continue;
            }
        };
        let align_start = match record.alignment_start() {
            Some(p) => usize::from(p) as u64 - 1,
            None => {
                count += 1;
                continue;
            }
        };
        let read_len = record.sequence().len() as u64;
        let align_end = align_start + read_len;

        let safe = rid < min_rid || (rid == min_rid && align_end <= min_pos);
        if safe {
            count += 1;
        } else {
            // Records are ordered; once one is unsafe, all later ones may be too.
            break;
        }
    }
    count
}

// ---------------------------------------------------------------------------
// Read selection
// ---------------------------------------------------------------------------

/// Select `n` read indices to modify, respecting strand balance.
fn select_reads_for_modification(
    records: &[RecordBuf],
    overlapping: &[usize],
    n: usize,
    rng: &mut StdRng,
) -> Vec<usize> {
    if n == 0 || overlapping.is_empty() {
        return Vec::new();
    }

    let n = n.min(overlapping.len());

    // Separate forward and reverse reads.
    let mut fwd: Vec<usize> = Vec::new();
    let mut rev: Vec<usize> = Vec::new();
    for &idx in overlapping {
        if records[idx].flags().is_reverse_complemented() {
            rev.push(idx);
        } else {
            fwd.push(idx);
        }
    }

    // Proportion of forward/reverse reads.
    let total = overlapping.len();
    let n_fwd = ((n as f64 * fwd.len() as f64 / total as f64).round() as usize).min(fwd.len());
    let n_rev = (n - n_fwd).min(rev.len());
    let n_fwd = n - n_rev; // re-balance in case n_rev was capped

    let selected_fwd = reservoir_sample(&fwd, n_fwd, rng);
    let selected_rev = reservoir_sample(&rev, n_rev, rng);

    let mut selected: Vec<usize> = selected_fwd.into_iter().chain(selected_rev).collect();
    selected.sort_unstable();
    selected
}

/// Reservoir sampling without replacement.
fn reservoir_sample(items: &[usize], k: usize, rng: &mut StdRng) -> Vec<usize> {
    if k == 0 || items.is_empty() {
        return Vec::new();
    }
    let k = k.min(items.len());
    let mut result: Vec<usize> = items[..k].to_vec();
    for (i, &item) in items[k..].iter().enumerate() {
        let j = rng.gen_range(0..=(i + k));
        if j < k {
            result[j] = item;
        }
    }
    result
}

// ---------------------------------------------------------------------------
// BAM I/O helpers (kept for tests)
// ---------------------------------------------------------------------------

/// Load all records from a BAM file into memory.
///
/// Used in tests. Production code uses the streaming path in `BamEditor::run`.
#[cfg(test)]
pub fn load_bam(path: &Path) -> Result<(sam::Header, Vec<RecordBuf>)> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("failed to open BAM file: {}", path.display()))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader
        .read_header()
        .with_context(|| format!("failed to read BAM header: {}", path.display()))?;

    let mut records = Vec::new();
    for result in reader.record_bufs(&header) {
        let record = result.with_context(|| "failed to read BAM record")?;
        records.push(record);
    }

    Ok((header, records))
}

/// Write records to a BAM file.
#[cfg(test)]
pub fn write_bam(path: &Path, header: &sam::Header, records: &[RecordBuf]) -> Result<()> {
    use noodles_sam::alignment::io::Write as _;

    let file = std::fs::File::create(path)
        .with_context(|| format!("failed to create BAM file: {}", path.display()))?;
    let mut writer = bam::io::Writer::new(file);
    writer
        .write_header(header)
        .context("failed to write BAM header")?;

    for record in records {
        writer
            .write_alignment_record(header, record)
            .context("failed to write BAM record")?;
    }

    writer.try_finish().context("failed to finalize BAM file")?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Truth VCF writer
// ---------------------------------------------------------------------------

/// Write a truth VCF for the spiked variants.
fn write_truth_vcf(
    path: &Path,
    sample_name: &str,
    contigs: &[(String, u64)],
    spiked: &[SpikedVariant],
) -> Result<()> {
    use std::io::Write;

    let file = std::fs::File::create(path)
        .with_context(|| format!("failed to create truth VCF: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "##fileformat=VCFv4.3")?;
    writeln!(w, "##source=VarForge-edit")?;
    writeln!(
        w,
        r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth at variant site">"#
    )?;
    writeln!(
        w,
        r#"##INFO=<ID=AC,Number=1,Type=Integer,Description="Actual alt count after spike-in">"#
    )?;
    writeln!(
        w,
        r#"##INFO=<ID=AF,Number=1,Type=Float,Description="Actual allele frequency after spike-in">"#
    )?;
    writeln!(
        w,
        r#"##INFO=<ID=EXPECTED_VAF,Number=1,Type=Float,Description="Expected variant allele frequency">"#
    )?;
    writeln!(
        w,
        r#"##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Variant type">"#
    )?;
    for (name, length) in contigs {
        writeln!(w, "##contig=<ID={name},length={length}>")?;
    }
    writeln!(
        w,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
    )?;
    writeln!(
        w,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"
    )?;

    for sv in spiked {
        let v = &sv.variant;
        let (ref_allele, alt_allele) = variant_alleles(v);
        let pos_1based = v.pos() + 1;
        let vartype = v.vartype();
        let info = format!(
            "DP={dp};AC={ac};AF={af:.6};EXPECTED_VAF={evaf:.6};VARTYPE={vartype}",
            dp = sv.total_depth,
            ac = sv.alt_count,
            af = sv.actual_vaf,
            evaf = v.expected_vaf,
        );
        writeln!(
            w,
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t0/1",
            chrom = v.chrom,
            pos = pos_1based,
            ref = String::from_utf8_lossy(&ref_allele),
            alt = String::from_utf8_lossy(&alt_allele),
        )?;
    }

    w.flush().context("failed to flush truth VCF")?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Variant helpers
// ---------------------------------------------------------------------------

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
        MutationType::Sv { .. } => (b"N".to_vec(), b"<SV>".to_vec()),
    }
}

// ---------------------------------------------------------------------------
// UMI-aware editing helpers
// ---------------------------------------------------------------------------

/// Group record indices by UMI family (MI tag).
///
/// Returns a map from family_id to Vec<record_index>.
// Called only in tests; future UMI-collapse logic will use this in production.
#[allow(dead_code)]
pub fn group_by_umi_family(records: &[RecordBuf]) -> HashMap<i64, Vec<usize>> {
    let mut families: HashMap<i64, Vec<usize>> = HashMap::new();
    for (i, record) in records.iter().enumerate() {
        if let Some(mi) = record.data().get(&Tag::UMI_ID).and_then(|v| v.as_int()) {
            families.entry(mi).or_default().push(i);
        }
    }
    families
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position as NoodlesPosition;
    use noodles_sam::alignment::{
        record::Flags,
        record_buf::{Cigar, Data, QualityScores, Sequence},
        RecordBuf,
    };
    use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    use std::num::NonZeroUsize;
    use tempfile::NamedTempFile;

    // ------------------------------------------------------------------
    // Helpers to build minimal in-memory BAM data
    // ------------------------------------------------------------------

    fn build_header(refs: &[(&str, u64)]) -> sam::Header {
        let mut builder = sam::Header::builder();
        for (name, len) in refs {
            let nz_len = NonZeroUsize::new(*len as usize).unwrap();
            builder = builder.add_reference_sequence(*name, Map::<ReferenceSequence>::new(nz_len));
        }
        builder.build()
    }

    fn make_record(
        seq: &[u8],
        qual: &[u8],
        ref_id: usize,
        align_start: u64, // 0-based
        flags: Flags,
    ) -> RecordBuf {
        let cigar_str = format!("{}M", seq.len());
        let cigar_ops = crate::io::bam::parse_cigar(&cigar_str).unwrap();
        RecordBuf::builder()
            .set_flags(flags)
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(NoodlesPosition::new(align_start as usize + 1).unwrap())
            .set_mapping_quality(noodles_sam::alignment::record::MappingQuality::new(60).unwrap())
            .set_cigar(cigar_ops.into_iter().collect::<Cigar>())
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(qual.to_vec()))
            .build()
    }

    fn make_records_at_pos(
        n: usize,
        ref_id: usize,
        align_start: u64,
        seq: &[u8],
    ) -> Vec<RecordBuf> {
        (0..n)
            .map(|_| {
                make_record(
                    seq,
                    &vec![30u8; seq.len()],
                    ref_id,
                    align_start,
                    Flags::empty(),
                )
            })
            .collect()
    }

    // ------------------------------------------------------------------
    // Test 1: SNV spike-in, base changed, quality preserved
    // ------------------------------------------------------------------
    #[test]
    fn test_snv_spike_in_bam() {
        let header = build_header(&[("chr1", 1000)]);
        // 20 reads of "ACGTACGTACGTACGT" at position 100.
        let seq = b"ACGTACGTACGTACGTACGT";
        let records = make_records_at_pos(20, 0, 100, seq);

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        write_bam(tmp_in.path(), &header, &records).unwrap();

        let variant = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 102,
                ref_base: b'G',
                alt_base: b'T',
            },
            expected_vaf: 1.0, // 100% so all reads are modified
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let config = EditConfig {
            input_bam: tmp_in.path().to_path_buf(),
            output_bam: tmp_out.path().to_path_buf(),
            variants: vec![variant],
            seed: 42,
            purity: 1.0,
            truth_vcf: None,
            sample_name: "TEST".to_string(),
        };

        let mut editor = BamEditor::new(config);
        let spiked = editor.run().unwrap();

        assert_eq!(spiked.len(), 1);
        assert!(
            spiked[0].alt_count > 0,
            "at least some reads should be modified"
        );

        let (_, out_records) = load_bam(tmp_out.path()).unwrap();
        // At least some reads should have 'T' at offset 2.
        let alt_count = out_records
            .iter()
            .filter(|r| {
                let s: Vec<u8> = r.sequence().as_ref().to_vec();
                s.get(2) == Some(&b'T')
            })
            .count();
        assert!(alt_count > 0, "SNV should be applied to at least one read");
    }

    // ------------------------------------------------------------------
    // Test 2: Stochastic VAF - 50% VAF at 1000x depth should be ~50%
    // ------------------------------------------------------------------
    #[test]
    fn test_snv_stochastic_vaf() {
        let header = build_header(&[("chr1", 10000)]);
        let seq = b"ACGTACGTACGTACGTACGT"; // 20 bases
        let records = make_records_at_pos(1000, 0, 100, seq);

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        write_bam(tmp_in.path(), &header, &records).unwrap();

        let variant = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 102,
                ref_base: b'G',
                alt_base: b'T',
            },
            expected_vaf: 0.5,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let config = EditConfig {
            input_bam: tmp_in.path().to_path_buf(),
            output_bam: tmp_out.path().to_path_buf(),
            variants: vec![variant],
            seed: 42,
            purity: 1.0,
            truth_vcf: None,
            sample_name: "TEST".to_string(),
        };

        let mut editor = BamEditor::new(config);
        let spiked = editor.run().unwrap();

        assert_eq!(spiked.len(), 1);
        let sv = &spiked[0];
        // Should be roughly 50% +/- 10%.
        assert!(
            sv.actual_vaf >= 0.40 && sv.actual_vaf <= 0.60,
            "VAF at 50% target should be in [0.40, 0.60], got {}",
            sv.actual_vaf
        );
    }

    // ------------------------------------------------------------------
    // Test 3: Insertion spike-in
    // ------------------------------------------------------------------
    #[test]
    fn test_indel_insertion_bam() {
        let header = build_header(&[("chr1", 1000)]);
        let seq = b"ACGTACGTACGTACGT";
        let records = make_records_at_pos(10, 0, 100, seq);

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        write_bam(tmp_in.path(), &header, &records).unwrap();

        let variant = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Indel {
                pos: 102,
                ref_seq: b"G".to_vec(),
                alt_seq: b"GTT".to_vec(),
            },
            expected_vaf: 1.0,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let config = EditConfig {
            input_bam: tmp_in.path().to_path_buf(),
            output_bam: tmp_out.path().to_path_buf(),
            variants: vec![variant],
            seed: 42,
            purity: 1.0,
            truth_vcf: None,
            sample_name: "TEST".to_string(),
        };

        let mut editor = BamEditor::new(config);
        let spiked = editor.run().unwrap();

        assert!(spiked[0].alt_count > 0, "insertion should be applied");

        let (_, out_records) = load_bam(tmp_out.path()).unwrap();
        let has_insertion = out_records.iter().any(|r| {
            let ops: Vec<_> = r.cigar().as_ref().to_vec();
            ops.iter()
                .any(|op| op.kind() == noodles_sam::alignment::record::cigar::op::Kind::Insertion)
        });
        assert!(
            has_insertion,
            "at least one record should have an Insertion CIGAR op"
        );
    }

    // ------------------------------------------------------------------
    // Test 4: Deletion spike-in
    // ------------------------------------------------------------------
    #[test]
    fn test_indel_deletion_bam() {
        let header = build_header(&[("chr1", 1000)]);
        let seq = b"ACGTACGTACGTACGT";
        let records = make_records_at_pos(10, 0, 100, seq);

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        write_bam(tmp_in.path(), &header, &records).unwrap();

        let variant = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Indel {
                pos: 102,
                ref_seq: b"GT".to_vec(),
                alt_seq: b"G".to_vec(),
            },
            expected_vaf: 1.0,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let config = EditConfig {
            input_bam: tmp_in.path().to_path_buf(),
            output_bam: tmp_out.path().to_path_buf(),
            variants: vec![variant],
            seed: 42,
            purity: 1.0,
            truth_vcf: None,
            sample_name: "TEST".to_string(),
        };

        let mut editor = BamEditor::new(config);
        let spiked = editor.run().unwrap();
        assert!(spiked[0].alt_count > 0, "deletion should be applied");

        let (_, out_records) = load_bam(tmp_out.path()).unwrap();
        let has_deletion = out_records.iter().any(|r| {
            let ops: Vec<_> = r.cigar().as_ref().to_vec();
            ops.iter()
                .any(|op| op.kind() == noodles_sam::alignment::record::cigar::op::Kind::Deletion)
        });
        assert!(
            has_deletion,
            "at least one record should have a Deletion CIGAR op"
        );
    }

    // ------------------------------------------------------------------
    // Test 5: NM tag updated
    // ------------------------------------------------------------------
    #[test]
    fn test_nm_tag_updated() {
        use noodles_sam::alignment::record::data::field::Tag;
        let header = build_header(&[("chr1", 1000)]);
        let seq = b"ACGTACGT";
        // Build records with NM=0.
        let records: Vec<RecordBuf> = (0..5)
            .map(|_| {
                let cigar_ops = crate::io::bam::parse_cigar("8M").unwrap();
                let mut data = Data::default();
                data.insert(
                    Tag::EDIT_DISTANCE,
                    noodles_sam::alignment::record_buf::data::field::Value::from(0i32),
                );
                RecordBuf::builder()
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(NoodlesPosition::new(101).unwrap())
                    .set_mapping_quality(
                        noodles_sam::alignment::record::MappingQuality::new(60).unwrap(),
                    )
                    .set_cigar(cigar_ops.into_iter().collect::<Cigar>())
                    .set_sequence(Sequence::from(seq.as_ref()))
                    .set_quality_scores(QualityScores::from(vec![30u8; seq.len()]))
                    .set_data(data)
                    .build()
            })
            .collect();

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        write_bam(tmp_in.path(), &header, &records).unwrap();

        let variant = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 100,
                ref_base: b'A',
                alt_base: b'C',
            },
            expected_vaf: 1.0,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let config = EditConfig {
            input_bam: tmp_in.path().to_path_buf(),
            output_bam: tmp_out.path().to_path_buf(),
            variants: vec![variant],
            seed: 42,
            purity: 1.0,
            truth_vcf: None,
            sample_name: "TEST".to_string(),
        };

        let mut editor = BamEditor::new(config);
        editor.run().unwrap();

        let (_, out_records) = load_bam(tmp_out.path()).unwrap();
        // Any modified record should have NM=1.
        let modified_nm: Vec<i64> = out_records
            .iter()
            .filter_map(|r| {
                let s: Vec<u8> = r.sequence().as_ref().to_vec();
                if s.first() == Some(&b'C') {
                    r.data().get(&Tag::EDIT_DISTANCE).and_then(|v| v.as_int())
                } else {
                    None
                }
            })
            .collect();

        for nm in &modified_nm {
            assert_eq!(
                *nm, 1,
                "NM should be 1 after one SNV modification, got {nm}"
            );
        }
        assert!(!modified_nm.is_empty(), "some records should be modified");
    }

    // ------------------------------------------------------------------
    // Test 6: Strand balance
    // ------------------------------------------------------------------
    #[test]
    fn test_strand_balance() {
        let header = build_header(&[("chr1", 1000)]);
        let seq = b"ACGTACGTACGTACGT";
        let mut records = Vec::new();

        // 50 forward, 50 reverse reads at the same position.
        for _ in 0..50 {
            records.push(make_record(
                seq,
                &vec![30u8; seq.len()],
                0,
                100,
                Flags::empty(),
            ));
        }
        for _ in 0..50 {
            records.push(make_record(
                seq,
                &vec![30u8; seq.len()],
                0,
                100,
                Flags::REVERSE_COMPLEMENTED,
            ));
        }

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        write_bam(tmp_in.path(), &header, &records).unwrap();

        let variant = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 102,
                ref_base: b'G',
                alt_base: b'T',
            },
            expected_vaf: 0.5, // 50 of 100 reads
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let config = EditConfig {
            input_bam: tmp_in.path().to_path_buf(),
            output_bam: tmp_out.path().to_path_buf(),
            variants: vec![variant],
            seed: 42,
            purity: 1.0,
            truth_vcf: None,
            sample_name: "TEST".to_string(),
        };

        let mut editor = BamEditor::new(config);
        editor.run().unwrap();

        let (_, out_records) = load_bam(tmp_out.path()).unwrap();

        let mut fwd_alt = 0usize;
        let mut rev_alt = 0usize;
        for r in &out_records {
            let s: Vec<u8> = r.sequence().as_ref().to_vec();
            if s.get(2) == Some(&b'T') {
                if r.flags().is_reverse_complemented() {
                    rev_alt += 1;
                } else {
                    fwd_alt += 1;
                }
            }
        }

        // Both strands should have some alt reads.
        // With 50 of 100 reads modified at 50:50 strand ratio, we expect ~25 per strand.
        // Allow generous margin (at least 5 on each strand out of 50).
        assert!(
            fwd_alt >= 5,
            "forward strand should have some alt reads, got {fwd_alt}"
        );
        assert!(
            rev_alt >= 5,
            "reverse strand should have some alt reads, got {rev_alt}"
        );
    }

    // ------------------------------------------------------------------
    // Test 7: Quality preserved for SNV
    // ------------------------------------------------------------------
    #[test]
    fn test_quality_preserved() {
        let header = build_header(&[("chr1", 1000)]);
        let seq = b"ACGTACGT";
        let quals: Vec<u8> = (20u8..28).collect();
        let record = {
            let cigar_ops = crate::io::bam::parse_cigar("8M").unwrap();
            RecordBuf::builder()
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(NoodlesPosition::new(101).unwrap())
                .set_mapping_quality(
                    noodles_sam::alignment::record::MappingQuality::new(60).unwrap(),
                )
                .set_cigar(cigar_ops.into_iter().collect::<Cigar>())
                .set_sequence(Sequence::from(seq.as_ref()))
                .set_quality_scores(QualityScores::from(quals.clone()))
                .build()
        };

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        write_bam(tmp_in.path(), &header, &[record]).unwrap();

        let variant = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 102,
                ref_base: b'G',
                alt_base: b'T',
            },
            expected_vaf: 1.0,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let config = EditConfig {
            input_bam: tmp_in.path().to_path_buf(),
            output_bam: tmp_out.path().to_path_buf(),
            variants: vec![variant],
            seed: 42,
            purity: 1.0,
            truth_vcf: None,
            sample_name: "TEST".to_string(),
        };

        let mut editor = BamEditor::new(config);
        editor.run().unwrap();

        let (_, out_records) = load_bam(tmp_out.path()).unwrap();
        assert_eq!(out_records.len(), 1);
        let out_quals: Vec<u8> = out_records[0].quality_scores().as_ref().to_vec();
        assert_eq!(
            out_quals, quals,
            "quality scores should be identical after SNV"
        );
    }

    // ------------------------------------------------------------------
    // Test 8: Unmodified reads are byte-identical to input
    // ------------------------------------------------------------------
    #[test]
    fn test_unmodified_reads_unchanged() {
        let header = build_header(&[("chr1", 2000)]);
        let near_seq = b"ACGTACGT";
        let far_seq = b"TTTTTTTT";

        let near_record = make_record(near_seq, &[30u8; 8], 0, 100, Flags::empty());
        let far_record = make_record(far_seq, &[30u8; 8], 0, 500, Flags::empty()); // won't overlap pos 102

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        write_bam(tmp_in.path(), &header, &[near_record, far_record]).unwrap();

        let variant = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 102,
                ref_base: b'G',
                alt_base: b'T',
            },
            expected_vaf: 1.0,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let config = EditConfig {
            input_bam: tmp_in.path().to_path_buf(),
            output_bam: tmp_out.path().to_path_buf(),
            variants: vec![variant],
            seed: 42,
            purity: 1.0,
            truth_vcf: None,
            sample_name: "TEST".to_string(),
        };

        let mut editor = BamEditor::new(config);
        editor.run().unwrap();

        let (_, out_records) = load_bam(tmp_out.path()).unwrap();
        assert_eq!(out_records.len(), 2);

        // The far read (index 1) should be byte-for-byte identical.
        let out_far_seq: Vec<u8> = out_records[1].sequence().as_ref().to_vec();
        assert_eq!(
            out_far_seq,
            far_seq.to_vec(),
            "non-overlapping read must be unchanged"
        );
    }

    // ------------------------------------------------------------------
    // Test 9: Streaming - total record count preserved
    // ------------------------------------------------------------------
    #[test]
    fn test_streaming_record_count_preserved() {
        let header = build_header(&[("chr1", 5000)]);
        let seq = b"ACGTACGTACGTACGT"; // 16 bases

        // Records at varying positions to exercise the flush boundary logic.
        let mut records = Vec::new();
        for start in [0u64, 100, 500, 1000, 2000, 3000] {
            for _ in 0..5 {
                records.push(make_record(seq, &[30u8; 16], 0, start, Flags::empty()));
            }
        }
        let total = records.len();

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        write_bam(tmp_in.path(), &header, &records).unwrap();

        // Spike a variant in the middle.
        let variant = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 502,
                ref_base: b'G',
                alt_base: b'T',
            },
            expected_vaf: 0.5,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let config = EditConfig {
            input_bam: tmp_in.path().to_path_buf(),
            output_bam: tmp_out.path().to_path_buf(),
            variants: vec![variant],
            seed: 42,
            purity: 1.0,
            truth_vcf: None,
            sample_name: "TEST".to_string(),
        };

        let mut editor = BamEditor::new(config);
        editor.run().unwrap();

        let (_, out_records) = load_bam(tmp_out.path()).unwrap();
        assert_eq!(
            out_records.len(),
            total,
            "output must contain exactly as many records as input"
        );
    }

    // ------------------------------------------------------------------
    // Test 10: Multiple variants on the same chromosome
    // ------------------------------------------------------------------
    #[test]
    fn test_multiple_variants_same_chrom() {
        let header = build_header(&[("chr1", 5000)]);
        let seq = b"AAAAAAAAAAAAAAAAAAAA"; // 20 As

        // Two groups of reads at different positions.
        let mut records = Vec::new();
        for _ in 0..10 {
            records.push(make_record(seq, &[30u8; 20], 0, 100, Flags::empty()));
        }
        for _ in 0..10 {
            records.push(make_record(seq, &[30u8; 20], 0, 1000, Flags::empty()));
        }

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        write_bam(tmp_in.path(), &header, &records).unwrap();

        // Two variants, one at each group.
        let v1 = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 105,
                ref_base: b'A',
                alt_base: b'T',
            },
            expected_vaf: 1.0,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };
        let v2 = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 1005,
                ref_base: b'A',
                alt_base: b'C',
            },
            expected_vaf: 1.0,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let config = EditConfig {
            input_bam: tmp_in.path().to_path_buf(),
            output_bam: tmp_out.path().to_path_buf(),
            variants: vec![v1, v2],
            seed: 42,
            purity: 1.0,
            truth_vcf: None,
            sample_name: "TEST".to_string(),
        };

        let mut editor = BamEditor::new(config);
        let spiked = editor.run().unwrap();

        assert_eq!(spiked.len(), 2, "both variants should be spiked");
        assert!(spiked[0].alt_count > 0, "first variant should modify reads");
        assert!(
            spiked[1].alt_count > 0,
            "second variant should modify reads"
        );

        let (_, out_records) = load_bam(tmp_out.path()).unwrap();
        assert_eq!(out_records.len(), 20, "record count must be preserved");

        // First group: base at offset 5 should be T.
        let first_group_alt = out_records[..10]
            .iter()
            .filter(|r| r.sequence().as_ref().get(5) == Some(&b'T'))
            .count();
        assert!(
            first_group_alt > 0,
            "first variant should modify first group"
        );

        // Second group: base at offset 5 should be C.
        let second_group_alt = out_records[10..]
            .iter()
            .filter(|r| r.sequence().as_ref().get(5) == Some(&b'C'))
            .count();
        assert!(
            second_group_alt > 0,
            "second variant should modify second group"
        );
    }
}
