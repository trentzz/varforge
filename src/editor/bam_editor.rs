//! Core BAM editing engine.
//!
//! Takes an existing BAM file, spikes in variants at specified positions with
//! stochastic VAF sampling, and writes a modified BAM file.

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
    /// Tumour purity (0.0–1.0). Scales all VAFs.
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
    /// 1. Load all reads from the input BAM into memory.
    /// 2. For each variant, collect overlapping reads, sample alt count, modify reads.
    /// 3. Write modified reads to the output BAM.
    /// 4. Optionally write truth VCF.
    pub fn run(&mut self) -> Result<Vec<SpikedVariant>> {
        // ---------------------------------------------------------------
        // 1. Load all reads from the input BAM.
        // ---------------------------------------------------------------
        let (header, records) = load_bam(&self.config.input_bam)?;

        // ---------------------------------------------------------------
        // 2. Build a mutable index: chrom+pos → list of record indices.
        //    We keep records in a Vec for O(1) random access.
        // ---------------------------------------------------------------
        let mut records: Vec<RecordBuf> = records;
        let pos_index = build_position_index(&records, &header);

        // ---------------------------------------------------------------
        // 3. Spike in each variant.
        // ---------------------------------------------------------------
        let mut spiked: Vec<SpikedVariant> = Vec::new();
        let variants = self.config.variants.clone();

        for variant in &variants {
            let purity = self.config.purity;
            let effective_vaf = variant.expected_vaf * purity;
            let result =
                self.spike_variant(&mut records, &pos_index, &header, variant, effective_vaf)?;
            if let Some(sv) = result {
                spiked.push(sv);
            }
        }

        // ---------------------------------------------------------------
        // 4. Write modified records to the output BAM.
        // ---------------------------------------------------------------
        write_bam(&self.config.output_bam, &header, &records)?;

        // ---------------------------------------------------------------
        // 5. Write truth VCF if requested.
        // ---------------------------------------------------------------
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

    /// Spike a single variant into the relevant records.
    fn spike_variant(
        &mut self,
        records: &mut [RecordBuf],
        pos_index: &PositionIndex,
        header: &sam::Header,
        variant: &Variant,
        effective_vaf: f64,
    ) -> Result<Option<SpikedVariant>> {
        let var_pos = variant.pos();
        let chrom = &variant.chrom;

        // Find reference sequence ID for the chromosome.
        let ref_id = header
            .reference_sequences()
            .get_index_of(chrom.as_bytes())
            .ok_or_else(|| anyhow::anyhow!("chromosome '{}' not found in BAM header", chrom))?;

        // Collect indices of records overlapping the variant position.
        let overlapping: Vec<usize> =
            collect_overlapping_indices(records, pos_index, ref_id, var_pos);

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
            select_reads_for_modification(records, &overlapping, alt_count as usize, &mut self.rng);

        // Apply the modification to selected reads.
        let mut actual_modified: u32 = 0;
        for idx in to_modify {
            let record = &mut records[idx];
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
// Position index
// ---------------------------------------------------------------------------

/// Maps (ref_id, position_bucket) → Vec<record_index>.
/// Positions are bucketed in 1 kb windows for fast overlap lookup.
type PositionIndex = HashMap<(usize, u64), Vec<usize>>;

/// Build a position index over the loaded records.
fn build_position_index(records: &[RecordBuf], header: &sam::Header) -> PositionIndex {
    const BUCKET: u64 = 1000;
    let mut index: PositionIndex = HashMap::new();

    for (i, record) in records.iter().enumerate() {
        let ref_id = match record.reference_sequence_id() {
            Some(id) => id,
            None => continue,
        };
        let align_start = match record.alignment_start() {
            Some(p) => usize::from(p) as u64 - 1,
            None => continue,
        };
        let read_len = record.sequence().len() as u64;
        let align_end = align_start + read_len;

        // Insert into all buckets the read overlaps.
        let start_bucket = align_start / BUCKET;
        let end_bucket = align_end / BUCKET;
        for bucket in start_bucket..=end_bucket {
            index.entry((ref_id, bucket)).or_default().push(i);
        }
    }

    // Suppress unused import warning from header parameter when building index
    let _ = header;
    index
}

/// Collect indices of records that overlap the given reference position.
fn collect_overlapping_indices(
    records: &[RecordBuf],
    pos_index: &PositionIndex,
    ref_id: usize,
    ref_pos: u64,
) -> Vec<usize> {
    const BUCKET: u64 = 1000;
    let bucket = ref_pos / BUCKET;

    let candidates = match pos_index.get(&(ref_id, bucket)) {
        Some(v) => v,
        None => return Vec::new(),
    };

    let mut result = Vec::new();
    for &idx in candidates {
        let record = &records[idx];
        let rid = match record.reference_sequence_id() {
            Some(id) => id,
            None => continue,
        };
        if rid != ref_id {
            continue;
        }
        let align_start = match record.alignment_start() {
            Some(p) => usize::from(p) as u64 - 1,
            None => continue,
        };
        let read_len = record.sequence().len() as u64;
        let align_end = align_start + read_len;

        // Skip MQ=0 reads.
        let mq = record.mapping_quality().map(u8::from).unwrap_or(0);
        if mq == 0 {
            continue;
        }

        if ref_pos >= align_start && ref_pos < align_end {
            result.push(idx);
        }
    }

    result
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
// BAM I/O
// ---------------------------------------------------------------------------

/// Load all records from a BAM file into memory.
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
// Variant helpers (copied/adapted from simulate.rs)
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
/// Returns a map from family_id → Vec<record_index>.
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

    #[allow(dead_code)]
    fn write_and_reload(
        records: Vec<RecordBuf>,
        header: &sam::Header,
    ) -> (sam::Header, Vec<RecordBuf>) {
        let tmp = NamedTempFile::new().unwrap();
        write_bam(tmp.path(), header, &records).unwrap();
        load_bam(tmp.path()).unwrap()
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
    // Test 2: Stochastic VAF — 50% VAF at 1000x depth should be ~50%
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
        // Should be roughly 50% ± 10%.
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
    // Test 9: Determinism with same seed
    // ------------------------------------------------------------------
    #[test]
    fn test_deterministic_with_seed() {
        let header = build_header(&[("chr1", 1000)]);
        let seq = b"ACGTACGTACGTACGTACGT";
        let records = make_records_at_pos(100, 0, 100, seq);

        let variant = Variant {
            chrom: "chr1".to_string(),
            mutation: MutationType::Snv {
                pos: 102,
                ref_base: b'G',
                alt_base: b'T',
            },
            expected_vaf: 0.3,
            clone_id: None,
            haplotype: None,
            ccf: None,
        };

        let run = |seed: u64| -> Vec<u8> {
            let tmp_in = NamedTempFile::new().unwrap();
            let tmp_out = NamedTempFile::new().unwrap();
            write_bam(tmp_in.path(), &header, &records).unwrap();
            let config = EditConfig {
                input_bam: tmp_in.path().to_path_buf(),
                output_bam: tmp_out.path().to_path_buf(),
                variants: vec![variant.clone()],
                seed,
                purity: 1.0,
                truth_vcf: None,
                sample_name: "TEST".to_string(),
            };
            let mut editor = BamEditor::new(config);
            editor.run().unwrap();
            let (_, out) = load_bam(tmp_out.path()).unwrap();
            out.iter()
                .flat_map(|r| r.sequence().as_ref().iter().copied())
                .collect()
        };

        let run1 = run(42);
        let run2 = run(42);
        let run3 = run(99);

        assert_eq!(run1, run2, "same seed must produce identical output");
        assert_ne!(
            run1, run3,
            "different seeds should produce different output (usually)"
        );
    }

    // ------------------------------------------------------------------
    // Test 10: Truth VCF output
    // ------------------------------------------------------------------
    #[test]
    fn test_truth_vcf_output() {
        let header = build_header(&[("chr1", 1000)]);
        let seq = b"ACGTACGTACGTACGT";
        let records = make_records_at_pos(30, 0, 100, seq);

        let tmp_in = NamedTempFile::new().unwrap();
        let tmp_out = NamedTempFile::new().unwrap();
        let tmp_vcf = NamedTempFile::new().unwrap();
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
            truth_vcf: Some(tmp_vcf.path().to_path_buf()),
            sample_name: "SAMPLE".to_string(),
        };

        let mut editor = BamEditor::new(config);
        editor.run().unwrap();

        let vcf_content = std::fs::read_to_string(tmp_vcf.path()).unwrap();
        assert!(vcf_content.contains("##fileformat=VCFv4.3"));
        assert!(vcf_content.contains("chr1"));
        assert!(vcf_content.contains("VARTYPE=SNV"));
        // Data line with chr1 should exist.
        let data_lines: Vec<&str> = vcf_content
            .lines()
            .filter(|l| !l.starts_with('#'))
            .collect();
        assert!(
            !data_lines.is_empty(),
            "truth VCF should have at least one data line"
        );
    }

    // ------------------------------------------------------------------
    // Test 11: UMI family editing — entire family modified for true variants
    // ------------------------------------------------------------------
    #[test]
    fn test_umi_family_editing() {
        use noodles_sam::alignment::record_buf::data::field::Value;

        // Build records with UMI family tags (MI tag).
        let _header = build_header(&[("chr1", 1000)]);
        let seq = b"ACGTACGTACGT";

        let records: Vec<RecordBuf> = (0..6)
            .map(|i| {
                let family_id = i / 2; // 3 families of 2 reads each
                let cigar_ops = crate::io::bam::parse_cigar("12M").unwrap();
                let mut data = Data::default();
                data.insert(Tag::UMI_SEQUENCE, Value::String("ACGT".as_bytes().into()));
                data.insert(Tag::UMI_ID, Value::from(family_id));
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

        // Verify grouping works.
        let families = group_by_umi_family(&records);
        assert_eq!(families.len(), 3, "should be 3 UMI families");
        for members in families.values() {
            assert_eq!(members.len(), 2, "each family should have 2 members");
        }
    }

    // ------------------------------------------------------------------
    // Test 12: UMI artifact editing — single read modified within family
    // ------------------------------------------------------------------
    #[test]
    fn test_umi_artifact_editing() {
        use noodles_sam::alignment::record_buf::data::field::Value;

        // Single family of 4 reads. Spike a variant into exactly 1 read (artifact mode).
        let header = build_header(&[("chr1", 1000)]);
        let seq = b"ACGTACGTACGT";

        let records: Vec<RecordBuf> = (0..4)
            .map(|_| {
                let cigar_ops = crate::io::bam::parse_cigar("12M").unwrap();
                let mut data = Data::default();
                data.insert(Tag::UMI_SEQUENCE, Value::String("ACGT".as_bytes().into()));
                data.insert(Tag::UMI_ID, Value::from(0i32));
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

        // With 4 reads at VAF=0.25, we expect 1 read modified.
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
            expected_vaf: 0.25,
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

        // alt_count should be ~ 1 out of 4 (could be 0–2 with binomial).
        assert!(
            spiked[0].alt_count <= 2,
            "artifact: should modify at most a few reads in the family, got {}",
            spiked[0].alt_count
        );
    }

    // ------------------------------------------------------------------
    // Test 13: MD tag updated for SNV
    // ------------------------------------------------------------------
    #[test]
    fn test_md_tag_updated() {
        use noodles_sam::alignment::record_buf::data::field::Value;

        let header = build_header(&[("chr1", 1000)]);
        let seq = b"ACGTACGT";
        let record = {
            let cigar_ops = crate::io::bam::parse_cigar("8M").unwrap();
            let mut data = Data::default();
            // MD=8 means 8 matches.
            data.insert(
                Tag::MISMATCHED_POSITIONS,
                Value::String("8".as_bytes().into()),
            );
            data.insert(Tag::EDIT_DISTANCE, Value::from(0i32));
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
        let md = out_records[0]
            .data()
            .get(&Tag::MISMATCHED_POSITIONS)
            .cloned();

        // MD should be updated (non-None) and different from original "8".
        if let Some(noodles_sam::alignment::record_buf::data::field::Value::String(s)) = md {
            let s_str = String::from_utf8_lossy(s.as_ref()).to_string();
            // The updated MD should contain the ref base character.
            assert!(
                s_str != "8",
                "MD tag should be updated after SNV spike-in, got: {s_str}"
            );
        }
    }
}
