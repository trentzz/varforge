use std::fs::File;
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        record::{
            cigar::{op::Kind, Op},
            data::field::Tag,
            Flags, MappingQuality,
        },
        record_buf::{data::field::Value, Cigar, Data, QualityScores, Sequence},
        RecordBuf,
    },
    header::record::value::{
        map::{self, header::Version, ReadGroup, ReferenceSequence},
        Map,
    },
};

use crate::io::config::SampleConfig;

/// BAM writer that outputs aligned paired-end reads with proper SAM headers.
pub struct BamWriter {
    writer: bam::io::Writer<bgzf::Writer<File>>,
    header: sam::Header,
    sample_name: String,
    /// Mapping quality written to every record.
    mapq: u8,
}

impl BamWriter {
    /// Create a new BAM writer at `path`.
    ///
    /// The BAM header is constructed from:
    /// - `ref_sequences`: list of (name, length) pairs for @SQ lines
    /// - `sample_config`: used for @RG read group fields
    pub fn new(
        path: &Path,
        ref_sequences: &[(String, u64)],
        sample_config: &SampleConfig,
        mapq: u8,
    ) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("failed to create BAM file: {}", path.display()))?;

        let platform = sample_config
            .platform
            .as_deref()
            .unwrap_or("ILLUMINA")
            .to_string();

        let sample_name = sample_config.name.clone();

        // Build the @HD line with VN:1.6 SO:coordinate
        let hd = {
            let mut hd = Map::<map::Header>::new(Version::new(1, 6));
            // Set sort order to coordinate via other_fields
            let so_tag =
                noodles_sam::header::record::value::map::tag::Other::try_from([b'S', b'O'])
                    .map_err(|_| anyhow::anyhow!("invalid SO tag"))?;
            hd.other_fields_mut().insert(so_tag, "coordinate".into());
            hd
        };

        // Build header with @SQ lines
        let mut builder = sam::Header::builder().set_header(hd);

        for (name, length) in ref_sequences {
            let len = NonZeroUsize::try_from(*length as usize)
                .with_context(|| format!("reference sequence length 0 for {}", name))?;
            builder =
                builder.add_reference_sequence(name.as_str(), Map::<ReferenceSequence>::new(len));
        }

        // Build @RG with SM, PL, LB fields via other_fields
        let mut rg = Map::<ReadGroup>::default();
        {
            let sm_tag =
                noodles_sam::header::record::value::map::tag::Other::try_from([b'S', b'M'])
                    .map_err(|_| anyhow::anyhow!("invalid SM tag"))?;
            rg.other_fields_mut()
                .insert(sm_tag, sample_name.as_str().into());

            let pl_tag =
                noodles_sam::header::record::value::map::tag::Other::try_from([b'P', b'L'])
                    .map_err(|_| anyhow::anyhow!("invalid PL tag"))?;
            rg.other_fields_mut()
                .insert(pl_tag, platform.as_str().into());

            let lb_tag =
                noodles_sam::header::record::value::map::tag::Other::try_from([b'L', b'B'])
                    .map_err(|_| anyhow::anyhow!("invalid LB tag"))?;
            rg.other_fields_mut()
                .insert(lb_tag, sample_name.as_str().into());
        }

        let header = builder.add_read_group(sample_name.as_str(), rg).build();

        let mut writer = bam::io::Writer::new(file);
        writer
            .write_header(&header)
            .context("failed to write BAM header")?;

        Ok(Self {
            writer,
            header,
            sample_name,
            mapq,
        })
    }

    /// Write an aligned read pair to the BAM file.
    ///
    /// - `pair`: the [`ReadPair`] holding sequences and qualities
    /// - `ref_id`: 0-based reference sequence index (index into @SQ lines)
    /// - `pos`: 0-based alignment start position on the reference (converted to 1-based internally)
    /// - `cigar_r1`: CIGAR string for read 1 (e.g. `"150M"`)
    /// - `cigar_r2`: CIGAR string for read 2 (e.g. `"150M"`)
    pub fn write_pair(
        &mut self,
        pair: &crate::core::types::ReadPair,
        ref_id: usize,
        pos: u64,
        cigar_r1: &str,
        cigar_r2: &str,
    ) -> Result<()> {
        // NM=0: simulated reads with no injected variants have zero mismatches.
        let data_r1: Data = [
            (Tag::READ_GROUP, Value::from(self.sample_name.as_str())),
            (Tag::EDIT_DISTANCE, Value::from(0i32)),
        ]
        .into_iter()
        .collect();
        let data_r2: Data = [
            (Tag::READ_GROUP, Value::from(self.sample_name.as_str())),
            (Tag::EDIT_DISTANCE, Value::from(0i32)),
        ]
        .into_iter()
        .collect();
        self.write_pair_inner(pair, ref_id, pos, cigar_r1, cigar_r2, data_r1, data_r2)
    }

    /// Write a single aligned read to the BAM file.
    ///
    /// Used for long-read platforms where each molecule produces one record.
    /// The record uses the pair QNAME and read1 sequence/quality. No mate
    /// information is set; flags indicate a primary, unpaired, mapped read.
    ///
    /// - `pair`: the [`ReadPair`] whose `read1` field is written
    /// - `ref_id`: 0-based reference sequence index
    /// - `pos`: 0-based alignment start position on the reference
    /// - `cigar`: CIGAR string for the read (e.g. `"15000M"`)
    pub fn write_single_read(
        &mut self,
        pair: &crate::core::types::ReadPair,
        ref_id: usize,
        pos: u64,
        cigar: &str,
    ) -> Result<()> {
        let pos1 = Position::new(pos as usize + 1)
            .ok_or_else(|| anyhow::anyhow!("invalid alignment position"))?;

        let mapq = MappingQuality::new(60).expect("60 is a valid mapping quality");

        let data: Data = [(Tag::READ_GROUP, Value::from(self.sample_name.as_str()))]
            .into_iter()
            .collect();

        // Flag 0: primary, mapped, single-segment (not paired).
        let flags = Flags::empty();

        let cigar_ops =
            parse_cigar(cigar).with_context(|| format!("failed to parse CIGAR: {cigar}"))?;

        let record = RecordBuf::builder()
            .set_name(pair.name.as_str())
            .set_flags(flags)
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(pos1)
            .set_mapping_quality(mapq)
            .set_cigar(cigar_ops.into_iter().collect::<Cigar>())
            .set_sequence(Sequence::from(pair.read1.seq.as_slice()))
            .set_quality_scores(QualityScores::from(pair.read1.qual.clone()))
            .set_data(data)
            .build();

        self.emit(&record)
    }

    /// Write a read pair with optional UMI tags (RX and MI).
    ///
    /// - `umi`: UMI barcode sequence string (written as `RX:Z:...`)
    /// - `family_id`: molecular family identifier (written as `MI:i:...`)
    #[allow(dead_code)]
    #[allow(clippy::too_many_arguments)]
    pub fn write_pair_with_umi(
        &mut self,
        pair: &crate::core::types::ReadPair,
        ref_id: usize,
        pos: u64,
        cigar_r1: &str,
        cigar_r2: &str,
        umi: &str,
        family_id: i32,
    ) -> Result<()> {
        let data_r1: Data = [
            (Tag::READ_GROUP, Value::from(self.sample_name.as_str())),
            (Tag::UMI_SEQUENCE, Value::from(umi)),
            (Tag::UMI_ID, Value::from(family_id)),
            (Tag::EDIT_DISTANCE, Value::from(0i32)),
        ]
        .into_iter()
        .collect();
        let data_r2: Data = [
            (Tag::READ_GROUP, Value::from(self.sample_name.as_str())),
            (Tag::UMI_SEQUENCE, Value::from(umi)),
            (Tag::UMI_ID, Value::from(family_id)),
            (Tag::EDIT_DISTANCE, Value::from(0i32)),
        ]
        .into_iter()
        .collect();
        self.write_pair_inner(pair, ref_id, pos, cigar_r1, cigar_r2, data_r1, data_r2)
    }

    /// Build and emit a paired-end read pair with caller-supplied per-read data fields.
    ///
    /// This is the single place that sets alignment positions, flags, CIGAR, sequence,
    /// and quality scores for a pair. Both `write_pair` and `write_pair_with_umi`
    /// delegate here, passing only the data fields that differ between them.
    #[allow(clippy::too_many_arguments)]
    fn write_pair_inner(
        &mut self,
        pair: &crate::core::types::ReadPair,
        ref_id: usize,
        pos: u64,
        cigar_r1: &str,
        cigar_r2: &str,
        data_r1: Data,
        data_r2: Data,
    ) -> Result<()> {
        // Positions in noodles are 1-based.
        let pos1 = Position::new(pos as usize + 1)
            .ok_or_else(|| anyhow::anyhow!("invalid alignment position"))?;

        // R2 starts at fragment_start + fragment_length - r2_read_length (0-based).
        // This is the correct paired-end position formula: R2 ends at the far end of
        // the fragment and reads back toward R1.
        let r2_len = pair.read2.seq.len();
        let pos2_zero = (pos as usize + pair.fragment_length).saturating_sub(r2_len);
        let pos2 = Position::new(pos2_zero + 1)
            .ok_or_else(|| anyhow::anyhow!("invalid mate alignment position"))?;

        let mapq = MappingQuality::new(self.mapq).expect("mapq is a valid mapping quality");
        // Use a saturating cast; fragment_length > i32::MAX is possible for long reads.
        let template_len = i32::try_from(pair.fragment_length).unwrap_or(i32::MAX);

        // Flags for R1: paired, properly segmented, mate reverse, first segment.
        let flags_r1 = Flags::SEGMENTED
            | Flags::PROPERLY_SEGMENTED
            | Flags::MATE_REVERSE_COMPLEMENTED
            | Flags::FIRST_SEGMENT;

        // Flags for R2: paired, properly segmented, reverse complemented, last segment.
        let flags_r2 = Flags::SEGMENTED
            | Flags::PROPERLY_SEGMENTED
            | Flags::REVERSE_COMPLEMENTED
            | Flags::LAST_SEGMENT;

        let cigar_ops_r1 =
            parse_cigar(cigar_r1).with_context(|| format!("failed to parse CIGAR: {cigar_r1}"))?;
        let cigar_ops_r2 =
            parse_cigar(cigar_r2).with_context(|| format!("failed to parse CIGAR: {cigar_r2}"))?;

        // R2 sequence must be reverse-complemented: BAM stores the sequence as it
        // appears on the forward strand, but R2 aligns to the reverse strand.
        let r2_seq_rc = crate::seq_utils::reverse_complement(&pair.read2.seq);
        let mut r2_qual_rev = pair.read2.qual.clone();
        r2_qual_rev.reverse();

        // Compute MD strings. R1 compares forward read against forward reference.
        // R2 is stored as RC, so the reference must also be RC for comparison.
        let md_r1 = build_md_string(&pair.read1.seq, &pair.ref_seq_r1);
        let ref_r2_rc = crate::seq_utils::reverse_complement(&pair.ref_seq_r2);
        let md_r2 = build_md_string(&r2_seq_rc, &ref_r2_rc);

        // Insert MD tags into the auxiliary data for each read.
        let mut data_r1 = data_r1;
        data_r1.insert(Tag::MISMATCHED_POSITIONS, Value::from(md_r1.as_str()));

        let mut data_r2 = data_r2;
        data_r2.insert(Tag::MISMATCHED_POSITIONS, Value::from(md_r2.as_str()));

        let r1 = build_record(
            &pair.name,
            flags_r1,
            ref_id,
            pos1,
            pos2,
            mapq,
            cigar_ops_r1,
            template_len,
            &pair.read1.seq,
            pair.read1.qual.clone(),
            data_r1,
        );

        let r2 = build_record(
            &pair.name,
            flags_r2,
            ref_id,
            pos2,
            pos1,
            mapq,
            cigar_ops_r2,
            -template_len,
            &r2_seq_rc,
            r2_qual_rev,
            data_r2,
        );

        self.emit(&r1)?;
        self.emit(&r2)?;

        Ok(())
    }

    /// Finalize and flush the BAM file.
    pub fn finish(mut self) -> Result<()> {
        self.writer
            .try_finish()
            .context("failed to finalize BAM file")?;
        Ok(())
    }

    /// Return a reference to the SAM header.
    #[allow(dead_code)]
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Write a single alignment record to the BAM stream.
    fn emit(&mut self, record: &RecordBuf) -> Result<()> {
        use noodles_sam::alignment::io::Write as _;
        // Split borrows: &mut self.writer and &self.header are disjoint.
        let writer = &mut self.writer;
        let header = &self.header;
        writer
            .write_alignment_record(header, record)
            .context("failed to write BAM record")
    }
}

/// Build the MD tag string for a single read.
///
/// The MD string encodes positions where the read differs from the reference.
/// Matching bases are counted and emitted as integers. Each mismatch emits the
/// current match count (even if zero) followed by the reference base at that
/// position, then resets the count.
///
/// If `ref_seq` is empty, the function returns the read length as a decimal
/// string, indicating all bases match (no reference available).
///
/// # Examples
///
/// - All 150 bp match: `"150"`
/// - Mismatch at position 50 (ref base A): `"50A99"`
/// - Mismatch at position 0 (ref base T): `"0T149"`
pub fn build_md_string(read_seq: &[u8], ref_seq: &[u8]) -> String {
    if ref_seq.is_empty() {
        return read_seq.len().to_string();
    }

    let mut result = String::new();
    let mut match_count: usize = 0;

    for (&read_base, &ref_base) in read_seq.iter().zip(ref_seq.iter()) {
        // Compare uppercase to handle any lowercase reference bases.
        if read_base.eq_ignore_ascii_case(&ref_base) {
            match_count += 1;
        } else {
            // Emit accumulated match count, then the mismatched reference base.
            result.push_str(&match_count.to_string());
            result.push(ref_base.to_ascii_uppercase() as char);
            match_count = 0;
        }
    }

    // Emit any trailing matches.
    result.push_str(&match_count.to_string());
    result
}

/// Parse a CIGAR string like `"150M"` or `"5M2I143M"` into a list of [`Op`]s.
pub fn parse_cigar(cigar_str: &str) -> Result<Vec<Op>> {
    let mut ops = Vec::new();
    let mut num_buf = String::new();

    for ch in cigar_str.chars() {
        if ch.is_ascii_digit() {
            num_buf.push(ch);
        } else {
            let len: usize = num_buf
                .parse()
                .with_context(|| format!("invalid CIGAR length before '{ch}'"))?;
            num_buf.clear();
            let kind = match ch {
                'M' => Kind::Match,
                'I' => Kind::Insertion,
                'D' => Kind::Deletion,
                'N' => Kind::Skip,
                'S' => Kind::SoftClip,
                'H' => Kind::HardClip,
                'P' => Kind::Pad,
                '=' => Kind::SequenceMatch,
                'X' => Kind::SequenceMismatch,
                other => anyhow::bail!("unknown CIGAR operation: '{other}'"),
            };
            ops.push(Op::new(kind, len));
        }
    }

    if !num_buf.is_empty() {
        anyhow::bail!("trailing digits in CIGAR string with no operation character");
    }

    Ok(ops)
}

/// Calculate the reference span consumed by a CIGAR string.
#[allow(dead_code)]
fn cigar_ref_span(cigar_str: &str) -> usize {
    parse_cigar(cigar_str)
        .unwrap_or_default()
        .iter()
        .filter_map(|op| op.kind().consumes_reference().then_some(op.len()))
        .sum()
}

/// Build a single SAM/BAM alignment record.
///
/// Centralises all record construction so that `write_pair_inner` does not
/// repeat the builder chain twice. The caller is responsible for choosing the
/// correct flags, positions, CIGAR ops, and data fields for each read in the
/// pair.
///
/// - `name`: read name (QNAME)
/// - `flags`: SAM flags for this read
/// - `ref_id`: 0-based reference sequence index
/// - `align_start`: 1-based alignment start position (noodles convention)
/// - `mate_start`: 1-based mate alignment start position
/// - `mapq`: mapping quality
/// - `cigar_ops`: parsed CIGAR operations
/// - `template_len`: signed template length (negative for R2)
/// - `seq`: read sequence as raw bytes
/// - `qual`: base quality scores
/// - `data`: auxiliary data fields (RG, RX, MI, …)
#[allow(clippy::too_many_arguments)]
fn build_record(
    name: &str,
    flags: Flags,
    ref_id: usize,
    align_start: Position,
    mate_start: Position,
    mapq: MappingQuality,
    cigar_ops: Vec<Op>,
    template_len: i32,
    seq: &[u8],
    qual: Vec<u8>,
    data: Data,
) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(align_start)
        .set_mapping_quality(mapq)
        .set_cigar(cigar_ops.into_iter().collect::<Cigar>())
        .set_mate_reference_sequence_id(ref_id)
        .set_mate_alignment_start(mate_start)
        .set_template_length(template_len)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(qual))
        .set_data(data)
        .build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::types::{Read, ReadPair};
    use crate::io::config::SampleConfig;
    use tempfile::NamedTempFile;

    fn make_sample_config(name: &str) -> SampleConfig {
        SampleConfig {
            name: name.to_string(),
            read_length: 150,
            coverage: 30.0,
            platform: Some("ILLUMINA".to_string()),
        }
    }

    fn make_read_pair(name: &str, len: usize) -> ReadPair {
        let seq = vec![b'A'; len];
        let qual = vec![30u8; len];
        ReadPair {
            name: name.to_string(),
            read1: Read::new(seq.clone(), qual.clone()),
            read2: Read::new(seq, qual),
            fragment_start: 100,
            fragment_length: len + 50,
            chrom: "chr1".to_string(),
            variant_tags: Vec::new(),
            ref_seq_r1: Vec::new(),
            ref_seq_r2: Vec::new(),
        }
    }

    fn ref_seqs() -> Vec<(String, u64)> {
        vec![
            ("chr1".to_string(), 248_956_422),
            ("chr2".to_string(), 242_193_529),
        ]
    }

    fn read_back_records(path: &Path, header: &sam::Header) -> Vec<RecordBuf> {
        let mut reader = bam::io::Reader::new(std::fs::File::open(path).unwrap());
        let _hdr = reader.read_header().unwrap();
        let mut records = Vec::new();
        for result in reader.record_bufs(header) {
            records.push(result.unwrap());
        }
        records
    }

    #[test]
    fn test_write_header() {
        let tmp = NamedTempFile::new().unwrap();
        let refs = ref_seqs();
        let cfg = make_sample_config("SAMPLE");

        let writer = BamWriter::new(tmp.path(), &refs, &cfg, 60).unwrap();
        let header = writer.header().clone();
        writer.finish().unwrap();

        // @HD present
        assert!(header.header().is_some());

        // @SQ lines
        assert_eq!(header.reference_sequences().len(), 2);
        assert!(header.reference_sequences().contains_key(&b"chr1"[..]));
        assert!(header.reference_sequences().contains_key(&b"chr2"[..]));

        // @RG line
        assert_eq!(header.read_groups().len(), 1);
        assert!(header.read_groups().contains_key(&b"SAMPLE"[..]));
    }

    #[test]
    fn test_write_single_pair() {
        let tmp = NamedTempFile::new().unwrap();
        let refs = ref_seqs();
        let cfg = make_sample_config("SAMPLE");

        let pair = make_read_pair("read1", 150);

        let mut writer = BamWriter::new(tmp.path(), &refs, &cfg, 60).unwrap();
        writer.write_pair(&pair, 0, 1000, "150M", "150M").unwrap();
        let header = writer.header().clone();
        writer.finish().unwrap();

        let records = read_back_records(tmp.path(), &header);
        assert_eq!(records.len(), 2, "expected two records (R1 and R2)");

        let r1 = &records[0];
        assert_eq!(r1.name().unwrap(), b"read1".as_ref());
        assert_eq!(r1.reference_sequence_id(), Some(0));
        // Position is 1-based in noodles: we wrote pos=1000 (0-based) => 1001
        assert_eq!(usize::from(r1.alignment_start().unwrap()), 1001);
        assert!(!r1.sequence().is_empty());
    }

    #[test]
    fn test_mate_flags() {
        let tmp = NamedTempFile::new().unwrap();
        let refs = ref_seqs();
        let cfg = make_sample_config("SAMPLE");

        let pair = make_read_pair("flagtest", 100);

        let mut writer = BamWriter::new(tmp.path(), &refs, &cfg, 60).unwrap();
        writer.write_pair(&pair, 0, 500, "100M", "100M").unwrap();
        let header = writer.header().clone();
        writer.finish().unwrap();

        let records = read_back_records(tmp.path(), &header);
        let r1 = &records[0];
        let r2 = &records[1];

        // R1 flags
        assert!(r1.flags().is_segmented(), "R1 should be paired");
        assert!(
            r1.flags().is_properly_segmented(),
            "R1 should be proper pair"
        );
        assert!(
            r1.flags().is_mate_reverse_complemented(),
            "R1 mate should be reverse"
        );
        assert!(r1.flags().is_first_segment(), "R1 should be first in pair");
        assert!(
            !r1.flags().is_last_segment(),
            "R1 should not be last in pair"
        );

        // R2 flags
        assert!(r2.flags().is_segmented(), "R2 should be paired");
        assert!(
            r2.flags().is_properly_segmented(),
            "R2 should be proper pair"
        );
        assert!(
            r2.flags().is_reverse_complemented(),
            "R2 should be reverse complemented"
        );
        assert!(r2.flags().is_last_segment(), "R2 should be last in pair");
        assert!(
            !r2.flags().is_first_segment(),
            "R2 should not be first in pair"
        );
    }

    #[test]
    fn test_umi_tags() {
        let tmp = NamedTempFile::new().unwrap();
        let refs = ref_seqs();
        let cfg = make_sample_config("SAMPLE");

        let pair = make_read_pair("umipair", 100);

        let mut writer = BamWriter::new(tmp.path(), &refs, &cfg, 60).unwrap();
        writer
            .write_pair_with_umi(&pair, 0, 500, "100M", "100M", "ACGTACGT", 42)
            .unwrap();
        let header = writer.header().clone();
        writer.finish().unwrap();

        let records = read_back_records(tmp.path(), &header);
        assert_eq!(records.len(), 2);

        for record in &records {
            let data = record.data();

            // RX tag (UMI sequence)
            let rx = data.get(&Tag::UMI_SEQUENCE).expect("RX tag missing");
            if let Value::String(s) = rx {
                let s_bytes: &[u8] = s.as_ref();
                assert_eq!(s_bytes, b"ACGTACGT");
            } else {
                panic!("RX tag should be a string");
            }

            // MI tag (family ID)
            let mi = data.get(&Tag::UMI_ID).expect("MI tag missing");
            assert_eq!(mi.as_int(), Some(42));
        }
    }

    #[test]
    fn test_cigar_preserved() {
        let tmp = NamedTempFile::new().unwrap();
        let refs = ref_seqs();
        let cfg = make_sample_config("SAMPLE");

        // Use a more complex CIGAR: 5M2I143M
        let pair = make_read_pair("cigartest", 150);

        let mut writer = BamWriter::new(tmp.path(), &refs, &cfg, 60).unwrap();
        writer
            .write_pair(&pair, 0, 200, "5M2I143M", "150M")
            .unwrap();
        let header = writer.header().clone();
        writer.finish().unwrap();

        let records = read_back_records(tmp.path(), &header);
        let r1 = &records[0];

        let ops: Vec<Op> = r1.cigar().as_ref().to_vec();

        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0].kind(), Kind::Match);
        assert_eq!(ops[0].len(), 5);
        assert_eq!(ops[1].kind(), Kind::Insertion);
        assert_eq!(ops[1].len(), 2);
        assert_eq!(ops[2].kind(), Kind::Match);
        assert_eq!(ops[2].len(), 143);
    }

    // --- MD tag unit tests ---

    #[test]
    fn test_md_all_match() {
        // 50 bp read where every base matches the reference.
        let seq = vec![b'A'; 50];
        let ref_seq = vec![b'A'; 50];
        assert_eq!(build_md_string(&seq, &ref_seq), "50");
    }

    #[test]
    fn test_md_single_snv() {
        // Read is all A's. Reference has T at position 5 (0-based).
        // MD should be: 5 matches, ref base T, 5 more matches → "5T5".
        let read = b"AAAAAAAAAAA";
        let mut ref_seq = b"AAAAAAAAAAA".to_vec();
        ref_seq[5] = b'T';
        assert_eq!(build_md_string(read, &ref_seq), "5T5");
    }

    #[test]
    fn test_md_empty_ref() {
        // Empty reference: return read length as all-match string.
        let read = vec![b'A'; 150];
        assert_eq!(build_md_string(&read, &[]), "150");
    }

    #[test]
    fn test_md_mismatch_at_start() {
        // Mismatch at position 0: "0" + ref_base + remaining_matches.
        let read = b"ACGT";
        let ref_seq = b"TCGT";
        assert_eq!(build_md_string(read, ref_seq), "0T3");
    }

    #[test]
    fn test_md_mismatch_at_end() {
        // Mismatch at last position.
        let read = b"ACGT";
        let ref_seq = b"ACGA";
        assert_eq!(build_md_string(read, ref_seq), "3A0");
    }

    #[test]
    fn test_md_multiple_mismatches() {
        // Mismatches at positions 1 and 3.
        // read:  A C G T
        // ref:   A T G A
        // MD:    1T1A0
        let read = b"ACGT";
        let ref_seq = b"ATGA";
        assert_eq!(build_md_string(read, ref_seq), "1T1A0");
    }

    #[test]
    fn test_read_group() {
        let tmp = NamedTempFile::new().unwrap();
        let refs = ref_seqs();
        let cfg = make_sample_config("MY_SAMPLE");

        let pair = make_read_pair("rgtest", 100);

        let mut writer = BamWriter::new(tmp.path(), &refs, &cfg, 60).unwrap();
        writer.write_pair(&pair, 0, 0, "100M", "100M").unwrap();
        let header = writer.header().clone();
        writer.finish().unwrap();

        let records = read_back_records(tmp.path(), &header);
        assert_eq!(records.len(), 2);

        for record in &records {
            let data = record.data();
            let rg = data.get(&Tag::READ_GROUP).expect("RG tag missing");
            if let Value::String(s) = rg {
                let s_bytes: &[u8] = s.as_ref();
                assert_eq!(s_bytes, b"MY_SAMPLE");
            } else {
                panic!("RG tag should be a string, got {:?}", rg);
            }
        }

        // Also verify RG is in the header
        assert!(header.read_groups().contains_key(&b"MY_SAMPLE"[..]));
    }
}
