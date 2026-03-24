use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result};
use noodles_bgzf as bgzf;

use crate::core::types::{MutationType, Variant};
use crate::variants::structural::StructuralVariant;
use crate::variants::structural::{sv_vcf_alt, sv_vcf_info};

/// A buffered record waiting to be written.
///
/// Records are sorted by (chrom_index, pos) before output so that the
/// resulting VCF is coordinate-sorted and can be tabix-indexed without a
/// separate sort step.
struct PendingRecord {
    /// Index into the contig list — used as the primary sort key.
    chrom_index: usize,
    /// 1-based genomic position — used as the secondary sort key.
    pos_1based: u64,
    /// Fully formatted VCF line (without trailing newline).
    line: String,
}

/// Writer for truth VCF files produced by VarForge.
///
/// Outputs a bgzip-compressed VCF 4.3 file containing every variant spiked
/// into the simulated reads, annotated with expected VAF, clone assignment,
/// variant type, and cancer-cell fraction.  The file is intended as a
/// ground-truth reference for benchmarking variant callers.
///
/// Records are buffered in memory and written in coordinate order when
/// `finish` is called.  This guarantees a sorted output that can be
/// tabix-indexed directly with `tabix -p vcf output.vcf.gz`.
pub struct TruthVcfWriter {
    /// Destination path — the bgzip writer is opened in `finish`.
    path: PathBuf,
    /// Sample name written to the column header.
    sample_name: String,
    /// Contig list written as `##contig` meta-lines.  Also used to determine
    /// chromosome sort order: the index of a chromosome in this list is its
    /// primary sort key.
    contigs: Vec<(String, u64)>,
    /// Map from chromosome name to its index in `contigs`.
    chrom_index: HashMap<String, usize>,
    /// Pending records, accumulated until `finish` is called.
    pending: Vec<PendingRecord>,
}

impl TruthVcfWriter {
    /// Create a new buffered truth VCF writer targeting `path`.
    ///
    /// The header is written when `finish` is called alongside the sorted
    /// records.  `contigs` is a slice of `(name, length)` pairs that will be
    /// emitted as `##contig=<ID=…,length=…>` lines and determine chromosome
    /// sort order.
    pub fn new(
        path: &std::path::Path,
        sample_name: &str,
        contigs: &[(String, u64)],
    ) -> Result<Self> {
        let chrom_index: HashMap<String, usize> = contigs
            .iter()
            .enumerate()
            .map(|(i, (name, _))| (name.clone(), i))
            .collect();

        Ok(Self {
            path: path.to_path_buf(),
            sample_name: sample_name.to_string(),
            contigs: contigs.to_vec(),
            chrom_index,
            pending: Vec::new(),
        })
    }

    /// Write the VCF header to an open bgzip writer.
    fn write_header(
        writer: &mut bgzf::Writer<std::fs::File>,
        sample_name: &str,
        contigs: &[(String, u64)],
    ) -> Result<()> {
        // File-format and source
        writeln!(writer, "##fileformat=VCFv4.3")?;
        writeln!(writer, "##source=VarForge")?;

        // INFO meta-lines
        writeln!(
            writer,
            r#"##INFO=<ID=EXPECTED_VAF,Number=1,Type=Float,Description="Expected variant allele frequency">"#
        )?;
        writeln!(
            writer,
            r#"##INFO=<ID=CLONE,Number=1,Type=String,Description="Clone assignment">"#
        )?;
        writeln!(
            writer,
            r#"##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Variant type: SNV, INDEL, MNV, SV, CNV">"#
        )?;
        writeln!(
            writer,
            r#"##INFO=<ID=CCF,Number=1,Type=Float,Description="Cancer cell fraction of assigned clone">"#
        )?;
        writeln!(
            writer,
            r#"##INFO=<ID=N_ALT_MOL,Number=1,Type=Integer,Description="Number of read pairs carrying the alt allele">"#
        )?;
        writeln!(
            writer,
            r#"##INFO=<ID=N_DUPLEX_ALT,Number=1,Type=Integer,Description="Number of duplex AB+BA family pairs carrying the alt allele">"#
        )?;
        writeln!(
            writer,
            r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type (DEL, INS, INV, DUP, BND)">"#
        )?;
        writeln!(
            writer,
            r#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">"#
        )?;
        writeln!(
            writer,
            r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">"#
        )?;

        // Contig lines
        for (name, length) in contigs {
            writeln!(writer, "##contig=<ID={name},length={length}>")?;
        }

        // FORMAT meta-line
        writeln!(
            writer,
            r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
        )?;

        // Column header
        writeln!(
            writer,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"
        )?;

        Ok(())
    }

    /// Select the VCF genotype string for a variant.
    ///
    /// Returns `"1/1"` when the variant is homozygous (VAF >= 0.99 or
    /// clone_id contains "hom"), and `"0/1"` otherwise.
    fn genotype(variant: &Variant) -> &'static str {
        let is_hom = variant.expected_vaf >= 0.99
            || variant
                .clone_id
                .as_deref()
                .map(|id| id.contains("hom"))
                .unwrap_or(false);
        if is_hom {
            "1/1"
        } else {
            "0/1"
        }
    }

    /// Buffer a single variant record for sorted output.
    ///
    /// - `ref_allele` / `alt_allele` are the REF and ALT byte sequences.
    /// - `n_alt_mol` is the number of read pairs (molecule families) carrying the alt allele.
    /// - `n_duplex_alt` is the number of duplex AB+BA family pairs carrying the alt allele.
    /// - The position written is 1-based (VCF convention); `variant.pos()` is
    ///   assumed to be 0-based and is incremented by 1.
    /// - QUAL is `.`, FILTER is `PASS`.
    /// - GT is `1/1` for homozygous variants (VAF >= 0.99 or clone_id contains "hom"),
    ///   `0/1` otherwise.
    ///
    /// Records are not written immediately.  Call `finish` to sort and flush.
    pub fn write_variant(
        &mut self,
        variant: &Variant,
        ref_allele: &[u8],
        alt_allele: &[u8],
        n_alt_mol: u32,
        n_duplex_alt: u32,
    ) -> Result<()> {
        let pos_1based = variant.pos() + 1;
        let ref_str = std::str::from_utf8(ref_allele).context("REF allele is not valid UTF-8")?;
        let alt_str = std::str::from_utf8(alt_allele).context("ALT allele is not valid UTF-8")?;

        let vartype = variant.vartype();
        let clone_id = variant.clone_id.as_deref().unwrap_or(".");
        let vaf = variant.expected_vaf;
        // Use the true CCF from the clonal tree when available. Fall back to
        // expected_vaf as a proxy when no clonal tree model is in use.
        let ccf = variant.ccf.unwrap_or(vaf);
        let gt = Self::genotype(variant);

        let info = format!(
            "EXPECTED_VAF={vaf:.6};CLONE={clone_id};VARTYPE={vartype};CCF={ccf:.6};\
             N_ALT_MOL={n_alt_mol};N_DUPLEX_ALT={n_duplex_alt}"
        );

        let line = format!(
            "{chrom}\t{pos}\t.\t{ref_str}\t{alt_str}\t.\tPASS\t{info}\tGT\t{gt}",
            chrom = variant.chrom,
            pos = pos_1based,
        );

        // Chromosomes not found in the contig list sort after all known contigs.
        let chrom_index = self
            .chrom_index
            .get(&variant.chrom)
            .copied()
            .unwrap_or(usize::MAX);

        self.pending.push(PendingRecord {
            chrom_index,
            pos_1based,
            line,
        });

        Ok(())
    }

    /// Buffer a structural variant record for sorted output.
    ///
    /// The ALT uses symbolic alleles (`<DEL>`, `<INS>`, `<INV>`, `<DUP>`) or
    /// BND notation for translocations.  The INFO field includes SVTYPE, END,
    /// and SVLEN in standard VCF format.
    ///
    /// - `expected_vaf` is the target allele frequency.
    /// - `clone_id` is the clone assignment (used to select GT and annotate CLONE).
    ///
    /// Records are not written immediately.  Call `finish` to sort and flush.
    pub fn write_sv(
        &mut self,
        sv: &StructuralVariant,
        expected_vaf: f64,
        clone_id: Option<&str>,
    ) -> Result<()> {
        let chrom = sv.chrom();
        let pos_1based = sv.start() + 1;
        let _end_1based = sv.end();

        let alt = sv_vcf_alt(sv);
        let sv_info = sv_vcf_info(sv);
        let clone_str = clone_id.unwrap_or(".");

        // Homozygous if VAF >= 0.99 or clone_id contains "hom".
        let is_hom = expected_vaf >= 0.99 || clone_id.map(|id| id.contains("hom")).unwrap_or(false);
        let gt = if is_hom { "1/1" } else { "0/1" };

        let info = format!(
            "{sv_info};EXPECTED_VAF={expected_vaf:.6};CLONE={clone_str};\
             VARTYPE=SV;CCF={expected_vaf:.6};N_ALT_MOL=0;N_DUPLEX_ALT=0"
        );

        let line = format!("{chrom}\t{pos_1based}\t.\tN\t{alt}\t.\tPASS\t{info}\tGT\t{gt}",);

        // Chromosomes not found in the contig list sort after all known contigs.
        let chrom_index = self.chrom_index.get(chrom).copied().unwrap_or(usize::MAX);

        self.pending.push(PendingRecord {
            chrom_index,
            pos_1based,
            line,
        });

        Ok(())
    }

    /// Sort all buffered records by (chromosome index, position) and write
    /// them to the bgzip-compressed output file, then close the writer.
    ///
    /// Chromosome order follows the order of the `contigs` slice passed to
    /// `new`.  Variants on chromosomes absent from that list sort last.
    /// Within a chromosome, records are sorted by ascending 1-based position.
    ///
    /// To build a tabix index after generation, run: `tabix -p vcf output.vcf.gz`
    pub fn finish(mut self) -> Result<()> {
        // Sort by chromosome index, then by position.
        self.pending.sort_by_key(|r| (r.chrom_index, r.pos_1based));

        let file = std::fs::File::create(&self.path)
            .with_context(|| format!("failed to create truth VCF: {}", self.path.display()))?;
        let mut writer = bgzf::Writer::new(file);

        Self::write_header(&mut writer, &self.sample_name, &self.contigs)?;

        for record in &self.pending {
            writeln!(writer, "{}", record.line)?;
        }

        writer
            .try_finish()
            .context("failed to finish bgzip truth VCF")?;

        build_tabix_index(&self.path)
            .with_context(|| format!("failed to write tabix index for {}", self.path.display()))?;

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tabix index generation
// ---------------------------------------------------------------------------

/// Build a `.tbi` tabix index for the bgzip VCF at `vcf_path`.
///
/// Runs `tabix -p vcf <vcf_path>` as a subprocess. The index is written to
/// `<vcf_path>.tbi` by tabix itself.
///
/// If `tabix` is not on `PATH`, logs a warning and returns `Ok(())`. The
/// caller is responsible for indexing manually in that case.
fn build_tabix_index(vcf_path: &std::path::Path) -> Result<()> {
    use std::process::Command;

    let path_str = vcf_path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("VCF path is not valid UTF-8"))?;

    let output = Command::new("tabix").args(["-p", "vcf", path_str]).output();

    match output {
        Ok(result) => {
            if result.status.success() {
                tracing::debug!("tabix index written for {}", vcf_path.display());
            } else {
                let stderr = String::from_utf8_lossy(&result.stderr);
                anyhow::bail!(
                    "tabix exited with status {}: {}",
                    result.status,
                    stderr.trim()
                );
            }
        }
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            tracing::warn!(
                "tabix not found on PATH; skipping index for {}. \
                 Run `tabix -p vcf {}` manually.",
                vcf_path.display(),
                vcf_path.display()
            );
        }
        Err(e) => {
            return Err(anyhow::Error::from(e)
                .context(format!("failed to run tabix for {}", vcf_path.display())));
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Helper methods on Variant / MutationType
// ---------------------------------------------------------------------------

impl Variant {
    /// 0-based genomic position of the variant.
    pub fn pos(&self) -> u64 {
        match &self.mutation {
            MutationType::Snv { pos, .. } => *pos,
            MutationType::Indel { pos, .. } => *pos,
            MutationType::Mnv { pos, .. } => *pos,
            MutationType::Sv { start, .. } => *start,
        }
    }

    /// Variant type string for the VARTYPE INFO field.
    pub fn vartype(&self) -> &'static str {
        match &self.mutation {
            MutationType::Snv { .. } => "SNV",
            MutationType::Indel { .. } => "INDEL",
            MutationType::Mnv { .. } => "MNV",
            MutationType::Sv { .. } => "SV",
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    fn snv_variant(chrom: &str, pos: u64, vaf: f64, clone: Option<&str>) -> Variant {
        Variant {
            chrom: chrom.to_string(),
            mutation: MutationType::Snv {
                pos,
                ref_base: b'A',
                alt_base: b'T',
            },
            expected_vaf: vaf,
            clone_id: clone.map(|s| s.to_string()),
            haplotype: None,
            ccf: None,
        }
    }

    fn insertion_variant(chrom: &str, pos: u64, vaf: f64) -> Variant {
        Variant {
            chrom: chrom.to_string(),
            mutation: MutationType::Indel {
                pos,
                ref_seq: b"A".to_vec(),
                alt_seq: b"ATG".to_vec(),
            },
            expected_vaf: vaf,
            clone_id: None,
            haplotype: None,
            ccf: None,
        }
    }

    fn deletion_variant(chrom: &str, pos: u64, vaf: f64) -> Variant {
        Variant {
            chrom: chrom.to_string(),
            mutation: MutationType::Indel {
                pos,
                ref_seq: b"ATG".to_vec(),
                alt_seq: b"A".to_vec(),
            },
            expected_vaf: vaf,
            clone_id: None,
            haplotype: None,
            ccf: None,
        }
    }

    fn mnv_variant(chrom: &str, pos: u64, vaf: f64) -> Variant {
        Variant {
            chrom: chrom.to_string(),
            mutation: MutationType::Mnv {
                pos,
                ref_seq: b"AC".to_vec(),
                alt_seq: b"TG".to_vec(),
            },
            expected_vaf: vaf,
            clone_id: None,
            haplotype: None,
            ccf: None,
        }
    }

    fn default_contigs() -> Vec<(String, u64)> {
        vec![
            ("chr1".to_string(), 248956422),
            ("chr2".to_string(), 242193529),
        ]
    }

    /// Read and decompress a bgzip VCF for assertions.
    fn read_vcf(path: &std::path::Path) -> String {
        use std::io::Read;
        let file = std::fs::File::open(path).expect("failed to open VCF");
        let mut reader = bgzf::Reader::new(file);
        let mut buf = String::new();
        reader.read_to_string(&mut buf).expect("failed to read VCF");
        buf
    }

    // ------------------------------------------------------------------
    // test_write_header
    // ------------------------------------------------------------------
    #[test]
    fn test_write_header() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let writer = TruthVcfWriter::new(tmp.path(), "SAMPLE01", &contigs).unwrap();
        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());

        assert!(vcf.contains("##fileformat=VCFv4.3"), "missing fileformat");
        assert!(vcf.contains("##source=VarForge"), "missing source");
        assert!(
            vcf.contains("##INFO=<ID=EXPECTED_VAF"),
            "missing EXPECTED_VAF INFO"
        );
        assert!(vcf.contains("##INFO=<ID=CLONE"), "missing CLONE INFO");
        assert!(vcf.contains("##INFO=<ID=VARTYPE"), "missing VARTYPE INFO");
        assert!(vcf.contains("##INFO=<ID=CCF"), "missing CCF INFO");
        assert!(
            vcf.contains("##INFO=<ID=N_ALT_MOL"),
            "missing N_ALT_MOL INFO"
        );
        assert!(
            vcf.contains("##INFO=<ID=N_DUPLEX_ALT"),
            "missing N_DUPLEX_ALT INFO"
        );
        assert!(vcf.contains("##FORMAT=<ID=GT"), "missing GT FORMAT");
        // Column header must contain sample name
        assert!(
            vcf.contains("SAMPLE01"),
            "missing sample name in column header"
        );
    }

    // ------------------------------------------------------------------
    // test_write_snv
    // ------------------------------------------------------------------
    #[test]
    fn test_write_snv() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();

        let v = snv_variant("chr1", 999, 0.35, Some("clone_A"));
        writer.write_variant(&v, b"A", b"T", 7, 0).unwrap();
        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());
        // Check non-header lines
        let record: Vec<&str> = vcf
            .lines()
            .find(|l| !l.starts_with('#'))
            .expect("no data line")
            .split('\t')
            .collect();

        assert_eq!(record[0], "chr1", "CHROM");
        assert_eq!(record[1], "1000", "POS (1-based)");
        assert_eq!(record[2], ".", "ID");
        assert_eq!(record[3], "A", "REF");
        assert_eq!(record[4], "T", "ALT");
        assert_eq!(record[5], ".", "QUAL");
        assert_eq!(record[6], "PASS", "FILTER");
        assert!(record[7].contains("VARTYPE=SNV"), "VARTYPE");
        assert_eq!(record[8], "GT", "FORMAT");
        assert_eq!(record[9], "0/1", "SAMPLE genotype");
    }

    // ------------------------------------------------------------------
    // test_write_indel
    // ------------------------------------------------------------------
    #[test]
    fn test_write_indel() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();

        // Insertion
        let ins = insertion_variant("chr1", 100, 0.2);
        writer.write_variant(&ins, b"A", b"ATG", 4, 0).unwrap();

        // Deletion
        let del = deletion_variant("chr1", 200, 0.15);
        writer.write_variant(&del, b"ATG", b"A", 3, 0).unwrap();

        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());
        let data_lines: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 2, "expected 2 INDEL records");

        // Insertion
        let ins_fields: Vec<&str> = data_lines[0].split('\t').collect();
        assert_eq!(ins_fields[3], "A");
        assert_eq!(ins_fields[4], "ATG");
        assert!(ins_fields[7].contains("VARTYPE=INDEL"));

        // Deletion
        let del_fields: Vec<&str> = data_lines[1].split('\t').collect();
        assert_eq!(del_fields[3], "ATG");
        assert_eq!(del_fields[4], "A");
        assert!(del_fields[7].contains("VARTYPE=INDEL"));
    }

    // ------------------------------------------------------------------
    // test_write_mnv
    // ------------------------------------------------------------------
    #[test]
    fn test_write_mnv() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();

        let v = mnv_variant("chr2", 500, 0.4);
        writer.write_variant(&v, b"AC", b"TG", 2, 0).unwrap();
        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());
        let record: Vec<&str> = vcf
            .lines()
            .find(|l| !l.starts_with('#'))
            .unwrap()
            .split('\t')
            .collect();

        assert_eq!(record[0], "chr2");
        assert_eq!(record[3], "AC");
        assert_eq!(record[4], "TG");
        assert!(record[7].contains("VARTYPE=MNV"));
    }

    // ------------------------------------------------------------------
    // test_expected_vaf_in_info
    // ------------------------------------------------------------------
    #[test]
    fn test_expected_vaf_in_info() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();

        let v = snv_variant("chr1", 0, 0.123456, None);
        writer.write_variant(&v, b"C", b"G", 1, 0).unwrap();
        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());
        let info = vcf
            .lines()
            .find(|l| !l.starts_with('#'))
            .unwrap()
            .split('\t')
            .nth(7)
            .unwrap()
            .to_string();

        assert!(
            info.contains("EXPECTED_VAF=0.123456"),
            "EXPECTED_VAF not found or incorrect in INFO: {info}"
        );
    }

    // ------------------------------------------------------------------
    // test_clone_assignment
    // ------------------------------------------------------------------
    #[test]
    fn test_clone_assignment() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();

        let v = snv_variant("chr1", 0, 0.3, Some("subclone_B"));
        writer.write_variant(&v, b"G", b"A", 5, 0).unwrap();
        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());
        let info = vcf
            .lines()
            .find(|l| !l.starts_with('#'))
            .unwrap()
            .split('\t')
            .nth(7)
            .unwrap()
            .to_string();

        assert!(
            info.contains("CLONE=subclone_B"),
            "CLONE not found or incorrect in INFO: {info}"
        );
    }

    // ------------------------------------------------------------------
    // test_multiple_variants
    // ------------------------------------------------------------------
    #[test]
    fn test_multiple_variants() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();

        // Write in non-sorted order: chr2 before chr1, chr1 out of position order.
        let v1 = snv_variant("chr1", 100, 0.5, None);
        let v2 = snv_variant("chr1", 200, 0.3, None);
        let v3 = snv_variant("chr2", 50, 0.1, None);

        writer.write_variant(&v3, b"C", b"A", 1, 0).unwrap();
        writer.write_variant(&v2, b"G", b"T", 6, 0).unwrap();
        writer.write_variant(&v1, b"A", b"C", 10, 0).unwrap();
        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());
        let data_lines: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();

        assert_eq!(data_lines.len(), 3, "expected 3 data records");

        // Output must be sorted: chr1:101, chr1:201, chr2:51
        let chrom0 = data_lines[0].split('\t').next().unwrap();
        let pos0: u64 = data_lines[0].split('\t').nth(1).unwrap().parse().unwrap();
        let chrom1 = data_lines[1].split('\t').next().unwrap();
        let pos1: u64 = data_lines[1].split('\t').nth(1).unwrap().parse().unwrap();
        let chrom2 = data_lines[2].split('\t').next().unwrap();
        let pos2: u64 = data_lines[2].split('\t').nth(1).unwrap().parse().unwrap();

        assert_eq!(chrom0, "chr1");
        assert_eq!(pos0, 101);
        assert_eq!(chrom1, "chr1");
        assert_eq!(pos1, 201);
        assert_eq!(chrom2, "chr2");
        assert_eq!(pos2, 51);
    }

    // ------------------------------------------------------------------
    // test_sorted_output (T114)
    // ------------------------------------------------------------------
    /// Verify that records submitted in arbitrary order are written in
    /// coordinate-sorted order, matching tabix requirements.
    #[test]
    fn test_sorted_output() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();

        // Submit in reverse chromosome and position order.
        let v_chr2 = snv_variant("chr2", 1000, 0.3, None);
        let v_chr1_late = snv_variant("chr1", 5000, 0.2, None);
        let v_chr1_early = snv_variant("chr1", 100, 0.4, None);

        writer.write_variant(&v_chr2, b"A", b"T", 3, 0).unwrap();
        writer
            .write_variant(&v_chr1_late, b"C", b"G", 5, 0)
            .unwrap();
        writer
            .write_variant(&v_chr1_early, b"G", b"A", 4, 0)
            .unwrap();
        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());
        let data_lines: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();

        assert_eq!(data_lines.len(), 3);

        let chrom_pos: Vec<(&str, u64)> = data_lines
            .iter()
            .map(|l| {
                let mut cols = l.split('\t');
                let chrom = cols.next().unwrap();
                let pos: u64 = cols.next().unwrap().parse().unwrap();
                (chrom, pos)
            })
            .collect();

        assert_eq!(chrom_pos[0], ("chr1", 101), "first record must be chr1:101");
        assert_eq!(
            chrom_pos[1],
            ("chr1", 5001),
            "second record must be chr1:5001"
        );
        assert_eq!(
            chrom_pos[2],
            ("chr2", 1001),
            "third record must be chr2:1001"
        );
    }

    // ------------------------------------------------------------------
    // test_contig_headers
    // ------------------------------------------------------------------
    #[test]
    fn test_contig_headers() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = vec![
            ("chrX".to_string(), 156040895u64),
            ("chrY".to_string(), 57227415u64),
        ];
        let writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();
        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());

        assert!(
            vcf.contains("##contig=<ID=chrX,length=156040895>"),
            "missing chrX contig line"
        );
        assert!(
            vcf.contains("##contig=<ID=chrY,length=57227415>"),
            "missing chrY contig line"
        );
    }

    // ------------------------------------------------------------------
    // test_homozygous_genotype (T079)
    // ------------------------------------------------------------------
    #[test]
    fn test_homozygous_genotype() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();

        // Germline homozygous by clone_id
        let hom_by_clone = snv_variant("chr1", 100, 0.5, Some("germline_hom"));
        writer
            .write_variant(&hom_by_clone, b"A", b"T", 10, 0)
            .unwrap();

        // Homozygous by VAF >= 0.99
        let hom_by_vaf = snv_variant("chr1", 200, 1.0, None);
        writer
            .write_variant(&hom_by_vaf, b"C", b"G", 10, 0)
            .unwrap();

        // Somatic heterozygous
        let het = snv_variant("chr1", 300, 0.35, Some("clone_A"));
        writer.write_variant(&het, b"G", b"A", 10, 0).unwrap();

        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());
        let data_lines: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 3);

        let gt0 = data_lines[0].split('\t').nth(9).unwrap();
        let gt1 = data_lines[1].split('\t').nth(9).unwrap();
        let gt2 = data_lines[2].split('\t').nth(9).unwrap();

        assert_eq!(gt0, "1/1", "germline_hom clone should produce GT=1/1");
        assert_eq!(gt1, "1/1", "VAF=1.0 should produce GT=1/1");
        assert_eq!(gt2, "0/1", "somatic variant should produce GT=0/1");
    }

    // ------------------------------------------------------------------
    // test_sv_record (T082)
    // ------------------------------------------------------------------
    #[test]
    fn test_sv_record() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();

        // Deletion SV: chr1:1000-2000
        let del_sv = StructuralVariant::Deletion {
            chrom: "chr1".to_string(),
            start: 1000,
            end: 2000,
        };
        writer.write_sv(&del_sv, 0.45, Some("clone_A")).unwrap();

        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());
        let record: Vec<&str> = vcf
            .lines()
            .find(|l| !l.starts_with('#'))
            .expect("no data line")
            .split('\t')
            .collect();

        assert_eq!(record[0], "chr1", "CHROM");
        assert_eq!(record[1], "1001", "POS (1-based, start+1)");
        assert_eq!(record[3], "N", "REF for SV is N");
        assert_eq!(record[4], "<DEL>", "ALT for deletion is <DEL>");
        assert_eq!(record[6], "PASS", "FILTER");

        let info = record[7];
        assert!(info.contains("SVTYPE=DEL"), "INFO must contain SVTYPE=DEL");
        assert!(info.contains("END=2000"), "INFO must contain END=2000");
        assert!(
            info.contains("SVLEN=-1000"),
            "INFO must contain SVLEN=-1000"
        );

        assert_eq!(record[8], "GT", "FORMAT");
        assert_eq!(record[9], "0/1", "GT for VAF=0.45 is 0/1");
    }

    // ------------------------------------------------------------------
    // test_bgzip_output (T081)
    // ------------------------------------------------------------------
    #[test]
    fn test_bgzip_output() {
        use std::io::Read;

        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let contigs = default_contigs();
        let writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();
        writer.finish().unwrap();

        // Read raw bytes — bgzip files start with the gzip magic number 0x1f 0x8b.
        let mut raw = Vec::new();
        std::fs::File::open(tmp.path())
            .unwrap()
            .read_to_end(&mut raw)
            .unwrap();
        assert!(
            raw.starts_with(&[0x1f, 0x8b]),
            "output must start with gzip magic bytes"
        );
    }

    #[test]
    fn test_tbi_index_created() {
        // After finish(), a .tbi tabix index must exist if tabix is on PATH.
        // Skip this test when tabix is not installed.
        use std::process::Command;
        use tempfile::TempDir;

        if Command::new("tabix").arg("--version").output().is_err() {
            eprintln!("skipping test_tbi_index_created: tabix not found on PATH");
            return;
        }

        let dir = TempDir::new().unwrap();
        let vcf_path = dir.path().join("test.vcf.gz");
        let tbi_path = dir.path().join("test.vcf.gz.tbi");

        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(&vcf_path, "SAMPLE", &contigs).unwrap();

        let v = snv_variant("chr1", 999, 0.35, Some("clone_A"));
        writer.write_variant(&v, b"A", b"T", 7, 0).unwrap();
        writer.finish().unwrap();

        assert!(
            tbi_path.exists(),
            ".tbi index was not created at {:?}",
            tbi_path
        );
        assert!(
            tbi_path.metadata().unwrap().len() > 0,
            ".tbi index file is empty"
        );
    }
}
