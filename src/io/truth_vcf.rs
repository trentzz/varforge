use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::core::types::{MutationType, Variant};

/// Writer for truth VCF files produced by VarForge.
///
/// Outputs a VCF 4.3 file containing every variant spiked into the simulated
/// reads, annotated with expected VAF, clone assignment, variant type, and
/// cancer-cell fraction.  The file is intended as a ground-truth reference for
/// benchmarking variant callers.
pub struct TruthVcfWriter {
    writer: BufWriter<std::fs::File>,
}

impl TruthVcfWriter {
    /// Create a new truth VCF at `path`.
    ///
    /// Writes the complete VCF header immediately, including:
    /// - `##fileformat`, `##source`, INFO meta-lines, contig lines, FORMAT line.
    /// - A single-sample column named `sample_name`.
    ///
    /// `contigs` is a slice of `(name, length)` pairs that will be emitted as
    /// `##contig=<ID=…,length=…>` lines.
    pub fn new(path: &Path, sample_name: &str, contigs: &[(String, u64)]) -> Result<Self> {
        let file = std::fs::File::create(path)
            .with_context(|| format!("failed to create truth VCF: {}", path.display()))?;
        let mut writer = BufWriter::new(file);

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

        Ok(Self { writer })
    }

    /// Write a single variant record.
    ///
    /// - `ref_allele` / `alt_allele` are the REF and ALT byte sequences.
    /// - `n_alt_mol` is the number of read pairs (molecule families) carrying the alt allele.
    /// - `n_duplex_alt` is the number of duplex AB+BA family pairs carrying the alt allele.
    /// - The position written is 1-based (VCF convention); `variant.pos()` is
    ///   assumed to be 0-based and is incremented by 1.
    /// - QUAL is `.`, FILTER is `PASS`, FORMAT/SAMPLE is `GT\t0/1`.
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
        // CCF is not stored directly on Variant; use expected_vaf as a proxy when
        // no richer tumour model is present.  Callers that have a full ClonalTree
        // should populate a CCF field on Variant in a future iteration.
        let ccf = vaf;

        let info = format!(
            "EXPECTED_VAF={vaf:.6};CLONE={clone_id};VARTYPE={vartype};CCF={ccf:.6};\
             N_ALT_MOL={n_alt_mol};N_DUPLEX_ALT={n_duplex_alt}"
        );

        writeln!(
            self.writer,
            "{chrom}\t{pos}\t.\t{ref_str}\t{alt_str}\t.\tPASS\t{info}\tGT\t0/1",
            chrom = variant.chrom,
            pos = pos_1based,
        )?;

        Ok(())
    }

    /// Flush and close the writer.
    pub fn finish(mut self) -> Result<()> {
        self.writer.flush().context("failed to flush truth VCF")?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Helper methods on Variant / MutationType
// ---------------------------------------------------------------------------

impl Variant {
    /// 0-based genomic position of the variant.
    pub fn pos(&self) -> u64 {
        variant_pos(self)
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

/// Free-function equivalent of `Variant::pos()` for use in other modules.
pub fn variant_pos(v: &Variant) -> u64 {
    match &v.mutation {
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
        }
    }

    fn default_contigs() -> Vec<(String, u64)> {
        vec![
            ("chr1".to_string(), 248956422),
            ("chr2".to_string(), 242193529),
        ]
    }

    fn read_vcf(path: &std::path::Path) -> String {
        std::fs::read_to_string(path).expect("failed to read VCF file")
    }

    // ------------------------------------------------------------------
    // test_write_header
    // ------------------------------------------------------------------
    #[test]
    fn test_write_header() {
        let tmp = NamedTempFile::new().unwrap();
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
        let tmp = NamedTempFile::new().unwrap();
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
        let tmp = NamedTempFile::new().unwrap();
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
        let tmp = NamedTempFile::new().unwrap();
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
        let tmp = NamedTempFile::new().unwrap();
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
        let tmp = NamedTempFile::new().unwrap();
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
        let tmp = NamedTempFile::new().unwrap();
        let contigs = default_contigs();
        let mut writer = TruthVcfWriter::new(tmp.path(), "SAMPLE", &contigs).unwrap();

        // Write three variants; caller is responsible for sort order.
        let v1 = snv_variant("chr1", 100, 0.5, None);
        let v2 = snv_variant("chr1", 200, 0.3, None);
        let v3 = snv_variant("chr2", 50, 0.1, None);

        writer.write_variant(&v1, b"A", b"C", 10, 0).unwrap();
        writer.write_variant(&v2, b"G", b"T", 6, 0).unwrap();
        writer.write_variant(&v3, b"C", b"A", 1, 0).unwrap();
        writer.finish().unwrap();

        let vcf = read_vcf(tmp.path());
        let data_lines: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();

        assert_eq!(data_lines.len(), 3, "expected 3 data records");

        // Order preserved as written
        let pos0: u64 = data_lines[0].split('\t').nth(1).unwrap().parse().unwrap();
        let pos1: u64 = data_lines[1].split('\t').nth(1).unwrap().parse().unwrap();
        let pos2: u64 = data_lines[2].split('\t').nth(1).unwrap().parse().unwrap();
        assert_eq!(pos0, 101);
        assert_eq!(pos1, 201);
        assert_eq!(pos2, 51);
    }

    // ------------------------------------------------------------------
    // test_contig_headers
    // ------------------------------------------------------------------
    #[test]
    fn test_contig_headers() {
        let tmp = NamedTempFile::new().unwrap();
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
}
