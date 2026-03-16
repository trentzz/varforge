use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::{write::GzEncoder, Compression};

use crate::core::types::ReadPair;

/// Streaming gzip-compressed paired FASTQ writer.
///
/// Creates two output files (`{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`) in the
/// specified directory and writes [`ReadPair`] records to them.
pub struct FastqWriter {
    r1: GzEncoder<BufWriter<std::fs::File>>,
    r2: GzEncoder<BufWriter<std::fs::File>>,
}

impl FastqWriter {
    /// Create a new [`FastqWriter`] that writes to `{output_dir}/{sample_name}_R1.fastq.gz`
    /// and `{output_dir}/{sample_name}_R2.fastq.gz`.
    pub fn new(output_dir: &Path, sample_name: &str) -> Result<Self> {
        let r1_path = output_dir.join(format!("{}_R1.fastq.gz", sample_name));
        let r2_path = output_dir.join(format!("{}_R2.fastq.gz", sample_name));

        let r1_file = std::fs::File::create(&r1_path)
            .with_context(|| format!("failed to create R1 file: {}", r1_path.display()))?;
        let r2_file = std::fs::File::create(&r2_path)
            .with_context(|| format!("failed to create R2 file: {}", r2_path.display()))?;

        let r1 = GzEncoder::new(BufWriter::new(r1_file), Compression::default());
        let r2 = GzEncoder::new(BufWriter::new(r2_file), Compression::default());

        Ok(Self { r1, r2 })
    }

    /// Write one [`ReadPair`] with the given read name.
    ///
    /// Outputs records in FASTQ format:
    /// ```text
    /// @{read_name}/1
    /// {sequence}
    /// +
    /// {quality_as_phred33_ascii}
    /// ```
    pub fn write_pair(&mut self, pair: &ReadPair, read_name: &str) -> Result<()> {
        write_fastq_record(&mut self.r1, read_name, 1, &pair.read1.seq, &pair.read1.qual)?;
        write_fastq_record(&mut self.r2, read_name, 2, &pair.read2.seq, &pair.read2.qual)?;
        Ok(())
    }

    /// Write a batch of [`ReadPair`]s, naming each record `{name_prefix}:{index}` (1-based).
    #[allow(dead_code)]
    pub fn write_pairs(&mut self, pairs: &[ReadPair], name_prefix: &str) -> Result<()> {
        for (i, pair) in pairs.iter().enumerate() {
            let read_name = format!("{}:{}", name_prefix, i + 1);
            self.write_pair(pair, &read_name)?;
        }
        Ok(())
    }

    /// Flush and finalise both gzip streams.
    pub fn finish(self) -> Result<()> {
        self.r1
            .finish()
            .context("failed to finalise R1 gzip stream")?;
        self.r2
            .finish()
            .context("failed to finalise R2 gzip stream")?;
        Ok(())
    }
}

/// Write a single FASTQ record to `writer`.
///
/// Quality values are raw Phred scores (0–93); this function adds 33 to each byte to produce
/// the Phred+33 ASCII encoding required by the FASTQ specification.
fn write_fastq_record(
    writer: &mut dyn Write,
    read_name: &str,
    mate: u8,
    seq: &[u8],
    qual: &[u8],
) -> Result<()> {
    // Header line: @{read_name}/{mate}
    writer.write_all(b"@")?;
    writer.write_all(read_name.as_bytes())?;
    writer.write_all(b"/")?;
    writer.write_all(&[b'0' + mate])?;
    writer.write_all(b"\n")?;

    // Sequence
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;

    // Separator
    writer.write_all(b"+\n")?;

    // Quality — Phred+33 encoded (Fix H-7: write bytes directly, no allocation).
    for &q in qual {
        writer.write_all(&[q + 33])?;
    }
    writer.write_all(b"\n")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::types::{Read, ReadPair};
    use flate2::read::GzDecoder;
    use std::io::Read as IoRead;
    use tempfile::TempDir;

    fn make_read_pair(seq1: &[u8], qual1: &[u8], seq2: &[u8], qual2: &[u8]) -> ReadPair {
        ReadPair {
            name: "test".to_string(),
            read1: Read::new(seq1.to_vec(), qual1.to_vec()),
            read2: Read::new(seq2.to_vec(), qual2.to_vec()),
            fragment_start: 0,
            fragment_length: seq1.len(),
            chrom: "chr1".to_string(),
        }
    }

    fn decompress_file(path: &std::path::Path) -> String {
        let file = std::fs::File::open(path).expect("failed to open file");
        let mut gz = GzDecoder::new(file);
        let mut contents = String::new();
        gz.read_to_string(&mut contents).expect("failed to decompress");
        contents
    }

    #[test]
    fn test_write_single_pair() {
        let dir = TempDir::new().unwrap();
        let seq = b"ACGTACGT";
        let qual = vec![30u8; 8];

        let pair = make_read_pair(seq, &qual, seq, &qual);

        let mut writer = FastqWriter::new(dir.path(), "sample").unwrap();
        writer.write_pair(&pair, "sample:chr1:1").unwrap();
        writer.finish().unwrap();

        let r1 = decompress_file(&dir.path().join("sample_R1.fastq.gz"));
        let lines: Vec<&str> = r1.lines().collect();

        assert_eq!(lines[0], "@sample:chr1:1/1");
        assert_eq!(lines[1], "ACGTACGT");
        assert_eq!(lines[2], "+");
        // Q30 = 30 + 33 = 63 = '?'
        assert_eq!(lines[3], "????????");
    }

    #[test]
    fn test_write_multiple_pairs() {
        let dir = TempDir::new().unwrap();
        let seq = b"TTTT";
        let qual = vec![20u8; 4];

        let pairs: Vec<ReadPair> = (0..5)
            .map(|_| make_read_pair(seq, &qual, seq, &qual))
            .collect();

        let mut writer = FastqWriter::new(dir.path(), "multi").unwrap();
        writer.write_pairs(&pairs, "multi:chr1").unwrap();
        writer.finish().unwrap();

        let r1 = decompress_file(&dir.path().join("multi_R1.fastq.gz"));
        let r2 = decompress_file(&dir.path().join("multi_R2.fastq.gz"));

        // 5 pairs × 4 lines each = 20 lines per file
        assert_eq!(r1.lines().count(), 20);
        assert_eq!(r2.lines().count(), 20);

        // Check ordering: first record should be named multi:chr1:1/1
        let first_header = r1.lines().next().unwrap();
        assert_eq!(first_header, "@multi:chr1:1/1");
        let fifth_header = r1.lines().nth(16).unwrap(); // 4*(5-1) = 16
        assert_eq!(fifth_header, "@multi:chr1:5/1");
    }

    #[test]
    fn test_quality_encoding() {
        let dir = TempDir::new().unwrap();
        let seq = b"ACGT";
        // Raw Phred scores: Q30 and Q0
        let qual = vec![30u8, 0u8, 30u8, 0u8];

        let pair = make_read_pair(seq, &qual, seq, &qual);

        let mut writer = FastqWriter::new(dir.path(), "qual").unwrap();
        writer.write_pair(&pair, "qual:chr1:1").unwrap();
        writer.finish().unwrap();

        let r1 = decompress_file(&dir.path().join("qual_R1.fastq.gz"));
        let lines: Vec<&str> = r1.lines().collect();

        // Q30 + 33 = 63 = '?', Q0 + 33 = 33 = '!'
        assert_eq!(lines[3], "?!?!");
    }

    #[test]
    fn test_gzip_compressed() {
        let dir = TempDir::new().unwrap();
        let seq = b"AAAA";
        let qual = vec![20u8; 4];
        let pair = make_read_pair(seq, &qual, seq, &qual);

        let mut writer = FastqWriter::new(dir.path(), "gz").unwrap();
        writer.write_pair(&pair, "gz:chr1:1").unwrap();
        writer.finish().unwrap();

        // Read raw bytes and check gzip magic number (0x1f 0x8b)
        let r1_path = dir.path().join("gz_R1.fastq.gz");
        let raw = std::fs::read(&r1_path).unwrap();
        assert!(raw.len() >= 2, "file should not be empty");
        assert_eq!(raw[0], 0x1f, "expected gzip magic byte 0");
        assert_eq!(raw[1], 0x8b, "expected gzip magic byte 1");

        let r2_path = dir.path().join("gz_R2.fastq.gz");
        let raw2 = std::fs::read(&r2_path).unwrap();
        assert_eq!(raw2[0], 0x1f);
        assert_eq!(raw2[1], 0x8b);
    }

    #[test]
    fn test_paired_consistency() {
        let dir = TempDir::new().unwrap();
        let seq1 = b"ACGT";
        let seq2 = b"TGCA";
        let qual = vec![25u8; 4];

        let pairs: Vec<ReadPair> = (0..10)
            .map(|_| make_read_pair(seq1, &qual, seq2, &qual))
            .collect();

        let mut writer = FastqWriter::new(dir.path(), "paired").unwrap();
        writer.write_pairs(&pairs, "paired:chr1").unwrap();
        writer.finish().unwrap();

        let r1 = decompress_file(&dir.path().join("paired_R1.fastq.gz"));
        let r2 = decompress_file(&dir.path().join("paired_R2.fastq.gz"));

        assert_eq!(r1.lines().count(), r2.lines().count());
        // 10 pairs × 4 lines = 40
        assert_eq!(r1.lines().count(), 40);
    }

    #[test]
    fn test_empty_writes() {
        let dir = TempDir::new().unwrap();

        let mut writer = FastqWriter::new(dir.path(), "empty").unwrap();
        writer.write_pairs(&[], "empty").unwrap();
        writer.finish().unwrap();

        let r1 = decompress_file(&dir.path().join("empty_R1.fastq.gz"));
        let r2 = decompress_file(&dir.path().join("empty_R2.fastq.gz"));

        assert_eq!(r1, "");
        assert_eq!(r2, "");
    }

    #[test]
    fn test_output_paths() {
        let dir = TempDir::new().unwrap();
        let seq = b"GGGG";
        let qual = vec![30u8; 4];
        let pair = make_read_pair(seq, &qual, seq, &qual);

        let mut writer = FastqWriter::new(dir.path(), "mysample").unwrap();
        writer.write_pair(&pair, "mysample:chr1:1").unwrap();
        writer.finish().unwrap();

        assert!(dir.path().join("mysample_R1.fastq.gz").exists());
        assert!(dir.path().join("mysample_R2.fastq.gz").exists());
        // Make sure no other FASTQ files were created
        assert!(!dir.path().join("mysample_R1.fastq").exists());
        assert!(!dir.path().join("mysample_R2.fastq").exists());
    }
}
