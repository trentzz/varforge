//! Reference genome loading via indexed FASTA.
//!
//! Provides random-access sequence extraction from a FASTA file with an associated `.fai` index.
//! The index must exist at `{path}.fai` before calling [`ReferenceGenome::open`].

use std::collections::HashMap;
use std::fmt;
use std::path::Path;
use std::sync::Mutex;

use anyhow::{anyhow, bail, Context, Result};
use noodles_core::Position;
use noodles_fasta::io::indexed_reader::Builder;
use noodles_fasta::io::IndexedReader;

use crate::core::types::Region;

/// An indexed FASTA reference genome.
///
/// Wraps an [`IndexedReader`] and a precomputed chromosome-length map loaded from the `.fai`
/// index.  Because [`IndexedReader`] needs `&mut self` for `query`, we protect it with a
/// [`Mutex`] so that `ReferenceGenome` can remain `Send + Sync` while still allowing immutable
/// `&self` calls from multiple threads.
pub struct ReferenceGenome {
    /// Path to the FASTA file; stored so the struct can be cloned by reopening.
    path: std::path::PathBuf,
    reader: Mutex<IndexedReader<noodles_fasta::io::BufReader<std::fs::File>>>,
    chrom_lengths: HashMap<String, u64>,
}

impl fmt::Debug for ReferenceGenome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ReferenceGenome")
            .field("chrom_lengths", &self.chrom_lengths)
            .finish_non_exhaustive()
    }
}

impl ReferenceGenome {
    /// Open an indexed FASTA file.
    ///
    /// The index is expected at `{path}.fai`.  Both the FASTA and the index must exist.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - the FASTA file does not exist,
    /// - the `.fai` file does not exist or is malformed, or
    /// - any I/O error occurs.
    pub fn open(path: &Path) -> Result<Self> {
        // Verify the FASTA file itself exists to give a clear error message.
        if !path.exists() {
            bail!("FASTA file not found: {}", path.display());
        }

        // The builder auto-loads `{path}.fai`.  We check explicitly so we get a
        // descriptive error rather than a cryptic I/O failure.
        let index_path = {
            let mut p = path.as_os_str().to_owned();
            p.push(".fai");
            std::path::PathBuf::from(p)
        };
        if !index_path.exists() {
            bail!(
                "FASTA index not found: {} (run `samtools faidx {}` to create it)",
                index_path.display(),
                path.display()
            );
        }

        let reader = Builder::default()
            .build_from_path(path)
            .with_context(|| format!("failed to open indexed FASTA: {}", path.display()))?;

        // Build chromosome-length map from the index records.
        let chrom_lengths: HashMap<String, u64> = reader
            .index()
            .as_ref()
            .iter()
            .map(|rec| {
                let name = String::from_utf8_lossy(rec.name()).into_owned();
                (name, rec.length())
            })
            .collect();

        Ok(Self {
            path: path.to_path_buf(),
            reader: Mutex::new(reader),
            chrom_lengths,
        })
    }

    /// Extract the uppercase nucleotide sequence for `region`.
    ///
    /// `region` uses **0-based, half-open** coordinates (`start` is inclusive, `end` is
    /// exclusive), matching BED-style conventions.  The returned `Vec<u8>` contains only ASCII
    /// uppercase bytes (`A`, `C`, `G`, `T`, `N`, …).
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - the chromosome does not exist in the index,
    /// - `region` extends beyond the chromosome length, or
    /// - an I/O error occurs.
    pub fn sequence(&self, region: &Region) -> Result<Vec<u8>> {
        let chrom_len = self
            .chrom_lengths
            .get(&region.chrom)
            .ok_or_else(|| anyhow!("chromosome '{}' not found in reference index", region.chrom))?;

        if region.end > *chrom_len {
            bail!(
                "region {}:{}-{} is out of bounds (chromosome length = {})",
                region.chrom,
                region.start,
                region.end,
                chrom_len
            );
        }

        if region.is_empty() {
            return Ok(Vec::new());
        }

        // Convert from 0-based half-open [start, end) to 1-based inclusive [start+1, end].
        // We add 1 to start before converting so that the conversion catches overflow.
        let start1: usize = usize::try_from(region.start + 1).map_err(|e| {
            anyhow!(
                "start position {} is too large for this platform: {}",
                region.start,
                e
            )
        })?;
        let end1: usize = usize::try_from(region.end).map_err(|e| {
            anyhow!(
                "end position {} is too large for this platform: {}",
                region.end,
                e
            )
        })?;

        let noodles_start = Position::try_from(start1)
            .map_err(|e| anyhow!("invalid start position {}: {}", region.start, e))?;
        let noodles_end = Position::try_from(end1)
            .map_err(|e| anyhow!("invalid end position {}: {}", region.end, e))?;

        let noodles_region =
            noodles_core::Region::new(region.chrom.as_bytes(), noodles_start..=noodles_end);

        let record = self
            .reader
            .lock()
            .map_err(|e| anyhow!("reference reader lock poisoned: {}", e))?
            .query(&noodles_region)
            .with_context(|| {
                format!(
                    "failed to query region {}:{}-{}",
                    region.chrom, region.start, region.end
                )
            })?;

        let seq: Vec<u8> = record
            .sequence()
            .as_ref()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();

        Ok(seq)
    }

    /// Returns a reference to the chromosome-name → length map built from the `.fai` index.
    pub fn chromosome_lengths(&self) -> &HashMap<String, u64> {
        &self.chrom_lengths
    }

    /// Returns the length of the named chromosome, or `None` if it is not in the index.
    pub fn chrom_len(&self, name: &str) -> Option<u64> {
        self.chrom_lengths.get(name).copied()
    }

    /// Returns `true` if the reference contains a chromosome named `name`.
    // Called only in tests; production code iterates known_chromosomes directly.
    #[allow(dead_code)]
    pub fn contains_chromosome(&self, name: &str) -> bool {
        self.chrom_lengths.contains_key(name)
    }

    /// Returns the total size of the genome (sum of all chromosome lengths).
    // Called only in tests; production code sums chrom_lengths inline.
    #[allow(dead_code)]
    pub fn genome_size(&self) -> u64 {
        self.chrom_lengths.values().sum()
    }
}

/// Clone by reopening the underlying FASTA file.
///
/// This is used by the parallel simulation pipeline so that each rayon worker
/// thread can own an independent file handle rather than contending on a
/// shared [`Mutex`].
///
/// # Panics
///
/// Panics if the FASTA file (or its `.fai` index) at the original path can no
/// longer be opened. The file must remain accessible for the lifetime of any
/// clone operation. In practice this is safe because the reference genome is
/// an immutable input that exists before the pipeline starts.
impl Clone for ReferenceGenome {
    fn clone(&self) -> Self {
        Self::open(&self.path).expect("failed to reopen reference genome for parallel worker")
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    /// Write a minimal FASTA + matching `.fai` into a temporary directory and
    /// return the path to the FASTA file.
    ///
    /// Test FASTA contents:
    /// ```
    /// >chr1
    /// ACGTACGTACGTACGT  (16 bases)
    /// >chr2
    /// NNNNacgtNNNN      (12 bases, mixed-case)
    /// ```
    fn write_test_fasta(dir: &TempDir) -> std::path::PathBuf {
        let fa_path = dir.path().join("test.fa");
        let fai_path = dir.path().join("test.fa.fai");

        // FASTA – each sequence on a single line for simplicity.
        let fa_content = b">chr1\nACGTACGTACGTACGT\n>chr2\nNNNNacgtNNNN\n";
        std::fs::write(&fa_path, fa_content).unwrap();

        // FAI – columns: name, length, offset, line_bases, line_width
        //   chr1: 16 bases, offset 6  (">chr1\n" = 6 bytes), line_bases=16, line_width=17 (\n)
        //   chr2: 12 bases, offset 29 (6 + 16 + 1 ["\n"] + 6 [">chr2\n"]), line_bases=12, line_width=13
        let fai_content = b"chr1\t16\t6\t16\t17\nchr2\t12\t29\t12\t13\n";
        std::fs::write(&fai_path, fai_content).unwrap();

        fa_path
    }

    // -----------------------------------------------------------------------
    // open() error cases
    // -----------------------------------------------------------------------

    #[test]
    fn test_open_missing_file() {
        let result = ReferenceGenome::open(Path::new("/nonexistent/path/genome.fa"));
        assert!(result.is_err(), "expected error for missing FASTA file");
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("not found"),
            "error message should mention 'not found': {msg}"
        );
    }

    #[test]
    fn test_open_missing_index() {
        let dir = TempDir::new().unwrap();
        let fa_path = dir.path().join("genome.fa");
        // Create the FASTA but NOT the .fai
        std::fs::write(&fa_path, b">chr1\nACGT\n").unwrap();

        let result = ReferenceGenome::open(&fa_path);
        assert!(result.is_err(), "expected error when .fai is absent");
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("index not found") || msg.contains("not found"),
            "error should mention missing index: {msg}"
        );
    }

    // -----------------------------------------------------------------------
    // Sequence extraction
    // -----------------------------------------------------------------------

    #[test]
    fn test_sequence_extraction() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let genome = ReferenceGenome::open(&fa).unwrap();

        // First 4 bases of chr1: "ACGT"
        let seq = genome.sequence(&Region::new("chr1", 0, 4)).unwrap();
        assert_eq!(seq, b"ACGT");

        // Middle 4 bases of chr1 starting at offset 4: "ACGT"
        let seq2 = genome.sequence(&Region::new("chr1", 4, 8)).unwrap();
        assert_eq!(seq2, b"ACGT");

        // All 16 bases of chr1
        let seq3 = genome.sequence(&Region::new("chr1", 0, 16)).unwrap();
        assert_eq!(seq3, b"ACGTACGTACGTACGT");
    }

    #[test]
    fn test_sequence_uppercase() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let genome = ReferenceGenome::open(&fa).unwrap();

        // chr2 is "NNNNacgtNNNN" – the middle 4 lowercase bases should be uppercased.
        let seq = genome.sequence(&Region::new("chr2", 4, 8)).unwrap();
        assert_eq!(seq, b"ACGT", "lowercase bases must be uppercased");
    }

    #[test]
    fn test_sequence_with_n_bases() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let genome = ReferenceGenome::open(&fa).unwrap();

        // First 4 bases of chr2 are 'N'
        let seq = genome.sequence(&Region::new("chr2", 0, 4)).unwrap();
        assert_eq!(seq, b"NNNN");
    }

    // -----------------------------------------------------------------------
    // Chromosome lengths
    // -----------------------------------------------------------------------

    #[test]
    fn test_chromosome_lengths() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let genome = ReferenceGenome::open(&fa).unwrap();

        let lengths = genome.chromosome_lengths();
        assert_eq!(lengths.get("chr1"), Some(&16u64));
        assert_eq!(lengths.get("chr2"), Some(&12u64));
    }

    // -----------------------------------------------------------------------
    // Out-of-bounds
    // -----------------------------------------------------------------------

    #[test]
    fn test_out_of_bounds() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let genome = ReferenceGenome::open(&fa).unwrap();

        // chr1 is 16 bases long; asking for [0, 17) should fail.
        let result = genome.sequence(&Region::new("chr1", 0, 17));
        assert!(result.is_err(), "expected out-of-bounds error");
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("out of bounds"),
            "error should mention out-of-bounds: {msg}"
        );
    }

    // -----------------------------------------------------------------------
    // contains_chromosome
    // -----------------------------------------------------------------------

    #[test]
    fn test_contains_chromosome() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let genome = ReferenceGenome::open(&fa).unwrap();

        assert!(genome.contains_chromosome("chr1"));
        assert!(genome.contains_chromosome("chr2"));
        assert!(!genome.contains_chromosome("chr3"));
        assert!(!genome.contains_chromosome("1")); // no "1"-style names in fixture
    }

    // -----------------------------------------------------------------------
    // genome_size
    // -----------------------------------------------------------------------

    #[test]
    fn test_genome_size() {
        let dir = TempDir::new().unwrap();
        let fa = write_test_fasta(&dir);
        let genome = ReferenceGenome::open(&fa).unwrap();

        // chr1 (16) + chr2 (12) = 28
        assert_eq!(genome.genome_size(), 28);
    }
}
