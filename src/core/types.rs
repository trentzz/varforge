/// A genomic region defined by chromosome, start, and end positions.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Region {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl Region {
    pub fn new(chrom: impl Into<String>, start: u64, end: u64) -> Self {
        Self {
            chrom: chrom.into(),
            start,
            end,
        }
    }

    #[must_use]
    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }
}

/// Records which variant a read pair carries, for FASTQ name annotation and sidecar output.
// Fields are populated for future use; not yet consumed by any output path.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VariantTag {
    pub chrom: String,
    pub pos: u64,
    pub vartype: String,
    pub vaf: f64,
    pub clone_id: Option<String>,
}

/// A simulated read pair.
#[derive(Debug, Clone)]
pub struct ReadPair {
    pub name: String,
    pub read1: Read,
    pub read2: Read,
    pub fragment_start: u64,
    pub fragment_length: usize,
    pub chrom: String,
    /// Variants carried by this read pair, populated during spike-in.
    pub variant_tags: Vec<VariantTag>,
    /// Reference sequence underlying read1, used to build the MD tag.
    ///
    /// Sliced from the unmodified reference before variant application, so
    /// the MD tag reflects true mismatches between the read and the genome.
    pub ref_seq_r1: Vec<u8>,
    /// Reference sequence underlying read2, used to build the MD tag.
    pub ref_seq_r2: Vec<u8>,
}

/// A single read with sequence and qualities.
#[derive(Debug, Clone)]
pub struct Read {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

impl Read {
    pub fn new(seq: Vec<u8>, qual: Vec<u8>) -> Self {
        assert_eq!(
            seq.len(),
            qual.len(),
            "sequence and quality length mismatch"
        );
        Self { seq, qual }
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    #[must_use]
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// Mutation types that can be spiked into reads.
#[derive(Debug, Clone, PartialEq)]
pub enum MutationType {
    Snv {
        pos: u64,
        ref_base: u8,
        alt_base: u8,
    },
    Indel {
        pos: u64,
        ref_seq: Vec<u8>,
        alt_seq: Vec<u8>,
    },
    Mnv {
        pos: u64,
        ref_seq: Vec<u8>,
        alt_seq: Vec<u8>,
    },
    /// Structural variant (>50 bp) with breakend notation for truth VCF.
    #[allow(dead_code)]
    Sv {
        sv_type: SvType,
        chrom: String,
        start: u64,
        end: u64,
    },
}

#[allow(dead_code)]
impl MutationType {
    /// Return the primary 0-based genomic position of this mutation.
    pub fn position(&self) -> u64 {
        match self {
            Self::Snv { pos, .. } => *pos,
            Self::Indel { pos, .. } => *pos,
            Self::Mnv { pos, .. } => *pos,
            Self::Sv { start, .. } => *start,
        }
    }

    /// Return a short string label for this mutation type.
    pub fn vartype_str(&self) -> &'static str {
        match self {
            Self::Snv { .. } => "SNV",
            Self::Indel { .. } => "INDEL",
            Self::Mnv { .. } => "MNV",
            Self::Sv { .. } => "SV",
        }
    }
}

/// Structural variant type classifications.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SvType {
    Deletion,
    Insertion,
    Inversion,
    Duplication,
    Translocation,
}

/// A variant with its genomic location and target VAF.
#[derive(Debug, Clone)]
pub struct Variant {
    pub chrom: String,
    pub mutation: MutationType,
    pub expected_vaf: f64,
    pub clone_id: Option<String>,
    /// Optional haplotype assignment (0 or 1).
    ///
    /// When set, the variant is only spiked into fragments sampled from the
    /// matching haplotype. When `None`, the variant is applied regardless.
    pub haplotype: Option<u8>,
    /// Cancer cell fraction of the clone carrying this variant.
    ///
    /// Populated when a `ClonalTree` is constructed from config.  `None`
    /// indicates no clonal-tree model was used; callers fall back to
    /// `expected_vaf` as a proxy.
    #[allow(dead_code)]
    pub ccf: Option<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_region_len() {
        let r = Region::new("chr1", 100, 200);
        assert_eq!(r.len(), 100);
        assert!(!r.is_empty());
    }

    #[test]
    fn test_region_empty() {
        let r = Region::new("chr1", 200, 200);
        assert!(r.is_empty());
    }

    #[test]
    fn test_read_construction() {
        let read = Read::new(vec![b'A', b'C', b'G'], vec![30, 30, 30]);
        assert_eq!(read.len(), 3);
    }

    #[test]
    #[should_panic(expected = "sequence and quality length mismatch")]
    fn test_read_length_mismatch() {
        Read::new(vec![b'A', b'C'], vec![30]);
    }
}
