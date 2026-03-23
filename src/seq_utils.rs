//! Sequence utility functions shared across modules.

/// Return the DNA complement of a single base.
///
/// Handles both uppercase and lowercase IUPAC bases. Unknown bases map to `N`.
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        _ => b'N',
    }
}

/// Return the reverse complement of a DNA sequence.
///
/// Unknown bases map to `N`. The input is not modified.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement_basic() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
    }

    #[test]
    fn test_complement_lowercase() {
        assert_eq!(complement(b'a'), b't');
        assert_eq!(complement(b't'), b'a');
        assert_eq!(complement(b'c'), b'g');
        assert_eq!(complement(b'g'), b'c');
    }

    #[test]
    fn test_complement_unknown() {
        assert_eq!(complement(b'N'), b'N');
        assert_eq!(complement(b'X'), b'N');
    }

    #[test]
    fn test_reverse_complement_basic() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GCGC"), b"GCGC");
    }

    #[test]
    fn test_reverse_complement_asymmetric() {
        // ATCG → reverse is GCTA, complement is CGAT
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
    }

    #[test]
    fn test_reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn test_reverse_complement_roundtrip() {
        let seq = b"ACGTACGTNNNN";
        let rc = reverse_complement(seq);
        let rc_rc = reverse_complement(&rc);
        assert_eq!(&rc_rc, seq);
    }
}
