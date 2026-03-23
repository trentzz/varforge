//! COSMIC mutational signatures for SBS96 trinucleotide context weighting.
//!
//! Each signature is a probability vector of length 96 over the SBS96 channels.
//! Channels are ordered: 6 substitution types × 16 trinucleotide contexts.
//!
//! The 6 substitution types follow the COSMIC pyrimidine convention (central
//! base is always C or T):
//!   0 = C>A, 1 = C>G, 2 = C>T, 3 = T>A, 4 = T>C, 5 = T>G
//!
//! The 16 trinucleotide contexts are indexed by:
//!   ctx_idx = base_to_idx(5p_base) * 4 + base_to_idx(3p_base)
//! where A=0, C=1, G=2, T=3.

// Public API items are provided for use by downstream callers (e.g. output
// formatters, summary tools). Not all are currently wired into the binary.
#![allow(dead_code)]

use crate::seq_utils::complement;

/// The SBS96 substitution types in channel order.
pub const SUBSTITUTION_TYPES: [&str; 6] = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"];

/// The 16 trinucleotide contexts (5'-base, 3'-base) in channel order.
pub const CONTEXTS: [[u8; 2]; 16] = [
    [b'A', b'A'],
    [b'A', b'C'],
    [b'A', b'G'],
    [b'A', b'T'],
    [b'C', b'A'],
    [b'C', b'C'],
    [b'C', b'G'],
    [b'C', b'T'],
    [b'G', b'A'],
    [b'G', b'C'],
    [b'G', b'G'],
    [b'G', b'T'],
    [b'T', b'A'],
    [b'T', b'C'],
    [b'T', b'G'],
    [b'T', b'T'],
];

/// Return the SBS96 channel index for a substitution type and trinucleotide context.
///
/// `sub_type`: 0=C>A, 1=C>G, 2=C>T, 3=T>A, 4=T>C, 5=T>G
/// `context_5p`: the 5' flanking base (A/C/G/T)
/// `context_3p`: the 3' flanking base (A/C/G/T)
pub fn channel_index(sub_type: usize, context_5p: u8, context_3p: u8) -> usize {
    let ctx_idx = base_to_idx(context_5p) * 4 + base_to_idx(context_3p);
    sub_type * 16 + ctx_idx
}

fn base_to_idx(base: u8) -> usize {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

/// Parse the SBS96 channel index from a trinucleotide context and ref/alt bases.
///
/// Follows the COSMIC SBS96 pyrimidine convention: the central base is always
/// C or T. If the ref base is A or G, take the reverse complement of the
/// trinucleotide and swap the substitution accordingly.
pub fn sbs96_index(context_5p: u8, ref_base: u8, alt_base: u8, context_3p: u8) -> Option<usize> {
    let (sub_type, ctx_5p, ctx_3p) =
        normalise_to_pyrimidine(context_5p, ref_base, alt_base, context_3p)?;
    Some(channel_index(sub_type, ctx_5p, ctx_3p))
}

/// Normalise to the pyrimidine convention for SBS96 indexing.
///
/// Returns (sub_type_idx, context_5p, context_3p) where sub_type_idx is
/// 0=C>A, 1=C>G, 2=C>T, 3=T>A, 4=T>C, 5=T>G.
fn normalise_to_pyrimidine(
    ctx_5p: u8,
    ref_base: u8,
    alt_base: u8,
    ctx_3p: u8,
) -> Option<(usize, u8, u8)> {
    let r = ref_base.to_ascii_uppercase();
    let a = alt_base.to_ascii_uppercase();

    if matches!(r, b'C' | b'T') {
        // Already in pyrimidine convention.
        let sub_type = match (r, a) {
            (b'C', b'A') => 0,
            (b'C', b'G') => 1,
            (b'C', b'T') => 2,
            (b'T', b'A') => 3,
            (b'T', b'C') => 4,
            (b'T', b'G') => 5,
            _ => return None,
        };
        Some((
            sub_type,
            ctx_5p.to_ascii_uppercase(),
            ctx_3p.to_ascii_uppercase(),
        ))
    } else {
        // Take reverse complement.
        let r_rc = complement(r);
        let a_rc = complement(a);
        let sub_type = match (r_rc, a_rc) {
            (b'C', b'A') => 0,
            (b'C', b'G') => 1,
            (b'C', b'T') => 2,
            (b'T', b'A') => 3,
            (b'T', b'C') => 4,
            (b'T', b'G') => 5,
            _ => return None,
        };
        // Reverse complement swaps 5'/3' positions.
        Some((
            sub_type,
            complement(ctx_3p).to_ascii_uppercase(),
            complement(ctx_5p).to_ascii_uppercase(),
        ))
    }
}

/// Extract the trinucleotide context at `pos_in_region` within a pre-fetched region sequence.
///
/// Returns (5p_base, central_base, 3p_base) or None if at the sequence boundary.
///
/// `region_seq`: the sequence of the pre-fetched region (0-based positions).
/// `pos_in_region`: the position of the variant within `region_seq`.
pub fn trinucleotide_context(region_seq: &[u8], pos_in_region: usize) -> Option<(u8, u8, u8)> {
    if pos_in_region == 0 || pos_in_region + 1 >= region_seq.len() {
        return None;
    }
    Some((
        region_seq[pos_in_region - 1],
        region_seq[pos_in_region],
        region_seq[pos_in_region + 1],
    ))
}

/// Return the SBS96 channel index for a variant at `pos_in_region` in `region_seq`.
///
/// Returns None if at the sequence boundary or the substitution is not a recognised SNV type.
pub fn variant_sbs96_index(
    region_seq: &[u8],
    pos_in_region: usize,
    ref_base: u8,
    alt_base: u8,
) -> Option<usize> {
    let (ctx_5p, _central, ctx_3p) = trinucleotide_context(region_seq, pos_in_region)?;
    sbs96_index(ctx_5p, ref_base, alt_base, ctx_3p)
}

/// A built-in COSMIC SBS signature as a 96-channel probability vector.
///
/// All values are non-negative and sum to 1.0.
pub type SignatureVec = [f64; 96];

/// Retrieve the built-in probability vector for a named COSMIC SBS signature.
///
/// Returns `None` if the signature name is not recognised.
pub fn builtin_signature(name: &str) -> Option<&'static SignatureVec> {
    match name {
        "SBS1" => Some(&SBS1),
        "SBS2" => Some(&SBS2),
        "SBS3" => Some(&SBS3),
        "SBS4" => Some(&SBS4),
        "SBS5" => Some(&SBS5),
        "SBS7a" => Some(&SBS7A),
        "SBS13" => Some(&SBS13),
        "SBS17a" => Some(&SBS17A),
        _ => None,
    }
}

/// List all built-in signature names.
pub fn builtin_signature_names() -> &'static [&'static str] {
    &[
        "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS7a", "SBS13", "SBS17a",
    ]
}

// ---------------------------------------------------------------------------
// SBS96 signature probability vectors
//
// Channel layout: sub_type * 16 + ctx_idx
//   sub_type: 0=C>A, 1=C>G, 2=C>T, 3=T>A, 4=T>C, 5=T>G
//   ctx_idx:  5p_base_idx * 4 + 3p_base_idx  (A=0, C=1, G=2, T=3)
//
// Channel indices for key contexts used below:
//   C>T in ACG: sub=2, 5p=A(0), 3p=G(2) => 2*16 + 0*4+2 = 34
//   C>T in CCG: sub=2, 5p=C(1), 3p=G(2) => 2*16 + 1*4+2 = 38
//   C>T in GCG: sub=2, 5p=G(2), 3p=G(2) => 2*16 + 2*4+2 = 42
//   C>T in TCG: sub=2, 5p=T(3), 3p=G(2) => 2*16 + 3*4+2 = 46
//   C>T in TCA: sub=2, 5p=T(3), 3p=A(0) => 2*16 + 3*4+0 = 44
//   C>T in TCT: sub=2, 5p=T(3), 3p=T(3) => 2*16 + 3*4+3 = 47
//   C>G in TCA: sub=1, 5p=T(3), 3p=A(0) => 1*16 + 3*4+0 = 28
//   C>G in TCT: sub=1, 5p=T(3), 3p=T(3) => 1*16 + 3*4+3 = 31
//   T>G in CTT: sub=5, 5p=C(1), 3p=T(3) => 5*16 + 1*4+3 = 87
// ---------------------------------------------------------------------------

/// SBS1: spontaneous deamination at CpG sites (C>T in xCG context).
///
/// Concentrated at C>T changes where the 3' neighbour is G (CpG context).
/// ~10% each for ACG, CCG, GCG, TCG contexts; remaining weight spread
/// uniformly across all other channels.
pub static SBS1: SignatureVec = {
    // Four CpG channels each get 0.10; remaining 0.60 spread over 92 channels.
    const CG_WEIGHT: f64 = 0.10;
    const REMAINDER: f64 = 1.0 - 4.0 * CG_WEIGHT;
    const OTHER: f64 = REMAINDER / 92.0;

    let mut v = [OTHER; 96];
    // C>T in ACG: channel 34
    v[34] = CG_WEIGHT;
    // C>T in CCG: channel 38
    v[38] = CG_WEIGHT;
    // C>T in GCG: channel 42
    v[42] = CG_WEIGHT;
    // C>T in TCG: channel 46
    v[46] = CG_WEIGHT;
    v
};

/// SBS2: APOBEC signature (C>T at TCA and TCT contexts).
///
/// ~40% at C>T in TCA, ~40% at C>T in TCT, remaining 20% spread uniformly.
pub static SBS2: SignatureVec = {
    const TCA_WEIGHT: f64 = 0.40;
    const TCT_WEIGHT: f64 = 0.40;
    const REMAINDER: f64 = 1.0 - TCA_WEIGHT - TCT_WEIGHT;
    const OTHER: f64 = REMAINDER / 94.0;

    let mut v = [OTHER; 96];
    // C>T in TCA: channel 44
    v[44] = TCA_WEIGHT;
    // C>T in TCT: channel 47
    v[47] = TCT_WEIGHT;
    v
};

/// SBS13: APOBEC signature (C>G at TCA and TCT contexts).
///
/// ~40% at C>G in TCA, ~40% at C>G in TCT, remaining 20% spread uniformly.
pub static SBS13: SignatureVec = {
    const TCA_WEIGHT: f64 = 0.40;
    const TCT_WEIGHT: f64 = 0.40;
    const REMAINDER: f64 = 1.0 - TCA_WEIGHT - TCT_WEIGHT;
    const OTHER: f64 = REMAINDER / 94.0;

    let mut v = [OTHER; 96];
    // C>G in TCA: channel 28
    v[28] = TCA_WEIGHT;
    // C>G in TCT: channel 31
    v[31] = TCT_WEIGHT;
    v
};

/// SBS4: tobacco smoking signature (C>A in various GCN contexts).
///
/// ~60% concentrated at C>A in GCA, GCC, GCG, GCT contexts (15% each),
/// remaining 40% spread uniformly.
pub static SBS4: SignatureVec = {
    // C>A channels where 5p=G: sub=0, 5p=G(2)
    //   GCA: 0*16 + 2*4+0 = 8
    //   GCC: 0*16 + 2*4+1 = 9
    //   GCG: 0*16 + 2*4+2 = 10
    //   GCT: 0*16 + 2*4+3 = 11
    const GCN_WEIGHT: f64 = 0.15;
    const REMAINDER: f64 = 1.0 - 4.0 * GCN_WEIGHT;
    const OTHER: f64 = REMAINDER / 92.0;

    let mut v = [OTHER; 96];
    v[8] = GCN_WEIGHT; // C>A in GCA
    v[9] = GCN_WEIGHT; // C>A in GCC
    v[10] = GCN_WEIGHT; // C>A in GCG
    v[11] = GCN_WEIGHT; // C>A in GCT
    v
};

/// SBS5: clock-like ageing signature.
///
/// Relatively flat across all 96 channels; uniform distribution.
pub static SBS5: SignatureVec = [1.0 / 96.0; 96];

/// SBS7a: UV radiation signature (C>T in CCN contexts).
///
/// ~15% each at C>T in CCA, CCC, CCG, CCT contexts (60% total),
/// remaining 40% spread uniformly.
pub static SBS7A: SignatureVec = {
    // C>T channels where 5p=C: sub=2, 5p=C(1)
    //   CCA: 2*16 + 1*4+0 = 36
    //   CCC: 2*16 + 1*4+1 = 37
    //   CCG: 2*16 + 1*4+2 = 38
    //   CCT: 2*16 + 1*4+3 = 39
    const CC_WEIGHT: f64 = 0.15;
    const REMAINDER: f64 = 1.0 - 4.0 * CC_WEIGHT;
    const OTHER: f64 = REMAINDER / 92.0;

    let mut v = [OTHER; 96];
    v[36] = CC_WEIGHT; // C>T in CCA
    v[37] = CC_WEIGHT; // C>T in CCC
    v[38] = CC_WEIGHT; // C>T in CCG
    v[39] = CC_WEIGHT; // C>T in CCT
    v
};

/// SBS17a: 5-FU / unknown aetiology (T>G in CTT context).
///
/// ~30% at T>G in CTT, ~40% spread across other T>G contexts, ~30% across
/// remaining channels. Values normalised to sum to 1.0.
pub static SBS17A: SignatureVec = {
    // T>G channels: sub=5, ctx_idx = 5p*4 + 3p
    //   CTT: 5*16 + 1*4+3 = 87
    // Other T>G channels (15 of them) each get 40%/15 ≈ 0.0267
    // Remaining 46 non-T>G channels each get 30%/46 ≈ 0.0065
    const CTT_WEIGHT: f64 = 0.30;
    const OTHER_TG_WEIGHT: f64 = 0.40 / 15.0;
    const BASE_WEIGHT: f64 = 0.30 / 80.0;

    let mut v = [BASE_WEIGHT; 96];
    // Set all T>G channels to OTHER_TG_WEIGHT first.
    let mut i = 80; // 5 * 16
    while i < 96 {
        v[i] = OTHER_TG_WEIGHT;
        i += 1;
    }
    // Override the primary channel.
    v[87] = CTT_WEIGHT;
    v
};

/// SBS3: BRCA / homologous recombination deficiency signature.
///
/// Spread across all 96 channels with a modest elevation at C>T changes,
/// reflecting the genomic instability phenotype seen in HRD tumours.
pub static SBS3: SignatureVec = {
    // Base weight slightly elevated at C>T (channels 32–47).
    const CT_WEIGHT: f64 = 0.013;
    const BASE_WEIGHT: f64 = (1.0 - 16.0 * CT_WEIGHT) / 80.0;

    let mut v = [BASE_WEIGHT; 96];
    let mut i = 32; // 2 * 16 (C>T block)
    while i < 48 {
        v[i] = CT_WEIGHT;
        i += 1;
    }
    v
};

#[cfg(test)]
mod tests {
    use super::*;

    // Tolerance for floating-point sum checks.
    const TOL: f64 = 1e-10;

    fn sum_sig(sig: &SignatureVec) -> f64 {
        sig.iter().sum()
    }

    #[test]
    fn test_sbs1_sums_to_one() {
        let s = sum_sig(&SBS1);
        assert!((s - 1.0).abs() < TOL, "SBS1 sums to {}", s);
    }

    #[test]
    fn test_sbs2_sums_to_one() {
        let s = sum_sig(&SBS2);
        assert!((s - 1.0).abs() < TOL, "SBS2 sums to {}", s);
    }

    #[test]
    fn test_sbs3_sums_to_one() {
        let s = sum_sig(&SBS3);
        assert!((s - 1.0).abs() < TOL, "SBS3 sums to {}", s);
    }

    #[test]
    fn test_sbs4_sums_to_one() {
        let s = sum_sig(&SBS4);
        assert!((s - 1.0).abs() < TOL, "SBS4 sums to {}", s);
    }

    #[test]
    fn test_sbs5_sums_to_one() {
        let s = sum_sig(&SBS5);
        assert!((s - 1.0).abs() < TOL, "SBS5 sums to {}", s);
    }

    #[test]
    fn test_sbs7a_sums_to_one() {
        let s = sum_sig(&SBS7A);
        assert!((s - 1.0).abs() < TOL, "SBS7a sums to {}", s);
    }

    #[test]
    fn test_sbs13_sums_to_one() {
        let s = sum_sig(&SBS13);
        assert!((s - 1.0).abs() < TOL, "SBS13 sums to {}", s);
    }

    #[test]
    fn test_sbs17a_sums_to_one() {
        let s = sum_sig(&SBS17A);
        assert!((s - 1.0).abs() < TOL, "SBS17a sums to {}", s);
    }

    #[test]
    fn test_all_signatures_non_negative() {
        let sigs: &[&SignatureVec] = &[&SBS1, &SBS2, &SBS3, &SBS4, &SBS5, &SBS7A, &SBS13, &SBS17A];
        for sig in sigs {
            for &v in sig.iter() {
                assert!(v >= 0.0, "negative probability {}", v);
            }
        }
    }

    #[test]
    fn test_channel_index_ct_tca() {
        // C>T (sub=2) in TCA context: 5p=T, 3p=A
        // Expected: 2*16 + base_to_idx(T)*4 + base_to_idx(A) = 32 + 3*4 + 0 = 44
        let idx = channel_index(2, b'T', b'A');
        assert_eq!(idx, 44, "C>T in TCA should be channel 44");
    }

    #[test]
    fn test_channel_index_ct_tcg() {
        // C>T (sub=2) in TCG context: 5p=T, 3p=G
        // Expected: 2*16 + 3*4 + 2 = 46
        let idx = channel_index(2, b'T', b'G');
        assert_eq!(idx, 46, "C>T in TCG should be channel 46");
    }

    #[test]
    fn test_channel_index_tg_ctt() {
        // T>G (sub=5) in CTT context: 5p=C, 3p=T
        // Expected: 5*16 + 1*4 + 3 = 87
        let idx = channel_index(5, b'C', b'T');
        assert_eq!(idx, 87, "T>G in CTT should be channel 87");
    }

    #[test]
    fn test_sbs96_index_pyrimidine_direct() {
        // C>T at position with 5p=T, 3p=A (TCA context) => channel 44
        let idx = sbs96_index(b'T', b'C', b'T', b'A');
        assert_eq!(idx, Some(44));
    }

    #[test]
    fn test_sbs96_index_purine_complement() {
        // G>A with 5p=T, 3p=A is the complement of C>T in TCG context.
        // Reverse complement: 5p_rc = complement(A) = T, 3p_rc = complement(T) = A
        // ref_rc = complement(G) = C, alt_rc = complement(A) = T => C>T
        // So this maps to C>T in T_C_A context (5p=T, 3p=A) = channel 44.
        let idx = sbs96_index(b'T', b'G', b'A', b'A');
        // 5p=T, ref=G, alt=A, 3p=A
        // r_rc=C, a_rc=T => C>T (sub=2)
        // ctx swapped: new 5p=complement(3p)=complement(A)=T, new 3p=complement(5p)=complement(T)=A
        // channel = 2*16 + base_to_idx(T)*4 + base_to_idx(A) = 32 + 12 + 0 = 44
        assert_eq!(idx, Some(44));
    }

    #[test]
    fn test_sbs96_index_invalid_substitution() {
        // C>C is not a valid SNV.
        let idx = sbs96_index(b'A', b'C', b'C', b'T');
        assert_eq!(idx, None);
    }

    #[test]
    fn test_trinucleotide_context_middle() {
        let seq = b"ACGT";
        assert_eq!(trinucleotide_context(seq, 1), Some((b'A', b'C', b'G')));
        assert_eq!(trinucleotide_context(seq, 2), Some((b'C', b'G', b'T')));
    }

    #[test]
    fn test_trinucleotide_context_boundary() {
        let seq = b"ACGT";
        // Position 0 has no 5' neighbour.
        assert_eq!(trinucleotide_context(seq, 0), None);
        // Last position has no 3' neighbour.
        assert_eq!(trinucleotide_context(seq, 3), None);
    }

    #[test]
    fn test_variant_sbs96_index_valid() {
        // Sequence: TACG; variant at pos 1 (ref=A, but we use ref_base directly).
        // Trinucleotide: 5p=T, central=A, 3p=C.
        // A>T: ref=A => purine, take RC: ref_rc=T, alt_rc=A => T>A (sub=3)
        // ctx_rc: 5p=complement(C)=G, 3p=complement(T)=A
        // channel = 3*16 + base_to_idx(G)*4 + base_to_idx(A) = 48 + 2*4 + 0 = 56
        let seq = b"TACG";
        let idx = variant_sbs96_index(seq, 1, b'A', b'T');
        assert!(idx.is_some());
    }

    #[test]
    fn test_sbs1_cpg_channels_dominant() {
        // CpG C>T channels (34, 38, 42, 46) should each be 0.10.
        assert!((SBS1[34] - 0.10).abs() < 1e-12);
        assert!((SBS1[38] - 0.10).abs() < 1e-12);
        assert!((SBS1[42] - 0.10).abs() < 1e-12);
        assert!((SBS1[46] - 0.10).abs() < 1e-12);
    }

    #[test]
    fn test_sbs2_tca_tct_dominant() {
        // C>T in TCA (44) and TCT (47) should be 0.40 each.
        assert!((SBS2[44] - 0.40).abs() < 1e-12);
        assert!((SBS2[47] - 0.40).abs() < 1e-12);
    }

    #[test]
    fn test_sbs13_tca_tct_dominant() {
        // C>G in TCA (28) and TCT (31) should be 0.40 each.
        assert!((SBS13[28] - 0.40).abs() < 1e-12);
        assert!((SBS13[31] - 0.40).abs() < 1e-12);
    }

    #[test]
    fn test_sbs17a_ctt_dominant() {
        // T>G in CTT (channel 87) should be 0.30.
        assert!((SBS17A[87] - 0.30).abs() < 1e-12);
    }

    #[test]
    fn test_builtin_signature_lookup() {
        assert!(builtin_signature("SBS1").is_some());
        assert!(builtin_signature("SBS5").is_some());
        assert!(builtin_signature("UNKNOWN").is_none());
    }

    #[test]
    fn test_builtin_signature_names_complete() {
        let names = builtin_signature_names();
        assert_eq!(names.len(), 8);
        for name in names {
            assert!(builtin_signature(name).is_some(), "name {} not found", name);
        }
    }
}
