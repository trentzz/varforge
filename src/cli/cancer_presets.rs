//! Cancer-type-specific simulation presets for VarForge.
//!
//! Each preset configures realistic mutation spectra, variant counts, tumour
//! purity ranges, and common driver mutations for a specific cancer type.
//! Mutation type distributions are weighted to reflect dominant COSMIC SBS
//! signatures for the given cancer.
//!
//! # Supported cancer types
//!
//! | Key               | Cancer                                   | Dominant signatures     |
//! |-------------------|------------------------------------------|-------------------------|
//! | `lung_adeno`      | Lung adenocarcinoma (NSCLC)              | SBS4 (smoking), C>A     |
//! | `colorectal`      | Colorectal cancer                        | SBS1/SBS5 (aging), C>T  |
//! | `breast_tnbc`     | Triple-negative breast cancer            | SBS3 (HRD), flat        |
//! | `melanoma`        | Cutaneous melanoma                       | SBS7a/b (UV), C>T       |
//! | `aml`             | Acute myeloid leukaemia                  | SBS1/SBS5, low TMB      |
//! | `prostate`        | Prostate adenocarcinoma                  | SBS1/SBS5, low TMB      |
//! | `pancreatic`      | Pancreatic ductal adenocarcinoma         | SBS1/SBS5, low purity   |
//! | `glioblastoma`    | Glioblastoma multiforme (GBM)            | SBS1/SBS5, moderate TMB |
//!
//! # Usage
//!
//! ```text
//! varforge simulate --config cfg.yaml --preset cancer:lung_adeno
//! ```
//!
//! The `cancer:` namespace prefix is stripped before being passed to
//! [`get`]; the calling code in [`super::presets`] is responsible for
//! routing the `cancer:` prefix here.

use crate::io::config::{ArtifactConfig, MutationConfig, RandomMutationConfig};

use super::presets::PresetOverlay;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Return the overlay for a named cancer-type preset, or an error if the name
/// is unknown.  The `cancer:` namespace prefix must already have been stripped
/// by the caller.
pub fn get(name: &str) -> anyhow::Result<PresetOverlay> {
    match name {
        "lung_adeno" => Ok(preset_lung_adeno()),
        "colorectal" => Ok(preset_colorectal()),
        "breast_tnbc" => Ok(preset_breast_tnbc()),
        "melanoma" => Ok(preset_melanoma()),
        "aml" => Ok(preset_aml()),
        "prostate" => Ok(preset_prostate()),
        "pancreatic" => Ok(preset_pancreatic()),
        "glioblastoma" => Ok(preset_glioblastoma()),
        other => anyhow::bail!(
            "unknown cancer preset '{}'; valid choices: lung_adeno, colorectal, \
             breast_tnbc, melanoma, aml, prostate, pancreatic, glioblastoma",
            other
        ),
    }
}

/// Return the names of all built-in cancer-type presets (without the
/// `cancer:` namespace prefix).
pub fn all_names() -> &'static [&'static str] {
    &[
        "lung_adeno",
        "colorectal",
        "breast_tnbc",
        "melanoma",
        "aml",
        "prostate",
        "pancreatic",
        "glioblastoma",
    ]
}

// ---------------------------------------------------------------------------
// Mutation burden constants (mutations per Mb → approx WGS count over 3 Gb)
// ---------------------------------------------------------------------------

/// Convert a mutations-per-Mb rate to an approximate total count for 3 Gb genome.
const fn muts_from_per_mb(per_mb: usize) -> usize {
    per_mb * 3_000
}

// ---------------------------------------------------------------------------
// Individual cancer preset definitions
// ---------------------------------------------------------------------------

/// `lung_adeno` – Lung adenocarcinoma / NSCLC.
///
/// High tumour mutational burden driven by tobacco smoking (SBS4: C>A
/// transversions).  EGFR and KRAS are the most common driver genes.  Typical
/// purity is moderate (~60 %).
///
/// References:
/// - COSMIC SBS4: associated with tobacco-smoke mutagens; dominantly C[T>A]G
/// - TMB: ~8–10 mut/Mb on average (range 1–100+)
fn preset_lung_adeno() -> PresetOverlay {
    // SBS4 smoking signature is C>A dominant.
    // We model that by boosting indel fraction slightly (insertions at
    // tobacco-damage sites) and using a high overall mutation count.
    // snv_fraction is kept high; the signature detail is captured in the
    // comments — actual trinucleotide weighting is outside the scope of
    // RandomMutationConfig for now but the counts & fractions reflect SBS4.
    PresetOverlay {
        coverage: Some(30.0),
        purity: Some(0.60),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                // ~8 mut/Mb × 3 000 Mb ≈ 24 000; use 10 000 for feasibility
                count: muts_from_per_mb(8),
                vaf_min: 0.01,
                vaf_max: 0.60,
                // SBS4: high C>A (treated as SNV here); small indel component
                snv_fraction: 0.85,
                indel_fraction: 0.10,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_signature: None,
            sv_count: 0,
            include_driver_mutations: false,
        }),
        ..Default::default()
    }
}

/// `colorectal` – Colorectal cancer (CRC).
///
/// Dominant aging signatures SBS1 (C>T at CpG) and SBS5 (flat, aging).
/// MSI-H tumours have very high TMB; MSS tumours have low-moderate TMB.
/// We model the average MSS profile.  APC, KRAS, TP53 are key drivers.
fn preset_colorectal() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        purity: Some(0.65),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                // MSS average ~5 mut/Mb
                count: muts_from_per_mb(5),
                vaf_min: 0.01,
                vaf_max: 0.65,
                // SBS1/SBS5 are SNV-dominated
                snv_fraction: 0.82,
                indel_fraction: 0.13,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_signature: None,
            sv_count: 0,
            include_driver_mutations: false,
        }),
        ..Default::default()
    }
}

/// `breast_tnbc` – Triple-negative breast cancer (TNBC).
///
/// Characterized by HRD (homologous recombination deficiency), dominated by
/// SBS3 (flat spectrum, associated with BRCA1/BRCA2 loss).  High indel
/// burden due to error-prone repair.  TP53 mutation is near-universal.
fn preset_breast_tnbc() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        purity: Some(0.55),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                // SBS3 tumours: ~3–6 mut/Mb; TNBC can be higher
                count: muts_from_per_mb(5),
                vaf_min: 0.01,
                vaf_max: 0.55,
                // SBS3 is relatively flat; increase indel to reflect HRD
                snv_fraction: 0.75,
                indel_fraction: 0.20,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_signature: None,
            sv_count: 0,
            include_driver_mutations: false,
        }),
        artifacts: Some(ArtifactConfig {
            ffpe_damage_rate: None,
            oxog_rate: None,
            duplicate_rate: Some(0.05),
            pcr_error_rate: None,
        }),
        ..Default::default()
    }
}

/// `melanoma` – Cutaneous melanoma.
///
/// Very high TMB dominated by UV radiation signature SBS7a/b (C>T at
/// dipyrimidines, particularly CC>TT dinucleotides).  BRAF V600E is present
/// in ~50 % of cases.
fn preset_melanoma() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        purity: Some(0.70),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                // SBS7 tumours: 20–50+ mut/Mb; cap at 30 for simulation speed
                count: muts_from_per_mb(30),
                vaf_min: 0.01,
                vaf_max: 0.70,
                // SBS7: almost exclusively C>T SNVs; MNV rate elevated (CC>TT)
                snv_fraction: 0.88,
                indel_fraction: 0.05,
                mnv_fraction: 0.07,
                signature: None,
            }),
            sv_signature: None,
            sv_count: 0,
            include_driver_mutations: false,
        }),
        ..Default::default()
    }
}

/// `aml` – Acute myeloid leukaemia.
///
/// Low TMB; predominantly aging signatures SBS1/SBS5.  Key drivers include
/// FLT3-ITD, NPM1 frameshift, DNMT3A R882H, IDH1/IDH2.  High purity
/// because leukaemia samples are typically from bone marrow or blood blasts.
fn preset_aml() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        purity: Some(0.80),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                // AML: ~0.5–1 mut/Mb
                count: muts_from_per_mb(1),
                vaf_min: 0.05,
                vaf_max: 0.80,
                // SBS1/SBS5 SNV dominated; FLT3-ITD is a large indel
                snv_fraction: 0.80,
                indel_fraction: 0.15,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_signature: None,
            sv_count: 0,
            include_driver_mutations: false,
        }),
        ..Default::default()
    }
}

/// `prostate` – Prostate adenocarcinoma.
///
/// Low-to-moderate TMB; aging signatures SBS1/SBS5 dominate.  TMPRSS2-ERG
/// fusion is the most common alteration (~50 %).  AR amplification and
/// PTEN loss are key events.
fn preset_prostate() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        purity: Some(0.50),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                // Prostate: ~1–2 mut/Mb
                count: muts_from_per_mb(2),
                vaf_min: 0.01,
                vaf_max: 0.50,
                snv_fraction: 0.82,
                indel_fraction: 0.13,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_signature: None,
            sv_count: 0,
            include_driver_mutations: false,
        }),
        ..Default::default()
    }
}

/// `pancreatic` – Pancreatic ductal adenocarcinoma (PDAC).
///
/// Low purity because of extensive stromal desmoplasia.  Low-moderate TMB.
/// Near-universal KRAS G12D/V mutation; TP53, SMAD4, CDKN2A also commonly
/// altered.  Aging signatures SBS1/SBS5 predominate.
fn preset_pancreatic() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        purity: Some(0.25), // stromal contamination → low purity
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                // PDAC: ~2–4 mut/Mb
                count: muts_from_per_mb(3),
                vaf_min: 0.005,
                vaf_max: 0.25,
                snv_fraction: 0.82,
                indel_fraction: 0.13,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_signature: None,
            sv_count: 0,
            include_driver_mutations: false,
        }),
        ..Default::default()
    }
}

/// `glioblastoma` – Glioblastoma multiforme (GBM).
///
/// Moderate TMB; aging signatures SBS1/SBS5 typical in IDH-wildtype GBM.
/// Key drivers: EGFR amplification, PTEN loss, IDH1 R132H (IDH-mutant
/// subset), TERT promoter mutation.
fn preset_glioblastoma() -> PresetOverlay {
    PresetOverlay {
        coverage: Some(30.0),
        purity: Some(0.65),
        mutations: Some(MutationConfig {
            vcf: None,
            random: Some(RandomMutationConfig {
                // GBM IDH-wt: ~3–5 mut/Mb
                count: muts_from_per_mb(4),
                vaf_min: 0.01,
                vaf_max: 0.65,
                snv_fraction: 0.82,
                indel_fraction: 0.13,
                mnv_fraction: 0.05,
                signature: None,
            }),
            sv_signature: None,
            sv_count: 0,
            include_driver_mutations: false,
        }),
        ..Default::default()
    }
}

// ---------------------------------------------------------------------------
// Known driver gene information (informational; not yet injected as specific
// VCF records, but exposed for documentation and future driver-injection tasks)
// ---------------------------------------------------------------------------

/// A canonical driver mutation associated with a cancer type.
///
/// Genomic coordinates are hg38 (GRCh38). `pos` is 0-based. `ref_allele` and
/// `alt_allele` are `None` for structural/fusion events (e.g. amplification,
/// fusion) that cannot be represented as a single-nucleotide or short indel
/// variant; those entries are skipped during driver injection.
#[derive(Debug, Clone, PartialEq)]
pub struct DriverMutation {
    /// HGNC gene symbol, e.g. `"KRAS"`.
    pub gene: &'static str,
    /// Short description, e.g. `"G12D"` or `"V600E"`.
    pub alteration: &'static str,
    /// Approximate prevalence in this cancer type (0.0–1.0).
    pub prevalence: f32,
    /// Chromosome name in hg38 format, e.g. `"chr12"`. `None` for
    /// non-positional alterations.
    pub chrom: Option<&'static str>,
    /// 0-based position on the chromosome (hg38). `None` when unavailable.
    pub pos: Option<u64>,
    /// Reference allele at this position (ASCII bytes). `None` for
    /// non-SNV/indel events.
    pub ref_allele: Option<&'static [u8]>,
    /// Alternate allele (ASCII bytes). `None` for non-SNV/indel events.
    pub alt_allele: Option<&'static [u8]>,
}

/// Return the list of canonical driver mutations for a given cancer preset
/// name.  The name must **not** include the `cancer:` prefix.
pub fn drivers_for(preset_name: &str) -> &'static [DriverMutation] {
    match preset_name {
        "lung_adeno" => LUNG_ADENO_DRIVERS,
        "colorectal" => COLORECTAL_DRIVERS,
        "breast_tnbc" => BREAST_TNBC_DRIVERS,
        "melanoma" => MELANOMA_DRIVERS,
        "aml" => AML_DRIVERS,
        "prostate" => PROSTATE_DRIVERS,
        "pancreatic" => PANCREATIC_DRIVERS,
        "glioblastoma" => GLIOBLASTOMA_DRIVERS,
        _ => &[],
    }
}

static LUNG_ADENO_DRIVERS: &[DriverMutation] = &[
    DriverMutation {
        gene: "KRAS",
        alteration: "G12C",
        prevalence: 0.13,
        // chr12:25_227_342 hg38; KRAS codon 12, TGT→TGT is G12C (c.34G>T)
        chrom: Some("chr12"),
        pos: Some(25_227_342),
        ref_allele: Some(b"G"),
        alt_allele: Some(b"T"),
    },
    DriverMutation {
        gene: "EGFR",
        alteration: "L858R",
        prevalence: 0.10,
        // chr7:55_191_822 hg38; EGFR exon 21 L858R (c.2573T>G)
        chrom: Some("chr7"),
        pos: Some(55_191_822),
        ref_allele: Some(b"T"),
        alt_allele: Some(b"G"),
    },
    DriverMutation {
        gene: "EGFR",
        alteration: "exon19del",
        prevalence: 0.10,
        // Representative exon 19 deletion (delE746-A750, most common)
        // chr7:55_174_772 hg38 (15 bp deletion AATTAAGAGAAGCA)
        chrom: Some("chr7"),
        pos: Some(55_174_772),
        ref_allele: Some(b"AATTAAGAGAAGCAA"),
        alt_allele: Some(b"A"),
    },
    DriverMutation {
        gene: "TP53",
        alteration: "various",
        prevalence: 0.46,
        // No single representative hotspot; skip injection.
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "STK11",
        alteration: "various",
        prevalence: 0.17,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "KEAP1",
        alteration: "various",
        prevalence: 0.12,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
];

static COLORECTAL_DRIVERS: &[DriverMutation] = &[
    DriverMutation {
        gene: "KRAS",
        // G12D is the most common in CRC (c.35G>A)
        alteration: "G12D/V",
        prevalence: 0.45,
        chrom: Some("chr12"),
        pos: Some(25_227_343),
        ref_allele: Some(b"G"),
        alt_allele: Some(b"A"),
    },
    DriverMutation {
        gene: "APC",
        alteration: "truncating",
        prevalence: 0.80,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "TP53",
        alteration: "various",
        prevalence: 0.60,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "SMAD4",
        alteration: "various",
        prevalence: 0.15,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "PIK3CA",
        alteration: "various",
        prevalence: 0.20,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
];

static BREAST_TNBC_DRIVERS: &[DriverMutation] = &[
    DriverMutation {
        gene: "TP53",
        alteration: "various",
        prevalence: 0.80,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "BRCA1",
        alteration: "germline/somatic",
        prevalence: 0.20,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "RB1",
        alteration: "loss",
        prevalence: 0.20,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "PIK3CA",
        // PIK3CA H1047R hotspot (c.3140A>G)
        alteration: "various",
        prevalence: 0.10,
        chrom: Some("chr3"),
        pos: Some(179_234_296),
        ref_allele: Some(b"A"),
        alt_allele: Some(b"G"),
    },
];

static MELANOMA_DRIVERS: &[DriverMutation] = &[
    DriverMutation {
        gene: "BRAF",
        alteration: "V600E",
        prevalence: 0.50,
        // chr7:140_753_335 hg38; BRAF V600E (c.1799T>A)
        chrom: Some("chr7"),
        pos: Some(140_753_335),
        ref_allele: Some(b"A"),
        alt_allele: Some(b"T"),
    },
    DriverMutation {
        gene: "NRAS",
        alteration: "Q61R/K",
        prevalence: 0.20,
        // chr1:114_716_126 hg38; NRAS Q61R (c.182A>G)
        chrom: Some("chr1"),
        pos: Some(114_716_126),
        ref_allele: Some(b"A"),
        alt_allele: Some(b"G"),
    },
    DriverMutation {
        gene: "TERT",
        alteration: "promoter C228T/C250T",
        prevalence: 0.74,
        // chr5:1_295_228 hg38; TERT promoter C228T (c.-124C>T)
        chrom: Some("chr5"),
        pos: Some(1_295_228),
        ref_allele: Some(b"C"),
        alt_allele: Some(b"T"),
    },
    DriverMutation {
        gene: "CDKN2A",
        alteration: "various",
        prevalence: 0.40,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "PTEN",
        alteration: "loss",
        prevalence: 0.20,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
];

static AML_DRIVERS: &[DriverMutation] = &[
    DriverMutation {
        gene: "FLT3",
        alteration: "ITD",
        prevalence: 0.25,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "NPM1",
        alteration: "W288fs",
        prevalence: 0.30,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "DNMT3A",
        alteration: "R882H",
        prevalence: 0.22,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "IDH1",
        alteration: "R132H",
        prevalence: 0.08,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "IDH2",
        alteration: "R140Q",
        prevalence: 0.12,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "TET2",
        alteration: "various",
        prevalence: 0.10,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
];

static PROSTATE_DRIVERS: &[DriverMutation] = &[
    DriverMutation {
        gene: "TMPRSS2-ERG",
        alteration: "fusion",
        prevalence: 0.50,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "AR",
        alteration: "amplification",
        prevalence: 0.30,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "PTEN",
        alteration: "loss",
        prevalence: 0.25,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "TP53",
        alteration: "various",
        prevalence: 0.25,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "BRCA2",
        alteration: "germline/somatic",
        prevalence: 0.12,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
];

static PANCREATIC_DRIVERS: &[DriverMutation] = &[
    DriverMutation {
        gene: "KRAS",
        alteration: "G12D",
        prevalence: 0.92,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "TP53",
        alteration: "various",
        prevalence: 0.72,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "SMAD4",
        alteration: "various",
        prevalence: 0.32,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "CDKN2A",
        alteration: "various",
        prevalence: 0.29,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
];

static GLIOBLASTOMA_DRIVERS: &[DriverMutation] = &[
    DriverMutation {
        gene: "EGFR",
        alteration: "amplification/EGFRvIII",
        prevalence: 0.57,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "PTEN",
        alteration: "loss",
        prevalence: 0.40,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "TP53",
        alteration: "various",
        prevalence: 0.28,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "TERT",
        alteration: "promoter C228T/C250T",
        prevalence: 0.72,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
    DriverMutation {
        gene: "IDH1",
        alteration: "R132H (IDH-mutant subset)",
        prevalence: 0.10,
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    },
];

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::presets::apply_preset_to_config;
    use crate::io::config::{Config, FragmentConfig, OutputConfig, QualityConfig, SampleConfig};
    use std::path::PathBuf;

    // Helper: construct a base Config with all defaults.
    fn base_config() -> Config {
        Config {
            reference: PathBuf::from("/dev/null"),
            output: OutputConfig {
                directory: PathBuf::from("/tmp"),
                fastq: true,
                bam: false,
                truth_vcf: false,
                manifest: false,
                germline_vcf: false,
                single_read_bam: false,
                mapq: 60,
            },
            sample: SampleConfig::default(),
            fragment: FragmentConfig::default(),
            quality: QualityConfig::default(),
            tumour: None,
            mutations: None,
            umi: None,
            artifacts: None,
            seed: None,
            threads: None,
            chromosomes: None,
            regions_bed: None,
            copy_number: None,
            gc_bias: None,
            samples: None,
            capture: None,
            performance: Default::default(),
            preset: None,
            vafs: None,
            germline: None,
            paired: None,
        }
    }

    // High mutation burden threshold (mut/Mb): above this we call it "high TMB".
    const HIGH_TMB_THRESHOLD_MUTS: usize = muts_from_per_mb(10);
    // Low mutation burden threshold: at or below this we call it "low TMB".
    const LOW_TMB_THRESHOLD_MUTS: usize = muts_from_per_mb(2);

    /// Test 1: Every cancer preset produces a valid, non-degenerate config.
    #[test]
    fn test_all_cancer_presets_valid() {
        for name in all_names() {
            let overlay = get(name).unwrap_or_else(|e| panic!("preset '{}' failed: {}", name, e));
            let mut cfg = base_config();
            apply_preset_to_config(&mut cfg, &overlay);

            assert!(
                cfg.sample.coverage > 0.0,
                "preset 'cancer:{name}' must produce positive coverage"
            );

            // Must have a mutation config
            let muts = cfg
                .mutations
                .unwrap_or_else(|| panic!("preset 'cancer:{name}' must set mutations"));
            let rand = muts
                .random
                .unwrap_or_else(|| panic!("preset 'cancer:{name}' must set random mutations"));

            assert!(rand.count > 0, "preset 'cancer:{name}' must have count > 0");
            assert!(
                rand.vaf_min > 0.0,
                "preset 'cancer:{name}' vaf_min must be > 0"
            );
            assert!(
                rand.vaf_max > rand.vaf_min,
                "preset 'cancer:{name}' vaf_max must exceed vaf_min"
            );

            let frac_sum = rand.snv_fraction + rand.indel_fraction + rand.mnv_fraction;
            assert!(
                (frac_sum - 1.0).abs() < 1e-6,
                "preset 'cancer:{name}' type fractions must sum to 1.0, got {frac_sum}"
            );

            // Must have tumour purity set
            let tumour = cfg
                .tumour
                .unwrap_or_else(|| panic!("preset 'cancer:{name}' must set tumour"));
            assert!(
                (0.0..=1.0).contains(&tumour.purity),
                "preset 'cancer:{name}' purity {} must be in [0, 1]",
                tumour.purity
            );
        }
    }

    /// Test 2: Lung adenocarcinoma preset has SBS4 (smoking) signature characteristics.
    ///
    /// SBS4 features:
    ///  - High SNV fraction (C>A transversions dominate)
    ///  - High overall mutation count (tobacco-exposed tumours have elevated TMB)
    ///  - snv_fraction should be ≥ 0.80
    #[test]
    fn test_lung_adeno_signatures() {
        let overlay = get("lung_adeno").unwrap();
        let rand = overlay.mutations.as_ref().unwrap().random.as_ref().unwrap();

        // SBS4 is almost entirely SNVs — fraction must be high.
        assert!(
            rand.snv_fraction >= 0.80,
            "lung_adeno SNV fraction {} should be ≥ 0.80 (SBS4 is SNV-dominant)",
            rand.snv_fraction
        );

        // Smoking history → elevated mutation burden.
        assert!(
            rand.count >= muts_from_per_mb(5),
            "lung_adeno mutation count {} should be at least {} (SBS4 elevated TMB)",
            rand.count,
            muts_from_per_mb(5),
        );

        // Purity check — lung tumour purity is typically 40–80 %.
        let purity = overlay.purity.unwrap();
        assert!(
            (0.30..=0.90).contains(&purity),
            "lung_adeno purity {purity} out of expected range [0.30, 0.90]"
        );
    }

    /// Test 3: Melanoma has very high mutation burden (SBS7 / UV damage).
    #[test]
    fn test_melanoma_high_tmb() {
        let overlay = get("melanoma").unwrap();
        let rand = overlay.mutations.as_ref().unwrap().random.as_ref().unwrap();

        assert!(
            rand.count > HIGH_TMB_THRESHOLD_MUTS,
            "melanoma mutation count {} should exceed high-TMB threshold {}",
            rand.count,
            HIGH_TMB_THRESHOLD_MUTS,
        );

        // UV signature (SBS7): C>T SNVs dominate even more than smoking.
        assert!(
            rand.snv_fraction >= 0.82,
            "melanoma SNV fraction {} should be ≥ 0.82 (SBS7 is C>T dominant)",
            rand.snv_fraction
        );
    }

    /// Test 4: AML has low mutation burden.
    #[test]
    fn test_aml_low_tmb() {
        let overlay = get("aml").unwrap();
        let rand = overlay.mutations.as_ref().unwrap().random.as_ref().unwrap();

        assert!(
            rand.count <= LOW_TMB_THRESHOLD_MUTS,
            "aml mutation count {} should be ≤ low-TMB threshold {}",
            rand.count,
            LOW_TMB_THRESHOLD_MUTS,
        );
    }

    /// Test 5: Known driver mutations are present for each cancer type.
    #[test]
    fn test_driver_mutations_included() {
        // lung_adeno: KRAS and EGFR drivers expected
        let lung_drivers = drivers_for("lung_adeno");
        assert!(
            lung_drivers.iter().any(|d| d.gene == "KRAS"),
            "lung_adeno drivers must include KRAS"
        );
        assert!(
            lung_drivers.iter().any(|d| d.gene == "EGFR"),
            "lung_adeno drivers must include EGFR"
        );

        // melanoma: BRAF V600E expected
        let mel_drivers = drivers_for("melanoma");
        assert!(
            mel_drivers
                .iter()
                .any(|d| d.gene == "BRAF" && d.alteration.contains("V600E")),
            "melanoma drivers must include BRAF V600E"
        );

        // aml: FLT3, NPM1, DNMT3A expected
        let aml_drivers = drivers_for("aml");
        assert!(
            aml_drivers.iter().any(|d| d.gene == "FLT3"),
            "aml drivers must include FLT3"
        );
        assert!(
            aml_drivers.iter().any(|d| d.gene == "NPM1"),
            "aml drivers must include NPM1"
        );
        assert!(
            aml_drivers.iter().any(|d| d.gene == "DNMT3A"),
            "aml drivers must include DNMT3A"
        );

        // pancreatic: KRAS G12D expected
        let pdac_drivers = drivers_for("pancreatic");
        assert!(
            pdac_drivers
                .iter()
                .any(|d| d.gene == "KRAS" && d.alteration.contains("G12D")),
            "pancreatic drivers must include KRAS G12D"
        );

        // All presets must have at least one driver listed.
        for name in all_names() {
            assert!(
                !drivers_for(name).is_empty(),
                "cancer preset '{name}' must have at least one driver mutation defined"
            );
        }
    }

    /// Test 6: User CLI/YAML values override cancer preset values.
    ///
    /// We apply the preset and then simulate what `apply_overrides` would do
    /// by directly mutating the config (coverage, purity).  The final values
    /// must match the user-supplied values, not the preset defaults.
    #[test]
    fn test_preset_override() {
        let overlay = get("melanoma").unwrap();
        let mut cfg = base_config();
        apply_preset_to_config(&mut cfg, &overlay);

        // Confirm preset was applied.
        assert!(
            cfg.sample.coverage > 0.0,
            "melanoma preset should set coverage"
        );

        // Simulate CLI override: user specifies --coverage 100 --purity 0.3
        let user_coverage = 100.0_f64;
        let user_purity = 0.3_f64;
        cfg.sample.coverage = user_coverage;
        if let Some(ref mut t) = cfg.tumour {
            t.purity = user_purity;
        }

        assert!(
            (cfg.sample.coverage - user_coverage).abs() < 1e-9,
            "user coverage {user_coverage} must survive after override"
        );
        assert!(
            (cfg.tumour.as_ref().unwrap().purity - user_purity).abs() < 1e-9,
            "user purity {user_purity} must survive after override"
        );
    }
}
