//! VarForge: synthetic cancer sequencing data generator.
//!
//! Generates realistic FASTQ and BAM files with controlled mutations, tumour parameters,
//! UMI tags, and cfDNA fragment profiles for benchmarking bioinformatics tools.

pub mod artifacts;
pub mod cli;
pub mod core;
pub mod editor;
pub mod io;
pub mod seq_utils;
pub mod tumour;
pub mod umi;
pub mod variants;
