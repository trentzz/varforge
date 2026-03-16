//! BAM editing engine for spiking variants into existing sequencing data.
//!
//! The editor reads a real BAM file, spikes in SNVs and indels at specified
//! positions using stochastic VAF sampling (binomial), and writes a modified
//! BAM. Unlike BAMSurgeon, every variant at a given VAF produces realistic
//! read-count scatter rather than a fixed deterministic count.

pub mod bam_editor;
pub mod read_modifier;
