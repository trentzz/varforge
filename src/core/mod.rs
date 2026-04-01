//! Core simulation primitives: types, coverage, fragment sampling, quality models, and the read engine.

pub mod capture;
pub mod coverage;
pub mod end_motifs;
pub mod engine;
pub mod error_orchestrator;
pub mod error_profile;
pub mod fragment;
pub mod gc_bias;
pub mod multi_sample;
pub mod profile_learner;
pub mod quality;
pub mod seq_errors;
pub mod types;
