# EPIC-SIGNATURES: COSMIC Mutational Signatures

## Goal
Enable mutation generation weighted by COSMIC SBS96 mutational signatures, so simulated variants reflect the trinucleotide context patterns of real cancer types.

## Motivation
Benchmarking somatic callers on biologically realistic mutation spectra is more informative than uniform random mutations. COSMIC signatures are the standard reference for cancer mutation patterns.

## Scope
- In scope: embedded SBS96 matrix (COSMIC v3.4), trinucleotide context lookup, signature-weighted alt base selection.
- Out of scope: DBS and ID signature channels, signature fitting or deconvolution, per-clone signature assignment.

## Tasks
- T029: Embed SBS96 signature matrix as compressed static asset
- T030: Implement trinucleotide context lookup from reference sequence
- T031: Signature-weighted alt base selection in random mutation generator
