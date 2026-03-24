# VarForge: Deferred Features

Features explicitly out of scope for v1. Do not implement or document as supported.

## Native long-read error models

Basic read-length support for PacBio and Nanopore exists. Deep platform-specific error modelling (kinetics, homopolymer errors, native basecalling artefacts) is deferred. Users who need this should use NanoSim or PBSIM3.

## Contamination simulation

The `ContaminationConfig` struct exists in the config schema but no simulation code backs it. Remove from the v1 config entirely. Reintroduce in a future version with a real implementation.

## FragmentModel::Custom

Dead code. Remove for v1. Users who need a custom fragment model can use the existing parametric options. A proper custom model can be added post-release.

## Chromosome name remapping

No aliasing between `chr1` and `1` naming conventions. Users must ensure the reference FASTA, input VCF, and config all use the same naming convention. Document this as a user requirement.

## Complex structural variants

Chromothripsis, templated insertions, and complex rearrangements are out of scope. Basic SVs (deletions, inversions, duplications, translocations) are in scope.

## Inline UMI mode

The config field exists but the feature is not implemented. Document explicitly as not supported in v1. Remove the config field or make it error on use.
