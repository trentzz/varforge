# EPIC-CONFIG-VALIDATION: Config Validation Hardening

## Goal

Add missing validation checks so that misconfigured runs fail fast with actionable error
messages rather than producing silently wrong output.

## Motivation

Several combinations of config values produce nonsense output without any error: reads longer
than fragments, UMIs longer than reads, malformed copy number regions, missing VCF paths.
These all pass validate() today.

## Scope

- In scope: cross-field checks, file existence checks, range checks, warning for suspicious
  combinations.
- Out of scope: full semantic validation of VCF contents, FASTA format checking.

## Tasks

- T083: Validate read_length <= fragment mean; error if not
- T084: Validate capture BED path exists when capture mode is set
- T085: Validate UMI length < read_length
- T086: Validate copy_number region strings are parseable at config load time
- T087: Warn when purity=1.0 with cfDNA fragment model
- T088: Validate mutation VCF path exists at config load time
