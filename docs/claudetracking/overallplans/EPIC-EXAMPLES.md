# EPIC-EXAMPLES: Example Config Accuracy

## Goal

Make all example YAML configs in examples/ runnable out-of-the-box and free of misleading
comments.

## Motivation

Every example hardcodes /data/ref/hg38.fa, which will not exist for any user. Misleading
comments about VCF compression requirements and mode descriptions reduce trust in the
examples and produce confusing error messages.

## Scope

- In scope: path templating, comment corrections, mode description accuracy.
- Out of scope: adding new example configs (separate task if desired).

## Tasks

- T100: Replace hardcoded reference paths in all examples with ${reference} placeholder
- T101: Fix misleading VCF compression comment in custom_mutations.yaml
- T102: Fix tumor_normal.yaml mode description (samples vs paired mode)
