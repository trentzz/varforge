# Migrate YAML parser: serde_yml to serde_yaml

**Category**: infra
**Related epic**: EPIC-RELEASE-GATE

## Context

`serde_yml` (v0.0.12) has two security advisories:
- RUSTSEC-2025-0068: serde_yml is unsound and unmaintained.
- RUSTSEC-2025-0067: libyml (transitive dep) is unsound and unmaintained.

This blocks any release that claims to be production-quality.

## Options

1. **Migrate to `serde_yaml` (v0.9)**: the canonical Rust YAML crate, actively maintained. API is nearly identical. Migration is mechanical.
2. **Migrate to `serde_yml` v0.1+**: if a new maintainer picks it up. Currently no sign of this.
3. **Switch to TOML**: would require rewriting all example configs and documentation. Not practical for a project that advertises YAML.

## Recommendation

Option 1. `serde_yaml` 0.9 is the standard choice. The migration is a find-and-replace of import paths.
