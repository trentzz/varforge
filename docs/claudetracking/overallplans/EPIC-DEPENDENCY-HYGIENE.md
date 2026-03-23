# EPIC-DEPENDENCY-HYGIENE: Dependency Cleanup

## Goal

Remove unused dependencies and migrate off unmaintained crates to reduce compile time and
supply-chain risk.

## Motivation

serde_yaml 0.9 is unmaintained and its successor (serde_yml) is a drop-in replacement.
Several noodles sub-crates may be pulling in unused features. Bloated Cargo.toml increases
cold compile times and widens the attack surface.

## Scope

- In scope: unused crate removal, serde_yaml migration, unused feature flags.
- Out of scope: upgrading all crates to latest (too risky without a test suite running CI).

## Tasks

- T093: Audit and remove unused dependencies (run cargo-machete or manual audit)
- T094: Migrate from serde_yaml 0.9 to serde_yml
