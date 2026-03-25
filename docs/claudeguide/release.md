# Release Configuration

## Registry

- **Primary**: crates.io (`cargo publish`)
- **Secondary**: GitHub Releases (tag + release notes)

## Publish command

```bash
cargo publish
```

## Auth

Requires `cargo login` with a crates.io API token before publishing.

## Pre-publish checks

```bash
cargo fmt -- --check
cargo clippy -- -D warnings
cargo test
```

## Post-publish steps

1. Create a GitHub release from the tag with release notes.
2. Release notes should summarise features, fixes, and breaking changes.

## Package exclusions

Managed in `Cargo.toml` `[package] exclude`. Non-package files (docs, paper,
benchmarks, research, scripts, CI config) must be excluded.

## Version scheme

`v0.Y.Z` per CLAUDE.md. Z bump for most releases. Y bump for breaking changes.
