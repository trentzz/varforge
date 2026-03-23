# EPIC-CLI-ERGONOMICS: CLI Usability

## Goal

Make common workflows easier to discover and execute from the command line.

## Motivation

Users cannot discover available presets without reading source code. The validate subcommand
only checks config structure, not file existence. Misspelled preset names give a terse error
with no suggestions.

## Scope

- In scope: --list-presets, enhanced validate, fuzzy preset name matching.
- Out of scope: interactive mode, shell completion scripts.

## Tasks

- T095: Add --list-presets flag to simulate subcommand (lists base and cancer presets)
- T096: Extend validate subcommand to check reference readability, BED validity, VCF parse
- T097: Suggest closest preset name when an unknown preset is specified
