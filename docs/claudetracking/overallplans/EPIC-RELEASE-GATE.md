# EPIC-RELEASE-GATE: Release Blockers

## Goal

Address all issues that must be resolved before v0.1.0 can ship. This includes
security advisories in dependencies, silent config misconfigurations, and the
primary quality gate (downstream variant caller validation).

## Tasks

- T133: Migrate off serde_yml (security advisory)
- T134: Reject unsupported inline UMI config
- T135: Run downstream variant caller validation (Mutect2)
- T136: Validate stochastic VAF with chi-squared test

## Priority

All tasks are high priority. The dependency migration (T133) and config guard
(T134) are quick wins. The downstream validation (T135) is the primary release
gate per the vision docs.
