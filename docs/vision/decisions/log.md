# VarForge: Decision Log

Append-only. Add new entries at the bottom. Do not edit past entries.

---

## 2026-03-24: Pre-release direction set

Decision: VarForge is feature-complete for v1. Focus shifts to stability, testing, and downstream tool compatibility. The paper releases alongside the tool.

Context: 106 tasks completed. The core simulation pipeline works end-to-end. Three in-progress features (BAM editor T066, ClonalTree integration T067, cancer presets T068) must land. Critical bugs found in UMI BAM output, truth VCF sorting, and a panic on zero standard deviation.

Key choices:
- ContaminationConfig and FragmentModel::Custom deferred to post-v1. Remove from the v1 config schema.
- Inline UMI mode deferred. Document as unsupported in v1.
- Chromosome name remapping deferred. Document the requirement for consistent naming as a user responsibility.
- Downstream tool validation (Mutect2, fgbio, hap.py) is the primary quality gate for release. No release until these experiments pass.
