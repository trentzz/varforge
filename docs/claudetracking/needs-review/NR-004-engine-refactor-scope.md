# Engine refactor: scope and timing

**Category**: architecture
**Related epic**: EPIC-POLISH

## Context

`generate_region` in `engine.rs` is ~500 lines. It handles fragment sampling,
variant spike-in (small variants and SVs), quality generation, UMI family
expansion, artifact injection, and read pair construction. This makes it hard
to test individual stages in isolation and difficult to modify one stage
without risk to others.

## Options

1. **Full decomposition now (T140)**: extract 4-5 helper methods before release.
   Lower risk of regressions if done carefully, but adds work before release.

2. **Defer to post-release**: ship v0.1 as-is, refactor for v0.2. The code
   works and passes all tests. The risk is that post-release features (clonal
   evolution, contamination) will further enlarge the function.

3. **Minimal extraction**: extract only the SV spike-in block (~80 lines) and
   UMI expansion block (~120 lines), which are the most self-contained. Leaves
   the core fragment loop intact.

## Recommendation

Option 3 for now. Extract the two most independent blocks (SV and UMI) to
reduce the function to ~300 lines. Defer the full decomposition to post-release.
