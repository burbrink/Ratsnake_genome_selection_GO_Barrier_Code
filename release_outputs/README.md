# Release Outputs

This folder contains a lightweight, GitHub-safe subset of outputs from the unified ratsnake workflow.

It is intended to provide:
- manuscript-facing summary tables
- a compact figure set
- key ZW/global-vs-local interpretation tables

It intentionally excludes:
- full window-level assignment tables
- full selected-window and selected-gene master tables
- large intermediate outputs used only for local reruns

## Layout

- `tables/`
  - summary TSV/CSV outputs for selection tiers, chromosome summaries, barrier comparisons, GO, gene overlap, and immune analyses
- `figures/`
  - a compact numbered PDF figure set
- `metadata/`
  - manifest files describing the release bundle

## Key interpretation hierarchy

1. Broad genome-state map
   - `GMM/RF + continuous-distance agreement`
2. Global chromosome-level signal
   - chromosome summaries and genome-wide contrasts
3. Local significant windows
   - chromosome-aware empirical `FDR <= 0.05`
4. Extreme candidate loci
   - distance-margin tiers, `Top 1%`, and `Top 500`

This distinction is especially important for `ZW`, which remains broadly shifted toward selected-like states while contributing very few windows to the strict genome-wide extreme tail.
