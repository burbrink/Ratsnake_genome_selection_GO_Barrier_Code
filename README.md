# Ratsnake Unified Workflow

This repository combines the current ratsnake selection-calling, barrier, gene/GO, immune, QC, and figure-numbering workflow into one modular codebase.

## Layout

- `config/`
  - local path defaults and runtime settings
- `scripts/01_selection/`
  - master selection pipeline
  - explicit-neutral continuous-distance QC/FDR helper
- `scripts/02_barrier/`
  - connected barrier, GO/gene-map, decision-table, QC, and figure-catalog scripts
- `outputs/selection/`
  - selection pipeline outputs
- `outputs/barrier_analysis_outputs/`
  - downstream barrier, GO, immune, QC, and figure outputs

## Main driver

Run the full workflow from the repository root:

```r
Rscript run_full_workflow.R
```

Supported modes:

```r
Rscript run_full_workflow.R selection
Rscript run_full_workflow.R downstream
Rscript run_full_workflow.R full
```

Fresh rerun from an empty output tree:

```r
PIPELINE_CLEAN_OUTPUTS=1 Rscript run_full_workflow.R full
```

## Current defaults

The default input `.rdata`, GTF, output locations, and object name are defined in:

- `config/default_paths.R`

Update that file if the raw data or annotation locations move.

## Important interpretation hierarchy

The workflow now keeps these layers separate:

1. Broad genome-state map
   - `GMM/RF + continuous-distance agreement`
2. Global chromosome-level signal
   - chromosome summaries and genome-wide contrasts
3. Local significant windows
   - chromosome-aware empirical `FDR <= 0.05`
4. Extreme candidate loci
   - `distance-margin`, `Top 1%`, and `Top 500`

This is especially important for `ZW`, which can be broadly shifted toward selected-like states while contributing few windows to the whole-genome extreme tail.

## Numbered figures

The unified figure catalog is written to:

- `outputs/selection/final_figures/figure_manifest.tsv`

and the numbered PDFs are copied into:

- `outputs/selection/final_figures/`
