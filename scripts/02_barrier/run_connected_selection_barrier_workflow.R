#!/usr/bin/env Rscript

script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
selection_output_dir <- file.path(repo_root, "outputs", "selection")
input_rdata <- Sys.getenv("PIPELINE_INPUT_RDATA", unset = "")
object_name <- "sd3_new_af"
outputs_root <- file.path(repo_root, "outputs", "barrier_analysis_outputs")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && nzchar(args[[1]])) selection_output_dir <- args[[1]]
if (length(args) >= 2L && nzchar(args[[2]])) input_rdata <- args[[2]]
if (length(args) >= 3L && nzchar(args[[3]])) object_name <- args[[3]]
if (length(args) >= 4L && nzchar(args[[4]])) outputs_root <- args[[4]]

label_rdata <- file.path(selection_output_dir, "barrier_label_objects.rdata")
handoff_manifest <- file.path(selection_output_dir, "Selection_Barrier_Handoff_Manifest.tsv")

if (!file.exists(label_rdata)) {
  stop(
    "Missing selection-to-barrier handoff label file:\n  ", label_rdata,
    "\nRun the updated selection pipeline first, then rerun this script.",
    call. = FALSE
  )
}

if (!nzchar(input_rdata)) {
  stop("Missing input RData path. Pass it as arg 2 or set PIPELINE_INPUT_RDATA.", call. = FALSE)
}

if (file.exists(handoff_manifest)) {
  message("Using selection handoff manifest: ", handoff_manifest)
} else {
  message("Selection handoff manifest was not found, but label object exists: ", label_rdata)
}

rules <- data.frame(
  rule = c(
    "majority_2of3",
    "gmm_rf",
    "best_continuous_3of3",
    "continuous_distance_3class",
    "gmm_rf_contdist",
    "empirical_fdr_selected",
    "genomewide_fdr_selected",
    "global_shift_selected",
    "dist_margin_025",
    "dist_margin_050",
    "dist_margin_100",
    "score_margin_025",
    "score_margin_050",
    "score_margin_100",
    "dist_top1pct",
    "dist_top500",
    "margin_025",
    "margin_050",
    "margin_100",
    "top1pct",
    "tail500"
  ),
  min_agreement = c(2L, 2L, 3L, 1L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L),
  output_dir = file.path(
    outputs_root,
    c(
      "barrier_recalibrated_outputs_majority_2of3",
      "barrier_recalibrated_outputs_gmm_rf",
      "barrier_recalibrated_outputs_best_continuous_3of3",
      "barrier_recalibrated_outputs_continuous_distance_3class",
      "barrier_recalibrated_outputs_gmm_rf_contdist",
      "barrier_recalibrated_outputs_empirical_fdr_selected",
      "barrier_recalibrated_outputs_genomewide_fdr_selected",
      "barrier_recalibrated_outputs_global_shift_selected",
      "barrier_recalibrated_outputs_dist_margin_025",
      "barrier_recalibrated_outputs_dist_margin_050",
      "barrier_recalibrated_outputs_dist_margin_100",
      "barrier_recalibrated_outputs_score_margin_025",
      "barrier_recalibrated_outputs_score_margin_050",
      "barrier_recalibrated_outputs_score_margin_100",
      "barrier_recalibrated_outputs_dist_top1pct",
      "barrier_recalibrated_outputs_dist_top500",
      "barrier_recalibrated_outputs_margin_025",
      "barrier_recalibrated_outputs_margin_050",
      "barrier_recalibrated_outputs_margin_100",
      "barrier_recalibrated_outputs_top1pct",
      "barrier_recalibrated_outputs_tail500"
    )
  ),
  stringsAsFactors = FALSE
)

recalibration_script <- file.path(script_dir, "barrier_recalibration_followup_sd3_new_af.R")
comparison_script <- file.path(script_dir, "compare_barrier_2of3_vs_triple.R")
go_gene_map_script <- file.path(script_dir, "run_connected_go_gene_map_workflow.R")
decision_table_update_script <- file.path(script_dir, "update_connected_decision_table.R")
qc_summary_script <- file.path(script_dir, "run_connected_qc_summary_workflow.R")

for (i in seq_len(nrow(rules))) {
  message("\n=== Running barrier recalibration rule: ", rules$rule[[i]], " ===")
  dir.create(rules$output_dir[[i]], recursive = TRUE, showWarnings = FALSE)
  args <- c(
    recalibration_script,
    input_rdata,
    object_name,
    rules$output_dir[[i]],
    label_rdata,
    as.character(rules$min_agreement[[i]]),
    "__none__",
    "__none__",
    "__none__",
    rules$rule[[i]]
  )
  status <- system2("Rscript", args = args)
  if (!identical(status, 0L)) {
    stop("Barrier recalibration failed for rule: ", rules$rule[[i]], call. = FALSE)
  }
}

message("\n=== Running agreement-tier comparison plots ===")
status <- system2(
  "Rscript",
  args = c(
    comparison_script,
    label_rdata,
    file.path(outputs_root, "barrier_comparison_outputs")
  )
)
if (!identical(status, 0L)) stop("Agreement-tier comparison failed.", call. = FALSE)

message("\n=== Running connected GO / gene-map workflow ===")
status <- system2(
  "Rscript",
  args = c(
    go_gene_map_script,
    selection_output_dir,
    file.path(outputs_root, "go_gene_map_outputs")
  )
)
if (!identical(status, 0L)) stop("Connected GO / gene-map workflow failed.", call. = FALSE)

message("\n=== Updating consolidated decision table ===")
status <- system2(
  "Rscript",
  args = c(
    decision_table_update_script,
    selection_output_dir,
    outputs_root
  )
)
if (!identical(status, 0L)) stop("Decision-table update failed.", call. = FALSE)

message("\n=== Writing connected QC summary outputs ===")
status <- system2(
  "Rscript",
  args = c(
    qc_summary_script,
    selection_output_dir,
    outputs_root
  )
)
if (!identical(status, 0L)) stop("Connected QC summary workflow failed.", call. = FALSE)

write.csv(
  rules,
  file.path(outputs_root, "connected_selection_barrier_rules.csv"),
  row.names = FALSE
)

message("\nConnected workflow complete. Outputs written under: ", outputs_root)
