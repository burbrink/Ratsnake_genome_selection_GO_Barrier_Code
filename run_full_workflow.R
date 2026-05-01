#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tools)
})

script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(script_dir, mustWork = FALSE)

source(file.path(repo_root, "config", "default_paths.R"))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L && nzchar(args[[1]])) args[[1]] else "full"
clean_outputs <- identical(Sys.getenv("PIPELINE_CLEAN_OUTPUTS", unset = "0"), "1")

selection_output_dir <- normalizePath(file.path(repo_root, workflow_defaults$selection_output_dir), mustWork = FALSE)
barrier_output_dir <- normalizePath(file.path(repo_root, workflow_defaults$barrier_output_dir), mustWork = FALSE)

selection_script <- file.path(repo_root, "scripts", "01_selection", "master_pipeline_UPDATED_v2026-02-07-02_ranknormRec_QCfix_Fig17.R")
selection_qc_script <- file.path(repo_root, "scripts", "01_selection", "selection_qc_and_fdr_helper.R")
downstream_script <- file.path(repo_root, "scripts", "02_barrier", "run_connected_selection_barrier_workflow.R")
qc_summary_script <- file.path(repo_root, "scripts", "02_barrier", "run_connected_qc_summary_workflow.R")

run_step <- function(cmd, args, workdir = repo_root, env = character()) {
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(workdir)
  status <- system2(cmd, args = args, env = env, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("Step failed: ", cmd, " ", paste(args, collapse = " "), call. = FALSE)
  }
}

message("Unified workflow repo: ", repo_root)
message("Mode: ", mode)
message("Selection outputs: ", selection_output_dir)
message("Barrier outputs: ", barrier_output_dir)
message("Clean outputs: ", clean_outputs)

if (clean_outputs) {
  if (dir.exists(selection_output_dir)) unlink(selection_output_dir, recursive = TRUE, force = TRUE)
  if (dir.exists(barrier_output_dir)) unlink(barrier_output_dir, recursive = TRUE, force = TRUE)
}

dir.create(selection_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(barrier_output_dir, recursive = TRUE, showWarnings = FALSE)

if (mode %in% c("selection", "full")) {
  message("\n=== Running selection master pipeline ===")
  run_step(
    "Rscript",
    c(selection_script),
    workdir = selection_output_dir,
    env = c(
      paste0("PIPELINE_INPUT_RDATA=", workflow_defaults$input_rdata),
      paste0("PIPELINE_N_CORES=", workflow_defaults$n_cores)
    )
  )

  message("\n=== Running selection QC/FDR helper ===")
  run_step(
    "Rscript",
    c(selection_qc_script, selection_output_dir),
    workdir = selection_output_dir
  )
}

if (mode %in% c("downstream", "full")) {
  message("\n=== Running connected downstream workflow ===")
  run_step(
    "Rscript",
    c(
      downstream_script,
      selection_output_dir,
      workflow_defaults$input_rdata,
      workflow_defaults$object_name,
      barrier_output_dir
    ),
    workdir = repo_root,
    env = c(paste0("PIPELINE_GTF_FILE=", workflow_defaults$gtf_file))
  )

  message("\n=== Refreshing numbered figure catalog ===")
  run_step(
    "Rscript",
    c(qc_summary_script, selection_output_dir, barrier_output_dir),
    workdir = repo_root
  )
}

message("\nUnified workflow complete.")
