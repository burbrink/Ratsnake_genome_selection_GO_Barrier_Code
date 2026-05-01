#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
source(file.path(script_dir, "analysis_sanity_checks.R"))

cfg <- list(
  input_rdata = Sys.getenv("PIPELINE_INPUT_RDATA", unset = file.path(repo_root, "data", "input_placeholder.rdata")),
  object_name = "sd3_new_af",
  prior_barrier_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_gmm_rf"),
  prior_hybrid_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "hybrid_zone_followup_outputs"),
  prior_driver_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_driver_outputs"),
  prior_domain_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_domain_outputs"),
  prior_moran_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_moransI_outputs"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "gene_density_followup_outputs"),
  primary_high_quantile = 0.95
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_rdata <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$object_name <- args[[2]]
if (length(args) >= 3 && nzchar(args[[3]])) cfg$output_dir <- args[[3]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
sanity <- init_sanity("gene_density_barrier_followup_sd3_new_af.R", cfg$output_dir)

# Biological framing:
# Gene density matters because barrier accumulation could be stronger in
# functionally dense regions with more potential targets of divergent selection.
# A positive association would support target-rich selection effects, whereas a
# negative association would suggest recombination deserts or structural effects
# dominate. Weak or mixed results would imply barrier domains are driven more by
# recombination and sweep history than by gene density alone.

find_col <- function(df, candidates = NULL, regex = NULL) {
  nms <- names(df)
  if (!is.null(candidates)) {
    hit <- intersect(candidates, nms)
    if (length(hit)) return(hit[[1]])
  }
  if (!is.null(regex)) {
    hit <- nms[grepl(regex, nms, ignore.case = TRUE, perl = TRUE)]
    if (length(hit)) return(hit[[1]])
  }
  NA_character_
}

find_all_cols <- function(df, regex) {
  names(df)[grepl(regex, names(df), ignore.case = TRUE, perl = TRUE)]
}

order_chromosomes <- function(chr) {
  tibble(chromosome = unique(as.character(chr))) %>%
    mutate(
      is_zw = grepl("ZW", chromosome, ignore.case = TRUE),
      chr_num = suppressWarnings(as.integer(gsub("^Chromosome_", "", chromosome))),
      chr_sort = dplyr::case_when(
        is_zw ~ 1e9,
        is.finite(chr_num) ~ chr_num,
        TRUE ~ 1e8 + rank(chromosome, ties.method = "first")
      )
    ) %>%
    arrange(chr_sort, chromosome) %>%
    pull(chromosome)
}

load(cfg$input_rdata)
assert_true(exists(cfg$object_name, inherits = FALSE), paste("Object not found:", cfg$object_name), sanity, "object_exists")
dat_raw <- get(cfg$object_name, inherits = FALSE)
assert_true(nrow(as.data.frame(dat_raw)) > 0L, "Loaded dataset is empty.", sanity, "nonempty_input_dataset")

chr_col <- find_col(dat_raw, c("Chromosome_Names", "Chromosome", "chr", "chromosome"))
start_col <- find_col(dat_raw, c(".start_bp", "window_start", "start_bp", "start"))
end_col <- find_col(dat_raw, c(".end_bp", "window_end", "end_bp", "end"))
mid_col <- find_col(dat_raw, c("midpoint", ".mid_bp", "window_mid", "mid_bp", "mid"))
recomb_col <- find_col(dat_raw, c("recombination", "rec", "rho"))
fst_col <- find_col(dat_raw, c("avg_wc_fst", "wc_fst", "fst"))
dxy_col <- find_col(dat_raw, c("avg_dxy", "dxy"))
pi_col <- find_col(dat_raw, c("pi_average", "avg_pi", "pi"))
taj_col <- find_col(dat_raw, c("TajimaD", "tajima_d"))
delta_col <- find_col(dat_raw, c("mean_abs_delta_p", "abs_delta_p"))
f3_col <- find_col(dat_raw, c("f3", "F3"))
snp_col <- find_col(dat_raw, c("n_snps_af", "N_SNPS", "no_snps", "n_variants", "snp_density"))
label_col <- find_col(dat_raw, c("sweep_group", "cont_label_h", "classification", "label"))

gene_candidates <- find_all_cols(
  dat_raw,
  "(^gene_density$|genes?_per_window|genes?_per_10kb|gene_count|n_genes|cds_density|exon_density|gene_.*density|density_gene|genecontent)"
)

detected_columns <- tibble(
  role = c("chromosome", "start", "end", "mid", "recombination", "avg_wc_fst", "avg_dxy", "pi_average", "TajimaD", "abs_delta_p", "f3", "snp_density", "label"),
  detected_column = c(chr_col, start_col, end_col, mid_col, recomb_col, fst_col, dxy_col, pi_col, taj_col, delta_col, f3_col, snp_col, label_col)
)
if (length(gene_candidates)) {
  detected_columns <- bind_rows(
    detected_columns,
    tibble(role = paste0("gene_density_candidate_", seq_along(gene_candidates)), detected_column = gene_candidates)
  )
}

dat <- as.data.frame(dat_raw) %>%
  mutate(
    Chromosome_Names = as.character(.data[[chr_col]]),
    .start_bp = as.numeric(.data[[start_col]]),
    .end_bp = if (!is.na(end_col)) as.numeric(.data[[end_col]]) else .start_bp + 9999,
    window_mid_bp = if (!is.na(mid_col)) as.numeric(.data[[mid_col]]) else (.start_bp + .end_bp) / 2,
    recombination = if (!is.na(recomb_col)) as.numeric(.data[[recomb_col]]) else NA_real_,
    avg_wc_fst = if (!is.na(fst_col)) as.numeric(.data[[fst_col]]) else NA_real_,
    avg_dxy = if (!is.na(dxy_col)) as.numeric(.data[[dxy_col]]) else NA_real_,
    pi_average = if (!is.na(pi_col)) as.numeric(.data[[pi_col]]) else NA_real_,
    TajimaD = if (!is.na(taj_col)) as.numeric(.data[[taj_col]]) else NA_real_,
    abs_delta_p = if (!is.na(delta_col)) as.numeric(.data[[delta_col]]) else NA_real_,
    f3 = if (!is.na(f3_col)) as.numeric(.data[[f3_col]]) else NA_real_,
    snp_density = if (!is.na(snp_col)) as.numeric(.data[[snp_col]]) else NA_real_,
    sweep_group = if (!is.na(label_col)) as.character(.data[[label_col]]) else NA_character_
  )

assert_has_cols(dat, c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"), sanity, "gene_density_input_data")
assert_unique_keys(dat, c("Chromosome_Names", ".start_bp", ".end_bp"), sanity, "gene_density_input_data")

prior_local <- read.csv(file.path(cfg$prior_barrier_dir, "barrier_recalibrated_local_summary.csv"), stringsAsFactors = FALSE)
n_before_join <- nrow(dat)
dat <- dat %>% left_join(prior_local, by = c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"))
assert_join_rowcount(n_before_join, dat, sanity, "gene_density_prior_local_join")
if ("abs_delta_p.x" %in% names(dat) || "abs_delta_p.y" %in% names(dat)) dat$abs_delta_p <- dplyr::coalesce(dat$abs_delta_p.x, dat$abs_delta_p.y)
if ("f3.x" %in% names(dat) || "f3.y" %in% names(dat)) dat$f3 <- dplyr::coalesce(dat$f3.x, dat$f3.y)
if ("sweep_group.x" %in% names(dat) || "sweep_group.y" %in% names(dat)) dat$sweep_group <- dplyr::coalesce(dat$sweep_group.y, dat$sweep_group.x)

domain_table <- NULL
domain_table_path <- file.path(cfg$prior_domain_dir, "k_local_domain_table.csv")
if (file.exists(domain_table_path)) domain_table <- read.csv(domain_table_path, stringsAsFactors = FALSE)

domain_summary <- NULL
domain_summary_path <- file.path(cfg$prior_domain_dir, "k_local_domain_summary.csv")
if (file.exists(domain_summary_path)) domain_summary <- read.csv(domain_summary_path, stringsAsFactors = FALSE)

gene_project_candidates <- c(
  "ratsnake_GFF/gene_eggnog_by_windows.txt",
  "ratsnake_GFF/intersected_genes.txt",
  "ratsnake_GFF/intersected_genes.gtf",
  "ratsnake_GFF/braker_passqc_eggnog_genes.gtf",
  "ratsnake_GFF/windows_10kb.bed"
)

project_file_tbl <- tibble(
  file = gene_project_candidates,
  exists = file.exists(gene_project_candidates),
  size_bytes = ifelse(file.exists(gene_project_candidates), file.info(gene_project_candidates)$size, NA_real_)
)

selected_gene_col <- NA_character_
if (length(gene_candidates)) {
  finite_counts <- vapply(gene_candidates, function(nm) sum(is.finite(suppressWarnings(as.numeric(dat_raw[[nm]])))), numeric(1))
  ord <- order(-finite_counts, gene_candidates)
  if (length(ord) && finite_counts[ord[1]] > 0) selected_gene_col <- gene_candidates[ord[1]]
}

domain_gene_info <- tibble(
  source = c("k_local_domain_table.mean_gene_density", "k_local_domain_table.mean_snp_density"),
  n_finite = c(
    if (!is.null(domain_table) && "mean_gene_density" %in% names(domain_table)) sum(is.finite(domain_table$mean_gene_density)) else NA_real_,
    if (!is.null(domain_table) && "mean_snp_density" %in% names(domain_table)) sum(is.finite(domain_table$mean_snp_density)) else NA_real_
  )
)

detection_report <- tibble(
  dataset_gene_candidates = if (length(gene_candidates)) paste(gene_candidates, collapse = "; ") else "<none>",
  selected_gene_density_column = ifelse(is.na(selected_gene_col), "<none>", selected_gene_col),
  selected_column_finite_n = ifelse(is.na(selected_gene_col), 0, sum(is.finite(suppressWarnings(as.numeric(dat_raw[[selected_gene_col]]))))),
  domain_table_has_mean_gene_density = !is.null(domain_table) && "mean_gene_density" %in% names(domain_table),
  domain_table_mean_gene_density_finite_n = if (!is.null(domain_table) && "mean_gene_density" %in% names(domain_table)) sum(is.finite(domain_table$mean_gene_density)) else 0
)

write.csv(detected_columns, file.path(cfg$output_dir, "detected_columns.csv"), row.names = FALSE)
write.csv(project_file_tbl, file.path(cfg$output_dir, "gene_density_project_file_check.csv"), row.names = FALSE)
write.csv(domain_gene_info, file.path(cfg$output_dir, "gene_density_domain_field_check.csv"), row.names = FALSE)
write.csv(detection_report, file.path(cfg$output_dir, "gene_density_detection_report.csv"), row.names = FALSE)

if (is.na(selected_gene_col)) {
  message("\nGene-density follow-up could not proceed.")
  message("No usable gene-density or gene-count column was found in sd3_new_af.")
  message("Dataset gene-density candidates: ", if (length(gene_candidates)) paste(gene_candidates, collapse = ", ") else "<none>")
  message("Project-side gene-density candidate files checked:")
  print(project_file_tbl)
  if (!is.null(domain_table) && "mean_gene_density" %in% names(domain_table)) {
    message("Prior domain table contains mean_gene_density but finite values = ", sum(is.finite(domain_table$mean_gene_density)), ".")
  }
  message("Because no real per-window gene-density variable is available in this workspace, the script is stopping gracefully after writing detection reports.")
  finalize_sanity(
    sanity,
    files = c(
      "detected_columns.csv",
      "gene_density_project_file_check.csv",
      "gene_density_domain_field_check.csv",
      "gene_density_detection_report.csv"
    )
  )
  quit(save = "no", status = 0)
}

# If a real gene-density column is ever present, the selected field will be
# mapped here and the downstream analysis can proceed.
dat$gene_density <- suppressWarnings(as.numeric(dat[[selected_gene_col]]))
dat$Chromosome_Names <- factor(dat$Chromosome_Names, levels = order_chromosomes(dat$Chromosome_Names))

message("Gene-density column detected: ", selected_gene_col)
message("A usable gene-density field is now available, so this script can be extended to run the full modeling workflow.")

finalize_sanity(
  sanity,
  files = c(
    "detected_columns.csv",
    "gene_density_project_file_check.csv",
    "gene_density_domain_field_check.csv",
    "gene_density_detection_report.csv"
  )
)
