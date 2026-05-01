#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

selection_output_dir <- if (length(args) >= 1L && nzchar(args[[1]])) args[[1]] else file.path(repo_root, "outputs", "selection")
barrier_outputs_root <- if (length(args) >= 2L && nzchar(args[[2]])) args[[2]] else file.path(repo_root, "outputs", "barrier_analysis_outputs")

selection_qc_dir <- file.path(selection_output_dir, "QC")
barrier_qc_dir <- file.path(barrier_outputs_root, "QC")
final_figures_dir <- file.path(selection_output_dir, "final_figures")
dir.create(selection_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(barrier_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_figures_dir, recursive = TRUE, showWarnings = FALSE)

go_dir <- file.path(barrier_outputs_root, "go_gene_map_outputs")

read_if_exists <- function(path) if (file.exists(path)) as.data.table(fread(path)) else data.table()
write_both <- function(dt, filename) {
  fwrite(dt, file.path(selection_qc_dir, filename), sep = "\t")
  fwrite(dt, file.path(barrier_qc_dir, filename), sep = "\t")
}
ordered_chr_levels <- function(x) {
  x <- unique(as.character(x))
  nn <- suppressWarnings(as.integer(sub("^Chromosome_([0-9]+)$", "\\1", x)))
  out <- x[order(ifelse(is.na(nn), Inf, nn), x == "Chromosome_ZW", x)]
  if ("Chromosome_ZW" %in% out) out <- c(setdiff(out, "Chromosome_ZW"), "Chromosome_ZW")
  out
}
copy_if_exists <- function(src, dst) {
  if (!file.exists(src)) return(FALSE)
  file.copy(src, dst, overwrite = TRUE)
}
path_slug <- function(path) {
  gsub("[^A-Za-z0-9]+", "_", sub("^_+|_+$", "", path))
}

barrier_models <- read_if_exists(file.path(barrier_outputs_root, "revised_barrier_models_by_tier.tsv"))
bgs_by_tier <- read_if_exists(file.path(barrier_outputs_root, "revised_BGS_by_tier.tsv"))
gene_counts <- read_if_exists(file.path(go_dir, "selected_gene_counts_by_tier.tsv"))
go_counts <- read_if_exists(file.path(go_dir, "GO_term_counts_by_tier.tsv"))
go_results <- read_if_exists(file.path(go_dir, "GO_gprofiler_results_by_tier.tsv"))
immune_fisher <- read_if_exists(file.path(go_dir, "immune_gene_enrichment_fisher_by_tier.tsv"))
immune_window <- read_if_exists(file.path(go_dir, "immune_window_class_comparisons.tsv"))
immune_terms <- read_if_exists(file.path(go_dir, "immune_related_GO_terms_by_tier.tsv"))
selected_windows <- read_if_exists(file.path(go_dir, "selected_windows_by_tier.tsv"))
gene_overlap_null <- read_if_exists(file.path(go_dir, "gene_overlap_matched_null_by_tier.tsv"))
go_null_summary <- read_if_exists(file.path(go_dir, "GO_matched_null_summary.tsv"))
immune_null <- read_if_exists(file.path(go_dir, "immune_matched_null_by_tier.tsv"))
global_chr_summary <- read_if_exists(file.path(selection_qc_dir, "QC_global_chromosome_summary.tsv"))
zw_vs_autosomes <- read_if_exists(file.path(selection_qc_dir, "QC_ZW_vs_autosomes.tsv"))
local_outlier_by_chr <- read_if_exists(file.path(selection_qc_dir, "QC_local_outlier_by_chromosome.tsv"))
window_fdr <- read_if_exists(file.path(selection_qc_dir, "QC_empirical_FDR_selected_windows.tsv"))

write_both(barrier_models, "QC_barrier_models.tsv")
write_both(bgs_by_tier, "QC_BGS_by_tier.tsv")

if (nrow(gene_counts)) {
  gene_qc <- copy(gene_counts)
  gene_qc[, genes_per_overlap_row := unique_gene_ids / pmax(overlapping_cds_rows, 1)]
  if (nrow(gene_overlap_null)) {
    gene_qc <- merge(gene_qc, gene_overlap_null, by = c("tier", "class"), all.x = TRUE)
  }
  write_both(gene_qc, "QC_gene_overlap_by_tier.tsv")
}

if (nrow(go_counts)) {
  broad_terms <- c("cytoplasm", "protein binding", "developmental process", "anatomical structure development", "multicellular organism development")
  go_summary <- merge(
    go_counts,
    gene_counts[, .(tier, class, unique_gene_ids, unique_named_genes)],
    by = c("tier", "class"),
    all.x = TRUE
  )
  if (nrow(go_results)) {
    broad_flag <- go_results[
      is.finite(p_value),
      .(
        n_sig_terms = sum(p_value <= 0.05, na.rm = TRUE),
        n_broad_generic_terms = sum(tolower(term_name) %chin% broad_terms & p_value <= 0.05, na.rm = TRUE),
        min_p_value = min(p_value, na.rm = TRUE)
      ),
      by = .(tier, class, subset)
    ]
    go_summary <- merge(go_summary, broad_flag, by = c("tier", "class", "subset"), all.x = TRUE)
  }
  if (nrow(go_null_summary)) {
    go_summary <- merge(go_summary, go_null_summary, by = c("tier", "class"), all.x = TRUE)
  }
  write_both(go_summary, "QC_GO_summary.tsv")
}

if (nrow(immune_fisher) || nrow(immune_window) || nrow(immune_terms)) {
  immune_summary <- merge(
    immune_fisher,
    immune_window,
    by = "tier",
    all = TRUE,
    suffixes = c("_gene", "_window")
  )
  if (nrow(immune_terms)) {
    immune_term_counts <- immune_terms[is.finite(p_value), .(n_immune_terms = .N, best_immune_term_p = min(p_value, na.rm = TRUE)), by = .(tier, class)]
    immune_summary <- merge(immune_summary, immune_term_counts, by = c("tier", "class"), all = TRUE)
  }
  if (nrow(immune_null)) {
    immune_summary <- merge(immune_summary, immune_null, by = c("tier", "class"), all = TRUE)
  }
  write_both(immune_summary, "QC_immune_summary.tsv")
}

if (nrow(selected_windows)) {
  selected_windows[, `:=`(
    Chromosome_Names = as.character(Chromosome_Names),
    window_start = as.integer(window_start),
    window_end = as.integer(window_end)
  )]
  spatial_qc <- selected_windows[
    selection_class %chin% c("Geographic sweep", "Balancing selection"),
    .(
      n_windows = .N,
      n_chromosomes = uniqueN(Chromosome_Names),
      mean_start = mean(window_start, na.rm = TRUE),
      mean_end = mean(window_end, na.rm = TRUE)
    ),
    by = .(tier, selection_class, Chromosome_Names)
  ]
  write_both(spatial_qc, "QC_chromosome_spatial_summary.tsv")
}

if (nrow(barrier_models)) {
  long_r2 <- melt(
    barrier_models[!is.na(tier)],
    id.vars = "tier",
    measure.vars = intersect(c("abs_delta_p_R2", "f3_R2"), names(barrier_models)),
    variable.name = "response",
    value.name = "R2"
  )
  p_r2 <- ggplot(long_r2[is.finite(R2)], aes(tier, R2, fill = response)) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(x = NULL, y = expression(R^2), fill = "Response", title = "Barrier-model support by tier")
  ggsave(file.path(selection_qc_dir, "barrier_model_R2_by_tier.pdf"), p_r2, width = 9, height = 5.5)
  ggsave(file.path(barrier_qc_dir, "barrier_model_R2_by_tier.pdf"), p_r2, width = 9, height = 5.5)
}

if (nrow(gene_counts)) {
  p_gene <- ggplot(gene_counts, aes(tier, unique_gene_ids, fill = class)) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(x = NULL, y = "Unique genes", fill = "Class", title = "Gene overlap by tier")
  ggsave(file.path(selection_qc_dir, "gene_overlap_by_tier.pdf"), p_gene, width = 9, height = 6)
  ggsave(file.path(barrier_qc_dir, "gene_overlap_by_tier.pdf"), p_gene, width = 9, height = 6)
}

if (nrow(go_counts)) {
  p_go <- ggplot(go_counts[subset == "whole_genome"], aes(tier, n_go_terms, fill = class)) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(x = NULL, y = "Significant GO terms", fill = "Class", title = "GO-term summary by tier")
  ggsave(file.path(selection_qc_dir, "GO_term_summary.pdf"), p_go, width = 9, height = 6)
  ggsave(file.path(barrier_qc_dir, "GO_term_summary.pdf"), p_go, width = 9, height = 6)
}

if (nrow(immune_window)) {
  p_immune <- ggplot(immune_window, aes(tier, odds_ratio, color = interaction(class_a, class_b))) +
    geom_point(size = 2) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(x = NULL, y = "Odds ratio", color = "Comparison", title = "Immune enrichment by tier")
  ggsave(file.path(selection_qc_dir, "immune_enrichment_by_tier.pdf"), p_immune, width = 9, height = 5.5)
  ggsave(file.path(barrier_qc_dir, "immune_enrichment_by_tier.pdf"), p_immune, width = 9, height = 5.5)
}

if (nrow(global_chr_summary) && nrow(local_outlier_by_chr) && nrow(window_fdr)) {
  chr_levels <- ordered_chr_levels(global_chr_summary$Chromosome_Names)
  global_chr_summary[, Chromosome_Names := factor(as.character(Chromosome_Names), levels = chr_levels)]
  local_outlier_by_chr[, Chromosome_Names := factor(as.character(Chromosome_Names), levels = chr_levels)]
  window_fdr[, Chromosome_Names := factor(as.character(Chromosome_Names), levels = chr_levels)]

  auto_mean_sweep <- global_chr_summary[Chromosome_Names != "Chromosome_ZW", mean(mean_sweep_strength, na.rm = TRUE)]
  p_a <- ggplot(global_chr_summary, aes(x = Chromosome_Names, y = mean_sweep_strength)) +
    geom_hline(yintercept = auto_mean_sweep, linetype = 2, color = "grey35", linewidth = 0.5) +
    geom_col(aes(fill = Chromosome_Names == "Chromosome_ZW"), color = "grey20", linewidth = 0.2) +
    scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "grey75"), guide = "none") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(
      x = NULL,
      y = "Mean sweep_strength",
      title = "A. Global chromosome-level shift"
    )

  tail_dt <- copy(window_fdr[is.finite(sweep_strength), .(
    group = fifelse(as.character(Chromosome_Names) == "Chromosome_ZW", "ZW", "Autosomes"),
    sweep_strength = as.numeric(sweep_strength)
  )])
  threshold_grid <- sort(unique(as.numeric(quantile(tail_dt$sweep_strength, probs = seq(0.5, 0.999, length.out = 220), na.rm = TRUE))))
  if (!length(threshold_grid)) threshold_grid <- sort(unique(tail_dt$sweep_strength))
  tail_curve <- rbindlist(lapply(c("Autosomes", "ZW"), function(gr) {
    vals <- tail_dt[group == gr, sweep_strength]
    if (!length(vals)) return(NULL)
    data.table(
      group = gr,
      threshold = threshold_grid,
      prop_above = vapply(threshold_grid, function(th) mean(vals >= th, na.rm = TRUE), numeric(1))
    )
  }), use.names = TRUE, fill = TRUE)
  p_b <- ggplot(tail_curve, aes(x = threshold, y = prop_above, color = group)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c("Autosomes" = "grey35", "ZW" = "#b2182b")) +
    theme_bw(base_size = 10) +
    labs(
      x = "Sweep-strength threshold",
      y = "Proportion of windows above threshold",
      color = NULL,
      title = "B. Global shift vs extreme tail"
    )

  c_dt <- melt(
    local_outlier_by_chr[, .(Chromosome_Names, sweep_fdr_genomewide, balancing_fdr_genomewide)],
    id.vars = "Chromosome_Names",
    variable.name = "class",
    value.name = "n_hits"
  )
  c_dt[, class := fifelse(class == "sweep_fdr_genomewide", "Sweep", "Balancing")]
  p_c <- ggplot(c_dt, aes(x = Chromosome_Names, y = n_hits, fill = class)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.62, color = "grey20", linewidth = 0.15) +
    scale_fill_manual(values = c("Sweep" = "#b2182b", "Balancing" = "#2166ac")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = NULL,
      y = "Genome-wide FDR hits",
      fill = NULL,
      title = "C. Genome-wide FDR <= 0.05 by chromosome"
    )

  d_dt <- melt(
    local_outlier_by_chr[, .(Chromosome_Names, sweep_fdr_chr_aware, balancing_fdr_chr_aware)],
    id.vars = "Chromosome_Names",
    variable.name = "class",
    value.name = "n_hits"
  )
  d_dt[, class := fifelse(class == "sweep_fdr_chr_aware", "Sweep", "Balancing")]
  p_d <- ggplot(d_dt, aes(x = Chromosome_Names, y = n_hits, fill = class)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.62, color = "grey20", linewidth = 0.15) +
    scale_fill_manual(values = c("Sweep" = "#b2182b", "Balancing" = "#2166ac")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = NULL,
      y = "Chromosome-aware FDR hits",
      fill = NULL,
      title = "D. Local chromosome-aware FDR <= 0.05"
    )

  fig05_pdf <- file.path(selection_qc_dir, "ZW_global_shift_vs_extreme_tail.pdf")
  fig05_png <- file.path(selection_qc_dir, "Fig05_ZW_global_shift_vs_extreme_tail.png")
  fig05_caption <- file.path(selection_qc_dir, "Fig05_ZW_global_shift_vs_extreme_tail_caption.txt")
  cap_text <- paste(
    "ZW shows a broad chromosome-wide shift in population-genetic profile,",
    "but contributes few windows to the genomewide extreme tail.",
    "Genomewide FDR is therefore dominated by highly localized autosomal outliers,",
    "whereas ZW signal is more diffuse and chromosome-wide."
  )
  writeLines(cap_text, fig05_caption)

  render_zw_figure <- function(device_fun, filename, width, height, res = NULL) {
    if (is.null(res)) {
      device_fun(filename, width = width, height = height)
    } else {
      device_fun(filename, width = width, height = height, units = "in", res = res)
    }
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2)))
    print(p_a, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p_b, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(p_c, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(p_d, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
    dev.off()
  }
  render_zw_figure(pdf, fig05_pdf, width = 13, height = 10)
  render_zw_figure(png, fig05_png, width = 13, height = 10, res = 300)
  file.copy(fig05_pdf, file.path(barrier_qc_dir, basename(fig05_pdf)), overwrite = TRUE)
  file.copy(fig05_png, file.path(barrier_qc_dir, basename(fig05_png)), overwrite = TRUE)
  file.copy(fig05_caption, file.path(barrier_qc_dir, basename(fig05_caption)), overwrite = TRUE)
}

old_fig_files <- list.files(final_figures_dir, full.names = TRUE, recursive = FALSE)
if (length(old_fig_files)) unlink(old_fig_files, recursive = TRUE, force = TRUE)

all_pdf_sources <- sort(unique(c(
  list.files(selection_output_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE),
  list.files(barrier_outputs_root, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
)))
exclude_patterns <- c(
  "/final_figures/",
  "/go_gene_map_outputs_revised_fdr_qc_",
  "/QC/go_gene_map_outputs_revised_fdr_qc_"
)
keep_pdf <- !Reduce(`|`, lapply(exclude_patterns, function(pat) grepl(pat, all_pdf_sources, fixed = TRUE)))
all_pdf_sources <- all_pdf_sources[keep_pdf]
if (length(all_pdf_sources)) {
  pdf_hash <- unname(tools::md5sum(all_pdf_sources))
  keep_unique <- !duplicated(pdf_hash)
  all_pdf_sources <- all_pdf_sources[keep_unique]
}

source_roots <- c(selection_output_dir, barrier_outputs_root)
source_relative <- vapply(all_pdf_sources, function(src) {
  rel <- src
  for (root in source_roots) {
    prefix <- paste0(normalizePath(root, winslash = "/", mustWork = FALSE), "/")
    nsrc <- normalizePath(src, winslash = "/", mustWork = FALSE)
    if (startsWith(nsrc, prefix)) {
      rel <- substring(nsrc, nchar(prefix) + 1L)
      break
    }
  }
  rel
}, character(1))

figure_specs <- data.table(
  figure_id = sprintf("Fig%03d", seq_along(all_pdf_sources)),
  source_path = all_pdf_sources,
  source_relative = source_relative
)
figure_specs[, filename := sprintf("%s_%s.pdf", figure_id, path_slug(sub("\\.pdf$", "", source_relative, ignore.case = TRUE)))]
figure_specs[, source_script := fifelse(
  grepl("go_gene_map_outputs|run_connected_go_gene_map_workflow", source_relative),
  "run_connected_go_gene_map_workflow.R",
  fifelse(
    grepl("barrier_recalibrated_outputs|barrier_model|plot_", source_relative),
    "barrier_recalibration_followup_sd3_new_af.R",
    fifelse(
      grepl("^QC/|selection_qc_and_fdr_helper", source_relative),
      "selection_qc_and_fdr_helper.R",
      "workflow_output"
    )
  )
)]
figure_specs[, interpretation := fifelse(
  grepl("ZW_global_shift_vs_extreme_tail", source_relative),
  "ZW global shift versus extreme-tail contrast",
  fifelse(
    grepl("immune", source_relative, ignore.case = TRUE),
    "Immune-focused summary",
    fifelse(
      grepl("GO|gprofiler", source_relative),
      "GO summary",
      fifelse(
        grepl("chromosome", source_relative, ignore.case = TRUE),
        "Chromosome-level figure",
        "Workflow figure"
      )
    )
  )
)]
figure_specs[, main_or_supplement := "supplement"]
figure_specs[grepl("GMM_BIC_silhouette_ARI|tier_population_genetic_profiles|chromosome_selected_window_tracks|plot_B_no_delta_p_vs_abs_delta_p|ZW_global_shift_vs_extreme_tail|immune_window_odds_ratio_comparisons|gene_overlap_by_tier|GO_term_summary", source_relative), main_or_supplement := "main"]
figure_specs[, status := vapply(seq_len(.N), function(i) {
  src <- source_path[[i]]
  dst <- file.path(final_figures_dir, filename[[i]])
  if (copy_if_exists(src, dst)) "copied" else "missing"
}, character(1))]
fwrite(figure_specs, file.path(final_figures_dir, "figure_manifest.tsv"), sep = "\t")
fig05_png_src <- file.path(selection_qc_dir, "Fig05_ZW_global_shift_vs_extreme_tail.png")
fig05_txt_src <- file.path(selection_qc_dir, "Fig05_ZW_global_shift_vs_extreme_tail_caption.txt")
copy_if_exists(fig05_png_src, file.path(final_figures_dir, "Fig005_ZW_global_shift_vs_extreme_tail.png"))
copy_if_exists(fig05_txt_src, file.path(final_figures_dir, "Fig005_ZW_global_shift_vs_extreme_tail_caption.txt"))

message("Connected QC summary workflow complete.")
