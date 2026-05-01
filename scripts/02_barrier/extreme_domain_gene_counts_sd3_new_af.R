#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})
script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

cfg <- list(
  extreme_domain_csv = file.path(repo_root, "outputs", "barrier_analysis_outputs", "two_scale_barrier_outputs", "two_scale_barrier_extreme_domain_table.csv"),
  gene_window_tsv = file.path(repo_root, "outputs", "selection", "BGS_GeneDensity_Windows.tsv"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "extreme_domain_gene_counts_outputs")
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$extreme_domain_csv <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$gene_window_tsv <- args[[2]]
if (length(args) >= 3 && nzchar(args[[3]])) cfg$output_dir <- args[[3]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cfg$extreme_domain_csv)) stop("Extreme domain CSV not found: ", cfg$extreme_domain_csv, call. = FALSE)
if (!file.exists(cfg$gene_window_tsv)) stop("Gene window TSV not found: ", cfg$gene_window_tsv, call. = FALSE)

dom <- read.csv(cfg$extreme_domain_csv, stringsAsFactors = FALSE, check.names = FALSE)
gw <- read.delim(cfg$gene_window_tsv, stringsAsFactors = FALSE, check.names = FALSE)

required_dom <- c("Chromosome_Names", "threshold_prob", "domain_id", "start_bp", "end_bp", "n_windows", "domain_length_bp")
required_gw <- c("Chromosome_Names", "window_start", "window_end", "gene_n", "gene_density_perMb")
if (length(setdiff(required_dom, names(dom)))) stop("Extreme domain CSV missing required columns.", call. = FALSE)
if (length(setdiff(required_gw, names(gw)))) stop("Gene window TSV missing required columns.", call. = FALSE)

dom <- dom %>%
  mutate(
    Chromosome_Names = as.character(Chromosome_Names),
    threshold_prob = as.numeric(threshold_prob),
    start_bp = as.numeric(start_bp),
    end_bp = as.numeric(end_bp),
    n_windows = as.numeric(n_windows),
    domain_length_bp = as.numeric(domain_length_bp)
  )

gw <- gw %>%
  transmute(
    Chromosome_Names = as.character(Chromosome_Names),
    window_start = as.numeric(window_start),
    window_end = as.numeric(window_end),
    gene_n = as.numeric(gene_n),
    gene_density_perMb = as.numeric(gene_density_perMb)
  )

gw_split <- split(gw, gw$Chromosome_Names)

annotated <- bind_rows(lapply(seq_len(nrow(dom)), function(i) {
  d <- dom[i, , drop = FALSE]
  chr_df <- gw_split[[d$Chromosome_Names]]
  if (is.null(chr_df) || !nrow(chr_df)) {
    d$gene_count_sum <- NA_real_
    d$mean_gene_n_per_window <- NA_real_
    d$mean_gene_density_perMb <- NA_real_
    return(d)
  }
  hits <- chr_df %>% filter(!(window_end < d$start_bp | window_start > d$end_bp))
  if (!nrow(hits)) {
    d$gene_count_sum <- 0
    d$mean_gene_n_per_window <- 0
    d$mean_gene_density_perMb <- 0
    return(d)
  }
  d$gene_count_sum <- sum(hits$gene_n, na.rm = TRUE)
  d$mean_gene_n_per_window <- mean(hits$gene_n, na.rm = TRUE)
  d$mean_gene_density_perMb <- mean(hits$gene_density_perMb, na.rm = TRUE)
  d
}))

summary_by_threshold <- annotated %>%
  group_by(threshold_prob) %>%
  summarise(
    n_domains = n(),
    mean_genes_per_domain = mean(gene_count_sum, na.rm = TRUE),
    median_genes_per_domain = median(gene_count_sum, na.rm = TRUE),
    min_genes_per_domain = min(gene_count_sum, na.rm = TRUE),
    max_genes_per_domain = max(gene_count_sum, na.rm = TRUE),
    total_genes_across_domains = sum(gene_count_sum, na.rm = TRUE),
    total_windows_across_domains = sum(n_windows, na.rm = TRUE),
    mean_genes_per_window = sum(gene_count_sum, na.rm = TRUE) / sum(n_windows, na.rm = TRUE),
    mean_gene_density_perMb = mean(mean_gene_density_perMb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(threshold_prob)

summary_primary_vs_background <- {
  primary <- annotated %>% filter(abs(threshold_prob - 0.99) < 1e-12)
  bg <- gw %>%
    transmute(
      group = "All windows",
      gene_n = gene_n,
      gene_density_perMb = gene_density_perMb
    )
  bind_rows(
    primary %>% transmute(group = "Top 1% domains", gene_n = mean_gene_n_per_window, gene_density_perMb = mean_gene_density_perMb),
    bg
  ) %>%
    group_by(group) %>%
    summarise(
      mean_gene_n = mean(gene_n, na.rm = TRUE),
      median_gene_n = median(gene_n, na.rm = TRUE),
      mean_gene_density_perMb = mean(gene_density_perMb, na.rm = TRUE),
      .groups = "drop"
    )
}

write.csv(annotated, file.path(cfg$output_dir, "extreme_domain_gene_counts_per_domain.csv"), row.names = FALSE)
write.csv(summary_by_threshold, file.path(cfg$output_dir, "extreme_domain_gene_counts_summary_by_threshold.csv"), row.names = FALSE)
write.csv(summary_primary_vs_background, file.path(cfg$output_dir, "extreme_domain_gene_counts_top1_vs_background.csv"), row.names = FALSE)

p <- ggplot(annotated, aes(x = factor(threshold_prob), y = gene_count_sum)) +
  geom_violin(fill = "#9ecae1", color = NA, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, outlier.alpha = 0.15) +
  labs(
    title = "Genes per Extreme Barrier Domain",
    subtitle = "Distribution across extreme-K thresholds",
    x = "Threshold probability",
    y = "Gene count per domain"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(cfg$output_dir, "plot_extreme_domain_gene_counts_by_threshold.pdf"), p, width = 7, height = 5.2)

top1 <- summary_by_threshold %>% filter(abs(threshold_prob - 0.99) < 1e-12) %>% slice_head(n = 1)
if (!nrow(top1)) top1 <- summary_by_threshold %>% slice_head(n = 1)

cat("\n================ Extreme Domain Gene Counts ================\n")
for (i in seq_len(nrow(summary_by_threshold))) {
  row <- summary_by_threshold[i, ]
  cat("Threshold q", row$threshold_prob,
      ": n_domains=", row$n_domains,
      ", mean genes/domain=", signif(row$mean_genes_per_domain, 6),
      ", median=", signif(row$median_genes_per_domain, 6),
      ", range=[", signif(row$min_genes_per_domain, 6), ", ", signif(row$max_genes_per_domain, 6), "]",
      ", mean genes/window=", signif(row$mean_genes_per_window, 6),
      "\n", sep = "")
}
cat("Primary top 1% result: mean genes/domain=", signif(top1$mean_genes_per_domain, 6),
    ", median=", signif(top1$median_genes_per_domain, 6),
    ", range=[", signif(top1$min_genes_per_domain, 6), ", ", signif(top1$max_genes_per_domain, 6), "]\n", sep = "")
cat("Outputs written to: ", cfg$output_dir, "\n", sep = "")
