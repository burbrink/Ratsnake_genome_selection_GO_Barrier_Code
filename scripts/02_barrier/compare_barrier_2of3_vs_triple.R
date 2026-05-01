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
  label_rdata = file.path(repo_root, "outputs", "selection", "barrier_label_objects.rdata"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_comparison_outputs"),
  tier_output_dirs = c(
    "Majority 2-of-3" = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_majority_2of3"),
    "GMM/RF" = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_gmm_rf"),
    "Best continuous 3-of-3" = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_best_continuous_3of3"),
    "Margin >= 0.25" = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_margin_025"),
    "Margin >= 0.50" = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_margin_050"),
    "Margin >= 1.00" = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_margin_100"),
    "Top 1%" = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_top1pct"),
    "Top 500" = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_tail500")
  )
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$label_rdata <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$output_dir <- args[[2]]

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
sanity <- init_sanity("compare_barrier_2of3_vs_triple.R", cfg$output_dir)
assert_true(file.exists(cfg$label_rdata), paste("Missing label rdata:", cfg$label_rdata), sanity, "comparison_label_rdata_exists")

read_csv_base <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
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

harmonize_model_label <- function(x) {
  x_chr <- as.character(x)
  dplyr::case_when(
    x_chr %in% c("Neutral/other", "Neutral", "Neutral / other") ~ "Neutral / equilibrium",
    x_chr %in% c("Balancing-like", "Balancing", "Balancing_like") ~ "Balancing selection",
    x_chr %in% c("Sweep-like", "Sweep", "Sweep_like") ~ "Geographic sweep",
    x_chr %in% c("DGF-like", "DGF", "DGF_like") ~ "Divergence with gene flow",
    x_chr %in% c("Allopatry-like", "Allopatry", "Allopatry_like") ~ "Allopatry-like",
    TRUE ~ x_chr
  )
}

majority_label <- function(gmm_label, rf_label, cont_label) {
  out <- rep(NA_character_, length(gmm_label))
  for (i in seq_along(out)) {
    vals <- c(gmm_label[i], rf_label[i], cont_label[i])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (!length(vals)) next
    tt <- sort(table(vals), decreasing = TRUE)
    if (as.integer(tt[1]) >= 2L) out[i] <- names(tt)[1]
  }
  out
}

flag_label <- function(flag, label) {
  flag <- !is.na(flag) & as.logical(flag)
  ifelse(flag, label, NA_character_)
}

extract_validation_summary <- function(path, rule_label) {
  x <- read_csv_base(path)
  if (is.null(x)) return(tibble())
  cor_rows <- x %>%
    filter(is.na(term), analysis %in% c(
      "abs_delta_p ~ B_no_delta_p",
      "abs_delta_p ~ K_local_mean",
      "recombination ~ B_no_delta_p",
      "f3 ~ K_local_mean"
    )) %>%
    transmute(rule = rule_label, analysis, metric = model_type, estimate, p_value, n)

  lm_rows <- x %>%
    filter(model_type == "lm", term %in% c("B_no_delta_p", "K_local_mean"), analysis %in% c(
      "abs_delta_p ~ B_no_delta_p",
      "abs_delta_p ~ K_local_mean",
      "abs_delta_p ~ B_no_delta_p + K_local_mean"
    )) %>%
    transmute(rule = rule_label, analysis, metric = paste0("lm_", term), estimate, p_value, n, r_squared, adj_r_squared, aic)

  bind_rows(cor_rows, lm_rows)
}

extract_best_models <- function(path, rule_label) {
  x <- read_csv_base(path)
  if (is.null(x)) return(tibble())
  x %>%
    filter(model_type == "lm") %>%
    group_by(analysis) %>%
    summarise(
      rule = rule_label,
      r_squared = max(r_squared, na.rm = TRUE),
      adj_r_squared = max(adj_r_squared, na.rm = TRUE),
      aic = min(aic, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(grepl("^(abs_delta_p|f3)", analysis))
}

extract_chromosomes <- function(path, rule_label) {
  x <- read_csv_base(path)
  if (is.null(x)) return(tibble())
  x %>% mutate(rule = rule_label)
}

load(cfg$label_rdata)
if (!exists("agree_dat")) stop("Label file does not contain agree_dat.", call. = FALSE)
agree_dat <- as.data.frame(agree_dat)

agree_dat <- agree_dat %>%
  mutate(
    gmm_label_h = harmonize_model_label(gmm_label),
    rf_label_h = harmonize_model_label(rf_label),
    cont_label_h2 = harmonize_model_label(cont_label_h),
    cont_label_top1pct_h2 = harmonize_model_label(cont_label_top1pct_h),
    cont_label_tail500_h2 = harmonize_model_label(cont_label_tail500_h)
  )

labels_long <- bind_rows(
  tibble(rule = "Majority 2-of-3", label = majority_label(agree_dat$gmm_label_h, agree_dat$rf_label_h, agree_dat$cont_label_h2)),
  tibble(rule = "GMM/RF", label = flag_label(agree_dat$gmm_rf_match, agree_dat$rf_label_h)),
  tibble(rule = "Best continuous 3-of-3", label = flag_label(agree_dat$triple_match, agree_dat$cont_label_h2)),
  tibble(rule = "Margin >= 0.25", label = flag_label(agree_dat$triple_margin_025, agree_dat$cont_label_h2)),
  tibble(rule = "Margin >= 0.50", label = flag_label(agree_dat$triple_margin_050, agree_dat$cont_label_h2)),
  tibble(rule = "Margin >= 1.00", label = flag_label(agree_dat$triple_margin_100, agree_dat$cont_label_h2)),
  tibble(rule = "Top 1%", label = flag_label(agree_dat$triple_top1pct, agree_dat$cont_label_top1pct_h2)),
  tibble(rule = "Top 500", label = flag_label(agree_dat$triple_tail500, agree_dat$cont_label_tail500_h2))
)

rule_levels <- names(cfg$tier_output_dirs)
labels_long$rule <- factor(labels_long$rule, levels = rule_levels)

label_counts <- labels_long %>%
  mutate(label = ifelse(is.na(label) | !nzchar(label), "No consensus", label)) %>%
  count(rule, label, name = "n") %>%
  group_by(rule) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

existing_dirs <- cfg$tier_output_dirs[file.exists(cfg$tier_output_dirs)]
validation_summary <- bind_rows(lapply(names(existing_dirs), function(rule) {
  extract_validation_summary(file.path(existing_dirs[[rule]], "barrier_validation_models.csv"), rule)
}))
best_models <- bind_rows(lapply(names(existing_dirs), function(rule) {
  extract_best_models(file.path(existing_dirs[[rule]], "barrier_validation_models.csv"), rule)
}))
chrom_summary <- bind_rows(lapply(names(existing_dirs), function(rule) {
  extract_chromosomes(file.path(existing_dirs[[rule]], "barrier_recalibrated_chromosome_summary.csv"), rule)
}))

write.csv(label_counts, file.path(cfg$output_dir, "comparison_label_counts.csv"), row.names = FALSE)
write.csv(validation_summary, file.path(cfg$output_dir, "comparison_validation_metrics.csv"), row.names = FALSE)
write.csv(best_models, file.path(cfg$output_dir, "comparison_best_models.csv"), row.names = FALSE)
if (nrow(chrom_summary) > 0L) {
  chrom_summary$Chromosome_Names <- factor(chrom_summary$Chromosome_Names, levels = order_chromosomes(chrom_summary$Chromosome_Names))
  write.csv(chrom_summary %>% arrange(rule, Chromosome_Names), file.path(cfg$output_dir, "comparison_chromosome_summary.csv"), row.names = FALSE)
} else {
  write.csv(chrom_summary, file.path(cfg$output_dir, "comparison_chromosome_summary.csv"), row.names = FALSE)
}

plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

class_cols <- c(
  "Balancing selection" = "#1f77b4",
  "Geographic sweep" = "#d62728",
  "Divergence with gene flow" = "#ff7f0e",
  "Allopatry-like" = "#9467bd",
  "Neutral / equilibrium" = "grey60",
  "No consensus" = "grey90"
)

p1 <- ggplot(label_counts, aes(rule, n, fill = label)) +
  geom_col(color = "grey30", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = class_cols, drop = FALSE) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Agreement-tier label counts", x = NULL, y = "Windows", fill = "Class") +
  plot_theme
ggsave(file.path(cfg$output_dir, "figure_label_counts_by_agreement_tier.png"), p1, width = 10, height = 6.2, dpi = 300)

p2 <- ggplot(label_counts, aes(rule, frac, fill = label)) +
  geom_col(color = "grey30", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = class_cols, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Agreement-tier label composition", x = NULL, y = "Fraction", fill = "Class") +
  plot_theme
ggsave(file.path(cfg$output_dir, "figure_label_fractions_by_agreement_tier.png"), p2, width = 10, height = 6.2, dpi = 300)

if (nrow(best_models) > 0L) {
  p3 <- best_models %>%
    filter(grepl("^abs_delta_p", analysis)) %>%
    ggplot(aes(analysis, adj_r_squared, fill = rule)) +
    geom_col(position = "dodge", color = "grey30") +
    labs(title = "abs_delta_p model fit across agreement tiers", x = NULL, y = "Adjusted R^2", fill = "Agreement tier") +
    plot_theme
  ggsave(file.path(cfg$output_dir, "figure_abs_delta_model_fit_by_agreement_tier.png"), p3, width = 11, height = 6, dpi = 300)
}

if (nrow(chrom_summary) > 0L) {
  p4 <- ggplot(chrom_summary, aes(Chromosome_Names, mean_B_chr, fill = rule)) +
    geom_col(position = "dodge", color = "grey30", linewidth = 0.15) +
    labs(title = "Chromosome mean barrier score across agreement tiers", x = "Chromosome", y = "mean_B_chr", fill = "Agreement tier") +
    plot_theme
  ggsave(file.path(cfg$output_dir, "figure_chromosome_mean_B_by_agreement_tier.png"), p4, width = 12, height = 6.5, dpi = 300)
}

message("Comparison outputs written to: ", normalizePath(cfg$output_dir))
finalize_sanity(
  sanity,
  files = c(
    "comparison_label_counts.csv",
    "comparison_validation_metrics.csv",
    "comparison_best_models.csv",
    "comparison_chromosome_summary.csv",
    "figure_label_counts_by_agreement_tier.png",
    "figure_label_fractions_by_agreement_tier.png"
  )
)
