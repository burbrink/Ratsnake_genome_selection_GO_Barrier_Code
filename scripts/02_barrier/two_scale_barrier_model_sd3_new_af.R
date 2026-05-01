#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})
script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

# Two-scale barrier model
# -----------------------
# Local scale:
#   P_local = exp(-alpha_local * K_local)
# calibrated empirically from independent mixing proxies.
#
# Global anti-collapse scale:
#   species persistence can be maintained by a limited number of very strong
#   barrier blocks, even if average compatibility is not extremely low.
# We therefore summarise extreme high-K domains rather than relying only on
# a genome-wide mean-K compatibility equation.

cfg <- list(
  local_summary_csv = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_2of3_corrected", "barrier_recalibrated_local_summary.csv"),
  empirical_alpha_csv = file.path(repo_root, "outputs", "barrier_analysis_outputs", "alpha_empirical_calibration_outputs", "empirical_alpha_estimates.csv"),
  empirical_alpha_proxy = "compat_composite",
  empirical_alpha_method = "direct_mse",
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "two_scale_barrier_outputs"),
  extreme_probs = c(0.99, 0.995, 0.999),
  top_plot_n_chromosomes = 8L
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$local_summary_csv <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$output_dir <- args[[2]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cfg$local_summary_csv)) stop("Local summary CSV not found: ", cfg$local_summary_csv, call. = FALSE)
dat0 <- read.csv(cfg$local_summary_csv, stringsAsFactors = FALSE, check.names = FALSE)

required <- c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp", "K_local_mean")
missing_req <- setdiff(required, names(dat0))
if (length(missing_req)) stop("Missing required columns: ", paste(missing_req, collapse = ", "), call. = FALSE)

dat <- dat0 %>%
  transmute(
    Chromosome_Names = as.character(Chromosome_Names),
    .start_bp = as.numeric(.start_bp),
    .end_bp = as.numeric(.end_bp),
    window_mid_bp = as.numeric(window_mid_bp),
    K_local_mean = as.numeric(K_local_mean),
    abs_delta_p = if ("abs_delta_p" %in% names(dat0)) as.numeric(abs_delta_p) else NA_real_,
    f3 = if ("f3" %in% names(dat0)) as.numeric(f3) else NA_real_,
    recombination = if ("recombination" %in% names(dat0)) as.numeric(recombination) else NA_real_,
    sweep_group = if ("sweep_group" %in% names(dat0)) as.character(sweep_group) else NA_character_
  ) %>%
  arrange(Chromosome_Names, .start_bp)

if (file.exists(cfg$empirical_alpha_csv)) {
  alpha_tbl <- read.csv(cfg$empirical_alpha_csv, stringsAsFactors = FALSE, check.names = FALSE)
  alpha_row <- alpha_tbl %>%
    filter(proxy == cfg$empirical_alpha_proxy, method == cfg$empirical_alpha_method) %>%
    slice_head(n = 1)
  if (nrow(alpha_row) == 0L) stop("Could not find requested empirical alpha row in ", cfg$empirical_alpha_csv, call. = FALSE)
  alpha_local <- alpha_row$alpha_estimate[[1]]
} else {
  stop("Empirical alpha estimate file not found: ", cfg$empirical_alpha_csv, call. = FALSE)
}

dat <- dat %>%
  mutate(
    P_local_empirical = ifelse(is.finite(K_local_mean), exp(-alpha_local * K_local_mean), NA_real_),
    log10_P_local_empirical = ifelse(is.finite(P_local_empirical), log10(P_local_empirical), NA_real_)
  )

get_runs <- function(flag) {
  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  tibble(run_value = r$values, run_start = starts, run_end = ends, n_windows = r$lengths)
}

build_domain_table <- function(df, prob) {
  cutoff <- unname(stats::quantile(df$K_local_mean, probs = prob, na.rm = TRUE))
  flagged <- df %>%
    mutate(extreme_flag = is.finite(K_local_mean) & K_local_mean >= cutoff) %>%
    group_by(Chromosome_Names) %>%
    group_modify(~{
      runs <- get_runs(.x$extreme_flag)
      runs <- runs %>% filter(run_value)
      if (nrow(runs) == 0L) return(tibble())
      bind_rows(lapply(seq_len(nrow(runs)), function(i) {
        rr <- runs[i, ]
        seg <- .x[rr$run_start:rr$run_end, , drop = FALSE]
        tibble(
          threshold_prob = prob,
          K_cutoff = cutoff,
          Chromosome_Names = seg$Chromosome_Names[[1]],
          domain_id = paste0(seg$Chromosome_Names[[1]], "_q", gsub("\\.", "", sprintf("%.3f", prob)), "_dom_", i),
          start_bp = min(seg$.start_bp, na.rm = TRUE),
          end_bp = max(seg$.end_bp, na.rm = TRUE),
          n_windows = nrow(seg),
          domain_length_bp = max(seg$.end_bp, na.rm = TRUE) - min(seg$.start_bp, na.rm = TRUE) + 1,
          mean_K_local = mean(seg$K_local_mean, na.rm = TRUE),
          max_K_local = max(seg$K_local_mean, na.rm = TRUE),
          mean_P_local_empirical = mean(seg$P_local_empirical, na.rm = TRUE),
          mean_abs_delta_p = mean(seg$abs_delta_p, na.rm = TRUE),
          mean_f3 = mean(seg$f3, na.rm = TRUE),
          mean_recombination = mean(seg$recombination, na.rm = TRUE),
          dominant_sweep_group = {
            sw <- seg$sweep_group[!is.na(seg$sweep_group) & nzchar(seg$sweep_group)]
            if (!length(sw)) NA_character_ else names(sort(table(sw), decreasing = TRUE))[1]
          }
        )
      }))
    }) %>%
    ungroup()
  list(cutoff = cutoff, domains = flagged)
}

domain_tables <- lapply(cfg$extreme_probs, function(prob) build_domain_table(dat, prob))
names(domain_tables) <- paste0("q", cfg$extreme_probs)

domains_all <- bind_rows(lapply(domain_tables, `[[`, "domains"))

domain_summary <- domains_all %>%
  group_by(threshold_prob, K_cutoff) %>%
  summarise(
    n_domains = n(),
    n_chromosomes_with_domains = n_distinct(Chromosome_Names),
    total_domain_windows = sum(n_windows, na.rm = TRUE),
    total_domain_bp = sum(domain_length_bp, na.rm = TRUE),
    mean_domain_windows = mean(n_windows, na.rm = TRUE),
    median_domain_windows = median(n_windows, na.rm = TRUE),
    max_domain_windows = max(n_windows, na.rm = TRUE),
    mean_domain_length_bp = mean(domain_length_bp, na.rm = TRUE),
    max_domain_length_bp = max(domain_length_bp, na.rm = TRUE),
    mean_K_local = mean(mean_K_local, na.rm = TRUE),
    max_K_local = max(max_K_local, na.rm = TRUE),
    mean_P_local_empirical = mean(mean_P_local_empirical, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(threshold_prob)

chr_summary_primary <- domains_all %>%
  filter(abs(threshold_prob - 0.995) < 1e-12) %>%
  group_by(Chromosome_Names) %>%
  summarise(
    n_extreme_domains = n(),
    total_extreme_windows = sum(n_windows, na.rm = TRUE),
    total_extreme_bp = sum(domain_length_bp, na.rm = TRUE),
    max_domain_windows = max(n_windows, na.rm = TRUE),
    max_domain_length_bp = max(domain_length_bp, na.rm = TRUE),
    mean_domain_K = mean(mean_K_local, na.rm = TRUE),
    min_domain_P_local = min(mean_P_local_empirical, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_extreme_domains), desc(max_domain_windows))

top_chr <- chr_summary_primary %>% slice_head(n = cfg$top_plot_n_chromosomes) %>% pull(Chromosome_Names)

genome_summary <- tibble(
  alpha_local = alpha_local,
  n_windows_with_K = sum(is.finite(dat$K_local_mean)),
  K_genome_mean = mean(dat$K_local_mean, na.rm = TRUE),
  P_overall_mean_empirical = exp(-alpha_local * mean(dat$K_local_mean, na.rm = TRUE)),
  log10_P_overall_mean_empirical = log10(exp(-alpha_local * mean(dat$K_local_mean, na.rm = TRUE))),
  n_domains_q99 = sum(abs(domains_all$threshold_prob - 0.99) < 1e-12, na.rm = TRUE),
  n_domains_q995 = sum(abs(domains_all$threshold_prob - 0.995) < 1e-12, na.rm = TRUE),
  n_domains_q999 = sum(abs(domains_all$threshold_prob - 0.999) < 1e-12, na.rm = TRUE)
)

write.csv(genome_summary, file.path(cfg$output_dir, "two_scale_barrier_genome_summary.csv"), row.names = FALSE)
write.csv(domains_all, file.path(cfg$output_dir, "two_scale_barrier_extreme_domain_table.csv"), row.names = FALSE)
write.csv(domain_summary, file.path(cfg$output_dir, "two_scale_barrier_extreme_domain_summary.csv"), row.names = FALSE)
write.csv(chr_summary_primary, file.path(cfg$output_dir, "two_scale_barrier_extreme_domain_chromosomes.csv"), row.names = FALSE)

p1 <- ggplot(dat %>% filter(is.finite(K_local_mean), is.finite(P_local_empirical)),
             aes(x = K_local_mean, y = P_local_empirical)) +
  geom_point(alpha = 0.08, size = 0.45, color = "#2b6cb0") +
  geom_smooth(method = "loess", se = FALSE, color = "#c53030", linewidth = 1) +
  labs(
    title = "Local Permeability Scale",
    subtitle = paste0("P_local = exp(-alpha_local * K_local), alpha_local = ", signif(alpha_local, 6)),
    x = "K_local_mean",
    y = "Empirical local compatibility"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(cfg$output_dir, "plot_two_scale_local_permeability.pdf"), p1, width = 7, height = 5.2)

plot_df <- dat %>%
  filter(Chromosome_Names %in% top_chr, is.finite(K_local_mean)) %>%
  left_join(
    domains_all %>%
      filter(abs(threshold_prob - 0.995) < 1e-12) %>%
      transmute(Chromosome_Names, start_bp, end_bp, domain_type = "Extreme high-K domain"),
    by = "Chromosome_Names"
  ) %>%
  mutate(in_extreme_domain = ifelse(is.na(start_bp), FALSE, .start_bp >= start_bp & .end_bp <= end_bp))

p2 <- ggplot(plot_df, aes(x = window_mid_bp / 1e6, y = K_local_mean)) +
  geom_point(aes(color = in_extreme_domain), alpha = 0.3, size = 0.45) +
  facet_wrap(~ Chromosome_Names, scales = "free_x") +
  scale_color_manual(values = c("FALSE" = "grey55", "TRUE" = "#c53030")) +
  labs(
    title = "Extreme High-K Domains Along Chromosomes",
    subtitle = "Primary threshold: top 0.5% of K_local_mean",
    x = "Genomic position (Mb)",
    y = "K_local_mean",
    color = "Extreme domain"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(cfg$output_dir, "plot_two_scale_extreme_domains_top_chromosomes.pdf"), p2, width = 10, height = 7)

p3 <- ggplot(domain_summary, aes(x = factor(threshold_prob), y = n_domains)) +
  geom_col(fill = "#8c2d04") +
  labs(
    title = "Number of Extreme High-K Domains",
    subtitle = "Genome-wide count of independent extreme-K areas",
    x = "Threshold quantile",
    y = "Number of domains"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(cfg$output_dir, "plot_two_scale_extreme_domain_counts.pdf"), p3, width = 6.5, height = 4.8)

primary <- domain_summary %>% filter(abs(threshold_prob - 0.995) < 1e-12) %>% slice_head(n = 1)
if (!nrow(primary)) primary <- domain_summary %>% slice_head(n = 1)

cat("\n================ Two-Scale Barrier Model ================\n")
cat("Empirical alpha_local: ", signif(alpha_local, 6), "\n", sep = "")
cat("Genome-wide mean K: ", signif(genome_summary$K_genome_mean, 6), "\n", sep = "")
cat("Empirical mean-scale P_overall: ", signif(genome_summary$P_overall_mean_empirical, 6), "\n", sep = "")
cat("Extreme high-K areas counted as consecutive windows above a global threshold within chromosomes.\n")
for (i in seq_len(nrow(domain_summary))) {
  row <- domain_summary[i, ]
  cat("Threshold q", row$threshold_prob, " (K >= ", signif(row$K_cutoff, 6), "): ",
      row$n_domains, " independent extreme domains across ",
      row$n_chromosomes_with_domains, " chromosomes; max run = ",
      row$max_domain_windows, " windows (", signif(row$max_domain_length_bp / 1e6, 4), " Mb)\n",
      sep = "")
}
cat("Primary extreme threshold (top 0.5%): ", primary$n_domains, " independent high-K areas across the genome.\n", sep = "")
cat("These areas are counted genome-wide, but domain boundaries are constrained within chromosomes so non-adjacent chromosomes are not merged.\n")
cat("Outputs written to: ", cfg$output_dir, "\n", sep = "")
