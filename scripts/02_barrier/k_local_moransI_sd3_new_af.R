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
  prior_driver_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_driver_outputs"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_moransI_outputs"),
  n_permutations = 199L,
  lags = c(1L, 2L, 5L, 10L, 20L, 50L, 100L)
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_rdata <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$object_name <- args[[2]]
if (length(args) >= 3 && nzchar(args[[3]])) cfg$output_dir <- args[[3]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
sanity <- init_sanity("k_local_moransI_sd3_new_af.R", cfg$output_dir)

# Biological meaning:
# Spatial autocorrelation asks whether neighboring windows have similar barrier
# scores more often than expected by chance. Positive Moran's I and slow decay
# of lag correlations support the idea that barrier loci cluster into domains
# rather than occurring independently one window at a time.

find_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
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

morans_i_line <- function(x) {
  keep <- is.finite(x)
  x <- x[keep]
  n <- length(x)
  if (n < 3L) {
    return(tibble(n = n, W = NA_real_, morans_I = NA_real_, expected_I = NA_real_))
  }
  z <- x - mean(x)
  denom <- sum(z^2)
  if (!is.finite(denom) || denom == 0) {
    return(tibble(n = n, W = 2 * (n - 1), morans_I = NA_real_, expected_I = -1 / (n - 1)))
  }
  num <- 2 * sum(z[-n] * z[-1])
  W <- 2 * (n - 1)
  tibble(n = n, W = W, morans_I = (n / W) * (num / denom), expected_I = -1 / (n - 1))
}

morans_i_perm <- function(x, n_perm = 199L) {
  obs <- morans_i_line(x)
  x <- x[is.finite(x)]
  if (nrow(obs) == 0L || !is.finite(obs$morans_I) || length(x) < 3L) {
    return(bind_cols(obs, tibble(p_value = NA_real_, p_greater = NA_real_, p_less = NA_real_, perm_mean = NA_real_, perm_sd = NA_real_)))
  }
  perm <- replicate(n_perm, morans_i_line(sample(x, replace = FALSE))$morans_I)
  tibble(
    n = obs$n,
    W = obs$W,
    morans_I = obs$morans_I,
    expected_I = obs$expected_I,
    p_value = mean(abs(perm) >= abs(obs$morans_I)),
    p_greater = mean(perm >= obs$morans_I),
    p_less = mean(perm <= obs$morans_I),
    perm_mean = mean(perm),
    perm_sd = stats::sd(perm)
  )
}

lag_corr_one <- function(x, lag_k) {
  keep <- is.finite(x)
  x <- x[keep]
  n <- length(x)
  if (n <= lag_k) return(tibble(n_pairs = 0L, correlation = NA_real_))
  x1 <- x[seq_len(n - lag_k)]
  x2 <- x[(lag_k + 1L):n]
  keep2 <- is.finite(x1) & is.finite(x2)
  if (sum(keep2) < 3L) return(tibble(n_pairs = sum(keep2), correlation = NA_real_))
  tibble(n_pairs = sum(keep2), correlation = suppressWarnings(stats::cor(x1[keep2], x2[keep2])))
}

lag_corr_perm <- function(x, lag_k, n_perm = 199L) {
  obs <- lag_corr_one(x, lag_k)
  x <- x[is.finite(x)]
  if (obs$n_pairs < 3L || !is.finite(obs$correlation)) {
    return(bind_cols(obs, tibble(p_value = NA_real_, perm_mean = NA_real_, perm_sd = NA_real_)))
  }
  perm <- replicate(n_perm, lag_corr_one(sample(x, replace = FALSE), lag_k)$correlation)
  bind_cols(obs, tibble(
    p_value = mean(abs(perm) >= abs(obs$correlation)),
    perm_mean = mean(perm, na.rm = TRUE),
    perm_sd = stats::sd(perm, na.rm = TRUE)
  ))
}

global_morans_blockdiag <- function(df, value_col, n_perm = 199L) {
  df2 <- df %>% filter(is.finite(.data[[value_col]])) %>% arrange(Chromosome_Names, window_index)
  split_vals <- split(df2[[value_col]], df2$Chromosome_Names)
  split_vals <- split_vals[lengths(split_vals) >= 3L]
  if (!length(split_vals)) {
    return(tibble(metric = value_col, n = 0L, W = NA_real_, morans_I = NA_real_, expected_I = NA_real_, p_value = NA_real_, perm_mean = NA_real_, perm_sd = NA_real_))
  }
  x <- unlist(split_vals, use.names = FALSE)
  n <- length(x)
  z <- x - mean(x)
  denom <- sum(z^2)
  edge_pairs <- lapply(split_vals, function(v) {
    v <- as.numeric(v)
    cbind(v[-length(v)], v[-1])
  })
  pairs_mat <- do.call(rbind, edge_pairs)
  num <- 2 * sum((pairs_mat[, 1] - mean(x)) * (pairs_mat[, 2] - mean(x)))
  W <- 2 * sum(vapply(split_vals, function(v) length(v) - 1L, integer(1)))
  obs_I <- (n / W) * (num / denom)
  perm <- replicate(n_perm, {
    perm_vals <- unlist(lapply(split_vals, sample, replace = FALSE), use.names = FALSE)
    z_perm <- perm_vals - mean(perm_vals)
    denom_perm <- sum(z_perm^2)
    idx <- 0L
    perm_num <- 0
    for (v in split_vals) {
      m <- length(v)
      block <- perm_vals[(idx + 1L):(idx + m)]
      idx <- idx + m
      perm_num <- perm_num + 2 * sum((block[-m] - mean(perm_vals)) * (block[-1] - mean(perm_vals)))
    }
    (n / W) * (perm_num / denom_perm)
  })
  tibble(
    metric = value_col,
    n = n,
    W = W,
    morans_I = obs_I,
    expected_I = -1 / (n - 1),
    p_value = mean(abs(perm) >= abs(obs_I)),
    perm_mean = mean(perm),
    perm_sd = stats::sd(perm)
  )
}

local_summary_path <- file.path(cfg$prior_barrier_dir, "barrier_recalibrated_local_summary.csv")
if (!file.exists(local_summary_path)) stop("Missing prior local summary: ", local_summary_path)
local_df <- read.csv(local_summary_path, stringsAsFactors = FALSE)
assert_true(nrow(local_df) > 0L, "Prior local summary is empty.", sanity, "nonempty_prior_local_summary")

load(cfg$input_rdata)
if (!exists(cfg$object_name, inherits = FALSE)) stop("Object not found: ", cfg$object_name)
dat_raw <- get(cfg$object_name, inherits = FALSE)
assert_true(nrow(as.data.frame(dat_raw)) > 0L, "Loaded dataset is empty.", sanity, "nonempty_input_dataset")

chr_col <- find_col(dat_raw, c("Chromosome_Names", "Chromosome", "chr", "chromosome"))
start_col <- find_col(dat_raw, c(".start_bp", "window_start", "start_bp", "start"))
end_col <- find_col(dat_raw, c(".end_bp", "window_end", "end_bp", "end"))
mid_col <- find_col(dat_raw, c("midpoint", ".mid_bp", "window_mid", "mid_bp", "mid"))

detected_columns <- tibble(
  role = c("chromosome", "start", "end", "mid"),
  detected_column = c(chr_col, start_col, end_col, mid_col)
)

dat <- as.data.frame(dat_raw) %>%
  mutate(
    Chromosome_Names = as.character(.data[[chr_col]]),
    .start_bp = as.numeric(.data[[start_col]]),
    .end_bp = if (!is.na(end_col)) as.numeric(.data[[end_col]]) else .start_bp + 9999,
    window_mid_bp = if (!is.na(mid_col)) as.numeric(.data[[mid_col]]) else (.start_bp + .end_bp) / 2
  ) %>%
  select(Chromosome_Names, .start_bp, .end_bp, window_mid_bp)
assert_has_cols(dat, c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"), sanity, "morans_input_data")
assert_unique_keys(dat, c("Chromosome_Names", ".start_bp", ".end_bp"), sanity, "morans_input_data")

dat <- dat %>%
  left_join(local_df %>% select(Chromosome_Names, .start_bp, .end_bp, window_mid_bp, B_no_delta_p, K_local_mean), by = c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp")) %>%
  mutate(Chromosome_Names = factor(Chromosome_Names, levels = order_chromosomes(Chromosome_Names))) %>%
  arrange(Chromosome_Names, .start_bp) %>%
  group_by(Chromosome_Names) %>%
  mutate(window_index = row_number()) %>%
  ungroup()
assert_missing_fraction(dat$B_no_delta_p, 0.1, sanity, "B_no_delta_p")

if (!"B_no_delta_p" %in% names(dat) || !any(is.finite(dat$B_no_delta_p))) stop("B_no_delta_p missing after merge.")
if (!"K_local_mean" %in% names(dat)) dat$K_local_mean <- NA_real_

morans_chr_B <- bind_rows(lapply(split(dat, dat$Chromosome_Names), function(df_chr) {
  res <- morans_i_perm(df_chr$B_no_delta_p, cfg$n_permutations)
  bind_cols(tibble(Chromosome_Names = as.character(df_chr$Chromosome_Names[1]), metric = "B_no_delta_p"), res)
}))

morans_chr_K <- bind_rows(lapply(split(dat, dat$Chromosome_Names), function(df_chr) {
  res <- morans_i_perm(df_chr$K_local_mean, cfg$n_permutations)
  bind_cols(tibble(Chromosome_Names = as.character(df_chr$Chromosome_Names[1]), metric = "K_local_mean"), res)
}))

morans_chr <- bind_rows(morans_chr_B, morans_chr_K)

morans_global <- bind_rows(
  global_morans_blockdiag(dat, "B_no_delta_p", cfg$n_permutations),
  global_morans_blockdiag(dat, "K_local_mean", cfg$n_permutations)
)

lags_tbl <- bind_rows(lapply(cfg$lags, function(lag_k) {
  per_chr_B <- bind_rows(lapply(split(dat, dat$Chromosome_Names), function(df_chr) {
    bind_cols(
      tibble(Chromosome_Names = as.character(df_chr$Chromosome_Names[1]), metric = "B_no_delta_p", lag = lag_k),
      lag_corr_perm(df_chr$B_no_delta_p, lag_k, cfg$n_permutations)
    )
  }))
  per_chr_K <- bind_rows(lapply(split(dat, dat$Chromosome_Names), function(df_chr) {
    bind_cols(
      tibble(Chromosome_Names = as.character(df_chr$Chromosome_Names[1]), metric = "K_local_mean", lag = lag_k),
      lag_corr_perm(df_chr$K_local_mean, lag_k, cfg$n_permutations)
    )
  }))
  pooled_B_pairs <- bind_rows(lapply(split(dat, dat$Chromosome_Names), function(df_chr) {
    v <- df_chr$B_no_delta_p[is.finite(df_chr$B_no_delta_p)]
    n <- length(v)
    if (n <= lag_k) return(tibble(x = numeric(), y = numeric()))
    tibble(x = v[seq_len(n - lag_k)], y = v[(lag_k + 1L):n])
  }))
  pooled_K_pairs <- bind_rows(lapply(split(dat, dat$Chromosome_Names), function(df_chr) {
    v <- df_chr$K_local_mean[is.finite(df_chr$K_local_mean)]
    n <- length(v)
    if (n <= lag_k) return(tibble(x = numeric(), y = numeric()))
    tibble(x = v[seq_len(n - lag_k)], y = v[(lag_k + 1L):n])
  }))
  pooled_B <- tibble(
    Chromosome_Names = "Genome_pooled",
    metric = "B_no_delta_p",
    lag = lag_k,
    n_pairs = nrow(pooled_B_pairs),
    correlation = if (nrow(pooled_B_pairs) >= 3L) cor(pooled_B_pairs$x, pooled_B_pairs$y) else NA_real_,
    p_value = NA_real_,
    perm_mean = NA_real_,
    perm_sd = NA_real_
  )
  pooled_K <- tibble(
    Chromosome_Names = "Genome_pooled",
    metric = "K_local_mean",
    lag = lag_k,
    n_pairs = nrow(pooled_K_pairs),
    correlation = if (nrow(pooled_K_pairs) >= 3L) cor(pooled_K_pairs$x, pooled_K_pairs$y) else NA_real_,
    p_value = NA_real_,
    perm_mean = NA_real_,
    perm_sd = NA_real_
  )
  bind_rows(per_chr_B, per_chr_K, pooled_B, pooled_K)
}))

domain_run_path <- file.path(cfg$prior_driver_dir, "k_local_runlength_summary.csv")
domain_run_tbl <- if (file.exists(domain_run_path)) read.csv(domain_run_path, stringsAsFactors = FALSE) else tibble()

lag_scale_summary <- lags_tbl %>%
  filter(Chromosome_Names == "Genome_pooled") %>%
  group_by(metric) %>%
  summarise(
    max_positive_lag = if (any(correlation > 0, na.rm = TRUE)) max(lag[correlation > 0], na.rm = TRUE) else NA_real_,
    first_lag_below_0.1 = if (any(correlation < 0.1, na.rm = TRUE)) min(lag[correlation < 0.1], na.rm = TRUE) else NA_real_,
    first_lag_below_0.05 = if (any(correlation < 0.05, na.rm = TRUE)) min(lag[correlation < 0.05], na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

domain_compare <- tibble()
if (nrow(domain_run_tbl) > 0L) {
  domain_compare <- domain_run_tbl %>%
    filter(cutoff_quantile == 0.95) %>%
    summarise(
      mean_high_domain_run = mean(mean_run_length, na.rm = TRUE),
      median_high_domain_run = median(mean_run_length, na.rm = TRUE),
      max_high_domain_run = max(max_run_length, na.rm = TRUE)
    ) %>%
    crossing(lag_scale_summary)
}

morans_plot <- morans_chr %>%
  filter(metric == "B_no_delta_p") %>%
  mutate(Chromosome_Names = factor(Chromosome_Names, levels = order_chromosomes(Chromosome_Names))) %>%
  ggplot(aes(x = Chromosome_Names, y = morans_I, fill = p_value < 0.05)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = unique(morans_chr$expected_I[morans_chr$metric == "B_no_delta_p"])[1], linetype = "dotted", color = "grey60") +
  scale_fill_manual(values = c("TRUE" = "red3", "FALSE" = "grey70"), guide = "none") +
  labs(x = NULL, y = "Moran's I for B_no_delta_p") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(cfg$output_dir, "plot_morans_I_by_chromosome.png"), morans_plot, width = 11, height = 5, dpi = 300)

lag_plot <- lags_tbl %>%
  filter(Chromosome_Names == "Genome_pooled") %>%
  ggplot(aes(x = lag, y = correlation, color = metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(values = c("B_no_delta_p" = "red3", "K_local_mean" = "dodgerblue3")) +
  scale_x_continuous(breaks = cfg$lags) +
  labs(x = "Lag (windows)", y = "Autocorrelation", color = "Metric") +
  theme_bw(base_size = 11)
ggsave(file.path(cfg$output_dir, "plot_autocorrelation_vs_distance.png"), lag_plot, width = 7, height = 5, dpi = 300)

lag_chr_plot <- lags_tbl %>%
  filter(metric == "B_no_delta_p", Chromosome_Names != "Genome_pooled") %>%
  mutate(Chromosome_Names = factor(Chromosome_Names, levels = order_chromosomes(Chromosome_Names))) %>%
  ggplot(aes(x = lag, y = correlation, group = Chromosome_Names)) +
  geom_line(alpha = 0.25, color = "grey55") +
  geom_line(data = subset(lags_tbl, metric == "B_no_delta_p" & Chromosome_Names == "Genome_pooled"), aes(x = lag, y = correlation), inherit.aes = FALSE, color = "red3", linewidth = 1.1) +
  geom_point(data = subset(lags_tbl, metric == "B_no_delta_p" & Chromosome_Names == "Genome_pooled"), aes(x = lag, y = correlation), inherit.aes = FALSE, color = "red3", size = 2) +
  scale_x_continuous(breaks = cfg$lags) +
  labs(x = "Lag (windows)", y = "Autocorrelation for B_no_delta_p") +
  theme_bw(base_size = 11)
ggsave(file.path(cfg$output_dir, "plot_autocorrelation_decay_by_chromosome.png"), lag_chr_plot, width = 7.5, height = 5, dpi = 300)

if (nrow(domain_compare) > 0L) {
  compare_plot <- domain_run_tbl %>%
    filter(cutoff_quantile == 0.95) %>%
    ggplot(aes(x = Chromosome_Names, y = max_run_length)) +
    geom_col(fill = "grey65") +
    geom_hline(yintercept = lag_scale_summary$max_positive_lag[lag_scale_summary$metric == "B_no_delta_p"], color = "red3", linetype = "dashed") +
    geom_hline(yintercept = lag_scale_summary$max_positive_lag[lag_scale_summary$metric == "K_local_mean"], color = "dodgerblue3", linetype = "dashed") +
    labs(x = NULL, y = "Max high-K domain run length (windows)") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(cfg$output_dir, "plot_domain_runlength_vs_autocorrelation_scale.png"), compare_plot, width = 11, height = 5, dpi = 300)
}

write.csv(detected_columns, file.path(cfg$output_dir, "detected_columns.csv"), row.names = FALSE)
write.csv(morans_chr, file.path(cfg$output_dir, "k_local_morans_I_per_chromosome.csv"), row.names = FALSE)
write.csv(morans_global, file.path(cfg$output_dir, "k_local_morans_I_global.csv"), row.names = FALSE)
write.csv(lags_tbl, file.path(cfg$output_dir, "k_local_B_autocorrelation_by_distance.csv"), row.names = FALSE)
write.csv(domain_compare, file.path(cfg$output_dir, "k_local_autocorrelation_domain_comparison.csv"), row.names = FALSE)

b_global <- morans_global %>% filter(metric == "B_no_delta_p") %>% slice_head(n = 1)
k_global <- morans_global %>% filter(metric == "K_local_mean") %>% slice_head(n = 1)
b_scale <- lag_scale_summary %>% filter(metric == "B_no_delta_p") %>% slice_head(n = 1)
k_scale <- lag_scale_summary %>% filter(metric == "K_local_mean") %>% slice_head(n = 1)
top_chr <- morans_chr %>% filter(metric == "B_no_delta_p") %>% arrange(desc(morans_I)) %>% slice_head(n = min(5L, nrow(filter(morans_chr, metric == "B_no_delta_p")))) %>% pull(Chromosome_Names)

cat("\nSpatial autocorrelation summary\n")
cat("Detected columns:\n")
apply(detected_columns, 1, function(x) cat(" -", x[["role"]], "=>", x[["detected_column"]], "\n"))
cat("\nIs B_no_delta_p significantly spatially autocorrelated? ",
    ifelse(nrow(b_global) == 1L && is.finite(b_global$morans_I) && b_global$p_value < 0.05 && b_global$morans_I > 0, "Yes", "No"),
    "; global Moran's I=", signif(b_global$morans_I, 4),
    ", p=", signif(b_global$p_value, 3), "\n", sep = "")
cat("At what distances does autocorrelation persist? ",
    "B_no_delta_p stays positive through lag ", b_scale$max_positive_lag,
    " and first drops below 0.1 at lag ", b_scale$first_lag_below_0.1, ".\n", sep = "")
if (nrow(domain_compare) == 1L) {
  cat("Does the autocorrelation scale match observed domain sizes? ",
      "High-K mean run length=", round(domain_compare$mean_high_domain_run, 1),
      " windows, max run length=", round(domain_compare$max_high_domain_run, 1),
      ", versus positive B autocorrelation through lag ", b_scale$max_positive_lag, ".\n", sep = "")
}
cat("Does K_local show stronger autocorrelation than B_i? ",
    ifelse(nrow(k_global) == 1L && is.finite(k_global$morans_I) && k_global$morans_I > b_global$morans_I, "Yes", "No"),
    "; Moran's I(B)=", signif(b_global$morans_I, 4),
    ", Moran's I(K_local)=", signif(k_global$morans_I, 4), "\n", sep = "")
cat("Chromosomes with strongest B_no_delta_p spatial clustering:", paste(top_chr, collapse = ", "), "\n")
cat("Saved outputs to:", cfg$output_dir, "\n")
finalize_sanity(
  sanity,
  files = c(
    "detected_columns.csv",
    "k_local_morans_I_per_chromosome.csv",
    "k_local_morans_I_global.csv",
    "k_local_B_autocorrelation_by_distance.csv",
    "k_local_autocorrelation_domain_comparison.csv",
    "plot_morans_I_by_chromosome.png",
    "plot_autocorrelation_vs_distance.png",
    "plot_autocorrelation_decay_by_chromosome.png"
  )
)
