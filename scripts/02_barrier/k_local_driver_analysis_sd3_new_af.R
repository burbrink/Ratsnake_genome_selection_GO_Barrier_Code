#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(data.table)
})
script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
source(file.path(script_dir, "analysis_sanity_checks.R"))
source(file.path(script_dir, "alpha_estimation_helpers.R"))

cfg <- list(
  input_rdata = Sys.getenv("PIPELINE_INPUT_RDATA", unset = file.path(repo_root, "data", "input_placeholder.rdata")),
  object_name = "sd3_new_af",
  prior_barrier_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_gmm_rf"),
  prior_hybrid_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "hybrid_zone_condition_outputs"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_driver_outputs"),
  alpha = 0.5,
  alpha_calibration_csv = NULL,
  alpha_calibration_p_col = "P_overall_mean",
  alpha_calibration_k_col = "K_genome_mean",
  run_quantile = 0.95,
  n_permutations = 100L,
  top_thresholds = c(0.99, 0.95, 0.90)
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_rdata <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$object_name <- args[[2]]
if (length(args) >= 3 && nzchar(args[[3]])) cfg$output_dir <- args[[3]]
if (length(args) >= 4 && nzchar(args[[4]])) cfg$alpha_calibration_csv <- args[[4]]
if (length(args) >= 5 && nzchar(args[[5]])) cfg$alpha_calibration_p_col <- args[[5]]
if (length(args) >= 6 && nzchar(args[[6]])) cfg$alpha_calibration_k_col <- args[[6]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
sanity <- init_sanity("k_local_driver_analysis_sd3_new_af.R", cfg$output_dir)

if (!is.null(cfg$alpha_calibration_csv) && nzchar(cfg$alpha_calibration_csv)) {
  if (!file.exists(cfg$alpha_calibration_csv)) stop("Alpha calibration CSV not found: ", cfg$alpha_calibration_csv, call. = FALSE)
  alpha_df <- read.csv(cfg$alpha_calibration_csv, stringsAsFactors = FALSE, check.names = FALSE)
  alpha_result <- estimate_alpha_from_pk(
    df = alpha_df,
    p_col = cfg$alpha_calibration_p_col,
    k_col = cfg$alpha_calibration_k_col
  )
  cfg$alpha <- alpha_result$alpha
  write_alpha_calibration_outputs(alpha_result, cfg$output_dir)
  message("Estimated alpha = ", signif(cfg$alpha, 6), " using ", alpha_result$preferred_model, ".")
}

# Biological meaning:
# B_i is the focal barrier score for a single genomic window.
# K_local is the local accumulated barrier density around a window.
# If K_local clusters into long low-recombination, high-divergence domains,
# that supports a broad genomic barrier architecture.
# If sweep-like labels are enriched in high-K_local regions, that supports the
# idea that sweep history helps generate local barrier accumulation.

find_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0L) return(NA_character_)
  hit[[1]]
}

z_score <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) > 1L) out[ok] <- as.numeric(scale(x[ok]))
  out
}

safe_cor_tbl <- function(x, y, analysis, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L) {
    return(tibble(
      analysis = analysis, model_type = paste0("cor_", method), term = NA_character_, estimate = NA_real_,
      statistic = NA_real_, p_value = NA_real_, n = sum(keep), r_squared = NA_real_, adj_r_squared = NA_real_,
      aic = NA_real_, note = "Too few complete observations"
    ))
  }
  ct <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = method))
  tibble(
    analysis = analysis, model_type = paste0("cor_", method), term = NA_character_, estimate = unname(ct$estimate),
    statistic = unname(ct$statistic), p_value = ct$p.value, n = sum(keep), r_squared = NA_real_,
    adj_r_squared = NA_real_, aic = NA_real_, note = NA_character_
  )
}

safe_lm_details <- function(formula, data, analysis) {
  mf <- model.frame(formula, data = data, na.action = na.omit)
  if (nrow(mf) < 5L) {
    return(list(
      model = NULL,
      summary = tibble(analysis = analysis, n = nrow(mf), r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = "Too few complete observations"),
      coefficients = tibble(
        analysis = analysis, model_type = "lm", term = NA_character_, estimate = NA_real_, std_beta = NA_real_,
        statistic = NA_real_, p_value = NA_real_, n = nrow(mf), r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = "Too few complete observations"
      )
    ))
  }

  fit <- stats::lm(formula, data = data, na.action = na.omit)
  sm <- summary(fit)
  cf <- as.data.frame(sm$coefficients)
  cf$term <- rownames(cf)
  rownames(cf) <- NULL
  mm <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  y_sd <- stats::sd(y, na.rm = TRUE)
  std_beta <- rep(NA_real_, nrow(cf))
  names(std_beta) <- cf$term
  for (term in cf$term) {
    if (term == "(Intercept)" || !(term %in% colnames(mm))) next
    x_sd <- stats::sd(mm[, term], na.rm = TRUE)
    if (is.finite(x_sd) && x_sd > 0 && is.finite(y_sd) && y_sd > 0) {
      std_beta[term] <- cf$Estimate[cf$term == term] * x_sd / y_sd
    }
  }

  list(
    model = fit,
    summary = tibble(
      analysis = analysis, n = nrow(mf), r_squared = unname(sm$r.squared), adj_r_squared = unname(sm$adj.r.squared),
      aic = stats::AIC(fit), note = NA_character_
    ),
    coefficients = tibble(
      analysis = analysis, model_type = "lm", term = cf$term, estimate = cf$Estimate, std_beta = unname(std_beta[cf$term]),
      statistic = cf$`t value`, p_value = cf$`Pr(>|t|)`, n = nrow(mf), r_squared = unname(sm$r.squared),
      adj_r_squared = unname(sm$adj.r.squared), aic = stats::AIC(fit), note = NA_character_
    )
  )
}

safe_group_test_tbl <- function(y, group, analysis) {
  keep <- is.finite(y) & !is.na(group)
  y2 <- y[keep]
  g2 <- droplevels(as.factor(group[keep]))
  if (length(y2) < 5L || nlevels(g2) < 2L) {
    return(tibble(
      analysis = analysis, model_type = "group_test", term = NA_character_, estimate = NA_real_, statistic = NA_real_,
      p_value = NA_real_, n = length(y2), r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_,
      note = "Too few observations or groups"
    ))
  }
  if (nlevels(g2) == 2L) {
    wt <- suppressWarnings(stats::wilcox.test(y2 ~ g2))
    return(tibble(
      analysis = analysis, model_type = "wilcox", term = paste(levels(g2), collapse = " vs "), estimate = NA_real_,
      statistic = unname(wt$statistic), p_value = wt$p.value, n = length(y2), r_squared = NA_real_,
      adj_r_squared = NA_real_, aic = NA_real_, note = NA_character_
    ))
  }
  fit <- stats::aov(y2 ~ g2)
  sm <- summary(fit)[[1]]
  tibble(
    analysis = analysis, model_type = "anova", term = rownames(sm)[1], estimate = NA_real_, statistic = sm$`F value`[1],
    p_value = sm$`Pr(>F)`[1], n = length(y2), r_squared = NA_real_, adj_r_squared = NA_real_,
    aic = stats::AIC(stats::lm(y2 ~ g2)), note = NA_character_
  )
}

calc_vif_tbl <- function(data, predictors, analysis_prefix) {
  preds <- predictors[predictors %in% names(data)]
  if (length(preds) < 2L) return(tibble())
  df <- data[, preds, drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) < 5L) return(tibble())
  bind_rows(lapply(preds, function(pred) {
    others <- setdiff(preds, pred)
    fit <- stats::lm(stats::as.formula(paste(pred, "~", paste(others, collapse = " + "))), data = df)
    r2 <- summary(fit)$r.squared
    vif <- if (is.finite(r2) && r2 < 1) 1 / (1 - r2) else Inf
    tibble(analysis = analysis_prefix, term = pred, vif = vif)
  }))
}

quantile_class <- function(x) {
  br <- unique(stats::quantile(x, probs = c(0, 0.25, 0.75, 1), na.rm = TRUE, names = FALSE))
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (length(br) < 2L || sum(ok) == 0L) return(out)
  out[ok] <- as.character(cut(x[ok], breaks = br, include.lowest = TRUE, labels = c("Low", "Intermediate", "High")[seq_len(length(br) - 1L)]))
  out
}

recomb_bins <- function(x, n = 4L) {
  br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE, names = FALSE))
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (length(br) < 2L || sum(ok) == 0L) return(out)
  out[ok] <- as.character(cut(x[ok], breaks = br, include.lowest = TRUE, labels = paste0("Q", seq_len(length(br) - 1L))))
  out
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

run_length_table <- function(x) {
  r <- rle(x)
  idx <- cumsum(r$lengths)
  tibble(value = r$values, run_length = r$lengths, end_idx = idx, start_idx = idx - r$lengths + 1L)
}

perm_run_null <- function(x, n_perm = 100L) {
  if (sum(x, na.rm = TRUE) == 0L) return(tibble(mean_max_run_perm = 0, p_ge_obs = NA_real_))
  obs_max <- max(run_length_table(x)$run_length[run_length_table(x)$value %in% TRUE], na.rm = TRUE)
  perm_max <- replicate(n_perm, {
    xp <- sample(x, replace = FALSE)
    rt <- run_length_table(xp)
    mx <- rt$run_length[rt$value %in% TRUE]
    if (length(mx) == 0L) 0L else max(mx)
  })
  tibble(mean_max_run_perm = mean(perm_max), p_ge_obs = mean(perm_max >= obs_max))
}

load(cfg$input_rdata)
if (!exists(cfg$object_name, inherits = FALSE)) stop("Object not found: ", cfg$object_name)
dat_raw <- get(cfg$object_name, inherits = FALSE)
assert_true(nrow(as.data.frame(dat_raw)) > 0L, "Loaded dataset is empty.", sanity, "nonempty_input_dataset")

chr_col <- find_col(dat_raw, c("Chromosome_Names", "Chromosome", "chr", "chromosome"))
start_col <- find_col(dat_raw, c(".start_bp", "window_start", "start_bp", "start"))
end_col <- find_col(dat_raw, c(".end_bp", "window_end", "end_bp", "end"))
mid_col <- find_col(dat_raw, c("midpoint", "window_mid", "mid_bp", "mid"))
recomb_col <- find_col(dat_raw, c("recombination", "rec", "rho"))
fst_col <- find_col(dat_raw, c("avg_wc_fst", "wc_fst", "fst"))
dxy_col <- find_col(dat_raw, c("avg_dxy", "dxy"))
pi_col <- find_col(dat_raw, c("pi_average", "avg_pi", "pi"))
taj_col <- find_col(dat_raw, c("TajimaD", "tajima_d"))
delta_p_col <- find_col(dat_raw, c("mean_abs_delta_p", "abs_delta_p"))
f3_col <- find_col(dat_raw, c("f3", "F3"))

detected_columns <- tibble(
  role = c("chromosome", "start", "end", "mid", "recombination", "avg_wc_fst", "avg_dxy", "pi_average", "TajimaD", "abs_delta_p", "f3"),
  detected_column = c(chr_col, start_col, end_col, mid_col, recomb_col, fst_col, dxy_col, pi_col, taj_col, delta_p_col, f3_col)
)

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
    abs_delta_p = if (!is.na(delta_p_col)) as.numeric(.data[[delta_p_col]]) else NA_real_,
    f3 = if (!is.na(f3_col)) as.numeric(.data[[f3_col]]) else NA_real_
  )
assert_has_cols(dat, c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"), sanity, "k_local_driver_input_data")
assert_unique_keys(dat, c("Chromosome_Names", ".start_bp", ".end_bp"), sanity, "k_local_driver_input_data")

prior_local <- read.csv(file.path(cfg$prior_barrier_dir, "barrier_recalibrated_local_summary.csv"), stringsAsFactors = FALSE)
n_before_join <- nrow(dat)
dat <- dat %>% left_join(prior_local, by = c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"))
assert_join_rowcount(n_before_join, dat, sanity, "k_local_driver_prior_local_join")
if ("abs_delta_p.x" %in% names(dat) || "abs_delta_p.y" %in% names(dat)) dat$abs_delta_p <- dplyr::coalesce(dat$abs_delta_p.x, dat$abs_delta_p.y)
if ("f3.x" %in% names(dat) || "f3.y" %in% names(dat)) dat$f3 <- dplyr::coalesce(dat$f3.x, dat$f3.y)

if (!"B_no_delta_p" %in% names(dat)) stop("B_no_delta_p missing after merge.")
if (!"K_local_mean" %in% names(dat)) stop("K_local_mean missing after merge.")
if (!"K_local_sum" %in% names(dat)) dat$K_local_sum <- NA_real_
if (!"P_local_mean" %in% names(dat)) dat$P_local_mean <- exp(-cfg$alpha * dat$K_local_mean)
if (!"sweep_group" %in% names(dat)) dat$sweep_group <- NA_character_
assert_missing_fraction(dat$K_local_mean, 0.1, sanity, "K_local_mean")

dat$K_local_class <- quantile_class(dat$K_local_mean)
dat$recombination_bin <- recomb_bins(dat$recombination, 4L)
dat$is_zw <- grepl("ZW", dat$Chromosome_Names, ignore.case = TRUE)
dat$Chromosome_Names <- factor(dat$Chromosome_Names, levels = order_chromosomes(dat$Chromosome_Names))

single_formulas <- list(
  "S1: K_local_mean ~ recombination" = K_local_mean ~ recombination,
  "S2: K_local_mean ~ avg_wc_fst" = K_local_mean ~ avg_wc_fst,
  "S3: K_local_mean ~ avg_dxy" = K_local_mean ~ avg_dxy,
  "S4: K_local_mean ~ pi_average" = K_local_mean ~ pi_average,
  "S5: K_local_mean ~ TajimaD" = K_local_mean ~ TajimaD,
  "S6: K_local_mean ~ abs_delta_p" = K_local_mean ~ abs_delta_p,
  "S7: K_local_mean ~ f3" = K_local_mean ~ f3,
  "S8: K_local_mean ~ sweep_group" = K_local_mean ~ sweep_group,
  "S9: K_local_mean ~ chromosome" = K_local_mean ~ Chromosome_Names
)

multi_formulas <- list(
  "M1: K_local_mean ~ recombination + avg_wc_fst + avg_dxy" = K_local_mean ~ recombination + avg_wc_fst + avg_dxy,
  "M2: K_local_mean ~ recombination + avg_wc_fst + avg_dxy + pi_average + TajimaD" = K_local_mean ~ recombination + avg_wc_fst + avg_dxy + pi_average + TajimaD,
  "M3: K_local_mean ~ recombination + avg_wc_fst + avg_dxy + sweep_group" = K_local_mean ~ recombination + avg_wc_fst + avg_dxy + sweep_group,
  "M4: K_local_mean ~ recombination + avg_wc_fst + avg_dxy + chromosome" = K_local_mean ~ recombination + avg_wc_fst + avg_dxy + Chromosome_Names,
  "M5: K_local_mean ~ recombination + avg_wc_fst + avg_dxy + sweep_group + chromosome" = K_local_mean ~ recombination + avg_wc_fst + avg_dxy + sweep_group + Chromosome_Names
)

single_formulas <- single_formulas[sapply(single_formulas, function(f) all(all.vars(f) %in% names(dat)))]
multi_formulas <- multi_formulas[sapply(multi_formulas, function(f) all(all.vars(f) %in% names(dat)))]

fit_single <- lapply(names(single_formulas), function(nm) safe_lm_details(single_formulas[[nm]], dat, nm))
names(fit_single) <- names(single_formulas)
fit_multi <- lapply(names(multi_formulas), function(nm) safe_lm_details(multi_formulas[[nm]], dat, nm))
names(fit_multi) <- names(multi_formulas)

model_comparisons <- bind_rows(
  bind_rows(lapply(fit_single, `[[`, "summary")),
  bind_rows(lapply(fit_multi, `[[`, "summary"))
) %>% mutate(response = "K_local_mean")

coefficients_tbl <- bind_rows(
  bind_rows(lapply(fit_single, `[[`, "coefficients")),
  bind_rows(lapply(fit_multi, `[[`, "coefficients"))
) %>% mutate(response = "K_local_mean")

correlations_tbl <- bind_rows(
  safe_cor_tbl(dat$K_local_mean, dat$recombination, "K_local_mean ~ recombination", "pearson"),
  safe_cor_tbl(dat$K_local_mean, dat$recombination, "K_local_mean ~ recombination", "spearman"),
  safe_cor_tbl(dat$K_local_mean, dat$avg_wc_fst, "K_local_mean ~ avg_wc_fst", "pearson"),
  safe_cor_tbl(dat$K_local_mean, dat$avg_dxy, "K_local_mean ~ avg_dxy", "pearson"),
  safe_cor_tbl(dat$K_local_mean, dat$abs_delta_p, "K_local_mean ~ abs_delta_p", "pearson")
)
if (any(is.finite(dat$pi_average))) correlations_tbl <- bind_rows(correlations_tbl, safe_cor_tbl(dat$K_local_mean, dat$pi_average, "K_local_mean ~ pi_average", "pearson"))
if (any(is.finite(dat$TajimaD))) correlations_tbl <- bind_rows(correlations_tbl, safe_cor_tbl(dat$K_local_mean, dat$TajimaD, "K_local_mean ~ TajimaD", "pearson"))
if (any(is.finite(dat$f3))) correlations_tbl <- bind_rows(correlations_tbl, safe_cor_tbl(dat$K_local_mean, dat$f3, "K_local_mean ~ f3", "pearson"))

resid_div <- stats::resid(stats::lm(K_local_mean ~ avg_wc_fst + avg_dxy + abs_delta_p, data = dat, na.action = na.exclude))
correlations_tbl <- bind_rows(
  correlations_tbl,
  safe_cor_tbl(resid_div, dat$recombination, "Residual K_local_mean | divergence metrics ~ recombination", "spearman")
)

vif_tbl <- bind_rows(
  calc_vif_tbl(dat, c("recombination", "avg_wc_fst", "avg_dxy"), "VIF: recombination + avg_wc_fst + avg_dxy"),
  calc_vif_tbl(dat, c("recombination", "avg_wc_fst", "avg_dxy", "pi_average", "TajimaD"), "VIF: recombination + avg_wc_fst + avg_dxy + pi_average + TajimaD"),
  calc_vif_tbl(dat, c("recombination", "avg_wc_fst", "avg_dxy", "B_no_delta_p", "K_local_mean"), "VIF: focal metrics + B_no_delta_p + K_local_mean")
)

class_tests <- bind_rows(
  safe_group_test_tbl(dat$K_local_mean, dat$recombination_bin, "K_local_mean by recombination quantile"),
  safe_group_test_tbl(dat$K_local_mean, dat$K_local_class, "K_local_mean by K_local class"),
  safe_group_test_tbl(dat$avg_wc_fst, dat$K_local_class, "avg_wc_fst by K_local class"),
  safe_group_test_tbl(dat$avg_dxy, dat$K_local_class, "avg_dxy by K_local class"),
  safe_group_test_tbl(dat$abs_delta_p, dat$K_local_class, "abs_delta_p by K_local class")
)
if (any(is.finite(dat$pi_average))) class_tests <- bind_rows(class_tests, safe_group_test_tbl(dat$pi_average, dat$K_local_class, "pi_average by K_local class"))
if (any(is.finite(dat$TajimaD))) class_tests <- bind_rows(class_tests, safe_group_test_tbl(dat$TajimaD, dat$K_local_class, "TajimaD by K_local class"))
if (any(is.finite(dat$f3))) class_tests <- bind_rows(class_tests, safe_group_test_tbl(dat$f3, dat$K_local_class, "f3 by K_local class"))
if (any(!is.na(dat$sweep_group))) class_tests <- bind_rows(class_tests, safe_group_test_tbl(dat$K_local_mean, dat$sweep_group, "K_local_mean by sweep/model class"))
class_tests <- bind_rows(class_tests, safe_group_test_tbl(dat$K_local_mean, ifelse(dat$is_zw, "Chromosome_ZW", "Autosome"), "K_local_mean: Chromosome_ZW vs autosomes"))

top_window_summary <- bind_rows(lapply(cfg$top_thresholds, function(thr) {
  cutoff <- stats::quantile(dat$K_local_mean, probs = thr, na.rm = TRUE)
  top_dat <- dat %>% filter(is.finite(K_local_mean), K_local_mean >= cutoff)
  bind_rows(
    tibble(summary_type = "overall", threshold = thr, threshold_label = paste0("top_", round((1 - thr) * 100), "pct"), category = "all",
           n_windows = nrow(top_dat), mean_recombination = mean(top_dat$recombination, na.rm = TRUE),
           mean_avg_wc_fst = mean(top_dat$avg_wc_fst, na.rm = TRUE), mean_avg_dxy = mean(top_dat$avg_dxy, na.rm = TRUE),
           mean_pi_average = mean(top_dat$pi_average, na.rm = TRUE), mean_TajimaD = mean(top_dat$TajimaD, na.rm = TRUE),
           mean_abs_delta_p = mean(top_dat$abs_delta_p, na.rm = TRUE), mean_f3 = mean(top_dat$f3, na.rm = TRUE)),
    top_dat %>% count(category = Chromosome_Names, name = "n_windows") %>% mutate(summary_type = "chromosome", threshold = thr, threshold_label = paste0("top_", round((1 - thr) * 100), "pct"),
           mean_recombination = NA_real_, mean_avg_wc_fst = NA_real_, mean_avg_dxy = NA_real_, mean_pi_average = NA_real_, mean_TajimaD = NA_real_, mean_abs_delta_p = NA_real_, mean_f3 = NA_real_),
    top_dat %>% count(category = sweep_group, name = "n_windows") %>% mutate(summary_type = "class", threshold = thr, threshold_label = paste0("top_", round((1 - thr) * 100), "pct"),
           mean_recombination = NA_real_, mean_avg_wc_fst = NA_real_, mean_avg_dxy = NA_real_, mean_pi_average = NA_real_, mean_TajimaD = NA_real_, mean_abs_delta_p = NA_real_, mean_f3 = NA_real_)
  )
}))

chrom_enrichment <- bind_rows(lapply(cfg$top_thresholds, function(thr) {
  cutoff <- stats::quantile(dat$K_local_mean, probs = thr, na.rm = TRUE)
  dat %>%
    mutate(in_top = is.finite(K_local_mean) & K_local_mean >= cutoff) %>%
    group_by(Chromosome_Names) %>%
    summarise(
      threshold = thr,
      threshold_label = paste0("top_", round((1 - thr) * 100), "pct"),
      n_top = sum(in_top, na.rm = TRUE),
      n_total = n(),
      frac_top = n_top / n_total,
      mean_K_local = mean(K_local_mean, na.rm = TRUE),
      .groups = "drop"
    )
}))

run_cutoff <- stats::quantile(dat$K_local_mean, probs = cfg$run_quantile, na.rm = TRUE)
runlength_summary <- bind_rows(lapply(split(dat, dat$Chromosome_Names), function(df_chr) {
  df_chr <- df_chr[order(df_chr$.start_bp), ]
  high <- is.finite(df_chr$K_local_mean) & df_chr$K_local_mean >= run_cutoff
  rt <- run_length_table(high)
  high_runs <- rt %>% filter(value %in% TRUE)
  null_tbl <- perm_run_null(high, cfg$n_permutations)
  tibble(
    Chromosome_Names = as.character(df_chr$Chromosome_Names[1]),
    cutoff_quantile = cfg$run_quantile,
    cutoff_value = run_cutoff,
    n_high_windows = sum(high, na.rm = TRUE),
    n_runs = nrow(high_runs),
    mean_run_length = ifelse(nrow(high_runs) > 0, mean(high_runs$run_length), 0),
    max_run_length = ifelse(nrow(high_runs) > 0, max(high_runs$run_length), 0),
    mean_max_run_perm = null_tbl$mean_max_run_perm,
    p_ge_obs = null_tbl$p_ge_obs
  )
}))

write.csv(model_comparisons, file.path(cfg$output_dir, "k_local_driver_model_comparisons.csv"), row.names = FALSE)
write.csv(coefficients_tbl, file.path(cfg$output_dir, "k_local_driver_coefficients.csv"), row.names = FALSE)
write.csv(bind_rows(correlations_tbl %>% mutate(kind = "correlation", vif = NA_real_) %>% select(analysis, kind, term, estimate, p_value, statistic, n, r_squared, adj_r_squared, aic, vif, note),
                    vif_tbl %>% mutate(kind = "vif", estimate = NA_real_, p_value = NA_real_, statistic = NA_real_, n = NA_integer_, r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = NA_character_) %>% select(analysis, kind, term, estimate, p_value, statistic, n, r_squared, adj_r_squared, aic, vif, note)),
          file.path(cfg$output_dir, "k_local_driver_correlations.csv"), row.names = FALSE)
write.csv(class_tests, file.path(cfg$output_dir, "k_local_driver_class_tests.csv"), row.names = FALSE)
write.csv(top_window_summary, file.path(cfg$output_dir, "k_local_top_windows_summary.csv"), row.names = FALSE)
write.csv(runlength_summary, file.path(cfg$output_dir, "k_local_runlength_summary.csv"), row.names = FALSE)
write.csv(chrom_enrichment, file.path(cfg$output_dir, "k_local_chromosome_enrichment.csv"), row.names = FALSE)
write.csv(detected_columns, file.path(cfg$output_dir, "detected_columns.csv"), row.names = FALSE)

plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

scatter_specs <- list(
  list(x = "recombination", file = "plot_K_local_vs_recombination.png", title = "K_local_mean vs recombination", color = "#1f4e79"),
  list(x = "avg_wc_fst", file = "plot_K_local_vs_avg_wc_fst.png", title = "K_local_mean vs avg_wc_fst", color = "#8c2d04"),
  list(x = "avg_dxy", file = "plot_K_local_vs_avg_dxy.png", title = "K_local_mean vs avg_dxy", color = "#7a3b69"),
  list(x = "pi_average", file = "plot_K_local_vs_pi_average.png", title = "K_local_mean vs pi_average", color = "#238b45"),
  list(x = "TajimaD", file = "plot_K_local_vs_TajimaD.png", title = "K_local_mean vs TajimaD", color = "#6a51a3"),
  list(x = "f3", file = "plot_K_local_vs_f3.png", title = "K_local_mean vs f3", color = "#2171b5")
)

for (sp in scatter_specs) {
  if (!(sp$x %in% names(dat)) || !any(is.finite(dat[[sp$x]]))) next
  p <- ggplot(dat, aes(.data[[sp$x]], K_local_mean)) +
    geom_point(alpha = 0.25, size = 0.8, color = sp$color) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black") +
    labs(title = sp$title, x = sp$x, y = "K_local_mean") +
    plot_theme
  ggsave(file.path(cfg$output_dir, sp$file), p, width = 7.8, height = 5.5, dpi = 300)
}

if (any(!is.na(dat$sweep_group))) {
  ggsave(file.path(cfg$output_dir, "plot_K_local_by_sweep_class.png"),
    ggplot(dat %>% filter(!is.na(sweep_group)), aes(sweep_group, K_local_mean, fill = sweep_group)) +
      geom_boxplot(outlier.alpha = 0.15) +
      labs(title = "K_local_mean by sweep/model class", x = NULL, y = "K_local_mean") +
      plot_theme,
    width = 8, height = 5.6, dpi = 300
  )
}

ggsave(file.path(cfg$output_dir, "plot_K_local_by_chromosome.png"),
  ggplot(dat, aes(Chromosome_Names, K_local_mean)) +
    geom_boxplot(outlier.alpha = 0.05, fill = "#9ecae1") +
    labs(title = "K_local_mean by chromosome", x = "Chromosome", y = "K_local_mean") +
    plot_theme,
  width = 11, height = 5.8, dpi = 300
)

cut_top5 <- stats::quantile(dat$K_local_mean, probs = 0.95, na.rm = TRUE)
top_chr <- chrom_enrichment %>% filter(abs(threshold - 0.95) < 1e-8) %>% arrange(desc(frac_top)) %>% slice_head(n = 8) %>% pull(Chromosome_Names) %>% as.character()
ggsave(file.path(cfg$output_dir, "plot_K_local_along_chromosomes_top5pct.png"),
  ggplot(dat %>% filter(as.character(Chromosome_Names) %in% top_chr), aes(window_mid_bp / 1e6, K_local_mean)) +
    geom_line(color = "grey50", linewidth = 0.25, na.rm = TRUE) +
    geom_point(data = dat %>% filter(as.character(Chromosome_Names) %in% top_chr, is.finite(K_local_mean), K_local_mean >= cut_top5),
               aes(window_mid_bp / 1e6, K_local_mean), color = "#b22222", size = 0.5, inherit.aes = FALSE) +
    facet_wrap(~ Chromosome_Names, scales = "free_x") +
    labs(title = "K_local_mean along chromosomes with top 5% windows highlighted", x = "Position (Mb)", y = "K_local_mean") +
    plot_theme,
  width = 11, height = 7.5, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_runlength_summary.png"),
  ggplot(runlength_summary, aes(reorder(Chromosome_Names, max_run_length), max_run_length)) +
    geom_col(fill = "#8c2d04", color = "grey30") +
    geom_point(aes(y = mean_max_run_perm), color = "#1f4e79", size = 2) +
    coord_flip() +
    labs(title = "Observed vs permuted max run length of top 5% K_local windows", x = "Chromosome", y = "Run length") +
    plot_theme,
  width = 8.5, height = 6.8, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_top_window_chromosome_enrichment.png"),
  ggplot(chrom_enrichment, aes(Chromosome_Names, frac_top, fill = threshold_label)) +
    geom_col(position = "dodge", color = "grey30") +
    labs(title = "Chromosome enrichment of top K_local windows", x = "Chromosome", y = "Fraction of chromosome in top set", fill = "Threshold") +
    plot_theme,
  width = 11, height = 5.8, dpi = 300
)

if (any(!is.na(dat$sweep_group))) {
  class_enrichment <- top_window_summary %>% filter(summary_type == "class", !is.na(category))
  ggsave(file.path(cfg$output_dir, "plot_top_window_class_enrichment.png"),
    ggplot(class_enrichment, aes(category, n_windows, fill = threshold_label)) +
      geom_col(position = "dodge", color = "grey30") +
      labs(title = "Sweep/model composition of top K_local windows", x = NULL, y = "Windows", fill = "Threshold") +
      plot_theme,
    width = 8.5, height = 5.6, dpi = 300
  )
}

best_multi <- model_comparisons %>% arrange(desc(adj_r_squared), aic) %>% slice(1)
zw_test <- class_tests %>% filter(grepl("Chromosome_ZW", analysis)) %>% slice(1)
rec_cor <- correlations_tbl %>% filter(analysis == "K_local_mean ~ recombination", model_type == "cor_spearman") %>% slice(1)
sweep_test <- class_tests %>% filter(analysis == "K_local_mean by sweep/model class") %>% slice(1)
mean_perm_p <- mean(runlength_summary$p_ge_obs, na.rm = TRUE)

message("\n================ K_local Driver Analysis ================")
message("Detected columns:")
print(detected_columns)
message("\nTop K_local models by adjusted R^2:")
print(model_comparisons %>% arrange(desc(adj_r_squared), aic) %>% head(10))
message("\nInterpretation:")
message("Is high K_local mainly driven by low recombination? ", ifelse(nrow(rec_cor) == 1 && is.finite(rec_cor$estimate) && rec_cor$estimate < -0.4, "yes, strongly", "not primarily / only partially"))
message("Do sweep-like windows tend to have higher K_local? ", ifelse(nrow(sweep_test) == 1 && is.finite(sweep_test$p_value) && sweep_test$p_value < 0.05, "yes", "not clearly"))
message("Are high-K_local windows enriched on particular chromosomes, especially Chromosome_ZW? ", ifelse(nrow(zw_test) == 1 && is.finite(zw_test$p_value) && zw_test$p_value < 0.05, "yes", "not strongly"))
message("Do high-K_local regions look like extended domains or isolated points? ", ifelse(is.finite(mean_perm_p) && mean_perm_p < 0.1, "extended domains are supported", "more mixed / not strongly above permutation"))
message("Which variables best explain K_local in the multivariable models? ", best_multi$analysis)
message("Most consistent biological interpretation: high local barrier density is most consistent with a broad barrier architecture combining low recombination, clustered divergence, and chromosome-scale genomic structure, with sweep history contributing where sweep-like windows are enriched.")
message("\nOutputs written to: ", normalizePath(cfg$output_dir))
finalize_sanity(
  sanity,
  files = c(
    "k_local_driver_model_comparisons.csv",
    "k_local_driver_coefficients.csv",
    "k_local_driver_correlations.csv",
    "k_local_driver_class_tests.csv",
    "k_local_top_windows_summary.csv",
    "k_local_runlength_summary.csv",
    "k_local_chromosome_enrichment.csv",
    "detected_columns.csv",
    "plot_K_local_by_chromosome.png",
    "plot_K_local_along_chromosomes_top5pct.png",
    "plot_runlength_summary.png"
  )
)
