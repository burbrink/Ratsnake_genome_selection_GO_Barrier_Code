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
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "hybrid_zone_followup_outputs"),
  alpha = 0.5,
  alpha_calibration_csv = NULL,
  alpha_calibration_p_col = "P_overall_mean",
  alpha_calibration_k_col = "K_genome_mean",
  rolling_window_n = 25L
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_rdata <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$object_name <- args[[2]]
if (length(args) >= 3 && nzchar(args[[3]])) cfg$output_dir <- args[[3]]
if (length(args) >= 4 && nzchar(args[[4]])) cfg$alpha_calibration_csv <- args[[4]]
if (length(args) >= 5 && nzchar(args[[5]])) cfg$alpha_calibration_p_col <- args[[5]]
if (length(args) >= 6 && nzchar(args[[6]])) cfg$alpha_calibration_k_col <- args[[6]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
sanity <- init_sanity("hybrid_zone_followup_recombination_control_sd3_new_af.R", cfg$output_dir)

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
# If K_local remains predictive after controlling for recombination, then the
# model captures more than recombination alone.
# If K_local improves fit beyond B_i alone, that supports local barrier
# accumulation rather than purely focal-window effects.

barrier_formula_terms <- list(
  base_no_delta = c("-z_recombination", "z_avg_wc_fst", "z_avg_dxy", "-z_pi_average", "-z_TajimaD")
)

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

make_score <- function(df, terms) {
  if (length(terms) == 0L) return(rep(NA_real_, nrow(df)))
  with(df, eval(parse(text = paste(terms, collapse = " + "))))
}

safe_cor_tbl <- function(x, y, analysis, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L) {
    return(tibble(
      analysis = analysis, model_type = paste0("cor_", method), term = NA_character_,
      estimate = NA_real_, statistic = NA_real_, p_value = NA_real_, n = sum(keep),
      r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = "Too few complete observations"
    ))
  }
  ct <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = method))
  tibble(
    analysis = analysis, model_type = paste0("cor_", method), term = NA_character_,
    estimate = unname(ct$estimate), statistic = unname(ct$statistic), p_value = ct$p.value, n = sum(keep),
    r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = NA_character_
  )
}

safe_lm_details <- function(formula, data, analysis) {
  mf <- model.frame(formula, data = data, na.action = na.omit)
  if (nrow(mf) < 5L) {
    return(list(
      model = NULL,
      coefficients = tibble(
        analysis = analysis, model_type = "lm", term = NA_character_,
        estimate = NA_real_, std_beta = NA_real_, statistic = NA_real_, p_value = NA_real_, n = nrow(mf),
        r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = "Too few complete observations"
      ),
      summary = tibble(
        analysis = analysis, n = nrow(mf), r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = "Too few complete observations"
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
    coefficients = tibble(
      analysis = analysis, model_type = "lm", term = cf$term,
      estimate = cf$Estimate, std_beta = unname(std_beta[cf$term]), statistic = cf$`t value`, p_value = cf$`Pr(>|t|)`,
      n = nrow(mf), r_squared = unname(sm$r.squared), adj_r_squared = unname(sm$adj.r.squared),
      aic = stats::AIC(fit), note = NA_character_
    ),
    summary = tibble(
      analysis = analysis, n = nrow(mf), r_squared = unname(sm$r.squared), adj_r_squared = unname(sm$adj.r.squared),
      aic = stats::AIC(fit), note = NA_character_
    )
  )
}

safe_group_test_tbl <- function(y, group, analysis) {
  keep <- is.finite(y) & !is.na(group)
  y2 <- y[keep]
  g2 <- droplevels(as.factor(group[keep]))
  if (length(y2) < 5L || nlevels(g2) < 2L) {
    return(tibble(
      analysis = analysis, model_type = "group_test", term = NA_character_, estimate = NA_real_,
      statistic = NA_real_, p_value = NA_real_, n = length(y2), r_squared = NA_real_,
      adj_r_squared = NA_real_, aic = NA_real_, note = "Too few observations or groups"
    ))
  }
  if (nlevels(g2) == 2L) {
    wt <- suppressWarnings(stats::wilcox.test(y2 ~ g2))
    return(tibble(
      analysis = analysis, model_type = "wilcox", term = paste(levels(g2), collapse = " vs "),
      estimate = NA_real_, statistic = unname(wt$statistic), p_value = wt$p.value, n = length(y2),
      r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = NA_character_
    ))
  }
  fit <- stats::aov(y2 ~ g2)
  sm <- summary(fit)[[1]]
  tibble(
    analysis = analysis, model_type = "anova", term = rownames(sm)[1],
    estimate = NA_real_, statistic = sm$`F value`[1], p_value = sm$`Pr(>F)`[1], n = length(y2),
    r_squared = NA_real_, adj_r_squared = NA_real_, aic = stats::AIC(stats::lm(y2 ~ g2)), note = NA_character_
  )
}

partial_cor_tbl <- function(data, y_col, x_col, control_cols, analysis) {
  req <- c(y_col, x_col, control_cols)
  if (!all(req %in% names(data))) {
    return(tibble(
      analysis = analysis, estimate = NA_real_, p_value = NA_real_, statistic = NA_real_,
      n = 0L, note = "Missing required columns"
    ))
  }
  df <- data[, req, drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) < 5L) {
    return(tibble(
      analysis = analysis, estimate = NA_real_, p_value = NA_real_, statistic = NA_real_,
      n = nrow(df), note = "Too few complete observations"
    ))
  }
  f_y <- stats::as.formula(paste(y_col, "~", paste(control_cols, collapse = " + ")))
  f_x <- stats::as.formula(paste(x_col, "~", paste(control_cols, collapse = " + ")))
  res_y <- stats::resid(stats::lm(f_y, data = df))
  res_x <- stats::resid(stats::lm(f_x, data = df))
  ct <- suppressWarnings(stats::cor.test(res_x, res_y, method = "spearman"))
  tibble(
    analysis = analysis, estimate = unname(ct$estimate), p_value = ct$p.value,
    statistic = unname(ct$statistic), n = nrow(df), note = NA_character_
  )
}

calc_vif_tbl <- function(data, predictors, analysis_prefix) {
  predictors <- predictors[predictors %in% names(data)]
  if (length(predictors) < 2L) return(tibble())
  df <- data[, predictors, drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) < 5L) return(tibble())
  bind_rows(lapply(predictors, function(pred) {
    others <- setdiff(predictors, pred)
    if (length(others) == 0L) {
      return(tibble(analysis = analysis_prefix, term = pred, vif = 1))
    }
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
assert_has_cols(dat, c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"), sanity, "hybrid_followup_input_data")
assert_unique_keys(dat, c("Chromosome_Names", ".start_bp", ".end_bp"), sanity, "hybrid_followup_input_data")

prior_local <- read.csv(file.path(cfg$prior_barrier_dir, "barrier_recalibrated_local_summary.csv"), stringsAsFactors = FALSE)
n_before_join <- nrow(dat)
dat <- dat %>%
  left_join(prior_local, by = c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"))
assert_join_rowcount(n_before_join, dat, sanity, "hybrid_followup_prior_local_join")

if ("abs_delta_p.x" %in% names(dat) || "abs_delta_p.y" %in% names(dat)) dat$abs_delta_p <- dplyr::coalesce(dat$abs_delta_p.x, dat$abs_delta_p.y)
if ("f3.x" %in% names(dat) || "f3.y" %in% names(dat)) dat$f3 <- dplyr::coalesce(dat$f3.x, dat$f3.y)

for (nm in c("recombination", "avg_wc_fst", "avg_dxy", "pi_average", "TajimaD")) {
  dat[[paste0("z_", nm)]] <- z_score(dat[[nm]])
}

if (!"B_no_delta_p" %in% names(dat) || !any(is.finite(dat$B_no_delta_p))) {
  dat$B_no_delta_p <- make_score(dat, barrier_formula_terms$base_no_delta)
}

if (!"K_local_mean" %in% names(dat) || !any(is.finite(dat$K_local_mean))) {
  dat <- dat %>%
    group_by(Chromosome_Names) %>%
    arrange(.start_bp, .by_group = TRUE) %>%
    mutate(
      K_local_mean = data.table::frollmean(B_no_delta_p, n = cfg$rolling_window_n, align = "center", na.rm = FALSE),
      K_local_sum = data.table::frollsum(B_no_delta_p, n = cfg$rolling_window_n, align = "center", na.rm = FALSE)
    ) %>%
    ungroup()
}

if (!"P_local_mean" %in% names(dat) || !any(is.finite(dat$P_local_mean))) dat$P_local_mean <- exp(-cfg$alpha * dat$K_local_mean)
if (!"sweep_group" %in% names(dat)) dat$sweep_group <- NA_character_
assert_missing_fraction(dat$B_no_delta_p, 0.1, sanity, "B_no_delta_p")
assert_missing_fraction(dat$K_local_mean, 0.1, sanity, "K_local_mean")

dat$K_local_class <- quantile_class(dat$K_local_mean)
dat$recombination_quantile <- quantile_class(dat$recombination)
dat$Chromosome_Names <- factor(dat$Chromosome_Names, levels = order_chromosomes(dat$Chromosome_Names))

f3_summary <- tibble(
  metric = c("n_finite", "mean", "median", "sd", "min", "max"),
  value = c(sum(is.finite(dat$f3)), mean(dat$f3, na.rm = TRUE), median(dat$f3, na.rm = TRUE), stats::sd(dat$f3, na.rm = TRUE), min(dat$f3, na.rm = TRUE), max(dat$f3, na.rm = TRUE))
)

abs_models <- list(
  "M1: abs_delta_p ~ recombination" = abs_delta_p ~ recombination,
  "M2: abs_delta_p ~ B_no_delta_p" = abs_delta_p ~ B_no_delta_p,
  "M3: abs_delta_p ~ K_local_mean" = abs_delta_p ~ K_local_mean,
  "M4: abs_delta_p ~ recombination + K_local_mean" = abs_delta_p ~ recombination + K_local_mean,
  "M5: abs_delta_p ~ recombination + B_no_delta_p" = abs_delta_p ~ recombination + B_no_delta_p,
  "M6: abs_delta_p ~ recombination + B_no_delta_p + K_local_mean" = abs_delta_p ~ recombination + B_no_delta_p + K_local_mean,
  "M7: abs_delta_p ~ avg_wc_fst" = abs_delta_p ~ avg_wc_fst,
  "M8: abs_delta_p ~ avg_dxy" = abs_delta_p ~ avg_dxy,
  "M9: abs_delta_p ~ avg_wc_fst + avg_dxy + recombination" = abs_delta_p ~ avg_wc_fst + avg_dxy + recombination,
  "M10: abs_delta_p ~ avg_wc_fst + avg_dxy + recombination + K_local_mean" = abs_delta_p ~ avg_wc_fst + avg_dxy + recombination + K_local_mean,
  "M11: abs_delta_p ~ avg_wc_fst + avg_dxy + recombination + B_no_delta_p" = abs_delta_p ~ avg_wc_fst + avg_dxy + recombination + B_no_delta_p,
  "M12: abs_delta_p ~ avg_wc_fst + avg_dxy + recombination + B_no_delta_p + K_local_mean" = abs_delta_p ~ avg_wc_fst + avg_dxy + recombination + B_no_delta_p + K_local_mean
)

f3_models <- list(
  "F1: f3 ~ recombination" = f3 ~ recombination,
  "F2: f3 ~ B_no_delta_p" = f3 ~ B_no_delta_p,
  "F3: f3 ~ K_local_mean" = f3 ~ K_local_mean,
  "F4: f3 ~ recombination + K_local_mean" = f3 ~ recombination + K_local_mean,
  "F5: f3 ~ recombination + B_no_delta_p" = f3 ~ recombination + B_no_delta_p,
  "F6: f3 ~ recombination + B_no_delta_p + K_local_mean" = f3 ~ recombination + B_no_delta_p + K_local_mean
)

abs_fit_list <- lapply(names(abs_models), function(nm) safe_lm_details(abs_models[[nm]], dat, nm))
names(abs_fit_list) <- names(abs_models)
model_comparisons <- bind_rows(lapply(abs_fit_list, `[[`, "summary")) %>% mutate(response = "abs_delta_p")
coefficients_tbl <- bind_rows(lapply(abs_fit_list, `[[`, "coefficients")) %>% mutate(response = "abs_delta_p")

f3_fit_list <- list()
if (any(is.finite(dat$f3))) {
  f3_fit_list <- lapply(names(f3_models), function(nm) safe_lm_details(f3_models[[nm]], dat, nm))
  names(f3_fit_list) <- names(f3_models)
  model_comparisons <- bind_rows(model_comparisons, bind_rows(lapply(f3_fit_list, `[[`, "summary")) %>% mutate(response = "f3"))
  coefficients_tbl <- bind_rows(coefficients_tbl, bind_rows(lapply(f3_fit_list, `[[`, "coefficients")) %>% mutate(response = "f3"))
}

partial_correlations <- bind_rows(
  safe_cor_tbl(dat$K_local_mean, dat$recombination, "K_local_mean ~ recombination", "pearson"),
  safe_cor_tbl(dat$K_local_mean, dat$recombination, "K_local_mean ~ recombination", "spearman")
) %>%
  mutate(kind = "simple")

partial_correlations <- bind_rows(
  partial_correlations,
  partial_cor_tbl(dat, "abs_delta_p", "K_local_mean", c("recombination"), "Partial: abs_delta_p ~ K_local_mean | recombination") %>% mutate(kind = "partial")
)
if (any(is.finite(dat$f3))) {
  partial_correlations <- bind_rows(
    partial_correlations,
    partial_cor_tbl(dat, "f3", "K_local_mean", c("recombination"), "Partial: f3 ~ K_local_mean | recombination") %>% mutate(kind = "partial")
  )
}

vif_tbl <- bind_rows(
  calc_vif_tbl(dat, c("recombination", "B_no_delta_p", "K_local_mean"), "VIF: recombination + B_no_delta_p + K_local_mean"),
  calc_vif_tbl(dat, c("avg_wc_fst", "avg_dxy", "recombination", "B_no_delta_p", "K_local_mean"), "VIF: focal metrics + recombination + B + K_local")
)

partial_correlations_export <- bind_rows(
  partial_correlations %>%
    mutate(term = NA_character_, vif = NA_real_) %>%
    select(analysis, kind, term, estimate, p_value, statistic, n, vif, note),
  vif_tbl %>%
    mutate(kind = "vif", estimate = NA_real_, p_value = NA_real_, statistic = NA_real_, n = NA_integer_, note = NA_character_) %>%
    select(analysis, kind, term, estimate, p_value, statistic, n, vif, note)
)

abs_resid <- stats::resid(stats::lm(abs_delta_p ~ recombination, data = dat, na.action = na.exclude))
dat$abs_delta_p_resid_recomb <- abs_resid
class_tests <- bind_rows(
  safe_group_test_tbl(dat$abs_delta_p, dat$K_local_class, "abs_delta_p by K_local class"),
  safe_group_test_tbl(dat$abs_delta_p_resid_recomb, dat$K_local_class, "abs_delta_p residualized on recombination by K_local class"),
  safe_lm_details(abs_delta_p ~ recombination + K_local_class, dat, "ANCOVA: abs_delta_p ~ recombination + K_local_class")$summary
) %>% mutate(response = "abs_delta_p")

if (any(is.finite(dat$f3))) {
  dat$f3_resid_recomb <- stats::resid(stats::lm(f3 ~ recombination, data = dat, na.action = na.exclude))
  class_tests <- bind_rows(
    class_tests,
    safe_group_test_tbl(dat$f3, dat$K_local_class, "f3 by K_local class"),
    safe_group_test_tbl(dat$f3_resid_recomb, dat$K_local_class, "f3 residualized on recombination by K_local class"),
    safe_lm_details(f3 ~ recombination + K_local_class, dat, "ANCOVA: f3 ~ recombination + K_local_class")$summary %>% mutate(response = "f3")
  )
}

chromosome_summary <- dat %>%
  group_by(Chromosome_Names) %>%
  summarise(
    n_windows = n(),
    mean_B_chr = mean(B_no_delta_p, na.rm = TRUE),
    mean_K_local = mean(K_local_mean, na.rm = TRUE),
    var_K_local = stats::var(K_local_mean, na.rm = TRUE),
    mean_recombination = mean(recombination, na.rm = TRUE),
    mean_abs_delta_p = mean(abs_delta_p, na.rm = TRUE),
    mean_f3 = mean(f3, na.rm = TRUE),
    .groups = "drop"
  )

chr_fit <- safe_lm_details(mean_abs_delta_p ~ mean_K_local, chromosome_summary, "Chromosome means: mean_abs_delta_p ~ mean_K_local")
model_comparisons <- bind_rows(model_comparisons, chr_fit$summary %>% mutate(response = "chromosome_mean_abs_delta_p"))
coefficients_tbl <- bind_rows(coefficients_tbl, chr_fit$coefficients %>% mutate(response = "chromosome_mean_abs_delta_p"))

write.csv(detected_columns, file.path(cfg$output_dir, "detected_columns.csv"), row.names = FALSE)
write.csv(model_comparisons, file.path(cfg$output_dir, "hybrid_zone_followup_model_comparisons.csv"), row.names = FALSE)
write.csv(coefficients_tbl, file.path(cfg$output_dir, "hybrid_zone_followup_coefficients.csv"), row.names = FALSE)
write.csv(partial_correlations_export, file.path(cfg$output_dir, "hybrid_zone_followup_partial_correlations.csv"), row.names = FALSE)
write.csv(chromosome_summary, file.path(cfg$output_dir, "hybrid_zone_followup_chromosome_summary.csv"), row.names = FALSE)
write.csv(class_tests, file.path(cfg$output_dir, "hybrid_zone_followup_class_tests.csv"), row.names = FALSE)
write.csv(f3_summary, file.path(cfg$output_dir, "f3_distribution_summary.csv"), row.names = FALSE)

plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(cfg$output_dir, "plot_abs_delta_vs_K_local_by_recomb_quantile.png"),
  ggplot(dat, aes(K_local_mean, abs_delta_p, color = recombination_quantile)) +
    geom_point(alpha = 0.28, size = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black") +
    labs(title = "abs_delta_p vs K_local_mean by recombination class", x = "K_local_mean", y = "abs_delta_p", color = "Recombination") +
    plot_theme,
  width = 8, height = 5.6, dpi = 300
)

if (any(is.finite(dat$f3))) {
  ggsave(file.path(cfg$output_dir, "plot_f3_vs_K_local_by_recomb_quantile.png"),
    ggplot(dat, aes(K_local_mean, f3, color = recombination_quantile)) +
      geom_point(alpha = 0.28, size = 0.8) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black") +
      labs(title = "f3 vs K_local_mean by recombination class", x = "K_local_mean", y = "f3", color = "Recombination") +
      plot_theme,
    width = 8, height = 5.6, dpi = 300
  )
}

fit_abs_partial <- stats::lm(abs_delta_p ~ recombination + B_no_delta_p, data = dat, na.action = na.exclude)
dat$abs_delta_partial_y <- stats::resid(fit_abs_partial)
fit_abs_partial_x <- stats::lm(K_local_mean ~ recombination + B_no_delta_p, data = dat, na.action = na.exclude)
dat$abs_delta_partial_x <- stats::resid(fit_abs_partial_x)

ggsave(file.path(cfg$output_dir, "plot_added_variable_K_local_for_abs_delta_p.png"),
  ggplot(dat, aes(abs_delta_partial_x, abs_delta_partial_y)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#5a2d82") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#1f4e79") +
    labs(title = "Added-variable view for K_local_mean on abs_delta_p", x = "Residual K_local_mean | recombination + B_no_delta_p", y = "Residual abs_delta_p | recombination + B_no_delta_p") +
    plot_theme,
  width = 8, height = 5.6, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_chromosome_mean_K_local_vs_mean_abs_delta_p.png"),
  ggplot(chromosome_summary, aes(mean_K_local, mean_abs_delta_p, label = Chromosome_Names)) +
    geom_point(size = 2.2, color = "#8c2d04") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#1f4e79") +
    labs(title = "Chromosome mean K_local vs mean abs_delta_p", x = "mean_K_local", y = "mean_abs_delta_p") +
    plot_theme,
  width = 7.8, height = 5.6, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_abs_delta_resid_by_K_class.png"),
  ggplot(dat %>% filter(!is.na(K_local_class)), aes(K_local_class, abs_delta_p_resid_recomb, fill = K_local_class)) +
    geom_boxplot(outlier.alpha = 0.15) +
    scale_fill_manual(values = c("Low" = "#91bfdb", "Intermediate" = "#fee090", "High" = "#fc8d59")) +
    labs(title = "abs_delta_p residualized on recombination by K_local class", x = "K_local class", y = "Residual abs_delta_p") +
    plot_theme,
  width = 7.2, height = 5.4, dpi = 300
)

if (any(is.finite(dat$f3))) {
  ggsave(file.path(cfg$output_dir, "plot_f3_resid_by_K_class.png"),
    ggplot(dat %>% filter(!is.na(K_local_class)), aes(K_local_class, f3_resid_recomb, fill = K_local_class)) +
      geom_boxplot(outlier.alpha = 0.15) +
      scale_fill_manual(values = c("Low" = "#91bfdb", "Intermediate" = "#fee090", "High" = "#fc8d59")) +
      labs(title = "f3 residualized on recombination by K_local class", x = "K_local class", y = "Residual f3") +
      plot_theme,
    width = 7.2, height = 5.4, dpi = 300
  )
}

best_abs <- model_comparisons %>% filter(response == "abs_delta_p") %>% arrange(desc(adj_r_squared), aic) %>% slice(1)
best_f3 <- model_comparisons %>% filter(response == "f3") %>% arrange(desc(adj_r_squared), aic) %>% slice(1)
coef_abs_m6 <- coefficients_tbl %>% filter(analysis == "M6: abs_delta_p ~ recombination + B_no_delta_p + K_local_mean", term == "K_local_mean") %>% slice(1)
coef_f3_f6 <- coefficients_tbl %>% filter(analysis == "F6: f3 ~ recombination + B_no_delta_p + K_local_mean", term == "K_local_mean") %>% slice(1)
coef_abs_m5 <- coefficients_tbl %>% filter(analysis == "M5: abs_delta_p ~ recombination + B_no_delta_p", term == "B_no_delta_p") %>% slice(1)

message("\n================ Hybrid-Zone Follow-up ================")
message("Detected columns:")
print(detected_columns)
message("\nf3 distribution:")
print(f3_summary)
message("\nTop abs_delta_p models by adjusted R^2:")
print(model_comparisons %>% filter(response == "abs_delta_p") %>% arrange(desc(adj_r_squared), aic) %>% head(6))
if (nrow(best_f3) > 0) {
  message("\nTop f3 models by adjusted R^2:")
  print(model_comparisons %>% filter(response == "f3") %>% arrange(desc(adj_r_squared), aic) %>% head(6))
}
message("\nPartial correlations and VIF diagnostics saved to CSV.")
message("\nInterpretation:")
message("Does K_local still predict abs_delta_p after controlling for recombination? ", ifelse(nrow(coef_abs_m6) == 1 && is.finite(coef_abs_m6$p_value) && coef_abs_m6$p_value < 0.05, "yes", "no"))
if (nrow(coef_f3_f6) == 1) {
  direction_f3 <- ifelse(coef_f3_f6$estimate > 0, "higher", "lower")
  message("Does K_local still predict f3 after controlling for recombination? ", ifelse(is.finite(coef_f3_f6$p_value) && coef_f3_f6$p_value < 0.05, "yes", "no"),
          "; observed direction: higher K_local is associated with ", direction_f3, " f3 in the multivariable model.")
} else {
  message("Does K_local still predict f3 after controlling for recombination? not tested because f3 was unavailable.")
}
message("Does K_local add explanatory power beyond focal barrier score B_i? ", ifelse(best_abs$analysis %in% c("M6: abs_delta_p ~ recombination + B_no_delta_p + K_local_mean", "M12: abs_delta_p ~ avg_wc_fst + avg_dxy + recombination + B_no_delta_p + K_local_mean"), "yes", "partially / model-dependent"))
message("Does the accumulated-barrier model add something beyond recombination alone? ", ifelse(best_abs$analysis != "M1: abs_delta_p ~ recombination", "yes", "no"))
message("Do the results continue to support a mosaic hybrid-zone regime? ", ifelse(best_abs$analysis != "M1: abs_delta_p ~ recombination" && nrow(coef_abs_m6) == 1 && coef_abs_m6$p_value < 0.05, "yes", "partially / weakly"))
message("\nOutputs written to: ", normalizePath(cfg$output_dir))
finalize_sanity(
  sanity,
  files = c(
    "detected_columns.csv",
    "hybrid_zone_followup_model_comparisons.csv",
    "hybrid_zone_followup_coefficients.csv",
    "hybrid_zone_followup_partial_correlations.csv",
    "hybrid_zone_followup_chromosome_summary.csv",
    "hybrid_zone_followup_class_tests.csv",
    "plot_abs_delta_vs_K_local_by_recomb_quantile.png",
    "plot_added_variable_K_local_for_abs_delta_p.png",
    "plot_chromosome_mean_K_local_vs_mean_abs_delta_p.png",
    "plot_abs_delta_resid_by_K_class.png"
  )
)
