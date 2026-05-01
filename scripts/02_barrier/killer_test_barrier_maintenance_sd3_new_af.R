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
  prior_domain_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_domain_outputs"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "killer_test_outputs"),
  n_bins = 4L,
  primary_high_quantile = 0.95,
  edge_window_span = 2L,
  nearby_outside_span = 2L
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_rdata <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$object_name <- args[[2]]
if (length(args) >= 3 && nzchar(args[[3]])) cfg$output_dir <- args[[3]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
sanity <- init_sanity("killer_test_barrier_maintenance_sd3_new_af.R", cfg$output_dir)

# Biological framing:
# B_i is the focal barrier score of a single window, while K_local is the
# accumulated local barrier density across neighboring windows. If K_local keeps
# improving fit after focal divergence, recombination, and chromosome effects
# are already included, and if that gain generalizes out of sample, then local
# barrier architecture is doing more than simply describing existing patterns.

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

qbin <- function(x, n = 4L) {
  if (all(!is.finite(x))) return(rep(NA_character_, length(x)))
  probs <- seq(0, 1, length.out = n + 1L)
  qs <- unique(stats::quantile(x, probs = probs, na.rm = TRUE, type = 7))
  if (length(qs) <= 2L) return(rep("Q1", length(x)))
  cut(x, breaks = qs, include.lowest = TRUE, ordered_result = TRUE, labels = paste0("Q", seq_len(length(qs) - 1L)))
}

rmse <- function(obs, pred) {
  keep <- is.finite(obs) & is.finite(pred)
  if (!sum(keep)) return(NA_real_)
  sqrt(mean((obs[keep] - pred[keep])^2))
}

test_r2 <- function(obs, pred) {
  keep <- is.finite(obs) & is.finite(pred)
  if (sum(keep) < 3L) return(NA_real_)
  1 - sum((obs[keep] - pred[keep])^2) / sum((obs[keep] - mean(obs[keep]))^2)
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(c(estimate = NA_real_, p_value = NA_real_))
  ct <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = "spearman"))
  c(estimate = unname(ct$estimate), p_value = ct$p.value)
}

dominant_label <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

build_formula <- function(response, include_k = FALSE, include_chr = TRUE, dat) {
  terms <- character()
  if ("recombination" %in% names(dat) && any(is.finite(dat$recombination))) terms <- c(terms, "recombination")
  if ("avg_wc_fst" %in% names(dat) && any(is.finite(dat$avg_wc_fst))) terms <- c(terms, "avg_wc_fst")
  if ("avg_dxy" %in% names(dat) && any(is.finite(dat$avg_dxy))) terms <- c(terms, "avg_dxy")
  if ("pi_average" %in% names(dat) && any(is.finite(dat$pi_average))) terms <- c(terms, "pi_average")
  if ("TajimaD" %in% names(dat) && any(is.finite(dat$TajimaD))) terms <- c(terms, "TajimaD")
  if ("sweep_group" %in% names(dat) && any(!is.na(dat$sweep_group))) terms <- c(terms, "sweep_group")
  if (include_chr && "Chromosome_Names" %in% names(dat)) terms <- c(terms, "Chromosome_Names")
  if (include_k && "K_local_mean" %in% names(dat) && any(is.finite(dat$K_local_mean))) terms <- c(terms, "K_local_mean")
  stats::as.formula(paste(response, "~", paste(unique(terms), collapse = " + ")))
}

fit_model_summary <- function(dat, response, include_k = FALSE, include_chr = TRUE, model_name) {
  form <- build_formula(response, include_k = include_k, include_chr = include_chr, dat = dat)
  vars_needed <- all.vars(form)
  plot_vars <- unique(c(vars_needed, if ("K_local_mean" %in% names(dat)) "K_local_mean"))
  df_plot <- dat[, plot_vars, drop = FALSE]
  df <- dat[, vars_needed, drop = FALSE]
  keep <- stats::complete.cases(df)
  df <- df[keep, , drop = FALSE]
  df_plot <- df_plot[keep, , drop = FALSE]
  if (nrow(df) < 50L) {
    return(list(
      comparison = tibble(
        outcome = response, model_name = model_name, n = nrow(df), formula = paste(deparse(form), collapse = " "),
        r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_,
        k_local_estimate = NA_real_, k_local_p_value = NA_real_, delta_r2 = NA_real_,
        note = "Too few complete cases"
      ),
      coefficients = tibble(),
      model = NULL,
      data = df_plot
    ))
  }
  fit <- stats::lm(form, data = df)
  sm <- summary(fit)
  coef_tbl <- as.data.frame(sm$coefficients)
  coef_tbl$term <- rownames(coef_tbl)
  rownames(coef_tbl) <- NULL
  names(coef_tbl)[1:4] <- c("estimate", "std_error", "statistic", "p_value")
  k_row <- coef_tbl %>% filter(term == "K_local_mean")
  base_r2 <- NA_real_
  if (include_k) {
    base_fit <- stats::lm(build_formula(response, include_k = FALSE, include_chr = include_chr, dat = df), data = df)
    base_r2 <- summary(base_fit)$r.squared
  }
  list(
    comparison = tibble(
      outcome = response,
      model_name = model_name,
      n = nrow(df),
      formula = paste(deparse(form), collapse = " "),
      r_squared = sm$r.squared,
      adj_r_squared = sm$adj.r.squared,
      aic = AIC(fit),
      k_local_estimate = if (nrow(k_row)) k_row$estimate[[1]] else NA_real_,
      k_local_p_value = if (nrow(k_row)) k_row$p_value[[1]] else NA_real_,
      delta_r2 = if (include_k) sm$r.squared - base_r2 else NA_real_,
      note = NA_character_
    ),
    coefficients = coef_tbl %>% mutate(outcome = response, model_name = model_name, formula = paste(deparse(form), collapse = " ")),
    model = fit,
    data = df_plot
  )
}

leave_one_chr_out <- function(dat, response) {
  chrs <- unique(as.character(dat$Chromosome_Names))
  bind_rows(lapply(chrs, function(chr_holdout) {
    train <- dat %>% filter(as.character(Chromosome_Names) != chr_holdout)
    test <- dat %>% filter(as.character(Chromosome_Names) == chr_holdout)
    form_base <- build_formula(response, include_k = FALSE, include_chr = FALSE, dat = train)
    form_k <- build_formula(response, include_k = TRUE, include_chr = FALSE, dat = train)
    base_vars <- all.vars(form_base)
    k_vars <- all.vars(form_k)
    train_base <- train[stats::complete.cases(train[, base_vars, drop = FALSE]), base_vars, drop = FALSE]
    test_base <- test[stats::complete.cases(test[, base_vars, drop = FALSE]), base_vars, drop = FALSE]
    train_k <- train[stats::complete.cases(train[, k_vars, drop = FALSE]), k_vars, drop = FALSE]
    test_k <- test[stats::complete.cases(test[, k_vars, drop = FALSE]), k_vars, drop = FALSE]
    out <- list()
    if (nrow(train_base) >= 50L && nrow(test_base) >= 20L) {
      fit_b <- stats::lm(form_base, data = train_base)
      pred_b <- stats::predict(fit_b, newdata = test_base)
      sp_b <- safe_spearman(test_base[[response]], pred_b)
      out[[1]] <- tibble(
        outcome = response, held_out_chromosome = chr_holdout, model_name = "baseline",
        n_train = nrow(train_base), n_test = nrow(test_base),
        test_r_squared = test_r2(test_base[[response]], pred_b),
        rmse = rmse(test_base[[response]], pred_b),
        spearman = sp_b[["estimate"]], spearman_p_value = sp_b[["p_value"]]
      )
    }
    if (nrow(train_k) >= 50L && nrow(test_k) >= 20L) {
      fit_k <- stats::lm(form_k, data = train_k)
      pred_k <- stats::predict(fit_k, newdata = test_k)
      sp_k <- safe_spearman(test_k[[response]], pred_k)
      out[[2]] <- tibble(
        outcome = response, held_out_chromosome = chr_holdout, model_name = "barrier",
        n_train = nrow(train_k), n_test = nrow(test_k),
        test_r_squared = test_r2(test_k[[response]], pred_k),
        rmse = rmse(test_k[[response]], pred_k),
        spearman = sp_k[["estimate"]], spearman_p_value = sp_k[["p_value"]]
      )
    }
    bind_rows(out)
  }))
}

run_table <- function(flag) {
  r <- rle(flag)
  idx <- cumsum(r$lengths)
  tibble(value = r$values, run_length = r$lengths, end_idx = idx, start_idx = idx - r$lengths + 1L)
}

classify_domain_context <- function(df_chr, high_flag, edge_span = 2L, outside_span = 2L) {
  n <- nrow(df_chr)
  cls <- rep(NA_character_, n)
  rt <- run_table(high_flag)
  high_runs <- rt %>% filter(value %in% TRUE)
  if (!nrow(high_runs)) return(cls)
  for (i in seq_len(nrow(high_runs))) {
    s <- high_runs$start_idx[i]
    e <- high_runs$end_idx[i]
    inside_idx <- s:e
    edge_idx <- unique(c(seq.int(s, min(e, s + edge_span - 1L)), seq.int(max(s, e - edge_span + 1L), e)))
    edge_idx <- edge_idx[edge_idx >= s & edge_idx <= e]
    core_idx <- setdiff(inside_idx, edge_idx)
    outside_idx <- unique(c(seq.int(max(1L, s - outside_span), s - 1L), seq.int(e + 1L, min(n, e + outside_span))))
    outside_idx <- outside_idx[!high_flag[outside_idx]]
    cls[core_idx] <- "Inside"
    cls[edge_idx] <- "Edge"
    cls[outside_idx[is.na(cls[outside_idx])]] <- "Nearby outside"
  }
  cls
}

safe_group_test <- function(y, group, analysis) {
  keep <- is.finite(y) & !is.na(group)
  y2 <- y[keep]
  g2 <- droplevels(as.factor(group[keep]))
  if (length(y2) < 5L || nlevels(g2) < 2L) {
    return(tibble(analysis = analysis, n = length(y2), statistic = NA_real_, p_value = NA_real_, effect = NA_real_, note = "Too few complete cases"))
  }
  fit <- stats::aov(y2 ~ g2)
  sm <- summary(fit)[[1]]
  tibble(
    analysis = analysis,
    n = length(y2),
    statistic = sm$`F value`[1],
    p_value = sm$`Pr(>F)`[1],
    effect = stats::sd(tapply(y2, g2, mean, na.rm = TRUE), na.rm = TRUE),
    note = NA_character_
  )
}

raw <- local({
  load(cfg$input_rdata)
  if (!exists(cfg$object_name, inherits = FALSE)) stop("Object not found: ", cfg$object_name)
  get(cfg$object_name, inherits = FALSE)
})
assert_true(nrow(as.data.frame(raw)) > 0L, "Loaded dataset is empty.", sanity, "nonempty_input_dataset")

chr_col <- find_col(raw, c("Chromosome_Names", "Chromosome", "chr", "chromosome"))
start_col <- find_col(raw, c(".start_bp", "window_start", "start_bp", "start"))
end_col <- find_col(raw, c(".end_bp", "window_end", "end_bp", "end"))
mid_col <- find_col(raw, c("midpoint", ".mid_bp", "window_mid", "mid_bp", "mid"))
recomb_col <- find_col(raw, c("recombination", "rec", "rho"))
fst_col <- find_col(raw, c("avg_wc_fst", "wc_fst", "fst"))
dxy_col <- find_col(raw, c("avg_dxy", "dxy"))
pi_col <- find_col(raw, c("pi_average", "avg_pi", "pi"))
taj_col <- find_col(raw, c("TajimaD", "tajima_d"))
delta_col <- find_col(raw, c("mean_abs_delta_p", "abs_delta_p"))
f3_col <- find_col(raw, c("f3", "F3"))

detected_columns <- tibble(
  role = c("chromosome", "start", "end", "mid", "recombination", "avg_wc_fst", "avg_dxy", "pi_average", "TajimaD", "abs_delta_p", "f3"),
  detected_column = c(chr_col, start_col, end_col, mid_col, recomb_col, fst_col, dxy_col, pi_col, taj_col, delta_col, f3_col)
)

dat <- as.data.frame(raw) %>%
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
    f3 = if (!is.na(f3_col)) as.numeric(.data[[f3_col]]) else NA_real_
  )
assert_has_cols(dat, c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"), sanity, "killer_test_input_data")
assert_unique_keys(dat, c("Chromosome_Names", ".start_bp", ".end_bp"), sanity, "killer_test_input_data")

prior_local <- read.csv(file.path(cfg$prior_barrier_dir, "barrier_recalibrated_local_summary.csv"), stringsAsFactors = FALSE)
n_before_join <- nrow(dat)
dat <- dat %>%
  left_join(prior_local, by = c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"))
assert_join_rowcount(n_before_join, dat, sanity, "killer_test_prior_local_join")

if ("abs_delta_p.x" %in% names(dat) || "abs_delta_p.y" %in% names(dat)) dat$abs_delta_p <- dplyr::coalesce(dat$abs_delta_p.x, dat$abs_delta_p.y)
if ("f3.x" %in% names(dat) || "f3.y" %in% names(dat)) dat$f3 <- dplyr::coalesce(dat$f3.x, dat$f3.y)
if (!"sweep_group" %in% names(dat)) dat$sweep_group <- NA_character_
if ("sweep_group.x" %in% names(dat) || "sweep_group.y" %in% names(dat)) dat$sweep_group <- dplyr::coalesce(dat$sweep_group.y, dat$sweep_group.x)

dat$Chromosome_Names <- factor(dat$Chromosome_Names, levels = order_chromosomes(dat$Chromosome_Names))
dat <- dat %>% arrange(Chromosome_Names, .start_bp) %>% group_by(Chromosome_Names) %>% mutate(window_index = row_number()) %>% ungroup()

if (!"B_no_delta_p" %in% names(dat)) stop("B_no_delta_p missing after merge.")
if (!"K_local_mean" %in% names(dat)) stop("K_local_mean missing after merge.")
assert_missing_fraction(dat$B_no_delta_p, 0.1, sanity, "B_no_delta_p")
assert_missing_fraction(dat$K_local_mean, 0.1, sanity, "K_local_mean")

res_abs_base <- fit_model_summary(dat, "abs_delta_p", include_k = FALSE, include_chr = TRUE, model_name = "baseline")
res_abs_bar <- fit_model_summary(dat, "abs_delta_p", include_k = TRUE, include_chr = TRUE, model_name = "barrier")
model_comp <- bind_rows(res_abs_base$comparison, res_abs_bar$comparison)
coef_tbl <- bind_rows(res_abs_base$coefficients, res_abs_bar$coefficients)

if (any(is.finite(dat$f3))) {
  res_f3_base <- fit_model_summary(dat, "f3", include_k = FALSE, include_chr = TRUE, model_name = "baseline")
  res_f3_bar <- fit_model_summary(dat, "f3", include_k = TRUE, include_chr = TRUE, model_name = "barrier")
  model_comp <- bind_rows(model_comp, res_f3_base$comparison, res_f3_bar$comparison)
  coef_tbl <- bind_rows(coef_tbl, res_f3_base$coefficients, res_f3_bar$coefficients)
} else {
  res_f3_base <- list(model = NULL)
  res_f3_bar <- list(model = NULL)
}

loco_tbl <- bind_rows(
  leave_one_chr_out(dat, "abs_delta_p"),
  if (any(is.finite(dat$f3))) leave_one_chr_out(dat, "f3") else tibble()
)

within_bin_df <- dat %>%
  mutate(
    recomb_bin = qbin(recombination, cfg$n_bins),
    fst_bin = qbin(avg_wc_fst, cfg$n_bins),
    dxy_bin = qbin(avg_dxy, cfg$n_bins),
    sweep_group = ifelse(is.na(sweep_group), "Unknown", sweep_group),
    match_bin = interaction(recomb_bin, fst_bin, dxy_bin, sweep_group, drop = TRUE, lex.order = TRUE)
  )

within_bin_results <- bind_rows(lapply(split(within_bin_df, within_bin_df$match_bin), function(df_bin) {
  out <- list()
  if (nrow(df_bin) >= 10L && dplyr::n_distinct(df_bin$K_local_mean[is.finite(df_bin$K_local_mean)]) > 2L) {
    sp_abs <- safe_spearman(df_bin$K_local_mean, df_bin$abs_delta_p)
    out[[1]] <- tibble(
      outcome = "abs_delta_p", match_bin = as.character(df_bin$match_bin[1]), n = nrow(df_bin),
      spearman = sp_abs[["estimate"]], p_value = sp_abs[["p_value"]],
      mean_k = mean(df_bin$K_local_mean, na.rm = TRUE), dominant_sweep_group = dominant_label(df_bin$sweep_group)
    )
    if (any(is.finite(df_bin$f3))) {
      sp_f3 <- safe_spearman(df_bin$K_local_mean, df_bin$f3)
      out[[2]] <- tibble(
        outcome = "f3", match_bin = as.character(df_bin$match_bin[1]), n = nrow(df_bin),
        spearman = sp_f3[["estimate"]], p_value = sp_f3[["p_value"]],
        mean_k = mean(df_bin$K_local_mean, na.rm = TRUE), dominant_sweep_group = dominant_label(df_bin$sweep_group)
      )
    }
  }
  bind_rows(out)
}))

resid_plot_df <- tibble()
if (!is.null(res_abs_base$model)) {
  df0 <- res_abs_base$data
  df0$residual_outcome <- stats::residuals(res_abs_base$model)
  df0$outcome <- "abs_delta_p"
  resid_plot_df <- bind_rows(resid_plot_df, df0 %>% select(outcome, K_local_mean, residual_outcome))
}
if (!is.null(res_f3_base$model)) {
  df0 <- res_f3_base$data
  df0$residual_outcome <- stats::residuals(res_f3_base$model)
  df0$outcome <- "f3"
  resid_plot_df <- bind_rows(resid_plot_df, df0 %>% select(outcome, K_local_mean, residual_outcome))
}

q_high <- stats::quantile(dat$K_local_mean, probs = cfg$primary_high_quantile, na.rm = TRUE)
dat$highK_5 <- is.finite(dat$K_local_mean) & dat$K_local_mean >= q_high
dat$domain_context <- unlist(lapply(split(dat, dat$Chromosome_Names), function(df_chr) {
  classify_domain_context(df_chr, df_chr$highK_5, cfg$edge_window_span, cfg$nearby_outside_span)
}), use.names = FALSE)

edge_results <- bind_rows(
  safe_group_test(dat$abs_delta_p, dat$domain_context, "abs_delta_p by domain context"),
  if (any(is.finite(dat$f3))) safe_group_test(dat$f3, dat$domain_context, "f3 by domain context") else tibble()
) %>%
  bind_rows(
    dat %>%
      filter(!is.na(domain_context)) %>%
      group_by(domain_context) %>%
      summarise(
        analysis = "context_means",
        mean_abs_delta_p = mean(abs_delta_p, na.rm = TRUE),
        mean_f3 = mean(f3, na.rm = TRUE),
        mean_K_local = mean(K_local_mean, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
  )

obs_pred_abs <- loco_tbl %>%
  filter(outcome == "abs_delta_p") %>%
  left_join(dat %>% select(Chromosome_Names, abs_delta_p) %>% group_by(Chromosome_Names) %>% summarise(mean_abs_delta = mean(abs_delta_p, na.rm = TRUE), .groups = "drop"), by = c("held_out_chromosome" = "Chromosome_Names"))

loco_summary <- loco_tbl %>%
  group_by(outcome, model_name) %>%
  summarise(
    mean_test_r2 = mean(test_r_squared, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_spearman = mean(spearman, na.rm = TRUE),
    .groups = "drop"
  )

coef_plot_df <- coef_tbl %>%
  filter(term == "K_local_mean") %>%
  mutate(
    conf_low = estimate - 1.96 * std_error,
    conf_high = estimate + 1.96 * std_error,
    label = paste(outcome, model_name, sep = ": ")
  )

write.csv(model_comp, file.path(cfg$output_dir, "killer_test_model_comparisons.csv"), row.names = FALSE)
write.csv(loco_tbl, file.path(cfg$output_dir, "killer_test_leave_one_chromosome_out.csv"), row.names = FALSE)
write.csv(within_bin_results, file.path(cfg$output_dir, "killer_test_within_bin_results.csv"), row.names = FALSE)
write.csv(edge_results, file.path(cfg$output_dir, "killer_test_domain_edge_results.csv"), row.names = FALSE)
write.csv(coef_tbl, file.path(cfg$output_dir, "killer_test_coefficients_full.csv"), row.names = FALSE)
write.csv(detected_columns, file.path(cfg$output_dir, "detected_columns.csv"), row.names = FALSE)

if (nrow(loco_tbl)) {
  p_loco <- loco_tbl %>%
    ggplot(aes(x = held_out_chromosome, y = test_r_squared, fill = model_name)) +
    geom_col(position = "dodge") +
    facet_wrap(~ outcome, scales = "free_y") +
    labs(x = "Held-out chromosome", y = "Out-of-sample test R^2", fill = "Model") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(cfg$output_dir, "plot_heldout_chromosome_r2.png"), p_loco, width = 11, height = 6, dpi = 300)
}

if (nrow(coef_plot_df)) {
  p_coef <- ggplot(coef_plot_df, aes(x = label, y = estimate, ymin = conf_low, ymax = conf_high, color = outcome)) +
    geom_pointrange() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = NULL, y = "K_local_mean coefficient") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(cfg$output_dir, "plot_k_local_coefficient_residual_models.png"), p_coef, width = 8, height = 5, dpi = 300)
}

if (nrow(resid_plot_df)) {
  resid_plot_df2 <- bind_rows(lapply(split(resid_plot_df, resid_plot_df$outcome), function(df0) {
    if (nrow(df0) > 30000L) df0 <- df0[sample(seq_len(nrow(df0)), 30000L), , drop = FALSE]
    df0
  }))
  p_resid <- ggplot(resid_plot_df2, aes(x = K_local_mean, y = residual_outcome)) +
    geom_point(alpha = 0.18, size = 0.5, color = "grey35") +
    geom_smooth(method = "lm", se = FALSE, color = "red3") +
    facet_wrap(~ outcome, scales = "free_y") +
    labs(x = "K_local_mean", y = "Residual outcome after baseline controls") +
    theme_bw(base_size = 10)
  ggsave(file.path(cfg$output_dir, "plot_within_bin_residual_relationships.png"), p_resid, width = 8, height = 5, dpi = 300)
}

if (any(!is.na(dat$domain_context))) {
  edge_plot_df <- dat %>% filter(!is.na(domain_context))
  p_edge_abs <- ggplot(edge_plot_df, aes(x = domain_context, y = abs_delta_p, fill = domain_context)) +
    geom_boxplot(outlier.size = 0.25) +
    labs(x = NULL, y = "abs_delta_p") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
  ggsave(file.path(cfg$output_dir, "plot_domain_edge_abs_delta_p.png"), p_edge_abs, width = 7, height = 5, dpi = 300)
  if (any(is.finite(edge_plot_df$f3))) {
    p_edge_f3 <- ggplot(edge_plot_df, aes(x = domain_context, y = f3, fill = domain_context)) +
      geom_boxplot(outlier.size = 0.25) +
      labs(x = NULL, y = "f3") +
      theme_bw(base_size = 10) +
      theme(legend.position = "none")
    ggsave(file.path(cfg$output_dir, "plot_domain_edge_f3.png"), p_edge_f3, width = 7, height = 5, dpi = 300)
  }
}

abs_bar <- model_comp %>% filter(outcome == "abs_delta_p", model_name == "barrier") %>% slice_head(n = 1)
abs_base <- model_comp %>% filter(outcome == "abs_delta_p", model_name == "baseline") %>% slice_head(n = 1)
f3_bar <- model_comp %>% filter(outcome == "f3", model_name == "barrier") %>% slice_head(n = 1)
f3_base <- model_comp %>% filter(outcome == "f3", model_name == "baseline") %>% slice_head(n = 1)
abs_loco <- loco_summary %>% filter(outcome == "abs_delta_p")
f3_loco <- loco_summary %>% filter(outcome == "f3")
bin_abs <- within_bin_results %>% filter(outcome == "abs_delta_p")
bin_f3 <- within_bin_results %>% filter(outcome == "f3")
context_means <- edge_results %>% filter(analysis == "context_means")

cat("\nKiller test summary\n")
cat("Detected columns:\n")
apply(detected_columns, 1, function(x) cat(" -", x[["role"]], "=>", x[["detected_column"]], "\n"))
cat("\nDoes K_local_mean remain significant after controlling for focal predictors? ",
    ifelse(nrow(abs_bar) == 1L && is.finite(abs_bar$k_local_p_value) && abs_bar$k_local_p_value < 0.05, "Yes", "No"),
    "; abs_delta_p baseline R^2=", signif(abs_base$r_squared, 4),
    ", barrier R^2=", signif(abs_bar$r_squared, 4),
    ", delta R^2=", signif(abs_bar$delta_r2, 4),
    ", K_local p=", signif(abs_bar$k_local_p_value, 3), "\n", sep = "")
if (nrow(f3_bar) == 1L) {
  cat("Does K_local_mean remain significant for f3 after controls? ",
      ifelse(is.finite(f3_bar$k_local_p_value) && f3_bar$k_local_p_value < 0.05, "Yes", "No"),
      "; f3 baseline R^2=", signif(f3_base$r_squared, 4),
      ", barrier R^2=", signif(f3_bar$r_squared, 4),
      ", delta R^2=", signif(f3_bar$delta_r2, 4),
      ", K_local p=", signif(f3_bar$k_local_p_value, 3), "\n", sep = "")
}
if (nrow(abs_loco) == 2L) {
  cat("Does K_local_mean improve out-of-sample prediction on held-out chromosomes? ",
      ifelse(abs_loco$mean_test_r2[abs_loco$model_name == "barrier"] > abs_loco$mean_test_r2[abs_loco$model_name == "baseline"], "Yes", "No"),
      "; abs_delta_p mean held-out R^2 baseline=", signif(abs_loco$mean_test_r2[abs_loco$model_name == "baseline"], 4),
      ", barrier=", signif(abs_loco$mean_test_r2[abs_loco$model_name == "barrier"], 4), "\n", sep = "")
}
if (nrow(f3_loco) == 2L) {
  cat("Does K_local_mean improve held-out prediction for f3? ",
      ifelse(f3_loco$mean_test_r2[f3_loco$model_name == "barrier"] > f3_loco$mean_test_r2[f3_loco$model_name == "baseline"], "Yes", "No"),
      "; baseline=", signif(f3_loco$mean_test_r2[f3_loco$model_name == "baseline"], 4),
      ", barrier=", signif(f3_loco$mean_test_r2[f3_loco$model_name == "barrier"], 4), "\n", sep = "")
}
if (nrow(bin_abs)) {
  cat("Does K_local_mean still matter within matched recombination/divergence bins? ",
      ifelse(mean(bin_abs$spearman > 0, na.rm = TRUE) > 0.5, "Yes", "Mixed"),
      "; mean within-bin Spearman for abs_delta_p=", signif(mean(bin_abs$spearman, na.rm = TRUE), 4), "\n", sep = "")
}
if (nrow(bin_f3)) {
  cat("Within matched bins for f3, the mean Spearman with K_local_mean is ", signif(mean(bin_f3$spearman, na.rm = TRUE), 4), ".\n", sep = "")
}
if (nrow(context_means)) {
  inside_abs <- context_means$mean_abs_delta_p[context_means$domain_context == "Inside"]
  outside_abs <- context_means$mean_abs_delta_p[context_means$domain_context == "Nearby outside"]
  inside_f3 <- context_means$mean_f3[context_means$domain_context == "Inside"]
  outside_f3 <- context_means$mean_f3[context_means$domain_context == "Nearby outside"]
  cat("Do high-K_local domain interiors show stronger maintained differentiation and lower introgression than nearby outside regions? ",
      ifelse(length(inside_abs) && length(outside_abs) && inside_abs > outside_abs, "Yes", "Mixed"),
      "; inside abs_delta_p=", signif(inside_abs, 4),
      ", outside abs_delta_p=", signif(outside_abs, 4),
      ifelse(length(inside_f3) && length(outside_f3), paste0("; inside f3=", signif(inside_f3, 4), ", outside f3=", signif(outside_f3, 4)), ""),
      "\n", sep = "")
}
cat("Do the results support the claim that accumulated local barrier architecture helps maintain species differences? ",
    ifelse(nrow(abs_bar) == 1L && is.finite(abs_bar$k_local_p_value) && abs_bar$k_local_p_value < 0.05, "Yes", "Mixed"),
    ".\n", sep = "")
cat("Saved outputs to:", cfg$output_dir, "\n")
finalize_sanity(
  sanity,
  files = c(
    "killer_test_model_comparisons.csv",
    "killer_test_leave_one_chromosome_out.csv",
    "killer_test_within_bin_results.csv",
    "killer_test_domain_edge_results.csv",
    "killer_test_coefficients_full.csv",
    "detected_columns.csv",
    "plot_heldout_chromosome_r2.png",
    "plot_k_local_coefficient_residual_models.png",
    "plot_within_bin_residual_relationships.png",
    "plot_domain_edge_abs_delta_p.png"
  )
)
