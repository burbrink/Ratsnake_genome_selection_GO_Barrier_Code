#!/usr/bin/env Rscript

# Empirical, non-circular alpha calibration for the sd3_new_af barrier framework.
# -------------------------------------------------------------------------------
# This script calibrates alpha in:
#   P_hat = exp(-alpha * K_local_mean)
#
# using independent empirical proxies of compatibility / introgression rather than
# P values that were themselves generated from a chosen alpha.
#
# Primary logic:
# - lower abs_delta_p implies higher compatibility
# - higher f3 implies higher compatibility
# - a composite proxy can combine both
#
# Output files:
# - empirical_alpha_estimates.csv
# - empirical_alpha_loco_cv.csv
# - empirical_alpha_summary.txt
# - empirical_alpha_proxy_vs_K.pdf
# - empirical_alpha_cv.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})
script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

cfg <- list(
  input_csv = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_2of3_corrected", "barrier_recalibrated_local_summary.csv"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "alpha_empirical_calibration_outputs"),
  k_col = "K_local_mean",
  chr_col = "Chromosome_Names",
  abs_col = "abs_delta_p",
  f3_col = "f3",
  alpha_bounds = c(0.001, 2),
  eps = 1e-6
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_csv <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$output_dir <- args[[2]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

estimate_path <- file.path(cfg$output_dir, "empirical_alpha_estimates.csv")
cv_path <- file.path(cfg$output_dir, "empirical_alpha_loco_cv.csv")
summary_path <- file.path(cfg$output_dir, "empirical_alpha_summary.txt")
plot_path <- file.path(cfg$output_dir, "empirical_alpha_proxy_vs_K.pdf")
cv_plot_path <- file.path(cfg$output_dir, "empirical_alpha_cv.pdf")

if (!file.exists(cfg$input_csv)) stop("Input CSV not found: ", cfg$input_csv, call. = FALSE)
dat0 <- read.csv(cfg$input_csv, stringsAsFactors = FALSE, check.names = FALSE)

required <- c(cfg$k_col, cfg$chr_col)
missing_req <- setdiff(required, names(dat0))
if (length(missing_req)) stop("Missing required columns: ", paste(missing_req, collapse = ", "), call. = FALSE)

if (!(cfg$abs_col %in% names(dat0)) && !(cfg$f3_col %in% names(dat0))) {
  stop("Need at least one empirical proxy column: abs_delta_p and/or f3.", call. = FALSE)
}

safe_midrank01 <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n <- sum(ok)
  if (n == 0L) return(out)
  r <- rank(x[ok], ties.method = "average")
  out[ok] <- (r - 0.5) / n
  out
}

clip01 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)

fit_alpha_direct <- function(df, proxy_col, k_col, lower, upper) {
  dd <- df[, c(proxy_col, k_col), drop = FALSE]
  names(dd) <- c("proxy", "K")
  dd <- dd[is.finite(dd$proxy) & is.finite(dd$K), , drop = FALSE]
  if (nrow(dd) < 10L) return(NULL)

  obj <- function(alpha) {
    mean((dd$proxy - exp(-alpha * dd$K))^2)
  }
  opt <- optimize(obj, interval = c(lower, upper))
  alpha_hat <- opt$minimum
  pred <- exp(-alpha_hat * dd$K)
  rmse <- sqrt(mean((dd$proxy - pred)^2))
  spearman <- suppressWarnings(cor(dd$proxy, pred, method = "spearman"))
  pearson <- suppressWarnings(cor(dd$proxy, pred, method = "pearson"))
  tibble(
    method = "direct_mse",
    proxy = proxy_col,
    n = nrow(dd),
    alpha_estimate = alpha_hat,
    rmse = rmse,
    spearman_proxy_pred = spearman,
    pearson_proxy_pred = pearson
  )
}

fit_alpha_loglinear <- function(df, proxy_col, k_col, with_intercept = FALSE) {
  dd <- df[, c(proxy_col, k_col), drop = FALSE]
  names(dd) <- c("proxy", "K")
  dd <- dd[is.finite(dd$proxy) & is.finite(dd$K) & dd$proxy > 0, , drop = FALSE]
  if (nrow(dd) < 10L) return(NULL)
  dd$log_proxy <- log(dd$proxy)
  fit <- if (with_intercept) lm(log_proxy ~ K, data = dd) else lm(log_proxy ~ 0 + K, data = dd)
  sm <- summary(fit)
  cf <- sm$coefficients
  slope_name <- if ("K" %in% rownames(cf)) "K" else tail(rownames(cf), 1)
  alpha_hat <- -unname(cf[slope_name, "Estimate"])
  se <- unname(cf[slope_name, "Std. Error"])
  ci <- tryCatch(confint(fit, parm = slope_name), error = function(e) c(NA_real_, NA_real_))
  ci_low <- -ci[2]
  ci_high <- -ci[1]
  if (ci_low > ci_high) {
    tmp <- ci_low
    ci_low <- ci_high
    ci_high <- tmp
  }
  pred <- exp(predict(fit, newdata = dd))
  tibble(
    method = if (with_intercept) "loglinear_intercept" else "loglinear_no_intercept",
    proxy = proxy_col,
    n = nrow(dd),
    alpha_estimate = alpha_hat,
    std_error = se,
    conf_low = ci_low,
    conf_high = ci_high,
    r_squared = unname(sm$r.squared),
    adj_r_squared = unname(sm$adj.r.squared),
    aic = AIC(fit),
    intercept = if ("(Intercept)" %in% rownames(cf)) unname(cf["(Intercept)", "Estimate"]) else 0,
    rmse = sqrt(mean((dd$proxy - pred)^2)),
    spearman_proxy_pred = suppressWarnings(cor(dd$proxy, pred, method = "spearman")),
    pearson_proxy_pred = suppressWarnings(cor(dd$proxy, pred, method = "pearson"))
  )
}

loco_alpha_cv <- function(df, proxy_col, k_col, chr_col, lower, upper) {
  dd <- df[, c(proxy_col, k_col, chr_col), drop = FALSE]
  names(dd) <- c("proxy", "K", "chr")
  dd <- dd[is.finite(dd$proxy) & is.finite(dd$K) & !is.na(dd$chr), , drop = FALSE]
  chrs <- unique(dd$chr)
  out <- lapply(chrs, function(chr) {
    train <- dd[dd$chr != chr, , drop = FALSE]
    test <- dd[dd$chr == chr, , drop = FALSE]
    if (nrow(train) < 10L || nrow(test) < 5L) return(NULL)
    opt <- optimize(function(alpha) mean((train$proxy - exp(-alpha * train$K))^2), interval = c(lower, upper))
    alpha_hat <- opt$minimum
    pred <- exp(-alpha_hat * test$K)
    tibble(
      proxy = proxy_col,
      heldout_chr = chr,
      n_train = nrow(train),
      n_test = nrow(test),
      alpha_estimate = alpha_hat,
      test_rmse = sqrt(mean((test$proxy - pred)^2)),
      test_spearman = suppressWarnings(cor(test$proxy, pred, method = "spearman")),
      test_pearson = suppressWarnings(cor(test$proxy, pred, method = "pearson"))
    )
  })
  bind_rows(out)
}

dat <- dat0 %>%
  transmute(
    Chromosome_Names = .data[[cfg$chr_col]],
    K_local_mean = as.numeric(.data[[cfg$k_col]]),
    abs_delta_p = if (cfg$abs_col %in% names(dat0)) as.numeric(.data[[cfg$abs_col]]) else NA_real_,
    f3 = if (cfg$f3_col %in% names(dat0)) as.numeric(.data[[cfg$f3_col]]) else NA_real_
  ) %>%
  mutate(
    compat_from_abs_delta = clip01(1 - safe_midrank01(abs_delta_p), cfg$eps),
    compat_from_f3 = clip01(safe_midrank01(f3), cfg$eps)
  )

if (all(!is.finite(dat$compat_from_abs_delta)) && all(!is.finite(dat$compat_from_f3))) {
  stop("No usable empirical compatibility proxies could be constructed.", call. = FALSE)
}

if (any(is.finite(dat$compat_from_abs_delta)) && any(is.finite(dat$compat_from_f3))) {
  dat$compat_composite <- rowMeans(cbind(dat$compat_from_abs_delta, dat$compat_from_f3), na.rm = TRUE)
} else if (any(is.finite(dat$compat_from_f3))) {
  dat$compat_composite <- dat$compat_from_f3
} else {
  dat$compat_composite <- dat$compat_from_abs_delta
}
dat$compat_composite <- clip01(dat$compat_composite, cfg$eps)

proxy_cols <- c("compat_from_abs_delta", "compat_from_f3", "compat_composite")
proxy_cols <- proxy_cols[vapply(proxy_cols, function(x) any(is.finite(dat[[x]])), logical(1))]

estimate_rows <- bind_rows(lapply(proxy_cols, function(proxy_col) {
  bind_rows(
    fit_alpha_direct(dat, proxy_col, "K_local_mean", cfg$alpha_bounds[1], cfg$alpha_bounds[2]),
    fit_alpha_loglinear(dat, proxy_col, "K_local_mean", with_intercept = FALSE),
    fit_alpha_loglinear(dat, proxy_col, "K_local_mean", with_intercept = TRUE)
  )
}))

cv_rows <- bind_rows(lapply(proxy_cols, function(proxy_col) {
  loco_alpha_cv(dat, proxy_col, "K_local_mean", "Chromosome_Names", cfg$alpha_bounds[1], cfg$alpha_bounds[2])
}))

write.csv(estimate_rows, estimate_path, row.names = FALSE)
write.csv(cv_rows, cv_path, row.names = FALSE)

best_proxy <- "compat_composite"
if (!(best_proxy %in% proxy_cols)) best_proxy <- proxy_cols[[1]]
best_row <- estimate_rows %>%
  filter(proxy == best_proxy, method == "direct_mse") %>%
  slice_head(n = 1)

plot_df <- dat %>%
  select(Chromosome_Names, K_local_mean, all_of(proxy_cols)) %>%
  tidyr::pivot_longer(cols = all_of(proxy_cols), names_to = "proxy", values_to = "compatibility") %>%
  filter(is.finite(K_local_mean), is.finite(compatibility))

curve_df <- estimate_rows %>%
  filter(method == "direct_mse") %>%
  select(proxy, alpha_estimate) %>%
  distinct() %>%
  tidyr::expand_grid(K_local_mean = seq(min(plot_df$K_local_mean, na.rm = TRUE), max(plot_df$K_local_mean, na.rm = TRUE), length.out = 300)) %>%
  mutate(predicted_compatibility = exp(-alpha_estimate * K_local_mean))

p1 <- ggplot(plot_df, aes(x = K_local_mean, y = compatibility)) +
  geom_point(alpha = 0.08, size = 0.45, color = "#2b6cb0") +
  geom_line(data = curve_df, aes(y = predicted_compatibility), color = "#c53030", linewidth = 1) +
  facet_wrap(~ proxy, scales = "free_y") +
  labs(
    title = "Empirical Compatibility Proxies vs K_local_mean",
    subtitle = "Red curve: P_hat = exp(-alpha * K_local_mean), alpha fitted by direct MSE",
    x = "K_local_mean",
    y = "Empirical compatibility proxy"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(plot_path, p1, width = 10, height = 7, units = "in")

cv_summary <- cv_rows %>%
  group_by(proxy) %>%
  summarise(
    mean_alpha = mean(alpha_estimate, na.rm = TRUE),
    sd_alpha = sd(alpha_estimate, na.rm = TRUE),
    mean_test_rmse = mean(test_rmse, na.rm = TRUE),
    mean_test_spearman = mean(test_spearman, na.rm = TRUE),
    mean_test_pearson = mean(test_pearson, na.rm = TRUE),
    .groups = "drop"
  )

p2 <- ggplot(cv_rows, aes(x = proxy, y = alpha_estimate, fill = proxy)) +
  geom_violin(alpha = 0.5, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.1) +
  labs(
    title = "Leave-One-Chromosome-Out Alpha Estimates",
    subtitle = "Alpha fitted by direct MSE on training chromosomes",
    x = NULL,
    y = "Estimated alpha"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none", panel.grid.minor = element_blank())
ggsave(cv_plot_path, p2, width = 8, height = 5.5, units = "in")

summary_lines <- c(
  "Empirical alpha calibration summary",
  "===================================",
  "",
  paste0("Input file: ", cfg$input_csv),
  paste0("K column: ", cfg$k_col),
  paste0("Chromosome column: ", cfg$chr_col),
  paste0("Compatibility proxies fitted: ", paste(proxy_cols, collapse = ", ")),
  ""
)

for (proxy_col in proxy_cols) {
  est_sub <- estimate_rows %>% filter(proxy == proxy_col)
  cv_sub <- cv_summary %>% filter(proxy == proxy_col)
  summary_lines <- c(summary_lines, paste0("Proxy: ", proxy_col), strrep("-", 7 + nchar(proxy_col)))
  for (i in seq_len(nrow(est_sub))) {
    row <- est_sub[i, ]
    summary_lines <- c(
      summary_lines,
      paste0("Method: ", row$method),
      paste0("  alpha: ", formatC(row$alpha_estimate, digits = 6, format = "fg", flag = "#")),
      if ("std_error" %in% names(row) && is.finite(row$std_error)) paste0("  SE: ", formatC(row$std_error, digits = 6, format = "fg", flag = "#")) else NULL,
      if ("conf_low" %in% names(row) && is.finite(row$conf_low)) paste0("  95% CI: [", formatC(row$conf_low, digits = 6, format = "fg", flag = "#"), ", ", formatC(row$conf_high, digits = 6, format = "fg", flag = "#"), "]") else NULL,
      if ("r_squared" %in% names(row) && is.finite(row$r_squared)) paste0("  R^2: ", formatC(row$r_squared, digits = 6, format = "fg", flag = "#")) else NULL,
      paste0("  RMSE(proxy vs fitted): ", formatC(row$rmse, digits = 6, format = "fg", flag = "#")),
      paste0("  Spearman(proxy, fitted): ", formatC(row$spearman_proxy_pred, digits = 6, format = "fg", flag = "#")),
      paste0("  Pearson(proxy, fitted): ", formatC(row$pearson_proxy_pred, digits = 6, format = "fg", flag = "#")),
      ""
    )
  }
  if (nrow(cv_sub)) {
    summary_lines <- c(
      summary_lines,
      "LOCO CV (direct MSE alpha):",
      paste0("  mean alpha: ", formatC(cv_sub$mean_alpha, digits = 6, format = "fg", flag = "#")),
      paste0("  sd alpha: ", formatC(cv_sub$sd_alpha, digits = 6, format = "fg", flag = "#")),
      paste0("  mean test RMSE: ", formatC(cv_sub$mean_test_rmse, digits = 6, format = "fg", flag = "#")),
      paste0("  mean test Spearman: ", formatC(cv_sub$mean_test_spearman, digits = 6, format = "fg", flag = "#")),
      paste0("  mean test Pearson: ", formatC(cv_sub$mean_test_pearson, digits = 6, format = "fg", flag = "#")),
      ""
    )
  }
}

if (nrow(best_row) == 1) {
  summary_lines <- c(
    summary_lines,
    "Recommended working alpha",
    "-------------------------",
    paste0("Proxy: ", best_proxy),
    paste0("Method: direct_mse"),
    paste0("Alpha: ", formatC(best_row$alpha_estimate, digits = 6, format = "fg", flag = "#")),
    "This is an empirical compatibility calibration using independent introgression/differentiation proxies, not a back-calculation from previously derived P values.",
    ""
  )
}

writeLines(summary_lines, summary_path)

cat("\n================ Empirical Alpha Calibration ================\n")
cat("Input file: ", cfg$input_csv, "\n", sep = "")
for (proxy_col in proxy_cols) {
  row <- estimate_rows %>% filter(proxy == proxy_col, method == "direct_mse") %>% slice_head(n = 1)
  cv_row <- cv_summary %>% filter(proxy == proxy_col)
  if (nrow(row) == 1) {
    cat("Proxy ", proxy_col, ": alpha=", signif(row$alpha_estimate, 6),
        ", RMSE=", signif(row$rmse, 4),
        ", Spearman(proxy,fitted)=", signif(row$spearman_proxy_pred, 4),
        if (nrow(cv_row) == 1) paste0(", LOCO mean alpha=", signif(cv_row$mean_alpha, 6),
                                       ", LOCO mean test Spearman=", signif(cv_row$mean_test_spearman, 4)) else "",
        "\n", sep = "")
  }
}
if (nrow(best_row) == 1) {
  cat("Recommended working alpha: ", signif(best_row$alpha_estimate, 6), " (proxy=", best_proxy, ", method=direct_mse)\n", sep = "")
}
cat("Outputs written to: ", cfg$output_dir, "\n", sep = "")
