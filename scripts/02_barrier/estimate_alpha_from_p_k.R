#!/usr/bin/env Rscript

# Estimate the scaling parameter alpha from:
#   P = exp(-alpha * K)
#
# This script supports two use cases:
# 1. Single supplied values of P_overall_mean and K_genome_mean
# 2. A data frame / CSV with multiple observations of P and K
#
# Outputs:
# - alpha_estimates.csv
# - alpha_model_summary.txt
# - alpha_fit_plot.pdf
# - alpha_log_fit_plot.pdf
#
# Notes:
# - alpha is estimated from log(P) = -alpha * K
# - the no-intercept fit enforces the model above directly
# - the intercept fit allows a baseline offset and can reveal misspecification

suppressPackageStartupMessages({
  library(ggplot2)
})

cfg <- list(
  # Optional single-value inputs.
  single_P = NA_real_,
  single_K = NA_real_,

  # Optional table input.
  input_csv = NA_character_,
  p_col = "P_overall_mean",
  k_col = "K_genome_mean",

  # Output directory.
  output_dir = getwd()
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && identical(args[[1]], "--single")) {
  if (length(args) < 3) stop("For --single mode, provide: --single P K [output_dir]", call. = FALSE)
  cfg$single_P <- as.numeric(args[[2]])
  cfg$single_K <- as.numeric(args[[3]])
  if (length(args) >= 4 && nzchar(args[[4]])) cfg$output_dir <- args[[4]]
} else {
  if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_csv <- args[[1]]
  if (length(args) >= 2 && nzchar(args[[2]])) cfg$p_col <- args[[2]]
  if (length(args) >= 3 && nzchar(args[[3]])) cfg$k_col <- args[[3]]
  if (length(args) >= 4 && nzchar(args[[4]])) cfg$output_dir <- args[[4]]
}

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

csv_path <- file.path(cfg$output_dir, "alpha_estimates.csv")
txt_path <- file.path(cfg$output_dir, "alpha_model_summary.txt")
fit_plot_path <- file.path(cfg$output_dir, "alpha_fit_plot.pdf")
log_plot_path <- file.path(cfg$output_dir, "alpha_log_fit_plot.pdf")

warn_msgs <- character()

safe_fmt <- function(x, digits = 6) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  formatC(x, digits = digits, format = "fg", flag = "#")
}

extract_model_row <- function(fit, model_name) {
  sm <- summary(fit)
  cf <- sm$coefficients
  slope_name <- if ("K" %in% rownames(cf)) "K" else tail(rownames(cf), 1)
  slope <- unname(cf[slope_name, "Estimate"])
  se <- unname(cf[slope_name, "Std. Error"])
  pval <- unname(cf[slope_name, grep("^Pr", colnames(cf))])
  alpha_hat <- -slope
  ci <- tryCatch(confint(fit, parm = slope_name), error = function(e) c(NA_real_, NA_real_))
  alpha_ci_low <- -ci[2]
  alpha_ci_high <- -ci[1]
  if (alpha_ci_low > alpha_ci_high) {
    tmp <- alpha_ci_low
    alpha_ci_low <- alpha_ci_high
    alpha_ci_high <- tmp
  }
  data.frame(
    model = model_name,
    n = nobs(fit),
    intercept = if ("(Intercept)" %in% rownames(cf)) unname(cf["(Intercept)", "Estimate"]) else 0,
    slope_K = slope,
    alpha_estimate = alpha_hat,
    std_error = se,
    conf_low = alpha_ci_low,
    conf_high = alpha_ci_high,
    r_squared = unname(sm$r.squared),
    adj_r_squared = unname(sm$adj.r.squared),
    aic = AIC(fit),
    p_value = pval,
    stringsAsFactors = FALSE
  )
}

clean_pk_data <- function(df, p_col, k_col) {
  keep_cols <- unique(c(p_col, k_col))
  dat <- df[, keep_cols, drop = FALSE]
  names(dat) <- c("P", "K")

  original_n <- nrow(dat)
  dat <- dat[is.finite(dat$P) & is.finite(dat$K), , drop = FALSE]
  removed_nonfinite <- original_n - nrow(dat)

  if (any(dat$P > 1, na.rm = TRUE)) {
    warn_msgs <<- c(warn_msgs, "Some P values were > 1. These were retained, but interpreted with caution.")
    warning("Some P values were > 1. These were retained, but interpreted with caution.", call. = FALSE)
  }

  removed_nonpositive <- sum(dat$P <= 0, na.rm = TRUE)
  dat <- dat[dat$P > 0, , drop = FALSE]
  dat$logP <- log(dat$P)

  list(
    data = dat,
    removed_nonfinite = removed_nonfinite,
    removed_nonpositive = removed_nonpositive
  )
}

results_list <- list()
summary_lines <- c(
  "Alpha estimation summary",
  "========================",
  ""
)

single_mode <- is.finite(cfg$single_P) && is.finite(cfg$single_K)
table_mode <- is.character(cfg$input_csv) && length(cfg$input_csv) == 1 && !is.na(cfg$input_csv) && nzchar(cfg$input_csv)

if (!single_mode && !table_mode) {
  summary_lines <- c(
    summary_lines,
    "No valid single-value inputs or input CSV were supplied.",
    "Edit cfg$single_P and cfg$single_K, or run:",
    "Rscript estimate_alpha_from_p_k.R --single P K output_dir",
    "or",
    "Rscript estimate_alpha_from_p_k.R data.csv P_col K_col output_dir",
    ""
  )
}

if (single_mode) {
  if (cfg$single_P <= 0) stop("single_P must be > 0.", call. = FALSE)
  if (cfg$single_K == 0) stop("single_K must be non-zero.", call. = FALSE)
  if (cfg$single_P > 1) {
    warn_msgs <- c(warn_msgs, "single_P > 1. This is unusual for a probability-like quantity.")
    warning("single_P > 1. This is unusual for a probability-like quantity.", call. = FALSE)
  }

  alpha_single <- -log(cfg$single_P) / cfg$single_K
  results_list[["single_value"]] <- data.frame(
    model = "single_value",
    n = 1,
    intercept = NA_real_,
    slope_K = -alpha_single,
    alpha_estimate = alpha_single,
    std_error = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    r_squared = NA_real_,
    adj_r_squared = NA_real_,
    aic = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  summary_lines <- c(
    summary_lines,
    "Single-value estimate",
    "---------------------",
    paste0("P = ", safe_fmt(cfg$single_P)),
    paste0("K = ", safe_fmt(cfg$single_K)),
    paste0("alpha = -log(P) / K = ", safe_fmt(alpha_single)),
    ""
  )
}

if (table_mode) {
  if (!file.exists(cfg$input_csv)) stop("Input CSV not found: ", cfg$input_csv, call. = FALSE)
  raw_df <- read.csv(cfg$input_csv, stringsAsFactors = FALSE, check.names = FALSE)

  if (!(cfg$p_col %in% names(raw_df))) stop("P column not found: ", cfg$p_col, call. = FALSE)
  if (!(cfg$k_col %in% names(raw_df))) stop("K column not found: ", cfg$k_col, call. = FALSE)

  cleaned <- clean_pk_data(raw_df, cfg$p_col, cfg$k_col)
  dat <- cleaned$data

  if (nrow(dat) < 3) stop("Need at least 3 finite observations with P > 0 to fit models.", call. = FALSE)

  fit_no_intercept <- lm(logP ~ 0 + K, data = dat)
  fit_intercept <- lm(logP ~ K, data = dat)

  results_list[["no_intercept"]] <- extract_model_row(fit_no_intercept, "log(P) ~ 0 + K")
  results_list[["with_intercept"]] <- extract_model_row(fit_intercept, "log(P) ~ K")

  preferred_model <- fit_no_intercept
  preferred_name <- "log(P) ~ 0 + K"
  preferred_row <- results_list[["no_intercept"]]

  summary_lines <- c(
    summary_lines,
    "Multi-observation regression estimates",
    "-------------------------------------",
    paste0("Input file: ", normalizePath(cfg$input_csv, winslash = "/", mustWork = FALSE)),
    paste0("P column: ", cfg$p_col),
    paste0("K column: ", cfg$k_col),
    paste0("Rows retained: ", nrow(dat)),
    paste0("Rows removed for non-finite values: ", cleaned$removed_nonfinite),
    paste0("Rows removed for P <= 0: ", cleaned$removed_nonpositive),
    ""
  )

  for (nm in c("no_intercept", "with_intercept")) {
    row <- results_list[[nm]]
    summary_lines <- c(
      summary_lines,
      paste0("Model: ", row$model),
      paste0("  alpha estimate: ", safe_fmt(row$alpha_estimate)),
      paste0("  standard error: ", safe_fmt(row$std_error)),
      paste0("  95% CI: [", safe_fmt(row$conf_low), ", ", safe_fmt(row$conf_high), "]"),
      paste0("  R^2: ", safe_fmt(row$r_squared)),
      paste0("  adjusted R^2: ", safe_fmt(row$adj_r_squared)),
      paste0("  AIC: ", safe_fmt(row$aic)),
      paste0("  p-value for slope: ", safe_fmt(row$p_value)),
      ""
    )
  }

  # Publication-quality plot on original scale.
  grid_K <- seq(min(dat$K, na.rm = TRUE), max(dat$K, na.rm = TRUE), length.out = 300)
  pred_df <- data.frame(
    K = grid_K,
    P_no_intercept = exp(predict(fit_no_intercept, newdata = data.frame(K = grid_K))),
    P_with_intercept = exp(predict(fit_intercept, newdata = data.frame(K = grid_K)))
  )

  p_fit <- ggplot(dat, aes(x = K, y = P)) +
    geom_point(color = "#2b6cb0", alpha = 0.65, size = 2) +
    geom_line(data = pred_df, aes(x = K, y = P_no_intercept), color = "#c53030", linewidth = 1.1) +
    geom_line(data = pred_df, aes(x = K, y = P_with_intercept), color = "#2f855a", linewidth = 1.0, linetype = "22") +
    labs(
      title = "Exponential Fit of P vs K",
      subtitle = "Red: log(P) ~ 0 + K; Green dashed: log(P) ~ K",
      x = "K",
      y = "P"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  ggsave(fit_plot_path, p_fit, width = 7.2, height = 5.4, units = "in")

  # Log-scale fit plus diagnostics.
  p_log <- ggplot(dat, aes(x = K, y = logP)) +
    geom_point(color = "#2b6cb0", alpha = 0.65, size = 2) +
    geom_abline(
      intercept = 0,
      slope = coef(fit_no_intercept)[["K"]],
      color = "#c53030",
      linewidth = 1.1
    ) +
    geom_abline(
      intercept = coef(fit_intercept)[["(Intercept)"]],
      slope = coef(fit_intercept)[["K"]],
      color = "#2f855a",
      linewidth = 1.0,
      linetype = "22"
    ) +
    labs(
      title = "Linearized Fit of log(P) vs K",
      subtitle = "Red: no intercept; Green dashed: with intercept",
      x = "K",
      y = "log(P)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  pdf(log_plot_path, width = 7.2, height = 5.4)
  print(p_log)
  par(mfrow = c(2, 2))
  plot(preferred_model, which = 1:4)
  dev.off()

  summary_lines <- c(
    summary_lines,
    paste0("Diagnostic plots in ", basename(log_plot_path), " use the preferred model: ", preferred_name),
    ""
  )
} else {
  # Empty placeholder PDFs if only single-value mode was used.
  if (single_mode) {
    pdf(fit_plot_path, width = 7.2, height = 5.4)
    plot.new()
    title("No multi-observation data supplied")
    dev.off()

    pdf(log_plot_path, width = 7.2, height = 5.4)
    plot.new()
    title("No multi-observation data supplied")
    dev.off()
  }
}

if (length(results_list)) {
  out_df <- do.call(rbind, results_list)
  write.csv(out_df, csv_path, row.names = FALSE)
} else {
  write.csv(data.frame(), csv_path, row.names = FALSE)
}

if (length(warn_msgs)) {
  summary_lines <- c(summary_lines, "Warnings", "--------", unique(warn_msgs), "")
}

writeLines(summary_lines, con = txt_path)

message("Wrote: ", csv_path)
message("Wrote: ", txt_path)
message("Wrote: ", fit_plot_path)
message("Wrote: ", log_plot_path)
