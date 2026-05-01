estimate_alpha_from_pk <- function(df = NULL,
                                   p_col = "P_overall_mean",
                                   k_col = "K_genome_mean",
                                   single_P = NA_real_,
                                   single_K = NA_real_,
                                   prefer_no_intercept = TRUE) {
  warn_msgs <- character()

  clean_pk_data <- function(df, p_col, k_col) {
    dat <- df[, unique(c(p_col, k_col)), drop = FALSE]
    names(dat) <- c("P", "K")
    original_n <- nrow(dat)
    dat <- dat[is.finite(dat$P) & is.finite(dat$K), , drop = FALSE]
    removed_nonfinite <- original_n - nrow(dat)
    if (any(dat$P > 1, na.rm = TRUE)) {
      warn_msgs <<- c(warn_msgs, "Some P values were > 1. These were retained, but interpreted with caution.")
    }
    removed_nonpositive <- sum(dat$P <= 0, na.rm = TRUE)
    dat <- dat[dat$P > 0, , drop = FALSE]
    dat$logP <- log(dat$P)
    list(data = dat, removed_nonfinite = removed_nonfinite, removed_nonpositive = removed_nonpositive)
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

  if (is.finite(single_P) && is.finite(single_K)) {
    if (single_P <= 0) stop("single_P must be > 0.", call. = FALSE)
    if (single_K == 0) stop("single_K must be non-zero.", call. = FALSE)
    if (single_P > 1) warn_msgs <- c(warn_msgs, "single_P > 1. This is unusual for a probability-like quantity.")
    alpha_single <- -log(single_P) / single_K
    out <- data.frame(
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
    return(list(
      alpha = alpha_single,
      preferred_model = "single_value",
      estimates = out,
      cleaned = NULL,
      warnings = unique(warn_msgs)
    ))
  }

  if (is.null(df)) stop("Provide either df or single_P/single_K.", call. = FALSE)
  if (!(p_col %in% names(df))) stop("P column not found: ", p_col, call. = FALSE)
  if (!(k_col %in% names(df))) stop("K column not found: ", k_col, call. = FALSE)

  cleaned <- clean_pk_data(df, p_col, k_col)
  dat <- cleaned$data
  if (nrow(dat) == 0) stop("No finite observations with P > 0 were available for alpha calibration.", call. = FALSE)
  if (nrow(dat) == 1) {
    return(estimate_alpha_from_pk(
      single_P = dat$P[[1]],
      single_K = dat$K[[1]]
    ))
  }

  fit_no_intercept <- lm(logP ~ 0 + K, data = dat)
  fit_intercept <- lm(logP ~ K, data = dat)

  est_no_intercept <- extract_model_row(fit_no_intercept, "log(P) ~ 0 + K")
  est_intercept <- extract_model_row(fit_intercept, "log(P) ~ K")
  estimates <- rbind(est_no_intercept, est_intercept)
  preferred_idx <- if (prefer_no_intercept) 1 else 2

  list(
    alpha = estimates$alpha_estimate[preferred_idx],
    preferred_model = estimates$model[preferred_idx],
    estimates = estimates,
    cleaned = cleaned,
    warnings = unique(warn_msgs)
  )
}

write_alpha_calibration_outputs <- function(alpha_result, output_dir, prefix = "alpha_calibration") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  csv_path <- file.path(output_dir, paste0(prefix, "_estimates.csv"))
  txt_path <- file.path(output_dir, paste0(prefix, "_summary.txt"))

  write.csv(alpha_result$estimates, csv_path, row.names = FALSE)

  lines <- c(
    "Alpha calibration summary",
    "=========================",
    paste0("Preferred model: ", alpha_result$preferred_model),
    paste0("Preferred alpha: ", formatC(alpha_result$alpha, digits = 6, format = "fg", flag = "#")),
    ""
  )
  if (!is.null(alpha_result$cleaned)) {
    lines <- c(
      lines,
      paste0("Rows removed for non-finite values: ", alpha_result$cleaned$removed_nonfinite),
      paste0("Rows removed for P <= 0: ", alpha_result$cleaned$removed_nonpositive),
      ""
    )
  }
  if (length(alpha_result$warnings)) {
    lines <- c(lines, "Warnings", "--------", alpha_result$warnings, "")
  }
  writeLines(lines, txt_path)
  invisible(list(csv = csv_path, txt = txt_path))
}
