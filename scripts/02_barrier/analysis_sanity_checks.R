init_sanity <- function(script_name, output_dir) {
  tracker <- new.env(parent = emptyenv())
  tracker$script_name <- script_name
  tracker$output_dir <- output_dir
  tracker$records <- data.frame(
    script = character(),
    check_name = character(),
    status = character(),
    details = character(),
    stringsAsFactors = FALSE
  )
  tracker
}

append_check <- function(tracker, check_name, status, details = "") {
  if (is.null(tracker) || !is.environment(tracker)) return(invisible(NULL))
  tracker$records <- rbind(
    tracker$records,
    data.frame(
      script = tracker$script_name,
      check_name = check_name,
      status = status,
      details = as.character(details),
      stringsAsFactors = FALSE
    )
  )
  invisible(NULL)
}

write_sanity_report <- function(tracker, filename = "sanity_report.csv") {
  if (is.null(tracker) || !is.environment(tracker)) return(invisible(NULL))
  if (!dir.exists(tracker$output_dir)) dir.create(tracker$output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tracker$records, file.path(tracker$output_dir, filename), row.names = FALSE)
  invisible(NULL)
}

assert_true <- function(cond, msg, tracker = NULL, check_name = msg) {
  if (!isTRUE(cond)) {
    append_check(tracker, check_name, "fail", msg)
    write_sanity_report(tracker)
    stop(msg, call. = FALSE)
  }
  append_check(tracker, check_name, "pass", msg)
  invisible(TRUE)
}

warn_if_false <- function(cond, msg, tracker = NULL, check_name = msg) {
  if (!isTRUE(cond)) {
    append_check(tracker, check_name, "warn", msg)
    warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }
  append_check(tracker, check_name, "pass", msg)
  invisible(TRUE)
}

assert_has_cols <- function(df, cols, tracker = NULL, df_name = "data") {
  miss <- setdiff(cols, names(df))
  assert_true(length(miss) == 0L,
              paste0(df_name, " is missing columns: ", paste(miss, collapse = ", ")),
              tracker,
              paste0("required_columns_", df_name))
}

assert_unique_keys <- function(df, cols, tracker = NULL, df_name = "data") {
  dup_n <- sum(duplicated(df[, cols, drop = FALSE]))
  assert_true(dup_n == 0L,
              paste0(df_name, " has duplicated keys for: ", paste(cols, collapse = ", "), " (n_dup=", dup_n, ")"),
              tracker,
              paste0("unique_keys_", df_name))
}

assert_join_rowcount <- function(n_before, df_after, tracker = NULL, label = "join") {
  assert_true(nrow(df_after) == n_before,
              paste0(label, " changed row count from ", n_before, " to ", nrow(df_after)),
              tracker,
              paste0("join_rowcount_", label))
}

assert_missing_fraction <- function(x, max_frac, tracker = NULL, label = deparse(substitute(x))) {
  frac <- mean(is.na(x))
  assert_true(is.finite(frac) && frac <= max_frac,
              paste0(label, " missing fraction ", signif(frac, 4), " exceeds ", max_frac),
              tracker,
              paste0("missing_fraction_", label))
}

warn_missing_fraction <- function(x, max_frac, tracker = NULL, label = deparse(substitute(x))) {
  frac <- mean(is.na(x))
  warn_if_false(is.finite(frac) && frac <= max_frac,
                paste0(label, " missing fraction ", signif(frac, 4), " exceeds ", max_frac),
                tracker,
                paste0("missing_fraction_warn_", label))
}

assert_quantile_fraction <- function(flag, target, tol = 0.01, tracker = NULL, label = deparse(substitute(flag))) {
  frac <- mean(flag, na.rm = TRUE)
  assert_true(is.finite(frac) && abs(frac - target) <= tol,
              paste0(label, " fraction ", signif(frac, 4), " is not within ", tol, " of target ", target),
              tracker,
              paste0("quantile_fraction_", label))
}

assert_probability_bounds <- function(x, lower = 0, upper = 1, tracker = NULL, label = deparse(substitute(x))) {
  keep <- is.finite(x)
  ok <- all(x[keep] > lower & x[keep] <= upper)
  assert_true(ok,
              paste0(label, " contains values outside (", lower, ", ", upper, "]"),
              tracker,
              paste0("probability_bounds_", label))
}

check_complete_cases <- function(df, cols, min_n = 50L, tracker = NULL, label = "model_complete_cases") {
  n_cc <- sum(stats::complete.cases(df[, cols, drop = FALSE]))
  assert_true(n_cc >= min_n,
              paste0(label, " has too few complete cases: ", n_cc, " < ", min_n),
              tracker,
              label)
}

check_output_files <- function(tracker, files) {
  for (f in files) {
    path <- if (grepl("^/", f)) f else file.path(tracker$output_dir, f)
    assert_true(file.exists(path), paste0("Missing output file: ", path), tracker, paste0("output_exists_", basename(path)))
    assert_true(is.finite(file.info(path)$size) && file.info(path)$size > 0,
                paste0("Empty output file: ", path),
                tracker,
                paste0("output_nonempty_", basename(path)))
  }
  invisible(TRUE)
}

finalize_sanity <- function(tracker, files = character(), filename = "sanity_report.csv") {
  if (length(files)) check_output_files(tracker, files)
  write_sanity_report(tracker, filename = filename)
  invisible(tracker)
}
