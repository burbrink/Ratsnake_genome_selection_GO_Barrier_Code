#!/usr/bin/env Rscript

# Barrier accumulation / hybrid-zone analysis for sd3_new_af
# -----------------------------------------------------------
# This script inspects the loaded dataset, detects the relevant columns,
# builds modular barrier scores, converts them into comparative reversal
# probabilities, and summarizes barrier accumulation at window, local,
# chromosome, and genome-wide scales.
#
# Important interpretation note:
# p_i = exp(-alpha * B_scaled), P_chr = exp(-alpha * K_chr_scaled),
# and P_overall = exp(-alpha * K_genome_scaled) are comparative indices
# of reversal / compatibility, not biologically calibrated probabilities,
# unless alpha and the score construction are empirically calibrated.

required_packages <- c("dplyr", "ggplot2", "tibble", "tidyr", "purrr", "data.table")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them first, for example:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(data.table)
})

options(stringsAsFactors = FALSE)
script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
source(file.path(script_dir, "analysis_sanity_checks.R"))
source(file.path(script_dir, "alpha_estimation_helpers.R"))

# ----------------------------- #
# User-editable parameters
# ----------------------------- #
cfg <- list(
  input_rdata = Sys.getenv("PIPELINE_INPUT_RDATA", unset = ""),
  object_name = "sd3_new_af",
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_model_outputs"),
  alpha = 1,
  alpha_calibration_csv = NULL,
  alpha_calibration_p_col = "P_overall_mean",
  alpha_calibration_k_col = "K_genome_mean",
  rolling_window_n = 25L,
  bin_n = 10L,
  top_n_chromosomes_for_local_plot = 8L,
  include_pi_in_barrier = TRUE,
  include_tajima_in_barrier = TRUE,
  include_sweep_evidence = TRUE,
  use_shifted_scaling = TRUE,
  min_label_agreement = 2L
)
if (!nzchar(cfg$input_rdata)) cfg$input_rdata <- file.path(repo_root, "data", "input_placeholder.rdata")

# Core formulas are intentionally kept near the top for easy swapping.
# Biological default logic here:
# - lower recombination -> higher barrier
# - higher Fst and higher dxy -> higher barrier
# - lower pi may reflect stronger linked barriers / selective depletion
# - more negative TajimaD may reflect sweep-like or skewed allele spectra
barrier_terms <- list(
  base_no_delta = c("-z_recombination", "z_avg_wc_fst", "z_avg_dxy"),
  optional_pi = "-z_pi_average",
  optional_taj = "-z_TajimaD",
  delta_term = "z_abs_delta_p"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_rdata <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$object_name <- args[[2]]
if (length(args) >= 3 && nzchar(args[[3]])) cfg$output_dir <- args[[3]]
if (length(args) >= 4 && nzchar(args[[4]])) cfg$alpha <- as.numeric(args[[4]])
if (length(args) >= 5 && nzchar(args[[5]])) cfg$alpha_calibration_csv <- args[[5]]
if (length(args) >= 6 && nzchar(args[[6]])) cfg$alpha_calibration_p_col <- args[[6]]
if (length(args) >= 7 && nzchar(args[[7]])) cfg$alpha_calibration_k_col <- args[[7]]

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
sanity <- init_sanity("barrier_score_analysis_sd3_new_af.R", cfg$output_dir)

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

pick_first <- function(x) if (length(x) == 0L) NA_character_ else x[[1]]

find_col <- function(df, candidates = NULL, regex = NULL) {
  nms <- names(df)
  if (!is.null(candidates)) {
    hit <- intersect(candidates, nms)
    if (length(hit) > 0L) return(hit[[1]])
  }
  if (!is.null(regex)) {
    hit <- nms[grepl(regex, nms, ignore.case = TRUE, perl = TRUE)]
    if (length(hit) > 0L) return(hit[[1]])
  }
  NA_character_
}

find_all_cols <- function(df, regex) {
  names(df)[grepl(regex, names(df), ignore.case = TRUE, perl = TRUE)]
}

z_score <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) < 2L) return(out)
  mu <- mean(x[ok], na.rm = TRUE)
  sdv <- stats::sd(x[ok], na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(out)
  out[ok] <- (x[ok] - mu) / sdv
  out
}

shift_scale_nonnegative <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  xmin <- min(x[ok], na.rm = TRUE)
  xshift <- x[ok] - xmin
  out[ok] <- xshift
  out
}

quantile_bin <- function(x, n = 10L) {
  out <- rep(NA_integer_, length(x))
  ok <- is.finite(x)
  if (sum(ok) < n) return(out)
  probs <- seq(0, 1, length.out = n + 1L)
  cuts <- unique(stats::quantile(x[ok], probs = probs, na.rm = TRUE, names = FALSE))
  if (length(cuts) < 2L) return(out)
  out[ok] <- cut(x[ok], breaks = cuts, include.lowest = TRUE, labels = FALSE)
  out
}

safe_cor <- function(x, y, method, label) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L) {
    return(tibble(
      analysis = label,
      model_type = paste0("cor_", method),
      term = NA_character_,
      estimate = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      n = sum(keep),
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      aic = NA_real_,
      note = "Too few complete observations"
    ))
  }
  ct <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = method))
  tibble(
    analysis = label,
    model_type = paste0("cor_", method),
    term = NA_character_,
    estimate = unname(ct$estimate),
    statistic = unname(ct$statistic),
    p_value = ct$p.value,
    n = sum(keep),
    r_squared = NA_real_,
    adj_r_squared = NA_real_,
    aic = NA_real_,
    note = NA_character_
  )
}

safe_lm <- function(formula, data, label) {
  mf <- model.frame(formula, data = data, na.action = na.omit)
  if (nrow(mf) < 5L) {
    return(tibble(
      analysis = label,
      model_type = "lm",
      term = NA_character_,
      estimate = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      n = nrow(mf),
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      aic = NA_real_,
      note = "Too few complete observations"
    ))
  }
  fit <- stats::lm(formula, data = data, na.action = na.omit)
  sm <- summary(fit)
  coef_df <- as.data.frame(sm$coefficients)
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  tibble(
    analysis = label,
    model_type = "lm",
    term = coef_df$term,
    estimate = coef_df$Estimate,
    statistic = coef_df$`t value`,
    p_value = coef_df$`Pr(>|t|)`,
    n = nrow(mf),
    r_squared = unname(sm$r.squared),
    adj_r_squared = unname(sm$adj.r.squared),
    aic = stats::AIC(fit),
    note = NA_character_
  )
}

safe_group_test <- function(y, group, label) {
  keep <- is.finite(y) & !is.na(group)
  y2 <- y[keep]
  g2 <- droplevels(as.factor(group[keep]))
  if (nlevels(g2) < 2L) {
    return(tibble(
      analysis = label,
      model_type = "group_test",
      term = NA_character_,
      estimate = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      n = length(y2),
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      aic = NA_real_,
      note = "Fewer than 2 groups"
    ))
  }

  if (nlevels(g2) == 2L) {
    wt <- stats::wilcox.test(y2 ~ g2)
    return(tibble(
      analysis = label,
      model_type = "wilcox",
      term = paste(levels(g2), collapse = " vs "),
      estimate = NA_real_,
      statistic = unname(wt$statistic),
      p_value = wt$p.value,
      n = length(y2),
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      aic = NA_real_,
      note = NA_character_
    ))
  }

  fit <- stats::aov(y2 ~ g2)
  sm <- summary(fit)[[1]]
  tibble(
    analysis = label,
    model_type = "anova",
    term = rownames(sm)[1],
    estimate = NA_real_,
    statistic = sm$`F value`[1],
    p_value = sm$`Pr(>F)`[1],
    n = length(y2),
    r_squared = NA_real_,
    adj_r_squared = NA_real_,
    aic = stats::AIC(stats::lm(y2 ~ g2)),
    note = NA_character_
  )
}

load_dataset_bundle <- function(input_rdata = NULL, object_name = "sd3_new_af") {
  if (exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
    message("Using object already present in the global environment: ", object_name)
    return(list(data = get(object_name, envir = .GlobalEnv, inherits = FALSE), env = .GlobalEnv))
  }

  if (is.null(input_rdata)) {
    stop(
      "Object ", object_name, " was not found in the global environment.\n",
      "Run this script after creating ", object_name, ", or pass an .rdata/.RData file path as argument 1."
    )
  }

  e <- new.env(parent = emptyenv())
  load(input_rdata, envir = e)
  if (!exists(object_name, envir = e, inherits = FALSE)) {
    stop("Loaded file does not contain object named ", object_name)
  }
  list(data = get(object_name, envir = e, inherits = FALSE), env = e)
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

derive_consensus_label <- function(gmm_label, rf_label, cont_label_h, min_agreement = 2L) {
  n <- length(gmm_label)
  out <- rep(NA_character_, n)
  agree_n <- rep(0L, n)

  for (i in seq_len(n)) {
    vals <- c(gmm_label[i], rf_label[i], cont_label_h[i])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(vals) == 0L) next
    tab <- sort(table(vals), decreasing = TRUE)
    best_label <- names(tab)[1]
    best_n <- as.integer(tab[1])
    if (best_n >= min_agreement) {
      out[i] <- best_label
      agree_n[i] <- best_n
    }
  }

  tibble(consensus_label = out, agreement_n = agree_n)
}

join_label_data <- function(df, env_ref, min_agreement = 2L) {
  out <- df

  gmm_col <- find_col(out, c("gmm_label", "gmm_model_raw", "gmm_model", "gmm_class"))
  rf_col <- find_col(out, c("rf_label", "rf_model", "rf_class", "rf_model_raw"))
  cont_col <- find_col(out, c("cont_label_h", "cont_label", "continuous_label", "model_cont_raw", "cont_model_raw"))
  triple_col <- find_col(out, c("triple_match", "triple_sil"))

  have_labels_now <- sum(!is.na(c(gmm_col, rf_col, cont_col))) >= 2L

  if (!have_labels_now && exists("sd3_new_af2", envir = env_ref, inherits = FALSE)) {
    lab_df <- as.data.frame(get("sd3_new_af2", envir = env_ref, inherits = FALSE))
    join_cols <- intersect(c("Chromosome_Names", ".start_bp", "gmm_label", "gmm_model_raw", "rf_label", "rf_model", "cont_label_h", "cont_label", "triple_match"), names(lab_df))
    if (all(c("Chromosome_Names", ".start_bp") %in% join_cols)) {
      out <- out %>% left_join(lab_df[, join_cols, drop = FALSE], by = c("Chromosome_Names", ".start_bp"))
    }
  }

  gmm_col <- find_col(out, c("gmm_label", "gmm_model_raw", "gmm_model", "gmm_class"))
  rf_col <- find_col(out, c("rf_label", "rf_model", "rf_class", "rf_model_raw"))
  cont_col <- find_col(out, c("cont_label_h", "cont_label", "continuous_label", "model_cont_raw", "cont_model_raw"))
  triple_col <- find_col(out, c("triple_match", "triple_sil"))
  have_labels_now <- sum(!is.na(c(gmm_col, rf_col, cont_col))) >= 2L

  if (!have_labels_now && exists("agree_dat", envir = env_ref, inherits = FALSE)) {
    lab_df <- as.data.frame(get("agree_dat", envir = env_ref, inherits = FALSE))
    join_cols <- intersect(c("Chromosome_Names", ".start_bp", "gmm_label", "rf_label", "cont_label_h", "triple_match", "triple_sil"), names(lab_df))
    if (all(c("Chromosome_Names", ".start_bp") %in% join_cols)) {
      out <- out %>% left_join(lab_df[, join_cols, drop = FALSE], by = c("Chromosome_Names", ".start_bp"))
    }
  }

  gmm_col <- find_col(out, c("gmm_label", "gmm_model_raw", "gmm_model", "gmm_class"))
  rf_col <- find_col(out, c("rf_label", "rf_model", "rf_class", "rf_model_raw"))
  cont_col <- find_col(out, c("cont_label_h", "cont_label", "continuous_label", "model_cont_raw", "cont_model_raw"))
  triple_col <- find_col(out, c("triple_match", "triple_sil"))

  if (!is.na(gmm_col)) out$gmm_label_h <- harmonize_model_label(out[[gmm_col]])
  if (!is.na(rf_col)) out$rf_label_h <- harmonize_model_label(out[[rf_col]])
  if (!is.na(cont_col)) out$cont_label_h2 <- harmonize_model_label(out[[cont_col]])
  if (!is.na(triple_col)) out$triple_match_available <- as.logical(out[[triple_col]])

  consensus_tbl <- derive_consensus_label(
    gmm_label = if ("gmm_label_h" %in% names(out)) out$gmm_label_h else rep(NA_character_, nrow(out)),
    rf_label = if ("rf_label_h" %in% names(out)) out$rf_label_h else rep(NA_character_, nrow(out)),
    cont_label_h = if ("cont_label_h2" %in% names(out)) out$cont_label_h2 else rep(NA_character_, nrow(out)),
    min_agreement = min_agreement
  )

  bind_cols(out, consensus_tbl) %>%
    mutate(
      sweep_model_consensus = case_when(
        consensus_label == "Geographic sweep" ~ "Sweep-like",
        consensus_label == "Balancing selection" ~ "Balancing-like",
        consensus_label == "Neutral / equilibrium" ~ "Neutral-like",
        TRUE ~ NA_character_
      )
    )
}

order_chromosomes <- function(chr) {
  chr_chr <- as.character(chr)
  base_num <- suppressWarnings(as.integer(gsub("^Chromosome_", "", chr_chr)))
  is_zw <- grepl("ZW", chr_chr, ignore.case = TRUE)
  ord_tbl <- tibble(chromosome = unique(chr_chr)) %>%
    mutate(
      is_zw = grepl("ZW", chromosome, ignore.case = TRUE),
      chr_num = suppressWarnings(as.integer(gsub("^Chromosome_", "", chromosome))),
      chr_sort = dplyr::case_when(
        is_zw ~ 1e9,
        is.finite(chr_num) ~ chr_num,
        TRUE ~ 1e8 + rank(chromosome, ties.method = "first")
      )
    ) %>%
    arrange(chr_sort, chromosome)
  ord_tbl$chromosome
}

make_score_from_terms <- function(df, term_strings) {
  if (length(term_strings) == 0L) return(rep(NA_real_, nrow(df)))
  expr <- paste(term_strings, collapse = " + ")
  with(df, eval(parse(text = expr)))
}

bundle <- load_dataset_bundle(cfg$input_rdata, cfg$object_name)
dat_raw <- as.data.frame(bundle$data)
loaded_env <- bundle$env
assert_true(nrow(dat_raw) > 0L, "Loaded dataset is empty.", sanity, "nonempty_input_dataset")

message("Rows: ", nrow(dat_raw))
message("Columns: ", ncol(dat_raw))
message("Column names:")
message(paste(names(dat_raw), collapse = ", "))
message("\nStructure:")
utils::str(dat_raw)

# ----------------------------- #
# Detect columns using the live dataset
# ----------------------------- #
chr_col <- find_col(dat_raw, c("Chromosome_Names", "Chromosomes", "Chromosome", "chr", "chromosome", "scaffold"))
start_col <- find_col(dat_raw, c(".start_bp", "window_start", "start_bp", "start", "win_start"))
end_col <- find_col(dat_raw, c(".end_bp", "window_end", "end_bp", "end", "win_end"))
mid_col <- find_col(dat_raw, c("midpoint", ".mid_bp", "window_mid", "mid_bp", "mid"))

recomb_col <- find_col(dat_raw, c("recombination", "rec", "rho", "r"))
fst_col <- find_col(dat_raw, c("avg_wc_fst", "wc_fst", "fst"))
dxy_col <- find_col(dat_raw, c("avg_dxy", "dxy"))
pi_col <- find_col(dat_raw, c("pi_average", "avg_pi", "pi"))
taj_col <- find_col(dat_raw, c("TajimaD", "tajima_d", "tajd"))
delta_p_col <- find_col(dat_raw, c("abs_delta_p", "mean_abs_delta_p", "delta_p_abs", "mean_delta_p_abs"))
model_label_col <- find_col(
  dat_raw,
  c("cont_label_h", "cont_label", "continuous_label", "gmm_model_raw", "rf_model", "model", "class"),
  regex = "(cont_label|continuous_label|gmm_model|rf_model|model|class)"
)
gmm_label_col_raw <- find_col(dat_raw, c("gmm_label", "gmm_model_raw", "gmm_model", "gmm_class"))
rf_label_col_raw <- find_col(dat_raw, c("rf_label", "rf_model", "rf_class", "rf_model_raw"))
cont_label_col_raw <- find_col(dat_raw, c("cont_label_h", "cont_label", "continuous_label", "model_cont_raw", "cont_model_raw"))

sweep_score_cols <- find_all_cols(
  dat_raw,
  "(^sweep_.*score$|^sweep_.*support|sweep.*prob|prob.*sweep|sweep.*posterior|posterior.*sweep|sweep_pmax|p_sweep)"
)

detected_columns <- tibble(
  role = c(
    "chromosome", "start", "end", "mid", "recombination", "avg_wc_fst",
    "avg_dxy", "pi_average", "TajimaD", "abs_delta_p", "model_label",
    "gmm_label", "rf_label", "cont_label"
  ),
  detected_column = c(
    chr_col, start_col, end_col, mid_col, recomb_col, fst_col,
    dxy_col, pi_col, taj_col, delta_p_col, model_label_col,
    gmm_label_col_raw, rf_label_col_raw, cont_label_col_raw
  )
)

if (length(sweep_score_cols) > 0L) {
  detected_columns <- bind_rows(
    detected_columns,
    tibble(role = paste0("sweep_score_", seq_along(sweep_score_cols)), detected_column = sweep_score_cols)
  )
}

message("\nDetected columns:")
print(detected_columns)

required_cols <- c(recomb_col, fst_col, dxy_col, chr_col, start_col)
if (any(is.na(required_cols))) {
  stop(
    "Missing one or more required columns.\n",
    paste(capture.output(print(detected_columns)), collapse = "\n")
  )
}

dat <- dat_raw %>%
  mutate(.row_id = row_number())

dat$chromosome <- as.character(dat[[chr_col]])
dat$window_start_bp <- suppressWarnings(as.numeric(dat[[start_col]]))
dat$window_end_bp <- if (!is.na(end_col)) suppressWarnings(as.numeric(dat[[end_col]])) else dat$window_start_bp + 9999
dat$window_mid_bp <- if (!is.na(mid_col)) suppressWarnings(as.numeric(dat[[mid_col]])) else (dat$window_start_bp + dat$window_end_bp) / 2

dat$recombination <- suppressWarnings(as.numeric(dat[[recomb_col]]))
dat$avg_wc_fst <- suppressWarnings(as.numeric(dat[[fst_col]]))
dat$avg_dxy <- suppressWarnings(as.numeric(dat[[dxy_col]]))
if (!is.na(pi_col)) dat$pi_average <- suppressWarnings(as.numeric(dat[[pi_col]]))
if (!is.na(taj_col)) dat$TajimaD <- suppressWarnings(as.numeric(dat[[taj_col]]))
if (!is.na(delta_p_col)) dat$abs_delta_p <- suppressWarnings(as.numeric(dat[[delta_p_col]]))
if (!is.na(model_label_col)) dat$model_label <- as.character(dat[[model_label_col]])
assert_has_cols(dat, c("chromosome", "window_start_bp", "window_end_bp", "window_mid_bp"), sanity, "barrier_input_data")
assert_unique_keys(dat, c("chromosome", "window_start_bp", "window_end_bp"), sanity, "barrier_input_data")

dat <- join_label_data(dat, loaded_env, min_agreement = cfg$min_label_agreement)

if ("consensus_label" %in% names(dat)) {
  message("\n2-of-3 consensus label counts:")
  print(table(dat$consensus_label, useNA = "ifany"))
  if (all(is.na(dat$consensus_label))) {
    message("No 2-of-3 consensus labels were available in the loaded objects. ",
            "To activate sweep/balancing comparisons, load an .rdata that also contains ",
            "agree_dat or sd3_new_af2, or run this script in the same R session after the master pipeline.")
  }
}

if (length(sweep_score_cols) > 0L) {
  for (i in seq_along(sweep_score_cols)) {
    dat[[paste0("sweep_score_", i)]] <- suppressWarnings(as.numeric(dat[[sweep_score_cols[[i]]]]))
  }
}

missingness_tbl <- tibble(
  column = names(dat_raw),
  class = vapply(dat_raw, function(x) paste(class(x), collapse = "/"), character(1)),
  n_missing = vapply(dat_raw, function(x) sum(is.na(x)), numeric(1)),
  frac_missing = vapply(dat_raw, function(x) mean(is.na(x)), numeric(1))
) %>%
  arrange(desc(frac_missing), column)

message("\nTop missingness columns:")
print(utils::head(missingness_tbl, 20))

# ----------------------------- #
# Standardize numeric predictors
# ----------------------------- #
standardize_if_present <- function(df, src, zname) {
  if (src %in% names(df)) df[[zname]] <- z_score(df[[src]])
  df
}

dat <- dat %>%
  standardize_if_present("recombination", "z_recombination") %>%
  standardize_if_present("avg_wc_fst", "z_avg_wc_fst") %>%
  standardize_if_present("avg_dxy", "z_avg_dxy")

if ("pi_average" %in% names(dat)) dat <- standardize_if_present(dat, "pi_average", "z_pi_average")
if ("TajimaD" %in% names(dat)) dat <- standardize_if_present(dat, "TajimaD", "z_TajimaD")
if ("abs_delta_p" %in% names(dat)) dat <- standardize_if_present(dat, "abs_delta_p", "z_abs_delta_p")

sweep_z_cols <- character(0)
if (length(sweep_score_cols) > 0L) {
  for (i in seq_along(sweep_score_cols)) {
    src <- paste0("sweep_score_", i)
    znm <- paste0("z_sweep_score_", i)
    dat <- standardize_if_present(dat, src, znm)
    if (znm %in% names(dat)) sweep_z_cols <- c(sweep_z_cols, znm)
  }
}

dat$sweep_evidence_numeric <- NA_real_
if (cfg$include_sweep_evidence && length(sweep_z_cols) > 0L) {
  dat$sweep_evidence_numeric <- rowMeans(dat[, sweep_z_cols, drop = FALSE], na.rm = TRUE)
  dat$sweep_evidence_numeric[!is.finite(dat$sweep_evidence_numeric)] <- NA_real_
}

dat$sweep_neutral_group <- NA_character_
dat$sweep_class_binary <- NA_real_
label_for_group <- NULL
if ("consensus_label" %in% names(dat) && any(!is.na(dat$consensus_label))) {
  label_for_group <- dat$consensus_label
} else if ("model_label" %in% names(dat)) {
  label_for_group <- dat$model_label
}

if (!is.null(label_for_group)) {
  lab <- tolower(label_for_group)
  is_sweep <- grepl("sweep", lab)
  is_neutral <- grepl("neutral|equilibrium|other", lab)
  is_balancing <- grepl("balancing", lab)
  dat$sweep_neutral_group <- case_when(
    is_sweep ~ "Sweep-like",
    is_balancing ~ "Balancing-like",
    is_neutral ~ "Neutral-like",
    TRUE ~ NA_character_
  )
  dat$sweep_class_binary <- case_when(
    is_sweep ~ 1,
    is_neutral ~ 0,
    TRUE ~ NA_real_
  )
}

# ----------------------------- #
# Barrier score construction
# ----------------------------- #
terms_no_delta <- barrier_terms$base_no_delta
terms_full <- c(terms_no_delta, barrier_terms$delta_term)

if (cfg$include_pi_in_barrier && "z_pi_average" %in% names(dat)) {
  terms_no_delta <- c(terms_no_delta, barrier_terms$optional_pi)
  terms_full <- c(terms_full, barrier_terms$optional_pi)
}

if (cfg$include_tajima_in_barrier && "z_TajimaD" %in% names(dat)) {
  terms_no_delta <- c(terms_no_delta, barrier_terms$optional_taj)
  terms_full <- c(terms_full, barrier_terms$optional_taj)
}

if (cfg$include_sweep_evidence) {
  if (any(is.finite(dat$sweep_evidence_numeric))) {
    terms_no_delta <- c(terms_no_delta, "sweep_evidence_numeric")
    terms_full <- c(terms_full, "sweep_evidence_numeric")
  } else if (any(is.finite(dat$sweep_class_binary))) {
    terms_no_delta <- c(terms_no_delta, "sweep_class_binary")
    terms_full <- c(terms_full, "sweep_class_binary")
  }
}

dat$B_no_delta_p <- make_score_from_terms(dat, terms_no_delta)
dat$B_full <- if ("z_abs_delta_p" %in% names(dat)) make_score_from_terms(dat, terms_full) else dat$B_no_delta_p

dat$B_no_delta_p_scaled <- shift_scale_nonnegative(dat$B_no_delta_p)
dat$B_full_scaled <- shift_scale_nonnegative(dat$B_full)
dat$B_no_delta_p_scaled_sd <- dat$B_no_delta_p_scaled / stats::sd(dat$B_no_delta_p_scaled, na.rm = TRUE)
dat$B_full_scaled_sd <- dat$B_full_scaled / stats::sd(dat$B_full_scaled, na.rm = TRUE)

dat$p_i_no_delta <- exp(-cfg$alpha * dat$B_no_delta_p_scaled)
dat$p_i_full <- exp(-cfg$alpha * dat$B_full_scaled)

# ----------------------------- #
# Rolling local barrier density
# ----------------------------- #
chromosome_levels <- order_chromosomes(dat$chromosome)
dat$chromosome <- factor(dat$chromosome, levels = chromosome_levels)

dt <- as.data.table(dat)
setorder(dt, chromosome, window_start_bp, .row_id)

dt[, K_local_mean_no_delta := data.table::frollmean(B_no_delta_p, n = cfg$rolling_window_n, align = "center", na.rm = TRUE), by = chromosome]
dt[, K_local_sum_no_delta := data.table::frollsum(B_no_delta_p, n = cfg$rolling_window_n, align = "center", na.rm = TRUE), by = chromosome]
dt[, K_local_mean_full := data.table::frollmean(B_full, n = cfg$rolling_window_n, align = "center", na.rm = TRUE), by = chromosome]
dt[, K_local_sum_full := data.table::frollsum(B_full, n = cfg$rolling_window_n, align = "center", na.rm = TRUE), by = chromosome]

dt[, K_local_mean_no_delta_scaled := shift_scale_nonnegative(K_local_mean_no_delta), by = chromosome]
dt[, K_local_sum_no_delta_scaled := shift_scale_nonnegative(K_local_sum_no_delta), by = chromosome]
dt[, K_local_mean_full_scaled := shift_scale_nonnegative(K_local_mean_full), by = chromosome]
dt[, K_local_sum_full_scaled := shift_scale_nonnegative(K_local_sum_full), by = chromosome]

dt[, P_local_mean_no_delta := exp(-cfg$alpha * K_local_mean_no_delta_scaled)]
dt[, P_local_sum_no_delta := exp(-cfg$alpha * K_local_sum_no_delta_scaled)]
dt[, P_local_mean_full := exp(-cfg$alpha * K_local_mean_full_scaled)]
dt[, P_local_sum_full := exp(-cfg$alpha * K_local_sum_full_scaled)]

dt[, log10_P_local_mean_no_delta := -cfg$alpha * K_local_mean_no_delta_scaled / log(10)]
dt[, log10_P_local_sum_no_delta := -cfg$alpha * K_local_sum_no_delta_scaled / log(10)]
dt[, log10_P_local_mean_full := -cfg$alpha * K_local_mean_full_scaled / log(10)]
dt[, log10_P_local_sum_full := -cfg$alpha * K_local_sum_full_scaled / log(10)]

dat <- as.data.frame(dt)

# ----------------------------- #
# Chromosome and genome accumulation
# ----------------------------- #
chromosome_barrier_summary <- dat %>%
  group_by(chromosome) %>%
  summarise(
    n_windows = n(),
    mean_B_no_delta_p = mean(B_no_delta_p, na.rm = TRUE),
    mean_B_full = mean(B_full, na.rm = TRUE),
    median_B_no_delta_p = median(B_no_delta_p, na.rm = TRUE),
    median_B_full = median(B_full, na.rm = TRUE),
    K_chr_no_delta_raw = sum(B_no_delta_p, na.rm = TRUE),
    K_chr_full_raw = sum(B_full, na.rm = TRUE),
    mean_recombination = mean(recombination, na.rm = TRUE),
    mean_abs_delta_p = if ("abs_delta_p" %in% names(dat)) mean(abs_delta_p, na.rm = TRUE) else NA_real_,
    mean_K_local_mean_no_delta = mean(K_local_mean_no_delta, na.rm = TRUE),
    mean_K_local_sum_no_delta = mean(K_local_sum_no_delta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    K_chr_no_delta_scaled = K_chr_no_delta_raw - min(K_chr_no_delta_raw, na.rm = TRUE),
    K_chr_full_scaled = K_chr_full_raw - min(K_chr_full_raw, na.rm = TRUE),
    P_chr_no_delta = exp(-cfg$alpha * K_chr_no_delta_scaled),
    P_chr_full = exp(-cfg$alpha * K_chr_full_scaled),
    log10_P_chr_no_delta = -cfg$alpha * K_chr_no_delta_scaled / log(10),
    log10_P_chr_full = -cfg$alpha * K_chr_full_scaled / log(10)
  ) %>%
  arrange(desc(K_chr_no_delta_raw))

K_genome_no_delta_raw <- sum(dat$B_no_delta_p, na.rm = TRUE)
K_genome_full_raw <- sum(dat$B_full, na.rm = TRUE)
K_genome_no_delta_scaled <- 0
K_genome_full_scaled <- 0

P_overall_no_delta <- exp(-cfg$alpha * K_genome_no_delta_scaled)
P_overall_full <- exp(-cfg$alpha * K_genome_full_scaled)
log_P_overall_no_delta <- -cfg$alpha * K_genome_no_delta_raw
log_P_overall_full <- -cfg$alpha * K_genome_full_raw
log10_P_overall_no_delta <- -cfg$alpha * K_genome_no_delta_raw / log(10)
log10_P_overall_full <- -cfg$alpha * K_genome_full_raw / log(10)

# The raw K sums are the accumulated barrier densities requested by the user.
# The exp(-K) values for chromosome/genome are reported as comparative indices.
# For numerical stability and interpretability, scaled versions are used locally
# and per-window; for the accumulated raw totals we report both raw log-space
# quantities and the directly exponentiated values, which may underflow to ~0.
P_overall_no_delta_raw <- exp(log_P_overall_no_delta)
P_overall_full_raw <- exp(log_P_overall_full)

dat$barrier_bin <- quantile_bin(dat$B_no_delta_p, cfg$bin_n)

local_barrier_summary <- dat %>%
  filter(is.finite(K_local_sum_no_delta) | is.finite(K_local_mean_no_delta)) %>%
  transmute(
    chromosome,
    window_start_bp,
    window_end_bp,
    window_mid_bp,
    B_no_delta_p,
    B_full,
    p_i_no_delta,
    p_i_full,
    K_local_mean_no_delta,
    K_local_sum_no_delta,
    P_local_mean_no_delta,
    P_local_sum_no_delta,
    log10_P_local_mean_no_delta,
    log10_P_local_sum_no_delta,
    abs_delta_p = if ("abs_delta_p" %in% names(dat)) abs_delta_p else NA_real_
  )

# ----------------------------- #
# Main tests
# ----------------------------- #
model_summary <- bind_rows(
  safe_cor(dat$B_no_delta_p, dat$abs_delta_p, "pearson", "B_no_delta_p predicts abs_delta_p"),
  safe_cor(dat$B_no_delta_p, dat$abs_delta_p, "spearman", "B_no_delta_p predicts abs_delta_p"),
  safe_lm(abs_delta_p ~ B_no_delta_p, dat, "abs_delta_p ~ B_no_delta_p"),
  safe_cor(dat$recombination, dat$B_no_delta_p, "pearson", "recombination vs B_no_delta_p"),
  safe_cor(dat$recombination, dat$B_no_delta_p, "spearman", "recombination vs B_no_delta_p"),
  safe_lm(B_no_delta_p ~ recombination, dat, "B_no_delta_p ~ recombination"),
  safe_lm(abs_delta_p ~ K_local_sum_no_delta, dat, "abs_delta_p ~ K_local_sum_no_delta"),
  safe_lm(abs_delta_p ~ B_no_delta_p + K_local_sum_no_delta, dat, "abs_delta_p ~ B_no_delta_p + K_local_sum_no_delta")
)

if (any(!is.na(dat$sweep_neutral_group))) {
  model_summary <- bind_rows(
    model_summary,
    safe_group_test(dat$B_no_delta_p, dat$sweep_neutral_group, "B_no_delta_p by sweep-like vs neutral-like")
  )
}

bin_summary <- dat %>%
  filter(!is.na(barrier_bin), is.finite(abs_delta_p)) %>%
  group_by(barrier_bin) %>%
  summarise(
    n = n(),
    mean_B_no_delta_p = mean(B_no_delta_p, na.rm = TRUE),
    mean_abs_delta_p = mean(abs_delta_p, na.rm = TRUE),
    median_abs_delta_p = median(abs_delta_p, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------- #
# Output tables
# ----------------------------- #
write.csv(detected_columns, file.path(cfg$output_dir, "detected_columns.csv"), row.names = FALSE)
write.csv(missingness_tbl, file.path(cfg$output_dir, "missingness_summary.csv"), row.names = FALSE)
write.csv(chromosome_barrier_summary, file.path(cfg$output_dir, "chromosome_barrier_summary.csv"), row.names = FALSE)
write.csv(local_barrier_summary, file.path(cfg$output_dir, "local_barrier_summary.csv"), row.names = FALSE)
write.csv(model_summary, file.path(cfg$output_dir, "model_summary.csv"), row.names = FALSE)
write.csv(bin_summary, file.path(cfg$output_dir, "barrier_bin_summary.csv"), row.names = FALSE)
write.csv(dat, file.path(cfg$output_dir, "sd3_new_af_barrier_scored_windows.csv"), row.names = FALSE)

# ----------------------------- #
# Plots
# ----------------------------- #
plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p1 <- ggplot(dat, aes(x = B_no_delta_p, y = abs_delta_p)) +
  geom_point(alpha = 0.3, size = 0.9, color = "#2b6c8f") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#b22222") +
  labs(
    title = "Barrier score excluding delta-p predicts ancestry differentiation",
    x = "B_no_delta_p",
    y = "abs_delta_p"
  ) +
  plot_theme

ggsave(file.path(cfg$output_dir, "plot_B_no_delta_p_vs_abs_delta_p.png"), p1, width = 7.5, height = 5.5, dpi = 300)

p2 <- ggplot(dat, aes(x = recombination, y = B_no_delta_p)) +
  geom_point(alpha = 0.3, size = 0.9, color = "#7a3b69") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#1f4e79") +
  labs(
    title = "Lower recombination is associated with higher barrier score",
    x = "Recombination",
    y = "B_no_delta_p"
  ) +
  plot_theme

ggsave(file.path(cfg$output_dir, "plot_recombination_vs_B_no_delta_p.png"), p2, width = 7.5, height = 5.5, dpi = 300)

if (any(!is.na(dat$sweep_neutral_group))) {
  p3 <- dat %>%
    filter(!is.na(sweep_neutral_group)) %>%
    ggplot(aes(x = sweep_neutral_group, y = B_no_delta_p, fill = sweep_neutral_group)) +
    geom_violin(alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = c("Sweep-like" = "#d7301f", "Balancing-like" = "#1b7837", "Neutral-like" = "#3182bd")) +
    labs(
      title = "Barrier score by 2-of-3 consensus window class",
      x = NULL,
      y = "B_no_delta_p"
    ) +
    plot_theme +
    theme(legend.position = "none")

  ggsave(file.path(cfg$output_dir, "plot_sweep_vs_neutral_B_no_delta_p.png"), p3, width = 6.5, height = 5.5, dpi = 300)
}

p4 <- chromosome_barrier_summary %>%
  mutate(chromosome = factor(chromosome, levels = chromosome_levels)) %>%
  ggplot(aes(x = chromosome, y = K_chr_no_delta_raw, fill = mean_B_no_delta_p)) +
  geom_col(color = "grey25") +
  scale_fill_gradient(low = "#d9f0ea", high = "#01665e") +
  labs(
    title = "Chromosome-level accumulated barrier density",
    x = "Chromosome",
    y = "K_chr = sum(B_no_delta_p)",
    fill = "Mean B"
  ) +
  plot_theme

ggsave(file.path(cfg$output_dir, "plot_chromosome_barrier_density.png"), p4, width = 9, height = 5.5, dpi = 300)

keep_chr_n <- min(cfg$top_n_chromosomes_for_local_plot, nrow(chromosome_barrier_summary))
keep_chr <- chromosome_barrier_summary %>%
  slice_max(order_by = K_chr_no_delta_raw, n = keep_chr_n) %>%
  pull(chromosome) %>%
  as.character()

p5 <- dat %>%
  filter(as.character(chromosome) %in% keep_chr) %>%
  ggplot(aes(x = window_mid_bp / 1e6, y = K_local_sum_no_delta)) +
  geom_line(color = "#8c2d04", linewidth = 0.4, na.rm = TRUE) +
  facet_wrap(~ chromosome, scales = "free_x") +
  labs(
    title = "Rolling local barrier density",
    subtitle = paste0("Centered rolling sum across ", cfg$rolling_window_n, " windows"),
    x = "Position (Mb)",
    y = "K_local"
  ) +
  plot_theme

ggsave(file.path(cfg$output_dir, "plot_local_barrier_density.png"), p5, width = 11, height = 7.5, dpi = 300)

p6 <- bin_summary %>%
  ggplot(aes(x = barrier_bin, y = mean_abs_delta_p)) +
  geom_line(color = "#08519c", linewidth = 0.8) +
  geom_point(color = "#08519c", size = 2) +
  labs(
    title = "Quantile-bin summary: barrier score vs ancestry differentiation",
    x = "Barrier score quantile bin",
    y = "Mean abs_delta_p"
  ) +
  plot_theme

ggsave(file.path(cfg$output_dir, "plot_barrier_bin_summary.png"), p6, width = 7, height = 5, dpi = 300)

# ----------------------------- #
# Console interpretation
# ----------------------------- #
top_chr <- chromosome_barrier_summary %>%
  select(chromosome, K_chr_no_delta_raw, P_chr_no_delta, log10_P_chr_no_delta) %>%
  head(5)

low_reversal_chr <- chromosome_barrier_summary %>%
  arrange(P_chr_no_delta) %>%
  select(chromosome, K_chr_no_delta_raw, P_chr_no_delta, log10_P_chr_no_delta) %>%
  head(5)

best_models <- model_summary %>%
  filter(model_type == "lm") %>%
  group_by(analysis) %>%
  summarise(
    n = max(n, na.rm = TRUE),
    r_squared = max(r_squared, na.rm = TRUE),
    adj_r_squared = max(adj_r_squared, na.rm = TRUE),
    aic = min(aic, na.rm = TRUE),
    .groups = "drop"
  )

message("\n================ Barrier Summary ================")
message("Highest chromosome-level barrier density (top 5):")
print(top_chr)

message("\nLowest chromosome-level reversal probability (top 5):")
print(low_reversal_chr)

message("\nModel summaries:")
print(best_models)

if ("abs_delta_p" %in% names(dat)) {
  rho_b <- safe_cor(dat$B_no_delta_p, dat$abs_delta_p, "spearman", "tmp")$estimate[1]
  rho_k <- safe_cor(dat$K_local_sum_no_delta, dat$abs_delta_p, "spearman", "tmp")$estimate[1]
  message("\nPrediction summary:")
  message("B_no_delta_p vs abs_delta_p Spearman rho = ", signif(rho_b, 4))
  message("K_local_sum_no_delta vs abs_delta_p Spearman rho = ", signif(rho_k, 4))
}

rho_rec <- safe_cor(dat$recombination, dat$B_no_delta_p, "spearman", "tmp")$estimate[1]
message("Recombination vs B_no_delta_p Spearman rho = ", signif(rho_rec, 4))

message("\nGenome-wide accumulation:")
message("K_genome_no_delta_raw = ", signif(K_genome_no_delta_raw, 6))
message("K_genome_full_raw = ", signif(K_genome_full_raw, 6))
message("P_overall_no_delta_raw = ", format(P_overall_no_delta_raw, scientific = TRUE))
message("P_overall_full_raw = ", format(P_overall_full_raw, scientific = TRUE))
message("log_P_overall_no_delta = ", signif(log_P_overall_no_delta, 6))
message("log10_P_overall_no_delta = ", signif(log10_P_overall_no_delta, 6))
message("log_P_overall_full = ", signif(log_P_overall_full, 6))
message("log10_P_overall_full = ", signif(log10_P_overall_full, 6))

message("\nInterpretation in terms of reversal probability")
message("- High B_i => low p_i.")
message("- High K => low P_overall.")
message("- If K_local is high where abs_delta_p is high, that supports the idea that accumulated barrier density restricts ancestry mixing.")

writeLines(
  c(
    "Interpretation in terms of reversal probability",
    "- High B_i => low p_i.",
    "- High K => low P_overall.",
    "- If K_local is high where abs_delta_p is high, that supports the idea that accumulated barrier density restricts ancestry mixing."
  ),
  con = file.path(cfg$output_dir, "interpretation_reversal_probability.txt")
)

message("\nAnalysis complete. Outputs written to: ", normalizePath(cfg$output_dir))
finalize_sanity(
  sanity,
  files = c(
    "detected_columns.csv",
    "missingness_summary.csv",
    "chromosome_barrier_summary.csv",
    "local_barrier_summary.csv",
    "model_summary.csv",
    "barrier_bin_summary.csv",
    "sd3_new_af_barrier_scored_windows.csv",
    "plot_B_no_delta_p_vs_abs_delta_p.png",
    "plot_recombination_vs_B_no_delta_p.png",
    "plot_chromosome_barrier_density.png",
    "plot_local_barrier_density.png",
    "plot_barrier_bin_summary.png",
    "interpretation_reversal_probability.txt"
  )
)
