#!/usr/bin/env Rscript

# Follow-up barrier recalibration and local-barrier validation for sd3_new_af
# --------------------------------------------------------------------------
# This script:
# 1. inspects the current dataset and prior barrier-analysis outputs
# 2. rebuilds modular barrier scores
# 3. recalibrates comparative compatibility metrics so they are interpretable
# 4. tests whether stronger local barrier density predicts reduced mixing
#
# Compatibility quantities here are comparative indices, not literal biological
# probabilities, unless alpha and the barrier score are biologically calibrated.

required_packages <- c("dplyr", "ggplot2", "tibble", "tidyr", "purrr", "data.table")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them first, for example:\ninstall.packages(c(",
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

cfg <- list(
  input_rdata = NULL,
  label_rdata = NULL,
  object_name = "sd3_new_af",
  output_dir = "barrier_recalibrated_outputs",
  prior_output_dir = "barrier_model_outputs",
  alpha_default = 0.1,
  alpha_calibration_csv = NULL,
  alpha_calibration_p_col = "P_overall_mean",
  alpha_calibration_k_col = "K_genome_mean",
  alpha_grid = c(0.01, 0.05, 0.1, 0.5, 1.0),
  rolling_window_n = 25L,
  top_n_chromosomes_for_plots = 8L,
  include_pi_in_barrier = TRUE,
  include_tajima_in_barrier = TRUE,
  include_f3_in_barrier = FALSE,
  include_sweep_evidence = TRUE,
  min_label_agreement = 2L,
  agreement_rule = "gmm_rf"
)

# Keep formulas easy to edit.
barrier_terms <- list(
  base_no_delta = c("-z_recombination", "z_avg_wc_fst", "z_avg_dxy"),
  delta_term = "z_abs_delta_p",
  optional_pi = "-z_pi_average",
  optional_taj = "-z_TajimaD",
  optional_f3 = "-z_f3"
)

args <- commandArgs(trailingOnly = TRUE)
is_arg_value <- function(x) nzchar(x) && !tolower(x) %in% c("na", "null", "none", "__none__")
if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_rdata <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$object_name <- args[[2]]
if (length(args) >= 3 && nzchar(args[[3]])) cfg$output_dir <- args[[3]]
if (length(args) >= 4 && nzchar(args[[4]])) cfg$label_rdata <- args[[4]]
if (length(args) >= 5 && nzchar(args[[5]])) cfg$min_label_agreement <- as.integer(args[[5]])
if (length(args) >= 6 && is_arg_value(args[[6]])) cfg$alpha_calibration_csv <- args[[6]]
if (length(args) >= 7 && is_arg_value(args[[7]])) cfg$alpha_calibration_p_col <- args[[7]]
if (length(args) >= 8 && is_arg_value(args[[8]])) cfg$alpha_calibration_k_col <- args[[8]]
if (length(args) >= 9 && nzchar(args[[9]])) cfg$agreement_rule <- args[[9]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
sanity <- init_sanity("barrier_recalibration_followup_sd3_new_af.R", cfg$output_dir)

if (!is.null(cfg$alpha_calibration_csv) && nzchar(cfg$alpha_calibration_csv)) {
  if (!file.exists(cfg$alpha_calibration_csv)) stop("Alpha calibration CSV not found: ", cfg$alpha_calibration_csv, call. = FALSE)
  alpha_df <- read.csv(cfg$alpha_calibration_csv, stringsAsFactors = FALSE, check.names = FALSE)
  alpha_result <- estimate_alpha_from_pk(
    df = alpha_df,
    p_col = cfg$alpha_calibration_p_col,
    k_col = cfg$alpha_calibration_k_col
  )
  cfg$alpha_default <- alpha_result$alpha
  cfg$alpha_grid <- sort(unique(c(alpha_result$alpha, cfg$alpha_grid)))
  write_alpha_calibration_outputs(alpha_result, cfg$output_dir)
  message("Estimated alpha_default = ", signif(cfg$alpha_default, 6), " using ", alpha_result$preferred_model, ".")
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
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) < 2L) return(out)
  sdv <- stats::sd(x[ok], na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(out)
  out[ok] <- (x[ok] - mean(x[ok], na.rm = TRUE)) / sdv
  out
}

shift_nonnegative <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  out[ok] <- x[ok] - min(x[ok], na.rm = TRUE)
  out
}

positive_part_centered <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  mu <- mean(x[ok], na.rm = TRUE)
  out[ok] <- pmax(x[ok] - mu, 0)
  out
}

safe_empirical_r2 <- function(formula, data, response, chr_col = "Chromosome_Names", n_perm = 200L) {
  mf <- model.frame(formula, data = data, na.action = na.omit)
  if (nrow(mf) < 20L) {
    return(tibble(observed_r2 = NA_real_, null_mean_r2 = NA_real_, null_q95_r2 = NA_real_, empirical_p = NA_real_))
  }
  fit <- stats::lm(formula, data = mf)
  obs_r2 <- summary(fit)$r.squared
  chr_vals <- mf[[chr_col]]
  y <- mf[[response]]
  null_r2 <- replicate(n_perm, {
    y_perm <- y
    for (cc in unique(chr_vals)) {
      idx <- which(chr_vals == cc)
      y_perm[idx] <- sample(y_perm[idx], length(idx), replace = FALSE)
    }
    mf_perm <- mf
    mf_perm[[response]] <- y_perm
    suppressWarnings(summary(stats::lm(formula, data = mf_perm))$r.squared)
  })
  tibble(
    observed_r2 = obs_r2,
    null_mean_r2 = mean(null_r2, na.rm = TRUE),
    null_q95_r2 = stats::quantile(null_r2, probs = 0.95, na.rm = TRUE, names = FALSE),
    empirical_p = (1 + sum(null_r2 >= obs_r2, na.rm = TRUE)) / (1 + sum(is.finite(null_r2)))
  )
}

safe_loco_r2 <- function(formula, data, response, chr_col = "Chromosome_Names") {
  loco_formula <- if (chr_col %in% all.vars(formula)) {
    stats::update.formula(formula, paste(". ~ . -", chr_col))
  } else {
    formula
  }
  mf <- model.frame(stats::update.formula(loco_formula, paste(". ~ . +", chr_col)), data = data, na.action = na.omit)
  if (nrow(mf) < 20L || !(chr_col %in% names(mf))) return(tibble())
  char_cols <- names(mf)[vapply(mf, is.character, logical(1))]
  for (cc in char_cols) mf[[cc]] <- factor(mf[[cc]])
  if (!is.factor(mf[[chr_col]])) mf[[chr_col]] <- factor(mf[[chr_col]])
  chr_vals <- as.character(mf[[chr_col]])
  rows <- lapply(unique(chr_vals), function(cc) {
    test_idx <- chr_vals == cc
    train <- mf[!test_idx, , drop = FALSE]
    test <- mf[test_idx, , drop = FALSE]
    if (nrow(train) < 20L || nrow(test) < 5L) return(NULL)
    factor_cols <- names(train)[vapply(train, is.factor, logical(1))]
    for (ff in factor_cols) {
      test[[ff]] <- factor(as.character(test[[ff]]), levels = levels(train[[ff]]))
    }
    fit <- stats::lm(loco_formula, data = train)
    keep_test <- stats::complete.cases(test)
    if (sum(keep_test) < 5L) return(NULL)
    pred <- stats::predict(fit, newdata = test[keep_test, , drop = FALSE])
    obs <- test[[response]][keep_test]
    ss_res <- sum((obs - pred)^2, na.rm = TRUE)
    ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
    tibble(heldout_chromosome = cc, n_test = sum(keep_test), holdout_r2 = ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_))
  })
  bind_rows(rows)
}

safe_primary_model_diagnostics <- function(formula, data, response, model_name, output_dir, chr_col = "Chromosome_Names") {
  mf <- model.frame(formula, data = data, na.action = na.omit)
  if (nrow(mf) < 20L) return(NULL)
  fit <- stats::lm(formula, data = mf)
  sm <- summary(fit)
  coef_dt <- tibble::rownames_to_column(as.data.frame(sm$coefficients), "term")
  names(coef_dt) <- c("term", "estimate", "std_error", "t_value", "p_value")
  coef_dt$model <- model_name

  num_terms <- setdiff(names(mf), c(response, chr_col))
  num_terms <- num_terms[vapply(mf[, num_terms, drop = FALSE], is.numeric, logical(1))]
  std_coef_dt <- tibble(model = model_name, term = character(), std_estimate = numeric())
  if (length(num_terms) >= 1L) {
    std_df <- mf
    std_df[[response]] <- as.numeric(scale(std_df[[response]]))
    for (tt in num_terms) std_df[[tt]] <- as.numeric(scale(std_df[[tt]]))
    std_fit <- stats::lm(formula, data = std_df)
    std_coef <- coef(std_fit)
    std_coef_dt <- tibble(model = model_name, term = names(std_coef), std_estimate = as.numeric(std_coef))
  }

  vif_dt <- tibble(model = model_name, term = character(), vif = numeric())
  if (length(num_terms) >= 2L) {
    vif_rows <- lapply(num_terms, function(tt) {
      others <- setdiff(num_terms, tt)
      r2 <- summary(stats::lm(stats::as.formula(paste(tt, "~", paste(others, collapse = " + "))), data = mf))$r.squared
      tibble(model = model_name, term = tt, vif = 1 / max(1e-8, 1 - r2))
    })
    vif_dt <- bind_rows(vif_rows)
  }

  loco_dt <- safe_loco_r2(formula, mf, response = response, chr_col = chr_col)
  perm_dt <- safe_empirical_r2(formula, mf, response = response, chr_col = chr_col, n_perm = ifelse(nrow(mf) > 50000L, 100L, 250L))
  perm_dt$model <- model_name

  resid_df <- tibble(
    fitted = fitted(fit),
    residual = resid(fit),
    std_resid = rstandard(fit)
  )
  p_resid <- ggplot(resid_df, aes(fitted, residual)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#2b6c8f") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    theme_bw(base_size = 10) +
    labs(x = "Fitted", y = "Residual", title = paste("Residuals vs fitted:", model_name))
  p_qq <- ggplot(resid_df, aes(sample = std_resid)) +
    stat_qq(alpha = 0.25, size = 0.6, color = "#7a0177") +
    stat_qq_line(color = "#252525") +
    theme_bw(base_size = 10) +
    labs(title = paste("QQ plot:", model_name))
  ggsave(file.path(output_dir, paste0("plot_", gsub("[^A-Za-z0-9]+", "_", model_name), "_residuals.pdf")), p_resid, width = 6.8, height = 5.1)
  ggsave(file.path(output_dir, paste0("plot_", gsub("[^A-Za-z0-9]+", "_", model_name), "_qq.pdf")), p_qq, width = 6.8, height = 5.1)

  list(
    coefficients = coef_dt,
    standardized = std_coef_dt,
    vif = vif_dt,
    loco = loco_dt,
    permutation = perm_dt
  )
}

safe_cor <- function(x, y, method, analysis) {
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

safe_lm <- function(formula, data, analysis) {
  mf <- model.frame(formula, data = data, na.action = na.omit)
  if (nrow(mf) < 5L) {
    return(tibble(
      analysis = analysis, model_type = "lm", term = NA_character_,
      estimate = NA_real_, statistic = NA_real_, p_value = NA_real_, n = nrow(mf),
      r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = "Too few complete observations"
    ))
  }
  fit <- stats::lm(formula, data = data, na.action = na.omit)
  sm <- summary(fit)
  cf <- as.data.frame(sm$coefficients)
  cf$term <- rownames(cf)
  rownames(cf) <- NULL
  tibble(
    analysis = analysis, model_type = "lm", term = cf$term,
    estimate = cf$Estimate, statistic = cf$`t value`, p_value = cf$`Pr(>|t|)`,
    n = nrow(mf), r_squared = unname(sm$r.squared), adj_r_squared = unname(sm$adj.r.squared),
    aic = stats::AIC(fit), note = NA_character_
  )
}

safe_group_test <- function(y, group, analysis) {
  keep <- is.finite(y) & !is.na(group)
  y2 <- y[keep]
  g2 <- droplevels(as.factor(group[keep]))
  if (nlevels(g2) < 2L) {
    return(tibble(
      analysis = analysis, model_type = "group_test", term = NA_character_,
      estimate = NA_real_, statistic = NA_real_, p_value = NA_real_, n = length(y2),
      r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = "Fewer than 2 groups"
    ))
  }
  if (nlevels(g2) == 2L) {
    wt <- stats::wilcox.test(y2 ~ g2)
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

safe_cor_by_group <- function(data, x_col, y_col, group_col, method, analysis_prefix) {
  req <- c(x_col, y_col, group_col)
  if (!all(req %in% names(data))) return(tibble())
  split_df <- split(data, data[[group_col]])
  purrr::map_dfr(names(split_df), function(grp) {
    df_g <- split_df[[grp]]
    x <- df_g[[x_col]]
    y <- df_g[[y_col]]
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 5L) {
      return(tibble(
        analysis = paste0(analysis_prefix, " [", grp, "]"),
        model_type = paste0("cor_", method),
        term = grp,
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
      analysis = paste0(analysis_prefix, " [", grp, "]"),
      model_type = paste0("cor_", method),
      term = grp,
      estimate = unname(ct$estimate),
      statistic = unname(ct$statistic),
      p_value = ct$p.value,
      n = sum(keep),
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      aic = NA_real_,
      note = NA_character_
    )
  })
}

quantile_bin <- function(x, n = 10L) {
  out <- rep(NA_integer_, length(x))
  ok <- is.finite(x)
  if (sum(ok) < n) return(out)
  br <- unique(stats::quantile(x[ok], probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE, names = FALSE))
  if (length(br) < 2L) return(out)
  out[ok] <- cut(x[ok], breaks = br, include.lowest = TRUE, labels = FALSE)
  out
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
    tt <- sort(table(vals), decreasing = TRUE)
    if (as.integer(tt[1]) >= min_agreement) {
      out[i] <- names(tt)[1]
      agree_n[i] <- as.integer(tt[1])
    }
  }
  tibble(consensus_label = out, agreement_n = agree_n)
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

make_score_from_terms <- function(df, terms) {
  if (length(terms) == 0L) return(rep(NA_real_, nrow(df)))
  with(df, eval(parse(text = paste(terms, collapse = " + "))))
}

load_bundle <- function(input_rdata, object_name) {
  if (exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
    return(list(data = get(object_name, envir = .GlobalEnv, inherits = FALSE), env = .GlobalEnv))
  }
  if (is.null(input_rdata)) {
    stop("Provide an .rdata path or load ", object_name, " in the current R session.")
  }
  e <- new.env(parent = emptyenv())
  load(input_rdata, envir = e)
  if (!exists(object_name, envir = e, inherits = FALSE)) {
    stop("Loaded file does not contain object named ", object_name)
  }
  list(data = get(object_name, envir = e, inherits = FALSE), env = e)
}

load_optional_label_env <- function(label_rdata = NULL) {
  default_candidate <- file.path(repo_root, "outputs", "selection", "barrier_label_objects.rdata")
  legacy_candidate <- file.path(repo_root, "barrier_label_objects.rdata")
  label_path <- label_rdata
  if (is.null(label_path) && file.exists(default_candidate)) label_path <- default_candidate
  if (is.null(label_path) && file.exists(legacy_candidate)) label_path <- legacy_candidate
  if (is.null(label_path)) return(NULL)
  if (!file.exists(label_path)) {
    message("Optional label .rdata not found: ", label_path)
    return(NULL)
  }
  e <- new.env(parent = emptyenv())
  load(label_path, envir = e)
  message("Loaded optional label objects from: ", label_path)
  message("Objects in label file: ", paste(ls(e), collapse = ", "))
  e
}

label_from_rule <- function(out, rule = "gmm_rf", min_agreement = 2L) {
  rule <- tolower(rule)

  consensus_tbl <- derive_consensus_label(
    gmm_label = if ("gmm_label_h" %in% names(out)) out$gmm_label_h else rep(NA_character_, nrow(out)),
    rf_label = if ("rf_label_h" %in% names(out)) out$rf_label_h else rep(NA_character_, nrow(out)),
    cont_label_h = if ("cont_label_h2" %in% names(out)) out$cont_label_h2 else rep(NA_character_, nrow(out)),
    min_agreement = min_agreement
  )

  flag_label <- function(flag_col, label_col = "rf_label_h") {
    flag <- if (flag_col %in% names(out)) !is.na(out[[flag_col]]) & as.logical(out[[flag_col]]) else rep(FALSE, nrow(out))
    lab <- if (label_col %in% names(out)) out[[label_col]] else if ("rf_label_h" %in% names(out)) out$rf_label_h else rep(NA_character_, nrow(out))
    tibble(consensus_label = ifelse(flag, lab, NA_character_), agreement_n = ifelse(flag & !is.na(lab), 3L, 0L))
  }

  direct_label <- function(label_col) {
    lab <- if (label_col %in% names(out)) as.character(out[[label_col]]) else rep(NA_character_, nrow(out))
    tibble(consensus_label = lab, agreement_n = ifelse(!is.na(lab) & nzchar(lab), 1L, 0L))
  }

  if (rule %in% c("majority_2of3", "2of3", "legacy_2of3")) return(consensus_tbl)
  if (rule %in% c("gmm_rf", "gmmrf", "primary")) return(flag_label("gmm_rf_match", "rf_label_h"))
  if (rule %in% c("best_continuous_3of3", "best_continuous", "triple", "triple_best")) return(flag_label("triple_match", "cont_label_h2"))
  if (rule %in% c("continuous_distance_3class", "contdist3", "distance3class")) return(direct_label("continuous_distance_label_empirical_h2"))
  if (rule %in% c("gmm_rf_contdist", "gmmrf_contdist", "main_3class_agreement", "broad_agreement_selected")) return(flag_label("contdist_match", "rf_label_h"))
  if (rule %in% c("empirical_fdr_selected", "empirical_fdr", "primary_fdr_selected")) return(direct_label("empirical_fdr_selected_label"))
  if (rule %in% c("global_shift_selected")) return(direct_label("global_shift_selected_label"))
  if (rule %in% c("genomewide_fdr_selected", "genomewide_fdr_all", "genomewide_fdr_primary", "genomewide_fdr_selected_full")) return(direct_label("genomewide_fdr_selected_label"))
  if (rule %in% c("dist_margin_025", "distance_margin_025")) return(flag_label("dist_margin_025", "rf_label_h"))
  if (rule %in% c("dist_margin_050", "distance_margin_050", "primary_selected")) return(flag_label("dist_margin_050", "rf_label_h"))
  if (rule %in% c("dist_margin_100", "distance_margin_100")) return(flag_label("dist_margin_100", "rf_label_h"))
  if (rule %in% c("score_margin_025")) return(flag_label("score_margin_025", "rf_label_h"))
  if (rule %in% c("score_margin_050")) return(flag_label("score_margin_050", "rf_label_h"))
  if (rule %in% c("score_margin_100")) return(flag_label("score_margin_100", "rf_label_h"))
  if (rule %in% c("dist_top1pct", "revised_top1pct")) return(flag_label("dist_top1pct", "rf_label_h"))
  if (rule %in% c("dist_top500", "revised_top500")) return(flag_label("dist_top500", "rf_label_h"))
  if (rule %in% c("margin_025", "margin025", "triple_margin_025")) return(flag_label("triple_margin_025", "cont_label_h2"))
  if (rule %in% c("margin_050", "margin050", "triple_margin_050")) return(flag_label("triple_margin_050", "cont_label_h2"))
  if (rule %in% c("margin_100", "margin100", "triple_margin_100")) return(flag_label("triple_margin_100", "cont_label_h2"))
  if (rule %in% c("top1pct", "top_1pct", "top1", "triple_top1pct")) return(flag_label("triple_top1pct", "cont_label_top1pct_h2"))
  if (rule %in% c("tail500", "top500", "strict_top500", "triple_tail500")) return(flag_label("triple_tail500", "cont_label_tail500_h2"))

  stop(
    "Unknown agreement_rule: ", rule,
    "\nSupported rules: gmm_rf, best_continuous_3of3, continuous_distance_3class, gmm_rf_contdist, empirical_fdr_selected, dist_margin_025, dist_margin_050, dist_margin_100, score_margin_025, score_margin_050, score_margin_100, dist_top1pct, dist_top500, margin_025, margin_050, margin_100, top1pct, tail500, majority_2of3",
    call. = FALSE
  )
}

join_label_data <- function(df, env_ref, min_agreement = 2L) {
  out <- df
  have_two_labels <- function(x) sum(!is.na(c(
    find_col(x, c("gmm_label", "gmm_model_raw", "gmm_model", "gmm_class")),
    find_col(x, c("rf_label", "rf_model", "rf_class", "rf_model_raw")),
    find_col(x, c("cont_label_h", "cont_label", "continuous_label", "model_cont_raw", "cont_model_raw"))
  ))) >= 2L

  if (!have_two_labels(out) && exists("agree_dat", envir = env_ref, inherits = FALSE)) {
    lab_df <- as.data.frame(get("agree_dat", envir = env_ref, inherits = FALSE))
    join_cols <- intersect(c(
      "Chromosome_Names", ".start_bp", "gmm_label", "rf_label", "cont_label_h",
      "cont_label_tail500_h", "cont_label_top1pct_h", "cont_margin", "cont_score_quantile",
      "gmm_rf_match", "triple_match", "triple_margin_025", "triple_margin_050",
      "triple_margin_100", "triple_top1pct", "triple_tail500", "triple_sil",
      "continuous_distance_label_origin_h", "continuous_distance_label_empirical_h",
      "continuous_distance_margin_origin", "continuous_distance_margin_empirical",
      "continuous_margin_label_025", "continuous_margin_label_050", "continuous_margin_label_100",
      "contdist_match", "dist_margin_025", "dist_margin_050", "dist_margin_100",
      "score_margin_025", "score_margin_050", "score_margin_100",
      "dist_top1pct", "dist_top500",
      "empirical_fdr_selected_label", "sweep_q_value", "balancing_q_value"
    ), names(lab_df))
    if (all(c("Chromosome_Names", ".start_bp") %in% join_cols)) out <- out %>% left_join(lab_df[, join_cols, drop = FALSE], by = c("Chromosome_Names", ".start_bp"))
  }

  if (!have_two_labels(out) && exists("sd3_new_af2", envir = env_ref, inherits = FALSE)) {
    lab_df <- as.data.frame(get("sd3_new_af2", envir = env_ref, inherits = FALSE))
    join_cols <- intersect(c(
      "Chromosome_Names", ".start_bp", "gmm_label", "gmm_model_raw", "rf_label", "rf_model",
      "cont_label_h", "cont_label", "cont_label_tail500_h", "cont_label_top1pct_h",
      "cont_margin", "cont_score_quantile", "gmm_rf_match", "triple_match",
      "triple_margin_025", "triple_margin_050", "triple_margin_100",
      "triple_top1pct", "triple_tail500",
      "continuous_distance_label_origin_h", "continuous_distance_label_empirical_h",
      "continuous_distance_margin_origin", "continuous_distance_margin_empirical",
      "continuous_margin_label_025", "continuous_margin_label_050", "continuous_margin_label_100",
      "contdist_match", "dist_margin_025", "dist_margin_050", "dist_margin_100",
      "score_margin_025", "score_margin_050", "score_margin_100",
      "dist_top1pct", "dist_top500",
      "empirical_fdr_selected_label", "sweep_q_value", "balancing_q_value"
    ), names(lab_df))
    if (all(c("Chromosome_Names", ".start_bp") %in% join_cols)) out <- out %>% left_join(lab_df[, join_cols, drop = FALSE], by = c("Chromosome_Names", ".start_bp"))
  }

  gmm_col <- find_col(out, c("gmm_label", "gmm_model_raw", "gmm_model", "gmm_class"))
  rf_col <- find_col(out, c("rf_label", "rf_model", "rf_class", "rf_model_raw"))
  cont_col <- find_col(out, c("cont_label_h", "cont_label", "continuous_label", "model_cont_raw", "cont_model_raw"))

  if (!is.na(gmm_col)) out$gmm_label_h <- harmonize_model_label(out[[gmm_col]])
  if (!is.na(rf_col)) out$rf_label_h <- harmonize_model_label(out[[rf_col]])
  if (!is.na(cont_col)) out$cont_label_h2 <- harmonize_model_label(out[[cont_col]])
  if ("cont_label_tail500_h" %in% names(out)) out$cont_label_tail500_h2 <- harmonize_model_label(out$cont_label_tail500_h)
  if ("cont_label_top1pct_h" %in% names(out)) out$cont_label_top1pct_h2 <- harmonize_model_label(out$cont_label_top1pct_h)
  if ("continuous_distance_label_origin_h" %in% names(out)) out$continuous_distance_label_origin_h2 <- harmonize_model_label(out$continuous_distance_label_origin_h)
  if ("continuous_distance_label_empirical_h" %in% names(out)) out$continuous_distance_label_empirical_h2 <- harmonize_model_label(out$continuous_distance_label_empirical_h)

  consensus_tbl <- label_from_rule(out, rule = cfg$agreement_rule, min_agreement = min_agreement)
  out <- bind_cols(out, consensus_tbl)

  triple_col <- find_col(out, c("triple_match", "triple_sil"))
  if (cfg$agreement_rule %in% c("majority_2of3", "2of3", "legacy_2of3") && min_agreement >= 3L && !is.na(triple_col)) {
    triple_idx <- isTRUE(FALSE)
    triple_idx <- !is.na(out[[triple_col]]) & as.logical(out[[triple_col]])
    fallback_label <- dplyr::coalesce(
      if ("cont_label_h2" %in% names(out)) out$cont_label_h2 else rep(NA_character_, nrow(out)),
      if ("gmm_label_h" %in% names(out)) out$gmm_label_h else rep(NA_character_, nrow(out)),
      if ("rf_label_h" %in% names(out)) out$rf_label_h else rep(NA_character_, nrow(out))
    )
    out$consensus_label[triple_idx] <- fallback_label[triple_idx]
    out$agreement_n[triple_idx & !is.na(out$consensus_label)] <- 3L
    out$consensus_label[!triple_idx] <- NA_character_
    out$agreement_n[!triple_idx] <- 0L
  }

  out$sweep_group <- dplyr::case_when(
    out$consensus_label == "Geographic sweep" ~ "Sweep-like",
    out$consensus_label == "Balancing selection" ~ "Balancing-like",
    out$consensus_label == "Neutral / equilibrium" ~ "Neutral-like",
    TRUE ~ NA_character_
  )

  out
}

bundle <- load_bundle(cfg$input_rdata, cfg$object_name)
dat_raw <- as.data.frame(bundle$data)
loaded_env <- bundle$env
assert_true(nrow(dat_raw) > 0L, "Loaded dataset is empty.", sanity, "nonempty_input_dataset")
label_env <- load_optional_label_env(cfg$label_rdata)
if (!is.null(label_env)) {
  parent.env(label_env) <- loaded_env
  loaded_env <- label_env
}

message("Rows: ", nrow(dat_raw))
message("Columns: ", ncol(dat_raw))
message("Column names:")
message(paste(names(dat_raw), collapse = ", "))
message("\nStructure:")
utils::str(dat_raw)

prior_files <- c(
  file.path(cfg$prior_output_dir, "detected_columns.csv"),
  file.path(cfg$prior_output_dir, "chromosome_barrier_summary.csv"),
  file.path(cfg$prior_output_dir, "model_summary.csv")
)
message("\nPrior output files present:")
print(tibble(file = prior_files, exists = file.exists(prior_files)))

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
f3_col <- find_col(dat_raw, c("f3", "F3", "f_3"))
sweep_score_cols <- find_all_cols(dat_raw, "(^sweep_.*score$|^sweep_.*support|sweep.*prob|prob.*sweep|sweep.*posterior|posterior.*sweep|sweep_pmax|p_sweep)")

detected_columns <- tibble(
  role = c("chromosome", "start", "end", "mid", "recombination", "avg_wc_fst", "avg_dxy", "pi_average", "TajimaD", "abs_delta_p", "f3"),
  detected_column = c(chr_col, start_col, end_col, mid_col, recomb_col, fst_col, dxy_col, pi_col, taj_col, delta_p_col, f3_col)
)
if (length(sweep_score_cols) > 0L) {
  detected_columns <- bind_rows(detected_columns, tibble(role = paste0("sweep_score_", seq_along(sweep_score_cols)), detected_column = sweep_score_cols))
}
message("\nDetected columns:")
print(detected_columns)

required_cols <- c(chr_col, start_col, recomb_col, fst_col, dxy_col)
if (any(is.na(required_cols))) {
  stop("Missing one or more required columns.\n", paste(capture.output(print(detected_columns)), collapse = "\n"))
}

dat <- dat_raw %>% mutate(.row_id = row_number())
dat$Chromosome_Names <- as.character(dat[[chr_col]])
dat$.start_bp <- suppressWarnings(as.numeric(dat[[start_col]]))
dat$.end_bp <- if (!is.na(end_col)) suppressWarnings(as.numeric(dat[[end_col]])) else dat$.start_bp + 9999
dat$window_mid_bp <- if (!is.na(mid_col)) suppressWarnings(as.numeric(dat[[mid_col]])) else (dat$.start_bp + dat$.end_bp) / 2
assert_has_cols(dat, c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"), sanity, "recalibration_input_data")
assert_unique_keys(dat, c("Chromosome_Names", ".start_bp", ".end_bp"), sanity, "recalibration_input_data")
dat$recombination <- suppressWarnings(as.numeric(dat[[recomb_col]]))
dat$avg_wc_fst <- suppressWarnings(as.numeric(dat[[fst_col]]))
dat$avg_dxy <- suppressWarnings(as.numeric(dat[[dxy_col]]))
if (!is.na(pi_col)) dat$pi_average <- suppressWarnings(as.numeric(dat[[pi_col]]))
if (!is.na(taj_col)) dat$TajimaD <- suppressWarnings(as.numeric(dat[[taj_col]]))
if (!is.na(delta_p_col)) dat$abs_delta_p <- suppressWarnings(as.numeric(dat[[delta_p_col]]))
if (!is.na(f3_col)) dat$f3 <- suppressWarnings(as.numeric(dat[[f3_col]]))
if (length(sweep_score_cols) > 0L) {
  for (i in seq_along(sweep_score_cols)) dat[[paste0("sweep_score_", i)]] <- suppressWarnings(as.numeric(dat[[sweep_score_cols[[i]]]]))
}

dat <- join_label_data(dat, loaded_env, min_agreement = cfg$min_label_agreement)
message("\nAgreement-rule label counts (", cfg$agreement_rule, "):")
print(table(dat$consensus_label, useNA = "ifany"))

missingness_tbl <- tibble(
  column = names(dat_raw),
  class = vapply(dat_raw, function(x) paste(class(x), collapse = "/"), character(1)),
  n_missing = vapply(dat_raw, function(x) sum(is.na(x)), numeric(1)),
  frac_missing = vapply(dat_raw, function(x) mean(is.na(x)), numeric(1))
) %>% arrange(desc(frac_missing), column)
message("\nTop missingness columns:")
print(utils::head(missingness_tbl, 20))

for (nm in c("recombination", "avg_wc_fst", "avg_dxy", "pi_average", "TajimaD", "abs_delta_p", "f3")) {
  if (nm %in% names(dat)) dat[[paste0("z_", nm)]] <- z_score(dat[[nm]])
}

sweep_z_cols <- character(0)
if (length(sweep_score_cols) > 0L) {
  for (i in seq_along(sweep_score_cols)) {
    src <- paste0("sweep_score_", i)
    znm <- paste0("z_sweep_score_", i)
    dat[[znm]] <- z_score(dat[[src]])
    sweep_z_cols <- c(sweep_z_cols, znm)
  }
}
dat$sweep_evidence_numeric <- if (length(sweep_z_cols) > 0L) rowMeans(dat[, sweep_z_cols, drop = FALSE], na.rm = TRUE) else NA_real_
dat$sweep_evidence_numeric[!is.finite(dat$sweep_evidence_numeric)] <- NA_real_
# Treat non-sweep labeled windows as zero sweep evidence so balancing windows remain in the score.
dat$sweep_class_binary <- ifelse(is.na(dat$sweep_group), NA_real_, ifelse(dat$sweep_group == "Sweep-like", 1, 0))

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
if (cfg$include_f3_in_barrier && "z_f3" %in% names(dat)) {
  terms_no_delta <- c(terms_no_delta, barrier_terms$optional_f3)
  terms_full <- c(terms_full, barrier_terms$optional_f3)
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

dat$B_no_delta_shift <- shift_nonnegative(dat$B_no_delta_p)
dat$B_full_shift <- shift_nonnegative(dat$B_full)
dat$B_no_delta_z <- z_score(dat$B_no_delta_p)
dat$B_full_z <- z_score(dat$B_full)

dat$p_i_shift <- exp(-cfg$alpha_default * dat$B_no_delta_shift)
dat$p_i_mean <- exp(-cfg$alpha_default * positive_part_centered(dat$B_no_delta_p))
dat$p_i_full_shift <- exp(-cfg$alpha_default * dat$B_full_shift)

chrom_levels <- order_chromosomes(dat$Chromosome_Names)
dat$Chromosome_Names <- factor(dat$Chromosome_Names, levels = chrom_levels)
dt <- as.data.table(dat)
setorder(dt, Chromosome_Names, .start_bp, .row_id)
dt[, K_local_mean := data.table::frollmean(B_no_delta_shift, n = cfg$rolling_window_n, align = "center", na.rm = TRUE), by = Chromosome_Names]
dt[, K_local_sum := data.table::frollsum(B_no_delta_shift, n = cfg$rolling_window_n, align = "center", na.rm = TRUE), by = Chromosome_Names]
dt[, K_local_mean_full := data.table::frollmean(B_full_shift, n = cfg$rolling_window_n, align = "center", na.rm = TRUE), by = Chromosome_Names]
dt[, K_local_sum_full := data.table::frollsum(B_full_shift, n = cfg$rolling_window_n, align = "center", na.rm = TRUE), by = Chromosome_Names]
dt[, P_local_mean := exp(-cfg$alpha_default * K_local_mean)]
dt[, P_local_sum := exp(-cfg$alpha_default * K_local_sum)]
dt[, log10_P_local_mean := -cfg$alpha_default * K_local_mean / log(10)]
dt[, log10_P_local_sum := -cfg$alpha_default * K_local_sum / log(10)]
dat <- as.data.frame(dt)

mean_B_genome <- mean(dat$B_no_delta_p, na.rm = TRUE)
mean_B_full_genome <- mean(dat$B_full, na.rm = TRUE)
K_genome_sum <- sum(dat$B_no_delta_shift, na.rm = TRUE)
K_genome_mean <- mean(dat$B_no_delta_shift, na.rm = TRUE)
K_genome_sum_full <- sum(dat$B_full_shift, na.rm = TRUE)
K_genome_mean_full <- mean(dat$B_full_shift, na.rm = TRUE)

barrier_recalibrated_genome_summary <- tibble(
  metric = c("mean_B_genome", "mean_B_full_genome", "K_genome_sum", "K_genome_mean", "K_genome_sum_full", "K_genome_mean_full"),
  value = c(mean_B_genome, mean_B_full_genome, K_genome_sum, K_genome_mean, K_genome_sum_full, K_genome_mean_full)
)

alpha_sensitivity <- purrr::map_dfr(cfg$alpha_grid, function(alpha_val) {
  tibble(
    alpha = alpha_val,
    P_overall_sum = exp(-alpha_val * K_genome_sum),
    log_P_overall_sum = -alpha_val * K_genome_sum,
    log10_P_overall_sum = -alpha_val * K_genome_sum / log(10),
    P_overall_mean = exp(-alpha_val * K_genome_mean),
    log_P_overall_mean = -alpha_val * K_genome_mean,
    log10_P_overall_mean = -alpha_val * K_genome_mean / log(10)
  )
})

chromosome_summary_base <- dat %>%
  group_by(Chromosome_Names) %>%
  summarise(
    n_windows = n(),
    mean_B_chr = mean(B_no_delta_p, na.rm = TRUE),
    mean_B_chr_full = mean(B_full, na.rm = TRUE),
    K_chr_sum = sum(B_no_delta_shift, na.rm = TRUE),
    K_chr_mean = mean(B_no_delta_shift, na.rm = TRUE),
    K_chr_sum_full = sum(B_full_shift, na.rm = TRUE),
    K_chr_mean_full = mean(B_full_shift, na.rm = TRUE),
    mean_recombination = mean(recombination, na.rm = TRUE),
    mean_abs_delta_p = if ("abs_delta_p" %in% names(dat)) mean(abs_delta_p, na.rm = TRUE) else NA_real_,
    mean_f3 = if ("f3" %in% names(dat)) mean(f3, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

chromosome_summary_alpha <- chromosome_summary_base %>%
  mutate(
    P_chr_sum = exp(-cfg$alpha_default * K_chr_sum),
    P_chr_mean = exp(-cfg$alpha_default * K_chr_mean),
    log10_P_chr_sum = -cfg$alpha_default * K_chr_sum / log(10),
    log10_P_chr_mean = -cfg$alpha_default * K_chr_mean / log(10)
  ) %>%
  arrange(desc(mean_B_chr))

alpha_chr_rank <- purrr::map_dfr(cfg$alpha_grid, function(alpha_val) {
  chromosome_summary_base %>%
    transmute(
      alpha = alpha_val,
      Chromosome_Names,
      mean_B_chr,
      K_chr_mean,
      P_chr_mean = exp(-alpha_val * K_chr_mean),
      log10_P_chr_mean = -alpha_val * K_chr_mean / log(10)
    ) %>%
    arrange(alpha, P_chr_mean) %>%
    group_by(alpha) %>%
    mutate(rank_low_compatibility = row_number()) %>%
    ungroup()
})

barrier_alpha_sensitivity <- alpha_sensitivity %>%
  left_join(
    alpha_chr_rank %>%
      group_by(alpha) %>%
      summarise(
        chr_P_chr_mean_min = min(P_chr_mean, na.rm = TRUE),
        chr_P_chr_mean_max = max(P_chr_mean, na.rm = TRUE),
        chr_log10_mean_range = diff(range(log10_P_chr_mean, na.rm = TRUE)),
        top_low_compatibility_chr = paste(head(Chromosome_Names[order(P_chr_mean)], 3), collapse = "; "),
        .groups = "drop"
      ),
    by = "alpha"
  ) %>%
  mutate(
    interpretable_spread = dplyr::case_when(
      is.finite(chr_log10_mean_range) & chr_log10_mean_range >= 0.5 & chr_log10_mean_range <= 5 ~ "good spread",
      is.finite(chr_log10_mean_range) & chr_log10_mean_range < 0.5 ~ "compressed",
      TRUE ~ "extreme / collapsed"
    )
  )

local_summary <- dat %>%
  transmute(
    Chromosome_Names,
    .start_bp,
    .end_bp,
    window_mid_bp,
    B_no_delta_p,
    B_full,
    B_no_delta_shift,
    B_full_shift,
    p_i_shift,
    p_i_mean,
    K_local_mean,
    K_local_sum,
    K_local_mean_full,
    K_local_sum_full,
    P_local_mean,
    P_local_sum,
    log10_P_local_mean,
    log10_P_local_sum,
    abs_delta_p = if ("abs_delta_p" %in% names(dat)) abs_delta_p else NA_real_,
    f3 = if ("f3" %in% names(dat)) f3 else NA_real_,
    sweep_group
  )

validation_models <- bind_rows(
  safe_cor(dat$B_no_delta_p, dat$abs_delta_p, "pearson", "abs_delta_p ~ B_no_delta_p"),
  safe_cor(dat$B_no_delta_p, dat$abs_delta_p, "spearman", "abs_delta_p ~ B_no_delta_p"),
  safe_lm(abs_delta_p ~ B_no_delta_p, dat, "abs_delta_p ~ B_no_delta_p"),
  safe_cor(dat$K_local_mean, dat$abs_delta_p, "pearson", "abs_delta_p ~ K_local_mean"),
  safe_cor(dat$K_local_mean, dat$abs_delta_p, "spearman", "abs_delta_p ~ K_local_mean"),
  safe_lm(abs_delta_p ~ K_local_mean, dat, "abs_delta_p ~ K_local_mean"),
  safe_lm(abs_delta_p ~ B_no_delta_p + K_local_mean, dat, "abs_delta_p ~ B_no_delta_p + K_local_mean"),
  safe_cor(dat$recombination, dat$B_no_delta_p, "pearson", "recombination ~ B_no_delta_p"),
  safe_cor(dat$recombination, dat$B_no_delta_p, "spearman", "recombination ~ B_no_delta_p"),
  safe_lm(B_no_delta_p ~ recombination, dat, "B_no_delta_p ~ recombination")
)

class_adjusted_models <- tibble()

f3_summary <- tibble()
if ("f3" %in% names(dat)) {
  f3_summary <- tibble(
    metric = c("n_finite", "mean", "median", "sd", "min", "max"),
    value = c(sum(is.finite(dat$f3)), mean(dat$f3, na.rm = TRUE), median(dat$f3, na.rm = TRUE), stats::sd(dat$f3, na.rm = TRUE), min(dat$f3, na.rm = TRUE), max(dat$f3, na.rm = TRUE))
  )

  validation_models <- bind_rows(
    validation_models,
    safe_cor(dat$B_no_delta_p, dat$f3, "pearson", "f3 ~ B_no_delta_p"),
    safe_cor(dat$B_no_delta_p, dat$f3, "spearman", "f3 ~ B_no_delta_p"),
    safe_lm(f3 ~ B_no_delta_p, dat, "f3 ~ B_no_delta_p"),
    safe_cor(dat$K_local_mean, dat$f3, "pearson", "f3 ~ K_local_mean"),
    safe_cor(dat$K_local_mean, dat$f3, "spearman", "f3 ~ K_local_mean"),
    safe_lm(f3 ~ K_local_mean, dat, "f3 ~ K_local_mean"),
    safe_lm(f3 ~ B_no_delta_p + K_local_mean, dat, "f3 ~ B_no_delta_p + K_local_mean")
  )
}

if (any(!is.na(dat$sweep_group))) {
  validation_models <- bind_rows(
    validation_models,
    safe_group_test(dat$B_no_delta_p, dat$sweep_group, "B_no_delta_p by sweep/model group"),
    safe_group_test(dat$K_local_mean, dat$sweep_group, "K_local_mean by sweep/model group"),
    safe_group_test(dat$P_local_mean, dat$sweep_group, "P_local_mean by sweep/model group")
  )

  class_adjusted_models <- bind_rows(
    safe_lm(abs_delta_p ~ B_no_delta_p + sweep_group, dat, "abs_delta_p ~ B_no_delta_p + sweep_group"),
    safe_lm(abs_delta_p ~ B_no_delta_p * sweep_group, dat, "abs_delta_p ~ B_no_delta_p * sweep_group"),
    safe_lm(abs_delta_p ~ B_no_delta_p + K_local_mean + sweep_group, dat, "abs_delta_p ~ B_no_delta_p + K_local_mean + sweep_group"),
    safe_lm(abs_delta_p ~ B_no_delta_p * sweep_group + K_local_mean, dat, "abs_delta_p ~ B_no_delta_p * sweep_group + K_local_mean"),
    safe_cor_by_group(dat, "B_no_delta_p", "abs_delta_p", "sweep_group", "spearman", "abs_delta_p ~ B_no_delta_p within sweep/model group"),
    safe_cor_by_group(dat, "K_local_mean", "abs_delta_p", "sweep_group", "spearman", "abs_delta_p ~ K_local_mean within sweep/model group")
  )

  if ("f3" %in% names(dat)) {
    class_adjusted_models <- bind_rows(
      class_adjusted_models,
      safe_lm(f3 ~ B_no_delta_p + sweep_group, dat, "f3 ~ B_no_delta_p + sweep_group"),
      safe_lm(f3 ~ B_no_delta_p * sweep_group, dat, "f3 ~ B_no_delta_p * sweep_group"),
      safe_lm(f3 ~ B_no_delta_p + K_local_mean + sweep_group, dat, "f3 ~ B_no_delta_p + K_local_mean + sweep_group"),
      safe_cor_by_group(dat, "B_no_delta_p", "f3", "sweep_group", "spearman", "f3 ~ B_no_delta_p within sweep/model group"),
      safe_cor_by_group(dat, "K_local_mean", "f3", "sweep_group", "spearman", "f3 ~ K_local_mean within sweep/model group")
    )
  }
}

primary_diag_abs <- safe_primary_model_diagnostics(
  abs_delta_p ~ B_no_delta_p + K_local_mean + Chromosome_Names,
  dat,
  response = "abs_delta_p",
  model_name = "abs_delta_p_primary_with_chr",
  output_dir = cfg$output_dir
)

primary_diag_f3 <- if ("f3" %in% names(dat)) {
  safe_primary_model_diagnostics(
    f3 ~ B_no_delta_p + K_local_mean + Chromosome_Names,
    dat,
    response = "f3",
    model_name = "f3_primary_with_chr",
    output_dir = cfg$output_dir
  )
} else {
  NULL
}

diag_coef <- bind_rows(
  if (!is.null(primary_diag_abs)) primary_diag_abs$coefficients,
  if (!is.null(primary_diag_f3)) primary_diag_f3$coefficients
)
diag_std <- bind_rows(
  if (!is.null(primary_diag_abs)) primary_diag_abs$standardized,
  if (!is.null(primary_diag_f3)) primary_diag_f3$standardized
)
diag_vif <- bind_rows(
  if (!is.null(primary_diag_abs)) primary_diag_abs$vif,
  if (!is.null(primary_diag_f3)) primary_diag_f3$vif
)
diag_loco <- bind_rows(
  if (!is.null(primary_diag_abs) && nrow(primary_diag_abs$loco) > 0) mutate(primary_diag_abs$loco, model = "abs_delta_p_primary_with_chr"),
  if (!is.null(primary_diag_f3) && nrow(primary_diag_f3$loco) > 0) mutate(primary_diag_f3$loco, model = "f3_primary_with_chr")
)
diag_perm <- bind_rows(
  if (!is.null(primary_diag_abs)) primary_diag_abs$permutation,
  if (!is.null(primary_diag_f3)) primary_diag_f3$permutation
)

write.csv(detected_columns, file.path(cfg$output_dir, "detected_columns.csv"), row.names = FALSE)
write.csv(missingness_tbl, file.path(cfg$output_dir, "missingness_summary.csv"), row.names = FALSE)
write.csv(barrier_recalibrated_genome_summary, file.path(cfg$output_dir, "barrier_recalibrated_genome_summary.csv"), row.names = FALSE)
write.csv(chromosome_summary_alpha, file.path(cfg$output_dir, "barrier_recalibrated_chromosome_summary.csv"), row.names = FALSE)
write.csv(local_summary, file.path(cfg$output_dir, "barrier_recalibrated_local_summary.csv"), row.names = FALSE)
write.csv(validation_models, file.path(cfg$output_dir, "barrier_validation_models.csv"), row.names = FALSE)
if (nrow(class_adjusted_models) > 0) write.csv(class_adjusted_models, file.path(cfg$output_dir, "barrier_class_adjusted_models.csv"), row.names = FALSE)
write.csv(barrier_alpha_sensitivity, file.path(cfg$output_dir, "barrier_alpha_sensitivity.csv"), row.names = FALSE)
write.csv(alpha_chr_rank, file.path(cfg$output_dir, "barrier_alpha_chromosome_rankings.csv"), row.names = FALSE)
if (nrow(f3_summary) > 0) write.csv(f3_summary, file.path(cfg$output_dir, "f3_distribution_summary.csv"), row.names = FALSE)
if (nrow(diag_coef) > 0) write.csv(diag_coef, file.path(cfg$output_dir, "barrier_primary_model_coefficients.csv"), row.names = FALSE)
if (nrow(diag_std) > 0) write.csv(diag_std, file.path(cfg$output_dir, "barrier_primary_model_standardized_coefficients.csv"), row.names = FALSE)
if (nrow(diag_vif) > 0) write.csv(diag_vif, file.path(cfg$output_dir, "barrier_primary_model_vif.csv"), row.names = FALSE)
if (nrow(diag_loco) > 0) write.csv(diag_loco, file.path(cfg$output_dir, "barrier_primary_model_loco.csv"), row.names = FALSE)
if (nrow(diag_perm) > 0) write.csv(diag_perm, file.path(cfg$output_dir, "barrier_primary_model_permutation_r2.csv"), row.names = FALSE)

plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

p1 <- ggplot(dat, aes(B_no_delta_p, abs_delta_p)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#2b6c8f") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#b22222") +
  labs(title = "Single-window barrier score vs ancestry differentiation", x = "B_no_delta_p", y = "abs_delta_p") +
  plot_theme
ggsave(file.path(cfg$output_dir, "plot_B_no_delta_p_vs_abs_delta_p.pdf"), p1, width = 7.5, height = 5.5)

p2 <- ggplot(dat, aes(K_local_mean, abs_delta_p)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#7a3b69") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#1f4e79") +
  labs(title = "Local barrier density vs ancestry differentiation", x = "K_local_mean", y = "abs_delta_p") +
  plot_theme
ggsave(file.path(cfg$output_dir, "plot_K_local_mean_vs_abs_delta_p.pdf"), p2, width = 7.5, height = 5.5)

if ("f3" %in% names(dat)) {
  p3 <- ggplot(dat, aes(K_local_mean, f3)) +
    geom_point(alpha = 0.3, size = 0.8, color = "#238b45") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#54278f") +
    labs(title = "Local barrier density vs f3", x = "K_local_mean", y = "f3") +
    plot_theme
  ggsave(file.path(cfg$output_dir, "plot_K_local_mean_vs_f3.pdf"), p3, width = 7.5, height = 5.5)
}

p4 <- ggplot(dat, aes(recombination, B_no_delta_p)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#7a3b69") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#1f4e79") +
  labs(title = "Recombination vs barrier score", x = "Recombination", y = "B_no_delta_p") +
  plot_theme
ggsave(file.path(cfg$output_dir, "plot_recombination_vs_B_no_delta_p.pdf"), p4, width = 7.5, height = 5.5)

p5 <- chromosome_summary_alpha %>%
  mutate(Chromosome_Names = factor(Chromosome_Names, levels = chrom_levels)) %>%
  ggplot(aes(Chromosome_Names, mean_B_chr, fill = mean_B_chr)) +
  geom_col(color = "grey25") +
  scale_fill_gradient(low = "#d9f0ea", high = "#01665e") +
  labs(title = "Chromosome-level mean barrier score", x = "Chromosome", y = "mean_B_chr") +
  plot_theme
ggsave(file.path(cfg$output_dir, "plot_chromosome_mean_B.pdf"), p5, width = 9, height = 5.5)

p6 <- chromosome_summary_alpha %>%
  mutate(Chromosome_Names = factor(Chromosome_Names, levels = chrom_levels)) %>%
  ggplot(aes(Chromosome_Names, log10_P_chr_mean, fill = log10_P_chr_mean)) +
  geom_col(color = "grey25") +
  scale_fill_gradient(low = "#fee8c8", high = "#e34a33") +
  labs(title = "Chromosome-level log10 compatibility index", x = "Chromosome", y = "log10_P_chr_mean") +
  plot_theme
ggsave(file.path(cfg$output_dir, "plot_chromosome_log10_compatibility.pdf"), p6, width = 9, height = 5.5)

keep_chr_n <- min(cfg$top_n_chromosomes_for_plots, nrow(chromosome_summary_alpha))
keep_chr <- chromosome_summary_alpha %>%
  slice_max(order_by = mean_B_chr, n = keep_chr_n) %>%
  pull(Chromosome_Names) %>% as.character()

p7 <- dat %>%
  filter(as.character(Chromosome_Names) %in% keep_chr) %>%
  ggplot(aes(window_mid_bp / 1e6, K_local_mean)) +
  geom_line(color = "#8c2d04", linewidth = 0.4, na.rm = TRUE) +
  facet_wrap(~ Chromosome_Names, scales = "free_x") +
  labs(title = "Rolling local barrier density across chromosomes", x = "Position (Mb)", y = "K_local_mean") +
  plot_theme
ggsave(file.path(cfg$output_dir, "plot_rolling_K_local_mean.pdf"), p7, width = 11, height = 7.5)

p8 <- dat %>%
  filter(as.character(Chromosome_Names) %in% keep_chr) %>%
  ggplot(aes(window_mid_bp / 1e6, P_local_mean)) +
  geom_line(color = "#08519c", linewidth = 0.4, na.rm = TRUE) +
  facet_wrap(~ Chromosome_Names, scales = "free_x") +
  labs(title = "Rolling local compatibility index across chromosomes", x = "Position (Mb)", y = "P_local_mean") +
  plot_theme
ggsave(file.path(cfg$output_dir, "plot_rolling_P_local_mean.pdf"), p8, width = 11, height = 7.5)

safe_lm_summary <- function(tbl, prefix) {
  req <- c("model_type", "analysis", "r_squared", "adj_r_squared", "aic")
  if (!nrow(tbl) || !all(req %in% names(tbl))) {
    return(tibble(
      analysis = character(),
      r_squared = numeric(),
      adj_r_squared = numeric(),
      aic = numeric()
    ))
  }
  tbl %>%
    filter(model_type == "lm", grepl(prefix, analysis)) %>%
    group_by(analysis) %>%
    summarise(
      r_squared = max(r_squared, na.rm = TRUE),
      adj_r_squared = max(adj_r_squared, na.rm = TRUE),
      aic = min(aic, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(adj_r_squared))
}

best_abs_models <- safe_lm_summary(validation_models, "^abs_delta_p")
best_f3_models <- safe_lm_summary(validation_models, "^f3")
best_class_adjusted <- safe_lm_summary(class_adjusted_models, "^abs_delta_p")

rho_rec <- safe_cor(dat$recombination, dat$B_no_delta_p, "spearman", "tmp")$estimate[1]
rho_abs_b <- safe_cor(dat$B_no_delta_p, dat$abs_delta_p, "spearman", "tmp")$estimate[1]
rho_abs_k <- safe_cor(dat$K_local_mean, dat$abs_delta_p, "spearman", "tmp")$estimate[1]
rho_f3_k <- if ("f3" %in% names(dat)) safe_cor(dat$K_local_mean, dat$f3, "spearman", "tmp")$estimate[1] else NA_real_

top_barrier_chr <- chromosome_summary_alpha %>% arrange(desc(mean_B_chr)) %>% select(Chromosome_Names, mean_B_chr, K_chr_mean) %>% head(5)
low_compat_chr <- chromosome_summary_alpha %>% arrange(P_chr_mean) %>% select(Chromosome_Names, P_chr_mean, log10_P_chr_mean) %>% head(5)
best_alpha_row <- barrier_alpha_sensitivity %>%
  mutate(score = case_when(interpretable_spread == "good spread" ~ 3, interpretable_spread == "compressed" ~ 2, TRUE ~ 1)) %>%
  arrange(desc(score), abs(alpha - cfg$alpha_default)) %>%
  slice(1)

message("\n================ Recalibrated Barrier Follow-up ================")
message("Detected abs_delta_p column: ", delta_p_col)
message("Detected f3 column: ", ifelse(is.na(f3_col), "none", f3_col))
message("Detected label consensus available: ", any(!is.na(dat$consensus_label)))
message("\nHighest barrier-density chromosomes:")
print(top_barrier_chr)
message("\nLowest compatibility chromosomes:")
print(low_compat_chr)
message("\nabs_delta_p model comparison:")
print(best_abs_models)
if (nrow(best_f3_models) > 0) {
  message("\nf3 model comparison:")
  print(best_f3_models)
}
if (nrow(best_class_adjusted) > 0) {
  message("\nClass-adjusted abs_delta_p model comparison:")
  print(best_class_adjusted)
}
if (nrow(f3_summary) > 0) {
  message("\nf3 distribution summary:")
  print(f3_summary)
}
message("\nAlpha sensitivity summary:")
print(barrier_alpha_sensitivity)
message("\nInterpretation:")
best_abs_label <- if (nrow(best_abs_models)) best_abs_models$analysis[[1]] else NA_character_
message("Local K predicts abs_delta_p better than single-window B_i: ", ifelse(is.na(best_abs_label), "insufficient data", ifelse(best_abs_label != "abs_delta_p ~ B_no_delta_p", "yes", "no")))
message("Spearman(B_no_delta_p, abs_delta_p) = ", signif(rho_abs_b, 4))
message("Spearman(K_local_mean, abs_delta_p) = ", signif(rho_abs_k, 4))
message("Spearman(recombination, B_no_delta_p) = ", signif(rho_rec, 4))
if (is.finite(rho_f3_k)) {
  direction <- ifelse(rho_f3_k > 0, "higher", "lower")
  message("Observed association with f3: higher local barrier density is associated with ", direction, " f3 (Spearman = ", signif(rho_f3_k, 4), ").")
} else {
  message("Observed association with f3: not tested because f3 was unavailable.")
}
message("Most interpretable alpha by spread heuristic: alpha = ", best_alpha_row$alpha, " (", best_alpha_row$interpretable_spread, ").")
message("These results support accumulated barrier density restricting ancestry mixing if local K remains positively associated with abs_delta_p and compatibility indices decline where local K is high.")
message("\nAnalysis complete. Outputs written to: ", normalizePath(cfg$output_dir))
finalize_sanity(
  sanity,
  files = c(
    "detected_columns.csv",
    "missingness_summary.csv",
    "barrier_recalibrated_genome_summary.csv",
    "barrier_recalibrated_chromosome_summary.csv",
    "barrier_recalibrated_local_summary.csv",
    "barrier_validation_models.csv",
    "barrier_alpha_sensitivity.csv",
    "barrier_alpha_chromosome_rankings.csv",
    "plot_B_no_delta_p_vs_abs_delta_p.pdf",
    "plot_K_local_mean_vs_abs_delta_p.pdf",
    "plot_recombination_vs_B_no_delta_p.pdf",
    "plot_chromosome_mean_B.pdf",
    "plot_chromosome_log10_compatibility.pdf",
    "plot_rolling_K_local_mean.pdf",
    "plot_rolling_P_local_mean.pdf"
  )
)
