#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(randomForest)
})

args <- commandArgs(trailingOnly = TRUE)

selection_output_dir <- if (length(args) >= 1L && nzchar(args[[1]])) {
  args[[1]]
} else {
  getwd()
}

script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()

label_rdata <- file.path(selection_output_dir, "barrier_label_objects.rdata")
qc_dir <- file.path(selection_output_dir, "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

sanity_script <- normalizePath(file.path(script_dir, "..", "02_barrier", "analysis_sanity_checks.R"), mustWork = FALSE)
if (file.exists(sanity_script)) source(sanity_script)
sanity <- if (exists("init_sanity")) init_sanity("selection_qc_and_fdr_helper.R", qc_dir) else NULL

stopifnot(file.exists(label_rdata))
load(label_rdata)

clean_names <- function(x) make.unique(ifelse(is.na(x) | !nzchar(x), "unnamed", x))
names(sd3_new_af2) <- clean_names(names(sd3_new_af2))
names(agree_dat) <- clean_names(names(agree_dat))

sd3 <- as.data.table(sd3_new_af2)
agd <- as.data.table(agree_dat)

safe_mean <- function(x) if (all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_sd <- function(x) if (sum(is.finite(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)
safe_quant <- function(x, p) if (all(!is.finite(x))) NA_real_ else as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 8))
safe_skew <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean(((x - m) / s)^3)
}
safe_wilcox_p <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  suppressWarnings(tryCatch(wilcox.test(x, y)$p.value, error = function(...) NA_real_))
}
mode_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}
block_permute <- function(x, block_n = 25L) {
  n <- length(x)
  if (n <= 1L) return(x)
  block_n <- max(1L, min(as.integer(block_n), n))
  idx <- split(seq_len(n), ceiling(seq_len(n) / block_n))
  idx <- idx[sample.int(length(idx))]
  x[unlist(idx, use.names = FALSE)]
}
build_empirical_p <- function(score, chr, n_iter = 1000L, block_n = 25L) {
  dt <- data.table(score = as.numeric(score), chr = as.character(chr))
  dt[, empirical_p := NA_real_]
  dt[, null_n := 0L]
  for (cc in unique(dt$chr)) {
    sub_idx <- which(dt$chr == cc & is.finite(dt$score))
    if (!length(sub_idx)) next
    sc <- dt$score[sub_idx]
    null_pool <- numeric(length(sc) * n_iter)
    pos <- 1L
    for (ii in seq_len(n_iter)) {
      perm <- block_permute(sc, block_n = block_n)
      null_pool[pos:(pos + length(sc) - 1L)] <- perm
      pos <- pos + length(sc)
    }
    null_pool <- sort(null_pool)
    null_n <- length(null_pool)
    n_ge <- null_n - findInterval(sc, null_pool, left.open = FALSE, rightmost.closed = FALSE)
    dt$empirical_p[sub_idx] <- (1 + n_ge) / (1 + null_n)
    dt$null_n[sub_idx] <- length(null_pool)
  }
  dt$empirical_p
}
build_reference_tail_p <- function(score, chr, reference_mask, min_ref_n = 100L) {
  dt <- data.table(
    score = as.numeric(score),
    chr = as.character(chr),
    reference_mask = as.logical(reference_mask)
  )
  dt[, empirical_p := NA_real_]
  global_ref <- sort(dt[reference_mask == TRUE & is.finite(score), score])
  global_n <- length(global_ref)
  for (cc in unique(dt$chr)) {
    obs_idx <- which(dt$chr == cc & is.finite(dt$score))
    if (!length(obs_idx)) next
    ref_chr <- sort(dt[chr == cc & reference_mask == TRUE & is.finite(score), score])
    ref_use <- if (length(ref_chr) >= min_ref_n) ref_chr else global_ref
    ref_n <- length(ref_use)
    if (!ref_n) next
    obs_scores <- dt$score[obs_idx]
    n_ge <- ref_n - findInterval(obs_scores, ref_use, left.open = FALSE, rightmost.closed = FALSE)
    dt$empirical_p[obs_idx] <- (1 + n_ge) / (1 + ref_n)
  }
  dt$empirical_p
}
ordered_chr_levels <- function(x) {
  x <- unique(as.character(x))
  nn <- suppressWarnings(as.integer(sub("^Chromosome_([0-9]+)$", "\\1", x)))
  out <- x[order(ifelse(is.na(nn), Inf, nn), x == "Chromosome_ZW", x)]
  if ("Chromosome_ZW" %in% out) out <- c(setdiff(out, "Chromosome_ZW"), "Chromosome_ZW")
  out
}
permute_chr_summary_p <- function(chr_summary, value_col, chr_col = "Chromosome_Names", zw_name = "Chromosome_ZW", n_iter = 5000L) {
  x <- copy(chr_summary)
  x <- x[is.finite(get(value_col))]
  if (!nrow(x) || !(zw_name %chin% x[[chr_col]]) || uniqueN(x[[chr_col]]) < 3L) return(NA_real_)
  obs_zw <- x[get(chr_col) == zw_name, get(value_col)][1]
  obs_auto <- x[get(chr_col) != zw_name, mean(get(value_col), na.rm = TRUE)]
  obs_diff <- obs_zw - obs_auto
  vals <- x[[value_col]]
  null_diff <- replicate(n_iter, {
    perm <- sample(vals, length(vals), replace = FALSE)
    zw_idx <- which(x[[chr_col]] == zw_name)[1]
    perm[zw_idx] - mean(perm[-zw_idx], na.rm = TRUE)
  })
  (1 + sum(abs(null_diff) >= abs(obs_diff), na.rm = TRUE)) / (1 + length(null_diff))
}

raw_vars <- c("avg_dxy", "avg_wc_fst", "pi_average", "f3", "TajimaD", "recombination", "abs_delta_p", "B_no_delta_p", "K_local_mean")
raw_vars_present <- raw_vars[raw_vars %chin% names(sd3)]
raw_long <- melt(
  sd3[, ..raw_vars_present],
  measure.vars = raw_vars_present,
  variable.name = "variable",
  value.name = "value"
)

qc_input_summary <- rbindlist(lapply(raw_vars_present, function(vv) {
  x <- as.numeric(sd3[[vv]])
  data.table(
    variable = vv,
    n = length(x),
    missing = sum(!is.finite(x)),
    mean = safe_mean(x),
    sd = safe_sd(x),
    min = safe_quant(x, 0),
    q01 = safe_quant(x, 0.01),
    median = safe_quant(x, 0.5),
    q99 = safe_quant(x, 0.99),
    max = safe_quant(x, 1)
  )
}), use.names = TRUE, fill = TRUE)

coord_checks <- data.table(
  variable = c("missing_chromosome", "missing_start", "missing_end", "duplicated_windows", "nonpositive_width_windows"),
  n = c(
    sum(is.na(sd3$Chromosome_Names) | !nzchar(as.character(sd3$Chromosome_Names))),
    sum(!is.finite(sd3$.start_bp)),
    sum(!is.finite(sd3$.end_bp)),
    sum(duplicated(sd3[, .(Chromosome_Names, .start_bp, .end_bp)])),
    sum(is.finite(sd3$.start_bp) & is.finite(sd3$.end_bp) & (sd3$.end_bp - sd3$.start_bp + 1) <= 0)
  )
)
fwrite(qc_input_summary, file.path(qc_dir, "QC_input_summary.tsv"), sep = "\t")
fwrite(coord_checks, file.path(qc_dir, "QC_input_coordinate_checks.tsv"), sep = "\t")
fwrite(sd3[, .N, by = Chromosome_Names][order(Chromosome_Names)], file.path(qc_dir, "QC_windows_per_chromosome.tsv"), sep = "\t")

if (length(raw_vars_present) >= 2L) {
  cor_mat <- suppressWarnings(cor(sd3[, ..raw_vars_present], use = "pairwise.complete.obs", method = "spearman"))
  qc_cor <- as.data.table(as.table(cor_mat))
  corr_name_map <- c(Var1 = "var1", Var2 = "var2", Freq = "spearman_rho", N = "spearman_rho")
  setnames(qc_cor, old = intersect(names(corr_name_map), names(qc_cor)), new = unname(corr_name_map[intersect(names(corr_name_map), names(qc_cor))]))
} else {
  qc_cor <- data.table(var1 = character(), var2 = character(), spearman_rho = numeric())
}
fwrite(qc_cor, file.path(qc_dir, "QC_correlations.tsv"), sep = "\t")

outlier_rows <- rbindlist(lapply(raw_vars_present, function(vv) {
  x <- as.numeric(sd3[[vv]])
  q <- quantile(x, probs = c(0.001, 0.01, 0.99, 0.999), na.rm = TRUE, names = FALSE, type = 8)
  data.table(
    variable = vv,
    threshold = c("below_0.1pct", "below_1pct", "above_99pct", "above_99.9pct"),
    cutoff = q,
    n_windows = c(sum(x <= q[1], na.rm = TRUE), sum(x <= q[2], na.rm = TRUE), sum(x >= q[3], na.rm = TRUE), sum(x >= q[4], na.rm = TRUE))
  )
}), use.names = TRUE, fill = TRUE)
fwrite(outlier_rows, file.path(qc_dir, "QC_input_outlier_report.tsv"), sep = "\t")

if (nrow(raw_long)) {
  p_raw_hist <- ggplot(raw_long[is.finite(value)], aes(value)) +
    geom_histogram(bins = 60, fill = "#729fcf", color = "grey25", linewidth = 0.2) +
    facet_wrap(~ variable, scales = "free", ncol = 3) +
    theme_bw(base_size = 10) +
    labs(x = "Raw value", y = "Windows", title = "Raw variable histograms")
  ggsave(file.path(qc_dir, "raw_variable_histograms.pdf"), p_raw_hist, width = 11, height = 8.5)
}

z_vars <- c("z_dxy", "z_fst", "z_pi", "z_f3", "z_taj", "z_rec")
z_present <- z_vars[z_vars %chin% names(agd)]
qc_ranknorm_summary <- rbindlist(lapply(z_present, function(vv) {
  x <- as.numeric(agd[[vv]])
  uq <- unique(x[is.finite(x)])
  data.table(
    variable = vv,
    n = length(x),
    missing = sum(!is.finite(x)),
    mean = safe_mean(x),
    sd = safe_sd(x),
    min = safe_quant(x, 0),
    max = safe_quant(x, 1),
    skew = safe_skew(x),
    n_unique = length(uq),
    tied_fraction = if (sum(is.finite(x)) > 0) 1 - (length(uq) / sum(is.finite(x))) else NA_real_
  )
}), use.names = TRUE, fill = TRUE)
fwrite(qc_ranknorm_summary, file.path(qc_dir, "QC_ranknorm_summary.tsv"), sep = "\t")

qq_dt <- rbindlist(lapply(z_present, function(vv) {
  x <- sort(as.numeric(agd[[vv]][is.finite(agd[[vv]])]))
  if (!length(x)) return(NULL)
  data.table(variable = vv, theoretical = qnorm(ppoints(length(x))), sample = x)
}), use.names = TRUE, fill = TRUE)
if (nrow(qq_dt)) {
  p_qq <- ggplot(qq_dt, aes(theoretical, sample)) +
    geom_point(alpha = 0.25, size = 0.4, color = "#2b8cbe") +
    geom_abline(intercept = 0, slope = 1, color = "#b30000", linewidth = 0.5) +
    facet_wrap(~ variable, scales = "free", ncol = 3) +
    theme_bw(base_size = 10) +
    labs(x = "Theoretical quantile", y = "Observed quantile", title = "QQ plots for rank-normalized variables")
  ggsave(file.path(qc_dir, "z_variable_QQplots.pdf"), p_qq, width = 11, height = 8.5)
}

gmm_model_selection_path <- file.path(selection_output_dir, "GMM_Gselection_StabilitySummary.tsv")
if (file.exists(gmm_model_selection_path)) {
  gsel <- as.data.table(fread(gmm_model_selection_path))
} else {
  gsel <- data.table(G = integer(), BIC = numeric(), sil_mean = numeric(), ari_mean = numeric(), ari_sd = numeric())
}
fwrite(gsel, file.path(qc_dir, "QC_GMM_model_selection.tsv"), sep = "\t")
if (nrow(gsel)) {
  gsel_long <- melt(gsel, id.vars = "G", measure.vars = intersect(c("BIC", "sil_mean", "ari_mean"), names(gsel)))
  p_gmm <- ggplot(gsel_long, aes(G, value, color = variable)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.4) +
    facet_wrap(~ variable, scales = "free_y", ncol = 1) +
    theme_bw(base_size = 10) +
    labs(x = "G", y = "Value", color = NULL, title = "GMM model selection diagnostics")
  ggsave(file.path(qc_dir, "GMM_BIC_silhouette_ARI.pdf"), p_gmm, width = 7.4, height = 8.4)
}

gmm_cent <- agd[
  !is.na(gmm_label),
  .(
    n_windows = .N,
    cluster_prop = .N / nrow(agd),
    mean_z_dxy = safe_mean(z_dxy),
    mean_z_fst = safe_mean(z_fst),
    mean_z_pi = safe_mean(z_pi),
    mean_z_f3 = safe_mean(z_f3),
    mean_z_taj = safe_mean(z_taj),
    dominant_chromosome = mode_chr(Chromosome_Names),
    dominant_chromosome_fraction = max(prop.table(table(Chromosome_Names)))
  ),
  by = .(gmm_label)
]
fwrite(gmm_cent, file.path(qc_dir, "QC_GMM_centroids.tsv"), sep = "\t")
if (nrow(gmm_cent)) {
  heat_dt <- melt(gmm_cent[, .(gmm_label, mean_z_dxy, mean_z_fst, mean_z_pi, mean_z_f3, mean_z_taj)],
                  id.vars = "gmm_label", variable.name = "metric", value.name = "mean_z")
  p_heat <- ggplot(heat_dt, aes(metric, gmm_label, fill = mean_z)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    theme_bw(base_size = 10) +
    labs(x = NULL, y = "GMM label", fill = "Mean z", title = "GMM centroid heatmap")
  ggsave(file.path(qc_dir, "GMM_centroid_heatmap.pdf"), p_heat, width = 7.8, height = 4.8)
}

rf_conf <- agd[, .N, by = .(gmm_label, rf_label)][order(gmm_label, rf_label)]
rf_agree <- data.table(
  metric = c("overall_agreement"),
  value = c(mean(agd$gmm_label == agd$rf_label, na.rm = TRUE))
)
fwrite(rbind(rf_conf, data.table(gmm_label = "SUMMARY", rf_label = rf_agree$metric, N = rf_agree$value), fill = TRUE),
       file.path(qc_dir, "QC_RF_confusion_matrix.tsv"), sep = "\t")

rf_qc_rows <- data.table()
rf_imp <- data.table()
rf_prob_summary <- data.table()
rf_fit_ok <- FALSE
rf_dat <- agd[complete.cases(agd[, ..z_present]) & !is.na(gmm_label), c("gmm_label", z_present), with = FALSE]
if (nrow(rf_dat) >= 100L && uniqueN(rf_dat$gmm_label) >= 2L) {
  set.seed(4201)
  rf_fit <- randomForest(as.factor(gmm_label) ~ ., data = as.data.frame(rf_dat), ntree = 400, importance = TRUE)
  rf_fit_ok <- TRUE
  rf_qc_rows <- data.table(
    metric = c("rf_qc_refit_used", "oob_error"),
    value = c(1, rf_fit$err.rate[rf_fit$ntree, "OOB"])
  )
  if (!is.null(rf_fit$confusion)) {
    conf_dt <- as.data.table(rf_fit$confusion, keep.rownames = "class")
    fwrite(conf_dt, file.path(qc_dir, "QC_RF_oob_confusion_matrix.tsv"), sep = "\t")
  }
  imp <- importance(rf_fit)
  rf_imp <- data.table(variable = rownames(imp), imp)
  prob <- predict(rf_fit, type = "prob")
  rf_prob_summary <- data.table(
    class = levels(as.factor(rf_dat$gmm_label)),
    mean_max_probability = tapply(apply(prob, 1, max), rf_dat$gmm_label, mean, na.rm = TRUE),
    frac_lt_0_6 = tapply(apply(prob, 1, max) < 0.6, rf_dat$gmm_label, mean, na.rm = TRUE),
    frac_lt_0_7 = tapply(apply(prob, 1, max) < 0.7, rf_dat$gmm_label, mean, na.rm = TRUE),
    frac_lt_0_8 = tapply(apply(prob, 1, max) < 0.8, rf_dat$gmm_label, mean, na.rm = TRUE)
  )
  imp_plot <- data.table(variable = rownames(imp), importance = imp[, 1])
  p_imp <- ggplot(imp_plot, aes(reorder(variable, importance), importance)) +
    geom_col(fill = "#6a51a3") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(x = NULL, y = colnames(imp)[1], title = "RF variable importance (QC refit)")
  ggsave(file.path(qc_dir, "RF_variable_importance.pdf"), p_imp, width = 7.2, height = 4.8)
  prob_long <- data.table(class = rf_dat$gmm_label, max_prob = apply(prob, 1, max))
  p_prob <- ggplot(prob_long, aes(max_prob, fill = class)) +
    geom_histogram(bins = 45, alpha = 0.6, position = "identity") +
    theme_bw(base_size = 10) +
    labs(x = "RF max class probability", y = "Windows", fill = "Class", title = "RF max-probability distributions (QC refit)")
  ggsave(file.path(qc_dir, "RF_probability_distributions.pdf"), p_prob, width = 8.4, height = 5.2)
}
if (!nrow(rf_qc_rows)) rf_qc_rows <- data.table(metric = "rf_qc_refit_used", value = 0)
fwrite(rf_qc_rows, file.path(qc_dir, "QC_RF_model_notes.tsv"), sep = "\t")
if (!nrow(rf_prob_summary)) rf_prob_summary <- data.table(class = "not_available", mean_max_probability = NA_real_, frac_lt_0_6 = NA_real_, frac_lt_0_7 = NA_real_, frac_lt_0_8 = NA_real_)
fwrite(rf_prob_summary, file.path(qc_dir, "QC_RF_probability_summary.tsv"), sep = "\t")
if (nrow(rf_imp)) fwrite(rf_imp, file.path(qc_dir, "QC_RF_variable_importance.tsv"), sep = "\t")

cont_cent_path <- file.path(selection_output_dir, "continuous_distance_centroids.tsv")
if (file.exists(cont_cent_path)) {
  cont_cent <- as.data.table(fread(cont_cent_path))
} else {
  cont_cent <- data.table()
}
fwrite(cont_cent, file.path(qc_dir, "QC_continuous_centroids.tsv"), sep = "\t")

cont_conf <- rbindlist(list(
  agd[, .(comparison = "GMM_vs_continuous", label_a = gmm_label, label_b = continuous_distance_label_empirical_h)][, .N, by = .(comparison, label_a, label_b)],
  agd[, .(comparison = "RF_vs_continuous", label_a = rf_label, label_b = continuous_distance_label_empirical_h)][, .N, by = .(comparison, label_a, label_b)],
  agd[, .(comparison = "GMMRF_vs_continuous", label_a = gmm_rf_label, label_b = continuous_distance_label_empirical_h)][, .N, by = .(comparison, label_a, label_b)]
), use.names = TRUE, fill = TRUE)
fwrite(cont_conf, file.path(qc_dir, "QC_continuous_confusion_matrices.tsv"), sep = "\t")

cont_diag <- agd[, .(
  continuous_distance_label_empirical_h,
  d_neutral_empirical, d_balancing, d_sweep,
  continuous_distance_margin_empirical,
  distance_ratio = fifelse(
    pmax(d_neutral_empirical, d_balancing, d_sweep, na.rm = TRUE) > 0,
    pmin(d_neutral_empirical, d_balancing, d_sweep, na.rm = TRUE) /
      pmax(
        pmin(pmax(d_neutral_empirical, d_balancing), pmax(d_neutral_empirical, d_sweep), na.rm = TRUE),
        pmin(pmax(d_balancing, d_sweep), pmax(d_neutral_empirical, d_sweep), na.rm = TRUE),
        na.rm = TRUE
      ),
    NA_real_
  )
)]
cont_long <- melt(cont_diag, id.vars = "continuous_distance_label_empirical_h", variable.name = "metric", value.name = "value")
if (nrow(cont_long)) {
  p_cont <- ggplot(cont_long[is.finite(value)], aes(value, fill = continuous_distance_label_empirical_h)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    facet_wrap(~ metric, scales = "free", ncol = 2) +
    theme_bw(base_size = 10) +
    labs(x = "Value", y = "Windows", fill = "Continuous label", title = "Continuous distance diagnostics")
  ggsave(file.path(qc_dir, "continuous_distance_margin_histograms.pdf"), p_cont, width = 10.5, height = 8)
}

pca_dt <- agd[complete.cases(agd[, ..z_present])]
if (nrow(pca_dt) >= 10L) {
  pca_obj <- prcomp(as.matrix(pca_dt[, ..z_present]), center = TRUE, scale. = FALSE)
  pca_plot_dt <- data.table(
    PC1 = pca_obj$x[, 1],
    PC2 = pca_obj$x[, 2],
    label = pca_dt$continuous_distance_label_empirical_h,
    margin = pca_dt$continuous_distance_margin_empirical
  )
  p_pca1 <- ggplot(pca_plot_dt, aes(PC1, PC2, color = label)) +
    geom_point(alpha = 0.35, size = 0.5) +
    theme_bw(base_size = 10) +
    labs(title = "PCA colored by continuous label", color = "Label")
  p_pca2 <- ggplot(pca_plot_dt, aes(PC1, PC2, color = margin)) +
    geom_point(alpha = 0.35, size = 0.5) +
    scale_color_viridis_c() +
    theme_bw(base_size = 10) +
    labs(title = "PCA colored by distance margin", color = "Margin")
  pdf(file.path(qc_dir, "PCA_or_UMAP_labels.pdf"), width = 11, height = 5.2)
  print(p_pca1)
  print(p_pca2)
  dev.off()
}

score_dir_path <- file.path(selection_output_dir, "revised_selection_tier_counts.tsv")
tier_counts <- if (file.exists(score_dir_path)) as.data.table(fread(score_dir_path)) else data.table()

n_iter <- if (identical(Sys.getenv("PIPELINE_FAST_TEST"), "1")) 100L else 1000L
block_n <- 25L

agd[, sweep_strength := d_neutral_empirical - d_sweep]
agd[, balancing_strength := d_neutral_empirical - d_balancing]
agd[, sweep_direction_score := z_fst + z_dxy - z_pi - z_f3]
agd[, balancing_direction_score := z_pi + z_taj - z_fst]
agd[, primary_sweep_candidate := gmm_rf_match & contdist_match & rf_label == "Geographic sweep" & sweep_strength > 0 & z_fst > 0 & z_pi < 0]
agd[, primary_balancing_candidate := gmm_rf_match & contdist_match & rf_label == "Balancing selection" & balancing_strength > 0 & z_pi > 0 & z_taj > 0]
agd[, neutral_reference_mask := gmm_rf_match & contdist_match & rf_label == "Neutral / equilibrium"]
agd[, sweep_neutralref_p := build_reference_tail_p(sweep_strength, Chromosome_Names, neutral_reference_mask)]
agd[, balancing_neutralref_p := build_reference_tail_p(balancing_strength, Chromosome_Names, neutral_reference_mask)]
agd[, sweep_empirical_p := build_empirical_p(sweep_strength, Chromosome_Names, n_iter = n_iter, block_n = block_n)]
agd[, balancing_empirical_p := build_empirical_p(balancing_strength, Chromosome_Names, n_iter = n_iter, block_n = block_n)]
agd[, sweep_empirical_p_genomewide := build_empirical_p(sweep_strength, rep("Genomewide", .N), n_iter = n_iter, block_n = block_n)]
agd[, balancing_empirical_p_genomewide := build_empirical_p(balancing_strength, rep("Genomewide", .N), n_iter = n_iter, block_n = block_n)]
agd[, sweep_q_value := NA_real_]
agd[, balancing_q_value := NA_real_]
if (sum(agd$primary_sweep_candidate, na.rm = TRUE) > 0) {
  agd[primary_sweep_candidate == TRUE, sweep_q_value := p.adjust(sweep_empirical_p, method = "BH"), by = Chromosome_Names]
}
if (sum(agd$primary_balancing_candidate, na.rm = TRUE) > 0) {
  agd[primary_balancing_candidate == TRUE, balancing_q_value := p.adjust(balancing_empirical_p, method = "BH"), by = Chromosome_Names]
}
agd[, sweep_q_value_genomewide := NA_real_]
agd[, balancing_q_value_genomewide := NA_real_]
if (sum(agd$primary_sweep_candidate, na.rm = TRUE) > 0) {
  idx <- which(agd$primary_sweep_candidate %in% TRUE)
  agd$sweep_q_value_genomewide[idx] <- p.adjust(agd$sweep_empirical_p_genomewide[idx], method = "BH")
}
if (sum(agd$primary_balancing_candidate, na.rm = TRUE) > 0) {
  idx <- which(agd$primary_balancing_candidate %in% TRUE)
  agd$balancing_q_value_genomewide[idx] <- p.adjust(agd$balancing_empirical_p_genomewide[idx], method = "BH")
}
agd[, sweep_q_value_chr_aware := sweep_q_value]
agd[, balancing_q_value_chr_aware := balancing_q_value]
agd[, empirical_fdr_selected_label := fifelse(
  primary_sweep_candidate == TRUE & is.finite(sweep_q_value) & sweep_q_value <= 0.05,
  "Geographic sweep",
  fifelse(
    primary_balancing_candidate == TRUE & is.finite(balancing_q_value) & balancing_q_value <= 0.05,
    "Balancing selection",
    "No consensus"
  )
)]
agd[, empirical_fdr_selected := empirical_fdr_selected_label %chin% c("Geographic sweep", "Balancing selection")]
agd[, genomewide_fdr_selected_label := fifelse(
  primary_sweep_candidate == TRUE & is.finite(sweep_q_value_genomewide) & sweep_q_value_genomewide <= 0.05,
  "Geographic sweep",
  fifelse(
    primary_balancing_candidate == TRUE & is.finite(balancing_q_value_genomewide) & balancing_q_value_genomewide <= 0.05,
    "Balancing selection",
    "No consensus"
  )
)]
agd[, global_shift_selected_label := fifelse(
  genomewide_fdr_selected_label == "Geographic sweep" & !(empirical_fdr_selected_label == "Geographic sweep"),
  "Geographic sweep",
  fifelse(
    genomewide_fdr_selected_label == "Balancing selection" & !(empirical_fdr_selected_label == "Balancing selection"),
    "Balancing selection",
    "No consensus"
  )
)]
if (!"gmm_rf_contdist_label" %in% names(agd)) {
  agd[, gmm_rf_contdist_label := fifelse(gmm_rf_match & contdist_match, rf_label, "No consensus")]
}
agd[, broad_agreement_selected_label := fifelse(
  gmm_rf_contdist_label %chin% c("Geographic sweep", "Balancing selection"),
  gmm_rf_contdist_label,
  "No consensus"
)]
agd[, interpretive_selection_class := fifelse(
  empirical_fdr_selected_label %chin% c("Geographic sweep", "Balancing selection"),
  "local_outlier",
  fifelse(
    global_shift_selected_label %chin% c("Geographic sweep", "Balancing selection"),
    "global_shift",
    fifelse(
      broad_agreement_selected_label %chin% c("Geographic sweep", "Balancing selection"),
      "broad_agreement",
      "no_consensus_or_neutral"
    )
  )
)]

fdr_out <- agd[, .(
  Chromosome_Names, .start_bp, .end_bp,
  gmm_label, rf_label, continuous_distance_label_empirical_h,
  gmm_rf_match, contdist_match,
  sweep_strength, balancing_strength,
  sweep_direction_score, balancing_direction_score,
  sweep_neutralref_p,
  sweep_empirical_p, sweep_q_value_chr_aware = sweep_q_value,
  sweep_empirical_p_genomewide, sweep_q_value_genomewide,
  balancing_neutralref_p,
  balancing_empirical_p, balancing_q_value_chr_aware = balancing_q_value,
  balancing_empirical_p_genomewide, balancing_q_value_genomewide,
  primary_sweep_candidate, primary_balancing_candidate,
  empirical_fdr_selected_label,
  genomewide_fdr_selected_label,
  global_shift_selected_label,
  broad_agreement_selected_label,
  interpretive_selection_class
)]
fwrite(fdr_out, file.path(qc_dir, "QC_empirical_FDR_selected_windows.tsv"), sep = "\t")

key_cols <- c("Chromosome_Names", ".start_bp")
setkeyv(agd, key_cols)
setkeyv(sd3, key_cols)
sd3[agd, empirical_fdr_selected_label := i.empirical_fdr_selected_label]
sd3[agd, genomewide_fdr_selected_label := i.genomewide_fdr_selected_label]
sd3[agd, global_shift_selected_label := i.global_shift_selected_label]
sd3[agd, broad_agreement_selected_label := i.broad_agreement_selected_label]
sd3[agd, interpretive_selection_class := i.interpretive_selection_class]
sd3[agd, sweep_q_value := i.sweep_q_value_chr_aware]
sd3[agd, balancing_q_value := i.balancing_q_value_chr_aware]
sd3[agd, sweep_q_value_genomewide := i.sweep_q_value_genomewide]
sd3[agd, balancing_q_value_genomewide := i.balancing_q_value_genomewide]
if ("rf_pmax" %in% names(agd)) sd3[agd, rf_pmax := i.rf_pmax]

assign_path <- file.path(selection_output_dir, "revised_selection_tier_window_assignments.tsv")
if (file.exists(assign_path)) {
  assign_dt <- as.data.table(fread(assign_path))
  setnames(assign_dt, clean_names(names(assign_dt)))
  setkeyv(assign_dt, key_cols)
  assign_dt[agd, empirical_fdr_selected_label := i.empirical_fdr_selected_label]
  assign_dt[agd, genomewide_fdr_selected_label := i.genomewide_fdr_selected_label]
  assign_dt[agd, global_shift_selected_label := i.global_shift_selected_label]
  assign_dt[agd, broad_agreement_selected_label := i.broad_agreement_selected_label]
  assign_dt[agd, interpretive_selection_class := i.interpretive_selection_class]
  assign_dt[agd, sweep_q_value := i.sweep_q_value_chr_aware]
  assign_dt[agd, balancing_q_value := i.balancing_q_value_chr_aware]
  assign_dt[agd, sweep_q_value_genomewide := i.sweep_q_value_genomewide]
  assign_dt[agd, balancing_q_value_genomewide := i.balancing_q_value_genomewide]
  fwrite(assign_dt, assign_path, sep = "\t")
}

summarise_label_tier <- function(df, tier_name, label_vec) {
  z <- copy(df)
  z[, tier_label := as.character(label_vec)]
  z[, tier := tier_name]
  z[is.na(tier_label) | !nzchar(tier_label), tier_label := "No consensus"]
  z[, window_bp := fifelse(is.finite(.end_bp) & is.finite(.start_bp), pmax(.end_bp - .start_bp + 1, 0), NA_real_)]
  abs_dp_vec <- if ("mean_abs_delta_p" %in% names(z)) z$mean_abs_delta_p else if ("abs_delta_p" %in% names(z)) z$abs_delta_p else rep(NA_real_, nrow(z))
  b_vec <- if ("B_no_delta_p" %in% names(z)) z$B_no_delta_p else rep(NA_real_, nrow(z))
  k_vec <- if ("K_local_mean" %in% names(z)) z$K_local_mean else rep(NA_real_, nrow(z))
  z[, .(
    n_windows = .N,
    fraction_of_genome = .N / nrow(z),
    bp_covered = sum(window_bp, na.rm = TRUE),
    mean_avg_wc_fst = safe_mean(avg_wc_fst),
    mean_avg_dxy = safe_mean(avg_dxy),
    mean_pi_average = safe_mean(pi_average),
    mean_TajimaD = safe_mean(TajimaD),
    mean_f3 = safe_mean(f3),
    mean_recombination = safe_mean(recombination),
    mean_abs_delta_p = safe_mean(abs_dp_vec),
    mean_B_no_delta_p = safe_mean(b_vec),
    mean_K_local_mean = safe_mean(k_vec)
  ), by = .(tier, class = tier_label)]
}

emp_fdr_tier <- summarise_label_tier(agd, "Empirical FDR <= 0.05", agd$empirical_fdr_selected_label)
genomewide_fdr_tier <- summarise_label_tier(agd, "Genome-wide FDR <= 0.05", agd$genomewide_fdr_selected_label)
global_shift_tier <- summarise_label_tier(agd, "Global shift selected windows", agd$global_shift_selected_label)
if (nrow(tier_counts)) {
  tier_counts <- rbind(tier_counts, emp_fdr_tier, genomewide_fdr_tier, global_shift_tier, use.names = TRUE, fill = TRUE)
  fwrite(tier_counts, file.path(selection_output_dir, "revised_selection_tier_counts.tsv"), sep = "\t")
}
fwrite(tier_counts, file.path(qc_dir, "QC_tier_summary_by_label.tsv"), sep = "\t")

decision_path <- file.path(selection_output_dir, "revised_model_decision_table.tsv")
if (file.exists(decision_path)) {
  dec <- as.data.table(fread(decision_path))
  if (!("Empirical FDR <= 0.05" %chin% dec$tier)) {
    dec <- rbind(
      dec,
      data.table(
        tier = "Empirical FDR <= 0.05",
        neutral_windows = 0L,
        balancing_windows = emp_fdr_tier[class == "Balancing selection", n_windows][1],
        sweep_windows = emp_fdr_tier[class == "Geographic sweep", n_windows][1],
        ambiguous_no_consensus_windows = nrow(agd) - sum(emp_fdr_tier$n_windows[emp_fdr_tier$class %chin% c("Balancing selection", "Geographic sweep")], na.rm = TRUE),
        balancing_genes = NA_real_,
        sweep_genes = NA_real_,
        BGS_like = NA_real_,
        sweep_like_excess = NA_real_,
        abs_delta_p_R2 = NA_real_,
        f3_R2 = NA_real_,
        number_of_significant_GO_terms = NA_real_,
        immune_Fisher_OR = NA_real_,
        immune_Fisher_p_value = NA_real_,
        immune_Fisher_OR_balancing = NA_real_,
        immune_Fisher_p_value_balancing = NA_real_,
        recommended_paper_use = "primary significant sweep/balancing calls"
      ),
      data.table(
        tier = "Genome-wide FDR <= 0.05",
        neutral_windows = 0L,
        balancing_windows = genomewide_fdr_tier[class == "Balancing selection", n_windows][1],
        sweep_windows = genomewide_fdr_tier[class == "Geographic sweep", n_windows][1],
        ambiguous_no_consensus_windows = nrow(agd) - sum(genomewide_fdr_tier$n_windows[genomewide_fdr_tier$class %chin% c("Balancing selection", "Geographic sweep")], na.rm = TRUE),
        balancing_genes = NA_real_, sweep_genes = NA_real_, BGS_like = NA_real_, sweep_like_excess = NA_real_,
        abs_delta_p_R2 = NA_real_, f3_R2 = NA_real_, number_of_significant_GO_terms = NA_real_,
        immune_Fisher_OR = NA_real_, immune_Fisher_p_value = NA_real_,
        immune_Fisher_OR_balancing = NA_real_, immune_Fisher_p_value_balancing = NA_real_,
        recommended_paper_use = "global chromosome-shift sensitivity"
      ),
      data.table(
        tier = "Global shift selected windows",
        neutral_windows = 0L,
        balancing_windows = global_shift_tier[class == "Balancing selection", n_windows][1],
        sweep_windows = global_shift_tier[class == "Geographic sweep", n_windows][1],
        ambiguous_no_consensus_windows = nrow(agd) - sum(global_shift_tier$n_windows[global_shift_tier$class %chin% c("Balancing selection", "Geographic sweep")], na.rm = TRUE),
        balancing_genes = NA_real_, sweep_genes = NA_real_, BGS_like = NA_real_, sweep_like_excess = NA_real_,
        abs_delta_p_R2 = NA_real_, f3_R2 = NA_real_, number_of_significant_GO_terms = NA_real_,
        immune_Fisher_OR = NA_real_, immune_Fisher_p_value = NA_real_,
        immune_Fisher_OR_balancing = NA_real_, immune_Fisher_p_value_balancing = NA_real_,
        recommended_paper_use = "genomewide-significant but not chromosome-local outliers"
      ),
      use.names = TRUE,
      fill = TRUE
    )
    fwrite(dec, decision_path, sep = "\t")
  }
}

global_chr_summary <- agd[, .(
  n_windows = .N,
  mean_z_dxy = safe_mean(z_dxy),
  mean_z_fst = safe_mean(z_fst),
  mean_z_pi = safe_mean(z_pi),
  mean_z_f3 = safe_mean(z_f3),
  mean_z_taj = safe_mean(z_taj),
  mean_z_rec = safe_mean(z_rec),
  mean_sweep_strength = safe_mean(sweep_strength),
  mean_balancing_strength = safe_mean(balancing_strength),
  frac_gmmrf_neutral = mean(gmm_rf_label == "Neutral / equilibrium", na.rm = TRUE),
  frac_gmmrf_balancing = mean(gmm_rf_label == "Balancing selection", na.rm = TRUE),
  frac_gmmrf_sweep = mean(gmm_rf_label == "Geographic sweep", na.rm = TRUE),
  frac_gmmrf_contdist_neutral = mean(gmm_rf_contdist_label == "Neutral / equilibrium", na.rm = TRUE),
  frac_gmmrf_contdist_balancing = mean(gmm_rf_contdist_label == "Balancing selection", na.rm = TRUE),
  frac_gmmrf_contdist_sweep = mean(gmm_rf_contdist_label == "Geographic sweep", na.rm = TRUE),
  frac_empirical_fdr_balancing = mean(empirical_fdr_selected_label == "Balancing selection", na.rm = TRUE),
  frac_empirical_fdr_sweep = mean(empirical_fdr_selected_label == "Geographic sweep", na.rm = TRUE),
  frac_genomewide_fdr_balancing = mean(genomewide_fdr_selected_label == "Balancing selection", na.rm = TRUE),
  frac_genomewide_fdr_sweep = mean(genomewide_fdr_selected_label == "Geographic sweep", na.rm = TRUE)
), by = Chromosome_Names][order(factor(Chromosome_Names, levels = ordered_chr_levels(Chromosome_Names)))]
fwrite(global_chr_summary, file.path(qc_dir, "QC_global_chromosome_summary.tsv"), sep = "\t")

auto_mask <- agd$Chromosome_Names != "Chromosome_ZW"
zw_vs_auto_rows <- rbindlist(lapply(list(
  list(metric = "z_fst", value = agd$z_fst),
  list(metric = "z_dxy", value = agd$z_dxy),
  list(metric = "z_pi", value = agd$z_pi),
  list(metric = "z_taj", value = agd$z_taj),
  list(metric = "z_f3", value = agd$z_f3),
  list(metric = "sweep_strength", value = agd$sweep_strength),
  list(metric = "balancing_strength", value = agd$balancing_strength),
  list(metric = "frac_broad_agreement_selected_like", value = as.numeric(agd$broad_agreement_selected_label %chin% c("Geographic sweep", "Balancing selection"))),
  list(metric = "frac_local_outlier_selected_like", value = as.numeric(agd$empirical_fdr_selected_label %chin% c("Geographic sweep", "Balancing selection"))),
  list(metric = "frac_genomewide_shift_selected_like", value = as.numeric(agd$genomewide_fdr_selected_label %chin% c("Geographic sweep", "Balancing selection")))
), function(obj) {
  chr_col <- paste0("mean_", obj$metric)
  chr_tmp <- agd[, .(metric_value = safe_mean(obj$value)), by = Chromosome_Names]
  setnames(chr_tmp, "metric_value", chr_col)
  data.table(
    metric = obj$metric,
    zw_mean = safe_mean(obj$value[agd$Chromosome_Names == "Chromosome_ZW"]),
    autosome_mean = safe_mean(obj$value[auto_mask]),
    wilcoxon_p = safe_wilcox_p(obj$value[agd$Chromosome_Names == "Chromosome_ZW"], obj$value[auto_mask]),
    chromosome_permutation_p = permute_chr_summary_p(chr_tmp, chr_col)
  )
}), use.names = TRUE, fill = TRUE)
fwrite(zw_vs_auto_rows, file.path(qc_dir, "QC_ZW_vs_autosomes.tsv"), sep = "\t")

local_chr_summary <- agd[, .(
  sweep_candidates_preFDR = sum(primary_sweep_candidate %in% TRUE, na.rm = TRUE),
  balancing_candidates_preFDR = sum(primary_balancing_candidate %in% TRUE, na.rm = TRUE),
  best_sweep_q_chr_aware = suppressWarnings(min(sweep_q_value_chr_aware[is.finite(sweep_q_value_chr_aware)], na.rm = TRUE)),
  best_balancing_q_chr_aware = suppressWarnings(min(balancing_q_value_chr_aware[is.finite(balancing_q_value_chr_aware)], na.rm = TRUE)),
  best_sweep_q_genomewide = suppressWarnings(min(sweep_q_value_genomewide[is.finite(sweep_q_value_genomewide)], na.rm = TRUE)),
  best_balancing_q_genomewide = suppressWarnings(min(balancing_q_value_genomewide[is.finite(balancing_q_value_genomewide)], na.rm = TRUE)),
  sweep_fdr_chr_aware = sum(empirical_fdr_selected_label == "Geographic sweep", na.rm = TRUE),
  balancing_fdr_chr_aware = sum(empirical_fdr_selected_label == "Balancing selection", na.rm = TRUE),
  sweep_fdr_genomewide = sum(genomewide_fdr_selected_label == "Geographic sweep", na.rm = TRUE),
  balancing_fdr_genomewide = sum(genomewide_fdr_selected_label == "Balancing selection", na.rm = TRUE)
), by = Chromosome_Names][order(factor(Chromosome_Names, levels = ordered_chr_levels(Chromosome_Names)))]
for (cc in c("best_sweep_q_chr_aware", "best_balancing_q_chr_aware", "best_sweep_q_genomewide", "best_balancing_q_genomewide")) {
  local_chr_summary[!is.finite(get(cc)), (cc) := NA_real_]
}
fwrite(local_chr_summary, file.path(qc_dir, "QC_local_outlier_by_chromosome.tsv"), sep = "\t")

if ("rf_pmax" %in% names(agd)) {
  agd[, rf_lowconf_lt_0_6 := rf_pmax < 0.6]
  agd[, rf_lowconf_lt_0_7 := rf_pmax < 0.7]
  agd[, rf_lowconf_lt_0_8 := rf_pmax < 0.8]
}

spatial_specs <- list(
  "GMM/RF + continuous-distance agreement" = agd$gmm_rf_contdist_label,
  "Distance margin >= 0.50" = agd$conservative_distance_margin_050_label,
  "Empirical FDR <= 0.05" = agd$empirical_fdr_selected_label,
  "Top 1%" = agd$top1pct_distance_margin_050_label
)
spatial_rows <- rbindlist(lapply(names(spatial_specs), function(tt) {
  lbl <- as.character(spatial_specs[[tt]])
  keep <- lbl %chin% c("Geographic sweep", "Balancing selection")
  if (!any(keep, na.rm = TRUE)) return(NULL)
  x <- agd[keep][order(Chromosome_Names, .start_bp)]
  gap_cutoff <- x[, median(diff(.start_bp), na.rm = TRUE)]
  if (!is.finite(gap_cutoff)) gap_cutoff <- 0
  x[, new_run := is.na(shift(.start_bp)) | Chromosome_Names != shift(Chromosome_Names) | (.start_bp - shift(.start_bp)) > (gap_cutoff * 1.5)]
  x[, run_id := cumsum(fifelse(is.na(new_run), TRUE, new_run))]
  rr <- x[, .(run_len = .N, bp = sum(.end_bp - .start_bp + 1, na.rm = TRUE)), by = .(Chromosome_Names, run_id)]
  chr_sum <- x[, .N, by = Chromosome_Names]
  data.table(
    tier = tt,
    n_selected_windows = nrow(x),
    n_chromosomes = uniqueN(x$Chromosome_Names),
    longest_run_windows = max(rr$run_len, na.rm = TRUE),
    mean_run_length = mean(rr$run_len, na.rm = TRUE),
    n_selected_domains = nrow(rr),
    top_chromosome = chr_sum[which.max(N), Chromosome_Names]
  )
}), use.names = TRUE, fill = TRUE)
fwrite(spatial_rows, file.path(qc_dir, "QC_chromosome_spatial_summary.tsv"), sep = "\t")

run_null_rows <- rbindlist(lapply(names(spatial_specs), function(tt) {
  lbl <- as.character(spatial_specs[[tt]])
  keep <- lbl %chin% c("Geographic sweep", "Balancing selection")
  if (!any(keep, na.rm = TRUE)) return(NULL)
  obs_mask <- as.integer(keep)
  obs_run <- rle(obs_mask)
  obs_longest <- max(obs_run$lengths[obs_run$values == 1], na.rm = TRUE)
  null_long <- replicate(200, {
    perm <- integer(length(obs_mask))
    for (cc in unique(agd$Chromosome_Names)) {
      idx <- which(agd$Chromosome_Names == cc)
      perm[idx] <- as.integer(block_permute(obs_mask[idx], block_n = 25L))
    }
    rr <- rle(perm)
    max(rr$lengths[rr$values == 1], na.rm = TRUE)
  })
  data.table(
    tier = tt,
    observed_longest_run = obs_longest,
    null_mean = mean(null_long, na.rm = TRUE),
    null_q95 = safe_quant(null_long, 0.95),
    empirical_p = (1 + sum(null_long >= obs_longest, na.rm = TRUE)) / (1 + length(null_long))
  )
}), use.names = TRUE, fill = TRUE)
fwrite(run_null_rows, file.path(qc_dir, "QC_run_length_nulls.tsv"), sep = "\t")
if (nrow(run_null_rows)) {
  p_run <- ggplot(run_null_rows, aes(tier, observed_longest_run)) +
    geom_col(fill = "#3182bd") +
    geom_point(aes(y = null_q95), color = "#cb181d", size = 2.2) +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(x = NULL, y = "Longest run (windows)", title = "Observed longest selected-window runs vs null 95th percentile")
  ggsave(file.path(qc_dir, "selected_window_run_length_nulls.pdf"), p_run, width = 8.2, height = 4.8)
}

if (nrow(tier_counts)) {
  prof_plot_dt <- tier_counts[class %chin% c("Geographic sweep", "Balancing selection") & tier %chin% c("GMM/RF", "GMM/RF + continuous-distance agreement", "Distance margin >= 0.50", "Empirical FDR <= 0.05", "Top 1%", "Top 500")]
  prof_long <- melt(
    prof_plot_dt,
    id.vars = c("tier", "class"),
    measure.vars = intersect(c("mean_avg_wc_fst", "mean_avg_dxy", "mean_pi_average", "mean_TajimaD", "mean_f3", "mean_recombination", "mean_abs_delta_p"), names(prof_plot_dt)),
    variable.name = "metric",
    value.name = "value"
  )
  if (nrow(prof_long)) {
    p_prof <- ggplot(prof_long, aes(tier, value, color = class, group = class)) +
      geom_point() +
      geom_line() +
      facet_wrap(~ metric, scales = "free_y", ncol = 3) +
      coord_flip() +
      theme_bw(base_size = 10) +
      labs(x = NULL, y = "Mean value", color = "Class", title = "Tier population-genetic profiles")
    ggsave(file.path(qc_dir, "tier_population_genetic_profiles.pdf"), p_prof, width = 11.5, height = 8.5)
  }
}

qc_bgs_path <- file.path(selection_output_dir, "revised_BGS_by_tier.tsv")
if (file.exists(qc_bgs_path)) {
  qc_bgs <- as.data.table(fread(qc_bgs_path))
} else {
  qc_bgs <- data.table(note = "downstream BGS summary not yet available")
}
fwrite(qc_bgs, file.path(qc_dir, "QC_BGS_by_tier.tsv"), sep = "\t")
if (nrow(qc_bgs) && "short_tier" %in% names(qc_bgs) && "fraction_sweep_windows_beyond_BGS" %in% names(qc_bgs)) {
  p_bgs <- ggplot(qc_bgs[!is.na(short_tier)], aes(short_tier, fraction_sweep_windows_beyond_BGS)) +
    geom_col(fill = "#dd8452") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(x = NULL, y = "Fraction beyond BGS", title = "Sweep-tier BGS profiles")
  ggsave(file.path(qc_dir, "BGS_sweep_tier_profiles.pdf"), p_bgs, width = 8.2, height = 5.2)
}

track_dt <- rbindlist(list(
  agd[gmm_rf_contdist_label %chin% c("Geographic sweep", "Balancing selection"), .(tier = "GMM/RF + contdist", class = gmm_rf_contdist_label, Chromosome_Names, midpoint = (.start_bp + .end_bp) / 2)],
  agd[empirical_fdr_selected_label %chin% c("Geographic sweep", "Balancing selection"), .(tier = "Empirical FDR <= 0.05", class = empirical_fdr_selected_label, Chromosome_Names, midpoint = (.start_bp + .end_bp) / 2)]
), use.names = TRUE, fill = TRUE)
if (nrow(track_dt)) {
  track_dt[, Chromosome_Names := factor(Chromosome_Names, levels = rev(ordered_chr_levels(Chromosome_Names)))]
  p_tracks <- ggplot(track_dt, aes(midpoint / 1e6, Chromosome_Names, color = class)) +
    geom_point(alpha = 0.5, size = 0.3) +
    facet_wrap(~ tier, ncol = 1) +
    theme_bw(base_size = 10) +
    labs(x = "Position (Mb)", y = "Chromosome", color = "Class", title = "Selected-window chromosome tracks")
  ggsave(file.path(qc_dir, "chromosome_selected_window_tracks.pdf"), p_tracks, width = 11.5, height = 8.5)
}

agree_dat <- agd
sd3_new_af2 <- sd3
save(agree_dat, sd3_new_af2, file = label_rdata)

if (!is.null(sanity) && exists("finalize_sanity")) {
  finalize_sanity(
    sanity,
    files = c(
      "QC_input_summary.tsv",
      "QC_correlations.tsv",
      "QC_ranknorm_summary.tsv",
      "QC_GMM_model_selection.tsv",
      "QC_GMM_centroids.tsv",
      "QC_RF_confusion_matrix.tsv",
      "QC_RF_probability_summary.tsv",
      "QC_continuous_centroids.tsv",
      "QC_continuous_confusion_matrices.tsv",
      "QC_tier_summary_by_label.tsv",
      "QC_empirical_FDR_selected_windows.tsv",
      "QC_BGS_by_tier.tsv",
      "QC_chromosome_spatial_summary.tsv"
    ),
    filename = "selection_qc_sanity_report.csv"
  )
}

message("Selection QC/FDR helper complete: ", qc_dir)
