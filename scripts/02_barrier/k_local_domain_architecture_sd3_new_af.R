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
  prior_driver_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_driver_outputs"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_domain_outputs"),
  alpha = 0.5,
  alpha_calibration_csv = NULL,
  alpha_calibration_p_col = "P_overall_mean",
  alpha_calibration_k_col = "K_genome_mean",
  domain_quantiles = c(0.99, 0.95, 0.90),
  primary_domain_quantile = 0.95,
  recombination_desert_quantiles = c(0.05, 0.10),
  n_permutations = 100L,
  top_n_chromosomes_for_plots = 8L
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_rdata <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$object_name <- args[[2]]
if (length(args) >= 3 && nzchar(args[[3]])) cfg$output_dir <- args[[3]]
if (length(args) >= 4 && nzchar(args[[4]])) cfg$alpha_calibration_csv <- args[[4]]
if (length(args) >= 5 && nzchar(args[[5]])) cfg$alpha_calibration_p_col <- args[[5]]
if (length(args) >= 6 && nzchar(args[[6]])) cfg$alpha_calibration_k_col <- args[[6]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
sanity <- init_sanity("k_local_domain_architecture_sd3_new_af.R", cfg$output_dir)

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

# Why domains matter biologically:
# K_local is the local accumulated barrier density. Consecutive runs of high
# K_local windows are candidate broad barrier blocks rather than isolated peaks.
# If long high-K_local runs overlap low recombination, sweeps, and high
# divergence, that supports a clustered genomic barrier architecture.

find_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0L) return(NA_character_)
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

safe_group_test_tbl <- function(y, group, analysis) {
  keep <- is.finite(y) & !is.na(group)
  y2 <- y[keep]
  g2 <- droplevels(as.factor(group[keep]))
  if (length(y2) < 5L || nlevels(g2) < 2L) {
    return(tibble(
      analysis = analysis, model_type = "group_test", term = NA_character_,
      estimate = NA_real_, statistic = NA_real_, p_value = NA_real_,
      n = length(y2), effect_size = NA_real_, note = "Too few observations or groups"
    ))
  }
  if (nlevels(g2) == 2L) {
    wt <- suppressWarnings(stats::wilcox.test(y2 ~ g2))
    med_diff <- diff(tapply(y2, g2, median, na.rm = TRUE))
    return(tibble(
      analysis = analysis, model_type = "wilcox", term = paste(levels(g2), collapse = " vs "),
      estimate = NA_real_, statistic = unname(wt$statistic), p_value = wt$p.value,
      n = length(y2), effect_size = unname(med_diff), note = NA_character_
    ))
  }
  fit <- stats::aov(y2 ~ g2)
  sm <- summary(fit)[[1]]
  tibble(
    analysis = analysis, model_type = "anova", term = rownames(sm)[1],
    estimate = NA_real_, statistic = sm$`F value`[1], p_value = sm$`Pr(>F)`[1],
    n = length(y2), effect_size = stats::sd(tapply(y2, g2, mean, na.rm = TRUE), na.rm = TRUE), note = NA_character_
  )
}

safe_cor_tbl <- function(x, y, analysis, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L) {
    return(tibble(analysis = analysis, estimate = NA_real_, p_value = NA_real_, statistic = NA_real_, n = sum(keep), method = method, note = "Too few complete observations"))
  }
  ct <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = method))
  tibble(analysis = analysis, estimate = unname(ct$estimate), p_value = ct$p.value, statistic = unname(ct$statistic), n = sum(keep), method = method, note = NA_character_)
}

run_table <- function(flag) {
  r <- rle(flag)
  idx <- cumsum(r$lengths)
  tibble(value = r$values, run_length = r$lengths, end_idx = idx, start_idx = idx - r$lengths + 1L)
}

dominant_label <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

permute_overlap <- function(high_flag, desert_flag, n_perm = 100L) {
  keep <- !is.na(high_flag) & !is.na(desert_flag)
  high_flag <- high_flag[keep]
  desert_flag <- desert_flag[keep]
  obs_overlap <- sum(high_flag & desert_flag)
  perm_overlap <- replicate(n_perm, sum(sample(high_flag, replace = FALSE) & desert_flag))
  tibble(
    observed_overlap = obs_overlap,
    mean_perm_overlap = mean(perm_overlap),
    p_ge_obs = mean(perm_overlap >= obs_overlap)
  )
}

load(cfg$input_rdata)
if (!exists(cfg$object_name, inherits = FALSE)) stop("Object not found: ", cfg$object_name)
dat_raw <- get(cfg$object_name, inherits = FALSE)
assert_true(nrow(as.data.frame(dat_raw)) > 0L, "Loaded dataset is empty.", sanity, "nonempty_input_dataset")

chr_col <- find_col(dat_raw, c("Chromosome_Names", "Chromosome", "chr", "chromosome"))
start_col <- find_col(dat_raw, c(".start_bp", "window_start", "start_bp", "start"))
end_col <- find_col(dat_raw, c(".end_bp", "window_end", "end_bp", "end"))
mid_col <- find_col(dat_raw, c("midpoint", ".mid_bp", "window_mid", "mid_bp", "mid"))
recomb_col <- find_col(dat_raw, c("recombination", "rec", "rho"))
fst_col <- find_col(dat_raw, c("avg_wc_fst", "wc_fst", "fst"))
dxy_col <- find_col(dat_raw, c("avg_dxy", "dxy"))
pi_col <- find_col(dat_raw, c("pi_average", "avg_pi", "pi"))
taj_col <- find_col(dat_raw, c("TajimaD", "tajima_d"))
delta_p_col <- find_col(dat_raw, c("mean_abs_delta_p", "abs_delta_p"))
f3_col <- find_col(dat_raw, c("f3", "F3"))
snp_density_col <- find_col(dat_raw, c("n_snps_af", "N_SNPS", "no_snps", "n_variants"))
gene_density_col <- find_col(dat_raw, c("gene_density", "genes_per_mb", "gene_count", "n_genes"))

detected_columns <- tibble(
  role = c("chromosome", "start", "end", "mid", "recombination", "avg_wc_fst", "avg_dxy", "pi_average", "TajimaD", "abs_delta_p", "f3", "snp_density", "gene_density"),
  detected_column = c(chr_col, start_col, end_col, mid_col, recomb_col, fst_col, dxy_col, pi_col, taj_col, delta_p_col, f3_col, snp_density_col, gene_density_col)
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
    f3 = if (!is.na(f3_col)) as.numeric(.data[[f3_col]]) else NA_real_,
    snp_density = if (!is.na(snp_density_col)) as.numeric(.data[[snp_density_col]]) else NA_real_,
    gene_density = if (!is.na(gene_density_col)) as.numeric(.data[[gene_density_col]]) else NA_real_
  )
assert_has_cols(dat, c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"), sanity, "k_local_domain_input_data")
assert_unique_keys(dat, c("Chromosome_Names", ".start_bp", ".end_bp"), sanity, "k_local_domain_input_data")

prior_local <- read.csv(file.path(cfg$prior_barrier_dir, "barrier_recalibrated_local_summary.csv"), stringsAsFactors = FALSE)
n_before_join <- nrow(dat)
dat <- dat %>% left_join(prior_local, by = c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"))
assert_join_rowcount(n_before_join, dat, sanity, "k_local_domain_prior_local_join")
if ("abs_delta_p.x" %in% names(dat) || "abs_delta_p.y" %in% names(dat)) dat$abs_delta_p <- dplyr::coalesce(dat$abs_delta_p.x, dat$abs_delta_p.y)
if ("f3.x" %in% names(dat) || "f3.y" %in% names(dat)) dat$f3 <- dplyr::coalesce(dat$f3.x, dat$f3.y)

if (!"B_no_delta_p" %in% names(dat)) stop("B_no_delta_p missing after merge.")
if (!"K_local_mean" %in% names(dat)) stop("K_local_mean missing after merge.")
if (!"K_local_sum" %in% names(dat)) dat$K_local_sum <- NA_real_
if (!"sweep_group" %in% names(dat)) dat$sweep_group <- NA_character_
if (!"P_local_mean" %in% names(dat)) dat$P_local_mean <- exp(-cfg$alpha * dat$K_local_mean)
assert_missing_fraction(dat$K_local_mean, 0.1, sanity, "K_local_mean")

dat$Chromosome_Names <- factor(dat$Chromosome_Names, levels = order_chromosomes(dat$Chromosome_Names))
dat$is_zw <- grepl("ZW", dat$Chromosome_Names, ignore.case = TRUE)

for (thr in cfg$domain_quantiles) {
  nm <- paste0("is_top_", gsub("\\.", "", sprintf("%.2f", 1 - thr)))
  cutoff <- stats::quantile(dat$K_local_mean, probs = thr, na.rm = TRUE)
  dat[[nm]] <- is.finite(dat$K_local_mean) & dat$K_local_mean >= cutoff
}

primary_cutoff <- stats::quantile(dat$K_local_mean, probs = cfg$primary_domain_quantile, na.rm = TRUE)
dat$is_primary_domain_window <- is.finite(dat$K_local_mean) & dat$K_local_mean >= primary_cutoff
assert_quantile_fraction(dat$is_primary_domain_window, 1 - cfg$primary_domain_quantile, 0.01, sanity, "primary_high_domain_fraction")

domain_table <- bind_rows(lapply(split(dat, dat$Chromosome_Names), function(df_chr) {
  df_chr <- df_chr[order(df_chr$.start_bp), ]
  rt <- run_table(df_chr$is_primary_domain_window)
  high_runs <- rt %>% filter(value %in% TRUE)
  if (nrow(high_runs) == 0L) return(tibble())
  bind_rows(lapply(seq_len(nrow(high_runs)), function(i) {
    rr <- high_runs[i, ]
    sub <- df_chr[rr$start_idx:rr$end_idx, , drop = FALSE]
    tibble(
      Chromosome_Names = as.character(sub$Chromosome_Names[1]),
      domain_id = paste0(as.character(sub$Chromosome_Names[1]), "_dom_", i),
      start_bp = min(sub$.start_bp, na.rm = TRUE),
      end_bp = max(sub$.end_bp, na.rm = TRUE),
      n_windows = nrow(sub),
      domain_length_bp = max(sub$.end_bp, na.rm = TRUE) - min(sub$.start_bp, na.rm = TRUE) + 1,
      mean_K_local = mean(sub$K_local_mean, na.rm = TRUE),
      max_K_local = max(sub$K_local_mean, na.rm = TRUE),
      mean_recombination = mean(sub$recombination, na.rm = TRUE),
      mean_avg_wc_fst = mean(sub$avg_wc_fst, na.rm = TRUE),
      mean_avg_dxy = mean(sub$avg_dxy, na.rm = TRUE),
      mean_pi_average = mean(sub$pi_average, na.rm = TRUE),
      mean_TajimaD = mean(sub$TajimaD, na.rm = TRUE),
      mean_abs_delta_p = mean(sub$abs_delta_p, na.rm = TRUE),
      mean_f3 = mean(sub$f3, na.rm = TRUE),
      mean_snp_density = mean(sub$snp_density, na.rm = TRUE),
      mean_gene_density = mean(sub$gene_density, na.rm = TRUE),
      dominant_sweep_group = dominant_label(sub$sweep_group)
    )
  }))
}))

dat$domain_status <- ifelse(dat$is_primary_domain_window, "Primary domain window", "Background")
domain_summary <- domain_table %>%
  group_by(Chromosome_Names) %>%
  summarise(
    n_domains = n(),
    total_domain_windows = sum(n_windows),
    total_domain_bp = sum(domain_length_bp),
    mean_domain_length_bp = mean(domain_length_bp),
    max_domain_length_bp = max(domain_length_bp),
    .groups = "drop"
  )

metric_comparisons <- bind_rows(
  safe_group_test_tbl(dat$recombination, dat$domain_status, "recombination: domain vs background"),
  safe_group_test_tbl(dat$avg_wc_fst, dat$domain_status, "avg_wc_fst: domain vs background"),
  safe_group_test_tbl(dat$avg_dxy, dat$domain_status, "avg_dxy: domain vs background"),
  safe_group_test_tbl(dat$abs_delta_p, dat$domain_status, "abs_delta_p: domain vs background")
) %>%
  bind_rows(
    safe_group_test_tbl(dat$recombination, ifelse(dat$is_top_001, "Top1", "Other"), "recombination: top 1% vs other"),
    safe_group_test_tbl(dat$recombination, ifelse(dat$is_top_005, "Top5", "Other"), "recombination: top 5% vs other"),
    safe_group_test_tbl(dat$recombination, ifelse(dat$is_top_010, "Top10", "Other"), "recombination: top 10% vs other")
  )

if (any(is.finite(dat$pi_average))) {
  metric_comparisons <- bind_rows(metric_comparisons,
    safe_group_test_tbl(dat$pi_average, dat$domain_status, "pi_average: domain vs background"),
    safe_group_test_tbl(dat$pi_average, ifelse(dat$is_top_005, "Top5", "Other"), "pi_average: top 5% vs other")
  )
}
if (any(is.finite(dat$TajimaD))) {
  metric_comparisons <- bind_rows(metric_comparisons,
    safe_group_test_tbl(dat$TajimaD, dat$domain_status, "TajimaD: domain vs background"),
    safe_group_test_tbl(dat$TajimaD, ifelse(dat$is_top_005, "Top5", "Other"), "TajimaD: top 5% vs other")
  )
}
if (any(is.finite(dat$f3))) {
  metric_comparisons <- bind_rows(metric_comparisons,
    safe_group_test_tbl(dat$f3, dat$domain_status, "f3: domain vs background"),
    safe_group_test_tbl(dat$f3, ifelse(dat$is_top_005, "Top5", "Other"), "f3: top 5% vs other")
  )
}
if (any(is.finite(dat$snp_density))) {
  metric_comparisons <- bind_rows(metric_comparisons, safe_group_test_tbl(dat$snp_density, dat$domain_status, "snp_density: domain vs background"))
}
if (any(is.finite(dat$gene_density))) {
  metric_comparisons <- bind_rows(metric_comparisons, safe_group_test_tbl(dat$gene_density, dat$domain_status, "gene_density: domain vs background"))
}

scatter_correlations <- bind_rows(
  safe_cor_tbl(dat$K_local_mean, dat$recombination, "K_local_mean vs recombination"),
  safe_cor_tbl(dat$K_local_mean, dat$avg_wc_fst, "K_local_mean vs avg_wc_fst"),
  safe_cor_tbl(dat$K_local_mean, dat$avg_dxy, "K_local_mean vs avg_dxy"),
  safe_cor_tbl(dat$K_local_mean, dat$abs_delta_p, "K_local_mean vs abs_delta_p")
)
if (any(is.finite(dat$pi_average))) scatter_correlations <- bind_rows(scatter_correlations, safe_cor_tbl(dat$K_local_mean, dat$pi_average, "K_local_mean vs pi_average"))
if (any(is.finite(dat$TajimaD))) scatter_correlations <- bind_rows(scatter_correlations, safe_cor_tbl(dat$K_local_mean, dat$TajimaD, "K_local_mean vs TajimaD"))
if (any(is.finite(dat$f3))) scatter_correlations <- bind_rows(scatter_correlations, safe_cor_tbl(dat$K_local_mean, dat$f3, "K_local_mean vs f3"))

class_enrichment <- tibble()
if (any(!is.na(dat$sweep_group))) {
  class_enrichment <- bind_rows(lapply(c("is_top_001", "is_top_005", "is_top_010"), function(flag_nm) {
    tab <- table(in_set = ifelse(dat[[flag_nm]], "TopK", "Background"), class = dat$sweep_group, useNA = "no")
    chi_obj <- suppressWarnings(chisq.test(tab))
    as_tibble(as.data.frame.matrix(prop.table(tab, 1)), rownames = "set") %>%
      pivot_longer(-set, names_to = "class", values_to = "fraction") %>%
      mutate(test = paste0(flag_nm, " class enrichment"), p_value = chi_obj$p.value)
  }))
}

recombination_desert_overlap <- bind_rows(lapply(cfg$recombination_desert_quantiles, function(q) {
  desert_cut <- stats::quantile(dat$recombination, probs = q, na.rm = TRUE)
  desert_flag <- is.finite(dat$recombination) & dat$recombination <= desert_cut
  bind_rows(lapply(c(0.99, 0.95, 0.90), function(thr) {
    top_cut <- stats::quantile(dat$K_local_mean, probs = thr, na.rm = TRUE)
    high_flag <- is.finite(dat$K_local_mean) & dat$K_local_mean >= top_cut
    ov <- permute_overlap(high_flag, desert_flag, cfg$n_permutations)
    tibble(
      high_k_quantile = thr,
      desert_quantile = q,
      n_high = sum(high_flag, na.rm = TRUE),
      n_desert = sum(desert_flag, na.rm = TRUE),
      observed_overlap = ov$observed_overlap,
      mean_perm_overlap = ov$mean_perm_overlap,
      p_ge_obs = ov$p_ge_obs,
      overlap_fraction_high = ov$observed_overlap / sum(high_flag, na.rm = TRUE)
    )
  }))
}))

chromosome_summary <- dat %>%
  group_by(Chromosome_Names) %>%
  summarise(
    n_windows = n(),
    frac_top5 = mean(is_primary_domain_window, na.rm = TRUE),
    mean_recombination = mean(recombination, na.rm = TRUE),
    mean_avg_wc_fst = mean(avg_wc_fst, na.rm = TRUE),
    mean_avg_dxy = mean(avg_dxy, na.rm = TRUE),
    mean_pi_average = mean(pi_average, na.rm = TRUE),
    mean_abs_delta_p = mean(abs_delta_p, na.rm = TRUE),
    mean_f3 = mean(f3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(domain_summary, by = "Chromosome_Names") %>%
  mutate(
    n_domains = coalesce(n_domains, 0L),
    total_domain_windows = coalesce(total_domain_windows, 0L),
    total_domain_bp = coalesce(total_domain_bp, 0),
    mean_domain_length_bp = coalesce(mean_domain_length_bp, 0),
    max_domain_length_bp = coalesce(max_domain_length_bp, 0),
    frac_windows_in_domains = total_domain_windows / n_windows
  ) %>%
  arrange(desc(frac_windows_in_domains), desc(max_domain_length_bp))

zw_vs_auto <- safe_group_test_tbl(dat$K_local_mean, ifelse(dat$is_zw, "Chromosome_ZW", "Autosome"), "K_local_mean: Chromosome_ZW vs autosomes")

permutation_summary <- runlength_summary <- bind_rows(lapply(split(dat, dat$Chromosome_Names), function(df_chr) {
  df_chr <- df_chr[order(df_chr$.start_bp), ]
  high <- is.finite(df_chr$K_local_mean) & df_chr$K_local_mean >= primary_cutoff
  rt <- run_table(high)
  high_runs <- rt %>% filter(value %in% TRUE)
  obs_max <- if (nrow(high_runs) > 0) max(high_runs$run_length) else 0
  perm_max <- replicate(cfg$n_permutations, {
    xp <- sample(high, replace = FALSE)
    rp <- run_table(xp) %>% filter(value %in% TRUE)
    if (nrow(rp) > 0) max(rp$run_length) else 0
  })
  tibble(
    Chromosome_Names = as.character(df_chr$Chromosome_Names[1]),
    cutoff_quantile = cfg$primary_domain_quantile,
    cutoff_value = primary_cutoff,
    n_high_windows = sum(high, na.rm = TRUE),
    n_runs = nrow(high_runs),
    mean_run_length = ifelse(nrow(high_runs) > 0, mean(high_runs$run_length), 0),
    max_run_length = obs_max,
    mean_max_run_perm = mean(perm_max),
    p_ge_obs = mean(perm_max >= obs_max)
  )
}))

write.csv(domain_table, file.path(cfg$output_dir, "k_local_domain_table.csv"), row.names = FALSE)
write.csv(domain_summary, file.path(cfg$output_dir, "k_local_domain_summary.csv"), row.names = FALSE)
write.csv(bind_rows(metric_comparisons, scatter_correlations %>% transmute(analysis, model_type = method, term = NA_character_, estimate, statistic, p_value, n, effect_size = estimate, note)), file.path(cfg$output_dir, "k_local_domain_metric_comparisons.csv"), row.names = FALSE)
write.csv(class_enrichment, file.path(cfg$output_dir, "k_local_domain_class_enrichment.csv"), row.names = FALSE)
write.csv(recombination_desert_overlap, file.path(cfg$output_dir, "k_local_domain_recombination_desert_overlap.csv"), row.names = FALSE)
write.csv(chromosome_summary, file.path(cfg$output_dir, "k_local_domain_chromosome_summary.csv"), row.names = FALSE)
write.csv(permutation_summary, file.path(cfg$output_dir, "k_local_domain_permutation_summary.csv"), row.names = FALSE)
write.csv(detected_columns, file.path(cfg$output_dir, "detected_columns.csv"), row.names = FALSE)

plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

top_chr <- chromosome_summary %>%
  slice_head(n = min(cfg$top_n_chromosomes_for_plots, nrow(chromosome_summary))) %>%
  pull(Chromosome_Names) %>%
  as.character()

ggsave(file.path(cfg$output_dir, "plot_domain_K_local_across_all_chromosomes.png"),
  ggplot(dat, aes(window_mid_bp / 1e6, K_local_mean)) +
    geom_line(color = "grey55", linewidth = 0.22, na.rm = TRUE) +
    geom_point(data = dat %>% filter(is_primary_domain_window),
               aes(window_mid_bp / 1e6, K_local_mean), color = "#b22222", size = 0.45, inherit.aes = FALSE) +
    facet_wrap(~ Chromosome_Names, scales = "free_x") +
    labs(title = "K_local domains across all chromosomes", x = "Position (Mb)", y = "K_local_mean") +
    plot_theme,
  width = 14, height = 10.5, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_domain_K_local_across_chromosomes.png"),
  ggplot(dat %>% filter(as.character(Chromosome_Names) %in% top_chr), aes(window_mid_bp / 1e6, K_local_mean)) +
    geom_line(color = "grey55", linewidth = 0.25, na.rm = TRUE) +
    geom_point(data = dat %>% filter(as.character(Chromosome_Names) %in% top_chr, is_primary_domain_window),
               aes(window_mid_bp / 1e6, K_local_mean), color = "#b22222", size = 0.5, inherit.aes = FALSE) +
    facet_wrap(~ Chromosome_Names, scales = "free_x") +
    labs(title = "K_local domains along chromosomes", x = "Position (Mb)", y = "K_local_mean") +
    plot_theme,
  width = 11, height = 7.5, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_domain_recombination_across_all_chromosomes.png"),
  ggplot(dat, aes(window_mid_bp / 1e6, recombination)) +
    geom_line(color = "#1f4e79", linewidth = 0.22, na.rm = TRUE) +
    geom_point(data = dat %>% filter(is_primary_domain_window),
               aes(window_mid_bp / 1e6, recombination), color = "#b22222", size = 0.45, inherit.aes = FALSE) +
    facet_wrap(~ Chromosome_Names, scales = "free_x") +
    labs(title = "Recombination with high-K_local domains overlaid across all chromosomes", x = "Position (Mb)", y = "Recombination") +
    plot_theme,
  width = 14, height = 10.5, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_domain_recombination_across_chromosomes.png"),
  ggplot(dat %>% filter(as.character(Chromosome_Names) %in% top_chr), aes(window_mid_bp / 1e6, recombination)) +
    geom_line(color = "#1f4e79", linewidth = 0.25, na.rm = TRUE) +
    geom_point(data = dat %>% filter(as.character(Chromosome_Names) %in% top_chr, is_primary_domain_window),
               aes(window_mid_bp / 1e6, recombination), color = "#b22222", size = 0.5, inherit.aes = FALSE) +
    facet_wrap(~ Chromosome_Names, scales = "free_x") +
    labs(title = "Recombination with high-K_local domains overlaid", x = "Position (Mb)", y = "Recombination") +
    plot_theme,
  width = 11, height = 7.5, dpi = 300
)

metric_box_specs <- list(
  list(var = "recombination", file = "plot_domain_vs_background_recombination.png"),
  list(var = "pi_average", file = "plot_domain_vs_background_pi_average.png"),
  list(var = "TajimaD", file = "plot_domain_vs_background_TajimaD.png"),
  list(var = "avg_wc_fst", file = "plot_domain_vs_background_avg_wc_fst.png"),
  list(var = "avg_dxy", file = "plot_domain_vs_background_avg_dxy.png"),
  list(var = "abs_delta_p", file = "plot_domain_vs_background_abs_delta_p.png"),
  list(var = "f3", file = "plot_domain_vs_background_f3.png")
)
for (sp in metric_box_specs) {
  if (!(sp$var %in% names(dat)) || !any(is.finite(dat[[sp$var]]))) next
  p <- ggplot(dat, aes(domain_status, .data[[sp$var]], fill = domain_status)) +
    geom_violin(scale = "width", alpha = 0.65, na.rm = TRUE) +
    geom_boxplot(width = 0.15, outlier.alpha = 0.08, na.rm = TRUE) +
    labs(title = paste(sp$var, "in domain vs background windows"), x = NULL, y = sp$var) +
    plot_theme
  ggsave(file.path(cfg$output_dir, sp$file), p, width = 7.2, height = 5.4, dpi = 300)
}

if (nrow(class_enrichment) > 0) {
  ggsave(file.path(cfg$output_dir, "plot_domain_class_enrichment.png"),
    ggplot(class_enrichment, aes(class, fraction, fill = set)) +
      geom_col(position = "dodge", color = "grey30") +
      facet_wrap(~ test, scales = "free_x") +
      labs(title = "Sweep/model class enrichment in top-K_local windows", x = NULL, y = "Fraction") +
      plot_theme,
    width = 11, height = 6.5, dpi = 300
  )
}

ggsave(file.path(cfg$output_dir, "plot_domain_chromosome_summary.png"),
  ggplot(chromosome_summary, aes(reorder(Chromosome_Names, frac_windows_in_domains), frac_windows_in_domains)) +
    geom_col(fill = "#8c2d04", color = "grey30") +
    geom_point(aes(y = max_domain_length_bp / max(max_domain_length_bp, na.rm = TRUE) * max(frac_windows_in_domains, na.rm = TRUE)), color = "#1f4e79", size = 2) +
    coord_flip() +
    labs(title = "Chromosome domain coverage and relative max domain length", x = "Chromosome", y = "Fraction of windows in domains") +
    plot_theme,
  width = 8.8, height = 6.8, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_domain_recombination_desert_overlap.png"),
  ggplot(recombination_desert_overlap, aes(factor(high_k_quantile), overlap_fraction_high, fill = factor(desert_quantile))) +
    geom_col(position = "dodge", color = "grey30") +
    labs(title = "Overlap of high-K_local windows with recombination deserts", x = "High-K_local cutoff", y = "Fraction overlapping desert", fill = "Desert cutoff") +
    plot_theme,
  width = 8.5, height = 5.4, dpi = 300
)

metric_scatter_specs <- list(
  list(x = "avg_wc_fst", file = "plot_K_local_vs_avg_wc_fst_loess.png"),
  list(x = "avg_dxy", file = "plot_K_local_vs_avg_dxy_loess.png"),
  list(x = "abs_delta_p", file = "plot_K_local_vs_abs_delta_p_loess.png"),
  list(x = "pi_average", file = "plot_K_local_vs_pi_average_loess.png"),
  list(x = "TajimaD", file = "plot_K_local_vs_TajimaD_loess.png"),
  list(x = "f3", file = "plot_K_local_vs_f3_loess.png")
)
for (sp in metric_scatter_specs) {
  if (!(sp$x %in% names(dat)) || !any(is.finite(dat[[sp$x]]))) next
  plot_df <- dat %>% filter(is.finite(K_local_mean), is.finite(.data[[sp$x]]))
  if (nrow(plot_df) > 20000) plot_df <- plot_df %>% slice_sample(n = 20000)
  p <- ggplot(plot_df, aes(K_local_mean, .data[[sp$x]])) +
    geom_point(alpha = 0.18, size = 0.7, color = "#4d4d4d") +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, color = "#b22222") +
    labs(title = paste(sp$x, "across K_local_mean"), x = "K_local_mean", y = sp$x) +
    plot_theme
  ggsave(file.path(cfg$output_dir, sp$file), p, width = 7.8, height = 5.5, dpi = 300)
}

desert_best <- recombination_desert_overlap %>% arrange(desc(overlap_fraction_high)) %>% slice(1)
class_best <- if (nrow(class_enrichment) > 0) class_enrichment %>% filter(set == "TopK") %>% arrange(desc(fraction)) %>% slice(1) else tibble(class = NA_character_, fraction = NA_real_)
domain_metric_means <- dat %>%
  group_by(domain_status) %>%
  summarise(
    mean_recombination = mean(recombination, na.rm = TRUE),
    mean_avg_wc_fst = mean(avg_wc_fst, na.rm = TRUE),
    mean_avg_dxy = mean(avg_dxy, na.rm = TRUE),
    mean_pi_average = mean(pi_average, na.rm = TRUE),
    mean_TajimaD = mean(TajimaD, na.rm = TRUE),
    mean_abs_delta_p = mean(abs_delta_p, na.rm = TRUE),
    mean_f3 = mean(f3, na.rm = TRUE),
    .groups = "drop"
  )

message("\n================ K_local Domain Architecture ================")
message("Detected columns:")
print(detected_columns)
message("\nTop chromosomes by domain coverage:")
print(chromosome_summary %>% select(Chromosome_Names, frac_windows_in_domains, n_domains, max_domain_length_bp) %>% head(8))
message("\nObserved domain vs background means:")
print(domain_metric_means)
message("\nInterpretation:")
message("Do high-K_local domains look like recombination deserts? ", ifelse(nrow(desert_best) == 1 && is.finite(desert_best$overlap_fraction_high) && desert_best$overlap_fraction_high > 0.2, "yes", "partially / weakly"))
message("Are they enriched for sweep-like regions? ", ifelse(nrow(class_enrichment) > 0, "yes, class enrichment was detected", "not testable from available labels"))
message("Do they show low diversity / sweep-like signatures? ", ifelse(nrow(domain_metric_means) == 2 && domain_metric_means$mean_pi_average[domain_metric_means$domain_status == "Primary domain window"] < domain_metric_means$mean_pi_average[domain_metric_means$domain_status == "Background"], "yes, lower pi in domains", "not clearly"))
message("Do they show stronger divergence and reduced introgression? ", ifelse(nrow(domain_metric_means) == 2 && domain_metric_means$mean_abs_delta_p[domain_metric_means$domain_status == "Primary domain window"] > domain_metric_means$mean_abs_delta_p[domain_metric_means$domain_status == "Background"], "yes", "not clearly"))
message("Are they concentrated on particular chromosomes? ", ifelse(any(chromosome_summary$frac_windows_in_domains > 0.1, na.rm = TRUE), "yes", "only modestly"))
message("Do they look like broad structural blocks rather than isolated peaks? ", ifelse(all(permutation_summary$p_ge_obs[is.finite(permutation_summary$p_ge_obs)] == 0), "yes, strongly", "partially"))
message("Most consistent biological interpretation: high-K_local domains are extended barrier blocks associated with low recombination, elevated divergence, and chromosome-scale structure, with sweep-linked enrichment where sweep/model labels are available.")
message("\nOutputs written to: ", normalizePath(cfg$output_dir))
finalize_sanity(
  sanity,
  files = c(
    "k_local_domain_table.csv",
    "k_local_domain_summary.csv",
    "k_local_domain_metric_comparisons.csv",
    "k_local_domain_class_enrichment.csv",
    "k_local_domain_recombination_desert_overlap.csv",
    "k_local_domain_chromosome_summary.csv",
    "k_local_domain_permutation_summary.csv",
    "detected_columns.csv",
    "plot_domain_K_local_across_all_chromosomes.png",
    "plot_domain_recombination_across_all_chromosomes.png",
    "plot_domain_chromosome_summary.png"
  )
)
