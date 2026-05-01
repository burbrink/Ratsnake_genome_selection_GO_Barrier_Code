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
  prior_output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_gmm_rf"),
  label_rdata = file.path(repo_root, "outputs", "selection", "barrier_label_objects.rdata"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "hybrid_zone_condition_outputs"),
  alpha = 0.5,
  alpha_calibration_csv = NULL,
  alpha_calibration_p_col = "P_overall_mean",
  alpha_calibration_k_col = "K_genome_mean",
  rolling_window_n = 25L,
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
sanity <- init_sanity("hybrid_zone_condition_sd3_new_af.R", cfg$output_dir)

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

roll_mean_center <- function(x, n) data.table::frollmean(x, n = n, align = "center", na.rm = FALSE)
roll_sum_center <- function(x, n) data.table::frollsum(x, n = n, align = "center", na.rm = FALSE)

quantile_class <- function(x, probs, labels) {
  br <- unique(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (length(br) < 2L || sum(ok) == 0L) return(out)
  out[ok] <- as.character(cut(x[ok], breaks = br, include.lowest = TRUE, labels = labels[seq_len(length(br) - 1L)]))
  out
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

safe_lm_tbl <- function(formula, data, analysis) {
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

safe_group_test_tbl <- function(y, group, analysis) {
  keep <- is.finite(y) & !is.na(group)
  y2 <- y[keep]
  g2 <- droplevels(as.factor(group[keep]))
  if (length(y2) < 5L || nlevels(g2) < 2L) {
    return(tibble(
      analysis = analysis, model_type = "group_test", term = NA_character_,
      estimate = NA_real_, statistic = NA_real_, p_value = NA_real_, n = length(y2),
      r_squared = NA_real_, adj_r_squared = NA_real_, aic = NA_real_, note = "Too few observations or groups"
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
assert_has_cols(dat, c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"), sanity, "hybrid_zone_input_data")
assert_unique_keys(dat, c("Chromosome_Names", ".start_bp", ".end_bp"), sanity, "hybrid_zone_input_data")

prior_local_path <- file.path(cfg$prior_output_dir, "barrier_recalibrated_local_summary.csv")
prior_chr_path <- file.path(cfg$prior_output_dir, "barrier_recalibrated_chromosome_summary.csv")
prior_genome_path <- file.path(cfg$prior_output_dir, "barrier_recalibrated_genome_summary.csv")

if (file.exists(prior_local_path)) {
  local_prev <- read.csv(prior_local_path, stringsAsFactors = FALSE)
  n_before_join <- nrow(dat)
  dat <- dat %>%
    left_join(local_prev, by = c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"))
  assert_join_rowcount(n_before_join, dat, sanity, "hybrid_zone_prior_local_join")
  if ("abs_delta_p.x" %in% names(dat) || "abs_delta_p.y" %in% names(dat)) {
    dat$abs_delta_p <- dplyr::coalesce(dat$abs_delta_p.x, dat$abs_delta_p.y)
  }
  if ("f3.x" %in% names(dat) || "f3.y" %in% names(dat)) {
    dat$f3 <- dplyr::coalesce(dat$f3.x, dat$f3.y)
  }
} else {
  dat <- dat %>%
    mutate(
      z_recombination = z_score(recombination),
      z_avg_wc_fst = z_score(avg_wc_fst),
      z_avg_dxy = z_score(avg_dxy),
      z_pi_average = z_score(pi_average),
      z_TajimaD = z_score(TajimaD),
      B_no_delta_p = -z_recombination + z_avg_wc_fst + z_avg_dxy - z_pi_average - z_TajimaD
    ) %>%
    group_by(Chromosome_Names) %>%
    arrange(.start_bp, .by_group = TRUE) %>%
    mutate(
      K_local_mean = roll_mean_center(B_no_delta_p, cfg$rolling_window_n),
      K_local_sum = roll_sum_center(B_no_delta_p, cfg$rolling_window_n)
    ) %>%
    ungroup() %>%
    mutate(P_local_mean = exp(-cfg$alpha * K_local_mean))
}

if (file.exists(cfg$label_rdata)) {
  load(cfg$label_rdata)
  if (exists("agree_dat", inherits = FALSE)) {
    maj <- apply(agree_dat[, c("gmm_label", "rf_label", "cont_label_h")], 1, function(v) {
      v <- v[!is.na(v) & nzchar(v)]
      if (!length(v)) return(NA_character_)
      tt <- sort(table(v), decreasing = TRUE)
      if (as.integer(tt[1]) >= 2L) names(tt)[1] else NA_character_
    })
    label_tbl <- tibble(
      Chromosome_Names = agree_dat$Chromosome_Names,
      .start_bp = agree_dat$.start_bp,
      consensus_label = maj
    ) %>%
      mutate(
        sweep_group = case_when(
          consensus_label == "Geographic sweep" ~ "Sweep-like",
          consensus_label == "Balancing selection" ~ "Balancing-like",
          consensus_label == "Neutral / equilibrium" ~ "Neutral-like",
          TRUE ~ NA_character_
        )
      )
    dat <- dat %>% left_join(label_tbl, by = c("Chromosome_Names", ".start_bp"))
  }
}

if (!"B_no_delta_p" %in% names(dat)) stop("B_no_delta_p not found or could not be recreated.")
if (!"K_local_mean" %in% names(dat)) stop("K_local_mean not found or could not be recreated.")
if (!"K_local_sum" %in% names(dat)) dat$K_local_sum <- NA_real_
if (!"P_local_mean" %in% names(dat)) dat$P_local_mean <- exp(-cfg$alpha * dat$K_local_mean)
dat$log10_P_local_mean <- log10(dat$P_local_mean)
assert_missing_fraction(dat$B_no_delta_p, 0.1, sanity, "B_no_delta_p")
assert_missing_fraction(dat$K_local_mean, 0.1, sanity, "K_local_mean")
assert_probability_bounds(dat$P_local_mean, tracker = sanity, label = "P_local_mean")

dat$K_local_quartile <- quantile_class(dat$K_local_mean, c(0, 0.25, 0.75, 1), c("Low", "Intermediate", "High"))
dat$K_local_tertile <- quantile_class(dat$K_local_mean, c(0, 1/3, 2/3, 1), c("Low", "Intermediate", "High"))

finite_B <- is.finite(dat$B_no_delta_p)
finite_K <- is.finite(dat$K_local_mean)
finite_P <- is.finite(dat$P_local_mean)

K_genome_mean <- mean(dat$K_local_mean, na.rm = TRUE)
K_genome_sum <- sum(pmax(dat$B_no_delta_p, 0), na.rm = TRUE)
P_overall_mean <- exp(-cfg$alpha * K_genome_mean)
log10_P_overall_mean <- log10(P_overall_mean)

genome_summary <- tibble(
  alpha = cfg$alpha,
  n_windows = nrow(dat),
  n_finite_B = sum(finite_B),
  n_finite_K_local = sum(finite_K),
  K_genome_mean = K_genome_mean,
  K_genome_sum = K_genome_sum,
  P_overall_mean = P_overall_mean,
  log10_P_overall_mean = log10_P_overall_mean,
  K_local_variance = stats::var(dat$K_local_mean, na.rm = TRUE),
  K_local_sd = stats::sd(dat$K_local_mean, na.rm = TRUE),
  K_local_IQR = stats::IQR(dat$K_local_mean, na.rm = TRUE),
  P_local_variance = stats::var(dat$P_local_mean, na.rm = TRUE)
)

class_summary <- bind_rows(
  dat %>%
    filter(!is.na(K_local_quartile)) %>%
    count(class_system = "quartile", class = K_local_quartile, name = "n_windows") %>%
    mutate(fraction = n_windows / sum(n_windows)),
  dat %>%
    filter(!is.na(K_local_tertile)) %>%
    count(class_system = "tertile", class = K_local_tertile, name = "n_windows") %>%
    mutate(fraction = n_windows / sum(n_windows))
)

chrom_levels <- order_chromosomes(dat$Chromosome_Names)
dat$Chromosome_Names <- factor(dat$Chromosome_Names, levels = chrom_levels)

chromosome_summary <- dat %>%
  group_by(Chromosome_Names) %>%
  summarise(
    n_windows = n(),
    mean_B_chr = mean(B_no_delta_p, na.rm = TRUE),
    mean_K_local = mean(K_local_mean, na.rm = TRUE),
    var_K_local = stats::var(K_local_mean, na.rm = TRUE),
    sd_K_local = stats::sd(K_local_mean, na.rm = TRUE),
    frac_low_quartile = mean(K_local_quartile == "Low", na.rm = TRUE),
    frac_intermediate_quartile = mean(K_local_quartile == "Intermediate", na.rm = TRUE),
    frac_high_quartile = mean(K_local_quartile == "High", na.rm = TRUE),
    mean_abs_delta_p = mean(abs_delta_p, na.rm = TRUE),
    mean_f3 = mean(f3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_B_chr))

model_tests <- bind_rows(
  safe_cor_tbl(dat$K_local_mean, dat$abs_delta_p, "abs_delta_p ~ K_local_mean", "pearson"),
  safe_cor_tbl(dat$K_local_mean, dat$abs_delta_p, "abs_delta_p ~ K_local_mean", "spearman"),
  safe_lm_tbl(abs_delta_p ~ K_local_mean, dat, "abs_delta_p ~ K_local_mean"),
  safe_group_test_tbl(dat$abs_delta_p, dat$K_local_quartile, "abs_delta_p by K_local quartile"),
  safe_group_test_tbl(dat$abs_delta_p, dat$K_local_tertile, "abs_delta_p by K_local tertile")
) 

if ("f3" %in% names(dat) && any(is.finite(dat$f3))) {
  model_tests <- bind_rows(
    model_tests,
    safe_cor_tbl(dat$K_local_mean, dat$f3, "f3 ~ K_local_mean", "pearson"),
    safe_cor_tbl(dat$K_local_mean, dat$f3, "f3 ~ K_local_mean", "spearman"),
    safe_lm_tbl(f3 ~ K_local_mean, dat, "f3 ~ K_local_mean"),
    safe_group_test_tbl(dat$f3, dat$K_local_quartile, "f3 by K_local quartile"),
    safe_group_test_tbl(dat$f3, dat$K_local_tertile, "f3 by K_local tertile")
  )
}

if ("sweep_group" %in% names(dat) && any(!is.na(dat$sweep_group))) {
  ct <- table(dat$sweep_group, dat$K_local_quartile, useNA = "no")
  if (all(dim(ct) > 1L)) {
    chi <- suppressWarnings(chisq.test(ct))
    model_tests <- bind_rows(
      model_tests,
      tibble(
        analysis = "sweep/model label enrichment across K_local quartiles",
        model_type = "chisq",
        term = "sweep_group x K_local_quartile",
        estimate = NA_real_,
        statistic = unname(chi$statistic),
        p_value = chi$p.value,
        n = sum(ct),
        r_squared = NA_real_,
        adj_r_squared = NA_real_,
        aic = NA_real_,
        note = NA_character_
      )
    )
  }
}

write.csv(detected_columns, file.path(cfg$output_dir, "detected_columns.csv"), row.names = FALSE)
write.csv(genome_summary, file.path(cfg$output_dir, "hybrid_zone_condition_genome_summary.csv"), row.names = FALSE)
write.csv(chromosome_summary, file.path(cfg$output_dir, "hybrid_zone_condition_chromosome_summary.csv"), row.names = FALSE)
write.csv(class_summary, file.path(cfg$output_dir, "hybrid_zone_condition_local_classes.csv"), row.names = FALSE)
write.csv(model_tests, file.path(cfg$output_dir, "hybrid_zone_condition_model_tests.csv"), row.names = FALSE)

plot_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(cfg$output_dir, "plot_hybrid_zone_K_local_distribution.png"),
  ggplot(dat, aes(K_local_mean)) +
    geom_histogram(bins = 60, fill = "#2b6c8f", color = "white", na.rm = TRUE) +
    labs(title = "Distribution of local barrier density", x = "K_local_mean", y = "Windows") +
    plot_theme,
  width = 7.5, height = 5.5, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_hybrid_zone_P_local_distribution.png"),
  ggplot(dat, aes(P_local_mean)) +
    geom_histogram(bins = 60, fill = "#b22222", color = "white", na.rm = TRUE) +
    labs(title = "Distribution of local compatibility index", x = "P_local_mean", y = "Windows") +
    plot_theme,
  width = 7.5, height = 5.5, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_hybrid_zone_abs_delta_vs_K_local.png"),
  ggplot(dat, aes(K_local_mean, abs_delta_p)) +
    geom_point(alpha = 0.25, size = 0.8, color = "#5a2d82") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#1f4e79") +
    labs(title = "Ancestry differentiation vs local barrier density", x = "K_local_mean", y = "abs_delta_p") +
    plot_theme,
  width = 7.5, height = 5.5, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_hybrid_zone_abs_delta_by_K_class.png"),
  ggplot(dat %>% filter(!is.na(K_local_quartile)), aes(K_local_quartile, abs_delta_p, fill = K_local_quartile)) +
    geom_boxplot(outlier.alpha = 0.15) +
    scale_fill_manual(values = c("Low" = "#91bfdb", "Intermediate" = "#fee090", "High" = "#fc8d59")) +
    labs(title = "Ancestry differentiation by local barrier class", x = "K_local quartile class", y = "abs_delta_p") +
    plot_theme,
  width = 7, height = 5.2, dpi = 300
)

if ("f3" %in% names(dat) && any(is.finite(dat$f3))) {
  ggsave(file.path(cfg$output_dir, "plot_hybrid_zone_f3_vs_K_local.png"),
    ggplot(dat, aes(K_local_mean, f3)) +
      geom_point(alpha = 0.25, size = 0.8, color = "#238b45") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#8c2d04") +
      labs(title = "f3 vs local barrier density", x = "K_local_mean", y = "f3") +
      plot_theme,
    width = 7.5, height = 5.5, dpi = 300
  )

  ggsave(file.path(cfg$output_dir, "plot_hybrid_zone_f3_by_K_class.png"),
    ggplot(dat %>% filter(!is.na(K_local_quartile)), aes(K_local_quartile, f3, fill = K_local_quartile)) +
      geom_boxplot(outlier.alpha = 0.15) +
      scale_fill_manual(values = c("Low" = "#91bfdb", "Intermediate" = "#fee090", "High" = "#fc8d59")) +
      labs(title = "f3 by local barrier class", x = "K_local quartile class", y = "f3") +
      plot_theme,
    width = 7, height = 5.2, dpi = 300
  )
}

top_chr <- chromosome_summary %>%
  slice_max(order_by = var_K_local, n = min(cfg$top_n_chromosomes_for_plots, nrow(chromosome_summary))) %>%
  pull(Chromosome_Names) %>%
  as.character()

ggsave(file.path(cfg$output_dir, "plot_hybrid_zone_K_local_across_chromosomes.png"),
  ggplot(dat %>% filter(as.character(Chromosome_Names) %in% top_chr), aes(window_mid_bp / 1e6, K_local_mean)) +
    geom_line(color = "#8c2d04", linewidth = 0.35, na.rm = TRUE) +
    facet_wrap(~ Chromosome_Names, scales = "free_x") +
    labs(title = "Local barrier density along top heterogeneous chromosomes", x = "Position (Mb)", y = "K_local_mean") +
    plot_theme,
  width = 11, height = 7.5, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_hybrid_zone_P_local_across_chromosomes.png"),
  ggplot(dat %>% filter(as.character(Chromosome_Names) %in% top_chr), aes(window_mid_bp / 1e6, P_local_mean)) +
    geom_line(color = "#08519c", linewidth = 0.35, na.rm = TRUE) +
    facet_wrap(~ Chromosome_Names, scales = "free_x") +
    labs(title = "Local compatibility along top heterogeneous chromosomes", x = "Position (Mb)", y = "P_local_mean") +
    plot_theme,
  width = 11, height = 7.5, dpi = 300
)

chrom_long <- chromosome_summary %>%
  select(Chromosome_Names, mean_K_local, var_K_local) %>%
  pivot_longer(cols = c(mean_K_local, var_K_local), names_to = "metric", values_to = "value")

ggsave(file.path(cfg$output_dir, "plot_hybrid_zone_chromosome_mean_and_variance.png"),
  ggplot(chrom_long, aes(Chromosome_Names, value, fill = metric)) +
    geom_col(position = "dodge", color = "grey30") +
    labs(title = "Chromosome-level mean and variance of local barrier density", x = "Chromosome", y = "Value") +
    plot_theme,
  width = 11, height = 5.8, dpi = 300
)

ggsave(file.path(cfg$output_dir, "plot_hybrid_zone_top_chromosome_segments.png"),
  ggplot(dat %>% filter(as.character(Chromosome_Names) %in% head(top_chr, 4), !is.na(K_local_quartile)),
         aes(window_mid_bp / 1e6, K_local_mean, color = K_local_quartile)) +
    geom_line(linewidth = 0.35, na.rm = TRUE) +
    facet_wrap(~ Chromosome_Names, scales = "free_x") +
    scale_color_manual(values = c("Low" = "#91bfdb", "Intermediate" = "#fee090", "High" = "#fc8d59")) +
    labs(title = "Low/intermediate/high local barrier segments", x = "Position (Mb)", y = "K_local_mean", color = "K_local class") +
    plot_theme,
  width = 11, height = 7, dpi = 300
)

low_frac <- class_summary %>% filter(class_system == "quartile", class == "Low") %>% pull(fraction)
high_frac <- class_summary %>% filter(class_system == "quartile", class == "High") %>% pull(fraction)
rho_abs_k <- model_tests %>% filter(analysis == "abs_delta_p ~ K_local_mean", model_type == "cor_spearman") %>% pull(estimate)
rho_f3_k <- model_tests %>% filter(analysis == "f3 ~ K_local_mean", model_type == "cor_spearman") %>% pull(estimate)

high_genome_barrier <- is.finite(log10_P_overall_mean) && log10_P_overall_mean < -1
heterogeneous_local <- is.finite(genome_summary$K_local_sd) && genome_summary$K_local_sd > 0.5
non_negligible_local_compatibility <- mean(dat$P_local_mean > 0.1, na.rm = TRUE) > 0.05
low_barrier_more_mixing <- is.finite(rho_abs_k) && rho_abs_k > 0
mosaic_regime <- high_genome_barrier && heterogeneous_local && non_negligible_local_compatibility && low_barrier_more_mixing

message("\n================ Hybrid-Zone Condition Summary ================")
message("Using alpha = ", cfg$alpha, " from prior barrier calibration.")
message("Detected columns:")
print(detected_columns)
message("\nGenome-wide condition:")
print(genome_summary)
message("\nTop chromosomes by mean barrier:")
print(head(chromosome_summary %>% select(Chromosome_Names, mean_B_chr, mean_K_local, var_K_local), 5))
message("\nMost heterogeneous chromosomes:")
print(chromosome_summary %>% arrange(desc(var_K_local)) %>% select(Chromosome_Names, var_K_local, frac_low_quartile, frac_high_quartile) %>% head(5))
message("\nPlain-language assessment:")
message("K_genome is large: ", ifelse(high_genome_barrier, "yes", "no"))
message("K_local is heterogeneous across the genome: ", ifelse(heterogeneous_local, "yes", "no"))
message("Some regions retain non-negligible local compatibility: ", ifelse(non_negligible_local_compatibility, "yes", "no"))
message("Low-K_local regions show more ancestry mixing: ", ifelse(low_barrier_more_mixing, "yes", "no"))
if (is.finite(rho_f3_k)) {
  message("Observed K_local vs f3 association (Spearman): ", signif(rho_f3_k, 4))
}
message("The data support a mosaic hybrid-zone regime: ", ifelse(mosaic_regime, "yes", "no"))
message("\nPaper-style interpretation:")
message(
  ifelse(
    mosaic_regime,
    "The data support a hybrid-zone regime in which genome-wide barrier density is high enough to make full panmixia unlikely, while local barrier density remains heterogeneous enough to permit partial ancestry mixing in subsets of the genome.",
    "The data do not fully support a mosaic hybrid-zone regime under the current barrier formulation, because one or more of the expected conditions on genome-wide barrier density, local heterogeneity, or local ancestry mixing were not met."
  )
)
message("\nOutputs written to: ", normalizePath(cfg$output_dir))
finalize_sanity(
  sanity,
  files = c(
    "detected_columns.csv",
    "hybrid_zone_condition_genome_summary.csv",
    "hybrid_zone_condition_chromosome_summary.csv",
    "hybrid_zone_condition_local_classes.csv",
    "hybrid_zone_condition_model_tests.csv",
    "plot_hybrid_zone_K_local_distribution.png",
    "plot_hybrid_zone_P_local_distribution.png",
    "plot_hybrid_zone_abs_delta_vs_K_local.png",
    "plot_hybrid_zone_abs_delta_by_K_class.png",
    "plot_hybrid_zone_K_local_across_chromosomes.png",
    "plot_hybrid_zone_P_local_across_chromosomes.png",
    "plot_hybrid_zone_chromosome_mean_and_variance.png"
  )
)
