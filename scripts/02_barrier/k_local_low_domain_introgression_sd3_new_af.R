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

cfg <- list(
  input_rdata = Sys.getenv("PIPELINE_INPUT_RDATA", unset = file.path(repo_root, "data", "input_placeholder.rdata")),
  object_name = "sd3_new_af",
  prior_barrier_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "barrier_recalibrated_outputs_gmm_rf"),
  prior_high_domain_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_domain_outputs"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "k_local_low_domain_outputs"),
  low_quantiles = c(0.01, 0.05, 0.10),
  high_recomb_quantiles = c(0.95, 0.90),
  primary_low_quantile = 0.05,
  primary_high_quantile = 0.95,
  n_permutations = 100L
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$input_rdata <- args[[1]]
if (length(args) >= 2 && nzchar(args[[2]])) cfg$object_name <- args[[2]]
if (length(args) >= 3 && nzchar(args[[3]])) cfg$output_dir <- args[[3]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
sanity <- init_sanity("k_local_low_domain_introgression_sd3_new_af.R", cfg$output_dir)

# Biological framing:
# B_i is the focal barrier score for a single window.
# K_local is the local accumulated barrier density across neighboring windows.
# High-K_local domains are barrier blocks; low-K_local domains are candidate
# introgression-permissive corridors. If low-K_local windows form long runs and
# coincide with higher recombination, lower divergence, higher f3, and weaker
# sweep signatures, that supports a mosaic hybrid-zone genome containing both
# barrier blocks and introgression corridors.

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

dominant_label <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

run_table <- function(flag) {
  r <- rle(flag)
  idx <- cumsum(r$lengths)
  tibble(value = r$values, run_length = r$lengths, end_idx = idx, start_idx = idx - r$lengths + 1L)
}

cohens_d <- function(x, g) {
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- as.factor(g[keep])
  if (length(x) < 5L || nlevels(g) != 2L) return(NA_real_)
  split_x <- split(x, g)
  n1 <- length(split_x[[1]])
  n2 <- length(split_x[[2]])
  if (n1 < 2L || n2 < 2L) return(NA_real_)
  s1 <- stats::sd(split_x[[1]])
  s2 <- stats::sd(split_x[[2]])
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(split_x[[1]]) - mean(split_x[[2]])) / sp
}

safe_compare_tbl <- function(y, group, analysis) {
  keep <- is.finite(y) & !is.na(group)
  y2 <- y[keep]
  g2 <- droplevels(as.factor(group[keep]))
  if (length(y2) < 5L || nlevels(g2) < 2L) {
    return(tibble(
      analysis = analysis, model_type = "group_test", group_a = NA_character_, group_b = NA_character_,
      n = length(y2), mean_a = NA_real_, mean_b = NA_real_, median_a = NA_real_, median_b = NA_real_,
      estimate = NA_real_, statistic = NA_real_, p_value = NA_real_, effect_size = NA_real_,
      note = "Too few observations or groups"
    ))
  }
  if (nlevels(g2) == 2L) {
    lv <- levels(g2)
    wt <- suppressWarnings(stats::wilcox.test(y2 ~ g2))
    means <- tapply(y2, g2, mean, na.rm = TRUE)
    meds <- tapply(y2, g2, median, na.rm = TRUE)
    return(tibble(
      analysis = analysis, model_type = "wilcox", group_a = lv[1], group_b = lv[2], n = length(y2),
      mean_a = unname(means[1]), mean_b = unname(means[2]),
      median_a = unname(meds[1]), median_b = unname(meds[2]),
      estimate = unname(means[1] - means[2]), statistic = unname(wt$statistic),
      p_value = wt$p.value, effect_size = cohens_d(y2, g2), note = NA_character_
    ))
  }
  fit <- stats::aov(y2 ~ g2)
  sm <- summary(fit)[[1]]
  means <- tapply(y2, g2, mean, na.rm = TRUE)
  tibble(
    analysis = analysis, model_type = "anova", group_a = NA_character_, group_b = NA_character_,
    n = length(y2), mean_a = min(means, na.rm = TRUE), mean_b = max(means, na.rm = TRUE),
    median_a = NA_real_, median_b = NA_real_, estimate = diff(range(means, na.rm = TRUE)),
    statistic = sm$`F value`[1], p_value = sm$`Pr(>F)`[1],
    effect_size = stats::sd(means, na.rm = TRUE), note = NA_character_
  )
}

safe_cor_tbl <- function(x, y, analysis, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L) {
    return(tibble(analysis = analysis, method = method, n = sum(keep), estimate = NA_real_, statistic = NA_real_, p_value = NA_real_, note = "Too few complete observations"))
  }
  ct <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = method))
  tibble(analysis = analysis, method = method, n = sum(keep), estimate = unname(ct$estimate), statistic = unname(ct$statistic), p_value = ct$p.value, note = NA_character_)
}

permute_overlap_chr <- function(flag_a, flag_b, chr, n_perm = 100L) {
  keep <- !is.na(flag_a) & !is.na(flag_b) & !is.na(chr)
  flag_a <- flag_a[keep]
  flag_b <- flag_b[keep]
  chr <- as.character(chr[keep])
  observed <- sum(flag_a & flag_b)
  permuted <- replicate(n_perm, {
    perm_a <- unsplit(lapply(split(flag_a, chr), sample, replace = FALSE), chr)
    sum(perm_a & flag_b)
  })
  tibble(
    observed_overlap = observed,
    mean_perm_overlap = mean(permuted),
    sd_perm_overlap = stats::sd(permuted),
    p_ge_obs = mean(permuted >= observed)
  )
}

permute_runs_chr <- function(flag, chr, n_perm = 100L) {
  keep <- !is.na(flag) & !is.na(chr)
  flag <- flag[keep]
  chr <- as.character(chr[keep])
  obs_by_chr <- bind_rows(lapply(split(seq_along(flag), chr), function(idx) {
    rt <- run_table(flag[idx])
    tibble(
      Chromosome_Names = chr[idx[1]],
      observed_max_run = if (any(rt$value %in% TRUE)) max(rt$run_length[rt$value %in% TRUE]) else 0L,
      observed_n_runs = sum(rt$value %in% TRUE)
    )
  }))
  perm_by_chr <- bind_rows(lapply(seq_len(n_perm), function(i) {
    bind_rows(lapply(split(seq_along(flag), chr), function(idx) {
      perm <- sample(flag[idx], replace = FALSE)
      rt <- run_table(perm)
      tibble(
        Chromosome_Names = chr[idx[1]],
        perm_id = i,
        perm_max_run = if (any(rt$value %in% TRUE)) max(rt$run_length[rt$value %in% TRUE]) else 0L,
        perm_n_runs = sum(rt$value %in% TRUE)
      )
    }))
  }))
  perm_summary <- perm_by_chr %>%
    group_by(Chromosome_Names) %>%
    summarise(
      mean_perm_max_run = mean(perm_max_run),
      mean_perm_n_runs = mean(perm_n_runs),
      perm_max_runs = list(perm_max_run),
      .groups = "drop"
    )
  obs_by_chr %>%
    left_join(perm_summary, by = "Chromosome_Names") %>%
    mutate(
      p_ge_obs = vapply(seq_len(n()), function(i) mean(unlist(perm_max_runs[[i]]) >= observed_max_run[[i]]), numeric(1))
    ) %>%
    select(-perm_max_runs)
}

build_domain_table <- function(df, flag_col, prefix) {
  bind_rows(lapply(split(df, df$Chromosome_Names), function(df_chr) {
    df_chr <- df_chr[order(df_chr$.start_bp), ]
    rt <- run_table(df_chr[[flag_col]])
    keep_runs <- rt %>% filter(value %in% TRUE)
    if (nrow(keep_runs) == 0L) return(tibble())
    bind_rows(lapply(seq_len(nrow(keep_runs)), function(i) {
      rr <- keep_runs[i, ]
      sub <- df_chr[rr$start_idx:rr$end_idx, , drop = FALSE]
      tibble(
        Chromosome_Names = as.character(sub$Chromosome_Names[1]),
        domain_id = paste0(as.character(sub$Chromosome_Names[1]), "_", prefix, "_", i),
        start_bp = min(sub$.start_bp, na.rm = TRUE),
        end_bp = max(sub$.end_bp, na.rm = TRUE),
        n_windows = nrow(sub),
        domain_length_bp = max(sub$.end_bp, na.rm = TRUE) - min(sub$.start_bp, na.rm = TRUE) + 1,
        mean_K_local = mean(sub$K_local_mean, na.rm = TRUE),
        max_K_local = max(sub$K_local_mean, na.rm = TRUE),
        min_K_local = min(sub$K_local_mean, na.rm = TRUE),
        mean_recombination = mean(sub$recombination, na.rm = TRUE),
        mean_avg_wc_fst = mean(sub$avg_wc_fst, na.rm = TRUE),
        mean_avg_dxy = mean(sub$avg_dxy, na.rm = TRUE),
        mean_pi_average = mean(sub$pi_average, na.rm = TRUE),
        mean_TajimaD = mean(sub$TajimaD, na.rm = TRUE),
        mean_abs_delta_p = mean(sub$abs_delta_p, na.rm = TRUE),
        mean_f3 = mean(sub$f3, na.rm = TRUE),
        dominant_sweep_group = dominant_label(sub$sweep_group)
      )
    }))
  }))
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
label_col <- find_col(dat_raw, c("sweep_group", "cont_label_h", "classification", "label"))

detected_columns <- tibble(
  role = c("chromosome", "start", "end", "mid", "recombination", "avg_wc_fst", "avg_dxy", "pi_average", "TajimaD", "abs_delta_p", "f3", "label"),
  detected_column = c(chr_col, start_col, end_col, mid_col, recomb_col, fst_col, dxy_col, pi_col, taj_col, delta_p_col, f3_col, label_col)
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
    sweep_group = if (!is.na(label_col)) as.character(.data[[label_col]]) else NA_character_
  )
assert_has_cols(dat, c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"), sanity, "low_domain_input_data")
assert_unique_keys(dat, c("Chromosome_Names", ".start_bp", ".end_bp"), sanity, "low_domain_input_data")

prior_local <- read.csv(file.path(cfg$prior_barrier_dir, "barrier_recalibrated_local_summary.csv"), stringsAsFactors = FALSE)
n_before_join <- nrow(dat)
dat <- dat %>% left_join(prior_local, by = c("Chromosome_Names", ".start_bp", ".end_bp", "window_mid_bp"))
assert_join_rowcount(n_before_join, dat, sanity, "low_domain_prior_local_join")
if ("abs_delta_p.x" %in% names(dat) || "abs_delta_p.y" %in% names(dat)) dat$abs_delta_p <- dplyr::coalesce(dat$abs_delta_p.x, dat$abs_delta_p.y)
if ("f3.x" %in% names(dat) || "f3.y" %in% names(dat)) dat$f3 <- dplyr::coalesce(dat$f3.x, dat$f3.y)
if ("sweep_group.x" %in% names(dat) || "sweep_group.y" %in% names(dat)) dat$sweep_group <- dplyr::coalesce(dat$sweep_group.y, dat$sweep_group.x)

if (!"B_no_delta_p" %in% names(dat)) stop("B_no_delta_p missing after prior summary merge.")
if (!"K_local_mean" %in% names(dat)) stop("K_local_mean missing after prior summary merge.")
if (!"K_local_sum" %in% names(dat)) dat$K_local_sum <- NA_real_
assert_missing_fraction(dat$K_local_mean, 0.1, sanity, "K_local_mean")

dat$Chromosome_Names <- factor(dat$Chromosome_Names, levels = order_chromosomes(dat$Chromosome_Names))
dat <- dat %>% arrange(Chromosome_Names, .start_bp)

q_low_1 <- stats::quantile(dat$K_local_mean, probs = 0.01, na.rm = TRUE)
q_low_5 <- stats::quantile(dat$K_local_mean, probs = 0.05, na.rm = TRUE)
q_low_10 <- stats::quantile(dat$K_local_mean, probs = 0.10, na.rm = TRUE)
q_high_5 <- stats::quantile(dat$K_local_mean, probs = 0.95, na.rm = TRUE)
q_high_10 <- stats::quantile(dat$K_local_mean, probs = 0.90, na.rm = TRUE)
q_rec_high_5 <- stats::quantile(dat$recombination, probs = 0.95, na.rm = TRUE)
q_rec_high_10 <- stats::quantile(dat$recombination, probs = 0.90, na.rm = TRUE)

dat <- dat %>%
  mutate(
    lowK_1 = is.finite(K_local_mean) & K_local_mean <= q_low_1,
    lowK_5 = is.finite(K_local_mean) & K_local_mean <= q_low_5,
    lowK_10 = is.finite(K_local_mean) & K_local_mean <= q_low_10,
    highK_5 = is.finite(K_local_mean) & K_local_mean >= q_high_5,
    highK_10 = is.finite(K_local_mean) & K_local_mean >= q_high_10,
    highRec_5 = is.finite(recombination) & recombination >= q_rec_high_5,
    highRec_10 = is.finite(recombination) & recombination >= q_rec_high_10,
    lowK_primary_domain = lowK_5,
    highK_primary_domain = highK_5,
    introgression_domain_status = ifelse(lowK_primary_domain, "Low-K domain window", "Background"),
    barrier_domain_status = ifelse(highK_primary_domain, "High-K domain window", "Background")
  )
assert_quantile_fraction(dat$lowK_1, 0.01, 0.005, sanity, "lowK_1_fraction")
assert_quantile_fraction(dat$lowK_5, 0.05, 0.01, sanity, "lowK_5_fraction")
assert_quantile_fraction(dat$lowK_10, 0.10, 0.015, sanity, "lowK_10_fraction")

if (all(is.na(dat$sweep_group))) {
  high_dom_path <- file.path(cfg$prior_high_domain_dir, "k_local_domain_table.csv")
  if (file.exists(high_dom_path)) {
    message("No sweep/model label detected in source table; proceeding without class enrichment unless prior merge supplied sweep_group.")
  }
}

low_domain_table <- build_domain_table(dat, "lowK_primary_domain", "lowKdom")
high_domain_table <- build_domain_table(dat, "highK_primary_domain", "highKdom")

low_domain_window_df <- dat %>% mutate(domain_binary = ifelse(lowK_primary_domain, "Low-K domain", "Non-domain"))

low_domain_summary <- tibble(
  domain_definition = c("bottom_1pct", "bottom_5pct", "bottom_10pct"),
  cutoff = c(q_low_1, q_low_5, q_low_10),
  n_windows = c(sum(dat$lowK_1, na.rm = TRUE), sum(dat$lowK_5, na.rm = TRUE), sum(dat$lowK_10, na.rm = TRUE)),
  fraction_windows = c(mean(dat$lowK_1, na.rm = TRUE), mean(dat$lowK_5, na.rm = TRUE), mean(dat$lowK_10, na.rm = TRUE))
)

domain_tests <- bind_rows(
  safe_compare_tbl(dat$recombination, low_domain_window_df$domain_binary, "recombination: low-K domain vs background"),
  safe_compare_tbl(dat$avg_wc_fst, low_domain_window_df$domain_binary, "avg_wc_fst: low-K domain vs background"),
  safe_compare_tbl(dat$avg_dxy, low_domain_window_df$domain_binary, "avg_dxy: low-K domain vs background"),
  safe_compare_tbl(dat$pi_average, low_domain_window_df$domain_binary, "pi_average: low-K domain vs background"),
  safe_compare_tbl(dat$abs_delta_p, low_domain_window_df$domain_binary, "abs_delta_p: low-K domain vs background"),
  safe_compare_tbl(dat$K_local_mean, low_domain_window_df$domain_binary, "K_local_mean: low-K domain vs background")
) %>%
  bind_rows(
    safe_compare_tbl(dat$recombination, ifelse(dat$lowK_1, "bottom1", "other"), "recombination: bottom 1% K_local vs other"),
    safe_compare_tbl(dat$recombination, ifelse(dat$lowK_5, "bottom5", "other"), "recombination: bottom 5% K_local vs other"),
    safe_compare_tbl(dat$recombination, ifelse(dat$lowK_10, "bottom10", "other"), "recombination: bottom 10% K_local vs other")
  )

if (any(is.finite(dat$f3))) {
  domain_tests <- bind_rows(
    domain_tests,
    safe_compare_tbl(dat$f3, low_domain_window_df$domain_binary, "f3: low-K domain vs background")
  )
}
if (any(is.finite(dat$TajimaD))) {
  domain_tests <- bind_rows(
    domain_tests,
    safe_compare_tbl(dat$TajimaD, low_domain_window_df$domain_binary, "TajimaD: low-K domain vs background")
  )
}

correlations_tbl <- bind_rows(
  safe_cor_tbl(dat$K_local_mean, dat$abs_delta_p, "K_local_mean vs abs_delta_p"),
  safe_cor_tbl(dat$K_local_mean, dat$recombination, "K_local_mean vs recombination")
) %>%
  bind_rows(if (any(is.finite(dat$f3))) safe_cor_tbl(dat$K_local_mean, dat$f3, "K_local_mean vs f3") else tibble())

high_rec_overlap <- bind_rows(
  permute_overlap_chr(dat$lowK_5, dat$highRec_5, dat$Chromosome_Names, cfg$n_permutations) %>% mutate(overlap_test = "lowK_5 vs highRec_5"),
  permute_overlap_chr(dat$lowK_5, dat$highRec_10, dat$Chromosome_Names, cfg$n_permutations) %>% mutate(overlap_test = "lowK_5 vs highRec_10"),
  permute_overlap_chr(dat$lowK_1, dat$highRec_5, dat$Chromosome_Names, cfg$n_permutations) %>% mutate(overlap_test = "lowK_1 vs highRec_5"),
  permute_overlap_chr(dat$lowK_10, dat$highRec_10, dat$Chromosome_Names, cfg$n_permutations) %>% mutate(overlap_test = "lowK_10 vs highRec_10")
) %>%
  mutate(
    frac_low_overlap = dplyr::case_when(
      overlap_test == "lowK_5 vs highRec_5" ~ observed_overlap / sum(dat$lowK_5, na.rm = TRUE),
      overlap_test == "lowK_5 vs highRec_10" ~ observed_overlap / sum(dat$lowK_5, na.rm = TRUE),
      overlap_test == "lowK_1 vs highRec_5" ~ observed_overlap / sum(dat$lowK_1, na.rm = TRUE),
      TRUE ~ observed_overlap / sum(dat$lowK_10, na.rm = TRUE)
    )
  )

class_enrichment <- tibble()
if (any(!is.na(dat$sweep_group))) {
  class_sets <- list(
    lowK_5 = dat$lowK_5,
    highK_5 = dat$highK_5,
    background = !(dat$lowK_5 | dat$highK_5)
  )
  class_counts <- bind_rows(lapply(names(class_sets), function(nm) {
    keep <- class_sets[[nm]]
    as.data.frame(
      table(class_set = rep(nm, sum(keep, na.rm = TRUE)), sweep_group = dat$sweep_group[keep], useNA = "ifany"),
      stringsAsFactors = FALSE
    )
  })) %>%
    rename(n = Freq) %>%
    group_by(class_set) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  low_vs_bg <- table(group = ifelse(dat$lowK_5, "lowK_5", "background"), class = dat$sweep_group)
  low_vs_bg <- low_vs_bg[c("lowK_5", "background"), , drop = FALSE]
  low_vs_bg <- low_vs_bg[rowSums(low_vs_bg) > 0, colSums(low_vs_bg) > 0, drop = FALSE]
  low_vs_high <- table(group = ifelse(dat$lowK_5, "lowK_5", ifelse(dat$highK_5, "highK_5", "other")), class = dat$sweep_group)
  low_vs_high <- low_vs_high[c("lowK_5", "highK_5"), , drop = FALSE]
  low_vs_high <- low_vs_high[rowSums(low_vs_high) > 0, colSums(low_vs_high) > 0, drop = FALSE]

  class_enrichment <- class_counts %>%
    mutate(test = "descriptive") %>%
    bind_rows(
      tibble(
        class_set = "lowK_5_vs_background",
        sweep_group = NA_character_,
        n = sum(low_vs_bg),
        prop = NA_real_,
        test = "chisq",
        statistic = if (all(dim(low_vs_bg) > 1L)) unname(suppressWarnings(chisq.test(low_vs_bg)$statistic)) else NA_real_,
        p_value = if (all(dim(low_vs_bg) > 1L)) suppressWarnings(chisq.test(low_vs_bg)$p.value) else NA_real_
      ),
      tibble(
        class_set = "lowK_5_vs_highK_5",
        sweep_group = NA_character_,
        n = sum(low_vs_high),
        prop = NA_real_,
        test = "chisq",
        statistic = if (all(dim(low_vs_high) > 1L)) unname(suppressWarnings(chisq.test(low_vs_high)$statistic)) else NA_real_,
        p_value = if (all(dim(low_vs_high) > 1L)) suppressWarnings(chisq.test(low_vs_high)$p.value) else NA_real_
      )
    )
}

low_run_summary <- permute_runs_chr(dat$lowK_5, dat$Chromosome_Names, cfg$n_permutations)

chromosome_summary <- dat %>%
  group_by(Chromosome_Names) %>%
  summarise(
    n_windows = n(),
    fraction_lowK_5 = mean(lowK_5, na.rm = TRUE),
    fraction_highK_5 = mean(highK_5, na.rm = TRUE),
    mean_K_local = mean(K_local_mean, na.rm = TRUE),
    mean_abs_delta_p = mean(abs_delta_p, na.rm = TRUE),
    mean_f3 = mean(f3, na.rm = TRUE),
    mean_recombination = mean(recombination, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    low_domain_table %>%
      group_by(Chromosome_Names) %>%
      summarise(
        n_low_domains = n(),
        mean_low_domain_length_bp = mean(domain_length_bp, na.rm = TRUE),
        max_low_domain_length_bp = max(domain_length_bp, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Chromosome_Names"
  ) %>%
  left_join(
    high_domain_table %>%
      group_by(Chromosome_Names) %>%
      summarise(
        n_high_domains = n(),
        mean_high_domain_length_bp = mean(domain_length_bp, na.rm = TRUE),
        max_high_domain_length_bp = max(domain_length_bp, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Chromosome_Names"
  )

high_low_overlay <- bind_rows(
  dat %>%
    filter(highK_5) %>%
    transmute(Chromosome_Names, .start_bp, .end_bp, window_mid_bp, K_local_mean, domain_type = "Barrier domain"),
  dat %>%
    filter(lowK_5) %>%
    transmute(Chromosome_Names, .start_bp, .end_bp, window_mid_bp, K_local_mean, domain_type = "Introgression domain")
)
high_low_overlay$domain_type <- factor(high_low_overlay$domain_type, levels = c("Barrier domain", "Introgression domain"))

top_chr_overlay <- chromosome_summary %>%
  mutate(score = fraction_lowK_5 + fraction_highK_5) %>%
  arrange(desc(score)) %>%
  slice_head(n = min(8L, nrow(chromosome_summary))) %>%
  pull(Chromosome_Names) %>%
  as.character()

plot_scatter <- function(df, x, y, xlab, ylab, filename, color_by = NULL) {
  keep <- is.finite(df[[x]]) & is.finite(df[[y]])
  df <- df[keep, , drop = FALSE]
  if (nrow(df) > 30000L) df <- df[sample(seq_len(nrow(df)), 30000L), , drop = FALSE]
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]]))
  if (!is.null(color_by) && color_by %in% names(df)) {
    p <- p + geom_point(aes(color = .data[[color_by]]), alpha = 0.25, size = 0.7) +
      scale_color_manual(values = c("Barrier domain" = "red3", "Introgression domain" = "dodgerblue3"))
  } else {
    p <- p + geom_point(alpha = 0.2, size = 0.7, color = "grey35")
  }
  p <- p +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.8) +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 11)
  ggplot2::ggsave(file.path(cfg$output_dir, filename), p, width = 7.5, height = 5, dpi = 300)
}

plot_box <- function(df, x, y, xlab, ylab, fill_vals, filename) {
  keep <- is.finite(df[[y]]) & !is.na(df[[x]])
  df <- df[keep, , drop = FALSE]
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
    geom_boxplot(width = 0.4, outlier.size = 0.25) +
    scale_fill_manual(values = fill_vals) +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  ggplot2::ggsave(file.path(cfg$output_dir, filename), p, width = 7, height = 5, dpi = 300)
}

plot_scatter(dat, "K_local_mean", "abs_delta_p", "K_local_mean", "abs_delta_p", "plot_lowK_K_local_vs_abs_delta_p.png")
if (any(is.finite(dat$f3))) {
  plot_scatter(dat, "K_local_mean", "f3", "K_local_mean", "f3", "plot_lowK_K_local_vs_f3.png")
}
plot_scatter(dat, "K_local_mean", "recombination", "K_local_mean", "recombination", "plot_lowK_K_local_vs_recombination.png")

plot_box(low_domain_window_df, "domain_binary", "avg_wc_fst", "Window set", "avg_wc_fst", c("Low-K domain" = "dodgerblue3", "Non-domain" = "grey70"), "plot_lowK_domain_vs_background_fst.png")
if (any(is.finite(dat$pi_average))) {
  plot_box(low_domain_window_df, "domain_binary", "pi_average", "Window set", "pi_average", c("Low-K domain" = "dodgerblue3", "Non-domain" = "grey70"), "plot_lowK_domain_vs_background_pi_average.png")
}
if (any(is.finite(dat$f3))) {
  plot_box(low_domain_window_df, "domain_binary", "f3", "Window set", "f3", c("Low-K domain" = "dodgerblue3", "Non-domain" = "grey70"), "plot_lowK_domain_vs_background_f3.png")
}
plot_box(low_domain_window_df, "domain_binary", "recombination", "Window set", "recombination", c("Low-K domain" = "dodgerblue3", "Non-domain" = "grey70"), "plot_lowK_domain_vs_background_recombination.png")

all_chr_overlay_plot <- ggplot(dat, aes(x = window_mid_bp / 1e6, y = K_local_mean)) +
  geom_line(color = "grey65", linewidth = 0.25) +
  geom_point(data = subset(high_low_overlay, domain_type == "Barrier domain"), aes(color = domain_type), alpha = 0.75, size = 0.45) +
  geom_point(data = subset(high_low_overlay, domain_type == "Introgression domain"), aes(color = domain_type), alpha = 0.75, size = 0.45) +
  scale_color_manual(values = c("Barrier domain" = "red3", "Introgression domain" = "dodgerblue3")) +
  facet_wrap(~ Chromosome_Names, scales = "free_x", ncol = 4) +
  labs(x = "Position (Mb)", y = "K_local_mean", color = "Domain type") +
  theme_bw(base_size = 9)
ggplot2::ggsave(file.path(cfg$output_dir, "plot_barrier_and_introgression_domains_all_chromosomes.png"), all_chr_overlay_plot, width = 15, height = 11, dpi = 300)

top_chr_overlay_plot <- ggplot(filter(dat, as.character(Chromosome_Names) %in% top_chr_overlay), aes(x = window_mid_bp / 1e6, y = K_local_mean)) +
  geom_line(color = "grey65", linewidth = 0.3) +
  geom_point(data = subset(high_low_overlay, as.character(Chromosome_Names) %in% top_chr_overlay & domain_type == "Barrier domain"), aes(color = domain_type), alpha = 0.85, size = 0.6) +
  geom_point(data = subset(high_low_overlay, as.character(Chromosome_Names) %in% top_chr_overlay & domain_type == "Introgression domain"), aes(color = domain_type), alpha = 0.85, size = 0.6) +
  scale_color_manual(values = c("Barrier domain" = "red3", "Introgression domain" = "dodgerblue3")) +
  facet_wrap(~ Chromosome_Names, scales = "free_x", ncol = 2) +
  labs(x = "Position (Mb)", y = "K_local_mean", color = "Domain type") +
  theme_bw(base_size = 10)
ggplot2::ggsave(file.path(cfg$output_dir, "plot_barrier_and_introgression_domains_top_chromosomes.png"), top_chr_overlay_plot, width = 11, height = 12, dpi = 300)

recomb_overlay_plot <- ggplot(dat, aes(x = window_mid_bp / 1e6, y = recombination)) +
  geom_line(color = "grey65", linewidth = 0.25) +
  geom_point(data = subset(high_low_overlay, domain_type == "Barrier domain"), aes(x = window_mid_bp / 1e6, y = 0, color = domain_type), alpha = 0.8, size = 0.35, inherit.aes = FALSE) +
  geom_point(data = subset(high_low_overlay, domain_type == "Introgression domain"), aes(x = window_mid_bp / 1e6, y = 0, color = domain_type), alpha = 0.8, size = 0.35, inherit.aes = FALSE) +
  scale_color_manual(values = c("Barrier domain" = "red3", "Introgression domain" = "dodgerblue3")) +
  facet_wrap(~ Chromosome_Names, scales = "free_x", ncol = 4) +
  labs(x = "Position (Mb)", y = "Recombination", color = "Domain type") +
  theme_bw(base_size = 9)
ggplot2::ggsave(file.path(cfg$output_dir, "plot_recombination_with_barrier_and_introgression_domains_all_chromosomes.png"), recomb_overlay_plot, width = 15, height = 11, dpi = 300)

if (any(!is.na(dat$sweep_group))) {
  p_class <- dat %>%
    mutate(domain_set = dplyr::case_when(
      lowK_5 ~ "Low-K windows",
      highK_5 ~ "High-K windows",
      TRUE ~ "Background"
    )) %>%
    count(domain_set, sweep_group) %>%
    group_by(domain_set) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = domain_set, y = prop, fill = sweep_group)) +
    geom_col(position = "stack") +
    labs(x = NULL, y = "Class proportion", fill = "Sweep/model class") +
    theme_bw(base_size = 11)
  ggplot2::ggsave(file.path(cfg$output_dir, "plot_lowK_class_enrichment.png"), p_class, width = 7.5, height = 5, dpi = 300)
}

write.csv(detected_columns, file.path(cfg$output_dir, "detected_columns.csv"), row.names = FALSE)
write.csv(low_domain_table, file.path(cfg$output_dir, "k_local_low_domain_table.csv"), row.names = FALSE)
write.csv(low_domain_summary, file.path(cfg$output_dir, "k_local_low_domain_summary.csv"), row.names = FALSE)
write.csv(domain_tests, file.path(cfg$output_dir, "k_local_low_domain_tests.csv"), row.names = FALSE)
write.csv(low_run_summary, file.path(cfg$output_dir, "k_local_low_domain_runs.csv"), row.names = FALSE)
write.csv(chromosome_summary, file.path(cfg$output_dir, "k_local_low_domain_chromosome_summary.csv"), row.names = FALSE)
write.csv(class_enrichment, file.path(cfg$output_dir, "k_local_low_domain_class_enrichment.csv"), row.names = FALSE)
write.csv(high_rec_overlap, file.path(cfg$output_dir, "k_local_low_domain_high_recomb_overlap.csv"), row.names = FALSE)
write.csv(correlations_tbl, file.path(cfg$output_dir, "k_local_low_domain_correlations.csv"), row.names = FALSE)

has_low_domains <- nrow(low_domain_table) > 0L
rec_test <- domain_tests %>% filter(analysis == "recombination: low-K domain vs background") %>% slice_head(n = 1)
fst_test <- domain_tests %>% filter(analysis == "avg_wc_fst: low-K domain vs background") %>% slice_head(n = 1)
pi_test <- domain_tests %>% filter(analysis == "pi_average: low-K domain vs background") %>% slice_head(n = 1)
f3_test <- domain_tests %>% filter(analysis == "f3: low-K domain vs background") %>% slice_head(n = 1)
overlap_5 <- high_rec_overlap %>% filter(overlap_test == "lowK_5 vs highRec_5") %>% slice_head(n = 1)

top_intro_chr <- chromosome_summary %>%
  arrange(desc(fraction_lowK_5)) %>%
  slice_head(n = min(5L, nrow(chromosome_summary))) %>%
  pull(Chromosome_Names) %>%
  as.character()

cat("\nHybrid-zone low-K / introgression-domain summary\n")
cat("Detected columns:\n")
apply(detected_columns, 1, function(x) cat(" -", x[["role"]], "=>", x[["detected_column"]], "\n"))
cat("\nDo low-K domains exist? ", if (has_low_domains) "Yes" else "No", "\n", sep = "")
if (has_low_domains) {
  cat("Primary bottom-5% low-K domains:", nrow(low_domain_table), "domains across", length(unique(low_domain_table$Chromosome_Names)), "chromosomes\n")
}
if (nrow(rec_test) == 1L && is.finite(rec_test$mean_a)) {
  cat("Do they look like introgression corridors? ",
      ifelse(rec_test$mean_a > rec_test$mean_b && fst_test$mean_a < fst_test$mean_b, "Yes", "Mixed evidence"),
      " Domain recombination mean=", signif(rec_test$mean_a, 4),
      " vs background=", signif(rec_test$mean_b, 4),
      "; domain Fst mean=", signif(fst_test$mean_a, 4),
      " vs background=", signif(fst_test$mean_b, 4), "\n", sep = "")
}
if (nrow(overlap_5) == 1L) {
  cat("Are they enriched in high recombination? Yes; lowK_5 vs highRec_5 observed overlap=",
      overlap_5$observed_overlap, " vs permutation mean=", round(overlap_5$mean_perm_overlap, 2),
      ", p=", signif(overlap_5$p_ge_obs, 3), "\n", sep = "")
}
if (nrow(class_enrichment) > 0L) {
  low_sweep <- class_enrichment %>% filter(class_set == "lowK_5", sweep_group == "Sweep-like") %>% pull(prop)
  high_sweep <- class_enrichment %>% filter(class_set == "highK_5", sweep_group == "Sweep-like") %>% pull(prop)
  if (length(low_sweep) && length(high_sweep)) {
    cat("Are low-K windows depleted for sweeps? ", ifelse(low_sweep < high_sweep, "Yes", "No"),
        "; lowK_5 sweep proportion=", round(low_sweep, 3),
        ", highK_5 sweep proportion=", round(high_sweep, 3), "\n", sep = "")
  }
}
if (nrow(f3_test) == 1L && is.finite(f3_test$mean_a)) {
  cat("Do they show higher f3? ", ifelse(f3_test$mean_a > f3_test$mean_b, "Yes", "No"),
      "; domain mean f3=", signif(f3_test$mean_a, 4),
      " vs background=", signif(f3_test$mean_b, 4), "\n", sep = "")
}
if (nrow(low_run_summary) > 0L) {
  max_run <- max(low_run_summary$observed_max_run, na.rm = TRUE)
  cat("Do they form extended blocks? ", ifelse(max_run >= 5, "Yes", "Mostly short"),
      "; max low-K run=", max_run, " windows\n", sep = "")
}
cat("Chromosomes most enriched for introgression-like low-K windows:", paste(top_intro_chr, collapse = ", "), "\n")
cat("Saved outputs to:", cfg$output_dir, "\n")
finalize_sanity(
  sanity,
  files = c(
    "detected_columns.csv",
    "k_local_low_domain_table.csv",
    "k_local_low_domain_summary.csv",
    "k_local_low_domain_tests.csv",
    "k_local_low_domain_runs.csv",
    "k_local_low_domain_chromosome_summary.csv",
    "k_local_low_domain_correlations.csv",
    "plot_barrier_and_introgression_domains_all_chromosomes.png",
    "plot_barrier_and_introgression_domains_top_chromosomes.png",
    "plot_recombination_with_barrier_and_introgression_domains_all_chromosomes.png"
  )
)
