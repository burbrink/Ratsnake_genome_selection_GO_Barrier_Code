#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

selection_output_dir <- if (length(args) >= 1 && nzchar(args[[1]])) {
  args[[1]]
} else {
  file.path(repo_root, "outputs", "selection")
}

barrier_outputs_root <- if (length(args) >= 2 && nzchar(args[[2]])) {
  args[[2]]
} else {
  file.path(repo_root, "outputs", "barrier_analysis_outputs")
}

selection_decision_path <- file.path(selection_output_dir, "revised_model_decision_table.tsv")
label_rdata <- file.path(selection_output_dir, "barrier_label_objects.rdata")
go_output_dir <- file.path(barrier_outputs_root, "go_gene_map_outputs")

must_exist <- c(
  selection_decision_path,
  label_rdata,
  file.path(go_output_dir, "selected_windows_by_tier.tsv"),
  file.path(go_output_dir, "selected_gene_counts_by_tier.tsv"),
  file.path(go_output_dir, "GO_gprofiler_results_by_tier.tsv"),
  file.path(go_output_dir, "immune_gene_enrichment_fisher_by_tier.tsv")
)

missing_files <- must_exist[!file.exists(must_exist)]
if (length(missing_files)) {
  stop(
    "Missing required files for decision-table update:\n",
    paste(sprintf("  - %s", missing_files), collapse = "\n"),
    call. = FALSE
  )
}

tier_map <- data.table(
  tier = c(
    "GMM/RF",
    "Old best continuous 3-of-3",
    "Continuous-distance 3-class",
    "GMM/RF + continuous-distance agreement",
    "Empirical FDR <= 0.05",
    "Global shift selected windows",
    "Genome-wide FDR <= 0.05",
    "Distance margin >= 0.25",
    "Distance margin >= 0.50",
    "Distance margin >= 1.00",
    "Score margin >= 0.25",
    "Score margin >= 0.50",
    "Score margin >= 1.00",
    "Top 1%",
    "Top 500"
  ),
  short_tier = c(
    "gmm_rf",
    "best_continuous_3of3",
    "continuous_distance_3class",
    "gmm_rf_contdist",
    "empirical_fdr_selected",
    "global_shift_selected",
    "genomewide_fdr_selected",
    "dist_margin_025",
    "dist_margin_050",
    "dist_margin_100",
    "score_margin_025",
    "score_margin_050",
    "score_margin_100",
    "dist_top1pct",
    "dist_top500"
  ),
  barrier_dir = file.path(
    barrier_outputs_root,
    c(
      "barrier_recalibrated_outputs_gmm_rf",
      "barrier_recalibrated_outputs_best_continuous_3of3",
      "barrier_recalibrated_outputs_continuous_distance_3class",
      "barrier_recalibrated_outputs_gmm_rf_contdist",
      "barrier_recalibrated_outputs_empirical_fdr_selected",
      "barrier_recalibrated_outputs_global_shift_selected",
      "barrier_recalibrated_outputs_genomewide_fdr_selected",
      "barrier_recalibrated_outputs_dist_margin_025",
      "barrier_recalibrated_outputs_dist_margin_050",
      "barrier_recalibrated_outputs_dist_margin_100",
      "barrier_recalibrated_outputs_score_margin_025",
      "barrier_recalibrated_outputs_score_margin_050",
      "barrier_recalibrated_outputs_score_margin_100",
      "barrier_recalibrated_outputs_dist_top1pct",
      "barrier_recalibrated_outputs_dist_top500"
    )
  )
)

decision <- as.data.table(fread(selection_decision_path))
decision <- decision[, .(
  tier,
  neutral_windows,
  balancing_windows,
  sweep_windows,
  ambiguous_no_consensus_windows,
  recommended_paper_use
)]
setkey(decision, tier)

selected_windows <- as.data.table(fread(file.path(go_output_dir, "selected_windows_by_tier.tsv")))
selected_genes <- as.data.table(fread(file.path(go_output_dir, "selected_gene_counts_by_tier.tsv")))
go_results <- as.data.table(fread(file.path(go_output_dir, "GO_gprofiler_results_by_tier.tsv")))
immune_results <- as.data.table(fread(file.path(go_output_dir, "immune_gene_enrichment_fisher_by_tier.tsv")))

selected_windows[, `:=`(
  Chromosome_Names = as.character(Chromosome_Names),
  window_start = as.integer(window_start),
  window_end = as.integer(window_end)
)]

selected_genes[, class := fifelse(class == "sweep", "Geographic sweep", "Balancing selection")]
gene_summary <- selected_genes[
  tier %chin% tier_map$short_tier,
  .(unique_gene_ids = max(unique_gene_ids, na.rm = TRUE)),
  by = .(short_tier = tier, class)
]

go_summary <- go_results[
  subset == "whole_genome" & is.finite(p_value) & p_value <= 0.05,
  .(number_of_significant_GO_terms = .N),
  by = .(short_tier = tier)
]

immune_summary <- immune_results[
  class == "sweep" & tier %chin% tier_map$short_tier,
  .(
    short_tier = tier,
    immune_Fisher_OR = fisher_odds_ratio,
    immune_Fisher_p_value = fisher_p_value
  )
]

immune_summary_balancing <- immune_results[
  class == "balancing" & tier %chin% tier_map$short_tier,
  .(
    short_tier = tier,
    immune_Fisher_OR_balancing = fisher_odds_ratio,
    immune_Fisher_p_value_balancing = fisher_p_value
  )
]

barrier_rows <- rbindlist(lapply(seq_len(nrow(tier_map)), function(i) {
  f <- file.path(tier_map$barrier_dir[[i]], "barrier_validation_models.csv")
  if (!file.exists(f)) return(NULL)
  x <- as.data.table(fread(f))
  x[, short_tier := tier_map$short_tier[[i]]]
  x
}), use.names = TRUE, fill = TRUE)

barrier_models_summary <- barrier_rows[
  analysis %chin% c(
    "abs_delta_p ~ B_no_delta_p + K_local_mean",
    "f3 ~ B_no_delta_p + K_local_mean"
  ) & model_type == "lm"
][
  ,
  .(
    n = max(n, na.rm = TRUE),
    r_squared = max(r_squared, na.rm = TRUE),
    adj_r_squared = max(adj_r_squared, na.rm = TRUE),
    aic = max(aic, na.rm = TRUE)
  ),
  by = .(short_tier, analysis)
]

barrier_summary_wide <- merge(
  barrier_models_summary[analysis == "abs_delta_p ~ B_no_delta_p + K_local_mean", .(
    short_tier,
    abs_delta_p_n = n,
    abs_delta_p_R2 = r_squared,
    abs_delta_p_adj_R2 = adj_r_squared,
    abs_delta_p_AIC = aic
  )],
  barrier_models_summary[analysis == "f3 ~ B_no_delta_p + K_local_mean", .(
    short_tier,
    f3_n = n,
    f3_R2 = r_squared,
    f3_adj_R2 = adj_r_squared,
    f3_AIC = aic
  )],
  by = "short_tier",
  all = TRUE
)

load(label_rdata)
if (!exists("sd3_new_af2")) {
  stop("Expected sd3_new_af2 in ", label_rdata, call. = FALSE)
}

bgs_abs_delta_col <- if ("abs_delta_p" %in% names(sd3_new_af2)) "abs_delta_p" else "mean_abs_delta_p"

bgs_dt <- as.data.table(sd3_new_af2)[
  ,
  .(
    Chromosome_Names = as.character(Chromosome_Names),
    window_start = as.integer(.start_bp),
    window_end = as.integer(.end_bp),
    recombination = as.numeric(recombination),
    pi_average = as.numeric(pi_average),
    avg_wc_fst = as.numeric(avg_wc_fst),
    avg_dxy = as.numeric(avg_dxy),
    TajimaD = as.numeric(TajimaD),
    abs_delta_p = as.numeric(get(bgs_abs_delta_col)),
    f3 = as.numeric(f3)
  )
]

bgs_dt <- bgs_dt[
  is.finite(recombination) &
    is.finite(pi_average) &
    is.finite(avg_wc_fst) &
    is.finite(avg_dxy)
]

eps_rec <- max(min(bgs_dt$recombination[bgs_dt$recombination > 0], na.rm = TRUE) * 0.5, 1e-12)
eps_pi <- max(min(bgs_dt$pi_average[bgs_dt$pi_average > 0], na.rm = TRUE) * 0.5, 1e-12)
eps_dxy <- max(min(bgs_dt$avg_dxy[bgs_dt$avg_dxy > 0], na.rm = TRUE) * 0.5, 1e-12)

bgs_dt[, `:=`(
  log_rec = log10(recombination + eps_rec),
  log_pi = log10(pi_average + eps_pi),
  log_dxy = log10(avg_dxy + eps_dxy)
)]

m_pi2 <- lm(log_pi ~ poly(log_rec, 2, raw = TRUE), data = bgs_dt)
bgs_dt[, pi_resid_from_rec := log_pi - predict(m_pi2, newdata = bgs_dt)]

q_pi_low <- as.numeric(quantile(bgs_dt$pi_resid_from_rec, 0.10, na.rm = TRUE))
q_fst_hi <- as.numeric(quantile(bgs_dt$avg_wc_fst, 0.90, na.rm = TRUE))
q_dxy_low <- as.numeric(quantile(bgs_dt$avg_dxy, 0.10, na.rm = TRUE))
q_taj_low <- as.numeric(quantile(bgs_dt$TajimaD, 0.10, na.rm = TRUE))

bgs_dt[, `:=`(
  low_pi_resid = pi_resid_from_rec <= q_pi_low,
  high_fst = avg_wc_fst >= q_fst_hi,
  low_dxy = avg_dxy <= q_dxy_low,
  low_taj = fifelse(is.finite(TajimaD), TajimaD <= q_taj_low, FALSE)
)]

bgs_dt[, BGS_class := fifelse(
  low_pi_resid & (high_fst | low_taj),
  "sweep_like_excess",
  fifelse(low_pi_resid & !high_fst & low_dxy, "BGS_like", "other")
)]

setkey(bgs_dt, Chromosome_Names, window_start)

sweep_windows <- selected_windows[selection_class == "Geographic sweep"]
bgs_by_tier <- rbindlist(lapply(tier_map$short_tier, function(st) {
  sw <- sweep_windows[tier == st]
  if (!nrow(sw)) {
    return(data.table(
      short_tier = st,
      n_sweep_windows = 0L,
      BGS_like = 0L,
      sweep_like_excess = 0L,
      other = 0L,
      sweep_windows_beyond_BGS = 0L,
      fraction_sweep_windows_beyond_BGS = NA_real_,
      mean_recombination = NA_real_,
      mean_pi_average = NA_real_,
      mean_avg_wc_fst = NA_real_,
      mean_abs_delta_p = NA_real_,
      mean_f3 = NA_real_
    ))
  }
  sw <- unique(sw[, .(Chromosome_Names = as.character(Chromosome_Names), window_start = as.integer(window_start))])
  joined <- bgs_dt[sw, on = .(Chromosome_Names, window_start), nomatch = 0L]
  data.table(
    short_tier = st,
    n_sweep_windows = nrow(joined),
    BGS_like = sum(joined$BGS_class == "BGS_like", na.rm = TRUE),
    sweep_like_excess = sum(joined$BGS_class == "sweep_like_excess", na.rm = TRUE),
    other = sum(joined$BGS_class == "other", na.rm = TRUE),
    sweep_windows_beyond_BGS = sum(joined$BGS_class == "sweep_like_excess", na.rm = TRUE),
    fraction_sweep_windows_beyond_BGS = mean(joined$BGS_class == "sweep_like_excess", na.rm = TRUE),
    mean_recombination = mean(joined$recombination, na.rm = TRUE),
    mean_pi_average = mean(joined$pi_average, na.rm = TRUE),
    mean_avg_wc_fst = mean(joined$avg_wc_fst, na.rm = TRUE),
    mean_abs_delta_p = mean(joined$abs_delta_p, na.rm = TRUE),
    mean_f3 = mean(joined$f3, na.rm = TRUE)
  )
}), use.names = TRUE, fill = TRUE)

decision_enriched <- merge(decision, tier_map, by = "tier", all.x = TRUE)
decision_enriched <- merge(
  decision_enriched,
  gene_summary[class == "Balancing selection", .(short_tier, balancing_genes = unique_gene_ids)],
  by = "short_tier",
  all.x = TRUE
)
decision_enriched <- merge(
  decision_enriched,
  gene_summary[class == "Geographic sweep", .(short_tier, sweep_genes = unique_gene_ids)],
  by = "short_tier",
  all.x = TRUE
)
decision_enriched <- merge(decision_enriched, bgs_by_tier[, .(short_tier, BGS_like, sweep_like_excess)], by = "short_tier", all.x = TRUE)
decision_enriched <- merge(decision_enriched, barrier_summary_wide[, .(short_tier, abs_delta_p_R2, f3_R2)], by = "short_tier", all.x = TRUE)
decision_enriched <- merge(decision_enriched, go_summary, by = "short_tier", all.x = TRUE)
decision_enriched <- merge(decision_enriched, immune_summary, by = "short_tier", all.x = TRUE)
decision_enriched <- merge(decision_enriched, immune_summary_balancing, by = "short_tier", all.x = TRUE)

for (col in c("balancing_genes", "sweep_genes", "BGS_like", "sweep_like_excess", "number_of_significant_GO_terms")) {
  set(
    decision_enriched,
    i = which(is.na(decision_enriched[[col]])),
    j = col,
    value = 0
  )
}

decision_enriched <- decision_enriched[
  ,
  .(
    tier,
    neutral_windows,
    balancing_windows,
    sweep_windows,
    ambiguous_no_consensus_windows,
    balancing_genes = as.integer(balancing_genes),
    sweep_genes = as.integer(sweep_genes),
    BGS_like = as.integer(BGS_like),
    sweep_like_excess = as.integer(sweep_like_excess),
    abs_delta_p_R2,
    f3_R2,
    number_of_significant_GO_terms = as.integer(number_of_significant_GO_terms),
    immune_Fisher_OR,
    immune_Fisher_p_value,
    immune_Fisher_OR_balancing,
    immune_Fisher_p_value_balancing,
    recommended_paper_use
  )
]

revised_barrier_models_by_tier <- merge(
  tier_map[, .(tier, short_tier)],
  barrier_summary_wide,
  by = "short_tier",
  all.x = TRUE
)[, .(
  tier,
  abs_delta_p_n,
  abs_delta_p_R2,
  abs_delta_p_adj_R2,
  abs_delta_p_AIC,
  f3_n,
  f3_R2,
  f3_adj_R2,
  f3_AIC
)]

revised_go_by_tier <- merge(
  tier_map[, .(tier, short_tier)],
  go_summary,
  by = "short_tier",
  all.x = TRUE
)[, .(
  tier,
  number_of_significant_GO_terms = fifelse(is.na(number_of_significant_GO_terms), 0L, as.integer(number_of_significant_GO_terms))
)]

revised_immune_fisher_by_tier <- merge(
  tier_map[, .(tier, short_tier)],
  immune_summary,
  by = "short_tier",
  all.x = TRUE
)
revised_immune_fisher_by_tier <- merge(
  revised_immune_fisher_by_tier,
  immune_summary_balancing,
  by = "short_tier",
  all.x = TRUE
)

selection_bgs_out <- file.path(selection_output_dir, "revised_BGS_by_tier.tsv")
selection_barrier_out <- file.path(selection_output_dir, "revised_barrier_models_by_tier.tsv")
selection_go_out <- file.path(selection_output_dir, "revised_GO_by_tier.tsv")
selection_immune_out <- file.path(selection_output_dir, "revised_immune_fisher_by_tier.tsv")

barrier_bgs_out <- file.path(barrier_outputs_root, "revised_BGS_by_tier.tsv")
barrier_barrier_out <- file.path(barrier_outputs_root, "revised_barrier_models_by_tier.tsv")
barrier_go_out <- file.path(barrier_outputs_root, "revised_GO_by_tier.tsv")
barrier_immune_out <- file.path(barrier_outputs_root, "revised_immune_fisher_by_tier.tsv")
barrier_decision_out <- file.path(barrier_outputs_root, "revised_model_decision_table.tsv")

fwrite(bgs_by_tier, selection_bgs_out, sep = "\t")
fwrite(revised_barrier_models_by_tier, selection_barrier_out, sep = "\t")
fwrite(revised_go_by_tier, selection_go_out, sep = "\t")
fwrite(revised_immune_fisher_by_tier, selection_immune_out, sep = "\t")
fwrite(decision_enriched, selection_decision_path, sep = "\t")

fwrite(bgs_by_tier, barrier_bgs_out, sep = "\t")
fwrite(revised_barrier_models_by_tier, barrier_barrier_out, sep = "\t")
fwrite(revised_go_by_tier, barrier_go_out, sep = "\t")
fwrite(revised_immune_fisher_by_tier, barrier_immune_out, sep = "\t")
fwrite(decision_enriched, barrier_decision_out, sep = "\t")

message("Updated consolidated decision table: ", selection_decision_path)
message("Copied downstream summaries to: ", barrier_outputs_root)
