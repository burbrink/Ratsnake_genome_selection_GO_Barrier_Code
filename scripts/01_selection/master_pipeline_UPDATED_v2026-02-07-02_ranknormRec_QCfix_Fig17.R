############################################################
## Speciation model testing 2 (Option 4b: recombination EXCLUDED from clustering/labeling; tested downstream)
## UPDATED MASTER PIPELINE v2026-02-06-01
##
## Major changes vs your pasted v2025-12-15-01:
##   (1) GMM clustering is now 5D: z_dxy, z_fst, z_pi, z_f3, z_taj
##       - recombination (z_rec) is NOT used for clustering (non-circular).
##   (2) Model assignment scores now explicitly use Tajima's D and f3, but NOT recombination:
##       - TajimaD: negative tail supports Sweep; positive tail supports Balancing; positive tail penalizes Allopatry.
##       - f3: your f3 = mean((pt-ps1)*(pt-ps2)); negative f3 supports admixture -> boosts DGF only.
##   (3) RF now includes z_f3 consistently (and still excludes recombination).
##   (4) A new module tests whether recombination predicts model labels (and proposes an expected ordering),
##       using nonparametric tests + (multi)logistic models with effect sizes.
##   (5) [14H] genome tick plot now uses Chromosome_Names (no i.Chromosome_Names dependency).
##
## Assumes input data frame: sd3_new_af (or set sd3t <- your object)
## Required columns (minimum):
##   avg_dxy, avg_wc_fst, pi_average, f3, TajimaD, recombination,
##   Chromosome_Names (e.g. "Chromosome_1"), Chromosome_Size,
##   window_start, window_end
## Optional/extra columns okay.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mclust)
  library(randomForest)
  library(cluster)
  library(pROC)
  library(tibble)
  library(purrr)
  library(parallel)
  library(data.table)
})

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

find_existing_path <- function(candidates, label) {
  hits <- unique(normalizePath(candidates[file.exists(candidates)], winslash = "/", mustWork = FALSE))
  if (length(hits) == 0) {
    stop(label, " not found. Checked:\n", paste(sprintf("  - %s", unique(candidates)), collapse = "\n"))
  }
  hits[[1]]
}

show_plot <- function(p) {
  if (interactive()) print(p)
  invisible(p)
}

script_dir <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
data_rdata <- Sys.getenv("PIPELINE_INPUT_RDATA", unset = "")
if (!nzchar(data_rdata)) {
  data_rdata <- find_existing_path(
    c(
      file.path(getwd(), "east-west_rats_fst_dxy_pi_recombination_da_GC_div_model_f3_tajD_mds_het_Fhet_freq_p_parentals.rdata"),
      file.path(script_dir, "east-west_rats_fst_dxy_pi_recombination_da_GC_div_model_f3_tajD_mds_het_Fhet_freq_p_parentals.rdata"),
      file.path(project_dir, "east-west_rats_fst_dxy_pi_recombination_da_GC_div_model_f3_tajD_mds_het_Fhet_freq_p_parentals.rdata")
    ),
    "Pipeline input RData"
  )
}
cat("Loading pipeline input from:", data_rdata, "\n")
load(data_rdata) -> xx
## cores for parallel loops
n_cores <- {
  env_cores <- suppressWarnings(as.integer(Sys.getenv("PIPELINE_N_CORES", unset = NA_character_)))
  if (is.finite(env_cores) && !is.na(env_cores) && env_cores >= 1L) {
    out <- env_cores
  } else {
    nc <- suppressWarnings(parallel::detectCores())
    if (!is.finite(nc) || nc < 1) nc <- 1
    out <- max(1L, as.integer(nc) - 1L)
  }
  if (!is.finite(out) || out < 1L) out <- 1L
  out
}
cat("Using", n_cores, "cores for parallel operations.\n")

## RF confidence threshold for strict neutral override
conf_thresh <- 0.95

## Optional fast mode for smoke testing (set env PIPELINE_FAST_TEST=1)
fast_mode <- identical(Sys.getenv("PIPELINE_FAST_TEST"), "1")
if (fast_mode) cat("FAST mode enabled: reduced bootstrap/permutation settings for quick validation.\n")

## Neutral centroid radius threshold (in 5D z-space)
neutral_r0 <- 1.4  # slightly tighter in 5D than your 6D 1.5, but tune as needed

## DGF centroid gating based on f3: allow DGF only if mean z_f3 is sufficiently negative
## (negative f3 is consistent with admixture / gene flow signal)
dgf_f3_gate <- -0.25  # tune as needed (e.g., -0.2 to -0.6)

############################################################
## Helpers
############################################################

## Rank-based inverse-normal transform: maps any continuous distribution to ~N(0,1)
ranknorm <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  n_ok <- sum(ok)
  if (n_ok < 2L) return(out)
  r <- rank(x[ok], ties.method = "average")
  u <- (r - 0.5) / n_ok
  out[ok] <- qnorm(u)
  out
}

## tail helpers for TajD and f3 logic
pos_tail <- function(z) pmax(0,  z)
neg_tail <- function(z) pmax(0, -z)

harmonize_continuous_model_label <- function(x) {
  dplyr::case_when(
    x %in% c("Neutral/other", "Neutral", "Neutral / other") ~ "Neutral / equilibrium",
    x %in% c("Balancing-like", "Balancing", "Balancing_like") ~ "Balancing selection",
    x %in% c("Sweep-like", "Sweep", "Sweep_like") ~ "Geographic sweep",
    x %in% c("DGF-like", "DGF", "DGF_like") ~ "Divergence with gene flow",
    x %in% c("Allopatry-like", "Allopatry", "Allopatry_like") ~ "Allopatry-like",
    TRUE ~ x
  )
}

## robust key column finders (used later in agreement/joining)
.pick_first <- function(x) if (length(x) == 0L) NA_character_ else x[[1]]
.find_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  .pick_first(hit)
}
.find_chr_col <- function(df) {
  candidates <- c("Chromosome_Names","Chromosomes","Chromosome","chrom","chr","scaffold","contig")
  .find_col(df, candidates)
}
.find_start_col <- function(df) {
  candidates <- c(".start_bp","window_start","start_bp","start","win_start","pos_start")
  .find_col(df, candidates)
}
.find_end_col <- function(df) {
  candidates <- c(".end_bp","window_end","end_bp","end","win_end","pos_end")
  .find_col(df, candidates)
}
.rename_keys_to_standard <- function(df, chr_col, start_col, end_col = NA_character_) {
  out <- df
  if (!is.na(chr_col)   && chr_col   != "Chromosome_Names") out <- dplyr::rename(out, Chromosome_Names = !!chr_col)
  if (!is.na(start_col) && start_col != ".start_bp")        out <- dplyr::rename(out, .start_bp = !!start_col)
  if (!is.na(end_col)   && end_col   != ".end_bp")          out <- dplyr::rename(out, .end_bp = !!end_col)
  out
}

############################################################
## [0] BUILD sd3t FROM sd3_new_af (RANK-NORMALIZED z METRICS)
############################################################

sd3t <- sd3_new_af

num_cols <- c("avg_dxy", "avg_wc_fst", "pi_average",
              "recombination", "f3", "TajimaD")

sd3t <- sd3t %>%
  mutate(across(all_of(num_cols), as.numeric))

## Genomic coordinates
if (!".start_bp" %in% names(sd3t) && "window_start" %in% names(sd3t)) {
  sd3t <- sd3t %>% mutate(.start_bp = as.numeric(window_start))
}
if (!".end_bp" %in% names(sd3t) && "window_end" %in% names(sd3t)) {
  sd3t <- sd3t %>% mutate(.end_bp = as.numeric(window_end))
}
if (!".mid_bp" %in% names(sd3t)) {
  sd3t <- sd3t %>% mutate(.mid_bp = (.start_bp + .end_bp) / 2)
}

## Rank-normalized ~N(0,1) “z” for all metrics
sd3t <- sd3t %>%
  mutate(
    z_dxy = ranknorm(avg_dxy),
    z_fst = ranknorm(avg_wc_fst),
    z_pi  = ranknorm(pi_average),
    z_f3  = ranknorm(f3),
    z_taj = ranknorm(TajimaD),
    z_rec = ranknorm(recombination)
  )

cat("\n========== [0] sd3t constructed (rank-normalized z for all metrics) ==========\n")
cat("Rows:", nrow(sd3t), "\n")
cat("Missing counts (z_*):\n")
print(colSums(!is.finite(as.matrix(sd3t %>% dplyr::select(z_dxy,z_fst,z_pi,z_f3,z_taj,z_rec)))))

############################################################
## [1] 5D MATRIX FOR GMM (EXCLUDES RECOMBINATION)
##     IMPORTANT: Xz already standardized-ish (ranknorm), so NO scale()
############################################################

valid <- with(sd3t,
              is.finite(z_dxy) & is.finite(z_fst) &
                is.finite(z_pi)  & is.finite(z_f3)  &
                is.finite(z_taj))

Xz <- cbind(
  z_dxy = sd3t$z_dxy[valid],
  z_fst = sd3t$z_fst[valid],
  z_pi  = sd3t$z_pi [valid],
  z_f3  = sd3t$z_f3 [valid],
  z_taj = sd3t$z_taj[valid]
)

############################################################
## [1A] G-SELECTION BY STABILITY (VEV)
############################################################

cat("\n================ [1A] G-SELECTION BY STABILITY (VEV) ==================\n")

model_name <- "VEV"
n_boot_G   <- if (fast_mode) 5 else 50
n_sil_G    <- min(5000L, nrow(Xz))

## [1A.1] G = 1: BIC baseline
cat("\n--- Evaluating G = 1 (BIC only) ---\n")
gmm_G1 <- tryCatch(
  Mclust(Xz, G = 1, modelNames = model_name),
  error = function(e) {
    cat("  [FAIL] Mclust G=1 failed:", conditionMessage(e), "\n")
    NULL
  }
)
BIC_G1 <- if (!is.null(gmm_G1)) gmm_G1$bic else NA_real_

bic_G1_df <- tibble(
  G        = 1L,
  BIC      = BIC_G1,
  sil_mean = NA_real_,
  ari_mean = NA_real_,
  ari_sd   = NA_real_
)

## [1A.2] G = 2–9: silhouette + bootstrap ARI
G_candidates <- if (fast_mode) 2:4 else 2:10

eval_one_G <- function(G) {
  cat("\n--- Evaluating G =", G, " ---\n")
  set.seed(100 + G)
  gmm_try <- tryCatch(
    Mclust(Xz, G = G, modelNames = model_name),
    error = function(e) {
      cat("  [FAIL] Mclust failed for G =", G, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(gmm_try)) {
    return(tibble(G = G, BIC = NA_real_, sil_mean = NA_real_, ari_mean = NA_real_, ari_sd = NA_real_))
  }

  labs    <- gmm_try$classification
  BIC_val <- gmm_try$bic

  if (length(unique(labs)) < 2L) {
    cat("  [WARN] All points in one effective cluster for G =", G, "\n")
    return(tibble(G = G, BIC = BIC_val, sil_mean = NA_real_, ari_mean = NA_real_, ari_sd = NA_real_))
  }

  ## Subsampled silhouette
  set.seed(200 + G)
  idx_sil <- sample(seq_len(nrow(Xz)), n_sil_G)
  X_sub   <- Xz[idx_sil, ]
  lab_sub <- labs[idx_sil]
  sil     <- silhouette(lab_sub, dist(X_sub))
  sil_mean <- summary(sil)$avg.width
  cat("  Mean silhouette:", round(sil_mean, 4), "\n")

  ## Bootstrap ARI stability
  cat("  Bootstrap ARI (", n_boot_G, " reps)...\n", sep = "")
  ari_vals <- unlist(
    mclapply(seq_len(n_boot_G), function(b) {
      set.seed(300 + G*1000 + b)
      out <- tryCatch({
        idx_boot  <- sample(seq_len(nrow(Xz)), replace = TRUE)
        X_boot    <- Xz[idx_boot, ]
        gmm_boot  <- Mclust(X_boot, G = G, modelNames = model_name)
        labs_boot <- gmm_boot$classification
        adjustedRandIndex(labs[idx_boot], labs_boot)
      }, error = function(e) NA_real_)
      out
    }, mc.cores = n_cores)
  )

  ari_vals <- ari_vals[is.finite(ari_vals)]
  ari_mean <- if (length(ari_vals) > 0) mean(ari_vals) else NA_real_
  ari_sd   <- if (length(ari_vals) > 0) sd(ari_vals)   else NA_real_
  cat("  ARI mean:", round(ari_mean,4), " SD:", round(ari_sd,4), "\n")

  tibble(G = G, BIC = BIC_val, sil_mean = sil_mean, ari_mean = ari_mean, ari_sd = ari_sd)
}

stab_list_2plus <- lapply(G_candidates, eval_one_G)
stab_summary_2plus <- bind_rows(stab_list_2plus) %>% arrange(G)
stab_summary <- bind_rows(bic_G1_df, stab_summary_2plus) %>% arrange(G)

cat("\n[1A] Stability summary over G (including G=1 BIC baseline):\n")
print(as.data.frame(stab_summary))
data.table::fwrite(as.data.table(stab_summary), "GMM_Gselection_StabilitySummary.tsv", sep = "\t")

## [1A.3] Elbow plots
p_G_elbow_sil_ari <- stab_summary %>%
  filter(G >= 2) %>%
  pivot_longer(cols = c("sil_mean","ari_mean"),
               names_to = "metric",
               values_to = "value") %>%
  ggplot(aes(x = G, y = value, color = metric)) +
  geom_line() + geom_point() +
  theme_bw() +
  labs(title = "G-selection: silhouette & ARI vs G (VEV, G ≥ 2)",
       x = "Number of components (G)", y = "Score", color = "Metric")
show_plot(p_G_elbow_sil_ari)
ggsave("GMM_Gselection_Silhouette_ARI_vs_G.pdf", p_G_elbow_sil_ari, width = 7.4, height = 5.4)
ggsave("GMM_Gselection_Silhouette_ARI_vs_G.png", p_G_elbow_sil_ari, width = 7.4, height = 5.4, dpi = 300)

p_G_elbow_BIC <- stab_summary %>%
  ggplot(aes(x = G, y = BIC)) +
  geom_line() + geom_point() +
  theme_bw() +
  labs(title = "BIC vs G (including G=1 baseline)",
       x = "Number of components (G)", y = "BIC")
show_plot(p_G_elbow_BIC)
ggsave("GMM_Gselection_BIC_vs_G.pdf", p_G_elbow_BIC, width = 7.4, height = 5.4)
ggsave("GMM_Gselection_BIC_vs_G.png", p_G_elbow_BIC, width = 7.4, height = 5.4, dpi = 300)

## [1A.4] Choose bestG among G ≥ 3
ari_floor <- 0.60
sil_floor <- 0.10

cat("\n[1A.4] G-selection rule (for bestG ≥ 3):\n")
cat("  - Consider G = 3..10\n")
cat("  - Keep those with ari_mean ≥", ari_floor, "and sil_mean ≥", sil_floor, "\n")
cat("  - Choose largest G among passers; else max ARI among G≥3; else fallback to 2\n\n")

stab_3plus <- stab_summary_2plus %>% filter(G >= 3)
stab_stable_3plus <- stab_3plus %>%
  filter(!is.na(ari_mean), ari_mean >= ari_floor, sil_mean >= sil_floor)

if (nrow(stab_stable_3plus) > 0) {
  stab_best <- stab_stable_3plus %>% arrange(desc(G), desc(ari_mean), desc(sil_mean)) %>% slice(1)
  cat("[1A.4] Stable G candidates (G≥3):\n")
  print(as.data.frame(stab_stable_3plus))
} else if (nrow(stab_3plus %>% filter(!is.na(ari_mean))) > 0) {
  stab_best <- stab_3plus %>% filter(!is.na(ari_mean)) %>% arrange(desc(ari_mean), desc(sil_mean)) %>% slice(1)
  cat("[1A.4] No G≥3 passed floors; falling back to max-ARI G≥3.\n")
} else {
  cat("[1A.4] No usable G ≥ 3 found; bestG will be set to 2.\n")
  stab_best <- tibble(G = 2, BIC = NA_real_, sil_mean = NA_real_, ari_mean = NA_real_, ari_sd = NA_real_)
}

bestG <- stab_best$G[1]
cat("\n[1A] Chosen bestG (prefer G≥3, else 2):", bestG, "\n")
cat("Row for chosen bestG:\n")
print(as.data.frame(stab_best))

############################################################
## Helper: SAVE RF DIAGNOSTICS
############################################################

save_rf_diagnostics <- function(rf, level_tag, out_dir = ".") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  file_stub <- paste0("rf_", level_tag)

  err_rate <- rf$err.rate
  if (is.null(err_rate)) {
    err_rate <- matrix(NA_real_, nrow = 1L, ncol = 1L,
                       dimnames = list(NULL, "OOB"))
  }
  err_rate_df <- as.data.frame(err_rate)
  err_rate_df$tree <- seq_len(nrow(err_rate_df))
  err_rate_df <- err_rate_df %>% select(tree, everything())

  final_err <- err_rate_df[nrow(err_rate_df), , drop = FALSE]
  oob_col <- if ("OOB" %in% names(final_err)) "OOB" else names(final_err)[2]

  rf_summary <- tibble(
    level_tag = level_tag,
    ntree = rf$ntree,
    final_oob_error = suppressWarnings(as.numeric(final_err[[oob_col]]))
  )

  class_cols <- setdiff(names(final_err), c("tree", "OOB"))
  if (length(class_cols) > 0) {
    class_df <- final_err %>%
      select(all_of(class_cols)) %>%
      pivot_longer(cols = everything(),
                   names_to = "class",
                   values_to = "class_error")
    rf_summary$class_specific_error <- paste0(class_df$class, "=", signif(class_df$class_error, 4), collapse = "; ")
  } else {
    rf_summary$class_specific_error <- NA_character_
  }

  imp_raw <- importance(rf)
  imp_df <- as.data.frame(imp_raw)
  imp_df$variable <- rownames(imp_df)
  rownames(imp_df) <- NULL

  if (!"MeanDecreaseAccuracy" %in% names(imp_df)) imp_df$MeanDecreaseAccuracy <- NA_real_
  if (!"MeanDecreaseGini" %in% names(imp_df)) imp_df$MeanDecreaseGini <- NA_real_

  imp_cols <- c("variable", "MeanDecreaseAccuracy", "MeanDecreaseGini",
                setdiff(names(imp_df), c("variable", "MeanDecreaseAccuracy", "MeanDecreaseGini")))
  imp_df <- imp_df[, imp_cols[imp_cols %in% names(imp_df)], drop = FALSE]

  rank_metric <- if (all(is.na(imp_df$MeanDecreaseAccuracy))) "MeanDecreaseGini" else "MeanDecreaseAccuracy"
  imp_top6 <- imp_df %>%
    arrange(desc(.data[[rank_metric]]), desc(MeanDecreaseGini), variable) %>%
    slice_head(n = 6)

  cat("\n[", level_tag, " RF] Final OOB error: ", signif(rf_summary$final_oob_error[1], 5), "\n", sep = "")
  cat("[", level_tag, " RF] Top 6 variables (", rank_metric, "):\n", sep = "")
  print(as.data.frame(imp_top6[, c("variable", "MeanDecreaseAccuracy", "MeanDecreaseGini")]))

  imp_plot_df <- imp_df %>%
    select(variable, MeanDecreaseAccuracy, MeanDecreaseGini) %>%
    pivot_longer(cols = c("MeanDecreaseAccuracy", "MeanDecreaseGini"),
                 names_to = "metric",
                 values_to = "importance") %>%
    filter(is.finite(importance))

  plot_pdf <- file.path(out_dir, paste0(file_stub, "_variable_importance_plot.pdf"))
  plot_png <- file.path(out_dir, paste0(file_stub, "_variable_importance_plot.png"))

  if (nrow(imp_plot_df) > 0) {
    metric_order <- imp_plot_df %>%
      filter(metric == rank_metric) %>%
      arrange(importance) %>%
      pull(variable)
    if (length(metric_order) == 0) {
      metric_order <- imp_plot_df %>%
        group_by(variable) %>%
        summarise(max_importance = max(importance, na.rm = TRUE), .groups = "drop") %>%
        arrange(max_importance) %>%
        pull(variable)
    }
    imp_plot_df$variable <- factor(imp_plot_df$variable, levels = metric_order)

    p_imp <- ggplot(imp_plot_df, aes(x = variable, y = importance, fill = metric)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.68, color = "grey25") +
      coord_flip() +
      scale_fill_manual(values = c("MeanDecreaseAccuracy" = "#2C7FB8",
                                   "MeanDecreaseGini" = "#D95F0E")) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "top",
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold")
      ) +
      labs(
        title = paste0("Random Forest variable importance (", level_tag, ")"),
        x = NULL,
        y = "Importance",
        fill = NULL
      )

    ggsave(plot_pdf, p_imp, width = 7.5, height = 4.8)
    ggsave(plot_png, p_imp, width = 7.5, height = 4.8, dpi = 300)
  }

  write.csv(rf_summary,
            file.path(out_dir, paste0(file_stub, "_oob_error_summary.csv")),
            row.names = FALSE)
  write.csv(err_rate_df,
            file.path(out_dir, paste0(file_stub, "_error_by_tree.csv")),
            row.names = FALSE)
  write.csv(imp_df,
            file.path(out_dir, paste0(file_stub, "_variable_importance.csv")),
            row.names = FALSE)

  cat("[", level_tag, " RF] Wrote files:\n", sep = "")
  cat("  ", file.path(out_dir, paste0(file_stub, "_oob_error_summary.csv")), "\n", sep = "")
  cat("  ", file.path(out_dir, paste0(file_stub, "_error_by_tree.csv")), "\n", sep = "")
  cat("  ", file.path(out_dir, paste0(file_stub, "_variable_importance.csv")), "\n", sep = "")
  if (nrow(imp_plot_df) > 0) {
    cat("  ", plot_pdf, "\n", sep = "")
    cat("  ", plot_png, "\n", sep = "")
  }

  invisible(list(
    summary = rf_summary,
    error_by_tree = err_rate_df,
    variable_importance = imp_df
  ))
}

############################################################
## Helper: RUN FULL PIPELINE FOR A GIVEN G
############################################################

run_level_pipeline <- function(level_tag, G_target) {

  # robustly display model_name if defined earlier
  model_name_local <- if (exists("model_name", inherits = TRUE)) get("model_name", inherits = TRUE) else "(unknown)"

  cat("
============================================================
")
  cat("RUNNING FULL PIPELINE FOR LEVEL: ", level_tag, "
", sep="")
  cat("GMM model: ", model_name_local, "   G: ", G_target, "
", sep="")
  cat("============================================================
")

  cat("\n\n")

  sd3_clusters <- sd3t

  ##########################################################
  ## [1B] FINAL GMM FIT (VEV, G_target) in 5D
  ##########################################################
  set.seed(123 + G_target)
  gmm <- Mclust(Xz, G = G_target, modelNames = model_name)

  cat("\n================ [", level_tag, " 1B] FINAL GMM FIT ==================\n", sep = "")
  cat("Model:", gmm$modelName, "  G:", gmm$G, "\n")
  cat("BIC:", gmm$bic, "\n\n")

  labs   <- gmm$classification
  post   <- gmm$z
  p_max  <- apply(post, 1, max)

  cluster_full   <- rep(NA_integer_, nrow(sd3_clusters))
  p_cluster_full <- rep(NA_real_,    nrow(sd3_clusters))
  cluster_full[valid]   <- labs
  p_cluster_full[valid] <- p_max

  sd3_clusters <- sd3_clusters %>%
    mutate(
      cluster   = factor(cluster_full),
      p_cluster = p_cluster_full
    )

  ##########################################################
  ## [1C] GMM CENTROIDS IN 5D FIT SPACE
  ##########################################################
  mu_mat <- t(gmm$parameters$mean)
  colnames(mu_mat) <- colnames(Xz)
  centroids_fit5 <- as.data.frame(mu_mat)
  centroids_fit5$cluster <- factor(seq_len(nrow(centroids_fit5)))

  cat("\n================ [", level_tag, " 1C] GMM CENTROIDS (fit 5D z-space) ==================\n", sep = "")
  print(as.data.frame(centroids_fit5))

  ##########################################################
  ## [1D] GMM SILHOUETTE (5D)
  ##########################################################
  cat("\n================ [", level_tag, " 1D] GMM CLUSTER SILHOUETTE ==================\n", sep = "")

  labs_gmm <- as.integer(sd3_clusters$cluster[valid])

  set.seed(42 + G_target)
  n_sil <- min(5000L, nrow(Xz))
  sub_idx <- sample(seq_len(nrow(Xz)), n_sil)

  Xz_sub   <- Xz[sub_idx, ]
  labs_sub <- labs_gmm[sub_idx]

  sil <- silhouette(labs_sub, dist(Xz_sub))
  sil_sum <- summary(sil)

  cat("Mean silhouette (GMM, ", level_tag, "): ", sil_sum$avg.width, "\n", sep = "")
  cat("Silhouette width by cluster:\n")
  print(sil_sum$clus.avg.widths)

  gmm_sil_full <- rep(NA_real_, nrow(sd3_clusters))
  gmm_sil_full[valid][sub_idx] <- sil[, "sil_width"]

  sd3_clusters <- sd3_clusters %>%
    mutate(gmm_silhouette = gmm_sil_full)

  ##########################################################
  ## [1E] PERMUTATION NULL FOR GMM SILHOUETTE (5D)
  ##########################################################
  cat("\n================ [", level_tag, " 1E] PERMUTATION NULL: GMM SILHOUETTE ==================\n", sep = "")

  set.seed(2025 + G_target)
  n_perm <- if (fast_mode) 20 else 100

  sil_null <- unlist(
    mclapply(seq_len(n_perm), function(b) {
      X_perm <- apply(Xz, 2, sample, replace = FALSE)
      out <- tryCatch({
        gmm_perm  <- Mclust(X_perm, G = G_target, modelNames = gmm$modelName)
        labs_perm <- gmm_perm$classification
        n_sub2    <- min(3000L, nrow(X_perm))
        idx_sub2  <- sample(seq_len(nrow(X_perm)), n_sub2)
        sil_perm  <- silhouette(labs_perm[idx_sub2], dist(X_perm[idx_sub2, ]))
        summary(sil_perm)$avg.width
      }, error = function(e) NA_real_)
      out
    }, mc.cores = n_cores)
  )

  cat("\nObserved GMM mean silhouette (", level_tag, "): ", sil_sum$avg.width, "\n", sep = "")
  cat("Null silhouette summary (non-NA):\n")
  print(summary(sil_null[is.finite(sil_null)]))

  p_sil_perm <- ggplot(
    data.frame(
      type  = c(rep("null", sum(is.finite(sil_null))), "observed"),
      value = c(sil_null[is.finite(sil_null)], sil_sum$avg.width)
    ),
    aes(x = value, fill = type)
  ) +
    geom_histogram(data = ~ subset(., type == "null"), bins = 20, alpha = 0.6) +
    geom_vline(xintercept = sil_sum$avg.width, linetype = "dashed", linewidth = 1) +
    theme_bw() +
    labs(title = paste("Permutation null for GMM silhouette (", level_tag, ")", sep = ""),
         x = "Mean silhouette", y = "Count")
  show_plot(p_sil_perm)

  ##########################################################
  ## [1F] MAHALANOBIS DISTANCES & CENTROID SEPARATION (5D)
  ##########################################################
  cat("\n================ [", level_tag, " 1F] MAHALANOBIS + CENTROID DISTANCES ==================\n", sep = "")

  mu_list    <- lapply(seq_len(G_target), function(k) gmm$parameters$mean[, k])
  Sigma_list <- lapply(seq_len(G_target), function(k) gmm$parameters$variance$sigma[, , k])

  mahal_gmm <- function(x, k) {
    mu <- mu_list[[k]]
    S  <- Sigma_list[[k]]
    md <- mahalanobis(x, center = mu, cov = S)
    sqrt(md)
  }

  cluster_idx_full <- as.integer(sd3_clusters$cluster[valid])
  dists <- numeric(length(cluster_idx_full))
  for (k in seq_len(G_target)) {
    idx_k <- which(cluster_idx_full == k)
    if (length(idx_k) > 0) dists[idx_k] <- mahal_gmm(Xz[idx_k, ], k)
  }

  mahal_full <- rep(NA_real_, nrow(sd3_clusters))
  mahal_full[valid] <- dists

  sd3_clusters <- sd3_clusters %>%
    mutate(mahal_cluster = mahal_full)

  mahal_summary <- sd3_clusters %>%
    filter(!is.na(mahal_cluster), !is.na(cluster)) %>%
    group_by(cluster) %>%
    summarise(
      n       = n(),
      mean_md = mean(mahal_cluster),
      sd_md   = sd(mahal_cluster),
      q50_md  = median(mahal_cluster),
      q90_md  = quantile(mahal_cluster, 0.9),
      q99_md  = quantile(mahal_cluster, 0.99),
      .groups = "drop"
    ) %>%
    arrange(as.integer(cluster))

  cat("\nMahalanobis distance summary per cluster:\n")
  print(as.data.frame(mahal_summary))

  cent_mat  <- as.matrix(centroids_fit5[, colnames(Xz)])
  dist_cent <- as.matrix(dist(cent_mat))
  cat("\nBetween-centroid Euclidean distances (5D z-space):\n")
  print(round(dist_cent, 3))

  ##########################################################
  ## [2A] z_* metrics are already rank-normalized
  ##########################################################
  cat("\n================ [", level_tag, " 2A] z_* METRICS (ALREADY RANK-NORMALIZED) ==================\n", sep = "")
  cat("Using columns: z_dxy, z_fst, z_pi, z_f3, z_taj (no further scaling)\n")

  ##########################################################
  ## [2B] EMPIRICAL CENTROIDS (means per cluster) in 5D
  ##########################################################
  cat("\n================ [", level_tag, " 2B] GMM CENTROIDS IN z-SPACE (EMPIRICAL) ==================\n", sep = "")

  centroids_z5 <- sd3_clusters %>%
    filter(!is.na(cluster)) %>%
    group_by(cluster) %>%
    summarise(
      z_dxy = mean(z_dxy, na.rm = TRUE),
      z_fst = mean(z_fst, na.rm = TRUE),
      z_pi  = mean(z_pi,  na.rm = TRUE),
      z_f3  = mean(z_f3,  na.rm = TRUE),
      z_taj = mean(z_taj, na.rm = TRUE),
      .groups = "drop"
    )

  print(as.data.frame(centroids_z5))

  ##########################################################
  ## [2C] SCORE-BASED MODEL LABELS (USES TajD + f3; EXCLUDES recombination)
  ##########################################################
  cat("\n================ [", level_tag, " 2C] SCORE-BASED MODEL LABELS FOR CENTROIDS ==================\n", sep = "")

  ## f3: your definition => negative f3 supports admixture, so use negative tail
  centroids_z5 <- centroids_z5 %>%
    mutate(
      admix = neg_tail(z_f3),  # only negative tail counts (admixture-like)
      sweep_score_c = ( +1.2*z_fst -1.0*z_dxy -1.2*z_pi -1.0*neg_tail(z_taj) -0.3*admix ),
      bal_score_c   = ( -1.0*z_fst +1.3*z_pi +1.2*pos_tail(z_taj) -0.3*abs(z_dxy) ),
      dgf_score_c   = ( +1.1*z_fst +1.1*z_dxy +0.8*admix -0.4*z_pi ),
      allop_score_c = ( +1.0*z_fst -1.0*z_pi -0.8*abs(z_dxy) -0.6*pos_tail(z_taj) -0.3*admix ),
      r_centroid    = sqrt(z_dxy^2 + z_fst^2 + z_pi^2 + z_f3^2 + z_taj^2)
    )

  assign_model_from_centroid <- function(sw, ba, dg, al, f3,
                                       neutral_thresh = 0.5,
                                       dgf_gate = dgf_f3_gate) {

    scores <- c(
      "Geographic sweep"           = sw,
      "Balancing selection"        = ba,
      "Divergence with gene flow"  = dg,
      "Allopatry-like"             = al
    )

    # Ensure names are plain character (prevents weird NSE/backtick issues)
    names(scores) <- as.character(names(scores))

    ## Gate DGF: only allow if centroid mean z_f3 is sufficiently negative
    ## (negative f3 is consistent with admixture / gene flow signal)
    if (!is.finite(f3) || f3 > dgf_gate) {
      scores["Divergence with gene flow"] <- -Inf
    }


    best_name <- names(which.max(scores))
    best_val  <- max(scores)
    if (best_val < neutral_thresh) "Neutral / equilibrium" else best_name
  }

  
.tag <- if (exists("level_tag")) level_tag else "MAIN"
cat("[", .tag, " 2C] DGF gating: require centroid mean z_f3 <= ", dgf_f3_gate, "
", sep = "")
cat("[", .tag, " 2C] Centroids eligible for DGF: ",
    sum(is.finite(centroids_z5$z_f3) & centroids_z5$z_f3 <= dgf_f3_gate),
    " / ", nrow(centroids_z5), "
", sep = "")
centroids_z5$model <- mapply(
    assign_model_from_centroid,
    sw = centroids_z5$sweep_score_c,
    ba = centroids_z5$bal_score_c,
    dg = centroids_z5$dgf_score_c,
    al = centroids_z5$allop_score_c,
    f3 = centroids_z5$z_f3
  )

  cat("\n[", level_tag, " 2C+] Neutral centroid rule: exactly ONE centroid forced to 'Neutral / equilibrium'.\n", sep = "")
  cat("    Using radius threshold neutral_r0 =", neutral_r0, "\n")

  neutral_candidates <- which(centroids_z5$r_centroid < neutral_r0)
  if (length(neutral_candidates) >= 1) {
    k_neutral <- neutral_candidates[which.min(centroids_z5$r_centroid[neutral_candidates])]
  } else {
    k_neutral <- which.min(centroids_z5$r_centroid)
  }

  centroids_z5$model <- as.character(centroids_z5$model)
  centroids_z5$model[k_neutral] <- "Neutral / equilibrium"

  print(as.data.frame(centroids_z5))

  ## Map centroid → model for each window
  cluster_to_model <- setNames(centroids_z5$model, centroids_z5$cluster)

  sd3_clusters <- sd3_clusters %>%
    mutate(
      gmm_model_raw = cluster_to_model[as.character(cluster)],
      gmm_model_raw = ifelse(is.na(gmm_model_raw), "Neutral / equilibrium", gmm_model_raw),
      gmm_model_raw = as.character(gmm_model_raw)
    )

  cat("\n[", level_tag, " 2D] Table of GMM model labels (gmm_model_raw):\n", sep = "")
  model_summary_gmm <- sd3_clusters %>%
    dplyr::count(gmm_model_raw, name = "n_windows") %>%
    arrange(desc(n_windows))
  print(as.data.frame(model_summary_gmm))

  cat("\n[", level_tag, " 2E] CHECKING DIRECTIONS PER GMM MODEL (means of z_*) ==================\n", sep = "")
  gmm_direction_check <- sd3_clusters %>%
    group_by(gmm_model_raw) %>%
    summarise(
      n          = n(),
      mean_z_dxy = mean(z_dxy, na.rm = TRUE),
      mean_z_fst = mean(z_fst, na.rm = TRUE),
      mean_z_pi  = mean(z_pi,  na.rm = TRUE),
      mean_z_f3  = mean(z_f3,  na.rm = TRUE),
      mean_z_taj = mean(z_taj, na.rm = TRUE),
      mean_z_rec = mean(z_rec, na.rm = TRUE),
      .groups = "drop"
    )
  print(as.data.frame(gmm_direction_check))

  ##########################################################
  ## [2F] ACTIVE MODEL SET
  ##########################################################
  models_found_gmm  <- sort(unique(centroids_z5$model))
  models_nonneutral <- setdiff(models_found_gmm, "Neutral / equilibrium")
  active_models <- c(models_nonneutral, "Neutral / equilibrium")

  cat("\n[", level_tag, " 2F] Models found by GMM (centroids):\n", sep = "")
  print(models_found_gmm)
  cat("\n[", level_tag, " 2F] Active models downstream (non-neutral + Neutral):\n", sep = "")
  print(active_models)

  ##########################################################
  ## [3] RANDOM FOREST (uses z_*; includes z_f3; excludes recombination)
  ##########################################################
  palette_models <- c(
    "Geographic sweep"          = "#E24A33",
    "Balancing selection"       = "#1B7837",
    "Divergence with gene flow" = "gold",
    "Allopatry-like"            = "#984EA3",
    "Neutral / equilibrium"     = "grey80"
  )

  cat("\n================ [", level_tag, " 3A] RF: CORE TRAINING SET (ACTIVE MODELS ONLY) ==================\n", sep = "")

  core <- sd3_clusters %>%
    filter(!is.na(gmm_model_raw),
           gmm_model_raw %in% active_models,
           p_cluster >= 0.9) %>%
    mutate(model = factor(gmm_model_raw, levels = active_models)) %>%
    select(model, z_dxy, z_fst, z_pi, z_f3, z_taj) %>%
    filter(is.finite(z_dxy), is.finite(z_fst), is.finite(z_pi), is.finite(z_f3), is.finite(z_taj)) %>%
    mutate(model = droplevels(model))

  cat("Core training set size (", level_tag, "): ", nrow(core), " windows\n", sep = "")
  cat("Training classes present:\n")
  print(table(core$model))

  rf_available <- (length(unique(core$model)) >= 2)

  if (!rf_available) {
    cat("\n[", level_tag, " 3B] RF skipped: only one class present in core.\n", sep = "")
    sd3_clusters <- sd3_clusters %>%
      mutate(
        rf_model      = factor(NA_character_, levels = active_models),
        rf_pmax       = NA_real_,
        rf_silhouette = NA_real_
      )
  } else {
    cat("\n================ [", level_tag, " 3B] RF: FIT ==================\n", sep = "")
    set.seed(777 + G_target)

    rf <- randomForest(
      model ~ z_dxy + z_fst + z_pi + z_f3 + z_taj,
      data = core,
      ntree = 500,
      importance = TRUE
    )

    print(rf)
    cat("\nOOB error rate (", level_tag, "): ", rf$err.rate[rf$ntree, "OOB"], "\n", sep = "")

    cat("\n================ [", level_tag, " 3C] RF CONFUSION MATRIX (CORE) ==================\n", sep = "")
    pred_core <- predict(rf, core)
    print(table(predicted = pred_core, actual = core$model))

    cat("\n================ [", level_tag, " 3D] RF VARIABLE IMPORTANCE ==================\n", sep = "")
    print(importance(rf))

    cat("\n================ [", level_tag, " 3D+] RF DIAGNOSTICS EXPORT ==================\n", sep = "")
    rf_diag <- save_rf_diagnostics(rf = rf, level_tag = level_tag, out_dir = ".")

    cat("\n================ [", level_tag, " 3E] RF PER-CLASS AUC (CORE) ==================\n", sep = "")
    prob_core <- predict(rf, core, type = "prob")
    for (mod in colnames(prob_core)) {
      auc_val <- auc(core$model == mod, prob_core[, mod])
      cat(mod, " : AUC =", as.numeric(auc_val), "\n")
    }

    ## [3F] RF prediction genome-wide + strict neutral override
    cat("\n================ [", level_tag, " 3F] RF: PREDICT ALL WINDOWS ==================\n", sep = "")

    pred_dat <- sd3_clusters %>%
      select(z_dxy, z_fst, z_pi, z_f3, z_taj) %>%
      mutate(across(everything(), as.numeric))

    valid_pred <- with(pred_dat,
                       is.finite(z_dxy) & is.finite(z_fst) &
                         is.finite(z_pi)  & is.finite(z_f3)  &
                         is.finite(z_taj))

    rf_model_raw <- rep(NA_character_, nrow(pred_dat))
    rf_pmax      <- rep(NA_real_,      nrow(pred_dat))

    prob_all <- predict(rf, newdata = pred_dat[valid_pred, ], type = "prob")
    rf_model_raw[valid_pred] <- colnames(prob_all)[max.col(prob_all)]
    rf_pmax[valid_pred]      <- apply(prob_all, 1, max)

    rf_model_strict <- rf_model_raw
    rf_model_strict[is.na(rf_model_strict)] <- "Neutral / equilibrium"
    rf_model_strict[rf_pmax < conf_thresh]  <- "Neutral / equilibrium"

    sd3_clusters <- sd3_clusters %>%
      mutate(
        rf_model = factor(rf_model_strict, levels = active_models),
        rf_pmax  = rf_pmax
      )

    rf_model_props <- sd3_clusters %>%
      dplyr::count(rf_model, name = "n_windows") %>%
      mutate(total = sum(n_windows, na.rm = TRUE),
             prop  = n_windows / total) %>%
      arrange(desc(prop))

    cat("\n================ [", level_tag, " 3G] RF MODEL PROPORTIONS (STRICT) ==================\n", sep = "")
    print(as.data.frame(rf_model_props))

    ## [3H] RF silhouette
    cat("\n================ [", level_tag, " 3H] RF SILHOUETTE (SUBSAMPLED) ==================\n", sep = "")

    rf_pca_raw <- sd3_clusters %>%
      select(z_dxy, z_fst, z_pi, z_f3, z_taj, rf_model)

    idx_rf <- with(rf_pca_raw,
                   is.finite(z_dxy) & is.finite(z_fst) &
                     is.finite(z_pi)  & is.finite(z_f3)  &
                     is.finite(z_taj) & !is.na(rf_model))

    X_rf    <- as.matrix(rf_pca_raw[idx_rf, c("z_dxy","z_fst","z_pi","z_f3","z_taj")])
    labs_rf <- as.integer(factor(rf_pca_raw$rf_model[idx_rf], levels = active_models))

    set.seed(101 + G_target)
    n_sil_rf <- min(5000L, nrow(X_rf))
    sub_idx_rf <- sample(seq_len(nrow(X_rf)), n_sil_rf)

    X_rf_sub    <- X_rf[sub_idx_rf, ]
    labs_rf_sub <- labs_rf[sub_idx_rf]

    sil_rf <- silhouette(labs_rf_sub, dist(X_rf_sub))
    sil_rf_sum <- summary(sil_rf)

    cat("Mean silhouette (RF labels, ", level_tag, "): ", sil_rf_sum$avg.width, "\n", sep = "")
    cat("Silhouette width by RF model:\n")
    print(sil_rf_sum$clus.avg.widths)

    rf_silhouette_full <- rep(NA_real_, nrow(sd3_clusters))
    rf_valid_idx_full  <- which(idx_rf)[sub_idx_rf]
    rf_silhouette_full[rf_valid_idx_full] <- sil_rf[, "sil_width"]

    sd3_clusters <- sd3_clusters %>%
      mutate(rf_silhouette = rf_silhouette_full)
  }

  ##########################################################
  ## [4] PCA OF z-SCORES (GMM vs RF) in 5D
  ##########################################################
  cat("\n================ [", level_tag, " 4] PCA OF z-SCORES (GMM vs RF) ==================\n", sep = "")

  gmm_levels <- sort(unique(sd3_clusters$gmm_model_raw))

  pca_raw <- sd3_clusters %>%
    select(z_dxy, z_fst, z_pi, z_f3, z_taj,
           gmm_model = gmm_model_raw,
           rf_model)

  idx_pca <- with(pca_raw,
                  is.finite(z_dxy) & is.finite(z_fst) &
                    is.finite(z_pi)  & is.finite(z_f3)  &
                    is.finite(z_taj))

  X_res   <- as.matrix(pca_raw[idx_pca, c("z_dxy","z_fst","z_pi","z_f3","z_taj")])
  pca_res <- prcomp(X_res, center = TRUE, scale. = FALSE)

  scores <- as.data.frame(pca_res$x[, 1:3])
  scores$gmm_model <- factor(pca_raw$gmm_model[idx_pca], levels = gmm_levels)
  scores$rf_model  <- pca_raw$rf_model[idx_pca]

  cat("PCA (", level_tag, ") done on ", nrow(scores), " windows.\n", sep = "")

  set.seed(1234 + G_target)
  n_plot <- min(30000, nrow(scores))
  sub_idx <- sample(seq_len(nrow(scores)), size = n_plot, replace = FALSE)
  scores_sub <- scores[sub_idx, ]

  p_pca_gmm <- ggplot(scores_sub, aes(x = PC1, y = PC2, color = gmm_model)) +
    geom_point(alpha = 0.25, size = 0.7) +
    stat_ellipse(aes(group = gmm_model), type = "norm", level = 0.68,
                 linewidth = 0.5, alpha = 0.8, na.rm = TRUE) +
    scale_color_manual(values = palette_models[gmm_levels], drop = FALSE) +
    theme_bw() +
    labs(title = paste("PCA PC1–PC2 colored by GMM model (", level_tag, ")", sep = ""),
         x = "PC1", y = "PC2", color = "GMM model")
  show_plot(p_pca_gmm)
  ggsave(paste0("PCA_GMM_modelspace_", level_tag, ".pdf"), p_pca_gmm, width = 7.2, height = 5.8)
  ggsave(paste0("PCA_GMM_modelspace_", level_tag, ".png"), p_pca_gmm, width = 7.2, height = 5.8, dpi = 300)

  if (rf_available) {
    p_pca_rf <- ggplot(scores_sub, aes(x = PC1, y = PC2, color = rf_model)) +
      geom_point(alpha = 0.25, size = 0.7) +
      stat_ellipse(aes(group = rf_model), type = "norm", level = 0.68,
                   linewidth = 0.5, alpha = 0.8, na.rm = TRUE) +
      scale_color_manual(values = palette_models[levels(sd3_clusters$rf_model)], drop = FALSE) +
      theme_bw() +
      labs(title = paste("PCA PC1–PC2 colored by RF model (", level_tag, ")", sep = ""),
           x = "PC1", y = "PC2", color = "RF model")
    show_plot(p_pca_rf)
    ggsave(paste0("PCA_RF_modelspace_", level_tag, ".pdf"), p_pca_rf, width = 7.2, height = 5.8)
    ggsave(paste0("PCA_RF_modelspace_", level_tag, ".png"), p_pca_rf, width = 7.2, height = 5.8, dpi = 300)
  }

  ##########################################################
  ## [4B] ENHANCED G=4 VISUALIZATION (3D + MULTI-PANEL)
  ##########################################################
  if (identical(level_tag, "G4")) {
    cat("\n================ [G4 4B] ENHANCED PCA/UMAP VISUALIZATION ==================\n")

    pca_pairs <- dplyr::bind_rows(
      scores_sub %>% transmute(x = PC1, y = PC2, gmm_model, panel = "PC1 vs PC2"),
      scores_sub %>% transmute(x = PC1, y = PC3, gmm_model, panel = "PC1 vs PC3"),
      scores_sub %>% transmute(x = PC2, y = PC3, gmm_model, panel = "PC2 vs PC3")
    )

    p_pca_tripanel <- ggplot(pca_pairs, aes(x = x, y = y, color = gmm_model)) +
      geom_point(alpha = 0.22, size = 0.6) +
      stat_ellipse(aes(group = gmm_model), type = "norm", level = 0.68,
                   linewidth = 0.45, alpha = 0.85, na.rm = TRUE) +
      facet_wrap(~ panel, scales = "free") +
      scale_color_manual(values = palette_models[gmm_levels], drop = FALSE) +
      theme_bw() +
      labs(title = "G4 PCA model-space (all 3 PC pair views)",
           x = "Principal component", y = "Principal component", color = "GMM model")
    show_plot(p_pca_tripanel)
    ggsave("PCA_GMM_modelspace_G4_tripanel.pdf", p_pca_tripanel, width = 11, height = 4.2)
    ggsave("PCA_GMM_modelspace_G4_tripanel.png", p_pca_tripanel, width = 11, height = 4.2, dpi = 300)

    if (requireNamespace("plotly", quietly = TRUE) &&
        requireNamespace("htmlwidgets", quietly = TRUE)) {
      p3 <- plotly::plot_ly(
        data = scores_sub,
        x = ~PC1, y = ~PC2, z = ~PC3,
        color = ~gmm_model,
        colors = unname(palette_models[gmm_levels]),
        type = "scatter3d", mode = "markers",
        marker = list(size = 2.2, opacity = 0.38),
        text = ~paste0("Model: ", gmm_model),
        hoverinfo = "text"
      )

      cent3 <- scores_sub %>%
        group_by(gmm_model) %>%
        summarise(PC1 = mean(PC1, na.rm = TRUE),
                  PC2 = mean(PC2, na.rm = TRUE),
                  PC3 = mean(PC3, na.rm = TRUE),
                  .groups = "drop")

      p3 <- plotly::add_trace(
        p3, data = cent3, inherit = FALSE,
        x = ~PC1, y = ~PC2, z = ~PC3,
        type = "scatter3d", mode = "markers+text",
        marker = list(size = 6, symbol = "diamond", color = "black"),
        text = ~as.character(gmm_model), textposition = "top center",
        hoverinfo = "text", showlegend = FALSE
      ) %>%
        plotly::layout(
          title = "G4 PCA 3D model-space (interactive)",
          scene = list(xaxis = list(title = "PC1"),
                       yaxis = list(title = "PC2"),
                       zaxis = list(title = "PC3"))
        )

      tryCatch(
        htmlwidgets::saveWidget(
          plotly::as_widget(p3),
          file = "PCA_GMM_modelspace_G4_3D_interactive.html",
          selfcontained = TRUE
        ),
        error = function(e) {
          htmlwidgets::saveWidget(
            plotly::as_widget(p3),
            file = "PCA_GMM_modelspace_G4_3D_interactive.html",
            selfcontained = FALSE
          )
        }
      )
      cat("[G4 4B] Wrote: PCA_GMM_modelspace_G4_3D_interactive.html\n")
    } else {
      cat("[G4 4B] SKIP 3D PCA HTML: requires packages 'plotly' + 'htmlwidgets'.\n")
    }

    if (requireNamespace("uwot", quietly = TRUE)) {
      set.seed(4404)
      um3 <- uwot::umap(
        X_res,
        n_components = 3,
        n_neighbors = 30,
        min_dist = 0.15,
        metric = "euclidean",
        scale = FALSE,
        verbose = FALSE
      )

      um_df <- data.frame(
        UMAP1 = um3[, 1], UMAP2 = um3[, 2], UMAP3 = um3[, 3],
        gmm_model = factor(pca_raw$gmm_model[idx_pca], levels = gmm_levels)
      )

      set.seed(4405)
      um_n <- min(30000, nrow(um_df))
      um_sub <- um_df[sample(seq_len(nrow(um_df)), size = um_n), , drop = FALSE]

      p_um2 <- ggplot(um_sub, aes(x = UMAP1, y = UMAP2, color = gmm_model)) +
        geom_point(alpha = 0.25, size = 0.7) +
        stat_ellipse(aes(group = gmm_model), type = "norm", level = 0.68,
                     linewidth = 0.5, alpha = 0.8, na.rm = TRUE) +
        scale_color_manual(values = palette_models[gmm_levels], drop = FALSE) +
        theme_bw() +
        labs(title = "G4 UMAP 2D colored by GMM model", color = "GMM model")
      show_plot(p_um2)
      ggsave("UMAP_GMM_modelspace_G4_2D.pdf", p_um2, width = 7.2, height = 5.8)
      ggsave("UMAP_GMM_modelspace_G4_2D.png", p_um2, width = 7.2, height = 5.8, dpi = 300)

      if (requireNamespace("plotly", quietly = TRUE) &&
          requireNamespace("htmlwidgets", quietly = TRUE)) {
        p_um3 <- plotly::plot_ly(
          data = um_sub,
          x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
          color = ~gmm_model,
          colors = unname(palette_models[gmm_levels]),
          type = "scatter3d", mode = "markers",
          marker = list(size = 2.2, opacity = 0.38),
          text = ~paste0("Model: ", gmm_model),
          hoverinfo = "text"
        ) %>%
          plotly::layout(
            title = "G4 UMAP 3D model-space (interactive)",
            scene = list(xaxis = list(title = "UMAP1"),
                         yaxis = list(title = "UMAP2"),
                         zaxis = list(title = "UMAP3"))
          )
        tryCatch(
          htmlwidgets::saveWidget(
            plotly::as_widget(p_um3),
            file = "UMAP_GMM_modelspace_G4_3D_interactive.html",
            selfcontained = TRUE
          ),
          error = function(e) {
            htmlwidgets::saveWidget(
              plotly::as_widget(p_um3),
              file = "UMAP_GMM_modelspace_G4_3D_interactive.html",
              selfcontained = FALSE
            )
          }
        )
        cat("[G4 4B] Wrote: UMAP_GMM_modelspace_G4_3D_interactive.html\n")
      }
    } else {
      cat("[G4 4B] SKIP UMAP: package 'uwot' not installed.\n")
    }

    if (requireNamespace("GGally", quietly = TRUE)) {
      p_pairs <- GGally::ggpairs(
        scores_sub %>% select(PC1, PC2, PC3, gmm_model),
        columns = 1:3,
        mapping = ggplot2::aes(color = gmm_model, alpha = 0.2),
        progress = FALSE
      )
      grDevices::pdf("PCA_GMM_modelspace_G4_pairwise_matrix.pdf", width = 9.5, height = 9)
      show_plot(p_pairs)
      grDevices::dev.off()
      cat("[G4 4B] Wrote: PCA_GMM_modelspace_G4_pairwise_matrix.pdf\n")
    } else {
      cat("[G4 4B] SKIP pairwise matrix: package 'GGally' not installed.\n")
    }
  }

  ##########################################################
  ## [6] CONTINUOUS MODEL SCORES (USES TajD + f3; EXCLUDES recombination)
  ##########################################################
  cat("\n================ [", level_tag, " 6] CONTINUOUS MODEL SCORES ==================\n", sep = "")

  sd3_cont <- sd3_clusters %>%
    mutate(
      admix = neg_tail(z_f3),
      sweep_score = ( +1.2*z_fst -1.0*z_dxy -1.2*z_pi -1.0*neg_tail(z_taj) -0.3*admix ),
      bal_score   = ( -1.0*z_fst +1.3*z_pi +1.2*pos_tail(z_taj) -0.3*abs(z_dxy) ),
      dgf_score   = ( +1.1*z_fst +1.1*z_dxy +0.8*admix -0.4*z_pi ),
      allop_score = ( +1.0*z_fst -1.0*z_pi -0.8*abs(z_dxy) -0.6*pos_tail(z_taj) -0.3*admix ),
      mid_Mb = .mid_bp / 1e6
    )

  feature_cols_6d <- c("z_dxy", "z_fst", "z_pi", "z_f3", "z_taj", "z_rec")
  selection_classes <- c("Neutral / equilibrium", "Balancing selection", "Geographic sweep")

  centroid_from_subset <- function(df, label_name) {
    out <- rep(NA_real_, length(feature_cols_6d))
    names(out) <- feature_cols_6d
    if (nrow(df) > 0) {
      mat <- as.matrix(df[, ..feature_cols_6d])
      storage.mode(mat) <- "double"
      ok <- rowSums(!is.finite(mat)) == 0L
      if (any(ok)) out <- colMeans(mat[ok, , drop = FALSE], na.rm = TRUE)
    }
    tibble::tibble(
      centroid_type = label_name,
      z_dxy = out[["z_dxy"]],
      z_fst = out[["z_fst"]],
      z_pi = out[["z_pi"]],
      z_f3 = out[["z_f3"]],
      z_taj = out[["z_taj"]],
      z_rec = out[["z_rec"]]
    )
  }

  euclid_to_centroid <- function(df, centroid_vec) {
    mat <- as.matrix(df[, ..feature_cols_6d])
    storage.mode(mat) <- "double"
    out <- rep(NA_real_, nrow(mat))
    centroid_vec <- as.numeric(centroid_vec[feature_cols_6d])
    if (length(centroid_vec) != length(feature_cols_6d) || any(!is.finite(centroid_vec))) return(out)
    ok <- rowSums(!is.finite(mat)) == 0L
    if (!any(ok)) return(out)
    diffs <- sweep(mat[ok, , drop = FALSE], 2, centroid_vec, FUN = "-")
    out[ok] <- sqrt(rowSums(diffs^2))
    out
  }

  nearest_distance_label <- function(dist_mat, label_names) {
    n <- nrow(dist_mat)
    best_label <- rep(NA_character_, n)
    margin <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      vals <- as.numeric(dist_mat[i, ])
      if (!all(is.finite(vals))) next
      ord <- order(vals, decreasing = FALSE)
      best_label[i] <- label_names[ord[1]]
      margin[i] <- if (length(ord) >= 2L) vals[ord[2]] - vals[ord[1]] else NA_real_
    }
    list(label = best_label, margin = margin)
  }

  assign_margin_label <- function(bal_score, sweep_score, threshold) {
    dplyr::case_when(
      !is.finite(bal_score) | !is.finite(sweep_score) ~ NA_character_,
      abs(bal_score - sweep_score) < threshold ~ "Ambiguous",
      bal_score > sweep_score ~ "Balancing selection",
      sweep_score > bal_score ~ "Geographic sweep",
      TRUE ~ "Ambiguous"
    )
  }

  p_hist_scores <- sd3_cont %>%
    select(sweep_score, bal_score, dgf_score, allop_score) %>%
    pivot_longer(everything(), names_to = "score_type", values_to = "value") %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 60) +
    facet_wrap(~ score_type, scales = "free") +
    theme_bw() +
    labs(title = paste("Distributions of continuous model scores (", level_tag, ")", sep = ""),
         x = "Score", y = "Count")
  show_plot(p_hist_scores)

  ##########################################################
  ## [7] TOP-N CONTINUOUS MODEL CLASSES
  ##########################################################
  cat("\n================ [", level_tag, " 7] TOP-N CONTINUOUS MODEL CLASSES ==================\n", sep = "")

  n_top <- 500

  scores_sweep <- sd3_cont$sweep_score
  scores_bal   <- sd3_cont$bal_score
  scores_dgf   <- sd3_cont$dgf_score
  scores_allop <- sd3_cont$allop_score

  rank_sweep <- rank_bal <- rank_dgf <- rank_allop <- rep(NA_integer_, nrow(sd3_cont))

  idx_sweep <- which(is.finite(scores_sweep))
  idx_bal   <- which(is.finite(scores_bal))
  idx_dgf   <- which(is.finite(scores_dgf))
  idx_allop <- which(is.finite(scores_allop))

  rank_sweep[idx_sweep] <- rank(-scores_sweep[idx_sweep], ties.method = "first")
  rank_bal[idx_bal]     <- rank(-scores_bal[idx_bal],     ties.method = "first")
  rank_dgf[idx_dgf]     <- rank(-scores_dgf[idx_dgf],     ties.method = "first")
  rank_allop[idx_allop] <- rank(-scores_allop[idx_allop], ties.method = "first")

  allow_sweep <- "Geographic sweep"          %in% models_nonneutral
  allow_bal   <- "Balancing selection"       %in% models_nonneutral
  allow_dgf   <- "Divergence with gene flow" %in% models_nonneutral
  allow_allop <- "Allopatry-like"            %in% models_nonneutral

  sd3_cont <- sd3_cont %>%
    mutate(
      rank_sweep = rank_sweep,
      rank_bal   = rank_bal,
      rank_dgf   = rank_dgf,
      rank_allop = rank_allop,
      is_sweep_top_raw = !is.na(rank_sweep) & rank_sweep <= n_top,
      is_bal_top_raw   = !is.na(rank_bal)   & rank_bal   <= n_top,
      is_dgf_top_raw   = !is.na(rank_dgf)   & rank_dgf   <= n_top,
      is_allop_top_raw = !is.na(rank_allop) & rank_allop <= n_top,
      is_sweep_top = is_sweep_top_raw & allow_sweep,
      is_bal_top   = is_bal_top_raw   & allow_bal,
      is_dgf_top   = is_dgf_top_raw   & allow_dgf,
      is_allop_top = is_allop_top_raw & allow_allop,
      model_cont = dplyr::case_when(
        is_sweep_top  ~ "Sweep-like",
        is_bal_top    ~ "Balancing-like",
        is_dgf_top    ~ "DGF-like",
        is_allop_top  ~ "Allopatry-like",
        TRUE          ~ "Neutral/other"
      ),
      model_cont_tail500 = model_cont
    )

  score_mat <- as.matrix(sd3_cont[, c("sweep_score", "bal_score", "dgf_score", "allop_score")])
  colnames(score_mat) <- c("Sweep-like", "Balancing-like", "DGF-like", "Allopatry-like")
  allowed_cont <- c(
    "Sweep-like" = allow_sweep,
    "Balancing-like" = allow_bal,
    "DGF-like" = allow_dgf,
    "Allopatry-like" = allow_allop
  )
  score_mat_allowed <- score_mat
  score_mat_allowed[, !allowed_cont[colnames(score_mat_allowed)]] <- NA_real_

  best_cont_label <- rep(NA_character_, nrow(score_mat_allowed))
  best_cont_score <- rep(NA_real_, nrow(score_mat_allowed))
  second_cont_score <- rep(NA_real_, nrow(score_mat_allowed))
  for (i in seq_len(nrow(score_mat_allowed))) {
    vals <- score_mat_allowed[i, ]
    finite_vals <- is.finite(vals)
    if (any(finite_vals)) {
      ord <- order(vals, decreasing = TRUE, na.last = NA)
      best_cont_label[i] <- names(vals)[ord[1]]
      best_cont_score[i] <- vals[ord[1]]
      second_cont_score[i] <- if (length(ord) >= 2L) vals[ord[2]] else NA_real_
    }
  }

  score_quantile <- function(x) {
    ok <- is.finite(x)
    out <- rep(NA_real_, length(x))
    if (sum(ok) > 0L) out[ok] <- rank(x[ok], ties.method = "average") / sum(ok)
    out
  }
  q_sweep <- score_quantile(sd3_cont$sweep_score)
  q_bal <- score_quantile(sd3_cont$bal_score)
  q_dgf <- score_quantile(sd3_cont$dgf_score)
  q_allop <- score_quantile(sd3_cont$allop_score)
  best_cont_score_quantile <- dplyr::case_when(
    best_cont_label == "Sweep-like" ~ q_sweep,
    best_cont_label == "Balancing-like" ~ q_bal,
    best_cont_label == "DGF-like" ~ q_dgf,
    best_cont_label == "Allopatry-like" ~ q_allop,
    TRUE ~ NA_real_
  )

  sd3_cont <- sd3_cont %>%
    mutate(
      best_cont_label = best_cont_label,
      best_cont_score = best_cont_score,
      second_cont_score = second_cont_score,
      cont_margin = best_cont_score - second_cont_score,
      cont_score_quantile = best_cont_score_quantile,
      cont_top_90 = is.finite(cont_score_quantile) & cont_score_quantile >= 0.90,
      cont_top_95 = is.finite(cont_score_quantile) & cont_score_quantile >= 0.95,
      cont_top_99 = is.finite(cont_score_quantile) & cont_score_quantile >= 0.99,
      model_cont_top1pct = if_else(cont_top_99, best_cont_label, "Neutral/other"),
      model_cont_top5pct = if_else(cont_top_95, best_cont_label, "Neutral/other"),
      model_cont_top10pct = if_else(cont_top_90, best_cont_label, "Neutral/other")
    )

  gmm_rf_seed <- sd3_cont %>%
    filter(!is.na(gmm_model_raw), !is.na(rf_model), gmm_model_raw == rf_model)

  centroid_neutral_origin <- setNames(rep(0, length(feature_cols_6d)), feature_cols_6d)
  centroid_neutral_empirical <- centroid_from_subset(
    gmm_rf_seed %>% filter(gmm_model_raw == "Neutral / equilibrium"),
    "neutral_empirical"
  )
  centroid_balancing <- centroid_from_subset(
    gmm_rf_seed %>% filter(gmm_model_raw == "Balancing selection"),
    "balancing"
  )
  centroid_sweep <- centroid_from_subset(
    gmm_rf_seed %>% filter(gmm_model_raw == "Geographic sweep"),
    "sweep"
  )

  continuous_distance_centroids <- dplyr::bind_rows(
    tibble::tibble(
      centroid_type = "neutral_origin",
      z_dxy = 0, z_fst = 0, z_pi = 0, z_f3 = 0, z_taj = 0, z_rec = 0
    ),
    centroid_neutral_empirical,
    centroid_balancing,
    centroid_sweep
  )

  centroid_vec_neutral_empirical <- stats::setNames(
    as.numeric(centroid_neutral_empirical[1, feature_cols_6d]),
    feature_cols_6d
  )
  centroid_vec_balancing <- stats::setNames(
    as.numeric(centroid_balancing[1, feature_cols_6d]),
    feature_cols_6d
  )
  centroid_vec_sweep <- stats::setNames(
    as.numeric(centroid_sweep[1, feature_cols_6d]),
    feature_cols_6d
  )

  d_neutral_origin <- euclid_to_centroid(sd3_cont, centroid_neutral_origin)
  d_neutral_empirical <- euclid_to_centroid(sd3_cont, centroid_vec_neutral_empirical)
  d_balancing <- euclid_to_centroid(sd3_cont, centroid_vec_balancing)
  d_sweep <- euclid_to_centroid(sd3_cont, centroid_vec_sweep)

  nearest_origin <- nearest_distance_label(
    cbind(d_neutral_origin, d_balancing, d_sweep),
    c("Neutral / equilibrium", "Balancing selection", "Geographic sweep")
  )
  nearest_empirical <- nearest_distance_label(
    cbind(d_neutral_empirical, d_balancing, d_sweep),
    c("Neutral / equilibrium", "Balancing selection", "Geographic sweep")
  )

  continuous_score_margin <- abs(sd3_cont$bal_score - sd3_cont$sweep_score)

  sd3_cont <- sd3_cont %>%
    mutate(
      d_neutral_origin = d_neutral_origin,
      d_neutral_empirical = d_neutral_empirical,
      d_balancing = d_balancing,
      d_sweep = d_sweep,
      continuous_distance_label_origin = nearest_origin$label,
      continuous_distance_label_empirical = nearest_empirical$label,
      continuous_distance_margin_origin = nearest_origin$margin,
      continuous_distance_margin_empirical = nearest_empirical$margin,
      continuous_score_margin = continuous_score_margin,
      continuous_margin_label_025 = assign_margin_label(bal_score, sweep_score, 0.25),
      continuous_margin_label_050 = assign_margin_label(bal_score, sweep_score, 0.50),
      continuous_margin_label_100 = assign_margin_label(bal_score, sweep_score, 1.00),
      continuous_margin_label = continuous_margin_label_050
    )

  cat("[", level_tag, " 7] Top-N classification (N =", n_top, ") respecting active models:\n", sep = "")
  print(table(sd3_cont$model_cont, useNA = "ifany"))
  cat("[", level_tag, " 7] Best-supported continuous label classification (all windows):\n", sep = "")
  print(table(sd3_cont$best_cont_label, useNA = "ifany"))
  cat("[", level_tag, " 7] Continuous distance label classification (empirical neutral centroid):\n", sep = "")
  print(table(sd3_cont$continuous_distance_label_empirical, useNA = "ifany"))
  cat("[", level_tag, " 7] Continuous margin summary:\n", sep = "")
  print(summary(sd3_cont$cont_margin))
  cat("[", level_tag, " 7] Continuous distance-margin summary (empirical neutral centroid):\n", sep = "")
  print(summary(sd3_cont$continuous_distance_margin_empirical))

  ##########################################################
  ## [8] GENOME TRACK (CONTINUOUS TOP-N) - uses Chromosome_Names
  ##########################################################
  cat("\n================ [", level_tag, " 8] GENOME TRACK (CONTINUOUS TOP-N) ==================\n", sep = "")

  chrom_sizes_cont <- sd3_cont %>%
    distinct(Chromosome_Names, Chromosome_Size) %>%
    mutate(
      Chromosome_Names = as.character(Chromosome_Names),
      is_ZW = grepl("ZW", Chromosome_Names, ignore.case = TRUE),
      chr_index = if_else(is_ZW, 99L, suppressWarnings(as.integer(gsub("Chromosome_", "", Chromosome_Names))))
    ) %>%
    arrange(chr_index) %>%
    mutate(
      Chromosome_Names = factor(Chromosome_Names, levels = Chromosome_Names),
      chr_len_bp = Chromosome_Size
    )

  final_levels_cont <- levels(chrom_sizes_cont$Chromosome_Names)

  sd3_cont <- sd3_cont %>%
    mutate(Chromosome_Names = factor(as.character(Chromosome_Names), levels = final_levels_cont))

  chrom_bg_cont <- chrom_sizes_cont %>% mutate(y_base = seq_len(n()))
  lane_half_bg <- 0.2

  lane_rects_cont <- chrom_bg_cont %>%
    transmute(
      Chromosome_Names,
      xmin = 0,
      xmax = chr_len_bp / 1e6,
      ymin = y_base - lane_half_bg,
      ymax = y_base + lane_half_bg
    )

  palette_cont <- c(
    "Sweep-like"      = "#E24A33",
    "Balancing-like"  = "#1B7837",
    "DGF-like"        = "#348ABD",
    "Allopatry-like"  = "#984EA3"
  )

  model_cont_possible <- c("Sweep-like","Balancing-like","DGF-like","Allopatry-like")
  models_sel <- intersect(model_cont_possible, unique(sd3_cont$model_cont))

  sel_points <- sd3_cont %>%
    filter(model_cont %in% models_sel, is.finite(mid_Mb)) %>%
    inner_join(chrom_bg_cont %>% select(Chromosome_Names, y_base), by = "Chromosome_Names") %>%
    mutate(y = y_base + 0.45)

  p_genome_points <- ggplot() +
    geom_rect(data = lane_rects_cont,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey90", color = NA) +
    geom_point(data = sel_points,
               aes(x = mid_Mb, y = y, color = model_cont),
               size = 1.4, alpha = 0.9) +
    scale_y_reverse(
      breaks = chrom_bg_cont$y_base,
      labels = levels(chrom_bg_cont$Chromosome_Names),
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    scale_color_manual(values = palette_cont[models_sel], drop = FALSE) +
    labs(x = "Position (Mb)", y = "Chromosome",
         color = paste0("Top-", n_top, " windows"),
         title = paste("Continuous model windows (", level_tag, ")", sep = "")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(hjust = 1),
          legend.position = "bottom")
  show_plot(p_genome_points)
  ggsave(paste0("Chromosome_continuous_models_", level_tag, ".pdf"), p_genome_points, width = 12, height = 7)
  ggsave(paste0("Chromosome_continuous_models_", level_tag, ".png"), p_genome_points, width = 12, height = 7, dpi = 300)

  ##########################################################
  ## [9] RECOMB DESERTS (defined on rank-normalized recombination z_rec; preserves quantiles)
  ##########################################################
  cat("\n================ [", level_tag, " 9] OVERLAP WITH RECOMBINATION DESERTS ==================\n", sep = "")

  sd3_cont <- sd3_cont %>%
    mutate(
      rec_quantile = ntile(z_rec, 10),
      is_rec_desert = rec_quantile == 1
    )

  cat("Fraction of genome in recombination desert (lowest 10%, ", level_tag, "):\n", sep = "")
  print(mean(sd3_cont$is_rec_desert, na.rm = TRUE))

  overlap_tab <- sd3_cont %>%
    filter(model_cont %in% c("Sweep-like","Balancing-like","DGF-like","Allopatry-like","Neutral/other")) %>%
    group_by(model_cont) %>%
    summarise(
      n_total  = n(),
      n_desert = sum(is_rec_desert, na.rm = TRUE),
      frac_desert = n_desert / n_total,
      .groups = "drop"
    )
  cat("\nOverlap of continuous model classes with recombination deserts (", level_tag, "):\n", sep = "")
  print(as.data.frame(overlap_tab))

  ##########################################################
  ## [11] HEADER MODEL SUMMARY (GMM vs RF vs CONTINUOUS)
  ##########################################################
  cat("\n================ [", level_tag, " 11] HEADER MODEL SUMMARY ==================\n", sep = "")

  header_gmm <- sd3_clusters %>%
    dplyr::count(gmm_model_raw, name = "n_windows") %>%
    rename(model = gmm_model_raw) %>%
    mutate(method = "GMM")

  header_rf <- sd3_clusters %>%
    dplyr::count(rf_model, name = "n_windows") %>%
    rename(model = rf_model) %>%
    mutate(method = "RF (strict)")

  header_cont <- sd3_cont %>%
    dplyr::count(model_cont, name = "n_windows") %>%
    rename(model = model_cont) %>%
    mutate(method = "Continuous (top-N)")

  header_cont_best <- sd3_cont %>%
    dplyr::count(best_cont_label, name = "n_windows") %>%
    rename(model = best_cont_label) %>%
    mutate(method = "Continuous (best label)")

  header_cont_top1 <- sd3_cont %>%
    dplyr::count(model_cont_top1pct, name = "n_windows") %>%
    rename(model = model_cont_top1pct) %>%
    mutate(method = "Continuous (top 1%)")

  header_models <- bind_rows(header_gmm, header_rf, header_cont, header_cont_best, header_cont_top1) %>%
    group_by(method) %>%
    mutate(total = sum(n_windows), prop = n_windows / total) %>%
    ungroup() %>%
    arrange(method, desc(prop))

  print(as.data.frame(header_models))

  ##########################################################
  ## [12] AGREEMENT: GMM vs RF vs CONTINUOUS (kept; optional use later)
  ##########################################################
  cat("\n================ [", level_tag, " 12] AGREEMENT: GMM vs RF vs CONTINUOUS ==================\n", sep = "")

  agree_dat <- sd3_cont %>%
    select(Chromosome_Names, .start_bp, .mid_bp,
           z_dxy, z_fst, z_pi, z_f3, z_taj, z_rec,
           gmm_model = gmm_model_raw,
           sweep_score, bal_score,
           model_cont, model_cont_tail500, best_cont_label, best_cont_score,
           second_cont_score, cont_margin, cont_score_quantile,
           model_cont_top1pct, model_cont_top5pct, model_cont_top10pct,
           continuous_distance_label_origin, continuous_distance_label_empirical,
           continuous_distance_margin_origin, continuous_distance_margin_empirical,
           d_neutral_origin, d_neutral_empirical, d_balancing, d_sweep,
           continuous_score_margin,
           continuous_margin_label_025, continuous_margin_label_050,
           continuous_margin_label_100, continuous_margin_label) %>%
    left_join(
      sd3_clusters %>% select(Chromosome_Names, .start_bp, rf_model, rf_pmax, gmm_silhouette, rf_silhouette),
      by = c("Chromosome_Names", ".start_bp")
    )

  gmm_sil_thr <- quantile(agree_dat$gmm_silhouette, probs = 0.75, na.rm = TRUE)
  rf_sil_thr  <- quantile(agree_dat$rf_silhouette,  probs = 0.75, na.rm = TRUE)

  cat("  GMM silhouette threshold (75th pct):", gmm_sil_thr, "\n")
  cat("  RF  silhouette threshold (75th pct):", rf_sil_thr,  "\n")

  agree_dat <- agree_dat %>%
    mutate(
      high_gmm_sil = !is.na(gmm_silhouette) & gmm_silhouette >= gmm_sil_thr,
      high_rf_sil  = !is.na(rf_silhouette)  & rf_silhouette  >= rf_sil_thr,
      high_rf_p    = !is.na(rf_pmax)        & rf_pmax        >= conf_thresh,
      high_conf    = high_gmm_sil & high_rf_sil & high_rf_p
    )

  cat("Number of high-confidence windows (silhouette + RF pmax):\n")
  print(table(agree_dat$high_conf, useNA = "ifany"))

  agree_dat <- agree_dat %>%
    mutate(
      gmm_label  = as.character(gmm_model),
      rf_label   = as.character(rf_model),
      cont_label = as.character(best_cont_label),
      cont_label_h = harmonize_continuous_model_label(cont_label),
      cont_label_tail500 = as.character(model_cont_tail500),
      cont_label_tail500_h = harmonize_continuous_model_label(cont_label_tail500),
      cont_label_top1pct = as.character(model_cont_top1pct),
      cont_label_top1pct_h = harmonize_continuous_model_label(cont_label_top1pct),
      cont_distance_label_origin = as.character(continuous_distance_label_origin),
      cont_distance_label_origin_h = harmonize_continuous_model_label(cont_distance_label_origin),
      cont_distance_label_empirical = as.character(continuous_distance_label_empirical),
      cont_distance_label_empirical_h = harmonize_continuous_model_label(cont_distance_label_empirical),
      agree_gmm_rf = !is.na(gmm_label) & !is.na(rf_label) & (gmm_label == rf_label),
      agree_gmm_cont = !is.na(gmm_label) & !is.na(cont_label_h) & (gmm_label == cont_label_h),
      agree_rf_cont = !is.na(rf_label) & !is.na(cont_label_h) & (rf_label == cont_label_h),
      agree_all_three = agree_gmm_rf & agree_gmm_cont & agree_rf_cont,
      agree_gmm_rf_contdist = agree_gmm_rf & !is.na(cont_distance_label_empirical_h) & (rf_label == cont_distance_label_empirical_h),
      agree_all_three_top1pct = agree_gmm_rf & !is.na(cont_label_top1pct_h) &
        cont_label_top1pct != "Neutral/other" & rf_label == cont_label_top1pct_h,
      agree_all_three_tail500 = agree_gmm_rf & !is.na(cont_label_tail500_h) &
        cont_label_tail500 != "Neutral/other" & rf_label == cont_label_tail500_h,
      agree_all_three_margin_025 = agree_all_three & is.finite(cont_margin) & cont_margin >= 0.25,
      agree_all_three_margin_050 = agree_all_three & is.finite(cont_margin) & cont_margin >= 0.50,
      agree_all_three_margin_100 = agree_all_three & is.finite(cont_margin) & cont_margin >= 1.00,
      conservative_distance_margin_025 = agree_gmm_rf_contdist &
        rf_label %in% c("Geographic sweep", "Balancing selection") &
        is.finite(continuous_distance_margin_empirical) &
        continuous_distance_margin_empirical >= 0.25,
      conservative_distance_margin_050 = agree_gmm_rf_contdist &
        rf_label %in% c("Geographic sweep", "Balancing selection") &
        is.finite(continuous_distance_margin_empirical) &
        continuous_distance_margin_empirical >= 0.50,
      conservative_distance_margin_100 = agree_gmm_rf_contdist &
        rf_label %in% c("Geographic sweep", "Balancing selection") &
        is.finite(continuous_distance_margin_empirical) &
        continuous_distance_margin_empirical >= 1.00,
      score_margin_selected_025 = agree_gmm_rf &
        rf_label %in% c("Geographic sweep", "Balancing selection") &
        !is.na(continuous_margin_label_025) & continuous_margin_label_025 == rf_label,
      score_margin_selected_050 = agree_gmm_rf &
        rf_label %in% c("Geographic sweep", "Balancing selection") &
        !is.na(continuous_margin_label_050) & continuous_margin_label_050 == rf_label,
      score_margin_selected_100 = agree_gmm_rf &
        rf_label %in% c("Geographic sweep", "Balancing selection") &
        !is.na(continuous_margin_label_100) & continuous_margin_label_100 == rf_label
    )

  agree_summary_high <- agree_dat %>%
    filter(high_conf) %>%
    summarise(
      n_windows      = n(),
      frac_gmm_rf    = mean(agree_gmm_rf,    na.rm = TRUE),
      frac_gmm_cont  = mean(agree_gmm_cont,  na.rm = TRUE),
      frac_rf_cont   = mean(agree_rf_cont,   na.rm = TRUE),
      frac_all_three = mean(agree_all_three, na.rm = TRUE)
    )

  cat("\n[", level_tag, " 12A] Agreement for high-confidence windows:\n", sep = "")
  print(as.data.frame(agree_summary_high))

  agreement_tier_counts <- c(
    sum(agree_dat$agree_gmm_rf, na.rm = TRUE),
    sum(agree_dat$agree_all_three, na.rm = TRUE),
    sum(agree_dat$agree_all_three_margin_025, na.rm = TRUE),
    sum(agree_dat$agree_all_three_margin_050, na.rm = TRUE),
    sum(agree_dat$agree_all_three_margin_100, na.rm = TRUE),
    sum(agree_dat$agree_all_three_top1pct, na.rm = TRUE),
    sum(agree_dat$agree_all_three_tail500, na.rm = TRUE)
  )
  agreement_tiers <- tibble::tibble(
    level_tag = level_tag,
    tier = c(
      "GMM/RF agreed windows",
      "GMM/RF + continuous best-label agreed windows",
      "GMM/RF + continuous best-label + margin >= 0.25",
      "GMM/RF + continuous best-label + margin >= 0.50",
      "GMM/RF + continuous best-label + margin >= 1.00",
      "GMM/RF + continuous top-1% agreed windows",
      "Current strict top-500 triple-agreed windows"
    ),
    n_windows = agreement_tier_counts,
    frac_windows = agreement_tier_counts / nrow(agree_dat)
  )
  data.table::fwrite(agreement_tiers, paste0("Continuous_Agreement_Tiers_", level_tag, ".tsv"), sep = "\t")
  cat("\n[", level_tag, " 12B] Continuous agreement tiers:\n", sep = "")
  print(as.data.frame(agreement_tiers))

  agreement_by_class <- agree_dat %>%
    filter(agree_gmm_rf, !is.na(rf_label)) %>%
    group_by(level_tag = level_tag, gmm_rf_label = rf_label) %>%
    summarise(
      n_gmm_rf = n(),
      n_best_cont = sum(agree_all_three, na.rm = TRUE),
      n_best_margin_025 = sum(agree_all_three_margin_025, na.rm = TRUE),
      n_best_margin_050 = sum(agree_all_three_margin_050, na.rm = TRUE),
      n_best_margin_100 = sum(agree_all_three_margin_100, na.rm = TRUE),
      n_top1pct = sum(agree_all_three_top1pct, na.rm = TRUE),
      n_tail500 = sum(agree_all_three_tail500, na.rm = TRUE),
      frac_best_cont = n_best_cont / n_gmm_rf,
      frac_top1pct = n_top1pct / n_gmm_rf,
      frac_tail500 = n_tail500 / n_gmm_rf,
      .groups = "drop"
    )
  data.table::fwrite(agreement_by_class, paste0("Continuous_Agreement_ByClass_", level_tag, ".tsv"), sep = "\t")
  cat("\n[", level_tag, " 12C] Continuous agreement by GMM/RF agreed class:\n", sep = "")
  print(as.data.frame(agreement_by_class))

  cat("\n========== PIPELINE COMPLETE FOR LEVEL ", level_tag, " ==========\n\n", sep = "")

  invisible(list(
    sd3_clusters   = sd3_clusters,
    sd3_cont       = sd3_cont,
    centroids_z5   = centroids_z5,
    centroids_fit5 = centroids_fit5,
    gmm            = gmm
  ))

}

############################################################
## RUN BOTH LEVELS: G = 2 AND G = bestG
############################################################

res_G2 <- run_level_pipeline("G2", 2)

if (!is.na(bestG) && bestG > 2) {
  res_Gbest <- run_level_pipeline(paste0("G", bestG), bestG)
  sd3_clusters <- res_Gbest$sd3_clusters
  sd3_cont     <- res_Gbest$sd3_cont
  tag_for_diag <- paste0("G", bestG)
} else {
  cat("\n[MAIN] bestG is 2 or NA; skipping separate Gbest run.\n")
  sd3_clusters <- res_G2$sd3_clusters
  sd3_cont     <- res_G2$sd3_cont
  tag_for_diag <- "G2"
}

cat("\n========== UPDATED MASTER PIPELINE COMPLETE v2026-02-06-01 ==========\n")
cat("Chosen bestG from stability-based selection (prefer ≥3):", bestG, "\n")
cat("Diagnostics below use level:", tag_for_diag, "\n")

############################################################
## [15] RECOMBINATION AS A PREDICTOR OF MODEL (NON-CIRCULAR TESTS)
############################################################
## Hypothesized ordering (general expectation; test empirically):
##   Lowest recombination:   Sweep-like and Allopatry-like
##     - both reflect linked selection / long haplotypes / reduced diversity islands, often enriched in low-recomb.
##   Higher recombination:   DGF-like
##     - gene flow more likely to persist where recombination decouples barrier and background; introgression tends to be
##       enriched in higher recombination regions in many systems.
##   Highest recombination:  Balancing selection (often maintained as a local signal rather than a long linked block)
##     - balancing loci can occur anywhere, but detectability and persistence of linked diversity typically improves with recombination.
##
## Test plan implemented below:
##   A) Nonparametric: Kruskal-Wallis + pairwise Wilcoxon on recombination ~ GMM model
##   B) Binary logits: each model vs Neutral (or vs others), with recombination + optional covariates
##   C) Multinomial logit: model ~ recombination (+ covariates) using nnet::multinom if available
##   D) Effect sizes: standardized odds ratios and McFadden pseudo-R2 for nested models
############################################################

cat("\n================ [15] RECOMBINATION PREDICTS MODEL? ==================\n")

rec_dat <- sd3_clusters %>%
  select(gmm_model_raw, rf_model, recombination, z_rec, z_fst, z_dxy, z_pi, z_f3, z_taj, p_cluster) %>%
  filter(is.finite(z_rec), !is.na(gmm_model_raw)) %>%
  mutate(
    gmm_model_raw = factor(gmm_model_raw),
    rf_model = as.character(rf_model)
  )

cat("Rows with finite rank-normalized recombination (z_rec) and GMM model:", nrow(rec_dat), "\n")
cat("GMM model counts:\n")
print(table(rec_dat$gmm_model_raw, useNA="ifany"))

## A) Kruskal-Wallis (GMM labels)
if (length(unique(rec_dat$gmm_model_raw)) > 1) {
  cat("\n[15A] Kruskal-Wallis: rank-normalized recombination (z_rec) ~ gmm_model_raw\n")
  print(kruskal.test(z_rec ~ gmm_model_raw, data = rec_dat))

  cat("\n[15A2] Pairwise Wilcoxon (BH-adjusted):\n")
  pw <- pairwise.wilcox.test(rec_dat$z_rec, rec_dat$gmm_model_raw, p.adjust.method = "BH")
  print(pw)
}

## B) Binary logits vs Neutral (if Neutral present)
##    Uses rank-normalized recombination (z_rec); ORs are per +1 SD by construction
if ("Neutral / equilibrium" %in% levels(rec_dat$gmm_model_raw)) {
  cat("\n[15B] Binary logits vs Neutral (one-vs-rest among non-neutral models)\n")

  for (m in setdiff(levels(rec_dat$gmm_model_raw), "Neutral / equilibrium")) {
    dfm <- rec_dat %>%
      filter(gmm_model_raw %in% c("Neutral / equilibrium", m)) %>%
      mutate(y = as.integer(gmm_model_raw == m))

    if (length(unique(dfm$y)) < 2) next

    fit1 <- glm(y ~ z_rec, data = dfm, family = binomial())
    fit0 <- glm(y ~ 1,    data = dfm, family = binomial())

    ## McFadden pseudo-R2
    R2m <- 1 - (as.numeric(logLik(fit1)) / as.numeric(logLik(fit0)))

    OR  <- exp(coef(fit1)["z_rec"])
    CI  <- tryCatch(exp(confint(fit1)["z_rec", ]), error = function(e) c(NA, NA))

    cat("\nModel:", m, " vs Neutral\n")
    cat("  OR per +1 SD recombination (z_rec):", OR, "\n")
    cat("  95% CI:", CI[1], "-", CI[2], "\n")
    cat("  McFadden R2:", R2m, "\n")
    cat("  Summary:\n")
    print(summary(fit1))
  }
} else {
  cat("\n[15B] Neutral / equilibrium not present among GMM labels; skipping vs-neutral logits.\n")
}

## C) Multinomial logit: model ~ z_rec (+ optional covariates)
##    Uses nnet::multinom if installed.
if (requireNamespace("nnet", quietly = TRUE)) {
  cat("\n[15C] Multinomial logit: gmm_model_raw ~ z_rec (+ z_fst + z_dxy + z_pi)\n")
  suppressPackageStartupMessages(library(nnet))

  dfm <- rec_dat %>%
    mutate(gmm_model_raw = droplevels(gmm_model_raw)) %>%
    filter(!is.na(gmm_model_raw), is.finite(z_rec),
           is.finite(z_fst), is.finite(z_dxy), is.finite(z_pi))

  if (nrow(dfm) > 0 && length(unique(dfm$gmm_model_raw)) >= 3) {

    ## base vs full to get pseudo-R2 style deviance change
    m0 <- nnet::multinom(gmm_model_raw ~ 1, data = dfm, trace = FALSE)
    m1 <- nnet::multinom(gmm_model_raw ~ z_rec, data = dfm, trace = FALSE)
    m2 <- nnet::multinom(gmm_model_raw ~ z_rec + z_fst + z_dxy + z_pi, data = dfm, trace = FALSE)

    dev0 <- 2 * as.numeric(logLik(m0))
    dev1 <- 2 * as.numeric(logLik(m1))
    dev2 <- 2 * as.numeric(logLik(m2))

    ## McFadden-like R2 against null
    R2_1 <- 1 - (as.numeric(logLik(m1)) / as.numeric(logLik(m0)))
    R2_2 <- 1 - (as.numeric(logLik(m2)) / as.numeric(logLik(m0)))

    cat("  McFadden R2 (z_rec only):", R2_1, "\n")
    cat("  McFadden R2 (z_rec + covariates):", R2_2, "\n")

    cat("\n  Coefficients (z_rec model):\n")
    print(summary(m1))

    cat("\n  Coefficients (controlled model):\n")
    print(summary(m2))

  } else {
    cat("  Not enough rows or classes for multinomial model.\n")
  }
} else {
  cat("\n[15C] Package 'nnet' not installed; skipping multinomial.\n")
}

## D) Visual: recombination distributions by model (GMM labels)
p_rec_by_model <- ggplot(rec_dat, aes(x = gmm_model_raw, y = z_rec, fill = gmm_model_raw)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = NULL, y = "Recombination (rank-normalized z)", title = "Recombination by GMM-inferred model (non-circular)")
show_plot(p_rec_by_model)

############################################################
## [14] AGREEMENT MODULE (kept; includes updated [14H] tick plot)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
})

alpha <- 0.05

## silhouette filtering controls
sil_filter_mode <- "quantile"  # "quantile" or "hard" or "none"
sil_thr_hard    <- 0.7         # used only if sil_filter_mode == "hard"
sil_quantile    <- 0.75        # used only if sil_filter_mode == "quantile" (top 25%)

## continuous label harmonization
harmonize_continuous_labels <- TRUE

## ---- [14A] Standardize keys in sd3_clusters ----
if (!exists("sd3_clusters")) stop("sd3_clusters not found in environment.")

chr_col   <- .find_chr_col(sd3_clusters)
start_col <- .find_start_col(sd3_clusters)
end_col   <- .find_end_col(sd3_clusters)

if (is.na(chr_col) || is.na(start_col)) {
  stop(
    "Could not identify key columns in sd3_clusters.\n",
    "Found chr_col = ", chr_col, " ; start_col = ", start_col, "\n",
    "Available columns (first 80): ", paste(head(names(sd3_clusters), 80), collapse = ", ")
  )
}

sd3c <- sd3_clusters %>%
  .rename_keys_to_standard(chr_col, start_col, end_col) %>%
  mutate(
    Chromosome_Names = as.character(Chromosome_Names),
    .start_bp = as.numeric(.start_bp),
    .end_bp   = if (".end_bp" %in% names(.)) as.numeric(.end_bp) else NA_real_
  )

## ---- [14B] Bring in continuous labels (prefer sd3_cont) ----
cont_candidates <- c("best_cont_label","model_cont_best","cont_best_label",
                     "model_cont","model_continuous","cont_model","continuous_model",
                     "cont_label","continuous_label","model_cont_raw","cont_model_raw")
cont_extra_candidates <- c(
  "sweep_score", "bal_score",
  "model_cont", "model_cont_tail500", "best_cont_label", "best_cont_score",
  "second_cont_score", "cont_margin", "cont_score_quantile",
  "model_cont_top1pct", "model_cont_top5pct", "model_cont_top10pct",
  "continuous_distance_label_origin", "continuous_distance_label_empirical",
  "continuous_distance_margin_origin", "continuous_distance_margin_empirical",
  "d_neutral_origin", "d_neutral_empirical", "d_balancing", "d_sweep",
  "continuous_score_margin",
  "continuous_margin_label_025", "continuous_margin_label_050",
  "continuous_margin_label_100", "continuous_margin_label"
)

cont_col_in_clusters <- .find_col(sd3c, cont_candidates)

if (is.na(cont_col_in_clusters)) {
  if (!exists("sd3_cont")) stop("No continuous label column found in sd3_clusters AND sd3_cont not found.")

  chr2   <- .find_chr_col(sd3_cont)
  start2 <- .find_start_col(sd3_cont)
  end2   <- .find_end_col(sd3_cont)

  sd3t2 <- sd3_cont %>%
    .rename_keys_to_standard(chr2, start2, end2) %>%
    mutate(
      Chromosome_Names = as.character(Chromosome_Names),
      .start_bp = as.numeric(.start_bp),
      .end_bp   = if (".end_bp" %in% names(.)) as.numeric(.end_bp) else NA_real_
    )

  cont_col_in_cont <- .find_col(sd3t2, cont_candidates)
  if (is.na(cont_col_in_cont)) stop("Could not find continuous label in sd3_cont.")

  cont_join_cols <- unique(c(cont_col_in_cont, intersect(cont_extra_candidates, names(sd3t2))))

  sd3c <- sd3c %>%
    left_join(
      sd3t2 %>%
        select(Chromosome_Names, .start_bp, all_of(cont_join_cols)) %>%
        rename(cont_label_raw = !!rlang::sym(cont_col_in_cont)),
      by = c("Chromosome_Names", ".start_bp")
    )

  cont_col_in_clusters <- "cont_label_raw"
}

cont_tail_col <- .find_col(sd3c, c("model_cont_tail500", "model_cont"))
cont_top1_col <- .find_col(sd3c, c("model_cont_top1pct"))
cont_margin_col <- .find_col(sd3c, c("cont_margin"))
cont_quantile_col <- .find_col(sd3c, c("cont_score_quantile"))
cont_dist_origin_col <- .find_col(sd3c, c("continuous_distance_label_origin"))
cont_dist_emp_col <- .find_col(sd3c, c("continuous_distance_label_empirical"))
cont_dist_margin_origin_col <- .find_col(sd3c, c("continuous_distance_margin_origin"))
cont_dist_margin_emp_col <- .find_col(sd3c, c("continuous_distance_margin_empirical"))
dist_neutral_origin_col <- .find_col(sd3c, c("d_neutral_origin"))
dist_neutral_emp_col <- .find_col(sd3c, c("d_neutral_empirical"))
dist_bal_col <- .find_col(sd3c, c("d_balancing"))
dist_sweep_col <- .find_col(sd3c, c("d_sweep"))
score_margin_col <- .find_col(sd3c, c("continuous_score_margin"))
margin_lab_025_col <- .find_col(sd3c, c("continuous_margin_label_025"))
margin_lab_050_col <- .find_col(sd3c, c("continuous_margin_label_050", "continuous_margin_label"))
margin_lab_100_col <- .find_col(sd3c, c("continuous_margin_label_100"))

## ---- [14C] Find label + silhouette columns ----
gmm_label_col <- .find_col(sd3c, c("gmm_model_raw","gmm_label","gmm_model","gmm_class"))
rf_label_col  <- .find_col(sd3c, c("rf_model","rf_label","rf_class","rf_model_raw"))

if (is.na(gmm_label_col) || is.na(rf_label_col)) {
  stop("Could not find required label columns in sd3_clusters.")
}

gmm_sil_col <- .find_col(sd3c, c("gmm_silhouette","sil_gmm","silhouette_gmm","gmm_sil"))
rf_sil_col  <- .find_col(sd3c, c("rf_silhouette","sil_rf","silhouette_rf","rf_sil"))

## ---- [14D] Build agreement flags (+ harmonized cont labels) ----
agree_dat <- sd3c %>%
  mutate(
    gmm_label  = as.character(.data[[gmm_label_col]]),
    rf_label   = as.character(.data[[rf_label_col]]),
    cont_label = as.character(.data[[cont_col_in_clusters]]),

    cont_label_h = dplyr::case_when(
      !harmonize_continuous_labels ~ cont_label,
      TRUE ~ harmonize_continuous_model_label(cont_label)
    ),
    cont_label_tail500 = if (!is.na(cont_tail_col)) as.character(.data[[cont_tail_col]]) else NA_character_,
    cont_label_tail500_h = dplyr::case_when(
      !harmonize_continuous_labels ~ cont_label_tail500,
      TRUE ~ harmonize_continuous_model_label(cont_label_tail500)
    ),
    cont_label_top1pct = if (!is.na(cont_top1_col)) as.character(.data[[cont_top1_col]]) else NA_character_,
    cont_label_top1pct_h = dplyr::case_when(
      !harmonize_continuous_labels ~ cont_label_top1pct,
      TRUE ~ harmonize_continuous_model_label(cont_label_top1pct)
    ),
    cont_margin = if (!is.na(cont_margin_col)) as.numeric(.data[[cont_margin_col]]) else NA_real_,
    cont_score_quantile = if (!is.na(cont_quantile_col)) as.numeric(.data[[cont_quantile_col]]) else NA_real_,
    continuous_distance_label_origin = if (!is.na(cont_dist_origin_col)) as.character(.data[[cont_dist_origin_col]]) else NA_character_,
    continuous_distance_label_origin_h = harmonize_continuous_model_label(continuous_distance_label_origin),
    continuous_distance_label_empirical = if (!is.na(cont_dist_emp_col)) as.character(.data[[cont_dist_emp_col]]) else NA_character_,
    continuous_distance_label_empirical_h = harmonize_continuous_model_label(continuous_distance_label_empirical),
    continuous_distance_margin_origin = if (!is.na(cont_dist_margin_origin_col)) as.numeric(.data[[cont_dist_margin_origin_col]]) else NA_real_,
    continuous_distance_margin_empirical = if (!is.na(cont_dist_margin_emp_col)) as.numeric(.data[[cont_dist_margin_emp_col]]) else NA_real_,
    d_neutral_origin = if (!is.na(dist_neutral_origin_col)) as.numeric(.data[[dist_neutral_origin_col]]) else NA_real_,
    d_neutral_empirical = if (!is.na(dist_neutral_emp_col)) as.numeric(.data[[dist_neutral_emp_col]]) else NA_real_,
    d_balancing = if (!is.na(dist_bal_col)) as.numeric(.data[[dist_bal_col]]) else NA_real_,
    d_sweep = if (!is.na(dist_sweep_col)) as.numeric(.data[[dist_sweep_col]]) else NA_real_,
    continuous_score_margin = if (!is.na(score_margin_col)) as.numeric(.data[[score_margin_col]]) else NA_real_,
    continuous_margin_label_025 = if (!is.na(margin_lab_025_col)) as.character(.data[[margin_lab_025_col]]) else NA_character_,
    continuous_margin_label_050 = if (!is.na(margin_lab_050_col)) as.character(.data[[margin_lab_050_col]]) else NA_character_,
    continuous_margin_label_100 = if (!is.na(margin_lab_100_col)) as.character(.data[[margin_lab_100_col]]) else NA_character_,

    triple_match =
      !is.na(gmm_label) & !is.na(rf_label) & !is.na(cont_label_h) &
      (gmm_label == rf_label) & (rf_label == cont_label_h),
    gmm_rf_match =
      !is.na(gmm_label) & !is.na(rf_label) &
      (gmm_label == rf_label),
    triple_tail500 =
      gmm_rf_match & !is.na(cont_label_tail500_h) &
      cont_label_tail500 != "Neutral/other" &
      (rf_label == cont_label_tail500_h),
    triple_top1pct =
      gmm_rf_match & !is.na(cont_label_top1pct_h) &
      cont_label_top1pct != "Neutral/other" &
      (rf_label == cont_label_top1pct_h),
    triple_margin_025 = triple_match & is.finite(cont_margin) & cont_margin >= 0.25,
    triple_margin_050 = triple_match & is.finite(cont_margin) & cont_margin >= 0.50,
    triple_margin_100 = triple_match & is.finite(cont_margin) & cont_margin >= 1.00,
    contdist_match =
      gmm_rf_match & !is.na(continuous_distance_label_empirical_h) &
      (rf_label == continuous_distance_label_empirical_h),
    dist_margin_025 =
      contdist_match &
      rf_label %in% c("Geographic sweep", "Balancing selection") &
      is.finite(continuous_distance_margin_empirical) &
      continuous_distance_margin_empirical >= 0.25,
    dist_margin_050 =
      contdist_match &
      rf_label %in% c("Geographic sweep", "Balancing selection") &
      is.finite(continuous_distance_margin_empirical) &
      continuous_distance_margin_empirical >= 0.50,
    dist_margin_100 =
      contdist_match &
      rf_label %in% c("Geographic sweep", "Balancing selection") &
      is.finite(continuous_distance_margin_empirical) &
      continuous_distance_margin_empirical >= 1.00,
    score_margin_025 =
      gmm_rf_match &
      rf_label %in% c("Geographic sweep", "Balancing selection") &
      !is.na(continuous_margin_label_025) &
      continuous_margin_label_025 == rf_label,
    score_margin_050 =
      gmm_rf_match &
      rf_label %in% c("Geographic sweep", "Balancing selection") &
      !is.na(continuous_margin_label_050) &
      continuous_margin_label_050 == rf_label,
    score_margin_100 =
      gmm_rf_match &
      rf_label %in% c("Geographic sweep", "Balancing selection") &
      !is.na(continuous_margin_label_100) &
      continuous_margin_label_100 == rf_label,

    gmm_silhouette = if (!is.na(gmm_sil_col)) as.numeric(.data[[gmm_sil_col]]) else NA_real_,
    rf_silhouette  = if (!is.na(rf_sil_col))  as.numeric(.data[[rf_sil_col]])  else NA_real_,

    min_sil = pmin(gmm_silhouette, rf_silhouette,na.rm=T),
    min_sil = ifelse(is.infinite(min_sil), NA_real_, min_sil),

    sil_thr_effective = dplyr::case_when(
      sil_filter_mode == "hard"     ~ sil_thr_hard,
      sil_filter_mode == "quantile" ~ as.numeric(quantile(min_sil, probs = sil_quantile, na.rm = TRUE)),
      TRUE                          ~ NA_real_
    ),

    sil_ok = dplyr::case_when(
      sil_filter_mode == "none" ~ TRUE,
      is.na(min_sil) ~ TRUE,  # <-- treat missing silhouette as "pass"
      is.finite(min_sil) & is.finite(sil_thr_effective) ~ (min_sil >= sil_thr_effective),
      TRUE ~ FALSE
    ),

    triple_sil = triple_match & sil_ok
  ) %>%
  dplyr::group_by(rf_label) %>%
  dplyr::mutate(
    dist_top1pct = dplyr::case_when(
      rf_label %in% c("Geographic sweep", "Balancing selection") & dist_margin_050 ~ {
        evidence <- dplyr::if_else(
          rf_label == "Geographic sweep",
          sweep_score,
          dplyr::if_else(rf_label == "Balancing selection", bal_score, NA_real_)
        )
        rank_frac <- rank(-evidence, ties.method = "average", na.last = "keep") / sum(is.finite(evidence))
        is.finite(rank_frac) & rank_frac <= 0.01
      },
      TRUE ~ FALSE
    ),
    dist_top500 = dplyr::case_when(
      rf_label %in% c("Geographic sweep", "Balancing selection") & dist_margin_050 ~ {
        evidence <- dplyr::if_else(
          rf_label == "Geographic sweep",
          sweep_score,
          dplyr::if_else(rf_label == "Balancing selection", bal_score, NA_real_)
        )
        rk <- rank(-evidence, ties.method = "first", na.last = "keep")
        is.finite(rk) & rk <= 500
      },
      TRUE ~ FALSE
    )
  ) %>%
  dplyr::ungroup() %>%
  select(
    Chromosome_Names, .start_bp, dplyr::any_of(".end_bp"),
    dplyr::any_of(c(
      "z_dxy", "z_fst", "z_pi", "z_f3", "z_taj", "z_rec",
      "avg_dxy", "avg_wc_fst", "pi_average", "f3", "TajimaD", "recombination",
      "abs_delta_p", "mean_abs_delta_p", "B_no_delta_p", "K_local_mean"
    )),
    gmm_label, rf_label, cont_label, cont_label_h,
    dplyr::any_of(c("sweep_score", "bal_score")),
    cont_label_tail500, cont_label_tail500_h,
    cont_label_top1pct, cont_label_top1pct_h,
    cont_margin, cont_score_quantile, best_cont_score,
    continuous_distance_label_origin, continuous_distance_label_origin_h,
    continuous_distance_label_empirical, continuous_distance_label_empirical_h,
    continuous_distance_margin_origin, continuous_distance_margin_empirical,
    d_neutral_origin, d_neutral_empirical, d_balancing, d_sweep,
    continuous_score_margin,
    continuous_margin_label_025, continuous_margin_label_050, continuous_margin_label_100,
    dplyr::any_of("rf_pmax"),
    gmm_silhouette, rf_silhouette, min_sil, sil_thr_effective,
    gmm_rf_match, triple_match, triple_tail500, triple_top1pct,
    triple_margin_025, triple_margin_050, triple_margin_100,
    contdist_match, dist_margin_025, dist_margin_050, dist_margin_100,
    score_margin_025, score_margin_050, score_margin_100,
    dist_top1pct, dist_top500,
    sil_ok, triple_sil
  )

data.table::fwrite(agree_dat, "Continuous_Agreement_Windows_Global.tsv", sep = "\t")
cat("\n[14] Wrote window-level agreement table: Continuous_Agreement_Windows_Global.tsv\n")

cat("\n[14] Counts:\n")
print(as.data.frame(agree_dat %>% summarise(
  n_total = n(),
  n_gmm_rf_match = sum(gmm_rf_match, na.rm = TRUE),
  n_triple_match = sum(triple_match, na.rm = TRUE),
  n_triple_margin_025 = sum(triple_margin_025, na.rm = TRUE),
  n_triple_margin_050 = sum(triple_margin_050, na.rm = TRUE),
  n_triple_margin_100 = sum(triple_margin_100, na.rm = TRUE),
  n_triple_top1pct = sum(triple_top1pct, na.rm = TRUE),
  n_triple_tail500 = sum(triple_tail500, na.rm = TRUE),
  n_sil_ok       = sum(sil_ok, na.rm = TRUE),
  n_triple_sil   = sum(triple_sil, na.rm = TRUE),
  sil_filter_mode = sil_filter_mode,
  sil_thr_effective = unique(na.omit(sil_thr_effective))[1]
)))

agreement_tiers_global_counts <- c(
  sum(agree_dat$gmm_rf_match, na.rm = TRUE),
  sum(agree_dat$triple_match, na.rm = TRUE),
  sum(agree_dat$triple_margin_025, na.rm = TRUE),
  sum(agree_dat$triple_margin_050, na.rm = TRUE),
  sum(agree_dat$triple_margin_100, na.rm = TRUE),
  sum(agree_dat$triple_top1pct, na.rm = TRUE),
  sum(agree_dat$triple_tail500, na.rm = TRUE)
)
agreement_tiers_global <- tibble::tibble(
  tier = c(
    "GMM/RF agreed windows",
    "GMM/RF + continuous best-label agreed windows",
    "GMM/RF + continuous best-label + margin >= 0.25",
    "GMM/RF + continuous best-label + margin >= 0.50",
    "GMM/RF + continuous best-label + margin >= 1.00",
    "GMM/RF + continuous top-1% agreed windows",
    "Current strict top-500 triple-agreed windows"
  ),
  n_windows = agreement_tiers_global_counts,
  frac_windows = agreement_tiers_global_counts / nrow(agree_dat)
)
data.table::fwrite(agreement_tiers_global, "Continuous_Agreement_Tiers_Global.tsv", sep = "\t")
cat("\n[14] Continuous agreement tiers (global bestG object):\n")
print(as.data.frame(agreement_tiers_global))

agreement_by_class_global <- agree_dat %>%
  filter(gmm_rf_match, !is.na(rf_label)) %>%
  group_by(gmm_rf_label = rf_label) %>%
  summarise(
    n_gmm_rf = n(),
    n_best_cont = sum(triple_match, na.rm = TRUE),
    n_best_margin_025 = sum(triple_margin_025, na.rm = TRUE),
    n_best_margin_050 = sum(triple_margin_050, na.rm = TRUE),
    n_best_margin_100 = sum(triple_margin_100, na.rm = TRUE),
    n_top1pct = sum(triple_top1pct, na.rm = TRUE),
    n_tail500 = sum(triple_tail500, na.rm = TRUE),
    frac_best_cont = n_best_cont / n_gmm_rf,
    frac_top1pct = n_top1pct / n_gmm_rf,
    frac_tail500 = n_tail500 / n_gmm_rf,
    .groups = "drop"
  )
data.table::fwrite(agreement_by_class_global, "Continuous_Agreement_ByClass_Global.tsv", sep = "\t")
cat("\n[14] Continuous agreement by class (global bestG object):\n")
print(as.data.frame(agreement_by_class_global))

cat("\n[14] Top triple combinations (labels):\n")
print(as.data.frame(
  agree_dat %>%
    filter(!is.na(gmm_label), !is.na(rf_label), !is.na(cont_label_h)) %>%
    dplyr::count(gmm_label, rf_label, cont_label_h, sort = TRUE) %>%
    head(20)
))

cat("\n[14] Pairwise agreement rates (using cont_label_h):\n")
print(as.data.frame(agree_dat %>% summarise(
  frac_gmm_rf   = mean(gmm_rf_match, na.rm = TRUE),
  frac_gmm_cont = mean(gmm_label == cont_label_h, na.rm = TRUE),
  frac_rf_cont  = mean(rf_label  == cont_label_h, na.rm = TRUE),
  frac_triple   = mean(triple_match, na.rm = TRUE)
)))

############################################################
## [14G] REVISED EXPLICIT-NEUTRAL TIER TABLE + OUTPUTS
############################################################

metric_pick <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) return(rep(NA_real_, nrow(df)))
  as.numeric(df[[hit[[1]]]])
}

safe_mean_num <- function(x) if (all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_median_num <- function(x) if (all(!is.finite(x))) NA_real_ else median(x, na.rm = TRUE)

agree_dat <- agree_dat %>%
  mutate(
    window_bp = dplyr::if_else(is.finite(.end_bp), pmax(.end_bp - .start_bp + 1, 0), NA_real_),
    gmm_rf_label = dplyr::if_else(gmm_rf_match, rf_label, "No consensus"),
    old_best_continuous_3of3_label = dplyr::if_else(triple_match, rf_label, "No consensus"),
    continuous_distance_3class_label = dplyr::if_else(
      !is.na(continuous_distance_label_empirical_h),
      continuous_distance_label_empirical_h,
      "No consensus"
    ),
    gmm_rf_contdist_label = dplyr::if_else(contdist_match, rf_label, "No consensus"),
    conservative_distance_margin_025_label = dplyr::if_else(dist_margin_025, rf_label, "No consensus"),
    conservative_distance_margin_050_label = dplyr::if_else(dist_margin_050, rf_label, "No consensus"),
    conservative_distance_margin_100_label = dplyr::if_else(dist_margin_100, rf_label, "No consensus"),
    score_margin_selected_025_label = dplyr::if_else(score_margin_025, rf_label, "No consensus"),
    score_margin_selected_050_label = dplyr::if_else(score_margin_050, rf_label, "No consensus"),
    score_margin_selected_100_label = dplyr::if_else(score_margin_100, rf_label, "No consensus"),
    top1pct_distance_margin_050_label = dplyr::if_else(dist_top1pct, rf_label, "No consensus"),
    top500_distance_margin_050_label = dplyr::if_else(dist_top500, rf_label, "No consensus")
  )

continuous_distance_centroids_global <- dplyr::bind_rows(
  tibble::tibble(
    centroid_type = "neutral_origin",
    z_dxy = 0, z_fst = 0, z_pi = 0, z_f3 = 0, z_taj = 0, z_rec = 0
  ),
  agree_dat %>%
    filter(gmm_rf_match, rf_label == "Neutral / equilibrium") %>%
    summarise(
      centroid_type = "neutral_empirical",
      z_dxy = safe_mean_num(z_dxy), z_fst = safe_mean_num(z_fst), z_pi = safe_mean_num(z_pi),
      z_f3 = safe_mean_num(z_f3), z_taj = safe_mean_num(z_taj), z_rec = safe_mean_num(z_rec)
    ),
  agree_dat %>%
    filter(gmm_rf_match, rf_label == "Balancing selection") %>%
    summarise(
      centroid_type = "balancing",
      z_dxy = safe_mean_num(z_dxy), z_fst = safe_mean_num(z_fst), z_pi = safe_mean_num(z_pi),
      z_f3 = safe_mean_num(z_f3), z_taj = safe_mean_num(z_taj), z_rec = safe_mean_num(z_rec)
    ),
  agree_dat %>%
    filter(gmm_rf_match, rf_label == "Geographic sweep") %>%
    summarise(
      centroid_type = "sweep",
      z_dxy = safe_mean_num(z_dxy), z_fst = safe_mean_num(z_fst), z_pi = safe_mean_num(z_pi),
      z_f3 = safe_mean_num(z_f3), z_taj = safe_mean_num(z_taj), z_rec = safe_mean_num(z_rec)
    )
)
data.table::fwrite(continuous_distance_centroids_global, "continuous_distance_centroids.tsv", sep = "\t")

continuous_distance_window_scores <- agree_dat %>%
  select(
    Chromosome_Names, .start_bp, dplyr::any_of(".end_bp"),
    z_dxy, z_fst, z_pi, z_f3, z_taj, z_rec,
    d_neutral_origin, d_neutral_empirical, d_balancing, d_sweep,
    continuous_distance_label_origin, continuous_distance_label_empirical,
    continuous_distance_margin_origin, continuous_distance_margin_empirical
  )
data.table::fwrite(continuous_distance_window_scores, "continuous_distance_window_scores.tsv", sep = "\t")

continuous_margin_rejection_labels <- agree_dat %>%
  select(
    Chromosome_Names, .start_bp, dplyr::any_of(".end_bp"),
    dplyr::any_of(c("bal_score", "sweep_score")),
    continuous_score_margin,
    continuous_margin_label_025, continuous_margin_label_050, continuous_margin_label_100
  )
data.table::fwrite(continuous_margin_rejection_labels, "continuous_margin_rejection_labels.tsv", sep = "\t")

revised_selection_assignments <- agree_dat %>%
  select(
    Chromosome_Names, .start_bp, dplyr::any_of(".end_bp"), window_bp,
    gmm_label, rf_label,
    gmm_rf_label,
    old_best_continuous_3of3_label,
    continuous_distance_3class_label,
    gmm_rf_contdist_label,
    conservative_distance_margin_025_label,
    conservative_distance_margin_050_label,
    conservative_distance_margin_100_label,
    continuous_margin_label_025, continuous_margin_label_050, continuous_margin_label_100,
    score_margin_selected_025_label,
    score_margin_selected_050_label,
    score_margin_selected_100_label,
    top1pct_distance_margin_050_label,
    top500_distance_margin_050_label,
    continuous_distance_margin_empirical,
    continuous_score_margin
  )
data.table::fwrite(revised_selection_assignments, "revised_selection_tier_window_assignments.tsv", sep = "\t")

tier_specs_revised <- tibble::tribble(
  ~tier, ~label_col, ~recommended_paper_use,
  "GMM/RF", "gmm_rf_label", "broad genomic state",
  "Old best continuous 3-of-3", "old_best_continuous_3of3_label", "legacy comparison only",
  "Continuous-distance 3-class", "continuous_distance_3class_label", "continuous model with explicit neutral state",
  "GMM/RF + continuous-distance agreement", "gmm_rf_contdist_label", "main 3-class agreement model",
  "Distance margin >= 0.25", "conservative_distance_margin_025_label", "sensitivity analysis",
  "Distance margin >= 0.50", "conservative_distance_margin_050_label", "primary conservative selected-window map",
  "Distance margin >= 1.00", "conservative_distance_margin_100_label", "strict sensitivity analysis",
  "Score margin >= 0.25", "score_margin_selected_025_label", "alternative selection-only sensitivity",
  "Score margin >= 0.50", "score_margin_selected_050_label", "alternative selection-only map",
  "Score margin >= 1.00", "score_margin_selected_100_label", "strict alternative sensitivity",
  "Top 1%", "top1pct_distance_margin_050_label", "high-confidence sensitivity",
  "Top 500", "top500_distance_margin_050_label", "extreme candidate loci only"
)

summarise_tier_labels <- function(df, tier_name, label_col) {
  x <- df %>%
    mutate(
      tier = tier_name,
      tier_label = as.character(.data[[label_col]]),
      tier_label = ifelse(is.na(tier_label) | !nzchar(tier_label), "No consensus", tier_label),
      avg_dxy_metric = metric_pick(., c("avg_dxy")),
      avg_wc_fst_metric = metric_pick(., c("avg_wc_fst")),
      pi_average_metric = metric_pick(., c("pi_average")),
      f3_metric = metric_pick(., c("f3")),
      TajimaD_metric = metric_pick(., c("TajimaD")),
      recombination_metric = metric_pick(., c("recombination")),
      abs_delta_p_metric = metric_pick(., c("abs_delta_p", "mean_abs_delta_p")),
      B_no_delta_p_metric = metric_pick(., c("B_no_delta_p")),
      K_local_mean_metric = metric_pick(., c("K_local_mean"))
    )
  x %>%
    group_by(tier, class = tier_label) %>%
    summarise(
      n_windows = n(),
      bp_covered = sum(window_bp, na.rm = TRUE),
      n_chromosomes = dplyr::n_distinct(Chromosome_Names),
      mean_avg_dxy = safe_mean_num(avg_dxy_metric),
      median_avg_dxy = safe_median_num(avg_dxy_metric),
      mean_avg_wc_fst = safe_mean_num(avg_wc_fst_metric),
      median_avg_wc_fst = safe_median_num(avg_wc_fst_metric),
      mean_pi_average = safe_mean_num(pi_average_metric),
      median_pi_average = safe_median_num(pi_average_metric),
      mean_f3 = safe_mean_num(f3_metric),
      median_f3 = safe_median_num(f3_metric),
      mean_TajimaD = safe_mean_num(TajimaD_metric),
      median_TajimaD = safe_median_num(TajimaD_metric),
      mean_recombination = safe_mean_num(recombination_metric),
      median_recombination = safe_median_num(recombination_metric),
      mean_abs_delta_p = safe_mean_num(abs_delta_p_metric),
      median_abs_delta_p = safe_median_num(abs_delta_p_metric),
      mean_B_no_delta_p = safe_mean_num(B_no_delta_p_metric),
      median_B_no_delta_p = safe_median_num(B_no_delta_p_metric),
      mean_K_local_mean = safe_mean_num(K_local_mean_metric),
      median_K_local_mean = safe_median_num(K_local_mean_metric),
      .groups = "drop"
    )
}

summarise_tier_chromosomes <- function(df, tier_name, label_col) {
  df %>%
    mutate(
      tier = tier_name,
      tier_label = as.character(.data[[label_col]]),
      tier_label = ifelse(is.na(tier_label) | !nzchar(tier_label), "No consensus", tier_label)
    ) %>%
    group_by(tier, Chromosome_Names, class = tier_label) %>%
    summarise(
      n_windows = n(),
      bp_covered = sum(window_bp, na.rm = TRUE),
      .groups = "drop"
    )
}

revised_selection_tier_counts <- purrr::map2_dfr(
  tier_specs_revised$tier,
  tier_specs_revised$label_col,
  ~ summarise_tier_labels(agree_dat, .x, .y)
)
data.table::fwrite(revised_selection_tier_counts, "revised_selection_tier_counts.tsv", sep = "\t")

revised_selection_tier_chromosomes <- purrr::map2_dfr(
  tier_specs_revised$tier,
  tier_specs_revised$label_col,
  ~ summarise_tier_chromosomes(agree_dat, .x, .y)
)
data.table::fwrite(revised_selection_tier_chromosomes, "revised_selection_tier_chromosome_counts.tsv", sep = "\t")

decision_table_core <- tier_specs_revised %>%
  rowwise() %>%
  do({
    tt <- summarise_tier_labels(agree_dat, .$tier, .$label_col)
    get_class_n <- function(cls) {
      val <- tt$n_windows[tt$class == cls]
      if (length(val) == 0L) 0L else val[[1]]
    }
    tibble::tibble(
      tier = .$tier,
      neutral_windows = get_class_n("Neutral / equilibrium"),
      balancing_windows = get_class_n("Balancing selection"),
      sweep_windows = get_class_n("Geographic sweep"),
      ambiguous_no_consensus_windows = sum(tt$n_windows[tt$class %in% c("Ambiguous", "No consensus")], na.rm = TRUE),
      balancing_genes = NA_real_,
      sweep_genes = NA_real_,
      BGS_like = NA_real_,
      sweep_like_excess = NA_real_,
      abs_delta_p_R2 = NA_real_,
      f3_R2 = NA_real_,
      number_of_significant_GO_terms = NA_real_,
      immune_Fisher_OR = NA_real_,
      immune_Fisher_p_value = NA_real_,
      recommended_paper_use = .$recommended_paper_use
    )
  }) %>%
  ungroup()
data.table::fwrite(decision_table_core, "revised_model_decision_table.tsv", sep = "\t")

model_count_plot_dt <- revised_selection_tier_counts %>%
  filter(tier %in% c(
    "GMM/RF", "Old best continuous 3-of-3", "Continuous-distance 3-class",
    "GMM/RF + continuous-distance agreement", "Distance margin >= 0.50",
    "Top 1%", "Top 500"
  ))

revised_class_palette <- c(
  "Neutral / equilibrium" = "grey60",
  "Balancing selection" = "#1f77b4",
  "Geographic sweep" = "#d62728",
  "Ambiguous" = "#7f7f7f",
  "No consensus" = "grey85"
)

p_model_count_revised <- ggplot2::ggplot(
  model_count_plot_dt,
  ggplot2::aes(x = tier, y = n_windows, fill = class)
) +
  ggplot2::geom_col(color = "grey25", linewidth = 0.2) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = revised_class_palette, drop = FALSE) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "bottom") +
  ggplot2::labs(x = NULL, y = "Windows", fill = "Class", title = "Revised model counts across selection tiers")
ggsave("revised_model_count_barplot.pdf", p_model_count_revised, width = 9.5, height = 6.2)
ggsave("revised_model_count_barplot.png", p_model_count_revised, width = 9.5, height = 6.2, dpi = 300)

dist_diag_dt <- agree_dat %>%
  select(
    d_neutral_empirical, d_balancing, d_sweep, continuous_distance_margin_empirical,
    continuous_distance_label_empirical_h
  ) %>%
  pivot_longer(
    cols = c(d_neutral_empirical, d_balancing, d_sweep, continuous_distance_margin_empirical),
    names_to = "metric",
    values_to = "value"
  )
p_distance_diag <- ggplot2::ggplot(
  dist_diag_dt,
  ggplot2::aes(x = value, fill = continuous_distance_label_empirical_h)
) +
  ggplot2::geom_histogram(bins = 50, alpha = 0.65, position = "identity") +
  ggplot2::facet_wrap(~ metric, scales = "free") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::labs(
    x = "Value",
    y = "Windows",
    fill = "Empirical distance label",
    title = "Continuous distance-to-centroid diagnostics"
  )
ggsave("continuous_distance_diagnostics.pdf", p_distance_diag, width = 11, height = 7)
ggsave("continuous_distance_diagnostics.png", p_distance_diag, width = 11, height = 7, dpi = 300)

############################################################
## [14H] GENOME TRACK: triple agreement tick plot (UPDATED: uses Chromosome_Names)
############################################################

drop_neutral <- TRUE
neutral_label <- "Neutral / equilibrium"

tick_half  <- 0.35
tick_alpha <- 0.90
tick_lwd   <- 0.5

## ---- chromosome order + lengths based on Chromosome_Names ----
chrom_df <- sd3c %>%
  dplyr::distinct(Chromosome_Names, Chromosome_Size) %>%
  dplyr::mutate(
    Chromosome_Names = as.character(Chromosome_Names),
    is_ZW = grepl("ZW", Chromosome_Names, ignore.case = TRUE),
    chr_index = dplyr::if_else(
      is_ZW, 99L,
      suppressWarnings(as.integer(gsub("Chromosome_", "", Chromosome_Names)))
    )
  ) %>%
  dplyr::arrange(chr_index) %>%
  dplyr::mutate(
    Chromosome_Names = factor(Chromosome_Names, levels = Chromosome_Names),
    chr_len_Mb = Chromosome_Size / 1e6,
    y = as.numeric(Chromosome_Names)
  )

chrom_lines <- chrom_df %>%
  dplyr::transmute(y, x = 0, xend = chr_len_Mb)

## ---- triple-agreement windows ----
triple_df <- sd3c %>%
  dplyr::select(Chromosome_Names, .start_bp, .end_bp, gmm_model_raw) %>%
  dplyr::left_join(agree_dat %>% dplyr::select(Chromosome_Names, .start_bp, triple_match),
                   by = c("Chromosome_Names", ".start_bp")) %>%
  dplyr::mutate(
    Chromosome_Names = as.character(Chromosome_Names),
    .start_bp = as.numeric(.start_bp),
    .end_bp   = ifelse(is.finite(.end_bp), as.numeric(.end_bp), NA_real_),
    mid_Mb    = ifelse(is.finite(.end_bp),
                       (.start_bp + .end_bp)/2/1e6,
                       .start_bp/1e6),
    triple_match = as.logical(triple_match),
    label = as.character(gmm_model_raw)
  ) %>%
  dplyr::filter(triple_match, !is.na(Chromosome_Names), is.finite(mid_Mb), !is.na(label)) %>%
  dplyr::mutate(
    Chromosome_Names = factor(Chromosome_Names, levels = levels(chrom_df$Chromosome_Names))
  ) %>%
  dplyr::left_join(chrom_df %>% dplyr::select(Chromosome_Names, y), by = "Chromosome_Names")

if (isTRUE(drop_neutral)) {
  triple_df <- triple_df %>% dplyr::filter(label != neutral_label)
}

class_colors_all <- c(
  "Balancing selection"     = "#1f77b4",
  "Geographic sweep"        = "#d62728",
  "Divergence with gene flow" = "#ff7f0e",
  "Allopatry-like"          = "#9467bd",
  "Neutral / equilibrium"   = "grey60"
)
class_colors <- class_colors_all
class_colors <- class_colors[names(class_colors) %in% unique(triple_df$label)]

p_triple <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = chrom_lines,
    ggplot2::aes(x = x, xend = xend, y = y, yend = y),
    linewidth = 2.2,
    color = "grey90",
    lineend = "butt"
  ) +
  ggplot2::geom_linerange(
    data = triple_df,
    ggplot2::aes(
      x = mid_Mb,
      ymin = y - tick_half,
      ymax = y + tick_half,
      color = label
    ),
    alpha = tick_alpha,
    linewidth = tick_lwd
  ) +
  ggplot2::scale_y_reverse(
    breaks = chrom_df$y,
    labels = as.character(chrom_df$Chromosome_Names),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_color_manual(values = class_colors, name = "Triple-agreed class\n(neutral excluded)") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = "bottom"
  ) +
  ggplot2::labs(
    x = "Position (Mb)",
    y = "Chromosome",
    title = "Genome-wide triple-agreement windows",
    subtitle = "Neutral / equilibrium is absent here because the continuous best-label agreement does not include a neutral class"
  )

show_plot(p_triple)
ggsave("Chromosome_triple_agreement_models.pdf", p_triple, width = 12, height = 7)
ggsave("Chromosome_triple_agreement_models.png", p_triple, width = 12, height = 7, dpi = 300)

############################################################
## [14I] GENOME TRACKS: tiered GMM/RF/continuous agreement
############################################################

tier_order <- c(
  "GMM/RF",
  "Best continuous",
  "Margin >= 0.25",
  "Margin >= 0.50",
  "Margin >= 1.00",
  "Top 1%",
  "Top 500"
)

tier_long <- sd3c %>%
  dplyr::select(Chromosome_Names, .start_bp, .end_bp, gmm_model_raw) %>%
  dplyr::left_join(
    agree_dat %>%
      dplyr::select(
        Chromosome_Names, .start_bp, rf_label,
        gmm_rf_match, triple_match, triple_margin_025,
        triple_margin_050, triple_margin_100,
        triple_top1pct, triple_tail500
      ),
    by = c("Chromosome_Names", ".start_bp")
  ) %>%
  dplyr::mutate(
    Chromosome_Names = as.character(Chromosome_Names),
    .start_bp = as.numeric(.start_bp),
    .end_bp = ifelse(is.finite(.end_bp), as.numeric(.end_bp), NA_real_),
    mid_Mb = ifelse(is.finite(.end_bp), (.start_bp + .end_bp) / 2 / 1e6, .start_bp / 1e6),
    label = as.character(rf_label)
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      gmm_rf_match, triple_match, triple_margin_025,
      triple_margin_050, triple_margin_100,
      triple_top1pct, triple_tail500
    ),
    names_to = "tier_key",
    values_to = "is_selected"
  ) %>%
  dplyr::mutate(
    is_selected = as.logical(is_selected),
    tier = dplyr::recode(
      tier_key,
      gmm_rf_match = "GMM/RF",
      triple_match = "Best continuous",
      triple_margin_025 = "Margin >= 0.25",
      triple_margin_050 = "Margin >= 0.50",
      triple_margin_100 = "Margin >= 1.00",
      triple_top1pct = "Top 1%",
      triple_tail500 = "Top 500"
    ),
    tier = factor(tier, levels = tier_order),
    Chromosome_Names = factor(Chromosome_Names, levels = levels(chrom_df$Chromosome_Names))
  ) %>%
  dplyr::filter(is_selected, !is.na(Chromosome_Names), is.finite(mid_Mb), !is.na(label)) %>%
  dplyr::left_join(chrom_df %>% dplyr::select(Chromosome_Names, y), by = "Chromosome_Names")

tier_class_counts <- tier_long %>%
  dplyr::count(tier, label, name = "n_windows") %>%
  dplyr::group_by(tier) %>%
  dplyr::mutate(frac_tier = n_windows / sum(n_windows)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(tier, dplyr::desc(n_windows))
data.table::fwrite(tier_class_counts, "Agreement_Tier_ClassCounts_Global.tsv", sep = "\t")

plot_tier_chromosomes <- function(dat, out_stub, title, include_neutral = TRUE, subtitle = NULL) {
  dat_plot <- dat
  if (!include_neutral) {
    dat_plot <- dat_plot %>% dplyr::filter(label != neutral_label)
  }
  colors_use <- class_colors_all[names(class_colors_all) %in% unique(dat_plot$label)]

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = chrom_lines,
      ggplot2::aes(x = x, xend = xend, y = y, yend = y),
      linewidth = 1.6,
      color = "grey90",
      lineend = "butt"
    ) +
    ggplot2::geom_linerange(
      data = dat_plot,
      ggplot2::aes(
        x = mid_Mb,
        ymin = y - 0.28,
        ymax = y + 0.28,
        color = label
      ),
      alpha = 0.65,
      linewidth = 0.35
    ) +
    ggplot2::facet_wrap(~ tier, ncol = 1) +
    ggplot2::scale_y_reverse(
      breaks = chrom_df$y,
      labels = as.character(chrom_df$Chromosome_Names),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_color_manual(
      values = colors_use,
      name = if (isTRUE(include_neutral)) "Class" else "Class\n(neutral excluded)"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 7)
    ) +
    ggplot2::labs(
      x = "Position (Mb)",
      y = "Chromosome",
      title = title,
      subtitle = subtitle
    )

  show_plot(p)
  ggsave(paste0(out_stub, ".pdf"), p, width = 12.5, height = 18)
  ggsave(paste0(out_stub, ".png"), p, width = 12.5, height = 18, dpi = 300)
  invisible(p)
}

plot_tier_chromosomes(
  tier_long,
  out_stub = "Chromosome_agreement_tiers_all_models",
  title = "Genome-wide agreement tiers (including neutral)",
  include_neutral = TRUE,
  subtitle = "Neutral / equilibrium is shown only in the GMM/RF tier"
)
plot_tier_chromosomes(
  tier_long,
  out_stub = "Chromosome_agreement_tiers_nonneutral_models",
  title = "Genome-wide agreement tiers (neutral excluded)",
  include_neutral = FALSE,
  subtitle = "Best continuous, margin, Top 1%, and Top 500 tiers contain only sweep or balancing by construction"
)

p_tier_counts <- tier_class_counts %>%
  ggplot2::ggplot(ggplot2::aes(x = tier, y = n_windows, fill = label)) +
  ggplot2::geom_col(color = "grey25", linewidth = 0.2) +
  ggplot2::scale_fill_manual(values = class_colors_all, name = "Class") +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::coord_flip() +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                 legend.position = "bottom") +
  ggplot2::labs(x = NULL, y = "Windows", title = "Agreement-tier window counts by class")
show_plot(p_tier_counts)
ggsave("Agreement_Tier_ClassCounts.pdf", p_tier_counts, width = 9, height = 5.8)
ggsave("Agreement_Tier_ClassCounts.png", p_tier_counts, width = 9, height = 5.8, dpi = 300)

p_tier_frac <- tier_class_counts %>%
  ggplot2::ggplot(ggplot2::aes(x = tier, y = frac_tier, fill = label)) +
  ggplot2::geom_col(color = "grey25", linewidth = 0.2) +
  ggplot2::scale_fill_manual(values = class_colors_all, name = "Class") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::coord_flip() +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                 legend.position = "bottom") +
  ggplot2::labs(x = NULL, y = "Fraction within tier", title = "Agreement-tier class composition")
show_plot(p_tier_frac)
ggsave("Agreement_Tier_ClassFractions.pdf", p_tier_frac, width = 9, height = 5.8)
ggsave("Agreement_Tier_ClassFractions.png", p_tier_frac, width = 9, height = 5.8, dpi = 300)

############################################################
## Join agreement labels back onto your raw table (UPDATED key handling)
############################################################
## Your data now has:
##   Chromosomes (scaffold IDs) AND Chromosome_Names (pretty chromosome names).
## Agreement tables are keyed by Chromosome_Names + .start_bp.
## This join uses Chromosome_Names + window_start.
############################################################

raw_chr_col <- .find_col(sd3_new_af, c("Chromosome_Names","Chromosomes","Chromosome","chr"))
raw_start_col <- .find_col(sd3_new_af, c("window_start",".start_bp","start_bp","start"))
if (is.na(raw_chr_col) || is.na(raw_start_col)) {
  stop("Could not find Chromosome/start columns in sd3_new_af for joining agreement labels.")
}

sd3_new_af2 <- sd3_new_af

## --- rename chromosome column only if needed ---
if (raw_chr_col != "Chromosome_Names") {
  sd3_new_af2 <- sd3_new_af2 %>%
    dplyr::rename(Chromosome_Names = !!raw_chr_col)
}

## --- create .start_bp safely (avoid duplicate names) ---
if (raw_start_col != ".start_bp") {

  ## if .start_bp already exists, remove it first
  if (".start_bp" %in% names(sd3_new_af2)) {
    sd3_new_af2 <- sd3_new_af2 %>% dplyr::select(-.start_bp)
  }

  sd3_new_af2 <- sd3_new_af2 %>%
    dplyr::rename(.start_bp = !!raw_start_col)
}

## --- enforce types ---
sd3_new_af2 <- sd3_new_af2 %>%
  mutate(
    Chromosome_Names = as.character(Chromosome_Names),
    .start_bp = as.numeric(.start_bp)
  )

## --- final join ---
sd3_new_af2 <- sd3_new_af2 %>%
  left_join(
    agree_dat %>%
      select(
        Chromosome_Names, .start_bp,
        gmm_rf_match, triple_match, triple_tail500, triple_top1pct,
        triple_margin_025, triple_margin_050, triple_margin_100,
        contdist_match, dist_margin_025, dist_margin_050, dist_margin_100,
        score_margin_025, score_margin_050, score_margin_100,
        dist_top1pct, dist_top500,
        cont_label_h, cont_label_tail500_h, cont_label_top1pct_h,
        continuous_distance_label_origin_h, continuous_distance_label_empirical_h,
        continuous_distance_margin_origin, continuous_distance_margin_empirical,
        continuous_margin_label_025, continuous_margin_label_050, continuous_margin_label_100,
        d_neutral_origin, d_neutral_empirical, d_balancing, d_sweep,
        continuous_score_margin, dplyr::any_of("rf_pmax"),
        cont_margin, cont_score_quantile
      ),
    by = c("Chromosome_Names", ".start_bp")
  )

cat("\n[JOIN] Added revised agreement columns to sd3_new_af2: GMM/RF, old triple, explicit-neutral distance tiers, score-margin tiers, strict tails, and continuous diagnostics\n")
print(table(sd3_new_af2$triple_match, useNA = "ifany"))

save(agree_dat, sd3_new_af2, file = "barrier_label_objects.rdata")
cat("\n[JOIN] Saved barrier label objects for downstream barrier-loci repo: barrier_label_objects.rdata\n")

handoff_manifest <- tibble::tibble(
  artifact = c(
    "barrier_label_objects.rdata",
    "Continuous_Agreement_Windows_Global.tsv",
    "Continuous_Agreement_Tiers_Global.tsv",
    "Continuous_Agreement_ByClass_Global.tsv",
    "Agreement_Tier_ClassCounts_Global.tsv",
    "revised_selection_tier_counts.tsv",
    "revised_selection_tier_window_assignments.tsv",
    "continuous_distance_centroids.tsv",
    "continuous_distance_window_scores.tsv",
    "continuous_margin_rejection_labels.tsv",
    "revised_model_decision_table.tsv",
    "Chromosome_agreement_tiers_all_models.pdf",
    "Chromosome_agreement_tiers_nonneutral_models.pdf"
  ),
  description = c(
    "RData containing agree_dat and sd3_new_af2 for barrier-loci downstream scripts",
    "Window-level GMM/RF/continuous agreement flags used for downstream tiers",
    "Global counts for agreement tiers",
    "Agreement-tier counts by GMM/RF class",
    "Plot-ready agreement-tier class counts",
    "Revised explicit-neutral tier summaries with window counts, bp, and metric summaries",
    "Window-level assignments for revised tiers",
    "Explicit-neutral continuous-distance centroids in 6D z-space",
    "Window-level continuous-distance scores and labels",
    "Continuous score-margin rejection labels",
    "Decision table for comparing revised selection tiers",
    "Genome-wide chromosome plot for all agreement tiers including neutral",
    "Genome-wide chromosome plot for non-neutral agreement tiers"
  ),
  exists = file.exists(artifact),
  generated_at = as.character(Sys.time()),
  source_script = "master_pipeline_UPDATED_v2026-02-07-02_ranknormRec_QCfix_Fig17.R"
)
data.table::fwrite(handoff_manifest, "Selection_Barrier_Handoff_Manifest.tsv", sep = "\t")
cat("\n[JOIN] Wrote downstream handoff manifest: Selection_Barrier_Handoff_Manifest.tsv\n")


############################################################
## 


############################################################
## [15] BALANCING vs SWEEP TESTS (triple-agreement) with Ho/He + recombination
##      + Wahlund-effect diagnostics
##
## INPUT: sd3_new_af2  (raw table joined with agreement labels)
##
## Requires columns (typical): east_Ho, west_Ho, admixed_Ho,
##                            east_He, west_He, admixed_He,
##                            total_Ho, total_He, total_n_sites,
##                            delta_Fhet_Wahlund, Wahlund_supported,
##                            recombination
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

keep_models <- c("Balancing selection", "Geographic sweep")
require_recombination <- TRUE

if (!exists("sd3_new_af2")) {
  cat("\n[15] NOTE: sd3_new_af2 not found; skipping Wahlund/Balancing-vs-Sweep module.\n")
} else {

  cat("\n================ [15] BALANCING vs SWEEP (TRIPLE AGREEMENT) =================\n")

  cat("\ntriple_match class:\n")
  print(class(sd3_new_af2$triple_match))

  cat("\ntriple_match table (incl NA):\n")
  print(table(sd3_new_af2$triple_match, useNA = "ifany"))

  cat("\nUnique cont_label_h values (first 50):\n")
  print(head(sort(unique(sd3_new_af2$cont_label_h)), 50))

  sd3_new_af2_fix <- sd3_new_af2 %>%
    mutate(
      ## robust logical coercion
      triple_match2 = case_when(
        triple_match %in% c(TRUE, 1, "1", "TRUE", "T", "true", "t") ~ TRUE,
        triple_match %in% c(FALSE, 0, "0", "FALSE", "F", "false", "f") ~ FALSE,
        TRUE ~ FALSE
      ),
      ## robust label cleaning
      cont_label_h2 = str_trim(as.character(cont_label_h))
    )

  cat("\nAfter coercion: triple_match2 table:\n")
  print(table(sd3_new_af2_fix$triple_match2, useNA = "ifany"))

  cat("\nUnique cont_label_h2 values among triple_match2==TRUE (first 50):\n")
  print(head(sort(unique(sd3_new_af2_fix$cont_label_h2[sd3_new_af2_fix$triple_match2])), 50))

  ## Column presence check (fail gracefully with a helpful message)
  needed_cols <- c("east_Ho","west_Ho","admixed_Ho","east_He","west_He","admixed_He",
                   "total_Ho","total_He","total_n_sites","delta_Fhet_Wahlund","Wahlund_supported")
  miss_cols <- setdiff(needed_cols, names(sd3_new_af2_fix))
  if (length(miss_cols) > 0) {
    cat("\n[15] SKIP: Missing required columns in sd3_new_af2 for Wahlund module:\n")
    print(miss_cols)
  } else {

    ## ---- Filter to triple-agreement + target classes; compute derived metrics ----
    dat <- sd3_new_af2_fix %>%
      filter(triple_match2) %>%
      filter(cont_label_h2 %in% keep_models) %>%
      mutate(
        model = factor(cont_label_h2, levels = c("Geographic sweep", "Balancing selection")),

        ## rank-normalized recombination (kept separate from raw)
        z_rec = ranknorm(recombination),

        ## within-pop means (core “true balancing” signal)
        mean_within_Ho = rowMeans(cbind(east_Ho, west_Ho, admixed_Ho), na.rm = TRUE),
        mean_within_He = rowMeans(cbind(east_He, west_He, admixed_He), na.rm = TRUE),

        ## pooling inflation diagnostics (Wahlund-style)
        excess_pooling_Ho = total_Ho - mean_within_Ho,
        excess_pooling_He = total_He - mean_within_He,

        ## split balancing into Wahlund vs non-Wahlund
        Wahlund_class = case_when(
          model == "Balancing selection" & isTRUE(Wahlund_supported) ~ "Balancing + Wahlund",
          model == "Balancing selection" & !isTRUE(Wahlund_supported) ~ "Balancing (non-Wahlund)",
          model == "Geographic sweep" ~ "Sweep",
          TRUE ~ NA_character_
        ),
        Wahlund_class = factor(
          Wahlund_class,
          levels = c("Sweep", "Balancing (non-Wahlund)", "Balancing + Wahlund")
        )
      )

    cat("\nN windows (triple & model-filtered):", nrow(dat), "\n")
    print(table(dat$model, useNA = "ifany"))
    cat("\nWahlund_class counts:\n")
    print(table(dat$Wahlund_class, useNA = "ifany"))

    ## fail loudly if something breaks upstream again
    stopifnot(nrow(dat) > 0)
    stopifnot(length(unique(dat$model)) == 2)

    if (isTRUE(require_recombination) && !"recombination" %in% names(dat)) {
      stop("recombination column not found in dat (set require_recombination <- FALSE to bypass).")
    }

    ## ---- Summary diagnostics tables ----
    diag_model <- dat %>%
      group_by(model) %>%
      summarise(
        n = n(),
        mean_within_Ho = mean(mean_within_Ho, na.rm = TRUE),
        mean_total_Ho  = mean(total_Ho, na.rm = TRUE),
        mean_within_He = mean(mean_within_He, na.rm = TRUE),
        mean_total_He  = mean(total_He, na.rm = TRUE),
        mean_excess_pooling_Ho = mean(excess_pooling_Ho, na.rm = TRUE),
        mean_excess_pooling_He = mean(excess_pooling_He, na.rm = TRUE),
        mean_delta_Wahlund = mean(delta_Fhet_Wahlund, na.rm = TRUE),
        frac_Wahlund_supported = mean(Wahlund_supported, na.rm = TRUE),
        mean_z_rec = if ("recombination" %in% names(dat)) mean(z_rec, na.rm = TRUE) else NA_real_,
        median_z_rec = if ("recombination" %in% names(dat)) median(z_rec, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )

    diag_wahlund <- dat %>%
      group_by(Wahlund_class) %>%
      summarise(
        n = n(),
        mean_within_Ho = mean(mean_within_Ho, na.rm = TRUE),
        mean_total_Ho  = mean(total_Ho, na.rm = TRUE),
        mean_within_He = mean(mean_within_He, na.rm = TRUE),
        mean_total_He  = mean(total_He, na.rm = TRUE),
        mean_excess_pooling_Ho = mean(excess_pooling_Ho, na.rm = TRUE),
        mean_excess_pooling_He = mean(excess_pooling_He, na.rm = TRUE),
        mean_delta_Wahlund = mean(delta_Fhet_Wahlund, na.rm = TRUE),
        frac_Wahlund_supported = mean(Wahlund_supported, na.rm = TRUE),
        mean_z_rec = if ("recombination" %in% names(dat)) mean(z_rec, na.rm = TRUE) else NA_real_,
        median_z_rec = if ("recombination" %in% names(dat)) median(z_rec, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )

    cat("\n=== Model-level diagnostics ===\n")
    print(diag_model)
    cat("\n=== Wahlund-class diagnostics ===\n")
    print(diag_wahlund)

    ## ---- Key plots (Ho/He within + pooling inflation) ----
    p_within_Ho <- ggplot(dat, aes(x = model, y = mean_within_Ho, fill = model)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.15, outlier.shape = NA) +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(x = NULL, y = "Mean within-population Ho",
           title = "Within-population Ho by model (triple agreement)")

    p_within_He <- ggplot(dat, aes(x = model, y = mean_within_He, fill = model)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.15, outlier.shape = NA) +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(x = NULL, y = "Mean within-population He",
           title = "Within-population He by model (triple agreement)")

    p_poolHo <- ggplot(dat, aes(x = model, y = excess_pooling_Ho, fill = model)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.15, outlier.shape = NA) +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(x = NULL, y = "total_Ho − mean_within_Ho",
           title = "Pooling inflation of Ho (Wahlund-style) by model")

    p_poolHe <- ggplot(dat, aes(x = model, y = excess_pooling_He, fill = model)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.15, outlier.shape = NA) +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "none") +
      labs(x = NULL, y = "total_He − mean_within_He",
           title = "Pooling inflation of He (Wahlund-style) by model")

    show_plot(p_within_Ho); show_plot(p_within_He); show_plot(p_poolHo); show_plot(p_poolHe)

    ## ---- Wilcoxon tests (Balancing vs Sweep) ----
    cat("\n=== Wilcoxon tests (Balancing vs Sweep) ===\n")
    print(wilcox.test(mean_within_Ho ~ model, data = dat))
    print(wilcox.test(mean_within_He ~ model, data = dat))
    print(wilcox.test(excess_pooling_Ho ~ model, data = dat))
    print(wilcox.test(excess_pooling_He ~ model, data = dat))
    print(wilcox.test(delta_Fhet_Wahlund ~ model, data = dat))

    ## ---- Recombination: distribution + logistic prediction of Balancing vs Sweep ----
    if ("recombination" %in% names(dat)) {

      dat_rec <- dat %>% mutate(z_rec = ranknorm(recombination)) %>% filter(is.finite(z_rec))

      p_rec <- ggplot(dat_rec, aes(x = model, y = z_rec, fill = model)) +
        geom_violin(trim = FALSE, alpha = 0.6) +
        geom_boxplot(width = 0.15, outlier.shape = NA) +
        theme_bw() +
        theme(panel.grid = element_blank(), legend.position = "none") +
        labs(x = NULL, y = "Recombination (rank-normalized z)",
             title = "Recombination (rank-normalized) by model (triple agreement)")
      show_plot(p_rec)

      cat("\n=== Wilcoxon: rank-normalized recombination (z_rec) ~ model ===\n")
      print(wilcox.test(z_rec ~ model, data = dat_rec))

      ## Logistic: Balancing=1 vs Sweep=0
      dat_rec <- dat_rec %>% mutate(is_balancing = as.integer(model == "Balancing selection"))

      m_rec1 <- glm(is_balancing ~ z_rec, data = dat_rec, family = binomial())
      m_rec2 <- glm(is_balancing ~ z_rec + total_n_sites + delta_Fhet_Wahlund,
                    data = dat_rec, family = binomial())
      m_rec3 <- glm(is_balancing ~ z_rec + total_n_sites + delta_Fhet_Wahlund + excess_pooling_Ho,
                    data = dat_rec, family = binomial())

      cat("\n=== Logistic models predicting Balancing vs Sweep ===\n")
      cat("\n--- m_rec1: is_balancing ~ z_rec ---\n"); print(summary(m_rec1))
      cat("\n--- m_rec2: + total_n_sites + delta_Fhet_Wahlund ---\n"); print(summary(m_rec2))
      cat("\n--- m_rec3: + excess_pooling_Ho ---\n"); print(summary(m_rec3))

      ## Odds ratio for recombination in m_rec2
      OR <- exp(coef(m_rec2)["z_rec"])
      CI <- exp(suppressMessages(confint(m_rec2))["z_rec", ])
      cat("\nOdds ratio for +1 SD rank-normalized recombination (z_rec) (m_rec2):", OR, "\n")
      cat("95% CI:", CI[1], "-", CI[2], "\n")
    }

    ## ---- Linear models: do balancing effects persist after controls? ----
    m_totalHo <- lm(
      total_Ho ~ model + mean_within_Ho + delta_Fhet_Wahlund + total_n_sites +
        (if ("recombination" %in% names(dat)) z_rec else 0),
      data = dat
    )
    cat("\n=== Linear model: total_Ho controlled ===\n")
    print(summary(m_totalHo))

    m_totalHe <- lm(
      total_He ~ model + mean_within_He + delta_Fhet_Wahlund + total_n_sites +
        (if ("recombination" %in% names(dat)) z_rec else 0),
      data = dat
    )
    cat("\n=== Linear model: total_He controlled ===\n")
    print(summary(m_totalHe))

    ## Optional exports:
    # write.csv(diag_model,   "diag_model_bal_vs_sweep.csv",   row.names = FALSE)
    # write.csv(diag_wahlund, "diag_wahlund_bal_vs_sweep.csv", row.names = FALSE)
  }
}

############################################################
## [16] MODEL VALIDATION: Sweeps vs Balancing Selection (NO 0.9 cutoff)
## Uses triple-agreement windows; compares sweep vs balancing
## with Ho/He, Wahlund metrics, recombination, and mean_abs_delta_p (continuous)
## Composite score = ORTHOGONAL: z(|Δp|) − z(He)  (NO recombination inside)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

if (!exists("sd3_new_af2")) {
  cat("\n[16] NOTE: sd3_new_af2 not found; skipping validation module.\n")
} else if (!all(c("mean_abs_delta_p","recombination","east_He","west_He","admixed_He",
                  "east_Ho","west_Ho","admixed_Ho","total_He","total_n_sites") %in% names(sd3_new_af2))) {
  cat("\n[16] NOTE: Required columns missing for validation module; skipping.\n")
} else {

  model_colors <- c(
    "Geographic sweep"    = "#F8766D",
    "Balancing selection" = "#00BFC4"
  )

  dat <- sd3_new_af2 %>%
    filter(triple_match) %>%
    filter(cont_label_h %in% c("Geographic sweep", "Balancing selection")) %>%
    mutate(
      model = factor(cont_label_h, levels = c("Geographic sweep", "Balancing selection")),

      mean_within_Ho = rowMeans(cbind(east_Ho, west_Ho, admixed_Ho), na.rm = TRUE),
      mean_within_He = rowMeans(cbind(east_He, west_He, admixed_He), na.rm = TRUE),

      ## “Pooling effect” (positive = within > pooled, pooled depressed)
      excess_pooling_Ho = mean_within_Ho - total_Ho,
      excess_pooling_He = mean_within_He - total_He,
      z_rec = ranknorm(recombination)
    ) %>%
    filter(
      is.finite(mean_abs_delta_p),
      is.finite(z_rec),
      is.finite(mean_within_He),
      is.finite(mean_within_Ho),
      is.finite(total_He),
      is.finite(total_n_sites)
    )

  cat("\n================ [16] VALIDATION DATASET =================\n")
  cat("N windows:", nrow(dat), "\n")
  print(table(dat$model, useNA = "ifany"))

  if (nrow(dat) > 0 && length(unique(dat$model)) == 2) {

    dat_score <- dat %>%
      mutate(
        z_delta_p = as.numeric(scale(mean_abs_delta_p)),
        z_He      = as.numeric(scale(mean_within_He)),
        sweep_support_score = z_delta_p - z_He
      )

    cat("\n================ [16] NONPARAMETRIC TESTS (NO Δp cutoff) =================\n")
    cat("\n[1] mean_abs_delta_p (expect: Sweep > Balancing)\n")
    print(wilcox.test(mean_abs_delta_p ~ model, data = dat_score))

    cat("\n[2] recombination (expect: Sweep < Balancing)\n")
    print(wilcox.test(z_rec ~ model, data = dat_score))

    cat("\n[3] mean within-population He (expect: Balancing > Sweep)\n")
    print(wilcox.test(mean_within_He ~ model, data = dat_score))

    cat("\n[4] mean within-population Ho (expect: Balancing > Sweep)\n")
    print(wilcox.test(mean_within_Ho ~ model, data = dat_score))

    cat("\n[5] pooling depression of Ho = (within − pooled) (expect: Sweep > Balancing)\n")
    print(wilcox.test(excess_pooling_Ho ~ model, data = dat_score))

    cat("\n[6] Delta Fhet Wahlund (expect: Sweep > Balancing)\n")
    if ("delta_Fhet_Wahlund" %in% names(dat_score)) {
      print(wilcox.test(delta_Fhet_Wahlund ~ model, data = dat_score))
    } else {
      cat("  (delta_Fhet_Wahlund not present; skipped)\n")
    }

    cat("\n[7] ORTHOGONAL sweep-support score = z(|Δp|) − z(He) (expect: Sweep > Balancing)\n")
    print(wilcox.test(sweep_support_score ~ model, data = dat_score))

    cat("\n================ [16] REGRESSION: mean_abs_delta_p ~ model * recombination =================\n")
    lm_dp <- lm(mean_abs_delta_p ~ model * z_rec, data = dat_score)
    print(summary(lm_dp))

    cat("\n================ [16] REGRESSION: ORTHOGONAL score ~ model * recombination =================\n")
    lm_score <- lm(sweep_support_score ~ model * z_rec, data = dat_score)
    print(summary(lm_score))

    cat("\n================ [16] LOGISTIC: Predicting Sweep (NO Δp cutoff) =================\n")
    dat_score$IsSweep <- as.integer(dat_score$model == "Geographic sweep")
    m1 <- glm(IsSweep ~ mean_abs_delta_p, data = dat_score, family = binomial)
    m2 <- glm(IsSweep ~ mean_abs_delta_p + z_rec + total_n_sites, data = dat_score, family = binomial)
    m3 <- glm(IsSweep ~ mean_abs_delta_p * z_rec + total_n_sites, data = dat_score, family = binomial)

    cat("\n[Model 1] Sweep ~ |Δp|\n"); print(summary(m1))
    cat("\n[Model 2] Sweep ~ |Δp| + z_rec + total_n_sites\n"); print(summary(m2))
    cat("\n[Model 3] Sweep ~ |Δp| × z_rec + total_n_sites\n"); print(summary(m3))

    ## Plots
    p_score <- ggplot(dat_score, aes(model, sweep_support_score, fill = model)) +
      geom_violin(trim = FALSE, alpha = 0.8, color = "grey25") +
      geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
      scale_fill_manual(values = model_colors) +
      theme_bw() +
      labs(
        title = "Orthogonal sweep-support score by model",
        subtitle = "Score = z(|Δp|) − z(He) (recombination excluded)",
        y = "Orthogonal sweep-support score",
        x = NULL
      ) +
      theme(legend.position = "none")
    show_plot(p_score)

    p_dp_rec <- ggplot(dat_score, aes(z_rec, mean_abs_delta_p, color = model)) +
      geom_point(alpha = 0.45, size = 1.6) +
      geom_smooth(se = FALSE, linewidth = 1) +
      scale_color_manual(values = model_colors) +
      theme_bw() +
      labs(
        title = "Parental differentiation vs rank-normalized recombination by model",
        subtitle = "No thresholding; raw mean |Δp|; x-axis uses rank-normalized recombination (z_rec)",
        x = "Recombination (rank-normalized z)",
        y = "mean |Δp|"
      )
    show_plot(p_dp_rec)

    p_score_rec <- ggplot(dat_score, aes(z_rec, sweep_support_score, color = model)) +
      geom_point(alpha = 0.35, size = 1.5) +
      geom_smooth(se = FALSE, linewidth = 1) +
      scale_color_manual(values = model_colors) +
      theme_bw() +
      labs(
        title = "Orthogonal sweep-support score vs rank-normalized recombination",
        subtitle = "Score excludes recombination: z(|Δp|) − z(He)",
        x = "Recombination (rank-normalized z)",
        y = "Orthogonal sweep-support score"
      )
    show_plot(p_score_rec)

  } else {
    cat("\n[16] NOTE: Validation dataset empty or only one class; skipping.\n")
  }
}


# END-OF-PIPELINE QC HARNESS (updated for 5D clustering)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

cat("\n==================== RUNNING QC HARNESS ====================\n")

QC_RESULTS <- list()

# 1) z-score summaries
zcols <- c("z_dxy","z_fst","z_pi","z_f3","z_taj","z_rec")
QC_RESULTS$z_summary <- summary(sd3t %>% select(all_of(zcols)))

QC_RESULTS$z_sd <- sd3_clusters %>%
  summarise(across(all_of(zcols), ~sd(.x, na.rm = TRUE)))

QC_RESULTS$z_sd_flagged <- {
  bad <- names(QC_RESULTS$z_sd)[abs(unlist(QC_RESULTS$z_sd) - 1) > 0.05]
  if (length(bad) == 0) "OK (all SD ≈ 1)" else bad
}

# 2) p_cluster diagnostics
valid5 <- with(sd3_clusters,
               is.finite(z_dxy) & is.finite(z_fst) &
                 is.finite(z_pi)  & is.finite(z_f3)  &
                 is.finite(z_taj))

QC_RESULTS$p_cluster_range <- range(sd3_clusters$p_cluster[valid5], na.rm = TRUE)
QC_RESULTS$p_cluster_summary <- summary(sd3_clusters$p_cluster[valid5])
QC_RESULTS$cluster_counts <- table(sd3_clusters$cluster, useNA = "ifany")

# 3) Empirical centroids (z-space) + scores
QC_RESULTS$centroids_z5 <- sd3_clusters %>%
  filter(valid5, !is.na(cluster)) %>%
  group_by(cluster) %>%
  summarise(across(all_of(c("z_dxy","z_fst","z_pi","z_f3","z_taj")), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

QC_RESULTS$centroid_scores <- QC_RESULTS$centroids_z5 %>%
  mutate(
    admix = neg_tail(z_f3),
    sweep_score_c = ( +1.2*z_fst -1.0*z_dxy -1.2*z_pi -1.0*neg_tail(z_taj) -0.3*admix ),
    bal_score_c   = ( -1.0*z_fst +1.3*z_pi +1.2*pos_tail(z_taj) -0.3*abs(z_dxy) ),
    dgf_score_c   = ( +1.1*z_fst +1.1*z_dxy +0.8*admix -0.4*z_pi ),
    allop_score_c = ( +1.0*z_fst -1.0*z_pi -0.8*abs(z_dxy) -0.6*pos_tail(z_taj) -0.3*admix ),
    r_centroid    = sqrt(z_dxy^2 + z_fst^2 + z_pi^2 + z_f3^2 + z_taj^2)
  ) %>%
  arrange(r_centroid)

# 4) Continuous top-N checks
QC_RESULTS$model_cont_counts <- table(sd3_cont$model_cont, useNA = "ifany")

QC_RESULTS$topN_overlap <- sd3_cont %>%
  summarise(
    n_sweep_top = sum(is_sweep_top, na.rm = TRUE),
    n_bal_top   = sum(is_bal_top,   na.rm = TRUE),
    n_dgf_top   = sum(is_dgf_top,   na.rm = TRUE),
    n_allop_top = sum(is_allop_top, na.rm = TRUE),
    n_overlap_sweep_bal = sum(is_sweep_top & is_bal_top, na.rm = TRUE),
    n_overlap_any = sum(
      (is_sweep_top + is_bal_top + is_dgf_top + is_allop_top) > 1,
      na.rm = TRUE
    )
  )

# 5) Recombination-by-model sanity (GMM labels)
QC_RESULTS$recomb_by_gmm <- sd3_clusters %>%
  group_by(gmm_model_raw) %>%
  summarise(
    n = n(),
    mean_z_rec = mean(z_rec, na.rm = TRUE),
    med_z_rec  = median(z_rec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_z_rec)

# FINAL REPORT
cat("\n==================== QC SUMMARY ====================\n")
for (nm in names(QC_RESULTS)) {
  cat("\n----", nm, "----\n")
  print(QC_RESULTS[[nm]])
}
cat("\n==================== QC COMPLETE ====================\n")



 

############################################################
## [17] PAPER FIGURE PANEL SET: Sweeps vs Recombination + π residual
## Produces separate panel files (A–E) and (if patchwork available) a combined multi-panel PDF.
## Requires: sd3_new_af2 with at least recombination, pi_average, cont_label_h
############################################################

make_sweep_recomb_panels <- function(sd3_new_af2,
                                    out_prefix = "Figure_Sweep_Recomb",
                                    out_dir = ".",
                                    use_ranknorm = TRUE) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })

  if (!is.function(ranknorm)) {
    stop("ranknorm() not found in environment; please ensure it's defined earlier in the script.")
  }

  req_cols <- c("recombination", "pi_average", "cont_label_h")
  miss <- setdiff(req_cols, names(sd3_new_af2))
  if (length(miss) > 0) {
    stop("make_sweep_recomb_panels(): Missing required columns in sd3_new_af2: ",
         paste(miss, collapse = ", "))
  }

  # Build dat4
  dat4 <- sd3_new_af2 %>%
    mutate(
      z_rec   = ranknorm(recombination),
      z_pi    = ranknorm(pi_average),
      is_sweep = as.integer(cont_label_h == "Geographic sweep")
    ) %>%
    filter(is.finite(z_rec), is.finite(z_pi), !is.na(is_sweep))

  if (nrow(dat4) < 1000) {
    warning("make_sweep_recomb_panels(): dat4 has only ", nrow(dat4), " rows after filtering.")
  }

  # Residualize pi ~ recombination
  lm_pi <- lm(z_pi ~ z_rec, data = dat4)
  dat4$pi_resid <- resid(lm_pi)

  # Logistic model: sweep ~ pi_resid + recombination
  m_res <- glm(is_sweep ~ pi_resid + z_rec, data = dat4, family = binomial())

  # Odds ratios + 95% CI
  or_ci <- function(fit) {
    b  <- coef(fit)
    ci <- suppressMessages(confint(fit))   # profile CI
    tibble(
      term = names(b),
      estimate = b,
      OR = exp(b),
      lo = exp(ci[,1]),
      hi = exp(ci[,2])
    ) %>% filter(term != "(Intercept)")
  }
  or_tab <- or_ci(m_res) %>%
    mutate(term = recode(term,
                         "pi_resid" = "\u03C0 residual",
                         "z_rec"    = "Recombination (z_rec)"))

  # PANEL A: sweep fraction by recombination decile
  dec_tab <- dat4 %>%
    mutate(rec_decile = ntile(z_rec, 10)) %>%
    group_by(rec_decile) %>%
    summarise(
      n = n(),
      sweep_frac = mean(is_sweep, na.rm = TRUE),
      .groups = "drop"
    )

  pA <- ggplot(dec_tab, aes(x = rec_decile, y = sweep_frac)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      x = "Recombination decile (z_rec; 1 = lowest recombination)",
      y = "Fraction sweep",
      title = "A. Sweeps enriched in low recombination"
    )

  # PANEL B: partial effect P(sweep) vs z_rec holding pi_resid=0
  zseq <- seq(quantile(dat4$z_rec, 0.01, na.rm = TRUE),
              quantile(dat4$z_rec, 0.99, na.rm = TRUE),
              length.out = 250)

  predB <- tibble(
    z_rec = zseq,
    pi_resid = 0
  ) %>%
    mutate(p = predict(m_res, newdata = ., type = "response"))

  pB <- ggplot(predB, aes(x = z_rec, y = p)) +
    geom_line(linewidth = 0.8) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      x = "Recombination (z_rec; rank-normalized)",
      y = "P(Sweep)",
      title = "B. Partial effect of recombination (\u03C0 residual = 0)"
    )

  # PANEL C: pi_resid by sweep vs non-sweep
  dat4$SweepClass <- factor(dat4$is_sweep, levels = c(0,1), labels = c("Non-sweep", "Sweep"))

  pC <- ggplot(dat4, aes(x = SweepClass, y = pi_resid)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      x = NULL,
      y = "\u03C0 residual (after removing recombination effect)",
      title = "C. Sweeps show unusually low \u03C0 beyond recombination"
    )

  # PANEL D: joint prediction surface (red–blue high-contrast)
  rec_grid  <- seq(quantile(dat4$z_rec, 0.02, na.rm = TRUE),
                   quantile(dat4$z_rec, 0.98, na.rm = TRUE),
                   length.out = 140)
  pi_grid   <- seq(quantile(dat4$pi_resid, 0.02, na.rm = TRUE),
                   quantile(dat4$pi_resid, 0.98, na.rm = TRUE),
                   length.out = 140)

  gridD <- expand.grid(z_rec = rec_grid, pi_resid = pi_grid) %>%
    as_tibble() %>%
    mutate(p = predict(m_res, newdata = ., type = "response"))

  pD <- ggplot(gridD, aes(x = z_rec, y = pi_resid, fill = p)) +
    geom_raster() +
    scale_fill_gradient2(
      low  = "#2b6cb0",   # deep blue  → low sweep probability
      mid  = "white",     # midpoint
      high = "#c53030",   # strong red → high sweep probability
      midpoint = 0.5,
      limits = c(0, 1),
      name = "P(Sweep)"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      x = "Recombination (z_rec)",
      y = "\u03C0 residual",
      title = "D. Joint prediction surface"
    )

  # PANEL E: forest plot of ORs
  pE <- ggplot(or_tab, aes(x = term, y = OR, ymin = lo, ymax = hi)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_pointrange() +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      x = NULL,
      y = "Odds ratio (per +1 SD)",
      title = "E. Effect sizes (logistic model)"
    )

  # Write separate panel files
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, paste0(out_prefix, "_PanelA.pdf")), pA, width = 6.5, height = 4.8)
  ggsave(file.path(out_dir, paste0(out_prefix, "_PanelB.pdf")), pB, width = 6.5, height = 4.8)
  ggsave(file.path(out_dir, paste0(out_prefix, "_PanelC.pdf")), pC, width = 6.5, height = 4.8)
  ggsave(file.path(out_dir, paste0(out_prefix, "_PanelD.pdf")), pD, width = 6.5, height = 5.2)
  ggsave(file.path(out_dir, paste0(out_prefix, "_PanelE.pdf")), pE, width = 6.5, height = 4.0)

  # Combined panel set if patchwork is available
  if (requireNamespace("patchwork", quietly = TRUE)) {
    suppressPackageStartupMessages(library(patchwork))
    fig <- (pA | pB) / (pC | pD) / pE +
      plot_annotation(
        title = "Selective sweeps are enriched in low recombination and show reduced diversity beyond recombination expectations",
        subtitle = "Recombination and diversity were rank-normalized; \u03C0 residuals from z_pi ~ z_rec."
      )
    ggsave(file.path(out_dir, paste0(out_prefix, "_Panels_Combined.pdf")),
           fig, width = 12, height = 13)
  } else {
    message("NOTE: 'patchwork' not installed; wrote separate panel PDFs only (A–E).")
  }

  invisible(list(dat4 = dat4, lm_pi = lm_pi, m_res = m_res,
                 panelA = pA, panelB = pB, panelC = pC, panelD = pD, panelE = pE))
}

## Auto-run (safe): only if sd3_new_af2 exists and contains required columns
if (exists("sd3_new_af2")) {
  if (all(c("recombination","pi_average","cont_label_h") %in% names(sd3_new_af2))) {
    cat("\n================ [17] MAKING PAPER FIGURE PANELS =================\n")
    FIG17 <- make_sweep_recomb_panels(
      sd3_new_af2 = sd3_new_af2,
      out_prefix = "Figure_Sweep_Recomb",
      out_dir = "."
    )
    cat("Wrote panels: Figure_Sweep_Recomb_Panel[A–E].pdf\n")
    if (file.exists("Figure_Sweep_Recomb_Panels_Combined.pdf")) {
      cat("Wrote combined: Figure_Sweep_Recomb_Panels_Combined.pdf\n")
    }
  } else {
    cat("\n[17] SKIP: sd3_new_af2 missing required columns for figure panels.\n")
  }
} else {
  cat("\n[17] SKIP: sd3_new_af2 not found; figure panels not generated.\n")
}
# Show the panels in the plot window (in addition to saving PDFs)
if (exists("FIG17")) {
  show_plot(FIG17$panelA)
  show_plot(FIG17$panelB)
  show_plot(FIG17$panelC)
  show_plot(FIG17$panelD)
  show_plot(FIG17$panelE)
  
  # If patchwork combined exists (it doesn't get returned), re-create it quickly:
  if (requireNamespace("patchwork", quietly = TRUE)) {
    suppressPackageStartupMessages(library(patchwork))
    fig_screen <- (FIG17$panelA | FIG17$panelB) / (FIG17$panelC | FIG17$panelD) / FIG17$panelE +
      plot_annotation(
        title = "Selective sweeps are enriched in low recombination and show reduced diversity beyond recombination expectations",
        subtitle = "Recombination and diversity were rank-normalized; π residuals from z_pi ~ z_rec."
      )
    show_plot(fig_screen)
  }
}


####sweep enrichment and snps density vs recombination
library(dplyr)
library(ggplot2)
library(scales)

if (exists("chr_ok") &&
    all(c("Chromosome_Names", "median_rec", "median_snp_density", "chr_type", "n_windows", "sweep_frac") %in% names(chr_ok))) {
  # short labels: remove leading "Chromosome_" (keeps "ZW" as "ZW")
  chr_ok2 <- chr_ok %>%
    mutate(chr_lab = gsub("^Chromosome_", "", Chromosome_Names))
  
  base_theme <- theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  # -------- Figure A: SNP density vs recombination --------
  pA <- ggplot(chr_ok2, aes(x = log10(median_rec), y = median_snp_density)) +
    geom_point(aes(shape = chr_type, size = n_windows), alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_text(aes(label = chr_lab), size = 3.2, vjust = -0.8, check_overlap = FALSE) +
    scale_size_continuous(name = "Windows (n)", range = c(2.5, 6)) +
    labs(
      x = "Chromosome median recombination (log10 raw)",
      y = "Median SNP density (n_snps_af / no_sites)",
      title = "SNP density vs recombination"
    ) +
    base_theme
  
  print(pA)
  ggsave("chr_SNPdensity_vs_recombination_labeled.pdf", pA, width = 7.5, height = 6)
  ggsave("chr_SNPdensity_vs_recombination_labeled.png", pA, width = 7.5, height = 6, dpi = 300)
  
  # -------- Figure B: Sweep enrichment vs recombination --------
  pB <- ggplot(chr_ok2, aes(x = log10(median_rec), y = sweep_frac)) +
    geom_point(aes(shape = chr_type, size = n_windows), alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
    geom_text(aes(label = chr_lab), size = 3.2, vjust = -0.8, check_overlap = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_size_continuous(name = "Windows (n)", range = c(2.5, 6)) +
    labs(
      x = "Chromosome median recombination (log10 raw)",
      y = "Sweep windows (fraction)",
      title = "Sweep enrichment vs recombination"
    ) +
    base_theme
  
  print(pB)
  ggsave("chr_SweepEnrichment_vs_recombination_labeled.pdf", pB, width = 7.5, height = 6)
  ggsave("chr_SweepEnrichment_vs_recombination_labeled.png", pB, width = 7.5, height = 6, dpi = 300)
} else {
  cat("\n[18A] SKIP: chr_ok missing required columns for chromosome recombination summary plots.\n")
}

############################################################
## [18] OPTIONAL QC FIGURES: threshold sensitivity + bootstrap stability
############################################################
make_selection_qc_figures <- function(sd3_clusters, agree_dat,
                                      out_prefix = "Selection_QC",
                                      out_dir = ".",
                                      rf_grid = seq(0.80, 0.95, by = 0.05),
                                      gmm_grid = seq(0.80, 0.95, by = 0.05),
                                      rf_base = 0.90,
                                      gmm_base = 0.90,
                                      classes = c("Geographic sweep", "Balancing selection"),
                                      n_boot = 200) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    library(data.table)
  })

  req_sd3 <- c("Chromosome_Names", ".start_bp", "p_cluster", "rf_pmax")
  req_agd <- c("Chromosome_Names", ".start_bp", "gmm_label", "rf_label", "cont_label_h")
  if (!all(req_sd3 %in% names(sd3_clusters))) {
    miss <- setdiff(req_sd3, names(sd3_clusters))
    stop("make_selection_qc_figures(): missing in sd3_clusters: ", paste(miss, collapse = ", "))
  }
  if (!all(req_agd %in% names(agree_dat))) {
    miss <- setdiff(req_agd, names(agree_dat))
    stop("make_selection_qc_figures(): missing in agree_dat: ", paste(miss, collapse = ", "))
  }

  dat <- sd3_clusters %>%
    select(Chromosome_Names, .start_bp, p_cluster, rf_pmax, dplyr::any_of(".end_bp")) %>%
    left_join(
      agree_dat %>% select(Chromosome_Names, .start_bp, gmm_label, rf_label, cont_label_h),
      by = c("Chromosome_Names", ".start_bp")
    ) %>%
    mutate(
      window_id = paste(Chromosome_Names, .start_bp, sep = ":"),
      agree_k = (gmm_label == rf_label) + (gmm_label == cont_label_h) + (rf_label == cont_label_h),
      is_class = cont_label_h %in% classes
    )

  select_windows <- function(df, rf_thr, gmm_thr) {
    df %>%
      filter(
        is_class,
        agree_k >= 2,
        is.finite(rf_pmax), rf_pmax >= rf_thr,
        is.finite(p_cluster), p_cluster >= gmm_thr
      )
  }

  baseline <- select_windows(dat, rf_base, gmm_base)
  base_ids <- unique(baseline$window_id)
  if (length(base_ids) == 0) {
    warning("make_selection_qc_figures(): no windows selected at baseline thresholds.")
  }

  grid <- tidyr::expand_grid(rf_thr = rf_grid, gmm_thr = gmm_grid) %>%
    rowwise() %>%
    do({
      sel <- select_windows(dat, .$rf_thr, .$gmm_thr)
      ids <- unique(sel$window_id)
      inter <- length(intersect(ids, base_ids))
      uni <- length(unique(c(ids, base_ids)))
      tibble(
        rf_thr = .$rf_thr,
        gmm_thr = .$gmm_thr,
        n_total = nrow(sel),
        n_sweep = sum(sel$cont_label_h == "Geographic sweep", na.rm = TRUE),
        n_balancing = sum(sel$cont_label_h == "Balancing selection", na.rm = TRUE),
        jaccard_vs_base = if (uni > 0) inter / uni else NA_real_
      )
    }) %>%
    ungroup()

  p_heat_count <- ggplot(grid, aes(x = rf_thr, y = gmm_thr, fill = n_total)) +
    geom_tile() +
    geom_text(aes(label = n_total), size = 3) +
    scale_fill_viridis_c(option = "D", name = "Selected windows") +
    scale_x_continuous(breaks = rf_grid) +
    scale_y_continuous(breaks = gmm_grid) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      x = "RF confidence threshold",
      y = "GMM confidence threshold",
      title = "Threshold sensitivity: selected-window counts"
    )

  p_heat_jacc <- ggplot(grid, aes(x = rf_thr, y = gmm_thr, fill = jaccard_vs_base)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", jaccard_vs_base)), size = 3) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1), name = "Jaccard vs baseline") +
    scale_x_continuous(breaks = rf_grid) +
    scale_y_continuous(breaks = gmm_grid) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      x = "RF confidence threshold",
      y = "GMM confidence threshold",
      title = "Threshold sensitivity: overlap with baseline"
    )

  set.seed(20260305)
  boot_stats <- lapply(seq_len(n_boot), function(b) {
    idx <- sample.int(nrow(dat), nrow(dat), replace = TRUE)
    db <- dat[idx, , drop = FALSE]
    selb <- select_windows(db, rf_base, gmm_base)
    n_total <- nrow(selb)
    if (n_total == 0) {
      return(data.frame(boot = b, n_total = 0, frac_sweep = NA_real_, frac_balancing = NA_real_))
    }
    data.frame(
      boot = b,
      n_total = n_total,
      frac_sweep = mean(selb$cont_label_h == "Geographic sweep", na.rm = TRUE),
      frac_balancing = mean(selb$cont_label_h == "Balancing selection", na.rm = TRUE)
    )
  }) %>% bind_rows()

  boot_long <- boot_stats %>%
    select(boot, frac_sweep, frac_balancing) %>%
    pivot_longer(cols = c(frac_sweep, frac_balancing), names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
                           frac_sweep = "Sweep fraction",
                           frac_balancing = "Balancing fraction"))

  p_boot_counts <- ggplot(boot_stats, aes(x = n_total)) +
    geom_histogram(bins = 25, fill = "grey50", color = "white") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      x = "Selected windows per bootstrap",
      y = "Count",
      title = "Bootstrap stability: selected-window totals"
    )

  p_boot_frac <- ggplot(boot_long, aes(x = metric, y = value, fill = metric)) +
    geom_violin(trim = TRUE, alpha = 0.8, color = "grey30") +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    labs(
      x = NULL,
      y = "Fraction",
      title = "Bootstrap stability: class composition"
    )

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ThresholdCountHeatmap.pdf")), p_heat_count, width = 7, height = 5.5)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ThresholdCountHeatmap.png")), p_heat_count, width = 7, height = 5.5, dpi = 300)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ThresholdJaccardHeatmap.pdf")), p_heat_jacc, width = 7, height = 5.5)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ThresholdJaccardHeatmap.png")), p_heat_jacc, width = 7, height = 5.5, dpi = 300)
  ggsave(file.path(out_dir, paste0(out_prefix, "_BootstrapCounts.pdf")), p_boot_counts, width = 7, height = 5.5)
  ggsave(file.path(out_dir, paste0(out_prefix, "_BootstrapCounts.png")), p_boot_counts, width = 7, height = 5.5, dpi = 300)
  ggsave(file.path(out_dir, paste0(out_prefix, "_BootstrapFractions.pdf")), p_boot_frac, width = 7, height = 5.5)
  ggsave(file.path(out_dir, paste0(out_prefix, "_BootstrapFractions.png")), p_boot_frac, width = 7, height = 5.5, dpi = 300)

  invisible(list(grid = grid, boot = boot_stats,
                 p_heat_count = p_heat_count, p_heat_jacc = p_heat_jacc,
                 p_boot_counts = p_boot_counts, p_boot_frac = p_boot_frac))
}

if (exists("sd3_clusters") && exists("agree_dat")) {
  if (all(c("Chromosome_Names", ".start_bp", "p_cluster", "rf_pmax") %in% names(sd3_clusters)) &&
      all(c("Chromosome_Names", ".start_bp", "gmm_label", "rf_label", "cont_label_h") %in% names(agree_dat))) {
    cat("\n================ [18] MAKING QC SENSITIVITY FIGURES =================\n")
    QC_FIGS <- make_selection_qc_figures(
      sd3_clusters = sd3_clusters,
      agree_dat = agree_dat,
      out_prefix = "Selection_QC",
      out_dir = "."
    )
    cat("Wrote QC figures: Selection_QC_* (threshold and bootstrap diagnostics)\n")
  } else {
    cat("\n[18] SKIP: missing required columns in sd3_clusters/agree_dat for QC figure module.\n")
  }
} else {
  cat("\n[18] SKIP: sd3_clusters or agree_dat not found; QC figure module not run.\n")
}

############################################################
## [19] SWEEP SUBTYPE COMPARISON (G4): compare sweep components
############################################################

make_sweep_subtype_g4 <- function(sd3_g4,
                                  out_prefix = "SweepSubtype_Comparison_G4",
                                  out_dir = ".") {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    library(data.table)
  })

  need <- c("cluster", "gmm_model_raw", "Chromosome_Names", ".start_bp",
            "z_fst", "z_dxy", "z_pi", "z_taj", "z_f3", "z_rec")
  miss <- setdiff(need, names(sd3_g4))
  if (length(miss) > 0) {
    stop("make_sweep_subtype_g4(): missing columns: ", paste(miss, collapse = ", "))
  }

  sw <- sd3_g4 %>%
    mutate(cluster = as.factor(cluster)) %>%
    filter(gmm_model_raw == "Geographic sweep", !is.na(cluster))

  if (nrow(sw) == 0) {
    stop("make_sweep_subtype_g4(): no sweep windows found.")
  }

  sweep_clusters <- sort(unique(as.character(sw$cluster)))
  if (length(sweep_clusters) < 2) {
    stop("make_sweep_subtype_g4(): fewer than 2 sweep clusters found.")
  }
  # Keep top two by size if >2 (defensive)
  top2 <- sw %>%
    group_by(cluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n)) %>%
    slice_head(n = 2) %>%
    pull(cluster) %>%
    as.character()
  sw <- sw %>% filter(as.character(cluster) %in% top2) %>%
    mutate(cluster = factor(as.character(cluster), levels = sort(unique(as.character(cluster)))))

  prof <- sw %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      mean_z_fst = mean(z_fst, na.rm = TRUE),
      mean_z_dxy = mean(z_dxy, na.rm = TRUE),
      mean_z_pi  = mean(z_pi, na.rm = TRUE),
      mean_z_taj = mean(z_taj, na.rm = TRUE),
      mean_z_f3  = mean(z_f3, na.rm = TRUE),
      mean_z_rec = mean(z_rec, na.rm = TRUE),
      .groups = "drop"
    )

  # Wilcoxon tests for key contrasts
  wilcox_stat <- function(v, metric_name) {
    x <- sw %>% filter(cluster == levels(sw$cluster)[1]) %>% pull(dplyr::all_of(v))
    y <- sw %>% filter(cluster == levels(sw$cluster)[2]) %>% pull(dplyr::all_of(v))
    wt <- wilcox.test(x, y)
    tibble(
      metric = metric_name,
      p_value = wt$p.value,
      median_cluster1 = median(x, na.rm = TRUE),
      median_cluster2 = median(y, na.rm = TRUE),
      delta_median_c1_minus_c2 = median(x, na.rm = TRUE) - median(y, na.rm = TRUE)
    )
  }

  tests <- bind_rows(
    wilcox_stat("z_fst", "z_fst"),
    wilcox_stat("z_pi",  "z_pi"),
    wilcox_stat("z_dxy", "z_dxy"),
    wilcox_stat("z_taj", "z_taj"),
    wilcox_stat("z_f3",  "z_f3"),
    wilcox_stat("z_rec", "z_rec")
  ) %>% arrange(p_value)

  # Chromosome composition by subtype
  chr_comp <- sw %>%
    dplyr::count(cluster, Chromosome_Names, name = "n_windows") %>%
    group_by(cluster) %>%
    mutate(frac = n_windows / sum(n_windows)) %>%
    ungroup()

  # Profile heatmap
  prof_long <- prof %>%
    select(cluster, starts_with("mean_z_")) %>%
    pivot_longer(-cluster, names_to = "metric", values_to = "value")

  p_profile <- ggplot(prof_long, aes(x = metric, y = cluster, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c53030", midpoint = 0) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(x = NULL, y = "Sweep subtype (cluster)", fill = "Mean z", title = "G4 sweep subtype profile")

  # Distribution panels
  sw_long <- sw %>%
    select(cluster, z_fst, z_pi, z_dxy, z_taj, z_f3, z_rec) %>%
    pivot_longer(-cluster, names_to = "metric", values_to = "value")

  p_dist <- ggplot(sw_long, aes(x = cluster, y = value, fill = cluster)) +
    geom_violin(trim = TRUE, alpha = 0.75) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    labs(x = "Sweep subtype (cluster)", y = "Value", title = "Sweep subtype distributions (G4)")

  # Chromosome composition
  p_chr <- ggplot(chr_comp, aes(x = Chromosome_Names, y = frac, fill = cluster)) +
    geom_col(position = "dodge") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(x = "Chromosome", y = "Fraction of subtype windows", fill = "Subtype",
         title = "Chromosome composition of sweep subtypes (G4)")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  data.table::fwrite(prof, file.path(out_dir, paste0(out_prefix, "_CentroidProfile.tsv")), sep = "\t")
  data.table::fwrite(tests, file.path(out_dir, paste0(out_prefix, "_WilcoxonTests.tsv")), sep = "\t")
  data.table::fwrite(chr_comp, file.path(out_dir, paste0(out_prefix, "_ChromosomeComposition.tsv")), sep = "\t")

  ggsave(file.path(out_dir, paste0(out_prefix, "_ProfileHeatmap.pdf")), p_profile, width = 8.0, height = 4.2)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ProfileHeatmap.png")), p_profile, width = 8.0, height = 4.2, dpi = 300)
  ggsave(file.path(out_dir, paste0(out_prefix, "_DistributionPanels.pdf")), p_dist, width = 11.0, height = 7.0)
  ggsave(file.path(out_dir, paste0(out_prefix, "_DistributionPanels.png")), p_dist, width = 11.0, height = 7.0, dpi = 300)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ChromosomeComposition.pdf")), p_chr, width = 12.0, height = 5.5)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ChromosomeComposition.png")), p_chr, width = 12.0, height = 5.5, dpi = 300)

  invisible(list(profile = prof, tests = tests, chr_comp = chr_comp))
}

if (exists("bestG") && is.finite(bestG) && bestG == 4 && exists("res_Gbest")) {
  cat("\n================ [19] MAKING SWEEP SUBTYPE COMPARISON (G4) =================\n")
  SWEEP_SUBTYPE_G4 <- tryCatch(
    make_sweep_subtype_g4(
      sd3_g4 = res_Gbest$sd3_clusters,
      out_prefix = "SweepSubtype_Comparison_G4",
      out_dir = "."
    ),
    error = function(e) {
      cat("[19] SKIP:", conditionMessage(e), "\n")
      NULL
    }
  )
} else {
  cat("\n[19] SKIP: bestG is not 4 or res_Gbest not available; sweep subtype module not run.\n")
}

############################################################
## [20] BACKGROUND SELECTION (BGS) VS SWEEP-EXCESS MODULE
############################################################

make_bgs_vs_sweep_report <- function(sd3_raw,
                                     out_prefix = "BGS",
                                     out_dir = ".") {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(data.table)
  })

  chr_col <- .find_col(sd3_raw, c("Chromosome_Names", "Chromosomes", "Chromosome", "chr"))
  rec_col <- .find_col(sd3_raw, c("recombination", "rec", "rho", "r"))
  pi_col  <- .find_col(sd3_raw, c("pi_average", "avg_pi", "pi", "avg_pi_east"))
  fst_col <- .find_col(sd3_raw, c("avg_wc_fst", "fst", "FST"))
  dxy_col <- .find_col(sd3_raw, c("avg_dxy", "dxy", "DXY"))
  taj_col <- .find_col(sd3_raw, c("TajimaD", "tajD", "tajima_d", "tajimaD"))
  ws_col  <- .find_col(sd3_raw, c("window_start", ".start_bp", "start", "start_bp"))
  we_col  <- .find_col(sd3_raw, c("window_end", ".end_bp", "end", "end_bp"))
  mid_col <- .find_col(sd3_raw, c("midpoint", ".mid_bp", "mid_bp", "window_mid"))
  lbl_col <- .find_col(sd3_raw, c("cont_label_h", "gmm_model_raw", "rf_model"))

  req <- c(chr_col, rec_col, pi_col, fst_col, dxy_col)
  if (any(is.na(req))) {
    stop("make_bgs_vs_sweep_report(): missing required columns for BGS module.")
  }

  dt2 <- sd3_raw %>%
    mutate(
      Chromosome_Names = as.character(.data[[chr_col]]),
      recombination = as.numeric(.data[[rec_col]]),
      pi_average = as.numeric(.data[[pi_col]]),
      avg_wc_fst = as.numeric(.data[[fst_col]]),
      avg_dxy = as.numeric(.data[[dxy_col]]),
      TajimaD = if (!is.na(taj_col)) as.numeric(.data[[taj_col]]) else NA_real_,
      window_start = if (!is.na(ws_col)) as.numeric(.data[[ws_col]]) else NA_real_,
      window_end = if (!is.na(we_col)) as.numeric(.data[[we_col]]) else NA_real_,
      midpoint = if (!is.na(mid_col)) as.numeric(.data[[mid_col]]) else (window_start + window_end) / 2,
      sweep_label = if (!is.na(lbl_col)) as.character(.data[[lbl_col]]) else NA_character_
    ) %>%
    filter(
      !is.na(Chromosome_Names),
      is.finite(recombination), is.finite(pi_average),
      is.finite(avg_wc_fst), is.finite(avg_dxy)
    )

  eps_rec <- max(min(dt2$recombination[dt2$recombination > 0], na.rm = TRUE) * 0.5, 1e-12)
  eps_pi  <- max(min(dt2$pi_average[dt2$pi_average > 0], na.rm = TRUE) * 0.5, 1e-12)
  eps_dxy <- max(min(dt2$avg_dxy[dt2$avg_dxy > 0], na.rm = TRUE) * 0.5, 1e-12)

  dt2 <- dt2 %>%
    mutate(
      log_rec = log10(recombination + eps_rec),
      log_pi  = log10(pi_average + eps_pi),
      log_dxy = log10(avg_dxy + eps_dxy)
    )

  m_pi1  <- lm(log_pi ~ log_rec, data = dt2)
  m_pi2  <- lm(log_pi ~ poly(log_rec, 2, raw = TRUE), data = dt2)
  m_fst1 <- lm(avg_wc_fst ~ log_rec, data = dt2)
  m_fst2 <- lm(avg_wc_fst ~ log_rec + log_pi, data = dt2)
  m_fst3 <- lm(avg_wc_fst ~ log_rec + log_pi + log_dxy, data = dt2)
  m_dxy1 <- lm(log_dxy ~ log_rec, data = dt2)
  m_dxy2 <- lm(log_dxy ~ log_rec + log_pi, data = dt2)

  dt2 <- dt2 %>%
    mutate(
      pi_pred_from_rec = predict(m_pi2, newdata = dt2),
      pi_resid_from_rec = log_pi - pi_pred_from_rec
    )

  q_pi_low  <- as.numeric(quantile(dt2$pi_resid_from_rec, 0.10, na.rm = TRUE))
  q_fst_hi  <- as.numeric(quantile(dt2$avg_wc_fst, 0.90, na.rm = TRUE))
  q_dxy_low <- as.numeric(quantile(dt2$avg_dxy, 0.10, na.rm = TRUE))
  q_taj_low <- as.numeric(quantile(dt2$TajimaD, 0.10, na.rm = TRUE))

  dt2 <- dt2 %>%
    mutate(
      low_pi_resid = pi_resid_from_rec <= q_pi_low,
      high_fst = avg_wc_fst >= q_fst_hi,
      low_dxy = avg_dxy <= q_dxy_low,
      low_taj = ifelse(is.finite(TajimaD), TajimaD <= q_taj_low, FALSE),
      BGS_class = case_when(
        low_pi_resid & (high_fst | low_taj) ~ "sweep_like_excess",
        low_pi_resid & !high_fst & low_dxy ~ "BGS_like",
        TRUE ~ "other"
      ),
      is_sweep_model = !is.na(sweep_label) & grepl("Geographic sweep|Sweep-like", sweep_label, ignore.case = TRUE),
      is_sweep_beyond_bgs = is_sweep_model & BGS_class == "sweep_like_excess"
    )

  dt2 <- dt2 %>%
    mutate(rec_bin = ntile(log_rec, 10))

  bin_sum <- dt2 %>%
    group_by(rec_bin) %>%
    summarise(
      n = n(),
      med_rec = median(recombination, na.rm = TRUE),
      med_pi = median(pi_average, na.rm = TRUE),
      med_fst = median(avg_wc_fst, na.rm = TRUE),
      med_dxy = median(avg_dxy, na.rm = TRUE),
      frac_bgs_like = mean(BGS_class == "BGS_like", na.rm = TRUE),
      frac_sweep_excess = mean(BGS_class == "sweep_like_excess", na.rm = TRUE),
      .groups = "drop"
    )

  run_summary <- tibble(
    n_windows = nrow(dt2),
    n_sweep_models = sum(dt2$is_sweep_model, na.rm = TRUE),
    n_bgs_like = sum(dt2$BGS_class == "BGS_like", na.rm = TRUE),
    n_sweep_like_excess = sum(dt2$BGS_class == "sweep_like_excess", na.rm = TRUE),
    n_sweep_beyond_bgs = sum(dt2$is_sweep_beyond_bgs, na.rm = TRUE),
    frac_sweep_beyond_bgs_of_sweeps = mean(dt2$is_sweep_beyond_bgs[dt2$is_sweep_model], na.rm = TRUE)
  )

  have_ppcor <- requireNamespace("ppcor", quietly = TRUE)
  have_mgcv  <- requireNamespace("mgcv", quietly = TRUE)
  have_broom <- requireNamespace("broom", quietly = TRUE)

  pc_pi <- pc_fst <- pc_dxy <- "ppcor package not installed; partial correlations skipped."
  if (have_ppcor) {
    ctrl_pi <- dt2 %>% select(avg_wc_fst, log_dxy) %>% as.data.frame()
    ctrl_fst <- dt2 %>% select(log_pi, log_dxy) %>% as.data.frame()
    ctrl_dxy <- dt2 %>% select(log_pi, avg_wc_fst) %>% as.data.frame()
    pc_pi <- ppcor::pcor.test(dt2$log_rec, dt2$log_pi, ctrl_pi, method = "spearman")
    pc_fst <- ppcor::pcor.test(dt2$log_rec, dt2$avg_wc_fst, ctrl_fst, method = "spearman")
    pc_dxy <- ppcor::pcor.test(dt2$log_rec, dt2$log_dxy, ctrl_dxy, method = "spearman")
  }

  g_pi <- g_fst <- g_dxy <- "mgcv package not installed; GAM summaries skipped."
  if (have_mgcv) {
    g_pi <- mgcv::gam(log_pi ~ s(log_rec, k = 6), data = dt2, method = "REML")
    g_fst <- mgcv::gam(avg_wc_fst ~ s(log_rec, k = 6) + s(log_pi, k = 6), data = dt2, method = "REML")
    g_dxy <- mgcv::gam(log_dxy ~ s(log_rec, k = 6), data = dt2, method = "REML")
  }

  model_coefs <- function(m) {
    if (have_broom) return(broom::tidy(m))
    as.data.frame(summary(m)$coefficients) %>% tibble::rownames_to_column("term")
  }

  top_sweep <- dt2 %>%
    filter(BGS_class == "sweep_like_excess") %>%
    arrange(pi_resid_from_rec) %>%
    select(Chromosome_Names, midpoint, pi_average, avg_wc_fst, avg_dxy, TajimaD, recombination, pi_resid_from_rec) %>%
    slice_head(n = 10)

  top_bgs <- dt2 %>%
    filter(BGS_class == "BGS_like") %>%
    arrange(pi_resid_from_rec) %>%
    select(Chromosome_Names, midpoint, pi_average, avg_wc_fst, avg_dxy, TajimaD, recombination, pi_resid_from_rec) %>%
    slice_head(n = 10)

  chr_sizes <- dt2 %>%
    group_by(Chromosome_Names) %>%
    summarise(chr_len_bp = max(window_end, midpoint, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      chr_index = suppressWarnings(as.integer(gsub("^Chromosome_", "", Chromosome_Names))),
      chr_index = ifelse(is.na(chr_index), 9999L, chr_index)
    ) %>%
    arrange(chr_index, Chromosome_Names) %>%
    mutate(
      Chromosome_Names = factor(Chromosome_Names, levels = Chromosome_Names),
      y = as.numeric(Chromosome_Names)
    )

  lane_half <- 0.22
  p_chr <- ggplot() +
    geom_rect(
      data = chr_sizes,
      aes(xmin = 0, xmax = chr_len_bp / 1e6, ymin = y - lane_half, ymax = y + lane_half),
      fill = "grey92", color = "grey80", linewidth = 0.2
    ) +
    geom_point(
      data = dt2 %>% filter(is_sweep_model, is.finite(midpoint)) %>%
        mutate(Chromosome_Names = factor(Chromosome_Names, levels = levels(chr_sizes$Chromosome_Names))) %>%
        left_join(chr_sizes %>% select(Chromosome_Names, y), by = "Chromosome_Names"),
      aes(x = midpoint / 1e6, y = y, color = "All sweep windows"),
      alpha = 0.22, size = 0.7
    ) +
    geom_point(
      data = dt2 %>% filter(is_sweep_beyond_bgs, is.finite(midpoint)) %>%
        mutate(Chromosome_Names = factor(Chromosome_Names, levels = levels(chr_sizes$Chromosome_Names))) %>%
        left_join(chr_sizes %>% select(Chromosome_Names, y), by = "Chromosome_Names"),
      aes(x = midpoint / 1e6, y = y, color = "Sweep beyond BGS"),
      alpha = 0.95, size = 1.3
    ) +
    scale_color_manual(values = c("All sweep windows" = "grey35", "Sweep beyond BGS" = "#cc1f1a")) +
    scale_y_reverse(
      breaks = chr_sizes$y,
      labels = as.character(chr_sizes$Chromosome_Names),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    labs(
      x = "Position (Mb)",
      y = "Chromosome",
      color = NULL,
      title = "Sweep-model windows and the subset beyond the BGS baseline",
      subtitle = "Grey = all GMM sweep windows; red = sweep-like excess after the BGS screen, not the top-N continuous-tail plot"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "bottom")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_report <- file.path(out_dir, paste0(out_prefix, "_FullStatReport.txt"))
  out_plot_pdf <- file.path(out_dir, paste0(out_prefix, "_ChromosomeSweepBeyondBGS.pdf"))
  out_plot_png <- file.path(out_dir, paste0(out_prefix, "_ChromosomeSweepBeyondBGS.png"))

  capture <- capture.output({
    cat("\n\n==================== REPORT START ====================\n")
    cat("\n--- DATA SUMMARY ---\n")
    cat("Number of windows:", nrow(dt2), "\n")
    cat("Number of chromosomes:", dplyr::n_distinct(dt2$Chromosome_Names), "\n")

    cat("\n--- SPEARMAN CORRELATIONS ---\n")
    print(cor.test(dt2$log_rec, dt2$log_pi, method = "spearman"))
    print(cor.test(dt2$log_rec, dt2$avg_wc_fst, method = "spearman"))
    print(cor.test(dt2$log_rec, dt2$log_dxy, method = "spearman"))
    if (all(is.finite(dt2$TajimaD))) print(cor.test(dt2$log_rec, dt2$TajimaD, method = "spearman"))

    cat("\n--- PARTIAL CORRELATIONS ---\n")
    print(pc_pi); print(pc_fst); print(pc_dxy)

    cat("\n--- MODEL R-SQUARED ---\n")
    print(data.frame(
      model = c("m_pi1","m_pi2","m_fst1","m_fst2","m_fst3","m_dxy1","m_dxy2"),
      r2 = c(summary(m_pi1)$r.squared,
             summary(m_pi2)$r.squared,
             summary(m_fst1)$r.squared,
             summary(m_fst2)$r.squared,
             summary(m_fst3)$r.squared,
             summary(m_dxy1)$r.squared,
             summary(m_dxy2)$r.squared)
    ))

    cat("\n--- MODEL COEFFICIENTS ---\n")
    cat("\n[pi model]\n"); print(model_coefs(m_pi2))
    cat("\n[FST model]\n"); print(model_coefs(m_fst3))
    cat("\n[dXY model]\n"); print(model_coefs(m_dxy2))

    cat("\n--- GAM SUMMARIES ---\n")
    cat("\n[GAM pi ~ rec]\n"); print(if (is.character(g_pi)) g_pi else summary(g_pi))
    cat("\n[GAM fst ~ rec + pi]\n"); print(if (is.character(g_fst)) g_fst else summary(g_fst))
    cat("\n[GAM dxy ~ rec]\n"); print(if (is.character(g_dxy)) g_dxy else summary(g_dxy))

    cat("\n--- KEY QUANTILES ---\n")
    print(list(
      pi = quantile(dt2$pi_average, probs = c(0, .1, .5, .9, 1), na.rm = TRUE),
      fst = quantile(dt2$avg_wc_fst, probs = c(0, .1, .5, .9, 1), na.rm = TRUE),
      dxy = quantile(dt2$avg_dxy, probs = c(0, .1, .5, .9, 1), na.rm = TRUE),
      rec = quantile(dt2$recombination, probs = c(0, .1, .5, .9, 1), na.rm = TRUE),
      tajd = quantile(dt2$TajimaD, probs = c(0, .1, .5, .9, 1), na.rm = TRUE)
    ))

    cat("\n--- BGS / SWEEP CLASS COUNTS ---\n")
    print(table(dt2$BGS_class))

    cat("\n--- PI RESIDUAL SUMMARY ---\n")
    print(summary(dt2$pi_resid_from_rec))

    cat("\n--- RECOMBINATION BIN SUMMARY ---\n")
    print(bin_sum)

    cat("\n--- BGS RUN SUMMARY ---\n")
    print(run_summary)

    cat("\n--- TOP SWEEP-LIKE WINDOWS ---\n")
    print(top_sweep)

    cat("\n--- TOP BGS-LIKE WINDOWS ---\n")
    print(top_bgs)

    cat("\n==================== REPORT END ====================\n")
  })
  writeLines(capture, con = out_report)

  data.table::fwrite(dt2 %>% dplyr::count(BGS_class, name = "n_windows"),
                     file.path(out_dir, paste0(out_prefix, "_ClassCounts.tsv")), sep = "\t")
  data.table::fwrite(bin_sum, file.path(out_dir, paste0(out_prefix, "_RecombinationBinSummary.tsv")), sep = "\t")
  data.table::fwrite(run_summary, file.path(out_dir, paste0(out_prefix, "_RunSummary.tsv")), sep = "\t")
  data.table::fwrite(top_sweep, file.path(out_dir, paste0(out_prefix, "_TopSweepLike.tsv")), sep = "\t")
  data.table::fwrite(top_bgs, file.path(out_dir, paste0(out_prefix, "_TopBGSLike.tsv")), sep = "\t")
  data.table::fwrite(dt2 %>% filter(is_sweep_beyond_bgs),
                     file.path(out_dir, paste0(out_prefix, "_SweepBeyondBGS_Windows.tsv")), sep = "\t")

  ggsave(out_plot_pdf, p_chr, width = 12, height = 7)
  ggsave(out_plot_png, p_chr, width = 12, height = 7, dpi = 300)

  invisible(list(
    dt2 = dt2,
    run_summary = run_summary,
    bin_sum = bin_sum,
    report_file = out_report,
    plot_file_pdf = out_plot_pdf,
    plot_file_png = out_plot_png
  ))
}

if (exists("sd3_new_af2")) {
  cat("\n================ [20] BACKGROUND SELECTION VS SWEEP-EXCESS =================\n")
  BGS_REPORT <- tryCatch(
    make_bgs_vs_sweep_report(
      sd3_raw = sd3_new_af2,
      out_prefix = "BGS",
      out_dir = "."
    ),
    error = function(e) {
      cat("[20] SKIP:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(BGS_REPORT)) {
    cat("[20] Wrote: BGS_FullStatReport.txt, BGS_* summary tables, and BGS_ChromosomeSweepBeyondBGS.*\n")
  }
} else {
  cat("\n[20] SKIP: sd3_new_af2 not found; BGS module not run.\n")
}

############################################################
## [21] GENE DENSITY ASSOCIATION FOR BGS / SWEEP CLASSES
############################################################

make_bgs_gene_density_module <- function(bgs_dt,
                                         sd3_clusters,
                                         out_prefix = "BGS_GeneDensity",
                                         out_dir = ".",
                                         gtf_path = NULL) {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
  })

  need_bgs <- c("Chromosome_Names", "window_start", "window_end", "recombination", "BGS_class")
  miss_bgs <- setdiff(need_bgs, names(bgs_dt))
  if (length(miss_bgs) > 0) {
    stop("make_bgs_gene_density_module(): missing bgs_dt columns: ", paste(miss_bgs, collapse = ", "))
  }
  need_sd3 <- c("Chromosomes", "Chromosome_Names", "window_start", "window_end")
  miss_sd3 <- setdiff(need_sd3, names(sd3_clusters))
  if (length(miss_sd3) > 0) {
    stop("make_bgs_gene_density_module(): missing sd3_clusters columns: ", paste(miss_sd3, collapse = ", "))
  }

  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE) ||
      !requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("make_bgs_gene_density_module(): requires packages GenomicRanges, IRanges, rtracklayer.")
  }

  if (is.null(gtf_path) || !nzchar(gtf_path)) {
    cand <- c(
      file.path(getwd(), "Pantherophis_alleghaniensis_obsoletus.gtf"),
      file.path(script_dir, "Pantherophis_alleghaniensis_obsoletus.gtf"),
      file.path(project_dir, "Pantherophis_alleghaniensis_obsoletus.gtf"),
      file.path(project_dir, "..", "ratsnake_GFF", "Pantherophis_alleghaniensis_obsoletus.gtf"),
      file.path(project_dir, "..", "..", "ratsnake_GFF", "Pantherophis_alleghaniensis_obsoletus.gtf"),
      Sys.getenv("PIPELINE_GTF_PATH", unset = "")
    )
    gtf_path <- cand[nzchar(cand) & file.exists(cand)][1]
  }
  if (is.na(gtf_path) || !file.exists(gtf_path)) {
    stop("make_bgs_gene_density_module(): GTF file not found.")
  }

  bgs <- as.data.table(bgs_dt)[, .(
    Chromosome_Names = as.character(Chromosome_Names),
    window_start = as.integer(window_start),
    window_end = as.integer(window_end),
    recombination = as.numeric(recombination),
    BGS_class = as.character(BGS_class)
  )]

  win <- as.data.table(sd3_clusters)[, .(
    Chromosomes = as.character(Chromosomes),
    Chromosome_Names = as.character(Chromosome_Names),
    window_start = as.integer(window_start),
    window_end = as.integer(window_end)
  )]
  setkey(win, Chromosome_Names, window_start, window_end)
  setkey(bgs, Chromosome_Names, window_start, window_end)
  win <- win[bgs, nomatch = 0L]
  if (nrow(win) == 0) stop("No overlapping windows between BGS table and sd3_clusters.")

  ann <- rtracklayer::import(gtf_path)
  ann <- ann[as.character(S4Vectors::mcols(ann)$type) %in% c("gene", "CDS")]
  is_gene <- as.character(S4Vectors::mcols(ann)$type) == "gene"
  feat <- if (sum(is_gene, na.rm = TRUE) > 0) ann[is_gene] else ann
  gene_id <- S4Vectors::mcols(feat)$gene_id
  if (is.null(gene_id)) gene_id <- S4Vectors::mcols(feat)$transcript_id
  if (is.null(gene_id)) gene_id <- paste0("g", seq_along(feat))
  S4Vectors::mcols(feat)$gene_id_fix <- as.character(gene_id)

  win_gr <- GenomicRanges::GRanges(
    seqnames = win$Chromosomes,
    ranges = IRanges::IRanges(start = pmax(1L, win$window_start),
                              end = pmax(win$window_start, win$window_end))
  )
  hits <- GenomicRanges::findOverlaps(win_gr, feat, ignore.strand = TRUE)
  hit_dt <- data.table(
    win_i = S4Vectors::queryHits(hits),
    gene_id = S4Vectors::mcols(feat)$gene_id_fix[S4Vectors::subjectHits(hits)]
  )

  win[, win_i := .I]
  gene_win <- hit_dt[, .(gene_n = uniqueN(gene_id)), by = win_i]
  win <- merge(win, gene_win, by = "win_i", all.x = TRUE)
  win[is.na(gene_n), gene_n := 0L]
  win[, window_bp := pmax(1L, window_end - window_start + 1L)]
  win[, gene_density_perMb := gene_n / (window_bp / 1e6)]

  join_dt <- merge(
    bgs,
    win[, .(Chromosome_Names, window_start, window_end, gene_n, gene_density_perMb)],
    by = c("Chromosome_Names", "window_start", "window_end"),
    all.x = TRUE
  )
  join_dt <- join_dt[is.finite(gene_density_perMb)]
  join_dt[, BGS_class := factor(BGS_class, levels = c("BGS_like", "other", "sweep_like_excess"))]

  class_sum <- join_dt[, .(
    n = .N,
    mean_gene_density_perMb = mean(gene_density_perMb, na.rm = TRUE),
    median_gene_density_perMb = median(gene_density_perMb, na.rm = TRUE),
    q25_gene_density_perMb = quantile(gene_density_perMb, 0.25, na.rm = TRUE),
    q75_gene_density_perMb = quantile(gene_density_perMb, 0.75, na.rm = TRUE)
  ), by = BGS_class][order(BGS_class)]

  kw <- kruskal.test(gene_density_perMb ~ BGS_class, data = as.data.frame(join_dt))
  pw <- pairwise.wilcox.test(join_dt$gene_density_perMb, join_dt$BGS_class, p.adjust.method = "BH")
  sub2 <- join_dt[BGS_class %in% c("sweep_like_excess", "BGS_like")]
  wil2 <- wilcox.test(gene_density_perMb ~ BGS_class, data = as.data.frame(sub2))
  eps <- max(min(sub2$recombination[sub2$recombination > 0], na.rm = TRUE) * 0.5, 1e-12)
  sub2[, log_rec := log10(recombination + eps)]
  sub2[, is_sweep_excess := as.integer(BGS_class == "sweep_like_excess")]
  mod <- glm(is_sweep_excess ~ scale(gene_density_perMb) + scale(log_rec),
             data = as.data.frame(sub2), family = binomial())

  # Clear class-level plots
  join_dt[, gene_n_cat := ifelse(gene_n >= 3, "3+ genes/window", paste0(gene_n, " gene/window"))]
  join_dt[!(gene_n_cat %in% c("0 gene/window","1 gene/window","2 gene/window","3+ genes/window")),
          gene_n_cat := "3+ genes/window"]
  join_dt[, gene_n_cat := factor(gene_n_cat, levels = c("0 gene/window","1 gene/window","2 gene/window","3+ genes/window"))]

  p1dat <- join_dt[, .N, by = .(BGS_class, gene_n_cat)]
  p1 <- ggplot(p1dat, aes(x = BGS_class, y = N, fill = gene_n_cat)) +
    geom_col(position = "fill", color = "white", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_brewer(palette = "YlGnBu", direction = 1) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(x = "BGS class", y = "Within-class proportion", fill = "Gene count category",
         title = "Gene density composition by BGS class")

  sumdt <- join_dt[, .(n = .N, mean_gene_n = mean(gene_n), se = sd(gene_n) / sqrt(.N)), by = BGS_class]
  sumdt[, `:=`(lo = mean_gene_n - 1.96 * se, hi = mean_gene_n + 1.96 * se, label = paste0("n=", n))]
  p2 <- ggplot(sumdt, aes(x = BGS_class, y = mean_gene_n)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12) +
    geom_text(aes(label = label, y = hi + 0.08), size = 3) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(x = "BGS class", y = "Mean genes per 10 kb window",
         title = "Mean gene count by BGS class (95% CI)")

  set.seed(20260318)
  sub_idx <- sample(seq_len(nrow(join_dt)), min(20000, nrow(join_dt)))
  p3 <- ggplot(join_dt[sub_idx], aes(x = BGS_class, y = gene_n, color = BGS_class)) +
    geom_jitter(width = 0.22, alpha = 0.15, size = 0.6, show.legend = FALSE) +
    geom_boxplot(width = 0.28, outlier.shape = NA, alpha = 0.2, color = "black", fill = "white") +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(x = "BGS class", y = "Genes per 10 kb window",
         title = "Per-window gene count by BGS class")

  # Chromosome-stratified tests
  sum_chr <- join_dt[, .(
    n = .N,
    mean_gene_n = mean(gene_n, na.rm = TRUE),
    med_gene_n = median(gene_n, na.rm = TRUE),
    mean_gene_density_perMb = mean(gene_density_perMb, na.rm = TRUE),
    med_gene_density_perMb = median(gene_density_perMb, na.rm = TRUE)
  ), by = .(Chromosome_Names, BGS_class)]

  chr_levels <- sort(unique(join_dt$Chromosome_Names))
  kw_list <- wb_list <- wo_list <- vector("list", length(chr_levels))
  for (i in seq_along(chr_levels)) {
    ch <- chr_levels[i]
    sub <- join_dt[Chromosome_Names == ch]
    cnt <- sub[, .N, by = BGS_class]
    valid_classes <- cnt[N >= 5, as.character(BGS_class)]
    subk <- sub[as.character(BGS_class) %in% valid_classes]
    if (length(unique(subk$BGS_class)) >= 2) {
      kwc <- kruskal.test(gene_n ~ BGS_class, data = as.data.frame(subk))
      kw_list[[i]] <- data.table(
        Chromosome_Names = ch, n_windows = nrow(subk),
        n_classes = length(unique(subk$BGS_class)),
        kw_statistic = as.numeric(kwc$statistic),
        kw_df = as.integer(kwc$parameter), kw_p = kwc$p.value
      )
    }
    sub_wb <- sub[as.character(BGS_class) %in% c("sweep_like_excess", "BGS_like")]
    if (length(unique(sub_wb$BGS_class)) == 2 && all(sub_wb[, .N, by = BGS_class]$N >= 5)) {
      wt <- wilcox.test(gene_n ~ BGS_class, data = as.data.frame(sub_wb))
      med <- sub_wb[, .(med = median(gene_n, na.rm = TRUE)), by = BGS_class]
      wb_list[[i]] <- data.table(
        Chromosome_Names = ch,
        n_sweep = sub_wb[BGS_class == "sweep_like_excess", .N],
        n_bgs = sub_wb[BGS_class == "BGS_like", .N],
        med_sweep = med[BGS_class == "sweep_like_excess", med],
        med_bgs = med[BGS_class == "BGS_like", med],
        delta_median = med[BGS_class == "sweep_like_excess", med] - med[BGS_class == "BGS_like", med],
        wilcox_p = wt$p.value
      )
    }
    sub_wo <- sub[as.character(BGS_class) %in% c("sweep_like_excess", "other")]
    if (length(unique(sub_wo$BGS_class)) == 2 && all(sub_wo[, .N, by = BGS_class]$N >= 5)) {
      wt2 <- wilcox.test(gene_n ~ BGS_class, data = as.data.frame(sub_wo))
      med2 <- sub_wo[, .(med = median(gene_n, na.rm = TRUE)), by = BGS_class]
      wo_list[[i]] <- data.table(
        Chromosome_Names = ch,
        n_sweep = sub_wo[BGS_class == "sweep_like_excess", .N],
        n_other = sub_wo[BGS_class == "other", .N],
        med_sweep = med2[BGS_class == "sweep_like_excess", med],
        med_other = med2[BGS_class == "other", med],
        delta_median = med2[BGS_class == "sweep_like_excess", med] - med2[BGS_class == "other", med],
        wilcox_p = wt2$p.value
      )
    }
  }
  kw_dt <- rbindlist(kw_list, fill = TRUE)
  wb_dt <- rbindlist(wb_list, fill = TRUE)
  wo_dt <- rbindlist(wo_list, fill = TRUE)
  if (nrow(kw_dt) > 0) kw_dt[, kw_p_adj_BH := p.adjust(kw_p, method = "BH")]
  if (nrow(wb_dt) > 0) wb_dt[, wilcox_p_adj_BH := p.adjust(wilcox_p, method = "BH")]
  if (nrow(wo_dt) > 0) wo_dt[, wilcox_p_adj_BH := p.adjust(wilcox_p, method = "BH")]

  ord <- sum_chr[, .(mean_gene_n = mean(mean_gene_n, na.rm = TRUE)), by = Chromosome_Names][order(mean_gene_n), Chromosome_Names]
  sum_chr[, Chromosome_Names := factor(Chromosome_Names, levels = ord)]
  p_heat <- ggplot(sum_chr, aes(x = BGS_class, y = Chromosome_Names, fill = mean_gene_n)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_viridis_c(option = "C") +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(x = "BGS class", y = "Chromosome", fill = "Mean genes/window",
         title = "Gene density by class across chromosomes")

  sum_chr_small <- sum_chr[n >= 5]
  p_lines <- ggplot(sum_chr_small, aes(x = BGS_class, y = mean_gene_n, group = Chromosome_Names)) +
    geom_line(alpha = 0.35, color = "grey40") +
    geom_point(aes(color = BGS_class), alpha = 0.75, size = 1.2) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(x = "BGS class", y = "Mean genes/window",
         title = "Per-chromosome class shift in gene density")

  p_delta <- NULL
  if (nrow(wb_dt) > 0) {
    wb_plot <- copy(wb_dt)
    wb_plot[, significant := wilcox_p_adj_BH < 0.05]
    wb_plot <- wb_plot[order(delta_median)]
    wb_plot[, Chromosome_Names := factor(Chromosome_Names, levels = Chromosome_Names)]
    p_delta <- ggplot(wb_plot, aes(x = Chromosome_Names, y = delta_median, color = significant)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
      geom_segment(aes(xend = Chromosome_Names, y = 0, yend = delta_median), color = "grey70", linewidth = 0.6) +
      geom_point(size = 2.6) +
      scale_color_manual(values = c("FALSE" = "grey35", "TRUE" = "#b2182b"),
                         labels = c("FALSE" = "Not significant", "TRUE" = "BH < 0.05")) +
      coord_flip() + theme_bw() + theme(panel.grid = element_blank()) +
      labs(x = "Chromosome", y = "Delta genes/window (sweep_like_excess - BGS_like)",
           color = NULL,
           title = "Chromosome-specific gene-density delta: sweep_like_excess vs BGS_like")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Main outputs
  data.table::fwrite(class_sum, file.path(out_dir, paste0(out_prefix, "_Association.tsv")), sep = "\t")
  data.table::fwrite(join_dt, file.path(out_dir, paste0(out_prefix, "_Windows.tsv")), sep = "\t")
  ggsave(file.path(out_dir, paste0(out_prefix, "_byClass_Violin.pdf")),
         ggplot(join_dt, aes(x = BGS_class, y = gene_density_perMb, fill = BGS_class)) +
           geom_violin(trim = TRUE, alpha = 0.8, color = "grey30", linewidth = 0.3) +
           stat_summary(fun = median, geom = "point", shape = 21, size = 2.2, fill = "white", color = "black") +
           theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
           labs(x = NULL, y = "Gene density (genes/Mb)", title = "Gene density by BGS class"),
         width = 8, height = 5)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ClassComposition_Stacked.pdf")), p1, width = 8.5, height = 5)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ClassMeanCI.pdf")), p2, width = 7.5, height = 5)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ClassJitterBox.pdf")), p3, width = 8.5, height = 5.5)

  # Chromosome outputs
  data.table::fwrite(sum_chr, file.path(out_dir, paste0(out_prefix, "_ByChromosome_Summary.tsv")), sep = "\t")
  data.table::fwrite(kw_dt, file.path(out_dir, paste0(out_prefix, "_ByChromosome_Kruskal.tsv")), sep = "\t")
  data.table::fwrite(wb_dt, file.path(out_dir, paste0(out_prefix, "_ByChromosome_SweepVsBGS_Wilcoxon.tsv")), sep = "\t")
  data.table::fwrite(wo_dt, file.path(out_dir, paste0(out_prefix, "_ByChromosome_SweepVsOther_Wilcoxon.tsv")), sep = "\t")
  ggsave(file.path(out_dir, paste0(out_prefix, "_ByChromosome_Heatmap.pdf")), p_heat, width = 8.5, height = 9.5)
  ggsave(file.path(out_dir, paste0(out_prefix, "_ByChromosome_Lines.pdf")), p_lines, width = 8.2, height = 5.2)
  if (!is.null(p_delta)) {
    ggsave(file.path(out_dir, paste0(out_prefix, "_ByChromosome_Delta_SweepMinusBGS.pdf")), p_delta, width = 8.6, height = 6.4)
  }

  report <- capture.output({
    cat("BGS vs Gene Density Report\n==========================\n\n")
    cat("Input windows:", nrow(join_dt), "\n\n")
    cat("Class summary:\n"); print(class_sum)
    cat("\nKruskal-Wallis:\n"); print(kw)
    cat("\nPairwise Wilcoxon (BH-adjusted):\n"); print(pw)
    cat("\nSweep_like_excess vs BGS_like Wilcoxon:\n"); print(wil2)
    cat("\nLogistic model (sweep_like_excess ~ gene_density + recombination):\n"); print(summary(mod))
    cat("\nChromosome-stratified top Kruskal (BH):\n")
    if (nrow(kw_dt) > 0) print(kw_dt[order(kw_p_adj_BH)][1:min(20, .N)])
  })
  writeLines(report, file.path(out_dir, paste0(out_prefix, "_Report.txt")))

  invisible(list(join_dt = join_dt, class_sum = class_sum, kw_dt = kw_dt, wb_dt = wb_dt))
}

if (exists("BGS_REPORT") && !is.null(BGS_REPORT) && exists("sd3_clusters")) {
  cat("\n================ [21] BGS CLASS VS GENE DENSITY =================\n")
  BGS_GENE_DENSITY <- tryCatch(
    make_bgs_gene_density_module(
      bgs_dt = BGS_REPORT$dt2,
      sd3_clusters = sd3_clusters,
      out_prefix = "BGS_GeneDensity",
      out_dir = "."
    ),
    error = function(e) {
      cat("[21] SKIP:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(BGS_GENE_DENSITY)) {
    cat("[21] Wrote: BGS_GeneDensity_* class and chromosome association outputs.\n")
  }
} else {
  cat("\n[21] SKIP: requires BGS_REPORT and sd3_clusters in memory.\n")
}

############################################################
## [22] SELECTION QC + EMPIRICAL-FDR PRIMARY CALLS
############################################################

selection_qc_helper <- file.path(script_dir, "selection_qc_and_fdr_helper.R")
if (file.exists(selection_qc_helper)) {
  cat("\n================ [22] SELECTION QC + EMPIRICAL-FDR PRIMARY CALLS =================\n")
  qc_status <- tryCatch(
    system2("Rscript", args = c(selection_qc_helper, getwd())),
    error = function(e) {
      cat("[22] QC helper failed to launch:", conditionMessage(e), "\n")
      NA_integer_
    }
  )
  if (identical(qc_status, 0L)) {
    cat("[22] Wrote QC folder, empirical-FDR selected windows, and updated barrier_label_objects.rdata.\n")
  } else {
    cat("[22] SKIP/WARN: QC helper returned status ", qc_status, "\n", sep = "")
  }
} else {
  cat("\n[22] SKIP: selection QC helper script not found at ", selection_qc_helper, "\n", sep = "")
}
