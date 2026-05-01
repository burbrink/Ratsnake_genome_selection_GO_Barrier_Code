#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})
script_file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file_arg)) dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

cfg <- list(
  empirical_alpha_csv = file.path(repo_root, "outputs", "barrier_analysis_outputs", "alpha_empirical_calibration_outputs", "empirical_alpha_estimates.csv"),
  two_scale_genome_csv = file.path(repo_root, "outputs", "barrier_analysis_outputs", "two_scale_barrier_outputs", "two_scale_barrier_genome_summary.csv"),
  two_scale_domains_csv = file.path(repo_root, "outputs", "barrier_analysis_outputs", "two_scale_barrier_outputs", "two_scale_barrier_extreme_domain_summary.csv"),
  two_scale_chr_csv = file.path(repo_root, "outputs", "barrier_analysis_outputs", "two_scale_barrier_outputs", "two_scale_barrier_extreme_domain_chromosomes.csv"),
  output_dir = file.path(repo_root, "outputs", "barrier_analysis_outputs", "two_scale_barrier_summary_outputs")
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) cfg$output_dir <- args[[1]]
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

alpha_tbl <- read.csv(cfg$empirical_alpha_csv, stringsAsFactors = FALSE, check.names = FALSE)
genome_tbl <- read.csv(cfg$two_scale_genome_csv, stringsAsFactors = FALSE, check.names = FALSE)
domain_tbl <- read.csv(cfg$two_scale_domains_csv, stringsAsFactors = FALSE, check.names = FALSE)
chr_tbl <- read.csv(cfg$two_scale_chr_csv, stringsAsFactors = FALSE, check.names = FALSE)

best_alpha <- alpha_tbl %>%
  filter(proxy == "compat_composite", method == "direct_mse") %>%
  slice_head(n = 1)

primary <- domain_tbl %>%
  filter(abs(threshold_prob - 0.995) < 1e-12) %>%
  slice_head(n = 1)
if (!nrow(primary)) primary <- domain_tbl %>% slice_head(n = 1)

top_chr <- chr_tbl %>% slice_head(n = 5)

summary_tbl <- tibble(
  alpha_local = best_alpha$alpha_estimate[[1]],
  alpha_proxy = best_alpha$proxy[[1]],
  alpha_method = best_alpha$method[[1]],
  K_genome_mean = genome_tbl$K_genome_mean[[1]],
  P_overall_mean_empirical = genome_tbl$P_overall_mean_empirical[[1]],
  log10_P_overall_mean_empirical = genome_tbl$log10_P_overall_mean_empirical[[1]],
  n_extreme_domains_primary = primary$n_domains[[1]],
  n_chromosomes_with_extreme_domains = primary$n_chromosomes_with_domains[[1]],
  total_extreme_bp = primary$total_domain_bp[[1]],
  max_domain_bp = primary$max_domain_length_bp[[1]]
)
write.csv(summary_tbl, file.path(cfg$output_dir, "two_scale_barrier_summary.csv"), row.names = FALSE)

caption_lines <- c(
  "Two-scale barrier interpretation for sd3_new_af",
  "===============================================",
  "",
  paste0("Empirical local alpha (composite proxy, direct fit): ", signif(best_alpha$alpha_estimate[[1]], 6)),
  paste0("Genome-wide mean K_local: ", signif(genome_tbl$K_genome_mean[[1]], 6)),
  paste0("Mean-scale compatibility under empirical alpha: ", signif(genome_tbl$P_overall_mean_empirical[[1]], 6),
         " (log10=", signif(genome_tbl$log10_P_overall_mean_empirical[[1]], 6), ")"),
  "",
  paste0("At the primary extreme threshold (top 0.5% of K_local; K >= ", signif(primary$K_cutoff[[1]], 6), "), the genome contains ",
         primary$n_domains[[1]], " independent extreme high-K domains across ", primary$n_chromosomes_with_domains[[1]],
         " chromosomes, spanning a total of ", signif(primary$total_domain_bp[[1]] / 1e6, 4), " Mb."),
  paste0("The largest extreme domain spans ", primary$max_domain_windows[[1]], " windows (",
         signif(primary$max_domain_length_bp[[1]] / 1e6, 4), " Mb)."),
  "",
  "Interpretation:",
  "The empirical local calibration indicates that average compatibility is not vanishingly small.",
  "However, species persistence need not depend on the genome-wide mean alone.",
  "Instead, the presence of multiple long, extremely high-K barrier blocks implies a chromosome-scale anti-collapse architecture.",
  "Under this view, local permeability remains mosaic, while a limited number of strong barrier domains may be sufficient to prevent genome-wide collapse.",
  "",
  "Top chromosomes by number of extreme high-K domains:"
)

for (i in seq_len(nrow(top_chr))) {
  row <- top_chr[i, ]
  caption_lines <- c(
    caption_lines,
    paste0("- ", row$Chromosome_Names, ": ", row$n_extreme_domains, " domains; max length ",
           signif(row$max_domain_length_bp / 1e6, 4), " Mb; total extreme coverage ",
           signif(row$total_extreme_bp / 1e6, 4), " Mb")
  )
}

caption_lines <- c(
  caption_lines,
  "",
  "Paper-style summary:",
  "The barrier landscape is best interpreted as a two-scale system: a locally permeable mosaic calibrated by alpha_local, superimposed on a smaller set of extreme barrier domains that likely contribute disproportionately to long-term maintenance of species differences."
)

writeLines(caption_lines, con = file.path(cfg$output_dir, "two_scale_barrier_summary.txt"))

p <- ggplot(chr_tbl, aes(x = reorder(Chromosome_Names, n_extreme_domains), y = n_extreme_domains)) +
  geom_col(fill = "#8c2d04") +
  coord_flip() +
  labs(
    title = "Chromosome Distribution of Extreme High-K Domains",
    subtitle = "Primary threshold: top 0.5% of K_local_mean",
    x = NULL,
    y = "Number of extreme domains"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(cfg$output_dir, "plot_two_scale_extreme_domains_by_chromosome.pdf"), p, width = 7.5, height = 5.5)

cat("\n================ Two-Scale Barrier Summary ================\n")
cat("Empirical alpha_local: ", signif(best_alpha$alpha_estimate[[1]], 6), "\n", sep = "")
cat("Mean-scale P_overall: ", signif(genome_tbl$P_overall_mean_empirical[[1]], 6), "\n", sep = "")
cat("Primary extreme threshold (top 0.5%): ", primary$n_domains[[1]], " independent extreme high-K areas across ",
    primary$n_chromosomes_with_domains[[1]], " chromosomes.\n", sep = "")
cat("Largest extreme block: ", primary$max_domain_windows[[1]], " windows (",
    signif(primary$max_domain_length_bp[[1]] / 1e6, 4), " Mb)\n", sep = "")
cat("Outputs written to: ", cfg$output_dir, "\n", sep = "")
