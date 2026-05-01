#!/usr/bin/env Rscript

# Connect the updated selection agreement tiers to the downstream GTF/GO/gene-map
# outputs. This keeps the GO and chromosome-gene summaries tied to the same
# handoff bundle used by the barrier workflow.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
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

output_dir <- if (length(args) >= 2 && nzchar(args[[2]])) {
  args[[2]]
} else {
  file.path(repo_root, "outputs", "barrier_analysis_outputs", "go_gene_map_outputs")
}

gtf_file <- if (length(args) >= 3 && nzchar(args[[3]])) {
  args[[3]]
} else {
  Sys.getenv("PIPELINE_GTF_FILE", unset = file.path(repo_root, "data", "braker_passqc_eggnog_genes.gtf"))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

label_rdata <- file.path(selection_output_dir, "barrier_label_objects.rdata")
if (!file.exists(label_rdata)) {
  stop("Missing selection handoff RData: ", label_rdata, call. = FALSE)
}
if (!file.exists(gtf_file)) {
  stop("Missing GTF file: ", gtf_file, call. = FALSE)
}

message("Loading selection handoff: ", label_rdata)
load(label_rdata)
if (!exists("sd3_new_af2") || !exists("agree_dat")) {
  stop("Expected objects sd3_new_af2 and agree_dat in ", label_rdata, call. = FALSE)
}

sd3 <- as.data.table(sd3_new_af2)
agd <- as.data.table(agree_dat)

required_sd3 <- c("Chromosome_Names", ".start_bp", ".end_bp", "Chromosomes")
missing_sd3 <- setdiff(required_sd3, names(sd3))
if (length(missing_sd3)) {
  stop("sd3_new_af2 is missing required columns: ", paste(missing_sd3, collapse = ", "), call. = FALSE)
}

required_agd <- c(
  "Chromosome_Names", ".start_bp", ".end_bp", "gmm_label", "rf_label",
  "cont_label_h", "cont_label_top1pct_h", "cont_label_tail500_h",
  "gmm_rf_match", "triple_match", "triple_top1pct", "triple_tail500",
  "triple_margin_025", "triple_margin_050", "triple_margin_100",
  "continuous_distance_label_empirical_h", "contdist_match",
  "dist_margin_025", "dist_margin_050", "dist_margin_100",
  "score_margin_025", "score_margin_050", "score_margin_100",
  "dist_top1pct", "dist_top500",
  "empirical_fdr_selected_label", "global_shift_selected_label",
  "genomewide_fdr_selected_label", "broad_agreement_selected_label",
  "interpretive_selection_class"
)
missing_agd <- setdiff(required_agd, names(agd))
if (length(missing_agd)) {
  stop("agree_dat is missing required columns: ", paste(missing_agd, collapse = ", "), call. = FALSE)
}

sd3[, `:=`(Chromosome_Names = as.character(Chromosome_Names), .start_bp = as.integer(.start_bp))]
agd[, `:=`(Chromosome_Names = as.character(Chromosome_Names), .start_bp = as.integer(.start_bp))]
setkey(sd3, Chromosome_Names, .start_bp)
setkey(agd, Chromosome_Names, .start_bp)

dat <- sd3[agd, on = .(Chromosome_Names, .start_bp), nomatch = 0L]
if ("Chromosomes.x" %in% names(dat) && !"Chromosomes" %in% names(dat)) {
  setnames(dat, "Chromosomes.x", "Chromosomes")
}
if (".end_bp.x" %in% names(dat) && !".end_bp" %in% names(dat)) {
  setnames(dat, ".end_bp.x", ".end_bp")
}

selection_classes <- c("Geographic sweep", "Balancing selection")
class_to_short <- c("Geographic sweep" = "sweep", "Balancing selection" = "balancing")

ordered_chromosome_levels <- function(x) {
  x <- unique(as.character(x))
  chr_num <- suppressWarnings(as.integer(sub("^Chromosome_([0-9]+)$", "\\1", x)))
  out <- x[order(
    ifelse(is.na(chr_num), Inf, chr_num),
    x == "Chromosome_ZW",
    x
  )]
  if ("Chromosome_ZW" %in% out) {
    out <- c(setdiff(out, "Chromosome_ZW"), "Chromosome_ZW")
  }
  out
}

sample_matched_window_keys <- function(windows_dt, background_dt) {
  req <- windows_dt[, .N, by = Chromosome_Names]
  out <- vector("list", nrow(req))
  for (i in seq_len(nrow(req))) {
    chr_i <- req$Chromosome_Names[[i]]
    n_i <- req$N[[i]]
    pool <- background_dt[Chromosome_Names == chr_i, window_key]
    if (!length(pool)) next
    if (length(pool) >= n_i) {
      out[[i]] <- sample(pool, n_i, replace = FALSE)
    } else {
      out[[i]] <- sample(pool, n_i, replace = TRUE)
    }
  }
  unique(unlist(out, use.names = FALSE))
}

sample_matched_genes <- function(counts_dt, background_genes_dt) {
  out <- vector("list", nrow(counts_dt))
  for (i in seq_len(nrow(counts_dt))) {
    chr_i <- counts_dt$chr[[i]]
    n_i <- counts_dt$N[[i]]
    pool <- background_genes_dt[chr == chr_i & !is.na(eggnog_gene_name) & eggnog_gene_name != "-", unique(eggnog_gene_name)]
    if (!length(pool)) next
    if (length(pool) >= n_i) {
      out[[i]] <- sample(pool, n_i, replace = FALSE)
    } else {
      out[[i]] <- sample(pool, n_i, replace = TRUE)
    }
  }
  sort(unique(unlist(out, use.names = FALSE)))
}

chromosome_lengths_dt <- sd3[, .(chr_len = max(.end_bp, na.rm = TRUE)), by = Chromosome_Names]
chromosome_lengths_dt[, Chromosome_Names := as.character(Chromosome_Names)]

full_class_tier_specs <- data.table(
  tier = c("gmm_rf", "continuous_distance_3class", "gmm_rf_contdist"),
  label_col = c("gmm_rf_label", "continuous_distance_3class_label", "gmm_rf_contdist_label")
)

tier_specs <- data.table(
  tier = c(
    "gmm_rf", "best_continuous_3of3", "continuous_distance_3class", "gmm_rf_contdist",
    "empirical_fdr_selected",
    "genomewide_fdr_selected",
    "global_shift_selected",
    "dist_margin_025", "dist_margin_050", "dist_margin_100",
    "score_margin_025", "score_margin_050", "score_margin_100",
    "dist_top1pct", "dist_top500",
    "margin_025", "margin_050", "margin_100", "top1pct", "tail500"
  ),
  flag_col = c(
    "gmm_rf_match", "triple_match", "__direct__", "contdist_match",
    "__direct__",
    "__direct__",
    "__direct__",
    "dist_margin_025", "dist_margin_050", "dist_margin_100",
    "score_margin_025", "score_margin_050", "score_margin_100",
    "dist_top1pct", "dist_top500",
    "triple_margin_025", "triple_margin_050", "triple_margin_100", "triple_top1pct", "triple_tail500"
  ),
  label_col = c(
    "gmm_label", "cont_label_h", "continuous_distance_label_empirical_h", "rf_label",
    "empirical_fdr_selected_label",
    "genomewide_fdr_selected_label",
    "global_shift_selected_label",
    "rf_label", "rf_label", "rf_label",
    "rf_label", "rf_label", "rf_label",
    "rf_label", "rf_label",
    "cont_label_h", "cont_label_h", "cont_label_h", "cont_label_top1pct_h", "cont_label_tail500_h"
  )
)

get_tier_windows <- function(tier, flag_col, label_col) {
  x <- copy(if (identical(flag_col, "__direct__")) {
    dat[get(label_col) %in% selection_classes]
  } else {
    dat[get(flag_col) %in% TRUE & get(label_col) %in% selection_classes]
  })
  if (!nrow(x)) {
    return(data.table())
  }
  x[, `:=`(
    tier = tier,
    selection_class = get(label_col),
    selection_short = unname(class_to_short[get(label_col)]),
    chr = as.character(Chromosomes),
    window_start = as.integer(.start_bp),
    window_end = as.integer(.end_bp),
    window_midpoint = as.numeric((.start_bp + .end_bp) / 2)
  )]
  x[, .(
    tier, selection_class, selection_short,
    Chromosome_Names, chr, window_start, window_end, window_midpoint,
    gmm_label, rf_label, cont_label_h, cont_label_top1pct_h, cont_label_tail500_h,
    continuous_distance_label_empirical_h,
    gmm_rf_match, triple_match, triple_margin_025, triple_margin_050,
    triple_margin_100, triple_top1pct, triple_tail500,
    contdist_match, dist_margin_025, dist_margin_050, dist_margin_100,
    score_margin_025, score_margin_050, score_margin_100,
    dist_top1pct, dist_top500
  )]
}

selected_windows <- rbindlist(
  lapply(seq_len(nrow(tier_specs)), function(i) {
    get_tier_windows(tier_specs$tier[[i]], tier_specs$flag_col[[i]], tier_specs$label_col[[i]])
  }),
  use.names = TRUE,
  fill = TRUE
)

selected_windows[, window_key := paste(Chromosome_Names, window_start, window_end, sep = ":")]

fwrite(selected_windows, file.path(output_dir, "selected_windows_by_tier.tsv"), sep = "\t")

get_full_class_windows <- function(tier_name, label_col) {
  x <- copy(dat)
  x[, tier_label := as.character(get(label_col))]
  x[is.na(tier_label) | !nzchar(tier_label), tier_label := "No consensus"]
  x[, class_short := fcase(
    tier_label == "Geographic sweep", "sweep",
    tier_label == "Balancing selection", "balancing",
    tier_label == "Neutral / equilibrium", "neutral",
    tier_label == "Ambiguous", "ambiguous",
    default = "no_consensus"
  )]
  x[, .(
    tier = tier_name,
    tier_label,
    class_short,
    Chromosome_Names,
    chr = as.character(Chromosomes),
    window_start = as.integer(.start_bp),
    window_end = as.integer(.end_bp),
    window_midpoint = as.numeric((.start_bp + .end_bp) / 2)
  )]
}

full_class_windows <- rbindlist(
  lapply(seq_len(nrow(full_class_tier_specs)), function(i) {
    get_full_class_windows(full_class_tier_specs$tier[[i]], full_class_tier_specs$label_col[[i]])
  }),
  use.names = TRUE,
  fill = TRUE
)

all_windows_base <- unique(dat[, .(
  Chromosome_Names,
  chr = as.character(Chromosomes),
  window_start = as.integer(.start_bp),
  window_end = as.integer(.end_bp),
  window_midpoint = as.numeric((.start_bp + .end_bp) / 2)
)])
all_windows_base[, window_key := paste(Chromosome_Names, window_start, window_end, sep = ":")]

message("Importing GTF/CDS features: ", gtf_file)
gtf <- import(gtf_file)
cds_feat <- gtf[gtf$type == "CDS"]
cds_gene_dt <- unique(data.table(
  chr = as.character(seqnames(cds_feat)),
  gene_id = as.character(mcols(cds_feat)$gene_id),
  transcript_id = as.character(mcols(cds_feat)$transcript_id),
  eggnog_gene_name = as.character(mcols(cds_feat)$eggnog_gene_name)
))

get_hits_tbl <- function(windows_dt) {
  if (!nrow(windows_dt)) {
    return(data.table())
  }
  region <- GRanges(
    seqnames = windows_dt$chr,
    ranges = IRanges(start = windows_dt$window_start, end = windows_dt$window_end),
    tier = windows_dt$tier,
    selection_class = windows_dt$selection_class,
    selection_short = windows_dt$selection_short,
    Chromosome_Names = windows_dt$Chromosome_Names,
    window_start = windows_dt$window_start,
    window_end = windows_dt$window_end,
    window_midpoint = windows_dt$window_midpoint
  )
  hits <- findOverlaps(cds_feat, region, ignore.strand = TRUE)
  if (!length(hits)) {
    return(data.table())
  }
  out <- data.table(
    tier = mcols(region)$tier[subjectHits(hits)],
    class = mcols(region)$selection_short[subjectHits(hits)],
    selection_class = mcols(region)$selection_class[subjectHits(hits)],
    chr = as.character(seqnames(region)[subjectHits(hits)]),
    Chromosome_Names = mcols(region)$Chromosome_Names[subjectHits(hits)],
    window_start = mcols(region)$window_start[subjectHits(hits)],
    window_end = mcols(region)$window_end[subjectHits(hits)],
    window_midpoint = mcols(region)$window_midpoint[subjectHits(hits)],
    window_key = paste(
      mcols(region)$Chromosome_Names[subjectHits(hits)],
      mcols(region)$window_start[subjectHits(hits)],
      mcols(region)$window_end[subjectHits(hits)],
      sep = ":"
    ),
    gene_id = as.character(mcols(cds_feat)$gene_id[queryHits(hits)]),
    transcript_id = as.character(mcols(cds_feat)$transcript_id[queryHits(hits)]),
    eggnog_gene_name = as.character(mcols(cds_feat)$eggnog_gene_name[queryHits(hits)])
  )
  unique(out)
}

selected_genes <- get_hits_tbl(selected_windows)
fwrite(selected_genes, file.path(output_dir, "selected_genes_from_windows_by_tier.tsv"), sep = "\t")
all_window_hits <- get_hits_tbl(
  all_windows_base[, `:=`(tier = "background", selection_class = "background", selection_short = "background")]
)
all_window_gene_summary <- all_window_hits[
  ,
  .(
    overlapping_cds_rows = .N,
    unique_gene_ids = uniqueN(gene_id),
    unique_named_genes = uniqueN(eggnog_gene_name[!is.na(eggnog_gene_name) & eggnog_gene_name != "-"]),
    mean_cds_overlap_length = NA_real_
  ),
  by = .(window_key, Chromosome_Names)
]

get_hits_tbl_full <- function(windows_dt) {
  if (!nrow(windows_dt)) return(data.table())
  region <- GRanges(
    seqnames = windows_dt$chr,
    ranges = IRanges(start = windows_dt$window_start, end = windows_dt$window_end),
    tier = windows_dt$tier,
    tier_label = windows_dt$tier_label,
    class_short = windows_dt$class_short,
    Chromosome_Names = windows_dt$Chromosome_Names,
    window_start = windows_dt$window_start,
    window_end = windows_dt$window_end,
    window_midpoint = windows_dt$window_midpoint
  )
  hits <- findOverlaps(cds_feat, region, ignore.strand = TRUE)
  if (!length(hits)) return(data.table())
  out <- data.table(
    tier = mcols(region)$tier[subjectHits(hits)],
    tier_label = mcols(region)$tier_label[subjectHits(hits)],
    class_short = mcols(region)$class_short[subjectHits(hits)],
    chr = as.character(seqnames(region)[subjectHits(hits)]),
    Chromosome_Names = mcols(region)$Chromosome_Names[subjectHits(hits)],
    window_start = mcols(region)$window_start[subjectHits(hits)],
    window_end = mcols(region)$window_end[subjectHits(hits)],
    window_midpoint = mcols(region)$window_midpoint[subjectHits(hits)],
    gene_id = as.character(mcols(cds_feat)$gene_id[queryHits(hits)]),
    transcript_id = as.character(mcols(cds_feat)$transcript_id[queryHits(hits)]),
    eggnog_gene_name = as.character(mcols(cds_feat)$eggnog_gene_name[queryHits(hits)])
  )
  unique(out)
}

full_class_genes <- get_hits_tbl_full(full_class_windows)

for (tier_name in tier_specs$tier) {
  tier_genes <- selected_genes[tier == tier_name]
  fwrite(tier_genes[class == "sweep"], file.path(output_dir, sprintf("sweep_genes_from_windows_%s.tsv", tier_name)), sep = "\t")
  fwrite(tier_genes[class == "balancing"], file.path(output_dir, sprintf("balancing_genes_from_windows_%s.tsv", tier_name)), sep = "\t")
  fwrite(tier_genes, file.path(output_dir, sprintf("sweep_balancing_genes_from_windows_%s.tsv", tier_name)), sep = "\t")
}

if ("gmm_rf" %in% selected_genes$tier) {
  fwrite(selected_genes[tier == "gmm_rf" & class == "sweep"], file.path(output_dir, "sweep_genes_from_windows.tsv"), sep = "\t")
  fwrite(selected_genes[tier == "gmm_rf" & class == "balancing"], file.path(output_dir, "balancing_genes_from_windows.tsv"), sep = "\t")
  fwrite(selected_genes[tier == "gmm_rf"], file.path(output_dir, "sweep_balancing_genes_from_windows.tsv"), sep = "\t")
}

format_gost <- function(res, tier_name, class_name, subset_name) {
  if (is.null(res) || is.null(res$result) || !nrow(res$result)) {
    return(data.table(
      tier = tier_name, class = class_name, subset = subset_name,
      source = character(), term_id = character(), term_name = character(),
      p_value = numeric(), intersection_size = integer(), term_size = integer(),
      effective_domain_size = integer()
    ))
  }
  as.data.table(res$result)[order(p_value), .(
    tier = tier_name,
    class = class_name,
    subset = subset_name,
    source,
    term_id,
    term_name,
    p_value,
    intersection_size,
    term_size,
    effective_domain_size
  )]
}

run_gost_safe <- function(genes) {
  genes <- sort(unique(genes[!is.na(genes) & genes != "-" & nzchar(genes)]))
  if (length(genes) < 2L || !requireNamespace("gprofiler2", quietly = TRUE)) {
    return(NULL)
  }
  tryCatch(
    gprofiler2::gost(query = genes, organism = "hsapiens", sources = c("GO:BP", "GO:MF", "GO:CC")),
    error = function(e) {
      message("gProfiler query failed: ", conditionMessage(e))
      NULL
    }
  )
}

go_rows <- list()
go_i <- 1L
for (tier_name in c("gmm_rf", "continuous_distance_3class", "gmm_rf_contdist", "empirical_fdr_selected", "genomewide_fdr_selected", "global_shift_selected", "dist_margin_050", "dist_top1pct", "dist_top500")) {
  for (class_name in c("sweep", "balancing")) {
    genes <- selected_genes[tier == tier_name & class == class_name, eggnog_gene_name]
    go_rows[[go_i]] <- format_gost(run_gost_safe(genes), tier_name, class_name, "whole_genome")
    go_i <- go_i + 1L
  }
}

zw_chr <- "NC_045558.1_RagTag"
for (class_name in c("sweep", "balancing")) {
  genes <- selected_genes[tier == "gmm_rf" & class == class_name & chr == zw_chr, eggnog_gene_name]
  go_rows[[go_i]] <- format_gost(run_gost_safe(genes), "gmm_rf", class_name, "ZW")
  go_i <- go_i + 1L
}

go_results <- rbindlist(go_rows, use.names = TRUE, fill = TRUE)
fwrite(go_results, file.path(output_dir, "GO_gprofiler_results_by_tier.tsv"), sep = "\t")
fwrite(go_results[subset == "ZW"], file.path(output_dir, "ZW_GO_gprofiler_results.tsv"), sep = "\t")

immune_pattern <- paste(
  c(
    "^HLA", "^MR1", "^KLR", "^CLEC", "^IL[0-9]", "^IL1R", "^IL10R", "^IL12R",
    "^IL13R", "^IL17R", "^IL18R", "^IL20R", "^IL23R", "^IFN", "^IFNGR",
    "^TNFRSF", "^TNFSF", "^TNFAIP", "^NFKB", "^IKB", "^C2$", "^C3$",
    "^C4A$", "^C5$", "^C6$", "^C7$", "^C8A$", "^C9$"
  ),
  collapse = "|"
)

immune_go_keyword_pattern <- "(immune|antigen|mhc|interleukin|interferon|t cell|b cell|cytokine|chemokine|complement|defense|inflamm|leukocyte|lymphocyte|pathogen|viral|bacterial)"

gene_universe <- unique(cds_gene_dt[!is.na(gene_id) & nzchar(gene_id), .(gene_id, eggnog_gene_name)])
gene_universe[, immune_like := grepl(immune_pattern, eggnog_gene_name, ignore.case = TRUE)]

fisher_rows <- rbindlist(lapply(tier_specs$tier, function(tier_name) {
  rbindlist(lapply(c("sweep", "balancing"), function(class_name) {
    selected_ids <- unique(selected_genes[tier == tier_name & class == class_name, gene_id])
    in_set <- gene_universe$gene_id %in% selected_ids
    tab <- table(selected = in_set, immune_like = gene_universe$immune_like)
    pval <- odds <- NA_real_
    if (all(dim(tab) == c(2L, 2L))) {
      ft <- fisher.test(tab)
      pval <- ft$p.value
      odds <- unname(ft$estimate)
    }
    data.table(
      tier = tier_name,
      class = class_name,
      selected_genes = sum(in_set),
      immune_selected_genes = sum(in_set & gene_universe$immune_like),
      universe_genes = nrow(gene_universe),
      immune_universe_genes = sum(gene_universe$immune_like),
      fisher_odds_ratio = odds,
      fisher_p_value = pval
    )
  }))
}), use.names = TRUE)
fwrite(fisher_rows, file.path(output_dir, "immune_gene_enrichment_fisher_by_tier.tsv"), sep = "\t")

immune_go_terms <- go_results[
  grepl(immune_go_keyword_pattern, term_name, ignore.case = TRUE),
  .(tier, class, subset, source, term_id, term_name, p_value, intersection_size, term_size, effective_domain_size)
][order(tier, class, subset, p_value)]
fwrite(immune_go_terms, file.path(output_dir, "immune_related_GO_terms_by_tier.tsv"), sep = "\t")

immune_named_genes <- sort(unique(gene_universe[immune_like == TRUE & !is.na(eggnog_gene_name) & eggnog_gene_name != "-", eggnog_gene_name]))
immune_gene_sets_by_class <- full_class_genes[
  eggnog_gene_name %in% immune_named_genes,
  .(
    immune_gene_rows = .N,
    unique_immune_genes = uniqueN(gene_id),
    unique_immune_named_genes = uniqueN(eggnog_gene_name)
  ),
  by = .(tier, tier_label, class_short)
][order(tier, tier_label)]
fwrite(immune_gene_sets_by_class, file.path(output_dir, "immune_gene_counts_fullclass_by_tier.tsv"), sep = "\t")

immune_windows_by_class <- full_class_genes[
  ,
  .(
    has_immune_gene = any(eggnog_gene_name %in% immune_named_genes),
    immune_named_gene_n = uniqueN(eggnog_gene_name[eggnog_gene_name %in% immune_named_genes])
  ),
  by = .(tier, tier_label, class_short, Chromosome_Names, window_start, window_end, window_midpoint)
]

full_window_counts <- full_class_windows[
  ,
  .N,
  by = .(tier, tier_label, class_short)
]
setnames(full_window_counts, "N", "n_windows_total")

immune_window_summary <- merge(
  full_window_counts,
  immune_windows_by_class[
    ,
    .(
      n_windows_with_immune_gene = sum(has_immune_gene, na.rm = TRUE),
      mean_immune_named_gene_n = mean(immune_named_gene_n, na.rm = TRUE),
      median_immune_named_gene_n = median(immune_named_gene_n, na.rm = TRUE)
    ),
    by = .(tier, tier_label, class_short)
  ],
  by = c("tier", "tier_label", "class_short"),
  all.x = TRUE
)
immune_window_summary[is.na(n_windows_with_immune_gene), n_windows_with_immune_gene := 0L]
immune_window_summary[, frac_windows_with_immune_gene := fifelse(n_windows_total > 0, n_windows_with_immune_gene / n_windows_total, NA_real_)]
fwrite(immune_window_summary, file.path(output_dir, "immune_window_summary_fullclass_by_tier.tsv"), sep = "\t")

compare_classes_window_fisher <- function(dt, tier_name, class_a, class_b) {
  xa <- dt[tier == tier_name & class_short == class_a]
  xb <- dt[tier == tier_name & class_short == class_b]
  if (!nrow(xa) || !nrow(xb)) {
    return(data.table(
      tier = tier_name,
      class_a = class_a,
      class_b = class_b,
      n_a = nrow(xa),
      n_b = nrow(xb),
      immune_windows_a = sum(xa$has_immune_gene, na.rm = TRUE),
      immune_windows_b = sum(xb$has_immune_gene, na.rm = TRUE),
      odds_ratio = NA_real_,
      p_value = NA_real_
    ))
  }
  mat <- matrix(
    c(
      sum(xa$has_immune_gene, na.rm = TRUE),
      nrow(xa) - sum(xa$has_immune_gene, na.rm = TRUE),
      sum(xb$has_immune_gene, na.rm = TRUE),
      nrow(xb) - sum(xb$has_immune_gene, na.rm = TRUE)
    ),
    nrow = 2,
    byrow = TRUE
  )
  ft <- fisher.test(mat)
  data.table(
    tier = tier_name,
    class_a = class_a,
    class_b = class_b,
    n_a = nrow(xa),
    n_b = nrow(xb),
    immune_windows_a = sum(xa$has_immune_gene, na.rm = TRUE),
    immune_windows_b = sum(xb$has_immune_gene, na.rm = TRUE),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

immune_window_class_comparisons <- rbindlist(lapply(full_class_tier_specs$tier, function(tier_name) {
  rbindlist(list(
    compare_classes_window_fisher(immune_windows_by_class, tier_name, "balancing", "neutral"),
    compare_classes_window_fisher(immune_windows_by_class, tier_name, "balancing", "sweep")
  ), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)
fwrite(immune_window_class_comparisons, file.path(output_dir, "immune_window_class_comparisons.tsv"), sep = "\t")

plot_immune_fraction <- function() {
  x <- copy(immune_window_summary[class_short %in% c("balancing", "neutral", "sweep")])
  if (!nrow(x)) return(invisible(NULL))
  x[, class_label := factor(
    class_short,
    levels = c("balancing", "neutral", "sweep"),
    labels = c("Balancing", "Neutral", "Sweep")
  )]
  x[, tier := factor(tier, levels = c("gmm_rf", "continuous_distance_3class", "gmm_rf_contdist"))]
  p <- ggplot(x, aes(class_label, frac_windows_with_immune_gene, fill = class_label)) +
    geom_col(width = 0.72, color = "white") +
    facet_wrap(~ tier, scales = "free_y") +
    scale_fill_manual(values = c("Balancing" = "#2A7A68", "Neutral" = "grey65", "Sweep" = "#C64A2E")) +
    labs(
      x = NULL,
      y = "Fraction of windows with immune-like genes",
      title = "Immune-like gene content across balancing, neutral, and sweep windows"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(), legend.position = "none")
  ggsave(file.path(output_dir, "immune_window_fraction_by_class.pdf"), p, width = 10, height = 5.5)
}

plot_immune_or <- function() {
  x <- copy(immune_window_class_comparisons[is.finite(odds_ratio)])
  if (!nrow(x)) return(invisible(NULL))
  x[, comparison := factor(
    paste0("Balancing vs ", class_b),
    levels = c("Balancing vs neutral", "Balancing vs sweep")
  )]
  x[, tier := factor(tier, levels = c("gmm_rf", "continuous_distance_3class", "gmm_rf_contdist"))]
  p <- ggplot(x, aes(tier, odds_ratio, color = comparison, group = comparison)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
    geom_point(size = 2.8) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Balancing vs neutral" = "#1B5E20", "Balancing vs sweep" = "#8E3B1F")) +
    labs(
      x = NULL,
      y = "Odds ratio for immune-like windows",
      color = NULL,
      title = "Balancing windows are more likely to contain immune-like genes"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(), legend.position = "bottom")
  ggsave(file.path(output_dir, "immune_window_odds_ratio_comparisons.pdf"), p, width = 9.5, height = 5.5)
}

plot_balancing_immune_chromosomes <- function(tier_name) {
  x <- full_class_genes[
    tier == tier_name &
      class_short == "balancing" &
      eggnog_gene_name %in% immune_named_genes,
    .(immune_n = uniqueN(gene_id)),
    by = .(Chromosome_Names, window_start, window_end, window_midpoint)
  ]
  if (!nrow(x)) return(invisible(NULL))
  chr_levels <- ordered_chromosome_levels(x$Chromosome_Names)
  chr_df <- chromosome_lengths_dt[Chromosome_Names %in% chr_levels]
  chr_df[, Chromosome_Names := factor(Chromosome_Names, levels = chr_levels)]
  chr_df <- chr_df[order(Chromosome_Names)]
  chr_df[, y := seq_len(.N)]
  chr_df[, Chromosome_Names := factor(as.character(Chromosome_Names), levels = chr_levels)]
  x <- merge(x, chr_df[, .(Chromosome_Names, y)], by = "Chromosome_Names", all.x = TRUE)
  x[, Chromosome_Names := factor(Chromosome_Names, levels = chr_levels)]
  x[, mid_Mb := window_midpoint / 1e6]
  tick_half <- 0.24

  p_points <- ggplot() +
    geom_segment(
      data = chr_df,
      aes(x = 0, xend = chr_len / 1e6, y = y, yend = y),
      linewidth = 2, color = "grey90"
    ) +
    geom_linerange(
      data = x,
      aes(x = mid_Mb, ymin = y - tick_half, ymax = y + tick_half),
      linewidth = 0.7, alpha = 0.8, color = "#1f77b4"
    ) +
    geom_point(
      data = x,
      aes(x = mid_Mb, y = y, size = immune_n),
      alpha = 0.85, color = "#1f77b4"
    ) +
    scale_y_reverse(breaks = chr_df$y, labels = as.character(chr_df$Chromosome_Names)) +
    scale_size_continuous(range = c(2, 6.8), breaks = 1:4, name = "# immune genes\nin window") +
    labs(
      x = "Position (Mb)",
      y = "Chromosome/linkage group",
      title = paste("Balancing windows with immune-like genes:", tier_name)
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
  ggsave(file.path(output_dir, sprintf("balancing_immune_windows_%s_points.pdf", tier_name)), p_points, width = 12, height = 7)

  x[, immune_cat := factor(pmin(immune_n, 4L), levels = 1:4, labels = c("1", "2", "3", "4+"))]
  p_ticks <- ggplot() +
    geom_segment(
      data = chr_df,
      aes(x = 0, xend = chr_len / 1e6, y = y, yend = y),
      linewidth = 2, color = "grey90"
    ) +
    geom_linerange(
      data = x,
      aes(x = mid_Mb, ymin = y - tick_half, ymax = y + tick_half, color = immune_cat),
      linewidth = 1.1, alpha = 0.9
    ) +
    scale_y_reverse(breaks = chr_df$y, labels = as.character(chr_df$Chromosome_Names)) +
    scale_color_brewer(palette = "Set1", name = "# immune genes\nin window") +
    labs(
      x = "Position (Mb)",
      y = "Chromosome/linkage group",
      title = paste("Balancing windows with immune-like genes:", tier_name)
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
  ggsave(file.path(output_dir, sprintf("balancing_immune_windows_%s_ticks.pdf", tier_name)), p_ticks, width = 12, height = 7)
}

plot_immune_fraction()
plot_immune_or()
for (tier_name in full_class_tier_specs$tier) {
  plot_balancing_immune_chromosomes(tier_name)
}

plot_windows <- function(tier_name) {
  x <- selected_windows[tier == tier_name]
  if (!nrow(x)) return(invisible(NULL))
  chr_levels <- ordered_chromosome_levels(x$Chromosome_Names)
  chr_df <- chromosome_lengths_dt[Chromosome_Names %in% chr_levels]
  chr_df[, Chromosome_Names := factor(Chromosome_Names, levels = chr_levels)]
  chr_df <- chr_df[order(Chromosome_Names)]
  chr_df[, y := seq_len(.N)]
  x <- merge(x, chr_df[, .(Chromosome_Names, y)], by = "Chromosome_Names", all.x = TRUE)
  x[, mid_Mb := window_midpoint / 1e6]
  tick_half <- 0.22
  p <- ggplot() +
    geom_segment(
      data = chr_df,
      aes(x = 0, xend = chr_len / 1e6, y = y, yend = y),
      linewidth = 2, color = "grey88"
    ) +
    geom_linerange(
      data = x,
      aes(x = mid_Mb, ymin = y - tick_half, ymax = y + tick_half, color = selection_class),
      linewidth = 0.95, alpha = 0.85
    ) +
    scale_y_reverse(breaks = chr_df$y, labels = as.character(chr_df$Chromosome_Names)) +
    scale_color_manual(values = c("Geographic sweep" = "#C64A2E", "Balancing selection" = "#2A7A68")) +
    labs(
      x = "Position (Mb)",
      y = "Chromosome/linkage group",
      color = "Selection class",
      title = paste("Selected windows:", tier_name)
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
  ggsave(file.path(output_dir, sprintf("chromosome_selected_windows_%s.pdf", tier_name)), p, width = 12, height = 8)
}

plot_genes <- function(tier_name) {
  x <- unique(selected_genes[tier == tier_name & !is.na(window_midpoint), .(
    tier, class, selection_class, Chromosome_Names, window_midpoint, gene_id
  )])
  if (!nrow(x)) return(invisible(NULL))
  chr_levels <- ordered_chromosome_levels(x$Chromosome_Names)
  chr_df <- chromosome_lengths_dt[Chromosome_Names %in% chr_levels]
  chr_df[, Chromosome_Names := factor(Chromosome_Names, levels = chr_levels)]
  chr_df <- chr_df[order(Chromosome_Names)]
  chr_df[, y := seq_len(.N)]
  x <- merge(x, chr_df[, .(Chromosome_Names, y)], by = "Chromosome_Names", all.x = TRUE)
  x[, mid_Mb := window_midpoint / 1e6]
  tick_half <- 0.18
  p <- ggplot() +
    geom_segment(
      data = chr_df,
      aes(x = 0, xend = chr_len / 1e6, y = y, yend = y),
      linewidth = 2, color = "grey88"
    ) +
    geom_linerange(
      data = x,
      aes(x = mid_Mb, ymin = y - tick_half, ymax = y + tick_half, color = selection_class),
      linewidth = 0.8, alpha = 0.8
    ) +
    scale_y_reverse(breaks = chr_df$y, labels = as.character(chr_df$Chromosome_Names)) +
    scale_color_manual(values = c("Geographic sweep" = "#C64A2E", "Balancing selection" = "#2A7A68")) +
    labs(
      x = "Position (Mb)",
      y = "Chromosome/linkage group",
      color = "Selection class",
      title = paste("Genes overlapping selected windows:", tier_name)
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
  ggsave(file.path(output_dir, sprintf("chromosome_selected_genes_%s.pdf", tier_name)), p, width = 12, height = 8)
}

for (tier_name in c("gmm_rf", "continuous_distance_3class", "gmm_rf_contdist", "empirical_fdr_selected", "genomewide_fdr_selected", "global_shift_selected", "dist_margin_050", "dist_top1pct", "dist_top500")) {
  plot_windows(tier_name)
  plot_genes(tier_name)
}

null_iter <- if (identical(Sys.getenv("PIPELINE_FAST_TEST"), "1")) 25L else 100L
gene_overlap_null <- rbindlist(lapply(unique(selected_windows$tier), function(tier_name) {
  rbindlist(lapply(c("balancing", "sweep"), function(class_name) {
    obs_windows <- unique(selected_windows[tier == tier_name & selection_short == class_name, .(window_key, Chromosome_Names)])
    if (!nrow(obs_windows)) return(NULL)
    obs_hits <- selected_genes[tier == tier_name & class == class_name]
    obs_unique_genes <- uniqueN(obs_hits$gene_id)
    obs_genes_per_mb <- obs_unique_genes / max(sum((unique(obs_windows[, .(window_key, width_bp = 1L + 0L)])$width_bp), na.rm = TRUE), 1)
    null_vals <- replicate(null_iter, {
      sampled_keys <- sample_matched_window_keys(obs_windows, all_windows_base)
      samp_hits <- all_window_hits[window_key %chin% sampled_keys]
      c(unique_genes = uniqueN(samp_hits$gene_id), overlap_rows = nrow(samp_hits))
    })
    null_unique <- as.numeric(null_vals["unique_genes", ])
    data.table(
      tier = tier_name,
      class = class_name,
      observed_unique_genes = obs_unique_genes,
      null_mean_unique_genes = mean(null_unique, na.rm = TRUE),
      null_q95_unique_genes = as.numeric(quantile(null_unique, 0.95, na.rm = TRUE)),
      empirical_p_unique_genes = (1 + sum(null_unique >= obs_unique_genes, na.rm = TRUE)) / (1 + length(null_unique)),
      enrichment_ratio_unique_genes = obs_unique_genes / pmax(mean(null_unique, na.rm = TRUE), 1)
    )
  }), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)
fwrite(gene_overlap_null, file.path(output_dir, "gene_overlap_matched_null_by_tier.tsv"), sep = "\t")

background_genes_chr <- unique(all_window_hits[!is.na(gene_id) & nzchar(gene_id), .(chr, gene_id, eggnog_gene_name)])
go_null_summary <- rbindlist(lapply(c("gmm_rf", "continuous_distance_3class", "gmm_rf_contdist", "empirical_fdr_selected", "genomewide_fdr_selected", "global_shift_selected", "dist_margin_050"), function(tier_name) {
  rbindlist(lapply(c("balancing", "sweep"), function(class_name) {
    obs_gene_dt <- unique(selected_genes[tier == tier_name & class == class_name & !is.na(gene_id) & nzchar(gene_id), .(chr, gene_id, eggnog_gene_name)])
    obs_named <- sort(unique(obs_gene_dt[!is.na(eggnog_gene_name) & eggnog_gene_name != "-", eggnog_gene_name]))
    if (length(obs_named) < 2L) {
      return(data.table(
        tier = tier_name, class = class_name,
        observed_sig_go_terms = 0L, null_mean_sig_go_terms = NA_real_,
        null_q95_sig_go_terms = NA_real_, empirical_p_sig_go_terms = NA_real_,
        observed_mapped_genes = length(obs_named), null_iterations = 0L
      ))
    }
    chr_counts <- obs_gene_dt[!is.na(eggnog_gene_name) & eggnog_gene_name != "-", .N, by = chr]
    obs_sig <- nrow(go_results[tier == tier_name & class == class_name & subset == "whole_genome" & is.finite(p_value)])
    null_sig <- replicate(if (identical(Sys.getenv("PIPELINE_FAST_TEST"), "1")) 5L else 20L, {
      null_genes <- sample_matched_genes(chr_counts, background_genes_chr)
      res <- run_gost_safe(null_genes)
      if (is.null(res) || is.null(res$result) || !nrow(res$result)) 0L else nrow(res$result)
    })
    data.table(
      tier = tier_name,
      class = class_name,
      observed_sig_go_terms = obs_sig,
      null_mean_sig_go_terms = mean(null_sig, na.rm = TRUE),
      null_q95_sig_go_terms = as.numeric(quantile(null_sig, 0.95, na.rm = TRUE)),
      empirical_p_sig_go_terms = (1 + sum(null_sig >= obs_sig, na.rm = TRUE)) / (1 + length(null_sig)),
      observed_mapped_genes = length(obs_named),
      null_iterations = length(null_sig)
    )
  }), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)
fwrite(go_null_summary, file.path(output_dir, "GO_matched_null_summary.tsv"), sep = "\t")

immune_matched_null <- rbindlist(lapply(full_class_tier_specs$tier, function(tier_name) {
  rbindlist(lapply(c("balancing", "neutral", "sweep"), function(class_name) {
    obs_gene_dt <- unique(full_class_genes[tier == tier_name & class_short == class_name & !is.na(gene_id) & nzchar(gene_id), .(chr, gene_id, eggnog_gene_name)])
    obs_named <- sort(unique(obs_gene_dt[!is.na(eggnog_gene_name) & eggnog_gene_name != "-", eggnog_gene_name]))
    if (!length(obs_named)) return(NULL)
    obs_immune_n <- sum(obs_named %in% immune_named_genes)
    chr_counts <- obs_gene_dt[!is.na(eggnog_gene_name) & eggnog_gene_name != "-", .N, by = chr]
    null_immune <- replicate(if (identical(Sys.getenv("PIPELINE_FAST_TEST"), "1")) 100L else 500L, {
      sum(sample_matched_genes(chr_counts, background_genes_chr) %in% immune_named_genes)
    })
    top_chr <- obs_gene_dt[eggnog_gene_name %in% immune_named_genes, .N, by = chr][order(-N)][1]
    data.table(
      tier = tier_name,
      class = class_name,
      observed_immune_named_genes = obs_immune_n,
      null_mean_immune_named_genes = mean(null_immune, na.rm = TRUE),
      null_q95_immune_named_genes = as.numeric(quantile(null_immune, 0.95, na.rm = TRUE)),
      empirical_p_immune_named_genes = (1 + sum(null_immune >= obs_immune_n, na.rm = TRUE)) / (1 + length(null_immune)),
      enrichment_ratio_immune_named_genes = obs_immune_n / pmax(mean(null_immune, na.rm = TRUE), 1),
      top_immune_chr = if (nrow(top_chr)) top_chr$chr[[1]] else NA_character_
    )
  }), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)
fwrite(immune_matched_null, file.path(output_dir, "immune_matched_null_by_tier.tsv"), sep = "\t")

window_counts <- selected_windows[, .N, by = .(tier, selection_class)][order(tier, selection_class)]
gene_counts <- selected_genes[, .(
  overlapping_cds_rows = .N,
  unique_gene_ids = uniqueN(gene_id),
  unique_named_genes = uniqueN(eggnog_gene_name[!is.na(eggnog_gene_name) & eggnog_gene_name != "-"])
), by = .(tier, class)][order(tier, class)]
go_counts <- go_results[, .(n_go_terms = .N), by = .(tier, class, subset)][order(tier, class, subset)]

fwrite(window_counts, file.path(output_dir, "selected_window_counts_by_tier.tsv"), sep = "\t")
fwrite(gene_counts, file.path(output_dir, "selected_gene_counts_by_tier.tsv"), sep = "\t")
fwrite(go_counts, file.path(output_dir, "GO_term_counts_by_tier.tsv"), sep = "\t")

manifest <- data.table(
  file = list.files(output_dir, recursive = TRUE, full.names = FALSE),
  size_bytes = file.info(list.files(output_dir, recursive = TRUE, full.names = TRUE))$size
)
fwrite(manifest, file.path(output_dir, "go_gene_map_workflow_manifest.tsv"), sep = "\t")

report_file <- file.path(output_dir, "GO_GENE_MAP_WORKFLOW_REPORT.txt")
sink(report_file)
cat("Connected GO/gene-map workflow\n")
cat("Generated:", format(Sys.time()), "\n")
cat("Selection handoff:", label_rdata, "\n")
cat("GTF:", gtf_file, "\n")
cat("Output directory:", output_dir, "\n\n")
cat("Selected windows by tier/class:\n")
print(window_counts)
cat("\nSelected genes by tier/class:\n")
print(gene_counts)
cat("\nGO term counts by tier/class/subset:\n")
print(go_counts)
cat("\nImmune-related GO term counts by tier/class/subset:\n")
print(immune_go_terms[, .N, by = .(tier, class, subset)][order(tier, class, subset)])
cat("\nImmune-like gene Fisher tests:\n")
print(fisher_rows)
cat("\nImmune-like gene counts by full class and tier:\n")
print(immune_gene_sets_by_class)
cat("\nImmune-window summary by full class and tier:\n")
print(immune_window_summary)
cat("\nImmune-window class comparisons:\n")
print(immune_window_class_comparisons)
sink()

message("Connected GO/gene-map workflow complete: ", output_dir)
