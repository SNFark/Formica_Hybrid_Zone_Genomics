#!/usr/bin/env Rscript

# ============================================================
# 05_test_FST_centromere_selfalign_enrichment.R
#
# Role in pipeline:
#   Test whether highly differentiated genomic windows are enriched in
#   CentieR-predicted centromeres and/or high-density self-align regions.
#
# Important distinction from ROI analyses:
#   This is an independent genome-wide FST architecture test. It does not define
#   ROIs and it does not use the curated ROI table. The goal is to ask whether
#   regions of high differentiation are non-randomly associated with genomic
#   architecture.
#
# Main inputs:
#   - windowed or site-level Weir and Cockerham FST table
#   - CentieR/self-align centromere table
#
# Main outputs:
#   - results/FST_architecture/FST_windows_with_centromere_selfalign.tsv
#   - results/FST_architecture/high_FST_centromere_selfalign_tests.tsv
#   - results/FST_architecture/high_FST_*_summary.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# -------------------------------
# User settings: edit paths here if your files are elsewhere
# -------------------------------

fst_file <- "path/to/weir_cockerham_fst.tsv"
centromere_file <- "path/to/CentieR_and_selfalign_chromosomes.tsv"
output_dir <- "results/FST_architecture"

window_size <- 10000
high_fst_threshold <- 0.8
exclude_chromosomes <- character(0) # e.g. c("chr3") if needed

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# Helper functions
# -------------------------------

normalise_chr <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("FsiP_PB_v5_scf", "chr", x, fixed = TRUE)
  x <- ifelse(grepl("^\\d+$", x), paste0("chr", x), x)
  x <- gsub("^chr0+", "chr", x)
  x
}

read_fst <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE)
  names(x) <- gsub("WEIR_AND_COCKERHAM_FST", "fst", names(x))
  names(x) <- gsub("^CHROM$", "chrom_raw", names(x))
  names(x) <- gsub("^POS$", "pos", names(x))
  if (!all(c("chrom_raw", "pos", "fst") %in% names(x))) {
    stop("FST file must contain CHROM, POS, and WEIR_AND_COCKERHAM_FST columns.")
  }
  x %>% transmute(CHROM = normalise_chr(chrom_raw), POS = as.numeric(pos), fst = as.numeric(fst)) %>%
    filter(is.finite(POS), is.finite(fst))
}

read_centromeres <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE)
  required <- c("chr_fullname", "centieR_pred_start", "centieR_pred_stop",
                "High_density_selfalign_start", "High_density_selfalign_stop")
  if (!all(required %in% names(x))) stop("Centromere file is missing required columns.")
  x %>% transmute(
    CHROM = normalise_chr(chr_fullname),
    cent_start = as.numeric(centieR_pred_start),
    cent_end = as.numeric(centieR_pred_stop),
    self_start = as.numeric(High_density_selfalign_start),
    self_end = as.numeric(High_density_selfalign_stop)
  ) %>% filter(is.finite(cent_start), is.finite(cent_end), is.finite(self_start), is.finite(self_end))
}

run_feature_tests <- function(df, feature_col) {
  tab <- table(df$high_fst, df[[feature_col]])
  fisher <- fisher.test(tab)

  df_cmh <- df %>%
    group_by(CHROM) %>%
    filter(n_distinct(high_fst) > 1, n_distinct(.data[[feature_col]]) > 1) %>%
    ungroup()

  cmh_p <- NA_real_; cmh_or <- NA_real_
  if (nrow(df_cmh) > 0) {
    tab3d <- xtabs(as.formula(paste0("~ high_fst + ", feature_col, " + CHROM")), data = df_cmh)
    cmh <- mantelhaen.test(tab3d)
    cmh_p <- cmh$p.value
    cmh_or <- as.numeric(cmh$estimate)
  }

  summary <- df %>% group_by(.data[[feature_col]]) %>% summarise(
    n_windows = n(),
    prop_high_fst = mean(high_fst, na.rm = TRUE),
    mean_fst = mean(mean_fst, na.rm = TRUE),
    .groups = "drop"
  ) %>% rename(feature_state = 1)

  list(
    contingency_table = as.data.frame.matrix(tab),
    fisher = tibble(feature = feature_col,
                    odds_ratio = as.numeric(fisher$estimate),
                    p_value = fisher$p.value,
                    conf_low = fisher$conf.int[1],
                    conf_high = fisher$conf.int[2],
                    cmh_odds_ratio = cmh_or,
                    cmh_p_value = cmh_p),
    summary = summary
  )
}

# -------------------------------
# Load FST values, bin into fixed windows, and annotate architecture
# -------------------------------

fst <- read_fst(fst_file) %>% filter(!CHROM %in% exclude_chromosomes)
cent <- read_centromeres(centromere_file) %>% filter(!CHROM %in% exclude_chromosomes)

fst_win <- fst %>%
  mutate(win = floor(POS / window_size), start = win * window_size, end = start + window_size, mid = (start + end) / 2) %>%
  group_by(CHROM, start, end, mid) %>%
  summarise(mean_fst = mean(fst, na.rm = TRUE), n_snps = n(), .groups = "drop") %>%
  left_join(cent, by = "CHROM") %>%
  mutate(
    in_centromere = !is.na(cent_start) & end >= cent_start & start <= cent_end,
    in_selfalign = !is.na(self_start) & end >= self_start & start <= self_end,
    high_fst = mean_fst >= high_fst_threshold
  )

write_tsv(fst_win, file.path(output_dir, "FST_windows_with_centromere_selfalign.tsv"))

cent_tests <- run_feature_tests(fst_win, "in_centromere")
self_tests <- run_feature_tests(fst_win, "in_selfalign")

write_tsv(bind_rows(cent_tests$fisher, self_tests$fisher),
          file.path(output_dir, "high_FST_centromere_selfalign_tests.tsv"))
write_tsv(cent_tests$summary, file.path(output_dir, "high_FST_centromere_summary.tsv"))
write_tsv(self_tests$summary, file.path(output_dir, "high_FST_selfalign_summary.tsv"))

loo <- bind_rows(lapply(unique(fst_win$CHROM), function(chr) {
  tmp <- fst_win %>% filter(CHROM != chr)
  ft_cent <- fisher.test(table(tmp$high_fst, tmp$in_centromere))
  ft_self <- fisher.test(table(tmp$high_fst, tmp$in_selfalign))
  tibble(left_out = chr,
         feature = c("in_centromere", "in_selfalign"),
         odds_ratio = c(as.numeric(ft_cent$estimate), as.numeric(ft_self$estimate)),
         p_value = c(ft_cent$p.value, ft_self$p.value))
}))
write_tsv(loo, file.path(output_dir, "high_FST_leave_one_chromosome_out.tsv"))

chr_summary <- fst_win %>% group_by(CHROM) %>% summarise(
  n_centromere = sum(in_centromere, na.rm = TRUE),
  n_noncentromere = sum(!in_centromere, na.rm = TRUE),
  high_centromere = sum(high_fst & in_centromere, na.rm = TRUE),
  high_noncentromere = sum(high_fst & !in_centromere, na.rm = TRUE),
  prop_centromere_high = high_centromere / n_centromere,
  prop_noncentromere_high = high_noncentromere / n_noncentromere,
  enrichment_ratio = prop_centromere_high / prop_noncentromere_high,
  .groups = "drop"
)
write_tsv(chr_summary, file.path(output_dir, "high_FST_chromosome_summary.tsv"))

message("[INFO] Done.")
