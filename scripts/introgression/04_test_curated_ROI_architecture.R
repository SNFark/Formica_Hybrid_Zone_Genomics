#!/usr/bin/env Rscript

# ============================================================
# 04_test_curated_ROI_architecture.R
#
# Role in pipeline:
#   Test whether the manually curated female-defined ROIs are associated with
#   genome architecture, especially centromeric/self-align regions and high-
#   confidence inversions.
#
# Critical input distinction:
#   This script uses the final curated ROI table. It does not use the automatic
#   candidate ROI output from Script 01.
#
# Analyses performed:
#   A. Classify curated ROIs as centromeric/self-align when at least 50% of
#      their length overlaps the selected feature track.
#   B. Test whether centromeric/self-align ROIs are longer, and whether long
#      ROIs are enriched in centromeric/self-align regions.
#   C. Use chromosome-wise circular permutations to preserve chromosome identity
#      and ROI length while randomising genomic position.
#   D. Compare windowed carrier introgression inside centromeric/self-align
#      regions with matched adjacent chromosome-arm segments.
#   E. Annotate curated ROIs for overlap with high-confidence inversions.
#
# Main inputs:
#   - results/ROIs/curated_ROIs_female_defined.tsv
#   - results/introgression_windows/windowed_introgression_carriers_female_DI80_HI075_100kb.tsv
#   - CentieR/self-align centromere table
#   - high-confidence inversion table
#
# Main outputs:
#   - results/ROI_architecture/*.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
})

# -------------------------------
# User settings: edit paths here if your files are elsewhere
# -------------------------------

roi_file <- "results/ROIs/curated_ROIs_female_defined.tsv"
windowed_introgression_file <- "results/introgression_windows/windowed_introgression_carriers_female_DI80_HI075_100kb.tsv"
centromere_file <- "path/to/CentieR_and_selfalign_chromosomes.tsv"
inversion_file <- "path/to/high_confidence_inversions.tsv"

output_dir <- "results/ROI_architecture"
centromere_track <- "selfalign"  # "selfalign" or "centieR"
exclude_chromosomes <- c("chr3")

centromere_overlap_threshold <- 0.50
long_threshold_bp <- 250000
introgression_threshold <- 0.03
n_perm <- 2400

set.seed(1)
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

interval_overlap_bp <- function(start1, end1, start2, end2) {
  pmax(0, pmin(end1, end2) - pmax(start1, start2))
}

read_centromeres <- function(path, track = c("selfalign", "centieR")) {
  track <- match.arg(track)
  x <- read_tsv(path, show_col_types = FALSE)
  start_col <- if (track == "selfalign") "High_density_selfalign_start" else "centieR_pred_start"
  end_col <- if (track == "selfalign") "High_density_selfalign_stop" else "centieR_pred_stop"

  x %>%
    transmute(
      chr = normalise_chr(chr_fullname),
      start_bp = as.numeric(.data[[start_col]]),
      end_bp = as.numeric(.data[[end_col]])
    ) %>%
    filter(is.finite(start_bp), is.finite(end_bp), end_bp > start_bp) %>%
    group_by(chr) %>%
    summarise(start_bp = min(start_bp), end_bp = max(end_bp), .groups = "drop")
}

read_inversions <- function(path) {
  inv <- read_tsv(path, show_col_types = FALSE)
  nms <- names(inv)
  chr_col <- nms[grepl("^(chr|chrom|chromosome|scaffold|CHROM|Scaffold)$", nms)][1]
  start_col <- nms[grepl("start|START|POS|Start", nms)][1]
  end_col <- nms[grepl("end|END|stop|STOP|End", nms)][1]
  if (any(is.na(c(chr_col, start_col, end_col)))) stop("Could not infer inversion columns.")

  inv %>%
    transmute(
      chr = normalise_chr(.data[[chr_col]]),
      start_bp = as.numeric(.data[[start_col]]),
      end_bp = as.numeric(.data[[end_col]]),
      confidence = if ("confidence" %in% nms) as.character(.data[["confidence"]]) else NA_character_
    ) %>%
    mutate(
      start_bp = ifelse(start_bp < 1000, start_bp * 1e6, start_bp),
      end_bp = ifelse(end_bp < 1000, end_bp * 1e6, end_bp)
    ) %>%
    filter(is.finite(start_bp), is.finite(end_bp), end_bp > start_bp)
}

annotate_overlap <- function(rois, features, feature_name, frac_threshold = 0.50) {
  out <- rois %>%
    mutate(feature = feature_name,
           feature_overlap_bp = 0,
           feature_fraction = 0,
           overlaps_feature = FALSE)

  for (i in seq_len(nrow(out))) {
    f <- features %>% filter(chr == out$chr[i])
    if (!nrow(f)) next
    ov <- interval_overlap_bp(out$start_bp[i], out$end_bp[i], f$start_bp, f$end_bp)
    total_ov <- sum(ov, na.rm = TRUE)
    len <- out$end_bp[i] - out$start_bp[i]
    out$feature_overlap_bp[i] <- total_ov
    out$feature_fraction[i] <- ifelse(len > 0, total_ov / len, 0)
    out$overlaps_feature[i] <- out$feature_fraction[i] >= frac_threshold
  }
  out
}

fisher_haldane <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat)
  mat_ha <- mat + 0.5
  or_ha <- (mat_ha[1,1] * mat_ha[2,2]) / (mat_ha[1,2] * mat_ha[2,1])
  tibble(odds_ratio_fisher = as.numeric(ft$estimate),
         odds_ratio_HA = as.numeric(or_ha),
         p_value = ft$p.value,
         conf_low = ft$conf.int[1],
         conf_high = ft$conf.int[2])
}

circular_shift_rois <- function(rois, chr_lengths) {
  rois %>%
    group_by(chr) %>%
    group_modify(function(df, key) {
      L <- chr_lengths$chr_len[match(key$chr, chr_lengths$chr)]
      if (!is.finite(L) || L <= 0) return(df)
      shift <- sample.int(floor(L), 1)
      len <- df$end_bp - df$start_bp
      new_start <- ((df$start_bp + shift - 1) %% L) + 1
      new_end <- pmin(new_start + len, L)
      df$start_bp <- new_start
      df$end_bp <- new_end
      df$length_bp <- df$end_bp - df$start_bp
      df
    }) %>% ungroup()
}

# -------------------------------
# Load data
# -------------------------------

rois <- read_tsv(roi_file, show_col_types = FALSE) %>%
  mutate(chr = normalise_chr(chr)) %>%
  filter(!chr %in% exclude_chromosomes) %>%
  rename_with(~"start_bp", matches("^start$|^start_bp$")) %>%
  rename_with(~"end_bp", matches("^end$|^end_bp$")) %>%
  mutate(start_bp = as.numeric(start_bp), end_bp = as.numeric(end_bp), length_bp = end_bp - start_bp) %>%
  filter(is.finite(start_bp), is.finite(end_bp), end_bp > start_bp)

win <- read_tsv(windowed_introgression_file, show_col_types = FALSE) %>%
  mutate(chr = normalise_chr(chr),
         start = as.numeric(start), end = as.numeric(end),
         freq_use = if ("freq_sm" %in% names(.)) freq_sm else freq)

centromeres <- read_centromeres(centromere_file, track = centromere_track) %>%
  filter(!chr %in% exclude_chromosomes)
inversions <- read_inversions(inversion_file) %>% filter(!chr %in% exclude_chromosomes)

chr_lengths <- win %>% group_by(chr) %>% summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop")

# -------------------------------
# A. Curated ROI persistence in centromeric/self-align regions
# -------------------------------

rois_cent <- annotate_overlap(rois, centromeres, centromere_track, centromere_overlap_threshold) %>%
  mutate(is_long = length_bp >= long_threshold_bp)

write_tsv(rois_cent, file.path(output_dir, "ROIs_with_centromere_overlap.tsv"))

length_test <- wilcox.test(length_bp ~ overlaps_feature, data = rois_cent, exact = FALSE)
length_summary <- rois_cent %>%
  group_by(overlaps_feature) %>%
  summarise(n = n(), median_length_bp = median(length_bp), mean_length_bp = mean(length_bp), .groups = "drop")

delta_median <- with(length_summary,
  median_length_bp[overlaps_feature == TRUE] - median_length_bp[overlaps_feature == FALSE])

# Long ROI enrichment.
a <- sum(rois_cent$is_long & rois_cent$overlaps_feature)
b <- sum(rois_cent$is_long & !rois_cent$overlaps_feature)
c <- sum(!rois_cent$is_long & rois_cent$overlaps_feature)
d <- sum(!rois_cent$is_long & !rois_cent$overlaps_feature)

long_enrichment <- fisher_haldane(a, b, c, d) %>%
  mutate(long_threshold_bp = long_threshold_bp,
         n_long_centromere = a, n_long_arm = b,
         n_short_centromere = c, n_short_arm = d)

# Circular permutations.
message("[INFO] Running circular ROI permutations: ", n_perm)
obs_delta <- delta_median
obs_or <- long_enrichment$odds_ratio_HA[1]

perm_stats <- map_dfr(seq_len(n_perm), function(i) {
  rp <- circular_shift_rois(rois, chr_lengths)
  rp_cent <- annotate_overlap(rp, centromeres, centromere_track, centromere_overlap_threshold) %>%
    mutate(is_long = length_bp >= long_threshold_bp)

  med_tab <- rp_cent %>% group_by(overlaps_feature) %>% summarise(med = median(length_bp), .groups = "drop")
  delta <- with(med_tab, med[overlaps_feature == TRUE] - med[overlaps_feature == FALSE])
  if (length(delta) == 0) delta <- NA_real_

  a <- sum(rp_cent$is_long & rp_cent$overlaps_feature)
  b <- sum(rp_cent$is_long & !rp_cent$overlaps_feature)
  c <- sum(!rp_cent$is_long & rp_cent$overlaps_feature)
  d <- sum(!rp_cent$is_long & !rp_cent$overlaps_feature)
  or_ha <- fisher_haldane(a, b, c, d)$odds_ratio_HA[1]
  tibble(perm = i, delta_median_bp = delta, odds_ratio_HA = or_ha)
})

perm_p <- tibble(
  statistic = c("delta_median_bp", "long_ROI_odds_ratio_HA"),
  observed = c(obs_delta, obs_or),
  empirical_p_high = c(
    mean(perm_stats$delta_median_bp >= obs_delta, na.rm = TRUE),
    mean(perm_stats$odds_ratio_HA >= obs_or, na.rm = TRUE)
  )
)

write_tsv(length_summary, file.path(output_dir, "ROI_length_summary_by_centromere_overlap.tsv"))
write_tsv(tibble(test = "ROI length: centromere/self-align vs arms", delta_median_bp = delta_median, wilcox_p = length_test$p.value),
          file.path(output_dir, "ROI_length_wilcoxon.tsv"))
write_tsv(long_enrichment, file.path(output_dir, "long_ROI_centromere_enrichment.tsv"))
write_tsv(perm_stats, file.path(output_dir, "ROI_centromere_permutation_stats.tsv"))
write_tsv(perm_p, file.path(output_dir, "ROI_centromere_permutation_pvalues.tsv"))

# -------------------------------
# B. Centromere-associated depletion/enrichment of windowed carrier introgression
# -------------------------------

# Mark windows that overlap the selected centromere/self-align interval.
win_cent <- win %>%
  left_join(centromeres, by = "chr", suffix = c("", "_cent")) %>%
  mutate(
    overlap_bp = pmax(0, pmin(end, end_bp) - pmax(start, start_bp)),
    in_centromere = !is.na(overlap_bp) & overlap_bp > 0,
    introgressed_window = freq_use >= introgression_threshold
  )

# Build equal-length adjacent arm segments, choosing the side with more available sequence.
arm_segments <- centromeres %>%
  left_join(chr_lengths, by = "chr") %>%
  rowwise() %>%
  mutate(
    cen_len = end_bp - start_bp,
    left_len = max(0, start_bp - 1),
    right_len = max(0, chr_len - end_bp),
    use_right = right_len >= left_len,
    arm_start = ifelse(use_right, end_bp + 1, pmax(1, start_bp - cen_len)),
    arm_end = ifelse(use_right, pmin(chr_len, end_bp + cen_len), start_bp - 1)
  ) %>%
  ungroup() %>%
  transmute(chr, arm_start, arm_end) %>%
  filter(is.finite(arm_start), is.finite(arm_end), arm_end > arm_start)

win_arm <- win %>%
  left_join(arm_segments, by = "chr") %>%
  mutate(
    arm_overlap_bp = pmax(0, pmin(end, arm_end) - pmax(start, arm_start)),
    in_arm_match = !is.na(arm_overlap_bp) & arm_overlap_bp > 0,
    introgressed_window = freq_use >= introgression_threshold
  )

cent_chr <- win_cent %>% filter(in_centromere) %>% group_by(chr) %>%
  summarise(n_windows_cent = n(),
            intro_windows_cent = sum(introgressed_window, na.rm = TRUE),
            frac_intro_cent = mean(introgressed_window, na.rm = TRUE),
            mean_freq_cent = mean(freq_use, na.rm = TRUE), .groups = "drop")

arm_chr <- win_arm %>% filter(in_arm_match) %>% group_by(chr) %>%
  summarise(n_windows_arm = n(),
            intro_windows_arm = sum(introgressed_window, na.rm = TRUE),
            frac_intro_arm = mean(introgressed_window, na.rm = TRUE),
            mean_freq_arm = mean(freq_use, na.rm = TRUE), .groups = "drop")

cent_depletion <- centromeres %>%
  mutate(centromere_length_bp = end_bp - start_bp) %>%
  left_join(cent_chr, by = "chr") %>%
  left_join(arm_chr, by = "chr") %>%
  mutate(
    observed_minus_expected_intro_windows = intro_windows_cent - n_windows_cent * frac_intro_arm,
    zero_introgressed_centromere = intro_windows_cent == 0,
    binom_p_zero = ifelse(n_windows_cent > 0 & is.finite(frac_intro_arm),
                          dbinom(0, size = n_windows_cent, prob = frac_intro_arm), NA_real_),
    binom_p_zero_BH = p.adjust(binom_p_zero, method = "BH")
  )

paired_frac <- wilcox.test(cent_depletion$frac_intro_cent, cent_depletion$frac_intro_arm, paired = TRUE, exact = FALSE)
paired_mean <- wilcox.test(cent_depletion$mean_freq_cent, cent_depletion$mean_freq_arm, paired = TRUE, exact = FALSE)
size_lm <- lm(observed_minus_expected_intro_windows ~ centromere_length_bp, data = cent_depletion)

cent_depletion_tests <- tibble(
  test = c("paired_fraction_introgressed_windows", "paired_mean_introgression_frequency", "centromere_size_effect"),
  statistic = c(as.numeric(paired_frac$statistic), as.numeric(paired_mean$statistic), summary(size_lm)$r.squared),
  p_value = c(paired_frac$p.value, paired_mean$p.value, summary(size_lm)$coefficients[2,4])
)

write_tsv(cent_depletion, file.path(output_dir, "centromere_introgression_depletion_enrichment.tsv"))
write_tsv(cent_depletion_tests, file.path(output_dir, "centromere_depletion_tests.tsv"))

# -------------------------------
# C. Curated ROI overlap with high-confidence inversions
# -------------------------------

rois_inv <- annotate_overlap(rois, inversions, "inversion", frac_threshold = 0.001)
write_tsv(rois_inv, file.path(output_dir, "ROIs_with_inversion_overlap.tsv"))

message("[INFO] Done.")
