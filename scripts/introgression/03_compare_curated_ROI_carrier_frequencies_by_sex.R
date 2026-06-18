#!/usr/bin/env Rscript

# ============================================================
# 03_compare_curated_ROI_carrier_frequencies_by_sex.R
#
# Role in pipeline:
#   Compare female/worker and male carrier proportions within the FINAL,
#   manually curated ROI set.
#
# Critical input distinction:
#   This script must use the curated ROI table, not the automatic candidate
#   output from Script 01. Candidate ROIs were only a screening aid. The final
#   ROI set was curated after inspection of smoothed female introgression panels
#   so that biologically meaningful signals, including long low-frequency
#   contiguous tracts, could be retained.
#
# Why carrier proportion:
#   Females/workers are diploid and males are haploid. To compare sexes fairly,
#   this script classifies each individual as carrier/non-carrier within each
#   curated ROI and then compares the proportion of carriers between sexes.
#
# DIEM coding assumed here:
#   Females/workers: 0 = F. selysi homozygous, 1 = heterozygous,
#                    2 = F. cinerea homozygous
#   Males:           0 = F. selysi allele, 1/2 = non-selysi/cinerea-like state
#
# Main inputs:
#   - results/ROIs/curated_ROIs_female_defined.tsv
#   - filtered_geno_80_sex1.RData              object: genotypes, males
#   - filtered_geno_80_sex2.RData              object: genotypes, females/workers
#   - hybrid_idx_filtered_80_sex1.RData        object: h, males
#   - hybrid_idx_filtered_80_sex2.RData        object: h, females/workers
#   - markers.pos
#   - full_res.RData                           object: res$DI$DI
#   - samples_metadata.csv
#
# Main outputs:
#   - results/ROI_sex_comparison/*_all_ROIs.tsv
#   - results/ROI_sex_comparison/*_plot_ROIs.tsv/csv
#   - results/ROI_sex_comparison/*.png/pdf
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(readxl)
  library(ggplot2)
})

# -------------------------------
# User settings: edit paths here if your files are elsewhere
# -------------------------------

male_geno_rdata   <- "filtered_geno_80_sex1.RData"        # object: genotypes
female_geno_rdata <- "filtered_geno_80_sex2.RData"        # object: genotypes
male_hi_rdata     <- "hybrid_idx_filtered_80_sex1.RData"  # object: h
female_hi_rdata   <- "hybrid_idx_filtered_80_sex2.RData"  # object: h
metadata_file     <- "samples_metadata.csv"
markers_file      <- "markers.pos"
full_res_rdata    <- "full_res.RData"                     # object: res$DI$DI
roi_file          <- "results/ROIs/curated_ROIs_female_defined.tsv"

output_dir <- "results/ROI_sex_comparison"
output_prefix <- "ROI_worker_vs_male_carrier_proportions"

di_percentile <- 0.80       # top 20% diagnostic markers; used only to align marker coordinates
hi_threshold <- 0.75        # retain F. cinerea-background individuals

# Per-ROI carrier-classification thresholds.
coverage_min <- 0.60        # individual must have calls for >=60% of markers in the curated ROI
female_het_threshold <- 0.50 # female heterozygous carrier: >=50% genotype-1 sites
female_hom_threshold <- 0.75 # female homozygous introgressed carrier: >=75% genotype-0 sites
male_any_threshold <- 0.50   # male carrier: >=50% genotype-0 sites

# Manuscript-style plotting/export filters. These affect the plotted subset only.
exclude_chromosomes <- c("chr8")
exclude_artifact_chr <- "chr15"
exclude_artifact_start_mb <- 10.2
exclude_artifact_end_mb <- 10.5
exclude_artifact_tolerance_mb <- 0.02

# Significance labels.
star_p1 <- 0.051
star_p2 <- 0.01
star_p3 <- 0.001

# Plot output.
plot_width <- 14
plot_height <- 6
plot_dpi <- 300

# Set FALSE to skip manuscript-style plot filtering; all ROI results are always exported.
apply_final_plot_filters <- TRUE

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# Helper functions
# -------------------------------

normalise_chr <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("FsiP_PB_v5_scf", "chr", x, fixed = TRUE)
  x <- gsub("^scf", "chr", x)
  x <- ifelse(grepl("^\\d+$", x), paste0("chr", x), x)
  x <- gsub("^chr0+", "chr", x)
  x
}

chr_number <- function(x) suppressWarnings(as.integer(sub("^chr", "", as.character(x))))

read_diem_genotypes <- function(path) {
  load(path) # loads object named genotypes
  g <- as.matrix(genotypes)
  rm(genotypes)
  g[g == "_"] <- NA
  matrix(as.numeric(g), nrow = nrow(g), ncol = ncol(g))
}

read_metadata <- function(path) {
  # base R is robust to the extra empty columns in the metadata file.
  meta <- read.delim(path, sep = ",")
  if (!"genotype" %in% names(meta)) stop("Metadata file must contain a 'genotype' column.")
  meta$sex <- nchar(meta$genotype)
  meta
}

read_rois <- function(path) {
  ext <- tolower(tools::file_ext(path))
  x <- if (ext %in% c("xlsx", "xls")) {
    readxl::read_xlsx(path)
  } else {
    readr::read_tsv(path, show_col_types = FALSE)
  }

  nms <- names(x)
  chr_col <- nms[grepl("^(chr|chrom|chromosome|scaffold)$", nms, ignore.case = TRUE)][1]
  start_col <- nms[grepl("^(start_bp|start|mb_start|start_mb)$", nms, ignore.case = TRUE)][1]
  end_col <- nms[grepl("^(end_bp|end|mb_end|end_mb|stop)$", nms, ignore.case = TRUE)][1]

  if (any(is.na(c(chr_col, start_col, end_col)))) {
    stop("Could not infer ROI coordinate columns. Expected columns like chrom/chr, mb_start/start_bp, mb_end/end_bp.")
  }

  out <- x %>%
    transmute(
      chr = normalise_chr(.data[[chr_col]]),
      start_raw = as.numeric(.data[[start_col]]),
      end_raw = as.numeric(.data[[end_col]]),
      class = if ("class" %in% names(x)) as.character(.data[["class"]]) else NA_character_
    ) %>%
    mutate(
      # Treat values <1000 as Mb; larger values as bp.
      start_bp = ifelse(start_raw < 1000, start_raw * 1e6, start_raw),
      end_bp = ifelse(end_raw < 1000, end_raw * 1e6, end_raw)
    ) %>%
    select(chr, start_bp, end_bp, class) %>%
    filter(is.finite(start_bp), is.finite(end_bp), end_bp > start_bp)

  out
}

pick_test <- function(k_f_total, n_f, k_m_any, n_m) {
  if (any(is.na(c(k_f_total, n_f, k_m_any, n_m))) || n_f == 0 || n_m == 0) {
    return(list(p = NA_real_, method = NA_character_))
  }

  mat <- matrix(
    c(k_f_total, n_f - k_f_total,
      k_m_any,   n_m - k_m_any),
    nrow = 2,
    byrow = TRUE
  )

  if (min(mat) < 5) {
    list(p = fisher.test(mat)$p.value, method = "Fisher")
  } else {
    list(p = prop.test(c(k_f_total, k_m_any), c(n_f, n_m), correct = TRUE)$p.value,
         method = "prop.test")
  }
}

# -------------------------------
# Load curated ROIs and DIEM inputs
# -------------------------------

fm <- read_diem_genotypes(male_geno_rdata)
ff <- read_diem_genotypes(female_geno_rdata)

load(full_res_rdata) # object: res$DI$DI
if (!exists("res") || is.null(res$DI$DI)) stop("full_res_rdata must contain res$DI$DI.")

markers <- read.delim(markers_file, col.names = c("chr", "pos"))
keep_mark <- res$DI$DI > quantile(res$DI$DI, prob = di_percentile)
markers <- markers[keep_mark, , drop = FALSE]

# The genotype matrices are already DI80-filtered; do not subset fm/ff again.
stopifnot(ncol(fm) == nrow(markers), ncol(ff) == nrow(markers))

markers$chr <- normalise_chr(markers$chr)
markers$pos <- as.numeric(markers$pos)
markers$chr <- factor(markers$chr, levels = paste0("chr", 1:27))

ord <- order(markers$chr, markers$pos)
markers <- markers[ord, , drop = FALSE]
fm <- fm[, ord, drop = FALSE]
ff <- ff[, ord, drop = FALSE]

meta <- read_metadata(metadata_file)

load(female_hi_rdata); hf <- as.numeric(h); rm(h)
load(male_hi_rdata); hm <- as.numeric(h); rm(h)

idx_fem <- which(meta$sex == 2)
idx_mal <- which(meta$sex == 1)
if (length(hf) == nrow(meta)) hf <- hf[idx_fem]
if (length(hm) == nrow(meta)) hm <- hm[idx_mal]

stopifnot(length(hf) == nrow(ff), length(hm) == nrow(fm))

keep_fem <- which(hf >= hi_threshold & !is.na(hf))
keep_male <- which(hm >= hi_threshold & !is.na(hm))

message("[INFO] Cinerea-background workers/females: ", length(keep_fem), " / ", nrow(ff))
message("[INFO] Cinerea-background males:           ", length(keep_male), " / ", nrow(fm))

rois <- read_rois(roi_file) %>%
  filter(chr %in% unique(as.character(markers$chr))) %>%
  mutate(chr = factor(chr, levels = paste0("chr", 1:27))) %>%
  arrange(chr, start_bp)

# -------------------------------
# ROI carrier classification
# -------------------------------

eval_roi <- function(chr, start_bp, end_bp, class) {
  chr <- as.character(chr)
  idx <- which(as.character(markers$chr) == chr & markers$pos >= start_bp & markers$pos <= end_bp)

  if (!length(idx)) {
    return(tibble(
      chr = chr, start_mb = start_bp / 1e6, end_mb = end_bp / 1e6, class = class,
      n_markers = 0, n_f = 0, n_m = 0,
      k_f_het = 0, k_f_hom = 0, k_f_total = 0, k_m_any = 0,
      p_f_het = NA_real_, p_f_hom = NA_real_, p_f_total = NA_real_, p_m_any = NA_real_,
      diff_total = NA_real_, p_female_total_vs_male = NA_real_, test = NA_character_
    ))
  }

  gF <- ff[keep_fem, idx, drop = FALSE]
  gM <- fm[keep_male, idx, drop = FALSE]

  covF <- rowMeans(!is.na(gF))
  covM <- rowMeans(!is.na(gM))
  useF <- covF >= coverage_min
  useM <- covM >= coverage_min

  share0_F <- rowMeans(gF == 0, na.rm = TRUE)
  share1_F <- rowMeans(gF == 1, na.rm = TRUE)
  share0_M <- rowMeans(gM == 0, na.rm = TRUE)

  car_f_hom <- (share0_F >= female_hom_threshold) & useF
  car_f_het <- (share1_F >= female_het_threshold) & useF
  car_m_any <- (share0_M >= male_any_threshold) & useM

  n_f <- sum(useF)
  n_m <- sum(useM)

  k_f_hom <- sum(car_f_hom, na.rm = TRUE)
  k_f_het <- sum(car_f_het, na.rm = TRUE)
  k_f_total <- k_f_hom + k_f_het
  k_m_any <- sum(car_m_any, na.rm = TRUE)

  p_f_hom <- if (n_f > 0) k_f_hom / n_f else NA_real_
  p_f_het <- if (n_f > 0) k_f_het / n_f else NA_real_
  p_f_total <- if (n_f > 0) k_f_total / n_f else NA_real_
  p_m_any <- if (n_m > 0) k_m_any / n_m else NA_real_

  tst <- pick_test(k_f_total, n_f, k_m_any, n_m)

  tibble(
    chr = chr,
    start_mb = round(start_bp / 1e6, 2),
    end_mb = round(end_bp / 1e6, 2),
    class = class,
    n_markers = length(idx),
    n_f = n_f, n_m = n_m,
    k_f_het = k_f_het, k_f_hom = k_f_hom, k_f_total = k_f_total, k_m_any = k_m_any,
    p_f_het = p_f_het, p_f_hom = p_f_hom, p_f_total = p_f_total, p_m_any = p_m_any,
    diff_total = p_f_total - p_m_any,
    p_female_total_vs_male = tst$p,
    test = tst$method
  )
}

roi_results_all <- rois %>%
  mutate(row_id = row_number()) %>%
  pmap_dfr(function(chr, start_bp, end_bp, class, row_id) {
    eval_roi(chr, start_bp, end_bp, class)
  }) %>%
  mutate(
    p_adj_BH = p.adjust(p_female_total_vs_male, method = "BH"),
    star = case_when(
      is.na(p_female_total_vs_male) ~ "",
      p_female_total_vs_male <= star_p3 ~ "***",
      p_female_total_vs_male <= star_p2 ~ "**",
      p_female_total_vs_male <= star_p1 ~ "*",
      TRUE ~ ""
    ),
    chr = factor(chr, levels = paste0("chr", 1:27))
  ) %>%
  arrange(chr, start_mb)

roi_results_plot <- roi_results_all
if (apply_final_plot_filters) {
  roi_results_plot <- roi_results_plot %>%
    filter(!as.character(chr) %in% exclude_chromosomes) %>%
    filter(!(as.character(chr) == exclude_artifact_chr &
             abs(start_mb - exclude_artifact_start_mb) <= exclude_artifact_tolerance_mb &
             abs(end_mb - exclude_artifact_end_mb) <= exclude_artifact_tolerance_mb))
}

roi_results_plot <- roi_results_plot %>%
  mutate(
    roi_label = paste0(as.character(chr), ":",
                       formatC(start_mb, format = "f", digits = 1), "–",
                       formatC(end_mb, format = "f", digits = 1))
  ) %>%
  arrange(chr, start_mb)
roi_results_plot$roi_label <- factor(roi_results_plot$roi_label, levels = roi_results_plot$roi_label)

# -------------------------------
# Export results
# -------------------------------

write_tsv(roi_results_all, file.path(output_dir, paste0(output_prefix, "_all_ROIs.tsv")))
write_tsv(roi_results_plot, file.path(output_dir, paste0(output_prefix, "_plot_ROIs.tsv")))
write_csv(roi_results_plot, file.path(output_dir, paste0(output_prefix, "_plot_ROIs.csv")))

message("[INFO] Wrote ROI sex-comparison tables to: ", output_dir)

# -------------------------------
# Plot
# -------------------------------

p <- ggplot(roi_results_plot, aes(x = roi_label, y = p_f_total)) +
  geom_point(size = 3.2, color = "black") +
  geom_point(aes(y = p_m_any), shape = 1, size = 3.2, color = "navy") +
  geom_segment(aes(xend = roi_label, y = p_m_any, yend = p_f_total),
               alpha = 0.35, color = "gray40") +
  geom_text(aes(label = star), vjust = -0.9, size = 5.2, fontface = "bold", color = "black") +
  labs(x = "Chromosome region (Mb)",
       y = "Carrier proportion",
       title = "Introgressed regions of interest in workers versus males",
       subtitle = "Filled = workers/females total carriers; open = male carriers") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        panel.grid.minor = element_blank())

ggsave(file.path(output_dir, paste0(output_prefix, ".png")), p, width = plot_width, height = plot_height, dpi = plot_dpi)
ggsave(file.path(output_dir, paste0(output_prefix, ".pdf")), p, width = plot_width, height = plot_height)

message("[INFO] Done.")
