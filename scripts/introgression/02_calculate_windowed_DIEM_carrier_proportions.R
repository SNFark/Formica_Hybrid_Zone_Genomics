#!/usr/bin/env Rscript

# ============================================================
# 02_calculate_windowed_DIEM_carrier_proportions.R
#
# Role in pipeline:
#   Calculate genome-wide windowed carrier proportions from DIEM genotype
#   matrices, separately for females/workers and males.
#
# Why carrier proportion instead of allele frequency:
#   Females/workers are diploid whereas males are haploid. Allele-frequency
#   summaries are therefore not directly comparable across sexes. This script
#   classifies each individual as a carrier/non-carrier in each genomic window
#   and then calculates the proportion of carriers in that window.
#
# Interpretation:
#   These windowed tracks describe the genome-wide introgression landscape.
#   They are also used later for architecture plots and centromere-window
#   enrichment/depletion tests. 
#
# Main inputs:
#   - filtered_geno_80_sex1.RData              object: genotypes, males
#   - filtered_geno_80_sex2.RData              object: genotypes, females/workers
#   - hybrid_idx_filtered_80_sex1.RData        object: h, males
#   - hybrid_idx_filtered_80_sex2.RData        object: h, females/workers
#   - markers.pos
#   - full_res.RData                           object: res$DI$DI
#   - samples_metadata.csv
#
# Main outputs:
#   - results/introgression_windows/windowed_introgression_carriers_female_DI80_HI075_100kb.tsv
#   - results/introgression_windows/windowed_introgression_carriers_male_DI80_HI075_100kb.tsv
#   - results/introgression_windows/chromosome_mean_carrier_frequency_female_vs_male.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(zoo)
})

# -------------------------------
# User settings: edit paths here if your files are elsewhere
# -------------------------------

male_geno_rdata   <- "filtered_geno_80_sex1.RData"   # object: genotypes
female_geno_rdata <- "filtered_geno_80_sex2.RData"   # object: genotypes
male_hi_rdata     <- "hybrid_idx_filtered_80_sex1.RData" # object: h
female_hi_rdata   <- "hybrid_idx_filtered_80_sex2.RData" # object: h
metadata_file     <- "samples_metadata.csv"
markers_file      <- "markers.pos"
full_res_rdata    <- "full_res.RData" # object: res$DI$DI, used ONLY to filter marker table

output_dir <- "results/introgression_windows"

di_percentile <- 0.80       # retain the top 20% most diagnostic markers
hi_threshold <- 0.75        # retain F. cinerea-background individuals
window_size <- 100000       # 100 kb
smooth_kb <- 100            # set to 0 to skip smoothing
intro_state <- 0            # DIEM state corresponding to F. selysi ancestry
carrier_threshold <- 0.50   # carrier if >=50% informative markers carry F. selysi ancestry
cap <- 0.30                 # retained for plotting consistency; not used in tests

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

read_diem_genotypes <- function(path) {
  load(path) # loads object named genotypes
  g <- as.matrix(genotypes)
  rm(genotypes)
  g[g == "_"] <- NA
  matrix(as.numeric(g), nrow = nrow(g), ncol = ncol(g))
}

smooth_fun <- function(df, smooth_kb) {
  if (smooth_kb <= 0) {
    df$freq_sm <- df$freq
    return(df)
  }
  df %>%
    group_by(chr) %>%
    arrange(start, .by_group = TRUE) %>%
    mutate(
      wsize = median(end - start, na.rm = TRUE),
      k = max(1L, round((smooth_kb * 1000) / wsize)),
      freq_sm = zoo::rollapply(
        freq, width = k,
        FUN = function(x) mean(x, na.rm = TRUE),
        align = "center", fill = NA_real_
      )
    ) %>%
    ungroup()
}

make_carrier_windows <- function(geno, keep_idx, markers,
                                 window_size = 100000,
                                 intro_state = 0,
                                 carrier_threshold = 0.5) {
  chroms <- levels(markers$chr)

  map_dfr(chroms, function(chr) {
    sel <- markers$chr == chr
    pos <- markers$pos[sel]
    cols <- which(sel)
    if (!length(cols)) return(NULL)

    maxp <- max(pos, na.rm = TRUE)
    breaks <- seq(0, maxp + window_size, by = window_size)
    w_id <- cut(pos, breaks = breaks, include.lowest = TRUE,
                right = FALSE, labels = FALSE)
    uwin <- sort(unique(w_id))

    g_chr <- geno[keep_idx, cols, drop = FALSE]

    map_dfr(uwin, function(w) {
      cidx <- which(w_id == w)
      if (!length(cidx)) return(NULL)

      g_sub <- g_chr[, cidx, drop = FALSE]
      valid_rows <- rowSums(!is.na(g_sub)) > 0
      n_ind <- sum(valid_rows)

      if (n_ind == 0) {
        carrier_prop <- NA_real_
      } else {
        g_use <- g_sub[valid_rows, , drop = FALSE]
        frac_intro <- rowMeans(g_use == intro_state, na.rm = TRUE)
        carriers <- frac_intro >= carrier_threshold
        carrier_prop <- mean(carriers, na.rm = TRUE)
      }

      tibble(
        chr = as.character(chr),
        win = w,
        start = breaks[w],
        end = min(breaks[w] + window_size, maxp),
        n_markers = length(cidx),
        n_valid = n_ind,
        freq = carrier_prop,
        mid = (breaks[w] + min(breaks[w] + window_size, maxp)) / 2
      )
    })
  })
}

# -------------------------------
# Load DIEM genotypes, hybrid indices, metadata, and marker map
# -------------------------------

fm <- read_diem_genotypes(male_geno_rdata)
ff <- read_diem_genotypes(female_geno_rdata)

load(full_res_rdata) # res$DI$DI
markers <- read.delim(markers_file, col.names = c("chr", "pos"))
keep_mark <- res$DI$DI > quantile(res$DI$DI, prob = di_percentile)
markers <- markers[keep_mark, , drop = FALSE]

# IMPORTANT: filtered_geno_80_* matrices are already DI80; do not subset fm/ff again.
stopifnot(ncol(fm) == nrow(markers), ncol(ff) == nrow(markers))

markers$chr <- normalise_chr(markers$chr)
markers$pos <- as.numeric(markers$pos)
markers$chr <- factor(markers$chr, levels = unique(markers$chr))

meta <- read.delim(metadata_file, sep = ",")
meta$sex <- nchar(meta$genotype)

load(male_hi_rdata); hm <- as.numeric(h); rm(h)
idx_mal <- which(meta$sex == 1)
if (length(hm) == nrow(meta)) hm <- hm[idx_mal]
stopifnot(length(hm) == nrow(fm))

load(female_hi_rdata); hf <- as.numeric(h); rm(h)
idx_fem <- which(meta$sex == 2)
if (length(hf) == nrow(meta)) hf <- hf[idx_fem]
stopifnot(length(hf) == nrow(ff))

keep_male <- which(hm >= hi_threshold & !is.na(hm))
keep_fem <- which(hf >= hi_threshold & !is.na(hf))

message("[INFO] Cinerea-background males:   ", length(keep_male), " / ", nrow(fm))
message("[INFO] Cinerea-background females: ", length(keep_fem), " / ", nrow(ff))

# Order markers and genotype columns consistently.
ord <- order(markers$chr, markers$pos)
markers <- markers[ord, , drop = FALSE]
fm <- fm[, ord, drop = FALSE]
ff <- ff[, ord, drop = FALSE]

# -------------------------------
# Compute carrier-proportion windows
# -------------------------------

freq_female <- make_carrier_windows(
  geno = ff, keep_idx = keep_fem, markers = markers,
  window_size = window_size, intro_state = intro_state,
  carrier_threshold = carrier_threshold
) %>%
  smooth_fun(smooth_kb) %>%
  mutate(fillv = pmin(freq_sm, cap), sex = "female")

freq_male <- make_carrier_windows(
  geno = fm, keep_idx = keep_male, markers = markers,
  window_size = window_size, intro_state = intro_state,
  carrier_threshold = carrier_threshold
) %>%
  smooth_fun(smooth_kb) %>%
  mutate(fillv = pmin(freq_sm, cap), sex = "male")

write_tsv(freq_female, file.path(output_dir, paste0("windowed_introgression_carriers_female_DI80_HI075_", window_size/1000, "kb.tsv")))
write_tsv(freq_male, file.path(output_dir, paste0("windowed_introgression_carriers_male_DI80_HI075_", window_size/1000, "kb.tsv")))

chr_summary <- bind_rows(freq_female, freq_male) %>%
  group_by(sex, chr) %>%
  summarise(mean_carrier = mean(freq, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = sex, values_from = mean_carrier) %>%
  mutate(delta_female_minus_male = female - male) %>%
  arrange(desc(delta_female_minus_male))

write_tsv(chr_summary, file.path(output_dir, "chromosome_mean_carrier_frequency_female_vs_male.tsv"))

message("[INFO] Done.")
