#!/usr/bin/env Rscript

# ============================================================
# 01_generate_candidate_ROIs_from_female_DIEM.R
#
# Role in pipeline:
#   Generate a first-pass set of candidate introgressed regions from
#   female/worker DIEM output. This script does NOT define the final
#   ROI set used in downstream analyses.
#
# Biological/statistical rationale:
#   Candidate regions were screened from the smoothed female dosage-aware
#   introgression landscape in F. cinerea-background individuals.
#   The detector was intentionally designed to capture two different
#   signatures of introgression:
#     1. Localised peaks with relatively high introgression frequency.
#     2. Long contiguous tracts, including chromosome-arm or chromosome-end
#        signals, where introgression may occur at lower population frequency.
#
#   Long tracts were therefore allowed to pass a lower frequency threshold
#   than short peaks. This improves sensitivity to broad, biologically
#   interesting haplotypes that may be carried by only a subset of individuals.
#
# Important reproducibility note:
#   The output of this script is an AUTOMATED CANDIDATE set. Final manuscript
#   ROIs were curated after visual inspection of the smoothed frequency panels.
#   This was necessary because some biologically meaningful contiguous signals
#   can fall below the automated threshold. For example, the final curated ROI
#   set includes a long contiguous tract on chromosome 6 with mean frequency
#   around 0.023, which the automatic thresholding step missed.
#
# Downstream analyses:
#   Scripts 03 and 04 expect a manually curated ROI file, not the candidate
#   file produced here. Save the curated table as:
#       results/ROIs/curated_ROIs_female_defined.tsv
#   with at least chr, start_bp, and end_bp columns.
#
# Main inputs:
#   - filtered_geno_80_sex2.RData              object: genotypes
#   - hybrid_idx_filtered_80_sex2.RData        object: h
#   - markers.pos
#   - full_res.RData                           object: res$DI$DI
#   - samples_metadata.csv
#
# Main outputs:
#   - results/ROIs/female_dosage_aware_windowed_introgression.tsv
#   - results/ROIs/candidate_ROIs_from_female_DIEM.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(zoo)
})

# -------------------------------
# User settings: edit paths here if your files are elsewhere
# -------------------------------

# Repository root. Set this to the directory containing the input files.
base_dir <- "."

genotype_file_worker <- file.path(base_dir, "filtered_geno_80_sex2.RData")
hybrid_index_worker  <- file.path(base_dir, "hybrid_idx_filtered_80_sex2.RData")
metadata_file        <- file.path(base_dir, "samples_metadata.csv")
markers_file         <- file.path(base_dir, "markers.pos")
diem_results_file    <- file.path(base_dir, "full_res.RData")

output_dir <- file.path(base_dir, "results", "ROIs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hi_threshold <- 0.75
# DI quantile: keep markers with DI > 80th percentile, i.e. top 20%.
di_quantile  <- 0.80

window_size <- 100000   # 100 kb windows
smooth_kb   <- 100      # smoothing width in kb; set to 0 for no smoothing

# -------------------------------
# Automated candidate detector settings
# -------------------------------

# Broad-tract detector. Long contiguous regions can be biologically
# meaningful even at lower carrier frequency, so this threshold is lower
# than the threshold used for short peaks.
auto_long_min_mb <- 1.00
auto_mid_min     <- 0.04
auto_gap_windows <- 4
auto_edge_min    <- 0.03

# Local-peak detector. Short regions must exceed a higher threshold
# to avoid calling small stochastic fluctuations as candidates.
auto_short_max_mb <- 0.80
auto_high_min     <- 0.12

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

chr_number <- function(x) {
  suppressWarnings(as.integer(sub("^chr", "", as.character(x))))
}

read_diem_genotypes <- function(path) {
  load(path)  # expected object: genotypes
  if (!exists("genotypes")) {
    stop("File does not contain object named 'genotypes': ", path)
  }

  g <- as.matrix(genotypes)
  g[g == "_"] <- NA
  g <- matrix(as.numeric(g), nrow = nrow(g), ncol = ncol(g))

  rm(genotypes)
  g
}

read_hybrid_index <- function(path) {
  load(path)  # expected object: h
  if (!exists("h")) {
    stop("File does not contain object named 'h': ", path)
  }
  out <- as.numeric(h)
  rm(h)
  out
}

# Dosage-aware conversion for diploid workers/females.
# DIEM coding in this analysis:
#   0 = F. selysi state
#   1 = heterozygous
#   2 = F. cinerea state
selysi_dosage_worker <- function(x) {
  dplyr::case_when(
    x == 0 ~ 1.0,
    x == 1 ~ 0.5,
    x == 2 ~ 0.0,
    TRUE   ~ NA_real_
  )
}

empty_segs <- function() {
  tibble(
    chr = character(),
    start_bp = double(),
    end_bp = double(),
    length_bp = double(),
    mean_freq = double(),
    med_freq = double(),
    class = character()
  )
}

fill_small_gaps <- function(mask, max_gap = 0L) {
  if (!length(mask) || max_gap <= 0) {
    return(ifelse(is.na(mask), FALSE, mask))
  }

  m <- ifelse(is.na(mask), FALSE, mask)
  rl <- rle(m)
  ends <- cumsum(rl$lengths)
  starts <- c(1, head(ends, -1) + 1)

  for (i in seq_along(rl$values)) {
    if (!rl$values[i] && rl$lengths[i] <= max_gap) {
      m[starts[i]:ends[i]] <- TRUE
    }
  }

  m
}

runs_from_mask <- function(dfc, mask_raw) {
  mask <- ifelse(is.na(mask_raw), FALSE, mask_raw)
  if (!any(mask)) return(empty_segs())

  r <- rle(mask)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  keep <- which(r$values)

  out <- lapply(keep, function(i) {
    s <- starts[i]
    e <- ends[i]
    sl <- dfc[s:e, , drop = FALSE]

    tibble(
      chr       = as.character(sl$chr[1]),
      start_bp  = sl$start[1],
      end_bp    = sl$end[nrow(sl)],
      length_bp = sl$end[nrow(sl)] - sl$start[1],
      mean_freq = mean(sl$freq_sm, na.rm = TRUE),
      med_freq  = median(sl$freq_sm, na.rm = TRUE)
    )
  })

  bind_rows(out)
}

trim_edges <- function(dfc, segs, edge_min) {
  if (is.null(segs) || !nrow(segs)) return(empty_segs())

  out <- vector("list", nrow(segs))

  for (i in seq_len(nrow(segs))) {
    s <- segs$start_bp[i]
    e <- segs$end_bp[i]

    sl <- dfc %>%
      filter(.data$start >= s, .data$end <= e)

    keep <- which(!is.na(sl$freq_sm) & sl$freq_sm >= edge_min)
    if (!length(keep)) next

    s2 <- sl$start[min(keep)]
    e2 <- sl$end[max(keep)]

    out[[i]] <- segs[i, , drop = FALSE] %>%
      mutate(
        start_bp  = s2,
        end_bp    = e2,
        length_bp = e2 - s2,
        mean_freq = mean(sl$freq_sm[keep], na.rm = TRUE),
        med_freq  = median(sl$freq_sm[keep], na.rm = TRUE)
      )
  }

  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(empty_segs())

  bind_rows(out)
}

inside_any_long <- function(shorts, longs) {
  if (!nrow(shorts) || !nrow(longs)) return(rep(FALSE, nrow(shorts)))

  vapply(seq_len(nrow(shorts)), function(i) {
    any(longs$chr == shorts$chr[i] &
          shorts$start_bp[i] >= longs$start_bp &
          shorts$end_bp[i]   <= longs$end_bp)
  }, logical(1))
}

# -------------------------------
# Load female/worker DIEM genotypes
# -------------------------------

message("[INFO] Loading worker/female DIEM genotype matrix")
geno_worker <- read_diem_genotypes(genotype_file_worker)

# -------------------------------
# Load marker map and align it to the DI80 genotype matrix
# -------------------------------

message("[INFO] Loading marker positions and DIEM diagnostic index")
markers <- read_delim(
  markers_file,
  delim = "\t",
  col_names = c("chr", "pos"),
  show_col_types = FALSE
)

load(diem_results_file)  # expected object: res
if (!exists("res")) {
  stop("DIEM results file must contain object named 'res'.")
}

di <- res$DI$DI
keep_markers <- di > quantile(di, di_quantile, na.rm = TRUE)

# The genotype matrices used in the analysis were already DI80-filtered.
# If the matrix columns equal the number of top-DI markers, filter the marker map.
# If the marker map already equals the matrix columns, keep it as-is.
if (sum(keep_markers) == ncol(geno_worker)) {
  markers <- markers[keep_markers, , drop = FALSE]
} else if (nrow(markers) == ncol(geno_worker)) {
  warning("[WARN] Genotype matrix and marker map already have matching dimensions; keeping marker rows as-is.")
} else {
  stop("Marker/DI/genotype dimensions do not match. Check whether the matrix and markers are already DI80-filtered.")
}

markers <- markers %>%
  mutate(chr = normalise_chr(chr), pos = as.numeric(pos))

chr_levels <- paste0("chr", 1:27)
markers$chr <- factor(markers$chr, levels = chr_levels)

ord <- order(markers$chr, markers$pos)
markers <- markers[ord, , drop = FALSE]
geno_worker <- geno_worker[, ord, drop = FALSE]

stopifnot(ncol(geno_worker) == nrow(markers))

# -------------------------------
# Load metadata and hybrid index; retain F. cinerea-background females/workers
# -------------------------------

message("[INFO] Loading metadata and worker/female hybrid index")
meta <- read_csv(metadata_file, show_col_types = FALSE)
if (!"sex" %in% names(meta) && "genotype" %in% names(meta)) {
  meta <- meta %>% mutate(sex = nchar(genotype))
}

hi_worker <- read_hybrid_index(hybrid_index_worker)

# If h was stored for all individuals, subset to females/workers.
if (length(hi_worker) == nrow(meta) && "sex" %in% names(meta)) {
  hi_worker <- hi_worker[meta$sex == 2]
}

stopifnot(length(hi_worker) == nrow(geno_worker))

keep_worker <- which(hi_worker >= hi_threshold & !is.na(hi_worker))
message("[INFO] Workers/females kept (HI >= ", hi_threshold, "): ",
        length(keep_worker), " / ", nrow(geno_worker))

# -------------------------------
# Calculate the female dosage-aware windowed introgression profile
# -------------------------------

message("[INFO] Calculating dosage-aware worker/female windowed introgression")

freq_list <- vector("list", length(chr_levels))
names(freq_list) <- chr_levels

for (chr in chr_levels) {
  sel <- markers$chr == chr
  pos <- markers$pos[sel]
  cols <- which(sel)

  if (!length(cols)) next

  maxp <- max(pos, na.rm = TRUE)
  breaks <- seq(0, maxp + window_size, by = window_size)

  w_id <- cut(
    pos,
    breaks = breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = FALSE
  )

  windows <- sort(unique(w_id))

  out_chr <- lapply(windows, function(w) {
    cidx <- cols[w_id == w]
    if (!length(cidx)) return(NULL)

    g_sub <- geno_worker[keep_worker, cidx, drop = FALSE]
    dosage <- selysi_dosage_worker(as.numeric(g_sub))
    freq <- mean(dosage, na.rm = TRUE)

    tibble(
      chr = as.character(chr),
      win = w,
      start = breaks[w],
      end = min(breaks[w] + window_size, maxp),
      mid = (breaks[w] + min(breaks[w] + window_size, maxp)) / 2,
      n_markers = length(cidx),
      freq = freq
    )
  })

  freq_list[[chr]] <- bind_rows(out_chr)
}

freq_worker_dosage <- bind_rows(freq_list)

# Smooth windowed profile.
if (smooth_kb > 0) {
  freq_worker_dosage <- freq_worker_dosage %>%
    group_by(chr) %>%
    arrange(start, .by_group = TRUE) %>%
    mutate(
      wsize = median(end - start, na.rm = TRUE),
      k = max(1L, round((smooth_kb * 1000) / wsize)),
      freq_sm = zoo::rollapply(
        freq,
        width = k,
        FUN = function(x) mean(x, na.rm = TRUE),
        align = "center",
        fill = NA_real_
      )
    ) %>%
    ungroup()
} else {
  freq_worker_dosage$freq_sm <- freq_worker_dosage$freq
}

message("[INFO] Dosage-aware freq_sm summary:")
print(summary(freq_worker_dosage$freq_sm))

# -------------------------------
# Detect first-pass automatic candidate ROIs only
# -------------------------------

message("[INFO] Detecting first-pass automatic dosage-aware ROI candidates")

freq_df_use <- freq_worker_dosage %>%
  mutate(chr = normalise_chr(chr)) %>%
  filter(chr != "chr3", !is.na(freq_sm)) %>%
  arrange(factor(chr, levels = chr_levels), start)

seg_list <- lapply(split(freq_df_use, freq_df_use$chr), function(dfc) {
  dfc <- arrange(dfc, start)

  # Broad elevated regions: no upper frequency cap. In testing, an
  # upper cap fragmented broad signals such as chr16 and chr26. These
  # automated calls are candidates only and are checked visually later.
  long_mask <- dfc$freq_sm >= auto_mid_min
  long_mask <- fill_small_gaps(long_mask, max_gap = auto_gap_windows)

  seg_lm <- runs_from_mask(dfc, long_mask)

  if (nrow(seg_lm)) {
    seg_lm <- seg_lm %>%
      filter(length_bp >= auto_long_min_mb * 1e6)

    seg_lm <- trim_edges(dfc, seg_lm, edge_min = auto_edge_min)

    if (nrow(seg_lm)) {
      seg_lm <- seg_lm %>% mutate(class = "long_elevated")
    }
  } else {
    seg_lm <- empty_segs()
  }

  # Local-peak detector. Short regions must exceed a higher threshold
# to avoid calling small stochastic fluctuations as candidates..
  short_mask <- dfc$freq_sm >= auto_high_min
  seg_sh <- runs_from_mask(dfc, short_mask)

  if (nrow(seg_sh)) {
    seg_sh <- seg_sh %>%
      filter(length_bp <= auto_short_max_mb * 1e6)

    seg_sh <- trim_edges(dfc, seg_sh, edge_min = auto_high_min)

    if (nrow(seg_sh)) {
      seg_sh <- seg_sh %>% mutate(class = "short_high")
    }
  } else {
    seg_sh <- empty_segs()
  }

  # Avoid double-counting: remove short peaks fully nested inside broad tracts.
  if (nrow(seg_sh) && nrow(seg_lm)) {
    seg_sh <- seg_sh[!inside_any_long(seg_sh, seg_lm), , drop = FALSE]
  }

  bind_rows(seg_lm, seg_sh)
})

roi_candidates <- bind_rows(seg_list) %>%
  arrange(factor(chr, levels = chr_levels), start_bp) %>%
  mutate(
    start_mb = start_bp / 1e6,
    end_mb = end_bp / 1e6,
    length_mb = length_bp / 1e6
  )

roi_candidates_mb <- roi_candidates %>%
  transmute(
    chr,
    start_mb = round(start_mb, 2),
    end_mb = round(end_mb, 2),
    length_mb = round(length_mb, 2),
    mean_freq = round(mean_freq, 3),
    med_freq = round(med_freq, 3),
    class
  )

# -------------------------------
# Save automated screening outputs
# -------------------------------

freq_out <- file.path(output_dir, "female_dosage_aware_windowed_introgression.tsv")
roi_out  <- file.path(output_dir, "candidate_ROIs_from_female_DIEM.tsv")

write_tsv(freq_worker_dosage, freq_out)
write_tsv(roi_candidates_mb, roi_out)

message("[INFO] Saved dosage-aware window track: ", freq_out)
message("[INFO] Saved automatic ROI candidates: ", roi_out)
message("[INFO] Candidate ROI count: ", nrow(roi_candidates_mb))

print(roi_candidates_mb, n = Inf)

message("[INFO] Done.")
