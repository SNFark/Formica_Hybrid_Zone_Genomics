#!/usr/bin/env Rscript

# ============================================================
# 06_plot_introgression_architecture.R
#
# Role in pipeline:
#   Plot the female/worker windowed carrier-proportion introgression landscape
#   together with centromeric/self-align regions and high-confidence inversions.
#
# Interpretation:
#   This is a visualisation script. It does not define candidate ROIs, does not
#   curate final ROIs, and does not run statistical tests. It is intended to
#   show how introgression varies across chromosomes in relation to genome
#   architecture.
#
# Main inputs:
#   - results/introgression_windows/windowed_introgression_carriers_female_DI80_HI075_100kb.tsv
#   - CentieR/self-align centromere table
#   - high-confidence inversion table
#
# Main outputs:
#   - figures/introgression_carrier_architecture_female.pdf
#   - figures/introgression_carrier_architecture_female.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

# -------------------------------
# User settings: edit paths here if your files are elsewhere
# -------------------------------

windowed_introgression_file <- "results/introgression_windows/windowed_introgression_carriers_female_DI80_HI075_100kb.tsv"
centromere_file <- "path/to/CentieR_and_selfalign_chromosomes.tsv"
inversion_file <- "path/to/high_confidence_inversions.tsv"

output_dir <- "figures"
output_prefix <- "introgression_carrier_architecture_female"

cap <- 0.30
centromere_track <- "selfalign"  # "selfalign" or "centieR"
exclude_chromosomes <- c("chr3")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# Helper functions
# -------------------------------

normalise_chr <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("FsiP_PB_v5_scf", "chr", x, fixed = TRUE)
  x <- ifelse(grepl("^\\d+$", x), paste0("chr", x), x)
  x <- gsub("^scf", "chr", x)
  x <- gsub("^chr0+", "chr", x)
  x
}

chr_number <- function(x) suppressWarnings(as.integer(sub("^chr", "", as.character(x))))

read_centromeres <- function(path, track = c("selfalign", "centieR")) {
  track <- match.arg(track)
  x <- read_tsv(path, show_col_types = FALSE)
  start_col <- if (track == "selfalign") "High_density_selfalign_start" else "centieR_pred_start"
  end_col <- if (track == "selfalign") "High_density_selfalign_stop" else "centieR_pred_stop"
  x %>% transmute(chr = normalise_chr(chr_fullname),
                  start_mb = as.numeric(.data[[start_col]]) / 1e6,
                  stop_mb = as.numeric(.data[[end_col]]) / 1e6) %>%
    filter(is.finite(start_mb), is.finite(stop_mb), stop_mb > start_mb)
}

read_inversions <- function(path) {
  inv <- read_tsv(path, show_col_types = FALSE)
  nms <- names(inv)
  chr_col <- nms[grepl("^(chr|chrom|chromosome|scaffold|CHROM|Scaffold)$", nms)][1]
  start_col <- nms[grepl("start|START|POS|Start", nms)][1]
  end_col <- nms[grepl("end|END|stop|STOP|End", nms)][1]
  if (any(is.na(c(chr_col, start_col, end_col)))) stop("Could not infer inversion coordinate columns.")
  inv %>% transmute(chr = normalise_chr(.data[[chr_col]]),
                    start = as.numeric(.data[[start_col]]),
                    end = as.numeric(.data[[end_col]])) %>%
    mutate(start_mb = ifelse(start < 1000, start, start / 1e6),
           stop_mb = ifelse(end < 1000, end, end / 1e6)) %>%
    filter(is.finite(start_mb), is.finite(stop_mb), stop_mb > start_mb)
}

# -------------------------------
# Load windowed introgression and architecture tracks
# -------------------------------

freq_df <- read_tsv(windowed_introgression_file, show_col_types = FALSE) %>%
  mutate(chr = normalise_chr(chr),
         fillv = pmin(if ("freq_sm" %in% names(.)) freq_sm else freq, cap)) %>%
  filter(!chr %in% exclude_chromosomes)

lev <- unique(freq_df$chr)
lev <- lev[order(chr_number(lev), lev, na.last = TRUE)]
freq_df <- freq_df %>% mutate(chr = factor(chr, levels = lev))

chr_limits <- freq_df %>% group_by(chr) %>% summarise(xmax_mb = max(end, na.rm = TRUE) / 1e6, .groups = "drop")

cent <- read_centromeres(centromere_file, track = centromere_track) %>%
  filter(chr %in% lev) %>% mutate(chr = factor(chr, levels = lev))

inv <- read_inversions(inversion_file) %>%
  filter(chr %in% lev) %>% mutate(chr = factor(chr, levels = lev))

cent_edge_tiles <- bind_rows(cent %>% transmute(chr, x = start_mb),
                             cent %>% transmute(chr, x = stop_mb))

# -------------------------------
# Plot female carrier-proportion landscape with architecture overlays
# -------------------------------

band_h <- 0.85
cent_bar_w_mb <- 0.03

p <- ggplot() +
  geom_segment(data = chr_limits %>% mutate(x0 = 0),
               aes(x = x0, xend = xmax_mb, y = chr, yend = chr),
               linetype = "dotted", linewidth = 0.3, color = "grey55") +
  geom_tile(data = freq_df,
            aes(x = mid / 1e6, y = chr,
                width = (end - start) / 1e6,
                height = band_h,
                fill = fillv), color = NA) +
  scale_fill_gradientn(
    name = "Carrier freq",
    colours = c(alpha("darkorchid", 0), alpha("darkorchid", 0.9), "darkorchid4"),
    values = c(0, 2/3, 1), limits = c(0, cap), oob = scales::squish,
    breaks = c(0, 0.10, 0.20, 0.30), labels = label_percent(accuracy = 1)
  ) +
  geom_segment(data = inv,
               aes(x = start_mb, xend = stop_mb, y = chr, yend = chr),
               color = "firebrick3", linewidth = 0.9, alpha = 0.9, lineend = "butt") +
  geom_tile(data = cent_edge_tiles,
            aes(x = x, y = chr, width = cent_bar_w_mb, height = band_h),
            fill = "black", alpha = 1) +
  geom_segment(data = cent,
               aes(x = start_mb, xend = stop_mb, y = chr, yend = chr),
               color = "black", linetype = "22", linewidth = 0.9, alpha = 0.6) +
  scale_y_discrete(limits = rev(lev)) +
  scale_x_continuous(name = "Genomic position (Mb)") +
  labs(y = NULL,
       title = "Introgression carrier frequency, centromeric regions, and inversions",
       subtitle = "Carrier frequency capped for visualisation; chr3 excluded") +
  coord_cartesian(clip = "off", expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9),
        legend.position = "top",
        legend.key.width = unit(2.2, "lines"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, paste0(output_prefix, ".pdf")), p, width = 11, height = 8, units = "in")
ggsave(file.path(output_dir, paste0(output_prefix, ".png")), p, width = 11, height = 8, units = "in", dpi = 300)

message("[INFO] Done.")
