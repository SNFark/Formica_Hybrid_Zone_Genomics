#!/usr/bin/env Rscript
# ==============================================================================
# 05_plot_inversions.R
# ==============================================================================
# Plot large filtered inversion calls from a DELLY VCF.
#
# Input:
#   - Filtered inversion VCF from 02_filter_delly_inversions.sh
#   - Optional chromosome/scaffold length table with columns: CHROM, END_POS
#
# Output:
#   - PDF plot of inversion intervals across chromosomes/scaffolds
# ============================================================================== 

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# -------------------------------
# User settings
# -------------------------------

vcf_file <- "path/to/inversions_filtered_Q100_PE14_SR14_100kb.vcf"
scaffold_lengths_file <- "path/to/scaffold_lengths.tsv"  # optional; set to NA if unavailable
output_plot <- "path/to/inversion_plot.pdf"

min_length <- 100000
chrom_levels <- paste0("FsiP_PB_v5_scf", 1:27)

# -------------------------------
# Read VCF
# -------------------------------

vcf_raw <- read.table(
  vcf_file,
  comment.char = "#",
  header = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE
)

vcf <- vcf_raw[, 1:8]
colnames(vcf) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

large_inv <- vcf %>%
  mutate(
    END = as.numeric(sub(".*END=([0-9]+);?.*", "\\1", INFO)),
    LENGTH = END - POS + 1
  ) %>%
  filter(!is.na(END), LENGTH >= min_length) %>%
  filter(CHROM %in% chrom_levels) %>%
  mutate(CHROM = factor(CHROM, levels = chrom_levels)) %>%
  arrange(CHROM, POS)

# -------------------------------
# Optional scaffold end ticks
# -------------------------------

if (!is.na(scaffold_lengths_file) && file.exists(scaffold_lengths_file)) {
  scaffold_ends <- read_tsv(scaffold_lengths_file, show_col_types = FALSE) %>%
    filter(CHROM %in% levels(droplevels(large_inv$CHROM))) %>%
    mutate(
      CHROM = factor(CHROM, levels = levels(droplevels(large_inv$CHROM))),
      y_num = as.numeric(CHROM)
    )
} else {
  scaffold_ends <- tibble()
}

# -------------------------------
# Plot
# -------------------------------

p <- ggplot(large_inv, aes(x = POS, xend = END, y = CHROM, yend = CHROM)) +
  geom_segment(linewidth = 3, alpha = 0.55) +
  scale_x_continuous(labels = comma) +
  labs(
    x = "Genomic position (bp)",
    y = "Chromosome/scaffold",
    title = paste0("Large DELLY inversions (≥", comma(min_length), " bp)")
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))

if (nrow(scaffold_ends) > 0) {
  tick_height <- 0.3
  p <- p +
    geom_segment(
      data = scaffold_ends,
      aes(
        x = END_POS,
        xend = END_POS,
        y = y_num - tick_height / 2,
        yend = y_num + tick_height / 2
      ),
      inherit.aes = FALSE,
      linewidth = 0.5
    )
}

print(p)

ggsave(output_plot, p, width = 10, height = 7)
message("Plot written to: ", output_plot)
