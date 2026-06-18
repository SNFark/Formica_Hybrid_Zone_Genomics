# ==============================================================================
# F1 hybrid incompatibility analysis
# ==============================================================================
#
# Purpose
# -------
# This script estimates local F1 hybrid compatibility across the genome using
# windowed PC1 scores from parental and F1 individuals.
#
# For each genomic window, the script asks whether F1 individuals fall within the
# parental PC1 range for each parental species. Windows where few or no F1s fall
# within the parental range are treated as candidate regions of reduced hybrid
# compatibility.
#
# The script also:
#   1. compares compatibility estimates across alternative parental sampling schemes;
#   2. aggregates candidate incompatibility signals into genomic bins;
#   3. tests overlap with high-FST regions;
#   4. fits models controlling for FST and parental PC1 distance;
#   5. tests overlap with genomic features such as centromeres or inversions;
#   6. exports summary tables for figures and supplementary analyses.
#
# Input files
# -----------
# This script expects four processed input datasets:
#
# 1. Windowed PCA files
#    One table per chromosome, containing individual IDs and windowed PC1 scores.
#
# 2. Sample metadata
#    A table linking individual IDs to caste, locality, social form or genotype,
#    and genome-wide ancestry estimates.
#
# 3. Genome-wide FST estimates
#    Site-level or windowed Weir and Cockerham FST estimates between parental taxa.
#
# 4. Genome feature coordinates
#    Coordinates of genomic features used for overlap analyses, such as
#    centromeres, self-alignment blocks, or candidate inversions.
#
# Output files
# ------------
# The script writes CSV summary tables to the output directory defined below.
#
# Notes for reuse
# ---------------
# - File paths are defined in the "User settings" block.
# - Use project-relative paths where possible.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(readr)
  library(stringr)
  library(broom)
  library(ggrepel)
})

# ==============================================================================
# 0) User settings
# ==============================================================================
# Define paths to the processed input files used in the analysis.
# These can be absolute paths or paths relative to the working directory.

# Root directory of the project/repository. Change this if running the script
# from outside the repository root.
project_dir <- "."

windowed_pca_dir <- file.path(project_dir, "data", "windowed_PC1_files")
metadata_file    <- file.path(project_dir, "data", "sample_metadata_with_hybrid_indices.xlsx")
fst_file         <- file.path(project_dir, "data", "weir_cockerham_FST.tsv")
feature_file     <- file.path(project_dir, "data", "centromere_selfalignment_coordinates.tsv")
output_dir       <- file.path(project_dir, "results", "F1_compatibility")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis constants used throughout the script.
window_size <- 10000     # Window size used to aggregate FST and windowed PCA coordinates.
bin_size    <- 50000     # Bin size used for candidate incompatibility summaries.


# Hybrid-index thresholds used to classify workers for this analysis.
# HI = proportion of F. selysi ancestry inferred from genome-wide ancestry estimates.
#
# F1 hybrids are expected to have HI close to 0.5. Here, the F1 range is wider
# (0.38–0.69) because observed putative F1 workers did not fall exactly at 0.5.
# This is consistent with the manuscript results showing reduced interclass
# heterozygosity in putative F1s, likely due to historical introgression between
# the parental species and high but imperfect diagnostic differentiation.
#
# Pure parental thresholds are intentionally conservative. Individuals are treated
# as pure F. selysi only when HI > 0.985 and pure F. cinerea only when HI < 0.015,
# reducing the chance that late-generation backcrosses are included in the
# parental reference groups.
f1_min           <- 0.38
f1_max           <- 0.69
pure_selysi_min  <- 0.985
pure_cinerea_max <- 0.015

# SD multipliers used to define the parental PC1 range around each parental mean.
# For example, sd_multiplier = 2 tests whether an F1 is within mean ± 2 SD.
sd_multipliers <- c(1, 1.5, 2)

# Finges colonies used an older labelling scheme and therefore lack the
# standard SFIN prefix. These colony IDs are reassigned to SFIN so that
# locality information is standardised across samples.
finges_colonies <- c(
  "661", "CO54", "FF08", "FF14", "FF15", "FF23",
  "M08", "M13", "M36", "M52"
)

# Alternative parental-location schemes used for sensitivity analyses.
# The main manuscript analysis is "strict_all_F1_sites".
schemes <- list(
  strict_all_F1_sites = c("SMON", "SMAR", "SVIS"),
  no_Monthey = c("SMAR", "SVIS"),
  no_Martigny = c("SMON", "SVIS"),
  no_Visp = c("SMON", "SMAR"),
  Martigny_only = c("SMAR"),
  broad_original = c(
    "SMON", "SEVI", "SMAR", "SFUL", "SSAI", "SRID",
    "SARD", "SSIO", "SLEO", "SSIE", "SFIN", "SVIS", "SNAT"
  )
)

# Plotting/model defaults for the main figures generated below.
plot_scheme <- "strict_all_F1_sites"
plot_sd     <- 2

# Number of chromosome-preserving permutations used in feature-overlap tests.
n_permutations <- 1000

# Inversion set used for the primary feature-overlap analysis.
# Options are "main" for moderate-high/high/very-high confidence calls, or
# "sensitivity" for all moderate-or-higher calls. Low-confidence calls are
# excluded from both sets.
inversion_set <- "main"

# Save plots as PDF files in output_dir.
save_plots <- TRUE


# ==============================================================================
# 1) Inputs
# ==============================================================================
# Locate chromosome-level windowed PCA files and read metadata/FST files.

chrom_files <- list.files(
  windowed_pca_dir,
  pattern = "^chr(?:[0-9]+|X|Y).*\\.tsv\\.gz$",
  full.names = TRUE
)

if (length(chrom_files) == 0) {
  stop("No chromosome files found in: ", windowed_pca_dir)
}

hi <- readxl::read_excel(metadata_file)

fst_raw <- read_tsv(
  fst_file,
  show_col_types = FALSE
)

# Convert per-site or fine-window FST estimates into the same 10-kb coordinate
# system used by the windowed PCA files. Windows with mean FST >= 0.8 are marked
# as highly differentiated.
fst <- fst_raw %>%
  rename(
    chrom_raw = CHROM,
    pos = POS,
    fst = WEIR_AND_COCKERHAM_FST
  ) %>%
  mutate(
    chromosome = paste0("chr", chrom_raw),
    window_pos = floor(pos / window_size) * window_size
  ) %>%
  group_by(chromosome, window_pos) %>%
  summarise(
    mean_FST = mean(fst, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(high_FST = mean_FST >= 0.8)

# ==============================================================================
# 2) Metadata for all relevant workers
# ==============================================================================
# Build a worker-only metadata table with sample ID, colony ID, locality, and
# hybrid index. Locality is extracted from sample IDs, with manual correction for
# specified Finges colonies.


all_scheme_locs <- unique(unlist(schemes))

meta_all_workers <- hi %>%
  filter(Caste == "Worker") %>%
  mutate(
    id = ID,
    colony_id = str_remove(ID, "-.*$"),
    locality_raw = str_extract(ID, "^[A-Z]+"),
    locality = case_when(
      colony_id %in% finges_colonies ~ "SFIN",
      TRUE ~ locality_raw
    ),
    HI = `Selysi Ancestry`
  ) %>%
  filter(locality %in% all_scheme_locs) %>%
  select(id, colony_id, locality, HI, Caste, Location, `Supergene Genotype`)

# ==============================================================================
# 3) Load chromosome data
# ==============================================================================
# Reading all chromosome files.


all_chrom_long_all <- tibble()

for (chrom_file in chrom_files) {
  
  chrom_name <- str_match(
    basename(chrom_file),
    "^(chr(?:[0-9]+|X|Y))"
  )[,2]
  
  message("Loading ", chrom_name)
  
  chrom_long_tmp <- read_tsv(chrom_file, show_col_types = FALSE) %>%
    select(id, where(is.numeric)) %>%
    semi_join(meta_all_workers, by = "id") %>%
    pivot_longer(
      cols = -id,
      names_to = "window",
      values_to = "PC1"
    ) %>%
    mutate(
      window_pos_raw = parse_number(window),
      window_pos = floor(window_pos_raw / window_size) * window_size,
      chromosome = chrom_name
    ) %>%
    filter(!is.na(window_pos)) %>%
    left_join(meta_all_workers, by = "id") %>%
    select(chromosome, window_pos, id, locality, HI, PC1)
  
  all_chrom_long_all <- bind_rows(all_chrom_long_all, chrom_long_tmp)
}

dim(all_chrom_long_all)

object.size(all_chrom_long_all) / 1024^3

all_chrom_long_all %>%
  distinct(chromosome) %>%
  arrange(chromosome)

# ==============================================================================
# 4) F1 compatibility function
# ==============================================================================
# For each scheme, classify workers as F1, pure F. selysi, pure F. cinerea, or
# excluded. For every genomic window and parental species, calculate the fraction
# of F1 individuals whose PC1 value falls inside the parental mean ± SD range.


run_f1_compatibility_fast <- function(parent_locs, scheme_name) {
  
  message("Running scheme: ", scheme_name)
  
  meta_scheme <- meta_all_workers %>%
    mutate(
      class = case_when(
        locality %in% parent_locs & HI >= f1_min & HI <= f1_max ~ "F1",
        locality %in% parent_locs & HI > pure_selysi_min ~ "Selysi",
        locality %in% parent_locs & HI < pure_cinerea_max ~ "Cinerea",
        TRUE ~ "Exclude"
      )
    ) %>%
    filter(class != "Exclude") %>%
    select(id, colony_id, locality, HI, class, Caste, Location, `Supergene Genotype`)
  
  sample_check <- meta_scheme %>%
    count(locality, class) %>%
    pivot_wider(names_from = class, values_from = n, values_fill = 0) %>%
    mutate(scheme = scheme_name)
  
  print(sample_check, n = Inf)
  
  all_chrom_long <- all_chrom_long_all %>%
    semi_join(meta_scheme, by = "id") %>%
    left_join(meta_scheme %>% select(id, class), by = "id")
  
  parent_stats <- all_chrom_long %>%
    filter(class %in% c("Selysi", "Cinerea")) %>%
    group_by(chromosome, window_pos, class) %>%
    summarise(
      mean_PC1 = mean(PC1, na.rm = TRUE),
      sd_PC1 = sd(PC1, na.rm = TRUE),
      n_parent = sum(!is.na(PC1)),
      .groups = "drop"
    ) %>%
    rename(parent_species = class)
  
  f1s <- all_chrom_long %>%
    filter(class == "F1") %>%
    select(chromosome, window_pos, id, PC1)
  
  fraction_base <- parent_stats %>%
    left_join(
      f1s,
      by = c("chromosome", "window_pos"),
      relationship = "many-to-many"
    )
  
  fraction_realized <- map_dfr(sd_multipliers, function(sd_mult) {
    fraction_base %>%
      group_by(chromosome, window_pos, parent_species) %>%
      summarise(
        fraction_realized = mean(
          PC1 >= mean_PC1 - sd_mult * sd_PC1 &
            PC1 <= mean_PC1 + sd_mult * sd_PC1,
          na.rm = TRUE
        ),
        n_F1 = n_distinct(id),
        n_parent = first(n_parent),
        .groups = "drop"
      ) %>%
      rename(species = parent_species) %>%
      mutate(
        scheme = scheme_name,
        sd_multiplier = sd_mult
      )
  })
  
  local_window_means <- all_chrom_long %>%
    group_by(chromosome, window_pos, class) %>%
    summarise(
      mean_PC1 = mean(PC1, na.rm = TRUE),
      sd_PC1 = sd(PC1, na.rm = TRUE),
      n = sum(!is.na(PC1)),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = class,
      values_from = c(mean_PC1, sd_PC1, n)
    ) %>%
    mutate(
      parental_delta = mean_PC1_Selysi - mean_PC1_Cinerea,
      abs_parental_delta = abs(parental_delta),
      directional_score =
        (mean_PC1_F1 - mean_PC1_Cinerea) / parental_delta,
      scheme = scheme_name
    ) %>%
    filter(
      !is.na(directional_score),
      is.finite(directional_score),
      abs_parental_delta > quantile(abs_parental_delta, 0.25, na.rm = TRUE)
    ) %>%
    crossing(sd_multiplier = sd_multipliers)
  
  list(
    sample_check = sample_check,
    fraction_realized = fraction_realized,
    local_window_means = local_window_means
  )
}

# ==============================================================================
# 5) Run analyses once per parental-location scheme
# ==============================================================================


results <- imap(
  schemes,
  ~ run_f1_compatibility_fast(
    parent_locs = .x,
    scheme_name = .y
  )
)

sample_checks <- map_dfr(results, "sample_check")
all_fraction <- map_dfr(results, "fraction_realized")
all_local_means <- map_dfr(results, "local_window_means")

# ==============================================================================
# 6) Basic summaries
# ==============================================================================
# Summarise fraction realised by species and scheme, test species differences,
# and calculate whether F1s are shifted toward F. cinerea or F. selysi.


species_summary <- all_fraction %>%
  group_by(scheme, sd_multiplier, species) %>%
  summarise(
    n_windows = n(),
    mean_fraction_realized = mean(fraction_realized, na.rm = TRUE),
    median_fraction_realized = median(fraction_realized, na.rm = TRUE),
    prop_zero = mean(fraction_realized == 0, na.rm = TRUE),
    prop_one = mean(fraction_realized == 1, na.rm = TRUE),
    .groups = "drop"
  )

species_tests <- all_fraction %>%
  group_by(scheme, sd_multiplier) %>%
  summarise(
    mean_selysi = mean(fraction_realized[species == "Selysi"], na.rm = TRUE),
    mean_cinerea = mean(fraction_realized[species == "Cinerea"], na.rm = TRUE),
    delta_selysi_minus_cinerea = mean_selysi - mean_cinerea,
    wilcox_p = wilcox.test(
      fraction_realized[species == "Selysi"],
      fraction_realized[species == "Cinerea"]
    )$p.value,
    .groups = "drop"
  ) %>%
  arrange(sd_multiplier, scheme)

directional_summary <- all_local_means %>%
  group_by(scheme, sd_multiplier) %>%
  summarise(
    n_windows = n(),
    mean_directional_score = mean(directional_score, na.rm = TRUE),
    median_directional_score = median(directional_score, na.rm = TRUE),
    prop_toward_cinerea = mean(directional_score < 0.5, na.rm = TRUE),
    prop_toward_selysi = mean(directional_score > 0.5, na.rm = TRUE),
    .groups = "drop"
  )

sample_checks
species_summary
species_tests
directional_summary

# ==============================================================================
# 7) Strong candidate incompatibility bins
# ==============================================================================
# Convert window-level compatibility to 50-kb bins. A high incompatibility value
# means few or no F1s fall within the parental PC1 range for that species.


incompat_bins <- all_fraction %>%
  mutate(
    bin_pos = floor(window_pos / bin_size) * bin_size,
    zero_realized = fraction_realized == 0,
    incompatibility = 1 - fraction_realized
  ) %>%
  group_by(scheme, sd_multiplier, chromosome, bin_pos, species) %>%
  summarise(
    prop_zero_realized = mean(zero_realized, na.rm = TRUE),
    mean_incompatibility = mean(incompatibility, na.rm = TRUE),
    mean_fraction_realized = mean(fraction_realized, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  ) %>%
  group_by(scheme, sd_multiplier, species) %>%
  mutate(
    strong_top5_zero = prop_zero_realized >= quantile(prop_zero_realized, 0.95, na.rm = TRUE),
    strong_top10_zero = prop_zero_realized >= quantile(prop_zero_realized, 0.90, na.rm = TRUE),
    strong_top5_mean = mean_incompatibility >= quantile(mean_incompatibility, 0.95, na.rm = TRUE),
    strong_top10_mean = mean_incompatibility >= quantile(mean_incompatibility, 0.90, na.rm = TRUE)
  ) %>%
  ungroup()

strong_bin_summary <- incompat_bins %>%
  group_by(scheme, sd_multiplier, species) %>%
  summarise(
    n_bins = n(),
    n_top5_zero = sum(strong_top5_zero, na.rm = TRUE),
    n_top10_zero = sum(strong_top10_zero, na.rm = TRUE),
    n_top5_mean = sum(strong_top5_mean, na.rm = TRUE),
    n_top10_mean = sum(strong_top10_mean, na.rm = TRUE),
    .groups = "drop"
  )

strong_bin_summary

# ==============================================================================
# 8) FST binned to the same bin size
# ==============================================================================


fst_bins <- fst %>%
  mutate(
    bin_pos = floor(window_pos / bin_size) * bin_size
  ) %>%
  group_by(chromosome, bin_pos) %>%
  summarise(
    mean_FST_bin = mean(mean_FST, na.rm = TRUE),
    max_FST_bin = max(mean_FST, na.rm = TRUE),
    prop_high_FST = mean(high_FST, na.rm = TRUE),
    high_FST_bin = prop_high_FST > 0,
    .groups = "drop"
  )

incompat_fst_bins <- incompat_bins %>%
  left_join(fst_bins, by = c("chromosome", "bin_pos")) %>%
  filter(!is.na(mean_FST_bin))

# ==============================================================================
# 9) Overlap between candidate incompatibility and high-FST bins
# ==============================================================================


fst_overlap_top5_zero <- incompat_fst_bins %>%
  group_by(scheme, sd_multiplier, species) %>%
  summarise(
    n_strong = sum(strong_top5_zero, na.rm = TRUE),
    n_high_FST = sum(high_FST_bin, na.rm = TRUE),
    n_strong_high_FST = sum(strong_top5_zero & high_FST_bin, na.rm = TRUE),
    fisher_p = fisher.test(table(strong_top5_zero, high_FST_bin))$p.value,
    odds_ratio = as.numeric(fisher.test(table(strong_top5_zero, high_FST_bin))$estimate),
    .groups = "drop"
  )

fst_overlap_top10_zero <- incompat_fst_bins %>%
  group_by(scheme, sd_multiplier, species) %>%
  summarise(
    n_strong = sum(strong_top10_zero, na.rm = TRUE),
    n_high_FST = sum(high_FST_bin, na.rm = TRUE),
    n_strong_high_FST = sum(strong_top10_zero & high_FST_bin, na.rm = TRUE),
    fisher_p = fisher.test(table(strong_top10_zero, high_FST_bin))$p.value,
    odds_ratio = as.numeric(fisher.test(table(strong_top10_zero, high_FST_bin))$estimate),
    .groups = "drop"
  )

fst_overlap_top5_mean <- incompat_fst_bins %>%
  group_by(scheme, sd_multiplier, species) %>%
  summarise(
    n_strong = sum(strong_top5_mean, na.rm = TRUE),
    n_high_FST = sum(high_FST_bin, na.rm = TRUE),
    n_strong_high_FST = sum(strong_top5_mean & high_FST_bin, na.rm = TRUE),
    fisher_p = fisher.test(table(strong_top5_mean, high_FST_bin))$p.value,
    odds_ratio = as.numeric(fisher.test(table(strong_top5_mean, high_FST_bin))$estimate),
    .groups = "drop"
  )

fst_overlap_top5_zero
fst_overlap_top10_zero
fst_overlap_top5_mean

# ==============================================================================
# 10) Chromosome-level overlap between candidate incompatibility and high FST
# ==============================================================================


chrom_fst_overlap_top5_zero <- incompat_fst_bins %>%
  group_by(scheme, sd_multiplier, species, chromosome) %>%
  filter(
    n_distinct(strong_top5_zero) > 1,
    n_distinct(high_FST_bin) > 1
  ) %>%
  summarise(
    n_strong = sum(strong_top5_zero, na.rm = TRUE),
    n_high_FST = sum(high_FST_bin, na.rm = TRUE),
    n_strong_high_FST = sum(strong_top5_zero & high_FST_bin, na.rm = TRUE),
    fisher_p = fisher.test(table(strong_top5_zero, high_FST_bin))$p.value,
    odds_ratio = as.numeric(fisher.test(table(strong_top5_zero, high_FST_bin))$estimate),
    .groups = "drop"
  ) %>%
  group_by(scheme, sd_multiplier, species) %>%
  mutate(p_adj_BH = p.adjust(fisher_p, method = "BH")) %>%
  ungroup() %>%
  arrange(scheme, sd_multiplier, species, p_adj_BH)

chrom_fst_overlap_top5_zero

# Same chromosome-level test using the primary mean-incompatibility definition.
chrom_fst_overlap_top5_mean <- incompat_fst_bins %>%
  group_by(scheme, sd_multiplier, species, chromosome) %>%
  filter(
    n_distinct(strong_top5_mean) > 1,
    n_distinct(high_FST_bin) > 1
  ) %>%
  summarise(
    n_strong = sum(strong_top5_mean, na.rm = TRUE),
    n_high_FST = sum(high_FST_bin, na.rm = TRUE),
    n_strong_high_FST = sum(strong_top5_mean & high_FST_bin, na.rm = TRUE),
    fisher_p = fisher.test(table(strong_top5_mean, high_FST_bin))$p.value,
    odds_ratio = as.numeric(fisher.test(table(strong_top5_mean, high_FST_bin))$estimate),
    .groups = "drop"
  ) %>%
  group_by(scheme, sd_multiplier, species) %>%
  mutate(p_adj_BH = p.adjust(fisher_p, method = "BH")) %>%
  ungroup() %>%
  arrange(scheme, sd_multiplier, species, p_adj_BH)

chrom_fst_overlap_top5_mean

# ==============================================================================
# 11) Robust strong bins across sensitivity schemes
# ==============================================================================


robust_strong_bins <- incompat_bins %>%
  filter(sd_multiplier == 2, strong_top5_zero) %>%
  count(chromosome, bin_pos, species, name = "n_schemes_top5") %>%
  arrange(desc(n_schemes_top5), chromosome, bin_pos)

robust_strong_by_chrom <- robust_strong_bins %>%
  count(chromosome, species, wt = n_schemes_top5, name = "robust_score") %>%
  arrange(desc(robust_score))

robust_strong_bins
robust_strong_by_chrom

# ==============================================================================
# 12) Exploratory plots
# ==============================================================================


plot_scheme <- "strict_all_F1_sites"
plot_sd <- 2

plot_bins <- incompat_fst_bins %>%
  filter(
    scheme == plot_scheme,
    sd_multiplier == plot_sd
  )

p_incompat_fst <- ggplot() +
  geom_col(
    data = plot_bins,
    aes(
      x = bin_pos,
      y = mean_incompatibility,
      fill = species
    ),
    position = "identity",
    alpha = 0.45,
    width = bin_size
  ) +
  geom_point(
    data = fst_bins,
    aes(x = bin_pos, y = mean_FST_bin),
    color = "grey60",
    size = 0.25,
    alpha = 0.5
  ) +
  geom_point(
    data = fst_bins %>% filter(high_FST_bin),
    aes(x = bin_pos, y = mean_FST_bin),
    color = "red",
    size = 0.45,
    alpha = 0.8
  ) +
  facet_wrap(~chromosome, scales = "free_x", ncol = 7) +
  theme_minimal() +
  labs(
    x = "Genomic position",
    y = "Mean incompatibility / FST",
    fill = "Species",
    title = paste0(
      "Candidate F1 incompatibility regions and FST: ",
      plot_scheme,
      ", ±", plot_sd, " SD"
    )
  )

p_strong_bins <- ggplot(
  plot_bins,
  aes(x = bin_pos, y = prop_zero_realized, fill = species)
) +
  geom_col(position = "identity", alpha = 0.45, width = bin_size) +
  facet_wrap(~chromosome, scales = "free_x", ncol = 7) +
  theme_minimal() +
  labs(
    x = "Genomic position",
    y = "Proportion zero-realised windows",
    fill = "Species",
    title = "Strongest candidate incompatibility bins"
  )

p_directional <- ggplot(
  all_local_means,
  aes(x = directional_score)
) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  facet_grid(sd_multiplier ~ scheme) +
  theme_minimal() +
  labs(
    x = "F1 position between Cinerea (0) and Selysi (1)",
    y = "Number of windows"
  )

p_incompat_fst
p_strong_bins
p_directional

# ==============================================================================
# 13) Save tabular outputs
# ==============================================================================


write_csv(sample_checks, file.path(output_dir, "F1_sample_checks.csv"))
write_csv(species_summary, file.path(output_dir, "F1_species_summary.csv"))
write_csv(species_tests, file.path(output_dir, "F1_species_tests.csv"))
write_csv(directional_summary, file.path(output_dir, "F1_directional_summary.csv"))
write_csv(incompat_bins, file.path(output_dir, "F1_incompatibility_bins.csv"))
write_csv(strong_bin_summary, file.path(output_dir, "F1_strong_bin_summary.csv"))
write_csv(incompat_fst_bins, file.path(output_dir, "F1_incompatibility_FST_bins.csv"))
write_csv(fst_overlap_top5_zero, file.path(output_dir, "F1_FST_overlap_top5_zero.csv"))
write_csv(fst_overlap_top10_zero, file.path(output_dir, "F1_FST_overlap_top10_zero.csv"))
write_csv(fst_overlap_top5_mean, file.path(output_dir, "F1_FST_overlap_top5_mean.csv"))
write_csv(chrom_fst_overlap_top5_zero, file.path(output_dir, "F1_chrom_FST_overlap_top5_zero.csv"))
write_csv(chrom_fst_overlap_top5_mean, file.path(output_dir, "F1_chrom_FST_overlap_top5_mean.csv"))
write_csv(robust_strong_bins, file.path(output_dir, "F1_robust_strong_bins.csv"))
write_csv(robust_strong_by_chrom, file.path(output_dir, "F1_robust_strong_by_chrom.csv"))






# Use mean-based strongest incompatibility bins
# Recommended: strict F1 sites, ±2 SD first

fst_species_overlap <- incompat_fst_bins %>%
  filter(
    scheme == plot_scheme,
    sd_multiplier == plot_sd
  ) %>%
  group_by(species) %>%
  summarise(
    n_bins = n(),
    n_strong_incompat = sum(strong_top5_mean, na.rm = TRUE),
    n_high_FST = sum(high_FST_bin, na.rm = TRUE),
    n_overlap = sum(strong_top5_mean & high_FST_bin, na.rm = TRUE),
    prop_strong_overlapping_FST = n_overlap / n_strong_incompat,
    fisher_p = fisher.test(table(strong_top5_mean, high_FST_bin))$p.value,
    odds_ratio = as.numeric(fisher.test(table(strong_top5_mean, high_FST_bin))$estimate),
    .groups = "drop"
  )

print(fst_species_overlap, n = Inf)

fst_species_chrom_overlap <- incompat_fst_bins %>%
  filter(
    scheme == plot_scheme,
    sd_multiplier == plot_sd
  ) %>%
  group_by(species, chromosome) %>%
  filter(
    n_distinct(strong_top5_mean) > 1,
    n_distinct(high_FST_bin) > 1
  ) %>%
  summarise(
    n_bins = n(),
    n_strong_incompat = sum(strong_top5_mean, na.rm = TRUE),
    n_high_FST = sum(high_FST_bin, na.rm = TRUE),
    n_overlap = sum(strong_top5_mean & high_FST_bin, na.rm = TRUE),
    prop_strong_overlapping_FST = n_overlap / n_strong_incompat,
    fisher_p = fisher.test(table(strong_top5_mean, high_FST_bin))$p.value,
    odds_ratio = as.numeric(fisher.test(table(strong_top5_mean, high_FST_bin))$estimate),
    .groups = "drop"
  ) %>%
  group_by(species) %>%
  mutate(p_adj_BH = p.adjust(fisher_p, method = "BH")) %>%
  ungroup() %>%
  arrange(species, p_adj_BH)

print(fst_species_chrom_overlap, n = Inf)

write_csv(fst_species_overlap, file.path(output_dir, "F1_FST_species_overlap_top5_mean.csv"))
write_csv(fst_species_chrom_overlap, file.path(output_dir, "F1_FST_species_chrom_overlap_top5_mean.csv"))


test_control_delta <- all_local_means %>%
  filter(scheme == plot_scheme, sd_multiplier == plot_sd) %>%
  select(chromosome, window_pos, abs_parental_delta) %>%
  left_join(
    all_fraction %>%
      filter(scheme == plot_scheme, sd_multiplier == plot_sd) %>%
      mutate(incompatibility = 1 - fraction_realized),
    by = c("chromosome", "window_pos")
  ) %>%
  left_join(fst, by = c("chromosome", "window_pos"))

glm_fit <- glm(
  incompatibility ~ mean_FST + abs_parental_delta + species,
  data = test_control_delta,
  family = gaussian()
)

summary(glm_fit)



glm2 <- glm(
  incompatibility ~ mean_FST * species +
    abs_parental_delta * species,
  data = test_control_delta,
  family = gaussian()
)

summary(glm2)


test_control_delta %>%
  count(scheme, sd_multiplier)





# ===============================
# 1) Build model dataset
# ===============================

plot_scheme <- "strict_all_F1_sites"
plot_sd <- 2
bin_size <- 50000

test_control_delta <- all_local_means %>%
  filter(
    scheme == plot_scheme,
    sd_multiplier == plot_sd
  ) %>%
  select(chromosome, window_pos, abs_parental_delta) %>%
  distinct() %>%
  left_join(
    all_fraction %>%
      filter(
        scheme == plot_scheme,
        sd_multiplier == plot_sd
      ) %>%
      mutate(incompatibility = 1 - fraction_realized),
    by = c("chromosome", "window_pos")
  ) %>%
  left_join(fst, by = c("chromosome", "window_pos")) %>%
  filter(
    !is.na(incompatibility),
    !is.na(mean_FST),
    !is.na(abs_parental_delta),
    species %in% c("Cinerea", "Selysi")
  )

# ===============================
# 2) Fit interaction model
# ===============================

glm_interaction <- glm(
  incompatibility ~ mean_FST * species +
    abs_parental_delta * species,
  data = test_control_delta,
  family = gaussian()
)

summary(glm_interaction)

test_control_delta <- test_control_delta %>%
  mutate(
    predicted_incompatibility = predict(glm_interaction, newdata = .),
    residual_incompatibility = residuals(glm_interaction, type = "response")
  )

# ===============================
# 3) Bin model outputs to 50 kb
# ===============================

model_bins <- test_control_delta %>%
  mutate(
    bin_pos = floor(window_pos / bin_size) * bin_size
  ) %>%
  group_by(chromosome, bin_pos, species) %>%
  summarise(
    mean_incompatibility = mean(incompatibility, na.rm = TRUE),
    mean_predicted_incompatibility = mean(predicted_incompatibility, na.rm = TRUE),
    mean_residual_incompatibility = mean(residual_incompatibility, na.rm = TRUE),
    mean_FST_bin = mean(mean_FST, na.rm = TRUE),
    mean_parental_delta = mean(abs_parental_delta, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  ) %>%
  group_by(species) %>%
  mutate(
    strong_residual_top5 =
      mean_residual_incompatibility >= quantile(mean_residual_incompatibility, 0.95, na.rm = TRUE),
    strong_raw_top5 =
      mean_incompatibility >= quantile(mean_incompatibility, 0.95, na.rm = TRUE)
  ) %>%
  ungroup()

# chromosome ordering
chrom_levels <- paste0("chr", c(1:27))
model_bins <- model_bins %>%
  mutate(chromosome = factor(chromosome, levels = chrom_levels))

test_control_delta <- test_control_delta %>%
  mutate(chromosome = factor(chromosome, levels = chrom_levels))

# ===============================
# 4) Plot A: FST vs incompatibility with species-specific slopes
# ===============================

p_fst_incompat <- ggplot(
  test_control_delta,
  aes(x = mean_FST, y = incompatibility, color = species)
) +
  geom_point(alpha = 0.035, size = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1.1) +
  facet_wrap(~species) +
  theme_minimal(base_size = 12) +
  labs(
    x = expression("Mean " * F[ST]),
    y = "Incompatibility (1 - fraction realised)",
    color = "Parental species",
    title = "Differentiated regions show elevated candidate incompatibility",
    subtitle = paste0(plot_scheme, ", ±", plot_sd, " SD")
  )

p_fst_incompat





# ===============================
# 5) Plot B: FST vs incompatibility, coloured by parental divergence
# ===============================

p_fst_incompat_delta <- ggplot(
  test_control_delta,
  aes(
    x = mean_FST,
    y = incompatibility,
    color = abs_parental_delta
  )
) +
  geom_point(alpha = 0.08, size = 0.35) +
  geom_smooth(
    aes(group = species),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "black",
    linewidth = 0.9
  ) +
  facet_wrap(~species) +
  theme_minimal(base_size = 12) +
  labs(
    x = expression("Mean " * F[ST]),
    y = "Incompatibility (1 - fraction realised)",
    color = "Parental PC1 distance",
    title = "FST-incompatibility association after visualising parental divergence",
    subtitle = "Black lines show species-specific linear fits"
  )

p_fst_incompat_delta



# ===============================
# 6) Plot C: Raw incompatibility landscape + FST
# ===============================

p_raw_landscape <- ggplot() +
  geom_col(
    data = model_bins,
    aes(
      x = bin_pos,
      y = mean_incompatibility,
      fill = species
    ),
    position = "identity",
    alpha = 0.45,
    width = bin_size
  ) +
  geom_point(
    data = model_bins,
    aes(x = bin_pos, y = mean_FST_bin),
    color = "grey35",
    alpha = 0.45,
    size = 0.25
  ) +
  facet_wrap(~chromosome, scales = "free_x", ncol = 7) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Genomic position",
    y = expression("Mean incompatibility / mean " * F[ST]),
    fill = "Parental species",
    title = "Raw candidate incompatibility landscape",
    subtitle = paste0(plot_scheme, ", ±", plot_sd, " SD")
  )

p_raw_landscape


# ===============================
# 7) Plot D: Residual incompatibility landscape
# After accounting for FST, parental divergence, species, and interactions
# ===============================

p_residual_landscape <- ggplot(
  model_bins,
  aes(
    x = bin_pos,
    y = mean_residual_incompatibility,
    fill = species
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  geom_col(
    position = "identity",
    alpha = 0.45,
    width = bin_size
  ) +
  facet_wrap(~chromosome, scales = "free_x", ncol = 7) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Genomic position",
    y = "Residual incompatibility",
    fill = "Parental species",
    title = "Residual candidate incompatibility landscape",
    subtitle = "Residuals after accounting for FST, parental divergence, species, and interactions"
  )

p_residual_landscape



# ===============================
# 8) Plot E: Strong residual outlier bins
# ===============================

top_residual_bins <- model_bins %>%
  filter(strong_residual_top5) %>%
  arrange(species, desc(mean_residual_incompatibility))

p_residual_outliers <- ggplot(
  model_bins,
  aes(x = bin_pos, y = mean_residual_incompatibility)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  geom_col(
    aes(fill = species),
    alpha = 0.30,
    width = bin_size
  ) +
  geom_point(
    data = top_residual_bins,
    aes(x = bin_pos, y = mean_residual_incompatibility),
    size = 0.7
  ) +
  facet_wrap(~chromosome, scales = "free_x", ncol = 7) +
  theme_minimal(base_size = 10) +
  labs(
    x = "Genomic position",
    y = "Residual incompatibility",
    fill = "Parental species",
    title = "Top residual incompatibility bins",
    subtitle = "Points mark top 5% residual incompatibility bins within each species"
  )

p_residual_outliers





# ===============================
# 9) Plot F: Observed vs predicted incompatibility
# ===============================

p_observed_predicted <- ggplot(
  model_bins,
  aes(
    x = mean_predicted_incompatibility,
    y = mean_incompatibility,
    color = species
  )
) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~species) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Predicted incompatibility",
    y = "Observed incompatibility",
    color = "Parental species",
    title = "Observed versus model-predicted candidate incompatibility",
    subtitle = "Bins above the dashed line are more incompatible than expected"
  )

p_observed_predicted


# ===============================
# 10) Plot G: Label top residual outliers
# ===============================

label_bins <- model_bins %>%
  group_by(species) %>%
  slice_max(mean_residual_incompatibility, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    label = paste0(chromosome, ":", round(bin_pos / 1e6, 2), " Mb")
  )

p_observed_predicted_labeled <- ggplot(
  model_bins,
  aes(
    x = mean_predicted_incompatibility,
    y = mean_incompatibility,
    color = species
  )
) +
  geom_point(alpha = 0.20, size = 0.75) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = label_bins,
    aes(label = label),
    size = 2.4,
    max.overlaps = 50,
    show.legend = FALSE
  ) +
  facet_wrap(~species) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Predicted incompatibility",
    y = "Observed incompatibility",
    color = "Parental species",
    title = "Residual outlier candidate incompatibility regions",
    subtitle = "Labels indicate bins with highest positive residual incompatibility"
  )

p_observed_predicted_labeled

# ===============================
# 11) Plot H: Model coefficients
# ===============================

coef_df <- broom::tidy(glm_interaction, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_clean = recode(
      term,
      "mean_FST" = "FST",
      "speciesSelysi" = "Selysi main effect",
      "abs_parental_delta" = "Parental PC1 distance",
      "mean_FST:speciesSelysi" = "FST × Selysi",
      "speciesSelysi:abs_parental_delta" = "Parental distance × Selysi"
    )
  )

p_model_coefficients <- ggplot(
  coef_df,
  aes(x = estimate, y = reorder(term_clean, estimate))
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Model estimate",
    y = NULL,
    title = "Predictors of candidate incompatibility",
    subtitle = "Gaussian model with species interactions"
  )

p_model_coefficients

# Save model-level outputs and figures.
write_csv(test_control_delta, file.path(output_dir, "F1_model_window_level_data.csv"))
write_csv(model_bins, file.path(output_dir, "F1_model_50kb_bins.csv"))
write_csv(coef_df, file.path(output_dir, "F1_model_coefficients.csv"))

if (isTRUE(save_plots)) {
  ggsave(file.path(output_dir, "F1_plot_incompatibility_FST.pdf"),
         p_incompat_fst, width = 14, height = 8)
  ggsave(file.path(output_dir, "F1_plot_zero_realized_bins.pdf"),
         p_strong_bins, width = 14, height = 8)
  ggsave(file.path(output_dir, "F1_plot_directional_score.pdf"),
         p_directional, width = 12, height = 8)
  ggsave(file.path(output_dir, "F1_plot_FST_vs_incompatibility.pdf"),
         p_fst_incompat, width = 8, height = 4.5)
  ggsave(file.path(output_dir, "F1_plot_FST_vs_incompatibility_parental_delta.pdf"),
         p_fst_incompat_delta, width = 8, height = 4.5)
  ggsave(file.path(output_dir, "F1_plot_raw_incompatibility_landscape.pdf"),
         p_raw_landscape, width = 14, height = 8)
  ggsave(file.path(output_dir, "F1_plot_residual_incompatibility_landscape.pdf"),
         p_residual_landscape, width = 14, height = 8)
  ggsave(file.path(output_dir, "F1_plot_residual_outlier_bins.pdf"),
         p_residual_outliers, width = 14, height = 8)
  ggsave(file.path(output_dir, "F1_plot_observed_vs_predicted.pdf"),
         p_observed_predicted, width = 8, height = 4.5)
  ggsave(file.path(output_dir, "F1_plot_observed_vs_predicted_labeled.pdf"),
         p_observed_predicted_labeled, width = 8, height = 4.5)
  ggsave(file.path(output_dir, "F1_plot_model_coefficients.pdf"),
         p_model_coefficients, width = 7, height = 4.5)
}

# ==============================================================================
# 14) Genome-feature overlap analyses
# ==============================================================================
# This section tests whether candidate incompatibility bins overlap centromeres,
# high-density self-alignment blocks, or candidate inversions.

# ------------------------------------------------------------------------------
# 14.1 Read CentieR + self-alignment feature coordinates
# ------------------------------------------------------------------------------


centier_selfalign <- read_tsv(
  feature_file,
  show_col_types = FALSE
) %>%
  mutate(
    chromosome = paste0("chr", chr_number)
  )

centromere_intervals <- centier_selfalign %>%
  transmute(
    chromosome,
    start = centieR_pred_start,
    end = centieR_pred_stop,
    feature = "centromere"
  )

selfalign_intervals <- centier_selfalign %>%
  transmute(
    chromosome,
    start = High_density_selfalign_start,
    end = High_density_selfalign_stop,
    feature = "selfalign_block"
  )


# ------------------------------------------------------------------------------
# 14.2 Define candidate inversion intervals
# ------------------------------------------------------------------------------


## Main analysis:
## only robust inversions with strong support
## (Moderate-High, High, Very High)

inversion_main <- tribble(
  ~chromosome, ~start, ~end, ~feature,
  
  # chr8 large inversion
  "chr8",  8691765, 11243911, "inversion",
  
  # chr10 very high confidence
  "chr10", 6242232, 10864328, "inversion",
  
  # chr15 merged inversion block
  "chr15", 6775634, 10514107, "inversion"
)

## Sensitivity analysis:
## include all Moderate or stronger calls
## (Low-confidence calls excluded)

inversion_sensitivity <- tribble(
  ~chromosome, ~start, ~end, ~feature,
  
  # chr1 merged moderate calls
  "chr1",  2797512, 4338202, "inversion",
  
  # chr8 moderate calls
  "chr8",  6166092, 6467326, "inversion",
  "chr8",  7937535, 8106788, "inversion",
  
  # chr8 large moderate-high inversion
  "chr8",  8691765, 11243911, "inversion",
  
  # chr9 moderate
  "chr9",  8439398, 10065007, "inversion",
  
  # chr10 very high
  "chr10", 6242232, 10864328, "inversion",
  
  # chr15 merged moderate/high block
  "chr15", 6775634, 10514107, "inversion",
  
  # chr16 merged moderate calls
  "chr16", 3566041, 5440245, "inversion",
  
  # chr25 moderate
  "chr25", 2568006, 2733241, "inversion",
  
  # chr26 moderate
  "chr26",  955823, 2668657, "inversion"
)

## Choose which inversion set to analyse.
inversion_intervals <- switch(
  inversion_set,
  main = inversion_main,
  sensitivity = inversion_sensitivity,
  stop("inversion_set must be either 'main' or 'sensitivity'.")
)



# ------------------------------------------------------------------------------
# 14.3 Combine all genomic features
# ------------------------------------------------------------------------------


feature_intervals <- bind_rows(
  centromere_intervals,
  selfalign_intervals,
  inversion_intervals
)

# ------------------------------------------------------------------------------
# 14.4 Build candidate bins from model outputs
# ------------------------------------------------------------------------------


candidate_bins <- model_bins %>%
  mutate(
    bin_start = bin_pos,
    bin_end = bin_pos + bin_size
  ) %>%
  group_by(species) %>%
  mutate(
    raw_top5 = mean_incompatibility >= quantile(mean_incompatibility, 0.95, na.rm = TRUE),
    pos_resid_top5 = mean_residual_incompatibility >= quantile(mean_residual_incompatibility, 0.95, na.rm = TRUE),
    neg_resid_top5 = mean_residual_incompatibility <= quantile(mean_residual_incompatibility, 0.05, na.rm = TRUE)
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 14.5 Helper function to mark interval overlap
# ------------------------------------------------------------------------------


add_overlap <- function(bins, intervals, feature_name) {
  
  intervals_use <- intervals %>%
    filter(feature == feature_name)
  
  bins %>%
    rowwise() %>%
    mutate(
      overlap_feature = any(
        intervals_use$chromosome == chromosome &
          intervals_use$start < bin_end &
          intervals_use$end > bin_start
      )
    ) %>%
    ungroup() %>%
    mutate(feature = feature_name)
}

# ------------------------------------------------------------------------------
# 14.6 Fisher overlap tests
# ------------------------------------------------------------------------------


run_fisher_overlap <- function(feature_name) {
  
  x <- add_overlap(candidate_bins, feature_intervals, feature_name)
  
  x %>%
    group_by(species) %>%
    summarise(
      feature = feature_name,
      
      raw_overlap = sum(raw_top5 & overlap_feature, na.rm = TRUE),
      raw_total = sum(raw_top5, na.rm = TRUE),
      raw_p = fisher.test(table(raw_top5, overlap_feature))$p.value,
      raw_OR = as.numeric(fisher.test(table(raw_top5, overlap_feature))$estimate),
      
      pos_resid_overlap = sum(pos_resid_top5 & overlap_feature, na.rm = TRUE),
      pos_resid_total = sum(pos_resid_top5, na.rm = TRUE),
      pos_resid_p = fisher.test(table(pos_resid_top5, overlap_feature))$p.value,
      pos_resid_OR = as.numeric(fisher.test(table(pos_resid_top5, overlap_feature))$estimate),
      
      neg_resid_overlap = sum(neg_resid_top5 & overlap_feature, na.rm = TRUE),
      neg_resid_total = sum(neg_resid_top5, na.rm = TRUE),
      neg_resid_p = fisher.test(table(neg_resid_top5, overlap_feature))$p.value,
      neg_resid_OR = as.numeric(fisher.test(table(neg_resid_top5, overlap_feature))$estimate),
      
      .groups = "drop"
    )
}

feature_overlap_results <- map_dfr(
  c("centromere", "selfalign_block", "inversion"),
  run_fisher_overlap
)

print(feature_overlap_results, n = Inf)

# ------------------------------------------------------------------------------
# 14.7 Chromosome-level feature-overlap summaries
# ------------------------------------------------------------------------------


feature_overlap_by_chrom <- map_dfr(
  c("centromere", "selfalign_block", "inversion"),
  function(feature_name) {
    
    x <- add_overlap(candidate_bins, feature_intervals, feature_name)
    
    x %>%
      group_by(feature, species, chromosome) %>%
      summarise(
        n_bins = n(),
        n_feature_bins = sum(overlap_feature, na.rm = TRUE),
        raw_top5_in_feature = sum(raw_top5 & overlap_feature, na.rm = TRUE),
        pos_resid_top5_in_feature = sum(pos_resid_top5 & overlap_feature, na.rm = TRUE),
        neg_resid_top5_in_feature = sum(neg_resid_top5 & overlap_feature, na.rm = TRUE),
        .groups = "drop"
      )
  }
)

print(feature_overlap_by_chrom, n = Inf)

# ------------------------------------------------------------------------------
# 14.8 Chromosome-preserving permutation tests
# ------------------------------------------------------------------------------
# Candidate labels are permuted within chromosomes to preserve chromosome-level
# differences in window number and genomic architecture.


perm_overlap <- function(df, candidate_col, n_perm = 1000) {
  
  observed <- sum(df[[candidate_col]] & df$overlap_feature, na.rm = TRUE)
  
  null <- replicate(n_perm, {
    df %>%
      group_by(chromosome) %>%
      mutate(candidate_perm = sample(.data[[candidate_col]])) %>%
      ungroup() %>%
      summarise(overlap = sum(candidate_perm & overlap_feature, na.rm = TRUE)) %>%
      pull(overlap)
  })
  
  tibble(
    observed_overlap = observed,
    null_mean = mean(null),
    null_sd = sd(null),
    empirical_p_high = mean(null >= observed),
    empirical_p_low = mean(null <= observed)
  )
}

run_perm_feature <- function(feature_name, n_perm = 1000) {
  
  x <- add_overlap(candidate_bins, feature_intervals, feature_name)
  
  x %>%
    group_by(species) %>%
    group_modify(~ bind_rows(
      perm_overlap(.x, "raw_top5", n_perm) %>%
        mutate(test = "raw_top5"),
      perm_overlap(.x, "pos_resid_top5", n_perm) %>%
        mutate(test = "positive_residual_top5"),
      perm_overlap(.x, "neg_resid_top5", n_perm) %>%
        mutate(test = "negative_residual_top5")
    )) %>%
    ungroup() %>%
    mutate(feature = feature_name)
}

feature_perm_results <- map_dfr(
  c("centromere", "selfalign_block", "inversion"),
  run_perm_feature,
  n_perm = n_permutations
)

print(feature_perm_results, n = Inf)

# Save feature-overlap outputs.
write_csv(feature_overlap_results, file.path(output_dir, "F1_feature_overlap_fisher.csv"))
write_csv(feature_overlap_by_chrom, file.path(output_dir, "F1_feature_overlap_by_chromosome.csv"))
write_csv(feature_perm_results, file.path(output_dir, "F1_feature_overlap_permutations.csv"))
message("F1 compatibility analysis complete. Outputs written to: ", output_dir)
