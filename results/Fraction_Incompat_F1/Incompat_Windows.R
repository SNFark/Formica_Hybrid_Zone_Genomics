library(tidyverse)
library(readxl)
library(readr)
library(stringr)

# ===============================
# 1) Inputs
# ===============================

chrom_files <- list.files(
  "/Users/sfark/Desktop/Manuscript/data/tsv/",
  pattern = "^chr(?:[0-9]+|X|Y).*\\.tsv\\.gz$",
  full.names = TRUE
)

hi <- readxl::read_excel("/Users/sfark/Desktop/Table_1.xlsx")

fst_raw <- read_tsv(
  "/Users/sfark/Desktop/out.weir.workers.fst.txt",
  show_col_types = FALSE
)

fst <- fst_raw %>%
  rename(
    chrom_raw = CHROM,
    pos = POS,
    fst = WEIR_AND_COCKERHAM_FST
  ) %>%
  mutate(
    chromosome = paste0("chr", chrom_raw),
    window_pos = floor(pos / 10000) * 10000
  ) %>%
  group_by(chromosome, window_pos) %>%
  summarise(
    mean_FST = mean(fst, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(high_FST = mean_FST >= 0.8)

# ===============================
# 2) Parameters
# ===============================

schemes <- list(
  strict_all_F1_sites = c("SMON", "SMAR", "SVIS"),
  no_Monthey = c("SMAR", "SVIS"),
  no_Martigny = c("SMON", "SVIS"),
  no_Visp = c("SMON", "SMAR"),
  Martigny_only = c("SMAR"),
  broad_original = c("SVIS", "SNAT", "SMON", "SMAR", "SSIE", "SSIO")
)

sd_multipliers <- c(1, 1.5, 2)

f1_min <- 0.38
f1_max <- 0.69
pure_selysi_min <- 0.985
pure_cinerea_max <- 0.015

bin_size <- 50000

# ===============================
# 3) Metadata for all relevant workers
# ===============================

all_scheme_locs <- unique(unlist(schemes))

meta_all_workers <- hi %>%
  filter(Caste == "Worker") %>%
  mutate(
    id = ID,
    locality = str_extract(ID, "^[A-Z]+"),
    HI = `Selysi Ancestry`
  ) %>%
  filter(locality %in% all_scheme_locs) %>%
  select(id, locality, HI, Caste, Location, `Supergene Genotype`)

# ===============================
# 4) Load chromosome data ONCE
# ===============================

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
      window_pos = floor(window_pos_raw / 10000) * 10000,
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

# ===============================
# 5) Fast analysis function
# ===============================

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
    select(id, locality, HI, class, Caste, Location, `Supergene Genotype`)
  
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

# ===============================
# 6) Run analyses once per scheme
# ===============================

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


sample_checks <- map_dfr(results, "sample_check")
all_fraction <- map_dfr(results, "fraction_realized")
all_local_means <- map_dfr(results, "local_window_means")

# ===============================
# 5) Basic summaries
# ===============================

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

# ===============================
# 6) Strong incompatible bins
# ===============================

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

# ===============================
# 7) FST binned to same bin size
# ===============================

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

# ===============================
# 8) FST overlap tests
# ===============================

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

# ===============================
# 9) Chromosome-level overlap
# ===============================

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

# ===============================
# 10) Robust strong bins across schemes
# ===============================

robust_strong_bins <- incompat_bins %>%
  filter(sd_multiplier == 2, strong_top5_zero) %>%
  count(chromosome, bin_pos, species, name = "n_schemes_top5") %>%
  arrange(desc(n_schemes_top5), chromosome, bin_pos)

robust_strong_by_chrom <- robust_strong_bins %>%
  count(chromosome, species, wt = n_schemes_top5, name = "robust_score") %>%
  arrange(desc(robust_score))

robust_strong_bins
robust_strong_by_chrom

# ===============================
# 11) Plots
# ===============================

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

# ===============================
# 12) Save outputs
# ===============================

write_csv(sample_checks, "/Users/sfark/Desktop/F1_sample_checks.csv")
write_csv(species_summary, "/Users/sfark/Desktop/F1_species_summary.csv")
write_csv(species_tests, "/Users/sfark/Desktop/F1_species_tests.csv")
write_csv(directional_summary, "/Users/sfark/Desktop/F1_directional_summary.csv")
write_csv(incompat_bins, "/Users/sfark/Desktop/F1_incompatibility_bins.csv")
write_csv(strong_bin_summary, "/Users/sfark/Desktop/F1_strong_bin_summary.csv")
write_csv(incompat_fst_bins, "/Users/sfark/Desktop/F1_incompatibility_FST_bins.csv")
write_csv(fst_overlap_top5_zero, "/Users/sfark/Desktop/F1_FST_overlap_top5_zero.csv")
write_csv(fst_overlap_top10_zero, "/Users/sfark/Desktop/F1_FST_overlap_top10_zero.csv")
write_csv(fst_overlap_top5_mean, "/Users/sfark/Desktop/F1_FST_overlap_top5_mean.csv")
write_csv(chrom_fst_overlap_top5_zero, "/Users/sfark/Desktop/F1_chrom_FST_overlap_top5_zero.csv")
write_csv(robust_strong_bins, "/Users/sfark/Desktop/F1_robust_strong_bins.csv")
write_csv(robust_strong_by_chrom, "/Users/sfark/Desktop/F1_robust_strong_by_chrom.csv")






# Use mean-based strongest incompatibility bins
# Recommended: strict F1 sites, ±2 SD first

fst_species_overlap <- incompat_fst_bins %>%
  filter(
    scheme == "strict_all_F1_sites",
    sd_multiplier == 2
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
    scheme == "strict_all_F1_sites",
    sd_multiplier == 2
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



test_control_delta <- all_local_means %>%
  filter(scheme == "strict_all_F1_sites", sd_multiplier == 2) %>%
  select(chromosome, window_pos, abs_parental_delta) %>%
  left_join(
    all_fraction %>%
      filter(scheme == "strict_all_F1_sites", sd_multiplier == 2) %>%
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
