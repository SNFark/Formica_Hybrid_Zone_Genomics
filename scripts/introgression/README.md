# Introgression ROI and genome-architecture pipeline

This repository contains the R scripts used to analyse introgression from DIEM output and to test whether introgressed regions are associated with genome architecture.

## Pipeline order

1. `01_generate_candidate_ROIs_from_female_DIEM.R`  
   Generates an automated **candidate** ROI set from female/worker dosage-aware DIEM frequencies. This script is a screening step only.

2. Manual curation outside R  
   Candidate intervals are inspected together with the smoothed female introgression panels. Final ROIs are manually curated by refining boundaries, merging/removing candidates, and adding clear contiguous signals missed by the automatic thresholds. Save the final table as:

   `results/ROIs/curated_ROIs_female_defined.tsv`

   Required columns: `chr`, `start_bp`, `end_bp`. Additional annotation columns are allowed.

3. `02_calculate_windowed_DIEM_carrier_proportions.R`  
   Calculates windowed carrier proportions separately for females/workers and males.

4. `03_compare_curated_ROI_carrier_frequencies_by_sex.R`  
   Uses the curated ROI table to estimate and compare female and male carrier proportions per ROI.

5. `04_test_curated_ROI_architecture.R`  
   Tests the curated ROI set for centromere/self-align association, ROI length patterns, permutation support, centromere-window depletion/enrichment, and inversion overlap.

6. `05_test_FST_centromere_selfalign_enrichment.R`  
   Runs an independent genome-wide FST architecture enrichment test. This does not define or use ROIs.

7. `06_plot_introgression_architecture.R`  
   Plots female windowed carrier-proportion introgression with centromere/self-align and inversion overlays.

## Important note on ROI discovery

The automatic detector was designed to capture two complementary signatures: short localised peaks with high introgression and long contiguous tracts with lower tolerated frequencies. The output of Script 01 is therefore a candidate list, not the final ROI set. Visual curation was necessary because some potentially interesting biology can be missed by fixed thresholds. For example, the final curated set includes a long contiguous tract on chromosome 6 with mean frequency around 0.023, below the automated long-tract threshold.


## Reproducibility

The R scripts assume the working directory is the repository root. Edit the user settings at the top of each script if your files are stored elsewhere.
