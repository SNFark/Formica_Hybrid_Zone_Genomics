# Data directory

This directory contains the processed input files required to reproduce the introgression analyses presented in the manuscript.

## Contents

### `curated_ROIs_female_defined.tsv`
Final manually curated introgressed regions of interest (ROIs) used in downstream analyses. Coordinates are reported in megabases (Mb). These ROIs were derived from automated candidate detection followed by visual inspection of smoothed introgression profiles.

### `automatic_candidate_ROIs.tsv`
First-pass candidate introgressed regions generated from dosage-aware female DIEM output. These candidates served as a guide for manual curation and were not used directly in downstream analyses.

### `samples_metadata.csv`
Sample identifiers, species assignments, and social genotypes. Sex is inferred from genotype ploidy in the scripts.

### `markers.pos`
Genomic coordinates of diagnostic markers aligned to the DIEM genotype matrices.

### `weir_cockerham_fst.tsv`
Per-site Weir and Cockerham FST estimates used for downstream enrichment analyses.

### DIEM processed inputs
- filtered_geno_80_sex1.RData
- filtered_geno_80_sex2.RData
- hybrid_idx_filtered_80_sex1.RData
- hybrid_idx_filtered_80_sex2.RData
- full_res.RData

These files contain the processed DIEM outputs required to reproduce analyses from ROI discovery onward.

## External resources not redistributed

Centromere and self-alignment annotations originated from collaborator-generated resources used with permission and are therefore not redistributed in this repository.

## Reproducibility workflow

1. Female DIEM output → dosage-aware introgression frequencies.
2. Automated candidate ROI detection.
3. Manual curation of final ROIs.
4. Female and male carrier estimation.
5. Sex comparisons.
6. ROI architecture analyses.
7. FST enrichment analyses.
8. Final visualisation.
