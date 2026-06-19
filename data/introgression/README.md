# Data directory

This directory contains processed input files required to reproduce the introgression analyses presented in the manuscript.

Some larger input datasets exceeded GitHub file size limits and are therefore archived on Zenodo:

DOI: https://doi.org/10.5281/zenodo.20762532

## Contents

### curated_ROIs_female_defined.tsv
Final manually curated introgressed regions of interest (ROIs) used in downstream analyses. Coordinates are reported in megabases (Mb). These ROIs were derived from automated candidate detection followed by visual inspection of smoothed introgression profiles.

### samples_metadata.csv
Sample identifiers, species assignments, and social genotypes. Sex is inferred from genotype ploidy in the scripts.

## Additional input files archived on Zenodo

The following processed inputs are available from Zenodo (DOI: https://doi.org/10.5281/zenodo.20762532):

### Introgression inputs
Contained in:

`Formica_Hybrid_Zone_Genomics_Introgression_Inputs.zip`

- `markers.pos` – genomic coordinates of diagnostic markers aligned to the DIEM genotype matrices.
- `weir_cockerham_fst.tsv` – per-site Weir and Cockerham FST estimates used for downstream enrichment analyses.
- `filtered_geno_80_sex1.RData`
- `filtered_geno_80_sex2.RData`
- `hybrid_idx_filtered_80_sex1.RData`
- `hybrid_idx_filtered_80_sex2.RData`
- `full_res.RData`

These files contain the processed DIEM outputs and associated inputs required to reproduce analyses from ROI discovery onward.

## External resources not redistributed

Centromere and self-alignment annotations originated from collaborator-generated resources used with permission and are therefore not redistributed through either this repository or the Zenodo archive.

## Reproducibility workflow

1. Female DIEM output → dosage-aware introgression frequencies.
2. Automated candidate ROI detection.
3. Manual curation of final ROIs based on smoothed introgression profiles.
4. Female and male carrier estimation using curated ROIs.
5. Sex comparisons.
6. ROI architecture analyses.
7. FST enrichment analyses.
8. Final visualisation.

## Notes

The final curated ROI set was not taken directly from the automated candidate output. Automated detection was designed to identify both localised peaks of elevated introgression and extended introgressed tracts. Candidate regions were subsequently evaluated through visual inspection of smoothed introgression profiles to retain biologically meaningful regions that were missed by threshold-based approaches.
