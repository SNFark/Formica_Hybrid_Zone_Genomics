# Formica Hybrid Zone Genomics

This repository contains scripts, processed results, and associated data resources accompanying our study of introgression, hybrid incompatibilities, and structural variation in the *Formica selysi × Formica cinerea* hybrid zone.

The repository provides the workflows required to reproduce the analyses presented in the manuscript from processed analytical inputs onward. Large input datasets that exceeded GitHub file size limits are archived on Zenodo.

## Repository structure

### `data/`
Processed input files and documentation describing the datasets used in each analysis. Larger analytical inputs are archived on Zenodo and linked through the corresponding README files.

### `scripts/`
Analysis workflows used to generate the results presented in the manuscript.

- `introgression/` – identification and characterisation of introgressed regions of interest (ROIs) using DIEM-derived outputs.
- `F1_incompatibility/` – analyses quantifying realised parental variation in F₁ hybrids using local PCA scores.
- `structural_variants/` – identification and visualisation of candidate inversions.

### `results/`
Processed outputs underlying manuscript figures and statistical analyses, including summary tables and intermediate analytical results.

## Zenodo archive

Large processed datasets required to reproduce these analyses are available from Zenodo:

**DOI:** https://doi.org/10.5281/zenodo.20762532

The Zenodo archive includes:

- processed DIEM inputs used for introgression analyses;
- chromosome-specific windowed PCA files used in F₁ incompatibility analyses;
- additional large input files exceeding GitHub file size limits.

## Notes

This repository contains processed analytical inputs rather than raw sequencing data. Certain collaborator-generated resources (e.g. centromere and self-alignment annotations) were used with permission and are therefore not redistributed.

## Citation

If using these scripts or datasets, please cite both the associated manuscript and the Zenodo archive:

> DOI: https://doi.org/10.5281/zenodo.20762532
