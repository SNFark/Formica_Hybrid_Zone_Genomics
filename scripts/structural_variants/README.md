# Structural variant analysis scripts

These scripts reproduce the inversion-calling and inspection workflow used for the *Formica* hybrid zone manuscript.

## Overview

The analysis consisted of four steps:

1. Jointly call and genotype inversions with DELLY from high-coverage BAM files.
2. Filter candidate inversions using stringent evidence thresholds.
3. Inspect candidate inversions with Samplot in any representative high-coverage BAM and count read-pair orientations.
4. Check for coverage dips near inversion breakpoints and plot retained inversion intervals.

The original analysis used DELLY calls from 42 high-coverage BAM files: 30 males and 12 females. The scripts here are generalised so they can run on any number of BAMs listed in `bam_list.txt`.

Scaffold/chromosome 3 was excluded from the final inspection because it corresponds to the known social-supergene region.

## Scripts

### `01_delly_call_inversions.sh`

Jointly calls inversions with DELLY and genotypes calls across all BAMs in a user-provided `bam_list.txt`.

Required inputs:

- reference genome FASTA
- one sorted/indexed BAM path per line in `bam_list.txt`

Main outputs:

- `inversions_raw.bcf`
- `inversions_genotyped.bcf`
- `inversions_genotyped.vcf`

### `02_filter_delly_inversions.sh`

Filters DELLY calls using the criteria used in the manuscript:

- `SVTYPE = INV`
- `QUAL >= 100`
- `INFO/PE >= 14`
- `INFO/SR >= 14`
- inversion length `>= 100 kb`

Optional exclusion of scaffold 3 is built in.

### `03_samplot_inspection.py`

Generates a Samplot image for each filtered inversion and counts FR/RF/Other read-pair orientations in one representative BAM.

Original analysis used `SVIS10-M2_markdup.sorted.bam` as the representative high-coverage BAM.

Example run:

```bash
source ~/samplot_venv/bin/activate
python 03_samplot_inspection.py > inversion_summary.txt
```

### `04_breakpoint_coverage_dips.py`

Reads `inversion_summary.txt`, checks a 1 kb window at each inversion breakpoint, and reports whether mean coverage falls below a threshold. The original analysis used a coverage threshold of 5 reads.

Example run:

```bash
python 04_breakpoint_coverage_dips.py
```

### `05_plot_inversions.R`

Plots retained inversion intervals from the filtered VCF.

## Confidence interpretation

Final inversion confidence categories were based on combined evidence from:

- DELLY support after filtering;
- Samplot visual inspection;
- FR/RF/Other orientation counts in the representative BAM;
- breakpoint coverage dip summaries.

Low-confidence calls were not used in the main manuscript analyses. The main incompatibility-overlap analysis used only Moderate-High, High, and Very High confidence inversions, with Moderate calls used only in sensitivity analyses.

## Notes for reuse

All scripts use generic paths. Replace paths in the `User settings` section before running.

Large BAM, VCF, BCF, and FASTA files should not be committed to GitHub. Archive large datasets separately, for example on Zenodo or a sequencing archive.
