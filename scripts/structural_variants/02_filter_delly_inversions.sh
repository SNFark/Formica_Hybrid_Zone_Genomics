#!/usr/bin/env bash
#SBATCH --job-name=delly_filter
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=delly_filter_%j.log

# ==============================================================================
# 02_filter_delly_inversions.sh
# ==============================================================================
# Filter DELLY inversion calls using the high-confidence criteria used in the
# manuscript analysis.
#
# Original filtering criteria:
#   SVTYPE = INV
#   QUAL >= 100
#   INFO/PE >= 14
#   INFO/SR >= 14
#   inversion length >= 100 kb
#
# Scaffold/chromosome 3 can optionally be removed because it corresponds to the
# known social-supergene region in this study.
# ============================================================================== 

set -euo pipefail

# -------------------------------
# User settings
# -------------------------------

INPUT_VCF="path/to/inversions_genotyped.vcf"
OUTDIR="path/to/delly_output"
FILTERED_VCF="$OUTDIR/inversions_filtered_Q100_PE14_SR14_100kb.vcf"

# Set to "TRUE" to omit chromosome/scaffold 3 from the filtered output.
EXCLUDE_CHR3="TRUE"
CHR3_NAME="FsiP_PB_v5_scf3"

# Uncomment/adapt for your cluster environment if needed.
# module load bcftools

mkdir -p "$OUTDIR"

if [[ ! -f "$INPUT_VCF" ]]; then
  echo "[ERROR] Input VCF not found: $INPUT_VCF" >&2
  exit 1
fi

FILTER_EXPR='SVTYPE="INV" && QUAL>=100 && INFO/PE>=14 && INFO/SR>=14 && (INFO/END - POS + 1)>=100000'

if [[ "$EXCLUDE_CHR3" == "TRUE" ]]; then
  echo "[INFO] Filtering inversions and excluding $CHR3_NAME"
  bcftools view \
    -i "$FILTER_EXPR" \
    -T ^<(echo "$CHR3_NAME") \
    "$INPUT_VCF" \
    -Ov \
    -o "$FILTERED_VCF"
else
  echo "[INFO] Filtering inversions"
  bcftools view \
    -i "$FILTER_EXPR" \
    "$INPUT_VCF" \
    -Ov \
    -o "$FILTERED_VCF"
fi

echo "[INFO] Filtered VCF saved to: $FILTERED_VCF"

echo "[INFO] Inversion counts per chromosome/scaffold"
bcftools query -f '%CHROM\n' "$FILTERED_VCF" | sort | uniq -c | sort -nr
