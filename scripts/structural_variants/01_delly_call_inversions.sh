#!/usr/bin/env bash
#SBATCH --job-name=delly_inversions
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=delly_inversions_%j.log

# ==============================================================================
# 01_delly_call_inversions.sh
# ==============================================================================
# Jointly call and genotype inversions with DELLY from an arbitrary list of BAMs.
#
# Input:
#   - Reference genome FASTA
#   - Text file containing one sorted, indexed BAM path per line
#
# Output:
#   - Raw DELLY inversion BCF
#   - Genotyped inversion BCF
#   - Genotyped inversion VCF
#
# Notes:
# This script is intentionally generalised for reuse. The script will run on however many BAMs are 
# listed in BAM_LIST.
# ==============================================================================

set -euo pipefail

# -------------------------------
# User settings
# -------------------------------

REF="path/to/reference_genome.fasta"
BAM_LIST="path/to/bam_list.txt"
OUTDIR="path/to/delly_output"
DELLY="delly"

# Uncomment/adapt for your cluster environment if needed.
# module load samtools
# module load bcftools

mkdir -p "$OUTDIR"

RAW_BCF="$OUTDIR/inversions_raw.bcf"
GENOTYPED_BCF="$OUTDIR/inversions_genotyped.bcf"
GENOTYPED_VCF="$OUTDIR/inversions_genotyped.vcf"

# -------------------------------
# Checks
# -------------------------------

if [[ ! -f "$REF" ]]; then
  echo "[ERROR] Reference FASTA not found: $REF" >&2
  exit 1
fi

if [[ ! -f "$BAM_LIST" ]]; then
  echo "[ERROR] BAM list not found: $BAM_LIST" >&2
  exit 1
fi

mapfile -t BAMS < <(grep -v '^#' "$BAM_LIST" | sed '/^$/d')

if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "[ERROR] BAM list is empty: $BAM_LIST" >&2
  exit 1
fi

echo "[INFO] Number of BAM files: ${#BAMS[@]}"

for bam in "${BAMS[@]}"; do
  if [[ ! -f "$bam" ]]; then
    echo "[ERROR] BAM not found: $bam" >&2
    exit 1
  fi
  if [[ ! -f "${bam}.bai" ]]; then
    echo "[INFO] BAM index missing; indexing: $bam"
    samtools index "$bam"
  fi
done

# -------------------------------
# 1) Joint inversion calling
# -------------------------------

echo "[INFO] Calling inversions with DELLY"
"$DELLY" call \
  -t INV \
  -g "$REF" \
  -o "$RAW_BCF" \
  "${BAMS[@]}"

# -------------------------------
# 2) Genotype inversion calls
# -------------------------------

echo "[INFO] Genotyping inversion calls across samples"
"$DELLY" call \
  -t INV \
  -g "$REF" \
  -v "$RAW_BCF" \
  -o "$GENOTYPED_BCF" \
  "${BAMS[@]}"

# -------------------------------
# 3) Convert BCF to VCF
# -------------------------------

echo "[INFO] Converting genotyped BCF to VCF"
bcftools view "$GENOTYPED_BCF" -Ov -o "$GENOTYPED_VCF"

echo "[INFO] Done. Outputs written to: $OUTDIR"
