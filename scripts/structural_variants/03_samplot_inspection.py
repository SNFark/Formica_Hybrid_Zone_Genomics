#!/usr/bin/env python3
"""
03_samplot_inspection.py
========================

Generate Samplot images for filtered DELLY inversion calls and count read-pair
orientations in one representative BAM.

Original analysis:
    - DELLY calls came from 42 high-coverage BAMs.
    - Candidate inversions were inspected with Samplot using representative
      high-coverage BAMs.
    - Scaffold 3 was intentionally omitted because it corresponds to the known
      social-supergene region.

This script is generalised for GitHub use. Replace the paths in USER SETTINGS.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
import pysam

# ==============================================================================
# User settings
# ==============================================================================

bam_file = "path/to/representative_sample_markdup.sorted.bam"
fasta_file = "path/to/reference_genome.fasta"
vcf_file = "path/to/inversions_filtered_Q100_PE14_SR14_100kb.vcf"
output_dir = "path/to/samplot_output"
samplot_exe = "samplot"

exclude_chromosomes = {"FsiP_PB_v5_scf3"}

# ==============================================================================
# Helper functions
# ==============================================================================


def parse_info(info: str) -> dict[str, str]:
    """Parse the INFO column of a VCF row into a dictionary."""
    parsed: dict[str, str] = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
    return parsed


def count_orientations(bam_path: str, chrom: str, start: int, end: int) -> tuple[int, int, int]:
    """Count FR/RF/Other read-pair orientations in a candidate inversion interval."""
    fr = rf = other = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.mate_is_unmapped:
                other += 1
                continue

            if read.is_proper_pair:
                if read.is_reverse and not read.mate_is_reverse:
                    rf += 1
                elif not read.is_reverse and read.mate_is_reverse:
                    fr += 1
                else:
                    other += 1
            else:
                other += 1

    return fr, rf, other


def main() -> None:
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    for path in [bam_file, fasta_file, vcf_file]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Required input not found: {path}")

    with open(vcf_file, "r", encoding="utf-8") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom = cols[0]
            if chrom in exclude_chromosomes:
                continue

            start = int(cols[1])
            info = parse_info(cols[7])
            svtype = info.get("SVTYPE", "")
            end = int(info.get("END", start))

            if svtype != "INV":
                continue

            label = f"{chrom}_{start}_{end}"
            png_file = os.path.join(output_dir, f"{label}.png")

            print(f"Plotting inversion: {label}")

            subprocess.run(
                [
                    samplot_exe,
                    "plot",
                    "-b",
                    bam_file,
                    "-r",
                    fasta_file,
                    "-c",
                    chrom,
                    "-s",
                    str(start),
                    "-e",
                    str(end),
                    "-t",
                    "INV",
                    "-o",
                    png_file,
                ],
                check=True,
            )

            fr, rf, other = count_orientations(bam_file, chrom, start, end)
            print(f"Inversion {label}: FR={fr}, RF={rf}, Other={other}")


if __name__ == "__main__":
    main()
