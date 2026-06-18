#!/usr/bin/env python3
"""
04_breakpoint_coverage_dips.py
==============================

Check whether candidate inversions show coverage dips near their left and right
breakpoints in representative BAMs.

Input:
    - The text output produced by 03_samplot_inspection.py.
    - Representative BAM.

Output:
    - CSV table with Yes/No calls for coverage dips near each breakpoint.

The original analysis used a 1 kb breakpoint window and a mean coverage threshold
of 5 reads to mark a dip.
"""

from __future__ import annotations

import os
import pysam

# ==============================================================================
# User settings
# ==============================================================================

inversions_file = "path/to/inversion_summary.txt"
bam_file = "path/to/representative_sample_markdup.sorted.bam"
output_file = "path/to/coverage_dips_summary.csv"

coverage_threshold = 5
window_size = 1000

# ==============================================================================
# Main
# ==============================================================================


def mean_coverage(bam: pysam.AlignmentFile, chrom: str, start: int, end: int) -> float:
    """Return mean read depth across an interval."""
    if end <= start:
        return float("nan")
    cov = bam.count_coverage(chrom, start, end)
    return sum(map(sum, cov)) / (end - start)


def main() -> None:
    if not os.path.exists(inversions_file):
        raise FileNotFoundError(f"Inversion summary not found: {inversions_file}")
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM not found: {bam_file}")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    bam = pysam.AlignmentFile(bam_file, "rb")

    with open(output_file, "w", encoding="utf-8") as out:
        out.write(
            "Inversion,Chrom,Start,End,AvgCoverageStart,AvgCoverageEnd,"
            "CoverageDipStart,CoverageDipEnd\n"
        )

        with open(inversions_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line.startswith("Inversion"):
                    continue

                parts = line.split()
                inv_name = parts[1].rstrip(":")
                inv_parts = inv_name.split("_")

                chrom = "_".join(inv_parts[:-2])
                start = int(inv_parts[-2])
                end = int(inv_parts[-1])

                start_window_end = min(start + window_size, end)
                end_window_start = max(end - window_size, start)

                avg_cov_start = mean_coverage(bam, chrom, start, start_window_end)
                avg_cov_end = mean_coverage(bam, chrom, end_window_start, end)

                dip_start = "Yes" if avg_cov_start < coverage_threshold else "No"
                dip_end = "Yes" if avg_cov_end < coverage_threshold else "No"

                out.write(
                    f"{inv_name},{chrom},{start},{end},"
                    f"{avg_cov_start:.3f},{avg_cov_end:.3f},"
                    f"{dip_start},{dip_end}\n"
                )

    bam.close()
    print(f"Coverage dip summary written to: {output_file}")


if __name__ == "__main__":
    main()
