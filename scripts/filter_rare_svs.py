#!/usr/bin/env python3
"""
filter_rare_svs.py
Filter cohort-merged SV VCF to retain only rare SVs.

Usage:
    python3 filter_rare_svs.py \
        --vcf    cohort.merged.vcf.gz \
        --out    cohort.rare_svs.vcf.gz \
        --max_af 0.002 \
        --min_size 50

Arguments:
    --vcf      Input cohort VCF (SURVIVOR-merged, bgzipped)
    --out      Output VCF (bgzipped + tabix indexed)
    --max_af   Maximum allele frequency to keep (default: 0.002 = 2/1000)
    --min_size Minimum SV size in bp (default: 50)
    --stats    Optional: output stats TSV file

Author: Ayan Malakar
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Filter cohort SV VCF to retain rare SVs"
    )
    p.add_argument("--vcf",      required=True,  help="Input cohort VCF (.vcf.gz)")
    p.add_argument("--out",      required=True,  help="Output rare SV VCF (.vcf.gz)")
    p.add_argument("--max_af",   type=float, default=0.002,
                   help="Max allele frequency to retain (default: 0.002 = 2/1000)")
    p.add_argument("--min_size", type=int,   default=50,
                   help="Min SV size in bp (default: 50)")
    p.add_argument("--stats",    default=None,
                   help="Optional: output stats TSV file")
    return p.parse_args()


def count_samples(vcf: str) -> int:
    """Count number of samples in VCF header."""
    result = subprocess.run(
        ["bcftools", "query", "-l", vcf],
        capture_output=True, text=True, check=True
    )
    return len(result.stdout.strip().splitlines())


def get_sv_counts(vcf: str) -> dict:
    """Get SV counts by type from bcftools stats."""
    result = subprocess.run(
        ["bcftools", "stats", vcf],
        capture_output=True, text=True, check=True
    )
    counts = {}
    for line in result.stdout.splitlines():
        if line.startswith("SN"):
            parts = line.split("\t")
            if len(parts) >= 4 and "number of" in parts[2]:
                counts[parts[2].strip()] = parts[3].strip()
    return counts


def filter_rare_svs(args: argparse.Namespace) -> None:
    vcf_path = Path(args.vcf)
    if not vcf_path.exists():
        print(f"ERROR: Input VCF not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)

    # Count total samples to compute AF threshold
    n_samples = count_samples(args.vcf)
    max_ac = max(1, int(args.max_af * n_samples))  # max allele count

    print(f"=== Rare SV filter ===")
    print(f"Input:         {args.vcf}")
    print(f"Samples:       {n_samples}")
    print(f"Max AF:        {args.max_af} ({max_ac} samples)")
    print(f"Min SV size:   {args.min_size} bp")
    print(f"Output:        {args.out}")

    # Build bcftools filter expression
    # SUPP = number of samples carrying the SV (from SURVIVOR INFO field)
    # SVLEN may be negative for deletions — use abs via size filter
    filter_expr = (
        f'INFO/SUPP <= {max_ac} && '
        f'(INFO/SVLEN >= {args.min_size} || INFO/SVLEN <= -{args.min_size})'
    )

    out_path = Path(args.out)
    tmp_vcf = str(out_path).replace(".vcf.gz", ".tmp.vcf.gz")

    # Apply filter
    subprocess.run([
        "bcftools", "filter",
        "--include", filter_expr,
        "--output-type", "z",
        "--output", tmp_vcf,
        args.vcf
    ], check=True)

    # Index
    subprocess.run(["bcftools", "index", "--tbi", tmp_vcf], check=True)

    # Rename to final output
    os.rename(tmp_vcf, args.out)
    tbi_tmp = tmp_vcf + ".tbi"
    if os.path.exists(tbi_tmp):
        os.rename(tbi_tmp, args.out + ".tbi")
    else:
        subprocess.run(["bcftools", "index", "--tbi", args.out], check=True)

    # Stats
    before_counts = get_sv_counts(args.vcf)
    after_counts  = get_sv_counts(args.out)

    before_total = before_counts.get("number of records:", "?")
    after_total  = after_counts.get("number of records:", "?")

    print(f"\n=== Filter results ===")
    print(f"Before: {before_total} SVs")
    print(f"After:  {after_total} SVs (AF <= {args.max_af})")

    if args.stats:
        with open(args.stats, "w") as f:
            f.write("metric\tvalue\n")
            f.write(f"input_vcf\t{args.vcf}\n")
            f.write(f"n_samples\t{n_samples}\n")
            f.write(f"max_af\t{args.max_af}\n")
            f.write(f"max_ac\t{max_ac}\n")
            f.write(f"min_size_bp\t{args.min_size}\n")
            f.write(f"svs_before\t{before_total}\n")
            f.write(f"svs_after\t{after_total}\n")
        print(f"Stats written to: {args.stats}")


if __name__ == "__main__":
    filter_rare_svs(parse_args())
