#!/usr/bin/env bash
# =============================================================================
# merge_cohort.sh
# Merge per-sample svtyper-genotyped VCFs across a cohort using SURVIVOR
#
# Usage:
#   bash merge_cohort.sh \
#       --vcf_dir  <dir containing per-sample .vcf.gz files> \
#       --outdir   <output directory> \
#       --max_dist 1000 \
#       --min_sup  1 \
#       --min_size 50
#
# Output:
#   <outdir>/cohort.merged.vcf.gz
#   <outdir>/cohort.merged.vcf.gz.tbi
#   <outdir>/cohort_merge_stats.txt
#
# Notes:
#   - min_sup=1 here means any single sample having the call is enough to
#     include it in the cohort VCF. Rare SV filtering happens downstream
#     in filter_rare_svs.py.
#   - Input VCFs must be bgzipped and tabix-indexed (.vcf.gz + .tbi)
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
VCF_DIR=""
OUTDIR="cohort_merge_out"
MAX_DIST=1000
MIN_SUP=1
MIN_SIZE=50
SV_TYPE=1
SV_STRAND=1

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf_dir)  VCF_DIR="$2";  shift 2 ;;
        --outdir)   OUTDIR="$2";   shift 2 ;;
        --max_dist) MAX_DIST="$2"; shift 2 ;;
        --min_sup)  MIN_SUP="$2";  shift 2 ;;
        --min_size) MIN_SIZE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$VCF_DIR" ]]; then
    echo "ERROR: --vcf_dir is required"
    exit 1
fi

mkdir -p "$OUTDIR"

# ── Build VCF list ─────────────────────────────────────────────────────────────
VCF_LIST="${OUTDIR}/vcf_list.txt"
find "$VCF_DIR" -name "*.merged.svtyper.vcf.gz" | sort > "$VCF_LIST"

N_SAMPLES=$(wc -l < "$VCF_LIST")
echo "=== Cohort merge: ${N_SAMPLES} samples ==="
echo "Max breakpoint distance: ${MAX_DIST} bp"
echo "Min support:             ${MIN_SUP}"
echo "Min SV size:             ${MIN_SIZE} bp"
cat "$VCF_LIST"

if [[ "$N_SAMPLES" -eq 0 ]]; then
    echo "ERROR: No VCF files found in ${VCF_DIR}"
    exit 1
fi

# ── Run SURVIVOR cohort merge ─────────────────────────────────────────────────
SURVIVOR merge \
    "$VCF_LIST" \
    "$MAX_DIST" \
    "$MIN_SUP" \
    "$SV_TYPE" \
    "$SV_STRAND" \
    0 \
    "$MIN_SIZE" \
    "${OUTDIR}/cohort.merged.vcf"

bgzip "${OUTDIR}/cohort.merged.vcf"
tabix -p vcf "${OUTDIR}/cohort.merged.vcf.gz"

# ── Stats ─────────────────────────────────────────────────────────────────────
{
    echo "=== Cohort merge stats ==="
    echo "Samples merged: ${N_SAMPLES}"
    echo ""
    echo "--- SV counts by type ---"
    bcftools stats "${OUTDIR}/cohort.merged.vcf.gz" | grep "^SN"
    echo ""
    echo "--- SUPP distribution (how many samples carry each SV) ---"
    bcftools query -f '%INFO/SUPP\n' "${OUTDIR}/cohort.merged.vcf.gz" \
        | sort -n | uniq -c | sort -rn
} > "${OUTDIR}/cohort_merge_stats.txt"

cat "${OUTDIR}/cohort_merge_stats.txt"
echo ""
echo "Output: ${OUTDIR}/cohort.merged.vcf.gz"
