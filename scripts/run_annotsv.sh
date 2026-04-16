#!/usr/bin/env bash
# =============================================================================
# run_annotsv.sh
# Run AnnotSV on cohort rare SV VCF to produce annotated TSV
#
# Usage:
#   bash run_annotsv.sh \
#       --vcf          cohort.rare_svs.vcf.gz \
#       --outdir       annotsv_out \
#       --genome_build GRCh38 \
#       --hpo          "HP:0000707,HP:0001250"   # optional
#
# Output:
#   <outdir>/cohort.rare_svs.annotsv.tsv   — full AnnotSV annotation
#   <outdir>/cohort.rare_svs.pathogenic.tsv — filtered: class 4 + 5 only
#
# Notes:
#   - Runs with AnnotSV defaults — no custom thresholds needed
#   - HPO terms optional; if provided enables Exomiser phenotype scoring
#   - Requires AnnotSV to be installed and in PATH, or run via Docker:
#       docker run --rm \
#           -v $(pwd):/data \
#           ghcr.io/am5153/sv-pipeline-annotsv:3.4 \
#           bash /data/scripts/run_annotsv.sh --vcf /data/cohort.rare_svs.vcf.gz ...
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT_VCF=""
OUTDIR="annotsv_out"
GENOME_BUILD="GRCh38"
HPO_TERMS=""
MIN_SIZE=50

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)          INPUT_VCF="$2";    shift 2 ;;
        --outdir)       OUTDIR="$2";       shift 2 ;;
        --genome_build) GENOME_BUILD="$2"; shift 2 ;;
        --hpo)          HPO_TERMS="$2";    shift 2 ;;
        --min_size)     MIN_SIZE="$2";     shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT_VCF" ]]; then
    echo "ERROR: --vcf is required"
    exit 1
fi

SAMPLE_NAME=$(basename "$INPUT_VCF" .vcf.gz)
mkdir -p "$OUTDIR"

echo "=== AnnotSV annotation ==="
echo "Input:        ${INPUT_VCF}"
echo "Output dir:   ${OUTDIR}"
echo "Genome build: ${GENOME_BUILD}"
echo "HPO terms:    ${HPO_TERMS:-none}"

# ── Build AnnotSV command ──────────────────────────────────────────────────────
ANNOTSV_CMD=(
    AnnotSV
    -SVinputFile  "$INPUT_VCF"
    -outputDir    "$OUTDIR"
    -outputFile   "${SAMPLE_NAME}.annotsv"
    -genomeBuild  "$GENOME_BUILD"
    -SVminSize    "$MIN_SIZE"
    -overlap      70
    -reciprocal   yes
)

# Append HPO if provided
if [[ -n "$HPO_TERMS" ]]; then
    ANNOTSV_CMD+=(-hpo "$HPO_TERMS")
fi

# ── Run AnnotSV ───────────────────────────────────────────────────────────────
"${ANNOTSV_CMD[@]}"

FULL_TSV="${OUTDIR}/${SAMPLE_NAME}.annotsv.tsv"

if [[ ! -f "$FULL_TSV" ]]; then
    echo "ERROR: AnnotSV did not produce expected output: ${FULL_TSV}"
    exit 1
fi

echo ""
echo "=== AnnotSV complete ==="
echo "Full annotation: ${FULL_TSV}"

# ── Extract pathogenic / likely pathogenic (class 4 + 5) ─────────────────────
PATHOGENIC_TSV="${OUTDIR}/${SAMPLE_NAME}.pathogenic.tsv"

# AnnotSV ranking column is the last column; 5=Pathogenic, 4=Likely Pathogenic
head -1 "$FULL_TSV" > "$PATHOGENIC_TSV"
awk -F'\t' 'NR>1 && $2=="full" && ($NF=="5" || $NF=="4")' \
    "$FULL_TSV" >> "$PATHOGENIC_TSV"

N_PATH=$(awk 'NR>1' "$PATHOGENIC_TSV" | wc -l)
echo "Pathogenic / Likely Pathogenic SVs: ${N_PATH}"
echo "Pathogenic TSV: ${PATHOGENIC_TSV}"

# ── Summary stats ─────────────────────────────────────────────────────────────
echo ""
echo "=== AnnotSV class distribution ==="
awk -F'\t' 'NR>1 && $2=="full" {print $NF}' "$FULL_TSV" \
    | sort | uniq -c | sort -rn \
    | awk '{printf "  Class %s: %s SVs\n", $2, $1}'
