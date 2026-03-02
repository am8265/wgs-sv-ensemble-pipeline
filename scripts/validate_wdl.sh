#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# scripts/validate_wdl.sh
# Validate all WDL files using womtool
# Usage: bash scripts/validate_wdl.sh
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

WOMTOOL_URL="https://github.com/broadinstitute/cromwell/releases/download/85/womtool-85.jar"
WOMTOOL_JAR="/tmp/womtool.jar"

# Download womtool if not present
if [ ! -f "$WOMTOOL_JAR" ]; then
    echo "Downloading womtool..."
    wget -q "$WOMTOOL_URL" -O "$WOMTOOL_JAR"
fi

echo "=== Validating WDL files ==="
PASS=0
FAIL=0

validate() {
    local wdl="$1"
    if java -jar "$WOMTOOL_JAR" validate "$wdl" > /dev/null 2>&1; then
        echo "  ✓ $wdl"
        ((PASS++))
    else
        echo "  ✗ $wdl"
        java -jar "$WOMTOOL_JAR" validate "$wdl"
        ((FAIL++))
    fi
}

validate wdl/main.wdl
for wdl in wdl/tasks/*.wdl; do
    validate "$wdl"
done

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "$FAIL" -eq 0 ] || exit 1
