version 1.0

# ─────────────────────────────────────────────────────────────────────────────
# task: run_survivor
#
# Ensemble merging of per-caller filtered VCFs using SURVIVOR
#
# SURVIVOR merge parameters:
#   max_dist    — max breakpoint distance to merge (bp)
#   min_support — min callers required to support a call (ensemble threshold)
#   sv_type     — 1 = only merge same SVTYPE; 0 = merge across types
#   sv_strand   — 1 = consider strand; 0 = ignore strand
#   sv_estimate — 1 = estimate SV distance; 0 = use exact coordinates
#   min_size    — minimum SV size to retain (bp)
#
# DESIGN NOTES:
#   - min_support=2 balances sensitivity/specificity for 5-caller ensemble
#     (validated against GIAB HG002 — see docs/architecture.md)
#   - SUPP_VEC output field encodes which callers support each call
#     e.g. "10110" = Manta + Lumpy + BreakDancer but not Delly/CNVnator
#   - merge_stats logs SUPP_VEC distribution — useful QC metric
# ─────────────────────────────────────────────────────────────────────────────

task run_survivor {
    input {
        Array[File] filtered_vcfs
        String      sample_id
        Int         max_dist    = 1000
        Int         min_support = 2
        Int         sv_type     = 1
        Int         sv_strand   = 1
        Int         sv_estimate = 0
        Int         min_size    = 50
        Int         mem_gb      = 8
        String      docker      = "ghcr.io/am5153/sv-pipeline-survivor:1.0.7"
    }

    command <<<
        set -euo pipefail

        # Write VCF list file required by SURVIVOR
        printf '%s\n' ~{sep=' ' filtered_vcfs} > vcf_list.txt

        echo "=== SURVIVOR merge settings ==="
        echo "Callers: $(wc -l < vcf_list.txt)"
        echo "Max breakpoint distance: ~{max_dist} bp"
        echo "Min caller support:      ~{min_support}"
        echo "Min SV size:             ~{min_size} bp"
        cat vcf_list.txt

        SURVIVOR merge \
            vcf_list.txt \
            ~{max_dist} \
            ~{min_support} \
            ~{sv_type} \
            ~{sv_strand} \
            ~{sv_estimate} \
            ~{min_size} \
            ~{sample_id}.survivor.merged.vcf

        bgzip ~{sample_id}.survivor.merged.vcf
        tabix -p vcf ~{sample_id}.survivor.merged.vcf.gz

        # ── QC stats ──────────────────────────────────────────────────────
        {
            echo "=== Merged SV counts by SVTYPE ==="
            bcftools stats ~{sample_id}.survivor.merged.vcf.gz \
                | grep "^SN"

            echo ""
            echo "=== SUPP_VEC distribution (which callers support each call) ==="
            echo "Format: count  SUPP_VEC (1=supported, 0=not supported, order: Manta/Delly/Lumpy/CNVnator/BreakDancer)"
            bcftools query \
                -f '%INFO/SUPP_VEC\n' \
                ~{sample_id}.survivor.merged.vcf.gz \
                | sort | uniq -c | sort -rn

            echo ""
            echo "=== Calls supported by all 5 callers ==="
            bcftools view \
                --include 'INFO/SUPP=5' \
                ~{sample_id}.survivor.merged.vcf.gz \
                | bcftools stats | grep "^SN"
        } > merge_stats.txt

        cat merge_stats.txt
    >>>

    output {
        File merged_vcf     = "~{sample_id}.survivor.merged.vcf.gz"
        File merged_vcf_tbi = "~{sample_id}.survivor.merged.vcf.gz.tbi"
        File merge_stats    = "merge_stats.txt"
        File vcf_list       = "vcf_list.txt"
    }

    runtime {
        docker:     docker
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 50 SSD"
        maxRetries: 1
    }

    meta {
        description: "Ensemble merge of per-caller SV VCFs using SURVIVOR"
        author:      "Ayan Malakar"
    }

    parameter_meta {
        filtered_vcfs: {description: "Array of per-caller Duphold-filtered VCFs. Order determines SUPP_VEC bit position."}
        min_support:   {description: "Min callers required. 2 of 5 balances sensitivity/specificity."}
        max_dist:      {description: "Max breakpoint distance in bp to merge calls. 1000bp is standard for WGS."}
    }
}
