version 1.0

# ─────────────────────────────────────────────────────────────────────────────
# task: run_survivor
# Ensemble merging of SV calls from multiple callers using SURVIVOR
# Parameters:
#   max_dist:    max distance between breakpoints to merge (bp)
#   min_support: minimum number of callers supporting a call
#   sv_type:     1 = merge same SVTYPE only, 0 = merge all
#   sv_strand:   1 = consider strand, 0 = ignore
#   min_size:    minimum SV size to keep (bp)
# ─────────────────────────────────────────────────────────────────────────────

task run_survivor {
    input {
        Array[File] filtered_vcfs       # One per caller, post-duphold filter
        String      sample_id
        Int         max_dist    = 1000  # 1kb breakpoint window
        Int         min_support = 2     # At least 2 callers must agree
        Int         sv_type     = 1     # Merge same SVTYPE only
        Int         sv_strand   = 1     # Consider strand
        Int         sv_estimate = 0     # Do not estimate SV distance
        Int         min_size    = 50    # Minimum 50bp SV
        Int         mem_gb      = 8
        String      docker      = "ayan-malakar/survivor:1.0.7"
    }

    command <<<
        set -euo pipefail

        # Write VCF list file required by SURVIVOR
        printf '%s\n' ~{sep=' ' filtered_vcfs} > vcf_list.txt

        echo "=== SURVIVOR merge settings ===" 
        echo "Max distance:    ~{max_dist} bp"
        echo "Min support:     ~{min_support} callers"
        echo "Min SV size:     ~{min_size} bp"
        echo "Callers:"
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

        # Compress and index
        bgzip ~{sample_id}.survivor.merged.vcf
        tabix -p vcf ~{sample_id}.survivor.merged.vcf.gz

        # Summary stats
        echo "=== Merged SV counts by type ===" > merge_stats.txt
        bcftools stats ~{sample_id}.survivor.merged.vcf.gz \
            | grep "^SN" >> merge_stats.txt

        # Per-caller support summary
        echo "=== Per-caller support ===" >> merge_stats.txt
        bcftools query \
            -f '%SUPP_VEC\n' \
            ~{sample_id}.survivor.merged.vcf.gz \
            | sort | uniq -c | sort -rn >> merge_stats.txt
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
        description: "Merge SV calls from multiple callers using SURVIVOR ensemble approach"
        author:      "Ayan Malakar"
    }
}
