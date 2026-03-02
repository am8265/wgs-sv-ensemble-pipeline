version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_duphold
# Depth fold-change annotation + QC for SV calls
# Adds DHBFC (DEL filter) and DHFFC (DUP filter) annotations
# Recommended thresholds:
#   DEL: DHBFC < 0.7
#   DUP: DHFFC > 1.25
# ─────────────────────────────────────────────────────────────────────────────

task run_duphold {
    input {
        SampleInfo    sample
        SVCallerResult caller_result
        ReferenceData  reference
        Int            mem_gb = 8
        String         docker = "ayan-malakar/duphold:0.2.3"
    }

    String out_prefix = "~{sample.sample_id}.~{caller_result.caller}.duphold"

    command <<<
        set -euo pipefail

        duphold \
            --vcf ~{caller_result.vcf} \
            --bam ~{sample.bam} \
            --fasta ~{reference.fasta} \
            --threads 4 \
            --output ~{out_prefix}.vcf.gz

        bcftools index --tbi ~{out_prefix}.vcf.gz
    >>>

    output {
        DupholdResult result = object {
            caller:              caller_result.caller,
            annotated_vcf:       "~{out_prefix}.vcf.gz",
            annotated_vcf_index: "~{out_prefix}.vcf.gz.tbi"
        }
    }

    runtime {
        docker:     docker
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 50 SSD"
        maxRetries: 1
    }

    meta {
        description: "Annotate SV calls with Duphold depth fold-change for QC filtering"
        author:      "Ayan Malakar"
    }
}

# ─────────────────────────────────────────────────────────────────────────────
# task: filter_duphold
# Apply DHBFC/DHFFC thresholds to remove false positives
# ─────────────────────────────────────────────────────────────────────────────

task filter_duphold {
    input {
        DupholdResult duphold_result
        Float         del_dhbfc_threshold = 0.7    # DEL: keep DHBFC < threshold
        Float         dup_dhffc_threshold = 1.25   # DUP: keep DHFFC > threshold
        Int           mem_gb = 4
        String        docker = "ayan-malakar/duphold:0.2.3"
    }

    String out_prefix = "~{duphold_result.caller}.filtered"

    command <<<
        set -euo pipefail

        bcftools filter \
            --include '(SVTYPE="DEL" && FMT/DHBFC[0] < ~{del_dhbfc_threshold}) || \
                       (SVTYPE="DUP" && FMT/DHFFC[0] > ~{dup_dhffc_threshold}) || \
                       (SVTYPE!="DEL" && SVTYPE!="DUP")' \
            --output-type z \
            --output ~{out_prefix}.vcf.gz \
            ~{duphold_result.annotated_vcf}

        bcftools index --tbi ~{out_prefix}.vcf.gz

        # Log filter stats
        echo "=== Duphold filter stats: ~{duphold_result.caller} ===" \
            >> filter_stats.txt
        bcftools stats ~{out_prefix}.vcf.gz \
            | grep "^SN" >> filter_stats.txt
    >>>

    output {
        File filtered_vcf     = "~{out_prefix}.vcf.gz"
        File filtered_vcf_tbi = "~{out_prefix}.vcf.gz.tbi"
        File filter_stats     = "filter_stats.txt"
    }

    runtime {
        docker:     docker
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 50 SSD"
        maxRetries: 1
    }

    meta {
        description: "Filter SV calls using Duphold DHBFC/DHFFC thresholds"
        author:      "Ayan Malakar"
    }
}
