version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_delly
#
# Paired-end + split-read SV caller
# Detects: DEL, DUP, INV, TRA, INS
#
# DESIGN NOTES:
#   - exclude_regions is optional (File?) — sensible default is the Delly
#     curated exclude list for hg38 (removes satellite/centromeric regions)
#   - BCF → VCF conversion done inline so all callers emit consistent .vcf.gz
#   - FILTER=PASS prefilter applied here — Delly's own quality filter
# ─────────────────────────────────────────────────────────────────────────────

task run_delly {
    input {
        SampleInfo    sample
        ReferenceData reference
        File?         exclude_regions
        Int           cpu    = 4
        Int           mem_gb = 16
        String        docker = "ghcr.io/am5153/sv-pipeline-delly:1.1.6"
    }

    command <<<
        set -euo pipefail

        delly call \
            --genome ~{reference.fasta} \
            ~{"--exclude " + exclude_regions} \
            --outfile ~{sample.sample_id}.delly.bcf \
            ~{sample.bam}

        # Keep only PASS calls + convert BCF → VCF
        bcftools view \
            --apply-filters PASS \
            --output-type z \
            --output ~{sample.sample_id}.delly.vcf.gz \
            ~{sample.sample_id}.delly.bcf

        bcftools index --tbi ~{sample.sample_id}.delly.vcf.gz

        echo "Delly complete: $(bcftools stats ~{sample.sample_id}.delly.vcf.gz \
            | grep 'number of records' | cut -f4) SVs called (PASS only)"
    >>>

    output {
        SVCallerResult result = object {
            caller:    "delly",
            vcf:       "~{sample.sample_id}.delly.vcf.gz",
            vcf_index: "~{sample.sample_id}.delly.vcf.gz.tbi"
        }
    }

    runtime {
        docker:     docker
        cpu:        cpu
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 100 SSD"
        maxRetries: 1
    }

    meta {
        description: "Run Delly germline SV caller (paired-end + split-read)"
        author:      "Ayan Malakar"
    }

    parameter_meta {
        exclude_regions: {description: "Optional BED of regions to exclude (e.g. centromeres, satellites). Recommended: Delly hg38 exclude list."}
    }
}
