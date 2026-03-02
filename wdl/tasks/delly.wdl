version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_delly
# Paired-end + split-read SV caller
# Detects: DEL, DUP, INV, TRA, INS
# ─────────────────────────────────────────────────────────────────────────────

task run_delly {
    input {
        SampleInfo    sample
        ReferenceData reference
        File?         exclude_regions   # Optional mappability exclude BED
        Int           cpu    = 4
        Int           mem_gb = 16
        String        docker = "ayan-malakar/delly:1.1.6"
    }

    command <<<
        set -euo pipefail

        delly call \
            --genome ~{reference.fasta} \
            ~{"--exclude " + exclude_regions} \
            --outfile ~{sample.sample_id}.delly.bcf \
            ~{sample.bam}

        # Convert BCF → VCF and index
        bcftools view \
            --output-type z \
            --output ~{sample.sample_id}.delly.vcf.gz \
            ~{sample.sample_id}.delly.bcf

        bcftools index --tbi ~{sample.sample_id}.delly.vcf.gz
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
        description: "Run Delly germline SV caller"
        author:      "Ayan Malakar"
    }
}
