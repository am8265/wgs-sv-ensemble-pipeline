version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_manta
# Split-read + discordant pair SV caller (Illumina)
# Detects: DEL, DUP, INV, INS, BND
# ─────────────────────────────────────────────────────────────────────────────

task run_manta {
    input {
        SampleInfo    sample
        ReferenceData reference
        Int           cpu     = 8
        Int           mem_gb  = 16
        String        docker  = "ayan-malakar/manta:1.6.0"
    }

    String outdir = "manta_~{sample.sample_id}"

    command <<<
        set -euo pipefail

        # Configure Manta
        configManta.py \
            --bam ~{sample.bam} \
            --referenceFasta ~{reference.fasta} \
            --runDir ~{outdir}

        # Run Manta workflow
        ~{outdir}/runWorkflow.py \
            --mode local \
            --jobs ~{cpu} \
            --memGb ~{mem_gb}

        # Copy diploid SV VCF to working directory
        cp ~{outdir}/results/variants/diploidSV.vcf.gz \
            ~{sample.sample_id}.manta.vcf.gz
        cp ~{outdir}/results/variants/diploidSV.vcf.gz.tbi \
            ~{sample.sample_id}.manta.vcf.gz.tbi
    >>>

    output {
        SVCallerResult result = object {
            caller:    "manta",
            vcf:       "~{sample.sample_id}.manta.vcf.gz",
            vcf_index: "~{sample.sample_id}.manta.vcf.gz.tbi"
        }
    }

    runtime {
        docker: docker
        cpu:    cpu
        memory: "~{mem_gb} GB"
        disks:  "local-disk 100 SSD"
        # AWS HealthOmics / Terra compatible
        maxRetries: 1
    }

    meta {
        description: "Run Manta germline SV caller"
        author:      "Ayan Malakar"
    }
}
