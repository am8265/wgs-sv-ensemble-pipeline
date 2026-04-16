version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_manta
#
# Split-read + discordant pair SV caller (Illumina)
# Detects: DEL, DUP, INV, INS, BND
#
# DESIGN NOTES:
#   - Uses diploidSV.vcf.gz (germline) not candidateSV (unfiltered)
#   - --jobs tied to cpu input so it scales with cloud instance size
#   - Output is compressed + tabix indexed inline — downstream tasks
#     can assume .vcf.gz + .tbi always exist together
# ─────────────────────────────────────────────────────────────────────────────

task run_manta {
    input {
        SampleInfo    sample
        ReferenceData reference
        Int           cpu    = 8
        Int           mem_gb = 16
        String        docker = "ghcr.io/am5153/sv-pipeline-manta:1.6.0"
    }

    String outdir = "manta_~{sample.sample_id}"

    command <<<
        set -euo pipefail

        configManta.py \
            --bam ~{sample.bam} \
            --referenceFasta ~{reference.fasta} \
            --runDir ~{outdir}

        ~{outdir}/runWorkflow.py \
            --mode local \
            --jobs ~{cpu} \
            --memGb ~{mem_gb}

        cp ~{outdir}/results/variants/diploidSV.vcf.gz \
            ~{sample.sample_id}.manta.vcf.gz
        cp ~{outdir}/results/variants/diploidSV.vcf.gz.tbi \
            ~{sample.sample_id}.manta.vcf.gz.tbi

        echo "Manta complete: $(bcftools stats ~{sample.sample_id}.manta.vcf.gz \
            | grep 'number of records' | cut -f4) SVs called"
    >>>

    output {
        SVCallerResult result = object {
            caller:    "manta",
            vcf:       "~{sample.sample_id}.manta.vcf.gz",
            vcf_index: "~{sample.sample_id}.manta.vcf.gz.tbi"
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
        description: "Run Manta germline SV caller (split-read + discordant pair)"
        author:      "Ayan Malakar"
    }

    parameter_meta {
        sample:    {description: "Sample BAM + index + ID"}
        reference: {description: "Reference FASTA + indexes + genome build"}
        cpu:       {description: "Number of CPUs — passed to Manta --jobs"}
        mem_gb:    {description: "Memory in GB — passed to Manta --memGb"}
    }
}
