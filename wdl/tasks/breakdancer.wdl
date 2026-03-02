version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_breakdancer
# Discordant read pair SV caller
# Detects: DEL, INS, INV, ITX, CTX, Unknown
# ─────────────────────────────────────────────────────────────────────────────

task run_breakdancer {
    input {
        SampleInfo    sample
        ReferenceData reference
        Int           cpu    = 4
        Int           mem_gb = 8
        String        docker = "ayan-malakar/breakdancer:1.4.5"
    }

    command <<<
        set -euo pipefail

        # Generate config file from BAM
        bam2cfg.pl ~{sample.bam} \
            > ~{sample.sample_id}.breakdancer.cfg

        # Run BreakDancer
        breakdancer-max \
            ~{sample.sample_id}.breakdancer.cfg \
            > ~{sample.sample_id}.breakdancer.ctx

        # Convert to VCF
        breakdancer2vcf.py \
            --input  ~{sample.sample_id}.breakdancer.ctx \
            --output ~{sample.sample_id}.breakdancer.vcf \
            --sample ~{sample.sample_id} \
            --reference ~{reference.genome_build}

        bgzip ~{sample.sample_id}.breakdancer.vcf
        tabix -p vcf ~{sample.sample_id}.breakdancer.vcf.gz
    >>>

    output {
        SVCallerResult result = object {
            caller:    "breakdancer",
            vcf:       "~{sample.sample_id}.breakdancer.vcf.gz",
            vcf_index: "~{sample.sample_id}.breakdancer.vcf.gz.tbi"
        }
        File ctx_file = "~{sample.sample_id}.breakdancer.ctx"
    }

    runtime {
        docker:     docker
        cpu:        cpu
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 50 SSD"
        maxRetries: 1
    }

    meta {
        description: "Run BreakDancer discordant read pair SV caller"
        author:      "Ayan Malakar"
    }
}
