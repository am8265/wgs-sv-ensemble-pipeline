version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_lumpy
# Probabilistic SV caller using split-reads + discordant pairs
# Detects: DEL, DUP, INV, BND
# Note: svtyper used downstream for re-genotyping
# ─────────────────────────────────────────────────────────────────────────────

task run_lumpy {
    input {
        SampleInfo    sample
        ReferenceData reference
        Int           cpu    = 4
        Int           mem_gb = 16
        String        docker = "ayan-malakar/lumpy:0.3.1"
    }

    command <<<
        set -euo pipefail

        # Extract discordant reads
        samtools view -b -F 1294 ~{sample.bam} \
            | samtools sort -o ~{sample.sample_id}.discordants.bam
        samtools index ~{sample.sample_id}.discordants.bam

        # Extract split reads
        samtools view -h ~{sample.bam} \
            | extractSplitReads_BwaMem -i stdin \
            | samtools view -Sb \
            | samtools sort -o ~{sample.sample_id}.splitters.bam
        samtools index ~{sample.sample_id}.splitters.bam

        # Run lumpyexpress
        lumpyexpress \
            -B ~{sample.bam} \
            -S ~{sample.sample_id}.splitters.bam \
            -D ~{sample.sample_id}.discordants.bam \
            -o ~{sample.sample_id}.lumpy.vcf

        # Compress and index
        bgzip ~{sample.sample_id}.lumpy.vcf
        tabix -p vcf ~{sample.sample_id}.lumpy.vcf.gz
    >>>

    output {
        SVCallerResult result = object {
            caller:    "lumpy",
            vcf:       "~{sample.sample_id}.lumpy.vcf.gz",
            vcf_index: "~{sample.sample_id}.lumpy.vcf.gz.tbi"
        }
        # Expose BAMs for svtyper re-genotyping downstream
        File discordants_bam = "~{sample.sample_id}.discordants.bam"
        File splitters_bam   = "~{sample.sample_id}.splitters.bam"
    }

    runtime {
        docker:     docker
        cpu:        cpu
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 150 SSD"
        maxRetries: 1
    }

    meta {
        description: "Run Lumpy SV caller with split-read and discordant pair extraction"
        author:      "Ayan Malakar"
    }
}

# ─────────────────────────────────────────────────────────────────────────────
# task: run_svtyper
# Re-genotype Lumpy calls using svtyper
# ─────────────────────────────────────────────────────────────────────────────

task run_svtyper {
    input {
        SampleInfo sample
        File       lumpy_vcf
        File       discordants_bam
        File       splitters_bam
        Int        mem_gb = 8
        String     docker = "ayan-malakar/lumpy:0.3.1"
    }

    command <<<
        set -euo pipefail

        svtyper \
            --input_vcf ~{lumpy_vcf} \
            --bam ~{sample.bam} \
            --splitter_bam ~{splitters_bam} \
            --discordant_bam ~{discordants_bam} \
            --output_vcf ~{sample.sample_id}.lumpy.svtyper.vcf

        bgzip ~{sample.sample_id}.lumpy.svtyper.vcf
        tabix -p vcf ~{sample.sample_id}.lumpy.svtyper.vcf.gz
    >>>

    output {
        File genotyped_vcf     = "~{sample.sample_id}.lumpy.svtyper.vcf.gz"
        File genotyped_vcf_tbi = "~{sample.sample_id}.lumpy.svtyper.vcf.gz.tbi"
    }

    runtime {
        docker:     docker
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 50 SSD"
        maxRetries: 1
    }

    meta {
        description: "Re-genotype Lumpy SV calls with svtyper"
        author:      "Ayan Malakar"
    }
}
