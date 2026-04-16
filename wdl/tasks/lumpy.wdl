version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_lumpy
#
# Probabilistic SV caller — split-reads + discordant pairs
# Detects: DEL, DUP, INV, BND
#
# DESIGN NOTES:
#   - discordants + splitters extracted here and EXPOSED as outputs
#     so run_svtyper can consume them without re-extracting from BAM
#   - lumpyexpress wraps the full Lumpy workflow internally
#   - bgzip + tabix applied inline for consistency with other callers
# ─────────────────────────────────────────────────────────────────────────────

task run_lumpy {
    input {
        SampleInfo    sample
        ReferenceData reference
        Int           cpu    = 4
        Int           mem_gb = 16
        String        docker = "ghcr.io/am5153/sv-pipeline-lumpy:0.3.1"
    }

    command <<<
        set -euo pipefail

        # ── Extract discordant read pairs ─────────────────────────────────
        samtools view -b -F 1294 ~{sample.bam} \
            | samtools sort -o ~{sample.sample_id}.discordants.bam
        samtools index ~{sample.sample_id}.discordants.bam

        # ── Extract split reads ───────────────────────────────────────────
        samtools view -h ~{sample.bam} \
            | extractSplitReads_BwaMem -i stdin \
            | samtools view -Sb \
            | samtools sort -o ~{sample.sample_id}.splitters.bam
        samtools index ~{sample.sample_id}.splitters.bam

        # ── Run lumpyexpress ──────────────────────────────────────────────
        lumpyexpress \
            -B ~{sample.bam} \
            -S ~{sample.sample_id}.splitters.bam \
            -D ~{sample.sample_id}.discordants.bam \
            -o ~{sample.sample_id}.lumpy.vcf

        bgzip ~{sample.sample_id}.lumpy.vcf
        tabix -p vcf ~{sample.sample_id}.lumpy.vcf.gz

        echo "Lumpy complete: $(bcftools stats ~{sample.sample_id}.lumpy.vcf.gz \
            | grep 'number of records' | cut -f4) SVs called"
    >>>

    output {
        SVCallerResult result = object {
            caller:    "lumpy",
            vcf:       "~{sample.sample_id}.lumpy.vcf.gz",
            vcf_index: "~{sample.sample_id}.lumpy.vcf.gz.tbi"
        }
        # Expose BAMs for svtyper — avoids re-extraction downstream
        File discordants_bam     = "~{sample.sample_id}.discordants.bam"
        File discordants_bam_bai = "~{sample.sample_id}.discordants.bam.bai"
        File splitters_bam       = "~{sample.sample_id}.splitters.bam"
        File splitters_bam_bai   = "~{sample.sample_id}.splitters.bam.bai"
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
#
# Re-genotype SURVIVOR-merged SV calls using svtyper
#
# DESIGN NOTES:
#   - Uses Bayesian likelihood to assign GT, GQ, SQ per sample
#   - Discordant/splitter BAMs pre-extracted in run_lumpy — passed through
#     to avoid redundant BAM processing
#   - This produces the final genotyped output of the pipeline
# ─────────────────────────────────────────────────────────────────────────────

task run_svtyper {
    input {
        SampleInfo sample
        File       merged_vcf            # SURVIVOR-merged VCF (post-ensemble)
        File       discordants_bam
        File       discordants_bam_bai
        File       splitters_bam
        File       splitters_bam_bai
        Int        mem_gb = 8
        String     docker = "ghcr.io/am5153/sv-pipeline-lumpy:0.3.1"
    }

    command <<<
        set -euo pipefail
        

        svtyper \
            -i ~{merged_vcf} \
            -B ~{sample.bam} \
            -l ~{sample.bam}.json
            > ~{sample.sample_id}.merged.svtyper.vcf

        bgzip ~{sample.sample_id}.merged.svtyper.vcf
        tabix -p vcf ~{sample.sample_id}.merged.svtyper.vcf.gz

        echo "svtyper complete: $(bcftools stats ~{sample.sample_id}.merged.svtyper.vcf.gz \
            | grep 'number of records' | cut -f4) genotyped SVs"
    >>>

    output {
        File genotyped_vcf     = "~{sample.sample_id}.merged.svtyper.vcf.gz"
        File genotyped_vcf_tbi = "~{sample.sample_id}.merged.svtyper.vcf.gz.tbi"
    }

    runtime {
        docker:     docker
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 50 SSD"
        maxRetries: 1
    }

    meta {
        description: "Re-genotype SURVIVOR-merged SV calls with svtyper (Bayesian genotyping)"
        author:      "Ayan Malakar"
    }
}
