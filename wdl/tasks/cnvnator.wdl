version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_cnvnator
# Read-depth based CNV caller
# Detects: DEL, DUP (read-depth signal)
# Bin size 100bp recommended for WGS 30x
# ─────────────────────────────────────────────────────────────────────────────

task run_cnvnator {
    input {
        SampleInfo    sample
        ReferenceData reference
        File          reference_dir    # Directory containing per-chrom FASTAs
        Int           bin_size = 100   # Bin size in bp
        Int           cpu     = 4
        Int           mem_gb  = 16
        String        docker  = "ayan-malakar/cnvnator:0.4.1"
    }

    command <<<
        set -euo pipefail

        # Extract read mapping from BAM
        cnvnator \
            -root ~{sample.sample_id}.root \
            -tree ~{sample.bam} \
            -chrom $(seq 1 22 | sed 's/^/chr/') chrX chrY

        # Generate read depth histogram
        cnvnator \
            -root ~{sample.sample_id}.root \
            -his ~{bin_size} \
            -d ~{reference_dir}

        # Calculate statistics
        cnvnator \
            -root ~{sample.sample_id}.root \
            -stat ~{bin_size}

        # Partition into segments
        cnvnator \
            -root ~{sample.sample_id}.root \
            -partition ~{bin_size}

        # Call CNVs
        cnvnator \
            -root ~{sample.sample_id}.root \
            -call ~{bin_size} \
            > ~{sample.sample_id}.cnvnator.txt

        # Convert to VCF
        cnvnator2VCF.pl \
            ~{sample.sample_id}.cnvnator.txt \
            ~{reference.genome_build} \
            > ~{sample.sample_id}.cnvnator.vcf

        bgzip ~{sample.sample_id}.cnvnator.vcf
        tabix -p vcf ~{sample.sample_id}.cnvnator.vcf.gz
    >>>

    output {
        SVCallerResult result = object {
            caller:    "cnvnator",
            vcf:       "~{sample.sample_id}.cnvnator.vcf.gz",
            vcf_index: "~{sample.sample_id}.cnvnator.vcf.gz.tbi"
        }
        File root_file = "~{sample.sample_id}.root"
    }

    runtime {
        docker:     docker
        cpu:        cpu
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 100 SSD"
        maxRetries: 1
    }

    meta {
        description: "Run CNVnator read-depth CNV caller"
        author:      "Ayan Malakar"
    }
}
