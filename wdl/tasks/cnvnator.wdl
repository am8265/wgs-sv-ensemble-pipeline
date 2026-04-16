version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_cnvnator
#
# Read-depth CNV caller
# Detects: DEL, DUP (read-depth signal only)
#
# DESIGN NOTES:
#   - bin_size=100 is optimal for 30x WGS; use 500 for <15x coverage
#   - reference_dir must contain per-chromosome FASTAs named chr1.fa etc.
#   - cnvnator2VCF.pl handles the .txt → .vcf conversion
#   - ROOT file exposed as output — useful for rerunning with different
#     bin sizes without re-extracting read depth
# ─────────────────────────────────────────────────────────────────────────────

task run_cnvnator {
    input {
        SampleInfo    sample
        ReferenceData reference
        File          reference_dir
        Int           bin_size = 100
        Int           cpu      = 4
        Int           mem_gb   = 16
        String        docker   = "ghcr.io/am5153/sv-pipeline-cnvnator:0.4.1"
    }

    command <<<
        set -euo pipefail

        # ── Step 1: Extract read mapping from BAM into ROOT tree ──────────
        cnvnator \
            -root ~{sample.sample_id}.root \
            -tree ~{sample.bam} \
            -chrom $(seq 1 22 | sed 's/^/chr/') chrX chrY

        # ── Step 2: Generate read depth histogram ─────────────────────────
        cnvnator \
            -root ~{sample.sample_id}.root \
            -his ~{bin_size} \
            -d ~{reference_dir}

        # ── Step 3: Calculate mean and variance statistics ─────────────────
        cnvnator \
            -root ~{sample.sample_id}.root \
            -stat ~{bin_size}

        # ── Step 4: RD signal partitioning ────────────────────────────────
        cnvnator \
            -root ~{sample.sample_id}.root \
            -partition ~{bin_size}

        # ── Step 5: CNV calling ───────────────────────────────────────────
        cnvnator \
            -root ~{sample.sample_id}.root \
            -call ~{bin_size} \
            > ~{sample.sample_id}.cnvnator.txt

        # ── Step 6: Convert to VCF ────────────────────────────────────────
        cnvnator2VCF.pl \
            ~{sample.sample_id}.cnvnator.txt \
            ~{reference.genome_build} \
            > ~{sample.sample_id}.cnvnator.vcf

        bgzip ~{sample.sample_id}.cnvnator.vcf
        tabix -p vcf ~{sample.sample_id}.cnvnator.vcf.gz

        echo "CNVnator complete: $(bcftools stats ~{sample.sample_id}.cnvnator.vcf.gz \
            | grep 'number of records' | cut -f4) CNVs called (bin_size=~{bin_size}bp)"
    >>>

    output {
        SVCallerResult result = object {
            caller:    "cnvnator",
            vcf:       "~{sample.sample_id}.cnvnator.vcf.gz",
            vcf_index: "~{sample.sample_id}.cnvnator.vcf.gz.tbi"
        }
        # ROOT file exposed — can re-run calling with different bin sizes
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
        description: "Run CNVnator read-depth CNV caller (5-step pipeline)"
        author:      "Ayan Malakar"
    }

    parameter_meta {
        bin_size:      {description: "Bin size in bp. 100 for 30x WGS, 500 for <15x coverage."}
        reference_dir: {description: "Directory containing per-chromosome FASTAs (chr1.fa, chr2.fa ...)"}
    }
}
