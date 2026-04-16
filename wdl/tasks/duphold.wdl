version 1.0

import "../structs.wdl"

# ─────────────────────────────────────────────────────────────────────────────
# task: run_duphold
#
# Depth fold-change annotation for SV calls
# Adds per-sample FORMAT fields:
#   DHBFC — depth fold-change for the variant depth relative to bins in the genome with similar GC-content
#   DHFFC — depth fold-change vs flanking regions 
#

#   - Runs on each caller's VCF independently (scatter in main.wdl)
#   - Requires BAM + FASTA for depth extraction
#   - Does NOT filter — annotation only here, filtering is separate task
#     (separation of concerns: annotate first, filter second)
# ─────────────────────────────────────────────────────────────────────────────

task run_duphold {
    input {
        SampleInfo     sample
        SVCallerResult caller_result
        ReferenceData  reference
        Int            threads = 4
        Int            mem_gb  = 8
        String         docker  = "ghcr.io/am5153/sv-pipeline-duphold:0.2.3"
    }

    String out_prefix = "~{sample.sample_id}.~{caller_result.caller}.duphold"

    command <<<
        set -euo pipefail

        duphold \
            --vcf ~{caller_result.vcf} \
            --bam ~{sample.bam} \
            --fasta ~{reference.fasta} \
            --threads ~{threads} \
            --output ~{out_prefix}.vcf.gz

        bcftools index --tbi ~{out_prefix}.vcf.gz

        # Quick annotation summary
        echo "=== Duphold annotation complete: ~{caller_result.caller} ==="
        bcftools query \
            -f '%SVTYPE\t%INFO/SVLEN\t[%DHBFC]\t[%DHFFC]\n' \
            ~{out_prefix}.vcf.gz \
            | head -20
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
        cpu:        threads
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 50 SSD"
        maxRetries: 1
    }

    meta {
        description: "Annotate SV calls with Duphold depth fold-change (DHBFC/DHFFC)"
        author:      "Ayan Malakar"
    }
}

# ─────────────────────────────────────────────────────────────────────────────
# task: filter_duphold
#
# Apply depth fold-change thresholds to remove false positives
#
# FILTER LOGIC:
#
#   DUP: FILTER != LowQual AND (DHBFC > 1.3 OR DHFFC > 1.3)
#
#   DEL: FILTER != LowQual
#        AND (DHFFC < 0.7 OR DHBFC < 0.7)   
#        AND (DHGT[2] - DHGT[1] >= 0        ← hom_alt >= het count in region
#             OR DHGT[1] <= 5               ← very few het calls
#             OR GT == 1/1)                 ← homozygous deletion
#
#   Rationale for DEL genotype check:
#     Low DHBFC alone is insufficient for het DELs. In segmental duplication
#     regions, mapping artifacts produce low coverage without real deletions.
#     Checking that hom_alt >= het adds a genotype-level
#     confirmation that the depth signal is real.
# ─────────────────────────────────────────────────────────────────────────────

task filter_duphold {
    input {
        DupholdResult duphold_result
        # DEL thresholds — depth loss in flanks OR body
        Float         del_dhffc_threshold  = 0.7    # DHFFC < 0.7 (flank depth loss)
        Float         del_dhbfc_threshold  = 0.7    # DHBFC < 0.7 (body depth loss)
        # DUP thresholds — depth gain in flanks OR body
        Float         dup_dhbfc_threshold  = 1.3    # DHBFC > 1.3 (body depth gain)
        Float         dup_dhffc_threshold  = 1.3    # DHFFC > 1.3 (flank depth gain)
        Int           mem_gb = 4
        String        docker = "ghcr.io/am5153/sv-pipeline-duphold:0.2.3"
    }

    String out_prefix = "~{duphold_result.caller}.duphold_filtered"

    command <<<
        set -euo pipefail

        # ── Primary depth filter ──────────────────────────────────────────
        # DEL: DHBFC < threshold AND (homozygous OR has sufficient read support)
        # DUP: DHFFC > threshold
        # INV/BND/INS: pass through (no depth filter applicable)

        # ── DUP filter ────────────────────────────────────────────────────
        # Keep duplications that are:
        #   - Not LowQual AND
        #   - DHBFC > 1.3 OR DHFFC > 1.3 (depth gain in similar GC or flanks)
        bcftools filter \
            --include 'SVTYPE="DUP" && FILTER!="LowQual" && (FMT/DHBFC[0] > ~{dup_dhbfc_threshold} || FMT/DHFFC[0] > ~{dup_dhffc_threshold})' \
            --output-type z \
            ~{duphold_result.annotated_vcf} \
        | bcftools view --output-type z --output ~{out_prefix}.dup_filtered.vcf.gz
        bcftools index --tbi ~{out_prefix}.dup_filtered.vcf.gz

        # ── DEL filter ────────────────────────────────────────────────────
        # Keep deletions that are:
        #   - Not LowQual AND
        #   - DHFFC < 0.7 OR DHBFC < 0.7 (depth loss in flanks or body) AND
        #   - (DHGT[2] - DHGT[1] >= 0 OR DHGT[1] <= 5 OR GT == 1/1)
        #     i.e. hom_alt - het >= 0 (more hom than het calls in region)
        #          OR very few het calls (<=5) in region
        #          OR homozygous deletion call
        bcftools filter \
            --include 'SVTYPE="DEL" && FILTER!="LowQual" && (FMT/DHFFC[0] < ~{del_dhffc_threshold} || FMT/DHBFC[0] < ~{del_dhbfc_threshold}) && (FMT/DHGT[0:2] - FMT/DHGT[0:1] >= 0 || FMT/DHGT[0:1] <= 5 || FMT/GT[0]="1/1")' \
            --output-type z \
            ~{duphold_result.annotated_vcf} \
        | bcftools view --output-type z --output ~{out_prefix}.del_filtered.vcf.gz
        bcftools index --tbi ~{out_prefix}.del_filtered.vcf.gz

        # ── Pass through INV / BND / INS ──────────────────────────────────
        bcftools filter \
            --include 'SVTYPE!="DEL" && SVTYPE!="DUP"' \
            --output-type z \
            ~{duphold_result.annotated_vcf} \
        | bcftools view --output-type z --output ~{out_prefix}.other_filtered.vcf.gz
        bcftools index --tbi ~{out_prefix}.other_filtered.vcf.gz

        # ── Concatenate all filtered VCFs ─────────────────────────────────
        bcftools concat \
            --allow-overlaps \
            --output-type z \
            --output ~{out_prefix}.vcf.gz \
            ~{out_prefix}.dup_filtered.vcf.gz \
            ~{out_prefix}.del_filtered.vcf.gz \
            ~{out_prefix}.other_filtered.vcf.gz
        bcftools index --tbi ~{out_prefix}.vcf.gz

        bcftools index --tbi ~{out_prefix}.vcf.gz

        # ── Filter stats ──────────────────────────────────────────────────
        echo "=== Filter stats: ~{duphold_result.caller} ===" > filter_stats.txt

        echo "--- After Duphold depth filter ---" >> filter_stats.txt
        bcftools stats ~{out_prefix}.vcf.gz \
            | grep "^SN" >> filter_stats.txt

        cat filter_stats.txt
    >>>

    output {
        File filtered_vcf         = "~{out_prefix}.vcf.gz"
        File filtered_vcf_tbi     = "~{out_prefix}.vcf.gz.tbi"
        File filter_stats         = "filter_stats.txt"
    }

    runtime {
        docker:     docker
        memory:     "~{mem_gb} GB"
        disks:      "local-disk 50 SSD"
        maxRetries: 1
    }

    meta {
        description: "Filter SV calls using Duphold DHBFC/DHFFC thresholds with het DEL artifact detection"
        author:      "Ayan Malakar"
    }

    parameter_meta {
        del_dhffc_threshold: {description: "DEL: keep if DHFFC < threshold (depth loss in flanking region). Default 0.7."}
        del_dhbfc_threshold: {description: "DEL: keep if DHBFC < threshold (depth loss in GC bins). Default 0.7."}
        dup_dhbfc_threshold: {description: "DUP: keep if DHBFC > threshold (depth gain in GC bins). Default 1.3."}
        dup_dhffc_threshold: {description: "DUP: keep if DHFFC > threshold (depth gain in flanking region). Default 1.3."}
    }
}


