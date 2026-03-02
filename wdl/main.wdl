version 1.0

# ─────────────────────────────────────────────────────────────────────────────
# WGS Structural Variant Ensemble Pipeline
# Author:  Ayan Malakar
# Version: 1.0.0
#
# Five-caller ensemble SV detection pipeline:
#   1. Manta      — split-read + discordant pair
#   2. Delly      — paired-end + split-read
#   3. Lumpy      — probabilistic, split-read + discordant pair
#   4. CNVnator   — read-depth
#   5. BreakDancer — discordant read pair
#
# QC: Duphold depth fold-change annotation + filtering
# Merge: SURVIVOR ensemble merging (minSUPP=2)
#
# Compatible with:
#   - AWS HealthOmics
#   - Terra (Broad Institute)
#   - Local Cromwell
#   - HPC (via Cromwell SLURM backend)
# ─────────────────────────────────────────────────────────────────────────────

import "structs.wdl"
import "tasks/manta.wdl"       as manta_task
import "tasks/delly.wdl"       as delly_task
import "tasks/lumpy.wdl"       as lumpy_task
import "tasks/cnvnator.wdl"    as cnvnator_task
import "tasks/breakdancer.wdl" as breakdancer_task
import "tasks/duphold.wdl"     as duphold_task
import "tasks/survivor.wdl"    as survivor_task

workflow WGSSVEnsemblePipeline {

    input {
        # ── Sample inputs ──────────────────────────────────────────────────
        String sample_id
        File   bam
        File   bai

        # ── Reference ─────────────────────────────────────────────────────
        File   reference_fasta
        File   reference_fai
        File   reference_dict
        String genome_build = "hg38"

        # ── CNVnator reference dir (per-chrom FASTAs) ──────────────────────
        File   cnvnator_reference_dir

        # ── Optional: Delly exclude regions ───────────────────────────────
        File?  delly_exclude_regions

        # ── Duphold filter thresholds ─────────────────────────────────────
        Float  del_dhbfc_threshold = 0.7
        Float  dup_dhffc_threshold = 1.25

        # ── SURVIVOR merge parameters ─────────────────────────────────────
        Int    survivor_max_dist    = 1000
        Int    survivor_min_support = 2
        Int    survivor_min_size    = 50

        # ── Compute resources ─────────────────────────────────────────────
        Int    manta_cpu      = 8
        Int    manta_mem_gb   = 16
        Int    delly_cpu      = 4
        Int    delly_mem_gb   = 16
        Int    lumpy_cpu      = 4
        Int    lumpy_mem_gb   = 16
        Int    cnvnator_cpu   = 4
        Int    cnvnator_mem_gb = 16
        Int    breakdancer_cpu = 4
        Int    breakdancer_mem_gb = 8
    }

    # ── Build shared structs ───────────────────────────────────────────────────
    SampleInfo sample = object {
        sample_id: sample_id,
        bam:       bam,
        bai:       bai
    }

    ReferenceData reference = object {
        fasta:        reference_fasta,
        fai:          reference_fai,
        dict:         reference_dict,
        genome_build: genome_build
    }

    # ════════════════════════════════════════════════════════════════════════
    # STAGE 1: Run all 5 SV callers in parallel
    # ════════════════════════════════════════════════════════════════════════

    call manta_task.run_manta {
        input:
            sample    = sample,
            reference = reference,
            cpu       = manta_cpu,
            mem_gb    = manta_mem_gb
    }

    call delly_task.run_delly {
        input:
            sample          = sample,
            reference       = reference,
            exclude_regions = delly_exclude_regions,
            cpu             = delly_cpu,
            mem_gb          = delly_mem_gb
    }

    call lumpy_task.run_lumpy {
        input:
            sample    = sample,
            reference = reference,
            cpu       = lumpy_cpu,
            mem_gb    = lumpy_mem_gb
    }

    # svtyper re-genotyping of Lumpy calls
    call lumpy_task.run_svtyper {
        input:
            sample          = sample,
            lumpy_vcf       = run_lumpy.result.vcf,
            discordants_bam = run_lumpy.discordants_bam,
            splitters_bam   = run_lumpy.splitters_bam
    }

    call cnvnator_task.run_cnvnator {
        input:
            sample        = sample,
            reference     = reference,
            reference_dir = cnvnator_reference_dir,
            cpu           = cnvnator_cpu,
            mem_gb        = cnvnator_mem_gb
    }

    call breakdancer_task.run_breakdancer {
        input:
            sample    = sample,
            reference = reference,
            cpu       = breakdancer_cpu,
            mem_gb    = breakdancer_mem_gb
    }

    # ════════════════════════════════════════════════════════════════════════
    # STAGE 2: Duphold QC annotation for each caller
    # ════════════════════════════════════════════════════════════════════════

    # Replace Lumpy raw result with svtyper re-genotyped VCF
    SVCallerResult lumpy_regenotyped = object {
        caller:    "lumpy",
        vcf:       run_svtyper.genotyped_vcf,
        vcf_index: run_svtyper.genotyped_vcf_tbi
    }

    Array[SVCallerResult] caller_results = [
        run_manta.result,
        run_delly.result,
        lumpy_regenotyped,
        run_cnvnator.result,
        run_breakdancer.result
    ]

    scatter (caller_result in caller_results) {
        call duphold_task.run_duphold {
            input:
                sample        = sample,
                caller_result = caller_result,
                reference     = reference
        }

        call duphold_task.filter_duphold {
            input:
                duphold_result       = run_duphold.result,
                del_dhbfc_threshold  = del_dhbfc_threshold,
                dup_dhffc_threshold  = dup_dhffc_threshold
        }
    }

    # ════════════════════════════════════════════════════════════════════════
    # STAGE 3: SURVIVOR ensemble merge (minSUPP = 2)
    # ════════════════════════════════════════════════════════════════════════

    call survivor_task.run_survivor {
        input:
            filtered_vcfs = filter_duphold.filtered_vcf,
            sample_id     = sample_id,
            max_dist      = survivor_max_dist,
            min_support   = survivor_min_support,
            min_size      = survivor_min_size
    }

    # ════════════════════════════════════════════════════════════════════════
    # OUTPUTS
    # ════════════════════════════════════════════════════════════════════════

    output {
        # ── Per-caller raw VCFs ───────────────────────────────────────────
        File manta_vcf       = run_manta.result.vcf
        File delly_vcf       = run_delly.result.vcf
        File lumpy_vcf       = run_svtyper.genotyped_vcf
        File cnvnator_vcf    = run_cnvnator.result.vcf
        File breakdancer_vcf = run_breakdancer.result.vcf

        # ── Duphold-filtered VCFs ─────────────────────────────────────────
        Array[File] filtered_vcfs  = filter_duphold.filtered_vcf
        Array[File] filter_stats   = filter_duphold.filter_stats

        # ── SURVIVOR merged ensemble VCF ──────────────────────────────────
        File merged_vcf      = run_survivor.merged_vcf
        File merged_vcf_tbi  = run_survivor.merged_vcf_tbi
        File merge_stats     = run_survivor.merge_stats
    }

    meta {
        author:      "Ayan Malakar"
        description: "WGS Structural Variant Ensemble Pipeline — 5-caller + Duphold QC + SURVIVOR merge"
        version:     "1.0.0"
    }
}
