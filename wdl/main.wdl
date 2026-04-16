version 1.0

# ─────────────────────────────────────────────────────────────────────────────
# WGS Structural Variant Ensemble Pipeline — main.wdl
# adoption from Parliament2 ensemble approach.
# Author:  Ayan Malakar (github.com/am5153)
# Version: 1.0.0
#
# PIPELINE STAGES:
#   1. Run 5 SV callers in parallel (Manta, Delly, Lumpy, CNVnator, BreakDancer)
#   2. QC: Duphold depth fold-change + B-Allele annotation + filtering
#   4. SURVIVOR ensemble merge (minSUPP=2 of 5 callers)
#   5. svtyper re-genotyping on SURVIVOR-merged VCF
#
# SUPP_VEC caller order (preserved throughout):
#   Position 1 — Manta
#   Position 2 — Delly
#   Position 3 — Lumpy
#   Position 4 — CNVnator
#   Position 5 — BreakDancer
#
# COMPATIBLE WITH:
#   - Local Cromwell (configs/cromwell_local.conf)
#   - AWS HealthOmics (native WDL submission)
#   - Terra / Broad Institute
#   - HPC SGE/SLURM (via Cromwell SGE/SLURM backend)
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
        # ── Sample ─────────────────────────────────────────────────────────
        String sample_id
        File   bam
        File   bai

        # ── Reference ──────────────────────────────────────────────────────
        File   reference_fasta
        File   reference_fai
        File   reference_dict
        String genome_build = "hg38"

        # ── CNVnator: per-chrom FASTA directory (passed as path string) ───
        String cnvnator_reference_dir
        Int    cnvnator_bin_size = 100

        # ── Delly: optional exclude regions ────────────────────────────────
        File?  delly_exclude_regions

        # ── Duphold filter thresholds ──────────────────────────────────────
        Float  del_dhbfc_threshold = 0.7
        Float  dup_dhffc_threshold = 1.25

        # ── SURVIVOR merge parameters ──────────────────────────────────────
        Int    survivor_max_dist    = 1000
        Int    survivor_min_support = 2
        Int    survivor_min_size    = 50

        # ── Compute ────────────────────────────────────────────────────────
        Int    manta_cpu       = 8
        Int    manta_mem_gb    = 16
        Int    delly_cpu       = 4
        Int    delly_mem_gb    = 16
        Int    lumpy_cpu       = 4
        Int    lumpy_mem_gb    = 16
        Int    cnvnator_cpu    = 4
        Int    cnvnator_mem_gb = 16
        Int    breakdancer_cpu = 4
        Int    breakdancer_mem_gb = 8
    }

    # ── Build shared structs ───────────────────────────────────────────────
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

    # ══════════════════════════════════════════════════════════════════════
    # STAGE 1: Run all 5 callers in parallel
    # All tasks start simultaneously — Cromwell handles parallelism
    # ══════════════════════════════════════════════════════════════════════

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

    call cnvnator_task.run_cnvnator {
        input:
            sample        = sample,
            reference     = reference,
            reference_dir = cnvnator_reference_dir,
            bin_size      = cnvnator_bin_size,
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

    # ══════════════════════════════════════════════════════════════════════
    # STAGE 2: Build caller results array
    # IMPORTANT: Order here determines SUPP_VEC bit positions in SURVIVOR multi caller merge
    #   [0] Manta       → bit 1
    #   [1] Delly       → bit 2
    #   [2] Lumpy       → bit 3  (raw — svtyper runs AFTER merge)
    #   [3] CNVnator    → bit 4
    #   [4] BreakDancer → bit 5
    # ══════════════════════════════════════════════════════════════════════

    Array[SVCallerResult] caller_results = [
        run_manta.result,
        run_delly.result,
        run_lumpy.result,
        run_cnvnator.result,
        run_breakdancer.result
    ]

    # ══════════════════════════════════════════════════════════════════════
    # STAGE 3: Duphold annotation + filtering (scatter — parallel per caller)
    # ══════════════════════════════════════════════════════════════════════

    scatter (caller_result in caller_results) {

        call duphold_task.run_duphold {
            input:
                sample        = sample,
                caller_result = caller_result,
                reference     = reference
        }

        call duphold_task.filter_duphold {
            input:
                duphold_result      = run_duphold.result,
                del_dhbfc_threshold = del_dhbfc_threshold,
                dup_dhffc_threshold = dup_dhffc_threshold,
        }
    }

    # ══════════════════════════════════════════════════════════════════════
    # STAGE 4: SURVIVOR ensemble merge
    # Inputs are gathered from scatter — Array[File] of filtered VCFs
    # ══════════════════════════════════════════════════════════════════════

    call survivor_task.run_survivor {
        input:
            filtered_vcfs = filter_duphold.filtered_vcf,
            sample_id     = sample_id,
            max_dist      = survivor_max_dist,
            min_support   = survivor_min_support,
            min_size      = survivor_min_size
    }

    # ══════════════════════════════════════════════════════════════════════
    # STAGE 5: svtyper re-genotyping on SURVIVOR-merged VCF
    # Rationale: genotype the merged call set using all original BAM evidence
    # Input: merged multi-caller VCF (breakpoints resolved by SURVIVOR)
    # Output: genotyped VCF with GT, GQ, SQ per sample
    # ══════════════════════════════════════════════════════════════════════

    call lumpy_task.run_svtyper {
        input:
            sample               = sample,
            merged_vcf           = run_survivor.merged_vcf,
            discordants_bam      = run_lumpy.discordants_bam,
            discordants_bam_bai  = run_lumpy.discordants_bam_bai,
            splitters_bam        = run_lumpy.splitters_bam,
            splitters_bam_bai    = run_lumpy.splitters_bam_bai
    }

    # ══════════════════════════════════════════════════════════════════════
    # OUTPUTS
    # ══════════════════════════════════════════════════════════════════════

    output {
        # ── Per-caller raw VCFs ──────────────────────────────────────────
        File manta_vcf       = run_manta.result.vcf
        File delly_vcf       = run_delly.result.vcf
        File lumpy_vcf       = run_lumpy.result.vcf
        File cnvnator_vcf    = run_cnvnator.result.vcf
        File breakdancer_vcf = run_breakdancer.result.vcf

        # ── Duphold-filtered VCFs (one per caller) ───────────────────────
        Array[File] filtered_vcfs    = filter_duphold.filtered_vcf
        Array[File] per_caller_stats = filter_duphold.filter_stats

        # ── SURVIVOR ensemble merged VCF ─────────────────────────────────
        File merged_vcf     = run_survivor.merged_vcf
        File merged_vcf_tbi = run_survivor.merged_vcf_tbi
        File merge_stats    = run_survivor.merge_stats

        # ── svtyper re-genotyped final VCF ────────────────────────────────
        File genotyped_vcf     = run_svtyper.genotyped_vcf
        File genotyped_vcf_tbi = run_svtyper.genotyped_vcf_tbi

    }

    meta {
        author:      "Ayan Malakar"
        description: "WGS SV Ensemble Pipeline — 5-caller + Duphold QC + SURVIVOR merge"
        version:     "1.0.0"
    }
}
