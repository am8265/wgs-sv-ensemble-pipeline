version 1.0

# ─────────────────────────────────────────────────────────────────────────────
# structs.wdl
# Shared data structures for the WGS SV Ensemble Pipeline
# ─────────────────────────────────────────────────────────────────────────────

struct SampleInfo {
    String sample_id
    File   bam
    File   bai
}

struct ReferenceData {
    File   fasta
    File   fai
    File   dict
    String genome_build   # hg38 or hg19
}

struct SVCallerResult {
    String caller
    File   vcf
    File   vcf_index
}

struct DupholdResult {
    String caller
    File   annotated_vcf
    File   annotated_vcf_index
}
