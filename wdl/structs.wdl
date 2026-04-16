version 1.0

# ─────────────────────────────────────────────────────────────────────────────
# structs.wdl
# Shared data structures for the WGS SV Ensemble Pipeline
#
# DESIGN NOTE:
#   Using WDL structs instead of loose File/String pairs enforces a contract
#   between tasks — callers can't accidentally pass a Delly VCF to a Manta
#   task. 

struct SampleInfo {
    String sample_id
    File   bam
    File   bai
}

struct ReferenceData {
    File   fasta
    File   fai
    File   dict
    String genome_build   # "hg38" or "hg19"
}

struct SVCallerResult {
    String caller         # "manta" | "delly" | "lumpy" | "cnvnator" | "breakdancer"
    File   vcf
    File   vcf_index
}

struct DupholdResult {
    String caller
    File   annotated_vcf
    File   annotated_vcf_index
}
