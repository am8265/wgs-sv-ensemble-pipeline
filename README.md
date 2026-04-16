# WGS Structural Variant Ensemble Pipeline

[![CI — Validate & Lint](https://github.com/am8265/wgs-sv-ensemble-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/8265/wgs-sv-ensemble-pipeline/actions/workflows/ci.yml)
[![Docker — Build & Push](https://github.com/am8265/wgs-sv-ensemble-pipeline/actions/workflows/docker.yml/badge.svg)](https://github.com/am8265/wgs-sv-ensemble-pipeline/actions/workflows/docker.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![WDL 1.0](https://img.shields.io/badge/WDL-1.0-blue.svg)](https://openwdl.org/)

An almost production-grade, cloud-ready WDL pipeline for germline structural variant (SV) detection from whole genome sequencing (WGS) data. Five complementary SV callers are run in parallel, quality-filtered using Duphold depth fold-change thresholds, merged with SURVIVOR, and genotyped with svtyper. Cohort-level rare SV filtering and AnnotSV annotation are handled by a separate set of scripts designed to run after all per-sample processing is complete.

> This pipeline is based on SV calling work originally developed at Columbia University Medical Center (2018–2024). The original codebase was institutional. This repository reimplements the core logic as an open-source project.

---

## Pipeline Overview

The pipeline is split into two tiers:

### Tier 1 — Per-sample (`wdl/main.wdl`)

Run once per sample via Cromwell, AWS HealthOmics, or Terra:

```
WGS BAM (hg38 / hg19)
  │
  ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 1 — 5 SV CALLERS (all run in parallel)                │
│  Manta · Delly · Lumpy · CNVnator · BreakDancer              │
│                                                              │
│  Detection methods:                                          │
│    Split-read + discordant pair  — Manta, Delly, Lumpy       │
│    Read depth                    — CNVnator                  │
│    Discordant read pairs         — BreakDancer               │
└──────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 2 — DUPHOLD QC (scatter — one job per caller)         │
│                                                              │
│  Annotation:  DHBFC, DHFFC depth fold-change fields added    │
│                                                              │
│  DUP filter:  FILTER != LowQual                              │
│               AND (DHBFC > 1.3 OR DHFFC > 1.3)              │
│                                                              │
│  DEL filter:  FILTER != LowQual                              │
│               AND (DHFFC < 0.7 OR DHBFC < 0.7)              │
│               AND (DHGT[2]-DHGT[1] >= 0                      │
│                    OR DHGT[1] <= 5                           │
│                    OR GT == 1/1)                             │
└──────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 3 — SURVIVOR ENSEMBLE MERGE                           │
│  Merges 5 per-caller filtered VCFs into one call set         │
│  max_dist = 1000 bp · min_support = 2 · min_size = 50 bp     │
│  SUPP_VEC encodes per-caller support                         │
│  (bit order: Manta / Delly / Lumpy / CNVnator / BreakDancer) │
└──────────────────────────────────────────────────────────────┘
  │
  ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 4 — SVTYPER GENOTYPING                                │
│  Bayesian re-genotyping on SURVIVOR-merged call set          │
│  Uses pre-extracted discordant + splitter BAMs from Lumpy    │
│  Output: per-sample VCF with GT / GQ / SQ per variant        │
└──────────────────────────────────────────────────────────────┘
  │
  ▼
Per-sample filtered + genotyped SV VCF
```

---

### Tier 2 — Cohort-level (`scripts/`)

Run once after all per-sample VCFs have been generated:

```
All per-sample VCFs  (one per sample from Tier 1)
  │
  ▼  scripts/merge_cohort.sh
┌──────────────────────────────────────────────────────────────┐
│  SURVIVOR cohort merge                                       │
│  Produces multi-sample VCF across the full cohort            │
└──────────────────────────────────────────────────────────────┘
  │
  ▼  scripts/filter_rare_svs.py  (default: --max_af 0.002)
┌──────────────────────────────────────────────────────────────┐
│  Cohort allele frequency filter                              │
│  Retains SVs present in <= 2 / 1000 samples (configurable)  │
│  Output: multi-sample rare SV VCF                            │
└──────────────────────────────────────────────────────────────┘
  │
  ▼  scripts/run_annotsv.sh
┌──────────────────────────────────────────────────────────────┐
│  AnnotSV annotation on rare SV VCF                           │
│  Databases: ClinVar · ClinGen · gnomAD (pLI/LOEUF)          │
│             OMIM · DGV · GeneHancer                          │
│  Output: annotated TSV + pathogenic subset (class 4 + 5)     │
└──────────────────────────────────────────────────────────────┘
```

> **Why run AnnotSV only on the rare SV VCF?**
> AnnotSV runtime scales linearly with input size. Running it on the cohort-level rare SV VCF — after filtering out common SVs — reduces the annotation input from tens of thousands of calls to a few hundred, making the step fast and practical at scale.

---

## SV Callers

| Caller | Detection Method | SV Types | SUPP_VEC Bit |
|--------|-----------------|----------|-------------|
| [Manta](https://github.com/Illumina/manta) | Split-read + discordant pair | DEL, DUP, INV, INS, BND | 1 |
| [Delly](https://github.com/dellytools/delly) | Paired-end + split-read | DEL, DUP, INV, TRA, INS | 2 |
| [Lumpy](https://github.com/arq5x/lumpy-sv) | Probabilistic SR + DP | DEL, DUP, INV, BND | 3 |
| [CNVnator](https://github.com/abyzovlab/CNVnator) | Read depth | DEL, DUP | 4 |
| [BreakDancer](https://github.com/genome/breakdancer) | Discordant read pairs | DEL, INS, INV, CTX | 5 |

---

## Duphold Filter Logic

Production-validated thresholds applied per-caller before ensemble merging.

**Duplications** — keep if:
```
FILTER != LowQual AND (DHBFC > 1.3 OR DHFFC > 1.3)
```

**Deletions** — keep if:
```
FILTER != LowQual
AND (DHFFC < 0.7 OR DHBFC < 0.7)
AND (DHGT[2] - DHGT[1] >= 0   <- hom_alt >= het count in region
     OR DHGT[1] <= 5           <- very few het calls (artifact signal)
     OR GT == 1/1)             <- homozygous deletion
```

**INV / BND / INS** — pass through without depth filtering.

The genotype distribution check (`DHGT`) addresses a key limitation: in segmental duplication and low-complexity regions, low DHBFC alone is insufficient to confirm a real deletion. Confirming that homozygous calls dominate over heterozygous calls provides an independent genotype-level signal that the depth drop reflects a true variant.

---

## Quick Start

### Prerequisites

- [Cromwell](https://github.com/broadinstitute/cromwell/releases) >= 85
- Docker
- Java 11+

### Tier 1 — Run per-sample locally

```bash
git clone https://github.com/am8265/wgs-sv-ensemble-pipeline.git
cd wgs-sv-ensemble-pipeline

# Configure inputs
cp inputs/inputs.json inputs/my_sample.json
# Edit my_sample.json — update BAM path, reference paths, thresholds

# Run
java -Dconfig.file=configs/cromwell_local.conf \
     -jar cromwell.jar run \
     wdl/main.wdl \
     --inputs inputs/my_sample.json
```

### Tier 1 — Run on AWS HealthOmics

```bash
zip -r sv_pipeline.zip wdl/

aws omics create-workflow \
    --name sv-ensemble-pipeline \
    --definition-zip fileb://sv_pipeline.zip \
    --main wdl/main.wdl \
    --parameter-template file://inputs/inputs.json
```

### Tier 1 — Run on [Terra](https://terra.bio/)

Import this repo via **Workspace → Import Workflow → GitHub**, then submit via the Terra UI or `fiss` CLI.

### Tier 2 — Cohort processing

```bash
# 1. Merge all per-sample VCFs
bash scripts/merge_cohort.sh \
    --vcf_dir  /path/to/per_sample_vcfs/ \
    --outdir   cohort_out/

# 2. Filter for rare SVs
python3 scripts/filter_rare_svs.py \
    --vcf    cohort_out/cohort.merged.vcf.gz \
    --out    cohort_out/cohort.rare_svs.vcf.gz \
    --max_af 0.002 \
    --stats  cohort_out/filter_stats.tsv

# 3. Annotate rare SVs
bash scripts/run_annotsv.sh \
    --vcf    cohort_out/cohort.rare_svs.vcf.gz \
    --outdir cohort_out/annotsv/
```

---

## Key Parameters

### Tier 1 — `wdl/main.wdl`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `survivor_min_support` | 2 | Min callers required to retain a call |
| `survivor_max_dist` | 1000 bp | Max breakpoint distance to merge calls |
| `survivor_min_size` | 50 bp | Min SV size retained |
| `del_dhffc_threshold` | 0.7 | DEL: flanking region depth loss cutoff |
| `del_dhbfc_threshold` | 0.7 | DEL: SV body depth loss cutoff |
| `dup_dhbfc_threshold` | 1.3 | DUP: SV body depth gain cutoff |
| `dup_dhffc_threshold` | 1.3 | DUP: flanking region depth gain cutoff |
| `cnvnator_bin_size` | 100 bp | Read depth histogram bin size (100 for 30x WGS) |

### Tier 2 — `scripts/`

| Script | Parameter | Default | Description |
|--------|-----------|---------|-------------|
| `filter_rare_svs.py` | `--max_af` | 0.002 | Max cohort allele frequency |
| `filter_rare_svs.py` | `--min_size` | 50 bp | Min SV size |
| `merge_cohort.sh` | `--max_dist` | 1000 bp | Breakpoint merge distance |
| `merge_cohort.sh` | `--min_sup` | 1 | Min samples for cohort merge |

---

## Repository Structure

```
wgs-sv-ensemble-pipeline/
├── .github/
│   └── workflows/
│       ├── ci.yml          # WDL validation (womtool) + Dockerfile linting (hadolint)
│       └── docker.yml      # Build & push Docker images to GHCR on merge to main
├── wdl/
│   ├── main.wdl            # Top-level workflow — 4 stages, per-sample
│   ├── structs.wdl         # Typed data structures (SampleInfo, SVCallerResult, ...)
│   └── tasks/
│       ├── manta.wdl
│       ├── delly.wdl
│       ├── lumpy.wdl       # Includes svtyper task
│       ├── cnvnator.wdl
│       ├── breakdancer.wdl
│       ├── duphold.wdl     # Annotation task + production filter task
│       ├── survivor.wdl
│       └── annotsv.wdl     # Standalone task — used by Tier 2 scripts
├── docker/
│   ├── manta/Dockerfile
│   ├── delly/Dockerfile
│   ├── lumpy/Dockerfile    # Includes svtyper
│   ├── cnvnator/Dockerfile
│   ├── breakdancer/Dockerfile
│   ├── duphold/Dockerfile
│   └── survivor/Dockerfile
├── scripts/
│   ├── merge_cohort.sh     # SURVIVOR cohort merge
│   ├── filter_rare_svs.py  # Cohort AF filter — keeps SVs in <= N samples
│   └── run_annotsv.sh      # AnnotSV annotation on rare SV VCF
├── inputs/
│   └── inputs.json         # Example inputs for Tier 1
├── configs/
│   ├── cromwell_local.conf
│   └── cromwell_aws.conf
└── docs/
    └── architecture.md
```

---

## Cloud Compatibility

| Platform | Backend | Notes |
|----------|---------|-------|
| Local | Docker via Cromwell | `configs/cromwell_local.conf` |
| AWS HealthOmics | Native WDL | See Quick Start above |
| Terra / Broad | Google Cloud | Import via Terra UI |
| HPC (SGE/SLURM) | Cromwell SGE/SLURM backend | Extend `cromwell_local.conf` |

---

## Citation

If you use this pipeline, please cite the underlying tools:

- Manta: [Chen et al., Bioinformatics 2016](https://doi.org/10.1093/bioinformatics/btv710)
- Delly: [Rausch et al., Bioinformatics 2012](https://doi.org/10.1093/bioinformatics/bts378)
- Lumpy: [Layer et al., Genome Biology 2014](https://doi.org/10.1186/gb-2014-15-6-r84)
- CNVnator: [Abyzov et al., Genome Research 2011](https://doi.org/10.1101/gr.114876.110)
- BreakDancer: [Chen et al., Nature Methods 2009](https://doi.org/10.1038/nmeth.1363)
- Duphold: [Pedersen et al., GigaScience 2019](https://doi.org/10.1093/gigascience/giz040)
- SURVIVOR: [Jeffares et al., Nature Communications 2017](https://doi.org/10.1038/ncomms14061)
- svtyper: [Chiang et al., Nature Communications 2017](https://doi.org/10.1038/ncomms14061)
- AnnotSV: [Geoffroy et al., Bioinformatics 2018](https://doi.org/10.1093/bioinformatics/bty304)

- Parliament2 *(inspiration)*: [Zarate et al., GigaScience 2020](https://doi.org/10.1093/gigascience/giaa145)

---

## Author

**Ayan Malakar** · [github.com/am8265](https://github.com/am8265/)

---

## License

MIT — see [LICENSE](LICENSE)