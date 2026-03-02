# WGS Structural Variant Ensemble Pipeline

[![CI — Validate & Lint](https://github.com/ayan-malakar/wgs-sv-ensemble-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/ayan-malakar/wgs-sv-ensemble-pipeline/actions/workflows/ci.yml)
[![Docker — Build & Push](https://github.com/ayan-malakar/wgs-sv-ensemble-pipeline/actions/workflows/docker.yml/badge.svg)](https://github.com/ayan-malakar/wgs-sv-ensemble-pipeline/actions/workflows/docker.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![WDL 1.0](https://img.shields.io/badge/WDL-1.0-blue.svg)](https://openwdl.org/)

A production-grade, cloud-ready ensemble pipeline for germline structural variant (SV) detection from whole genome sequencing (WGS) data. Combines five complementary SV callers with Duphold QC annotation and SURVIVOR ensemble merging.

---

## Pipeline Overview

```
WGS BAM/CRAM
     │
     ▼
┌─────────────────────────────────────────────────────┐
│              5 SV CALLERS (parallel)                │
│  Manta  │  Delly  │  Lumpy  │ CNVnator │BreakDancer │
└─────────────────────────────────────────────────────┘
     │
     ▼ (svtyper re-genotyping for Lumpy)
┌─────────────────────────────────────────────────────┐
│           DUPHOLD QC ANNOTATION (scatter)           │
│     DHBFC (DEL < 0.7)  │  DHFFC (DUP > 1.25)      │
└─────────────────────────────────────────────────────┘
     │
     ▼
┌─────────────────────────────────────────────────────┐
│        SURVIVOR ENSEMBLE MERGE (minSUPP=2)          │
│   max_dist=1000bp  │  min_size=50bp                 │
└─────────────────────────────────────────────────────┘
     │
     ▼
Filtered Multi-Caller VCF
```

## SV Callers

| Caller | Method | SV Types | Reference |
|--------|--------|----------|-----------|
| [Manta](https://github.com/Illumina/manta) | Split-read + discordant pair | DEL, DUP, INV, INS, BND | Chen et al. 2016 |
| [Delly](https://github.com/dellytools/delly) | Paired-end + split-read | DEL, DUP, INV, TRA, INS | Rausch et al. 2012 |
| [Lumpy](https://github.com/arq5x/lumpy-sv) | Probabilistic, SR + DP | DEL, DUP, INV, BND | Layer et al. 2014 |
| [CNVnator](https://github.com/abyzovlab/CNVnator) | Read depth | DEL, DUP | Abyzov et al. 2011 |
| [BreakDancer](https://github.com/genome/breakdancer) | Discordant read pairs | DEL, INS, INV, CTX | Chen et al. 2009 |

## QC & Merging

- **[Duphold](https://github.com/brentp/duphold)** — depth fold-change annotation (DHBFC/DHFFC)
- **[SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)** — ensemble merging with breakpoint proximity window
- **[svtyper](https://github.com/hall-lab/svtyper)** — Lumpy call re-genotyping

---

## Quick Start

### Prerequisites

- [Cromwell](https://github.com/broadinstitute/cromwell/releases) ≥ 85
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/)
- Java 11+

### 1. Clone the repo

```bash
git clone https://github.com/ayan-malakar/wgs-sv-ensemble-pipeline.git
cd wgs-sv-ensemble-pipeline
```

### 2. Configure inputs

```bash
cp inputs/inputs.json inputs/my_sample.json
# Edit my_sample.json with your BAM, reference paths
```

### 3. Run locally

```bash
java -Dconfig.file=configs/cromwell_local.conf \
     -jar cromwell.jar run \
     wdl/main.wdl \
     --inputs inputs/my_sample.json
```

### 4. Run on Terra (Broad Institute)

1. Import this repo into Terra via **Workspace → Import Workflow → GitHub**
2. Upload your inputs JSON
3. Submit via Terra UI or `fiss` CLI

### 5. Run on AWS HealthOmics

```bash
# Zip the WDL directory
zip -r sv_pipeline.zip wdl/

# Create HealthOmics workflow
aws omics create-workflow \
    --name sv-ensemble-pipeline \
    --definition-zip fileb://sv_pipeline.zip \
    --main wdl/main.wdl \
    --parameter-template file://inputs/inputs.json
```

---

## Repository Structure

```
wgs-sv-ensemble-pipeline/
├── .github/workflows/
│   ├── ci.yml          # WDL validation + Dockerfile linting
│   └── docker.yml      # Build & push Docker images to GHCR
├── wdl/
│   ├── main.wdl        # Top-level workflow
│   ├── structs.wdl     # Shared data structures
│   └── tasks/
│       ├── manta.wdl
│       ├── delly.wdl
│       ├── lumpy.wdl   # Includes svtyper re-genotyping
│       ├── cnvnator.wdl
│       ├── breakdancer.wdl
│       ├── duphold.wdl # Annotation + filtering tasks
│       └── survivor.wdl
├── docker/             # One Dockerfile per tool
├── inputs/             # Example input JSONs
├── configs/
│   ├── cromwell_local.conf
│   └── cromwell_aws.conf
├── scripts/
│   └── validate_wdl.sh
├── tests/
│   └── test_inputs.json
└── docs/
    └── architecture.md
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `survivor_min_support` | 2 | Min callers supporting a call |
| `survivor_max_dist` | 1000 | Max breakpoint distance (bp) |
| `survivor_min_size` | 50 | Min SV size to retain (bp) |
| `del_dhbfc_threshold` | 0.7 | Duphold DEL filter threshold |
| `dup_dhffc_threshold` | 1.25 | Duphold DUP filter threshold |

---

## Cloud Compatibility

| Platform | Backend | Config |
|----------|---------|--------|
| Local | Docker | `configs/cromwell_local.conf` |
| AWS HealthOmics | Native WDL | See Quick Start |
| Terra (Broad) | Google Cloud | Import via Terra UI |
| HPC (SLURM) | Cromwell SLURM | Extend `cromwell_local.conf` |

---

## CI/CD

Every push triggers:
- **WDL validation** via `womtool`
- **Dockerfile linting** via `hadolint`
- **JSON validation** of all input files
- **Docker image builds** pushed to GitHub Container Registry on merge to `main`

---

## Citation

If you use this pipeline, please cite the individual tools:

- Manta: [Chen et al., Bioinformatics 2016](https://doi.org/10.1093/bioinformatics/btv710)
- Delly: [Rausch et al., Bioinformatics 2012](https://doi.org/10.1093/bioinformatics/bts378)
- Lumpy: [Layer et al., Genome Biology 2014](https://doi.org/10.1186/gb-2014-15-6-r84)
- CNVnator: [Abyzov et al., Genome Research 2011](https://doi.org/10.1101/gr.114876.110)
- BreakDancer: [Chen et al., Nature Methods 2009](https://doi.org/10.1038/nmeth.1363)
- Duphold: [Pedersen et al., GigaScience 2019](https://doi.org/10.1093/gigascience/giz040)
- SURVIVOR: [Jeffares et al., Nature Communications 2017](https://doi.org/10.1038/ncomms14061)

---

## Author

**Ayan Malakar**  
Bioinformatics Scientist  
[GitHub](https://github.com/ayan-malakar)

---

## License

MIT — see [LICENSE](LICENSE)
