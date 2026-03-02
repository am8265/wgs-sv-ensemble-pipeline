# Pipeline Architecture

## Design Principles

1. **Ensemble approach** — no single SV caller is optimal for all SV types. Five callers with complementary detection methods maximizes sensitivity while SURVIVOR merging controls specificity.
2. **Modular WDL tasks** — each caller is an independent task, making it easy to add/remove callers or swap implementations.
3. **Depth-based QC** — Duphold annotations provide an independent signal (read depth fold-change) to filter false positives before merging.
4. **Cloud-first** — all tasks use Docker containers and are compatible with AWS HealthOmics, Terra, and local Cromwell.

## Detection Methods

| Method | Callers | Strengths |
|--------|---------|-----------|
| Split-read | Manta, Delly, Lumpy | Precise breakpoints |
| Discordant pairs | Manta, Delly, Lumpy, BreakDancer | Sensitive for large SVs |
| Read depth | CNVnator | CNV detection, copy number |

## Duphold Thresholds

- `DHBFC < 0.7` for deletions — depth fold-change relative to bins flanking the call
- `DHFFC > 1.25` for duplications — fold-change relative to flanking sequence
- INV and BND calls pass through without depth filtering

## SURVIVOR Parameters

```
max_dist    = 1000   # Breakpoints within 1kb are merged
min_support = 2      # At least 2 callers must support the call
sv_type     = 1      # Only merge same SVTYPE (DEL with DEL, etc.)
sv_strand   = 1      # Consider strand orientation
min_size    = 50     # Minimum 50bp — filters noise
```

## AWS HealthOmics Setup

1. Create an S3 bucket for inputs/outputs
2. Create an IAM role with HealthOmics + S3 permissions
3. Zip the `wdl/` directory and upload as a workflow definition
4. Use `cromwell_aws.conf` for Cromwell-based submission, or submit natively via AWS CLI

## Terra Setup

1. Link your GitHub repo in Terra Workspace
2. Set the root entity as `sample`
3. Import `inputs/inputs.json` as the data table
4. Submit via Terra UI
