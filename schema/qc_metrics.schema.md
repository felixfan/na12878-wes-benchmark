# `qc_metrics` — sequencing quality table (Table 1 of the paper: three-platform quality comparison)

> One row = one (sample × reference) combination, holding raw-read, alignment and capture quality metrics.
> Reporting these is a requirement for a *Scientific Data* Data Descriptor.
>
> [`scripts/run_qc.sh`](../scripts/run_qc.sh) runs the tools on the server and produces the reports;
> [`scripts/aggregate_qc.py`](../scripts/aggregate_qc.py) collects them into this table.
> The authoritative column order is [`qc_metrics.header.csv`](qc_metrics.header.csv).

## Data flow

```
raw fastq  ── fastp ─────────────────▶ fastp.json          ┐
full-depth bam ── samtools flagstat ─▶ flagstat.txt        ├─▶ aggregate_qc.py ─▶ qc_metrics.csv
               ── Picard MarkDuplicates ─▶ dup_metrics.txt │
               ── GATK CollectHsMetrics ─▶ hs_metrics.txt  ┘  (needs V6 BED converted to interval_list)
```

> The full-depth BAM is the calibration-point BAM
> (`<WORKDIR>/bam/<sample>.<reference>.bwa-mem.full.s0/dedup.bam`). Raw-fastq QC does not depend on the
> reference build, so those values are repeated across the two build rows of a sample.

## Column dictionary

| Column | Source | Meaning |
|---|---|---|
| `sample` / `platform` / `reference` | identity | primary key (`platform` is 1:1 with `sample`; three platforms) |
| `raw_read_pairs` | fastp `before_filtering.total_reads` / 2 | raw read pairs |
| `raw_bases_gb` | fastp `total_bases` / 1e9 | raw base yield (target was ≥100 Gb) |
| `q30_pct` | fastp `q30_rate` × 100 | **Q30 fraction** — the headline base-quality metric |
| `gc_pct` | fastp `gc_content` × 100 | GC content |
| `mean_read_len` | fastp `read1_mean_length` | mean read length (PE150, so ~150) |
| `mapped_pct` | samtools flagstat | mapping rate (>99% argues against non-human contamination) |
| `properly_paired_pct` | samtools flagstat | properly paired fraction |
| `dup_pct` | Picard MarkDuplicates `PERCENT_DUPLICATION` × 100 | duplicate rate |
| `on_target_pct` | CollectHsMetrics `PCT_SELECTED_BASES` × 100 | **on-target enrichment** (on-bait + near-bait over aligned bases) |
| `mean_target_cov` | CollectHsMetrics `MEAN_TARGET_COVERAGE` | mean depth over the target at full depth |
| `fold_enrichment` | CollectHsMetrics `FOLD_ENRICHMENT` | capture enrichment factor |
| `fold80_penalty` | CollectHsMetrics `FOLD_80_BASE_PENALTY` | coverage uniformity (closer to 1 is more uniform) |
| `pct_target_10x` / `20x` / `30x` | CollectHsMetrics `PCT_TARGET_BASES_XX` × 100 | fraction of target bases at ≥10×/20×/30× |
| `zero_cvg_target_pct` | CollectHsMetrics `ZERO_CVG_TARGETS_PCT` × 100 | fraction of targets with zero coverage |
| `fastp_version` / `hsmetrics_tool` / `run_date` | provenance | fastp 1.3.3 / GATK 4.3.0.0 CollectHsMetrics |

> `CollectHsMetrics` requires an **interval_list**, not a BED:
> `gatk BedToIntervalList -I <V6.bed> -O <V6.<reference>.interval_list> -SD <reference.dict>`,
> once per build (`run_qc.sh` creates and caches it). Both bait and target use the V6 Covered design.
>
> `on_target_pct` and `dup_pct` are also written back into
> [`run_coverage`](../results/run_coverage_all.csv), where they are otherwise `NA`.
