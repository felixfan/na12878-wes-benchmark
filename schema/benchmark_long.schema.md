# `benchmark_long` — data dictionary for the aggregated results table

> **One row = one run × filter level × variant type × genotype class × region stratum.**
> This is the single table every downstream analysis consumes (saturation curves, platform comparison,
> CPU↔GPU concordance), together with the companion table `run_coverage`.
>
> The authoritative column order is [`benchmark_long.header.csv`](benchmark_long.header.csv);
> a synthetic example is in [`benchmark_long.example.csv`](benchmark_long.example.csv).
> If a column changes, change those two files and this document together.

---

## 0. Data flow

```
                        one row per run
   run_manifest.csv  ─────────────────────────┐  (run identity: platform / caller / seed / versions / BED / hap.py prefix)
                                              │
   <happy_prefix>.extended.csv ───────────────┤  (preferred; headline + het/hom + difficult regions, no run identity)
   <happy_prefix>.summary.csv   (fallback)    │  (headline only, used when extended.csv is absent)
                                              ▼
                          scripts/aggregate_happy.py    <- joins identity to metrics, normalises names
                                              │
                                              ▼
                            benchmark_long_all.csv      (this schema)
                                              │
   run_coverage.csv ──────────────────────────┤  (one row per run: on-target depth, ≥10x/20x fractions)
                                              ▼
                        plotting (pandas/matplotlib or R/ggplot)
```

**The key point:** hap.py output files contain **no** platform, caller, seed or tool-version information.
Those must be carried in by `run_manifest.csv`. Every run is therefore registered in the manifest first,
and the aggregation script joins identity to metrics afterwards. This join is the reproducibility anchor
of the whole pipeline.

---

## 1. Grain and primary key

A row is uniquely determined by this **composite key**:

```
run_id + benchmark_tool + filter + variant_type + genotype_class + region_stratum
```

- `run_id` identifies one "alignment + caller + depth + seed" analysis.
- The other four keys are **stratification dimensions**: a single run expands into many rows
  (ALL/PASS × SNV/INDEL × all/het/hom × all/difficult regions).

**`run_id` grammar (deterministic — do not improvise):**

```
{sample}.{reference}.{aligner}.{caller}.{hardware}.{X}x.s{seed}

  e.g.  ZY_Illumina_NA12878.grch38.bwa-mem.gatk4_hc.cpu.30x.s42
        ZY_BGI_NA12878.grch38.pbrun_fq2bam.pbrun_deepvariant.gpu.50x.s44

  Full-depth (not downsampled) runs use target_depth_x = full and seed = 0.
```

> The titration is performed in **on-target-X space**, not in gigabases (see
> [`EXPERIMENT-DESIGN.md`](../EXPERIMENT-DESIGN.md) §2): `fraction = target_X / X_full`, where `X_full`
> is measured during calibration. The depth field of `run_id` is therefore the **target** on-target depth
> (a round grid point); the realised depth is recorded in `mean_target_depth_x`.

---

## 2. Column dictionary

Legend: **required** = every row carries a value; `NA` = not applicable (e.g. `parabricks_version` on a CPU
row); *(blank)* = hap.py does not provide that metric for that stratum.

### 2.1 Run identity and comparison axes (from `run_manifest`)

| # | Column | Type | Allowed values / example | Req. | Notes |
|---|---|---|---|---|---|
| 1 | `run_id` | str | see §1 | yes | part of the composite key; also the join key to `run_coverage` |
| 2 | `sample` | enum | `ZY_Illumina_NA12878` \| `ZY_BGI_NA12878` \| `ZY_GeneMind_NA12878` | yes | the three same-DNA libraries |
| 3 | `platform` | enum | `illumina_novaseq6000` \| `bgi_dnbseq_t7` \| `genemind_surfseq5000` | yes | **primary comparison axis**; 1:1 with `sample` (redundant, but convenient for grouping) |
| 4 | `reference` | enum | `b37` \| `grch38` | yes | `b37` = GRCh37 legacy (`human_g1k_v37`, no `chr` prefix); GRCh38 is primary |
| 5 | `aligner` | enum | `bwa-mem` \| `pbrun_fq2bam` | yes | CPU vs GPU branch; see §5 |
| 6 | `caller` | enum | `gatk4_hc` \| `deepvariant` \| `pbrun_haplotypecaller` \| `pbrun_deepvariant` (also reserved: `gatk3`, `bcftools`, `freebayes`) | yes | this release contains the first four |
| 7 | `hardware` | enum | `cpu` \| `gpu` | yes | grouping key for CPU↔GPU concordance |
| 8 | `target_depth_x` | float \| `full` | `5`, `8`, `10`, `12`, `15`, `20`, `25`, `30`, `40`, `50`, `60`, `75`, `100`, `150`, `200`, `300`, `400`, `full` | yes | the downsampling **target** depth (the grid point); discrete and aligned across platforms; seed aggregation and platform matching use this column |
| 9 | `mean_target_depth_x` | float | e.g. `46.1` | yes | **the continuous axis**: measured mean depth over the V6 target (close to `target_depth_x`, varying slightly with seed). Saturation curves use this for x |
| 10 | `seed` | int | `42`, `43`, `44` (full depth = `0`; `1` marks excluded acceptance replicates) | yes | downsampling seed; three seeds per sub-full point give error bars |

### 2.2 Stratification dimensions (from hap.py output)

| # | Column | Type | Allowed values | Req. | Notes |
|---|---|---|---|---|---|
| 11 | `benchmark_tool` | enum | `happy` \| `vcfeval` | yes | evaluation tool; metric semantics are aligned between the two (§4) |
| 12 | `filter` | enum | `ALL` \| `PASS` | yes | **must be recorded**: `ALL` includes unfiltered calls, `PASS` only those passing filters. Precision differs substantially between them |
| 13 | `variant_type` | enum | `SNV` \| `INDEL` | yes | hap.py's `SNP` is normalised to `SNV` |
| 14 | `genotype_class` | enum | `all` \| `het` \| `hom` | yes | see §6; `hom` corresponds to hap.py's `homalt` |
| 15 | `region_stratum` | str | `all` \| `difficult_lowmappability` \| `segdup` \| `homopolymer` \| … | yes | `all` = the whole high-confidence region; difficult regions come from the GIAB genome stratifications |

### 2.3 Counts and metrics (from hap.py)

| # | Column | Type | Req. | hap.py source | Notes |
|---|---|---|---|---|---|
| 16 | `truth_total` | int | yes | `TRUTH.TOTAL` | variants in the truth set for this stratum (= TP + FN) |
| 17 | `tp` | int | yes | `TRUTH.TP` | true positives |
| 18 | `fn` | int | yes | `TRUTH.FN` | false negatives (missed) |
| 19 | `query_total` | int | yes | `QUERY.TOTAL` | calls made by the query VCF in this stratum |
| 20 | `fp` | int | yes | `QUERY.FP` | false positives |
| 21 | `unk` | int | opt | `QUERY.UNK` | calls outside the high-confidence region, ignored (absent for vcfeval) |
| 22 | `fp_gt` | int | opt | `FP.gt` | false positives with the right position but the wrong genotype (absent for vcfeval) |
| 23 | `precision` | float [0,1] | yes | `METRIC.Precision` | TP / (TP + FP) |
| 24 | `recall` | float [0,1] | yes | `METRIC.Recall` | TP / (TP + FN), i.e. sensitivity |
| 25 | `f1` | float [0,1] | yes | `METRIC.F1_Score` | 2·P·R / (P + R) |
| 26 | `frac_na` | float [0,1] | opt | `METRIC.Frac_NA` | UNK / QUERY.TOTAL, the fraction outside the high-confidence region (absent for vcfeval) |

> Required columns are present on **every** row, including het/hom and difficult-region rows.
> **het/hom rows carry full Precision/Recall/F1**, computed from the `.het` / `.homalt` count columns of
> `extended.csv` (see §4) — they are not left blank.

### 2.4 Provenance and pinned versions (from `run_manifest`)

| # | Column | Type | Req. | Notes |
|---|---|---|---|---|
| 27 | `aligner_version` | str | yes | e.g. `0.7.15` (BWA-MEM), `4.7.0` (Parabricks `fq2bam`) |
| 28 | `tool_version` | str | yes | **caller version**, e.g. `4.3.0.0` (GATK4), `1.9.0` (DeepVariant) |
| 29 | `parabricks_version` | str | yes | version on GPU runs; `NA` on CPU runs |
| 30 | `benchmark_tool_version` | str | yes | hap.py / rtg-tools version |
| 31 | `giab_version` | str | yes | e.g. `HG001_v4.2.1` |
| 32 | `target_bed` | str | yes | identifier/filename of the V6 ∩ GIAB high-confidence BED |
| 33 | `run_date` | date | yes | `YYYY-MM-DD` |

> Columns 27–31 exist because the CPU↔GPU comparison is only interpretable if Parabricks, GATK4 and the
> DeepVariant model versions are pinned and matched across the two branches.

---

## 3. Companion table `run_coverage` (one row per run)

Mean on-target depth, the ≥10×/≥20× covered fractions and similar quantities are **one value per run**;
storing them in the long table would repeat them on every stratum row, so they live in their own table,
joined on `run_id` (or on the BAM key `sample.reference.aligner.target_depth_x.seed`). See
[`run_coverage.example.csv`](run_coverage.example.csv).

Columns actually produced by `run_one.sh` via mosdepth over the V6 BED: `target_depth_x` (the grid point),
`mean_target_depth_x` (measured), `used_gb` (bases actually used by that run), and
`pct_target_10x` / `20x` / `30x` / `50x`. `median_target_depth_x` and the optional `on_target_pct` /
`dup_pct` (obtainable from Picard `CollectHsMetrics` / duplicate metrics) are left as `NA`.

> `mean_target_depth_x` is stored in both tables — redundantly in the long table so that saturation curves
> can be plotted without a join. The two must agree; `run_coverage` is authoritative.

---

## 4. Mapping from hap.py, and known limitations

> Column names and extraction rules were verified against real hap.py v0.3 output
> (`example/happy/expected-qfy.{summary,extended}.csv` in the hap.py source distribution).

- **Source priority.** `extended.csv` is a **superset** of `summary.csv` (its `Subtype=* & Subset=*` rows
  match summary cell for cell), so the script parses extended first and obtains all three stratifications
  at once. `summary.csv` is used only as a fallback, and then yields headline rows only.
- **Type.** hap.py `SNP` → `SNV`; `INDEL` → `INDEL`. Other rows are ignored (`Subtype≠*`, i.e. ti/tv and
  INDEL length classes, are not taken).
- **Filter.** hap.py emits both `ALL` and `PASS` rows for each type; both are kept.
- **Operating point.** In `extended.csv` the summary row of each group (no quality threshold) has `QQ='*'`;
  with `--roc` the same group also contains numeric-`QQ` ROC points. The script takes the `QQ='*'` row
  (`_operating_point()`).
- **Headline rows** (`region=all, genotype=all`) use the `METRIC.*` values hap.py reports directly.
- **het/hom rows.** hap.py does **not** split along the `Genotype` dimension (it only emits `*`); instead it
  provides `.het` / `.homalt` **count suffix columns** (`TRUTH.TP.het`, `QUERY.FP.homalt`, …). The script
  derives **both recall and precision** from these (`_suffix_metrics()`): recall = TP.gt / TRUTH.TOTAL.gt,
  precision = QUERY.TP.gt / (QUERY.TP.gt + QUERY.FP.gt). Note that `summary.csv` does **not** contain these
  `.het` columns (only `het_hom_ratio`), which is why het/hom stratification requires `extended.csv`.
- **Region strata** (`region≠all`) come from the `Subset` dimension of `extended.csv` and require
  `hap.py --stratification <giab-strat.tsv>`. `TS_boundary` and `TS_contained` are skipped (they are
  hap.py's own target-boundary subsets, not GIAB difficult regions). Implemented in
  `parse_happy_extended()`.
- **Granularity warning.** A complete GA4GH stratification set produces thousands of subsets, including
  per-gene exons such as `exons_ABCG1_1`. The script records whatever it is given; **which difficult
  regions appear is decided by the stratification set passed to hap.py** (this design uses only the
  top-level ones: segmental duplications, low mappability, homopolymers, extreme GC).
- **rtg vcfeval** as an alternative to hap.py: `True-pos-baseline` → `tp`, `False-neg` → `fn`,
  `False-pos` → `fp`, `Sensitivity` → `recall`, `Precision` → `precision`, `F-measure` → `f1`, taking the
  `None` (no-threshold) operating point. Such rows are marked `benchmark_tool=vcfeval`. Implemented in
  `aggregate_happy.py:parse_vcfeval()`, reading `<prefix>.snv.summary.txt` / `<prefix>.indel.summary.txt`
  (vcfeval's summary does not split by variant type, so it must be run once per type).

---

## 5. Deriving CPU↔GPU concordance from this table

No extra columns are needed: select the semantically equivalent caller pairs
(`gatk4_hc` ↔ `pbrun_haplotypecaller`, `deepvariant` ↔ `pbrun_deepvariant`), match on
`sample` / `reference` / `target_depth_x` / `seed`, pair by `hardware`, and compare F1, precision and
recall. This is exactly why the `hardware` column exists.

For the site-level comparison (which variants each pipeline called, independent of the truth set) see
`scripts/vcf_concordance.py` and the released `pair_concordance.tsv` / `pair_quality.tsv`.

---

## 6. Validation rules (enforced by `aggregate_happy.py`)

1. The composite key is unique — no duplicate rows.
2. `precision`, `recall`, `f1`, `frac_na` ∈ [0, 1] (or blank).
3. `f1` agrees with `2·P·R / (P + R)` to within 0.01 (floating-point tolerance).
4. `tp + fn == truth_total`.
5. Enum columns contain only allowed values.
6. `hardware=gpu` implies `parabricks_version ≠ NA`.
