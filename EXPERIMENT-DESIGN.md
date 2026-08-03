# Experimental design — NA12878 three-platform WES benchmark

> This document covers the **scientific design**: what is measured and why.
> For column-by-column definitions of the results table see [`schema/benchmark_long.schema.md`](schema/benchmark_long.schema.md);
> for how to run the pipeline see [`README.md`](README.md).
>
> The resource is a **Data Descriptor**: its value is a reusable, technically validated dataset.
> The depth titration and pipeline comparisons are **technical validation**, not claims of new biology.

---

## 1. Questions the design answers

Same-source NA12878 DNA, deep WES (PE150) on three platforms, evaluated throughout against the
GIAB HG001 v4.2.1 benchmark (Precision / Recall / F1):

1. **How sequencing depth affects variant-calling accuracy** — downsampling titration, producing
   F1-versus-on-target-depth saturation curves and their knee/plateau ("how much depth does WES need?").
2. **How analysis pipelines differ** — caller/pipeline, platform (Illumina / BGI / GeneMind),
   reference build (GRCh38 / b37) and CPU-versus-GPU (Parabricks).

**Four comparison axes:** platform × **on-target depth (X)** × caller/pipeline × reference build.

> **Benchmark scope (a limitation that the paper states explicitly).** The GIAB HG001 v4.2.1
> high-confidence set is the `1_22` release: it covers **autosomes 1–22 only, not X or Y**
> (NA12878 is female; each build's truth set contains ~3.89 M variants on contigs 1–22 only).
> All accuracy figures are therefore defined on autosomal exons (V6 ∩ GIAB 1–22 high-confidence
> regions); variants on X/Y exons are not evaluated.

---

## 2. Depth titration (the core of the design)

### 2.1 Why the axis is mean on-target depth (X), not gigabases

The downsampling knob is the fraction of reads retained (hence gigabases), but **accuracy tracks
on-target coverage depth**, and for this dataset the two are strongly non-linear:

- The V6 capture target is ≈ **60 Mb**; the raw data are ≈ **≥100 Gb per sample** (deliberately deep,
  to leave room for titration).
- mean on-target X ≈ `raw bases × on-target% × mapped% / target size` ≈ `100e9 × ~0.6 / 60e6`
  ≈ **~1000× at full depth**.
- So **~1 Gb ≈ ~10×**. The knee and plateau of WES accuracy fall at **~5–60×**, i.e. within the
  lowest ~1–6 Gb; anything above ~10 Gb (>100×) is essentially plateau.

A linear gigabase ladder (1, 2, … 100 Gb) would therefore spend most of its points on a redundant
plateau above 200× while compressing the entire knee into the first few gigabases — precisely the
region the "how much depth is enough" conclusion depends on. Worse, **equal gigabases is not equal
depth** across platforms (their on-target percentages differ), which would violate the requirement
that platforms be compared at matched on-target depth.

The titration is therefore performed directly in on-target-X space.

### 2.2 Two phases: calibrate, then titrate

| Phase | What it does | Output |
|---|---|---|
| **Phase 0 — calibration** | For each (sample, reference), align the full read set once and measure the target-region depth with mosdepth | `calibration.csv`: sample, reference, `x_full`, `gb_full` |
| **Phase 1 — titration** | For each target depth X_t: `fraction = X_t / X_full`; `seqtk sample -s<seed> <fraction>`, then the full pipeline | one VCF and one hap.py result set per run |

Downsampling is uniform, so the realised depth is ≈ `fraction × X_full ≈ X_t`; the actual
`mean_target_depth_x` is recorded per run and used as the continuous axis. Calibration is needed only
for GRCh38 (the titration runs there; b37 is only run at full depth). The GPU branch (`fq2bam`) reuses
the BWA-derived `X_full` — both embed BWA-MEM 0.7.15, so alignment and depth are near-identical and the
residual is far smaller than the titration step size (see §4).

**Measured values:** Illumina 627×, BGI 742×, GeneMind 687× at full depth.

### 2.3 Depth grid (on-target X, denser at the knee)

```
5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 75, 100, 150, 200, 300, 400, full
└──── dense at the knee (8 points ≤30×) ────┘└─ sparse plateau ─┘└ gap fill ┘ └ full (calibration point)
```

Eighteen points: eight at ≤30× to resolve the knee; a sparse 40–150× stretch to confirm saturation;
**200/300/400× added to fill the 150×-to-full gap**, which spans measured depths of ~204× to ~660× and
is the widest interval in the grid (SNV accuracy is already on a plateau there, but INDEL F1 still
gains ~0.007–0.009, so the points characterise the slow climb rather than interpolating a straight
line). `full` (measured 627–742×) is the ceiling and the calibration point.

### 2.4 Replication and error bars

Each sub-full depth point is run with **three fixed random seeds, 42 / 43 / 44** (`seqtk -s42/-s43/-s44`;
42 is the convention for single-seed comparisons). R1 and R2 use the same seed so pairing is preserved,
and seqtk's PRNG makes adjacent seeds produce fully independent subsamples. This gives error bars on F1,
which matter most at low depth where variance is largest. Each run records its seed in the `seed` column
of `run_manifest.csv` / `run_coverage.csv`, so any result is reproducible. The `full` point involves no
sampling and therefore no randomness: **one run only** (seed 0).

### 2.5 Platforms are aligned at matched on-target depth

All three platforms run the **same target-X grid**, so every cross-platform comparison is made at matched
on-target depth. Saturation curves for the three platforms share one x-axis, and Precision/Recall/F1
comparisons at a given depth are matched on `target_depth_x`.

### 2.6 Stratified saturation curves and knee/plateau readout

Saturation is not a single curve: the long-format results table supports stratification by
**variant type (SNV/INDEL) × genotype (het/hom) × region (all / difficult)**. Expected behaviour:

- SNVs saturate early (knee ~20–30×), INDELs later (~30–40×); heterozygous sites need more depth than
  homozygous ones.
- Difficult regions (homopolymers, segmental duplications, low mappability, extreme GC) saturate later
  and plateau lower.
- For each curve, report the **knee** (smallest depth reaching ~99% of plateau F1) and the **plateau F1**.

---

## 3. The pipeline behind every run

```
raw fastq (R1, R2)
  └─ [titration: seqtk sample -s<seed> <frac = X_t / X_full>]        # sub-full points only
     └─ align + mark duplicates:
          CPU  bwa-mem 0.7.15 | samtools sort | gatk MarkDuplicates | samtools index
          GPU  pbrun fq2bam
        └─ [BQSR: GATK family only] BaseRecalibrator -> ApplyBQSR -> recal.bam
           (DeepVariant skips BQSR and uses dedup.bam)
        └─ calling, restricted to the V6 target:
             gatk4 HaplotypeCaller / DeepVariant (WES model) / pbrun haplotypecaller / pbrun deepvariant
           └─ bcftools norm -m-any -f <reference>                    # left-align + split multi-allelics
              └─ hap.py against GIAB ∩ V6 (vcfeval engine, GA4GH stratification)
                 └─ aggregate_happy.py -> benchmark_long_all.csv -> plot_benchmark.py
```

Pinned tool versions are listed in [`README.md`](README.md) (Parabricks 4.7.0 / GATK 4.3.0.0 /
DeepVariant 1.9.0 / BWA-MEM 0.7.15).

---

## 4. The four test blocks

| Block | Purpose | Axes | Depth | Callers | Reference |
|---|---|---|---|---|---|
| **T1 — saturation** | depth versus accuracy (§2) | platform × depth × seed | **full X grid** | gatk4_hc + deepvariant (CPU) | GRCh38 |
| **T2 — legacy caller** (deferred) | which pipeline performs best | platform × caller | full | gatk3 / bcftools / freebayes — **not in this release**, to be added if reviewers request | GRCh38 |
| **T3 — CPU↔GPU** | GPU accelerates without losing accuracy | caller × hardware | 5/10/15/30/50/100/200× + full | gatk4_hc↔pbrun_hc, deepvariant↔pbrun_dv, bwa+Picard↔fq2bam | GRCh38 |
| **T4 — reference build** | GRCh38 primary, b37 legacy | reference | full | gatk4_hc + deepvariant | b37 |

**On T3.** The GPU branch reuses the BWA-derived `X_full`: `fq2bam` embeds the **same BWA-MEM 0.7.15**, so
alignment and depth agree far more closely than the titration step size — and the concordance analysis
quantifies exactly this. Versions must be pinned and the V6 BED shared, otherwise CPU-GPU differences are
not interpretable.

**On T4 (coordinate systems).** The legacy axis is **b37** (`human_g1k_v37`, no `chr` prefix), which aligns
natively with the GIAB HG001 GRCh37 truth set. The V6 BED ships in hg19 coordinates with `chr`, so the b37
copy has the prefix stripped (`sed 's/^chr//'`). On the GRCh38 route all three files (reference, truth set,
lifted-over V6) carry `chr`, so the issue does not arise.

**BQSR scope.** BQSR is applied to the **GATK family only** (gatk4_hc / pbrun_hc;
`BaseRecalibrator -> ApplyBQSR`, known-sites = dbSNP + Mills + 1000G/known indels, chosen per build).
**DeepVariant does not use it** — the DeepVariant authors advise against it and the CNN does not benefit.
Consequently `dedup.bam` is shared across callers while `recal.bam` is GATK-specific. This follows the
principle that each pipeline runs under its own best practice; the caller comparison therefore contains
BQSR as a covariate, which the paper states explicitly.

**Stratification** (passed to hap.py `--stratification`) uses only the top-level difficult regions:
`LowComplexity` (homopolymers and tandem repeats), `GCcontent`, `mappability` + `SegmentalDuplications`,
and `FunctionalRegions` (CDS). The 2000+ per-gene exon subsets are not used.

---

## 5. Run count (450)

| Block | Composition | Runs |
|---|---|---|
| T1 | 3 platforms × 2 callers × (17 sub-full points × 3 seeds + 1 full point) = 3×2×52 | **312** |
| T3 | 3 platforms × {pbrun_hc, pbrun_dv} × ({5,10,15,30,50,100,200}× × 3 seeds + full) = 3×2×22 | **132** |
| T4 | 3 platforms × {gatk4_hc, deepvariant} × b37 full | **6** |
| T2 | gatk3 legacy — **deferred** | — |
| | **Total** | **450** |

**Why not the full factorial (1,248).** The full factorial is 3 platforms × 52 depth-and-seed combinations
× 4 caller-and-hardware configurations × 2 reference builds = 1,248 cells, of which 450 (**36%**) were run.
The 798 cells deliberately left out are: (i) GPU × GRCh38 at the remaining 30 depth-and-seed combinations
(180) — equivalence is established at both extremes (5× and full), and intermediate depths interpolate;
(ii) b37 × CPU depth titration (306) — the build effect is measured at full depth (|ΔF1| ≤ 0.01) and a
depth × build interaction is not a question this dataset sets out to answer; (iii) b37 × GPU entirely (312),
the compound of the first two. A practical constraint pointed the same way: running the GPU branch over the
full grid would have required ~2.55 TB of additional BAM storage.

Note that "4 caller-and-hardware configurations" means **2 algorithms × 2 hardware implementations**
(`pbrun haplotypecaller` is the GPU implementation of GATK4 HaplotypeCaller on the same 4.3.0.0 engine),
not four independent callers — that correspondence is what makes the CPU↔GPU analysis interpretable.
See Supplementary Table S5 of the paper.

**Why T3 was widened to seven depth points.** The original {15, 30, 50} + full was too sparse. Since
CPU↔GPU disagreements were known to concentrate at low depth (median DP 10.5× at discordant versus 45.0× at
concordant sites), the low end is the informative stress test: 5× and 10× were added (if equivalence holds
there, the conclusion is considerably stronger; if it breaks, that is a finding that must be reported).
100× was added as a commonly used clinical WES depth and 200× to align with the new high-depth points.
Matched pairs increased from 60 to 132.

`scripts/gen_plan.sh` enumerates and de-duplicates the plan, emitting the exact counts and `jobs.list`.

---

## 6. Design assumptions that calibration resolved

The `X_full` figures in §2.1 (~1000×) were an estimate from 60 Mb of target, ~100 Gb of raw data and ~60%
on-target. Phase-0 calibration measured **627× (Illumina), 742× (BGI) and 687× (GeneMind)** — lower than the
estimate but ample for the grid, whose top sub-full point (400×) sits comfortably below the lowest platform's
full depth. The three platforms differ in `X_full` because their on-target percentages differ; because the
titration targets a depth rather than a read fraction, the grid points still align across platforms, and only
the per-platform sampling fraction differs.

---

*Companion files:* [`README.md`](README.md) (how to run) ·
[`scripts/gen_plan.sh`](scripts/gen_plan.sh) · [`scripts/run_one.sh`](scripts/run_one.sh) ·
[`scripts/calibrate.sh`](scripts/calibrate.sh) · [`schema/`](schema/) (data dictionary).
