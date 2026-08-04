# Deep WES of NA12878 (HG001) on three platforms — benchmarking framework and processed results

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21792674.svg)](https://doi.org/10.5281/zenodo.21792674)
[![License: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)

Analysis code and processed benchmarking outputs accompanying the *Scientific Data* Data Descriptor
**"Deep whole-exome sequencing of NA12878 (HG001) on three platforms for germline variant-calling benchmarking."**

The raw sequencing reads are deposited in the NCBI Sequence Read Archive under BioProject
**PRJNA1506792** (study SRP723789; BioSample SAMN62149305; Runs **SRR39986369** Illumina, **SRR39986368** BGI, **SRR39986367** GeneMind), held until publication.
This repository holds everything downstream of them: the pipeline that produced the 450 benchmarking runs,
and every table and figure reported in the paper. Nothing here is under embargo — the processed results,
QC tables and code are openly available now.

---

## What is here

| Path | Contents |
|---|---|
| `scripts/` | The analysis framework that produced the runs (shell + Python). Every run is fully determined by its `run_id`. |
| `scripts/jobs/` | The run lists actually executed (`jobs_*.list`). |
| `configs/` | `na12878.paths.template.sh` - data/resource locations. Copy, fill in the `/PATH/TO/...` placeholders, use. |
| `schema/` | Authoritative data dictionary for the processed result tables. |
| `results/` | Processed pipeline output: the full accuracy, coverage, timing, concordance and QC tables (see below). |
| `results/plan/` | `run_manifest.csv` - the registration record of every run that was executed, and the join key between the pipeline and the results. |
| `source_data/` | **One file per published figure and per published table** - the exact values behind Figures 2-7 and Tables 1-4, S1-S5. |
| `figure_code/` | **One script per published figure**, reading only the matching `source_data/` file. |
| `figures/` | The published figures `fig1.png` ... `fig7.png` (300 dpi). Figure 1 also ships as `fig1.pptx`, a fully editable PowerPoint file - it is a schematic, not a data plot. |
| `reference/` | Agilent SureSelect V6 Covered target BED, in GRCh38 and GRCh37/b37 coordinates. |
| `EXPERIMENT-DESIGN.md` | The experimental design: depth grid, blocks, and why the matrix is nested rather than full-factorial. |

### Key result files

| File | Rows | What it is |
|---|---|---|
| `results/benchmark_long_all.csv` | 5,448 | **The main table.** One row per run x filter x variant type x genotype class. 454 runs (450 analysed + 4 acceptance replicates with seed `s1`, excluded from all analyses). |
| `results/run_coverage_all.csv` | 454 | Per-run measured on-target depth and >=10/20/30/50x fractions (mosdepth). |
| `results/resource_usage_all.csv` | 454 | Per-run wall-clock time, one row per `run_id`. Seven runs were executed more than once; their timings are collapsed to the median. |
| `results/pair_concordance.tsv` | 402 | Site-by-site CPU-GPU comparison, 134 pairs x {SNV, INDEL, ALL}. |
| `results/pair_quality.tsv` | 536 | DP/GQ/QUAL by concordance category, 134 pairs x 4 categories. |
| `results/qc/qc_metrics.csv` | 3 | Per-library sequencing QC (fastp, samtools flagstat, Picard MarkDuplicates and CollectHsMetrics). |
| `results/qc/*` | - | Sample identity (somalier), contamination (VerifyBamID2) and the SNV hard-filter test. |
| `results/plan/run_manifest.csv` | 454 | **The execution record.** One row per `run_id`: platform, reference, aligner, caller, hardware, target and measured depth, seed, pinned tool versions, GIAB version, target BED. `aggregate_happy.py` joins it to the hap.py output to build the main table, so it is the reproducibility anchor of the whole chain. Runs executed more than once are recorded once (last registration). |
| `results/sra_file_manifest.tsv` | 6 | Raw FASTQ sizes (exact bytes) and MD5 checksums. |

Analyses in the paper exclude the four `s1` acceptance-replicate runs, giving **450 runs** and **132 CPU↔GPU pairs**.

---

## The `run_id` scheme

Every result traces back to its inputs through a deterministic identifier:

```
{sample}.{reference}.{aligner}.{caller}.{hardware}.{depth}.s{seed}
```

| Field | Allowed values |
|---|---|
| sample | `ZY_Illumina_NA12878`, `ZY_BGI_NA12878`, `ZY_GeneMind_NA12878` |
| reference | `grch38`, `b37` |
| aligner | `bwa-mem` (CPU), `pbrun_fq2bam` (GPU) |
| caller | `gatk4_hc`, `deepvariant`, `pbrun_haplotypecaller`, `pbrun_deepvariant` |
| hardware | `cpu`, `gpu` |
| depth | `5x` `8x` `10x` `12x` `15x` `20x` `25x` `30x` `40x` `50x` `60x` `75x` `100x` `150x` `200x` `300x` `400x` `full` |
| seed | `42` `43` `44` (sub-full), `0` (full depth, no sampling randomness), `1` (acceptance replicates, excluded) |

Example: `ZY_BGI_NA12878.grch38.bwa-mem.deepvariant.cpu.30x.s42`

---

## Reproducing

**1. Point the config at your data.**

```bash
cp configs/na12878.paths.template.sh configs/na12878.paths.sh
$EDITOR configs/na12878.paths.sh          # fill in every /PATH/TO/... placeholder
```

**2. Calibrate** (measures each platform's native on-target depth, needed to convert a target depth into a sampling fraction):

```bash
HAPPY_CONFIG=configs/na12878.paths.sh scripts/calibrate.sh
```

**3. Enumerate the run plan:**

```bash
HAPPY_CONFIG=configs/na12878.paths.sh scripts/gen_plan.sh      # writes jobs.list + run_plan.csv
```

**4. Execute runs** (each is idempotent; completed steps are skipped on re-run):

```bash
while read r; do HAPPY_CONFIG=configs/na12878.paths.sh scripts/run_one.sh "$r" </dev/null; done < jobs.list
```

`run_one.sh` does: seqtk downsampling → alignment + duplicate marking → mosdepth → BQSR (GATK family only)
→ variant calling → INDEL-only hard-filtering (GATK family) → `bcftools norm` → hap.py benchmarking,
then appends to the manifest, coverage and resource tables.

**5. Aggregate and plot:**

```bash
python scripts/aggregate_happy.py --manifest results/plan/run_manifest.csv \
       --happydir <workdir>/happy --out results/benchmark_long_all.csv
python scripts/make_source_data.py            # results/ -> source_data/
Rscript figure_code/Fig2_saturation.R         # source_data/ -> figures/  (one script per figure)
```

The `happy_prefix` column of the shipped manifest is relative (`happy/<run_id>`); point
`--happydir` at wherever your hap.py output lives.

### Pinned tool versions

BWA-MEM 0.7.15 · GATK 4.3.0.0 · DeepVariant 1.9.0 · NVIDIA Parabricks 4.7.0 · samtools/bcftools 1.23.1 ·
mosdepth 0.3.3 · seqtk 1.5 · hap.py 0.3.15 with rtg-tools 3.12.1 · fastp 1.3.3 · MultiQC 1.35 ·
VerifyBamID2 2.0.1 · somalier 0.2.19

Results are only comparable across runs when these versions match. The CPU↔GPU comparison in particular
relies on Parabricks embedding the same GATK 4.3.0.0 and DeepVariant 1.9.0 engines as the CPU branch.

---

## Source data and figures

Every published figure and table has its **own** data file under `source_data/`, and every figure has its
**own** script under `figure_code/`. Nothing is bundled: to check one number in one panel you open one CSV.

| Published item | Data file | Code |
|---|---|---|
| Figure 1 - study design | *(schematic, no data)* | *(no code - edit `figures/fig1.pptx` directly)* |
| Figure 2 - depth saturation | `source_data/Fig2_saturation.csv` | `figure_code/Fig2_saturation.R` |
| Figure 3 - platform comparison at 30x | `source_data/Fig3_platform.csv` | `figure_code/Fig3_platform.R` |
| Figure 4 - INDEL hard-filtering | `source_data/Fig4_hardfilter.csv` | `figure_code/Fig4_hardfilter.R` |
| Figure 5 - CPU-GPU concordance | `source_data/Fig5_concordance.csv` | `figure_code/Fig5_concordance.R` |
| Figure 6 - CPU-GPU acceleration | `source_data/Fig6_speedup.csv` | `figure_code/Fig6_speedup.R` |
| Figure 7 - reference build | `source_data/Fig7_reference.csv` | `figure_code/Fig7_reference.R` |
| Table 1 - sequencing QC | `source_data/Table1_QC.csv` | - |
| Table 2 - raw files, sizes, MD5 | `source_data/Table2_files.csv` | - |
| Table 3 - site-level concordance | `source_data/Table3_concordance.csv` | - |
| Table 4 - quality of concordant and discordant sites | `source_data/Table4_quality.csv` | - |
| Table S1 - per-pair quality (132 pairs) | `source_data/TableS1_pair_quality.csv` | - |
| Table S2 - identity and contamination | `source_data/TableS2_identity_contamination.csv` | - |
| Table S3 - SNV hard-filter test | `source_data/TableS3_snv_hardfilter.csv` | - |
| Table S4 - nominal versus measured depth | `source_data/TableS4_depth_mapping.csv` | - |
| Table S5 - design matrix | `source_data/TableS5_design_matrix.csv` | - |

**Redraw a figure** (writes a 300 dpi PNG into `figures/`):

```bash
Rscript figure_code/Fig2_saturation.R
```

The figure scripts are **R only** (`ggplot2` and `patchwork`; Fig 6 also needs `dplyr`). Each script
locates the repository relative to its own file, so it runs from any working directory, and writes
`figures/figN.png` at 300 dpi.

No title is drawn inside any figure: per the *Scientific Data* submission guidelines the title belongs in
the figure legend. Multi-panel figures are labelled (a), (b), ... Figure 1 is a schematic rather than a data
plot and has no code - edit `figures/fig1.pptx` in PowerPoint and re-export if you need to change it.

**Rebuild every source-data file from the processed results** (needs `pandas`):

```bash
python scripts/make_source_data.py
```

That regenerates all 15 files in `source_data/` from `results/`, so the chain
`results/ -> source_data/ -> figures/` is reproducible end to end.

---

## Notes for reusers

- **Compare at matched *measured* on-target depth, not at nominal depth or gigabases.** Downsampling removes
  proportionally fewer duplicate reads at smaller sampling fractions, so measured depth exceeds the nominal
  label by a factor that falls from ~1.5 at the lowest depths to 1.0 at full depth. Use
  `run_coverage_all.csv`; the mapping is tabulated in Supplementary Table S4 of the paper.
- **`PASS` means different things per caller family.** For the GATK family it is *raw SNVs plus
  hard-filtered INDELs* (SNV hard-filtering was measured to be net-negative and is deliberately not applied);
  DeepVariant uses its own internal model. For the GATK family the unfiltered set is under `filter=ALL`;
  for DeepVariant the `ALL` and `PASS` rows are identical, because the released call sets contain only
  records DeepVariant itself passed.
- **Benchmark scope.** All accuracy numbers are defined on the V6 Covered target ∩ GIAB HG001 v4.2.1
  high-confidence regions. Behaviour outside that intersection is not characterised.
- **Not covered by this dataset:** structural variants, copy-number variants, somatic calling, non-exonic
  regions, and anything requiring more than one individual.

---

## Licence

- **Code** (`scripts/`, `configs/`, `figure_code/`): MIT — see [`LICENSE`](LICENSE).
- **Data and documentation** (`results/`, `source_data/`, `figures/`, `schema/`, `reference/`,
  `EXPERIMENT-DESIGN.md`): CC BY 4.0 — see [`LICENSE-DATA`](LICENSE-DATA).

The Agilent SureSelect V6 target BED under `reference/` is redistributed for reproducibility of the
benchmark region; the capture design itself is Agilent's.

## Citing

See [`CITATION.cff`](CITATION.cff). Please cite the Data Descriptor, and the archived release if you
reuse the code or processed tables.

- **All versions (concept DOI):** [10.5281/zenodo.21792674](https://doi.org/10.5281/zenodo.21792674)
- **v1.0.1, the version described in the paper:** [10.5281/zenodo.21792675](https://doi.org/10.5281/zenodo.21792675)

The development repository is at <https://github.com/felixfan/na12878-wes-benchmark>.

## Third-party data

The GIAB HG001 v4.2.1 benchmark used as ground truth is from the Genome in a Bottle Consortium
(NIST/NCBI) and is **not** redistributed here — see the config template for the download location.
