# Deep WES of NA12878 (HG001) on three platforms — benchmarking framework and processed results

Analysis code and processed benchmarking outputs accompanying the *Scientific Data* Data Descriptor
**"Deep whole-exome sequencing of NA12878 (HG001) on three platforms for germline variant-calling benchmarking."**

The raw sequencing reads are in the NCBI Sequence Read Archive under BioProject **[PENDING-ACCESSION: PRJNA]**.
This repository holds everything downstream of them: the pipeline that produced the 450 benchmarking runs,
and every table and figure reported in the paper.

---

## What is here

| Path | Contents |
|---|---|
| `scripts/` | The complete analysis framework (shell + Python). Every run is fully determined by its `run_id`. |
| `scripts/jobs/` | The enumerated run lists actually executed (`jobs_*.list`). |
| `configs/` | `na12878.paths.template.sh` — data/resource locations. Copy, fill in the `/PATH/TO/...` placeholders, use. |
| `schema/` | Authoritative data dictionary for the results table (column-by-column definitions, validation rules). |
| `results/` | Processed benchmarking outputs (see below). |
| `results/qc/` | Sample identity (somalier), contamination (VerifyBamID2) and the SNV hard-filter test. |
| `results/figures/` | Figures as published (PNG, 300 dpi). |
| `reference/` | Agilent SureSelect V6 Covered target BED, in GRCh38 and GRCh37/b37 coordinates. |
| `source_data/` | `source_data.xlsx` (per-figure/table source values) and `figure_code/` (Python + R, one script per figure). |
| `EXPERIMENT-DESIGN.md` | The experimental design: depth grid, blocks, and why the matrix is nested rather than full-factorial. |

### Key result files

| File | Rows | What it is |
|---|---|---|
| `results/benchmark_long_all.csv` | 5,448 | **The main table.** One row per run × filter × variant type × genotype class. 454 runs (450 analysed + 4 acceptance replicates with seed `s1`, excluded from all analyses). |
| `results/run_coverage_all.csv` | — | Per-run measured on-target depth and ≥10/20/30/50× fractions (mosdepth). |
| `results/resource_usage_all.csv` | — | Per-run wall-clock time; the basis of the CPU-versus-GPU speed-up figures. |
| `results/pair_concordance.tsv` | 402 | Site-by-site CPU↔GPU comparison, 134 pairs × {SNV, INDEL, ALL}. |
| `results/pair_quality.tsv` | 536 | DP/GQ/QUAL by concordance category, 134 pairs × 4 categories. |
| `results/supplementary_cpu_gpu_quality.csv` | 134 | Supplementary Table S1 (per-pair quality, wide format). |

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
python scripts/aggregate_happy.py --manifest <workdir>/plan/run_manifest.csv --out benchmark_long_all.csv
python scripts/plot_benchmark.py --long benchmark_long_all.csv --outdir figures --resource resource_usage_all.csv
```

### Pinned tool versions

BWA-MEM 0.7.15 · GATK 4.3.0.0 · DeepVariant 1.9.0 · NVIDIA Parabricks 4.7.0 · samtools/bcftools 1.23.1 ·
mosdepth 0.3.3 · seqtk 1.5 · hap.py 0.3.15 with rtg-tools 3.12.1 · fastp 1.3.3 · MultiQC 1.35 ·
VerifyBamID2 2.0.1 · somalier 0.2.19

Results are only comparable across runs when these versions match. The CPU↔GPU comparison in particular
relies on Parabricks embedding the same GATK 4.3.0.0 and DeepVariant 1.9.0 engines as the CPU branch.

---

## Regenerating the figures

Each published figure has a standalone script in `source_data/figure_code/`, reading `source_data.xlsx`
and writing 300 dpi output. Python and R versions produce the same figure:

```bash
cd source_data/figure_code
python Fig2_saturation.py        # or: Rscript Fig2_saturation.R
```

Python needs `pandas matplotlib openpyxl`; R needs `readxl ggplot2 dplyr` (plus `tidyr` for Fig 3/4).

---

## Notes for reusers

- **Compare at matched *measured* on-target depth, not at nominal depth or gigabases.** Downsampling removes
  proportionally fewer duplicate reads at smaller sampling fractions, so measured depth exceeds the nominal
  label by a factor that falls from ~1.5 at the lowest depths to 1.0 at full depth. Use
  `run_coverage_all.csv`; the mapping is tabulated in Supplementary Table S4 of the paper.
- **`PASS` means different things per caller family.** For the GATK family it is *raw SNVs plus
  hard-filtered INDELs* (SNV hard-filtering was measured to be net-negative and is deliberately not applied);
  DeepVariant uses its own internal model. The unfiltered set is under `filter=ALL`.
- **Benchmark scope.** All accuracy numbers are defined on the V6 Covered target ∩ GIAB HG001 v4.2.1
  high-confidence regions. Behaviour outside that intersection is not characterised.
- **Not covered by this dataset:** structural variants, copy-number variants, somatic calling, non-exonic
  regions, and anything requiring more than one individual.

---

## Licence

- **Code** (`scripts/`, `configs/`, `source_data/figure_code/`): MIT — see [`LICENSE`](LICENSE).
- **Data and documentation** (`results/`, `schema/`, `reference/`, `source_data/source_data.xlsx`,
  `EXPERIMENT-DESIGN.md`): CC BY 4.0 — see [`LICENSE-DATA`](LICENSE-DATA).

The Agilent SureSelect V6 target BED under `reference/` is redistributed for reproducibility of the
benchmark region; the capture design itself is Agilent's.

## Citing

See [`CITATION.cff`](CITATION.cff). Please cite the Data Descriptor, and the archived release
(Zenodo DOI **[PENDING-ACCESSION]**) if you reuse the code or processed tables.

Repository: https://github.com/felixfan/na12878-wes-benchmark

## Third-party data

The GIAB HG001 v4.2.1 benchmark used as ground truth is from the Genome in a Bottle Consortium
(NIST/NCBI) and is **not** redistributed here — see the config template for the download location.
