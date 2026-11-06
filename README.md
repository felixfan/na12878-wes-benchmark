# NA12878 whole exome sequencing benchmark

Analysis code and processed results for NA12878/HG001 sequenced on Illumina NovaSeq 6000, BGI DNBSEQ-T7 and GeneMind SURFSeq 5000.

## Version 1.1.0

<!-- version-doi-start -->
**Current release: [v1.1.0 on Zenodo](https://doi.org/10.5281/zenodo.23984550) (DOI: 10.5281/zenodo.23984550).** The same release is available from [GitHub Releases](https://github.com/felixfan/na12878-wes-benchmark/releases/tag/v1.1.0).
<!-- version-doi-end -->

This release replaces the 66 formal GPU HaplotypeCaller runs with results from the 2026-11-05 rerun. Both GPU BQSR steps now use the same V6 capture target BED as the CPU workflow, with zero additional interval padding. The archive contains 454 run identifiers. The analyses use 450 runs and 132 matched CPU and GPU pairs. Four acceptance runs with seed s1 are retained as historical records and excluded from analyses.

The experimental run identifiers are stable across releases. Identify a result using both the release version and its run_id when comparing versions.

Raw reads are publicly available under [BioProject PRJNA1506792](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1506792), [SRA study SRP723789](https://www.ncbi.nlm.nih.gov/sra/?term=SRP723789) and [BioSample SAMN62149305](https://www.ncbi.nlm.nih.gov/biosample/SAMN62149305). The platform runs are [SRR39986369](https://www.ncbi.nlm.nih.gov/sra/?term=SRR39986369) for Illumina, [SRR39986368](https://www.ncbi.nlm.nih.gov/sra/?term=SRR39986368) for BGI and [SRR39986367](https://www.ncbi.nlm.nih.gov/sra/?term=SRR39986367) for GeneMind.

## Contents

| Path | Contents |
|---|---|
| results/benchmark_long_all.csv | 5,448 benchmark rows covering 454 archived run identifiers |
| results/plan/run_manifest.csv | Run identity, versions and benchmark output prefixes |
| results/run_coverage_all.csv | Measured on-target coverage for each archived run |
| results/resource_usage_all.csv | Calling and total wall-clock time for each archived run |
| results/pair_concordance.tsv | 134 archived CPU and GPU pairs, with three variant-type rows per pair |
| results/pair_quality.tsv | Four site categories per archived pair |
| results/revisions/v1.1.0/ | Affected run list, interval protocol, stage timings, previous affected records and old/new metric comparison |
| results/qc/ | Sequencing QC, identity, contamination and SNV filtering assessment |
| source_data/ | One CSV per figure or table, 15 files in total |
| figure_code/ | R scripts for Figures 2 to 7 |
| figures/ | 300 dpi PNG and vector PDF figures, plus the editable Figure 1 schematic |
| reference/ | SureSelect V6 Covered target BED files for GRCh38 and GRCh37/b37 |
| scripts/ | Pipeline, aggregation and source-data generation scripts |
| schema/ | Definitions of the benchmark and companion tables |
| SHA256SUMS | SHA-256 inventory of the package files |

The two archived acceptance pairs are excluded when generating the manuscript source data. Table4 reports 132 total matched pairs. Quality medians omit pairs with no sites in the relevant category. The discordant category has 91 nonempty pairs in this release.

## Results and interpretation

HaplotypeCaller genotype concordance is 99.999% for SNVs and 99.993% for INDELs. Its maximum absolute CPU/GPU F1 difference is 0.000340. The overall maximum across both caller pairs remains 0.004570. The median calling-time speedup is approximately 9.9 for HaplotypeCaller and 8.4 for DeepVariant.

The revised HaplotypeCaller GPU timings were measured in the November rerun. CPU and DeepVariant timings come from the earlier runs. Figure 6 reports the calling step alone. Total wall-clock time also includes processing and benchmarking, and depends on which upstream files were cached.

These are technical validation results for one reference sample. The libraries were prepared from the same DNA stock at different times. Any platform interpretation should account for the library preparation and storage history described in the manuscript.

## Reproduce tables and figures

From the repository root, with pandas installed, run

~~~bash
python scripts/make_source_data.py
~~~

This generates the 15 CSV files in source_data from the processed result tables. Regenerate the figures using

~~~bash
for script in figure_code/*.R; do Rscript "$script"; done
~~~

The figure scripts read source_data in this public package. The author's local manuscript folder uses source_data_v2. Figure 1 is an editable schematic in figures/fig1.pptx and has no numerical source-data CSV.

The November local rendering used Python 3.13, pandas 3.0.3, R 4.5.0, ggplot2 4.0.3, dplyr 1.2.1 and patchwork 1.3.2. The R installation supports Cairo PDF output.

## Reproduce the analysis

Use a fresh WORKDIR for version 1.1.0. Configure input files and tool wrappers before running the pipeline.

~~~bash
cp configs/na12878.paths.template.sh configs/na12878.paths.sh
~~~

Fill in the resource paths and a new output directory. Every input and output location must be visible inside the Docker containers. For example, set MOUNTS to include your home resources, input data and output filesystem when running scripts/install_wrappers.sh.

~~~bash
BIN="$HOME/bin-na12878-v110" MOUNTS="/home/USER /path/to/input /path/to/output" bash scripts/install_wrappers.sh
export WRAPPER_BIN="$HOME/bin-na12878-v110"
export PATH="$WRAPPER_BIN:$PATH"
export HAPPY_CONFIG="$PWD/configs/na12878.paths.sh"
bash scripts/calibrate.sh
bash scripts/gen_plan.sh
~~~

Review the generated plan before submitting jobs. It contains T1 with 312 CPU/GRCh38 runs, T3 with 132 GPU/GRCh38 runs and T4 with six CPU/b37 full-depth runs.

~~~bash
while read -r run_id; do
  THREADS=16 bash scripts/run_one.sh "$run_id" </dev/null
done < /path/to/WORKDIR/plan/jobs.list
~~~

GPU selection can be supplied using GPU_DEV. The default batch scheduler uses GPU indices starting at zero. Select resources appropriate to your server and avoid assigning jobs to occupied cards.

The driver records the target BED checksum, protocol and Parabricks version in a per-run protocol.txt file for GPU GATK runs. Its completion marker controls reuse of the interval-restricted recalibration BAM. Existing output files from an unmarked or incompatible protocol cause a failure with guidance to select a fresh WORKDIR.

The markers record the processing configuration. They do not replace review of the run logs and output files.

Aggregate completed hap.py outputs using

~~~bash
python scripts/aggregate_happy.py --manifest /path/to/WORKDIR/plan/run_manifest.csv --happydir /path/to/WORKDIR/happy --out /path/to/benchmark_long.csv
~~~

## Tool versions

BWA-MEM 0.7.15, GATK 4.3.0.0, DeepVariant 1.9.0, NVIDIA Parabricks 4.7.0, samtools and bcftools 1.23.1, mosdepth 0.3.3, seqtk 1.5, hap.py 0.3.15, rtg-tools 3.12.1, fastp 1.3.3 and MultiQC 1.35.

GATK-family calls use INDEL hard filtering. SNV records are exempted by the type guard. PASS selects records passing the applicable filters, and ALL includes records regardless of FILTER status. The benchmark evaluates regions shared by the V6 targets and the GIAB HG001 v4.2.1 high-confidence regions. The direct site-level CPU/GPU comparison uses the normalized called variants across the calling target regions.

## Citing this release

For version 1.1.0, cite the [Zenodo archive (DOI: 10.5281/zenodo.23984550)](https://doi.org/10.5281/zenodo.23984550). Machine-readable citation metadata are available in [`CITATION.cff`](CITATION.cff).

The [concept DOI](https://doi.org/10.5281/zenodo.21792674) identifies the resource across versions.

## Previous releases

The [historical v1.0.1 snapshot](https://doi.org/10.5281/zenodo.21792675) contains the earlier GPU BQSR results. The updated results described above are available in v1.1.0.

## Licenses

Project-authored analysis and plotting code is licensed under [MIT](LICENSE). Project-generated processed data, figures and documentation are licensed under [CC BY 4.0](LICENSE-DATA). Third-party materials retain any applicable original rights and licensing terms, as described in [LICENSE-DATA](LICENSE-DATA). The GIAB truth VCFs and raw reads are obtained from their original repositories and are not included in this repository.
