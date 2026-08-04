#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_source_data.py — build one CSV per published figure and per published table.

Every number that appears in a figure or a table of the Data Descriptor is written to its own
file under `source_data/`, derived from the processed result tables in `results/`. The figure
scripts in `figure_code/` read exactly these files, so a reader can check any published value
without re-running the pipeline.

    python scripts/make_source_data.py            # writes source_data/*.csv

Inputs (all under results/):
    benchmark_long_all.csv      accuracy, one row per run x filter x variant type x genotype class
    run_coverage_all.csv        per-run measured on-target depth
    resource_usage_all.csv      per-run wall-clock time
    pair_concordance.tsv        CPU vs GPU site-level counts
    pair_quality.tsv            CPU vs GPU DP/GQ/QUAL by concordance category
    sra_file_manifest.tsv       raw FASTQ sizes + MD5
    qc/qc_metrics.csv           sequencing QC
    qc/contamination_freemix.csv, qc/identity.pairs.tsv, qc/snv_hardfilter_test.csv

The four `s1` acceptance-replicate runs (and the two CPU/GPU pairs they form) are excluded
everywhere, giving 450 analysed runs and 132 matched pairs.
"""
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next((p for p in (HERE.parent, Path.cwd())
             if (p / "results" / "benchmark_long_all.csv").is_file()), HERE.parent)
RES = ROOT / "results"
OUT = ROOT / "source_data"
OUT.mkdir(exist_ok=True)

CPU_CALLERS = ("gatk4_hc", "deepvariant")
PAIR = {"pbrun_haplotypecaller": "gatk4_hc", "pbrun_deepvariant": "deepvariant"}
PLATFORM_LABEL = {
    "illumina_novaseq6000": "Illumina NovaSeq 6000",
    "bgi_dnbseq_t7": "BGI DNBSEQ-T7",
    "genemind_surfseq5000": "GeneMind SURFSeq 5000",
}
PLATFORM_ORDER = ["illumina_novaseq6000", "bgi_dnbseq_t7", "genemind_surfseq5000"]

written = []


def write(name, df):
    p = OUT / name
    df.to_csv(p, index=False)
    written.append((name, len(df), len(df.columns)))


# ---------------------------------------------------------------- accuracy tables
long = pd.read_csv(RES / "benchmark_long_all.csv")
long["target_depth_x"] = long["target_depth_x"].astype(str)
hl = long[(long.genotype_class == "all") & (long.region_stratum == "all")]
hl = hl[hl.seed != 1]                                     # drop acceptance replicates
PASS = hl[hl["filter"] == "PASS"]

# ---- Figure 2: depth-saturation, mean over the three seeds -----------------------
s = PASS[(PASS.reference == "grch38") & (PASS.caller.isin(CPU_CALLERS))]
write("Fig2_saturation.csv",
      s.groupby(["platform", "caller", "variant_type", "target_depth_x"])
       .agg(mean_on_target_depth=("mean_target_depth_x", "mean"),
            F1_mean=("f1", "mean"), F1_std=("f1", "std"),
            precision_mean=("precision", "mean"), recall_mean=("recall", "mean"),
            n_seeds=("f1", "size")).reset_index())

# ---- Figure 3: cross-platform comparison at nominal 30x --------------------------
s = PASS[(PASS.reference == "grch38") & (PASS.caller.isin(CPU_CALLERS)) & (PASS.target_depth_x == "30")]
write("Fig3_platform.csv",
      s.groupby(["platform", "caller", "variant_type"])
       .agg(precision=("precision", "mean"), recall=("recall", "mean"), F1=("f1", "mean"),
            mean_on_target_depth=("mean_target_depth_x", "mean"), n_seeds=("f1", "size")).reset_index())

# ---- Figure 4: GATK hard-filtering before/after, full depth ----------------------
s = hl[(hl.caller == "gatk4_hc") & (hl.reference == "grch38") & (hl.target_depth_x == "full")]
f4 = (s.groupby(["variant_type", "filter"])
        .agg(F1=("f1", "mean"), precision=("precision", "mean"), recall=("recall", "mean")).reset_index())
f4["filter"] = f4["filter"].map({"ALL": "raw (unfiltered)", "PASS": "hard-filtered"})
write("Fig4_hardfilter.csv", f4)

# ---- Figure 5: CPU vs GPU F1, one row per matched pair x variant type ------------
keys = ["sample", "reference", "target_depth_x", "seed", "variant_type"]
cpu = PASS[PASS.caller.isin(CPU_CALLERS)]
gpu = PASS[PASS.caller.isin(PAIR)].copy()
gpu["cpu_caller"] = gpu["caller"].map(PAIR)
m = gpu.merge(cpu, left_on=keys + ["cpu_caller"], right_on=keys + ["caller"], suffixes=("_gpu", "_cpu"))
f5 = (m[["sample", "reference", "cpu_caller", "target_depth_x", "seed", "variant_type", "f1_cpu", "f1_gpu"]]
      .rename(columns={"cpu_caller": "caller_pair", "f1_cpu": "CPU_F1", "f1_gpu": "GPU_F1"}))
f5["abs_delta_F1"] = (f5.CPU_F1 - f5.GPU_F1).abs()
write("Fig5_concordance.csv", f5)

# ---- Figure 6: CPU vs GPU calling wall-time and speed-up -------------------------
r = pd.read_csv(RES / "resource_usage_all.csv")
r["target_depth_x"] = r["target_depth_x"].astype(str)
num = ["call_sec", "total_sec"]
r = (r.groupby("run_id", as_index=False)          # a run executed more than once -> median
       .agg({**{c: "first" for c in r.columns if c not in num + ["run_id"]},
             **{c: "median" for c in num}}))
r = r[r.seed != 1]
rg = r[r.caller.isin(PAIR)].copy()
rg["cpu_caller"] = rg["caller"].map(PAIR)
kk = ["sample", "reference", "target_depth_x", "seed"]
mr = rg.merge(r[r.caller.isin(CPU_CALLERS)], left_on=kk + ["cpu_caller"], right_on=kk + ["caller"],
              suffixes=("_gpu", "_cpu"))
mr = mr[(mr.call_sec_cpu > 0) & (mr.call_sec_gpu > 0)]
mr["speedup"] = mr.call_sec_cpu / mr.call_sec_gpu
write("Fig6_speedup.csv",
      mr[["sample", "reference", "cpu_caller", "target_depth_x", "seed",
          "call_sec_cpu", "call_sec_gpu", "speedup"]].rename(columns={"cpu_caller": "caller_pair"}))

# ---- Figure 7: reference-build comparison, full depth ----------------------------
s = PASS[(PASS.caller.isin(CPU_CALLERS)) & (PASS.target_depth_x == "full")]
write("Fig7_reference.csv",
      s.groupby(["reference", "platform", "caller", "variant_type"])
       .agg(F1=("f1", "mean"), precision=("precision", "mean"), recall=("recall", "mean")).reset_index())

# ---------------------------------------------------------------- main tables
# ---- Table 1: sequencing QC ------------------------------------------------------
qc = pd.read_csv(RES / "qc" / "qc_metrics.csv")
mos = (hl[(hl.reference == "grch38") & (hl.target_depth_x == "full")]
       .groupby("platform").mean_target_depth_x.mean())
t1 = pd.DataFrame({
    "platform": qc.platform,
    "instrument": qc.platform.map(PLATFORM_LABEL),
    "q30_pct": qc.q30_pct.round(1),
    "on_target_pct": qc.on_target_pct.round(1),
    "duplicate_pct": qc.dup_pct.round(1),
    "mean_target_depth_mosdepth_x": qc.platform.map(mos).round(0),
    "mean_target_cov_hsmetrics_x": qc.mean_target_cov.round(0),
    "target_ge30x_pct": qc.pct_target_30x.round(1),
    "mapped_pct": qc.mapped_pct.round(2),
})
t1["_o"] = t1.platform.map({p: i for i, p in enumerate(PLATFORM_ORDER)})
write("Table1_QC.csv", t1.sort_values("_o").drop(columns="_o"))

# ---- Table 2: raw files, exact sizes and checksums -------------------------------
man = pd.read_csv(RES / "sra_file_manifest.tsv", sep="\t")
write("Table2_files.csv",
      man.rename(columns={"file": "file", "bytes": "size_bytes",
                          "gb_1e9": "size_gb_1e9", "gib": "size_gib", "md5": "md5"}))

# ---- Table 3 / Table 4 / S1: CPU vs GPU site-level comparison ---------------------
cnt = pd.read_csv(RES / "pair_concordance.tsv", sep="\t")
cnt = cnt[~cnt.label.str.endswith(".s1")]
cnt["caller"] = cnt.label.str.split(".").str[2]
t3 = (cnt[cnt.type != "ALL"].groupby(["caller", "type"])
        .agg(shared=("shared", "sum"), GT_concordant=("concordant", "sum"),
             GT_discordant=("discordant", "sum"), CPU_only=("cpu_only", "sum"),
             GPU_only=("gpu_only", "sum")).reset_index())
t3["GT_concordance_pct"] = (100 * t3.GT_concordant / t3.shared).round(3)
t3["GT_discordant_pct"] = (100 * t3.GT_discordant / t3.shared).round(3)
t3["CPU_only_pct"] = (100 * t3.CPU_only / t3.shared).round(3)
t3["GPU_only_pct"] = (100 * t3.GPU_only / t3.shared).round(3)
write("Table3_concordance.csv", t3)

qual = pd.read_csv(RES / "pair_quality.tsv", sep="\t")
qual = qual[~qual.label.str.endswith(".s1")]
write("Table4_quality.csv",
      qual.groupby("category").agg(n_pairs=("n", "size"), median_DP=("DP", "median"),
                                   median_GQ=("GQ", "median"), median_QUAL=("QUAL", "median"))
          .reindex(["concordant", "discordant", "cpu_only", "gpu_only"]).reset_index())

s1 = qual.pivot_table(index="label", columns="category", values=["n", "DP", "GQ", "QUAL"])
s1.columns = [f"{c}_{m}" for m, c in s1.columns]
order = [f"{c}_{m}" for c in ["concordant", "discordant", "cpu_only", "gpu_only"]
         for m in ["n", "DP", "GQ", "QUAL"]]
write("TableS1_pair_quality.csv",
      s1[[c for c in order if c in s1.columns]].reset_index().rename(columns={"label": "pair"}))

# ---- S2: sample identity and contamination ---------------------------------------
ident = pd.read_csv(RES / "qc" / "identity.pairs.tsv", sep="\t")
ident = ident.rename(columns={ident.columns[0]: "sample_a"})
cont = pd.read_csv(RES / "qc" / "contamination_freemix.csv")
s2 = pd.concat([
    pd.DataFrame({"metric": "somalier pairwise identity",
                  "unit": ident.sample_a + " vs " + ident.sample_b,
                  "value": ident.relatedness,
                  "detail": ("IBS0=" + ident.ibs0.astype(str)
                             + "; hom_concordance=" + ident.hom_concordance.astype(str)
                             + "; informative_sites=" + ident.n.astype(str))}),
    pd.DataFrame({"metric": "VerifyBamID2 contamination (FREEMIX)",
                  "unit": cont["sample"],          # note: cont.sample is the DataFrame method
                  "value": (100 * cont.freemix).round(3),
                  "detail": "percent; n_snp=" + cont.n_snp.astype(str) + "; " + cont.tool_source}),
], ignore_index=True)
write("TableS2_identity_contamination.csv", s2)

# ---- S3: SNV hard-filter test ----------------------------------------------------
write("TableS3_snv_hardfilter.csv", pd.read_csv(RES / "qc" / "snv_hardfilter_test.csv"))

# ---- S4: nominal versus measured on-target depth ---------------------------------
cov = hl[(hl.reference == "grch38") & (hl.aligner == "bwa-mem")].drop_duplicates("run_id")
s4 = cov.pivot_table(index="target_depth_x", columns="platform",
                     values="mean_target_depth_x", aggfunc="median")
s4 = s4.reindex([str(x) for x in [5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 75, 100, 150, 200, 300, 400]] + ["full"])
s4 = s4[[p for p in PLATFORM_ORDER if p in s4.columns]].round(1).reset_index()
s4.columns = ["nominal_depth_x"] + [f"measured_{c}_x" for c in s4.columns[1:]]
for c in s4.columns[1:]:
    s4[c.replace("measured_", "ratio_").replace("_x", "")] = [
        "" if d == "full" else round(v / float(d), 3) for d, v in zip(s4.nominal_depth_x, s4[c])]
write("TableS4_depth_mapping.csv", s4)

# ---- S5: design matrix -----------------------------------------------------------
runs = hl.drop_duplicates("run_id")
n_t1 = len(runs[(runs.reference == "grch38") & (runs.hardware == "cpu")])
n_t3 = len(runs[(runs.reference == "grch38") & (runs.hardware == "gpu")])
n_t4 = len(runs[runs.reference == "b37"])
n_depth = runs.target_depth_x.nunique()
cells = 3 * ((n_depth - 1) * 3 + 1) * 4 * 2
write("TableS5_design_matrix.csv", pd.DataFrame([
    {"block": "T1 depth titration", "runs": n_t1,
     "description": "GRCh38, CPU, 3 platforms x 2 callers x 17 sub-full depths x 3 seeds + full depth"},
    {"block": "T3 CPU-GPU pairing", "runs": n_t3,
     "description": "GRCh38, GPU (Parabricks), 3 platforms x 2 caller pairs x 7 depths x 3 seeds + full depth"},
    {"block": "T4 reference build", "runs": n_t4,
     "description": "GRCh37/b37, full depth, 3 platforms x 2 callers"},
    {"block": "Total analysed", "runs": n_t1 + n_t3 + n_t4,
     "description": f"{round(100 * (n_t1 + n_t3 + n_t4) / cells)}% of the {cells}-cell full factorial "
                    f"(3 platforms x {(n_depth - 1) * 3 + 1} depth-and-seed combinations x 4 caller/hardware "
                    f"configurations x 2 reference builds)"},
    {"block": "Excluded acceptance replicates", "runs": long[long.seed == 1].run_id.nunique(),
     "description": "Illumina 30x seed s1, CPU and GPU, both callers; recorded in the archive, excluded from analysis"},
]))

print(f"source_data -> {OUT}")
for n, r, c in written:
    print(f"  {n:38s} {r:5d} rows x {c:2d} cols")
