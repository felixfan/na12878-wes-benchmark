#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_benchmark.py -- generate the core figures from benchmark_long.csv.

    1) saturation     F1 versus mean on-target depth, one line per platform, seeds as error bars
    2) platform       Precision/Recall/F1 across platforms at a matched depth
    3) CPU/GPU        F1 agreement scatter for equivalent caller pairs (points on y = x mean no loss)

By default it reads schema/benchmark_long.example.csv and produces placeholder figures, which is
enough to verify the plotting pipeline. With real data:
    python plot_benchmark.py --long benchmark_long.csv --outdir figures

Requires: pandas, numpy, matplotlib.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")            # headless backend, writes files directly
import matplotlib.pyplot as plt

# font selection
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
from matplotlib.ticker import ScalarFormatter


def _logx_depth_ticks(ax):
    """Label the log x-axis at the actual depth grid points, rather than only at 10 and 100."""
    ticks = [8, 10, 15, 20, 30, 50, 75, 100, 200, 400, 700]
    ax.set_xticks(ticks)
    ax.set_xticks([], minor=True)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.tick_params(axis="x", labelsize=7.5)

PLATFORM_LABEL = {"illumina_novaseq6000": "Illumina NovaSeq6000",
                  "bgi_dnbseq_t7": "BGI DNBSEQ-T7",
                  "genemind_surfseq5000": "GeneMind SURFSeq5000"}
PLATFORM_COLOR = {"illumina_novaseq6000": "#4C72B0", "bgi_dnbseq_t7": "#DD8452",
                  "genemind_surfseq5000": "#55A868"}
PLATFORM_ORDER = ["illumina_novaseq6000", "bgi_dnbseq_t7", "genemind_surfseq5000"]
# semantically equivalent CPU <-> GPU caller pairs
GPU_PAIRS = [("gatk4_hc", "pbrun_haplotypecaller"), ("deepvariant", "pbrun_deepvariant")]
FOOT = "NA12878 · GRCh38 · GIAB HG001 v4.2.1 · hap.py vcfeval · seed 42/43/44"
# reference-build axis (T4): b37 legacy versus GRCh38
REF_LABEL = {"b37": "b37 (legacy)", "grch38": "GRCh38"}
REF_COLOR = {"b37": "#8C8C8C", "grch38": "#4C72B0"}
REF_ORDER = ["b37", "grch38"]


def headline(df, caller=None):
    """Select the headline subset: filter=PASS, genotype=all, region=all.
    Note the use of df['filter'] rather than df.filter, which is a DataFrame method."""
    m = ((df["filter"] == "PASS") & (df["genotype_class"] == "all")
         & (df["region_stratum"] == "all"))
    if caller is not None:
        m &= df["caller"] == caller
    return df[m].copy()


def _agg_seeds(sub, xcol="mean_target_depth_x", ycol="f1",
               by=("platform", "variant_type", "target_depth_x")):
    """Aggregate seeds by `by`: mean x, mean y, and std as yerr (0 for a single sample)."""
    g = sub.groupby(list(by), dropna=False)
    out = g.agg(x=(xcol, "mean"), y=(ycol, "mean"),
                yerr=(ycol, "std"), n=(ycol, "size")).reset_index()
    out["yerr"] = out["yerr"].fillna(0.0)
    return out


def fig_saturation(df, caller, outpath, reference="grch38"):
    sub = headline(df, caller)
    if reference:                       # restrict to one build, otherwise b37 full-depth points would
                                    # be mixed into the GRCh38 curve
        sub = sub[sub["reference"] == reference]
    sub = sub.dropna(subset=["f1", "mean_target_depth_x"])
    types = ["SNV", "INDEL"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, vt in zip(axes, types):
        a = _agg_seeds(sub[sub["variant_type"] == vt])
        for plat in PLATFORM_ORDER:
            grp = a[a["platform"] == plat].sort_values("x")
            if grp.empty:
                continue
            ax.errorbar(grp["x"], grp["y"], yerr=grp["yerr"], marker="o", capsize=3,
                        lw=1.8, color=PLATFORM_COLOR[plat], label=PLATFORM_LABEL[plat])
        ax.set_title(vt)
        ax.set_xlabel("On-target mean depth (X)")
        ax.set_xscale("log")
        _logx_depth_ticks(ax)
        ax.grid(alpha=.3, which="both")
    axes[0].set_ylabel("F1")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle(f"Depth saturation: F1 versus on-target depth   (caller={caller})")
    fig.text(0.99, 0.005, FOOT, ha="right", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def fig_platform(df, caller, depth, outpath, reference="grch38"):
    sub = headline(df, caller)
    if reference:
        sub = sub[sub["reference"] == reference]
    sub = sub[pd.to_numeric(sub["target_depth_x"], errors="coerce") == float(depth)]
    metrics = ["precision", "recall", "f1"]
    types = ["SNV", "INDEL"]
    plats = [p for p in PLATFORM_ORDER if p in set(sub["platform"])]
    x = np.arange(len(types))
    w = 0.35
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
    for ax, met in zip(axes, metrics):
        for i, plat in enumerate(plats):
            vals, errs = [], []
            for vt in types:
                v = pd.to_numeric(
                    sub[(sub["platform"] == plat) & (sub["variant_type"] == vt)][met],
                    errors="coerce")
                vals.append(v.mean() if len(v) else np.nan)
                errs.append(v.std() if v.notna().sum() > 1 else 0.0)
            off = (i - 0.5 * (len(plats) - 1)) * w
            ax.bar(x + off, vals, w, yerr=errs, capsize=3,
                   color=PLATFORM_COLOR[plat], label=PLATFORM_LABEL[plat])
        ax.set_xticks(x)
        ax.set_xticklabels(types)
        ax.set_title(met)
        ax.set_ylim(0.80, 1.005)
        ax.grid(alpha=.3, axis="y")
    if plats:
        axes[0].legend(fontsize=8, loc="lower left")
    mx = pd.to_numeric(sub["mean_target_depth_x"], errors="coerce").mean()
    fig.suptitle(f"Platform comparison at nominal {depth:g}x (measured ~{mx:.0f}x on-target)   (caller={caller})")
    fig.text(0.99, 0.005, FOOT, ha="right", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def fig_concordance(df, outpath):
    h = headline(df).dropna(subset=["f1"])
    keys = ["sample", "reference", "target_depth_x", "seed", "variant_type"]
    pts = []   # (pair_label, cpu_f1, gpu_f1)
    for cpu_caller, gpu_caller in GPU_PAIRS:
        c = h[h["caller"] == cpu_caller]
        g = h[h["caller"] == gpu_caller]
        if c.empty or g.empty:
            continue
        merged = c.merge(g, on=keys, suffixes=("_cpu", "_gpu"))
        for _, r in merged.iterrows():
            pts.append((f"{cpu_caller} vs {gpu_caller}",
                        float(r["f1_cpu"]), float(r["f1_gpu"])))
    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    if pts:
        for lab in sorted(set(p[0] for p in pts)):
            sel = [p for p in pts if p[0] == lab]
            ax.scatter([p[1] for p in sel], [p[2] for p in sel], s=70, label=lab, zorder=3)
        lo = min(min(p[1], p[2]) for p in pts) - 0.004
        hi = 1.001
        ax.plot([lo, hi], [lo, hi], "--", color="gray", lw=1, label="y = x (no loss)")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.legend(fontsize=8, loc="upper left")
    else:
        ax.text(0.5, 0.5, "no matched CPU/GPU pairs\n(need cpu and gpu runs at the same sample/depth/seed)",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_xlabel("CPU F1")
    ax.set_ylabel("GPU F1")
    ax.set_title("CPU vs GPU concordance (F1)")
    ax.grid(alpha=.3)
    fig.text(0.99, 0.005, FOOT, ha="right", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def fig_reference(df, outpath, callers=("gatk4_hc", "deepvariant")):
    """Reference-build comparison (T4): b37 legacy versus GRCh38 F1 at full depth.
    Layout: rows = SNV/INDEL, columns = caller; within each panel x = platform and the grouped bars
    are b37 versus grch38. Only full-depth points are used (single seed, so no error bars)."""
    h = headline(df)
    h = h[h["target_depth_x"].astype(str).str.lower() == "full"]
    types = ["SNV", "INDEL"]
    callers = [c for c in callers if c in set(h["caller"])]
    plats = [p for p in PLATFORM_ORDER if p in set(h["platform"])]
    if not callers or not plats:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "no full-depth data (needs b37/grch38 runs with target_depth_x=='full')",
                ha="center", va="center", transform=ax.transAxes)
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
        return
    x = np.arange(len(plats))
    w = 0.38
    nrow, ncol = len(types), len(callers)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.9 * nrow), squeeze=False)
    for ri, vt in enumerate(types):
        for ci, caller in enumerate(callers):
            ax = axes[ri][ci]
            for k, ref in enumerate(REF_ORDER):
                vals = []
                for plat in plats:
                    v = pd.to_numeric(
                        h[(h["platform"] == plat) & (h["caller"] == caller)
                          & (h["variant_type"] == vt) & (h["reference"] == ref)]["f1"],
                        errors="coerce")
                    vals.append(v.mean() if len(v) else np.nan)
                off = (k - 0.5) * w
                bars = ax.bar(x + off, vals, w, color=REF_COLOR[ref], label=REF_LABEL[ref])
                ax.bar_label(bars, fmt="%.3f", fontsize=7, padding=2)
            ax.set_xticks(x)
            ax.set_xticklabels([PLATFORM_LABEL[p].split()[0] for p in plats], fontsize=8)
            ax.set_title(f"{vt} · {caller}", fontsize=10)
            ax.set_ylim(0.75, 1.02)
            ax.grid(alpha=.3, axis="y")
            if ri == 0 and ci == 0:
                ax.legend(fontsize=8, loc="lower left")
    for ri in range(nrow):
        axes[ri][0].set_ylabel("F1")
    fig.suptitle("Reference build at full depth: b37 legacy versus GRCh38")
    fig.text(0.99, 0.005, "NA12878 · full-depth · GIAB HG001 v4.2.1 · hap.py vcfeval",
             ha="right", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def fig_speedup(resource_csv, outpath):
    """CPU/GPU speed-up: pair GPU and CPU runs on the same sample x reference x depth x seed and compare
    the wall-clock time of the calling step (call_sec).
    Left panel: CPU versus GPU time per pair (log-log) with iso-speed-up lines.
    Right panel: speed-up versus depth (the GPU advantage grows with depth).
    Source: resource_usage_all.csv, appended per run by run_one.sh (not benchmark_long)."""
    import csv as _csv
    PAIR = {"pbrun_haplotypecaller": "gatk4_hc", "pbrun_deepvariant": "deepvariant"}
    LAB = {"pbrun_haplotypecaller": "HaplotypeCaller (vs GATK4 CPU)", "pbrun_deepvariant": "DeepVariant (vs CPU)"}
    COL = {"pbrun_haplotypecaller": "#4C72B0", "pbrun_deepvariant": "#DD8452"}
    rows = {}
    for r in _csv.DictReader(open(resource_csv)):
        rows[r["run_id"]] = r
    cpu, gpu = {}, {}
    for r in rows.values():
        k = (r["sample"], r["reference"], r["target_depth_x"], r["seed"])
        (gpu if r["hardware"] == "gpu" else cpu)[(k, r["caller"])] = r
    P = []   # (gpu_caller, depth, cpu_sec, gpu_sec, speedup)
    for (k, gc), g in gpu.items():
        c = cpu.get((k, PAIR.get(gc)))
        if not c:
            continue
        try:
            cs, gs = float(c["call_sec"]), float(g["call_sec"])
        except (TypeError, ValueError):
            continue
        if cs > 0 and gs > 0:
            P.append((gc, k[2], cs, gs, cs / gs))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.6))
    if not P:
        ax1.text(.5, .5, "no paired timing data", ha="center", va="center", transform=ax1.transAxes)
        fig.savefig(outpath, dpi=150); plt.close(fig); return
    # left: time scatter with iso-speed-up reference lines
    allsec = [p[2] for p in P] + [p[3] for p in P]
    lo, hi = min(allsec) * 0.6, max(allsec) * 1.5
    for m, lab in [(5, "5×"), (10, "10×"), (20, "20×")]:
        ax1.plot([lo, hi], [lo / m, hi / m], "--", color="gray", lw=.8, alpha=.55)
        ax1.text(hi * 0.97, hi / m, lab, fontsize=7, color="gray", va="center", ha="right")
    for gc in COL:
        xs = [p[2] for p in P if p[0] == gc]; ys = [p[3] for p in P if p[0] == gc]
        ax1.scatter(xs, ys, s=42, color=COL[gc], label=LAB[gc], zorder=3, alpha=.85, edgecolor="white", lw=.5)
    ax1.set_xscale("log"); ax1.set_yscale("log"); ax1.set_xlim(lo, hi); ax1.set_ylim(lo / 25, hi / 3.5)
    ax1.set_xlabel("CPU calling wall-clock (s)"); ax1.set_ylabel("GPU calling wall-clock (s)")
    ax1.set_title("Per pair: CPU versus GPU time"); ax1.legend(fontsize=8, loc="lower right"); ax1.grid(alpha=.3, which="both")
    # right: speed-up versus nominal depth
    def dk(d): return (d == "full", 1e9 if d == "full" else float(d))
    for gc in COL:
        ds = sorted(set(p[1] for p in P if p[0] == gc), key=dk)
        med = [np.median([p[4] for p in P if p[0] == gc and p[1] == d]) for d in ds]
        ax2.plot(range(len(ds)), med, marker="o", color=COL[gc], lw=1.9, label=LAB[gc].split(" (")[0])
        ax2.set_xticks(range(len(ds))); ax2.set_xticklabels([str(d) + ("x" if d != "full" else "") for d in ds])
    ax2.set_xlabel("Nominal on-target depth"); ax2.set_ylabel("Speed-up (CPU / GPU calling time)")
    ax2.set_title("Speed-up grows with depth"); ax2.legend(fontsize=8, loc="upper left"); ax2.grid(alpha=.3, axis="y")
    ax2.axhline(1, color="gray", lw=.8, ls=":")
    fig.suptitle("GPU acceleration: Parabricks versus CPU (calling wall-clock, matched pairs)")
    fig.text(0.99, 0.005, "NA12878 - call_sec pairs (same sample, depth and seed) - 2x RTX PRO 6000 Blackwell",
             ha="right", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def fig_caller_compare(df, outpath, sample_kw="Illumina"):
    """Single-platform (Illumina by default) GATK4 HC versus DeepVariant saturation, one panel each for
    SNV and INDEL. Because headline rows use filter=PASS, the GATK INDEL curve is the hard-filtered
    one."""
    sub = df[(df["filter"] == "PASS") & (df["genotype_class"] == "all")
             & (df["region_stratum"] == "all") & df["sample"].str.contains(sample_kw)]
    CC = {"gatk4_hc": "#4C72B0", "deepvariant": "#DD8452"}
    CL = {"gatk4_hc": "GATK4 HC", "deepvariant": "DeepVariant"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, vt in zip(axes, ["SNV", "INDEL"]):
        s2 = sub[sub["variant_type"] == vt].dropna(subset=["f1", "mean_target_depth_x"])
        for c in ["gatk4_hc", "deepvariant"]:
            a = _agg_seeds(s2[s2["caller"] == c], by=("caller", "variant_type", "target_depth_x")).sort_values("x")
            if a.empty:
                continue
            ax.errorbar(a["x"], a["y"], yerr=a["yerr"], marker="o", capsize=3, lw=1.8,
                        color=CC[c], label=CL[c])
        ax.set_title(vt); ax.set_xlabel("On-target mean depth (X)")
        ax.set_xscale("log"); _logx_depth_ticks(ax); ax.grid(alpha=.3, which="both")
    axes[0].set_ylabel("F1"); axes[0].legend(fontsize=9, loc="lower right")
    fig.suptitle(f"GATK4 HC vs DeepVariant  ({sample_kw}, GRCh38)")
    fig.text(0.99, 0.005, FOOT, ha="right", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def fig_hardfilter(df, outpath, caller="gatk4_hc"):
    """GATK hard-filtering before and after: F1 and Precision at full depth (GRCh38, all three
    platforms), comparing raw (ALL) with hard-filtered (PASS). INDELs improve substantially through
    precision, while SNVs are unchanged because they are never filtered."""
    sub = df[(df["caller"] == caller) & (df["genotype_class"] == "all")
             & (df["region_stratum"] == "all") & (df["reference"] == "grch38")
             & (df["target_depth_x"].astype(str) == "full")]
    types = ["SNV", "INDEL"]
    FILT = [("ALL", "raw (unfiltered)", "#B0B0B0"), ("PASS", "hard-filtered", "#4C72B0")]
    x = np.arange(len(types)); w = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 5.4))
    for ax, met, mlab in zip(axes, ["f1", "precision"], ["F1", "Precision"]):
        for i, (fv, flab, col) in enumerate(FILT):
            vals = [pd.to_numeric(sub[(sub["variant_type"] == vt) & (sub["filter"] == fv)][met],
                                  errors="coerce").mean() for vt in types]
            bars = ax.bar(x + (i - 0.5) * w, vals, w, color=col, label=flab)
            ax.bar_label(bars, fmt="%.3f", fontsize=8, padding=2)
        ax.set_xticks(x); ax.set_xticklabels(types); ax.set_title(mlab)
        ax.set_ylim(0.6, 1.03); ax.grid(alpha=.3, axis="y")
    axes[0].legend(fontsize=9, loc="lower left")
    fig.suptitle(f"GATK hard-filtering, before and after ({caller}, full depth, three platforms, GRCh38)")
    fig.text(0.99, 0.005, "NA12878 - INDEL-only hard-filtering (GATK best practice); SNVs left unfiltered",
             ha="right", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Generate the core figures from benchmark_long")
    ap.add_argument("--long", default=os.path.join(here, "..", "schema", "benchmark_long.example.csv"),
                    help="benchmark_long.csv (defaults to the bundled example data)")
    ap.add_argument("--outdir", default=os.path.join(here, "..", "results", "figures"))
    ap.add_argument("--depth", default=30, type=float, help="which target_depth_x the platform-comparison figure uses")
    ap.add_argument("--resource", default=None, help="resource_usage_all.csv; when given, the CPU/GPU speed-up figure is produced")
    args = ap.parse_args()

    df = pd.read_csv(args.long)
    os.makedirs(args.outdir, exist_ok=True)

    # Produce the full figure set in one pass (fixed names, overwriting the output directory):
    # saturation and platform figures per caller, then concordance and the reference-build comparison
    jobs = []
    for c in ("gatk4_hc", "deepvariant"):
        jobs.append((f"sat3_{c}.png", lambda p, c=c: fig_saturation(df, c, p)))
        jobs.append((f"platform_{c}_{args.depth:g}x.png",
                     lambda p, c=c: fig_platform(df, c, args.depth, p)))
    jobs.append(("concordance_cpu_gpu.png", lambda p: fig_concordance(df, p)))
    jobs.append(("reference_compare.png", lambda p: fig_reference(df, p)))     # reference-build comparison
    jobs.append(("caller_compare.png", lambda p: fig_caller_compare(df, p)))   # GATK4 versus DeepVariant (Illumina)
    jobs.append(("hardfilter_gatk.png", lambda p: fig_hardfilter(df, p)))       # GATK hard-filtering before/after
    if args.resource and os.path.exists(args.resource):
        jobs.append(("speedup_cpu_gpu.png", lambda p: fig_speedup(args.resource, p)))   # CPU/GPU speed-up
    for name, fn in jobs:
        out = os.path.join(args.outdir, name)
        fn(out)
        print(f"  -> {out}")
    print(f"wrote {len(jobs)} figures (source: {os.path.basename(args.long)}, {len(df)} rows)")


if __name__ == "__main__":
    main()
