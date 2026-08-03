#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""plot_fig1_design.py — Figure 1: study-design schematic for the Scientific Data descriptor.
One NA12878 DNA aliquot → three platforms → depth titration × caller × hardware × reference build
→ GIAB benchmark → Precision/Recall/F1. Output: results/figures/fig1_design.png."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False

PLAT = [("Illumina NovaSeq 6000", "627×", "#4C72B0"),
        ("BGI DNBSEQ-T7", "742×", "#DD8452"),
        ("GeneMind SURFSeq 5000", "687×", "#55A868")]


def box(ax, x, y, w, h, text, fc="#F4F4F6", ec="#444", fs=10, bold=False, tc="#111"):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.12", fc=fc, ec=ec, lw=1.4))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", color=tc, linespacing=1.35)


def arrow(ax, x1, y1, x2, y2, color="#888"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                 mutation_scale=16, lw=1.4, color=color, shrinkA=1, shrinkB=1))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "..", "results", "figures", "fig1_design.png")
    fig, ax = plt.subplots(figsize=(8.8, 11))
    ax.set_xlim(0, 10); ax.set_ylim(0, 13); ax.axis("off")

    # 1 · source DNA
    box(ax, 5, 12.2, 6.4, 1.0,
        "NA12878 (HG001) genomic DNA — single aliquot\nAgilent SureSelect V6 capture · PE150",
        fc="#ECE7F6", ec="#7E57C2", fs=10.5, bold=True)
    for i in range(3):
        arrow(ax, 5, 11.68, 1.9 + i * 3.1, 10.78)
    # 2 · three platforms
    for i, (name, depth, col) in enumerate(PLAT):
        xx = 1.9 + i * 3.1
        box(ax, xx, 10.2, 2.9, 1.0, f"{name}\nmean on-target {depth}",
            fc="white", ec=col, fs=8.6, bold=True, tc=col)
        arrow(ax, xx, 9.68, 5, 9.08)
    # 3 · depth titration
    box(ax, 5, 8.5, 8.6, 1.1,
        "Depth titration — seqtk downsampling to an on-target-X grid\n"
        "5·8·10·12·15·20·25·30·40·50·60·75·100·150× and full   (seeds 42/43/44)",
        fc="#E8F5E9", ec="#43A047", fs=9.0, bold=True)
    arrow(ax, 5, 7.93, 5, 7.28)
    # 4 · analysis matrix
    mats = [("Align + dedup", "BWA-MEM (CPU)\nParabricks fq2bam (GPU)", "#4C72B0"),
            ("Variant calling", "GATK4 HaplotypeCaller\n+ DeepVariant\n(CPU & GPU)", "#C62828"),
            ("Reference build", "GRCh38 (primary)\n+ GRCh37 / b37", "#00838F")]
    for i, (t, d, col) in enumerate(mats):
        xx = 2.0 + i * 3.0
        box(ax, xx, 6.35, 2.8, 1.5, f"\n\n{d}", fc="white", ec=col, fs=8.3)
        ax.text(xx, 6.9, t, ha="center", va="center", fontsize=9.0, fontweight="bold", color=col)
    arrow(ax, 5, 5.55, 5, 4.98)
    # 5 · benchmark
    box(ax, 5, 4.4, 7.2, 1.0,
        "Benchmark — hap.py (vcfeval) vs GIAB HG001 v4.2.1\n"
        "restricted to V6 target ∩ high-confidence regions",
        fc="#FFF3E0", ec="#EF6C00", fs=9.2, bold=True)
    arrow(ax, 5, 3.87, 5, 3.28)
    # 6 · outputs
    box(ax, 5, 2.75, 7.4, 1.0,
        "Precision · Recall · F1   for   SNV / INDEL,  het / hom\n"
        "GATK-family INDELs hard-filtered (best-practice); SNVs unfiltered",
        fc="#F4F4F6", ec="#444", fs=9.0)
    # 324-runs badge
    box(ax, 8.9, 7.55, 1.7, 0.75, "324\nbenchmark runs",
        fc="#37474F", ec="#37474F", fs=8.5, bold=True, tc="white")

    fig.suptitle("Study design", fontsize=15, fontweight="bold", y=0.975)
    fig.text(0.99, 0.008, "NA12878 three-platform deep WES · GIAB HG001 v4.2.1",
             ha="right", fontsize=7.5, color="gray")
    fig.tight_layout(rect=[0, 0.01, 1, 0.955])
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("->", os.path.abspath(out))


if __name__ == "__main__":
    main()
