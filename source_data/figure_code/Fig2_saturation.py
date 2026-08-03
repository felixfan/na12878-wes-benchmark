#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure 2 — Depth-saturation curves (three platforms, DeepVariant). Reads source_data.xlsx."""
import pandas as pd, matplotlib.pyplot as plt

d = pd.read_excel('../source_data/source_data.xlsx', sheet_name='Fig2_saturation')
d = d[d.caller == 'deepvariant']
COL = {'illumina_novaseq6000': '#4C72B0', 'bgi_dnbseq_t7': '#DD8452', 'genemind_surfseq5000': '#55A868'}
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, vt in zip(axes, ['SNV', 'INDEL']):
    for plat, g in d[d.variant_type == vt].groupby('platform'):
        g = g.sort_values('mean_on_target_depth')
        ax.errorbar(g.mean_on_target_depth, g.F1_mean, yerr=g.F1_std, marker='o', capsize=3,
                    lw=1.8, color=COL.get(plat), label=plat)
    ax.set_xscale('log'); ax.set_xticks([8, 10, 15, 20, 30, 50, 75, 100, 200, 400, 700])
    ax.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.set_title(vt); ax.set_xlabel('On-target mean depth (X)'); ax.grid(alpha=.3, which='both')
axes[0].set_ylabel('F1'); axes[0].legend(fontsize=8, loc='lower right')
fig.suptitle('Figure 2 — Depth-saturation of variant-calling accuracy (DeepVariant)')
fig.tight_layout(); fig.savefig('Fig2_saturation.png', dpi=300); print('-> Fig2_saturation.png')
