#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure 3 — Cross-platform comparison at nominal 30x (DeepVariant). Reads source_data.xlsx."""
import pandas as pd, numpy as np, matplotlib.pyplot as plt

d = pd.read_excel('../source_data/source_data.xlsx', sheet_name='Fig3_platform')
d = d[d.caller == 'deepvariant']
COL = {'illumina_novaseq6000': '#4C72B0', 'bgi_dnbseq_t7': '#DD8452', 'genemind_surfseq5000': '#55A868'}
plats = list(COL); types = ['SNV', 'INDEL']; metrics = ['precision', 'recall', 'F1']
x = np.arange(len(types)); w = 0.25
fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
for ax, met in zip(axes, metrics):
    for i, plat in enumerate(plats):
        vals = [d[(d.platform == plat) & (d.variant_type == vt)][met].mean() for vt in types]
        ax.bar(x + (i - 1) * w, vals, w, color=COL[plat], label=plat)
    ax.set_xticks(x); ax.set_xticklabels(types); ax.set_title(met); ax.set_ylim(.8, 1.005); ax.grid(alpha=.3, axis='y')
axes[0].legend(fontsize=8, loc='lower left')
fig.suptitle('Figure 3 — Cross-platform comparison at nominal 30x (DeepVariant)')
fig.tight_layout(); fig.savefig('Fig3_platform.png', dpi=300); print('-> Fig3_platform.png')
