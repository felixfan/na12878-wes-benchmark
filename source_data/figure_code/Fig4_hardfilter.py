#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure 4 — GATK INDEL hard-filtering before/after (GATK4, full depth). Reads source_data.xlsx."""
import pandas as pd, numpy as np, matplotlib.pyplot as plt

d = pd.read_excel('../source_data/source_data.xlsx', sheet_name='Fig4_hardfilter')
types = ['SNV', 'INDEL']; filts = [('raw (unfiltered)', '#B0B0B0'), ('hard-filtered', '#4C72B0')]
x = np.arange(len(types)); w = 0.36
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
for ax, met in zip(axes, ['F1', 'precision']):
    for i, (fv, col) in enumerate(filts):
        vals = [d[(d.variant_type == vt) & (d['filter'] == fv)][met].values for vt in types]
        vals = [v[0] if len(v) else np.nan for v in vals]
        b = ax.bar(x + (i - .5) * w, vals, w, color=col, label=fv)
        ax.bar_label(b, fmt='%.3f', fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(types); ax.set_title(met); ax.set_ylim(.6, 1.03); ax.grid(alpha=.3, axis='y')
axes[0].legend(fontsize=9, loc='lower left')
fig.suptitle('Figure 4 — GATK INDEL hard-filtering before/after (full depth)')
fig.tight_layout(); fig.savefig('Fig4_hardfilter.png', dpi=300); print('-> Fig4_hardfilter.png')
