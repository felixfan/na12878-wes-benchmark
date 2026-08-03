#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure 7 — Reference-build comparison (b37 vs GRCh38, full depth). Reads source_data.xlsx."""
import pandas as pd, numpy as np, matplotlib.pyplot as plt

d = pd.read_excel('../source_data/source_data.xlsx', sheet_name='Fig7_reference')
callers = ['gatk4_hc', 'deepvariant']; types = ['SNV', 'INDEL']
REF = [('b37', '#8C8C8C'), ('grch38', '#4C72B0')]
plats = sorted(d.platform.unique()); x = np.arange(len(plats)); w = 0.38
fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.8))
for ri, vt in enumerate(types):
    for ci, caller in enumerate(callers):
        ax = axes[ri][ci]
        for k, (ref, col) in enumerate(REF):
            vals = [d[(d.platform == p) & (d.caller == caller) & (d.variant_type == vt) & (d.reference == ref)].F1.mean() for p in plats]
            b = ax.bar(x + (k - .5) * w, vals, w, color=col, label=ref)
            ax.bar_label(b, fmt='%.3f', fontsize=7)
        ax.set_xticks(x); ax.set_xticklabels([p.split('_')[0] for p in plats], fontsize=8)
        ax.set_title('%s . %s' % (vt, caller), fontsize=10); ax.set_ylim(.75, 1.02); ax.grid(alpha=.3, axis='y')
        if ri == 0 and ci == 0: ax.legend(fontsize=8, loc='lower left')
fig.suptitle('Figure 7 — Reference-build comparison (full depth)')
fig.tight_layout(); fig.savefig('Fig7_reference.png', dpi=300); print('-> Fig7_reference.png')
