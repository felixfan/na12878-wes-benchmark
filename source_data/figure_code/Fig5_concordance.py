#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure 5 — CPU vs GPU F1 concordance (60 matched pairs). Reads source_data.xlsx."""
import pandas as pd, matplotlib.pyplot as plt

d = pd.read_excel('../source_data/source_data.xlsx', sheet_name='Fig5_concordance')
fig, ax = plt.subplots(figsize=(5.4, 5.2))
for pair, g in d.groupby('caller_pair'):
    ax.scatter(g.CPU_F1, g.GPU_F1, s=70, label=pair, zorder=3)
lo = min(d.CPU_F1.min(), d.GPU_F1.min()) - .004
ax.plot([lo, 1.001], [lo, 1.001], '--', color='gray', lw=1, label='y = x')
ax.set_xlim(lo, 1.001); ax.set_ylim(lo, 1.001); ax.set_aspect('equal')
ax.set_xlabel('CPU F1'); ax.set_ylabel('GPU F1'); ax.set_title('Figure 5 — CPU vs GPU concordance (F1)')
ax.legend(fontsize=8, loc='upper left'); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig('Fig5_concordance.png', dpi=300); print('-> Fig5_concordance.png')
