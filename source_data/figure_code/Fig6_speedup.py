#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure 6 — CPU vs GPU acceleration (calling wall-time, 60 pairs). Reads source_data.xlsx."""
import pandas as pd, numpy as np, matplotlib.pyplot as plt

d = pd.read_excel('../source_data/source_data.xlsx', sheet_name='Fig6_speedup')
LAB = {'gatk4_hc': 'HaplotypeCaller', 'deepvariant': 'DeepVariant'}
COL = {'gatk4_hc': '#4C72B0', 'deepvariant': '#DD8452'}
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8))
lo, hi = d[['call_sec_cpu', 'call_sec_gpu']].min().min() * .6, d.call_sec_cpu.max() * 1.5
for m, lab in [(5, '5x'), (10, '10x'), (20, '20x')]:
    ax1.plot([lo, hi], [lo / m, hi / m], '--', color='gray', lw=.8, alpha=.5)
for c, g in d.groupby('caller_pair'):
    ax1.scatter(g.call_sec_cpu, g.call_sec_gpu, s=42, color=COL[c], label=LAB[c], alpha=.85)
ax1.set_xscale('log'); ax1.set_yscale('log'); ax1.set_xlabel('CPU calling (s)'); ax1.set_ylabel('GPU calling (s)')
ax1.set_title('per-pair time'); ax1.legend(fontsize=8, loc='lower right'); ax1.grid(alpha=.3, which='both')
def dk(x): return 1e9 if str(x) == 'full' else float(x)
for c, g in d.groupby('caller_pair'):
    med = g.groupby('target_depth_x').speedup.median()
    ds = sorted(med.index, key=dk)
    ax2.plot(range(len(ds)), [med[x] for x in ds], marker='o', color=COL[c], label=LAB[c])
    ax2.set_xticks(range(len(ds))); ax2.set_xticklabels([str(x) for x in ds])
ax2.set_xlabel('nominal depth'); ax2.set_ylabel('speedup (CPU/GPU)'); ax2.set_title('speedup vs depth')
ax2.legend(fontsize=8); ax2.grid(alpha=.3, axis='y')
fig.suptitle('Figure 6 — CPU vs GPU acceleration (calling step)')
fig.tight_layout(); fig.savefig('Fig6_speedup.png', dpi=300); print('-> Fig6_speedup.png')
