#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""export_source_data.py -- export the underlying data of every figure and table into one
multi-sheet xlsx. Each data figure (Fig 2-7) and each table (Table 1-4, S1) gets its own sheet;
the R and Python scripts in figure_code read this workbook to redraw the figures.
"""
import pandas as pd

BL = 'results/benchmark_long_all.csv'
RES = 'results/resource_usage_all.csv'
PC = 'results/pair_concordance.tsv'
PQ = 'results/pair_quality.tsv'
S1 = 'results/supplementary_cpu_gpu_quality.csv'
OUT = 'manuscript/source_data/source_data.xlsx'

df = pd.read_csv(BL)
# headline: genotype=all, region=all
hl = df[(df.genotype_class == 'all') & (df.region_stratum == 'all')].copy()
hl['target_depth_x'] = hl['target_depth_x'].astype(str)
PASS = hl[hl['filter'] == 'PASS']
gatk = ('gatk4_hc', 'deepvariant')

# ---- Fig 2: saturation (platform x caller x type x depth, F1 mean/std over seeds; PASS, GRCh38) ----
s = PASS[(PASS.reference == 'grch38') & (PASS.caller.isin(gatk)) & (PASS.seed != 1)]
fig2 = (s.groupby(['platform', 'caller', 'variant_type', 'target_depth_x'])
          .agg(mean_on_target_depth=('mean_target_depth_x', 'mean'),
               F1_mean=('f1', 'mean'), F1_std=('f1', 'std'),
               precision_mean=('precision', 'mean'), recall_mean=('recall', 'mean'),
               n_seeds=('f1', 'size')).reset_index())

# ---- Fig 3: platform comparison @ nominal 30x (PASS, GRCh38) ----
s = PASS[(PASS.reference == 'grch38') & (PASS.caller.isin(gatk)) & (PASS.target_depth_x == '30') & (PASS.seed != 1)]
fig3 = (s.groupby(['platform', 'caller', 'variant_type'])
          .agg(precision=('precision', 'mean'), recall=('recall', 'mean'),
               F1=('f1', 'mean'), mean_on_target_depth=('mean_target_depth_x', 'mean')).reset_index())

# ---- Fig 4: GATK hardfilter before/after (gatk4_hc, full, GRCh38, ALL vs PASS) ----
s = hl[(hl.caller == 'gatk4_hc') & (hl.reference == 'grch38') & (hl.target_depth_x == 'full')]
fig4 = (s.groupby(['variant_type', 'filter'])
          .agg(F1=('f1', 'mean'), precision=('precision', 'mean'), recall=('recall', 'mean')).reset_index())
fig4['filter'] = fig4['filter'].map({'ALL': 'raw (unfiltered)', 'PASS': 'hard-filtered'})

# ---- Fig 5: CPU vs GPU F1 pairs (60 pairs, PASS, exclude s1) ----
PAIR = {'pbrun_haplotypecaller': 'gatk4_hc', 'pbrun_deepvariant': 'deepvariant'}
keys = ['sample', 'reference', 'target_depth_x', 'seed', 'variant_type']
cpu = PASS[PASS.caller.isin(gatk)]
gpu = PASS[PASS.caller.isin(PAIR.keys())].copy()
gpu['cpu_caller'] = gpu['caller'].map(PAIR)
m = gpu.merge(cpu, left_on=keys + ['cpu_caller'], right_on=keys + ['caller'], suffixes=('_gpu', '_cpu'))
m = m[m.seed != 1]
fig5 = m[['sample', 'reference', 'cpu_caller', 'target_depth_x', 'seed', 'variant_type',
          'f1_cpu', 'f1_gpu']].rename(columns={'cpu_caller': 'caller_pair', 'f1_cpu': 'CPU_F1', 'f1_gpu': 'GPU_F1'})

# ---- Fig 6: speedup (resource_usage pairs, exclude s1) ----
r = pd.read_csv(RES)
r['run_id'] = r['run_id'].astype(str)
rg = r[r.caller.isin(PAIR.keys())].copy(); rg['cpu_caller'] = rg['caller'].map(PAIR)
rc = r[r.caller.isin(gatk)]
kk = ['sample', 'reference', 'target_depth_x', 'seed']
mr = rg.merge(rc, left_on=kk + ['cpu_caller'], right_on=kk + ['caller'], suffixes=('_gpu', '_cpu'))
mr = mr[(mr.seed != 1) & (mr.call_sec_cpu > 0) & (mr.call_sec_gpu > 0)]
mr['speedup'] = mr['call_sec_cpu'] / mr['call_sec_gpu']
fig6 = mr[['sample', 'reference', 'cpu_caller', 'target_depth_x', 'seed',
           'call_sec_cpu', 'call_sec_gpu', 'speedup']].rename(columns={'cpu_caller': 'caller_pair'})

# ---- Fig 7: reference build (full, PASS, b37 vs grch38) ----
s = PASS[(PASS.caller.isin(gatk)) & (PASS.target_depth_x == 'full')]
fig7 = (s.groupby(['reference', 'platform', 'caller', 'variant_type'])
          .agg(F1=('f1', 'mean'), precision=('precision', 'mean')).reset_index())

# ---- Table 1: QC ----
t1 = pd.DataFrame([
    ['Illumina NovaSeq 6000', '93%', '~89%', '~29%', '627x', '~97%', '≈100%'],
    ['BGI DNBSEQ-T7', '96%', '~89%', '~35%', '742x', '~97%', '99.8%'],
    ['GeneMind SURFSeq 5000', '95%', '~89%', '~33%', '687x', '~97%', '99.9%'],
], columns=['Platform', 'Q30', 'on-target', 'dup', 'mean cov', '>=30x', 'mapped'])

# ---- Table 2: MD5 ----
t2 = pd.DataFrame([
    ['ZY_Illumina_NA12878_R1.fq.gz', '26 GB', '98194a44e9fe98ecad619b49e8265bb7'],
    ['ZY_Illumina_NA12878_R2.fq.gz', '27 GB', '40da4905c53d79b0078b46e9eccb5c12'],
    ['ZY_BGI_NA12878_R1.fq.gz', '41 GB', '9e7afa31e020decc7a80cba4d7aa2965'],
    ['ZY_BGI_NA12878_R2.fq.gz', '41 GB', '1338fb0e957d3e93c3b0788610598fab'],
    ['ZY_GeneMind_NA12878_R1.fq.gz', '28 GB', 'ad92eac8fc79944e06f3c224b28b4757'],
    ['ZY_GeneMind_NA12878_R2.fq.gz', '28 GB', '48b5b57129bdcbc2c288b6b685290c75'],
], columns=['File', 'Size', 'MD5'])

# ---- Table 3 / Table 4 / S1: from pair_concordance/pair_quality (exclude s1 -> 60 pairs) ----
cnt = pd.read_csv(PC, sep='\t'); cnt = cnt[~cnt.label.str.endswith('.s1')]
cnt['caller'] = cnt.label.str.split('.').str[2]
t3 = (cnt[cnt.type != 'ALL'].groupby(['caller', 'type'])
        .agg(shared=('shared', 'sum'), GT_concordant=('concordant', 'sum'), GT_discordant=('discordant', 'sum'),
             CPU_only=('cpu_only', 'sum'), GPU_only=('gpu_only', 'sum')).reset_index())
t3['GT_concordance_pct'] = (100 * t3.GT_concordant / t3.shared).round(3)

qual = pd.read_csv(PQ, sep='\t'); qual = qual[~qual.label.str.endswith('.s1')]
t4 = (qual.groupby('category').agg(median_DP=('DP', 'median'), median_GQ=('GQ', 'median'),
                                   median_QUAL=('QUAL', 'median')).reindex(
      ['concordant', 'discordant', 'cpu_only', 'gpu_only']).reset_index())

s1 = pd.read_csv(S1); s1 = s1[~s1.pair.str.endswith('.s1')]

sheets = {
    'Fig2_saturation': fig2, 'Fig3_platform': fig3, 'Fig4_hardfilter': fig4,
    'Fig5_concordance': fig5, 'Fig6_speedup': fig6, 'Fig7_reference': fig7,
    'Table1_QC': t1, 'Table2_MD5': t2, 'Table3_concordance': t3,
    'Table4_quality': t4, 'TableS1_pair_quality': s1,
}
with pd.ExcelWriter(OUT, engine='openpyxl') as xw:
    for name, d in sheets.items():
        d.to_excel(xw, sheet_name=name, index=False)
print('saved', OUT)
for name, d in sheets.items():
    print('  %-22s %d rows x %d cols' % (name, len(d), len(d.columns)))
