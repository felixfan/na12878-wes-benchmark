#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""vcf_concordance.py -- site-by-site CPU versus GPU comparison, independent of the GIAB truth set,
with an optional quality breakdown.

For one data point (same sample x reference x depth x seed), the CPU and GPU caller outputs are
classified site by site:
  shared            called by both (identical CHROM:POS:REF:ALT)
    - concordant    ... and with an identical genotype
    - discordant    same site, different genotype
  cpu_only          called only by the CPU pipeline
  gpu_only          called only by the GPU pipeline
Results are stratified into SNP and INDEL, answering how many sites both pipelines call, how many of
those agree exactly, and how many are private to one pipeline.

--quality: additionally report the median quality of each category (concordant, discordant,
           cpu_only, gpu_only), using QUAL/DP/GQ for a raw VCF or QQ for a hap.py VCF. This tests
           whether disagreements concentrate at low-quality sites.

Input may be either:
  (a) a hap.py output VCF (with TRUTH and QUERY columns) -> the QUERY column is used
      (--source happy, the default)
  (b) a plain single-sample VCF such as calls.norm.vcf.gz -> the last sample column is used
      (--source raw)

Usage: vcf_concordance.py <cpu.vcf.gz> <gpu.vcf.gz> [--source happy|raw]
                          [--label NAME] [--tsv COUNTS.tsv] [--quality] [--qtsv QUALITY.tsv]
"""
import gzip, sys, argparse, statistics as st

QFIELDS = ['QUAL', 'DP', 'GQ', 'QQ']   # quality fields; whichever are present get used


def gt_norm(gt):
    """Normalise a genotype: drop phasing and sort the alleles, so 0/1 == 0|1 == 1/0."""
    g = gt.replace('|', '/').split('/')
    if any(a in ('.', '') for a in g):
        return None
    return '/'.join(sorted(g, key=lambda x: int(x)))


def load(path, source):
    """Return {(chrom,pos,ref,alt): (gt_norm, vtype, qual_dict)} for the variants this caller actually called."""
    calls = {}
    for line in gzip.open(path, 'rt'):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 10:
            continue
        chrom, pos, ref, alt = f[0], f[1], f[3], f[4]
        fmt = f[8].split(':')
        sample = f[10] if source == 'happy' and len(f) > 10 else f[-1]  # in hap.py output, QUERY is column 11
        d = dict(zip(fmt, sample.split(':')))
        gt = gt_norm(d.get('GT', '.'))
        if gt is None or gt == '0/0':
            continue
        bvt = d.get('BVT', '')
        if bvt == 'SNP':
            vt = 'SNV'
        elif bvt == 'INDEL':
            vt = 'INDEL'
        else:
            vt = 'SNV' if (len(ref) == 1 and len(alt) == 1) else 'INDEL'
        qual = {}
        try:
            qual['QUAL'] = float(f[5])
        except (ValueError, IndexError):
            pass
        for k in ('DP', 'GQ', 'QQ'):
            if k in d:
                try:
                    qual[k] = float(d[k])
                except ValueError:
                    pass
        calls[(chrom, pos, ref, alt)] = (gt, vt, qual)
    return calls


def compare(cpu, gpu):
    ck, gk = set(cpu), set(gpu)
    shared, conly, gonly = ck & gk, ck - gk, gk - ck
    rows = {}
    for vt in ('SNV', 'INDEL', 'ALL'):
        sel = (lambda k, d: True) if vt == 'ALL' else (lambda k, d: d[k][1] == vt)
        sh = [k for k in shared if sel(k, cpu)]
        conc = sum(1 for k in sh if cpu[k][0] == gpu[k][0])
        rows[vt] = dict(
            cpu_total=sum(1 for k in ck if sel(k, cpu)),
            gpu_total=sum(1 for k in gk if sel(k, gpu)),
            shared=len(sh), concordant=conc, discordant=len(sh) - conc,
            cpu_only=sum(1 for k in conly if sel(k, cpu)),
            gpu_only=sum(1 for k in gonly if sel(k, gpu)),
        )
    return rows


def quality_report(cpu, gpu):
    """Median quality per category. Returns (dict[cat] -> {n, field -> median}, list of usable fields)."""
    ck, gk = set(cpu), set(gpu); shared = ck & gk
    cats = {
        'concordant': [(k, cpu) for k in shared if cpu[k][0] == gpu[k][0]],
        'discordant': [(k, cpu) for k in shared if cpu[k][0] != gpu[k][0]],
        'cpu_only':   [(k, cpu) for k in ck - gk],
        'gpu_only':   [(k, gpu) for k in gk - ck],
    }
    out = {}
    avail = [f for f in QFIELDS
             if any(src[k][2].get(f) is not None for items in cats.values() for k, src in items)]
    for cat, items in cats.items():
        row = {'n': len(items)}
        for fld in QFIELDS:
            vs = [src[k][2][fld] for k, src in items if src[k][2].get(fld) is not None]
            row[fld] = st.median(vs) if vs else None
        out[cat] = row
    return out, avail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cpu'); ap.add_argument('gpu')
    ap.add_argument('--source', choices=['happy', 'raw'], default='happy')
    ap.add_argument('--label', default='')
    ap.add_argument('--tsv', help='append the counts to this TSV (for batch aggregation)')
    ap.add_argument('--quality', action='store_true', help='also report the median quality of each category')
    ap.add_argument('--qtsv', help='append the median qualities to this TSV')
    a = ap.parse_args()
    cpu = load(a.cpu, a.source); gpu = load(a.gpu, a.source)
    rows = compare(cpu, gpu)
    lab = a.label or 'CPU-vs-GPU'
    print(f"=== {lab}  (source={a.source}) ===")
    hdr = ['type', 'cpu_total', 'gpu_total', 'shared', 'concordant', 'discordant', 'cpu_only', 'gpu_only', 'GT_concord%']
    print('  '.join(f'{h:>11}' for h in hdr))
    for vt in ('SNV', 'INDEL', 'ALL'):
        r = rows[vt]
        pct = 100.0 * r['concordant'] / r['shared'] if r['shared'] else float('nan')
        vals = [vt, r['cpu_total'], r['gpu_total'], r['shared'], r['concordant'],
                r['discordant'], r['cpu_only'], r['gpu_only'], f'{pct:.3f}']
        print('  '.join(f'{str(v):>11}' for v in vals))
    if a.tsv:
        import os
        new = not os.path.exists(a.tsv)
        with open(a.tsv, 'a', encoding='utf-8') as o:
            if new:
                o.write('label\ttype\tcpu_total\tgpu_total\tshared\tconcordant\tdiscordant\tcpu_only\tgpu_only\n')
            for vt in ('SNV', 'INDEL', 'ALL'):
                r = rows[vt]
                o.write(f"{lab}\t{vt}\t{r['cpu_total']}\t{r['gpu_total']}\t{r['shared']}\t{r['concordant']}\t{r['discordant']}\t{r['cpu_only']}\t{r['gpu_only']}\n")

    if a.quality or a.qtsv:
        q, avail = quality_report(cpu, gpu)
        if avail:
            print('  --- median quality by category (are disagreements low-quality sites?) ---')
            print('  ' + '%-12s %8s' % ('category', 'n') + ''.join('%9s' % f for f in avail))
            for cat in ('concordant', 'discordant', 'cpu_only', 'gpu_only'):
                r = q[cat]
                cells = ''.join('%9s' % ('%.1f' % r[f] if r[f] is not None else 'n/a') for f in avail)
                print('  ' + '%-12s %8d' % (cat, r['n']) + cells)
        else:
            print('  (no usable quality field: a raw VCF should have QUAL/DP/GQ, a hap.py VCF has QQ)')
        if a.qtsv:
            import os
            new = not os.path.exists(a.qtsv)
            with open(a.qtsv, 'a', encoding='utf-8') as o:
                if new:
                    o.write('label\tcategory\tn\t' + '\t'.join(QFIELDS) + '\n')
                for cat in ('concordant', 'discordant', 'cpu_only', 'gpu_only'):
                    r = q[cat]
                    o.write(f"{lab}\t{cat}\t{r['n']}\t" +
                            '\t'.join('' if r[f] is None else f"{r[f]:.2f}" for f in QFIELDS) + '\n')


if __name__ == '__main__':
    main()
