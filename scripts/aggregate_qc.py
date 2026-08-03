#!/usr/bin/env python3
"""aggregate_qc.py -- collect the reports produced by run_qc.sh into qc_metrics.csv (Table 1).

Scans each <qcdir>/<sample>.<ref>/ subdirectory and parses:
  fastp.json (raw fastq) + flagstat.txt (alignment) + dup_metrics.txt (duplicates)
  + hs_metrics.txt (on-target capture metrics)
Column definitions are in schema/qc_metrics.schema.md.

Usage: python scripts/aggregate_qc.py --qcdir <WORKDIR>/qc [--out qc_metrics.csv]
"""
import argparse, csv, json, os, re, sys

HEADER = ("sample,platform,reference,raw_read_pairs,raw_bases_gb,q30_pct,gc_pct,mean_read_len,"
          "mapped_pct,properly_paired_pct,dup_pct,on_target_pct,mean_target_cov,fold_enrichment,"
          "fold80_penalty,pct_target_10x,pct_target_20x,pct_target_30x,zero_cvg_target_pct,"
          "fastp_version,hsmetrics_tool,run_date").split(",")

PLATFORM = {"ZY_Illumina": "illumina_novaseq6000", "ZY_BGI": "bgi_dnbseq_t7", "ZY_GeneMind": "genemind_surfseq5000"}


def platform_of(sample):
    return PLATFORM.get(sample.rsplit("_", 1)[0], "unknown")   # ZY_Illumina_NA12878 → ZY_Illumina


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def parse_fastp(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        j = json.load(f)
    bf = j.get("summary", {}).get("before_filtering", {})
    tr, tb = bf.get("total_reads"), bf.get("total_bases")
    q30, gc = bf.get("q30_rate"), bf.get("gc_content")
    return {
        "raw_read_pairs": int(tr / 2) if tr else None,
        "raw_bases_gb": round(tb / 1e9, 3) if tb else None,
        "q30_pct": round(q30 * 100, 3) if q30 is not None else None,
        "gc_pct": round(gc * 100, 3) if gc is not None else None,
        "mean_read_len": bf.get("read1_mean_length"),
        "fastp_version": j.get("summary", {}).get("fastp_version") or j.get("fastp_version"),
    }


def parse_flagstat(path):
    if not os.path.exists(path):
        return {}
    txt = open(path).read()
    d = {}
    m = re.search(r"mapped \(([\d.]+)%", txt)
    if m:
        d["mapped_pct"] = float(m.group(1))
    m = re.search(r"properly paired \(([\d.]+)%", txt)
    if m:
        d["properly_paired_pct"] = float(m.group(1))
    return d


def parse_picard(path):
    """Return a header->value dict for the first data row after '## METRICS CLASS'."""
    if not os.path.exists(path):
        return {}
    lines = open(path).read().splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("## METRICS CLASS") and i + 2 < len(lines):
            hdr = lines[i + 1].split("\t")
            for j in range(i + 2, len(lines)):
                if not lines[j].strip():
                    break
                return dict(zip(hdr, lines[j].split("\t")))
    return {}


def parse_dup(path):
    p = _num(parse_picard(path).get("PERCENT_DUPLICATION"))
    return {"dup_pct": round(p * 100, 3) if p is not None else None}


def parse_hs(path):
    m = parse_picard(path)
    pct = lambda k: (round(_num(m[k]) * 100, 3) if _num(m.get(k)) is not None else None)
    raw = lambda k: (round(_num(m[k]), 3) if _num(m.get(k)) is not None else None)
    return {
        "on_target_pct": pct("PCT_SELECTED_BASES"),
        "mean_target_cov": raw("MEAN_TARGET_COVERAGE"),
        "fold_enrichment": raw("FOLD_ENRICHMENT"),
        "fold80_penalty": raw("FOLD_80_BASE_PENALTY"),
        "pct_target_10x": pct("PCT_TARGET_BASES_10X"),
        "pct_target_20x": pct("PCT_TARGET_BASES_20X"),
        "pct_target_30x": pct("PCT_TARGET_BASES_30X"),
        "zero_cvg_target_pct": pct("ZERO_CVG_TARGETS_PCT"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qcdir", required=True, help="the $WORKDIR/qc root, containing <sample>.<ref>/ subdirectories")
    ap.add_argument("--out", default="qc_metrics.csv")
    ap.add_argument("--run-date", default="")
    a = ap.parse_args()

    rows = []
    for name in sorted(os.listdir(a.qcdir)):
        d = os.path.join(a.qcdir, name)
        if not os.path.isdir(d) or "." not in name:
            continue
        sample, ref = name.rsplit(".", 1)
        row = {k: "" for k in HEADER}
        row.update(sample=sample, platform=platform_of(sample), reference=ref,
                   run_date=a.run_date, hsmetrics_tool="gatk-CollectHsMetrics-4.3.0.0")
        for parser, fn in [(parse_fastp, "fastp.json"), (parse_flagstat, "flagstat.txt"),
                           (parse_dup, "dup_metrics.txt"), (parse_hs, "hs_metrics.txt")]:
            row.update({k: v for k, v in parser(os.path.join(d, fn)).items() if v is not None})
        rows.append(row)

    if not rows:
        sys.exit(f"no <sample>.<ref>/ subdirectories found in {a.qcdir}")
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {a.out}")


if __name__ == "__main__":
    main()
