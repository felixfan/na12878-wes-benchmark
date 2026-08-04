#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aggregate_happy.py -- collect the per-run hap.py results into one long-format benchmark table.

Data flow:
    run_manifest.csv         (one row per run: platform/caller/seed/versions/BED + hap.py output prefix)
        +  <happy_prefix>.extended.csv  (preferred; hap.py output, metrics only, no run identity)
        or <happy_prefix>.summary.csv   (fallback: headline rows only when extended is absent)
        =>  benchmark_long.csv          (columns defined by schema/benchmark_long.header.csv)

Usage:
    python aggregate_happy.py \
        --manifest schema/run_manifest.template.csv \
        --header   schema/benchmark_long.header.csv \
        --out      benchmark_long.csv

Design notes (column names and extraction rules verified against real hap.py v0.3 output):
    * Type=SNP/INDEL is normalised to SNV/INDEL.
    * extended.csv is a superset of summary.csv (its Subtype=* & Subset=* rows match summary
      cell for cell), so extended is parsed first and yields all three stratifications at once:
        - headline (region=all, genotype=all): uses the METRIC.* values hap.py reports;
        - het/hom  (region=all, genotype=het/hom): derived from the .het/.homalt count suffix
          columns, giving both recall and precision;
        - difficult regions (region=<Subset>, genotype=all): requires hap.py --stratification.
      The operating point is the group's summary row with QQ='*' (with --roc the same group
      also contains numeric-QQ ROC points).
    * summary.csv carries neither het/hom nor difficult regions -- it is used only when extended
      is missing (headline only). --no-extended forces this fallback.
    * rtg vcfeval is supported as an alternative to hap.py: for runs whose manifest says
      benchmark_tool=vcfeval, read <prefix>.snv.summary.txt / <prefix>.indel.summary.txt.
      See parse_vcfeval().

Requires: pandas.
"""

import argparse
import csv
import os
import sys

try:
    import pandas as pd
except ImportError:
    sys.exit("pandas is required: pip install pandas")


# ---- hap.py summary.csv columns -> our columns (headline rows, all metrics present) ----
SUMMARY_MAP = {
    "truth_total": "TRUTH.TOTAL",
    "tp":          "TRUTH.TP",
    "fn":          "TRUTH.FN",
    "query_total": "QUERY.TOTAL",
    "fp":          "QUERY.FP",
    "unk":         "QUERY.UNK",
    "fp_gt":       "FP.gt",
    "precision":   "METRIC.Precision",
    "recall":      "METRIC.Recall",
    "f1":          "METRIC.F1_Score",
    "frac_na":     "METRIC.Frac_NA",
}

# Run-identity fields copied verbatim from the manifest onto every row
IDENTITY_COLS = [
    "run_id", "sample", "platform", "reference", "aligner", "caller",
    "hardware", "target_depth_x", "mean_target_depth_x", "seed", "benchmark_tool",
    "aligner_version", "tool_version", "parabricks_version",
    "benchmark_tool_version", "giab_version", "target_bed", "run_date",
]

TYPE_MAP = {"SNP": "SNV", "INDEL": "INDEL"}

# our genotype_class -> extended.csv column suffix (het/hom live in suffix columns, not the Genotype dimension)
GENOTYPE_SUFFIX = {"het": "het", "hom": "homalt"}
# hap.py target-boundary subsets, not GIAB difficult regions -> skip
SKIP_SUBSETS = {"TS_contained", "TS_boundary"}

# extended.csv metric columns -> our columns (operating-point row, genotype=all, as reported by hap.py)
EXT_METRIC_MAP = {
    "truth_total": "TRUTH.TOTAL", "tp": "TRUTH.TP", "fn": "TRUTH.FN",
    "query_total": "QUERY.TOTAL", "fp": "QUERY.FP", "unk": "QUERY.UNK",
    "fp_gt": "FP.gt",
    "precision": "METRIC.Precision", "recall": "METRIC.Recall",
    "f1": "METRIC.F1_Score", "frac_na": "METRIC.Frac_NA",
}


def _num(v):
    """Safely convert a hap.py string to float; empty/nan returns None."""
    if v is None:
        return None
    s = str(v).strip()
    if s == "" or s.lower() in ("nan", "na", "."):
        return None
    try:
        f = float(s)
    except ValueError:
        return None
    return f if f == f else None  # filter out NaN


def parse_happy_summary(path):
    """Parse <prefix>.summary.csv -> headline rows [dict, ...] (region='all', genotype='all').

    NOTE: summary.csv holds only overall metrics (one row per Type x Filter). It has no het/hom
    split and no region stratification -- those are in extended.csv (see parse_happy_extended).
    This function therefore emits headline rows only, and is used solely as a fallback when
    extended.csv is unavailable.
    """
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    have = set(df.columns)
    rows = []
    for _, r in df.iterrows():
        vtype = TYPE_MAP.get(str(r.get("Type", "")).strip())
        filt = str(r.get("Filter", "")).strip()
        if vtype is None or filt not in ("ALL", "PASS"):
            continue
        row = {"filter": filt, "variant_type": vtype,
               "genotype_class": "all", "region_stratum": "all"}
        for col, src in SUMMARY_MAP.items():
            row[col] = _num(r[src]) if src in have else None
        rows.append(row)
    return rows


def _operating_point(rows):
    """Pick the operating point from a group of rows. In hap.py extended.csv the summary
    (no-threshold) row carries QQ='*' (non-numeric); with --roc the same group also contains
    numeric-QQ ROC points. Take the QQ='*' row (all calls included, i.e. the callset's own P/R);
    if no '*' row exists, fall back to the smallest numeric QQ.
    If a future hap.py version uses a different marker, only this function needs changing."""
    best, best_qq = None, None
    for r in rows:
        qq = _num(r.get("QQ"))               # '*' -> None
        if qq is None:                       # summary / operating-point row: use directly
            return r
        if best_qq is None or qq < best_qq:  # fallback: smallest numeric QQ
            best, best_qq = r, qq
    return best


def _suffix_metrics(r, suf, have):
    """Build a metric set from the .<suf> suffix columns (het/homalt) of an extended operating-point
    row. hap.py does not report per-genotype METRIC.* directly, but it does provide .het/.homalt
    count columns, from which both recall and precision can be derived."""
    def g(base):
        c = f"{base}.{suf}"
        return _num(r[c]) if c in have else None
    tt, tp, fn = g("TRUTH.TOTAL"), g("TRUTH.TP"), g("TRUTH.FN")
    qtp, fp, unk, qtot = g("QUERY.TP"), g("QUERY.FP"), g("QUERY.UNK"), g("QUERY.TOTAL")
    prec = qtp / (qtp + fp) if (qtp is not None and fp is not None and (qtp + fp) > 0) else None
    rec = tp / tt if (tp is not None and tt) else None
    f1 = 2 * prec * rec / (prec + rec) if (prec and rec and (prec + rec) > 0) else None
    frac_na = unk / qtot if (unk is not None and qtot) else None
    rd = lambda x: round(x, 6) if x is not None else None   # match the 6-digit precision hap.py reports
    return {"truth_total": tt, "tp": tp, "fn": fn, "query_total": qtot,
            "fp": fp, "unk": unk, "precision": rd(prec), "recall": rd(rec),
            "f1": rd(f1), "frac_na": rd(frac_na)}


def parse_happy_extended(path):
    """Parse the operating-point rows of <prefix>.extended.csv into stratified rows [dict, ...]:

      * headline           region='all', genotype='all'   -- the METRIC.* values hap.py reports
                                                             (identical to summary.csv)
      * het/hom            region='all', genotype=het/hom -- derived from the .het/.homalt suffix
                                                             columns (both recall and precision)
      * difficult regions  region=<Subset>, genotype='all' -- requires hap.py --stratification

    Only Subtype='*' (no ti/tv or INDEL-length breakdown) and Genotype='*' (genotype is handled via
    the suffix columns, not this dimension) are taken. hap.py's own target-boundary subsets
    TS_boundary / TS_contained are skipped, as they are not GIAB difficult regions.
    """
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    have = set(df.columns)
    groups = {}
    for _, r in df.iterrows():
        if str(r.get("Subtype", "*")).strip() != "*":
            continue
        if str(r.get("Genotype", "*")).strip() != "*":
            continue
        vtype = TYPE_MAP.get(str(r.get("Type", "")).strip())
        filt = str(r.get("Filter", "")).strip()
        subset = str(r.get("Subset", "*")).strip()
        if vtype is None or filt not in ("ALL", "PASS") or subset in SKIP_SUBSETS:
            continue
        groups.setdefault((vtype, filt, subset), []).append(r)

    rows = []
    for (vtype, filt, subset), grp in groups.items():
        r = _operating_point(grp)
        region = "all" if subset == "*" else subset
        base = {"filter": filt, "variant_type": vtype, "region_stratum": region}
        # genotype=all: use the metrics hap.py reports
        allrow = dict(base, genotype_class="all")
        for col, src in EXT_METRIC_MAP.items():
            allrow[col] = _num(r[src]) if src in have else None
        rows.append(allrow)
        # het/hom: only for the whole region (region=all), derived from the suffix columns
        if region == "all":
            for gt, suf in GENOTYPE_SUFFIX.items():
                m = _suffix_metrics(r, suf, have)
                if m["tp"] is not None:
                    rows.append(dict(base, genotype_class=gt, **m))
    return rows


def parse_vcfeval(path, variant_type, filt="PASS"):
    """Parse an rtg vcfeval summary.txt, taking the None (no-threshold) operating point -> one dict.

    Columns: Threshold, True-pos-baseline, True-pos-call, False-pos, False-neg,
             Precision, Sensitivity, F-measure.
    Mapping: True-pos-baseline->tp, False-neg->fn, False-pos->fp,
             Precision->precision, Sensitivity->recall, F-measure->f1.
    vcfeval evaluates only PASS records by default (--all-records includes everything), so filt
    defaults to PASS.
    Note: vcfeval's summary.txt does not split SNV/INDEL, so variant_type is supplied by the
    caller (run vcfeval once per variant type).
    """
    op = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            t = line.split()
            if len(t) < 8 or t[0] == "Threshold" or t[0].startswith("-"):
                continue
            op = t                       # keep the last data row as a fallback
            if t[0].lower() == "none":   # None = the no-threshold operating point
                break
    if op is None:
        raise ValueError(f"vcfeval summary has no data rows: {path}")
    tp, tp_call, fp, fn = _num(op[1]), _num(op[2]), _num(op[3]), _num(op[4])
    prec, sens, fm = _num(op[5]), _num(op[6]), _num(op[7])
    return {
        "filter": filt, "variant_type": variant_type,
        "genotype_class": "all", "region_stratum": "all",
        "truth_total": (tp + fn) if None not in (tp, fn) else None,
        "tp": tp, "fn": fn,
        "query_total": (tp_call + fp) if None not in (tp_call, fp) else None,
        "fp": fp, "precision": prec, "recall": sens, "f1": fm,
    }


def validate(row):
    """Return the list of warnings for a row (advisory only, never blocking). See schema section 6."""
    w = []
    for k in ("precision", "recall", "f1", "frac_na"):
        v = row.get(k)
        if v is not None and not (0.0 <= v <= 1.0):
            w.append(f"{k}={v} out of range [0,1]")
    tp, fn, tt = row.get("tp"), row.get("fn"), row.get("truth_total")
    if None not in (tp, fn, tt) and abs((tp + fn) - tt) > 0.5:
        w.append(f"tp+fn={tp + fn} != truth_total={tt}")
    p, rc, f1 = row.get("precision"), row.get("recall"), row.get("f1")
    if None not in (p, rc, f1) and (p + rc) > 0:
        f1c = 2 * p * rc / (p + rc)
        if abs(f1c - f1) > 0.01:
            w.append(f"f1={f1} disagrees with 2PR/(P+R)={f1c:.5f}")
    if row.get("hardware") == "gpu" and str(row.get("parabricks_version")) in ("NA", "", "None"):
        w.append("hardware=gpu but parabricks_version=NA")
    return w


def main():
    ap = argparse.ArgumentParser(description="Aggregate hap.py results into a long-format benchmark table")
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--manifest", required=True, help="path to run_manifest.csv")
    ap.add_argument("--header", default=os.path.join(here, "..", "schema", "benchmark_long.header.csv"),
                    help="authoritative column-order file (default: schema/benchmark_long.header.csv)")
    ap.add_argument("--out", default="benchmark_long.csv", help="output CSV")
    ap.add_argument("--happydir", default=None,
                    help="redirect the hap.py output directory: take the basename of happy_prefix and join it to this directory (useful when aggregating locally while the manifest holds server absolute paths)")
    ap.add_argument("--on-missing", choices=["skip", "error"], default="skip",
                    help="when a result file is missing: skip (default) or raise an error")
    ap.add_argument("--no-extended", action="store_true",
                    help="force use of summary.csv (headline rows only, no het/hom or difficult regions)")
    args = ap.parse_args()

    # Authoritative column order: a single source of truth, so script and schema cannot drift apart
    with open(args.header, encoding="utf-8") as fh:
        columns = next(csv.reader(fh))

    manifest = pd.read_csv(args.manifest, dtype=str, keep_default_na=False)
    # Acceptance runs and re-runs append the same run_id more than once -> de-duplicate by run_id,
    # keeping the most recent (last) entry
    n_dup = len(manifest) - manifest["run_id"].nunique()
    if n_dup:
        manifest = manifest.drop_duplicates(subset="run_id", keep="last").reset_index(drop=True)
        print(f"[dedup] dropped {n_dup} duplicate run_id rows from the manifest (kept the latest)", file=sys.stderr)
    out_rows, warns, n_runs, n_missing = [], 0, 0, 0

    for _, m in manifest.iterrows():
        prefix = m["happy_prefix"]
        if args.happydir:                                   # map server absolute path -> local happydir/<basename>
            prefix = os.path.join(args.happydir, os.path.basename(prefix))
        btool = str(m.get("benchmark_tool", "happy")).strip()
        prows = []

        if btool == "vcfeval":
            # look for <prefix>.snv.summary.txt / <prefix>.indel.summary.txt, one per variant type
            found = False
            for vt in ("SNV", "INDEL"):
                vpath = f"{prefix}.{vt.lower()}.summary.txt"
                if os.path.exists(vpath):
                    prows.append(parse_vcfeval(vpath, vt))
                    found = True
                else:
                    print(f"[missing] {vpath}", file=sys.stderr)
            if not found:
                n_missing += 1
                if args.on_missing == "error":
                    sys.exit(f"[missing] vcfeval summary for {m['run_id']}")
                continue
        else:  # hap.py
            # extended.csv is a superset of summary (headline rows match cell for cell), so prefer it:
            # it yields headline + het/hom + difficult regions at once. Fall back to summary
            # (headline only) when extended is absent.
            epath = prefix + ".extended.csv"
            spath = prefix + ".summary.csv"
            if not args.no_extended and os.path.exists(epath):
                prows = parse_happy_extended(epath)
            elif os.path.exists(spath):
                prows = parse_happy_summary(spath)
                print(f"[note] {m['run_id']} has no extended.csv: headline rows only, no het/hom or difficult regions",
                      file=sys.stderr)
            else:
                n_missing += 1
                msg = f"[missing] {epath} / {spath}"
                if args.on_missing == "error":
                    sys.exit(msg)
                print(msg, file=sys.stderr)
                continue

        identity = {c: m[c] for c in IDENTITY_COLS if c in manifest.columns}
        for prow in prows:
            row = dict(identity)
            row.update(prow)
            for wmsg in validate(row):
                print(f"[warn] {row['run_id']} {row.get('filter')}/{row.get('variant_type')}"
                      f"/{row.get('genotype_class')}/{row.get('region_stratum')}: {wmsg}",
                      file=sys.stderr)
                warns += 1
            out_rows.append(row)
        n_runs += 1

    if not out_rows:
        sys.exit("No rows were parsed -- check that happy_prefix in the manifest is correct.")

    out = pd.DataFrame(out_rows).reindex(columns=columns)
    # Use a nullable integer type for count columns so they print as 24850 rather than 24850.0,
    # leaving missing values blank
    for c in ("truth_total", "tp", "fn", "query_total", "fp", "unk", "fp_gt"):
        out[c] = pd.to_numeric(out[c], errors="coerce").round().astype("Int64")
    out.to_csv(args.out, index=False)

    # Composite-key uniqueness check (schema section 6, rule 1)
    key = ["run_id", "benchmark_tool", "filter", "variant_type", "genotype_class", "region_stratum"]
    dups = out.duplicated(subset=key).sum()
    if dups:
        print(f"[warn] {dups} rows share a duplicate composite key", file=sys.stderr)

    print(f"Done: {n_runs} runs, {len(out)} rows -> {args.out}"
          f"  ({n_missing} missing, {warns} warnings)")


if __name__ == "__main__":
    main()
