#!/usr/bin/env bash
# run_happy.sh -- benchmark one query VCF with hap.py, producing <out>.summary.csv and
#                 <out>.extended.csv, which feed directly into aggregate_happy.py.
#
# Restricted to the V6 exome target (-T) intersected with the GIAB high-confidence regions (-f),
# using the vcfeval engine and GA4GH --stratification. The output columns were verified against real
# hap.py output (see schema/benchmark_long.schema.md section 4).
#
# Usage:
#   HAPPY_CONFIG=configs/na12878.paths.sh \
#   scripts/run_happy.sh --query CALLS.vcf.gz --ref grch38 --out $WORKDIR/happy/<run_id>
#
# Afterwards register the run in run_manifest.csv (happy_prefix=<out>) and run aggregate_happy.py.
set -euo pipefail

# ---- configuration (data and resource paths) ----
CONFIG="${HAPPY_CONFIG:-$(dirname "$0")/../configs/na12878.paths.sh}"
[ -f "$CONFIG" ] || { echo "missing config: $CONFIG (copy configs/na12878.paths.template.sh and fill in the paths)"; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

usage() {
  cat <<'EOF'
Usage: run_happy.sh --query VCF --ref {grch38|b37} --out PREFIX [--threads N] [--engine vcfeval|xcmp]
  --query    VCF to evaluate (bgzip + tabix)
  --ref      reference build; selects the GIAB truth set, BED, stratification and SDF (see config)
  --out      output prefix; produces <PREFIX>.summary.csv and <PREFIX>.extended.csv
  --threads  thread count (default 8)
  --engine   comparison engine (default vcfeval; xcmp is the alternative)
EOF
}

THREADS=8; ENGINE=vcfeval; QUERY=""; REF=""; OUT=""
while [ $# -gt 0 ]; do
  case "$1" in
    --query)   QUERY=$2; shift 2;;
    --ref)     REF=$2;   shift 2;;
    --out)     OUT=$2;   shift 2;;
    --threads) THREADS=$2; shift 2;;
    --engine)  ENGINE=$2;  shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "unknown argument: $1" >&2; usage; exit 1;;
  esac
done
[ -n "$QUERY" ] && [ -n "$REF" ] && [ -n "$OUT" ] || { usage; exit 1; }

# ---- select truth set / regions / reference / SDF according to the build ----
case "$REF" in
  grch38)          TRUTH=$GIAB_VCF_GRCH38; CONF=$GIAB_BED_GRCH38; TARGET=$V6_BED_GRCH38
                   STRAT=$STRAT_GRCH38; FASTA=$REF_GRCH38; SDF=${SDF_GRCH38:-} ;;
  b37|grch37|hg19) TRUTH=$GIAB_VCF_B37;   CONF=$GIAB_BED_B37;   TARGET=$V6_BED_B37
                   STRAT=$STRAT_B37;    FASTA=$REF_B37;    SDF=${SDF_B37:-} ;;
  *) echo "unknown ref: $REF (use grch38 or b37)" >&2; exit 1;;
esac

# ---- existence checks (fail early, fail clearly) ----
for v in QUERY TRUTH CONF TARGET FASTA; do
  p=${!v}
  [ -n "$p" ] || { echo "config is missing $v" >&2; exit 1; }
  [ -e "$p" ] || { echo "cannot find $v: $p" >&2; exit 1; }
done

mkdir -p "$(dirname "$OUT")"
export HGREF="$FASTA"   # silences the hap.py "No reference file found at default locations" warning (-r is passed explicitly)

# ---- hap.py ----
# -f high-confidence regions (calls outside are counted as UNK); -T the V6 exome target;
# --stratification difficult regions; the vcfeval engine needs --engine-vcfeval-template pointing at
# the RTG SDF for this build.
hap.py "$TRUTH" "$QUERY" \
  -r "$FASTA" \
  -f "$CONF" \
  -T "$TARGET" \
  ${STRAT:+--stratification "$STRAT"} \
  --engine "$ENGINE" \
  ${SDF:+--engine-vcfeval-template "$SDF"} \
  --threads "$THREADS" \
  -o "$OUT"

echo "done: ${OUT}.summary.csv / ${OUT}.extended.csv"
echo "next: register this run in run_manifest.csv (happy_prefix=${OUT}), then run aggregate_happy.py"
