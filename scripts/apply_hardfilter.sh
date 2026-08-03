#!/usr/bin/env bash
# apply_hardfilter.sh -- retro-fit hard-filtering onto GATK-family runs that have ALREADY completed,
#                        without re-running variant calling.
#
#   runs/$RID/calls.vcf.gz  --hardfilter_vcf.sh-->  calls.filt.vcf.gz
#     --bcftools norm-->  calls.norm.vcf.gz  --run_happy.sh-->  overwrites happy/$RID.{summary,extended}.csv
#   Re-running aggregate_happy.py afterwards makes filter=PASS reflect the hard-filtered set, while
#   filter=ALL (the raw callset) is unchanged. DeepVariant runs are left alone: they filter internally.
#
# Usage:
#   scripts/apply_hardfilter.sh                 # every GATK-family run in the manifest, serially
#   scripts/apply_hardfilter.sh <run_id>        # a single run (for external parallel scheduling)
#   scripts/apply_hardfilter.sh --list          # print the GATK-family run_id list only
# Switches: DRY_RUN=1 | FORCE=1 (redo an existing filtered VCF) | THREADS=N | MANIFEST=<path>
# To parallelise when CPU is spare: scripts/apply_hardfilter.sh --list > /tmp/rids.txt
#   xargs -P8 -a /tmp/rids.txt -I{} sh -c 'scripts/apply_hardfilter.sh "$1" </dev/null' _ {}
#   NOTE: </dev/null must sit inside sh -c, redirecting the child's stdin so that docker -i cannot
#   consume it. Placing it at the end of the xargs line would redirect xargs' own stdin instead,
#   leaving it with nothing to read, so nothing would run.
set -euo pipefail
CONFIG="${HAPPY_CONFIG:-$(dirname "$0")/../configs/na12878.paths.sh}"
[ -f "$CONFIG" ] || { echo "missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
HERE="$(cd "$(dirname "$0")" && pwd)"
export PATH="${WRAPPER_BIN:-$HOME/bin}:$PATH"
DRY_RUN="${DRY_RUN:-0}"; FORCE="${FORCE:-0}"; THREADS="${THREADS:-8}"
MAN="${MANIFEST:-$WORKDIR/plan/run_manifest.csv}"   # for a local dry-run, point MANIFEST at a downloaded copy of the manifest
run(){ echo "+ $*"; [ "$DRY_RUN" = 1 ] || eval "$*"; }

is_gatk(){ case "$1" in gatk4_hc|gatk3|pbrun_haplotypecaller) return 0;; *) return 1;; esac; }

list_rids(){   # manifest columns: run_id=1 ... caller=7
  [ -f "$MAN" ] || { echo "missing manifest: $MAN" >&2; exit 1; }
  awk -F, 'NR>1 && ($7=="gatk4_hc"||$7=="gatk3"||$7=="pbrun_haplotypecaller"){print $1}' "$MAN" | sort -u
}

do_one(){
  local RID="$1" SAMPLE REF _A _C _HW _D _S FASTA OUT RAW FILT NORM HPFX
  IFS='.' read -r SAMPLE REF _A _C _HW _D _S <<<"$RID"
  is_gatk "$_C" || { echo "  skipping (not a GATK-family run): $RID"; return 0; }
  case "$REF" in
    grch38)          FASTA=$REF_GRCH38 ;;
    b37|grch37|hg19) FASTA=$REF_B37 ;;
    *) echo "  unknown reference: $REF ($RID)" >&2; return 1 ;;
  esac
  OUT="$WORKDIR/runs/$RID"; RAW="$OUT/calls.vcf.gz"; FILT="$OUT/calls.filt.vcf.gz"
  NORM="$OUT/calls.norm.vcf.gz"; HPFX="$WORKDIR/happy/$RID"
  echo "=== hardfilter $RID (ref=$REF caller=$_C) ==="
  if [ "$DRY_RUN" != 1 ] && [ ! -s "$RAW" ]; then
    echo "  !! raw VCF missing: $RAW -- skipping (this run needs calling to be re-run first)"; return 0
  fi
  if [ "$DRY_RUN" != 1 ] && [ -s "$FILT" ] && [ "$FORCE" != 1 ]; then
    echo "  hard-filtered VCF already exists (set FORCE=1 to redo), reusing $FILT"
  else
    run "DRY_RUN=$DRY_RUN THREADS=$THREADS bash '$HERE/hardfilter_vcf.sh' '$RAW' '$FASTA' '$FILT'"
  fi
  run "bcftools norm -m-any -f '$FASTA' '$FILT' -Oz -o '$NORM'"
  run "bcftools index -ft '$NORM'"
  run "HAPPY_CONFIG='$CONFIG' bash '$HERE/run_happy.sh' --query '$NORM' --ref '$REF' --out '$HPFX' --threads $THREADS"
  echo "  [ok] $RID hard-filtered and re-benchmarked -> $HPFX.summary.csv"
}

case "${1:-}" in
  --list) list_rids ;;
  "")
    mapfile -t RIDS < <(list_rids)
    echo ">>> GATK-family runs: ${#RIDS[@]} (serial; see the xargs example in the header to parallelise)"
    i=0; for RID in "${RIDS[@]}"; do i=$((i+1)); echo "[$i/${#RIDS[@]}]"; do_one "$RID"; done
    echo ">>> All done. Next: python3 scripts/aggregate_happy.py --manifest '$MAN' --happydir '$WORKDIR/happy'  (filter=PASS now means hard-filtered)"
    ;;
  *) do_one "$1" ;;
esac
