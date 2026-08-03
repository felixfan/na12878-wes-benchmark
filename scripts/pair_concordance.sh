#!/usr/bin/env bash
# pair_concordance.sh -- site-by-site CPU/GPU comparison, independent of the GIAB truth set,
#                        covering every T3 matched pair.
#
#   Walks the GPU callers in the manifest (pbrun_haplotypecaller / pbrun_deepvariant), derives the
#   CPU run for the same data point (same sample x reference x depth x seed, with bwa-mem and
#   gatk4_hc/deepvariant), runs vcf_concordance.py --source raw on the two calls.norm.vcf.gz files,
#   and appends shared/concordant/discordant/cpu_only/gpu_only counts to pair_concordance.tsv.
#   (calls.norm.vcf.gz is already restricted to the V6 target, so the comparison needs no further
#   restriction to the GIAB high-confidence regions.)
#
# Usage: [DRY_RUN=1] [MANIFEST=<path>] scripts/pair_concordance.sh
# Output: $WORKDIR/plan/pair_concordance.tsv
set -euo pipefail
CONFIG="${HAPPY_CONFIG:-$(dirname "$0")/../configs/na12878.paths.sh}"
[ -f "$CONFIG" ] || { echo "missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
HERE="$(cd "$(dirname "$0")" && pwd)"
DRY_RUN="${DRY_RUN:-0}"
MAN="${MANIFEST:-$WORKDIR/plan/run_manifest.csv}"
OUT="$WORKDIR/plan/pair_concordance.tsv"
QOUT="$WORKDIR/plan/pair_quality.tsv"                # median QUAL/DP/GQ per concordance category,
                                                     # to test whether disagreements are low-quality sites
run(){ echo "+ $*"; [ "$DRY_RUN" = 1 ] || eval "$*"; }
[ -f "$MAN" ] || { echo "missing manifest: $MAN" >&2; exit 1; }
[ "$DRY_RUN" = 1 ] || rm -f "$OUT" "$QOUT"

# GPU caller -> its CPU counterpart
declare -A CPUCALL=( [pbrun_haplotypecaller]=gatk4_hc [pbrun_deepvariant]=deepvariant )

# manifest columns: run_id=1 ... caller=7 ... hardware=9; select the GPU runs
mapfile -t GPUS < <(awk -F, 'NR>1 && ($7=="pbrun_haplotypecaller"||$7=="pbrun_deepvariant"){print $1}' "$MAN" | sort -u)
echo ">>> GPU runs to pair: ${#GPUS[@]}"
n=0; ok=0
for grid in "${GPUS[@]}"; do
  n=$((n+1))
  IFS='.' read -r S REF ALN GC HW DEP SEED <<<"$grid"
  ccall="${CPUCALL[$GC]:-}"; [ -n "$ccall" ] || { echo "  skipping (unknown GPU caller): $grid"; continue; }
  crid="$S.$REF.bwa-mem.$ccall.cpu.$DEP.$SEED"        # the CPU run_id for the same data point
  gvcf="$WORKDIR/runs/$grid/calls.norm.vcf.gz"
  cvcf="$WORKDIR/runs/$crid/calls.norm.vcf.gz"
  lab="$S.$REF.$ccall.$DEP.$SEED"
  echo "[$n/${#GPUS[@]}] $lab"
  if [ "$DRY_RUN" != 1 ] && { [ ! -s "$cvcf" ] || [ ! -s "$gvcf" ]; }; then
    echo "  !! missing VCF (cpu=$([ -s "$cvcf" ] && echo ok || echo MISS)  gpu=$([ -s "$gvcf" ] && echo ok || echo MISS)) -- skipping"
    continue
  fi
  run "python3 '$HERE/vcf_concordance.py' '$cvcf' '$gvcf' --source raw --label '$lab' --tsv '$OUT' --quality --qtsv '$QOUT'"
  ok=$((ok+1))
done
echo ">>> done: $ok pairs -> $OUT"
echo "    Summarise and plot pair_concordance.tsv together with pair_quality.tsv."
