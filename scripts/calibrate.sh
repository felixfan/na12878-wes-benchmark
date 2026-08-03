#!/usr/bin/env bash
# calibrate.sh -- Phase-0 calibration, a prerequisite for the depth titration.
#   See EXPERIMENT-DESIGN.md section 2.2. For each (sample, reference) it aligns the full fastq and
#   measures the mean on-target depth X_full with mosdepth, writing $WORKDIR/plan/calibration.csv.
#   run_one.sh then uses that to convert a target depth into a seqtk sampling fraction.
#
# Only combinations that are actually titrated need calibrating (grch38; b37 is run at full depth
# only and needs none). This reuses run_one.sh in ALIGN_ONLY mode (align + mosdepth, no calling),
# and the resulting BAM is reused by the later full-depth runs.
#
# Usage:
#   HAPPY_CONFIG=configs/na12878.paths.sh scripts/calibrate.sh              # all samples on grch38
#   HAPPY_CONFIG=... scripts/calibrate.sh --ref grch38 --samples ZY_Illumina_NA12878
#   DRY_RUN=1 scripts/calibrate.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

REF=grch38
SAMPLES="ZY_Illumina_NA12878 ZY_BGI_NA12878 ZY_GeneMind_NA12878"
while [ $# -gt 0 ]; do
  case "$1" in
    --ref) REF=$2; shift 2;;
    --samples) SAMPLES=$2; shift 2;;
    -h|--help) grep '^#' "$0" | sed 's/^# \?//'; exit 0;;
    *) echo "unknown argument: $1" >&2; exit 1;;
  esac
done

for s in $SAMPLES; do
  rid="$s.$REF.bwa-mem.gatk4_hc.cpu.full.s0"     # the caller field is irrelevant in ALIGN_ONLY mode; full depth means fraction=1
  echo ">>> calibrating $s ($REF): $rid"
  ALIGN_ONLY=1 bash "$HERE/run_one.sh" "$rid"
done

CONFIG="${HAPPY_CONFIG:-$HERE/../configs/na12878.paths.sh}"; source "$CONFIG"
echo
echo "calibration table: ${WORKDIR:-\$WORKDIR}/plan/calibration.csv"
[ "${DRY_RUN:-0}" = 1 ] || cat "$WORKDIR/plan/calibration.csv" 2>/dev/null || true
echo "next: scripts/gen_plan.sh to emit jobs.list, then run_one.sh per line (sub-full fractions can now be computed)"
