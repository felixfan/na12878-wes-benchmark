#!/usr/bin/env bash
# gen_plan.sh -- enumerate the T1-T4 test matrix, producing a deterministic run_id list (jobs.list)
#                and an identity table (run_plan.csv). Touches no data, only writes the plan, so it
#                can be run anywhere.
#
# Each run_id line is passed to scripts/run_one.sh (downsample -> align -> call -> norm -> coverage
# -> hap.py). The titration works in on-target-X space (see EXPERIMENT-DESIGN.md section 2), so
# sub-full points require calibrate.sh to have been run first.
# run_id grammar (schema section 1; deterministic, 7 dot-separated fields):
#   {sample}.{reference}.{aligner}.{caller}.{hardware}.{X}x.s{seed}   (full point: ...full.s0)
#   e.g. ZY_Illumina_NA12878.grch38.pbrun_fq2bam.pbrun_deepvariant.gpu.30x.s42
#
# Usage:
#   HAPPY_CONFIG=configs/na12878.paths.sh scripts/gen_plan.sh            # everything (T1+T2+T3+T4)
#   scripts/gen_plan.sh --tests T1,T3                                    # selected blocks only
#   scripts/gen_plan.sh --out ./plan                                     # output dir (default $WORKDIR/plan)
set -euo pipefail

CONFIG="${HAPPY_CONFIG:-$(dirname "$0")/../configs/na12878.paths.sh}"
[ -f "$CONFIG" ] && source "$CONFIG" || echo "note: $CONFIG not found; generating the plan with defaults (a config is still needed to run)" >&2

TESTS="T1,T3,T4"     # first-line callers = GATK4 + DeepVariant (CPU+GPU); gatk3/bcftools/freebayes are deferred
                     # (T2 = gatk3, request explicitly with --tests T2)
OUTDIR="${WORKDIR:+$WORKDIR/plan}"
while [ $# -gt 0 ]; do
  case "$1" in
    --tests) TESTS=$2; shift 2;;
    --out)   OUTDIR=$2; shift 2;;
    -h|--help) grep '^#' "$0" | sed 's/^# \?//'; exit 0;;
    *) echo "unknown argument: $1" >&2; exit 1;;
  esac
done
OUTDIR="${OUTDIR:-./plan}"
mkdir -p "$OUTDIR"

# ---- axis definitions (see EXPERIMENT-DESIGN.md sections 2-4) ----
SAMPLES=(ZY_Illumina_NA12878 ZY_BGI_NA12878 ZY_GeneMind_NA12878)   # three same-DNA NA12878 libraries
# Depth grid (on-target X): dense at the knee (8 points at <=30X), sparse across the plateau.
# 200/300/400 fill the 150-to-full gap (measured ~204x to ~660x, the widest interval in the grid;
# SNV accuracy has plateaued there but INDEL F1 still gains ~0.007-0.009).
DEPTHS_X=(5 8 10 12 15 20 25 30 40 50 60 75 100 150 200 300 400)
# GPU depth points for T3 (full depth is added separately). Widened from {15,30,50} to seven points:
#   5 and 10 were added because CPU/GPU disagreements concentrate at low depth (median DP 10.5x at
#   discordant versus 45.0x at concordant sites), making the low end the informative stress test;
#   100 is a common clinical WES depth and 200 aligns with the new high-depth points.
#   The remaining points (40/60/75/150/300/400) are omitted: their behaviour interpolates between the
#   two extremes, and running the GPU branch over the whole grid would need ~2.55 TB more storage.
DEPTHS_T3_X=(5 10 15 30 50 100 200)
SEEDS=(42 43 44)                          # downsampling seeds: 42 (the convention for single-seed comparisons), 43, 44; full depth uses seed 0

declare -A SEEN
JOBS="$OUTDIR/jobs.list"
PLAN="$OUTDIR/run_plan.csv"
: > "$JOBS"
echo "run_id,sample,platform,reference,aligner,caller,hardware,target_depth_x,seed,test" > "$PLAN"

platform_of() { case "$1" in ZY_Illumina_*) echo illumina_novaseq6000;; ZY_BGI_*) echo bgi_dnbseq_t7;; ZY_GeneMind_*) echo genemind_surfseq5000;; esac; }

emit() { # emit <sample> <ref> <aligner> <caller> <hw> <Xtoken> <seed> <test>   (Xtoken e.g. 30x|full)
  local rid="$1.$2.$3.$4.$5.$6.s$7"
  [ -n "${SEEN[$rid]:-}" ] && return 0        # de-duplicate (T1 and T3 depths can overlap)
  SEEN[$rid]=1
  echo "$rid" >> "$JOBS"
  echo "$rid,$1,$(platform_of "$1"),$2,$3,$4,$5,${6%x},$7,$8" >> "$PLAN"
}

has_test() { [[ ",$TESTS," == *",$1,"* ]]; }

# ---- T1 saturation: grch38, 3 platforms x {gatk4_hc, deepvariant} (cpu) x (17 depths x 3 seeds + full) ----
if has_test T1; then
  for s in "${SAMPLES[@]}"; do for c in gatk4_hc deepvariant; do
    for x in "${DEPTHS_X[@]}"; do for seed in "${SEEDS[@]}"; do
      emit "$s" grch38 bwa-mem "$c" cpu "${x}x" "$seed" T1
    done; done
    emit "$s" grch38 bwa-mem "$c" cpu full 0 T1              # full-depth point (no randomness, run once)
  done; done
fi

# ---- T2 legacy caller: grch38, full depth, gatk3 x seed 0 (gatk4/dv at full depth are already in T1) ----
if has_test T2; then
  for s in "${SAMPLES[@]}"; do
    emit "$s" grch38 bwa-mem gatk3 cpu full 0 T2
  done
fi

# ---- T3 CPU/GPU concordance: grch38, {pbrun_hc, pbrun_dv} (gpu, fq2bam) over the T3 depth points ----
# (the CPU counterparts gatk4_hc/deepvariant at the same depths are already in T1)
if has_test T3; then
  for s in "${SAMPLES[@]}"; do for c in pbrun_haplotypecaller pbrun_deepvariant; do
    for x in "${DEPTHS_T3_X[@]}"; do for seed in "${SEEDS[@]}"; do
      emit "$s" grch38 pbrun_fq2bam "$c" gpu "${x}x" "$seed" T3
    done; done
    emit "$s" grch38 pbrun_fq2bam "$c" gpu full 0 T3
  done; done
fi

# ---- T4 reference build: b37 legacy (no chr prefix), full depth, {gatk4_hc, deepvariant} (cpu), seed 0 ----
if has_test T4; then
  for s in "${SAMPLES[@]}"; do for c in gatk4_hc deepvariant; do
    emit "$s" b37 bwa-mem "$c" cpu full 0 T4
  done; done
fi

N=$(wc -l < "$JOBS")
echo "generated $N runs -> $JOBS"
echo "identity table -> $PLAN"
echo "runs per block:"; awk -F, 'NR>1{c[$NF]++} END{for(t in c) printf "  %s: %d\n", t, c[t]}' "$PLAN" | sort
echo
echo "next steps:"
echo "  0) calibrate first (needed to convert sub-full depths into fractions):  scripts/calibrate.sh"
echo "  1) run each line:  HAPPY_CONFIG=$CONFIG parallel -j2 scripts/run_one.sh :::: $JOBS   # keep GPU concurrency low"
echo "     or serially:  while read r; do scripts/run_one.sh \"\$r\" </dev/null; done < $JOBS"
echo "  note: run the full-depth points first -- they upsert the calibration table, and their BAMs are reused across callers"
