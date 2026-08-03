#!/usr/bin/env bash
# run_batch.sh <jobs.list> -- run a batch of run_one jobs, with resource monitoring and separate
#   CPU/GPU scheduling.
#   1) Starts scripts/monitor.sh in the background (sampling every MON_INT seconds) and stops it on
#      exit via a trap.
#   2) CPU jobs (run_id containing `.cpu.`) run in parallel through xargs -P $CPU_JOBS.
#   3) GPU jobs (`.gpu.`) are assigned to cards round-robin through GPU_DEV, $GPU_PER_CARD jobs per
#      card.
#
# NOTE: this is WES, far smaller than WGS, so start with conservative CPU_JOBS / GPU_PER_CARD, run
#   one batch, then tune using the monitor output (memory, GPU memory, utilisation).
# NOTE: CPU DeepVariant is the most memory-hungry step; hitting the memory wall slows everything.
# NOTE: if another process is holding GPU memory, free the cards before running the GPU jobs.
#
# Usage:
#   HAPPY_CONFIG=configs/na12878.paths.sh bash scripts/run_batch.sh $WORKDIR/plan/jobs.list
#   CPU_JOBS=8 GPU_PER_CARD=2 NGPU=2 RUN_WHAT=cpu bash scripts/run_batch.sh jobs.list   # CPU jobs only
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
JOBS="${1:?usage: run_batch.sh <jobs.list>}"; [ -f "$JOBS" ] || { echo "cannot find $JOBS" >&2; exit 1; }
CONFIG="${HAPPY_CONFIG:-$HERE/../configs/na12878.paths.sh}"; source "$CONFIG"

export THREADS="${THREADS:-16}"     # threads per run (bwa/gatk/dv/mosdepth); keep THREADS x CPU_JOBS <= core count
CPU_JOBS="${CPU_JOBS:-6}"           # concurrent CPU jobs (DeepVariant is memory-hungry, so start low)
NGPU="${NGPU:-2}"                   # number of GPUs
GPU_PER_CARD="${GPU_PER_CARD:-1}"   # concurrent jobs per card (start at 1, raise after checking GPU memory)
MON_INT="${MON_INT:-5}"             # monitor sampling interval in seconds
RUN_WHAT="${RUN_WHAT:-all}"         # all | cpu | gpu

# ---- 1) start resource monitoring (the trap guarantees it stops on exit) ----
MONCSV="${MONCSV:-$WORKDIR/plan/resource_timeseries.csv}"   # override this when running two batches
                                                            # concurrently, so the monitors do not
                                                            # interleave into one file
bash "$HERE/monitor.sh" "$MON_INT" "$MONCSV" >/dev/null 2>&1 &
MONPID=$!
trap 'kill "$MONPID" 2>/dev/null; echo ">>> monitoring stopped -> $MONCSV"' EXIT
echo ">>> resource monitoring started (every ${MON_INT}s -> $MONCSV, PID=$MONPID)"

CPU_LIST=$(grep '\.cpu\.' "$JOBS" || true)
GPU_LIST=$(grep '\.gpu\.' "$JOBS" || true)

# The BAM key excludes the caller, so all callers sharing (sample.ref.aligner.depth.seed) reuse one
# alignment. Running them in parallel would race to build the same BAM, so the batch is split into
# two phases: first build each unique BAM once (ALIGN_ONLY, keys are unique so there is no conflict),
# then run calling in parallel against the cached BAMs.
bamkeys() { awk -v c="$1" -v h="$2" -F. 'BEGIN{OFS="."}{k=$1"."$2"."$3"."$6"."$7; if(!(k in s)){s[k]=1; print $1,$2,$3,c,h,$6,$7}}'; }
sched_gpu() {  # <job list> <ALIGN_ONLY:1|0> : round-robin GPU_DEV across cards, concurrency = NGPU x GPU_PER_CARD
  # NOTE: only wait on the run_one PIDs started here. A bare `wait` or `wait -n` would also wait on
  #   the background monitor.sh, which samples forever, deadlocking between the two phases.
  local list="$1" ao="$2" i=0 r; local -a pids=()
  local cap=$(( NGPU * GPU_PER_CARD )); [ "$cap" -lt 1 ] && cap=1
  while IFS= read -r r; do [ -z "$r" ] && continue
    # < /dev/null is required: the pbrun wrapper uses `docker run -i` and wants stdin, so several
    # background jobs would compete for the terminal's stdin and the losers would fail outright.
    # (The CPU phase goes through xargs -P, whose children already get /dev/null as stdin.)
    GPU_DEV="$(( i % NGPU ))" ALIGN_ONLY="$ao" HAPPY_CONFIG="$CONFIG" bash "$HERE/run_one.sh" "$r" < /dev/null &
    pids+=("$!"); i=$((i+1))
    [ "${#pids[@]}" -ge "$cap" ] && { wait "${pids[0]}"; pids=("${pids[@]:1}"); }
  done <<< "$list"
  [ "${#pids[@]}" -gt 0 ] && wait "${pids[@]}"
}

# ---- 2) CPU: two phases ----
if [ "$RUN_WHAT" != gpu ] && [ -n "$CPU_LIST" ]; then
  A=$(printf '%s\n' "$CPU_LIST" | bamkeys gatk4_hc cpu)
  echo ">>> CPU phase A: building $(printf '%s\n' "$A"|grep -c .) BAMs (ALIGN_ONLY, -P $CPU_JOBS)"
  printf '%s\n' "$A" | HAPPY_CONFIG="$CONFIG" ALIGN_ONLY=1 xargs -P "$CPU_JOBS" -I{} bash "$HERE/run_one.sh" {}
  echo ">>> CPU phase B: calling for $(printf '%s\n' "$CPU_LIST"|grep -c .) runs (reusing BAMs, -P $CPU_JOBS)"
  printf '%s\n' "$CPU_LIST" | HAPPY_CONFIG="$CONFIG" xargs -P "$CPU_JOBS" -I{} bash "$HERE/run_one.sh" {}
fi

# ---- 3) GPU: two phases, spread across cards ----
if [ "$RUN_WHAT" != cpu ] && [ -n "$GPU_LIST" ]; then
  GA=$(printf '%s\n' "$GPU_LIST" | bamkeys pbrun_fq2bam gpu)
  echo ">>> GPU phase A: building $(printf '%s\n' "$GA"|grep -c .) fq2bam BAMs (ALIGN_ONLY, $NGPU cards)"
  sched_gpu "$GA" 1
  echo ">>> GPU phase B: calling for $(printf '%s\n' "$GPU_LIST"|grep -c .) runs ($NGPU x $GPU_PER_CARD concurrent)"
  sched_gpu "$GPU_LIST" 0
fi

echo ">>> batch complete ($JOBS)"
