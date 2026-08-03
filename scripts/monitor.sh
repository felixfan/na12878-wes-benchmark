#!/usr/bin/env bash
# monitor.sh [interval_seconds] [output_csv] -- sample CPU%, memory and per-GPU utilisation/memory
#   as a time series in the background. Start it before a batch and stop it afterwards. Useful for
#   tuning concurrency, for the resource curves behind the CPU/GPU speed-up analysis, and as
#   supplementary material.
#
# Usage:
#   nohup bash scripts/monitor.sh 5 /data/test_data/felix/output/plan/resource_timeseries.csv >/dev/null 2>&1 &
#   ...  run the calibrate / run_one batch  ...
#   pkill -f monitor.sh          # stop sampling
#
# Peak values: awk -F, 'NR>1{if($2>c)c=$2; if($5>g)g=$5} END{print "peak CPU%="c" peak GPU0 memory MB="g}' <csv>
set -uo pipefail
INT="${1:-5}"; OUT="${2:-resource_timeseries.csv}"
mkdir -p "$(dirname "$OUT")" 2>/dev/null || true
NGPU=$(nvidia-smi --query-gpu=count --format=csv,noheader 2>/dev/null | head -1 | tr -dc '0-9'); NGPU=${NGPU:-0}

hdr="timestamp,cpu_pct,mem_used_gb,mem_total_gb"
for ((g=0; g<NGPU; g++)); do hdr+=",gpu${g}_util_pct,gpu${g}_mem_used_mb,gpu${g}_mem_total_mb"; done
[ -f "$OUT" ] || echo "$hdr" > "$OUT"
echo "monitor: sampling every ${INT}s -> $OUT (GPUs=$NGPU); stop with: pkill -f monitor.sh" >&2

p_idle=0; p_tot=0
while true; do
  # CPU%: difference between two /proc/stat samples (the first line is the average since boot and is
  # ignored; subsequent lines are interval utilisation)
  read -r _ u n s idle iow irq soft _ < /proc/stat
  tot=$((u+n+s+idle+iow+irq+soft)); dt=$((tot-p_tot)); di=$((idle-p_idle)); p_tot=$tot; p_idle=$idle
  cpu=$([ "$dt" -gt 0 ] && echo $(( (100*(dt-di))/dt )) || echo 0)
  # memory in GB
  read -r mu mt < <(free -g | awk '/^Mem:/{print $3, $2}')
  line="$(date +%FT%T),$cpu,${mu:-NA},${mt:-NA}"
  # per-GPU utilisation and memory (--nounits strips units; paste joins the cards onto one line)
  if [ "$NGPU" -gt 0 ]; then
    g=$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null | tr -d ' ' | paste -sd,)
    line="$line,$g"
  fi
  echo "$line" >> "$OUT"
  sleep "$INT"
done
