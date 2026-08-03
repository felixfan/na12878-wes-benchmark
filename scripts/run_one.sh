#!/usr/bin/env bash
# run_one.sh -- execute *one* end-to-end benchmark run, fully determined by its run_id.
#
#   parse run_id -> seqtk downsampling (to the target on-target X) -> align + mark duplicates
#     (the BAM is cached and reused across callers) -> mosdepth coverage -> calling -> bcftools norm
#     -> run_happy.sh -> append rows to run_manifest.csv + run_coverage.csv
#
# The titration is done in on-target-X space (see EXPERIMENT-DESIGN.md section 2):
#   fraction = target_X / X_full, where X_full is measured by Phase-0 calibration and stored in
#   $WORKDIR/plan/calibration.csv. The full-depth point (fraction=1) needs no calibration, and the
#   mosdepth result of full + bwa-mem is upserted back into calibration.csv.
#
# Usage:
#   HAPPY_CONFIG=configs/na12878.paths.sh scripts/run_one.sh <run_id>
#   DRY_RUN=1  scripts/run_one.sh <run_id>     # print commands only
#   ALIGN_ONLY=1 scripts/run_one.sh <run_id>   # stop after mosdepth (calibration), no calling
#   RUN_HAPPY=0 scripts/run_one.sh <run_id>    # stop after norm, no benchmarking
#
# run_id format (7 dot-separated fields): {sample}.{ref}.{aligner}.{caller}.{hw}.{X}x.s{seed}
#   sub-full: ...cpu.30x.s42   |   full depth: ...cpu.full.s0
set -euo pipefail

# Detach inherited stdin: this script reads nothing from stdin, but child processes do
# (docker run -i, pbrun, ...). If the caller is `while read r; do run_one.sh "$r"; done < jobs.list`,
# a child would consume the remaining lines of jobs.list and the loop would stop after the first
# item. Pipes inside this script (seqtk etc.) get stdin from '|' and are unaffected.
exec 0</dev/null

CONFIG="${HAPPY_CONFIG:-$(dirname "$0")/../configs/na12878.paths.sh}"
[ -f "$CONFIG" ] || { echo "missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
HERE="$(cd "$(dirname "$0")" && pwd)"
export PATH="${WRAPPER_BIN:-$HOME/bin}:$PATH"   # the tool wrappers (seqtk/bwa/gatk/mosdepth/...) live in ~/bin; a non-interactive subshell does not
# read .bashrc, so PATH is extended explicitly -- otherwise "seqtk: command not found"

RID="${1:-}"; [ -n "$RID" ] || { echo "usage: run_one.sh <run_id>" >&2; exit 1; }
# Strip CR/LF and surrounding whitespace: if jobs.list was written on Windows (CRLF), `while read r`
# carries a trailing '\r' into the run_id, which corrupts the @RG line so MarkDuplicates fails with
# "@RG line missing SM tag" (and garbles terminal output).
RID="$(printf '%s' "$RID" | tr -d '\r\n' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
DRY_RUN="${DRY_RUN:-0}"; RUN_HAPPY="${RUN_HAPPY:-1}"; ALIGN_ONLY="${ALIGN_ONLY:-0}"
# Threads: default to using the machine, leaving 4 cores for the system. A fixed default of 8 badly
# under-uses the machine at the high-depth points (200-400X), where bwa aligns tens of GB of fastq.
# NOTE: when running several windows concurrently, set THREADS explicitly -- otherwise every window
# claims the whole machine and oversubscription makes things slower. For two CPU windows:
#   THREADS=$(( $(nproc)/2 - 2 ))
# None of the threaded steps changes results (bwa -t, samtools -@, mosdepth -t, pigz -p,
# hap.py --threads, and DeepVariant --num_shards which shards by region), so runs stay comparable.
if [ -z "${THREADS:-}" ]; then
  NCPU=$(nproc 2>/dev/null || echo 8)
  THREADS=$(( NCPU > 12 ? NCPU - 4 : NCPU ))
fi
# GATK (MarkDuplicates/BQSR) benefits from a larger heap: the default initial heap is only ~2.4G,
# grows during the run and spills to disk. This is a pure performance knob; results are unaffected.
# Divide by the number of concurrent windows when running several.
GATK_JAVA_OPTS="${GATK_JAVA_OPTS:--Xms4G -Xmx${GATK_MEM:-16}G}"
run() { echo "+ $*"; [ "$DRY_RUN" = 1 ] || eval "$*"; }
have() { [ "$DRY_RUN" != 1 ] && [ -s "$1" ]; }      # output exists (idempotent skip); always false in dry-run
# Decompression/compression: use pigz when available (downsampling reads the full ~600-740X fastq
# once, and single-threaded gzip decompression dominates); otherwise fall back to gzip. The
# decompressed byte stream is identical, so the seqtk subsample is unchanged -- this is purely faster.
DECOMP="gzip -dc"; COMPRESS="gzip"
command -v pigz >/dev/null 2>&1 && { DECOMP="pigz -dc"; COMPRESS="pigz -p ${THREADS}"; }

# ---- parse the run_id ----
IFS='.' read -r SAMPLE REF ALIGNER CALLER HW DEPTHTOK SEEDS <<<"$RID"
[ -n "${SEEDS:-}" ] || { echo "run_id has too few fields (7 required): $RID" >&2; exit 1; }
SEED="${SEEDS#s}"; _T0=$SECONDS   # start of wall-clock timing
if [ "$DEPTHTOK" = full ]; then TARGETX=full; else TARGETX="${DEPTHTOK%x}"; fi
case "$SAMPLE" in
  ZY_Illumina_*) PLATFORM=illumina_novaseq6000; PL=ILLUMINA; FQ_R1=$FQ_ILLUMINA_R1; FQ_R2=$FQ_ILLUMINA_R2 ;;
  ZY_BGI_*)   PLATFORM=bgi_dnbseq_t7;        PL=DNBSEQ;   FQ_R1=$FQ_BGI_R1;      FQ_R2=$FQ_BGI_R2 ;;
  ZY_GeneMind_*)       PLATFORM=genemind_surfseq5000; PL=OTHER;    FQ_R1=$FQ_GM_R1;       FQ_R2=$FQ_GM_R2 ;;
  *) echo "unknown sample: $SAMPLE" >&2; exit 1;;
esac
case "$REF" in
  grch38)          FASTA=$REF_GRCH38; TARGET=$V6_BED_GRCH38 ;;
  b37|grch37|hg19) FASTA=$REF_B37;    TARGET=$V6_BED_B37 ;;   # legacy axis = b37 (no chr prefix)
  *) echo "unknown reference: $REF" >&2; exit 1;;
esac
GIAB_VER_TAG="HG001_v4.2.1"

if [ "$DRY_RUN" != 1 ]; then
  for v in FQ_R1 FQ_R2 FASTA TARGET WORKDIR; do [ -n "${!v:-}" ] || { echo "config is missing $v" >&2; exit 1; }; done
fi
CALIB="$WORKDIR/plan/calibration.csv"
run "mkdir -p '$WORKDIR/bam' '$WORKDIR/runs/$RID' '$WORKDIR/happy' '$WORKDIR/plan'"
echo "=== $RID  (platform=$PLATFORM ref=$REF target=${TARGETX}X seed=$SEED hw=$HW align_only=$ALIGN_ONLY) ==="

# ---- 1) sampling fraction (in on-target-X space) ----
if [ "$TARGETX" = full ]; then
  FRAC=1
else
  XFULL=""
  if have "$CALIB"; then XFULL=$(awk -F, -v s="$SAMPLE" -v r="$REF" '$1==s&&$2==r{print $3}' "$CALIB" | head -1); fi
  if [ "$DRY_RUN" = 1 ]; then XFULL="${XFULL:-1000}"; fi   # dry-run uses a placeholder X_full to demonstrate the fraction
  [ -n "$XFULL" ] || { echo "no calibration: $CALIB has no x_full for $SAMPLE/$REF. Run Phase-0 first: scripts/calibrate.sh" >&2; exit 1; }
  FRAC=$(awk -v x="$TARGETX" -v xf="$XFULL" 'BEGIN{ if(xf+0<=0){print 1}else{f=x/xf; if(f>=1){print 1}else{printf "%.6f", f}} }')
fi

# ---- 2) align + mark duplicates (the BAM cache key excludes the caller, so gatk4 and dv share
#          one alignment) ----
BAMKEY="$SAMPLE.$REF.$ALIGNER.${DEPTHTOK}.s$SEED"
BAMDIR="$WORKDIR/bam/$BAMKEY"; DEDUP="$BAMDIR/dedup.bam"; COVPFX="$BAMDIR/cov"
run "mkdir -p '$BAMDIR'"
# 2a) subsample the fastq (full depth uses the originals; idempotent)
#   The sub-fastq cache key deliberately EXCLUDES the aligner, so the CPU (bwa) and GPU (fq2bam)
#   branches at the same depth x seed share one copy. This saves a full pass over the fastq and disk
#   space, and guarantees the CPU/GPU concordance comparison sees byte-identical input reads.
#   pigz is multi-threaded; seqtk -s<seed> is a deterministic single-pass Bernoulli sample; the
#   tmp+mv atomic write stops the CPU and GPU branches racing to create the same file.
if [ "$FRAC" = 1 ]; then
  SUB_R1=$FQ_R1; SUB_R2=$FQ_R2
else
  SUBDIR="$WORKDIR/bam/_sub/$SAMPLE.$REF.${DEPTHTOK}.s$SEED"
  SUB_R1="$SUBDIR/sub_R1.fq.gz"; SUB_R2="$SUBDIR/sub_R2.fq.gz"
  run "mkdir -p '$SUBDIR'"
  if have "$SUB_R1" && have "$SUB_R2"; then echo "  subsample already exists (shared), skipping (frac=$FRAC)"; else
    run "$DECOMP '$FQ_R1' | seqtk sample -s$SEED - $FRAC | $COMPRESS > '$SUB_R1.tmp.$$' && mv -f '$SUB_R1.tmp.$$' '$SUB_R1'"
    run "$DECOMP '$FQ_R2' | seqtk sample -s$SEED - $FRAC | $COMPRESS > '$SUB_R2.tmp.$$' && mv -f '$SUB_R2.tmp.$$' '$SUB_R2'"
  fi
fi
# 2b) alignment (idempotent: skipped when the BAM exists, so several callers reuse it)
if have "$DEDUP"; then
  echo "  alignment already exists, reusing $DEDUP"
else
  RG="@RG\\tID:$BAMKEY\\tSM:$SAMPLE\\tPL:$PL\\tPU:$SAMPLE\\tLB:$SAMPLE"   # PU is mandatory for pbrun fq2bam (bwa does not require it)
  case "$ALIGNER" in
    bwa-mem)
      run "bwa mem -t $THREADS -R '$RG' '$FASTA' '$SUB_R1' '$SUB_R2' | samtools sort -@ $THREADS -o '$BAMDIR/sorted.bam' -"
      run "gatk --java-options '$GATK_JAVA_OPTS' MarkDuplicates -I '$BAMDIR/sorted.bam' -O '$DEDUP' -M '$BAMDIR/dup_metrics.txt'"  # the GATK image bundles Picard 2.27.5, so Picard need not be installed separately
      run "samtools index '$DEDUP'"       # produces .bam.bai, accepted by both mosdepth and gatk (avoids Picard CREATE_INDEX .bai naming)
      run "rm -f '$BAMDIR/sorted.bam'"
      ;;
    pbrun_fq2bam)
      run "pbrun fq2bam --ref '$FASTA' --in-fq '$SUB_R1' '$SUB_R2' '$RG' --out-bam '$DEDUP'"
      run "samtools index '$DEDUP'"
      ;;
    *) echo "unknown aligner: $ALIGNER" >&2; exit 1;;
  esac
fi

# ---- 3) mosdepth coverage (on-target; idempotent) -> one run_coverage row ----
if ! have "$COVPFX.mosdepth.summary.txt"; then
  run "mosdepth -t $THREADS -n --by '$TARGET' --thresholds 1,10,20,30,50 '$COVPFX' '$DEDUP'"
fi
MEANX="NA"; MEDX="NA"; USEDGB="NA"; PCT10="NA"; PCT20="NA"; PCT30="NA"; PCT50="NA"
if have "$COVPFX.mosdepth.summary.txt"; then
  MEANX=$(awk '$1=="total_region"{print $4}' "$COVPFX.mosdepth.summary.txt")
  USEDGB=$(awk '$1=="total"{printf "%.2f", $3/1e9}' "$COVPFX.mosdepth.summary.txt")
  read -r PCT10 PCT20 PCT30 PCT50 < <(zcat "$COVPFX.thresholds.bed.gz" | awk '
    NR>1{len=$3-$2; tot+=len; c10+=$6; c20+=$7; c30+=$8; c50+=$9}
    END{ if(tot>0) printf "%.4f %.4f %.4f %.4f\n", c10/tot, c20/tot, c30/tot, c50/tot; else print "NA NA NA NA" }') || true
fi

# ---- 3b) full + bwa-mem: upsert the calibration table (x_full = measured on-target depth) ----
if [ "$TARGETX" = full ] && [ "$ALIGNER" = bwa-mem ] && [ "$DRY_RUN" != 1 ]; then
  [ -f "$CALIB" ] || echo "sample,reference,x_full,gb_full,run_date" > "$CALIB"
  if ! awk -F, -v s="$SAMPLE" -v r="$REF" 'NR>1&&$1==s&&$2==r{f=1} END{exit !f}' "$CALIB"; then
    echo "$SAMPLE,$REF,$MEANX,$USEDGB,$(date +%F)" >> "$CALIB"
    echo "  calibration written to $CALIB : $SAMPLE/$REF x_full=$MEANX gb_full=$USEDGB"
  fi
fi

append() { [ "$DRY_RUN" = 1 ] && { echo "  [row→$(basename "$1")] $3"; return; }; [ -f "$1" ] || echo "$2" > "$1"; echo "$3" >> "$1"; }
append "$WORKDIR/plan/run_coverage.csv" \
  "run_id,sample,platform,reference,target_depth_x,seed,mean_target_depth_x,median_target_depth_x,used_gb,pct_target_10x,pct_target_20x,pct_target_30x,pct_target_50x,coverage_tool,coverage_tool_version,run_date" \
  "$RID,$SAMPLE,$PLATFORM,$REF,$TARGETX,$SEED,$MEANX,$MEDX,$USEDGB,$PCT10,$PCT20,$PCT30,$PCT50,mosdepth,${MOSDEPTH_VERSION:-NA},$(date +%F)"

if [ "$ALIGN_ONLY" = 1 ]; then echo "=== align-only done: $RID (bam + coverage ready) ==="; exit 0; fi

# ---- 3c) BQSR (GATK family only; DeepVariant does not use it, see EXPERIMENT-DESIGN section 3)
#          -> recal.bam; DeepVariant keeps using dedup.bam ----
CALLBAM="$DEDUP"
case "$CALLER" in gatk4_hc|gatk3|pbrun_haplotypecaller) DO_BQSR=1;; *) DO_BQSR=0;; esac
if [ "$DO_BQSR" = 1 ]; then
  case "$REF" in
    grch38) KDBSNP=$KNOWN_DBSNP_GRCH38; KMILLS=$KNOWN_MILLS_GRCH38; KINDEL=$KNOWN_INDELS_GRCH38 ;;
    *)      KDBSNP=$KNOWN_DBSNP_B37;    KMILLS=$KNOWN_MILLS_B37;    KINDEL=$KNOWN_INDELS_B37 ;;
  esac
  [ "$DRY_RUN" = 1 ] || for k in "$KDBSNP" "$KMILLS" "$KINDEL"; do
    [ -n "$k" ] && [ -e "$k" ] || { echo "BQSR is missing a known-sites file: '$k' (set KNOWN_*_${REF^^} in the config)" >&2; exit 1; }; done
  RTABLE="$BAMDIR/recal.$CALLER.table"; RECAL="$BAMDIR/recal.$CALLER.bam"
  if have "$RECAL"; then echo "  BQSR recal.bam already exists, reusing $RECAL"; else
    if [ "$HW" = gpu ]; then
      run "pbrun bqsr --ref '$FASTA' --in-bam '$DEDUP' --knownSites '$KDBSNP' --knownSites '$KMILLS' --knownSites '$KINDEL' --out-recal-file '$RTABLE'"
      run "pbrun applybqsr --ref '$FASTA' --in-bam '$DEDUP' --in-recal-file '$RTABLE' --out-bam '$RECAL'"
    else
      run "gatk --java-options '$GATK_JAVA_OPTS' BaseRecalibrator -R '$FASTA' -I '$DEDUP' -L '$TARGET' --known-sites '$KDBSNP' --known-sites '$KMILLS' --known-sites '$KINDEL' -O '$RTABLE'"
      run "gatk --java-options '$GATK_JAVA_OPTS' ApplyBQSR -R '$FASTA' -I '$DEDUP' -L '$TARGET' --bqsr-recal-file '$RTABLE' -O '$RECAL'"
    fi
  fi
  CALLBAM="$RECAL"
fi

# ---- 4) calling (restricted to the V6 target; GATK family uses recal.bam, DeepVariant dedup.bam)
#          -> one directory per run ----
OUT="$WORKDIR/runs/$RID"; CALLS="$OUT/calls.vcf.gz"; NORM="$OUT/calls.norm.vcf.gz"; _TCALL=$SECONDS
case "$CALLER" in
  gatk4_hc)               run "gatk HaplotypeCaller -R '$FASTA' -I '$CALLBAM' -L '$TARGET' -O '$CALLS'" ;;
  gatk3)                  run "java -jar \$GATK3_JAR -T HaplotypeCaller -R '$FASTA' -I '$CALLBAM' -L '$TARGET' -o '$CALLS'" ;;
  deepvariant)            run "run_deepvariant --model_type=WES --ref='$FASTA' --reads='$CALLBAM' --regions='$TARGET' --output_vcf='$CALLS' --num_shards=$THREADS" ;;
  pbrun_haplotypecaller)  run "pbrun haplotypecaller --ref '$FASTA' --in-bam '$CALLBAM' --interval-file '$TARGET' --out-variants '$CALLS'" ;;
  pbrun_deepvariant)      run "pbrun deepvariant --mode shortread --use-wes-model --ref '$FASTA' --in-bam '$CALLBAM' --interval-file '$TARGET' --out-variants '$CALLS'" ;;   # --use-wes-model matches the CPU side's --model_type=WES; without it the WGS model is used, which
# misses heterozygous calls and breaks the concordance comparison
  *) echo "unknown caller: $CALLER" >&2; exit 1;;
esac

CALL_SEC=$(( SECONDS - _TCALL ))   # calling wall-clock seconds (basis of the CPU/GPU speed-up)

# ---- 4b) hard-filtering (GATK family only, and INDELs only; SNVs are left unfiltered because SNV
#          hard-filtering was measured to be net-negative. DeepVariant is skipped.) ----
#   VariantFiltration flags only INDELs, so hap.py filter=PASS is "raw SNVs + hard-filtered INDELs"
#   (the recommended set) while ALL remains the complete callset.
#   The logic lives in hardfilter_vcf.sh (shared with apply_hardfilter.sh, a single source of truth).
#   Idempotent: skipped when the filtered VCF already exists.
NORMIN="$CALLS"
case "$CALLER" in
  gatk4_hc|gatk3|pbrun_haplotypecaller)
    FILT="$OUT/calls.filt.vcf.gz"
    if have "$FILT"; then echo "  hard-filtered VCF already exists, reusing $FILT"; else
      run "DRY_RUN=$DRY_RUN THREADS=$THREADS bash '$HERE/hardfilter_vcf.sh' '$CALLS' '$FASTA' '$FILT'"
    fi
    NORMIN="$FILT"
    ;;
esac

# ---- 5) normalisation ----
run "bcftools norm -m-any -f '$FASTA' '$NORMIN' -Oz -o '$NORM'"
run "bcftools index -ft '$NORM'"   # produces .tbi (equivalent to tabix -p vcf; bcftools bundles htslib, so tabix is not needed)

# ---- 6) hap.py benchmark ----
HAPPY_PFX="$WORKDIR/happy/$RID"
[ "$RUN_HAPPY" = 1 ] && run "HAPPY_CONFIG='$CONFIG' bash '$HERE/run_happy.sh' --query '$NORM' --ref '$REF' --out '$HAPPY_PFX' --threads $THREADS"

# ---- 7) append the manifest row (in schema column order) ----
if [ "$HW" = gpu ]; then PBV=$PARABRICKS_VERSION; ALV=$PARABRICKS_VERSION; else PBV=NA; ALV=$BWA_VERSION; fi
case "$CALLER" in
  gatk4_hc|pbrun_haplotypecaller) TOOLV=$GATK4_VERSION;;
  deepvariant|pbrun_deepvariant)  TOOLV=$DEEPVARIANT_VERSION;;
  gatk3)                          TOOLV=$GATK3_VERSION;;
  *)                              TOOLV=NA;;
esac
append "$WORKDIR/plan/run_manifest.csv" \
  "run_id,sample,platform,reference,aligner,aligner_version,caller,tool_version,hardware,parabricks_version,target_depth_x,mean_target_depth_x,seed,benchmark_tool,benchmark_tool_version,giab_version,target_bed,happy_prefix,run_date" \
  "$RID,$SAMPLE,$PLATFORM,$REF,$ALIGNER,$ALV,$CALLER,$TOOLV,$HW,$PBV,$TARGETX,$MEANX,$SEED,happy,${HAPPY_VERSION:-NA},$GIAB_VER_TAG,$(basename "$TARGET"),$HAPPY_PFX,$(date +%F)"

# ---- 8) resource usage and timing (feeds the CPU/GPU speed-up analysis; peak host and GPU memory
#          are sampled separately by scripts/monitor.sh) ----
append "$WORKDIR/plan/resource_usage.csv" \
  "run_id,sample,platform,reference,caller,hardware,target_depth_x,seed,call_sec,total_sec,run_date" \
  "$RID,$SAMPLE,$PLATFORM,$REF,$CALLER,$HW,$TARGETX,$SEED,${CALL_SEC:-NA},$(( SECONDS - _T0 )),$(date +%F)"

echo "=== done: $RID -> manifest + coverage + resource rows appended; hap.py prefix $HAPPY_PFX ==="
echo "    aggregate with: python3 scripts/aggregate_happy.py --manifest '$WORKDIR/plan/run_manifest.csv'"
