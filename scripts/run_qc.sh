#!/usr/bin/env bash
# run_qc.sh <sample> <reference> -- sequencing quality QC (Table 1 of the paper).
#   Raw fastq via fastp, plus the full-depth BAM via samtools flagstat, Picard MarkDuplicates
#   (duplicate rate) and GATK CollectHsMetrics.
#   Reports are written to $WORKDIR/qc/<sample>.<ref>/ and collected by scripts/aggregate_qc.py.
#
# Prerequisite: the full-depth BAM already exists
# ($WORKDIR/bam/<sample>.<ref>.bwa-mem.full.s0/dedup.bam, built by calibrate.sh or a full-depth run).
# Usage:
#   HAPPY_CONFIG=configs/na12878.paths.sh scripts/run_qc.sh ZY_Illumina_NA12878 grch38
#   DRY_RUN=1 scripts/run_qc.sh ZY_GeneMind_NA12878 grch38
set -euo pipefail
CONFIG="${HAPPY_CONFIG:-$(dirname "$0")/../configs/na12878.paths.sh}"
[ -f "$CONFIG" ] || { echo "missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
DRY_RUN="${DRY_RUN:-0}"; THREADS="${THREADS:-8}"
run() { echo "+ $*"; [ "$DRY_RUN" = 1 ] || eval "$*"; }

SAMPLE="${1:?usage: run_qc.sh <sample> <reference>}"; REF="${2:?usage: run_qc.sh <sample> <reference>}"
case "$SAMPLE" in
  ZY_Illumina_*) PLATFORM=illumina_novaseq6000; FQ_R1=$FQ_ILLUMINA_R1; FQ_R2=$FQ_ILLUMINA_R2 ;;
  ZY_BGI_*)   PLATFORM=bgi_dnbseq_t7;        FQ_R1=$FQ_BGI_R1;      FQ_R2=$FQ_BGI_R2 ;;
  ZY_GeneMind_*)       PLATFORM=genemind_surfseq5000; FQ_R1=$FQ_GM_R1;       FQ_R2=$FQ_GM_R2 ;;
  *) echo "unknown sample: $SAMPLE" >&2; exit 1;;
esac
case "$REF" in
  grch38)          FASTA=$REF_GRCH38; TARGET=$V6_BED_GRCH38 ;;
  b37|grch37|hg19) FASTA=$REF_B37;    TARGET=$V6_BED_B37 ;;
  *) echo "unknown reference: $REF" >&2; exit 1;;
esac
DICT="${FASTA%.fasta}.dict"                       # human_g1k_v37.dict / Homo_sapiens_assembly38.dict
QCDIR="$WORKDIR/qc/$SAMPLE.$REF"
BAMDIR="$WORKDIR/bam/$SAMPLE.$REF.bwa-mem.full.s0"; BAM="$BAMDIR/dedup.bam"
run "mkdir -p '$QCDIR'"
echo "=== QC $SAMPLE ($PLATFORM) / $REF ==="

# 1) raw fastq QC with fastp ----
# Only the report is wanted, not trimmed reads, so -o/-O are omitted (fastp can emit a report alone;
# passing identical -o/-O paths would be an error). The statistics used come from before_filtering.
run "fastp -i '$FQ_R1' -I '$FQ_R2' -j '$QCDIR/fastp.json' -h '$QCDIR/fastp.html' -w $(( THREADS>16 ? 16 : THREADS ))"

# 2) full-depth BAM: mapping rate and duplicate rate ----
if [ "$DRY_RUN" != 1 ] && [ ! -s "$BAM" ]; then echo "full-depth BAM missing: $BAM (run calibrate.sh or the full-depth run first)" >&2; exit 1; fi
run "samtools flagstat '$BAM' > '$QCDIR/flagstat.txt'"
if [ "$DRY_RUN" = 1 ] || [ -f "$BAMDIR/dup_metrics.txt" ]; then run "cp '$BAMDIR/dup_metrics.txt' '$QCDIR/dup_metrics.txt'"; fi

# 3) V6 interval_list (CollectHsMetrics needs an interval_list, not a BED; built once per build and cached) ----
IL="$WORKDIR/qc/V6.$REF.interval_list"
if [ "$DRY_RUN" = 1 ] || [ ! -f "$IL" ]; then run "gatk BedToIntervalList -I '$TARGET' -O '$IL' -SD '$DICT'"; else echo "  interval_list already exists: $IL"; fi

# 4) GATK CollectHsMetrics (on-target fraction, enrichment and depth; bait = target = V6 Covered) ----
run "gatk CollectHsMetrics -I '$BAM' -O '$QCDIR/hs_metrics.txt' -R '$FASTA' -BI '$IL' -TI '$IL'"

echo "=== QC done for $SAMPLE/$REF -> $QCDIR (fastp.json / flagstat.txt / dup_metrics.txt / hs_metrics.txt) ==="
echo "    aggregate with: python3 scripts/aggregate_qc.py --qcdir '$WORKDIR/qc'"
