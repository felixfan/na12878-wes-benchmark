#!/usr/bin/env bash
# snv_hardfilter_test.sh -- evidence that SNV hard-filtering is net-negative (which is why the
#   production pipeline is INDEL-only).
#   In benchmark_long_all.csv the PASS and ALL SNV rows are identical, because SNVs are never
#   filtered -- so the released data contain no artefact backing the claim that SNV hard-filtering
#   lowers F1. This script takes the existing gatk4_hc raw VCFs, applies *SNV* hard-filtering,
#   re-runs hap.py, and reports SNV F1 raw (ALL) versus filtered (PASS) to demonstrate the loss.
#   It does not change the production pipeline (hardfilter_vcf.sh stays INDEL-only); this is a
#   one-off diagnostic.
#
# The GATK germline SNV thresholds:
#   QD<2.0 | FS>60.0 | MQ<40.0 | MQRankSum<-12.5 | ReadPosRankSum<-8.0 | SOR>3.0
#   The JEXL vc.isSNP() guard means only SNVs are flagged; INDELs are untouched and stay PASS. In the
#   hap.py SNP rows, PASS therefore means SNV-filtered and ALL means raw.
#
# Usage:
#   HAPPY_CONFIG=configs/na12878.paths.sh scripts/snv_hardfilter_test.sh [run_id ...]
#   The default run_id set is gatk4_hc.cpu.grch38 across the three platforms at {30x.s42, full.s0} (6 runs).
#   DRY_RUN=1 scripts/snv_hardfilter_test.sh
# Output: $WORKDIR/qc/snv_hardfilter_test.csv
#   (run_id,sample,platform,target_depth_x,snp_f1_raw,snp_f1_snvhf,delta_f1,snp_prec_raw,snp_prec_snvhf,snp_recall_raw,snp_recall_snvhf)
set -euo pipefail

CONFIG="${HAPPY_CONFIG:-$(dirname "$0")/../configs/na12878.paths.sh}"
[ -f "$CONFIG" ] || { echo "missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
HERE="$(cd "$(dirname "$0")" && pwd)"
export PATH="${WRAPPER_BIN:-$HOME/bin}:$PATH"
DRY_RUN="${DRY_RUN:-0}"; THREADS="${THREADS:-8}"
run(){ echo "+ $*"; [ "$DRY_RUN" = 1 ] || eval "$*"; }
have(){ [ "$DRY_RUN" != 1 ] && [ -s "$1" ]; }

# default run_id list
if [ $# -gt 0 ]; then RIDS=("$@"); else
  RIDS=(); for S in ZY_Illumina_NA12878 ZY_BGI_NA12878 ZY_GeneMind_NA12878; do
    RIDS+=("$S.grch38.bwa-mem.gatk4_hc.cpu.30x.s42" "$S.grch38.bwa-mem.gatk4_hc.cpu.full.s0"); done
fi

QCDIR="$WORKDIR/qc"; OUTCSV="$QCDIR/snv_hardfilter_test.csv"
run "mkdir -p '$QCDIR'"
[ "$DRY_RUN" = 1 ] || echo "run_id,sample,platform,target_depth_x,snp_f1_raw,snp_f1_snvhf,delta_f1,snp_prec_raw,snp_prec_snvhf,snp_recall_raw,snp_recall_snvhf,run_date" > "$OUTCSV"

# read one metric from the SNP row of a given Filter in hap.py summary.csv (located by header name)
metric(){ # $1=summary.csv $2=ALL|PASS $3=METRIC.F1_Score|METRIC.Precision|METRIC.Recall
  python3 - "$1" "$2" "$3" <<'PY'
import csv,sys
f,filt,col=sys.argv[1],sys.argv[2],sys.argv[3]
for r in csv.DictReader(open(f)):
    if r.get("Type")=="SNP" and r.get("Filter")==filt:
        print(r.get(col,"NA")); break
else: print("NA")
PY
}

for RID in "${RIDS[@]}"; do
  echo "=== SNV-hardfilter test: $RID ==="
  IFS='.' read -r SAMPLE REF ALIGNER CALLER HW DEPTHTOK SEEDS <<<"$RID"
  case "$SAMPLE" in ZY_Illumina_*) PLATFORM=illumina_novaseq6000;; ZY_BGI_*) PLATFORM=bgi_dnbseq_t7;; ZY_GeneMind_*) PLATFORM=genemind_surfseq5000;; *) echo "skipping unknown sample: $RID"; continue;; esac
  [ "$DEPTHTOK" = full ] && TARGETX=full || TARGETX="${DEPTHTOK%x}"
  RAW="$WORKDIR/runs/$RID/calls.vcf.gz"
  if [ "$DRY_RUN" != 1 ] && [ ! -s "$RAW" ]; then echo "  !! raw VCF missing: $RAW (run run_one.sh for this run_id first)" >&2; continue; fi

  SNVHF="$WORKDIR/runs/$RID/calls.snvhf.vcf.gz"; NORM="$WORKDIR/runs/$RID/calls.snvhf.norm.vcf.gz"
  # 1) hard-filter SNVs only (vc.isSNP() guard; INDELs untouched)
  if have "$SNVHF"; then echo "  SNV-filtered VCF already exists, reusing"; else
    run "gatk VariantFiltration -R '$REF_GRCH38' -V '$RAW' \
      --filter-expression 'vc.isSNP() && QD < 2.0'                --filter-name snp_QD2 \
      --filter-expression 'vc.isSNP() && FS > 60.0'               --filter-name snp_FS60 \
      --filter-expression 'vc.isSNP() && MQ < 40.0'               --filter-name snp_MQ40 \
      --filter-expression 'vc.isSNP() && MQRankSum < -12.5'       --filter-name snp_MQRankSum \
      --filter-expression 'vc.isSNP() && ReadPosRankSum < -8.0'   --filter-name snp_ReadPosRankSum \
      --filter-expression 'vc.isSNP() && SOR > 3.0'               --filter-name snp_SOR3 \
      --verbosity ERROR -O '$SNVHF'"
  fi
  # 2) normalise (same as the production pipeline)
  run "bcftools norm -m-any -f '$REF_GRCH38' '$SNVHF' -Oz -o '$NORM'"
  run "bcftools index -ft '$NORM'"
  # 3) hap.py via run_happy.sh, with exactly the same parameters as the main benchmark
  HPFX="$WORKDIR/happy/${RID}.snvhf"
  run "HAPPY_CONFIG='$CONFIG' bash '$HERE/run_happy.sh' --query '$NORM' --ref '$REF' --out '$HPFX' --threads $THREADS"

  if [ "$DRY_RUN" != 1 ]; then
    SUM="$HPFX.summary.csv"
    F1R=$(metric "$SUM" ALL  METRIC.F1_Score);  F1F=$(metric "$SUM" PASS METRIC.F1_Score)
    PR=$(metric "$SUM" ALL  METRIC.Precision); PF=$(metric "$SUM" PASS METRIC.Precision)
    RR=$(metric "$SUM" ALL  METRIC.Recall);    RF=$(metric "$SUM" PASS METRIC.Recall)
    DELTA=$(awk -v a="$F1F" -v b="$F1R" 'BEGIN{if(a=="NA"||b=="NA")print "NA"; else printf "%.4f", a-b}')
    echo "$RID,$SAMPLE,$PLATFORM,$TARGETX,$F1R,$F1F,$DELTA,$PR,$PF,$RR,$RF,$(date +%F)" >> "$OUTCSV"
    echo "  SNV F1 raw=$F1R -> SNV-hardfiltered=$F1F (delta=$DELTA; negative means a net loss, supporting INDEL-only)"
  fi
done
echo "=== done -> $OUTCSV  (a negative delta demonstrates that SNV hard-filtering is net-negative) ==="
