#!/usr/bin/env bash
# hardfilter_vcf.sh -- GATK best-practice hard-filtering (germline small variants, single-sample WES).
#
# INDELS ONLY -- SNVs are deliberately left unfiltered. This was calibrated empirically at two depths
# (30x and full):
#   INDEL filtering raises F1 by +0.02 to +0.08 (precision climbs sharply, recall is untouched), so it
#   is clearly worth applying.
#   SNV filtering lowers F1 by -0.005 to -0.013 (the recall lost to true variants exceeds the gain from
#   removing false positives; MQ<40 and MQRankSum are the main culprits and are depth-independent).
#   Raw SNV precision is already >0.98, so SNVs are left alone and all remain PASS.
#
# A single VariantFiltration pass with the JEXL guard !vc.isSNP(): only INDEL conditions exist, and
# SNV records short-circuit out, so they keep FILTER=PASS.
#   As a result hap.py's filter=PASS is the recommended set "raw SNVs + hard-filtered INDELs" (every
#   caller therefore uses PASS uniformly, and the aggregation and plotting code needs no changes),
#   while filter=ALL remains the complete unfiltered callset. No record is dropped and the output
#   stays coordinate-sorted, so no re-sort is needed.
#
# INDEL thresholds (the GATK germline recommendations):
#   QD<2.0 | QUAL<30 | FS>200.0 | ReadPosRankSum<-20.0 | SOR>10.0
#   Each condition gets its own filter name so the reason is traceable. JEXL '&&' short-circuits, so
#   SNVs and records missing an annotation are never mis-flagged.
#
# Usage: [DRY_RUN=1] scripts/hardfilter_vcf.sh <raw.vcf.gz> <ref.fasta> <out.filt.vcf.gz>
#   DeepVariant applies its own internal filtering and does not go through this script; it is used
#   only for gatk4_hc / gatk3 / pbrun_haplotypecaller.
set -euo pipefail
RAW="${1:?usage: hardfilter_vcf.sh <raw.vcf.gz> <ref.fasta> <out.vcf.gz>}"
FASTA="${2:?missing ref.fasta}"
OUT="${3:?missing out.vcf.gz}"
DRY_RUN="${DRY_RUN:-0}"
export PATH="${WRAPPER_BIN:-$HOME/bin}:$PATH"   # the gatk wrapper lives in ~/bin; extend PATH explicitly for non-interactive subshells
run(){ echo "+ $*"; [ "$DRY_RUN" = 1 ] || eval "$*"; }

# A single VariantFiltration pass; INDEL conditions carry the !vc.isSNP() guard.
# (For mixed sites vc.isSNP() is false, so they fall under the INDEL thresholds -- the conservative
# choice.)
run "gatk VariantFiltration -R '$FASTA' -V '$RAW' \
  --filter-expression '!vc.isSNP() && QD < 2.0'               --filter-name indel_QD2 \
  --filter-expression '!vc.isSNP() && QUAL < 30.0'            --filter-name indel_QUAL30 \
  --filter-expression '!vc.isSNP() && FS > 200.0'             --filter-name indel_FS200 \
  --filter-expression '!vc.isSNP() && ReadPosRankSum < -20.0' --filter-name indel_ReadPosRankSum \
  --filter-expression '!vc.isSNP() && SOR > 10.0'             --filter-name indel_SOR10 \
  --verbosity ERROR \
  -O '$OUT'"
# Only INDEL conditions exist (all guarded by !vc.isSNP()), so SNV records short-circuit out, are
# never flagged by any filter, and keep FILTER=PASS.
echo "  hard-filtering done -> $OUT"
