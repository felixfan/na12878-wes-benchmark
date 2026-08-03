#!/usr/bin/env bash
# qc_identity.sh -- quantify that the three platform libraries come from the same DNA (somalier relate).
#   Replaces a weak bound such as "r > 0.9" with exact genotype concordance and relatedness:
#   replicate libraries of one individual should give relatedness ~1 and homozygous concordance ~1.
#   Extracts directly from the three full-depth GRCh38 dedup.bam files.
#
# Prerequisites (not covered by the paths config):
#   1) A somalier binary or container: https://github.com/brentp/somalier/releases
#      docker: brentp/somalier:latest
#   2) A sites file matching the reference build. For GRCh38 use sites.hg38.vcf.gz (chr-prefixed, so
#      it matches assembly38). It is published as an issue attachment rather than a release asset:
#      https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
#
# Usage:
#   HAPPY_CONFIG=configs/na12878.paths.sh SOMALIER_SITES=/path/sites.hg38.vcf.gz \
#     scripts/qc_identity.sh
#   DRY_RUN=1 scripts/qc_identity.sh
#   SOM_BIN=somalier                                 # use a local binary instead of the container
# Output: $WORKDIR/qc/identity.pairs.tsv + identity.samples.tsv (somalier's native format).
#   The informative columns of pairs.tsv are relatedness, hom_concordance, shared_hets, ibs0, ibs2
#   and n, giving the exact agreement for each of the three platform pairs.
set -euo pipefail

CONFIG="${HAPPY_CONFIG:-$(dirname "$0")/../configs/na12878.paths.sh}"
[ -f "$CONFIG" ] || { echo "missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
export PATH="${WRAPPER_BIN:-$HOME/bin}:$PATH"
DRY_RUN="${DRY_RUN:-0}"
run(){ echo "+ $*"; [ "$DRY_RUN" = 1 ] || eval "$*"; }

: "${SOMALIER_SITES:?set SOMALIER_SITES=<sites.hg38.vcf.gz> (contig naming must match REF_GRCH38)}"
SOM_IMAGE="${SOM_IMAGE:-brentp/somalier:latest}"
SOM_BIN="${SOM_BIN:-}"   # non-empty: local binary; empty: docker run

QCDIR="$WORKDIR/qc"; EXTRACT="$QCDIR/somalier_extract"
run "mkdir -p '$EXTRACT'"
REFDIR="$(dirname "$REF_GRCH38")"; SITESDIR="$(dirname "$SOMALIER_SITES")"
som() {  # pass the somalier subcommand and arguments through; host paths resolve inside the container
  if [ -n "$SOM_BIN" ]; then echo "$SOM_BIN $*"
  else echo "docker run --rm -v '$WORKDIR:$WORKDIR' -v '$REFDIR:$REFDIR' -v '$SITESDIR:$SITESDIR' -w '$WORKDIR' '$SOM_IMAGE' somalier $*"; fi
}

# 1) extract per sample (the SM tag comes from the @RG line written by run_one.sh)
#    -> $EXTRACT/<SAMPLE>.somalier
for SAMPLE in ZY_Illumina_NA12878 ZY_BGI_NA12878 ZY_GeneMind_NA12878; do
  BAM="$WORKDIR/bam/$SAMPLE.grch38.bwa-mem.full.s0/dedup.bam"
  if [ "$DRY_RUN" != 1 ] && [ ! -s "$BAM" ]; then
    echo "  !! full-depth BAM missing: $BAM (run ALIGN_ONLY=1 run_one.sh $SAMPLE.grch38.bwa-mem.gatk4_hc.cpu.full.s0 first)" >&2; continue
  fi
  run "$(som extract -d "'$EXTRACT'" --sites "'$SOMALIER_SITES'" -f "'$REF_GRCH38'" "'$BAM'")"
done

# 2) relate -> pairs/samples tsv (all pairwise combinations of the three samples)
run "$(som relate -o "'$QCDIR/identity'" "'$EXTRACT'"/*.somalier)"
echo "=== done -> $QCDIR/identity.pairs.tsv (relatedness/hom_concordance/shared_hets) + identity.samples.tsv ==="
echo "    Same-DNA libraries are expected to give relatedness ~1 and hom_concordance ~1."
