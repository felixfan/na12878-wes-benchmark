#!/usr/bin/env bash
# qc_contamination.sh -- quantify cross-sample contamination per sample (VerifyBamID2 FREEMIX).
#   A high mapping rate only rules out non-human sequence; it says nothing about human-human
#   cross-contamination, which is what FREEMIX measures.
#   Runs against the cached full-depth GRCh38 dedup.bam of each sample.
#
# Prerequisites (not covered by the paths config):
#   1) A VerifyBamID2 binary. Using a local binary via VBID2_BIN is the least troublesome route:
#      - from source: cd VerifyBamID && mkdir -p build && cd build && cmake .. && make
#        -> the binary lands in <repo>/bin/VerifyBamID; set VBID2_BIN to it.
#      - or via conda: mamba install -c bioconda -c conda-forge verifybamid2, then VBID2_BIN=verifybamid2.
#      - a container also works, but check which biocontainer build tags currently exist at
#        https://quay.io/repository/biocontainers/verifybamid2?tab=tags -- do not guess a tag.
#   2) An SVD marker panel matching the reference build. For GRCh38 use the prefix
#      1000g.phase3.100k.b38.vcf.gz.dat, found under resource/ in the VerifyBamID repository.
#      NOTE on contig naming: assembly38 uses chr-prefixed names. A panel without the prefix produces
#      "no overlapping markers"; pick a panel matching the reference, and always check that #SNPS in
#      the .selfSM output is non-zero.
#
# Usage:
#   HAPPY_CONFIG=configs/na12878.paths.sh VBID2_SVD=/path/1000g.phase3.100k.b38.vcf.gz.dat \
#     scripts/qc_contamination.sh
#   DRY_RUN=1 scripts/qc_contamination.sh          # print commands only
#   VBID2_BIN=verifybamid2                          # use a local binary instead of the container
# Output: $WORKDIR/qc/contamination_freemix.csv  (sample,platform,freemix,n_snp,tool,run_date)
set -euo pipefail

CONFIG="${HAPPY_CONFIG:-$(dirname "$0")/../configs/na12878.paths.sh}"
[ -f "$CONFIG" ] || { echo "missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
export PATH="${WRAPPER_BIN:-$HOME/bin}:$PATH"
DRY_RUN="${DRY_RUN:-0}"; THREADS="${THREADS:-8}"
run(){ echo "+ $*"; [ "$DRY_RUN" = 1 ] || eval "$*"; }

: "${VBID2_SVD:?set VBID2_SVD=<SVD prefix> (GRCh38: .../1000g.phase3.100k.b38.vcf.gz.dat)}"
VBID2_IMAGE="${VBID2_IMAGE:-quay.io/biocontainers/verifybamid2:2.0.1--h468198e_6}"
VBID2_BIN="${VBID2_BIN:-}"   # non-empty: call a local binary; empty: fall back to docker run

QCDIR="$WORKDIR/qc"; OUTCSV="$QCDIR/contamination_freemix.csv"
run "mkdir -p '$QCDIR'"
[ "$DRY_RUN" = 1 ] || echo "sample,platform,freemix,n_snp,tool,tool_source,run_date" > "$OUTCSV"

# Invoke VerifyBamID2 either as a local binary or through docker (mounting WORKDIR, the reference
# directory and the SVD directory so that host paths resolve unchanged inside the container)
REFDIR="$(dirname "$REF_GRCH38")"; SVDDIR="$(dirname "$VBID2_SVD")"
vbid2() {  # $1=bam  $2=out_prefix
  if [ -n "$VBID2_BIN" ]; then
    echo "$VBID2_BIN --SVDPrefix '$VBID2_SVD' --Reference '$REF_GRCH38' --BamFile '$1' --Output '$2' --NumThread $THREADS"
  else
    echo "docker run --rm -v '$WORKDIR:$WORKDIR' -v '$REFDIR:$REFDIR' -v '$SVDDIR:$SVDDIR' -w '$WORKDIR' '$VBID2_IMAGE' \
verifybamid2 --SVDPrefix '$VBID2_SVD' --Reference '$REF_GRCH38' --BamFile '$1' --Output '$2' --NumThread $THREADS"
  fi
}

for SAMPLE in ZY_Illumina_NA12878 ZY_BGI_NA12878 ZY_GeneMind_NA12878; do
  case "$SAMPLE" in
    ZY_Illumina_*) PLATFORM=illumina_novaseq6000;;
    ZY_BGI_*)      PLATFORM=bgi_dnbseq_t7;;
    ZY_GeneMind_*) PLATFORM=genemind_surfseq5000;;
  esac
  BAM="$WORKDIR/bam/$SAMPLE.grch38.bwa-mem.full.s0/dedup.bam"
  if [ "$DRY_RUN" != 1 ] && [ ! -s "$BAM" ]; then
    echo "  !! full-depth BAM missing: $BAM  (run ALIGN_ONLY=1 run_one.sh $SAMPLE.grch38.bwa-mem.gatk4_hc.cpu.full.s0 first)" >&2
    continue
  fi
  PFX="$QCDIR/${SAMPLE}.vbid2"
  run "$(vbid2 "$BAM" "$PFX")"
  # .selfSM has a header row; FREEMIX is column 7 and #SNPS column 4 -- read the data row
  if [ "$DRY_RUN" != 1 ]; then
    FREEMIX=$(awk 'NR==2{print $7}' "$PFX.selfSM"); NSNP=$(awk 'NR==2{print $4}' "$PFX.selfSM")
    if [ -n "$VBID2_BIN" ]; then TOOLSRC="local:$(basename "$VBID2_BIN")"; else TOOLSRC="${VBID2_IMAGE##*/}"; fi
    echo "$SAMPLE,$PLATFORM,$FREEMIX,$NSNP,verifybamid2,$TOOLSRC,$(date +%F)" >> "$OUTCSV"
    echo "  $SAMPLE FREEMIX=$FREEMIX (#SNP=$NSNP)"
  fi
done
echo "=== done -> $OUTCSV ==="
