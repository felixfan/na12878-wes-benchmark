#!/usr/bin/env bash
# na12878.paths.template.sh -- data and resource locations for the benchmarking framework.
#
# Copy this file to na12878.paths.sh, replace every /PATH/TO/... placeholder with an absolute path on
# your system, then point HAPPY_CONFIG at it:
#     cp configs/na12878.paths.template.sh configs/na12878.paths.sh
#     $EDITOR configs/na12878.paths.sh
#     HAPPY_CONFIG=configs/na12878.paths.sh scripts/run_one.sh <run_id>
#
# The placeholder layout mirrors the system on which the released results were produced; any layout
# works as long as every variable resolves. Scripts source this file, so it must remain valid bash.

# ============ Raw fastq (three same-DNA NA12878 libraries, all Agilent V6, PE150) ============
# Naming convention: <sample>_R1.fq.gz / <sample>_R2.fq.gz
FQ_ILLUMINA_R1=/PATH/TO/FASTQ/ZY_Illumina_NA12878_R1.fq.gz
FQ_ILLUMINA_R2=/PATH/TO/FASTQ/ZY_Illumina_NA12878_R2.fq.gz
FQ_BGI_R1=/PATH/TO/FASTQ/ZY_BGI_NA12878_R1.fq.gz
FQ_BGI_R2=/PATH/TO/FASTQ/ZY_BGI_NA12878_R2.fq.gz
FQ_GM_R1=/PATH/TO/FASTQ/ZY_GeneMind_NA12878_R1.fq.gz
FQ_GM_R2=/PATH/TO/FASTQ/ZY_GeneMind_NA12878_R2.fq.gz

# ============ Reference genomes (fasta + .fai + .dict + BWA index) ============
# The legacy axis is b37 (no chr prefix), which aligns natively with the GIAB GRCh37 truth set.
# Note this is b37, not ucsc.hg19.
REF_GRCH38=/PATH/TO/REFERENCE_BUNDLE/GRCh38/Homo_sapiens_assembly38.fasta   # uncompressed, with .fai/.dict and BWA index
REF_B37=/PATH/TO/REFERENCE_BUNDLE/GRCh37/human_g1k_v37.fasta                # uncompressed, with .fai/.dict and BWA 0.7.15 index

# ============ RTG SDF (required by the vcfeval engine; build with: rtg format -o REF.sdf REF.fasta) ============
SDF_GRCH38=/PATH/TO/REFERENCE_BUNDLE/GRCh38/Homo_sapiens_assembly38.sdf
SDF_B37=/PATH/TO/REFERENCE_BUNDLE/GRCh37/human_g1k_v37.sdf

# ============ BQSR known-sites (GATK-family callers only; DeepVariant does not use BQSR) ============
# The GRCh38 Broad references ship with .tbi indexes. The b37 legacy bundles do not, so build them
# once with: gatk IndexFeatureFile -I <file>
KNOWN_DBSNP_GRCH38=/PATH/TO/REFERENCE_BUNDLE/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz
KNOWN_MILLS_GRCH38=/PATH/TO/REFERENCE_BUNDLE/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
KNOWN_INDELS_GRCH38=/PATH/TO/REFERENCE_BUNDLE/GRCh38/Homo_sapiens_assembly38.known_indels.vcf.gz
KNOWN_DBSNP_B37=/PATH/TO/REFERENCE_BUNDLE/GRCh37/dbsnp_138.b37.vcf.gz
KNOWN_MILLS_B37=/PATH/TO/REFERENCE_BUNDLE/GRCh37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
KNOWN_INDELS_B37=/PATH/TO/REFERENCE_BUNDLE/GRCh37/1000G_phase1.indels.b37.vcf.gz

# ============ GIAB HG001 v4.2.1 truth set (VCF + high-confidence BED, one pair per build) ============
# Download from:
#   https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/
# The files in the GRCh37 directory already use b37 coordinates (no chr prefix), matching REF_B37,
# so no renaming is needed.
GIAB_VCF_GRCH38=/PATH/TO/GIAB/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
GIAB_BED_GRCH38=/PATH/TO/GIAB/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
GIAB_VCF_B37=/PATH/TO/GIAB/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
GIAB_BED_B37=/PATH/TO/GIAB/HG001_GRCh37_1_22_v4.2.1_benchmark.bed

# ============ Agilent SureSelect V6 target BED (one per build) ============
# V6 = Agilent design S07604514 (SureSelect Human All Exon V6 r2), the "Target Regions" track.
# It must match the kit actually used for library preparation (not V6+UTR S07604624 or
# V6+COSMIC S07604715). SureDesign offers both builds: the hg38 version can be used directly (no
# liftOver), while the hg19 version needs the chr prefix stripped to match b37.
# Track choice: "Covered" (~60.5 Mb, the probe-covered space including splice edges), rather than
# "Regions" (~38 Mb core exons) or "Padded" (~109 Mb, too wide).
# Preparation: drop the track/browser header lines and cut -f1-3; for hg19 additionally
# sed 's/^chr//'. The same BED is used throughout: calling -L, hap.py -T and mosdepth --by.
# A prepared copy of both builds ships with this repository under reference/.
V6_BED_B37=/PATH/TO/REFERENCE_BUNDLE/SureSelectHumanAllExonV6r2/hg19/S07604514_Covered.GRCh37.bed
V6_BED_GRCH38=/PATH/TO/REFERENCE_BUNDLE/SureSelectHumanAllExonV6r2/hg38/S07604514_Covered.GRCh38.bed

# ============ GA4GH genome stratifications (.tsv index, one per build; optional) ============
# Only the top-level difficult regions are used: LowComplexity, GCcontent, mappability,
# SegmentalDuplications, FunctionalRegions. Leave empty to run hap.py without stratification.
STRAT_GRCH38=              # GRCh38 stratification .tsv
STRAT_B37=                 # GRCh37/b37 stratification .tsv (no-chr version)

# ============ Note on the depth titration (no data volume needs to be entered here) ============
# Downsampling works in on-target-X space: fraction = target_X / X_full.
# X_full is measured by Phase-0 calibration (scripts/calibrate.sh), written to
# $WORKDIR/plan/calibration.csv and read automatically by run_one.sh, so there is no need to enter
# raw gigabases here. See EXPERIMENT-DESIGN.md section 2.2.

# ============ GATK3 legacy jar (java -jar; only for the optional gatk3 caller) ============
GATK3_JAR=                 # absolute path to GenomeAnalysisTK.jar (3.8-1), if used

# ============ Output / working directory ============
WORKDIR=/PATH/TO/WORKDIR   # root for intermediate BAM/VCF files and hap.py output

# ============ Tool versions (recorded in the manifest for provenance) ============
PARABRICKS_VERSION=4.7.0   # nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1; bundles DeepVariant 1.9.0 + GATK 4.3.0.0
BWA_VERSION=0.7.15         # locally built binary
SAMTOOLS_VERSION=1.23.1    # staphb/samtools:1.23.1
GATK4_VERSION=4.3.0.0      # CPU broadinstitute/gatk:4.3.0.0, matching the GATK engine inside pbrun haplotypecaller
DEEPVARIANT_VERSION=1.9.0  # CPU google/deepvariant:1.9.0, matching the model inside pbrun deepvariant
GATK3_VERSION=3.8-1        # only if the optional gatk3 caller is used
PICARD_VERSION=2.27.5      # bundled in the GATK 4.3.0.0 image; MarkDuplicates runs as `gatk MarkDuplicates`
SEQTK_VERSION=1.5          # staphb/seqtk:1.5; downsampling via `sample -s<seed> <fraction>`
BCFTOOLS_VERSION=1.23.1    # staphb/bcftools:1.23.1; used for norm and index
HAPPY_VERSION=0.3.15       # mgibio/hap.py:v0.3.15; hap.py lives in /opt/hap.py/bin and is not on PATH,
                           # so the wrapper invokes it by full path
RTG_VERSION=3.12.1         # bundled inside the hap.py image (/opt/hap.py/libexec/rtg-tools-install/rtg)
MOSDEPTH_VERSION=0.3.3     # brentp/mosdepth:v0.3.3; run with --by <V6 BED> --thresholds
FASTP_VERSION=1.3.3        # staphb/fastp:1.3.3; raw fastq QC (Q30, GC, adapters, yield)
MULTIQC_VERSION=1.35       # multiqc/multiqc:v1.35; QC aggregation
