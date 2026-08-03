#!/usr/bin/env bash
# install_wrappers.sh -- expose the Docker-packaged tools as plain local commands, so that the bare
#   command names used by run_one.sh / run_happy.sh / calibrate.sh (samtools, gatk, bcftools,
#   mosdepth, seqtk, run_deepvariant, pbrun, hap.py) work unchanged. bwa is a local binary and is
#   simply symlinked.
#
# Each wrapper is one `docker run` layer, following the same conventions:
#   -i                     : attach stdin so pipes work (run_one has `bwa mem | samtools sort`)
#   -u $(id -u):$(id -g)   : run as the current user, so outputs are not owned by root
#   -v <mount>             : every path the container must see (reference, GIAB, fastq, WORKDIR)
#   -w $PWD                : keep the working directory consistent
#   --gpus all             : GPU tools only (pbrun)
#   pinned image tags      : reproducibility; these are the versions recorded in the manifest
#
# Usage (no root required; installs into ~/bin by default):
#   bash scripts/install_wrappers.sh
#   MOUNTS="/home/felix /data/test_data/felix /path/to/fastq" bash scripts/install_wrappers.sh
#   BIN=/somewhere/bin bash scripts/install_wrappers.sh
set -euo pipefail

BIN="${BIN:-$HOME/bin}"                                  # no root needed; make sure it is on PATH
# Mount points (space-separated): every top-level directory the containers must see, typically the
# one holding the reference/GIAB/bwa and the one holding WORKDIR.
# NOTE: if the fastq files live elsewhere, append their top-level directory to MOUNTS.
MOUNTS="${MOUNTS:-/home/felix /data/test_data/felix}"
BWA_BIN="${BWA_BIN:-/home/felix/bwa-0.7.15/bwa}"         # locally built bwa 0.7.15
mkdir -p "$BIN"

# assemble the -v arguments, one per mount point
VOPTS=""
for m in $MOUNTS; do VOPTS="$VOPTS -v \"$m:$m\""; done

# gen <command name> <extra docker args> <image:tag> <command inside the container>
gen() {
  local name=$1 extra=$2 image=$3 cmd=$4
  cat > "$BIN/$name" <<EOF
#!/usr/bin/env bash
exec docker run --rm -i -u "\$(id -u):\$(id -g)" $extra $VOPTS \\
  -w "\$PWD" $image $cmd "\$@"
EOF
  chmod +x "$BIN/$name"
  printf '  %-16s -> %s %s\n' "$name" "$image" "$cmd"
}

echo "installing wrappers into $BIN  (mounts: $MOUNTS):"

# ---- CPU tools ----
gen samtools        ""            staphb/samtools:1.23.1               samtools
gen bcftools        ""            staphb/bcftools:1.23.1               bcftools
gen gatk            ""            broadinstitute/gatk:4.3.0.0          gatk
gen mosdepth        ""            brentp/mosdepth:v0.3.3               ""    # the image entrypoint is already mosdepth, so the command is left empty
gen seqtk           ""            staphb/seqtk:1.5                     seqtk
gen run_deepvariant ""            google/deepvariant:1.9.0             run_deepvariant
gen fastp           ""            staphb/fastp:1.3.3                   fastp
gen multiqc         ""            multiqc/multiqc:v1.35                multiqc
gen hap.py          ""            mgibio/hap.py:v0.3.15                /opt/hap.py/bin/hap.py

# ---- GPU tool: Parabricks (GPU_DEV selects cards: empty = all, 0/1 = a single card,
#      which is how run_batch spreads jobs) ----
cat > "$BIN/pbrun" <<EOF
#!/usr/bin/env bash
GPUS="--gpus all"; [ -n "\${GPU_DEV:-}" ] && GPUS="--gpus device=\$GPU_DEV"
exec docker run --rm -i -u "\$(id -u):\$(id -g)" \$GPUS $VOPTS \\
  -w "\$PWD" nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1 pbrun "\$@"
EOF
chmod +x "$BIN/pbrun"; printf '  %-16s -> %s (GPU_DEV selects the card)\n' pbrun clara-parabricks:4.7.0-1

# ---- bwa: local binary, symlinked onto PATH ----
if [ -x "$BWA_BIN" ]; then
  ln -sf "$BWA_BIN" "$BIN/bwa"; printf '  %-16s -> %s (local binary)\n' bwa "$BWA_BIN"
else
  echo "  !! bwa not found ($BWA_BIN): skipping the symlink. Install it, then either run\n     ln -sf <bwa> $BIN/bwa manually or re-run this script with BWA_BIN set"
fi

cat <<EOF

Done. Make sure $BIN is on PATH:
  case ":\$PATH:" in *":$BIN:"*) ;; *) echo 'export PATH="$BIN:\$PATH"' >> ~/.bashrc && export PATH="$BIN:\$PATH";; esac
Quick self-check:
  bwa 2>&1 | head -3;               samtools --version | head -1
  bcftools --version | head -1;     gatk --version 2>&1 | grep -i 'GATK'
  mosdepth --version;               seqtk 2>&1 | grep -i version
  run_deepvariant --version 2>&1 | head -1;   hap.py --version
  pbrun version            # requires a visible GPU
Note: if pbrun reports insufficient shared memory, add --shm-size=... to the extra docker arguments
on the pbrun line of this script and reinstall.
EOF
