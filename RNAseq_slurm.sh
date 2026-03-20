#!/bin/bash

# MJL with Claude Sonnet 4.6 2026-03-09

#SBATCH --account=e33134
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --job-name=ONT_NanoPlot_EPI2ME
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

# =============================
# USER CONFIG (edit paths)
# =============================

FASTQ_FILE="/gpfs/home/mjl7892/DDPS/RNAseq_final_project/Coleus_amboinicus/Coleus_amboinicus_D11_03-09-2026_calls.fastq"

# C. barbatus reference (genome + annotation)
REF_GENOME="/gpfs/home/mjl7892/DDPS/RNAseq_final_project/Coleus_barbatus_test/PLBA_531000000.1.fasta"
REF_ANNOT_GFF="/gpfs/home/mjl7892/DDPS/RNAseq_final_project/Coleus_barbatus_test/PLBA_chr_aug_preds.gff"
REF_ANNOT_GTF=""    # leave empty to auto-convert GFF → GTF (recommended)

# ---- Overrides / options ----
# Force ONT protocol: export ONT_PROTOCOL_OVERRIDE=directRNA   # or cDNA
# Add extra minimap2 index flags (default tuned for ONT):
MM2_INDEX_OPTS="${MM2_INDEX_OPTS:--k 15}"

# Run EPI2ME DE analysis section? (needs sample sheet for multi-sample; OK single)
DO_DE="${DO_DE:-0}"   # 1 to enable: passes --de_analysis

# If you want to export matrices post-run via your attached prepDE.py:
RUN_PREPDE="${RUN_PREPDE:-1}"

# --- Unique output dir ---
RUN_TS=$(date +%Y-%m-%d_%H%M%S)
OUTDIR="$HOME/DDPS/RNAseq_final_project/${RUN_TS}_ont_nanoplot_epi2me"
THREADS=${SLURM_CPUS_PER_TASK:-16}

# =============================
# ENV / MODULES
# =============================
set -Eeuo pipefail
trap 'echo "[trap] Error on line $LINENO (exit=$?)" >&2' ERR

module purge
module unload nextflow 2>/dev/null || true
module load singularity || module load apptainer || true
module load gffread || true
module load pigz || true

mkdir -p "$OUTDIR"/{logs,ref,fastq,nanoplot,wf_transcriptomes,bin,singularity_cache}

# tee logs (overwrite, do not append)
exec > >(tee "$OUTDIR/logs/slurm.stdout") 2> >(tee "$OUTDIR/logs/slurm.stderr" >&2)
set -x

echo "RUN_ID: ${RUN_TS}  HOST: $(hostname)  THREADS=$THREADS"
echo "SCRIPT: $(readlink -f "$0")"
md5sum "$(readlink -f "$0")" || true

# -----------------------------
# Ensure Java 17+ for Nextflow
# -----------------------------
ensure_java17() {
  module load java/21 2>/dev/null || module load java/17 2>/dev/null || module load java/jdk17 2>/dev/null || true
  local maj; maj=$(java -version 2>&1 | awk -F\" '/version/ {split($2,a,"."); print a[1]}' || echo 0)
  if [[ -n "$maj" && "$maj" -ge 17 ]]; then
    export JAVA_HOME="$(dirname "$(dirname "$(readlink -f "$(command -v java)")")")"
    export JAVA_CMD="${JAVA_HOME}/bin/java"
    java -version
    return 0
  fi
  local JDK_DIR="$OUTDIR/jdk17"; mkdir -p "$JDK_DIR"
  local JDK_TAR="$JDK_DIR/temurin17.tar.gz"
  local URLS=(
    "https://api.adoptium.net/v3/binary/latest/17/ga/linux/x64/jdk/hotspot/normal/eclipse?project=jdk"
    "https://github.com/adoptium/temurin17-binaries/releases/latest/download/OpenJDK17U-jdk_x64_linux_hotspot.tar.gz"
  )
  rm -f "$JDK_TAR"
  for U in "${URLS[@]}"; do
    curl -fsSL -o "$JDK_TAR" "$U" || true
    [[ -s "$JDK_TAR" ]] && tar -tzf "$JDK_TAR" >/dev/null 2>&1 && break || rm -f "$JDK_TAR"
  done
  tar -xzf "$JDK_TAR" -C "$JDK_DIR"
  export JAVA_HOME="$(find "$JDK_DIR" -maxdepth 1 -type d -name 'jdk-*' | head -n1)"
  export PATH="$JAVA_HOME/bin:$PATH"
  export NXF_JAVA_HOME="$JAVA_HOME"
  java -version
}
ensure_java17

# --------------------------------------
# Ensure Nextflow ≥ 23.04.2
# --------------------------------------
NF_MIN_MAJOR=23; NF_MIN_MINOR=4
NFV=$(nextflow -version 2>/dev/null | sed -n 's/.*version \([0-9]*\)\.\([0-9]*\).*/\1 \2/p' || true)
if [[ -z "$NFV" ]]; then
  curl -s https://get.nextflow.io | bash
  mv nextflow "$OUTDIR/bin/nextflow"; export PATH="$OUTDIR/bin:$PATH"
else
  NF_MAJOR=$(echo "$NFV" | awk '{print $1}'); NF_MINOR=$(echo "$NFV" | awk '{print $2}')
  if (( NF_MAJOR < NF_MIN_MAJOR || (NF_MAJOR==NF_MIN_MAJOR && NF_MINOR < NF_MIN_MINOR) )); then
    curl -s https://get.nextflow.io | bash
    mv nextflow "$OUTDIR/bin/nextflow"; export PATH="$OUTDIR/bin:$PATH"
  fi
fi
NXF_BIN="$(command -v nextflow)"
"${NXF_BIN}" -version

# Singularity cache
export NXF_SINGULARITY_CACHEDIR="$OUTDIR/singularity_cache"

# =============================
# 1) Stage FASTQ
# =============================
if [[ ! -s "$FASTQ_FILE" ]]; then
  echo "ERROR: FASTQ not found: $FASTQ_FILE" >&2; exit 1
fi
FQ_BASENAME="$(basename "$FASTQ_FILE")"
STAGED_FASTQ="$OUTDIR/fastq/$FQ_BASENAME"
ln -s "$FASTQ_FILE" "$STAGED_FASTQ" 2>/dev/null || cp -a "$FASTQ_FILE" "$STAGED_FASTQ"
if [[ "$STAGED_FASTQ" != *.gz ]]; then
  pigz -p "$THREADS" -f "$STAGED_FASTQ"; STAGED_FASTQ="${STAGED_FASTQ}.gz"
fi
echo "FASTQ staged -> $STAGED_FASTQ"

# =============================
# 2) NanoPlot QC (ONT)
# =============================
run_nanoplot() {
  local fq="$1"; local out="$2"; local threads="$3"
  mkdir -p "$out"
  if command -v NanoPlot >/dev/null 2>&1; then
    NanoPlot --fastq "$fq" -o "$out" --threads "$threads" --tsv_stats --title "NanoPlot QC"
  else
    # BioContainers NanoPlot image via Singularity
    local IMG="docker://quay.io/biocontainers/nanoplot:1.46.2--pyhdfd78af_0"
    singularity exec "$IMG" NanoPlot --fastq "$fq" -o "$out" --threads "$threads" --tsv_stats --title "NanoPlot QC"
  fi
}
run_nanoplot "$STAGED_FASTQ" "$OUTDIR/nanoplot" "$THREADS"

# =============================
# 3) Convert GFF → GTF (patch AUGUSTUS GFF2 bare-ID format first)
# =============================
if [[ -z "$REF_ANNOT_GTF" ]]; then
  GTF_PATH="$OUTDIR/ref/$(basename "${REF_ANNOT_GFF%.*}").gtf"
  GFF_PATCHED="$OUTDIR/ref/$(basename "${REF_ANNOT_GFF%.*}").patched.gff"

  # AUGUSTUS GFF2 uses a bare transcript ID in column 9 (e.g. "g1.t1") with no
  # key=value pairs. gffread cannot parse this and silently drops all records,
  # producing an empty GTF. We rewrite col 9 into proper key=value format first.
  echo "Patching AUGUSTUS GFF2 bare-ID attributes → $GFF_PATCHED"
  awk -v OFS='\t' '
    BEGIN { FS=OFS="\t" }
    /^#/ { print; next }
    NF < 9 { print; next }
    {
      id = $9
      gsub(/^[ \t]+|[ \t]+$/, "", id)
      # Only rewrite if not already key=value style
      if (id !~ /=/ && id !~ /"/) {
        gene_id = id
        sub(/\.[^.]+$/, "", gene_id)
        $9 = "gene_id \"" gene_id "\"; transcript_id \"" id "\";"
      }
      print
    }
  ' "$REF_ANNOT_GFF" > "$GFF_PATCHED"

  echo "Converting patched GFF → GTF → $GTF_PATH"
  gffread -T "$GFF_PATCHED" -o "$GTF_PATH"

  # Sanity check: GTF must have transcript_id attributes
  TX_COUNT=$(grep -c 'transcript_id' "$GTF_PATH" 2>/dev/null || echo 0)
  if [[ "$TX_COUNT" -eq 0 ]]; then
    echo "ERROR: GTF has no transcript_id records after patching." >&2
    echo "       Inspect: $GFF_PATCHED" >&2
    exit 1
  fi
  echo "GTF sanity check passed: $TX_COUNT lines with transcript_id"
  echo "Sample GTF record (col 9):"
  grep -m1 $'\ttranscript\t' "$GTF_PATH" | cut -f9 || true
else
  GTF_PATH="$REF_ANNOT_GTF"
fi

# =============================
# 4) ONT protocol guess (can override)
# =============================
guess_ont_protocol() {
  local header
  header="$(zcat -f "$1" | head -n1 | tr '[:upper:]' '[:lower:]' || true)"
  if echo "$header" | grep -Eq 'rna00[2-9]|rna0[1-9]'; then echo "directRNA"; else echo "cDNA"; fi
}
ONT_PROTOCOL="${ONT_PROTOCOL_OVERRIDE:-$(guess_ont_protocol "$STAGED_FASTQ")}"
echo "ONT protocol: $ONT_PROTOCOL"

# =============================
# 5) Run EPI2ME Labs wf-transcriptomes (reference-guided)
# =============================
SAMPLE_NAME="${FQ_BASENAME%%.*}"
WF_OUT="$OUTDIR/wf_transcriptomes"
WF_REVISION="${WF_REVISION:-v1.7.2}"

# ---------------------------------------------------------------------------
# CONTAINER STRATEGY:
#
# The SHA-tagged image that v1.7.x's nextflow.config references (shaaaf20...)
# has a broken manifest on Docker Hub — Singularity cannot pull it via the
# SHA tag. However, 'ontresearch/wf-transcriptomes:latest' points to the
# IDENTICAL image (same full 888MB release with pysam/pychopper/etc.) and
# HAS a valid, pullable manifest.
#
# Fix: Pull via 'latest' tag (valid manifest, correct full image), but save
# the SIF under the SHA-tagged filename Nextflow expects. Nextflow finds the
# file as a cache hit and never tries to pull the broken tag itself.
#
# We first pull the Nextflow DSL, then read the exact SHA from nextflow.config
# so the SIF filename always matches, regardless of future version bumps.
# ---------------------------------------------------------------------------

SIF_DIR="$OUTDIR/singularity_cache"
export NXF_SINGULARITY_CACHEDIR="$SIF_DIR"

# Pull the Nextflow DSL first so we can read the exact SHA from its config
"${NXF_BIN}" pull epi2me-labs/wf-transcriptomes -r "$WF_REVISION"

# Read the exact container SHAs that this version's nextflow.config specifies
NXF_CFG="$HOME/.nextflow/assets/epi2me-labs/wf-transcriptomes/nextflow.config"
if [[ ! -f "$NXF_CFG" ]]; then
  echo "ERROR: Cannot find pulled nextflow.config at $NXF_CFG" >&2; exit 1
fi
echo "[Step 5] Container SHAs from $NXF_CFG:"
grep -E 'container_sha|common_sha' "$NXF_CFG" || true

WF_SHA=$(grep    'container_sha\b' "$NXF_CFG" | grep -v 'common' | grep -oP '(?<=")[^"]+' | head -1)
COMMON_SHA=$(grep 'common_sha\b'   "$NXF_CFG" | grep -oP '(?<=")[^"]+' | head -1)

# Fallback: grep for the sha pattern directly
[[ -z "$WF_SHA" ]]     && WF_SHA=$(grep     'container_sha' "$NXF_CFG" | grep -v common | grep -oE 'sha[a-f0-9]+' | head -1)
[[ -z "$COMMON_SHA" ]] && COMMON_SHA=$(grep  'common_sha'   "$NXF_CFG" | grep -oE 'sha[a-f0-9]+' | head -1)

if [[ -z "$WF_SHA" || -z "$COMMON_SHA" ]]; then
  echo "ERROR: Could not parse container SHAs. Relevant config lines:" >&2
  grep -E 'sha' "$NXF_CFG" | head -10 >&2; exit 1
fi
echo "[Step 5] wf_sha=$WF_SHA  common_sha=$COMMON_SHA"

# Nextflow SIF naming: replace / and : with -, append .img
sif_name() { local s="${1//\//-}"; echo "${s//:/-}.img"; }
SIF_WF="$SIF_DIR/$(sif_name "ontresearch/wf-transcriptomes:${WF_SHA}")"
SIF_COMMON="$SIF_DIR/$(sif_name "ontresearch/wf-common:${COMMON_SHA}")"
echo "[Step 5] Target SIF files:"
echo "  WF:     $SIF_WF"
echo "  Common: $SIF_COMMON"

build_sif() {
  local pull_tag="$1"  # tag to actually download (working manifest)
  local sif="$2"       # destination filename (may differ — uses SHA-tagged name)
  if [[ -s "$sif" ]]; then
    echo "[SIF] Already cached: $(basename "$sif")"; return 0
  fi

  # Primary: skopeo → local OCI dir → singularity build (avoids Docker Hub pull path entirely)
  if command -v skopeo >/dev/null 2>&1; then
    local oci_tmp; oci_tmp="$(mktemp -d "$SIF_DIR/oci_tmp_XXXXXX")"
    echo "[SIF] skopeo copy docker://$pull_tag → OCI → SIF"
    skopeo copy --override-os linux "docker://$pull_tag" "oci:${oci_tmp}:latest" \
      && singularity build --force "$sif" "oci:${oci_tmp}:latest" \
      && { rm -rf "$oci_tmp"; echo "[SIF] Built OK (skopeo): $(basename "$sif")"; return 0; }
    rm -rf "$oci_tmp"
    echo "[SIF] skopeo failed; trying --disable-cache..." >&2
  fi

  # Fallback: singularity build --disable-cache (fresh fetch, different code path)
  echo "[SIF] singularity build --disable-cache docker://$pull_tag"
  SINGULARITY_DISABLE_CACHE=1 singularity build --force "$sif" "docker://$pull_tag" \
    && { echo "[SIF] Built OK (--disable-cache): $(basename "$sif")"; return 0; }

  echo "[SIF] ERROR: all methods failed for $pull_tag" >&2
  return 1
}

echo "[Step 5] Building SIFs (pulling 'latest' tag to avoid broken SHA manifest)..."
# Pull 'latest' (valid manifest = same image as sha-tagged) → save under SHA filename
build_sif "ontresearch/wf-transcriptomes:latest"   "$SIF_WF"
build_sif "ontresearch/wf-common:${COMMON_SHA}"    "$SIF_COMMON"
echo "[Step 5] SIF build complete:"
ls -lh "$SIF_WF" "$SIF_COMMON"

# custom.config: point processes directly at pre-built SIFs by absolute path.
# This prevents Nextflow from ever invoking 'singularity pull docker://...'
CUSTOM_CFG="$OUTDIR/custom.config"
cat > "$CUSTOM_CFG" <<EOF
singularity.pullTimeout = '2 h'
singularity.enabled     = true
process {
  withLabel: 'isoforms' {
    container = '${SIF_WF}'
  }
  withLabel: 'wf_common' {
    container = '${SIF_COMMON}'
  }
}
EOF
echo "[Step 5] custom.config:"
cat "$CUSTOM_CFG"

NXF_OPTS='-Xms3g -Xmx8g'
WF_ARGS=(
  --fastq "$STAGED_FASTQ"
  --sample "$SAMPLE_NAME"
  --ref_genome "$REF_GENOME"
  --ref_annotation "$GTF_PATH"
  --transcriptome_source reference-guided
  --minimap2_index_opts "$MM2_INDEX_OPTS"
  --out_dir "$WF_OUT"
)
if [[ "$ONT_PROTOCOL" == "directRNA" ]]; then WF_ARGS+=( --direct_rna ); fi
if (( DO_DE )); then WF_ARGS+=( --de_analysis ); fi

"${NXF_BIN}" run epi2me-labs/wf-transcriptomes \
  -r "$WF_REVISION" \
  -c "$CUSTOM_CFG" \
  -profile singularity \
  -work-dir "$OUTDIR/work_wf" \
  -ansi-log false \
  "${WF_ARGS[@]}" \
  |& tee "$OUTDIR/logs/wf_transcriptomes.log"

# =============================
# 6) Optional: export count matrices from WF outputs (prepDE)
# =============================
if (( RUN_PREPDE )); then
  WF_DE_DIR="$WF_OUT/DE_final"
  if [[ -d "$WF_DE_DIR" ]]; then
    # If your prepDE.py is alongside the script or in CWD, try it.
    PREPDE="$(command -v python)"; PYTHON_BIN="${PREPDE:-python}"
    if [[ -f "./prepDE.py" ]]; then
      mkdir -p "$OUTDIR/postrun"
      $PYTHON_BIN ./prepDE.py --in_wf "$WF_OUT" --out_dir "$OUTDIR/postrun" --header_tsv "$OUTDIR/postrun/header.tsv" || true
      echo "prepDE outputs (if successful) in $OUTDIR/postrun"
    else
      echo "prepDE.py not found in current directory; skipping."
    fi
  else
    echo "DE_final/ not found under $WF_OUT (perhaps --de_analysis not run). Skipping prepDE."
  fi
fi

echo "DONE. Results:"
echo " - NanoPlot: $OUTDIR/nanoplot"
echo " - wf-transcriptomes: $WF_OUT"
