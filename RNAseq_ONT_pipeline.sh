#!/usr/bin/env bash
set -euo pipefail

# MJL with Claude Sonnet 4.6 2026-03-09

# Usage:
#   RNAseq_ONT_pipeline.sh <ont_fastq_dir> <outdir> <CPUs> <ref_genome.fa> <ref_annot.gtf|gff> <sample_sheet.csv> [--direct_rna] [--ref_transcriptome transcripts.fa] [--run_trinity]
#
# Notes:
#   - <sample_sheet.csv> for wf-transcriptomes must have columns: barcode,alias,condition
#   - Use --direct_rna if your library is ONT direct RNA; omit for cDNA
#   - If you provide --ref_transcriptome, wf-transcriptomes will also use it
#
# Tools assumed: toulligqc, nextflow(+singularity/docker), Rscript (for downstream)
# Trinity is optional and mostly useful for Illumina short reads; keep off for ONT-only runs.

if [[ $# -lt 6 ]]; then
  echo "Usage: $(basename "$0") <ont_fastq_dir> <outdir> <CPUs> <ref_genome.fa> <ref_annot.gtf|gff> <sample_sheet.csv> [--direct_rna] [--ref_transcriptome transcripts.fa] [--run_trinity]"
  exit 1
fi

ONT_FASTQ_DIR="$1"
OUTDIR="$2"
CPUS="$3"
REF_GENOME="$4"
REF_ANNOT="$5"
SAMPLE_SHEET="$6"
shift 6 || true

WF_LIBFLAG=""
REF_TRANSCRIPTOME=""
RUN_TRINITY=false

while (( "$#" )); do
  case "$1" in
    --direct_rna) WF_LIBFLAG="--direct_rna"; shift;;
    --ref_transcriptome) REF_TRANSCRIPTOME="$2"; shift 2;;
    --run_trinity) RUN_TRINITY=true; shift;;
    *) echo "Unknown flag: $1"; exit 1;;
  esac
done

mkdir -p "$OUTDIR"/{logs,toulligqc,wf-transcriptomes,downstream}

# 0) QC with ToulligQC (optional but recommended for ONT)  [4](https://github.com/GenomiqueENS/toulligQC)
if command -v toulligqc >/dev/null 2>&1; then
  echo "[ToulligQC] QC on $ONT_FASTQ_DIR"
  toulligqc --input "$ONT_FASTQ_DIR" --output "$OUTDIR/toulligqc" --threads "$CPUS" \
    || echo "[ToulligQC] QC returned non-zero status (continuing)."
else
  echo "[ToulligQC] not found; skipping."
fi

# 1) wf-transcriptomes: reference-guided assembly + DE  [1](https://epi2me.nanoporetech.com/epi2me-docs/workflows/wf-transcriptomes/)
nextflow pull epi2me-labs/wf-transcriptomes
MINIMAP2_OPTS="-k 15"
EXTRA_TRANSCRIPTOME_ARGS=()
if [[ -n "$REF_TRANSCRIPTOME" ]]; then
  EXTRA_TRANSCRIPTOME_ARGS+=(--ref_transcriptome "$REF_TRANSCRIPTOME")
fi

NXF_OPTS="-Xms2g -Xmx6g" nextflow run epi2me-labs/wf-transcriptomes \
  -profile singularity \
  --fastq "$ONT_FASTQ_DIR" \
  --ref_genome "$REF_GENOME" \
  --ref_annotation "$REF_ANNOT" \
  --sample_sheet "$SAMPLE_SHEET" \
  --minimap2_index_opts "$MINIMAP2_OPTS" \
  --de_analysis \
  $WF_LIBFLAG \
  "${EXTRA_TRANSCRIPTOME_ARGS[@]}" \
  --out_dir "$OUTDIR/wf-transcriptomes" \
  -work-dir "$OUTDIR/work" | tee "$OUTDIR/logs/wf-transcriptomes.log"

# 2) (Optional) Trinity de novo assembly (recommended for Illumina; not typical for ONT)  [8](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
if $RUN_TRINITY; then
  echo "[Trinity] WARNING: Trinity is optimized for Illumina; ensure inputs are suitable."
  mkdir -p "$OUTDIR/trinity"
  # Example: supply left/right reads if you have Illumina data
  # Trinity --seqType fq --left left.fq.gz --right right.fq.gz --CPU "$CPUS" --max_memory 100G --output "$OUTDIR/trinity"
fi

# 3) Downstream R analysis (volcano, GO, PCA/UMAP)
R_SCRIPT="$OUTDIR/downstream/longread_downstream.R"
cat > "$R_SCRIPT" <<'RSCRIPT'
# ===========================
# Downstream analysis for ONT
# Inputs: wf-transcriptomes outputs (DE_final & counts)
# Produces: volcano plots (EnhancedVolcano), GO enrichment (topGO), PCA/UMAP
# ===========================

install_if_missing <- function(pkg, bioc=FALSE){
  if (!require(pkg, character.only=TRUE, quietly=TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE)
    } else install.packages(pkg)
    library(pkg, character.only=TRUE)
  }
}

install_if_missing("BiocManager")
install_if_missing("EnhancedVolcano", bioc=TRUE)  # [9](https://bioconductor.org/packages//release/bioc/html/EnhancedVolcano.html)
install_if_missing("topGO", bioc=TRUE)            # [10](https://bioconductor.org/packages//release/bioc/vignettes/topGO/inst/doc/topGO_manual.html)
install_if_missing("ggplot2")
install_if_missing("data.table")
install_if_missing("Matrix")
install_if_missing("uwot")                        # for UMAP like refine.bio example  [11](https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html)

suppressPackageStartupMessages({
  library(EnhancedVolcano)
  library(topGO)
  library(ggplot2)
  library(data.table)
  library(Matrix)
  library(uwot)
})

# Auto-detect wf-transcriptomes outputs relative to script:
OUTDIR <- normalizePath(file.path("..","wf-transcriptomes"))
DE_DIR <- file.path(OUTDIR, "DE_final")
if (!dir.exists(DE_DIR)) stop("Cannot find DE_final under wf-transcriptomes output.")

# Expect edgeR-like tables in DE_final and count tables (gene_counts.tsv / transcript_counts.tsv)
# v0.3.0 notes mention count tables are included.  [13](https://zenodo.org/records/16756575)
gene_counts_file <- list.files(DE_DIR, pattern="gene_count.*\\.tsv$", full.names=TRUE)
tx_counts_file   <- list.files(DE_DIR, pattern="transcript_count.*\\.tsv$", full.names=TRUE)
if (length(gene_counts_file)==0) gene_counts_file <- list.files(DE_DIR, pattern="gene.*counts.*\\.tsv$", full.names=TRUE)
if (length(tx_counts_file)==0)   tx_counts_file   <- list.files(DE_DIR, pattern="transcript.*counts.*\\.tsv$", full.names=TRUE)

# Load count matrices if present
gene_counts <- tryCatch({ fread(gene_counts_file[1]) }, error=function(e) NULL)
tx_counts   <- tryCatch({ fread(tx_counts_file[1]) },   error=function(e) NULL)

# Load DE results (assume one or multiple contrasts)
de_tables <- list.files(DE_DIR, pattern="edgeR.*\\.tsv$|DE_results.*\\.tsv$|_DE_.*\\.tsv$", full.names=TRUE)
if (length(de_tables)==0) {
  message("Could not find DE tables; listing DE_final/ for inspection")
  print(list.files(DE_DIR))
} else {
  dir.create(file.path(OUTDIR, "plots"), showWarnings=FALSE)
  for (f in de_tables) {
    dt <- fread(f)
    # Heuristics: look for columns logFC and PValue or FDR
    xcol <- if ("logFC" %in% names(dt)) "logFC" else names(dt)[grep("logFC", names(dt))[1]]
    ycol <- if ("PValue" %in% names(dt)) "PValue" else names(dt)[grep("P", names(dt), ignore.case=TRUE)[1]]

    png(file.path(OUTDIR, "plots", paste0(basename(f), "_volcano.png")), width=1600, height=1200, res=200)
    EnhancedVolcano(dt,
      lab = dt[[1]],
      x = xcol, y = ycol,
      pCutoff = 0.05, FCcutoff = 1,
      title = paste0("Volcano: ", basename(f))
    )
    dev.off()
  }
}

# PCA + UMAP from normalized counts if available  (pattern inspired by refine.bio examples)  [11](https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html)
if (!is.null(gene_counts)) {
  mat <- as.matrix(gene_counts[,-1])
  rownames(mat) <- gene_counts[[1]]
  # simple logCPM-like transform
  mat_log <- log2(sweep(mat, 2, colSums(mat)/1e6, "/") + 1)
  mat_log <- mat_log[apply(mat_log,1,var)>0,]
  # PCA
  pca <- prcomp(t(mat_log), scale.=TRUE)
  df_pca <- as.data.frame(pca$x[,1:2])
  ggplot(df_pca, aes(PC1, PC2)) + geom_point() + ggtitle("PCA (logCPM)") +
    ggsave(filename = file.path(OUTDIR, "plots", "pca_logcpm.png"), width=8, height=6, dpi=200)

  # UMAP
  set.seed(1)
  emb <- umap(t(mat_log), n_neighbors=15, min_dist=0.2, metric="cosine")
  df_umap <- as.data.frame(emb); names(df_umap) <- c("UMAP1","UMAP2")
  ggplot(df_umap, aes(UMAP1, UMAP2)) + geom_point() + ggtitle("UMAP (logCPM)") +
    ggsave(filename = file.path(OUTDIR, "plots", "umap_logcpm.png"), width=8, height=6, dpi=200)
}

# GO enrichment with topGO (requires gene2GO mapping for your organism)  [10](https://bioconductor.org/packages//release/bioc/vignettes/topGO/inst/doc/topGO_manual.html)
# Provide a list of DE genes (example: from the first DE table)
if (length(de_tables)>0) {
  dt <- fread(de_tables[1])
  # Expect gene ID in first column; define 'interesting' by FDR < 0.05
  id_col <- names(dt)[1]
  fdr_col <- if ("FDR" %in% names(dt)) "FDR" else names(dt)[grep("FDR|adj", names(dt), ignore.case=TRUE)[1]]
  sig <- dt[[id_col]][ dt[[fdr_col]] < 0.05 ]

  # Load a gene2GO mapping file (tab: gene<tab>GO:xxxx;GO:yyyy)
  gene2go_path <- Sys.getenv("GENE2GO_TSV", unset = "")
  if (nzchar(gene2go_path) && file.exists(gene2go_path)) {
    m <- fread(gene2go_path, header=FALSE, sep="\t")
    allGenes <- unique(m[[1]])
    geneList <- factor(as.integer(allGenes %in% sig)); names(geneList) <- allGenes

    # build topGO object (BP ontology)
    GOdata <- new("topGOdata", ontology="BP", allGenes=geneList,
                  annot=annFUN.gene2GO, gene2GO = split(strsplit(m[[2]], ";"), seq_len(nrow(m))))
    resultFisher <- runTest(GOdata, algorithm="weight01", statistic="fisher")
    tab <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 30)
    fwrite(tab, file.path(OUTDIR, "plots", "topGO_BP_weight01.tsv"), sep="\t")
  } else {
    message("Set GENE2GO_TSV env var to a gene2GO mapping for GO enrichment.")
  }
}

message("Downstream analysis complete.")
RSCRIPT

echo "[R] Running downstream analysis..."
Rscript "$R_SCRIPT" | tee "$OUTDIR/logs/downstream_R.log"

echo "Pipeline complete. See: $OUTDIR"
