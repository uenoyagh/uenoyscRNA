suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

VERSION <- "v6.9.6.2"

BASE <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"
MON_RDS <- file.path(BASE, "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.5/objects/Mouse_MASH_Monocyte_state_module_scored_v6.9.5.rds")
WHOLE_RDS <- file.path(BASE, "Mouse_MASH_RDS/RDS3_annotation_visualization_v4.1.1/objects/RDS3_with_visualization_metadata_v4.1.1.rds")
DE_DIR <- file.path(BASE, "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.6/tables")
ANNOTATION_COL <- "celltype_for_R8plot_FIXED2"
OUTDIR <- file.path(BASE, paste0("Mouse_MASH_RDS/Mouse_MASH_Monocyte_", VERSION))
TABDIR <- file.path(OUTDIR, "tables")
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte genome-wide reference specificity\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(MON_RDS)) stop("Missing Monocyte RDS: ", MON_RDS)
if (!file.exists(WHOLE_RDS)) stop("Missing whole-cell RDS: ", WHOLE_RDS)

resolve_sample_col <- function(obj) {
  candidates <- c("sample", "sample_id", "Sample", "orig.ident")
  out <- candidates[candidates %in% colnames(obj@meta.data)][1]
  if (is.na(out)) stop("Could not resolve sample column.")
  out
}

get_counts <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  assay <- obj[["RNA"]]
  if (inherits(assay, "Assay5")) {
    count_layers <- grep("^counts", Layers(assay), value=TRUE)
    if (length(count_layers) == 0) stop("No RNA counts layer.")
    if (length(count_layers) > 1) {
      obj[["RNA"]] <- JoinLayers(obj[["RNA"]], layers=count_layers, new="counts")
    }
  }
  GetAssayData(obj, assay="RNA", layer="counts")
}

aggregate_logcpm <- function(mat, samples, genes, sample_order) {
  genes <- intersect(genes, rownames(mat))
  out <- matrix(NA_real_, nrow=length(genes), ncol=length(sample_order),
                dimnames=list(genes, sample_order))
  for (s in sample_order) {
    cells <- colnames(mat)[samples == s]
    if (length(cells) == 0) next
    pb <- Matrix::rowSums(mat[genes, cells, drop=FALSE])
    lib <- sum(Matrix::rowSums(mat[, cells, drop=FALSE]))
    out[genes, s] <- log2((1e6 * pb / max(lib, 1)) + 0.5)
  }
  out
}

row_median_na <- function(x) apply(x, 1, median, na.rm=TRUE)

samples_all6 <- c("STD_rep1", "CDHFD_rep1", "Sham1", "Sham20", "Tx17", "Tx5")
samples_shamtx <- c("Sham1", "Sham20", "Tx17", "Tx5")

mon <- readRDS(MON_RDS)
mon_counts <- get_counts(mon)
mon_sample_col <- resolve_sample_col(mon)
mon_sample <- as.character(mon@meta.data[[mon_sample_col]])
genes <- rownames(mon_counts)
cat("Genome-wide audit genes:", length(genes), "\n")

mon_all6 <- aggregate_logcpm(mon_counts, mon_sample, genes, samples_all6)
mon_shamtx <- aggregate_logcpm(mon_counts, mon_sample, genes, samples_shamtx)

whole <- readRDS(WHOLE_RDS)
if (!(ANNOTATION_COL %in% colnames(whole@meta.data))) stop("Missing annotation column: ", ANNOTATION_COL)
labels <- as.character(whole@meta.data[[ANNOTATION_COL]])
label_table <- sort(table(labels), decreasing=TRUE)
cat("\n=== AVAILABLE WHOLE-CELL LABELS ===\n")
print(label_table)
write.csv(data.frame(label=names(label_table), n_cells=as.integer(label_table), stringsAsFactors=FALSE),
          file.path(TABDIR, "Whole_cell_annotation_labels_v6.9.6.2.csv"), row.names=FALSE)

reference_labels <- c(
  Hepatocyte="Hepatocyte",
  Kupffer_Macrophage="Kupffer_Macrophage",
  HSC_Mesenchymal="HSC_Mesenchymal",
  LSEC="LSEC",
  Cholangiocyte="Cholangiocyte",
  Neutrophil="Neutrophil",
  Dendritic="Dendritic"
)
reference_labels <- reference_labels[reference_labels %in% unique(labels)]
if (length(reference_labels) == 0) stop("No reference lineages found.")
whole_sample_col <- resolve_sample_col(whole)

reference_all6 <- list()
reference_shamtx <- list()
reference_counts <- list()

for (ref_name in names(reference_labels)) {
  label <- reference_labels[[ref_name]]
  cells <- rownames(whole@meta.data)[labels == label]
  ref <- subset(whole, cells=cells)
  ref_counts <- get_counts(ref)
  ref_sample <- as.character(ref@meta.data[[whole_sample_col]])
  reference_all6[[ref_name]] <- aggregate_logcpm(ref_counts, ref_sample, genes, samples_all6)
  reference_shamtx[[ref_name]] <- aggregate_logcpm(ref_counts, ref_sample, genes, samples_shamtx)
  reference_counts[[ref_name]] <- data.frame(reference=ref_name, annotation_label=label, n_cells=ncol(ref), stringsAsFactors=FALSE)
  cat("Reference:", ref_name, "| cells:", ncol(ref), "\n")
  rm(ref, ref_counts); gc()
}

reference_counts <- do.call(rbind, reference_counts)
rownames(reference_counts) <- NULL
write.csv(reference_counts, file.path(TABDIR, "Reference_lineage_cell_counts_v6.9.6.2.csv"), row.names=FALSE)

comparison <- data.frame(
  gene=genes,
  Monocyte_median_logCPM_ShTx=row_median_na(mon_shamtx)[genes],
  Monocyte_median_logCPM_all6=row_median_na(mon_all6)[genes],
  stringsAsFactors=FALSE
)

for (ref_name in names(reference_labels)) {
  ref_st <- row_median_na(reference_shamtx[[ref_name]])
  ref_a6 <- row_median_na(reference_all6[[ref_name]])
  comparison[[paste0(ref_name, "_median_logCPM_ShTx")]] <- ref_st[comparison$gene]
  comparison[[paste0(ref_name, "_minus_Monocyte_ShTx")]] <- comparison[[paste0(ref_name, "_median_logCPM_ShTx")]] - comparison$Monocyte_median_logCPM_ShTx
  comparison[[paste0(ref_name, "_median_logCPM_all6")]] <- ref_a6[comparison$gene]
  comparison[[paste0(ref_name, "_minus_Monocyte_all6")]] <- comparison[[paste0(ref_name, "_median_logCPM_all6")]] - comparison$Monocyte_median_logCPM_all6
}

shamtx_delta_cols <- grep("_minus_Monocyte_ShTx$", colnames(comparison), value=TRUE)
all6_delta_cols <- grep("_minus_Monocyte_all6$", colnames(comparison), value=TRUE)
shamtx_delta_mat <- as.matrix(comparison[, shamtx_delta_cols, drop=FALSE])
all6_delta_mat <- as.matrix(comparison[, all6_delta_cols, drop=FALSE])

comparison$max_reference_delta_ShTx <- apply(shamtx_delta_mat, 1, max, na.rm=TRUE)
comparison$max_reference_lineage_ShTx <- sub("_minus_Monocyte_ShTx$", "", shamtx_delta_cols[apply(shamtx_delta_mat, 1, which.max)])
comparison$max_reference_delta_all6 <- apply(all6_delta_mat, 1, max, na.rm=TRUE)
comparison$max_reference_lineage_all6 <- sub("_minus_Monocyte_all6$", "", all6_delta_cols[apply(all6_delta_mat, 1, which.max)])

comparison$specificity_class_ShTx <- ifelse(
  comparison$max_reference_delta_ShTx > 2,
  "reference_lineage_dominant",
  ifelse(comparison$max_reference_delta_ShTx < -1, "Monocyte_enriched", "ambiguous_shared")
)
comparison$specificity_class_all6 <- ifelse(
  comparison$max_reference_delta_all6 > 2,
  "reference_lineage_dominant",
  ifelse(comparison$max_reference_delta_all6 < -1, "Monocyte_enriched", "ambiguous_shared")
)
comparison$keep_reference_aware <- comparison$max_reference_delta_ShTx <= 2
comparison$keep_strict_reference_sensitivity <- comparison$max_reference_delta_ShTx <= 1.5

de_files <- c(
  ALL=file.path(DE_DIR, "Monocyte_pseudobulk_ALL_Sham_vs_Tx_v6.9.6.csv"),
  NO_QCWATCH=file.path(DE_DIR, "Monocyte_pseudobulk_NO_QCWATCH_Sham_vs_Tx_v6.9.6.csv"),
  PRIMARY_CORE=file.path(DE_DIR, "Monocyte_pseudobulk_PRIMARY_CORE_Sham_vs_Tx_v6.9.6.csv")
)

for (fw in names(de_files)) {
  if (!file.exists(de_files[[fw]])) stop("Missing DE file: ", de_files[[fw]])
  de <- read.csv(de_files[[fw]], stringsAsFactors=FALSE)
  idx <- match(comparison$gene, de$gene)
  comparison[[paste0(fw, "_logFC")]] <- de$logFC[idx]
  comparison[[paste0(fw, "_PValue")]] <- de$PValue[idx]
  comparison[[paste0(fw, "_FDR")]] <- de$FDR[idx]
  comparison[[paste0(fw, "_replicate_direction")]] <- de$replicate_direction[idx]
}

write.csv(comparison, file.path(TABDIR, "Monocyte_genomewide_reference_specificity_v6.9.6.2.csv"), row.names=FALSE)

specificity_summary <- as.data.frame(table(specificity_class_ShTx=comparison$specificity_class_ShTx), stringsAsFactors=FALSE)
write.csv(specificity_summary, file.path(TABDIR, "Monocyte_genomewide_specificity_class_counts_v6.9.6.2.csv"), row.names=FALSE)

tested <- comparison[is.finite(comparison$PRIMARY_CORE_PValue), , drop=FALSE]
filter_summary <- data.frame(
  metric=c("PRIMARY_CORE_tested_genes", "REFERENCE_AWARE_keep_delta_le2", "STRICT_keep_delta_le1.5"),
  n_genes=c(
    nrow(tested),
    sum(tested$keep_reference_aware, na.rm=TRUE),
    sum(tested$keep_strict_reference_sensitivity, na.rm=TRUE)
  ),
  stringsAsFactors=FALSE
)
write.csv(filter_summary, file.path(TABDIR, "Monocyte_GSEA_filter_gene_counts_v6.9.6.2.csv"), row.names=FALSE)

key_genes <- c("Gdf15", "Thbs1", "Slc25a47", "Rgs1", "Cxcr4", "Ccn1")
key_table <- comparison[comparison$gene %in% key_genes, , drop=FALSE]
write.csv(key_table, file.path(TABDIR, "Monocyte_key_gene_genomewide_specificity_v6.9.6.2.csv"), row.names=FALSE)

top30 <- comparison[is.finite(comparison$PRIMARY_CORE_PValue), , drop=FALSE]
top30 <- top30[order(top30$PRIMARY_CORE_PValue), , drop=FALSE]
top30 <- head(top30, 30)
write.csv(top30, file.path(TABDIR, "Monocyte_PRIMARY_CORE_top30_genomewide_specificity_v6.9.6.2.csv"), row.names=FALSE)

cat("\n=== REFERENCE LINEAGE CELL COUNTS ===\n")
print(reference_counts, row.names=FALSE)

cat("\n=== GENOME-WIDE SPECIFICITY CLASS COUNTS ===\n")
print(specificity_summary, row.names=FALSE)

cat("\n=== GSEA FILTER GENE COUNTS ===\n")
print(filter_summary, row.names=FALSE)

cat("\n=== KEY GENE GENOME-WIDE SPECIFICITY ===\n")
key_cols <- c(
  "gene",
  "Monocyte_median_logCPM_ShTx",
  "max_reference_lineage_ShTx",
  "max_reference_delta_ShTx",
  "specificity_class_ShTx",
  "keep_reference_aware",
  "keep_strict_reference_sensitivity",
  "PRIMARY_CORE_logFC",
  "PRIMARY_CORE_FDR",
  "PRIMARY_CORE_replicate_direction"
)
print(key_table[, intersect(key_cols, colnames(key_table)), drop=FALSE], row.names=FALSE)

cat("\n=== TOP30 PRIMARY_CORE GENOME-WIDE SPECIFICITY ===\n")
print(top30[, intersect(key_cols, colnames(top30)), drop=FALSE], row.names=FALSE)

cat("\n====================================================\n")
cat("v6.9.6.2 COMPLETE\n")
cat("Genome-wide reference specificity audit complete\n")
cat("Primary specificity comparison: Sham1/Sham20/Tx17/Tx5\n")
cat("Reference-aware filter: max competing-lineage delta <= 2\n")
cat("Strict sensitivity filter: max competing-lineage delta <= 1.5\n")
cat("No genes removed from source data\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(OUTDIR, "sessionInfo_v6.9.6.2.txt")
)
