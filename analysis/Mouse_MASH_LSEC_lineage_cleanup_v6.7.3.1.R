suppressPackageStartupMessages({
  library(Seurat)
})

VERSION <- "v6.7.3.1"
RES_COL <- "LSEC_res0.6"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_v6.7.2/objects/",
  "Mouse_MASH_endothelial_res0.6_marker_audit_v6.7.2.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH LSEC lineage cleanup\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

obj <- readRDS(INPUT_RDS)

if (!RES_COL %in% colnames(obj@meta.data)) {
  stop("Missing cluster column: ", RES_COL)
}

cl <- as.character(
  obj@meta.data[[RES_COL]]
)

# --------------------------------------------------
# Fixed lineage decisions from v6.7.2 audit
# --------------------------------------------------

LSEC_CLUSTERS <- c(
  "0","1","2","3","4"
)

VASCULAR_CLUSTERS <- c(
  "6","7","12"
)

LYMPHATIC_CLUSTERS <- c(
  "15"
)

REMOVE_MAP <- c(
  "5"  = "Hepatocyte_like",
  "8"  = "Myeloid",
  "9"  = "HSC_Mesenchymal",
  "10" = "Cholangiocyte",
  "11" = "T_NK",
  "13" = "Neutrophil",
  "14" = "B_cell",
  "16" = "DC_Myeloid"
)

expected_clusters <- as.character(0:16)

observed_clusters <- sort(
  unique(cl)
)

missing_expected <- setdiff(
  expected_clusters,
  observed_clusters
)

if (length(missing_expected) > 0) {
  stop(
    "Expected res0.6 cluster(s) absent: ",
    paste(missing_expected, collapse=", ")
  )
}

audit_class <- rep(
  NA_character_,
  ncol(obj)
)

audit_class[
  cl %in% LSEC_CLUSTERS
] <- "LSEC_candidate"

audit_class[
  cl %in% VASCULAR_CLUSTERS
] <- "Vascular_EC"

audit_class[
  cl %in% LYMPHATIC_CLUSTERS
] <- "Lymphatic_EC"

for (x in names(REMOVE_MAP)) {
  audit_class[
    cl == x
  ] <- paste0(
    "Remove_",
    REMOVE_MAP[[x]]
  )
}

if (any(is.na(audit_class))) {
  bad <- unique(
    cl[is.na(audit_class)]
  )

  stop(
    "Unclassified cluster(s): ",
    paste(bad, collapse=", ")
  )
}

obj$endothelial_lineage_audit_v673 <-
  audit_class

cat("=== AUDIT CLASS COUNTS ===\n")
print(
  sort(
    table(
      obj$endothelial_lineage_audit_v673
    ),
    decreasing=TRUE
  )
)

# --------------------------------------------------
# Cluster-level decision table
# --------------------------------------------------

decision <- data.frame(
  cluster=expected_clusters,
  decision=NA_character_,
  biological_call=NA_character_,
  stringsAsFactors=FALSE
)

decision$decision[
  decision$cluster %in% LSEC_CLUSTERS
] <- "KEEP_LSEC"

decision$biological_call[
  decision$cluster %in% LSEC_CLUSTERS
] <- "LSEC_candidate"

decision$decision[
  decision$cluster %in% VASCULAR_CLUSTERS
] <- "KEEP_ENDOTHELIAL_NON_LSEC"

decision$biological_call[
  decision$cluster %in% VASCULAR_CLUSTERS
] <- "Vascular_EC"

decision$decision[
  decision$cluster %in% LYMPHATIC_CLUSTERS
] <- "KEEP_ENDOTHELIAL_NON_LSEC"

decision$biological_call[
  decision$cluster %in% LYMPHATIC_CLUSTERS
] <- "Lymphatic_EC"

for (x in names(REMOVE_MAP)) {

  decision$decision[
    decision$cluster == x
  ] <- "REMOVE_NON_ENDOTHELIAL"

  decision$biological_call[
    decision$cluster == x
  ] <- REMOVE_MAP[[x]]
}

cluster_n <- table(cl)

decision$n_cells <- as.integer(
  cluster_n[
    decision$cluster
  ]
)

write.csv(
  decision,
  file.path(
    TABDIR,
    "Endothelial_cluster_cleanup_decisions_v6.7.3.1.csv"
  ),
  row.names=FALSE
)

# --------------------------------------------------
# Counts before cleanup
# --------------------------------------------------

write.csv(
  as.data.frame(
    table(
      cluster=cl,
      sample=obj$sample
    )
  ),
  file.path(
    TABDIR,
    "Precleanup_cluster_by_sample_v6.7.3.1.csv"
  ),
  row.names=FALSE
)

write.csv(
  as.data.frame(
    table(
      cluster=cl,
      condition=obj$condition
    )
  ),
  file.path(
    TABDIR,
    "Precleanup_cluster_by_condition_v6.7.3.1.csv"
  ),
  row.names=FALSE
)

# --------------------------------------------------
# Clean endothelial object
# Keeps LSEC + vascular + lymphatic EC
# --------------------------------------------------

ENDO_KEEP_CLUSTERS <- c(
  LSEC_CLUSTERS,
  VASCULAR_CLUSTERS,
  LYMPHATIC_CLUSTERS
)

endo_cells <- colnames(obj)[
  cl %in% ENDO_KEEP_CLUSTERS
]

endo_clean <- subset(
  obj,
  cells=endo_cells
)

cat("\n=== CLEAN ENDOTHELIAL ===\n")
cat("Cells:", ncol(endo_clean), "\n")

print(
  table(
    endo_clean$endothelial_lineage_audit_v673
  )
)

if (ncol(endo_clean) != 11203) {
  stop(
    "Unexpected clean endothelial count: ",
    ncol(endo_clean),
    " ; expected 11203"
  )
}

saveRDS(
  endo_clean,
  file.path(
    OBJDIR,
    "Mouse_MASH_clean_endothelial_parent_v6.7.3.1.rds"
  ),
  compress=FALSE
)

write.csv(
  as.data.frame(
    table(
      sample=endo_clean$sample,
      class=endo_clean$endothelial_lineage_audit_v673
    )
  ),
  file.path(
    TABDIR,
    "Clean_endothelial_by_sample_v6.7.3.1.csv"
  ),
  row.names=FALSE
)

write.csv(
  as.data.frame(
    table(
      condition=endo_clean$condition,
      class=endo_clean$endothelial_lineage_audit_v673
    )
  ),
  file.path(
    TABDIR,
    "Clean_endothelial_by_condition_v6.7.3.1.csv"
  ),
  row.names=FALSE
)

# --------------------------------------------------
# Extract LSEC candidate cells
# --------------------------------------------------

lsec_cells <- colnames(obj)[
  cl %in% LSEC_CLUSTERS
]

lsec_parent <- subset(
  obj,
  cells=lsec_cells
)

cat("\n=== LSEC CANDIDATE FROM RES0.6 ===\n")
cat("Cells:", ncol(lsec_parent), "\n")

print(
  table(
    lsec_parent@meta.data[[RES_COL]]
  )
)

if (ncol(lsec_parent) != 9190) {
  stop(
    "Unexpected LSEC candidate count: ",
    ncol(lsec_parent),
    " ; expected 9190"
  )
}

cat("\nLSEC by sample:\n")
print(
  table(lsec_parent$sample)
)

cat("\nLSEC by condition:\n")
print(
  table(lsec_parent$condition)
)

write.csv(
  as.data.frame(
    table(
      sample=lsec_parent$sample,
      cluster=lsec_parent@meta.data[[RES_COL]]
    )
  ),
  file.path(
    TABDIR,
    "LSEC_candidate_by_sample_v6.7.3.1.csv"
  ),
  row.names=FALSE
)

write.csv(
  as.data.frame(
    table(
      condition=lsec_parent$condition,
      cluster=lsec_parent@meta.data[[RES_COL]]
    )
  ),
  file.path(
    TABDIR,
    "LSEC_candidate_by_condition_v6.7.3.1.csv"
  ),
  row.names=FALSE
)

# --------------------------------------------------
# Extract RAW RNA counts
# --------------------------------------------------

DefaultAssay(lsec_parent) <- "RNA"

cat("\n=== EXTRACT RAW COUNTS ===\n")

if (inherits(
  lsec_parent[["RNA"]],
  "Assay5"
)) {

  layers <- Layers(
    lsec_parent[["RNA"]]
  )

  cat(
    "RNA layers:",
    paste(layers, collapse=", "),
    "\n"
  )

  count_layers <- grep(
    "^counts($|\\.)",
    layers,
    value=TRUE
  )

  if (length(count_layers) == 0) {
    stop("No RNA counts layer found.")
  }

  if (
    length(count_layers) == 1 &&
    identical(count_layers, "counts")
  ) {

    counts <- LayerData(
      lsec_parent,
      assay="RNA",
      layer="counts"
    )

  } else {

    lsec_parent <- JoinLayers(
      lsec_parent,
      assay="RNA"
    )

    counts <- LayerData(
      lsec_parent,
      assay="RNA",
      layer="counts"
    )
  }

} else {

  counts <- GetAssayData(
    lsec_parent,
    assay="RNA",
    layer="counts"
  )
}

metadata <- lsec_parent@meta.data[
  colnames(counts),
  ,
  drop=FALSE
]

cat(
  "Counts:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "cells\n"
)

if (ncol(counts) != 9190) {
  stop(
    "Raw count extraction did not retain 9190 cells."
  )
}

# --------------------------------------------------
# Completely new LSEC object
# --------------------------------------------------

options(
  Seurat.object.assay.version="v3"
)

lsec_clean <- CreateSeuratObject(
  counts=counts,
  meta.data=metadata,
  project="Mouse_MASH_LSEC_clean"
)

rm(
  lsec_parent,
  counts
)

invisible(gc())

cat("\n=== CLEAN LSEC OBJECT ===\n")
cat("Cells:", ncol(lsec_clean), "\n")
cat("Features:", nrow(lsec_clean), "\n")
cat(
  "Assays:",
  paste(Assays(lsec_clean), collapse=", "),
  "\n"
)
cat(
  "RNA assay class:",
  class(lsec_clean[["RNA"]])[1],
  "\n"
)
cat(
  "Reductions:",
  length(Reductions(lsec_clean)),
  "\n"
)
cat(
  "Graphs:",
  length(lsec_clean@graphs),
  "\n"
)

if (!identical(Assays(lsec_clean), "RNA")) {
  stop(
    "Clean LSEC object contains unexpected assay."
  )
}

if (!inherits(
  lsec_clean[["RNA"]],
  "Assay"
)) {
  stop(
    "Clean LSEC RNA assay is not legacy Assay."
  )
}

if (length(Reductions(lsec_clean)) != 0) {
  stop(
    "Clean LSEC object unexpectedly has reductions."
  )
}

if (length(lsec_clean@graphs) != 0) {
  stop(
    "Clean LSEC object unexpectedly has graphs."
  )
}

if (ncol(lsec_clean) != 9190) {
  stop(
    "Clean LSEC object cell count != 9190."
  )
}

lsec_clean$LSEC_cleanup_status_v673 <-
  "LSEC_candidate_clean"

saveRDS(
  lsec_clean,
  file.path(
    OBJDIR,
    "Mouse_MASH_LSEC_raw_clean_v6.7.3.1.rds"
  ),
  compress=FALSE
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.7.3.1.txt"
  )
)

cat("\n====================================================\n")
cat("v6.7.3.1 COMPLETE\n")
cat("Original parent: 14402\n")
cat("Clean endothelial: ", ncol(endo_clean), "\n", sep="")
cat("Clean LSEC candidate: ", ncol(lsec_clean), "\n", sep="")
cat("Non-endothelial removed: 3199\n")
cat("No RPCA reclustering performed yet\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
