suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

VERSION <- "v6.9.2"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.1/objects/",
  "Mouse_MASH_Monocyte_RPCA_resolution_scan_v6.9.1.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
OBJDIR <- file.path(OUTDIR, "objects")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte lineage cleanup\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

CLUSTER_COL <- "monocyte_res0_4"

if (!(CLUSTER_COL %in% colnames(obj@meta.data))) {
  stop("Missing cluster column: ", CLUSTER_COL)
}

# ---------------------------------------------------------
# Final lineage-cleanup decision
#
# Remove:
#   cluster 7 = T/NK-like lymphoid contamination
#   cluster 8 = B-cell contamination
#
# Keep:
#   cluster 4 = Mmp8/Chil3/Sell/Cd177 monocyte-like state
#   cluster 3 = Ms4a7/Lpl/Mmp12/Gpr84 transitional state
#   cluster 5 = Nos2/Cxcl9/Saa3 inflammatory state
#   cluster 6 = interferon-response state
# ---------------------------------------------------------

remove_clusters <- c("7", "8")

cluster_vec <- as.character(
  obj@meta.data[[CLUSTER_COL]]
)

remove_cells <- rownames(obj@meta.data)[
  cluster_vec %in% remove_clusters
]

keep_cells <- setdiff(
  colnames(obj),
  remove_cells
)

cat("Input cells:", ncol(obj), "\n")
cat("Remove cells:", length(remove_cells), "\n")
cat("Expected retained:", length(keep_cells), "\n\n")

cat("=== REMOVED CLUSTERS ===\n")
print(
  table(
    cluster=cluster_vec[
      match(
        remove_cells,
        rownames(obj@meta.data)
      )
    ]
  )
)

# ---------------------------------------------------------
# Resolve sample column
# ---------------------------------------------------------

sample_candidates <- c(
  "sample",
  "sample_id",
  "Sample",
  "orig.ident"
)

sample_col <- sample_candidates[
  sample_candidates %in%
    colnames(obj@meta.data)
][1]

if (is.na(sample_col)) {
  stop("Could not resolve sample column.")
}

sample_vec <- as.character(
  obj@meta.data[[sample_col]]
)

# ---------------------------------------------------------
# Removal by sample
# ---------------------------------------------------------

sample_names <- sort(
  unique(sample_vec)
)

removal_summary <- do.call(
  rbind,
  lapply(
    sample_names,
    function(s) {

      all_cells <- rownames(obj@meta.data)[
        sample_vec == s
      ]

      removed <- intersect(
        all_cells,
        remove_cells
      )

      retained <- setdiff(
        all_cells,
        remove_cells
      )

      data.frame(
        sample=s,
        n_before=length(all_cells),
        n_removed=length(removed),
        n_retained=length(retained),
        removal_fraction=
          length(removed) /
          length(all_cells),
        stringsAsFactors=FALSE
      )
    }
  )
)

write.csv(
  removal_summary,
  file.path(
    TABDIR,
    "Monocyte_lineage_cleanup_by_sample_v6.9.2.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Raw RNA counts
# ---------------------------------------------------------

DefaultAssay(obj) <- "RNA"

get_raw_counts <- function(obj) {

  assay <- obj[["RNA"]]

  if (inherits(assay, "Assay5")) {

    layer_names <- Layers(assay)

    count_layers <- grep(
      "^counts",
      layer_names,
      value=TRUE
    )

    if (length(count_layers) == 0) {
      stop(
        "No RNA counts layer found: ",
        paste(layer_names, collapse=", ")
      )
    }

    if (length(count_layers) > 1) {

      obj[["RNA"]] <- JoinLayers(
        obj[["RNA"]],
        layers=count_layers,
        new="counts"
      )
    }

    counts <- GetAssayData(
      obj,
      assay="RNA",
      layer="counts"
    )

  } else {

    counts <- GetAssayData(
      obj,
      assay="RNA",
      layer="counts"
    )
  }

  counts
}

counts <- get_raw_counts(obj)

counts_clean <- counts[
  ,
  keep_cells,
  drop=FALSE
]

meta_clean <- obj@meta.data[
  keep_cells,
  ,
  drop=FALSE
]

# ---------------------------------------------------------
# Build raw-only lineage-clean object
# ---------------------------------------------------------

clean <- CreateSeuratObject(
  counts=counts_clean,
  meta.data=meta_clean,
  project="Mouse_MASH_Monocyte_lineage_clean_v6.9.2"
)

if (ncol(clean) != length(keep_cells)) {
  stop("Retained cell-count mismatch.")
}

saveRDS(
  clean,
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_lineage_clean_raw_v6.9.2.rds"
  )
)

# ---------------------------------------------------------
# Terminal output
# ---------------------------------------------------------

cat("\n=== LINEAGE CLEANUP BY SAMPLE ===\n")
print(
  removal_summary,
  row.names=FALSE
)

cat("\n=== FINAL CELL COUNTS ===\n")
cat("Before:", ncol(obj), "\n")
cat("Removed:", length(remove_cells), "\n")
cat("Retained:", ncol(clean), "\n")
cat("Features:", nrow(clean), "\n")

cat("\n=== CLEAN OBJECT ===\n")
cat(
  "Assays:",
  paste(Assays(clean), collapse=", "),
  "\n"
)
cat(
  "Reductions:",
  paste(Reductions(clean), collapse=", "),
  "\n"
)
cat(
  "Graphs:",
  paste(names(clean@graphs), collapse=", "),
  "\n"
)

cat("\n====================================================\n")
cat("v6.9.2 COMPLETE\n")
cat("Monocyte lineage cleanup complete\n")
cat("Removed res0.4 clusters: 7, 8\n")
cat("Interpretation: lymphoid + B-cell contamination\n")
cat("No other Monocyte state removed\n")
cat(
  "Clean raw object:",
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_lineage_clean_raw_v6.9.2.rds"
  ),
  "\n"
)
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.2.txt"
  )
)
