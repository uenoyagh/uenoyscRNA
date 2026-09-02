suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

set.seed(20260902)

VERSION <- "v6.8.3.2"

AUDIT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.2/objects/",
  "Mouse_MASH_Cholangiocyte_res0.3_audit_v6.8.2.rds"
)

RAW_PARENT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.0/objects/",
  "Mouse_MASH_Cholangiocyte_parent_raw_clean_v6.8.0.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte REVISED lineage cleanup\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(AUDIT_RDS)) {
  stop("Missing audit RDS: ", AUDIT_RDS)
}

if (!file.exists(RAW_PARENT_RDS)) {
  stop("Missing raw parent RDS: ", RAW_PARENT_RDS)
}

audit_obj <- readRDS(AUDIT_RDS)
raw_parent <- readRDS(RAW_PARENT_RDS)

if (!"Chol_res03_v682" %in% colnames(audit_obj@meta.data)) {
  stop("Missing Chol_res03_v682.")
}

if (!setequal(colnames(audit_obj), colnames(raw_parent))) {
  stop("Audit and raw-parent cell sets do not match.")
}

cat("=== INPUT ===\n")
cat("Audit cells:", ncol(audit_obj), "\n")
cat("Raw parent cells:", ncol(raw_parent), "\n")

# ---------------------------------------------------------
# Revised biological decision
#
# cluster 4 is restored because:
# - Biliary core >=2 hits in 100% of cells
# - strong STD/CDHFD enrichment
# - removing it caused disease-axis-specific cell loss
#
# It is retained as KEEP_QC_WATCH.
# ---------------------------------------------------------

decision <- data.frame(
  cluster=as.character(0:13),

  action=c(
    "KEEP",
    "KEEP",
    "KEEP",
    "REMOVE",
    "KEEP_QC_WATCH",
    "KEEP",
    "KEEP",
    "KEEP",
    "KEEP",
    "REMOVE",
    "REMOVE",
    "REMOVE",
    "REMOVE",
    "KEEP"
  ),

  interpretation=c(
    "Biliary_main",
    "Biliary_main",
    "Inflammatory_reactive_biliary_candidate",
    "Myeloid_contamination",
    "Disease_enriched_mixed_reactive_biliary_QC_watch",
    "Cycling_biliary_candidate",
    "Krt20_Cdh17_reactive_epithelial_candidate",
    "Ciliated_epithelial_candidate",
    "Dmbt1_Duox2_reactive_epithelial_candidate",
    "Vascular_mesenchymal_contamination",
    "Neutrophil_contamination",
    "Dendritic_cell_contamination",
    "B_cell_contamination",
    "Rare_Pou2f3_Trpm5_tuft_like_candidate"
  ),

  stringsAsFactors=FALSE
)

write.csv(
  decision,
  file.path(
    TABDIR,
    "Cholangiocyte_res0.3_revised_cleanup_decisions_v6.8.3.2.csv"
  ),
  row.names=FALSE
)

remove_clusters <- decision$cluster[
  decision$action == "REMOVE"
]

cat("\n=== REVISED CLEANUP DECISION ===\n")
cat(
  "REMOVE clusters:",
  paste(remove_clusters, collapse=", "),
  "\n"
)
cat("Cluster 4: KEEP_QC_WATCH\n")

cluster_vec <- as.character(
  audit_obj$Chol_res03_v682
)

keep_cells <- colnames(audit_obj)[
  !cluster_vec %in% remove_clusters
]

remove_cells <- colnames(audit_obj)[
  cluster_vec %in% remove_clusters
]

cat("\n=== CELL COUNTS ===\n")
cat("Before:", ncol(audit_obj), "\n")
cat("Retained:", length(keep_cells), "\n")
cat("Removed:", length(remove_cells), "\n")

if (length(keep_cells) != 15755) {
  warning(
    "Expected 15755 retained cells, observed ",
    length(keep_cells)
  )
}

# ---------------------------------------------------------
# Removed cells by cluster
# ---------------------------------------------------------

cat("\n=== REMOVED CELLS BY CLUSTER ===\n")

removed_cluster_table <- table(
  factor(
    cluster_vec[
      cluster_vec %in% remove_clusters
    ],
    levels=remove_clusters
  )
)

print(removed_cluster_table)

write.csv(
  data.frame(
    cluster=names(removed_cluster_table),
    removed_cells=as.integer(removed_cluster_table)
  ),
  file.path(
    TABDIR,
    "Cholangiocyte_removed_cells_by_cluster_v6.8.3.2.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Sample-level revised cleanup audit
# ---------------------------------------------------------

md_raw <- raw_parent@meta.data

sample_levels <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

sample_before <- table(
  factor(
    as.character(md_raw$sample),
    levels=sample_levels
  )
)

sample_retained <- table(
  factor(
    as.character(
      md_raw[
        keep_cells,
        "sample"
      ]
    ),
    levels=sample_levels
  )
)

sample_removed <-
  sample_before - sample_retained

sample_summary <- data.frame(
  sample=sample_levels,
  before=as.integer(sample_before),
  retained=as.integer(sample_retained),
  removed=as.integer(sample_removed),
  removed_fraction=
    as.integer(sample_removed) /
    as.integer(sample_before),
  stringsAsFactors=FALSE
)

write.csv(
  sample_summary,
  file.path(
    TABDIR,
    "Cholangiocyte_revised_cleanup_by_sample_v6.8.3.2.csv"
  ),
  row.names=FALSE
)

cat("\n=== REVISED CLEANUP BY SAMPLE ===\n")
print(sample_summary)

# ---------------------------------------------------------
# Condition-level revised cleanup audit
# ---------------------------------------------------------

condition_levels <- c(
  "STD",
  "CDHFD",
  "Sham",
  "Tx"
)

condition_before <- table(
  factor(
    as.character(md_raw$condition),
    levels=condition_levels
  )
)

condition_retained <- table(
  factor(
    as.character(
      md_raw[
        keep_cells,
        "condition"
      ]
    ),
    levels=condition_levels
  )
)

condition_removed <-
  condition_before - condition_retained

condition_summary <- data.frame(
  condition=condition_levels,
  before=as.integer(condition_before),
  retained=as.integer(condition_retained),
  removed=as.integer(condition_removed),
  removed_fraction=
    as.integer(condition_removed) /
    as.integer(condition_before),
  stringsAsFactors=FALSE
)

write.csv(
  condition_summary,
  file.path(
    TABDIR,
    "Cholangiocyte_revised_cleanup_by_condition_v6.8.3.2.csv"
  ),
  row.names=FALSE
)

cat("\n=== REVISED CLEANUP BY CONDITION ===\n")
print(condition_summary)

# ---------------------------------------------------------
# Cluster 4 QC-watch audit
# ---------------------------------------------------------

cluster4_cells <- colnames(audit_obj)[
  cluster_vec == "4"
]

cluster4_md <- audit_obj@meta.data[
  cluster4_cells,
  ,
  drop=FALSE
]

lineage_hit_cols <- c(
  "Biliary_core_hits_v682",
  "Hepatocyte_hits_v682",
  "Myeloid_hits_v682",
  "Vascular_endothelial_hits_v682",
  "LSEC_hits_v682",
  "HSC_mesenchymal_hits_v682",
  "Neutrophil_hits_v682",
  "Lymphoid_hits_v682"
)

cluster4_groups <- split(
  seq_len(nrow(cluster4_md)),
  as.character(cluster4_md$sample)
)

cluster4_summary <- do.call(
  rbind,
  lapply(
    cluster4_groups,
    function(idx) {

      x <- cluster4_md[
        idx,
        ,
        drop=FALSE
      ]

      out <- data.frame(
        sample=as.character(x$sample[1]),
        n_cells=nrow(x),
        median_nCount_RNA=
          median(x$nCount_RNA),
        median_nFeature_RNA=
          median(x$nFeature_RNA),
        stringsAsFactors=FALSE
      )

      for (col in lineage_hit_cols) {

        short <- sub(
          "_hits_v682$",
          "",
          col
        )

        out[[paste0(short, "_pct_ge2hits")]] <- mean(
          x[[col]] >= 2
        ) * 100
      }

      out
    }
  )
)

rownames(cluster4_summary) <- NULL

write.csv(
  cluster4_summary,
  file.path(
    TABDIR,
    "Cholangiocyte_cluster4_QC_watch_by_sample_v6.8.3.2.csv"
  ),
  row.names=FALSE
)

cat("\n=== CLUSTER 4 QC WATCH ===\n")
print(cluster4_summary)

# ---------------------------------------------------------
# Rebuild from RAW counts
# ---------------------------------------------------------

DefaultAssay(raw_parent) <- "RNA"

counts <- GetAssayData(
  raw_parent,
  assay="RNA",
  layer="counts"
)

counts_keep <- counts[
  ,
  keep_cells,
  drop=FALSE
]

meta_keep <- md_raw[
  keep_cells,
  ,
  drop=FALSE
]

meta_keep$Chol_res03_precleanup_v6832 <-
  cluster_vec[
    match(
      keep_cells,
      colnames(audit_obj)
    )
  ]

decision_map <- setNames(
  decision$action,
  decision$cluster
)

meta_keep$cleanup_status_v6832 <-
  unname(
    decision_map[
      meta_keep$Chol_res03_precleanup_v6832
    ]
  )

clean <- CreateSeuratObject(
  counts=counts_keep,
  assay="RNA",
  meta.data=meta_keep,
  project="Mouse_MASH_Cholangiocyte_clean_v6.8.3.2"
)

DefaultAssay(clean) <- "RNA"

cat("\n=== REVISED CLEAN OBJECT ===\n")
cat("Cells:", ncol(clean), "\n")
cat("Features:", nrow(clean), "\n")
cat("Assays:", paste(Assays(clean), collapse=", "), "\n")
cat("RNA assay class:", class(clean[["RNA"]])[1], "\n")
cat("Reductions:", length(Reductions(clean)), "\n")
cat("Graphs:", length(clean@graphs), "\n")

# ---------------------------------------------------------
# Hard validation
# ---------------------------------------------------------

stopifnot(
  ncol(clean) == 15755
)

stopifnot(
  identical(
    sort(colnames(clean)),
    sort(keep_cells)
  )
)

stopifnot(
  length(Reductions(clean)) == 0
)

stopifnot(
  length(clean@graphs) == 0
)

stopifnot(
  identical(
    Assays(clean),
    "RNA"
  )
)

cat("\nREVISED CLEAN CHOLANGIOCYTE CHECK: PASSED\n")

# ---------------------------------------------------------
# Save revised clean object
# ---------------------------------------------------------

saveRDS(
  clean,
  file.path(
    OBJDIR,
    "Mouse_MASH_Cholangiocyte_lineage_clean_raw_v6.8.3.2.rds"
  ),
  compress=FALSE
)

removed_table <- data.frame(
  cell=remove_cells,
  cluster=
    cluster_vec[
      match(
        remove_cells,
        colnames(audit_obj)
      )
    ],
  sample=
    as.character(
      audit_obj@meta.data[
        remove_cells,
        "sample"
      ]
    ),
  condition=
    as.character(
      audit_obj@meta.data[
        remove_cells,
        "condition"
      ]
    ),
  stringsAsFactors=FALSE
)

write.csv(
  removed_table,
  file.path(
    TABDIR,
    "Cholangiocyte_removed_cells_v6.8.3.2.csv"
  ),
  row.names=FALSE
)

summary_lines <- c(
  "# Mouse MASH Cholangiocyte revised lineage cleanup v6.8.3.2",
  "",
  "- v6.8.3 is superseded but retained for provenance.",
  paste0("- Parent cells: ", ncol(raw_parent)),
  paste0("- Retained cells: ", ncol(clean)),
  paste0("- Removed cells: ", length(remove_cells)),
  "- Removed res0.3 clusters: 3, 9, 10, 11, 12",
  "- Cluster 4 was restored and retained as KEEP_QC_WATCH.",
  "",
  "## Rationale for restoring cluster 4",
  "- 100% of cluster 4 cells had >=2 biliary-core markers in every sample.",
  "- Cluster 4 was strongly enriched in STD/CDHFD.",
  "- Removing cluster 4 caused disproportionate depletion of STD/CDHFD Cholangiocytes.",
  "- Retaining cluster 4 makes sample-level removal fractions substantially more balanced.",
  "- Cluster 4 remains a QC-watch population and will be re-evaluated after clean RPCA.",
  "",
  "## Confirmed removals",
  "- cluster 3: myeloid contamination",
  "- cluster 9: vascular/mesenchymal contamination",
  "- cluster 10: neutrophil contamination",
  "- cluster 11: dendritic-cell contamination",
  "- cluster 12: B-cell contamination",
  "",
  "- No final Cholangiocyte state annotation assigned.",
  "- Clean object rebuilt from raw RNA counts."
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_revised_lineage_cleanup_summary_v6.8.3.2.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.3.2.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.3.2 COMPLETE\n")
cat("Parent cells:", ncol(raw_parent), "\n")
cat("Retained cells:", ncol(clean), "\n")
cat("Removed cells:", length(remove_cells), "\n")
cat("Cluster 4 retained as KEEP_QC_WATCH\n")
cat("No final state annotation assigned\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
