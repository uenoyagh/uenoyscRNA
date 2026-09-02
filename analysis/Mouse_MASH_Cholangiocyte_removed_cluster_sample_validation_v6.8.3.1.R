suppressPackageStartupMessages({
  library(Seurat)
})

VERSION <- "v6.8.3.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.2/objects/",
  "Mouse_MASH_Cholangiocyte_res0.3_audit_v6.8.2.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)

remove_clusters <- c("3", "4", "9", "10", "11", "12")

required_cols <- c(
  "sample",
  "condition",
  "Chol_res03_v682",
  "Biliary_core_hits_v682",
  "Hepatocyte_hits_v682",
  "Myeloid_hits_v682",
  "Vascular_endothelial_hits_v682",
  "LSEC_hits_v682",
  "HSC_mesenchymal_hits_v682",
  "Neutrophil_hits_v682",
  "Lymphoid_hits_v682"
)

missing_cols <- setdiff(
  required_cols,
  colnames(obj@meta.data)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_cols, collapse=", ")
  )
}

md <- obj@meta.data
md$cluster <- as.character(md$Chol_res03_v682)

md <- md[
  md$cluster %in% remove_clusters,
  ,
  drop=FALSE
]

cat("====================================================\n")
cat("Removed-cluster sample-specific validation\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n")

cat("\n=== REMOVED CELLS BY CLUSTER x SAMPLE ===\n")

count_tab <- table(
  cluster=md$cluster,
  sample=md$sample
)

print(count_tab)

write.csv(
  as.data.frame.matrix(count_tab),
  file.path(
    TABDIR,
    "Removed_clusters_by_sample_counts_v6.8.3.1.csv"
  )
)

hit_cols <- c(
  Biliary_core="Biliary_core_hits_v682",
  Hepatocyte="Hepatocyte_hits_v682",
  Myeloid="Myeloid_hits_v682",
  Vascular_endothelial="Vascular_endothelial_hits_v682",
  LSEC="LSEC_hits_v682",
  HSC_mesenchymal="HSC_mesenchymal_hits_v682",
  Neutrophil="Neutrophil_hits_v682",
  Lymphoid="Lymphoid_hits_v682"
)

groups <- split(
  seq_len(nrow(md)),
  interaction(
    md$cluster,
    md$sample,
    drop=TRUE
  )
)

summary_list <- lapply(
  groups,
  function(idx) {

    x <- md[idx, , drop=FALSE]

    out <- data.frame(
      cluster=x$cluster[1],
      sample=as.character(x$sample[1]),
      condition=as.character(x$condition[1]),
      n_cells=nrow(x),
      stringsAsFactors=FALSE
    )

    for (nm in names(hit_cols)) {

      col <- hit_cols[[nm]]

      out[[paste0(nm, "_median_hits")]] <-
        median(x[[col]])

      out[[paste0(nm, "_pct_ge2hits")]] <-
        mean(x[[col]] >= 2) * 100
    }

    out
  }
)

summary_df <- do.call(
  rbind,
  summary_list
)

rownames(summary_df) <- NULL

summary_df <- summary_df[
  order(
    as.numeric(summary_df$cluster),
    summary_df$sample
  ),
  ,
  drop=FALSE
]

write.csv(
  summary_df,
  file.path(
    TABDIR,
    "Removed_clusters_lineage_by_cluster_sample_v6.8.3.1.csv"
  ),
  row.names=FALSE
)

cat("\n=== CLUSTER x SAMPLE LINEAGE SUMMARY ===\n")

print(
  summary_df[
    ,
    c(
      "cluster",
      "sample",
      "n_cells",
      "Biliary_core_pct_ge2hits",
      "Myeloid_pct_ge2hits",
      "Vascular_endothelial_pct_ge2hits",
      "LSEC_pct_ge2hits",
      "HSC_mesenchymal_pct_ge2hits",
      "Neutrophil_pct_ge2hits",
      "Lymphoid_pct_ge2hits"
    )
  ],
  row.names=FALSE
)

# ---------------------------------------------------------
# Sample-level contribution of each removed cluster
# ---------------------------------------------------------

sample_total <- table(
  obj$sample
)

removed_sample <- table(
  factor(
    md$sample,
    levels=names(sample_total)
  )
)

removed_fraction <- data.frame(
  sample=names(sample_total),
  total_parent=as.integer(sample_total),
  removed=as.integer(removed_sample),
  removed_fraction=
    as.integer(removed_sample) /
    as.integer(sample_total),
  stringsAsFactors=FALSE
)

write.csv(
  removed_fraction,
  file.path(
    TABDIR,
    "Removed_fraction_by_sample_v6.8.3.1.csv"
  ),
  row.names=FALSE
)

cat("\n=== TOTAL REMOVED FRACTION BY SAMPLE ===\n")
print(removed_fraction)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.3.1.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.3.1 COMPLETE\n")
cat("No cells changed\n")
cat("No object rewritten\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
