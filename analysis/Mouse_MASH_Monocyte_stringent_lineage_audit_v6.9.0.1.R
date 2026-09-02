suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

VERSION <- "v6.9.0.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_v6.9.0/objects/",
  "Mouse_MASH_Monocyte_parent_raw_clean_v6.9.0.rds"
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
cat("Mouse MASH Monocyte stringent lineage audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)
DefaultAssay(obj) <- "RNA"

if (inherits(obj[["RNA"]], "Assay5")) {
  counts <- GetAssayData(
    obj,
    assay="RNA",
    layer="counts"
  )
} else {
  counts <- GetAssayData(
    obj,
    assay="RNA",
    slot="counts"
  )
}

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
  sample_candidates %in% colnames(obj@meta.data)
][1]

if (is.na(sample_col)) {
  stop("Could not resolve sample column.")
}

sample_vec <- as.character(
  obj@meta.data[[sample_col]]
)

# ---------------------------------------------------------
# Marker definitions
#
# Deliberately exclude weak/non-specific markers:
#   Adgre1 from strict macrophage call
#   Itgax from strict DC call
#   Klf2 from strict endothelial call
#   S100a8/S100a9 from strict neutrophil call
# ---------------------------------------------------------

marker_sets <- list(

  Monocyte_identity=c(
    "Lyz2",
    "Ly6c2",
    "Ccr2",
    "Fcgr3",
    "Lst1",
    "Tyrobp",
    "Ctss"
  ),

  Neutrophil_strict=c(
    "Cxcr2",
    "Retnlg",
    "Camp",
    "Ngp",
    "Mmp8"
  ),

  Resident_Macrophage_strict=c(
    "Clec4f",
    "Timd4",
    "Vsig4",
    "Marco",
    "Cd5l"
  ),

  DC_strict=c(
    "Flt3",
    "Xcr1",
    "Clec10a",
    "Cd209a"
  ),

  Endothelial_strict=c(
    "Pecam1",
    "Kdr",
    "Emcn"
  ),

  Mesenchymal_strict=c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Des"
  ),

  Epithelial_strict=c(
    "Krt8",
    "Krt18",
    "Krt19",
    "Epcam"
  ),

  Lymphoid_strict=c(
    "Cd3d",
    "Cd3e",
    "Cd79a",
    "Ms4a1",
    "Nkg7"
  )
)

all_markers <- unique(
  unlist(marker_sets)
)

marker_presence <- data.frame(
  gene=all_markers,
  present=all_markers %in% rownames(counts),
  stringsAsFactors=FALSE
)

write.csv(
  marker_presence,
  file.path(
    TABDIR,
    "Monocyte_audit_marker_presence_v6.9.0.1.csv"
  ),
  row.names=FALSE
)

cat("=== MARKER PRESENCE ===\n")
print(
  marker_presence,
  row.names=FALSE
)

# ---------------------------------------------------------
# Count marker hits per cell
# ---------------------------------------------------------

get_hits <- function(genes) {

  genes <- intersect(
    genes,
    rownames(counts)
  )

  if (length(genes) == 0) {
    return(
      rep(
        0L,
        ncol(counts)
      )
    )
  }

  as.integer(
    Matrix::colSums(
      counts[
        genes,
        ,
        drop=FALSE
      ] > 0
    )
  )
}

for (nm in names(marker_sets)) {

  hits <- get_hits(
    marker_sets[[nm]]
  )

  obj[[paste0(
    nm,
    "_hits_v6.9.0.1"
  )]] <- hits
}

# ---------------------------------------------------------
# Stringent competing-lineage flags
# ---------------------------------------------------------

meta <- obj@meta.data

# Neutrophil:
# >=2 of Cxcr2/Retnlg/Camp/Ngp/Mmp8.
meta$Neutrophil_strict_flag_v6.9.0.1 <-
  meta$Neutrophil_strict_hits_v6.9.0.1 >= 2

# Resident macrophage:
# >=2 resident Kupffer markers.
meta$Resident_Macrophage_strict_flag_v6.9.0.1 <-
  meta$Resident_Macrophage_strict_hits_v6.9.0.1 >= 2

# DC:
# >=2 lineage-supportive DC markers.
meta$DC_strict_flag_v6.9.0.1 <-
  meta$DC_strict_hits_v6.9.0.1 >= 2

# Endothelial:
# >=2 core endothelial markers.
meta$Endothelial_strict_flag_v6.9.0.1 <-
  meta$Endothelial_strict_hits_v6.9.0.1 >= 2

# Mesenchymal / epithelial / lymphoid:
meta$Mesenchymal_strict_flag_v6.9.0.1 <-
  meta$Mesenchymal_strict_hits_v6.9.0.1 >= 2

meta$Epithelial_strict_flag_v6.9.0.1 <-
  meta$Epithelial_strict_hits_v6.9.0.1 >= 2

meta$Lymphoid_strict_flag_v6.9.0.1 <-
  meta$Lymphoid_strict_hits_v6.9.0.1 >= 2

flag_cols <- c(
  "Neutrophil_strict_flag_v6.9.0.1",
  "Resident_Macrophage_strict_flag_v6.9.0.1",
  "DC_strict_flag_v6.9.0.1",
  "Endothelial_strict_flag_v6.9.0.1",
  "Mesenchymal_strict_flag_v6.9.0.1",
  "Epithelial_strict_flag_v6.9.0.1",
  "Lymphoid_strict_flag_v6.9.0.1"
)

flag_matrix <- as.matrix(
  meta[
    ,
    flag_cols,
    drop=FALSE
  ]
)

meta$competing_lineage_flag_count_v6.9.0.1 <-
  rowSums(flag_matrix)

meta$any_competing_lineage_flag_v6.9.0.1 <-
  meta$competing_lineage_flag_count_v6.9.0.1 > 0

# ---------------------------------------------------------
# Sample-relative high-complexity audit
#
# Q3 + 3 IQR within each sample.
# Diagnostic only.
# ---------------------------------------------------------

meta$high_nCount_sample_relative_v6.9.0.1 <- FALSE
meta$high_nFeature_sample_relative_v6.9.0.1 <- FALSE

for (s in sort(unique(sample_vec))) {

  idx <- which(
    sample_vec == s
  )

  nc <- meta$nCount_RNA[idx]
  nf <- meta$nFeature_RNA[idx]

  nc_q3 <- as.numeric(
    quantile(
      nc,
      0.75,
      na.rm=TRUE
    )
  )

  nc_iqr <- IQR(
    nc,
    na.rm=TRUE
  )

  nf_q3 <- as.numeric(
    quantile(
      nf,
      0.75,
      na.rm=TRUE
    )
  )

  nf_iqr <- IQR(
    nf,
    na.rm=TRUE
  )

  nc_cut <- nc_q3 + 3 * nc_iqr
  nf_cut <- nf_q3 + 3 * nf_iqr

  meta$high_nCount_sample_relative_v6.9.0.1[idx] <-
    nc > nc_cut

  meta$high_nFeature_sample_relative_v6.9.0.1[idx] <-
    nf > nf_cut
}

meta$high_complexity_both_v6.9.0.1 <-
  meta$high_nCount_sample_relative_v6.9.0.1 &
  meta$high_nFeature_sample_relative_v6.9.0.1

obj@meta.data <- meta

# ---------------------------------------------------------
# Sample-level summary
# ---------------------------------------------------------

summary_rows <- list()

for (s in sort(unique(sample_vec))) {

  cells <- rownames(meta)[
    sample_vec == s
  ]

  md <- meta[
    cells,
    ,
    drop=FALSE
  ]

  summary_rows[[s]] <- data.frame(
    sample=s,
    n_cells=nrow(md),

    Monocyte_identity_median_hits=
      median(
        md$Monocyte_identity_hits_v6.9.0.1
      ),

    Monocyte_identity_fraction_ge2=
      mean(
        md$Monocyte_identity_hits_v6.9.0.1 >= 2
      ),

    Neutrophil_strict_n=
      sum(
        md$Neutrophil_strict_flag_v6.9.0.1
      ),

    Neutrophil_strict_fraction=
      mean(
        md$Neutrophil_strict_flag_v6.9.0.1
      ),

    Resident_Macrophage_strict_n=
      sum(
        md$Resident_Macrophage_strict_flag_v6.9.0.1
      ),

    Resident_Macrophage_strict_fraction=
      mean(
        md$Resident_Macrophage_strict_flag_v6.9.0.1
      ),

    DC_strict_n=
      sum(
        md$DC_strict_flag_v6.9.0.1
      ),

    DC_strict_fraction=
      mean(
        md$DC_strict_flag_v6.9.0.1
      ),

    Endothelial_strict_n=
      sum(
        md$Endothelial_strict_flag_v6.9.0.1
      ),

    Endothelial_strict_fraction=
      mean(
        md$Endothelial_strict_flag_v6.9.0.1
      ),

    Mesenchymal_strict_n=
      sum(
        md$Mesenchymal_strict_flag_v6.9.0.1
      ),

    Mesenchymal_strict_fraction=
      mean(
        md$Mesenchymal_strict_flag_v6.9.0.1
      ),

    Epithelial_strict_n=
      sum(
        md$Epithelial_strict_flag_v6.9.0.1
      ),

    Epithelial_strict_fraction=
      mean(
        md$Epithelial_strict_flag_v6.9.0.1
      ),

    Lymphoid_strict_n=
      sum(
        md$Lymphoid_strict_flag_v6.9.0.1
      ),

    Lymphoid_strict_fraction=
      mean(
        md$Lymphoid_strict_flag_v6.9.0.1
      ),

    Any_competing_lineage_n=
      sum(
        md$any_competing_lineage_flag_v6.9.0.1
      ),

    Any_competing_lineage_fraction=
      mean(
        md$any_competing_lineage_flag_v6.9.0.1
      ),

    High_complexity_both_n=
      sum(
        md$high_complexity_both_v6.9.0.1
      ),

    High_complexity_both_fraction=
      mean(
        md$high_complexity_both_v6.9.0.1
      ),

    Competing_and_high_complexity_n=
      sum(
        md$any_competing_lineage_flag_v6.9.0.1 &
        md$high_complexity_both_v6.9.0.1
      ),

    stringsAsFactors=FALSE
  )
}

summary_table <- do.call(
  rbind,
  summary_rows
)

rownames(summary_table) <- NULL

write.csv(
  summary_table,
  file.path(
    TABDIR,
    "Monocyte_stringent_lineage_summary_by_sample_v6.9.0.1.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Cross-tab: competing lineage flags
# ---------------------------------------------------------

flag_counts <- data.frame(
  lineage=c(
    "Neutrophil",
    "Resident_Macrophage",
    "DC",
    "Endothelial",
    "Mesenchymal",
    "Epithelial",
    "Lymphoid"
  ),

  n_flagged=c(
    sum(meta$Neutrophil_strict_flag_v6.9.0.1),
    sum(meta$Resident_Macrophage_strict_flag_v6.9.0.1),
    sum(meta$DC_strict_flag_v6.9.0.1),
    sum(meta$Endothelial_strict_flag_v6.9.0.1),
    sum(meta$Mesenchymal_strict_flag_v6.9.0.1),
    sum(meta$Epithelial_strict_flag_v6.9.0.1),
    sum(meta$Lymphoid_strict_flag_v6.9.0.1)
  ),

  fraction=c(
    mean(meta$Neutrophil_strict_flag_v6.9.0.1),
    mean(meta$Resident_Macrophage_strict_flag_v6.9.0.1),
    mean(meta$DC_strict_flag_v6.9.0.1),
    mean(meta$Endothelial_strict_flag_v6.9.0.1),
    mean(meta$Mesenchymal_strict_flag_v6.9.0.1),
    mean(meta$Epithelial_strict_flag_v6.9.0.1),
    mean(meta$Lymphoid_strict_flag_v6.9.0.1)
  ),

  stringsAsFactors=FALSE
)

write.csv(
  flag_counts,
  file.path(
    TABDIR,
    "Monocyte_stringent_lineage_flag_counts_v6.9.0.1.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Save diagnostic object
# ---------------------------------------------------------

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_stringent_lineage_audit_v6.9.0.1.rds"
  )
)

# ---------------------------------------------------------
# Terminal report
# ---------------------------------------------------------

cat("\n=== STRINGENT LINEAGE FLAGS OVERALL ===\n")
print(
  flag_counts,
  row.names=FALSE
)

cat("\n=== STRINGENT LINEAGE SUMMARY BY SAMPLE ===\n")
print(
  summary_table,
  row.names=FALSE
)

cat("\n=== OVERLAP WITH HIGH COMPLEXITY ===\n")

overall_overlap <- data.frame(
  metric=c(
    "Any_competing_lineage",
    "High_complexity_both",
    "Competing_and_high_complexity"
  ),
  n=c(
    sum(
      meta$any_competing_lineage_flag_v6.9.0.1
    ),
    sum(
      meta$high_complexity_both_v6.9.0.1
    ),
    sum(
      meta$any_competing_lineage_flag_v6.9.0.1 &
      meta$high_complexity_both_v6.9.0.1
    )
  ),
  fraction=c(
    mean(
      meta$any_competing_lineage_flag_v6.9.0.1
    ),
    mean(
      meta$high_complexity_both_v6.9.0.1
    ),
    mean(
      meta$any_competing_lineage_flag_v6.9.0.1 &
      meta$high_complexity_both_v6.9.0.1
    )
  ),
  stringsAsFactors=FALSE
)

print(
  overall_overlap,
  row.names=FALSE
)

cat("\n====================================================\n")
cat("v6.9.0.1 COMPLETE\n")
cat("Stringent lineage contamination audit complete\n")
cat("No cells removed\n")
cat("No annotation changed\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.0.1.txt"
  )
)
