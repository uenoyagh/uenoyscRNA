suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

VERSION <- "v6.9.0"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/RDS3_annotation_visualization_v4.1.1/objects/",
  "RDS3_with_visualization_metadata_v4.1.1.rds"
)

ANNOTATION_COL <- "celltype_for_R8plot_FIXED2"
TARGET_LABEL <- "Monocyte"

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Monocyte_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
OBJDIR <- file.path(OUTDIR, "objects")

dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH Monocyte parent audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Input RDS not found: ", INPUT_RDS)
}

whole <- readRDS(INPUT_RDS)

cat("Whole object:\n")
cat("Cells:", ncol(whole), "\n")
cat("Features:", nrow(whole), "\n\n")

if (!(ANNOTATION_COL %in% colnames(whole@meta.data))) {
  stop("Missing annotation column: ", ANNOTATION_COL)
}

# =========================================================
# Resolve sample column
# =========================================================

sample_candidates <- c(
  "sample",
  "sample_id",
  "Sample",
  "orig.ident"
)

sample_col <- sample_candidates[
  sample_candidates %in% colnames(whole@meta.data)
][1]

if (is.na(sample_col)) {
  stop(
    "Could not identify sample column. Available metadata columns:\n",
    paste(colnames(whole@meta.data), collapse=", ")
  )
}

cat("Sample column:", sample_col, "\n")

sample_vec <- as.character(
  whole@meta.data[[sample_col]]
)

cat(
  "Samples:",
  paste(sort(unique(sample_vec)), collapse=", "),
  "\n\n"
)

# =========================================================
# Parent Monocyte extraction
# =========================================================

annotation <- as.character(
  whole@meta.data[[ANNOTATION_COL]]
)

mon_cells <- colnames(whole)[
  annotation == TARGET_LABEL
]

cat("Annotated Monocyte cells:", length(mon_cells), "\n")

if (length(mon_cells) == 0) {
  cat("\nAvailable annotation labels:\n")
  print(sort(table(annotation), decreasing=TRUE))
  stop("No cells found with TARGET_LABEL = ", TARGET_LABEL)
}

mon_tmp <- subset(
  whole,
  cells=mon_cells
)

DefaultAssay(mon_tmp) <- "RNA"

# =========================================================
# Join RNA layers only within Monocyte subset
# =========================================================

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
        "No RNA counts layer found. Layers: ",
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
      slot="counts"
    )
  }

  counts
}

counts <- get_raw_counts(mon_tmp)

cat(
  "Raw counts matrix:",
  nrow(counts),
  "features x",
  ncol(counts),
  "cells\n\n"
)

# =========================================================
# Build clean raw-only Monocyte parent
# =========================================================

meta <- mon_tmp@meta.data[
  colnames(counts),
  ,
  drop=FALSE
]

mon_raw <- CreateSeuratObject(
  counts=counts,
  meta.data=meta,
  project="Mouse_MASH_Monocyte_v6.9.0"
)

saveRDS(
  mon_raw,
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_parent_raw_clean_v6.9.0.rds"
  )
)

# =========================================================
# Sample / condition mapping
# =========================================================

sample_raw <- as.character(
  mon_raw@meta.data[[sample_col]]
)

condition <- ifelse(
  grepl("^STD", sample_raw, ignore.case=TRUE),
  "STD",
  ifelse(
    grepl("CDHFD|CDAHFD", sample_raw, ignore.case=TRUE),
    "CDHFD",
    ifelse(
      grepl("^Sham", sample_raw, ignore.case=TRUE),
      "Sham",
      ifelse(
        grepl("^Tx", sample_raw, ignore.case=TRUE),
        "Tx",
        "Other"
      )
    )
  )
)

mon_raw$condition_v6.9.0 <- condition

sample_counts <- as.data.frame(
  table(
    sample=sample_raw
  ),
  stringsAsFactors=FALSE
)

condition_counts <- as.data.frame(
  table(
    condition=condition
  ),
  stringsAsFactors=FALSE
)

write.csv(
  sample_counts,
  file.path(
    TABDIR,
    "Monocyte_cell_counts_by_sample_v6.9.0.csv"
  ),
  row.names=FALSE
)

write.csv(
  condition_counts,
  file.path(
    TABDIR,
    "Monocyte_cell_counts_by_condition_v6.9.0.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Complexity by sample
# =========================================================

complexity <- do.call(
  rbind,
  lapply(
    sort(unique(sample_raw)),
    function(s) {

      cells <- colnames(mon_raw)[
        sample_raw == s
      ]

      md <- mon_raw@meta.data[
        cells,
        ,
        drop=FALSE
      ]

      data.frame(
        sample=s,
        n_cells=length(cells),

        median_nCount_RNA=
          median(
            md$nCount_RNA,
            na.rm=TRUE
          ),

        median_nFeature_RNA=
          median(
            md$nFeature_RNA,
            na.rm=TRUE
          ),

        q25_nCount_RNA=
          as.numeric(
            quantile(
              md$nCount_RNA,
              0.25,
              na.rm=TRUE
            )
          ),

        q75_nCount_RNA=
          as.numeric(
            quantile(
              md$nCount_RNA,
              0.75,
              na.rm=TRUE
            )
          ),

        q25_nFeature_RNA=
          as.numeric(
            quantile(
              md$nFeature_RNA,
              0.25,
              na.rm=TRUE
            )
          ),

        q75_nFeature_RNA=
          as.numeric(
            quantile(
              md$nFeature_RNA,
              0.75,
              na.rm=TRUE
            )
          ),

        stringsAsFactors=FALSE
      )
    }
  )
)

write.csv(
  complexity,
  file.path(
    TABDIR,
    "Monocyte_complexity_by_sample_v6.9.0.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Marker audit
# =========================================================

marker_sets <- list(

  Monocyte_core=c(
    "Lyz2",
    "Ly6c2",
    "Ccr2",
    "Fcgr3",
    "Lst1",
    "Tyrobp",
    "Ctss"
  ),

  Macrophage_Kupffer=c(
    "Adgre1",
    "Clec4f",
    "Timd4",
    "Vsig4",
    "Marco",
    "Cd5l",
    "Mertk"
  ),

  Neutrophil=c(
    "Cxcr2",
    "Retnlg",
    "Camp",
    "Ngp",
    "Mmp8",
    "S100a8",
    "S100a9"
  ),

  Dendritic=c(
    "Flt3",
    "Xcr1",
    "Clec10a",
    "Cd209a",
    "Itgax"
  ),

  Lymphoid=c(
    "Cd3d",
    "Cd3e",
    "Cd79a",
    "Ms4a1",
    "Nkg7"
  ),

  Endothelial=c(
    "Pecam1",
    "Kdr",
    "Emcn",
    "Klf2"
  ),

  Mesenchymal=c(
    "Col1a1",
    "Col1a2",
    "Col3a1",
    "Des"
  ),

  Epithelial=c(
    "Krt8",
    "Krt18",
    "Krt19",
    "Epcam"
  )
)

all_markers <- unique(
  unlist(
    marker_sets,
    use.names=FALSE
  )
)

present_markers <- intersect(
  all_markers,
  rownames(counts)
)

missing_markers <- setdiff(
  all_markers,
  rownames(counts)
)

cat(
  "Present audit markers:",
  length(present_markers),
  "/",
  length(all_markers),
  "\n"
)

if (length(missing_markers) > 0) {
  cat(
    "Missing markers:",
    paste(missing_markers, collapse=", "),
    "\n"
  )
}

# =========================================================
# Overall detection
# =========================================================

overall_rows <- list()

for (set_name in names(marker_sets)) {

  for (gene in marker_sets[[set_name]]) {

    if (!(gene %in% rownames(counts))) {
      next
    }

    x <- counts[
      gene,
      ,
      drop=TRUE
    ]

    overall_rows[[length(overall_rows) + 1]] <-
      data.frame(
        marker_set=set_name,
        gene=gene,
        detected_fraction=
          mean(x > 0),
        mean_raw_count=
          mean(x),
        median_raw_count=
          median(x),
        stringsAsFactors=FALSE
      )
  }
}

overall_marker_audit <- do.call(
  rbind,
  overall_rows
)

write.csv(
  overall_marker_audit,
  file.path(
    TABDIR,
    "Monocyte_marker_detection_overall_v6.9.0.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Detection by sample
# =========================================================

sample_marker_rows <- list()

for (s in sort(unique(sample_raw))) {

  cells <- colnames(mon_raw)[
    sample_raw == s
  ]

  for (set_name in names(marker_sets)) {

    for (gene in marker_sets[[set_name]]) {

      if (!(gene %in% rownames(counts))) {
        next
      }

      x <- counts[
        gene,
        cells,
        drop=TRUE
      ]

      sample_marker_rows[[length(sample_marker_rows) + 1]] <-
        data.frame(
          sample=s,
          marker_set=set_name,
          gene=gene,
          n_cells=length(cells),
          detected_fraction=
            mean(x > 0),
          mean_raw_count=
            mean(x),
          stringsAsFactors=FALSE
        )
    }
  }
}

sample_marker_audit <- do.call(
  rbind,
  sample_marker_rows
)

write.csv(
  sample_marker_audit,
  file.path(
    TABDIR,
    "Monocyte_marker_detection_by_sample_v6.9.0.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Simple marker-set hit burden per cell
#
# Diagnostic only.
# No cell exclusion.
# =========================================================

hit_sets <- c(
  "Monocyte_core",
  "Macrophage_Kupffer",
  "Neutrophil",
  "Dendritic",
  "Lymphoid",
  "Endothelial",
  "Mesenchymal",
  "Epithelial"
)

for (set_name in hit_sets) {

  genes <- intersect(
    marker_sets[[set_name]],
    rownames(counts)
  )

  if (length(genes) == 0) {
    mon_raw[[paste0(
      set_name,
      "_hits_v6.9.0"
    )]] <- 0
    next
  }

  hit_count <- Matrix::colSums(
    counts[
      genes,
      ,
      drop=FALSE
    ] > 0
  )

  mon_raw[[paste0(
    set_name,
    "_hits_v6.9.0"
  )]] <- hit_count[
    colnames(mon_raw)
  ]
}

saveRDS(
  mon_raw,
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_parent_raw_clean_v6.9.0.rds"
  )
)

# =========================================================
# Marker-set summary by sample
# =========================================================

hit_summary <- do.call(
  rbind,
  lapply(
    sort(unique(sample_raw)),
    function(s) {

      cells <- colnames(mon_raw)[
        sample_raw == s
      ]

      out <- data.frame(
        sample=s,
        n_cells=length(cells),
        stringsAsFactors=FALSE
      )

      for (set_name in hit_sets) {

        col_name <- paste0(
          set_name,
          "_hits_v6.9.0"
        )

        vals <- mon_raw@meta.data[
          cells,
          col_name
        ]

        out[[paste0(
          set_name,
          "_median_hits"
        )]] <- median(
          vals,
          na.rm=TRUE
        )

        out[[paste0(
          set_name,
          "_fraction_ge2"
        )]] <- mean(
          vals >= 2,
          na.rm=TRUE
        )
      }

      out
    }
  )
)

write.csv(
  hit_summary,
  file.path(
    TABDIR,
    "Monocyte_marker_set_hit_summary_by_sample_v6.9.0.csv"
  ),
  row.names=FALSE
)

# =========================================================
# Terminal summary
# =========================================================

cat("\n=== MONOCYTE CELLS BY SAMPLE ===\n")
print(
  sample_counts,
  row.names=FALSE
)

cat("\n=== MONOCYTE CELLS BY CONDITION ===\n")
print(
  condition_counts,
  row.names=FALSE
)

cat("\n=== COMPLEXITY BY SAMPLE ===\n")
print(
  complexity,
  row.names=FALSE
)

cat("\n=== MARKER DETECTION OVERALL ===\n")
print(
  overall_marker_audit,
  row.names=FALSE
)

cat("\n=== MARKER-SET HIT SUMMARY BY SAMPLE ===\n")
print(
  hit_summary,
  row.names=FALSE
)

cat("\n====================================================\n")
cat("v6.9.0 COMPLETE\n")
cat("Monocyte parent audit complete\n")
cat("No cells removed\n")
cat("No annotation changed\n")
cat(
  "Clean raw parent:",
  file.path(
    OBJDIR,
    "Mouse_MASH_Monocyte_parent_raw_clean_v6.9.0.rds"
  ),
  "\n"
)
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.9.0.txt"
  )
)
