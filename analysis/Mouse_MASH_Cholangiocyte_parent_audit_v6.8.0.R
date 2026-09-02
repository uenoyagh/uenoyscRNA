suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

set.seed(20260902)

VERSION <- "v6.8.0"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "RDS3_annotation_visualization_v4.1.1/objects/",
  "RDS3_with_visualization_metadata_v4.1.1.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

ANNO_COL <- "celltype_for_R8plot_FIXED2"
CHOL_LABEL <- "Cholangiocyte"
CYCLING_LABEL <- "Cycling"

cat("====================================================\n")
cat("Mouse MASH Cholangiocyte parent audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

# ---------------------------------------------------------
# Helper: raw RNA counts
# ---------------------------------------------------------

get_counts <- function(obj, assay="RNA") {

  DefaultAssay(obj) <- assay

  if (inherits(obj[[assay]], "Assay5")) {

    layers <- Layers(obj[[assay]])

    count_layers <- grep(
      "^counts($|\\.)",
      layers,
      value=TRUE
    )

    if (length(count_layers) == 0) {
      stop("No counts layer found.")
    }

    if (
      length(count_layers) == 1 &&
      identical(count_layers, "counts")
    ) {

      return(
        LayerData(
          obj,
          assay=assay,
          layer="counts"
        )
      )
    }

    obj <- JoinLayers(
      obj,
      assay=assay
    )

    return(
      LayerData(
        obj,
        assay=assay,
        layer="counts"
      )
    )

  } else {

    return(
      GetAssayData(
        obj,
        assay=assay,
        layer="counts"
      )
    )
  }
}

# ---------------------------------------------------------
# Read whole-cell RDS
# ---------------------------------------------------------

cat("Reading:\n", INPUT_RDS, "\n\n")

whole <- readRDS(INPUT_RDS)

cat("=== WHOLE-CELL OBJECT ===\n")
cat("Cells:", ncol(whole), "\n")
cat("Features:", nrow(whole), "\n")
cat("Assays:", paste(Assays(whole), collapse=", "), "\n\n")

required <- c(
  ANNO_COL,
  "sample",
  "condition"
)

missing <- setdiff(
  required,
  colnames(whole@meta.data)
)

if (length(missing) > 0) {
  stop(
    "Missing metadata: ",
    paste(missing, collapse=", ")
  )
}

cat("=== WHOLE-CELL ANNOTATIONS ===\n")
print(
  sort(
    table(
      whole@meta.data[[ANNO_COL]]
    ),
    decreasing=TRUE
  )
)

anno <- as.character(
  whole@meta.data[[ANNO_COL]]
)

if (!CHOL_LABEL %in% anno) {
  stop(
    "Exact Cholangiocyte label not found."
  )
}

# ---------------------------------------------------------
# Cholangiocyte parent
# ---------------------------------------------------------

chol_cells <- colnames(whole)[
  anno == CHOL_LABEL
]

cat("\n=== CHOLANGIOCYTE PARENT ===\n")
cat("Parent cells:", length(chol_cells), "\n")

chol_old <- subset(
  whole,
  cells=chol_cells
)

cat("\n=== CELLS BY SAMPLE ===\n")
print(
  table(chol_old$sample)
)

cat("\n=== CELLS BY CONDITION ===\n")
print(
  table(chol_old$condition)
)

write.csv(
  as.data.frame(
    table(
      sample=chol_old$sample
    )
  ),
  file.path(
    TABDIR,
    "Cholangiocyte_parent_cells_by_sample_v6.8.0.csv"
  ),
  row.names=FALSE
)

write.csv(
  as.data.frame(
    table(
      condition=chol_old$condition
    )
  ),
  file.path(
    TABDIR,
    "Cholangiocyte_parent_cells_by_condition_v6.8.0.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Extract raw RNA counts
# ---------------------------------------------------------

cat("\n=== EXTRACT RAW RNA COUNTS ===\n")

chol_counts <- get_counts(
  chol_old,
  assay="RNA"
)

chol_md <- chol_old@meta.data[
  colnames(chol_counts),
  ,
  drop=FALSE
]

cat(
  "Counts:",
  nrow(chol_counts),
  "genes x",
  ncol(chol_counts),
  "cells\n"
)

# ---------------------------------------------------------
# Build completely clean parent
# ---------------------------------------------------------

cat("\n=== CREATE CLEAN CHOLANGIOCYTE OBJECT ===\n")

options(
  Seurat.object.assay.version="v3"
)

chol <- CreateSeuratObject(
  counts=chol_counts,
  meta.data=chol_md,
  project="Mouse_MASH_Cholangiocyte"
)

rm(chol_old, chol_counts)
invisible(gc())

cat("Cells:", ncol(chol), "\n")
cat("Features:", nrow(chol), "\n")
cat("Assays:", paste(Assays(chol), collapse=", "), "\n")
cat("RNA assay class:", class(chol[["RNA"]])[1], "\n")
cat("Reductions:", length(Reductions(chol)), "\n")
cat("Graphs:", length(chol@graphs), "\n")

if (!identical(Assays(chol), "RNA")) {
  stop("Unexpected assay in clean Cholangiocyte object.")
}

if (!inherits(chol[["RNA"]], "Assay")) {
  stop("RNA assay is not legacy Assay.")
}

if (length(Reductions(chol)) != 0) {
  stop("Clean object unexpectedly contains reductions.")
}

if (length(chol@graphs) != 0) {
  stop("Clean object unexpectedly contains graphs.")
}

cat("\nCLEAN CHOLANGIOCYTE CHECK: PASSED\n")

# ---------------------------------------------------------
# Marker panels
# ---------------------------------------------------------

marker_sets <- list(

  Cholangiocyte_identity = c(
    "Krt19","Krt7","Krt8","Krt18",
    "Epcam","Sox9","Klf5","Muc1",
    "Hnf1b","Cftr","Slc4a2"
  ),

  Ductular_reactive = c(
    "Spp1","Mmp7","Krt19","Krt7",
    "Krt23","Sox9","Epcam",
    "Tacstd2","Prom1"
  ),

  Inflammatory = c(
    "Icam1","Vcam1",
    "Cxcl1","Cxcl2","Ccl2",
    "Il6","Nfkbia","Socs3"
  ),

  Fibrogenic_signaling = c(
    "Tgfb1","Tgfb2",
    "Ccn2","Jag1",
    "Serpine1","Thbs1"
  ),

  Senescence_stress = c(
    "Cdkn1a","Cdkn2a",
    "Trp53","Gadd45a",
    "Fos","Jun","Atf3"
  ),

  Cycling = c(
    "Mki67","Top2a","Birc5",
    "Ube2c","Cenpf","Pcna"
  ),

  Hepatocyte = c(
    "Alb","Ttr","Apoa1",
    "Fabp1","Cps1","Ass1",
    "Hnf4a"
  ),

  Myeloid = c(
    "Ptprc","Lyz2","Adgre1",
    "Tyrobp","Aif1"
  ),

  Endothelial = c(
    "Pecam1","Cdh5","Kdr",
    "Klf2","Stab2"
  ),

  HSC_mesenchymal = c(
    "Col1a1","Col3a1",
    "Pdgfra","Pdgfrb",
    "Lrat","Rbp1","Rgs5"
  ),

  Lymphoid = c(
    "Cd3d","Cd3e",
    "Nkg7","Ms4a1"
  ),

  Neutrophil = c(
    "S100a8","S100a9",
    "Retnlg","Mmp8"
  )
)

marker_present <- lapply(
  marker_sets,
  function(x) intersect(x, rownames(chol))
)

marker_missing <- lapply(
  marker_sets,
  function(x) setdiff(x, rownames(chol))
)

cat("\n=== MARKER AVAILABILITY ===\n")

for (nm in names(marker_sets)) {

  cat("\n[", nm, "]\n", sep="")

  cat(
    "Present: ",
    paste(
      marker_present[[nm]],
      collapse=", "
    ),
    "\n",
    sep=""
  )

  cat(
    "Missing: ",
    paste(
      marker_missing[[nm]],
      collapse=", "
    ),
    "\n",
    sep=""
  )
}

capture.output(
  {
    for (nm in names(marker_sets)) {
      cat("[", nm, "]\n", sep="")
      cat(
        "Present: ",
        paste(marker_present[[nm]], collapse=", "),
        "\n",
        sep=""
      )
      cat(
        "Missing: ",
        paste(marker_missing[[nm]], collapse=", "),
        "\n\n",
        sep=""
      )
    }
  },
  file=file.path(
    TABDIR,
    "Cholangiocyte_marker_availability_v6.8.0.txt"
  )
)

# ---------------------------------------------------------
# Normalize parent for descriptive marker audit
# ---------------------------------------------------------

DefaultAssay(chol) <- "RNA"

chol <- NormalizeData(
  chol,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

Idents(chol) <- factor(
  chol$sample
)

dot_genes <- unique(
  unlist(
    marker_present,
    use.names=FALSE
  )
)

p_dot <- DotPlot(
  chol,
  features=dot_genes,
  assay="RNA",
  dot.scale=6
) +
  RotatedAxis() +
  theme_classic(base_size=9) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
  ) +
  ggtitle(
    "Cholangiocyte parent marker audit by sample"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_parent_marker_DotPlot_by_sample_v6.8.0.pdf"
  ),
  p_dot,
  width=21,
  height=6.5
)

# ---------------------------------------------------------
# Average expression by sample
# ---------------------------------------------------------

avg <- AverageExpression(
  chol,
  assays="RNA",
  features=dot_genes,
  group.by="sample",
  slot="data",
  verbose=FALSE
)$RNA

write.csv(
  avg,
  file.path(
    TABDIR,
    "Cholangiocyte_parent_marker_AverageExpression_by_sample_v6.8.0.csv"
  )
)

# ---------------------------------------------------------
# Cycling-cell rescue audit
#
# Diagnostic only.
# No Cycling cells are added to the parent in v6.8.0.
# ---------------------------------------------------------

cat("\n=== CYCLING CHOLANGIOCYTE RESCUE AUDIT ===\n")

cycling_cells <- colnames(whole)[
  anno == CYCLING_LABEL
]

cat(
  "Whole-cell Cycling cells:",
  length(cycling_cells),
  "\n"
)

if (length(cycling_cells) > 0) {

  cyc <- subset(
    whole,
    cells=cycling_cells
  )

  cyc_counts <- get_counts(
    cyc,
    assay="RNA"
  )

  detected <- function(gene) {

    if (!gene %in% rownames(cyc_counts)) {
      return(
        rep(
          FALSE,
          ncol(cyc_counts)
        )
      )
    }

    as.vector(
      cyc_counts[gene, ] > 0
    )
  }

  biliary_specific_hits <-
    detected("Krt19") +
    detected("Epcam") +
    detected("Sox9") +
    detected("Krt7") +
    detected("Muc1")

  epithelial_support_hits <-
    detected("Krt8") +
    detected("Krt18")

  rescue_candidate <-
    biliary_specific_hits >= 1 &
    epithelial_support_hits >= 1

  cyc_audit <- data.frame(
    cell=colnames(cyc_counts),
    sample=as.character(
      cyc$sample[
        colnames(cyc_counts)
      ]
    ),
    condition=as.character(
      cyc$condition[
        colnames(cyc_counts)
      ]
    ),
    biliary_specific_hits=
      biliary_specific_hits,
    epithelial_support_hits=
      epithelial_support_hits,
    rescue_candidate=
      rescue_candidate,
    stringsAsFactors=FALSE
  )

  write.csv(
    cyc_audit,
    file.path(
      TABDIR,
      "Cycling_biliary_rescue_audit_per_cell_v6.8.0.csv"
    ),
    row.names=FALSE
  )

  rescue_by_sample <- as.data.frame(
    table(
      sample=cyc_audit$sample,
      rescue_candidate=
        cyc_audit$rescue_candidate
    )
  )

  write.csv(
    rescue_by_sample,
    file.path(
      TABDIR,
      "Cycling_biliary_rescue_candidates_by_sample_v6.8.0.csv"
    ),
    row.names=FALSE
  )

  cat(
    "Cycling cells meeting diagnostic biliary gate:",
    sum(rescue_candidate),
    "\n"
  )

  cat("\nBy sample:\n")

  print(
    table(
      cyc_audit$sample,
      cyc_audit$rescue_candidate
    )
  )

  rm(cyc, cyc_counts)
  invisible(gc())

} else {

  cat("No Cycling cells found.\n")
}

# ---------------------------------------------------------
# Save clean parent
# ---------------------------------------------------------

chol$Cholangiocyte_parent_status_v680 <-
  "Original_Cholangiocyte_annotation"

saveRDS(
  chol,
  file.path(
    OBJDIR,
    "Mouse_MASH_Cholangiocyte_parent_raw_clean_v6.8.0.rds"
  ),
  compress=FALSE
)

# ---------------------------------------------------------
# Summary
# ---------------------------------------------------------

summary_lines <- c(
  "# Mouse MASH Cholangiocyte parent audit v6.8.0",
  "",
  paste0("- Parent cells: ", ncol(chol)),
  paste0("- Features: ", nrow(chol)),
  "- Parent definition: original whole-cell annotation = Cholangiocyte",
  "- Object rebuilt from raw RNA counts",
  "- RNA assay only",
  "- No inherited reductions",
  "- No inherited graphs",
  "- Cycling cells were audited separately for possible biliary rescue",
  "- No Cycling cells were added to the parent in v6.8.0",
  "- No Cholangiocyte subclustering performed in v6.8.0"
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "Cholangiocyte_parent_audit_summary_v6.8.0.md"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.0.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.0 COMPLETE\n")
cat("Cholangiocyte parent cells:", ncol(chol), "\n")
cat("Features:", nrow(chol), "\n")
cat("Assays:", paste(Assays(chol), collapse=", "), "\n")
cat("Reductions:", length(Reductions(chol)), "\n")
cat("Graphs:", length(chol@graphs), "\n")
cat("No subclustering performed yet\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
