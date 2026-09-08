# ============================================================
# Human MASH v1.3.1
# Whole-liver annotation validation
#
# PURPOSE
#   Validate current REF12-derived whole-liver annotation
#   before freezing annotation and beginning disease DE.
#
# IMPORTANT
#   - Input canonical RDS is READ ONLY.
#   - No modified Seurat object is written.
#   - IL10-axis output in this script is feasibility/exploratory only.
#
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

input_rds <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "objects/",
  "Human_MASH_LT1to8_Canonical_v1.3.rds"
)

outdir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/",
  "Human_MASH_LT1to8_Canonical_v1.3/",
  "analysis/",
  "v1.3.1_whole_liver_annotation_validation"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

celltype_col <- "celltype_ref12_v13"
sample_col   <- "canonical_sample_v13"
nas_col      <- "NAS_score_v13"
nasgrp_col   <- "NAS_group_v13"

umap_native <- "umap"
umap_ref    <- "umapRPCA_ref12_v13"

# ------------------------------------------------------------
# 1. Load + invariant checks
# ------------------------------------------------------------

cat("============================================================\n")
cat("Human MASH v1.3.1 Whole-liver Annotation Validation\n")
cat("============================================================\n\n")

stopifnot(file.exists(input_rds))

obj <- readRDS(input_rds)

required_meta <- c(
  celltype_col,
  sample_col,
  nas_col,
  nasgrp_col
)

missing_meta <- setdiff(required_meta, colnames(obj@meta.data))

if (length(missing_meta) > 0) {
  stop(
    "Missing required metadata: ",
    paste(missing_meta, collapse = ", ")
  )
}

required_red <- c(umap_native, umap_ref)
missing_red <- setdiff(required_red, Reductions(obj))

if (length(missing_red) > 0) {
  stop(
    "Missing required reductions: ",
    paste(missing_red, collapse = ", ")
  )
}

cat("Cells      :", ncol(obj), "\n")
cat("RNA genes  :", nrow(obj[["RNA"]]), "\n")
cat("Reductions :", paste(Reductions(obj), collapse = " | "), "\n\n")

# ------------------------------------------------------------
# 2. Local vectors only
#    Do not alter metadata in the Seurat object
# ------------------------------------------------------------

celltype_vec <- as.character(obj[[celltype_col]][, 1])
sample_vec   <- as.character(obj[[sample_col]][, 1])
nas_vec      <- as.character(obj[[nas_col]][, 1])
nasgrp_vec   <- as.character(obj[[nasgrp_col]][, 1])

celltype_vec[
  is.na(celltype_vec) | celltype_vec == ""
] <- "NA"

sample_vec[
  is.na(sample_vec) | sample_vec == ""
] <- "NA"

nas_vec[
  is.na(nas_vec) | nas_vec == ""
] <- "NA"

nasgrp_vec[
  is.na(nasgrp_vec) | nasgrp_vec == ""
] <- "NA"

# Natural LT ordering when possible
sample_order_fun <- function(x) {

  u <- unique(x)

  if (
    length(u) > 0 &&
    all(grepl("^LT[0-9]+", u))
  ) {

    number <- suppressWarnings(
      as.integer(
        sub("^LT([0-9]+).*", "\\1", u)
      )
    )

    if (!anyNA(number)) {
      return(u[order(number)])
    }
  }

  sort(u)
}

sample_order <- sample_order_fun(sample_vec)

# ------------------------------------------------------------
# 3. Sample / NAS mapping
# ------------------------------------------------------------

sample_nas <- unique(
  data.frame(
    sample = sample_vec,
    NAS_score = nas_vec,
    NAS_group = nasgrp_vec,
    stringsAsFactors = FALSE
  )
)

sample_nas$.order <- match(
  sample_nas$sample,
  sample_order
)

sample_nas <- sample_nas[
  order(
    sample_nas$.order,
    sample_nas$NAS_score,
    sample_nas$NAS_group
  ),
]

sample_nas$.order <- NULL

write.csv(
  sample_nas,
  file.path(
    outdir,
    "01_sample_NAS_mapping.csv"
  ),
  row.names = FALSE
)

cat("=== SAMPLE / NAS MAPPING ===\n")
print(sample_nas, row.names = FALSE)
cat("\n")

# ------------------------------------------------------------
# 4. Current cell-type counts
# ------------------------------------------------------------

ct_counts <- as.data.frame(
  table(celltype = celltype_vec),
  stringsAsFactors = FALSE
)

names(ct_counts)[2] <- "n_cells"

ct_counts$percent <- round(
  100 * ct_counts$n_cells /
    sum(ct_counts$n_cells),
  4
)

ct_counts <- ct_counts[
  order(ct_counts$n_cells, decreasing = TRUE),
]

write.csv(
  ct_counts,
  file.path(
    outdir,
    "02_current_celltype_counts.csv"
  ),
  row.names = FALSE
)

cat("=== CURRENT CELL TYPES ===\n")
print(ct_counts, row.names = FALSE)
cat("\n")

# ------------------------------------------------------------
# 5. Cell type x sample composition
# ------------------------------------------------------------

ct_sample <- as.data.frame(
  table(
    celltype = celltype_vec,
    sample = sample_vec
  ),
  stringsAsFactors = FALSE
)

ct_sample$sample_total <- ave(
  ct_sample$Freq,
  ct_sample$sample,
  FUN = sum
)

ct_sample$celltype_total <- ave(
  ct_sample$Freq,
  ct_sample$celltype,
  FUN = sum
)

ct_sample$within_sample_pct <- ifelse(
  ct_sample$sample_total > 0,
  100 * ct_sample$Freq /
    ct_sample$sample_total,
  NA_real_
)

ct_sample$within_celltype_pct <- ifelse(
  ct_sample$celltype_total > 0,
  100 * ct_sample$Freq /
    ct_sample$celltype_total,
  NA_real_
)

ct_sample$sample <- factor(
  ct_sample$sample,
  levels = sample_order
)

write.csv(
  ct_sample,
  file.path(
    outdir,
    "03_celltype_by_sample_long.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Annotation provenance / crosswalk
# ------------------------------------------------------------

auto_col <- "auto_celltype_final_ref12_v13"
basis_col <- "celltype_annotation_basis_ref12_v13"

if (auto_col %in% colnames(obj@meta.data)) {

  auto_vec <- as.character(
    obj[[auto_col]][, 1]
  )

  auto_vec[
    is.na(auto_vec) | auto_vec == ""
  ] <- "NA"

  xw <- as.data.frame.matrix(
    table(
      Current = celltype_vec,
      Auto = auto_vec
    )
  )

  xw_out <- data.frame(
    celltype_ref12_v13 = rownames(xw),
    xw,
    check.names = FALSE,
    row.names = NULL
  )

  write.csv(
    xw_out,
    file.path(
      outdir,
      "04_current_vs_auto_annotation_crosswalk.csv"
    ),
    row.names = FALSE
  )
}

if (basis_col %in% colnames(obj@meta.data)) {

  basis_vec <- as.character(
    obj[[basis_col]][, 1]
  )

  basis_vec[
    is.na(basis_vec) | basis_vec == ""
  ] <- "NA"

  basis_tab <- as.data.frame.matrix(
    table(
      Current = celltype_vec,
      Basis = basis_vec
    )
  )

  basis_out <- data.frame(
    celltype_ref12_v13 = rownames(basis_tab),
    basis_tab,
    check.names = FALSE,
    row.names = NULL
  )

  write.csv(
    basis_out,
    file.path(
      outdir,
      "05_annotation_basis_by_celltype.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 7. Marker panel
#
# Includes lineages missing from current annotation
# to test whether they are genuinely absent or misassigned.
# ------------------------------------------------------------

marker_list <- list(

  Hepatocyte = c(
    "ALB",
    "APOA1",
    "TTR",
    "ASGR1",
    "CYP3A4",
    "CYP2E1"
  ),

  Immature_hepatocyte = c(
    "AFP",
    "DLK1",
    "PROM1",
    "EPCAM",
    "KRT8",
    "KRT18"
  ),

  Cholangiocyte = c(
    "KRT7",
    "KRT19",
    "SOX9",
    "MUC1",
    "CFTR"
  ),

  LSEC = c(
    "CLEC4G",
    "FCGR2B",
    "STAB2",
    "LYVE1",
    "PLVAP",
    "EMCN"
  ),

  Vascular_endothelial = c(
    "PECAM1",
    "VWF",
    "KDR",
    "ENG",
    "ACKR1",
    "RAMP2"
  ),

  Quiescent_HSC = c(
    "RBP1",
    "LRAT",
    "CYGB",
    "RELN",
    "RGS5",
    "PDGFRB"
  ),

  Activated_HSC_fibroblast = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "ACTA2"
  ),

  Monocyte_macrophage = c(
    "LST1",
    "TYROBP",
    "FCER1G",
    "CTSS",
    "LILRB1",
    "FCN1"
  ),

  Resident_macrophage = c(
    "C1QC",
    "C1QB",
    "MARCO",
    "VSIG4",
    "CD68",
    "MSR1"
  ),

  SPP1_TREM2_macrophage = c(
    "SPP1",
    "TREM2",
    "GPNMB",
    "LPL",
    "APOE",
    "CD9"
  ),

  T_NK = c(
    "CD3D",
    "CD3E",
    "TRBC1",
    "NKG7",
    "GNLY",
    "KLRD1"
  ),

  B_Plasma = c(
    "MS4A1",
    "CD79A",
    "CD37",
    "JCHAIN",
    "MZB1",
    "SDC1"
  ),

  Dendritic = c(
    "FCER1A",
    "CD1C",
    "CLEC10A",
    "HLA-DRA",
    "LILRA4"
  ),

  Mesothelial = c(
    "MSLN",
    "KRT17",
    "KRT23",
    "UPK3B",
    "WT1",
    "CAV1"
  ),

  Cycling = c(
    "MKI67",
    "TOP2A",
    "STMN1",
    "TUBA1B",
    "UBE2C",
    "CENPF"
  )
)

marker_availability <- do.call(
  rbind,
  lapply(
    names(marker_list),
    function(grp) {

      genes <- marker_list[[grp]]

      data.frame(
        marker_group = grp,
        gene = genes,
        present_in_RNA = genes %in%
          rownames(obj),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(marker_availability) <- NULL

write.csv(
  marker_availability,
  file.path(
    outdir,
    "06_marker_gene_availability.csv"
  ),
  row.names = FALSE
)

marker_list_present <- lapply(
  marker_list,
  function(x) {
    intersect(x, rownames(obj))
  }
)

marker_list_present <- marker_list_present[
  lengths(marker_list_present) > 0
]

# ------------------------------------------------------------
# 8. IL10-axis feasibility
# ------------------------------------------------------------

il10_genes <- c(
  "IL10",
  "IL10RA",
  "IL10RB",
  "JAK1",
  "TYK2",
  "STAT3",
  "SOCS3"
)

il10_availability <- data.frame(
  gene = il10_genes,
  present_in_RNA =
    il10_genes %in% rownames(obj),
  stringsAsFactors = FALSE
)

write.csv(
  il10_availability,
  file.path(
    outdir,
    "07_IL10_axis_gene_availability.csv"
  ),
  row.names = FALSE
)

cat("=== IL10 AXIS GENE AVAILABILITY ===\n")
print(il10_availability, row.names = FALSE)
cat("\n")

# ------------------------------------------------------------
# 9. High-saturation cell-type palette
# ------------------------------------------------------------

palette_known <- c(
  "Hepatocyte" =
    "#00A6A6",

  "Hepatoblast_or_immature_hepatocyte" =
    "#0072B2",

  "Cholangiocyte_or_ductal" =
    "#00C853",

  "LSEC" =
    "#00CFEA",

  "Endothelial" =
    "#0057B8",

  "HSC_quiescent" =
    "#F06292",

  "HSC_fibroblast_mesenchymal" =
    "#FF1493",

  "Fibroblast_mesenchymal" =
    "#B5179E",

  "Inflammatory_monocyte_macrophage" =
    "#E31A1C",

  "SPP1_TREM2_macrophage" =
    "#FF7F00",

  "Plasma_cell" =
    "#7B2CBF",

  "Cycling_cell" =
    "#FFD600",

  "Mesothelial" =
    "#8D6E63"
)

current_celltypes <- unique(celltype_vec)

missing_colors <- setdiff(
  current_celltypes,
  names(palette_known)
)

if (length(missing_colors) > 0) {

  extra_cols <- setNames(
    grDevices::hcl.colors(
      length(missing_colors),
      palette = "Dark 3"
    ),
    missing_colors
  )

  celltype_palette <- c(
    palette_known,
    extra_cols
  )

} else {

  celltype_palette <- palette_known
}

celltype_palette <- celltype_palette[
  current_celltypes
]

# ------------------------------------------------------------
# 10. UMAP comparison
# ------------------------------------------------------------

p_umap_native <- DimPlot(
  obj,
  reduction = umap_native,
  group.by = celltype_col,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  pt.size = 0.08,
  raster = TRUE,
  cols = celltype_palette
) +
  ggtitle("Current annotation on native UMAP") +
  theme_classic(base_size = 11)

p_umap_ref <- DimPlot(
  obj,
  reduction = umap_ref,
  group.by = celltype_col,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  pt.size = 0.08,
  raster = TRUE,
  cols = celltype_palette
) +
  ggtitle("Current annotation on RPCA reference UMAP") +
  theme_classic(base_size = 11)

p_umap_compare <- patchwork::wrap_plots(
  p_umap_native,
  p_umap_ref,
  ncol = 2
)

ggsave(
  file.path(
    outdir,
    "Fig01_CurrentAnnotation_native_vs_RPCA_UMAP.pdf"
  ),
  p_umap_compare,
  width = 17,
  height = 8
)

# ------------------------------------------------------------
# 11. Sample-colored UMAP
# ------------------------------------------------------------

p_sample_umap <- DimPlot(
  obj,
  reduction = umap_ref,
  group.by = sample_col,
  pt.size = 0.08,
  raster = TRUE,
  shuffle = TRUE
) +
  ggtitle("Canonical samples on RPCA reference UMAP") +
  theme_classic(base_size = 11)

ggsave(
  file.path(
    outdir,
    "Fig02_Sample_RPCA_UMAP.pdf"
  ),
  p_sample_umap,
  width = 10,
  height = 8
)

# ------------------------------------------------------------
# 12. Split UMAP by sample
# ------------------------------------------------------------

p_split <- DimPlot(
  obj,
  reduction = umap_ref,
  group.by = celltype_col,
  split.by = sample_col,
  pt.size = 0.08,
  raster = TRUE,
  cols = celltype_palette,
  ncol = 4
)

ggsave(
  file.path(
    outdir,
    "Fig03_Celltype_RPCA_UMAP_split_by_sample.pdf"
  ),
  p_split,
  width = 18,
  height = 12
)

# ------------------------------------------------------------
# 13. Whole-liver marker DotPlot
# ------------------------------------------------------------

p_dot <- DotPlot(
  obj,
  features = marker_list_present,
  assay = "RNA",
  group.by = celltype_col,
  dot.scale = 6
) +
  RotatedAxis() +
  scale_color_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  labs(
    title =
      "Whole-liver annotation validation marker panel",
    color = "Scaled\nexpression",
    size = "% expressed"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1,
        vjust = 1
      ),
    axis.text.y =
      element_text(size = 8)
  )

ggsave(
  file.path(
    outdir,
    "Fig04_WholeLiver_AnnotationValidation_DotPlot.pdf"
  ),
  p_dot,
  width = 24,
  height = 10
)

write.csv(
  p_dot$data,
  file.path(
    outdir,
    "08_WholeLiver_DotPlot_data.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Sentinel FeaturePlots
# ------------------------------------------------------------

sentinel_genes <- c(
  "ALB",
  "AFP",
  "EPCAM",
  "KRT19",
  "CLEC4G",
  "STAB2",
  "PECAM1",
  "RBP1",
  "COL1A1",
  "LST1",
  "C1QC",
  "SPP1",
  "TREM2",
  "CD3D",
  "NKG7",
  "MS4A1",
  "JCHAIN",
  "MKI67"
)

sentinel_present <- intersect(
  sentinel_genes,
  rownames(obj)
)

if (length(sentinel_present) > 0) {

  p_feature <- FeaturePlot(
    obj,
    features = sentinel_present,
    reduction = umap_ref,
    ncol = 4,
    pt.size = 0.08,
    raster = TRUE,
    order = FALSE,
    cols = c(
      "#F2F2F2",
      "#FF1A1A"
    )
  )

  ggsave(
    file.path(
      outdir,
      "Fig05_SentinelMarker_FeaturePlots.pdf"
    ),
    p_feature,
    width = 16,
    height = 18
  )
}

# ------------------------------------------------------------
# 15. Sample composition plot
# ------------------------------------------------------------

p_comp <- ggplot(
  ct_sample,
  aes(
    x = sample,
    y = within_sample_pct,
    fill = celltype
  )
) +
  geom_col(width = 0.85) +
  scale_fill_manual(
    values = celltype_palette
  ) +
  labs(
    title =
      "Whole-liver cell-type composition by sample",
    x = NULL,
    y = "% of cells within sample",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    outdir,
    "Fig06_CelltypeComposition_by_sample.pdf"
  ),
  p_comp,
  width = 12,
  height = 8
)

# ------------------------------------------------------------
# 16. Sample-bias heatmap
#     Percentage contribution of each sample to each cell type
# ------------------------------------------------------------

p_bias <- ggplot(
  ct_sample,
  aes(
    x = sample,
    y = celltype,
    fill = within_celltype_pct
  )
) +
  geom_tile(
    linewidth = 0.25,
    color = "grey80"
  ) +
  geom_text(
    aes(
      label = ifelse(
        Freq > 0,
        sprintf(
          "%.1f",
          within_celltype_pct
        ),
        ""
      )
    ),
    size = 2.5
  ) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#FF1A1A"
  ) +
  labs(
    title =
      "Sample contribution within each current cell type",
    subtitle =
      "High values identify sample-dominated annotations",
    x = NULL,
    y = NULL,
    fill = "%"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  file.path(
    outdir,
    "Fig07_SampleContribution_Celltype_QC.pdf"
  ),
  p_bias,
  width = 12,
  height = 8
)

# ------------------------------------------------------------
# 17. Exploratory IL10-axis DotPlot
#     Do NOT use for biological conclusion before annotation freeze.
# ------------------------------------------------------------

il10_present <- intersect(
  il10_genes,
  rownames(obj)
)

if (length(il10_present) > 0) {

  p_il10 <- DotPlot(
    obj,
    features = il10_present,
    group.by = celltype_col,
    assay = "RNA",
    dot.scale = 8
  ) +
    RotatedAxis() +
    scale_color_gradient2(
      low = "#0033FF",
      mid = "#FFFFFF",
      high = "#FF1A1A",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish
    ) +
    labs(
      title =
        "Exploratory IL10 / IL10R pathway",
      subtitle =
        "Interpret only after annotation validation",
      color = "Scaled\nexpression",
      size = "% expressed"
    ) +
    theme_bw(base_size = 10)

  ggsave(
    file.path(
      outdir,
      "Fig08_Exploratory_IL10_axis_DotPlot.pdf"
    ),
    p_il10,
    width = 10,
    height = 7
  )

  write.csv(
    p_il10$data,
    file.path(
      outdir,
      "09_Exploratory_IL10_axis_DotPlot_data.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 18. Automated QC flags
# ------------------------------------------------------------

flag_lines <- c(
  "Human MASH v1.3.1 annotation validation flags",
  "==============================================",
  "",
  paste0(
    "Cells: ",
    ncol(obj)
  ),
  paste0(
    "RNA genes: ",
    nrow(obj[["RNA"]])
  ),
  ""
)

small_types <- ct_counts[
  ct_counts$n_cells < 50,
]

if (nrow(small_types) > 0) {

  flag_lines <- c(
    flag_lines,
    "FLAG: current annotations with <50 cells:",
    paste0(
      "  ",
      small_types$celltype,
      " = ",
      small_types$n_cells,
      " cells (",
      small_types$percent,
      "%)"
    ),
    ""
  )
}

immature_row <- ct_counts[
  ct_counts$celltype ==
    "Hepatoblast_or_immature_hepatocyte",
]

if (
  nrow(immature_row) == 1 &&
  immature_row$percent > 5
) {

  flag_lines <- c(
    flag_lines,
    paste0(
      "FLAG: Hepatoblast_or_immature_hepatocyte = ",
      immature_row$n_cells,
      " cells (",
      immature_row$percent,
      "%); marker validation required."
    ),
    ""
  )
}

label_names <- as.character(
  ct_counts$celltype
)

if (
  !any(
    grepl(
      "T_cell|Tcell",
      label_names,
      ignore.case = TRUE
    )
  )
) {

  flag_lines <- c(
    flag_lines,
    "FLAG: no explicit T-cell annotation label.",
    ""
  )
}

if (
  !any(
    grepl(
      "NK|natural_killer",
      label_names,
      ignore.case = TRUE
    )
  )
) {

  flag_lines <- c(
    flag_lines,
    "FLAG: no explicit NK-cell annotation label.",
    ""
  )
}

if (
  !any(
    grepl(
      "B_cell|Bcell",
      label_names,
      ignore.case = TRUE
    )
  )
) {

  flag_lines <- c(
    flag_lines,
    "FLAG: no explicit B-cell annotation label (Plasma_cell is present).",
    ""
  )
}

if (
  !any(
    grepl(
      "dendritic",
      label_names,
      ignore.case = TRUE
    )
  )
) {

  flag_lines <- c(
    flag_lines,
    "FLAG: no explicit dendritic-cell annotation label.",
    ""
  )
}

# Sample dominance within each cell type
split_idx <- split(
  seq_len(nrow(ct_sample)),
  ct_sample$celltype
)

max_idx <- unlist(
  lapply(
    split_idx,
    function(idx) {
      idx[
        which.max(
          ct_sample$within_celltype_pct[idx]
        )
      ]
    }
  )
)

dominance <- ct_sample[
  max_idx,
  c(
    "celltype",
    "sample",
    "within_celltype_pct"
  )
]

dominance$sample <- as.character(
  dominance$sample
)

names(dominance)[3] <-
  "max_sample_contribution_pct"

dominance <- dominance[
  order(
    dominance$max_sample_contribution_pct,
    decreasing = TRUE
  ),
]

write.csv(
  dominance,
  file.path(
    outdir,
    "10_max_sample_contribution_by_celltype.csv"
  ),
  row.names = FALSE
)

dominant_types <- dominance[
  dominance$max_sample_contribution_pct >= 70,
]

if (nrow(dominant_types) > 0) {

  flag_lines <- c(
    flag_lines,
    "FLAG: cell types with >=70% contribution from one sample:",
    paste0(
      "  ",
      dominant_types$celltype,
      " : ",
      dominant_types$sample,
      " = ",
      round(
        dominant_types$max_sample_contribution_pct,
        1
      ),
      "%"
    ),
    ""
  )
}

flag_lines <- c(
  flag_lines,
  "NOTE:",
  "These are annotation/QC flags, not proof that a biological lineage is absent.",
  "Marker DotPlot and FeaturePlots must be reviewed before reannotation.",
  "",
  "IL10-axis output is exploratory only at v1.3.1.",
  "No modified RDS was written."
)

writeLines(
  flag_lines,
  file.path(
    outdir,
    "11_annotation_validation_flags.txt"
  )
)

cat(
  paste(
    flag_lines,
    collapse = "\n"
  )
)

cat("\n\n")

# ------------------------------------------------------------
# 19. Session info
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(
    outdir,
    "99_sessionInfo.txt"
  )
)

cat("============================================================\n")
cat("v1.3.1 COMPLETE\n")
cat("NO MODIFIED RDS WRITTEN\n")
cat("Output:\n")
cat(outdir, "\n")
cat("============================================================\n")

