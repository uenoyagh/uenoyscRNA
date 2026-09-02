suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

VERSION <- "v6.7.0"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "RDS3_annotation_visualization_v4.1.1/",
  "objects/",
  "RDS3_with_visualization_metadata_v4.1.1.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_MASH_LSEC_",
  VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("Mouse MASH LSEC endothelial parent audit\n")
cat("Version:", VERSION, "\n")
cat("============================================\n\n")

cat("Reading:\n", INPUT_RDS, "\n\n")
obj <- readRDS(INPUT_RDS)

anno_col <- "celltype_for_R8plot_FIXED2"
sample_col <- "sample"
condition_col <- "condition"

required_cols <- c(anno_col, sample_col, condition_col)

missing_cols <- setdiff(required_cols, colnames(obj@meta.data))

if (length(missing_cols) > 0) {
  stop(
    "Required metadata column(s) missing: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("Whole object cells:", ncol(obj), "\n")
cat("Whole object features:", nrow(obj), "\n")
cat("Annotation column:", anno_col, "\n\n")

anno <- as.character(obj@meta.data[[anno_col]])

cat("=== Endothelial-like labels in annotation ===\n")

endothelial_labels_found <- sort(
  unique(
    anno[
      grepl(
        "LSEC|Vascular|endothelial",
        anno,
        ignore.case = TRUE
      )
    ]
  )
)

print(endothelial_labels_found)

cat("\n=== Counts of endothelial-like labels ===\n")
print(
  sort(
    table(
      anno[
        grepl(
          "LSEC|Vascular|endothelial",
          anno,
          ignore.case = TRUE
        )
      ]
    ),
    decreasing = TRUE
  )
)

target_labels <- c(
  "LSEC",
  "Vascular_endothelial"
)

missing_target <- setdiff(
  target_labels,
  unique(anno)
)

if (length(missing_target) > 0) {
  stop(
    "Expected endothelial label(s) not found: ",
    paste(missing_target, collapse = ", ")
  )
}

keep_cells <- colnames(obj)[anno %in% target_labels]

cat("\n=== Endothelial parent ===\n")
cat("Parent cells:", length(keep_cells), "\n")

parent <- subset(
  obj,
  cells = keep_cells
)

parent$endothelial_parent_label_v670 <-
  as.character(parent@meta.data[[anno_col]])

cat("\n=== Parent label counts ===\n")
print(
  table(
    parent$endothelial_parent_label_v670,
    useNA = "ifany"
  )
)

cat("\n=== Parent by sample ===\n")
print(
  table(
    parent@meta.data[[sample_col]],
    parent$endothelial_parent_label_v670,
    useNA = "ifany"
  )
)

cat("\n=== Parent by condition ===\n")
print(
  table(
    parent@meta.data[[condition_col]],
    parent$endothelial_parent_label_v670,
    useNA = "ifany"
  )
)

sample_levels <- unique(
  as.character(obj@meta.data[[sample_col]])
)

whole_sample_n <- table(
  factor(
    as.character(obj@meta.data[[sample_col]]),
    levels = sample_levels
  )
)

parent_sample_n <- table(
  factor(
    as.character(parent@meta.data[[sample_col]]),
    levels = sample_levels
  )
)

sample_summary <- data.frame(
  sample = sample_levels,
  whole_cells = as.integer(whole_sample_n),
  endothelial_parent_cells = as.integer(parent_sample_n),
  stringsAsFactors = FALSE
)

sample_summary$endothelial_fraction <-
  sample_summary$endothelial_parent_cells /
  sample_summary$whole_cells

for (lab in target_labels) {

  n_lab <- table(
    factor(
      as.character(
        parent@meta.data[[sample_col]][
          parent$endothelial_parent_label_v670 == lab
        ]
      ),
      levels = sample_levels
    )
  )

  new_col_n <- paste0(
    gsub("[^A-Za-z0-9]+", "_", lab),
    "_cells"
  )

  new_col_fraction <- paste0(
    gsub("[^A-Za-z0-9]+", "_", lab),
    "_fraction_of_whole"
  )

  sample_summary[[new_col_n]] <- as.integer(n_lab)

  sample_summary[[new_col_fraction]] <-
    sample_summary[[new_col_n]] /
    sample_summary$whole_cells
}

write.csv(
  sample_summary,
  file.path(
    TABDIR,
    "Endothelial_parent_counts_and_fraction_by_sample_v6.7.0.csv"
  ),
  row.names = FALSE
)

condition_table <- as.data.frame(
  table(
    condition = as.character(
      parent@meta.data[[condition_col]]
    ),
    endothelial_class =
      parent$endothelial_parent_label_v670
  )
)

write.csv(
  condition_table,
  file.path(
    TABDIR,
    "Endothelial_parent_counts_by_condition_v6.7.0.csv"
  ),
  row.names = FALSE
)

cross_cols <- intersect(
  c(
    "celltype_auto_annotation",
    "celltype_for_R8plot_FIXED2",
    "vote_ueno_celltype_v34",
    "vote_ueno_summary_v40",
    "target_cellclass_v410"
  ),
  colnames(parent@meta.data)
)

cross_df <- parent@meta.data[, cross_cols, drop = FALSE]
cross_df$cell <- rownames(cross_df)

write.csv(
  cross_df,
  file.path(
    TABDIR,
    "Endothelial_parent_annotation_audit_v6.7.0.csv"
  ),
  row.names = FALSE
)

if (!"umapRPCA" %in% Reductions(obj)) {
  stop("Reduction 'umapRPCA' not found.")
}

um <- Embeddings(obj, "umapRPCA")

plot_df <- data.frame(
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  annotation = anno,
  stringsAsFactors = FALSE
)

plot_df$highlight <- "Other"
plot_df$highlight[
  plot_df$annotation == "LSEC"
] <- "LSEC"

plot_df$highlight[
  plot_df$annotation == "Vascular_endothelial"
] <- "Vascular_endothelial"

plot_df$highlight <- factor(
  plot_df$highlight,
  levels = c(
    "Other",
    "LSEC",
    "Vascular_endothelial"
  )
)

plot_df <- plot_df[
  order(plot_df$highlight != "Other"),
]

p_whole <- ggplot(
  plot_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = highlight
  )
) +
  geom_point(
    size = 0.20,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "Other" = "grey88",
      "LSEC" = "#00BFC4",
      "Vascular_endothelial" = "#0066FF"
    )
  ) +
  coord_equal() +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle(
    "Whole liver: endothelial parent"
  )

ggsave(
  file.path(
    FIGDIR,
    "Wholecell_UMAP_LSEC_Vascular_endothelial_highlight_v6.7.0.pdf"
  ),
  p_whole,
  width = 8.5,
  height = 7.0
)

parent_um <- Embeddings(
  parent,
  "umapRPCA"
)

parent_plot_df <- data.frame(
  UMAP_1 = parent_um[, 1],
  UMAP_2 = parent_um[, 2],
  endothelial_class =
    parent$endothelial_parent_label_v670,
  sample =
    as.character(
      parent@meta.data[[sample_col]]
    ),
  condition =
    as.character(
      parent@meta.data[[condition_col]]
    )
)

p_parent <- ggplot(
  parent_plot_df,
  aes(
    UMAP_1,
    UMAP_2,
    color = endothelial_class
  )
) +
  geom_point(
    size = 0.45,
    alpha = 0.90
  ) +
  scale_color_manual(
    values = c(
      "LSEC" = "#00BFC4",
      "Vascular_endothelial" = "#0066FF"
    )
  ) +
  coord_equal() +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle(
    "Endothelial parent on original RPCA UMAP"
  )

ggsave(
  file.path(
    FIGDIR,
    "Endothelial_parent_original_UMAP_v6.7.0.pdf"
  ),
  p_parent,
  width = 8.0,
  height = 7.0
)

marker_panel <- c(
  "Pecam1",
  "Cdh5",
  "Kdr",
  "Klf2",
  "Clec4g",
  "Stab1",
  "Stab2",
  "Lyve1",
  "Fcgr2b",
  "Mrc1",
  "Plvap",
  "Vwf",
  "Emcn",
  "Esm1",
  "Cd34",
  "Rgcc",
  "Car4",
  "Ackr1",
  "Ptprc",
  "Adgre1",
  "Lyz2",
  "Col1a1",
  "Col3a1",
  "Rgs5",
  "Alb",
  "Krt19"
)

marker_present <- marker_panel[
  marker_panel %in% rownames(parent)
]

marker_missing <- setdiff(
  marker_panel,
  marker_present
)

cat("\n=== Marker audit ===\n")
cat(
  "Present:",
  paste(marker_present, collapse = ", "),
  "\n"
)

cat(
  "Missing:",
  paste(marker_missing, collapse = ", "),
  "\n"
)

writeLines(
  c(
    paste(
      "Present:",
      paste(marker_present, collapse = ", ")
    ),
    paste(
      "Missing:",
      paste(marker_missing, collapse = ", ")
    )
  ),
  file.path(
    TABDIR,
    "Marker_availability_v6.7.0.txt"
  )
)

DefaultAssay(parent) <- "RNA"

if (length(marker_present) > 0) {

  Idents(parent) <-
    parent$endothelial_parent_label_v670

  p_dot <- DotPlot(
    parent,
    features = marker_present,
    assay = "RNA",
    dot.scale = 6
  ) +
    RotatedAxis() +
    theme_classic(base_size = 11) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    ggtitle(
      "Endothelial parent marker audit"
    )

  ggsave(
    file.path(
      FIGDIR,
      "Endothelial_parent_marker_DotPlot_v6.7.0.pdf"
    ),
    p_dot,
    width = 14,
    height = 5.5
  )
}

saveRDS(
  parent,
  file.path(
    OBJDIR,
    "Mouse_MASH_endothelial_parent_v6.7.0.rds"
  ),
  compress = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTDIR,
    "sessionInfo_v6.7.0.txt"
  )
)

cat("\n============================================\n")
cat("v6.7.0 COMPLETE\n")
cat("Endothelial parent cells:", ncol(parent), "\n")
cat("Output:", OUTDIR, "\n")
cat("============================================\n")
