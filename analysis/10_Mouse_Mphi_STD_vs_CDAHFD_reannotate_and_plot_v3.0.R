############################################################
## 10_Mouse_Mphi_STD_vs_CDAHFD_reannotate_and_plot_v3.0.R
##
## Full-mouse RDS -> macrophage extraction -> Layer2 reannotation
## -> STD/CDAHFD plots and summary tables
##
## Strategy
##   1. Extract Layer1 Kupffer_Macrophage + Monocyte cells
##   2. Reprocess this subset with RNA assay
##   3. Cluster macrophages de novo
##   4. Score five biologically defined marker modules
##   5. Assign one Layer2 label per cluster
##   6. Export integrated, condition-split, and sample-split UMAPs
##   7. Export module scores, subtype counts, M2/M1 ratio,
##      and SPP1/TREM2 percentage
##
## Palette
##   Ueno MASH Macro Palette v1.0
############################################################


############################################################
## 0. USER SETTINGS
############################################################

rds_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

output_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_Mphi_STD_vs_CDAHFD_v3.0"
)

## Metadata columns confirmed by diagnosis
layer1_col   <- "celltype_for_R8plot_FIXED2"
condition_col <- "condition"
sample_col    <- "sample"

## Layer1 values to include
mphi_layer1_values <- c(
  "Kupffer_Macrophage",
  "Monocyte"
)

## Conditions to analyze
conditions_keep <- c(
  "STD",
  "CDAHFD"
)

## Re-clustering settings
n_variable_features <- 3000L
n_pcs_use <- 30L
cluster_resolution <- 0.8
umap_min_dist <- 0.25
umap_n_neighbors <- 30L

## Plot settings
umap_point_size_integrated <- 0.85
umap_point_size_split <- 1.15
png_dpi <- 400
save_png <- TRUE

## Cluster-label confidence
## Top module score must exceed the second score by this margin.
minimum_score_margin <- 0.025

## If all five centered cluster scores are below this value,
## the cluster is classified as Other.
minimum_top_centered_score <- 0

random_seed <- 1234L
set.seed(random_seed)


############################################################
## 1. PACKAGES
############################################################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "tidyr",
  "data.table",
  "patchwork",
  "openxlsx",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing packages: ",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(patchwork)
  library(openxlsx)
  library(scales)
})


############################################################
## 2. OUTPUT DIRECTORIES
############################################################

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dir_rds    <- file.path(output_dir, "00_RDS")
dir_umap   <- file.path(output_dir, "01_UMAP")
dir_module <- file.path(output_dir, "02_ModuleScore")
dir_count  <- file.path(output_dir, "03_CellCounts")
dir_ratio  <- file.path(output_dir, "04_Ratios")
dir_table  <- file.path(output_dir, "05_Tables")
dir_marker <- file.path(output_dir, "06_ClusterMarkers")
dir_log    <- file.path(output_dir, "Logs")

for (x in c(
  dir_rds, dir_umap, dir_module, dir_count,
  dir_ratio, dir_table, dir_marker, dir_log
)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}


############################################################
## 3. LOGGING
############################################################

write_log <- function(text) {
  line <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    text
  )
  cat(line, "\n")
  cat(
    line,
    "\n",
    file = file.path(dir_log, "analysis.log"),
    append = TRUE
  )
}


############################################################
## 4. PALETTES
############################################################

palette_name <- "Ueno MASH Macro Palette v1.0"

ueno_mash_macro_palette_v1 <- c(
  "Resident Kupffer-like Mphi"      = "#00C853",
  "Monocyte-like Mphi"              = "#00A83B",
  "Inflammatory M1-like Mphi"       = "#FF1010",
  "Pro-resolution M2-like Mphi"     = "#F45BC1",
  "SPP1/TREM2 MASH-associated Mphi" = "#1261F3",
  "Other"                           = "#303030"
)

mphi_subtype_levels <- names(ueno_mash_macro_palette_v1)

condition_palette <- c(
  "STD" = "#3366F5",
  "CDAHFD" = "#FF4A4A"
)

m1_m2_palette <- c(
  "M2 cells" = "#3366F5",
  "M1 cells" = "#FF4A4A"
)


############################################################
## 5. MARKER SIGNATURES
############################################################

marker_sets <- list(

  Resident_Kupffer = c(
    "Clec4f", "Timd4", "Vsig4", "Marco",
    "Cd163", "Lyve1", "C1qa", "C1qb", "C1qc"
  ),

  Monocyte = c(
    "Ly6c2", "Ccr2", "Plac8", "S100a8",
    "S100a9", "Lyz2", "Ctss", "Fcgr3"
  ),

  Inflammatory_M1 = c(
    "Il1b", "Tnf", "Il6", "Il12b",
    "Il23a", "Nos2", "Ptgs2", "Cd80",
    "Cd86", "Cxcl9", "Cxcl10", "Ccl2"
  ),

  Pro_resolution_M2 = c(
    "Il10", "Mrc1", "Arg1", "Retnla",
    "Chil3", "Chil4", "Ccl17", "Ccl22",
    "Ccl24", "Maf"
  ),

  SPP1_TREM2_MASH = c(
    "Spp1", "Trem2", "Gpnmb", "Lgals3",
    "Fabp5", "Lpl", "Cd9", "Ctsb",
    "Ctsk", "Mmp12"
  )
)

label_lookup <- c(
  "Resident_Kupffer"   = "Resident Kupffer-like Mphi",
  "Monocyte"           = "Monocyte-like Mphi",
  "Inflammatory_M1"    = "Inflammatory M1-like Mphi",
  "Pro_resolution_M2"  = "Pro-resolution M2-like Mphi",
  "SPP1_TREM2_MASH"    = "SPP1/TREM2 MASH-associated Mphi"
)


############################################################
## 6. HELPER FUNCTIONS
############################################################

normalize_condition <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    grepl("CDAHFD|CDHFD|CDAH", x, ignore.case = TRUE) ~ "CDAHFD",
    grepl("^STD$|STD_|standard", x, ignore.case = TRUE) ~ "STD",
    TRUE ~ x
  )
}


theme_ueno_umap <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 2,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = base_size - 1,
        hjust = 0.5
      ),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size - 2),
      strip.background = element_rect(
        fill = "#F4F4F4",
        colour = "#888888",
        linewidth = 0.5
      ),
      strip.text = element_text(
        face = "bold",
        size = base_size
      ),
      panel.border = element_rect(
        colour = "#555555",
        fill = NA,
        linewidth = 0.5
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
}


theme_ueno_bar <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = base_size + 2
      ),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}


save_plot <- function(plot_object, filename, width, height) {

  ggsave(
    filename = paste0(filename, ".pdf"),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )

  if (isTRUE(save_png)) {
    ggsave(
      filename = paste0(filename, ".png"),
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      dpi = png_dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
}


prepare_signature_sets <- function(object, signatures) {

  available_features <- rownames(object)

  present <- lapply(
    signatures,
    function(g) intersect(g, available_features)
  )

  missing <- lapply(
    signatures,
    function(g) setdiff(g, available_features)
  )

  gene_status <- bind_rows(
    lapply(
      names(signatures),
      function(nm) {
        data.frame(
          signature = nm,
          gene = signatures[[nm]],
          present = signatures[[nm]] %in% available_features,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  fwrite(
    gene_status,
    file.path(dir_table, "Marker_gene_availability.csv")
  )

  too_short <- names(present)[
    lengths(present) < 3
  ]

  if (length(too_short) > 0) {
    stop(
      paste0(
        "Fewer than three usable genes in signature(s): ",
        paste(too_short, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  list(
    present = present,
    missing = missing
  )
}


make_umap_plot <- function(
    data,
    title,
    subtitle = NULL,
    point_size = 0.85,
    facet_col = NULL,
    show_legend = TRUE
) {

  p <- ggplot(
    data,
    aes(
      x = MphiUMAP_1,
      y = MphiUMAP_2,
      colour = Layer2_annotation
    )
  ) +
    geom_point(
      size = point_size,
      alpha = 1,
      stroke = 0,
      shape = 16
    ) +
    scale_colour_manual(
      values = ueno_mash_macro_palette_v1,
      limits = mphi_subtype_levels,
      drop = FALSE,
      name = "Mphi subtype"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "MphiUMAP_1",
      y = "MphiUMAP_2"
    ) +
    coord_fixed() +
    theme_ueno_umap(base_size = 12) +
    guides(
      colour = guide_legend(
        override.aes = list(
          size = 4.5,
          alpha = 1
        )
      )
    )

  if (!is.null(facet_col)) {
    p <- p +
      facet_wrap(
        stats::as.formula(paste0("~", facet_col)),
        nrow = 1
      )
  }

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  p
}


############################################################
## 7. LOAD FULL RDS
############################################################

write_log(paste0("Palette: ", palette_name))

if (!file.exists(rds_file)) {
  stop(
    paste0("RDS not found:\n", rds_file),
    call. = FALSE
  )
}

write_log(paste0("Loading RDS: ", rds_file))
obj <- readRDS(rds_file)

if (!inherits(obj, "Seurat")) {
  stop("Loaded object is not a Seurat object.", call. = FALSE)
}

required_metadata <- c(
  layer1_col,
  condition_col,
  sample_col
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(obj@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    paste0(
      "Missing metadata columns: ",
      paste(missing_metadata, collapse = ", ")
    ),
    call. = FALSE
  )
}

write_log(
  paste0(
    "Full object: ",
    format(ncol(obj), big.mark = ","),
    " cells; ",
    format(nrow(obj), big.mark = ","),
    " genes"
  )
)


############################################################
## 8. EXTRACT MACROPHAGES
############################################################

obj$Condition_v3 <- normalize_condition(
  obj@meta.data[[condition_col]]
)

layer1_values <- as.character(
  obj@meta.data[[layer1_col]]
)

cells_keep <- colnames(obj)[
  layer1_values %in% mphi_layer1_values &
    obj$Condition_v3 %in% conditions_keep
]

if (length(cells_keep) == 0) {
  stop(
    "No Kupffer_Macrophage/Monocyte cells found for STD/CDAHFD.",
    call. = FALSE
  )
}

mphi <- subset(
  obj,
  cells = cells_keep
)

mphi$Condition_v3 <- factor(
  normalize_condition(
    mphi@meta.data[[condition_col]]
  ),
  levels = conditions_keep
)

mphi$Sample_v3 <- as.character(
  mphi@meta.data[[sample_col]]
)

write_log(
  paste0(
    "Macrophages retained: ",
    format(ncol(mphi), big.mark = ",")
  )
)

fwrite(
  as.data.frame(
    table(
      Condition = mphi$Condition_v3,
      Layer1 = mphi@meta.data[[layer1_col]]
    )
  ),
  file.path(
    dir_table,
    "Extracted_cell_counts_by_condition_and_Layer1.csv"
  )
)


############################################################
## 9. RNA REPROCESSING
############################################################

if (!"RNA" %in% Assays(mphi)) {
  stop("RNA assay is absent.", call. = FALSE)
}

DefaultAssay(mphi) <- "RNA"

## Join Seurat v5 layers when necessary
mphi <- tryCatch(
  JoinLayers(mphi, assay = "RNA"),
  error = function(e) {
    write_log(
      paste0(
        "JoinLayers skipped: ",
        conditionMessage(e)
      )
    )
    mphi
  }
)

write_log("NormalizeData")
mphi <- NormalizeData(
  mphi,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

write_log("FindVariableFeatures")
mphi <- FindVariableFeatures(
  mphi,
  selection.method = "vst",
  nfeatures = n_variable_features,
  verbose = FALSE
)

write_log("ScaleData")
mphi <- ScaleData(
  mphi,
  features = VariableFeatures(mphi),
  verbose = FALSE
)

write_log("RunPCA")
mphi <- RunPCA(
  mphi,
  features = VariableFeatures(mphi),
  npcs = max(n_pcs_use, 30L),
  reduction.name = "mphi_pca",
  verbose = FALSE
)

write_log("FindNeighbors")
mphi <- FindNeighbors(
  mphi,
  reduction = "mphi_pca",
  dims = seq_len(n_pcs_use),
  graph.name = c(
    "mphi_nn",
    "mphi_snn"
  ),
  verbose = FALSE
)

write_log("FindClusters")
mphi <- FindClusters(
  mphi,
  graph.name = "mphi_snn",
  resolution = cluster_resolution,
  cluster.name = "mphi_cluster_v3",
  random.seed = random_seed,
  verbose = FALSE
)

write_log("RunUMAP")
mphi <- RunUMAP(
  mphi,
  reduction = "mphi_pca",
  dims = seq_len(n_pcs_use),
  reduction.name = "MphiUMAP",
  reduction.key = "MphiUMAP_",
  n.neighbors = umap_n_neighbors,
  min.dist = umap_min_dist,
  seed.use = random_seed,
  verbose = FALSE
)


############################################################
## 10. MODULE SCORES
############################################################

signature_info <- prepare_signature_sets(
  mphi,
  marker_sets
)

present_sets <- signature_info$present

for (signature_name in names(present_sets)) {

  score_prefix <- paste0(
    "Score_",
    signature_name
  )

  write_log(
    paste0(
      "AddModuleScore: ",
      signature_name,
      " (",
      length(present_sets[[signature_name]]),
      " genes)"
    )
  )

  mphi <- AddModuleScore(
    object = mphi,
    features = list(
      present_sets[[signature_name]]
    ),
    name = score_prefix,
    assay = "RNA",
    seed = random_seed
  )

  generated_col <- paste0(
    score_prefix,
    "1"
  )

  final_col <- paste0(
    "Module_",
    signature_name
  )

  mphi[[final_col]] <- mphi@meta.data[[generated_col]]
}

module_columns <- paste0(
  "Module_",
  names(marker_sets)
)


############################################################
## 11. CLUSTER-LEVEL LABEL ASSIGNMENT
############################################################

cluster_score_table <- mphi@meta.data %>%
  mutate(
    mphi_cluster_v3 = as.character(
      mphi_cluster_v3
    )
  ) %>%
  group_by(
    mphi_cluster_v3
  ) %>%
  summarise(
    n_cells = n(),
    across(
      all_of(module_columns),
      median,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

## Center each module across clusters so signatures are comparable.
for (col_now in module_columns) {
  centered_name <- paste0(
    col_now,
    "_centered"
  )

  cluster_score_table[[centered_name]] <-
    cluster_score_table[[col_now]] -
    median(
      cluster_score_table[[col_now]],
      na.rm = TRUE
    )
}

centered_columns <- paste0(
  module_columns,
  "_centered"
)

score_matrix <- as.matrix(
  cluster_score_table[
    ,
    centered_columns,
    drop = FALSE
  ]
)

top_index <- apply(
  score_matrix,
  1,
  which.max
)

top_score <- apply(
  score_matrix,
  1,
  max,
  na.rm = TRUE
)

second_score <- apply(
  score_matrix,
  1,
  function(x) {
    x_sorted <- sort(
      x,
      decreasing = TRUE,
      na.last = TRUE
    )
    if (length(x_sorted) < 2) {
      return(NA_real_)
    }
    x_sorted[2]
  }
)

top_signature <- names(marker_sets)[
  top_index
]

score_margin <- top_score - second_score

cluster_label <- unname(
  label_lookup[
    top_signature
  ]
)

cluster_label[
  is.na(top_score) |
    top_score < minimum_top_centered_score |
    is.na(score_margin) |
    score_margin < minimum_score_margin
] <- "Other"

cluster_score_table$top_signature <- top_signature
cluster_score_table$top_centered_score <- top_score
cluster_score_table$second_centered_score <- second_score
cluster_score_table$score_margin <- score_margin
cluster_score_table$Layer2_annotation <- cluster_label

cluster_score_table$Layer2_annotation <- factor(
  cluster_score_table$Layer2_annotation,
  levels = mphi_subtype_levels
)

fwrite(
  cluster_score_table,
  file.path(
    dir_table,
    "Cluster_module_scores_and_Layer2_assignment.csv"
  )
)

cluster_to_label <- setNames(
  as.character(
    cluster_score_table$Layer2_annotation
  ),
  cluster_score_table$mphi_cluster_v3
)

mphi$Layer2_annotation <- unname(
  cluster_to_label[
    as.character(
      mphi$mphi_cluster_v3
    )
  ]
)

mphi$Layer2_annotation <- factor(
  mphi$Layer2_annotation,
  levels = mphi_subtype_levels
)

write_log("Layer2 assignment completed")


############################################################
## 12. CLUSTER MARKERS FOR VALIDATION
############################################################

Idents(mphi) <- "mphi_cluster_v3"

write_log("FindAllMarkers for cluster validation")

cluster_markers <- FindAllMarkers(
  object = mphi,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

fwrite(
  cluster_markers,
  file.path(
    dir_marker,
    "Mphi_cluster_markers_all.csv"
  )
)

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(
    desc(avg_log2FC),
    .by_group = TRUE
  ) %>%
  slice_head(n = 30) %>%
  ungroup()

fwrite(
  top_markers,
  file.path(
    dir_marker,
    "Mphi_cluster_markers_top30.csv"
  )
)


############################################################
## 13. UMAP DATA
############################################################

emb <- Embeddings(
  mphi,
  reduction = "MphiUMAP"
)

umap_df <- data.frame(
  cell = rownames(emb),
  MphiUMAP_1 = emb[, 1],
  MphiUMAP_2 = emb[, 2],
  stringsAsFactors = FALSE
)

md <- mphi@meta.data
md$cell <- rownames(md)

umap_df <- left_join(
  umap_df,
  md,
  by = "cell"
)

umap_df$Condition_v3 <- factor(
  umap_df$Condition_v3,
  levels = conditions_keep
)

umap_df$Layer2_annotation <- factor(
  umap_df$Layer2_annotation,
  levels = mphi_subtype_levels
)

fwrite(
  umap_df,
  file.path(
    dir_table,
    "MphiUMAP_coordinates_and_metadata.csv"
  )
)


############################################################
## 14. INTEGRATED UMAP
############################################################

p_umap_integrated <- make_umap_plot(
  data = umap_df,
  title = "Mouse macrophage UMAP: integrated",
  subtitle = paste0(
    "STD and CDAHFD | ",
    palette_name
  ),
  point_size = umap_point_size_integrated,
  show_legend = TRUE
)

save_plot(
  p_umap_integrated,
  file.path(
    dir_umap,
    "UMAP_Mphi_integrated_v3"
  ),
  width = 10,
  height = 7.5
)


############################################################
## 15. CONDITION-SPLIT UMAP
############################################################

p_umap_condition <- make_umap_plot(
  data = umap_df,
  title = "Mouse macrophage UMAP by condition",
  subtitle = "Shared integrated coordinates",
  point_size = umap_point_size_split,
  facet_col = "Condition_v3",
  show_legend = TRUE
)

save_plot(
  p_umap_condition,
  file.path(
    dir_umap,
    "UMAP_Mphi_STD_vs_CDAHFD_shared_coordinates_v3"
  ),
  width = 14,
  height = 7
)


############################################################
## 16. INDIVIDUAL CONDITION UMAPS
############################################################

for (condition_now in conditions_keep) {

  df_now <- umap_df %>%
    filter(
      Condition_v3 == condition_now
    )

  p_now <- make_umap_plot(
    data = df_now,
    title = paste0(
      "Mouse macrophage UMAP: ",
      condition_now
    ),
    subtitle = "Shared integrated coordinates",
    point_size = umap_point_size_split,
    show_legend = TRUE
  )

  save_plot(
    p_now,
    file.path(
      dir_umap,
      paste0(
        "UMAP_Mphi_",
        condition_now,
        "_v3"
      )
    ),
    width = 9,
    height = 7
  )
}


############################################################
## 17. SAMPLE-SPLIT UMAP
############################################################

sample_levels <- unique(
  as.character(
    umap_df$Sample_v3
  )
)

umap_df$Sample_v3 <- factor(
  umap_df$Sample_v3,
  levels = sample_levels
)

n_samples <- length(sample_levels)
n_cols <- min(3L, n_samples)
n_rows <- ceiling(n_samples / n_cols)

p_umap_sample <- ggplot(
  umap_df,
  aes(
    x = MphiUMAP_1,
    y = MphiUMAP_2,
    colour = Layer2_annotation
  )
) +
  geom_point(
    size = umap_point_size_split,
    alpha = 1,
    stroke = 0
  ) +
  scale_colour_manual(
    values = ueno_mash_macro_palette_v1,
    limits = mphi_subtype_levels,
    drop = FALSE,
    name = "Mphi subtype"
  ) +
  facet_wrap(
    ~Sample_v3,
    ncol = n_cols
  ) +
  coord_fixed() +
  labs(
    title = "Mouse macrophage UMAP by sample",
    subtitle = "Shared integrated coordinates",
    x = "MphiUMAP_1",
    y = "MphiUMAP_2"
  ) +
  theme_ueno_umap(base_size = 11) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 4.5
      )
    )
  )

save_plot(
  p_umap_sample,
  file.path(
    dir_umap,
    "UMAP_Mphi_by_sample_v3"
  ),
  width = 14,
  height = max(6, n_rows * 4.6)
)


############################################################
## 18. CLUSTER-ID UMAP FOR VALIDATION
############################################################

p_cluster <- ggplot(
  umap_df,
  aes(
    x = MphiUMAP_1,
    y = MphiUMAP_2,
    colour = factor(mphi_cluster_v3)
  )
) +
  geom_point(
    size = 0.75,
    alpha = 1,
    stroke = 0
  ) +
  coord_fixed() +
  labs(
    title = "Mouse macrophage UMAP: de novo clusters",
    x = "MphiUMAP_1",
    y = "MphiUMAP_2",
    colour = "Cluster"
  ) +
  theme_ueno_umap(base_size = 12) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 4)
    )
  )

save_plot(
  p_cluster,
  file.path(
    dir_umap,
    "UMAP_Mphi_clusters_v3"
  ),
  width = 10,
  height = 7.5
)


############################################################
## 19. MODULE-SCORE VIOLIN PLOTS
############################################################

module_df <- mphi@meta.data %>%
  select(
    Condition_v3,
    all_of(module_columns)
  ) %>%
  pivot_longer(
    cols = all_of(module_columns),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  mutate(
    Module = recode(
      Module,
      "Module_Resident_Kupffer" =
        "Resident Kupffer module",
      "Module_Monocyte" =
        "Monocyte module",
      "Module_Inflammatory_M1" =
        "Inflammatory M1 module",
      "Module_Pro_resolution_M2" =
        "Pro-resolution M2 module",
      "Module_SPP1_TREM2_MASH" =
        "SPP1/TREM2 MASH module"
    )
  )

fwrite(
  module_df,
  file.path(
    dir_table,
    "Module_scores_cell_level.csv"
  )
)

p_module <- ggplot(
  module_df,
  aes(
    x = Condition_v3,
    y = Score,
    fill = Condition_v3
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4,
    colour = "#333333"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = NA,
    colour = "#202020",
    linewidth = 0.45
  ) +
  facet_wrap(
    ~Module,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = condition_palette,
    drop = FALSE
  ) +
  labs(
    title = "Macrophage marker modules",
    subtitle = "STD versus CDAHFD",
    x = NULL,
    y = "AddModuleScore"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(
      fill = "#F4F4F4",
      colour = "#777777"
    ),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black"
    )
  )

save_plot(
  p_module,
  file.path(
    dir_module,
    "ModuleScore_five_signatures_STD_vs_CDAHFD_v3"
  ),
  width = 12,
  height = 8
)


############################################################
## 20. SUBTYPE COUNTS AND FRACTIONS
############################################################

subtype_summary <- umap_df %>%
  count(
    Condition_v3,
    Layer2_annotation,
    name = "Cell_count"
  ) %>%
  group_by(
    Condition_v3
  ) %>%
  mutate(
    Total_Mphi = sum(Cell_count),
    Fraction = Cell_count / Total_Mphi,
    Percent = Fraction * 100
  ) %>%
  ungroup()

fwrite(
  subtype_summary,
  file.path(
    dir_table,
    "Mphi_subtype_counts_and_fractions.csv"
  )
)

p_subtype_fraction <- ggplot(
  subtype_summary,
  aes(
    x = Condition_v3,
    y = Fraction,
    fill = Layer2_annotation
  )
) +
  geom_col(
    width = 0.62
  ) +
  scale_fill_manual(
    values = ueno_mash_macro_palette_v1,
    limits = mphi_subtype_levels,
    drop = FALSE,
    name = "Mphi subtype"
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1)
  ) +
  labs(
    title = "Macrophage subtype composition",
    x = NULL,
    y = "Fraction of total macrophages"
  ) +
  theme_ueno_bar(base_size = 13)

save_plot(
  p_subtype_fraction,
  file.path(
    dir_count,
    "Mphi_subtype_fraction_STD_vs_CDAHFD_v3"
  ),
  width = 7.5,
  height = 7
)


############################################################
## 21. M1/M2 COUNTS AND RATIO
############################################################

m1_m2_counts <- umap_df %>%
  mutate(
    M1_M2_class = case_when(
      Layer2_annotation ==
        "Inflammatory M1-like Mphi" ~ "M1 cells",
      Layer2_annotation ==
        "Pro-resolution M2-like Mphi" ~ "M2 cells",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(M1_M2_class)
  ) %>%
  count(
    Condition_v3,
    M1_M2_class,
    name = "Cell_count"
  ) %>%
  complete(
    Condition_v3 = factor(
      conditions_keep,
      levels = conditions_keep
    ),
    M1_M2_class = c(
      "M2 cells",
      "M1 cells"
    ),
    fill = list(
      Cell_count = 0
    )
  )

m1_m2_counts$M1_M2_class <- factor(
  m1_m2_counts$M1_M2_class,
  levels = c(
    "M2 cells",
    "M1 cells"
  )
)

fwrite(
  m1_m2_counts,
  file.path(
    dir_table,
    "M1_M2_cell_counts.csv"
  )
)

p_m1_m2_count <- ggplot(
  m1_m2_counts,
  aes(
    x = Condition_v3,
    y = Cell_count,
    fill = M1_M2_class
  )
) +
  geom_col(width = 0.58) +
  scale_fill_manual(
    values = m1_m2_palette,
    drop = FALSE,
    name = "Cell type"
  ) +
  labs(
    title = "M1 and M2 macrophage cell counts",
    x = NULL,
    y = "Cell count"
  ) +
  theme_ueno_bar(base_size = 13)

save_plot(
  p_m1_m2_count,
  file.path(
    dir_count,
    "M1_M2_cell_counts_STD_vs_CDAHFD_v3"
  ),
  width = 7,
  height = 6.5
)

ratio_df <- m1_m2_counts %>%
  mutate(
    class_short = case_when(
      M1_M2_class == "M1 cells" ~ "M1",
      M1_M2_class == "M2 cells" ~ "M2"
    )
  ) %>%
  select(
    Condition_v3,
    class_short,
    Cell_count
  ) %>%
  pivot_wider(
    names_from = class_short,
    values_from = Cell_count,
    values_fill = 0
  ) %>%
  mutate(
    M2_M1_ratio = (M2 + 0.5) / (M1 + 0.5),
    Log2_M2_M1_ratio = log2(M2_M1_ratio)
  )

fwrite(
  ratio_df,
  file.path(
    dir_table,
    "M2_M1_ratio.csv"
  )
)

p_ratio <- ggplot(
  ratio_df,
  aes(
    x = Condition_v3,
    y = Log2_M2_M1_ratio,
    fill = Condition_v3
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "#888888",
    linewidth = 0.6
  ) +
  geom_col(
    width = 0.58,
    colour = "#333333",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Log2_M2_M1_ratio
      )
    ),
    vjust = ifelse(
      ratio_df$Log2_M2_M1_ratio >= 0,
      -0.4,
      1.4
    ),
    fontface = "bold",
    size = 4.2
  ) +
  scale_fill_manual(
    values = condition_palette
  ) +
  labs(
    title = "M2/M1 cell ratio",
    x = NULL,
    y = "Log2((M2 + 0.5)/(M1 + 0.5))"
  ) +
  theme_ueno_bar(base_size = 13) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_ratio,
  file.path(
    dir_ratio,
    "M2_M1_log2_ratio_STD_vs_CDAHFD_v3"
  ),
  width = 7,
  height = 6.5
)


############################################################
## 22. SPP1/TREM2 PERCENTAGE
############################################################

spp1_trem2_df <- umap_df %>%
  group_by(
    Condition_v3
  ) %>%
  summarise(
    Total_Mphi = n(),
    SPP1_TREM2_Mphi = sum(
      Layer2_annotation ==
        "SPP1/TREM2 MASH-associated Mphi",
      na.rm = TRUE
    ),
    Percent = 100 *
      SPP1_TREM2_Mphi /
      Total_Mphi,
    .groups = "drop"
  )

fwrite(
  spp1_trem2_df,
  file.path(
    dir_table,
    "SPP1_TREM2_Mphi_percentage.csv"
  )
)

p_spp1_trem2 <- ggplot(
  spp1_trem2_df,
  aes(
    x = Condition_v3,
    y = Percent,
    fill = Condition_v3
  )
) +
  geom_col(
    width = 0.58,
    colour = "#333333",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = paste0(
        sprintf("%.2f%%", Percent),
        "\n",
        SPP1_TREM2_Mphi,
        "/",
        Total_Mphi
      )
    ),
    vjust = -0.35,
    fontface = "bold",
    size = 4
  ) +
  scale_fill_manual(
    values = condition_palette
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.18)
    )
  ) +
  labs(
    title = "SPP1/TREM2 macrophages",
    x = NULL,
    y = "Percentage of total macrophages"
  ) +
  theme_ueno_bar(base_size = 13) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_spp1_trem2,
  file.path(
    dir_ratio,
    "SPP1_TREM2_Mphi_percentage_STD_vs_CDAHFD_v3"
  ),
  width = 7,
  height = 6.5
)


############################################################
## 23. DOTPLOT FOR VALIDATION
############################################################

dot_features <- unique(
  unlist(
    present_sets,
    use.names = FALSE
  )
)

Idents(mphi) <- "Layer2_annotation"

p_dot <- DotPlot(
  object = mphi,
  features = dot_features,
  assay = "RNA",
  group.by = "Layer2_annotation",
  dot.scale = 7
) +
  RotatedAxis() +
  scale_colour_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0
  ) +
  labs(
    title = "Marker validation by reannotated macrophage subtype",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    )
  )

save_plot(
  p_dot,
  file.path(
    dir_module,
    "DotPlot_marker_validation_v3"
  ),
  width = 16,
  height = 7
)


############################################################
## 24. COMBINED FIGURE
############################################################

combined_figure <-
  (
    p_umap_integrated |
      (
        p_umap_condition +
          theme(
            legend.position = "none"
          )
      )
  ) /
  (
    (
      p_module +
        theme(
          legend.position = "none"
        )
    ) |
      p_m1_m2_count |
      p_ratio |
      p_spp1_trem2
  ) +
  plot_layout(
    heights = c(1.25, 1)
  ) +
  plot_annotation(
    title = "Mouse MASH model liver macrophage analysis",
    subtitle = paste0(
      "STD versus CDAHFD | de novo Layer2 reannotation | ",
      palette_name
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 22,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 13,
        hjust = 0.5
      )
    )
  )

save_plot(
  combined_figure,
  file.path(
    output_dir,
    "Mouse_Mphi_STD_vs_CDAHFD_combined_figure_v3"
  ),
  width = 22,
  height = 13
)


############################################################
## 25. EXCEL SUMMARY
############################################################

settings_table <- data.frame(
  Item = c(
    "RDS file",
    "Layer1 column",
    "Condition column",
    "Sample column",
    "Layer1 values included",
    "Cells analyzed",
    "Variable features",
    "PCs",
    "Cluster resolution",
    "UMAP min.dist",
    "UMAP n.neighbors",
    "Minimum score margin",
    "Minimum top centered score",
    "Palette",
    "Analysis date"
  ),
  Value = c(
    rds_file,
    layer1_col,
    condition_col,
    sample_col,
    paste(
      mphi_layer1_values,
      collapse = ", "
    ),
    ncol(mphi),
    n_variable_features,
    n_pcs_use,
    cluster_resolution,
    umap_min_dist,
    umap_n_neighbors,
    minimum_score_margin,
    minimum_top_centered_score,
    palette_name,
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  stringsAsFactors = FALSE
)

palette_table <- data.frame(
  Palette_name = palette_name,
  Mphi_subtype = names(
    ueno_mash_macro_palette_v1
  ),
  HEX = unname(
    ueno_mash_macro_palette_v1
  ),
  stringsAsFactors = FALSE
)

marker_table <- bind_rows(
  lapply(
    names(marker_sets),
    function(nm) {
      data.frame(
        signature = nm,
        gene = marker_sets[[nm]],
        present = marker_sets[[nm]] %in%
          rownames(mphi),
        stringsAsFactors = FALSE
      )
    }
  )
)

wb <- createWorkbook()

addWorksheet(wb, "Settings")
writeData(wb, "Settings", settings_table)

addWorksheet(wb, "Palette")
writeData(wb, "Palette", palette_table)

addWorksheet(wb, "Marker_sets")
writeData(wb, "Marker_sets", marker_table)

addWorksheet(wb, "Cluster_assignment")
writeData(
  wb,
  "Cluster_assignment",
  cluster_score_table
)

addWorksheet(wb, "Subtype_summary")
writeData(
  wb,
  "Subtype_summary",
  subtype_summary
)

addWorksheet(wb, "M1_M2_counts")
writeData(
  wb,
  "M1_M2_counts",
  m1_m2_counts
)

addWorksheet(wb, "M2_M1_ratio")
writeData(
  wb,
  "M2_M1_ratio",
  ratio_df
)

addWorksheet(wb, "SPP1_TREM2")
writeData(
  wb,
  "SPP1_TREM2",
  spp1_trem2_df
)

saveWorkbook(
  wb,
  file.path(
    dir_table,
    "Mouse_Mphi_STD_vs_CDAHFD_summary_v3.xlsx"
  ),
  overwrite = TRUE
)


############################################################
## 26. SAVE REANNOTATED OBJECT AND PALETTE
############################################################

saveRDS(
  mphi,
  file.path(
    dir_rds,
    "Mouse_Mphi_STD_CDAHFD_reannotated_v3.rds"
  ),
  compress = FALSE
)

saveRDS(
  ueno_mash_macro_palette_v1,
  file.path(
    dir_table,
    "Ueno_MASH_Macro_Palette_v1.rds"
  )
)

palette_script <- c(
  '# Ueno MASH Macro Palette v1.0',
  '',
  'ueno_mash_macro_palette_v1 <- c(',
  '  "Resident Kupffer-like Mphi" = "#00C853",',
  '  "Monocyte-like Mphi" = "#00A83B",',
  '  "Inflammatory M1-like Mphi" = "#FF1010",',
  '  "Pro-resolution M2-like Mphi" = "#F45BC1",',
  '  "SPP1/TREM2 MASH-associated Mphi" = "#1261F3",',
  '  "Other" = "#303030"',
  ')'
)

writeLines(
  palette_script,
  file.path(
    dir_table,
    "Ueno_MASH_Macro_Palette_v1.R"
  )
)


############################################################
## 27. SESSION INFORMATION
############################################################

sink(
  file.path(
    dir_log,
    "sessionInfo.txt"
  )
)

cat("Analysis completed:\n")
cat(
  format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  ),
  "\n\n"
)

cat("RDS:\n", rds_file, "\n\n", sep = "")
cat("Cells analyzed:\n", ncol(mphi), "\n\n", sep = "")
cat("Palette:\n", palette_name, "\n\n", sep = "")
print(sessionInfo())

sink()


############################################################
## 28. FINISH
############################################################

write_log("Analysis completed successfully")

cat("\n")
cat("====================================================\n")
cat("Mouse Mphi reannotation and plotting v3.0 completed\n")
cat("====================================================\n")
cat("Cells analyzed :", ncol(mphi), "\n")
cat("Clusters       :", length(unique(mphi$mphi_cluster_v3)), "\n")
cat("Output         :", output_dir, "\n")
cat("Palette        :", palette_name, "\n")
cat("====================================================\n")
