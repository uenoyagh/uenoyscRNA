############################################################
## 10_Mouse_Mphi_STD_vs_CDAHFD_Figures.R
##
## Mouse macrophage analysis:
## STD versus CDAHFD
##
## Outputs
##   1. Integrated UMAP
##   2. Condition-split UMAP
##   3. STD-only UMAP
##   4. CDAHFD-only UMAP
##   5. Sample-split UMAP
##   6. M1/M2 module-score violin plots
##   7. M1/M2 cell counts
##   8. M2/M1 ratio
##   9. Profibrotic macrophage percentage
##  10. CSV / Excel summary
##
## Palette:
##   Ueno MASH Macro Palette v1.0
############################################################


############################################################
## 0. USER SETTINGS
############################################################

## マクロファージ解析用RDSのパス
## 現在使用するRDSに合わせて、ここだけ修正してください。
rds_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

## 出力先
output_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_Mphi_STD_vs_CDAHFD_Figures"
)

## 対象条件
conditions_keep <- c("STD", "CDAHFD")

## UMAPドットサイズ
umap_point_size <- 2.0

## UMAPのPNG解像度
png_dpi <- 400

## TRUEの場合、各図をPNGでも保存
save_png <- TRUE

## TRUEの場合、sample別UMAPも作成
make_sample_umap <- TRUE

## TRUEの場合、module scoreを新規計算
## 既存のM1/M2 score列が存在する場合は自動的に再利用します。
recalculate_module_scores <- FALSE

## Module scoreの乱数固定
random_seed <- 1234


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
      "以下のパッケージがありません:\n",
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

set.seed(random_seed)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir_umap <- file.path(output_dir, "01_UMAP")
dir_module <- file.path(output_dir, "02_ModuleScore")
dir_count <- file.path(output_dir, "03_CellCounts")
dir_ratio <- file.path(output_dir, "04_Ratios")
dir_table <- file.path(output_dir, "05_Tables")
dir_log <- file.path(output_dir, "Logs")

for (x in c(
  dir_umap,
  dir_module,
  dir_count,
  dir_ratio,
  dir_table,
  dir_log
)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}


############################################################
## 2. PALETTE DEFINITIONS
############################################################

palette_name <- "Ueno MASH Macro Palette v1.0"

## 添付図左上の発色を再現する標準色
ueno_mash_macro_palette_v1 <- c(
  "Resident Kupffer-like Mphi" =
    "#00C853",

  "Monocyte-like Mphi" =
    "#00A83B",

  "Inflammatory M1-like Mphi" =
    "#FF1010",

  "Pro-resolution M2-like Mphi" =
    "#F45BC1",

  "SPP1/TREM2 MASH-associated Mphi" =
    "#1261F3",

  "Other" =
    "#303030"
)

condition_palette <- c(
  "STD" = "#3366F5",
  "CDAHFD" = "#FF4A4A"
)

m1_m2_palette <- c(
  "M1 cells" = "#FF4A4A",
  "M2 cells" = "#3366F5"
)

## 表示順
mphi_subtype_levels <- c(
  "Resident Kupffer-like Mphi",
  "Monocyte-like Mphi",
  "Inflammatory M1-like Mphi",
  "Pro-resolution M2-like Mphi",
  "SPP1/TREM2 MASH-associated Mphi",
  "Other"
)


############################################################
## 3. MARKER GENE SETS
############################################################

## マウスM1-like macrophage module
m1_markers_mouse <- c(
  "Il1b",
  "Tnf",
  "Il6",
  "Il12b",
  "Il23a",
  "Ccl2",
  "Ccl3",
  "Ccl4",
  "Cxcl9",
  "Cxcl10",
  "Nos2",
  "Ptgs2",
  "Cd80",
  "Cd86"
)

## マウスM2 / resolution-like macrophage module
m2_markers_mouse <- c(
  "Il10",
  "Mrc1",
  "Arg1",
  "Retnla",
  "Chil3",
  "Chil4",
  "Ccl17",
  "Ccl22",
  "Ccl24",
  "Maf"
)

## Profibrotic macrophage module
profibrotic_markers_mouse <- c(
  "Spp1",
  "Trem2",
  "Gpnmb",
  "Lgals3",
  "Tgfb1",
  "Pdgfb",
  "Fn1",
  "Mmp9",
  "Mmp12",
  "Ctsb",
  "Ctsk"
)


############################################################
## 4. HELPER FUNCTIONS
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


detect_metadata_column <- function(
    object,
    candidates,
    label,
    required = TRUE
) {

  found <- candidates[
    candidates %in% colnames(object@meta.data)
  ]

  if (length(found) == 0) {

    if (required) {
      stop(
        paste0(
          label,
          "列を検出できませんでした。\n",
          "候補: ",
          paste(candidates, collapse = ", "),
          "\n\nMetadata columns:\n",
          paste(
            colnames(object@meta.data),
            collapse = "\n"
          )
        ),
        call. = FALSE
      )
    }

    return(NA_character_)
  }

  found[1]
}


detect_reduction <- function(object) {

  reductions_available <- Reductions(object)

  preferred <- c(
    "MphiUMAP",
    "mphi_umap",
    "mphiUMAP",
    "umap",
    "UMAP",
    "rpca.umap",
    "integrated.umap"
  )

  found <- preferred[
    preferred %in% reductions_available
  ]

  if (length(found) > 0) {
    return(found[1])
  }

  umap_like <- reductions_available[
    grepl(
      "umap",
      reductions_available,
      ignore.case = TRUE
    )
  ]

  if (length(umap_like) > 0) {
    return(umap_like[1])
  }

  stop(
    paste0(
      "UMAP reductionが見つかりません。\n",
      "Available reductions: ",
      paste(reductions_available, collapse = ", ")
    ),
    call. = FALSE
  )
}


normalize_condition <- function(x) {

  x <- as.character(x)

  x <- case_when(
    grepl(
      "CDAHFD|CDHFD|CDAH",
      x,
      ignore.case = TRUE
    ) ~ "CDAHFD",

    grepl(
      "STD|Standard|Control",
      x,
      ignore.case = TRUE
    ) ~ "STD",

    TRUE ~ x
  )

  x
}


normalize_mphi_annotation <- function(x) {

  x_original <- as.character(x)
  x_clean <- trimws(x_original)

  result <- case_when(

    grepl(
      "Resident.*Kupffer|Kupffer.*Resident|Resident Kupffer",
      x_clean,
      ignore.case = TRUE
    ) ~ "Resident Kupffer-like Mphi",

    grepl(
      "Monocyte",
      x_clean,
      ignore.case = TRUE
    ) ~ "Monocyte-like Mphi",

    grepl(
      "Inflammatory|M1-like|M1 like",
      x_clean,
      ignore.case = TRUE
    ) ~ "Inflammatory M1-like Mphi",

    grepl(
      "Pro-resolution|Pro resolution|M2-like|M2 like",
      x_clean,
      ignore.case = TRUE
    ) ~ "Pro-resolution M2-like Mphi",

    grepl(
      "SPP1|TREM2|MASH-associated|MASH associated",
      x_clean,
      ignore.case = TRUE
    ) ~ "SPP1/TREM2 MASH-associated Mphi",

    grepl(
      "^Other$|Unclassified|Unknown|NA",
      x_clean,
      ignore.case = TRUE
    ) ~ "Other",

    TRUE ~ "Other"
  )

  factor(
    result,
    levels = mphi_subtype_levels
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

      axis.title = element_text(
        face = "bold",
        size = base_size
      ),

      axis.text = element_text(
        size = base_size - 1,
        colour = "black"
      ),

      legend.title = element_text(
        face = "bold",
        size = base_size - 1
      ),

      legend.text = element_text(
        size = base_size - 2
      ),

      strip.background = element_rect(
        fill = "#F4F4F4",
        colour = "#AAAAAA",
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

      plot.margin = margin(
        8,
        8,
        8,
        8
      )
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

      axis.title = element_text(
        face = "bold"
      ),

      axis.text.x = element_text(
        size = base_size + 1,
        colour = "black"
      ),

      axis.text.y = element_text(
        colour = "black"
      ),

      legend.title = element_text(
        face = "bold"
      ),

      panel.grid = element_blank()
    )
}


save_plot_pdf_png <- function(
    plot_object,
    filename_without_extension,
    width,
    height
) {

  pdf_file <- paste0(
    filename_without_extension,
    ".pdf"
  )

  ggsave(
    filename = pdf_file,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )

  if (isTRUE(save_png)) {

    png_file <- paste0(
      filename_without_extension,
      ".png"
    )

    ggsave(
      filename = png_file,
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


prepare_umap_data <- function(
    object,
    reduction,
    annotation_col,
    condition_col,
    sample_col
) {

  emb <- Embeddings(
    object,
    reduction = reduction
  )

  if (ncol(emb) < 2) {
    stop(
      "UMAP reductionに2次元以上の座標がありません。",
      call. = FALSE
    )
  }

  df <- data.frame(
    cell = rownames(emb),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    stringsAsFactors = FALSE
  )

  md <- object@meta.data

  md$cell <- rownames(md)

  df <- left_join(
    df,
    md,
    by = "cell"
  )

  df$Mphi_subtype <- normalize_mphi_annotation(
    df[[annotation_col]]
  )

  df$Condition_plot <- normalize_condition(
    df[[condition_col]]
  )

  df$Sample_plot <- as.character(
    df[[sample_col]]
  )

  df
}


make_umap_plot <- function(
    data,
    title,
    subtitle = NULL,
    facet_variable = NULL,
    show_legend = TRUE,
    fixed_limits = NULL
) {

  p <- ggplot(
    data,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = Mphi_subtype
    )
  ) +

    geom_point(
      size = umap_point_size,
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
      x = paste0(reduction_use, "_1"),
      y = paste0(reduction_use, "_2")
    ) +

    coord_fixed() +

    theme_ueno_umap(base_size = 12) +

    guides(
      colour = guide_legend(
        override.aes = list(
          size = 4,
          alpha = 1
        )
      )
    )

  if (!is.null(facet_variable)) {

    p <- p +
      facet_wrap(
        stats::as.formula(
          paste0("~", facet_variable)
        ),
        nrow = 1
      )
  }

  if (!is.null(fixed_limits)) {

    p <- p +
      coord_fixed(
        xlim = fixed_limits$x,
        ylim = fixed_limits$y,
        expand = TRUE
      )
  }

  if (!show_legend) {

    p <- p +
      theme(
        legend.position = "none"
      )
  }

  p
}


find_existing_score_column <- function(
    metadata,
    candidates
) {

  found <- candidates[
    candidates %in% colnames(metadata)
  ]

  if (length(found) == 0) {
    return(NA_character_)
  }

  found[1]
}


calculate_or_find_module_scores <- function(object) {

  metadata <- object@meta.data

  m1_score_candidates <- c(
    "M1_score",
    "M1Score",
    "M1_ModuleScore",
    "M1_module_score",
    "M1_marker_module",
    "M1_marker_module1"
  )

  m2_score_candidates <- c(
    "M2_score",
    "M2Score",
    "M2_ModuleScore",
    "M2_module_score",
    "M2_marker_module",
    "M2_marker_module1"
  )

  m1_col <- find_existing_score_column(
    metadata,
    m1_score_candidates
  )

  m2_col <- find_existing_score_column(
    metadata,
    m2_score_candidates
  )

  use_existing <-
    !is.na(m1_col) &&
    !is.na(m2_col) &&
    !isTRUE(recalculate_module_scores)

  if (use_existing) {

    write_log(
      paste0(
        "既存Module scoreを使用: ",
        m1_col,
        ", ",
        m2_col
      )
    )

    object$M1_UenoScore <- object@meta.data[[m1_col]]
    object$M2_UenoScore <- object@meta.data[[m2_col]]

    return(object)
  }

  write_log("M1/M2 module scoreを計算します")

  DefaultAssay(object) <- "RNA"

  feature_names <- rownames(object)

  m1_present <- intersect(
    m1_markers_mouse,
    feature_names
  )

  m2_present <- intersect(
    m2_markers_mouse,
    feature_names
  )

  write_log(
    paste0(
      "M1 genes found: ",
      length(m1_present),
      "/",
      length(m1_markers_mouse),
      " | ",
      paste(m1_present, collapse = ", ")
    )
  )

  write_log(
    paste0(
      "M2 genes found: ",
      length(m2_present),
      "/",
      length(m2_markers_mouse),
      " | ",
      paste(m2_present, collapse = ", ")
    )
  )

  if (length(m1_present) < 3) {
    stop(
      "M1 module scoreに使用できる遺伝子が3個未満です。",
      call. = FALSE
    )
  }

  if (length(m2_present) < 3) {
    stop(
      "M2 module scoreに使用できる遺伝子が3個未満です。",
      call. = FALSE
    )
  }

  object <- AddModuleScore(
    object = object,
    features = list(m1_present),
    name = "M1_UenoScore",
    assay = "RNA",
    seed = random_seed
  )

  object <- AddModuleScore(
    object = object,
    features = list(m2_present),
    name = "M2_UenoScore",
    assay = "RNA",
    seed = random_seed
  )

  object$M1_UenoScore <-
    object@meta.data[["M1_UenoScore1"]]

  object$M2_UenoScore <-
    object@meta.data[["M2_UenoScore1"]]

  object
}


############################################################
## 5. LOAD RDS
############################################################

write_log(
  paste0(
    "Palette: ",
    palette_name
  )
)

if (!file.exists(rds_file)) {

  stop(
    paste0(
      "RDSが見つかりません:\n",
      rds_file
    ),
    call. = FALSE
  )
}

write_log(
  paste0(
    "Loading RDS: ",
    rds_file
  )
)

obj <- readRDS(rds_file)

if (!inherits(obj, "Seurat")) {

  stop(
    "読み込んだファイルはSeurat objectではありません。",
    call. = FALSE
  )
}

write_log(
  paste0(
    "Cells: ",
    format(ncol(obj), big.mark = ","),
    " | Features: ",
    format(nrow(obj), big.mark = ",")
  )
)


############################################################
## 6. DETECT METADATA AND REDUCTION
############################################################

condition_col <- detect_metadata_column(
  obj,
  candidates = c(
    "condition",
    "Condition",
    "condition_FIXED2",
    "group",
    "Group",
    "diet",
    "Diet"
  ),
  label = "Condition"
)

sample_col <- detect_metadata_column(
  obj,
  candidates = c(
    "sample",
    "Sample",
    "sample_for_R8plot_FIXED2",
    "sample_for_R8plot",
    "orig.ident",
    "sample_id",
    "SampleID"
  ),
  label = "Sample"
)

annotation_col <- detect_metadata_column(
  obj,
  candidates = c(
    "Layer2_annotation",
    "layer2_annotation",
    "Mphi_subtype",
    "Mphi_annotation",
    "celltype_manual",
    "celltype_annotation",
    "celltype_for_R8plot_FIXED2",
    "celltype_for_R8plot",
    "celltype_auto_annotation"
  ),
  label = "Macrophage annotation"
)

reduction_use <- detect_reduction(obj)

write_log(
  paste0(
    "Condition column: ",
    condition_col
  )
)

write_log(
  paste0(
    "Sample column: ",
    sample_col
  )
)

write_log(
  paste0(
    "Annotation column: ",
    annotation_col
  )
)

write_log(
  paste0(
    "Reduction: ",
    reduction_use
  )
)


############################################################
## 7. FILTER STD AND CDAHFD
############################################################

obj$Condition_plot <- normalize_condition(
  obj@meta.data[[condition_col]]
)

cells_keep <- rownames(obj@meta.data)[
  obj$Condition_plot %in% conditions_keep
]

if (length(cells_keep) == 0) {

  stop(
    paste0(
      "STD/CDAHFD細胞が検出されませんでした。\n",
      "Condition values:\n",
      paste(
        unique(obj@meta.data[[condition_col]]),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

obj_sub <- subset(
  obj,
  cells = cells_keep
)

obj_sub$Condition_plot <- factor(
  normalize_condition(
    obj_sub@meta.data[[condition_col]]
  ),
  levels = conditions_keep
)

write_log(
  paste0(
    "STD/CDAHFD cells retained: ",
    format(ncol(obj_sub), big.mark = ",")
  )
)


############################################################
## 8. PREPARE UMAP DATA
############################################################

umap_df <- prepare_umap_data(
  object = obj_sub,
  reduction = reduction_use,
  annotation_col = annotation_col,
  condition_col = condition_col,
  sample_col = sample_col
)

umap_df$Condition_plot <- factor(
  umap_df$Condition_plot,
  levels = conditions_keep
)

umap_limits <- list(
  x = range(
    umap_df$UMAP_1,
    na.rm = TRUE
  ),

  y = range(
    umap_df$UMAP_2,
    na.rm = TRUE
  )
)

data.table::fwrite(
  umap_df,
  file.path(
    dir_table,
    "UMAP_coordinates_STD_CDAHFD.csv"
  )
)


############################################################
## 9. INTEGRATED UMAP
############################################################

p_umap_integrated <- make_umap_plot(
  data = umap_df,
  title = "Mouse Mphi UMAP: integrated data",
  subtitle = paste0(
    "STD and CDAHFD / ",
    palette_name,
    " / dot size ×2"
  ),
  show_legend = TRUE,
  fixed_limits = umap_limits
)

save_plot_pdf_png(
  p_umap_integrated,
  file.path(
    dir_umap,
    "UMAP_Mphi_integrated_Ueno_MASH_Macro_Palette_v1"
  ),
  width = 10.5,
  height = 7.5
)


############################################################
## 10. CONDITION-SPLIT UMAP
############################################################

p_umap_condition <- make_umap_plot(
  data = umap_df,
  title = "Mouse Mphi UMAP by condition",
  subtitle = paste0(
    "Shared integrated coordinates / ",
    palette_name
  ),
  facet_variable = "Condition_plot",
  show_legend = TRUE,
  fixed_limits = umap_limits
)

save_plot_pdf_png(
  p_umap_condition,
  file.path(
    dir_umap,
    "UMAP_Mphi_STD_vs_CDAHFD_shared_coordinates"
  ),
  width = 14,
  height = 7
)


############################################################
## 11. INDIVIDUAL CONDITION UMAP
############################################################

for (condition_now in conditions_keep) {

  df_now <- umap_df %>%
    filter(
      Condition_plot == condition_now
    )

  p_now <- make_umap_plot(
    data = df_now,
    title = paste0(
      "Mouse Mphi UMAP: ",
      condition_now
    ),
    subtitle = "Shared integrated UMAP coordinates",
    show_legend = TRUE,
    fixed_limits = umap_limits
  )

  save_plot_pdf_png(
    p_now,
    file.path(
      dir_umap,
      paste0(
        "UMAP_Mphi_",
        condition_now,
        "_shared_coordinates"
      )
    ),
    width = 9,
    height = 7
  )
}


############################################################
## 12. SAMPLE-SPLIT UMAP
############################################################

if (isTRUE(make_sample_umap)) {

  n_samples <- length(
    unique(umap_df$Sample_plot)
  )

  sample_rows <- ceiling(
    n_samples / 3
  )

  p_umap_sample <- ggplot(
    umap_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = Mphi_subtype
    )
  ) +

    geom_point(
      size = umap_point_size,
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
      ~Sample_plot,
      ncol = 3
    ) +

    coord_fixed(
      xlim = umap_limits$x,
      ylim = umap_limits$y
    ) +

    labs(
      title = "Mouse Mphi UMAP by sample",
      subtitle = "Shared integrated coordinates",
      x = paste0(reduction_use, "_1"),
      y = paste0(reduction_use, "_2")
    ) +

    theme_ueno_umap(base_size = 11) +

    guides(
      colour = guide_legend(
        override.aes = list(
          size = 4
        )
      )
    )

  save_plot_pdf_png(
    p_umap_sample,
    file.path(
      dir_umap,
      "UMAP_Mphi_by_sample_shared_coordinates"
    ),
    width = 14,
    height = max(6, sample_rows * 4.5)
  )
}


############################################################
## 13. MODULE SCORES
############################################################

obj_sub <- calculate_or_find_module_scores(
  obj_sub
)

module_df <- obj_sub@meta.data %>%
  transmute(
    Condition = factor(
      Condition_plot,
      levels = conditions_keep
    ),

    M1_score = M1_UenoScore,
    M2_score = M2_UenoScore
  ) %>%

  pivot_longer(
    cols = c(
      M1_score,
      M2_score
    ),
    names_to = "Module",
    values_to = "Score"
  ) %>%

  mutate(
    Module = recode(
      Module,
      "M1_score" = "M1 marker module",
      "M2_score" = "M2 marker module"
    )
  )

data.table::fwrite(
  module_df,
  file.path(
    dir_table,
    "M1_M2_module_scores_cell_level.csv"
  )
)

p_module <- ggplot(
  module_df,
  aes(
    x = Condition,
    y = Score,
    fill = Condition
  )
) +

  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.45,
    colour = "#333333"
  ) +

  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    fill = NA,
    colour = "#202020",
    linewidth = 0.5
  ) +

  facet_wrap(
    ~Module,
    scales = "free_y",
    nrow = 1
  ) +

  scale_fill_manual(
    values = condition_palette,
    drop = FALSE
  ) +

  labs(
    title = "M1 and M2 marker modules in mouse Mphi",
    subtitle = "STD versus CDAHFD / all macrophages",
    x = NULL,
    y = "Module expression"
  ) +

  theme_classic(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    ),

    strip.background = element_rect(
      fill = "#F4F4F4",
      colour = "#777777"
    ),

    strip.text = element_text(
      face = "bold"
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 13,
      colour = "black"
    ),

    legend.position = "right"
  )

save_plot_pdf_png(
  p_module,
  file.path(
    dir_module,
    "ModuleScore_M1_M2_STD_vs_CDAHFD"
  ),
  width = 8.5,
  height = 6.5
)


############################################################
## 14. CELL COUNTS
############################################################

cell_class_df <- umap_df %>%
  mutate(
    M1_M2_class = case_when(

      Mphi_subtype ==
        "Inflammatory M1-like Mphi" ~
        "M1 cells",

      Mphi_subtype ==
        "Pro-resolution M2-like Mphi" ~
        "M2 cells",

      TRUE ~ NA_character_
    )
  ) %>%

  filter(
    !is.na(M1_M2_class)
  ) %>%

  count(
    Condition_plot,
    M1_M2_class,
    name = "Cell_count"
  ) %>%

  complete(
    Condition_plot = factor(
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

cell_class_df$M1_M2_class <- factor(
  cell_class_df$M1_M2_class,
  levels = c(
    "M2 cells",
    "M1 cells"
  )
)

data.table::fwrite(
  cell_class_df,
  file.path(
    dir_table,
    "M1_M2_cell_counts.csv"
  )
)

p_counts <- ggplot(
  cell_class_df,
  aes(
    x = Condition_plot,
    y = Cell_count,
    fill = M1_M2_class
  )
) +

  geom_col(
    width = 0.55
  ) +

  scale_fill_manual(
    values = m1_m2_palette[
      c(
        "M2 cells",
        "M1 cells"
      )
    ],
    drop = FALSE,
    name = "Cell type"
  ) +

  labs(
    title = "M1, M2 Mac cell counts",
    x = NULL,
    y = "Cell count"
  ) +

  theme_ueno_bar(base_size = 13)

save_plot_pdf_png(
  p_counts,
  file.path(
    dir_count,
    "M1_M2_cell_counts_STD_vs_CDAHFD"
  ),
  width = 6.5,
  height = 6.5
)


############################################################
## 15. M2/M1 RATIO
############################################################

ratio_df <- cell_class_df %>%
  mutate(
    class_short = case_when(
      M1_M2_class == "M1 cells" ~ "M1",
      M1_M2_class == "M2 cells" ~ "M2"
    )
  ) %>%

  select(
    Condition_plot,
    class_short,
    Cell_count
  ) %>%

  pivot_wider(
    names_from = class_short,
    values_from = Cell_count,
    values_fill = 0
  ) %>%

  mutate(
    M2_M1_ratio =
      (M2 + 0.5) /
      (M1 + 0.5),

    Log2_M2_M1_ratio =
      log2(M2_M1_ratio)
  )

data.table::fwrite(
  ratio_df,
  file.path(
    dir_table,
    "M2_M1_ratio_STD_vs_CDAHFD.csv"
  )
)

p_ratio <- ggplot(
  ratio_df,
  aes(
    x = Condition_plot,
    y = Log2_M2_M1_ratio,
    fill = Condition_plot
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "#888888",
    linewidth = 0.6
  ) +

  geom_col(
    width = 0.55,
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
    size = 4.3
  ) +

  scale_fill_manual(
    values = condition_palette,
    drop = FALSE
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.12,
        0.16
      )
    )
  ) +

  labs(
    title = "M2/M1 cell ratio in total Mac",
    x = NULL,
    y = "Log2 ratio"
  ) +

  theme_ueno_bar(base_size = 13) +

  theme(
    legend.position = "none"
  )

save_plot_pdf_png(
  p_ratio,
  file.path(
    dir_ratio,
    "M2_M1_log2_ratio_STD_vs_CDAHFD"
  ),
  width = 6.5,
  height = 6.5
)


############################################################
## 16. PROFIBROTIC MACROPHAGE PERCENTAGE
############################################################

## 今回はSPP1/TREM2 MASH-associated Mphiを
## Profibrotic macrophageとして集計
profibrotic_df <- umap_df %>%
  group_by(
    Condition_plot
  ) %>%

  summarise(
    Total_Mphi = n(),

    Profibrotic_Mphi = sum(
      Mphi_subtype ==
        "SPP1/TREM2 MASH-associated Mphi",
      na.rm = TRUE
    ),

    Profibrotic_percent =
      Profibrotic_Mphi /
      Total_Mphi *
      100,

    .groups = "drop"
  )

data.table::fwrite(
  profibrotic_df,
  file.path(
    dir_table,
    "Profibrotic_Mphi_percentage_STD_vs_CDAHFD.csv"
  )
)

p_profibrotic <- ggplot(
  profibrotic_df,
  aes(
    x = Condition_plot,
    y = Profibrotic_percent,
    fill = Condition_plot
  )
) +

  geom_col(
    width = 0.55,
    colour = "#333333",
    linewidth = 0.4
  ) +

  geom_text(
    aes(
      label = paste0(
        sprintf(
          "%.2f%%",
          Profibrotic_percent
        ),
        "\n",
        Profibrotic_Mphi,
        "/",
        Total_Mphi
      )
    ),
    vjust = -0.35,
    fontface = "bold",
    size = 4
  ) +

  scale_fill_manual(
    values = condition_palette,
    drop = FALSE
  ) +

  scale_y_continuous(
    labels = label_number(
      accuracy = 0.1
    ),

    expand = expansion(
      mult = c(
        0,
        0.18
      )
    )
  ) +

  labs(
    title = "Profibrotic Mac % (/total M\u03c6)",
    x = NULL,
    y = "Profibrotic macrophage / total Mphi (%)"
  ) +

  theme_ueno_bar(base_size = 13) +

  theme(
    legend.position = "none"
  )

save_plot_pdf_png(
  p_profibrotic,
  file.path(
    dir_ratio,
    "Profibrotic_Mphi_percentage_STD_vs_CDAHFD"
  ),
  width = 6.5,
  height = 6.5
)


############################################################
## 17. COMBINED FIGURE WITHOUT SLIDE DECORATION
############################################################

p_umap_integrated_small <-
  p_umap_integrated +

  theme(
    legend.position = "right"
  )

p_umap_condition_small <-
  p_umap_condition +

  theme(
    legend.position = "none"
  )

p_module_small <-
  p_module +

  theme(
    legend.position = "none"
  )

combined_figure <-
  (
    p_umap_integrated_small |
      p_umap_condition_small
  ) /

  (
    p_module_small |
      p_counts |
      p_ratio |
      p_profibrotic
  ) +

  plot_layout(
    heights = c(
      1.25,
      1
    )
  ) +

  plot_annotation(
    title =
      "Mouse MASH model liver macrophage analysis",

    subtitle =
      paste0(
        "STD versus CDAHFD | ",
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

save_plot_pdf_png(
  combined_figure,
  file.path(
    output_dir,
    "Mouse_Mphi_STD_vs_CDAHFD_combined_figure"
  ),
  width = 22,
  height = 13
)


############################################################
## 18. EXCEL SUMMARY
############################################################

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

metadata_summary <- data.frame(
  Item = c(
    "RDS file",
    "Condition column",
    "Sample column",
    "Annotation column",
    "UMAP reduction",
    "UMAP point size",
    "Palette",
    "Number of cells",
    "Analysis date"
  ),

  Value = c(
    rds_file,
    condition_col,
    sample_col,
    annotation_col,
    reduction_use,
    umap_point_size,
    palette_name,
    ncol(obj_sub),
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  stringsAsFactors = FALSE
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "AnalysisSettings"
)

writeData(
  wb,
  "AnalysisSettings",
  metadata_summary
)

addWorksheet(
  wb,
  "Palette"
)

writeData(
  wb,
  "Palette",
  palette_table
)

addWorksheet(
  wb,
  "M1_M2_counts"
)

writeData(
  wb,
  "M1_M2_counts",
  cell_class_df
)

addWorksheet(
  wb,
  "M2_M1_ratio"
)

writeData(
  wb,
  "M2_M1_ratio",
  ratio_df
)

addWorksheet(
  wb,
  "Profibrotic_percent"
)

writeData(
  wb,
  "Profibrotic_percent",
  profibrotic_df
)

addWorksheet(
  wb,
  "M1_markers"
)

writeData(
  wb,
  "M1_markers",
  data.frame(
    gene = m1_markers_mouse
  )
)

addWorksheet(
  wb,
  "M2_markers"
)

writeData(
  wb,
  "M2_markers",
  data.frame(
    gene = m2_markers_mouse
  )
)

saveWorkbook(
  wb,
  file.path(
    dir_table,
    "Mouse_Mphi_STD_vs_CDAHFD_summary.xlsx"
  ),
  overwrite = TRUE
)


############################################################
## 19. SAVE PALETTE AS RDS AND R SCRIPT
############################################################

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
  con = file.path(
    dir_table,
    "Ueno_MASH_Macro_Palette_v1.R"
  )
)


############################################################
## 20. SESSION INFORMATION
############################################################

sink(
  file.path(
    dir_log,
    "sessionInfo.txt"
  )
)

cat(
  "Analysis completed:\n",
  format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  ),
  "\n\n"
)

cat(
  "Palette:\n",
  palette_name,
  "\n\n"
)

cat(
  "RDS:\n",
  rds_file,
  "\n\n"
)

cat(
  "Condition column:\n",
  condition_col,
  "\n\n"
)

cat(
  "Sample column:\n",
  sample_col,
  "\n\n"
)

cat(
  "Annotation column:\n",
  annotation_col,
  "\n\n"
)

cat(
  "Reduction:\n",
  reduction_use,
  "\n\n"
)

print(
  sessionInfo()
)

sink()


############################################################
## 21. FINISH
############################################################

write_log("Analysis completed successfully")

cat("\n")
cat("====================================================\n")
cat("Mouse Mphi STD versus CDAHFD analysis completed\n")
cat("====================================================\n")
cat("RDS               :", rds_file, "\n")
cat("Condition column  :", condition_col, "\n")
cat("Sample column     :", sample_col, "\n")
cat("Annotation column :", annotation_col, "\n")
cat("Reduction         :", reduction_use, "\n")
cat("Palette           :", palette_name, "\n")
cat("Output directory  :", output_dir, "\n")
cat("====================================================\n")

