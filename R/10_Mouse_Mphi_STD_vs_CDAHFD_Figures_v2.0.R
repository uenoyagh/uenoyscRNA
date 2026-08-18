############################################################
## 10_Mouse_Mphi_STD_vs_CDAHFD_Figures_v2.0.R
##
## Mouse MASH liver macrophage analysis
## STD versus CDAHFD
##
## Major correction from v1:
##   - Macrophage cells are extracted explicitly
##   - Layer2 annotation column is detected by its contents
##   - Analysis stops if the six-class annotation is absent
##   - No silent conversion of all unmatched cells to "Other"
##
## Palette:
##   Ueno MASH Macro Palette v1.0
############################################################


############################################################
## 0. SETTINGS
############################################################

rds_file <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
)

output_dir <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_Mphi_STD_vs_CDAHFD_Figures_v2.0"
)

conditions_keep <- c(
  "STD",
  "CDAHFD"
)

umap_point_size <- 1.25
png_dpi <- 400
save_png <- TRUE
make_sample_umap <- TRUE
random_seed <- 1234

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
      "不足パッケージ:\n",
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

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir_umap <- file.path(
  output_dir,
  "01_UMAP"
)

dir_module <- file.path(
  output_dir,
  "02_ModuleScore"
)

dir_count <- file.path(
  output_dir,
  "03_CellCounts"
)

dir_ratio <- file.path(
  output_dir,
  "04_Ratios"
)

dir_table <- file.path(
  output_dir,
  "05_Tables"
)

dir_log <- file.path(
  output_dir,
  "Logs"
)

for (x in c(
  dir_umap,
  dir_module,
  dir_count,
  dir_ratio,
  dir_table,
  dir_log
)) {
  dir.create(
    x,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


############################################################
## 3. LOG FUNCTION
############################################################

write_log <- function(text) {

  line <- paste0(
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    " | ",
    text
  )

  cat(line, "\n")

  cat(
    line,
    "\n",
    file = file.path(
      dir_log,
      "analysis.log"
    ),
    append = TRUE
  )
}


############################################################
## 4. PALETTE
############################################################

palette_name <- "Ueno MASH Macro Palette v1.0"

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
  "M2 cells" = "#3366F5",
  "M1 cells" = "#FF4A4A"
)

mphi_subtype_levels <- c(
  "Resident Kupffer-like Mphi",
  "Monocyte-like Mphi",
  "Inflammatory M1-like Mphi",
  "Pro-resolution M2-like Mphi",
  "SPP1/TREM2 MASH-associated Mphi",
  "Other"
)


############################################################
## 5. MARKER MODULES
############################################################

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


############################################################
## 6. BASIC FUNCTIONS
############################################################

normalize_condition <- function(x) {

  x <- as.character(x)

  dplyr::case_when(
    grepl(
      "CDAHFD|CDHFD|CDAH",
      x,
      ignore.case = TRUE
    ) ~ "CDAHFD",

    grepl(
      "STD|standard",
      x,
      ignore.case = TRUE
    ) ~ "STD",

    TRUE ~ x
  )
}


detect_column <- function(
    object,
    candidates,
    label
) {

  found <- candidates[
    candidates %in%
      colnames(object@meta.data)
  ]

  if (length(found) == 0) {

    stop(
      paste0(
        label,
        "列を検出できませんでした。\n\n",
        paste(
          colnames(object@meta.data),
          collapse = "\n"
        )
      ),
      call. = FALSE
    )
  }

  found[1]
}


detect_reduction <- function(object) {

  reductions_available <- Reductions(object)

  preferred <- c(
    "MphiUMAP",
    "mphiUMAP",
    "mphi_umap",
    "umapRPCA",
    "umap",
    "UMAP"
  )

  found <- preferred[
    preferred %in% reductions_available
  ]

  if (length(found) > 0) {
    return(found[1])
  }

  found <- reductions_available[
    grepl(
      "umap",
      reductions_available,
      ignore.case = TRUE
    )
  ]

  if (length(found) == 0) {
    stop(
      "UMAP reductionが見つかりません。",
      call. = FALSE
    )
  }

  found[1]
}


############################################################
## 7. SIX-CLASS ANNOTATION NORMALIZATION
############################################################

normalize_layer2_annotation <- function(x) {

  x <- trimws(
    as.character(x)
  )

  result <- rep(
    NA_character_,
    length(x)
  )

  result[
    grepl(
      "Resident.*Kupffer|Kupffer.*Resident|Resident Kupffer",
      x,
      ignore.case = TRUE
    )
  ] <- "Resident Kupffer-like Mphi"

  result[
    is.na(result) &
      grepl(
        "Monocyte",
        x,
        ignore.case = TRUE
      )
  ] <- "Monocyte-like Mphi"

  result[
    is.na(result) &
      grepl(
        "Inflammatory|M1[- ]?like",
        x,
        ignore.case = TRUE
      )
  ] <- "Inflammatory M1-like Mphi"

  result[
    is.na(result) &
      grepl(
        "Pro[- ]?resolution|M2[- ]?like",
        x,
        ignore.case = TRUE
      )
  ] <- "Pro-resolution M2-like Mphi"

  result[
    is.na(result) &
      grepl(
        "SPP1|TREM2|MASH[- ]?associated",
        x,
        ignore.case = TRUE
      )
  ] <- "SPP1/TREM2 MASH-associated Mphi"

  result[
    is.na(result) &
      grepl(
        "^Other$|Unclassified|Unknown",
        x,
        ignore.case = TRUE
      )
  ] <- "Other"

  factor(
    result,
    levels = mphi_subtype_levels
  )
}


############################################################
## 8. FIND THE REAL LAYER2 COLUMN BY CONTENT
############################################################

score_layer2_column <- function(x) {

  x <- unique(
    na.omit(
      as.character(x)
    )
  )

  if (length(x) == 0) {
    return(0)
  }

  patterns <- c(
    "Resident.*Kupffer|Kupffer.*Resident",
    "Monocyte",
    "Inflammatory|M1[- ]?like",
    "Pro[- ]?resolution|M2[- ]?like",
    "SPP1|TREM2|MASH[- ]?associated"
  )

  sum(
    vapply(
      patterns,
      function(pattern_now) {
        any(
          grepl(
            pattern_now,
            x,
            ignore.case = TRUE
          )
        )
      },
      logical(1)
    )
  )
}


find_layer2_column <- function(object) {

  md <- object@meta.data

  preferred_candidates <- c(
    "Layer2_annotation",
    "layer2_annotation",
    "Layer2",
    "layer2",
    "Mphi_subtype",
    "Mphi_annotation",
    "Mphi_annotation_layer2",
    "celltype_layer2",
    "celltype_manual_layer2"
  )

  preferred_found <- preferred_candidates[
    preferred_candidates %in%
      colnames(md)
  ]

  if (length(preferred_found) > 0) {

    preferred_scores <- vapply(
      preferred_found,
      function(column_now) {
        score_layer2_column(
          md[[column_now]]
        )
      },
      numeric(1)
    )

    if (max(preferred_scores) >= 3) {

      selected <- preferred_found[
        which.max(preferred_scores)
      ]

      return(selected)
    }
  }

  character_columns <- colnames(md)[
    vapply(
      md,
      function(x) {
        is.character(x) ||
          is.factor(x)
      },
      logical(1)
    )
  ]

  all_scores <- vapply(
    character_columns,
    function(column_now) {
      score_layer2_column(
        md[[column_now]]
      )
    },
    numeric(1)
  )

  score_table <- data.frame(
    column = character_columns,
    score = all_scores,
    stringsAsFactors = FALSE
  ) %>%
    arrange(
      desc(score)
    )

  data.table::fwrite(
    score_table,
    file.path(
      dir_table,
      "Layer2_annotation_column_scores.csv"
    )
  )

  if (
    nrow(score_table) == 0 ||
    max(score_table$score) < 3
  ) {

    stop(
      paste0(
        "\n6分類Layer2 annotation列を検出できませんでした。\n",
        "全細胞RDSにLayer2分類が保存されていない可能性があります。\n\n",
        "候補列のスコアは次に保存しました:\n",
        file.path(
          dir_table,
          "Layer2_annotation_column_scores.csv"
        ),
        "\n\n",
        "この状態で解析を継続すると、M1/M2/SPP1-TREM2が",
        "全て0になるため、処理を停止しました。"
      ),
      call. = FALSE
    )
  }

  score_table$column[1]
}


############################################################
## 9. DETECT MACROPHAGE LAYER1 COLUMN
############################################################

score_macrophage_column <- function(x) {

  x <- as.character(x)

  sum(
    grepl(
      "Mphi|Macrophage|Kupffer|Monocyte",
      x,
      ignore.case = TRUE
    ),
    na.rm = TRUE
  )
}


find_macrophage_column <- function(
    object,
    exclude_column = NULL
) {

  md <- object@meta.data

  preferred_candidates <- c(
    "celltype_for_R8plot_FIXED2",
    "celltype_for_R8plot",
    "celltype_annotation",
    "celltype",
    "CellType",
    "Layer1_annotation",
    "layer1_annotation",
    "Layer1",
    "layer1"
  )

  preferred_candidates <- setdiff(
    preferred_candidates,
    exclude_column
  )

  found <- preferred_candidates[
    preferred_candidates %in%
      colnames(md)
  ]

  if (length(found) > 0) {

    scores <- vapply(
      found,
      function(column_now) {
        score_macrophage_column(
          md[[column_now]]
        )
      },
      numeric(1)
    )

    if (max(scores) > 0) {
      return(
        found[
          which.max(scores)
        ]
      )
    }
  }

  character_columns <- colnames(md)[
    vapply(
      md,
      function(x) {
        is.character(x) ||
          is.factor(x)
      },
      logical(1)
    )
  ]

  character_columns <- setdiff(
    character_columns,
    exclude_column
  )

  scores <- vapply(
    character_columns,
    function(column_now) {
      score_macrophage_column(
        md[[column_now]]
      )
    },
    numeric(1)
  )

  if (length(scores) == 0 || max(scores) == 0) {
    return(NA_character_)
  }

  character_columns[
    which.max(scores)
  ]
}


############################################################
## 10. THEMES AND EXPORT
############################################################

theme_ueno_umap <- function(base_size = 12) {

  theme_classic(
    base_size = base_size
  ) +

    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = base_size + 2
      ),

      plot.subtitle = element_text(
        hjust = 0.5,
        size = base_size - 1
      ),

      axis.title = element_text(
        face = "bold"
      ),

      axis.text = element_text(
        colour = "black"
      ),

      legend.title = element_text(
        face = "bold"
      ),

      strip.background = element_rect(
        fill = "#F4F4F4",
        colour = "#999999"
      ),

      strip.text = element_text(
        face = "bold"
      ),

      panel.border = element_rect(
        colour = "#555555",
        fill = NA,
        linewidth = 0.5
      )
    )
}


theme_ueno_bar <- function(base_size = 13) {

  theme_classic(
    base_size = base_size
  ) +

    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),

      axis.title = element_text(
        face = "bold"
      ),

      axis.text = element_text(
        colour = "black"
      ),

      legend.title = element_text(
        face = "bold"
      )
    )
}


save_plot <- function(
    plot_object,
    filename,
    width,
    height
) {

  ggsave(
    paste0(
      filename,
      ".pdf"
    ),
    plot_object,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    limitsize = FALSE
  )

  if (isTRUE(save_png)) {

    ggsave(
      paste0(
        filename,
        ".png"
      ),
      plot_object,
      width = width,
      height = height,
      units = "in",
      dpi = png_dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
}


############################################################
## 11. LOAD OBJECT
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

obj <- readRDS(
  rds_file
)

if (!inherits(obj, "Seurat")) {
  stop(
    "Seurat objectではありません。",
    call. = FALSE
  )
}

write_log(
  paste0(
    "Cells: ",
    format(
      ncol(obj),
      big.mark = ","
    ),
    " | Features: ",
    format(
      nrow(obj),
      big.mark = ","
    )
  )
)


############################################################
## 12. METADATA DETECTION
############################################################

condition_col <- detect_column(
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

sample_col <- detect_column(
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

reduction_use <- detect_reduction(
  obj
)

layer2_col <- find_layer2_column(
  obj
)

layer1_col <- find_macrophage_column(
  obj,
  exclude_column = layer2_col
)

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
    "Layer1 macrophage column: ",
    layer1_col
  )
)

write_log(
  paste0(
    "Layer2 annotation column: ",
    layer2_col
  )
)

write_log(
  paste0(
    "Reduction: ",
    reduction_use
  )
)


############################################################
## 13. DIAGNOSTIC TABLES
############################################################

layer2_values <- data.frame(
  value = as.character(
    obj@meta.data[[layer2_col]]
  ),
  stringsAsFactors = FALSE
) %>%
  count(
    value,
    name = "n",
    sort = TRUE
  )

data.table::fwrite(
  layer2_values,
  file.path(
    dir_table,
    "Layer2_annotation_original_values.csv"
  )
)

if (!is.na(layer1_col)) {

  layer1_values <- data.frame(
    value = as.character(
      obj@meta.data[[layer1_col]]
    ),
    stringsAsFactors = FALSE
  ) %>%
    count(
      value,
      name = "n",
      sort = TRUE
    )

  data.table::fwrite(
    layer1_values,
    file.path(
      dir_table,
      "Layer1_annotation_original_values.csv"
    )
  )
}


############################################################
## 14. DEFINE MΦ SUBTYPES
############################################################

obj$Mphi_subtype_v2 <- normalize_layer2_annotation(
  obj@meta.data[[layer2_col]]
)

mapping_check <- data.frame(
  original = as.character(
    obj@meta.data[[layer2_col]]
  ),

  normalized = as.character(
    obj$Mphi_subtype_v2
  ),

  stringsAsFactors = FALSE
) %>%
  count(
    original,
    normalized,
    name = "n",
    sort = TRUE
  )

data.table::fwrite(
  mapping_check,
  file.path(
    dir_table,
    "Layer2_annotation_mapping_check.csv"
  )
)


############################################################
## 15. FILTER MACROPHAGES AND CONDITIONS
############################################################

obj$Condition_plot <- normalize_condition(
  obj@meta.data[[condition_col]]
)

## Layer2分類が付与されている細胞をMφと定義
is_mphi_by_layer2 <- !is.na(
  obj$Mphi_subtype_v2
)

## Layer1列がある場合は補助的に利用
if (!is.na(layer1_col)) {

  is_mphi_by_layer1 <- grepl(
    "Mphi|Macrophage|Kupffer|Monocyte",
    as.character(
      obj@meta.data[[layer1_col]]
    ),
    ignore.case = TRUE
  )

} else {

  is_mphi_by_layer1 <- rep(
    FALSE,
    ncol(obj)
  )
}

cells_keep <- colnames(obj)[
  (
    is_mphi_by_layer2 |
      is_mphi_by_layer1
  ) &
    obj$Condition_plot %in%
    conditions_keep
]

if (length(cells_keep) == 0) {

  stop(
    paste0(
      "STD/CDAHFD macrophageが検出されませんでした。\n",
      "Layer1 column: ",
      layer1_col,
      "\nLayer2 column: ",
      layer2_col
    ),
    call. = FALSE
  )
}

obj_sub <- subset(
  obj,
  cells = cells_keep
)

## Layer1だけで拾われ、Layer2が空欄の細胞はOther
obj_sub$Mphi_subtype_v2 <- as.character(
  obj_sub$Mphi_subtype_v2
)

obj_sub$Mphi_subtype_v2[
  is.na(
    obj_sub$Mphi_subtype_v2
  ) |
    obj_sub$Mphi_subtype_v2 == ""
] <- "Other"

obj_sub$Mphi_subtype_v2 <- factor(
  obj_sub$Mphi_subtype_v2,
  levels = mphi_subtype_levels
)

obj_sub$Condition_plot <- factor(
  normalize_condition(
    obj_sub@meta.data[[condition_col]]
  ),
  levels = conditions_keep
)

write_log(
  paste0(
    "STD/CDAHFD macrophages retained: ",
    format(
      ncol(obj_sub),
      big.mark = ","
    )
  )
)


############################################################
## 16. VALIDATE CLASS COUNTS
############################################################

annotation_counts <- obj_sub@meta.data %>%
  count(
    Condition_plot,
    Mphi_subtype_v2,
    name = "Cell_count"
  ) %>%
  complete(
    Condition_plot = factor(
      conditions_keep,
      levels = conditions_keep
    ),

    Mphi_subtype_v2 = factor(
      mphi_subtype_levels,
      levels = mphi_subtype_levels
    ),

    fill = list(
      Cell_count = 0
    )
  )

data.table::fwrite(
  annotation_counts,
  file.path(
    dir_table,
    "Mphi_subtype_counts_STD_CDAHFD.csv"
  )
)

non_other_count <- sum(
  annotation_counts$Cell_count[
    annotation_counts$Mphi_subtype_v2 !=
      "Other"
  ]
)

if (non_other_count == 0) {

  stop(
    paste0(
      "\nLayer2分類の非Other細胞が0件です。\n",
      "誤ったmetadata列を選択している可能性があります。\n\n",
      "確認ファイル:\n",
      file.path(
        dir_table,
        "Layer2_annotation_mapping_check.csv"
      )
    ),
    call. = FALSE
  )
}

write_log(
  paste0(
    "Non-Other macrophages: ",
    non_other_count
  )
)


############################################################
## 17. UMAP DATA
############################################################

emb <- Embeddings(
  obj_sub,
  reduction = reduction_use
)

umap_df <- data.frame(
  cell = rownames(emb),
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2],
  stringsAsFactors = FALSE
)

md <- obj_sub@meta.data
md$cell <- rownames(md)

umap_df <- left_join(
  umap_df,
  md,
  by = "cell"
)

umap_df$Mphi_subtype <- factor(
  umap_df$Mphi_subtype_v2,
  levels = mphi_subtype_levels
)

umap_df$Condition_plot <- factor(
  umap_df$Condition_plot,
  levels = conditions_keep
)

umap_df$Sample_plot <- as.character(
  umap_df[[sample_col]]
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
    "UMAP_coordinates_macrophages_v2.csv"
  )
)


############################################################
## 18. UMAP FUNCTION
############################################################

make_umap <- function(
    data,
    title,
    subtitle = NULL,
    facet_col = NULL,
    show_legend = TRUE
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
      values =
        ueno_mash_macro_palette_v1,

      limits =
        mphi_subtype_levels,

      drop = FALSE,

      name =
        "Mphi subtype"
    ) +

    labs(
      title = title,
      subtitle = subtitle,
      x = paste0(
        reduction_use,
        "_1"
      ),
      y = paste0(
        reduction_use,
        "_2"
      )
    ) +

    coord_fixed(
      xlim = umap_limits$x,
      ylim = umap_limits$y
    ) +

    theme_ueno_umap(
      base_size = 12
    ) +

    guides(
      colour = guide_legend(
        override.aes = list(
          size = 4,
          alpha = 1
        )
      )
    )

  if (!is.null(facet_col)) {

    p <- p +
      facet_wrap(
        stats::as.formula(
          paste0(
            "~",
            facet_col
          )
        ),
        nrow = 1
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


############################################################
## 19. INTEGRATED UMAP
############################################################

p_umap_integrated <- make_umap(
  data = umap_df,
  title = "Mouse Mphi UMAP: integrated data",
  subtitle = paste0(
    "STD and CDAHFD / ",
    palette_name,
    " / macrophages only"
  ),
  show_legend = TRUE
)

save_plot(
  p_umap_integrated,
  file.path(
    dir_umap,
    "UMAP_Mphi_integrated_v2"
  ),
  width = 10,
  height = 7.5
)


############################################################
## 20. CONDITION-SPLIT UMAP
############################################################

p_umap_condition <- make_umap(
  data = umap_df,
  title = "Mouse Mphi UMAP by condition",
  subtitle =
    "STD versus CDAHFD / shared integrated coordinates",
  facet_col = "Condition_plot",
  show_legend = TRUE
)

save_plot(
  p_umap_condition,
  file.path(
    dir_umap,
    "UMAP_Mphi_STD_vs_CDAHFD_v2"
  ),
  width = 14,
  height = 7
)


############################################################
## 21. INDIVIDUAL CONDITION UMAP
############################################################

for (condition_now in conditions_keep) {

  df_now <- umap_df %>%
    filter(
      Condition_plot ==
        condition_now
    )

  p_now <- make_umap(
    data = df_now,
    title = paste0(
      "Mouse Mphi UMAP: ",
      condition_now
    ),
    subtitle =
      "Shared integrated coordinates",
    show_legend = TRUE
  )

  save_plot(
    p_now,
    file.path(
      dir_umap,
      paste0(
        "UMAP_Mphi_",
        condition_now,
        "_v2"
      )
    ),
    width = 9,
    height = 7
  )
}


############################################################
## 22. SAMPLE UMAP
############################################################

if (isTRUE(make_sample_umap)) {

  n_samples <- length(
    unique(
      umap_df$Sample_plot
    )
  )

  n_rows <- ceiling(
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
      values =
        ueno_mash_macro_palette_v1,

      limits =
        mphi_subtype_levels,

      drop = FALSE,

      name =
        "Mphi subtype"
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
      title =
        "Mouse Mphi UMAP by sample",

      subtitle =
        "Shared integrated coordinates",

      x = paste0(
        reduction_use,
        "_1"
      ),

      y = paste0(
        reduction_use,
        "_2"
      )
    ) +

    theme_ueno_umap(
      base_size = 11
    ) +

    guides(
      colour = guide_legend(
        override.aes = list(
          size = 4
        )
      )
    )

  save_plot(
    p_umap_sample,
    file.path(
      dir_umap,
      "UMAP_Mphi_by_sample_v2"
    ),
    width = 14,
    height = max(
      6,
      n_rows * 4.5
    )
  )
}


############################################################
## 23. MODULE SCORES
############################################################

DefaultAssay(
  obj_sub
) <- "RNA"

m1_present <- intersect(
  m1_markers_mouse,
  rownames(obj_sub)
)

m2_present <- intersect(
  m2_markers_mouse,
  rownames(obj_sub)
)

write_log(
  paste0(
    "M1 genes found: ",
    length(m1_present),
    "/",
    length(m1_markers_mouse),
    " | ",
    paste(
      m1_present,
      collapse = ", "
    )
  )
)

write_log(
  paste0(
    "M2 genes found: ",
    length(m2_present),
    "/",
    length(m2_markers_mouse),
    " | ",
    paste(
      m2_present,
      collapse = ", "
    )
  )
)

obj_sub <- AddModuleScore(
  object = obj_sub,
  features = list(
    m1_present
  ),
  name = "M1_UenoScore",
  assay = "RNA",
  seed = random_seed
)

obj_sub <- AddModuleScore(
  object = obj_sub,
  features = list(
    m2_present
  ),
  name = "M2_UenoScore",
  assay = "RNA",
  seed = random_seed
)

obj_sub$M1_score <-
  obj_sub$M1_UenoScore1

obj_sub$M2_score <-
  obj_sub$M2_UenoScore1

module_df <- obj_sub@meta.data %>%
  transmute(
    Condition = factor(
      Condition_plot,
      levels = conditions_keep
    ),

    M1_score = M1_score,
    M2_score = M2_score
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
      "M1_score" =
        "M1 marker module",
      "M2_score" =
        "M2 marker module"
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
    colour = "#333333",
    linewidth = 0.4
  ) +

  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    fill = NA,
    colour = "#202020"
  ) +

  facet_wrap(
    ~Module,
    scales = "free_y",
    nrow = 1
  ) +

  scale_fill_manual(
    values = condition_palette
  ) +

  labs(
    title =
      "M1 and M2 marker modules in mouse Mphi",

    subtitle =
      "STD versus CDAHFD / macrophages only",

    x = NULL,

    y =
      "Module expression"
  ) +

  theme_classic(
    base_size = 13
  ) +

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
      colour = "black"
    )
  )

save_plot(
  p_module,
  file.path(
    dir_module,
    "ModuleScore_M1_M2_v2"
  ),
  width = 8.5,
  height = 6.5
)


############################################################
## 24. M1 AND M2 COUNTS
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
    !is.na(
      M1_M2_class
    )
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
    "M1_M2_cell_counts_v2.csv"
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
    values = m1_m2_palette,
    name = "Cell type"
  ) +

  labs(
    title =
      "M1, M2 Mac cell counts",

    x = NULL,

    y = "Cell count"
  ) +

  theme_ueno_bar()

save_plot(
  p_counts,
  file.path(
    dir_count,
    "M1_M2_cell_counts_v2"
  ),
  width = 6.5,
  height = 6.5
)


############################################################
## 25. M2/M1 RATIO
############################################################

ratio_df <- cell_class_df %>%
  mutate(
    class_short = case_when(
      M1_M2_class ==
        "M1 cells" ~ "M1",

      M1_M2_class ==
        "M2 cells" ~ "M2"
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
      log2(
        M2_M1_ratio
      )
  )

data.table::fwrite(
  ratio_df,
  file.path(
    dir_table,
    "M2_M1_ratio_v2.csv"
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
    colour = "#888888"
  ) +

  geom_col(
    width = 0.55,
    colour = "#333333"
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
    title =
      "M2/M1 cell ratio in total Mac",

    x = NULL,

    y = "Log2 ratio"
  ) +

  theme_ueno_bar() +

  theme(
    legend.position = "none"
  )

save_plot(
  p_ratio,
  file.path(
    dir_ratio,
    "M2_M1_log2_ratio_v2"
  ),
  width = 6.5,
  height = 6.5
)


############################################################
## 26. PROFIBROTIC MACROPHAGE
############################################################

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
    "Profibrotic_Mphi_percentage_v2.csv"
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
    colour = "#333333"
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
    values = condition_palette
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.18
      )
    )
  ) +

  labs(
    title =
      "Profibrotic Mac % (/total M\u03c6)",

    x = NULL,

    y =
      "Profibrotic macrophage / total Mphi (%)"
  ) +

  theme_ueno_bar() +

  theme(
    legend.position = "none"
  )

save_plot(
  p_profibrotic,
  file.path(
    dir_ratio,
    "Profibrotic_Mphi_percentage_v2"
  ),
  width = 6.5,
  height = 6.5
)


############################################################
## 27. COMBINED FIGURE
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
        palette_name,
        " | v2.0"
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
    "Mouse_Mphi_STD_vs_CDAHFD_combined_figure_v2"
  ),
  width = 22,
  height = 13
)


############################################################
## 28. EXCEL OUTPUT
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

settings_table <- data.frame(
  Item = c(
    "RDS file",
    "Condition column",
    "Sample column",
    "Layer1 column",
    "Layer2 column",
    "Reduction",
    "Cells retained",
    "Palette",
    "Point size",
    "Analysis date"
  ),

  Value = c(
    rds_file,
    condition_col,
    sample_col,
    layer1_col,
    layer2_col,
    reduction_use,
    ncol(obj_sub),
    palette_name,
    umap_point_size,
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
  "Settings"
)

writeData(
  wb,
  "Settings",
  settings_table
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
  "Annotation_counts"
)

writeData(
  wb,
  "Annotation_counts",
  annotation_counts
)

addWorksheet(
  wb,
  "Annotation_mapping"
)

writeData(
  wb,
  "Annotation_mapping",
  mapping_check
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
  "Profibrotic"
)

writeData(
  wb,
  "Profibrotic",
  profibrotic_df
)

saveWorkbook(
  wb,
  file.path(
    dir_table,
    "Mouse_Mphi_STD_vs_CDAHFD_summary_v2.xlsx"
  ),
  overwrite = TRUE
)


############################################################
## 29. SAVE PALETTE
############################################################

saveRDS(
  ueno_mash_macro_palette_v1,
  file.path(
    dir_table,
    "Ueno_MASH_Macro_Palette_v1.rds"
  )
)


############################################################
## 30. SESSION INFO
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
  "RDS:\n",
  rds_file,
  "\n\n"
)

cat(
  "Layer1 column:\n",
  layer1_col,
  "\n\n"
)

cat(
  "Layer2 column:\n",
  layer2_col,
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
## 31. FINISH
############################################################

write_log(
  "Analysis completed successfully"
)

cat("\n")
cat("====================================================\n")
cat("Mouse Mphi analysis v2.0 completed\n")
cat("====================================================\n")
cat("Layer1 column :", layer1_col, "\n")
cat("Layer2 column :", layer2_col, "\n")
cat("Reduction     :", reduction_use, "\n")
cat("Cells         :", ncol(obj_sub), "\n")
cat("Output        :", output_dir, "\n")
cat("====================================================\n")
