############################################################
## 02_validate_RDS3_markers_v3.R
##
## RDS3 全クラスタ・マーカー遺伝子抽出
##
## 主な機能
##   1. RStudio Projectルートを自動認識
##   2. RDS3を自動検索
##   3. presto対応 Wilcoxon検定
##   4. クラスタ単位で逐次保存
##   5. 中断後の再開
##   6. FULL / FAST3000切り替え
##   7. Top10 / Top20 / Top30作成
##   8. Excel出力
##   9. ログ・解析条件・sessionInfo保存
############################################################


############################################################
## 0. ユーザー設定
############################################################

## 解析モード
##
## "FULL"
##   全細胞を使用する本解析
##
## "FAST3000"
##   各identity最大3,000細胞を用いる動作確認
##
analysis_mode <- "FULL"


## マーカー抽出条件
min_pct          <- 0.10
logfc_threshold  <- 0.25
only_positive    <- TRUE
random_seed      <- 1234


## 強制的に最初から再解析する場合だけ TRUE
##
## FALSE:
##   保存済みクラスタを飛ばして再開
##
## TRUE:
##   既存のクラスタ別CSVを削除して最初から実行
##
overwrite_existing <- FALSE


## 使用するアッセイ
assay_use <- "RNA"


## RDS3の正確なファイル名
rds3_filename <-
  "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"


############################################################
## 1. 必要パッケージ
############################################################

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "presto",
  "dplyr",
  "data.table",
  "openxlsx"
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
      paste(missing_packages, collapse = ", "),
      "\n\n必要なパッケージをインストールしてから再実行してください。"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(presto)
  library(dplyr)
  library(data.table)
  library(openxlsx)
})


############################################################
## 2. 基本設定
############################################################

set.seed(random_seed)

options(
  future.globals.maxSize = 64 * 1024^3,
  timeout = max(600, getOption("timeout"))
)

## data.tableが使用できるCPUスレッド数を設定
available_cores <- parallel::detectCores(logical = TRUE)

if (is.na(available_cores) || available_cores < 1) {
  available_cores <- 1
}

data.table::setDTthreads(
  threads = max(1, available_cores - 1)
)


############################################################
## 3. プロジェクトルートを認識
############################################################

find_project_root <- function(start_dir = getwd()) {

  current_dir <- normalizePath(
    start_dir,
    winslash = "/",
    mustWork = TRUE
  )

  repeat {

    rproj_files <- list.files(
      current_dir,
      pattern = "\\.Rproj$",
      full.names = TRUE
    )

    if (length(rproj_files) > 0) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)

    if (identical(parent_dir, current_dir)) {
      break
    }

    current_dir <- parent_dir
  }

  ## .Rprojが見つからない場合は現在の作業ディレクトリ
  normalizePath(
    start_dir,
    winslash = "/",
    mustWork = TRUE
  )
}


project_root <- find_project_root()

cat("\n")
cat("====================================================\n")
cat("RDS3 marker analysis v3\n")
cat("====================================================\n")
cat("Project root :", project_root, "\n")
cat("Working dir  :", getwd(), "\n")
cat("Mode         :", analysis_mode, "\n")
cat("presto       :", as.character(packageVersion("presto")), "\n")
cat("Seurat       :", as.character(packageVersion("Seurat")), "\n")
cat("CPU detected :", available_cores, "\n")
cat("====================================================\n\n")


############################################################
## 4. RDS3を自動検索
############################################################

preferred_paths <- c(
  file.path(project_root, rds3_filename),
  file.path(project_root, "data", rds3_filename),
  file.path(project_root, "Data", rds3_filename),
  file.path(project_root, "results", rds3_filename),
  file.path(project_root, "Results", rds3_filename),
  file.path(project_root, "rds", rds3_filename),
  file.path(project_root, "RDS", rds3_filename),
  file.path(project_root, "objects", rds3_filename)
)

preferred_paths <- preferred_paths[file.exists(preferred_paths)]

if (length(preferred_paths) > 0) {

  rds_file <- preferred_paths[1]

} else {

  cat("RDS3をプロジェクト内から検索しています...\n")

  rds_candidates <- list.files(
    path = project_root,
    pattern = paste0(
      "^",
      gsub("\\.", "\\\\.", rds3_filename),
      "$"
    ),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )

  ## 出力フォルダ内のコピーは除外
  rds_candidates <- rds_candidates[
    !grepl(
      "Phase2_Markers|AnnotationValidation|validation",
      rds_candidates,
      ignore.case = TRUE
    )
  ]

  if (length(rds_candidates) == 0) {

    stop(
      paste0(
        "\nRDS3が見つかりませんでした。\n",
        "検索したファイル名:\n",
        rds3_filename,
        "\n\n",
        "スクリプト上部の rds3_filename または",
        " rds_file の設定を確認してください。"
      ),
      call. = FALSE
    )
  }

  if (length(rds_candidates) > 1) {

    cat("\n複数のRDS3候補が見つかりました:\n")

    for (i in seq_along(rds_candidates)) {
      cat(sprintf("[%d] %s\n", i, rds_candidates[i]))
    }

    stop(
      paste0(
        "\nRDS3候補が複数あります。",
        "\n不要なコピーを移動するか、",
        "スクリプト内で rds_file を直接指定してください。"
      ),
      call. = FALSE
    )
  }

  rds_file <- rds_candidates[1]
}

rds_file <- normalizePath(
  rds_file,
  winslash = "/",
  mustWork = TRUE
)

cat("RDS3 file:\n", rds_file, "\n\n")


############################################################
## 5. 出力フォルダ
############################################################

mode_suffix <- switch(
  toupper(analysis_mode),
  "FULL"     = "FULL",
  "FAST3000" = "FAST3000",
  stop(
    "analysis_mode は 'FULL' または 'FAST3000' にしてください。",
    call. = FALSE
  )
)

outdir <- file.path(
  project_root,
  "results",
  "Mouse_MASH_RDS3_validation",
  paste0("Phase2_Markers_v3_", mode_suffix)
)

per_cluster_dir <- file.path(outdir, "PerCluster")
log_dir         <- file.path(outdir, "Logs")
checkpoint_dir  <- file.path(outdir, "Checkpoints")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(per_cluster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

progress_log <- file.path(
  log_dir,
  "FindMarkers_progress.log"
)

error_log <- file.path(
  log_dir,
  "FindMarkers_errors.log"
)


############################################################
## 6. 補助関数
############################################################

write_log <- function(text, file = progress_log) {

  message_text <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    text
  )

  cat(message_text, "\n")

  cat(
    message_text,
    "\n",
    file = file,
    append = TRUE
  )
}


safe_filename <- function(x) {

  x <- as.character(x)

  x <- gsub(
    "[^[:alnum:]_.-]+",
    "_",
    x
  )

  x <- gsub(
    "_+",
    "_",
    x
  )

  x
}


natural_cluster_order <- function(x) {

  x_character <- as.character(x)
  x_numeric   <- suppressWarnings(as.numeric(x_character))

  if (all(!is.na(x_numeric))) {
    return(x_character[order(x_numeric)])
  }

  sort(x_character)
}


get_fc_column <- function(df) {

  fc_candidates <- c(
    "avg_log2FC",
    "avg_logFC",
    "logFC"
  )

  found <- fc_candidates[
    fc_candidates %in% colnames(df)
  ]

  if (length(found) == 0) {
    return(NA_character_)
  }

  found[1]
}


############################################################
## 7. 既存結果の処理
############################################################

if (isTRUE(overwrite_existing)) {

  existing_csv <- list.files(
    per_cluster_dir,
    pattern = "\\.csv$",
    full.names = TRUE
  )

  existing_rds <- list.files(
    checkpoint_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )

  if (length(existing_csv) > 0) {
    file.remove(existing_csv)
  }

  if (length(existing_rds) > 0) {
    file.remove(existing_rds)
  }

  if (file.exists(progress_log)) {
    file.remove(progress_log)
  }

  if (file.exists(error_log)) {
    file.remove(error_log)
  }

  cat("既存結果を削除し、最初から解析します。\n")
}


############################################################
## 8. RDS3読み込み
############################################################

write_log("Loading RDS3")

load_start <- Sys.time()

obj <- readRDS(rds_file)

load_end <- Sys.time()

write_log(
  paste0(
    "RDS3 loaded in ",
    round(
      as.numeric(
        difftime(load_end, load_start, units = "mins")
      ),
      2
    ),
    " minutes"
  )
)

if (!inherits(obj, "Seurat")) {
  stop(
    "読み込んだオブジェクトはSeurat objectではありません。",
    call. = FALSE
  )
}

write_log(
  paste0(
    "Cells = ",
    format(ncol(obj), big.mark = ","),
    "; Features = ",
    format(nrow(obj), big.mark = ",")
  )
)


############################################################
## 9. RNA assayとLayer確認
############################################################

if (!assay_use %in% Assays(obj)) {

  stop(
    paste0(
      "指定した assay がありません: ",
      assay_use,
      "\nAvailable assays: ",
      paste(Assays(obj), collapse = ", ")
    ),
    call. = FALSE
  )
}

DefaultAssay(obj) <- assay_use

write_log(
  paste0(
    "DefaultAssay set to ",
    assay_use
  )
)


## Seurat v5の複数Layerを必要に応じて結合
layer_names <- tryCatch(
  SeuratObject::Layers(obj[[assay_use]]),
  error = function(e) character(0)
)

if (length(layer_names) > 0) {

  write_log(
    paste0(
      "RNA layers before JoinLayers: ",
      paste(layer_names, collapse = ", ")
    )
  )

  multiple_count_layers <- sum(
    grepl("^counts\\.", layer_names)
  ) > 1

  multiple_data_layers <- sum(
    grepl("^data\\.", layer_names)
  ) > 1

  if (multiple_count_layers || multiple_data_layers) {

    write_log("Joining RNA assay layers")

    obj <- tryCatch(
      {
        SeuratObject::JoinLayers(
          object = obj,
          assay = assay_use
        )
      },
      error = function(e1) {

        tryCatch(
          {
            Seurat::JoinLayers(
              object = obj,
              assay = assay_use
            )
          },
          error = function(e2) {

            stop(
              paste0(
                "JoinLayersに失敗しました。\n",
                "SeuratObject error: ",
                conditionMessage(e1),
                "\nSeurat error: ",
                conditionMessage(e2)
              ),
              call. = FALSE
            )
          }
        )
      }
    )

    layer_names_after <- tryCatch(
      SeuratObject::Layers(obj[[assay_use]]),
      error = function(e) character(0)
    )

    write_log(
      paste0(
        "RNA layers after JoinLayers: ",
        paste(layer_names_after, collapse = ", ")
      )
    )
  }
}


############################################################
## 10. クラスタ列を自動判定
############################################################

cluster_candidates <- c(
  "cluster_for_R8plot_FIXED2",
  "cluster_for_R8plot",
  "seurat_clusters",
  "RNA_snn_res.3",
  "integrated_snn_res.3"
)

cluster_col <- cluster_candidates[
  cluster_candidates %in% colnames(obj@meta.data)
][1]

if (is.na(cluster_col) || length(cluster_col) == 0) {

  stop(
    paste0(
      "クラスタ列を検出できませんでした。\n\n",
      "Metadata columns:\n",
      paste(
        colnames(obj@meta.data),
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}

cluster_values <- as.character(
  obj@meta.data[[cluster_col]]
)

if (anyNA(cluster_values)) {

  stop(
    paste0(
      "クラスタ列 ",
      cluster_col,
      " にNAがあります。"
    ),
    call. = FALSE
  )
}

cluster_levels <- natural_cluster_order(
  unique(cluster_values)
)

obj@meta.data[[cluster_col]] <- factor(
  cluster_values,
  levels = cluster_levels
)

Idents(obj) <- cluster_col

write_log(
  paste0(
    "Cluster column = ",
    cluster_col
  )
)

write_log(
  paste0(
    "Number of clusters = ",
    length(cluster_levels)
  )
)


############################################################
## 11. クラスタ細胞数
############################################################

cluster_cell_counts <- data.frame(
  cluster = names(table(Idents(obj))),
  n_cells = as.integer(table(Idents(obj))),
  stringsAsFactors = FALSE
)

cluster_cell_counts <- cluster_cell_counts[
  match(
    cluster_levels,
    cluster_cell_counts$cluster
  ),
  ,
  drop = FALSE
]

data.table::fwrite(
  cluster_cell_counts,
  file.path(
    outdir,
    "ClusterCellCounts.csv"
  )
)

print(cluster_cell_counts)


############################################################
## 12. 解析条件保存
############################################################

if (toupper(analysis_mode) == "FULL") {

  max_cells_per_ident <- Inf

} else {

  max_cells_per_ident <- 3000
}

analysis_settings <- data.frame(
  parameter = c(
    "project_root",
    "rds_file",
    "analysis_mode",
    "assay",
    "cluster_column",
    "number_of_cells",
    "number_of_features",
    "number_of_clusters",
    "test_use",
    "presto_version",
    "min_pct",
    "logfc_threshold",
    "only_positive",
    "max_cells_per_ident",
    "random_seed",
    "analysis_start"
  ),
  value = c(
    project_root,
    rds_file,
    analysis_mode,
    assay_use,
    cluster_col,
    ncol(obj),
    nrow(obj),
    length(cluster_levels),
    "wilcox",
    as.character(packageVersion("presto")),
    min_pct,
    logfc_threshold,
    only_positive,
    as.character(max_cells_per_ident),
    random_seed,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  analysis_settings,
  file.path(
    outdir,
    "AnalysisSettings.csv"
  )
)


############################################################
## 13. 完了済みクラスタを判定
############################################################

cluster_manifest <- data.frame(
  cluster = cluster_levels,
  safe_cluster = vapply(
    cluster_levels,
    safe_filename,
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)

cluster_manifest$output_file <- file.path(
  per_cluster_dir,
  paste0(
    "Cluster_",
    cluster_manifest$safe_cluster,
    ".csv"
  )
)

data.table::fwrite(
  cluster_manifest,
  file.path(
    outdir,
    "ClusterFileManifest.csv"
  )
)


is_completed <- vapply(
  cluster_manifest$output_file,
  function(x) {

    if (!file.exists(x)) {
      return(FALSE)
    }

    file_info <- file.info(x)

    if (is.na(file_info$size) || file_info$size == 0) {
      return(FALSE)
    }

    test_read <- tryCatch(
      data.table::fread(
        x,
        nrows = 2,
        showProgress = FALSE
      ),
      error = function(e) NULL
    )

    !is.null(test_read)
  },
  FUN.VALUE = logical(1)
)

completed_clusters <- cluster_manifest$cluster[
  is_completed
]

todo_clusters <- cluster_manifest$cluster[
  !is_completed
]

write_log(
  paste0(
    "Completed clusters detected = ",
    length(completed_clusters)
  )
)

write_log(
  paste0(
    "Clusters remaining = ",
    length(todo_clusters)
  )
)

if (length(completed_clusters) > 0) {

  cat(
    "\n完了済みクラスタ:\n",
    paste(completed_clusters, collapse = ", "),
    "\n"
  )
}

if (length(todo_clusters) > 0) {

  cat(
    "\n今回解析するクラスタ:\n",
    paste(todo_clusters, collapse = ", "),
    "\n\n"
  )
}


############################################################
## 14. クラスタ単位 FindMarkers
############################################################

if (length(todo_clusters) > 0) {

  for (i in seq_along(todo_clusters)) {

    current_cluster <- todo_clusters[i]

    safe_cluster <- cluster_manifest$safe_cluster[
      cluster_manifest$cluster == current_cluster
    ]

    output_csv <- cluster_manifest$output_file[
      cluster_manifest$cluster == current_cluster
    ]

    checkpoint_rds <- file.path(
      checkpoint_dir,
      paste0(
        "Cluster_",
        safe_cluster,
        ".rds"
      )
    )

    cluster_n <- sum(
      as.character(Idents(obj)) == current_cluster
    )

    write_log(
      paste0(
        "START cluster ",
        current_cluster,
        " [",
        i,
        "/",
        length(todo_clusters),
        "] cells=",
        format(cluster_n, big.mark = ",")
      )
    )

    cluster_start <- Sys.time()

    marker_result <- tryCatch(

      {
        Seurat::FindMarkers(
          object = obj,
          ident.1 = current_cluster,
          assay = assay_use,
          slot = "data",
          test.use = "wilcox",
          only.pos = only_positive,
          min.pct = min_pct,
          logfc.threshold = logfc_threshold,
          max.cells.per.ident = max_cells_per_ident,
          random.seed = random_seed,
          densify = FALSE,
          verbose = FALSE
        )
      },

      error = function(e) {

        error_message <- paste0(
          "ERROR cluster ",
          current_cluster,
          ": ",
          conditionMessage(e)
        )

        write_log(
          error_message,
          file = error_log
        )

        NULL
      }
    )

    if (is.null(marker_result)) {

      write_log(
        paste0(
          "Cluster ",
          current_cluster,
          " failed; proceeding to next cluster"
        )
      )

      next
    }

    marker_result <- as.data.frame(
      marker_result
    )

    marker_result$gene <- rownames(
      marker_result
    )

    rownames(marker_result) <- NULL

    marker_result$cluster <- current_cluster
    marker_result$n_cells_cluster <- cluster_n
    marker_result$analysis_mode <- analysis_mode
    marker_result$assay <- assay_use
    marker_result$test <- "wilcox_presto"

    preferred_first_columns <- c(
      "cluster",
      "gene",
      "n_cells_cluster",
      "analysis_mode",
      "assay",
      "test"
    )

    marker_result <- marker_result[
      ,
      c(
        preferred_first_columns,
        setdiff(
          colnames(marker_result),
          preferred_first_columns
        )
      ),
      drop = FALSE
    ]

    ## 一時RDS保存
    saveRDS(
      marker_result,
      checkpoint_rds,
      compress = FALSE
    )

    ## CSVは一時ファイルを経由して安全に保存
    temporary_csv <- paste0(
      output_csv,
      ".tmp"
    )

    data.table::fwrite(
      marker_result,
      temporary_csv
    )

    if (file.exists(output_csv)) {
      file.remove(output_csv)
    }

    rename_success <- file.rename(
      temporary_csv,
      output_csv
    )

    if (!rename_success) {

      stop(
        paste0(
          "クラスタ ",
          current_cluster,
          " のCSV保存に失敗しました:\n",
          output_csv
        ),
        call. = FALSE
      )
    }

    cluster_end <- Sys.time()

    elapsed_minutes <- round(
      as.numeric(
        difftime(
          cluster_end,
          cluster_start,
          units = "mins"
        )
      ),
      2
    )

    write_log(
      paste0(
        "FINISH cluster ",
        current_cluster,
        "; markers=",
        nrow(marker_result),
        "; elapsed=",
        elapsed_minutes,
        " min"
      )
    )

    rm(marker_result)
    invisible(gc(verbose = FALSE))
  }
}


############################################################
## 15. クラスタ別ファイルの完全性確認
############################################################

completed_after <- vapply(
  cluster_manifest$output_file,
  file.exists,
  FUN.VALUE = logical(1)
)

missing_after <- cluster_manifest$cluster[
  !completed_after
]

if (length(missing_after) > 0) {

  warning(
    paste0(
      "以下のクラスタは未完了です:\n",
      paste(missing_after, collapse = ", "),
      "\n\n",
      "エラーログを確認し、同じスクリプトを再実行してください:\n",
      error_log
    ),
    call. = FALSE
  )

} else {

  write_log("All clusters completed")
}


############################################################
## 16. クラスタ別結果を統合
############################################################

available_files <- cluster_manifest$output_file[
  file.exists(cluster_manifest$output_file)
]

if (length(available_files) == 0) {

  stop(
    "統合できるクラスタ別マーカーファイルがありません。",
    call. = FALSE
  )
}

write_log(
  paste0(
    "Merging ",
    length(available_files),
    " cluster files"
  )
)

all_markers_list <- lapply(
  available_files,
  function(x) {
    data.table::fread(
      x,
      showProgress = FALSE
    )
  }
)

all_markers <- data.table::rbindlist(
  all_markers_list,
  use.names = TRUE,
  fill = TRUE
)

all_markers <- as.data.frame(
  all_markers
)

all_markers$cluster <- factor(
  as.character(all_markers$cluster),
  levels = cluster_levels
)

fc_column <- get_fc_column(
  all_markers
)

if (is.na(fc_column)) {

  stop(
    paste0(
      "Fold-change列が見つかりません。\n",
      "Available columns:\n",
      paste(
        colnames(all_markers),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

write_log(
  paste0(
    "Fold-change column = ",
    fc_column
  )
)


############################################################
## 17. 並び替え
############################################################

all_markers <- all_markers %>%
  arrange(
    cluster,
    desc(.data[[fc_column]])
  )

all_markers$cluster <- as.character(
  all_markers$cluster
)

data.table::fwrite(
  all_markers,
  file.path(
    outdir,
    "AllClusterMarkers.csv"
  )
)

saveRDS(
  all_markers,
  file.path(
    outdir,
    "AllClusterMarkers.rds"
  ),
  compress = FALSE
)


############################################################
## 18. Top10 / Top20 / Top30
############################################################

make_top_markers <- function(marker_df, n_top, fc_col) {

  marker_df %>%
    group_by(cluster) %>%
    arrange(
      desc(.data[[fc_col]]),
      .by_group = TRUE
    ) %>%
    slice_head(n = n_top) %>%
    ungroup()
}

top10 <- make_top_markers(
  all_markers,
  n_top = 10,
  fc_col = fc_column
)

top20 <- make_top_markers(
  all_markers,
  n_top = 20,
  fc_col = fc_column
)

top30 <- make_top_markers(
  all_markers,
  n_top = 30,
  fc_col = fc_column
)

data.table::fwrite(
  top10,
  file.path(
    outdir,
    "Top10MarkersPerCluster.csv"
  )
)

data.table::fwrite(
  top20,
  file.path(
    outdir,
    "Top20MarkersPerCluster.csv"
  )
)

data.table::fwrite(
  top30,
  file.path(
    outdir,
    "Top30MarkersPerCluster.csv"
  )
)


############################################################
## 19. 横型Top marker表
############################################################

make_marker_wide <- function(top_df, n_top) {

  top_df %>%
    group_by(cluster) %>%
    mutate(rank = row_number()) %>%
    select(cluster, rank, gene) %>%
    tidyr::pivot_wider(
      names_from = rank,
      values_from = gene,
      names_prefix = "Rank_"
    ) %>%
    ungroup()
}

## tidyrが利用できる場合のみ作成
if (requireNamespace("tidyr", quietly = TRUE)) {

  top10_wide <- make_marker_wide(
    top10,
    10
  )

  top20_wide <- make_marker_wide(
    top20,
    20
  )

  top30_wide <- make_marker_wide(
    top30,
    30
  )

  data.table::fwrite(
    top10_wide,
    file.path(
      outdir,
      "Top10Markers_Wide.csv"
    )
  )

  data.table::fwrite(
    top20_wide,
    file.path(
      outdir,
      "Top20Markers_Wide.csv"
    )
  )

  data.table::fwrite(
    top30_wide,
    file.path(
      outdir,
      "Top30Markers_Wide.csv"
    )
  )

} else {

  top10_wide <- NULL
  top20_wide <- NULL
  top30_wide <- NULL

  write_log(
    "tidyr not installed: wide marker tables were skipped"
  )
}


############################################################
## 20. Excel workbook
############################################################

write_log("Creating Excel workbook")

excel_file <- file.path(
  outdir,
  paste0(
    "RDS3_ClusterMarkers_v3_",
    mode_suffix,
    ".xlsx"
  )
)

wb <- openxlsx::createWorkbook()

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom",
  wrapText = TRUE
)

add_excel_sheet <- function(
    workbook,
    sheet_name,
    data,
    freeze_first_row = TRUE
) {

  openxlsx::addWorksheet(
    workbook,
    sheet_name
  )

  openxlsx::writeData(
    workbook,
    sheet = sheet_name,
    x = data,
    headerStyle = header_style,
    withFilter = TRUE
  )

  if (freeze_first_row) {

    openxlsx::freezePane(
      workbook,
      sheet = sheet_name,
      firstRow = TRUE
    )
  }

  openxlsx::setColWidths(
    workbook,
    sheet = sheet_name,
    cols = seq_len(ncol(data)),
    widths = "auto"
  )
}


add_excel_sheet(
  wb,
  "AnalysisSettings",
  analysis_settings
)

add_excel_sheet(
  wb,
  "ClusterCellCounts",
  cluster_cell_counts
)

add_excel_sheet(
  wb,
  "Top10",
  as.data.frame(top10)
)

add_excel_sheet(
  wb,
  "Top20",
  as.data.frame(top20)
)

add_excel_sheet(
  wb,
  "Top30",
  as.data.frame(top30)
)

if (!is.null(top10_wide)) {

  add_excel_sheet(
    wb,
    "Top10_Wide",
    as.data.frame(top10_wide)
  )

  add_excel_sheet(
    wb,
    "Top20_Wide",
    as.data.frame(top20_wide)
  )

  add_excel_sheet(
    wb,
    "Top30_Wide",
    as.data.frame(top30_wide)
  )
}


## Excelの1シート上限を考慮
excel_row_limit <- 1048576

if (nrow(all_markers) + 1 < excel_row_limit) {

  add_excel_sheet(
    wb,
    "AllMarkers",
    all_markers
  )

} else {

  all_marker_note <- data.frame(
    message = paste0(
      "AllMarkersはExcelの行数上限を超えるため、",
      "CSVのみ保存しました。総行数: ",
      format(nrow(all_markers), big.mark = ",")
    ),
    csv_file = file.path(
      outdir,
      "AllClusterMarkers.csv"
    ),
    stringsAsFactors = FALSE
  )

  add_excel_sheet(
    wb,
    "AllMarkers_INFO",
    all_marker_note
  )
}


openxlsx::saveWorkbook(
  wb,
  excel_file,
  overwrite = TRUE
)


############################################################
## 21. sessionInfo保存
############################################################

session_info_file <- file.path(
  log_dir,
  "sessionInfo.txt"
)

sink(session_info_file)

cat("Analysis completed:\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("RDS file:\n")
cat(rds_file, "\n\n")

cat("Output directory:\n")
cat(outdir, "\n\n")

cat("presto version:\n")
print(packageVersion("presto"))

cat("\nSeurat version:\n")
print(packageVersion("Seurat"))

cat("\nSession information:\n")
print(sessionInfo())

sink()


############################################################
## 22. 完了メッセージ
############################################################

write_log(
  paste0(
    "ANALYSIS COMPLETE; total marker rows=",
    format(nrow(all_markers), big.mark = ",")
  )
)

cat("\n")
cat("====================================================\n")
cat("RDS3 marker analysis v3 completed\n")
cat("====================================================\n")
cat("Mode              :", analysis_mode, "\n")
cat("RDS3              :", rds_file, "\n")
cat("Cluster column    :", cluster_col, "\n")
cat("Number of clusters:", length(cluster_levels), "\n")
cat("Total marker rows :", format(nrow(all_markers), big.mark = ","), "\n")
cat("FC column         :", fc_column, "\n")
cat("Output directory  :", outdir, "\n")
cat("Excel file        :", excel_file, "\n")
cat("Progress log      :", progress_log, "\n")
cat("Error log         :", error_log, "\n")
cat("====================================================\n")
