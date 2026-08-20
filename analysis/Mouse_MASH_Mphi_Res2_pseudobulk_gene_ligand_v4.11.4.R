#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH scRNA-seq
# MΦ-only RPCA Res2.0
# SAFE pseudobulk rebuild + gene validation + ligand output
# Fixed manual annotation v4.8.4
# v4.11.4
# ==============================================================================
#
# CORE PRINCIPLE
#   Pseudobulk is rebuilt from scratch using the manually validated method:
#
#     cells <- rownames(meta)[sample == X & subtype == Y]
#     pb    <- Matrix::rowSums(counts[, cells, drop = FALSE])
#
#   No indirect aggregation helper is used.
#
# ANALYSES
#   1) Safe pseudobulk matrix for 5 fixed MΦ subtypes x biological samples
#   2) QC of each pseudobulk: n_cells, library size, marker counts / CPM
#   3) Descriptive STD vs CDAHFD log2FC
#   4) Descriptive Sham vs Tx log2FC
#   5) Optional exploratory edgeR for Sham vs Tx if edgeR is installed
#   6) Gene-level custom DotPlot from normalized RNA data
#   7) Sample-level pseudobulk heatmap
#   8) Candidate ligand-output DotPlot / heatmap / effect-size table
#
# IMPORTANT
#   - Res2.0 clustering unchanged
#   - v4.8.4 MΦ annotation unchanged
#   - STD and CDAHFD are descriptive only if each has n=1
#   - No claim of statistical significance from STD vs CDAHFD
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4114)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

PROJECT_DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
  PROJECT_DATA_ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_manual_annotation_v4.8.4",
  "RDS",
  "Mouse_Mphi_Res2_manual_class_annotated_v4.8.4.rds"
)

OUTPUT_DIR <- file.path(
  PROJECT_DATA_ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_pseudobulk_gene_ligand_v4.11.4"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
PB_DIR  <- file.path(OUTPUT_DIR, "Pseudobulk")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PB_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"

PREFERRED_SAMPLE_COL    <- "sample_4group"
PREFERRED_CONDITION_COL <- "condition_4group"
PREFERRED_CLASS_COL     <- "macrophage_class_Res2_v484"

CONDITION_ORDER <- c("STD", "CDAHFD", "Sham", "Tx")

MAJOR_CLASS_ORDER <- c(
  "Inflammatory-Mphi",
  "Anti-inflammatory-Mphi",
  "Fibrogenic-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

CLASS_LABELS <- c(
  "Inflammatory-Mphi"           = "Inflammatory-MΦ",
  "Anti-inflammatory-Mphi"      = "Anti-inflammatory-MΦ",
  "Fibrogenic-Mphi"             = "Fibrogenic-MΦ",
  "Repair/Resolution-Mphi"      = "Repair/Resolution-MΦ",
  "Lipid-associated/TREM2-Mphi" = "Lipid-associated/TREM2-MΦ"
)

CPM_PSEUDOCOUNT <- 1
MIN_CPM <- 1
MIN_SAMPLES_EXPRESSED <- 2L
TOP_N_EFFECT <- 20L

# ------------------------------------------------------------------------------
# 1. Validation genes and candidate ligands
# ------------------------------------------------------------------------------

VALIDATION_GENE_SETS <- list(
  Inflammatory = c(
    "Il1b","Tnf","Ccl2","Ccl3","Ccl4","Cxcl10","Nos2","Cd80","Cd86","Stat1"
  ),
  Anti_inflammatory = c(
    "Mrc1","Cd163","Il1rn","Retnla","Chil3","Arg1","Mertk","Igf1","Hmox1","Klf4","Maf"
  ),
  Fibrogenic = c(
    "Spp1","Tgfb1","Pdgfb","Thbs1","Lgals3","Gpnmb","Mmp12","Mmp14","Ctsb"
  ),
  Repair_Resolution = c(
    "Mertk","Axl","Mfge8","Gas6","Igf1","Hmox1","Mmp12","Mmp13","Mmp14","Plau"
  ),
  Lipid_TREM2 = c(
    "Trem2","Gpnmb","Cd9","Lpl","Apoe","Fabp5","Abca1","Plin2","Ctsd"
  ),
  IL10_STAT3 = c(
    "Il10","Il10ra","Il10rb","Jak1","Tyk2","Stat3","Socs3","Bcl3","Il1rn"
  )
)

LIGAND_GENE_SETS <- list(
  HSC_candidate_ligands = c(
    "Tgfb1","Pdgfb","Spp1","Osm","Il1b","Tnf","Ccl2","Ccl3","Ccl4","Igf1"
  ),
  LSEC_vascular_candidate_ligands = c(
    "Vegfa","Angptl4","Tgfb1","Spp1","Osm","Il1b","Tnf","Igf1"
  ),
  Hepatocyte_candidate_ligands = c(
    "Il1b","Tnf","Osm","Il6","Tgfb1","Igf1","Spp1"
  ),
  Repair_resolution_candidate_ligands = c(
    "Igf1","Gas6","Mfge8","Il10","Tgfb1","Areg"
  )
)

QC_GENES <- c(
  "Adgre1","Lyz2","Csf1r","Il1b","Tnf",
  "Mrc1","Cd163","Spp1","Trem2","Apoe"
)

# ------------------------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "Seurat","SeuratObject","Matrix",
  "dplyr","tidyr","tibble",
  "ggplot2","patchwork","scales","pheatmap"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

HAS_EDGER <- requireNamespace("edgeR", quietly = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(pheatmap)
})

if (HAS_EDGER) {
  suppressPackageStartupMessages(library(edgeR))
}

# ------------------------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------------------------

msg <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (length(hit) == 0L) return(NA_character_)
  hit[[1]]
}

canonical_condition <- function(condition, sample) {
  out <- as.character(condition)
  smp <- as.character(sample)

  out[grepl("^STD$", out, ignore.case = TRUE)] <- "STD"
  out[grepl("CDAHFD|CDHFD", out, ignore.case = TRUE)] <- "CDAHFD"
  out[grepl("^Sham", out, ignore.case = TRUE)] <- "Sham"
  out[grepl("^Tx", out, ignore.case = TRUE)] <- "Tx"

  out[!(out %in% CONDITION_ORDER)] <- NA_character_

  out[is.na(out) & grepl("^STD", smp, ignore.case = TRUE)] <- "STD"
  out[is.na(out) & grepl("CDAHFD|CDHFD", smp, ignore.case = TRUE)] <- "CDAHFD"
  out[is.na(out) & grepl("^Sham", smp, ignore.case = TRUE)] <- "Sham"
  out[is.na(out) & grepl("^Tx", smp, ignore.case = TRUE)] <- "Tx"

  factor(out, levels = CONDITION_ORDER)
}

normalize_class <- function(x) {
  y <- as.character(x)

  y[grepl("^Inflammatory", y, ignore.case = TRUE)] <- "Inflammatory-Mphi"
  y[grepl("Anti[-_ ]?inflammatory|Restrict[-_ ]?M2", y, ignore.case = TRUE)] <-
    "Anti-inflammatory-Mphi"
  y[grepl("^Fibrogenic", y, ignore.case = TRUE)] <- "Fibrogenic-Mphi"
  y[grepl("Repair.?/?Resolution", y, ignore.case = TRUE)] <-
    "Repair/Resolution-Mphi"
  y[grepl("Lipid.*TREM2", y, ignore.case = TRUE)] <-
    "Lipid-associated/TREM2-Mphi"

  y
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

save_pdf <- function(filename, plot, width, height) {
  ggsave(
    filename = file.path(FIG_DIR, filename),
    plot = plot,
    device = cairo_pdf,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
}

get_layer <- function(object, assay, layer) {
  x <- tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e) NULL
  )

  if (is.null(x)) {
    x <- tryCatch(
      GetAssayData(object, assay = assay, slot = layer),
      error = function(e) NULL
    )
  }

  x
}

calc_cpm <- function(pb_counts) {
  lib <- colSums(pb_counts)
  if (any(lib <= 0)) {
    stop("At least one pseudobulk library has zero total counts.")
  }

  sweep(
    pb_counts,
    2,
    lib / 1e6,
    "/"
  )
}

safe_row_z <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

make_custom_dotplot_data <- function(
  object,
  assay,
  group_vector,
  group_levels,
  features
) {
  features <- unique(intersect(features, rownames(object[[assay]])))

  expr <- get_layer(object, assay, "data")
  if (is.null(expr)) {
    stop("Normalized data layer was not found.")
  }

  expr <- expr[
    features,
    colnames(object),
    drop = FALSE
  ]

  group_vector <- factor(
    as.character(group_vector),
    levels = group_levels
  )

  res <- vector("list", length(group_levels))
  names(res) <- group_levels

  for (grp in group_levels) {
    idx <- which(group_vector == grp)
    if (length(idx) == 0L) next

    mat <- expr[, idx, drop = FALSE]

    res[[grp]] <- tibble(
      group = grp,
      gene = features,
      avg_expr = as.numeric(Matrix::rowMeans(mat)),
      pct_expr = as.numeric(Matrix::rowMeans(mat > 0)) * 100
    )
  }

  out <- bind_rows(res)

  out %>%
    group_by(gene) %>%
    mutate(
      avg_expr_scaled = {
        z <- as.numeric(scale(avg_expr))
        z[!is.finite(z)] <- 0
        pmax(pmin(z, 2.5), -2.5)
      }
    ) %>%
    ungroup() %>%
    mutate(
      group = factor(group, levels = group_levels),
      gene = factor(gene, levels = features)
    )
}

# ------------------------------------------------------------------------------
# 4. Load object and validate source layers
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
  stop("INPUT_RDS not found: ", INPUT_RDS)
}

msg("Loading: ", INPUT_RDS)
mphi <- readRDS(INPUT_RDS)

if (!inherits(mphi, "Seurat")) {
  stop("INPUT_RDS is not a Seurat object.")
}

if (!(ASSAY_USE %in% Assays(mphi))) {
  stop("RNA assay not found.")
}

DefaultAssay(mphi) <- ASSAY_USE

counts <- get_layer(mphi, ASSAY_USE, "counts")
data_norm <- get_layer(mphi, ASSAY_USE, "data")

if (is.null(counts)) stop("RNA counts layer not found.")
if (is.null(data_norm)) stop("RNA normalized data layer not found.")

if (!identical(colnames(counts), colnames(mphi))) {
  stop("Counts columns do not match Seurat cell order.")
}

if (!identical(colnames(mphi), rownames(mphi@meta.data))) {
  stop("Seurat cell order does not match metadata row order.")
}

msg(
  "Counts matrix: ",
  nrow(counts), " genes x ", ncol(counts), " cells"
)

# ------------------------------------------------------------------------------
# 5. Canonical metadata
# ------------------------------------------------------------------------------

meta_cols <- colnames(mphi@meta.data)

SAMPLE_COL <- first_existing(
  meta_cols,
  c(PREFERRED_SAMPLE_COL, "sample", "sample_id", "orig.ident")
)

CONDITION_COL <- first_existing(
  meta_cols,
  c(PREFERRED_CONDITION_COL, "condition", "group", "Condition")
)

CLASS_COL <- first_existing(
  meta_cols,
  c(
    PREFERRED_CLASS_COL,
    "macrophage_class_v484",
    "manual_class_v484",
    "macrophage_class",
    "Mphi_class"
  )
)

if (is.na(SAMPLE_COL)) stop("Sample column not found.")
if (is.na(CLASS_COL)) stop("v4.8.4 class column not found.")

if (is.na(CONDITION_COL)) {
  CONDITION_COL <- ".v4114_condition_placeholder"
  mphi@meta.data[[CONDITION_COL]] <- NA_character_
}

mphi$sample_v4114 <- as.character(mphi@meta.data[[SAMPLE_COL]])

mphi$condition_v4114 <- canonical_condition(
  mphi@meta.data[[CONDITION_COL]],
  mphi@meta.data[[SAMPLE_COL]]
)

mphi$class_v4114 <- normalize_class(
  mphi@meta.data[[CLASS_COL]]
)

if (anyNA(mphi$condition_v4114)) {
  stop(
    "Unresolved conditions: ",
    paste(
      unique(mphi$sample_v4114[is.na(mphi$condition_v4114)]),
      collapse = ", "
    )
  )
}

cell_meta <- mphi@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell = cell,
    sample = sample_v4114,
    condition = factor(condition_v4114, levels = CONDITION_ORDER),
    macrophage_class = class_v4114
  )

sample_condition <- cell_meta %>%
  distinct(sample, condition) %>%
  arrange(condition, sample)

write.csv(
  sample_condition,
  file.path(TAB_DIR, "00_sample_condition_map_v4.11.4.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. SAFE pseudobulk rebuild: explicit rowSums per sample x subtype
# ------------------------------------------------------------------------------

pb_meta <- cell_meta %>%
  filter(macrophage_class %in% MAJOR_CLASS_ORDER) %>%
  distinct(sample, condition, macrophage_class) %>%
  arrange(
    factor(macrophage_class, levels = MAJOR_CLASS_ORDER),
    condition,
    sample
  ) %>%
  mutate(
    pb_id = paste(macrophage_class, sample, sep = "__")
  )

pb_vectors <- vector("list", nrow(pb_meta))
names(pb_vectors) <- pb_meta$pb_id

qc_rows <- vector("list", nrow(pb_meta))

for (i in seq_len(nrow(pb_meta))) {

  smp <- pb_meta$sample[[i]]
  subtype <- pb_meta$macrophage_class[[i]]
  pb_id <- pb_meta$pb_id[[i]]

  cells <- rownames(mphi@meta.data)[
    mphi$sample_v4114 == smp &
      mphi$class_v4114 == subtype
  ]

  cells <- intersect(cells, colnames(counts))

  if (length(cells) == 0L) {
    stop(
      "No cells for pseudobulk: ",
      pb_id
    )
  }

  pb <- Matrix::rowSums(
    counts[, cells, drop = FALSE]
  )

  if (length(pb) != nrow(counts)) {
    stop(
      "Unexpected pseudobulk vector length for ",
      pb_id
    )
  }

  if (sum(pb) <= 0) {
    stop(
      "Zero-count pseudobulk library: ",
      pb_id
    )
  }

  pb_vectors[[pb_id]] <- pb

  qc_rows[[i]] <- tibble(
    pb_id = pb_id,
    sample = smp,
    condition = as.character(pb_meta$condition[[i]]),
    macrophage_class = subtype,
    n_cells = length(cells),
    library_size = sum(pb),
    mean_counts_per_cell = sum(pb) / length(cells)
  )
}

pb_counts <- do.call(
  cbind,
  pb_vectors
)

rownames(pb_counts) <- rownames(counts)
colnames(pb_counts) <- names(pb_vectors)

pb_qc <- bind_rows(qc_rows)

# ------------------------------------------------------------------------------
# 7. HARD QC: ensure pseudobulk columns are not accidentally identical
# ------------------------------------------------------------------------------

identical_pairs <- list()
pair_idx <- 1L

for (i in seq_len(ncol(pb_counts) - 1L)) {
  for (j in (i + 1L):ncol(pb_counts)) {

    if (identical(pb_counts[, i], pb_counts[, j])) {

      identical_pairs[[pair_idx]] <- tibble(
        pb1 = colnames(pb_counts)[[i]],
        pb2 = colnames(pb_counts)[[j]]
      )

      pair_idx <- pair_idx + 1L
    }
  }
}

identical_pairs_df <- bind_rows(identical_pairs)

if (nrow(identical_pairs_df) > 0L) {

  write.csv(
    identical_pairs_df,
    file.path(
      TAB_DIR,
      "ERROR_identical_pseudobulk_columns_v4.11.4.csv"
    ),
    row.names = FALSE
  )

  stop(
    "Identical pseudobulk columns were detected. ",
    "See ERROR_identical_pseudobulk_columns_v4.11.4.csv"
  )
}

msg("QC passed: no identical pseudobulk columns.")

# ------------------------------------------------------------------------------
# 8. CPM / log2CPM and QC marker table
# ------------------------------------------------------------------------------

pb_cpm <- calc_cpm(pb_counts)

pb_log2cpm <- log2(
  pb_cpm +
    CPM_PSEUDOCOUNT
)

qc_genes_use <- intersect(
  QC_GENES,
  rownames(pb_counts)
)

qc_gene_table <- as.data.frame(
  pb_counts[
    qc_genes_use,
    ,
    drop = FALSE
  ]
) %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "pb_id",
    values_to = "raw_count"
  ) %>%
  left_join(
    as.data.frame(
      pb_cpm[
        qc_genes_use,
        ,
        drop = FALSE
      ]
    ) %>%
      rownames_to_column("gene") %>%
      pivot_longer(
        cols = -gene,
        names_to = "pb_id",
        values_to = "CPM"
      ),
    by = c("gene", "pb_id")
  ) %>%
  left_join(
    pb_meta,
    by = "pb_id"
  )

write.csv(
  pb_qc,
  file.path(PB_DIR, "00_pseudobulk_QC_library_size_v4.11.4.csv"),
  row.names = FALSE
)

write.csv(
  qc_gene_table,
  file.path(PB_DIR, "01_pseudobulk_QC_marker_counts_CPM_v4.11.4.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = rownames(pb_counts),
    pb_counts,
    check.names = FALSE
  ),
  file.path(PB_DIR, "02_pseudobulk_raw_counts_v4.11.4.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = rownames(pb_cpm),
    pb_cpm,
    check.names = FALSE
  ),
  file.path(PB_DIR, "03_pseudobulk_CPM_v4.11.4.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = rownames(pb_log2cpm),
    pb_log2cpm,
    check.names = FALSE
  ),
  file.path(PB_DIR, "04_pseudobulk_log2CPM_v4.11.4.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Sample counts cross-check against metadata
# ------------------------------------------------------------------------------

metadata_crosscheck <- cell_meta %>%
  filter(macrophage_class %in% MAJOR_CLASS_ORDER) %>%
  count(
    sample,
    condition,
    macrophage_class,
    name = "metadata_n_cells"
  ) %>%
  mutate(
    pb_id = paste(macrophage_class, sample, sep = "__")
  ) %>%
  left_join(
    pb_qc %>%
      select(pb_id, n_cells, library_size),
    by = "pb_id"
  ) %>%
  mutate(
    cell_count_match = metadata_n_cells == n_cells
  )

write.csv(
  metadata_crosscheck,
  file.path(
    PB_DIR,
    "05_pseudobulk_metadata_cellcount_crosscheck_v4.11.4.csv"
  ),
  row.names = FALSE
)

if (!all(metadata_crosscheck$cell_count_match)) {
  stop("Pseudobulk cell-count crosscheck failed.")
}

msg("QC passed: metadata cell counts match pseudobulk cell selection.")

# ------------------------------------------------------------------------------
# 10. Descriptive disease and treatment effect sizes
# ------------------------------------------------------------------------------

mean_cpm_for <- function(subtype, condition_name) {

  cols <- pb_meta$pb_id[
    pb_meta$macrophage_class == subtype &
      as.character(pb_meta$condition) == condition_name
  ]

  cols <- intersect(cols, colnames(pb_cpm))

  if (length(cols) == 0L) {
    return(
      rep(
        NA_real_,
        nrow(pb_cpm)
      )
    )
  }

  if (length(cols) == 1L) {
    return(
      as.numeric(pb_cpm[, cols])
    )
  }

  rowMeans(
    pb_cpm[, cols, drop = FALSE]
  )
}

disease_results <- list()
treatment_results <- list()

for (subtype in MAJOR_CLASS_ORDER) {

  std_cpm <- mean_cpm_for(subtype, "STD")
  cdahfd_cpm <- mean_cpm_for(subtype, "CDAHFD")
  sham_cpm <- mean_cpm_for(subtype, "Sham")
  tx_cpm <- mean_cpm_for(subtype, "Tx")

  disease <- tibble(
    macrophage_class = subtype,
    gene = rownames(pb_cpm),
    mean_CPM_STD = std_cpm,
    mean_CPM_CDAHFD = cdahfd_cpm,
    log2FC_CDAHFD_vs_STD = log2(
      (cdahfd_cpm + CPM_PSEUDOCOUNT) /
        (std_cpm + CPM_PSEUDOCOUNT)
    )
  ) %>%
    arrange(
      desc(
        abs(log2FC_CDAHFD_vs_STD)
      )
    )

  treatment <- tibble(
    macrophage_class = subtype,
    gene = rownames(pb_cpm),
    mean_CPM_Sham = sham_cpm,
    mean_CPM_Tx = tx_cpm,
    log2FC_Tx_vs_Sham = log2(
      (tx_cpm + CPM_PSEUDOCOUNT) /
        (sham_cpm + CPM_PSEUDOCOUNT)
    )
  ) %>%
    arrange(
      desc(
        abs(log2FC_Tx_vs_Sham)
      )
    )

  disease_results[[subtype]] <- disease
  treatment_results[[subtype]] <- treatment

  write.csv(
    disease,
    file.path(
      PB_DIR,
      paste0(
        "06_Disease_",
        safe_filename(subtype),
        "_STD_vs_CDAHFD_v4.11.4.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    treatment,
    file.path(
      PB_DIR,
      paste0(
        "07_Treatment_",
        safe_filename(subtype),
        "_Sham_vs_Tx_v4.11.4.csv"
      )
    ),
    row.names = FALSE
  )
}

disease_combined <- bind_rows(disease_results)
treatment_combined <- bind_rows(treatment_results)

write.csv(
  disease_combined,
  file.path(
    PB_DIR,
    "08_Disease_all_subtypes_STD_vs_CDAHFD_v4.11.4.csv"
  ),
  row.names = FALSE
)

write.csv(
  treatment_combined,
  file.path(
    PB_DIR,
    "09_Treatment_all_subtypes_Sham_vs_Tx_v4.11.4.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. Anti-inflammatory MΦ manual-QC replication of known genes
# ------------------------------------------------------------------------------

anti_qc <- bind_rows(
  disease_results[["Anti-inflammatory-Mphi"]] %>%
    filter(gene %in% qc_genes_use)
)

write.csv(
  anti_qc,
  file.path(
    PB_DIR,
    "10_Anti_inflammatory_Mphi_manual_QC_known_genes_v4.11.4.csv"
  ),
  row.names = FALSE
)

# Expected diagnostic values should be non-zero and biologically variable.
if (
  all(
    abs(anti_qc$log2FC_CDAHFD_vs_STD) < 1e-8
  )
) {
  stop(
    "Anti-inflammatory-MΦ QC genes all have ~0 disease log2FC. ",
    "Pseudobulk rebuild failed."
  )
}

msg("QC passed: Anti-inflammatory-MΦ known genes show non-zero disease effects.")

# ------------------------------------------------------------------------------
# 12. Optional exploratory edgeR for Sham vs Tx
# ------------------------------------------------------------------------------

edger_all <- list()

if (HAS_EDGER) {

  for (subtype in MAJOR_CLASS_ORDER) {

    one_meta <- pb_meta %>%
      filter(
        macrophage_class == subtype,
        condition %in% c("Sham", "Tx")
      )

    if (
      sum(one_meta$condition == "Sham") >= 2L &&
      sum(one_meta$condition == "Tx") >= 2L
    ) {

      one_counts <- pb_counts[
        ,
        one_meta$pb_id,
        drop = FALSE
      ]

      group <- factor(
        as.character(one_meta$condition),
        levels = c("Sham", "Tx")
      )

      y <- edgeR::DGEList(
        counts = round(as.matrix(one_counts)),
        group = group
      )

      keep <- rowSums(edgeR::cpm(y) >= MIN_CPM) >= MIN_SAMPLES_EXPRESSED

      y <- y[
        keep,
        ,
        keep.lib.sizes = FALSE
      ]

      y <- edgeR::calcNormFactors(y)

      design <- model.matrix(~ group)

      y <- edgeR::estimateDisp(
        y,
        design,
        robust = TRUE
      )

      fit <- edgeR::glmQLFit(
        y,
        design,
        robust = TRUE
      )

      qlf <- edgeR::glmQLFTest(
        fit,
        coef = 2
      )

      tt <- edgeR::topTags(
        qlf,
        n = Inf,
        sort.by = "PValue"
      )$table %>%
        rownames_to_column("gene") %>%
        mutate(
          macrophage_class = subtype,
          .before = gene
        )

      edger_all[[subtype]] <- tt

      write.csv(
        tt,
        file.path(
          PB_DIR,
          paste0(
            "11_edgeR_",
            safe_filename(subtype),
            "_Sham_vs_Tx_exploratory_v4.11.4.csv"
          )
        ),
        row.names = FALSE
      )
    }
  }

  if (length(edger_all) > 0L) {
    write.csv(
      bind_rows(edger_all),
      file.path(
        PB_DIR,
        "12_edgeR_all_subtypes_Sham_vs_Tx_exploratory_v4.11.4.csv"
      ),
      row.names = FALSE
    )
  }
}

# ------------------------------------------------------------------------------
# 13. Effect-size figures
# ------------------------------------------------------------------------------

disease_top <- disease_combined %>%
  group_by(macrophage_class) %>%
  slice_max(
    order_by = abs(log2FC_CDAHFD_vs_STD),
    n = TOP_N_EFFECT,
    with_ties = FALSE
  ) %>%
  ungroup()

p_disease <- ggplot(
  disease_top,
  aes(
    x = log2FC_CDAHFD_vs_STD,
    y = reorder(gene, log2FC_CDAHFD_vs_STD)
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = log2FC_CDAHFD_vs_STD,
      yend = reorder(gene, log2FC_CDAHFD_vs_STD)
    ),
    linewidth = 0.5
  ) +
  geom_point(size = 2.2) +
  facet_wrap(
    ~ macrophage_class,
    scales = "free_y",
    ncol = 3,
    labeller = as_labeller(CLASS_LABELS)
  ) +
  labs(
    title = "Subtype-specific pseudobulk: STD → CDAHFD",
    subtitle = "Safe explicit rowSums pseudobulk | descriptive effect size only",
    x = "log2FC (CDAHFD / STD)",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

save_pdf(
  "01_SAFE_pseudobulk_Disease_top_log2FC_v4.11.4.pdf",
  p_disease,
  14,
  11
)

treatment_top <- treatment_combined %>%
  group_by(macrophage_class) %>%
  slice_max(
    order_by = abs(log2FC_Tx_vs_Sham),
    n = TOP_N_EFFECT,
    with_ties = FALSE
  ) %>%
  ungroup()

p_treatment <- ggplot(
  treatment_top,
  aes(
    x = log2FC_Tx_vs_Sham,
    y = reorder(gene, log2FC_Tx_vs_Sham)
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = log2FC_Tx_vs_Sham,
      yend = reorder(gene, log2FC_Tx_vs_Sham)
    ),
    linewidth = 0.5
  ) +
  geom_point(size = 2.2) +
  facet_wrap(
    ~ macrophage_class,
    scales = "free_y",
    ncol = 3,
    labeller = as_labeller(CLASS_LABELS)
  ) +
  labs(
    title = "Subtype-specific pseudobulk: Sham → Tx",
    subtitle = "Safe explicit rowSums pseudobulk",
    x = "log2FC (Tx / Sham)",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

save_pdf(
  "02_SAFE_pseudobulk_Treatment_top_log2FC_v4.11.4.pdf",
  p_treatment,
  14,
  11
)

# ------------------------------------------------------------------------------
# 14. Gene-level custom DotPlot
# ------------------------------------------------------------------------------

VALIDATION_USE <- lapply(
  VALIDATION_GENE_SETS,
  intersect,
  y = rownames(mphi[[ASSAY_USE]])
)

validation_features <- unique(
  unlist(
    VALIDATION_USE,
    use.names = FALSE
  )
)

mphi$subtype_condition_v4114 <- factor(
  paste(
    as.character(mphi$class_v4114),
    as.character(mphi$condition_v4114),
    sep = " | "
  ),
  levels = as.vector(
    outer(
      MAJOR_CLASS_ORDER,
      CONDITION_ORDER,
      function(a, b) {
        paste(a, b, sep = " | ")
      }
    )
  )
)

subtype_condition_levels <- levels(
  mphi$subtype_condition_v4114
)

marker_dot_df <- make_custom_dotplot_data(
  object = mphi,
  assay = ASSAY_USE,
  group_vector = mphi$subtype_condition_v4114,
  group_levels = subtype_condition_levels,
  features = validation_features
)

write.csv(
  marker_dot_df,
  file.path(
    TAB_DIR,
    "13_marker_DotPlot_numeric_data_v4.11.4.csv"
  ),
  row.names = FALSE
)

p_dot <- ggplot(
  marker_dot_df,
  aes(
    x = gene,
    y = group
  )
) +
  geom_point(
    aes(
      size = pct_expr,
      color = avg_expr_scaled
    )
  ) +
  scale_size_continuous(
    range = c(0.25, 7),
    limits = c(0, 100),
    name = "Percent\nexpressed"
  ) +
  scale_color_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Average\nexpression\n(z-score)"
  ) +
  labs(
    title = "Gene-level validation of fixed MΦ subtype programs",
    subtitle = "Custom DotPlot from normalized RNA data",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 7)
  )

save_pdf(
  "03_marker_gene_custom_DotPlot_v4.11.4.pdf",
  p_dot,
  22,
  10
)

# ------------------------------------------------------------------------------
# 15. Sample-level pseudobulk heatmap
# ------------------------------------------------------------------------------

heatmap_genes <- intersect(
  validation_features,
  rownames(pb_log2cpm)
)

hm <- pb_log2cpm[
  heatmap_genes,
  pb_meta$pb_id,
  drop = FALSE
]

hm_z <- safe_row_z(hm)

annotation_col <- pb_meta %>%
  select(
    pb_id,
    macrophage_class,
    condition,
    sample
  ) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$pb_id
annotation_col$pb_id <- NULL

grDevices::cairo_pdf(
  file.path(
    FIG_DIR,
    "04_SAFE_marker_gene_pseudobulk_sample_heatmap_v4.11.4.pdf"
  ),
  width = 18,
  height = 13
)

pheatmap::pheatmap(
  hm_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  color = grDevices::colorRampPalette(
    c("#0033FF", "#FFFFFF", "#FF1A1A")
  )(101),
  border_color = NA,
  annotation_col = annotation_col,
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 45,
  main = paste0(
    "SAFE marker-gene pseudobulk heatmap\n",
    "Explicit rowSums pseudobulk | gene-wise z-score"
  )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 16. Candidate ligand output
# ------------------------------------------------------------------------------

LIGAND_USE <- lapply(
  LIGAND_GENE_SETS,
  intersect,
  y = rownames(mphi[[ASSAY_USE]])
)

ligand_features <- unique(
  unlist(
    LIGAND_USE,
    use.names = FALSE
  )
)

ligand_dot_df <- make_custom_dotplot_data(
  object = mphi,
  assay = ASSAY_USE,
  group_vector = mphi$subtype_condition_v4114,
  group_levels = subtype_condition_levels,
  features = ligand_features
)

write.csv(
  ligand_dot_df,
  file.path(
    TAB_DIR,
    "14_ligand_DotPlot_numeric_data_v4.11.4.csv"
  ),
  row.names = FALSE
)

p_ligand_dot <- ggplot(
  ligand_dot_df,
  aes(
    x = gene,
    y = group
  )
) +
  geom_point(
    aes(
      size = pct_expr,
      color = avg_expr_scaled
    )
  ) +
  scale_size_continuous(
    range = c(0.25, 8),
    limits = c(0, 100),
    name = "Percent\nexpressed"
  ) +
  scale_color_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Average\nexpression\n(z-score)"
  ) +
  labs(
    title = "Candidate macrophage ligand-output map",
    subtitle = "Expression screen only",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 7)
  )

save_pdf(
  "05_candidate_ligand_custom_DotPlot_v4.11.4.pdf",
  p_ligand_dot,
  19,
  9.5
)

ligand_effect <- bind_rows(
  lapply(
    MAJOR_CLASS_ORDER,
    function(subtype) {

      dis <- disease_results[[subtype]] %>%
        filter(gene %in% ligand_features) %>%
        transmute(
          macrophage_class = subtype,
          comparison = "STD_to_CDAHFD",
          gene = gene,
          log2FC = log2FC_CDAHFD_vs_STD
        )

      trt <- treatment_results[[subtype]] %>%
        filter(gene %in% ligand_features) %>%
        transmute(
          macrophage_class = subtype,
          comparison = "Sham_to_Tx",
          gene = gene,
          log2FC = log2FC_Tx_vs_Sham
        )

      bind_rows(dis, trt)
    }
  )
)

ligand_membership <- bind_rows(
  lapply(
    names(LIGAND_USE),
    function(set_name) {
      tibble(
        ligand_set = set_name,
        gene = LIGAND_USE[[set_name]]
      )
    }
  )
)

ligand_effect <- ligand_effect %>%
  left_join(
    ligand_membership,
    by = "gene"
  )

write.csv(
  ligand_effect,
  file.path(
    TAB_DIR,
    "15_candidate_ligand_effect_sizes_v4.11.4.csv"
  ),
  row.names = FALSE
)

p_ligand_effect <- ggplot(
  ligand_effect,
  aes(
    x = log2FC,
    y = gene
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  geom_point(size = 2.1) +
  facet_grid(
    ligand_set ~ comparison,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    title = "Candidate macrophage ligand changes",
    subtitle = "Safe pseudobulk effect-size screen",
    x = "log2 fold-change",
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 8)
  )

save_pdf(
  "06_candidate_ligand_SAFE_effect_size_v4.11.4.pdf",
  p_ligand_effect,
  11,
  11
)

# ------------------------------------------------------------------------------
# 17. Summary
# ------------------------------------------------------------------------------

p_summary <- (
  p_disease /
  p_treatment /
  p_dot /
  p_ligand_dot /
  p_ligand_effect
) +
  patchwork::plot_layout(
    heights = c(1, 1, 0.95, 0.85, 0.85)
  ) +
  patchwork::plot_annotation(
    title = "Mouse MASH MΦ SAFE pseudobulk / gene / ligand analysis v4.11.4",
    subtitle = "Explicit sample × subtype rowSums pseudobulk",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 19
      )
    )
  )

save_pdf(
  "07_Mphi_SAFE_pseudobulk_gene_ligand_summary_v4.11.4.pdf",
  p_summary,
  22,
  37
)

# ------------------------------------------------------------------------------
# 18. README / session info
# ------------------------------------------------------------------------------

readme <- c(
  "Mouse MASH MΦ SAFE pseudobulk / gene / ligand analysis v4.11.4",
  "",
  paste0("Input: ", INPUT_RDS),
  "",
  "Pseudobulk method:",
  "  Explicit Matrix::rowSums(counts[, cells]) for every sample x subtype.",
  "  This is the same method manually validated before v4.11.4.",
  "",
  "Safety checks:",
  "  counts cell order == Seurat cell order",
  "  metadata row order == Seurat cell order",
  "  no zero-count pseudobulk libraries",
  "  no identical pseudobulk columns",
  "  sample x subtype cell counts match metadata",
  "  Anti-inflammatory-MΦ known genes must show non-zero STD/CDAHFD effects",
  "",
  "Disease comparison:",
  "  STD vs CDAHFD = descriptive effect size only.",
  "",
  "Treatment comparison:",
  "  Sham vs Tx = descriptive effect size always.",
  paste0(
    "  edgeR exploratory = ",
    ifelse(HAS_EDGER, "ON", "OFF")
  ),
  "",
  "Primary QC outputs:",
  "  00_pseudobulk_QC_library_size",
  "  01_pseudobulk_QC_marker_counts_CPM",
  "  05_pseudobulk_metadata_cellcount_crosscheck",
  "  10_Anti_inflammatory_Mphi_manual_QC_known_genes",
  "",
  "Primary figures:",
  "  01 SAFE disease log2FC",
  "  02 SAFE treatment log2FC",
  "  03 marker custom DotPlot",
  "  04 SAFE pseudobulk heatmap",
  "  05 candidate ligand DotPlot",
  "  06 candidate ligand SAFE effect-size",
  "  07 summary"
)

writeLines(
  readme,
  file.path(
    OUTPUT_DIR,
    "README_v4.11.4.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    LOG_DIR,
    "sessionInfo_v4.11.4.txt"
  )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)
msg("edgeR exploratory treatment analysis: ", ifelse(HAS_EDGER, "ON", "OFF"))
