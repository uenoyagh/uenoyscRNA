#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6630)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(nichenetr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# NicheNet validation of macrophage -> Hepatocyte communication
#
# Version: v6.6.3
#
# INPUTS
# ------
# 1) v6.6.0 interaction-ready Seurat object
#    5 macrophage subtypes + 5 primary Hepatocyte states
#
# 2) v6.6.2.1 refined CellChat high-confidence LR shortlist
#
# BIOLOGICAL REPLICATES
# ---------------------
#   Sham1
#   Sham20
#   Tx17
#   Tx5
#
# SENDERS
# -------
#   Mphi_Anti-inflammatory
#   Mphi_Inflammatory
#   Mphi_ECM-associated-inflammatory
#   Mphi_Repair-Resolution
#   Mphi_Lipid-associated-TREM2
#
# RECEIVERS
# ---------
#   Hep_Periportal
#   Hep_Pericentral
#   Hep_Injury-inflammatory
#   Hep_Intermediate
#   Hep_Cycling
#
# PRIMARY QUESTION
# ----------------
# Which macrophage-derived ligands can explain Tx-associated transcriptional
# changes in Hepatocytes, especially the Injury/inflammatory Hepatocyte state?
#
# IMPORTANT DESIGN FEATURES
# -------------------------
# - Receiver DE is sample-level pseudobulk, NOT cell-level pseudoreplication.
# - Sender DE is sample-level pseudobulk.
# - edgeR quasi-likelihood is used with biological n=2 Sham vs n=2 Tx.
# - Tx-up and Tx-down receiver programs are analyzed separately.
# - Mouse NicheNet v2 prior model is used directly.
# - NicheNet ligand activity is integrated with:
#     a) sender ligand expression fraction
#     b) sender pseudobulk ligand Tx/Sham logFC
#     c) v6.6.2.1 replicate-aware supported CellChat evidence
# - CellChat is NOT rerun.
# - NicheNet ligand activity does NOT prove ligand direction.
#
# PRIMARY MECHANISTIC FAMILIES
# ----------------------------
#   PDGF
#   FN1
#   SPP1
#   SEMA4
#   TNF
#   MIF
#   COLLAGEN
#   LAMININ
#
# INTERPRETATION
# --------------
# Tx_up_program:
#   NicheNet asks which ligands are capable of explaining genes higher in Tx.
#
# Tx_down_program:
#   NicheNet asks which ligands are capable of explaining genes lower in Tx.
#
# Directional concordance additionally requires:
#   - sender ligand pseudobulk direction matching the receiver program
#   - supported CellChat Tx-Sham direction matching the receiver program
#
# Because n=2/condition, all formal DE statistics remain exploratory.
# ==============================================================================


# ==============================================================================
# 1. Helpers
# ==============================================================================

msg <- function(...) {
  message(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste0(...)
  )
}


save_pdf <- function(
  p,
  file,
  width,
  height
) {
  grDevices::pdf(
    file,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  print(p)
  grDevices::dev.off()
}


safe_join_rna <- function(object) {

  if (
    "RNA" %in%
      Assays(object) &&
    length(
      Layers(
        object[["RNA"]]
      )
    ) > 1
  ) {

    object[["RNA"]] <- JoinLayers(
      object[["RNA"]]
    )
  }

  object
}


download_if_missing <- function(
  url,
  dest
) {

  if (
    file.exists(dest) &&
    file.info(dest)$size > 1000
  ) {

    msg(
      "Using cached model: ",
      dest
    )

    return(
      invisible(dest)
    )
  }

  msg(
    "Downloading NicheNet model: ",
    url
  )

  tmp <- paste0(
    dest,
    ".tmp"
  )

  if (
    file.exists(tmp)
  ) {
    unlink(tmp)
  }

  utils::download.file(
    url = url,
    destfile = tmp,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )

  if (
    !file.exists(tmp) ||
    file.info(tmp)$size <= 1000
  ) {
    stop(
      "Model download failed or file too small: ",
      url
    )
  }

  ok <- file.rename(
    tmp,
    dest
  )

  if (
    !ok
  ) {
    stop(
      "Could not move downloaded model to: ",
      dest
    )
  }

  invisible(dest)
}


aggregate_counts_by_sample <- function(
  counts,
  meta,
  celltype_col,
  sample_col,
  celltype_value,
  samples
) {

  out <- lapply(
    samples,
    function(sample_name) {

      cells <- rownames(meta)[
        as.character(
          meta[[celltype_col]]
        ) ==
          celltype_value &
        as.character(
          meta[[sample_col]]
        ) ==
          sample_name
      ]

      if (
        length(cells) == 0
      ) {
        stop(
          "No cells for ",
          celltype_value,
          " / ",
          sample_name
        )
      }

      Matrix::rowSums(
        counts[
          ,
          cells,
          drop = FALSE
        ]
      )
    }
  )

  mat <- do.call(
    cbind,
    out
  )

  rownames(mat) <- rownames(counts)
  colnames(mat) <- samples

  mat
}


run_edger_ql <- function(
  pb_counts,
  samples
) {

  condition <- factor(
    ifelse(
      grepl(
        "^Tx",
        samples
      ),
      "Tx",
      "Sham"
    ),
    levels = c(
      "Sham",
      "Tx"
    )
  )

  y <- edgeR::DGEList(
    counts = pb_counts,
    group = condition
  )

  keep <- edgeR::filterByExpr(
    y,
    group = condition
  )

  if (
    sum(keep) < 100
  ) {
    stop(
      "Too few genes passed edgeR filterByExpr: ",
      sum(keep)
    )
  }

  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  y <- edgeR::calcNormFactors(y)

  design <- model.matrix(
    ~ condition
  )

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

  edgeR::topTags(
    qlf,
    n = Inf,
    sort.by = "PValue"
  )$table %>%
    as.data.frame() %>%
    rownames_to_column(
      "gene"
    ) %>%
    as_tibble()
}


pct_expressed <- function(
  counts,
  cells
) {

  if (
    !length(cells)
  ) {
    return(
      numeric()
    )
  }

  x <- counts[
    ,
    cells,
    drop = FALSE
  ]

  Matrix::rowMeans(
    x > 0
  )
}


select_receiver_geneset <- function(
  de_table,
  direction,
  ligand_target_genes,
  min_genes = 20,
  max_genes = 200
) {

  if (
    direction == "Tx_up_program"
  ) {

    direction_primary <- de_table %>%
      filter(
        logFC >= 0.5
      )

    direction_fallback <- de_table %>%
      filter(
        logFC >= 0.25
      )

  } else if (
    direction == "Tx_down_program"
  ) {

    direction_primary <- de_table %>%
      filter(
        logFC <= -0.5
      )

    direction_fallback <- de_table %>%
      filter(
        logFC <= -0.25
      )

  } else {

    stop(
      "Unknown receiver program: ",
      direction
    )
  }

  primary <- direction_primary %>%
    filter(
      FDR <= 0.10,
      gene %in%
        ligand_target_genes
    ) %>%
    arrange(
      PValue
    )

  if (
    nrow(primary) >= min_genes
  ) {

    selected <- primary %>%
      slice_head(
        n = max_genes
      )

    method <-
      "FDR<=0.10_and_absLogFC>=0.5"

  } else {

    relaxed <- direction_primary %>%
      filter(
        PValue <= 0.05,
        gene %in%
          ligand_target_genes
      ) %>%
      arrange(
        PValue
      )

    if (
      nrow(relaxed) >= min_genes
    ) {

      selected <- relaxed %>%
        slice_head(
          n = max_genes
        )

      method <-
        "PValue<=0.05_and_absLogFC>=0.5_exploratory"

    } else {

      ranked <- direction_fallback %>%
        filter(
          gene %in%
            ligand_target_genes
        ) %>%
        arrange(
          PValue
        )

      n_take <- min(
        max_genes,
        max(
          min_genes,
          min(
            100,
            nrow(ranked)
          )
        )
      )

      selected <- ranked %>%
        slice_head(
          n = n_take
        )

      method <-
        "Ranked_fallback_absLogFC>=0.25_exploratory"
    }
  }

  list(
    genes =
      unique(
        selected$gene
      ),
    table =
      selected,
    method =
      method
  )
}


safe_predict_ligand_activities <- function(
  geneset,
  background,
  ligand_target_matrix,
  potential_ligands
) {

  geneset <- intersect(
    geneset,
    rownames(
      ligand_target_matrix
    )
  )

  background <- intersect(
    background,
    rownames(
      ligand_target_matrix
    )
  )

  geneset <- intersect(
    geneset,
    background
  )

  potential_ligands <- intersect(
    potential_ligands,
    colnames(
      ligand_target_matrix
    )
  )

  if (
    length(geneset) < 10
  ) {

    warning(
      "NicheNet skipped: fewer than 10 receiver genes in geneset."
    )

    return(
      tibble()
    )
  }

  if (
    length(potential_ligands) < 2
  ) {

    warning(
      "NicheNet skipped: fewer than 2 potential ligands."
    )

    return(
      tibble()
    )
  }

  nichenetr::predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  ) %>%
    as_tibble() %>%
    arrange(
      desc(
        aupr_corrected
      )
    ) %>%
    mutate(
      nichenet_rank =
        row_number(),

      nichenet_percentile =
        ifelse(
          n() > 1,
          1 -
            (
              nichenet_rank -
                1
            ) /
              (
                n() -
                  1
              ),
          1
        )
    )
}


extract_top_target_links <- function(
  ligand_activities,
  geneset_table,
  ligand_target_matrix,
  receiver,
  receiver_program,
  top_ligands = 20,
  top_targets_per_ligand = 20
) {

  if (
    !nrow(ligand_activities)
  ) {
    return(
      tibble()
    )
  }

  ligands_use <- ligand_activities %>%
    slice_head(
      n = top_ligands
    ) %>%
    pull(
      test_ligand
    )

  lapply(
    ligands_use,
    function(ligand_name) {

      if (
        !ligand_name %in%
          colnames(
            ligand_target_matrix
          )
      ) {
        return(
          tibble()
        )
      }

      targets <- intersect(
        geneset_table$gene,
        rownames(
          ligand_target_matrix
        )
      )

      if (
        !length(targets)
      ) {
        return(
          tibble()
        )
      }

      rp <- ligand_target_matrix[
        targets,
        ligand_name
      ]

      tibble(
        target =
          names(rp),
        regulatory_potential =
          as.numeric(rp)
      ) %>%
        filter(
          is.finite(
            regulatory_potential
          ),
          regulatory_potential > 0
        ) %>%
        arrange(
          desc(
            regulatory_potential
          )
        ) %>%
        slice_head(
          n =
            top_targets_per_ligand
        ) %>%
        left_join(
          geneset_table %>%
            select(
              gene,
              receiver_logFC =
                logFC,
              receiver_logCPM =
                logCPM,
              receiver_PValue =
                PValue,
              receiver_FDR =
                FDR
            ),
          by = c(
            "target" =
              "gene"
          )
        ) %>%
        mutate(
          ligand =
            ligand_name,
          receiver =
            receiver,
          receiver_program =
            receiver_program
        )
    }
  ) %>%
    bind_rows()
}


mechanism_family <- function(gene) {

  gene_chr <- as.character(gene)

  case_when(
    gene_chr == "Fn1" ~
      "FN1",

    gene_chr %in%
      c(
        "Pdgfa",
        "Pdgfb",
        "Pdgfc",
        "Pdgfd"
      ) ~
      "PDGF",

    gene_chr == "Spp1" ~
      "SPP1",

    grepl(
      "^Sema4",
      gene_chr
    ) ~
      "SEMA4",

    gene_chr == "Tnf" ~
      "TNF",

    gene_chr == "Mif" ~
      "MIF",

    grepl(
      "^Col[0-9]",
      gene_chr
    ) ~
      "COLLAGEN",

    grepl(
      "^Lam[abc]",
      gene_chr
    ) ~
      "LAMININ",

    gene_chr %in%
      c(
        "Tgfb1",
        "Tgfb2",
        "Tgfb3"
      ) ~
      "TGFB",

    gene_chr == "Gas6" ~
      "GAS6",

    gene_chr == "Hgf" ~
      "HGF",

    gene_chr == "Igf1" ~
      "IGF1",

    gene_chr == "Areg" ~
      "AREG",

    gene_chr == "Hbegf" ~
      "HBEGF",

    gene_chr == "Plau" ~
      "PLAU",

    gene_chr == "App" ~
      "APP",

    TRUE ~
      "Other"
  )
}


direction_matches_program <- function(
  receiver_program,
  value
) {

  case_when(
    is.na(value) ~
      FALSE,

    receiver_program ==
      "Tx_up_program" &
      value > 0 ~
      TRUE,

    receiver_program ==
      "Tx_down_program" &
      value < 0 ~
      TRUE,

    TRUE ~
      FALSE
  )
}


# ==============================================================================
# 2. Paths
# ==============================================================================

ROOT <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"


INPUT_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_interaction_ready_v6.6.0",
  "RDS",
  "Mouse_MASH_Mphi5_Hep5_interaction_ready_v6.6.0.rds"
)


CELLCHAT_HIGHCONF <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_CellChat_refine_v6.6.2.1",
  "Tables",
  "14_high_confidence_LR_shortlist_v6.6.2.1.csv"
)


CELLCHAT_REFINED_ALL <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_CellChat_refine_v6.6.2.1",
  "Tables",
  "03_LR_refined_replicate_comparison_v6.6.2.1.csv"
)


OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_Hep5_NicheNet_v6.6.3"
)


RDS_OUT <- file.path(
  OUT,
  "RDS"
)


TAB_OUT <- file.path(
  OUT,
  "Tables"
)


FIG_OUT <- file.path(
  OUT,
  "Figures"
)


LOG_OUT <- file.path(
  OUT,
  "Logs"
)


MODEL_DIR <- file.path(
  ROOT,
  "Reference_Models",
  "NicheNet_mouse_v2"
)


for (
  d in c(
    OUT,
    RDS_OUT,
    TAB_OUT,
    FIG_OUT,
    LOG_OUT,
    MODEL_DIR
  )
) {

  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 3. Settings
# ==============================================================================

GROUP_COL <-
  "interaction_celltype_v660"


SAMPLE_COL <-
  "sample_interaction_v660"


SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)


SENDERS <- c(
  "Mphi_Anti-inflammatory",
  "Mphi_Inflammatory",
  "Mphi_ECM-associated-inflammatory",
  "Mphi_Repair-Resolution",
  "Mphi_Lipid-associated-TREM2"
)


RECEIVERS <- c(
  "Hep_Periportal",
  "Hep_Pericentral",
  "Hep_Injury-inflammatory",
  "Hep_Intermediate",
  "Hep_Cycling"
)


PRIMARY_RECEIVER <-
  "Hep_Injury-inflammatory"


RECEIVER_PROGRAMS <- c(
  "Tx_up_program",
  "Tx_down_program"
)


EXPRESSION_PCT_THRESHOLD <-
  0.10


MIN_GENESET <-
  20


MAX_GENESET <-
  200


TOP_LIGANDS_FOR_TARGETS <-
  20


TOP_TARGETS_PER_LIGAND <-
  20


PRIMARY_MECHANISM_FAMILIES <- c(
  "PDGF",
  "FN1",
  "SPP1",
  "SEMA4",
  "TNF",
  "MIF",
  "COLLAGEN",
  "LAMININ"
)


# ==============================================================================
# 4. NicheNet mouse v2 model
# ==============================================================================

LR_URL <-
  "https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"


LTM_URL <-
  "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"


LR_RDS <- file.path(
  MODEL_DIR,
  "lr_network_mouse_21122021.rds"
)


LTM_RDS <- file.path(
  MODEL_DIR,
  "ligand_target_matrix_nsga2r_final_mouse.rds"
)


download_if_missing(
  LR_URL,
  LR_RDS
)


download_if_missing(
  LTM_URL,
  LTM_RDS
)


msg(
  "Loading NicheNet mouse v2 prior model..."
)


lr_network <- readRDS(
  LR_RDS
) %>%
  as_tibble() %>%
  distinct(
    from,
    to
  )


ligand_target_matrix <- readRDS(
  LTM_RDS
)


model_audit <- tibble(
  item = c(
    "lr_network_rows",
    "lr_network_ligands",
    "lr_network_receptors",
    "ligand_target_rows_targets",
    "ligand_target_cols_ligands"
  ),
  value = c(
    nrow(
      lr_network
    ),
    n_distinct(
      lr_network$from
    ),
    n_distinct(
      lr_network$to
    ),
    nrow(
      ligand_target_matrix
    ),
    ncol(
      ligand_target_matrix
    )
  )
)


write.csv(
  model_audit,
  file.path(
    TAB_OUT,
    "01_NicheNet_mouse_model_audit_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 5. Load v6.6.0 interaction-ready object
# ==============================================================================

if (
  !file.exists(
    INPUT_RDS
  )
) {

  stop(
    "Input v6.6.0 RDS not found: ",
    INPUT_RDS
  )
}


msg(
  "Loading v6.6.0 interaction-ready object..."
)


obj <- readRDS(
  INPUT_RDS
)


if (
  !"RNA" %in%
    Assays(
      obj
    )
) {

  stop(
    "RNA assay missing."
  )
}


DefaultAssay(
  obj
) <- "RNA"


obj <- safe_join_rna(
  obj
)


required_meta <- c(
  GROUP_COL,
  SAMPLE_COL
)


missing_meta <- setdiff(
  required_meta,
  colnames(
    obj@meta.data
  )
)


if (
  length(
    missing_meta
  )
) {

  stop(
    "Missing metadata: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}


present_groups <- unique(
  as.character(
    obj@meta.data[[
      GROUP_COL
    ]]
  )
)


missing_groups <- setdiff(
  c(
    SENDERS,
    RECEIVERS
  ),
  present_groups
)


if (
  length(
    missing_groups
  )
) {

  stop(
    "Missing required interaction states: ",
    paste(
      missing_groups,
      collapse = ", "
    )
  )
}


counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)


meta <- obj@meta.data


msg(
  "Cells: ",
  ncol(obj),
  " | Genes: ",
  nrow(obj)
)


# ==============================================================================
# 6. Load v6.6.2.1 CellChat evidence
# ==============================================================================

cellchat_highconf <- tibble()


if (
  file.exists(
    CELLCHAT_HIGHCONF
  )
) {

  cellchat_highconf <- read.csv(
    CELLCHAT_HIGHCONF,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()

  msg(
    "Loaded CellChat high-confidence LR rows: ",
    nrow(
      cellchat_highconf
    )
  )

} else {

  warning(
    "v6.6.2.1 high-confidence CellChat LR table not found. ",
    "NicheNet will run without high-confidence CellChat integration."
  )
}


cellchat_all <- tibble()


if (
  file.exists(
    CELLCHAT_REFINED_ALL
  )
) {

  cellchat_all <- read.csv(
    CELLCHAT_REFINED_ALL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()

  msg(
    "Loaded all refined CellChat LR rows: ",
    nrow(
      cellchat_all
    )
  )

} else {

  warning(
    "All refined CellChat LR table not found."
  )
}


# ==============================================================================
# 7. Expression fractions
# ==============================================================================

msg(
  "Computing expression fractions for 5 Mphi + 5 Hep states..."
)


expr_fraction_list <- list()


for (
  celltype_name in c(
    SENDERS,
    RECEIVERS
  )
) {

  cells <- rownames(meta)[
    as.character(
      meta[[
        GROUP_COL
      ]]
    ) ==
      celltype_name
  ]

  pct <- pct_expressed(
    counts,
    cells
  )

  expr_fraction_list[[
    celltype_name
  ]] <- tibble(
    gene =
      names(pct),
    celltype =
      celltype_name,
    pct_expressed =
      as.numeric(pct),
    expressed_ge_10pct =
      pct >=
        EXPRESSION_PCT_THRESHOLD
  )
}


expr_fraction <- bind_rows(
  expr_fraction_list
)


write.csv(
  expr_fraction,
  file.path(
    TAB_OUT,
    "02_expression_fraction_by_celltype_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Receiver pseudobulk DE
# ==============================================================================

msg(
  "Running Hepatocyte receiver pseudobulk DE..."
)


receiver_de_list <- list()


for (
  receiver_name in RECEIVERS
) {

  msg(
    "Receiver DE: ",
    receiver_name
  )

  pb <- aggregate_counts_by_sample(
    counts = counts,
    meta = meta,
    celltype_col = GROUP_COL,
    sample_col = SAMPLE_COL,
    celltype_value = receiver_name,
    samples = SAMPLES
  )

  de <- run_edger_ql(
    pb_counts = pb,
    samples = SAMPLES
  ) %>%
    mutate(
      receiver =
        receiver_name,
      comparison =
        "Tx_vs_Sham",
      interpretation =
        "positive_logFC_means_higher_in_Tx"
    )

  receiver_de_list[[
    receiver_name
  ]] <- de

  write.csv(
    de,
    file.path(
      TAB_OUT,
      paste0(
        "03_receiver_DE_",
        gsub(
          "[^A-Za-z0-9]+",
          "_",
          receiver_name
        ),
        "_Tx_vs_Sham_v6.6.3.csv"
      )
    ),
    row.names = FALSE
  )
}


receiver_de_all <- bind_rows(
  receiver_de_list
)


write.csv(
  receiver_de_all,
  file.path(
    TAB_OUT,
    "04_receiver_DE_all_Hep_states_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. Sender pseudobulk DE
# ==============================================================================

msg(
  "Running macrophage sender pseudobulk DE..."
)


sender_de_list <- list()


for (
  sender_name in SENDERS
) {

  msg(
    "Sender DE: ",
    sender_name
  )

  pb <- aggregate_counts_by_sample(
    counts = counts,
    meta = meta,
    celltype_col = GROUP_COL,
    sample_col = SAMPLE_COL,
    celltype_value = sender_name,
    samples = SAMPLES
  )

  de <- run_edger_ql(
    pb_counts = pb,
    samples = SAMPLES
  ) %>%
    mutate(
      sender =
        sender_name,
      comparison =
        "Tx_vs_Sham",
      interpretation =
        "positive_logFC_means_higher_in_Tx"
    )

  sender_de_list[[
    sender_name
  ]] <- de

  write.csv(
    de,
    file.path(
      TAB_OUT,
      paste0(
        "05_sender_DE_",
        gsub(
          "[^A-Za-z0-9]+",
          "_",
          sender_name
        ),
        "_Tx_vs_Sham_v6.6.3.csv"
      )
    ),
    row.names = FALSE
  )
}


sender_de_all <- bind_rows(
  sender_de_list
)


write.csv(
  sender_de_all,
  file.path(
    TAB_OUT,
    "06_sender_DE_all_Mphi_states_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. NicheNet per receiver x receiver program
# ==============================================================================

msg(
  "Running NicheNet ligand activity analysis..."
)


nichenet_activity_all <- list()

geneset_audit <- list()

target_link_all <- list()

sender_receiver_ligand_map_all <- list()


for (
  receiver_name in RECEIVERS
) {

  msg(
    "============================================================"
  )

  msg(
    "NicheNet receiver: ",
    receiver_name
  )

  receiver_expr <- expr_fraction %>%
    filter(
      celltype ==
        receiver_name,
      expressed_ge_10pct
    ) %>%
    pull(
      gene
    )

  expressed_receptors <- intersect(
    receiver_expr,
    unique(
      lr_network$to
    )
  )

  background <- intersect(
    receiver_expr,
    rownames(
      ligand_target_matrix
    )
  )

  sender_ligand_map <- lapply(
    SENDERS,
    function(sender_name) {

      sender_expr <- expr_fraction %>%
        filter(
          celltype ==
            sender_name,
          expressed_ge_10pct
        ) %>%
        pull(
          gene
        )

      lr_use <- lr_network %>%
        filter(
          from %in%
            sender_expr,
          to %in%
            expressed_receptors,
          from %in%
            colnames(
              ligand_target_matrix
            )
        )

      if (
        !nrow(lr_use)
      ) {

        return(
          tibble()
        )
      }

      lr_use %>%
        distinct(
          ligand =
            from,
          receptor =
            to
        ) %>%
        mutate(
          sender =
            sender_name,
          receiver =
            receiver_name
        )
    }
  ) %>%
    bind_rows()

  sender_receiver_ligand_map_all[[
    receiver_name
  ]] <-
    sender_ligand_map

  potential_ligands <- unique(
    sender_ligand_map$ligand
  )

  msg(
    "Potential ligands: ",
    length(
      potential_ligands
    )
  )

  for (
    receiver_program in RECEIVER_PROGRAMS
  ) {

    msg(
      "  Program: ",
      receiver_program
    )

    selected <- select_receiver_geneset(
      de_table =
        receiver_de_list[[
          receiver_name
        ]],
      direction =
        receiver_program,
      ligand_target_genes =
        rownames(
          ligand_target_matrix
        ),
      min_genes =
        MIN_GENESET,
      max_genes =
        MAX_GENESET
    )

    geneset <- selected$genes

    geneset_key <- paste(
      receiver_name,
      receiver_program,
      sep = "__"
    )

    geneset_audit[[
      geneset_key
    ]] <- tibble(
      receiver =
        receiver_name,
      receiver_program =
        receiver_program,
      geneset_method =
        selected$method,
      n_geneset =
        length(geneset),
      n_background =
        length(background),
      n_potential_ligands =
        length(
          potential_ligands
        )
    )

    write.csv(
      selected$table,
      file.path(
        TAB_OUT,
        paste0(
          "07_geneset_",
          gsub(
            "[^A-Za-z0-9]+",
            "_",
            receiver_name
          ),
          "_",
          receiver_program,
          "_v6.6.3.csv"
        )
      ),
      row.names = FALSE
    )

    ligand_activity <- safe_predict_ligand_activities(
      geneset =
        geneset,
      background =
        background,
      ligand_target_matrix =
        ligand_target_matrix,
      potential_ligands =
        potential_ligands
    )

    if (
      !nrow(ligand_activity)
    ) {
      next
    }

    ligand_activity <- ligand_activity %>%
      mutate(
        receiver =
          receiver_name,
        receiver_program =
          receiver_program,
        geneset_method =
          selected$method,
        n_geneset =
          length(geneset),
        inferred_activity_direction =
          ifelse(
            receiver_program ==
              "Tx_up_program",
            "Candidate_Tx_active_ligand",
            "Candidate_Sham_active_or_Tx_lost_ligand"
          )
      )

    nichenet_activity_all[[
      geneset_key
    ]] <-
      ligand_activity

    target_links <- extract_top_target_links(
      ligand_activities =
        ligand_activity,
      geneset_table =
        selected$table,
      ligand_target_matrix =
        ligand_target_matrix,
      receiver =
        receiver_name,
      receiver_program =
        receiver_program,
      top_ligands =
        TOP_LIGANDS_FOR_TARGETS,
      top_targets_per_ligand =
        TOP_TARGETS_PER_LIGAND
    )

    if (
      nrow(target_links)
    ) {

      target_link_all[[
        geneset_key
      ]] <-
        target_links
    }
  }
}


nichenet_activity <- bind_rows(
  nichenet_activity_all
)


geneset_audit_df <- bind_rows(
  geneset_audit
)


target_links <- bind_rows(
  target_link_all
)


sender_receiver_ligand_map <- bind_rows(
  sender_receiver_ligand_map_all
)


write.csv(
  geneset_audit_df,
  file.path(
    TAB_OUT,
    "08_NicheNet_geneset_audit_v6.6.3.csv"
  ),
  row.names = FALSE
)


write.csv(
  sender_receiver_ligand_map,
  file.path(
    TAB_OUT,
    "09_sender_receiver_potential_ligand_receptor_map_v6.6.3.csv"
  ),
  row.names = FALSE
)


write.csv(
  nichenet_activity,
  file.path(
    TAB_OUT,
    "10_NicheNet_ligand_activity_all_Hep_receivers_v6.6.3.csv"
  ),
  row.names = FALSE
)


write.csv(
  target_links,
  file.path(
    TAB_OUT,
    "11_NicheNet_top_ligand_target_links_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Expand NicheNet activities to sender-specific candidate rows
# ==============================================================================

if (
  !nrow(
    nichenet_activity
  )
) {

  stop(
    "No NicheNet ligand activity results were generated."
  )
}


integrated <- nichenet_activity %>%
  rename(
    ligand =
      test_ligand
  ) %>%
  inner_join(
    sender_receiver_ligand_map %>%
      distinct(
        sender,
        receiver,
        ligand
      ),
    by = c(
      "receiver",
      "ligand"
    )
  )


sender_ligand_de <- sender_de_all %>%
  select(
    sender,
    ligand =
      gene,
    sender_ligand_logFC =
      logFC,
    sender_ligand_logCPM =
      logCPM,
    sender_ligand_PValue =
      PValue,
    sender_ligand_FDR =
      FDR
  )


integrated <- integrated %>%
  left_join(
    sender_ligand_de,
    by = c(
      "sender",
      "ligand"
    )
  )


sender_expr_for_join <- expr_fraction %>%
  filter(
    celltype %in%
      SENDERS
  ) %>%
  select(
    sender =
      celltype,
    ligand =
      gene,
    sender_ligand_pct_expressed =
      pct_expressed
  )


integrated <- integrated %>%
  left_join(
    sender_expr_for_join,
    by = c(
      "sender",
      "ligand"
    )
  )


# ==============================================================================
# 12. Integrate v6.6.2.1 supported CellChat evidence
# ==============================================================================

if (
  nrow(
    cellchat_highconf
  )
) {

  required_cc <- c(
    "direction",
    "source",
    "target",
    "ligand",
    "receptor",
    "pathway_name",
    "supported_Tx_minus_Sham",
    "supported_pairwise_direction_consistency",
    "total_support_samples",
    "confidence_class"
  )

  missing_cc <- setdiff(
    required_cc,
    colnames(
      cellchat_highconf
    )
  )

  if (
    length(
      missing_cc
    )
  ) {

    stop(
      "v6.6.2.1 high-confidence CellChat table missing columns: ",
      paste(
        missing_cc,
        collapse = ", "
      )
    )
  }

  cellchat_summary <- cellchat_highconf %>%
    filter(
      direction ==
        "Mphi_to_Hep"
    ) %>%
    group_by(
      source,
      target,
      ligand
    ) %>%
    arrange(
      desc(
        abs(
          supported_Tx_minus_Sham
        )
      ),
      .by_group = TRUE
    ) %>%
    summarise(
      CellChat_highconf =
        TRUE,

      CellChat_best_supported_delta =
        first(
          supported_Tx_minus_Sham
        ),

      CellChat_direction_consistency =
        max(
          supported_pairwise_direction_consistency,
          na.rm = TRUE
        ),

      CellChat_support_samples =
        max(
          total_support_samples,
          na.rm = TRUE
        ),

      CellChat_confidence_class =
        first(
          confidence_class
        ),

      CellChat_pathways =
        paste(
          unique(
            pathway_name
          ),
          collapse = "|"
        ),

      CellChat_receptors =
        paste(
          unique(
            receptor
          ),
          collapse = "|"
        ),

      .groups = "drop"
    ) %>%
    rename(
      sender =
        source,
      receiver =
        target
    )

  integrated <- integrated %>%
    left_join(
      cellchat_summary,
      by = c(
        "sender",
        "receiver",
        "ligand"
      )
    )

} else {

  integrated$CellChat_highconf <-
    FALSE

  integrated$CellChat_best_supported_delta <-
    NA_real_

  integrated$CellChat_direction_consistency <-
    NA_real_

  integrated$CellChat_support_samples <-
    NA_real_

  integrated$CellChat_confidence_class <-
    NA_character_

  integrated$CellChat_pathways <-
    NA_character_

  integrated$CellChat_receptors <-
    NA_character_
}


integrated <- integrated %>%
  mutate(
    CellChat_highconf =
      replace_na(
        CellChat_highconf,
        FALSE
      ),

    mechanism_family =
      mechanism_family(
        ligand
      ),

    NicheNet_top20 =
      nichenet_rank <= 20,

    sender_expressed =
      !is.na(
        sender_ligand_pct_expressed
      ) &
      sender_ligand_pct_expressed >=
        EXPRESSION_PCT_THRESHOLD,

    sender_DE_direction_match =
      direction_matches_program(
        receiver_program,
        sender_ligand_logFC
      ),

    sender_DE_supported =
      sender_DE_direction_match &
      !is.na(
        sender_ligand_PValue
      ) &
      sender_ligand_PValue <=
        0.10,

    CellChat_direction_match =
      CellChat_highconf &
      direction_matches_program(
        receiver_program,
        CellChat_best_supported_delta
      ),

    CellChat_replicate_supported =
      CellChat_highconf &
      !is.na(
        CellChat_support_samples
      ) &
      CellChat_support_samples >=
        3 &
      !is.na(
        CellChat_direction_consistency
      ) &
      CellChat_direction_consistency >=
        0.75,

    CellChat_matching_supported =
      CellChat_direction_match &
      CellChat_replicate_supported,

    evidence_count =
      as.integer(
        NicheNet_top20
      ) +
      as.integer(
        sender_DE_supported
      ) +
      as.integer(
        CellChat_matching_supported
      ),

    evidence_class =
      case_when(
        NicheNet_top20 &
          sender_DE_supported &
          CellChat_matching_supported ~
          "3_way_supported",

        NicheNet_top20 &
          CellChat_matching_supported ~
          "NicheNet_plus_CellChat",

        NicheNet_top20 &
          sender_DE_supported ~
          "NicheNet_plus_senderDE",

        sender_DE_supported &
          CellChat_matching_supported ~
          "CellChat_plus_senderDE",

        CellChat_matching_supported ~
          "CellChat_only",

        NicheNet_top20 ~
          "NicheNet_only",

        sender_DE_supported ~
          "senderDE_only",

        TRUE ~
          "No_strict_concordance"
      )
  )


write.csv(
  integrated,
  file.path(
    TAB_OUT,
    "12_integrated_NicheNet_CellChat_sender_evidence_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Primary Injury-Hepatocyte shortlist
# ==============================================================================

injury_shortlist <- integrated %>%
  filter(
    receiver ==
      PRIMARY_RECEIVER
  ) %>%
  arrange(
    receiver_program,
    desc(
      evidence_count
    ),
    nichenet_rank,
    desc(
      abs(
        CellChat_best_supported_delta
      )
    )
  )


write.csv(
  injury_shortlist,
  file.path(
    TAB_OUT,
    "13_InjuryHep_integrated_Mphi_ligand_shortlist_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Primary mechanistic families
# ==============================================================================

core_mechanisms <- integrated %>%
  filter(
    mechanism_family %in%
      PRIMARY_MECHANISM_FAMILIES
  ) %>%
  arrange(
    receiver ==
      PRIMARY_RECEIVER,
    receiver_program,
    mechanism_family,
    desc(
      evidence_count
    ),
    nichenet_rank
  )


write.csv(
  core_mechanisms,
  file.path(
    TAB_OUT,
    "14_core_mechanism_family_integrated_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Strict concordance shortlist
# ==============================================================================

strict_shortlist <- integrated %>%
  filter(
    sender_expressed,
    NicheNet_top20 |
      sender_DE_supported |
      CellChat_matching_supported,
    evidence_count >=
      2
  ) %>%
  arrange(
    receiver ==
      PRIMARY_RECEIVER,
    desc(
      evidence_count
    ),
    nichenet_rank,
    desc(
      abs(
        CellChat_best_supported_delta
      )
    )
  )


write.csv(
  strict_shortlist,
  file.path(
    TAB_OUT,
    "15_strict_concordant_candidate_shortlist_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Injury-Hep strict shortlist
# ==============================================================================

injury_strict <- strict_shortlist %>%
  filter(
    receiver ==
      PRIMARY_RECEIVER
  )


write.csv(
  injury_strict,
  file.path(
    TAB_OUT,
    "16_InjuryHep_strict_concordant_candidates_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 17. Core ligand-target links for Injury-Hep
# ==============================================================================

injury_core_ligands <- injury_shortlist %>%
  filter(
    mechanism_family %in%
      PRIMARY_MECHANISM_FAMILIES,
    nichenet_rank <=
      30
  ) %>%
  pull(
    ligand
  ) %>%
  unique()


injury_core_target_links <- target_links %>%
  filter(
    receiver ==
      PRIMARY_RECEIVER,
    ligand %in%
      injury_core_ligands
  ) %>%
  left_join(
    injury_shortlist %>%
      select(
        receiver,
        receiver_program,
        ligand,
        nichenet_rank,
        aupr_corrected,
        sender,
        sender_ligand_logFC,
        sender_ligand_PValue,
        CellChat_best_supported_delta,
        CellChat_direction_consistency,
        CellChat_support_samples,
        evidence_class,
        evidence_count,
        mechanism_family
      ),
    by = c(
      "receiver",
      "receiver_program",
      "ligand"
    )
  )


write.csv(
  injury_core_target_links,
  file.path(
    TAB_OUT,
    "17_InjuryHep_core_ligand_target_links_v6.6.3.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 18. Focused CellChat axes audit
# ==============================================================================

if (
  nrow(
    cellchat_all
  )
) {

  focused_cellchat <- cellchat_all %>%
    filter(
      direction ==
        "Mphi_to_Hep",
      target ==
        PRIMARY_RECEIVER
    ) %>%
    mutate(
      mechanism_family =
        case_when(
          toupper(
            pathway_name
          ) == "PDGF" ~
            "PDGF",

          toupper(
            pathway_name
          ) == "FN1" ~
            "FN1",

          toupper(
            pathway_name
          ) == "SPP1" ~
            "SPP1",

          toupper(
            pathway_name
          ) == "SEMA4" ~
            "SEMA4",

          toupper(
            pathway_name
          ) == "TNF" ~
            "TNF",

          TRUE ~
            "Other"
        )
    ) %>%
    filter(
      mechanism_family !=
        "Other"
    )

  write.csv(
    focused_cellchat,
    file.path(
      TAB_OUT,
      "18_InjuryHep_focused_CellChat_axes_v6.6.3.csv"
    ),
    row.names = FALSE
  )
}


# ==============================================================================
# 19. Figure 1: receiver pseudobulk DE
# ==============================================================================

receiver_volcano <- receiver_de_all %>%
  mutate(
    neglog10FDR =
      -log10(
        pmax(
          FDR,
          1e-300
        )
      ),

    highlighted =
      abs(
        logFC
      ) >=
        0.5 &
      FDR <=
        0.10
  )


p_volcano <- ggplot(
  receiver_volcano,
  aes(
    x =
      logFC,
    y =
      neglog10FDR
  )
) +
  geom_point(
    aes(
      alpha =
        highlighted
    ),
    size =
      0.7
  ) +
  geom_vline(
    xintercept = c(
      -0.5,
      0.5
    ),
    linewidth =
      0.3
  ) +
  geom_hline(
    yintercept =
      -log10(
        0.10
      ),
    linewidth =
      0.3
  ) +
  facet_wrap(
    ~ receiver,
    nrow =
      1,
    scales =
      "free_y"
  ) +
  scale_alpha_manual(
    values = c(
      "TRUE" =
        0.9,
      "FALSE" =
        0.22
    ),
    guide =
      "none"
  ) +
  labs(
    title =
      "Hepatocyte-state pseudobulk DE | Tx vs Sham",
    subtitle =
      "edgeR QL; positive logFC = higher in Tx; biological n=2 vs 2",
    x =
      "log2FC Tx vs Sham",
    y =
      "-log10 FDR"
  ) +
  theme_classic(
    base_size =
      7.5
  )


save_pdf(
  p_volcano,
  file.path(
    FIG_OUT,
    "01_Hep5_pseudobulk_DE_volcano_v6.6.3.pdf"
  ),
  15,
  4.5
)


# ==============================================================================
# 20. Figure 2: NicheNet ligand activity heatmap
# ==============================================================================

if (
  nrow(
    nichenet_activity
  )
) {

  top_ligand_union <- nichenet_activity %>%
    group_by(
      receiver,
      receiver_program
    ) %>%
    slice_min(
      order_by =
        nichenet_rank,
      n =
        15,
      with_ties =
        FALSE
    ) %>%
    pull(
      test_ligand
    ) %>%
    unique()

  ligand_heat <- nichenet_activity %>%
    filter(
      test_ligand %in%
        top_ligand_union
    ) %>%
    mutate(
      receiver_program_label =
        paste0(
          receiver,
          " | ",
          receiver_program
        )
    )

  p_ligand <- ggplot(
    ligand_heat,
    aes(
      x =
        receiver_program_label,
      y =
        test_ligand,
      fill =
        aupr_corrected
    )
  ) +
    geom_tile(
      linewidth =
        0.2
    ) +
    scale_fill_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) +
    labs(
      title =
        "NicheNet ligand activity across Hepatocyte states",
      x =
        NULL,
      y =
        "Candidate ligand",
      fill =
        "AUPR corrected"
    ) +
    theme_classic(
      base_size =
        6.5
    ) +
    theme(
      axis.text.x =
        element_text(
          angle =
            55,
          hjust =
            1
        )
    )

  save_pdf(
    p_ligand,
    file.path(
      FIG_OUT,
      "02_NicheNet_ligand_activity_heatmap_v6.6.3.pdf"
    ),
    14,
    14
  )
}


# ==============================================================================
# 21. Figure 3: integrated Injury-Hep candidates
# ==============================================================================

injury_plot <- injury_shortlist %>%
  filter(
    nichenet_rank <=
      30 |
    evidence_count >=
      2
  ) %>%
  mutate(
    label =
      paste0(
        sender,
        " | ",
        ligand
      )
  )


if (
  nrow(
    injury_plot
  )
) {

  p_injury <- ggplot(
    injury_plot,
    aes(
      x =
        receiver_program,
      y =
        label,
      size =
        nichenet_percentile,
      fill =
        evidence_count
    )
  ) +
    geom_point(
      shape =
        21
    ) +
    scale_fill_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) +
    labs(
      title =
        "Integrated Mphi ligand -> Injury-Hep candidates",
      subtitle =
        "NicheNet + sender pseudobulk DE + replicate-supported CellChat",
      x =
        NULL,
      y =
        NULL,
      size =
        "NicheNet percentile",
      fill =
        "Evidence count"
    ) +
    theme_classic(
      base_size =
        7
    )

  save_pdf(
    p_injury,
    file.path(
      FIG_OUT,
      "03_integrated_Mphi_ligand_InjuryHep_candidates_v6.6.3.pdf"
    ),
    8,
    max(
      9,
      min(
        20,
        5 +
          0.20 *
            n_distinct(
              injury_plot$label
            )
      )
    )
  )
}


# ==============================================================================
# 22. Figure 4: primary mechanism family summary
# ==============================================================================

core_plot <- core_mechanisms %>%
  filter(
    receiver ==
      PRIMARY_RECEIVER,
    mechanism_family %in%
      PRIMARY_MECHANISM_FAMILIES,
    nichenet_rank <=
      50
  ) %>%
  mutate(
    label =
      paste0(
        sender,
        " | ",
        ligand
      ),
    cellchat_delta_plot =
      replace_na(
        CellChat_best_supported_delta,
        0
      )
  )


if (
  nrow(
    core_plot
  )
) {

  p_core <- ggplot(
    core_plot,
    aes(
      x =
        receiver_program,
      y =
        label,
      size =
        nichenet_percentile,
      fill =
        cellchat_delta_plot
    )
  ) +
    geom_point(
      shape =
        21
    ) +
    facet_wrap(
      ~ mechanism_family,
      scales =
        "free_y",
      ncol =
        4
    ) +
    scale_fill_gradient2(
      low =
        "#0033FF",
      mid =
        "#FFFFFF",
      high =
        "#FF1A1A",
      midpoint =
        0
    ) +
    labs(
      title =
        "Core Mphi -> Injury-Hep mechanistic families",
      subtitle =
        "Point size = NicheNet percentile; fill = supported CellChat Tx-Sham",
      x =
        NULL,
      y =
        NULL,
      size =
        "NicheNet percentile",
      fill =
        "CellChat\nTx-Sham"
    ) +
    theme_classic(
      base_size =
        6.5
    )

  save_pdf(
    p_core,
    file.path(
      FIG_OUT,
      "04_core_Mphi_InjuryHep_mechanisms_v6.6.3.pdf"
    ),
    13,
    10
  )
}


# ==============================================================================
# 23. Save result object
# ==============================================================================

results <- list(
  version =
    "v6.6.3",

  receiver_DE =
    receiver_de_all,

  sender_DE =
    sender_de_all,

  expression_fraction =
    expr_fraction,

  geneset_audit =
    geneset_audit_df,

  sender_receiver_ligand_map =
    sender_receiver_ligand_map,

  nichenet_activity =
    nichenet_activity,

  target_links =
    target_links,

  integrated =
    integrated,

  injury_shortlist =
    injury_shortlist,

  strict_shortlist =
    strict_shortlist,

  injury_core_target_links =
    injury_core_target_links
)


saveRDS(
  results,
  file.path(
    RDS_OUT,
    "Mouse_MASH_Mphi5_Hep5_NicheNet_results_v6.6.3.rds"
  ),
  compress = FALSE
)


# ==============================================================================
# 24. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "CellChat_highconf_input",
    "CellChat_all_refined_input",
    "edgeR_version",
    "nichenetr_version",
    "NicheNet_model",
    "expression_pct_threshold",
    "receiver_DE_design",
    "biological_replicates",
    "primary_receiver",
    "receiver_states",
    "sender_states",
    "primary_receiver_gene_selection",
    "fallback_gene_selection",
    "CellChat_rerun",
    "formal_group_inference_note"
  ),
  value = c(
    "v6.6.3",
    INPUT_RDS,
    CELLCHAT_HIGHCONF,
    CELLCHAT_REFINED_ALL,
    as.character(
      packageVersion(
        "edgeR"
      )
    ),
    as.character(
      packageVersion(
        "nichenetr"
      )
    ),
    "Mouse NicheNet v2 / Zenodo 7074291",
    as.character(
      EXPRESSION_PCT_THRESHOLD
    ),
    "sample-level pseudobulk edgeR QL; Tx vs Sham",
    "Sham n=2; Tx n=2",
    PRIMARY_RECEIVER,
    paste(
      RECEIVERS,
      collapse = " | "
    ),
    paste(
      SENDERS,
      collapse = " | "
    ),
    "FDR<=0.10 and abs(logFC)>=0.5",
    "P<=0.05 abs(logFC)>=0.5, then ranked abs(logFC)>=0.25 if needed",
    "FALSE",
    "Exploratory because biological n=2 per group"
  )
)


write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.6.3.csv"
  ),
  row.names = FALSE
)


capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.6.3.txt"
    )
)


# ==============================================================================
# 25. Console summary
# ==============================================================================

msg(
  "Geneset audit:"
)


print(
  geneset_audit_df
)


msg(
  "Top Injury-Hep integrated candidates:"
)


print(
  injury_shortlist %>%
    select(
      sender,
      receiver_program,
      ligand,
      mechanism_family,
      nichenet_rank,
      aupr_corrected,
      sender_ligand_logFC,
      sender_ligand_PValue,
      CellChat_best_supported_delta,
      CellChat_direction_consistency,
      CellChat_support_samples,
      evidence_count,
      evidence_class
    ) %>%
    arrange(
      desc(
        evidence_count
      ),
      nichenet_rank
    ) %>%
    head(
      40
    )
)


msg(
  "Strict Injury-Hep concordant candidates:"
)


print(
  injury_strict %>%
    select(
      sender,
      receiver_program,
      ligand,
      mechanism_family,
      nichenet_rank,
      sender_ligand_logFC,
      sender_ligand_PValue,
      CellChat_best_supported_delta,
      CellChat_direction_consistency,
      CellChat_support_samples,
      evidence_count,
      evidence_class
    ) %>%
    head(
      40
    )
)


msg(
  "DONE."
)


msg(
  "Output directory: ",
  OUT
)
