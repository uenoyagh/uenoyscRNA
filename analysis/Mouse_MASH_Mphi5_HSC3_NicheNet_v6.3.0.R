#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6300)

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
# NicheNet validation of macrophage -> HSC communication
#
# Version: v6.3.0
#
# INPUTS
#   1) v6.2.0 interaction-ready Seurat object
#      5 macrophage subtypes + 3 HSC states
#
#   2) v6.2.2 refined CellChat high-confidence LR shortlist
#
# ANALYTIC DESIGN
#
#   Biological replicates:
#     Sham1, Sham20, Tx17, Tx5
#
#   Sender populations:
#     Anti-inflammatory-Mphi
#     Inflammatory-Mphi
#     ECM-associated inflammatory-Mphi
#     Repair/Resolution-Mphi
#     Lipid-associated/TREM2-Mphi
#
#   Receiver populations:
#     qHSC
#     ECM-activated HSC
#     Contractile HSC
#
# PRIMARY QUESTION
#   Which macrophage-derived ligands can explain transcriptional changes in
#   each HSC receiver state after transplantation?
#
# IMPORTANT DESIGN FEATURES
#
#   - Receiver DE is sample-level pseudobulk, NOT cell-level pseudoreplication.
#   - edgeR QL is used with n=2 Sham and n=2 Tx; results are exploratory.
#   - Tx-up and Tx-down receiver gene programs are analyzed separately.
#   - NicheNet ligand activity is integrated with:
#       a) sender ligand expression
#       b) sender pseudobulk ligand Tx/Sham logFC
#       c) v6.2.2 CellChat replicate support
#   - No arbitrary CellChat log2 ratio is used.
#   - Mouse NicheNet v2 prior model is used directly.
#
# INTERPRETATION
#
#   Tx_up_program:
#     Ligands whose activity may explain genes increased in HSC after Tx.
#
#   Tx_down_program:
#     Ligands whose reduced activity in Tx / higher activity in Sham may explain
#     genes that decrease in HSC after Tx.
#
#   NicheNet ligand activity itself does NOT prove ligand direction.
#   Sender DE and CellChat direction are retained as independent evidence.
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
    file.exists(
      dest
    ) &&
    file.info(
      dest
    )$size >
      1000
  ) {
    msg(
      "Using cached model: ",
      dest
    )
    return(
      invisible(
        dest
      )
    )
  }

  msg(
    "Downloading: ",
    url
  )

  tmp <- paste0(
    dest,
    ".tmp"
  )

  if (
    file.exists(
      tmp
    )
  ) {
    unlink(
      tmp
    )
  }

  utils::download.file(
    url = url,
    destfile = tmp,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )

  if (
    !file.exists(
      tmp
    ) ||
    file.info(
      tmp
    )$size <=
      1000
  ) {
    stop(
      "Model download failed or file too small: ",
      url
    )
  }

  file.rename(
    tmp,
    dest
  )

  invisible(
    dest
  )
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
    function(smp) {

      cells <- rownames(
        meta
      )[
        as.character(
          meta[[
            celltype_col
          ]]
        ) ==
          celltype_value &
          as.character(
            meta[[
              sample_col
            ]]
          ) ==
            smp
      ]

      if (
        length(
          cells
        ) ==
          0
      ) {
        stop(
          "No cells for ",
          celltype_value,
          " / ",
          smp
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

  rownames(
    mat
  ) <- rownames(
    counts
  )

  colnames(
    mat
  ) <- samples

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
    sum(
      keep
    ) <
      100
  ) {
    stop(
      "Too few genes passed edgeR filterByExpr: ",
      sum(
        keep
      )
    )
  }

  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  y <- edgeR::calcNormFactors(
    y
  )

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

  tab <- edgeR::topTags(
    qlf,
    n = Inf,
    sort.by = "PValue"
  )$table %>%
    as.data.frame() %>%
    rownames_to_column(
      "gene"
    ) %>%
    as_tibble()

  tab
}

pct_expressed <- function(
  counts,
  cells
) {

  if (
    !length(
      cells
    )
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
    x >
      0
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
    direction ==
      "Tx_up_program"
  ) {

    direction_primary <- de_table %>%
      filter(
        logFC >=
          0.5
      )

    direction_fallback <- de_table %>%
      filter(
        logFC >=
          0.25
      )

  } else if (
    direction ==
      "Tx_down_program"
  ) {

    direction_primary <- de_table %>%
      filter(
        logFC <=
          -0.5
      )

    direction_fallback <- de_table %>%
      filter(
        logFC <=
          -0.25
      )

  } else {

    stop(
      "Unknown direction: ",
      direction
    )
  }

  primary <- direction_primary %>%
    filter(
      FDR <=
        0.10
    ) %>%
    filter(
      gene %in%
        ligand_target_genes
    ) %>%
    arrange(
      PValue
    )

  if (
    nrow(
      primary
    ) >=
      min_genes
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
        PValue <=
          0.05
      ) %>%
      filter(
        gene %in%
          ligand_target_genes
      ) %>%
      arrange(
        PValue
      )

    if (
      nrow(
        relaxed
      ) >=
        min_genes
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
        ) %>%
        slice_head(
          n = min(
            max_genes,
            max(
              min_genes,
              min(
                100,
                nrow(
                  direction_fallback
                )
              )
            )
          )
        )

      selected <- ranked

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
    length(
      geneset
    ) <
      10
  ) {
    warning(
      "NicheNet skipped: fewer than 10 receiver genes in geneset."
    )
    return(
      tibble()
    )
  }

  if (
    length(
      potential_ligands
    ) <
      2
  ) {
    warning(
      "NicheNet skipped: fewer than 2 potential ligands."
    )
    return(
      tibble()
    )
  }

  out <- nichenetr::predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix =
      ligand_target_matrix,
    potential_ligands =
      potential_ligands
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
          n() >
            1,
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

  out
}

extract_top_target_links <- function(
  ligand_activities,
  geneset_table,
  ligand_target_matrix,
  receiver,
  receiver_program,
  top_ligands = 15,
  top_targets_per_ligand = 20
) {

  if (
    !nrow(
      ligand_activities
    )
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

  out <- lapply(
    ligands_use,
    function(lig) {

      if (
        !lig %in%
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
        !length(
          targets
        )
      ) {
        return(
          tibble()
        )
      }

      rp <- ligand_target_matrix[
        targets,
        lig
      ]

      tmp <- tibble(
        target =
          names(
            rp
          ),
        regulatory_potential =
          as.numeric(
            rp
          )
      ) %>%
        filter(
          is.finite(
            regulatory_potential
          ),
          regulatory_potential >
            0
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
            lig,
          receiver =
            receiver,
          receiver_program =
            receiver_program
        )

      tmp
    }
  ) %>%
    bind_rows()

  out
}

mechanism_family <- function(gene) {

  case_when(
    gene ==
      "Fn1" ~
      "FN1",

    gene %in%
      c(
        "Pdgfa",
        "Pdgfb",
        "Pdgfc",
        "Pdgfd"
      ) ~
      "PDGF",

    gene %in%
      c(
        "Sema4a",
        "Sema4b",
        "Sema4c",
        "Sema4d",
        "Sema4f",
        "Sema4g"
      ) ~
      "SEMA4",

    gene ==
      "App" ~
      "APP",

    gene ==
      "Thbs1" ~
      "THBS",

    gene ==
      "Spp1" ~
      "SPP1",

    gene %in%
      c(
        "Tgfb1",
        "Tgfb2",
        "Tgfb3"
      ) ~
      "TGFB",

    gene %in%
      c(
        "Tnf"
      ) ~
      "TNF",

    TRUE ~
      "Other"
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
  "Mphi5_HSC3_interaction_ready_v6.2.0",
  "RDS",
  "Mouse_MASH_Mphi5_HSC3_interaction_ready_v6.2.0.rds"
)

CELLCHAT_HIGHCONF <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_CellChat_refined_v6.2.2",
  "Tables",
  "13_high_confidence_LR_candidate_shortlist_v6.2.2.csv"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_NicheNet_v6.3.0"
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
  "interaction_celltype_v620"

SAMPLE_COL <-
  "sample_interaction_v620"

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

SENDERS <- c(
  "Anti-inflammatory-Mphi",
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

RECEIVERS <- c(
  "qHSC",
  "ECM-activated HSC",
  "Contractile HSC"
)

RECEIVER_PROGRAMS <- c(
  "Tx_up_program",
  "Tx_down_program"
)

EXPRESSION_PCT_THRESHOLD <-
  0.10

MIN_GENESET <- 20

MAX_GENESET <- 200

TOP_LIGANDS_FOR_TARGETS <- 15

TOP_TARGETS_PER_LIGAND <- 20


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
  "Loading NicheNet mouse model..."
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
    "01_NicheNet_mouse_model_audit_v6.3.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 5. Load interaction-ready object
# ==============================================================================

if (
  !file.exists(
    INPUT_RDS
  )
) {
  stop(
    "Input v6.2.0 RDS not found: ",
    INPUT_RDS
  )
}

msg(
  "Loading v6.2.0 interaction-ready object..."
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

counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

meta <- obj@meta.data

msg(
  "Cells: ",
  ncol(
    obj
  )
)

msg(
  "Genes: ",
  nrow(
    obj
  )
)


# ==============================================================================
# 6. CellChat high-confidence evidence
# ==============================================================================

if (
  file.exists(
    CELLCHAT_HIGHCONF
  )
) {

  cellchat <- read.csv(
    CELLCHAT_HIGHCONF,
    check.names = FALSE
  ) %>%
    as_tibble()

  required_cc <- c(
    "source",
    "target",
    "ligand",
    "receptor",
    "delta_prob_supported",
    "pairwise_direction_consistency",
    "total_support_replicates",
    "confidence_class"
  )

  missing_cc <- setdiff(
    required_cc,
    colnames(
      cellchat
    )
  )

  if (
    length(
      missing_cc
    )
  ) {
    warning(
      "CellChat high-confidence table missing columns: ",
      paste(
        missing_cc,
        collapse = ", "
      ),
      ". CellChat integration will be partial."
    )
  }

} else {

  warning(
    "CellChat v6.2.2 high-confidence table not found. ",
    "NicheNet will run, but integrated evidence will lack CellChat support."
  )

  cellchat <- tibble()
}


# ==============================================================================
# 7. Expression fractions for all senders / receivers
# ==============================================================================

msg(
  "Computing expression fractions..."
)

expr_fraction_list <- list()

for (
  ct in c(
    SENDERS,
    RECEIVERS
  )
) {

  cells <- rownames(
    meta
  )[
    as.character(
      meta[[
        GROUP_COL
      ]]
    ) ==
      ct
  ]

  pct <- pct_expressed(
    counts,
    cells
  )

  expr_fraction_list[[
    ct
  ]] <- tibble(
    gene =
      names(
        pct
      ),
    celltype =
      ct,
    pct_expressed =
      as.numeric(
        pct
      ),
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
    "02_expression_fraction_by_celltype_v6.3.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Pseudobulk DE for HSC receivers
# ==============================================================================

msg(
  "Running HSC receiver pseudobulk DE..."
)

receiver_de_list <- list()

for (
  receiver in RECEIVERS
) {

  msg(
    "Receiver DE: ",
    receiver
  )

  pb <- aggregate_counts_by_sample(
    counts = counts,
    meta = meta,
    celltype_col = GROUP_COL,
    sample_col = SAMPLE_COL,
    celltype_value = receiver,
    samples = SAMPLES
  )

  de <- run_edger_ql(
    pb_counts = pb,
    samples = SAMPLES
  ) %>%
    mutate(
      receiver =
        receiver,
      comparison =
        "Tx_vs_Sham",
      interpretation =
        "positive_logFC_means_higher_in_Tx"
    )

  receiver_de_list[[
    receiver
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
          receiver
        ),
        "_Tx_vs_Sham_v6.3.0.csv"
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
    "04_receiver_DE_all_HSC_states_v6.3.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. Pseudobulk DE for macrophage senders
# ==============================================================================

msg(
  "Running macrophage sender pseudobulk DE..."
)

sender_de_list <- list()

for (
  sender in SENDERS
) {

  msg(
    "Sender DE: ",
    sender
  )

  pb <- aggregate_counts_by_sample(
    counts = counts,
    meta = meta,
    celltype_col = GROUP_COL,
    sample_col = SAMPLE_COL,
    celltype_value = sender,
    samples = SAMPLES
  )

  de <- run_edger_ql(
    pb_counts = pb,
    samples = SAMPLES
  ) %>%
    mutate(
      sender =
        sender,
      comparison =
        "Tx_vs_Sham",
      interpretation =
        "positive_logFC_means_higher_in_Tx"
    )

  sender_de_list[[
    sender
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
          sender
        ),
        "_Tx_vs_Sham_v6.3.0.csv"
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
    "06_sender_DE_all_Mphi_states_v6.3.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. NicheNet per receiver and receiver-program direction
# ==============================================================================

msg(
  "Running NicheNet ligand activity analysis..."
)

nichenet_activity_all <- list()

geneset_audit <- list()

target_link_all <- list()

sender_receiver_ligand_map_all <- list()

for (
  receiver in RECEIVERS
) {

  msg(
    "============================================================"
  )

  msg(
    "NicheNet receiver: ",
    receiver
  )

  receiver_expr <- expr_fraction %>%
    filter(
      celltype ==
        receiver,
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

  # --------------------------------------------------------------------------
  # Potential ligand map by macrophage sender
  # --------------------------------------------------------------------------

  sender_ligand_map <- lapply(
    SENDERS,
    function(sender) {

      sender_expr <- expr_fraction %>%
        filter(
          celltype ==
            sender,
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
        !nrow(
          lr_use
        )
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
            sender,
          receiver =
            receiver
        )
    }
  ) %>%
    bind_rows()

  sender_receiver_ligand_map_all[[
    receiver
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
          receiver
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

    geneset_audit[[
      paste(
        receiver,
        receiver_program,
        sep = "__"
      )
    ]] <- tibble(
      receiver =
        receiver,
      receiver_program =
        receiver_program,
      geneset_method =
        selected$method,
      n_geneset =
        length(
          geneset
        ),
      n_background =
        length(
          background
        ),
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
            receiver
          ),
          "_",
          receiver_program,
          "_v6.3.0.csv"
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
      !nrow(
        ligand_activity
      )
    ) {
      next
    }

    ligand_activity <- ligand_activity %>%
      mutate(
        receiver =
          receiver,
        receiver_program =
          receiver_program,
        geneset_method =
          selected$method,
        n_geneset =
          length(
            geneset
          ),
        inferred_activity_direction =
          ifelse(
            receiver_program ==
              "Tx_up_program",
            "Candidate_Tx_active_ligand",
            "Candidate_Sham_active_or_Tx_lost_ligand"
          )
      )

    nichenet_activity_all[[
      paste(
        receiver,
        receiver_program,
        sep = "__"
      )
    ]] <-
      ligand_activity

    links <- extract_top_target_links(
      ligand_activities =
        ligand_activity,
      geneset_table =
        selected$table,
      ligand_target_matrix =
        ligand_target_matrix,
      receiver =
        receiver,
      receiver_program =
        receiver_program,
      top_ligands =
        TOP_LIGANDS_FOR_TARGETS,
      top_targets_per_ligand =
        TOP_TARGETS_PER_LIGAND
    )

    if (
      nrow(
        links
      )
    ) {
      target_link_all[[
        paste(
          receiver,
          receiver_program,
          sep = "__"
        )
      ]] <-
        links
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
    "08_NicheNet_geneset_audit_v6.3.0.csv"
  ),
  row.names = FALSE
)

write.csv(
  sender_receiver_ligand_map,
  file.path(
    TAB_OUT,
    "09_sender_receiver_potential_ligand_receptor_map_v6.3.0.csv"
  ),
  row.names = FALSE
)

write.csv(
  nichenet_activity,
  file.path(
    TAB_OUT,
    "10_NicheNet_ligand_activity_all_receivers_v6.3.0.csv"
  ),
  row.names = FALSE
)

write.csv(
  target_links,
  file.path(
    TAB_OUT,
    "11_NicheNet_top_ligand_target_links_v6.3.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Expand NicheNet ligand activity to sender-specific candidates
# ==============================================================================

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

# Add sender pseudobulk ligand DE.
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

# Add expression fractions.
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
# 12. Add refined CellChat evidence
# ==============================================================================

if (
  nrow(
    cellchat
  ) &&
  all(
    c(
      "source",
      "target",
      "ligand"
    ) %in%
      colnames(
        cellchat
      )
  )
) {

  cellchat_summary <- cellchat %>%
    mutate(
      ligand =
        as.character(
          ligand
        )
    ) %>%
    group_by(
      source,
      target,
      ligand
    ) %>%
    arrange(
      desc(
        abs(
          delta_prob_supported
        )
      ),
      .by_group = TRUE
    ) %>%
    summarise(
      CellChat_supported =
        TRUE,

      CellChat_best_delta_prob =
        first(
          delta_prob_supported
        ),

      CellChat_direction_consistency =
        max(
          pairwise_direction_consistency,
          na.rm = TRUE
        ),

      CellChat_support_replicates =
        max(
          total_support_replicates,
          na.rm = TRUE
        ),

      CellChat_confidence_class =
        first(
          confidence_class
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

  if (
    "pathway_name" %in%
      colnames(
        cellchat
      )
  ) {

    cellchat_pathways <- cellchat %>%
      group_by(
        source,
        target,
        ligand
      ) %>%
      summarise(
        CellChat_pathways =
          paste(
            unique(
              pathway_name
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

    cellchat_summary <- cellchat_summary %>%
      left_join(
        cellchat_pathways,
        by = c(
          "sender",
          "receiver",
          "ligand"
        )
      )

  } else {

    cellchat_summary$CellChat_pathways <-
      NA_character_
  }

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

  integrated$CellChat_supported <-
    FALSE

  integrated$CellChat_best_delta_prob <-
    NA_real_

  integrated$CellChat_direction_consistency <-
    NA_real_

  integrated$CellChat_support_replicates <-
    NA_real_

  integrated$CellChat_confidence_class <-
    NA_character_

  integrated$CellChat_receptors <-
    NA_character_

  integrated$CellChat_pathways <-
    NA_character_
}

integrated <- integrated %>%
  mutate(
    CellChat_supported =
      replace_na(
        CellChat_supported,
        FALSE
      ),

    mechanism_family =
      mechanism_family(
        ligand
      ),

    expected_sender_direction =
      ifelse(
        receiver_program ==
          "Tx_up_program",
        "Tx_up",
        "Tx_down"
      ),

    sender_DE_direction_match =
      case_when(

        is.na(
          sender_ligand_logFC
        ) ~
          FALSE,

        receiver_program ==
          "Tx_up_program" &
          sender_ligand_logFC >
            0 ~
          TRUE,

        receiver_program ==
          "Tx_down_program" &
          sender_ligand_logFC <
            0 ~
          TRUE,

        TRUE ~
          FALSE
      ),

    CellChat_direction_match =
      case_when(

        !CellChat_supported |
          is.na(
            CellChat_best_delta_prob
          ) ~
          FALSE,

        receiver_program ==
          "Tx_up_program" &
          CellChat_best_delta_prob >
            0 ~
          TRUE,

        receiver_program ==
          "Tx_down_program" &
          CellChat_best_delta_prob <
            0 ~
          TRUE,

        TRUE ~
          FALSE
      ),

    NicheNet_top20 =
      nichenet_rank <=
        20,

    CellChat_replicate_supported =
      CellChat_supported &
        !is.na(
          CellChat_support_replicates
        ) &
        CellChat_support_replicates >=
          2 &
        !is.na(
          CellChat_direction_consistency
        ) &
        CellChat_direction_consistency >=
          0.75,

    sender_DE_supported_direction =
      sender_DE_direction_match &
        !is.na(
          sender_ligand_PValue
        ) &
        sender_ligand_PValue <=
          0.10,

    evidence_count =
      as.integer(
        NicheNet_top20
      ) +
        as.integer(
          CellChat_replicate_supported &
            CellChat_direction_match
        ) +
        as.integer(
          sender_DE_supported_direction
        )
  ) %>%
  arrange(
    receiver,
    receiver_program,
    desc(
      evidence_count
    ),
    nichenet_rank
  )

write.csv(
  integrated,
  file.path(
    TAB_OUT,
    "12_integrated_NicheNet_CellChat_sender_evidence_v6.3.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. Mechanistic shortlist
# ==============================================================================

mechanistic_shortlist <- integrated %>%
  filter(
    NicheNet_top20 |
      evidence_count >=
        2 |
      mechanism_family !=
        "Other"
  ) %>%
  mutate(
    primary_axis = case_when(

      mechanism_family %in%
        c(
          "PDGF",
          "FN1",
          "SEMA4",
          "APP"
        ) ~
        TRUE,

      TRUE ~
        FALSE
    )
  ) %>%
  arrange(
    desc(
      evidence_count
    ),
    desc(
      primary_axis
    ),
    nichenet_rank
  )

write.csv(
  mechanistic_shortlist,
  file.path(
    TAB_OUT,
    "13_mechanistic_ligand_shortlist_v6.3.0.csv"
  ),
  row.names = FALSE
)

primary_axis_only <- mechanistic_shortlist %>%
  filter(
    mechanism_family %in%
      c(
        "PDGF",
        "FN1",
        "SEMA4",
        "APP",
        "THBS",
        "SPP1",
        "TGFB",
        "TNF"
      )
  )

write.csv(
  primary_axis_only,
  file.path(
    TAB_OUT,
    "14_priority_axes_PDGF_FN1_SEMA4_APP_etc_v6.3.0.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Attach sender / CellChat evidence to ligand-target links
# ==============================================================================

if (
  nrow(
    target_links
  )
) {

  target_links_integrated <- target_links %>%
    inner_join(
      integrated %>%
        select(
          sender,
          receiver,
          receiver_program,
          ligand,
          nichenet_rank,
          aupr_corrected,
          evidence_count,
          mechanism_family,
          sender_ligand_logFC,
          sender_ligand_PValue,
          sender_ligand_FDR,
          CellChat_supported,
          CellChat_best_delta_prob,
          CellChat_direction_consistency,
          CellChat_support_replicates,
          CellChat_confidence_class
        ),
      by = c(
        "receiver",
        "receiver_program",
        "ligand"
      )
    ) %>%
    arrange(
      receiver,
      receiver_program,
      desc(
        evidence_count
      ),
      nichenet_rank,
      desc(
        regulatory_potential
      )
    )

  write.csv(
    target_links_integrated,
    file.path(
      TAB_OUT,
      "15_integrated_ligand_target_links_v6.3.0.csv"
    ),
    row.names = FALSE
  )

} else {

  target_links_integrated <-
    tibble()
}


# ==============================================================================
# 15. Figure 1: receiver pseudobulk DE
# ==============================================================================

volcano_df <- receiver_de_all %>%
  mutate(
    neglog10FDR =
      -log10(
        pmax(
          FDR,
          1e-300
        )
      ),

    DE_class =
      case_when(

        FDR <=
          0.10 &
          logFC >=
            0.5 ~
          "Tx_up",

        FDR <=
          0.10 &
          logFC <=
            -0.5 ~
          "Tx_down",

        TRUE ~
          "NS"
      )
  )

p_volcano <- ggplot(
  volcano_df,
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
        DE_class !=
          "NS"
    ),
    size = 0.8
  ) +
  facet_wrap(
    ~ receiver,
    ncol = 3,
    scales = "free_y"
  ) +
  geom_vline(
    xintercept = c(
      -0.5,
      0.5
    ),
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept =
      -log10(
        0.10
      ),
    linewidth = 0.3
  ) +
  scale_alpha_manual(
    values = c(
      "TRUE" = 0.9,
      "FALSE" = 0.25
    ),
    guide = "none"
  ) +
  labs(
    title =
      "HSC-state pseudobulk DE | Tx vs Sham",
    subtitle =
      "edgeR QL; positive logFC = higher in Tx; biological n=2 vs 2",
    x =
      "log2FC Tx vs Sham",
    y =
      "-log10 FDR"
  ) +
  theme_classic(
    base_size = 8
  )

save_pdf(
  p_volcano,
  file.path(
    FIG_OUT,
    "01_HSC3_pseudobulk_DE_volcano_v6.3.0.pdf"
  ),
  12,
  4.5
)


# ==============================================================================
# 16. Figure 2: NicheNet ligand activity heatmap
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
      n = 15,
      with_ties = FALSE
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
      linewidth = 0.2
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
        "NicheNet ligand activity across HSC states",
      x = NULL,
      y =
        "Candidate ligand",
      fill =
        "AUPR corrected"
    ) +
    theme_classic(
      base_size = 7
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        )
    )

  save_pdf(
    p_ligand,
    file.path(
      FIG_OUT,
      "02_NicheNet_ligand_activity_heatmap_v6.3.0.pdf"
    ),
    11,
    12
  )
}


# ==============================================================================
# 17. Figure 3: integrated mechanistic candidates
# ==============================================================================

integrated_plot_df <- mechanistic_shortlist %>%
  filter(
    nichenet_rank <=
      25 |
      evidence_count >=
        2
  ) %>%
  mutate(
    label =
      paste0(
        sender,
        " | ",
        ligand
      ),
    receiver_program_label =
      paste0(
        receiver,
        " | ",
        receiver_program
      )
  )

if (
  nrow(
    integrated_plot_df
  )
) {

  p_integrated <- ggplot(
    integrated_plot_df,
    aes(
      x =
        receiver_program_label,
      y =
        label,
      size =
        nichenet_percentile,
      fill =
        evidence_count
    )
  ) +
    geom_point(
      shape = 21
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
        "Integrated NicheNet + CellChat + sender-DE candidates",
      subtitle =
        "Evidence count: NicheNet top20 + replicate-supported matching CellChat + sender-DE direction",
      x = NULL,
      y = NULL,
      size =
        "NicheNet percentile",
      fill =
        "Evidence count"
    ) +
    theme_classic(
      base_size = 7
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        )
    )

  save_pdf(
    p_integrated,
    file.path(
      FIG_OUT,
      "03_integrated_Mphi_ligand_HSC_program_candidates_v6.3.0.pdf"
    ),
    12,
    max(
      8,
      min(
        20,
        4 +
          0.18 *
            n_distinct(
              integrated_plot_df$label
            )
      )
    )
  )
}


# ==============================================================================
# 18. Figure 4: priority-axis summary
# ==============================================================================

priority_plot_df <- primary_axis_only %>%
  filter(
    mechanism_family %in%
      c(
        "PDGF",
        "FN1",
        "SEMA4",
        "APP",
        "THBS",
        "SPP1",
        "TGFB",
        "TNF"
      )
  ) %>%
  mutate(
    label =
      paste0(
        sender,
        " -> ",
        receiver,
        " | ",
        ligand
      )
  )

if (
  nrow(
    priority_plot_df
  )
) {

  p_priority <- ggplot(
    priority_plot_df,
    aes(
      x =
        sender_ligand_logFC,
      y =
        reorder(
          label,
          sender_ligand_logFC
        ),
      size =
        nichenet_percentile,
      fill =
        evidence_count
    )
  ) +
    geom_point(
      shape = 21
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.3
    ) +
    facet_wrap(
      ~ receiver_program,
      scales = "free_y",
      ncol = 2
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
        "Priority mechanistic axes",
      subtitle =
        "Sender ligand pseudobulk logFC with NicheNet and CellChat evidence",
      x =
        "Sender ligand log2FC Tx vs Sham",
      y = NULL,
      size =
        "NicheNet percentile",
      fill =
        "Evidence count"
    ) +
    theme_classic(
      base_size = 7
    )

  save_pdf(
    p_priority,
    file.path(
      FIG_OUT,
      "04_priority_axes_PDGF_FN1_SEMA4_APP_v6.3.0.pdf"
    ),
    13,
    12
  )
}


# ==============================================================================
# 19. Save compact results RDS
# ==============================================================================

results <- list(
  receiver_DE =
    receiver_de_all,
  sender_DE =
    sender_de_all,
  geneset_audit =
    geneset_audit_df,
  potential_ligand_receptor_map =
    sender_receiver_ligand_map,
  ligand_activity =
    nichenet_activity,
  ligand_target_links =
    target_links,
  integrated_candidates =
    integrated,
  mechanistic_shortlist =
    mechanistic_shortlist,
  priority_axes =
    primary_axis_only
)

saveRDS(
  results,
  file.path(
    RDS_OUT,
    "Mouse_MASH_Mphi5_HSC3_NicheNet_results_v6.3.0.rds"
  ),
  compress = FALSE
)


# ==============================================================================
# 20. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "CellChat_input",
    "edgeR_version",
    "nichenetr_version",
    "NicheNet_model",
    "expression_pct_threshold",
    "receiver_DE_design",
    "biological_replicates",
    "primary_receiver_gene_selection",
    "fallback_gene_selection",
    "formal_group_inference_note"
  ),
  value = c(
    "v6.3.0",
    INPUT_RDS,
    CELLCHAT_HIGHCONF,
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
    "FDR<=0.10 and abs(logFC)>=0.5",
    "P<=0.05 abs(logFC)>=0.5, then ranked abs(logFC)>=0.25 if needed",
    "Exploratory because biological n=2 per group"
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.3.0.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.3.0.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
