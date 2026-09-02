#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(6320)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# Mouse MASH scRNA-seq
# Focused validation of the PDGFB macrophage -> HSC axis
#
# Version: v6.3.2
#
# PRIMARY AXIS
#   Repair/Resolution-Mphi
#       -> Pdgfb
#       -> Pdgfra / Pdgfrb
#       -> ECM-activated HSC
#
# INPUTS
#   1) v6.2.0 interaction-ready Seurat object
#   2) v6.3.0 sender / receiver pseudobulk DE tables
#   3) v6.3.1 strict 3-way concordance tables
#
# PURPOSE
#   Validate the leading mechanistic candidate directly from observed expression:
#
#   A. Population abundance
#      Is Repair/Resolution-Mphi more abundant in Tx?
#
#   B. Sender ligand expression
#      Is Pdgfb expression per Repair/Resolution-Mphi cell lower in Tx?
#
#   C. Receiver receptor availability
#      Are Pdgfra / Pdgfrb expressed in ECM-activated HSC?
#
#   D. Receiver target program
#      Do NicheNet-supported Pdgfb target genes show the expected Tx-down pattern
#      in ECM-activated HSC?
#
# IMPORTANT
#   - No CellChat rerun.
#   - No NicheNet rerun.
#   - No reclustering / reintegration / new UMAP.
#   - Biological n = 2 Sham vs n = 2 Tx.
#   - Formal inference is exploratory; sample-level values are shown explicitly.
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

mean_logexpr <- function(
  data_mat,
  genes,
  cells
) {
  genes <- intersect(
    genes,
    rownames(
      data_mat
    )
  )

  if (
    !length(
      genes
    ) ||
    !length(
      cells
    )
  ) {
    return(
      rep(
        NA_real_,
        length(
          genes
        )
      )
    )
  }

  Matrix::rowMeans(
    data_mat[
      genes,
      cells,
      drop = FALSE
    ]
  )
}

pct_expr <- function(
  counts,
  genes,
  cells
) {
  genes <- intersect(
    genes,
    rownames(
      counts
    )
  )

  if (
    !length(
      genes
    ) ||
    !length(
      cells
    )
  ) {
    return(
      rep(
        NA_real_,
        length(
          genes
        )
      )
    )
  }

  Matrix::rowMeans(
    counts[
      genes,
      cells,
      drop = FALSE
    ] >
      0
  )
}

aggregate_counts_vector <- function(
  counts,
  cells
) {
  Matrix::rowSums(
    counts[
      ,
      cells,
      drop = FALSE
    ]
  )
}

calc_cpm_gene <- function(
  pb_counts,
  gene
) {
  if (
    !gene %in%
      rownames(
        pb_counts
      )
  ) {
    return(
      rep(
        NA_real_,
        ncol(
          pb_counts
        )
      )
    )
  }

  lib <- Matrix::colSums(
    pb_counts
  )

  1e6 *
    as.numeric(
      pb_counts[
        gene,
        ,
        drop = TRUE
      ]
    ) /
    pmax(
      lib,
      1
    )
}

row_z <- function(x) {
  s <- stats::sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(
      s
    ) ||
    s ==
      0
  ) {
    return(
      rep(
        0,
        length(
          x
        )
      )
    )
  }

  (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) / s
}

canonical_condition <- function(sample) {
  ifelse(
    grepl(
      "^Tx",
      sample
    ),
    "Tx",
    "Sham"
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

V630 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_NicheNet_v6.3.0"
)

V631 <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "Mphi5_HSC3_NicheNet_concordance_v6.3.1"
)

SENDER_DE_FILE <- file.path(
  V630,
  "Tables",
  "06_sender_DE_all_Mphi_states_v6.3.0.csv"
)

RECEIVER_DE_FILE <- file.path(
  V630,
  "Tables",
  "04_receiver_DE_all_HSC_states_v6.3.0.csv"
)

THREEWAY_FILE <- file.path(
  V631,
  "Tables",
  "02_three_way_supported_ligands_v6.3.1.csv"
)

LRT_FILE <- file.path(
  V631,
  "Tables",
  "10_three_way_ligand_receptor_target_links_v6.3.1.csv"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_Interaction",
  "PDGFB_axis_validation_v6.3.2"
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

RDS_OUT <- file.path(
  OUT,
  "RDS"
)

for (
  d in c(
    OUT,
    TAB_OUT,
    FIG_OUT,
    LOG_OUT,
    RDS_OUT
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

CONDITION_COL <-
  "condition_interaction_v620"

WHOLE_UMAP <-
  "umapRPCA"

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

MPHI5 <- c(
  "Anti-inflammatory-Mphi",
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

HSC3 <- c(
  "qHSC",
  "ECM-activated HSC",
  "Contractile HSC"
)

FOCAL_SENDER <-
  "Repair/Resolution-Mphi"

FOCAL_RECEIVER <-
  "ECM-activated HSC"

LIGAND <-
  "Pdgfb"

RECEPTORS <- c(
  "Pdgfra",
  "Pdgfrb"
)

# These are the top v6.3.1 predicted targets, retained as a fixed validation set.
DEFAULT_TARGETS <- c(
  "Cd24a",
  "Aff3",
  "Sfrp1",
  "Ypel1",
  "Atp1b1",
  "Kcnma1",
  "Itga4",
  "Slfn4",
  "Has2",
  "Il7r"
)

# Higher-priority mesenchymal/HSC-relevant targets for focused display.
FOCUSED_TARGETS <- c(
  "Sfrp1",
  "Kcnma1",
  "Has2",
  "Itga4",
  "Atp1b1"
)

SAMPLE_COLORS <- c(
  "Sham1" = "#0072B2",
  "Sham20" = "#56B4E9",
  "Tx17" = "#D55E00",
  "Tx5" = "#E69F00"
)

CONDITION_COLORS <- c(
  "Sham" = "#0072B2",
  "Tx" = "#D55E00"
)


# ==============================================================================
# 4. Preflight
# ==============================================================================

for (
  f in c(
    INPUT_RDS,
    SENDER_DE_FILE,
    RECEIVER_DE_FILE,
    THREEWAY_FILE,
    LRT_FILE
  )
) {
  if (
    !file.exists(
      f
    )
  ) {
    stop(
      "Required input missing: ",
      f
    )
  }
}

msg(
  "Loading interaction-ready object..."
)

obj <- readRDS(
  INPUT_RDS
)

DefaultAssay(
  obj
) <- "RNA"

obj <- safe_join_rna(
  obj
)

required_meta <- c(
  GROUP_COL,
  SAMPLE_COL,
  CONDITION_COL
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

if (
  !WHOLE_UMAP %in%
    Reductions(
      obj
    )
) {
  stop(
    "UMAP reduction missing: ",
    WHOLE_UMAP
  )
}

counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

data_mat <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)

meta <- obj@meta.data

sender_de <- read.csv(
  SENDER_DE_FILE,
  check.names = FALSE
) %>%
  as_tibble()

receiver_de <- read.csv(
  RECEIVER_DE_FILE,
  check.names = FALSE
) %>%
  as_tibble()

threeway <- read.csv(
  THREEWAY_FILE,
  check.names = FALSE
) %>%
  as_tibble()

lrt <- read.csv(
  LRT_FILE,
  check.names = FALSE
) %>%
  as_tibble()


# ==============================================================================
# 5. Verify the focal 3-way-supported axis
# ==============================================================================

focal_threeway <- threeway %>%
  filter(
    sender ==
      FOCAL_SENDER,
    receiver ==
      FOCAL_RECEIVER,
    ligand ==
      LIGAND
  )

write.csv(
  focal_threeway,
  file.path(
    TAB_OUT,
    "01_PDGFB_threeway_supported_axis_audit_v6.3.2.csv"
  ),
  row.names = FALSE
)

if (
  !nrow(
    focal_threeway
  )
) {
  warning(
    "The focal Repair/Resolution-Mphi -> Pdgfb -> ECM-activated HSC axis ",
    "was not found in the v6.3.1 3-way-supported table."
  )
}


# ==============================================================================
# 6. Repair/Resolution-Mphi abundance by sample
# ==============================================================================

mphi_meta <- meta %>%
  as_tibble(
    rownames = "cell"
  ) %>%
  filter(
    .data[[
      GROUP_COL
    ]] %in%
      MPHI5
  )

mphi_abundance <- mphi_meta %>%
  count(
    .data[[
      SAMPLE_COL
    ]],
    .data[[
      GROUP_COL
    ]],
    name =
      "n_cells"
  ) %>%
  group_by(
    .data[[
      SAMPLE_COL
    ]]
  ) %>%
  mutate(
    total_Mphi5 =
      sum(
        n_cells
      ),
    fraction_within_Mphi5 =
      n_cells /
        total_Mphi5
  ) %>%
  ungroup() %>%
  mutate(
    condition =
      canonical_condition(
        as.character(
          .data[[
            SAMPLE_COL
          ]]
        )
      )
  )

write.csv(
  mphi_abundance,
  file.path(
    TAB_OUT,
    "02_Mphi5_abundance_by_sample_v6.3.2.csv"
  ),
  row.names = FALSE
)

repair_abundance <- mphi_abundance %>%
  filter(
    .data[[
      GROUP_COL
    ]] ==
      FOCAL_SENDER
  )

write.csv(
  repair_abundance,
  file.path(
    TAB_OUT,
    "03_RepairResolution_Mphi_abundance_by_sample_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Pdgfb expression by Mphi subtype and sample
# ==============================================================================

if (
  !LIGAND %in%
    rownames(
      counts
    )
) {
  stop(
    "Pdgfb not found in RNA assay."
  )
}

p_mphi_expr_list <- list()

for (
  smp in SAMPLES
) {

  for (
    ct in MPHI5
  ) {

    cells <- rownames(
      meta
    )[
      as.character(
        meta[[
          SAMPLE_COL
        ]]
      ) ==
        smp &
        as.character(
          meta[[
            GROUP_COL
          ]]
        ) ==
          ct
    ]

    if (
      !length(
        cells
      )
    ) {
      next
    }

    p_mphi_expr_list[[
      paste(
        smp,
        ct,
        sep = "__"
      )
    ]] <- tibble(
      sample =
        smp,
      condition =
        canonical_condition(
          smp
        ),
      Mphi_subtype =
        ct,
      n_cells =
        length(
          cells
        ),
      Pdgfb_pct_expressed =
        as.numeric(
          pct_expr(
            counts,
            LIGAND,
            cells
          )
        ),
      Pdgfb_mean_logexpr =
        as.numeric(
          mean_logexpr(
            data_mat,
            LIGAND,
            cells
          )
        )
    )
  }
}

mphi_pdgfb_expr <- bind_rows(
  p_mphi_expr_list
)

write.csv(
  mphi_pdgfb_expr,
  file.path(
    TAB_OUT,
    "04_Pdgfb_expression_by_Mphi_subtype_sample_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Sample-level pseudobulk CPM for Pdgfb
# ==============================================================================

pb_mphi_list <- list()

for (
  ct in MPHI5
) {

  mat_list <- lapply(
    SAMPLES,
    function(smp) {

      cells <- rownames(
        meta
      )[
        as.character(
          meta[[
            SAMPLE_COL
          ]]
        ) ==
          smp &
          as.character(
            meta[[
              GROUP_COL
            ]]
          ) ==
          ct
      ]

      aggregate_counts_vector(
        counts,
        cells
      )
    }
  )

  pb <- do.call(
    cbind,
    mat_list
  )

  rownames(
    pb
  ) <- rownames(
    counts
  )

  colnames(
    pb
  ) <- SAMPLES

  cpm <- calc_cpm_gene(
    pb,
    LIGAND
  )

  pb_mphi_list[[
    ct
  ]] <- tibble(
    Mphi_subtype =
      ct,
    sample =
      SAMPLES,
    condition =
      canonical_condition(
        SAMPLES
      ),
    Pdgfb_CPM =
      cpm,
    log2_Pdgfb_CPM1 =
      log2(
        cpm + 1
      )
  )
}

mphi_pdgfb_pb <- bind_rows(
  pb_mphi_list
)

write.csv(
  mphi_pdgfb_pb,
  file.path(
    TAB_OUT,
    "05_Pdgfb_pseudobulk_CPM_by_Mphi_subtype_sample_v6.3.2.csv"
  ),
  row.names = FALSE
)

repair_pdgfb_pb <- mphi_pdgfb_pb %>%
  filter(
    Mphi_subtype ==
      FOCAL_SENDER
  )

write.csv(
  repair_pdgfb_pb,
  file.path(
    TAB_OUT,
    "06_RepairResolution_Pdgfb_pseudobulk_CPM_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. Sender DE audit for Pdgfb
# ==============================================================================

sender_pdgfb_de <- sender_de %>%
  filter(
    sender ==
      FOCAL_SENDER,
    gene ==
      LIGAND
  )

write.csv(
  sender_pdgfb_de,
  file.path(
    TAB_OUT,
    "07_RepairResolution_Pdgfb_sender_DE_audit_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. Pdgfra / Pdgfrb expression by HSC state and sample
# ==============================================================================

present_receptors <- intersect(
  RECEPTORS,
  rownames(
    counts
  )
)

if (
  !length(
    present_receptors
  )
) {
  stop(
    "Neither Pdgfra nor Pdgfrb is present in RNA assay."
  )
}

hsc_receptor_list <- list()

for (
  smp in SAMPLES
) {

  for (
    ct in HSC3
  ) {

    cells <- rownames(
      meta
    )[
      as.character(
        meta[[
          SAMPLE_COL
        ]]
      ) ==
        smp &
        as.character(
          meta[[
            GROUP_COL
          ]]
        ) ==
          ct
    ]

    if (
      !length(
        cells
      )
    ) {
      next
    }

    pcts <- pct_expr(
      counts,
      present_receptors,
      cells
    )

    means <- mean_logexpr(
      data_mat,
      present_receptors,
      cells
    )

    hsc_receptor_list[[
      paste(
        smp,
        ct,
        sep = "__"
      )
    ]] <- tibble(
      sample =
        smp,
      condition =
        canonical_condition(
          smp
        ),
      HSC_state =
        ct,
      receptor =
        present_receptors,
      n_cells =
        length(
          cells
        ),
      pct_expressed =
        as.numeric(
          pcts
        ),
      mean_logexpr =
        as.numeric(
          means
        )
    )
  }
}

hsc_receptor_expr <- bind_rows(
  hsc_receptor_list
)

write.csv(
  hsc_receptor_expr,
  file.path(
    TAB_OUT,
    "08_Pdgfra_Pdgfrb_expression_by_HSC_state_sample_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. NicheNet target set from v6.3.1
# ==============================================================================

pdgfb_lrt <- lrt %>%
  filter(
    sender ==
      FOCAL_SENDER,
    receiver ==
      FOCAL_RECEIVER,
    ligand ==
      LIGAND
  ) %>%
  arrange(
    desc(
      LRT_score
    ),
    desc(
      regulatory_potential
    )
  )

write.csv(
  pdgfb_lrt,
  file.path(
    TAB_OUT,
    "09_PDGFB_ECM_HSC_ligand_receptor_target_links_v6.3.2.csv"
  ),
  row.names = FALSE
)

lrt_targets <- unique(
  pdgfb_lrt$target
)

TARGETS <- unique(
  c(
    intersect(
      DEFAULT_TARGETS,
      lrt_targets
    ),
    head(
      lrt_targets,
      10
    )
  )
)

TARGETS <- TARGETS[
  TARGETS %in%
    rownames(
      counts
    )
]

if (
  !length(
    TARGETS
  )
) {
  stop(
    "No PDGFB target genes available in RNA assay."
  )
}

write.csv(
  tibble(
    target =
      TARGETS,
    focused =
      TARGETS %in%
        FOCUSED_TARGETS
  ),
  file.path(
    TAB_OUT,
    "10_PDGFB_target_gene_set_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. ECM-activated HSC target expression by sample
# ==============================================================================

target_expr_list <- list()

for (
  smp in SAMPLES
) {

  cells <- rownames(
    meta
  )[
    as.character(
      meta[[
        SAMPLE_COL
      ]]
    ) ==
      smp &
      as.character(
        meta[[
          GROUP_COL
        ]]
      ) ==
        FOCAL_RECEIVER
  ]

  pcts <- pct_expr(
    counts,
    TARGETS,
    cells
  )

  means <- mean_logexpr(
    data_mat,
    TARGETS,
    cells
  )

  target_expr_list[[
    smp
  ]] <- tibble(
    sample =
      smp,
    condition =
      canonical_condition(
        smp
      ),
    target =
      TARGETS,
    pct_expressed =
      as.numeric(
        pcts
      ),
    mean_logexpr =
      as.numeric(
        means
      )
  )
}

target_expr <- bind_rows(
  target_expr_list
)

write.csv(
  target_expr,
  file.path(
    TAB_OUT,
    "11_ECM_HSC_PDGFB_target_expression_by_sample_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 13. ECM-activated HSC target pseudobulk CPM
# ==============================================================================

pb_hsc_list <- lapply(
  SAMPLES,
  function(smp) {

    cells <- rownames(
      meta
    )[
      as.character(
        meta[[
          SAMPLE_COL
        ]]
      ) ==
        smp &
        as.character(
          meta[[
            GROUP_COL
          ]]
        ) ==
          FOCAL_RECEIVER
    ]

    aggregate_counts_vector(
      counts,
      cells
    )
  }
)

pb_hsc <- do.call(
  cbind,
  pb_hsc_list
)

rownames(
  pb_hsc
) <- rownames(
  counts
)

colnames(
  pb_hsc
) <- SAMPLES

lib_hsc <- Matrix::colSums(
  pb_hsc
)

target_pb <- lapply(
  TARGETS,
  function(g) {

    cpm <- 1e6 *
      as.numeric(
        pb_hsc[
          g,
          ,
          drop = TRUE
        ]
      ) /
      pmax(
        lib_hsc,
        1
      )

    tibble(
      target =
        g,
      sample =
        SAMPLES,
      condition =
        canonical_condition(
          SAMPLES
        ),
      CPM =
        cpm,
      log2_CPM1 =
        log2(
          cpm + 1
        )
    )
  }
) %>%
  bind_rows()

write.csv(
  target_pb,
  file.path(
    TAB_OUT,
    "12_ECM_HSC_PDGFB_target_pseudobulk_CPM_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. Receiver DE audit for PDGFB targets
# ==============================================================================

target_receiver_de <- receiver_de %>%
  filter(
    receiver ==
      FOCAL_RECEIVER,
    gene %in%
      TARGETS
  ) %>%
  arrange(
    FDR,
    PValue
  )

write.csv(
  target_receiver_de,
  file.path(
    TAB_OUT,
    "13_ECM_HSC_PDGFB_target_receiver_DE_audit_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 15. Composite sample-level PDGFB axis table
# ==============================================================================

repair_sample_summary <- repair_abundance %>%
  transmute(
    sample =
      as.character(
        .data[[
          SAMPLE_COL
        ]]
      ),
    condition,
    RepairResolution_n =
      n_cells,
    RepairResolution_fraction_Mphi5 =
      fraction_within_Mphi5
  ) %>%
  left_join(
    mphi_pdgfb_expr %>%
      filter(
        Mphi_subtype ==
          FOCAL_SENDER
      ) %>%
      select(
        sample,
        Pdgfb_pct_expressed,
        Pdgfb_mean_logexpr
      ),
    by =
      "sample"
  ) %>%
  left_join(
    repair_pdgfb_pb %>%
      select(
        sample,
        Pdgfb_CPM,
        log2_Pdgfb_CPM1
      ),
    by =
      "sample"
  )

ecm_receptor_summary <- hsc_receptor_expr %>%
  filter(
    HSC_state ==
      FOCAL_RECEIVER
  ) %>%
  select(
    sample,
    receptor,
    receptor_pct_expressed =
      pct_expressed,
    receptor_mean_logexpr =
      mean_logexpr
  ) %>%
  pivot_wider(
    names_from =
      receptor,
    values_from = c(
      receptor_pct_expressed,
      receptor_mean_logexpr
    )
  )

axis_sample_summary <- repair_sample_summary %>%
  left_join(
    ecm_receptor_summary,
    by =
      "sample"
  )

write.csv(
  axis_sample_summary,
  file.path(
    TAB_OUT,
    "14_PDGFB_axis_sample_summary_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 16. Figure 1:
#     Repair/Resolution abundance vs Pdgfb expression
# ==============================================================================

p_abundance <- ggplot(
  repair_sample_summary,
  aes(
    x =
      sample,
    y =
      100 *
        RepairResolution_fraction_Mphi5,
    fill =
      condition
  )
) +
  geom_col(
    width = 0.7
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Repair/Resolution-Mphi abundance",
    x = NULL,
    y =
      "% of Mphi5",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p_pdgfb_pct <- ggplot(
  repair_sample_summary,
  aes(
    x =
      sample,
    y =
      100 *
        Pdgfb_pct_expressed,
    fill =
      condition
  )
) +
  geom_col(
    width = 0.7
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Pdgfb-positive Repair/Resolution-Mphi",
    x = NULL,
    y =
      "% Pdgfb+ cells",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p_pdgfb_mean <- ggplot(
  repair_sample_summary,
  aes(
    x =
      sample,
    y =
      Pdgfb_mean_logexpr,
    fill =
      condition
  )
) +
  geom_col(
    width = 0.7
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Pdgfb mean expression per Repair/Resolution-Mphi",
    x = NULL,
    y =
      "Mean log-normalized expression",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p_pdgfb_pb <- ggplot(
  repair_sample_summary,
  aes(
    x =
      sample,
    y =
      log2_Pdgfb_CPM1,
    fill =
      condition
  )
) +
  geom_col(
    width = 0.7
  ) +
  scale_fill_manual(
    values =
      CONDITION_COLORS
  ) +
  labs(
    title =
      "Repair/Resolution-Mphi Pdgfb pseudobulk",
    x = NULL,
    y =
      "log2(CPM+1)",
    fill =
      "Condition"
  ) +
  theme_classic(
    base_size = 9
  )

p1 <- (
  p_abundance +
    p_pdgfb_pct
) /
  (
    p_pdgfb_mean +
      p_pdgfb_pb
  ) +
  plot_annotation(
    title =
      "Repair/Resolution-Mphi abundance and Pdgfb output | v6.3.2"
  )

save_pdf(
  p1,
  file.path(
    FIG_OUT,
    "01_RepairResolution_abundance_vs_Pdgfb_output_v6.3.2.pdf"
  ),
  12,
  9
)


# ==============================================================================
# 17. Figure 2:
#     Pdgfb across all macrophage subtypes by sample
# ==============================================================================

mphi_plot_df <- mphi_pdgfb_expr %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      ),
    Mphi_subtype =
      factor(
        Mphi_subtype,
        levels =
          MPHI5
      )
  )

p2 <- ggplot(
  mphi_plot_df,
  aes(
    x =
      sample,
    y =
      Mphi_subtype,
    size =
      Pdgfb_pct_expressed,
    fill =
      Pdgfb_mean_logexpr
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
  scale_size(
    range = c(
      1,
      9
    )
  ) +
  labs(
    title =
      "Pdgfb expression across macrophage subtypes",
    subtitle =
      "Point size = fraction expressing; fill = mean log-normalized expression",
    x = NULL,
    y = NULL,
    size =
      "Fraction expressing",
    fill =
      "Mean expression"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p2,
  file.path(
    FIG_OUT,
    "02_Pdgfb_Mphi5_by_sample_DotPlot_v6.3.2.pdf"
  ),
  10,
  6
)


# ==============================================================================
# 18. Figure 3:
#     Pdgfra/Pdgfrb in HSC states by sample
# ==============================================================================

receptor_plot_df <- hsc_receptor_expr %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      ),
    HSC_state =
      factor(
        HSC_state,
        levels =
          HSC3
      )
  )

p3 <- ggplot(
  receptor_plot_df,
  aes(
    x =
      sample,
    y =
      HSC_state,
    size =
      pct_expressed,
    fill =
      mean_logexpr
  )
) +
  geom_point(
    shape = 21
  ) +
  facet_wrap(
    ~ receptor,
    ncol = 2
  ) +
  scale_fill_gradientn(
    colours = c(
      "#0033FF",
      "#FFFFFF",
      "#FF1A1A"
    )
  ) +
  scale_size(
    range = c(
      1,
      9
    )
  ) +
  labs(
    title =
      "PDGF receptor availability across HSC states",
    subtitle =
      "Pdgfra / Pdgfrb expression in Sham1, Sham20, Tx17, Tx5",
    x = NULL,
    y = NULL,
    size =
      "Fraction expressing",
    fill =
      "Mean expression"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p3,
  file.path(
    FIG_OUT,
    "03_Pdgfra_Pdgfrb_HSC3_by_sample_DotPlot_v6.3.2.pdf"
  ),
  11,
  6
)


# ==============================================================================
# 19. Figure 4:
#     PDGFB target pseudobulk heatmap in ECM-activated HSC
# ==============================================================================

target_heat <- target_pb %>%
  select(
    target,
    sample,
    log2_CPM1
  ) %>%
  pivot_wider(
    names_from =
      sample,
    values_from =
      log2_CPM1
  )

target_heat_mat <- as.matrix(
  target_heat[
    ,
    SAMPLES,
    drop = FALSE
  ]
)

rownames(
  target_heat_mat
) <- target_heat$target

target_z <- t(
  apply(
    target_heat_mat,
    1,
    row_z
  )
)

target_heat_long <- as.data.frame(
  target_z
) %>%
  rownames_to_column(
    "target"
  ) %>%
  pivot_longer(
    cols =
      all_of(
        SAMPLES
      ),
    names_to =
      "sample",
    values_to =
      "z"
  ) %>%
  mutate(
    sample =
      factor(
        sample,
        levels =
          SAMPLES
      ),
    target =
      factor(
        target,
        levels =
          rev(
            TARGETS
          )
      )
  )

p4 <- ggplot(
  target_heat_long,
  aes(
    x =
      sample,
    y =
      target,
    fill =
      z
  )
) +
  geom_tile(
    linewidth = 0.35
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0
  ) +
  labs(
    title =
      "ECM-activated HSC | predicted PDGFB target program",
    subtitle =
      "Sample-level pseudobulk log2(CPM+1), row z-score",
    x = NULL,
    y = NULL,
    fill =
      "Row z-score"
  ) +
  theme_classic(
    base_size = 9
  )

save_pdf(
  p4,
  file.path(
    FIG_OUT,
    "04_ECM_HSC_PDGFB_target_pseudobulk_heatmap_v6.3.2.pdf"
  ),
  8,
  7
)


# ==============================================================================
# 20. Figure 5:
#     Focused target sample-level pseudobulk values
# ==============================================================================

focused_present <- intersect(
  FOCUSED_TARGETS,
  TARGETS
)

if (
  length(
    focused_present
  )
) {

  focused_df <- target_pb %>%
    filter(
      target %in%
        focused_present
    ) %>%
    mutate(
      sample =
        factor(
          sample,
          levels =
            SAMPLES
        )
    )

  p5 <- ggplot(
    focused_df,
    aes(
      x =
        sample,
      y =
        log2_CPM1,
      group =
        1,
      fill =
        condition
    )
  ) +
    geom_col(
      width = 0.7
    ) +
    facet_wrap(
      ~ target,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(
      values =
        CONDITION_COLORS
    ) +
    labs(
      title =
        "Focused predicted PDGFB targets in ECM-activated HSC",
      subtitle =
        "Biological samples shown explicitly; no cell-level pseudoreplication",
      x = NULL,
      y =
        "log2(CPM+1)",
      fill =
        "Condition"
    ) +
    theme_classic(
      base_size = 9
    )

  save_pdf(
    p5,
    file.path(
      FIG_OUT,
      "05_ECM_HSC_focused_PDGFB_targets_by_sample_v6.3.2.pdf"
    ),
    11,
    7
  )
}


# ==============================================================================
# 21. Figure 6:
#     Pdgfb FeaturePlot on interaction-ready whole-cell UMAP
# ==============================================================================

if (
  LIGAND %in%
    rownames(
      obj
    )
) {

  p6 <- FeaturePlot(
    obj,
    features =
      LIGAND,
    reduction =
      WHOLE_UMAP,
    min.cutoff =
      "q05",
    max.cutoff =
      "q95",
    order =
      TRUE,
    pt.size =
      0.25,
    raster =
      FALSE
  ) +
    scale_color_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    ) +
    ggtitle(
      "Pdgfb | 5 Mphi + 3 HSC interaction-ready object"
    )

  save_pdf(
    p6,
    file.path(
      FIG_OUT,
      "06_Pdgfb_FeaturePlot_interaction_ready_v6.3.2.pdf"
    ),
    8,
    7
  )
}


# ==============================================================================
# 22. Figure 7:
#     Pdgfra/Pdgfrb FeaturePlot on HSC receiver cells
# ==============================================================================

hsc_cells <- rownames(
  meta
)[
  as.character(
    meta[[
      GROUP_COL
    ]]
  ) %in%
    HSC3
]

hsc_obj <- subset(
  obj,
  cells =
    hsc_cells
)

if (
  length(
    present_receptors
  )
) {

  p7 <- FeaturePlot(
    hsc_obj,
    features =
      present_receptors,
    reduction =
      WHOLE_UMAP,
    min.cutoff =
      "q05",
    max.cutoff =
      "q95",
    order =
      TRUE,
    pt.size =
      0.35,
    raster =
      FALSE,
    ncol =
      length(
        present_receptors
      )
  ) &
    scale_color_gradientn(
      colours = c(
        "#0033FF",
        "#FFFFFF",
        "#FF1A1A"
      )
    )

  save_pdf(
    p7,
    file.path(
      FIG_OUT,
      "07_Pdgfra_Pdgfrb_FeaturePlot_HSC3_v6.3.2.pdf"
    ),
    12,
    6
  )
}


# ==============================================================================
# 23. Summary table for manuscript-level interpretation
# ==============================================================================

sender_logFC <- if (
  nrow(
    sender_pdgfb_de
  )
) {
  sender_pdgfb_de$logFC[[1]]
} else {
  NA_real_
}

sender_P <- if (
  nrow(
    sender_pdgfb_de
  )
) {
  sender_pdgfb_de$PValue[[1]]
} else {
  NA_real_
}

sender_FDR <- if (
  nrow(
    sender_pdgfb_de
  )
) {
  sender_pdgfb_de$FDR[[1]]
} else {
  NA_real_
}

mean_sham_frac <- repair_sample_summary %>%
  filter(
    condition ==
      "Sham"
  ) %>%
  summarise(
    x =
      mean(
        RepairResolution_fraction_Mphi5
      )
  ) %>%
  pull(
    x
  )

mean_tx_frac <- repair_sample_summary %>%
  filter(
    condition ==
      "Tx"
  ) %>%
  summarise(
    x =
      mean(
        RepairResolution_fraction_Mphi5
      )
  ) %>%
  pull(
    x
  )

mean_sham_pdgfb <- repair_sample_summary %>%
  filter(
    condition ==
      "Sham"
  ) %>%
  summarise(
    x =
      mean(
        Pdgfb_mean_logexpr
      )
  ) %>%
  pull(
    x
  )

mean_tx_pdgfb <- repair_sample_summary %>%
  filter(
    condition ==
      "Tx"
  ) %>%
  summarise(
    x =
      mean(
        Pdgfb_mean_logexpr
      )
  ) %>%
  pull(
    x
  )

summary_table <- tibble(
  metric = c(
    "RepairResolution_Mphi_fraction_mean_Sham",
    "RepairResolution_Mphi_fraction_mean_Tx",
    "RepairResolution_fraction_Tx_minus_Sham",
    "RepairResolution_Pdgfb_mean_logexpr_Sham",
    "RepairResolution_Pdgfb_mean_logexpr_Tx",
    "RepairResolution_Pdgfb_mean_logexpr_Tx_minus_Sham",
    "RepairResolution_Pdgfb_sender_logFC_Tx_vs_Sham",
    "RepairResolution_Pdgfb_sender_PValue",
    "RepairResolution_Pdgfb_sender_FDR",
    "Pdgfra_present_in_RNA",
    "Pdgfrb_present_in_RNA",
    "n_PDGFB_targets_validated",
    "threeway_axis_present"
  ),
  value = c(
    mean_sham_frac,
    mean_tx_frac,
    mean_tx_frac -
      mean_sham_frac,
    mean_sham_pdgfb,
    mean_tx_pdgfb,
    mean_tx_pdgfb -
      mean_sham_pdgfb,
    sender_logFC,
    sender_P,
    sender_FDR,
    "Pdgfra" %in%
      rownames(
        counts
      ),
    "Pdgfrb" %in%
      rownames(
        counts
      ),
    length(
      TARGETS
    ),
    nrow(
      focal_threeway
    ) >
      0
  )
)

write.csv(
  summary_table,
  file.path(
    TAB_OUT,
    "15_PDGFB_axis_validation_summary_v6.3.2.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 24. Save compact RDS
# ==============================================================================

results <- list(
  focal_threeway =
    focal_threeway,
  RepairResolution_abundance =
    repair_abundance,
  Pdgfb_expression_Mphi5 =
    mphi_pdgfb_expr,
  Pdgfb_pseudobulk_Mphi5 =
    mphi_pdgfb_pb,
  sender_Pdgfb_DE =
    sender_pdgfb_de,
  HSC_receptors =
    hsc_receptor_expr,
  Pdgfb_LRT =
    pdgfb_lrt,
  targets =
    TARGETS,
  target_expression =
    target_expr,
  target_pseudobulk =
    target_pb,
  target_receiver_DE =
    target_receiver_de,
  axis_sample_summary =
    axis_sample_summary,
  validation_summary =
    summary_table
)

saveRDS(
  results,
  file.path(
    RDS_OUT,
    "Mouse_MASH_PDGFB_axis_validation_results_v6.3.2.rds"
  ),
  compress = FALSE
)


# ==============================================================================
# 25. Manifest
# ==============================================================================

manifest <- tibble(
  parameter = c(
    "version",
    "input_RDS",
    "focal_sender",
    "ligand",
    "focal_receiver",
    "receptors",
    "samples",
    "biological_replicates",
    "CellChat_recomputed",
    "NicheNet_recomputed",
    "formal_inference_note"
  ),
  value = c(
    "v6.3.2",
    INPUT_RDS,
    FOCAL_SENDER,
    LIGAND,
    FOCAL_RECEIVER,
    paste(
      RECEPTORS,
      collapse = ","
    ),
    paste(
      SAMPLES,
      collapse = ","
    ),
    "Sham n=2; Tx n=2",
    "FALSE",
    "FALSE",
    "Exploratory; sample-level values shown explicitly"
  )
)

write.csv(
  manifest,
  file.path(
    LOG_OUT,
    "analysis_manifest_v6.3.2.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file =
    file.path(
      LOG_OUT,
      "sessionInfo_v6.3.2.txt"
    )
)

msg(
  "DONE."
)

msg(
  "Output directory: ",
  OUT
)
