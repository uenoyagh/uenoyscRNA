#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(Matrix)
  library(ggplot2)
})

set.seed(20261119)

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HUMAN_MPHI_RDS <- file.path(
  ROOT, "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects", "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

HUMAN_HSC_RDS <- file.path(
  ROOT, "Human_MASH_RDS",
  "Human_MASH_HSC_Receiver_Primary_v1.0",
  "objects", "Human_MASH_HSC_Receiver_Primary_v1.0.rds"
)

MOUSE_MPHI_RDS <- file.path(
  ROOT, "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5", "RDS",
  "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

MOUSE_HSC_RDS <- file.path(
  ROOT, "Mouse_MASH_RDS",
  "Mouse_MASH_HSC_PrimaryClean_v7.0", "objects",
  "Mouse_MASH_HSC_PrimaryClean_v7.0.rds"
)

OUT <- file.path(
  ROOT, "Human_MASH_RDS",
  "Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.2"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG <- file.path(OUT, "figures")
DIR_OBJ <- file.path(OUT, "objects")
DIR_LOG <- file.path(OUT, "logs")
DIR_PROV <- file.path(OUT, "provenance")

for (d in c(DIR_TABLE, DIR_FIG, DIR_OBJ, DIR_LOG, DIR_PROV)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(
  DIR_LOG,
  "Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.2.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human-Mouse MASH Mphi-HSC PDGF cross-species v1.0.2\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")
cat("Patch note: v1.0.2 retains the v1.0.1 biological fixes and restores the canonical Human HSC input path/version (Human_MASH_HSC_Receiver_Primary_v1.0); outputs remain versioned separately.\n\n")

cat("Human myeloid source:\n", HUMAN_MPHI_RDS, "\n\n")
cat("Human HSC source:\n", HUMAN_HSC_RDS, "\n\n")
cat("Mouse Mphi source:\n", MOUSE_MPHI_RDS, "\n\n")
cat("Mouse HSC source:\n", MOUSE_HSC_RDS, "\n\n")

stopifnot(
  file.exists(HUMAN_MPHI_RDS),
  file.exists(HUMAN_HSC_RDS),
  file.exists(MOUSE_MPHI_RDS),
  file.exists(MOUSE_HSC_RDS)
)

get_feature_counts_v5 <- function(
  object,
  features,
  cells = NULL,
  assay_name = "RNA",
  object_label = "object"
) {
  assay_obj <- object[[assay_name]]
  count_layers <- Layers(assay_obj, search = "^counts")
  if (length(count_layers) == 0) stop("No counts layers in ", object_label)

  mats <- lapply(count_layers, function(layer_name) {
    x <- LayerData(assay_obj, layer = layer_name)
    feat <- intersect(features, rownames(x))
    if (length(feat) == 0) return(NULL)
    x <- x[feat, , drop = FALSE]
    if (!is.null(cells)) {
      cc <- intersect(cells, colnames(x))
      if (length(cc) == 0) return(NULL)
      x <- x[, cc, drop = FALSE]
    }
    x
  })

  keep <- !vapply(mats, is.null, logical(1))
  mats <- mats[keep]
  used_layers <- count_layers[keep]
  if (length(mats) == 0) stop("No requested data in ", object_label)

  common_features <- Reduce(intersect, lapply(mats, rownames))
  common_features <- features[features %in% common_features]
  if (length(common_features) == 0) stop("No common features in ", object_label)

  mats <- lapply(mats, function(x) x[common_features, , drop = FALSE])

  all_cells <- unlist(lapply(mats, colnames), use.names = FALSE)
  if (anyDuplicated(all_cells) > 0) {
    stop("Duplicated barcodes across count layers in ", object_label)
  }

  mat <- do.call(cbind, mats)

  if (!is.null(cells)) {
    missing_cells <- setdiff(cells, colnames(mat))
    if (length(missing_cells) > 0) {
      stop(object_label, ": requested cells missing; n=", length(missing_cells))
    }
    mat <- mat[, cells, drop = FALSE]
  }

  list(matrix = mat, layers = used_layers)
}

resolve_ncount <- function(object, cells, assay_name = "RNA") {
  md <- object[[]]
  candidates <- c(paste0("nCount_", assay_name), "nCount_RNA", "nCount")
  cc <- candidates[candidates %in% colnames(md)][1]
  if (length(cc) == 0 || is.na(cc)) return(NULL)
  out <- as.numeric(md[cells, cc, drop = TRUE])
  names(out) <- cells
  out
}

make_cp10k <- function(count_mat, libsize) {
  if (is.null(libsize)) stop("nCount metadata unavailable.")
  libsize <- libsize[colnames(count_mat)]
  if (any(!is.finite(libsize)) || any(libsize <= 0)) stop("Invalid nCount.")
  t(t(count_mat) / libsize) * 1e4
}

program_score <- function(cp10k, genes) {
  g <- intersect(genes, rownames(cp10k))
  if (length(g) == 0) return(rep(NA_real_, ncol(cp10k)))
  Matrix::colMeans(cp10k[g, , drop = FALSE])
}

safe_spearman <- function(x, y, human = FALSE) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) {
    return(tibble(n = sum(ok), rho = NA_real_, p_nominal = NA_real_))
  }
  rho <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  p <- if (human) {
    suppressWarnings(
      tryCatch(
        cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    )
  } else NA_real_
  tibble(n = sum(ok), rho = rho, p_nominal = p)
}

summarise_state_ligand <- function(md, cp10k, ligand_genes) {
  sample_totals <- md %>% count(sample, name = "sample_cells")

  bind_rows(lapply(ligand_genes, function(gene_name) {
    expr <- if (gene_name %in% rownames(cp10k)) {
      as.numeric(cp10k[gene_name, md$cell])
    } else {
      rep(NA_real_, nrow(md))
    }

    md %>%
      mutate(
        gene = gene_name,
        expr = expr,
        detected = expr > 0
      ) %>%
      group_by(sample, sender_state, gene) %>%
      summarise(
        state_cells = n(),
        mean_CP10k = mean(expr, na.rm = TRUE),
        detection_fraction = mean(detected, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(sample_totals, by = "sample") %>%
      mutate(
        state_fraction = state_cells / sample_cells,
        ligand_burden = state_fraction * mean_CP10k
      )
  }))
}

summarise_hsc <- function(
  md, cp10k, receptor_genes,
  core_genes, remodel_genes, cycling_genes
) {

  rec_long <- bind_rows(lapply(receptor_genes, function(gene_name) {
    expr <- if (gene_name %in% rownames(cp10k)) {
      as.numeric(cp10k[gene_name, md$cell])
    } else {
      rep(NA_real_, nrow(md))
    }

    md %>%
      mutate(
        receptor = gene_name,
        expr = expr,
        detected = expr > 0
      ) %>%
      group_by(sample, receiver_state, receptor) %>%
      summarise(
        cells = n(),
        mean_CP10k = mean(expr, na.rm = TRUE),
        detection_fraction = mean(detected, na.rm = TRUE),
        .groups = "drop"
      )
  }))

  md2 <- md %>%
    mutate(
      core_collagen = program_score(cp10k, core_genes),
      ECM_remodeling = program_score(cp10k, remodel_genes),
      cycling = program_score(cp10k, cycling_genes)
    )

  state_outputs <- md2 %>%
    group_by(sample, receiver_state) %>%
    summarise(
      cells = n(),
      core_collagen = mean(core_collagen, na.rm = TRUE),
      ECM_remodeling = mean(ECM_remodeling, na.rm = TRUE),
      cycling = mean(cycling, na.rm = TRUE),
      .groups = "drop"
    )

  state_fractions <- md2 %>%
    count(sample, receiver_state, name = "cells") %>%
    group_by(sample) %>%
    mutate(receiver_fraction = cells / sum(cells)) %>%
    ungroup()

  whole_hsc <- md2 %>%
    group_by(sample) %>%
    summarise(
      HSC_cells = n(),
      wholeHSC_core_collagen = mean(core_collagen, na.rm = TRUE),
      wholeHSC_ECM_remodeling = mean(ECM_remodeling, na.rm = TRUE),
      wholeHSC_cycling = mean(cycling, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    receptors = rec_long,
    state_outputs = state_outputs,
    state_fractions = state_fractions,
    whole_hsc = whole_hsc
  )
}

HUMAN_LIGANDS <- c("PDGFA", "PDGFB")
HUMAN_RECEPTORS <- c("PDGFRA", "PDGFRB")
HUMAN_CORE <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2")
HUMAN_REMODEL <- c(
  "FN1", "SPARC", "THBS2", "LOXL2", "TIMP1",
  "TGFBI", "CTHRC1", "SERPINE1", "POSTN"
)
HUMAN_CYCLING <- c("MKI67", "TOP2A", "CDK1", "CCNB1", "PCNA")

MOUSE_LIGANDS <- c("Pdgfa", "Pdgfb")
MOUSE_RECEPTORS <- c("Pdgfra", "Pdgfrb")
MOUSE_CORE <- c("Col1a1", "Col1a2", "Col3a1", "Col5a1", "Col5a2")
MOUSE_REMODEL <- c(
  "Fn1", "Sparc", "Thbs2", "Loxl2", "Timp1",
  "Tgfbi", "Cthrc1", "Serpine1", "Postn"
)
MOUSE_CYCLING <- c("Mki67", "Top2a", "Cdk1", "Ccnb1", "Pcna")

human_mphi <- readRDS(HUMAN_MPHI_RDS)
human_hsc <- readRDS(HUMAN_HSC_RDS)
mouse_mphi <- readRDS(MOUSE_MPHI_RDS)
mouse_hsc <- readRDS(MOUSE_HSC_RDS)

DefaultAssay(human_mphi) <- "RNA"
DefaultAssay(human_hsc) <- "RNA"
DefaultAssay(mouse_mphi) <- "RNA"
DefaultAssay(mouse_hsc) <- "RNA"

# ---------------- HUMAN ----------------

human_mphi_md <- human_mphi[[]] %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample = as.character(sample_v17),
    sender_state = as.character(myeloid_state_FINAL_v19)
  ) %>%
  filter(grepl("^LT[1-8]$", sample))

human_hsc_md0 <- human_hsc[[]] %>%
  rownames_to_column("cell")

hsc_celltype_candidates <- c(
  "celltype",
  "celltype_for_R8plot_FIXED2",
  "celltype_internal"
)

human_hsc_celltype_col <- hsc_celltype_candidates[
  hsc_celltype_candidates %in% colnames(human_hsc_md0)
][1]

if (length(human_hsc_celltype_col) == 0 || is.na(human_hsc_celltype_col)) {
  stop("Human HSC celltype column not found.")
}

human_hsc_md <- human_hsc_md0 %>%
  transmute(
    cell,
    sample = as.character(sample),
    source_celltype =
      as.character(
        .data[[human_hsc_celltype_col]]
      ),
    receiver_state = case_when(
      source_celltype == "HSC_quiescent" ~
        "qHSC",
      source_celltype == "HSC_fibroblast_mesenchymal" ~
        "aHSC_mesenchymal",
      TRUE ~
        NA_character_
    )
  ) %>%
  filter(
    grepl("^LT[1-8]$", sample),
    !is.na(receiver_state)
  )

human_hsc_state_audit <- human_hsc_md %>%
  count(
    sample,
    receiver_state,
    name = "cells"
  ) %>%
  group_by(sample) %>%
  mutate(
    state_fraction =
      cells / sum(cells)
  ) %>%
  ungroup()

write_tsv(
  human_hsc_state_audit,
  file.path(
    DIR_TABLE,
    "00_Human_HSC_state_reconstruction_audit.tsv"
  )
)

if (
  !all(
    c(
      "qHSC",
      "aHSC_mesenchymal"
    ) %in%
      unique(
        human_hsc_md$receiver_state
      )
  )
) {
  stop(
    "Human HSC state reconstruction did not recover both qHSC and aHSC_mesenchymal."
  )
}

human_samples <- intersect(
  sort(unique(human_mphi_md$sample)),
  sort(unique(human_hsc_md$sample))
)

human_mphi_md <- human_mphi_md %>% filter(sample %in% human_samples)
human_hsc_md <- human_hsc_md %>% filter(sample %in% human_samples)

human_mphi_counts <- get_feature_counts_v5(
  human_mphi,
  HUMAN_LIGANDS,
  human_mphi_md$cell,
  object_label = "human_mphi"
)$matrix

human_hsc_counts <- get_feature_counts_v5(
  human_hsc,
  unique(c(HUMAN_RECEPTORS, HUMAN_CORE, HUMAN_REMODEL, HUMAN_CYCLING)),
  human_hsc_md$cell,
  object_label = "human_hsc"
)$matrix

human_mphi_cp10k <- make_cp10k(
  human_mphi_counts,
  resolve_ncount(human_mphi, human_mphi_md$cell)
)

human_hsc_cp10k <- make_cp10k(
  human_hsc_counts,
  resolve_ncount(human_hsc, human_hsc_md$cell)
)

human_ligand_state <- summarise_state_ligand(
  human_mphi_md, human_mphi_cp10k, HUMAN_LIGANDS
)

human_hsc_sum <- summarise_hsc(
  human_hsc_md, human_hsc_cp10k,
  HUMAN_RECEPTORS, HUMAN_CORE, HUMAN_REMODEL, HUMAN_CYCLING
)

write_tsv(
  human_ligand_state,
  file.path(DIR_TABLE, "01_Human_myeloid_state_PDGF_ligand_burden.tsv")
)

write_tsv(
  human_hsc_sum$receptors,
  file.path(DIR_TABLE, "02_Human_HSC_PDGF_receptors.tsv")
)

write_tsv(
  human_hsc_sum$state_outputs,
  file.path(DIR_TABLE, "03_Human_HSC_state_outputs.tsv")
)

write_tsv(
  human_hsc_sum$state_fractions,
  file.path(DIR_TABLE, "04_Human_HSC_state_fractions.tsv")
)

human_trem2_states <- intersect(
  "Lipid-associated/TREM2_Mphi",
  sort(
    unique(
      human_mphi_md$sender_state
    )
  )
)

if (length(human_trem2_states) != 1) {
  stop(
    "Exact Human TREM2 state 'Lipid-associated/TREM2_Mphi' was not found."
  )
}

human_trem2_state_summary <- human_ligand_state %>%
  filter(
    sender_state ==
      "Lipid-associated/TREM2_Mphi"
  )

write_tsv(
  human_trem2_state_summary,
  file.path(DIR_TABLE, "05_Human_TREM2_associated_PDGF_ligand_burden.tsv")
)

human_trem2_sparsity <- human_trem2_state_summary %>%
  group_by(gene) %>%
  summarise(
    samples_n = n_distinct(sample),
    samples_nonzero_burden =
      sum(
        ligand_burden > 0,
        na.rm = TRUE
      ),
    median_burden =
      median(
        ligand_burden,
        na.rm = TRUE
      ),
    max_burden =
      max(
        ligand_burden,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write_tsv(
  human_trem2_sparsity,
  file.path(
    DIR_TABLE,
    "05b_Human_TREM2_PDGF_sparsity_audit.tsv"
  )
)

human_total_ligand <- human_ligand_state %>%
  group_by(sample, gene) %>%
  summarise(
    total_myeloid_ligand_burden = sum(ligand_burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = gene,
    values_from = total_myeloid_ligand_burden,
    names_prefix = "total_"
  )

human_trem2_ligand <- human_trem2_state_summary %>%
  group_by(sample, gene) %>%
  summarise(
    TREM2_ligand_burden = sum(ligand_burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = gene,
    values_from = TREM2_ligand_burden,
    names_prefix = "TREM2_"
  )

human_receptor_whole <- human_hsc_sum$receptors %>%
  group_by(sample, receptor) %>%
  summarise(
    mean_CP10k = weighted.mean(mean_CP10k, w = cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(sample, receptor, mean_CP10k) %>%
  pivot_wider(
    names_from = receptor,
    values_from = mean_CP10k,
    names_prefix = "HSC_"
  )

human_ahsc <- human_hsc_sum$state_fractions %>%
  filter(receiver_state == "aHSC_mesenchymal") %>%
  select(sample, aHSC_fraction = receiver_fraction)

human_ahsc_outputs <- human_hsc_sum$state_outputs %>%
  filter(receiver_state == "aHSC_mesenchymal") %>%
  select(
    sample,
    aHSC_core_collagen = core_collagen,
    aHSC_ECM_remodeling = ECM_remodeling,
    aHSC_cycling = cycling
  )

human_sample_axis <- human_total_ligand %>%
  full_join(human_trem2_ligand, by = "sample") %>%
  left_join(human_receptor_whole, by = "sample") %>%
  left_join(human_ahsc, by = "sample") %>%
  left_join(human_ahsc_outputs, by = "sample") %>%
  left_join(human_hsc_sum$whole_hsc, by = "sample") %>%
  mutate(
    total_PDGFA_PDGFRA_product = total_PDGFA * HSC_PDGFRA,
    total_PDGFB_PDGFRB_product = total_PDGFB * HSC_PDGFRB,
    TREM2_PDGFA_PDGFRA_product = TREM2_PDGFA * HSC_PDGFRA,
    TREM2_PDGFB_PDGFRB_product = TREM2_PDGFB * HSC_PDGFRB
  )

write_tsv(
  human_sample_axis,
  file.path(DIR_TABLE, "06_Human_sample_PDGF_axis_summary.tsv")
)

if (
  any(
    human_sample_axis$aHSC_fraction <= 0 |
      human_sample_axis$aHSC_fraction >= 1
  )
) {
  warning(
    "One or more Human aHSC fractions are at a boundary (0 or 1); inspect 00_Human_HSC_state_reconstruction_audit.tsv."
  )
}

if (
  length(
    unique(
      human_sample_axis$aHSC_fraction
    )
  ) < 2
) {
  stop(
    "Human aHSC_fraction has no between-sample variation; HSC state reconstruction failed."
  )
}

human_predictors <- c(
  "total_PDGFA",
  "total_PDGFB",
  "TREM2_PDGFA",
  "TREM2_PDGFB",
  "total_PDGFA_PDGFRA_product",
  "total_PDGFB_PDGFRB_product",
  "TREM2_PDGFA_PDGFRA_product",
  "TREM2_PDGFB_PDGFRB_product"
)

human_endpoints <- c(
  "aHSC_fraction",
  "aHSC_core_collagen",
  "aHSC_ECM_remodeling",
  "aHSC_cycling",
  "wholeHSC_core_collagen",
  "wholeHSC_ECM_remodeling",
  "wholeHSC_cycling"
)

human_cor <- bind_rows(lapply(human_predictors, function(pred) {
  bind_rows(lapply(human_endpoints, function(endp) {
    res <- safe_spearman(
      human_sample_axis[[pred]],
      human_sample_axis[[endp]],
      human = TRUE
    )
    tibble(
      species = "Human",
      predictor = pred,
      endpoint = endp,
      n = res$n,
      predictor_nonzero_n =
        sum(
          is.finite(
            human_sample_axis[[pred]]
          ) &
            human_sample_axis[[pred]] > 0
        ),
      endpoint_nonzero_n =
        sum(
          is.finite(
            human_sample_axis[[endp]]
          ) &
            human_sample_axis[[endp]] > 0
        ),
      rho = res$rho,
      p_nominal = res$p_nominal,
      inference = "exploratory_sample_level; nominal_p_unadjusted"
    )
  }))
})) %>%
  mutate(
    p_BH_all56 =
      p.adjust(
        p_nominal,
        method = "BH"
      )
  )

write_tsv(
  human_cor,
  file.path(DIR_TABLE, "07_Human_PDGF_HSC_correlations.tsv")
)

# ---------------- MOUSE ----------------

mouse_mphi_md0 <- mouse_mphi[[]] %>%
  rownames_to_column("cell")

mouse_mphi_sample_col <- c(
  "sample_4group", "sample", "sample_for_annotation"
)
mouse_mphi_sample_col <- mouse_mphi_sample_col[
  mouse_mphi_sample_col %in% colnames(mouse_mphi_md0)
][1]

mouse_mphi_state_col <- c(
  "macrophage_class_Res2_FINAL_v4145_char",
  "macrophage_class_Res2_FINAL_v4145"
)
mouse_mphi_state_col <- mouse_mphi_state_col[
  mouse_mphi_state_col %in% colnames(mouse_mphi_md0)
][1]

if (
  length(mouse_mphi_sample_col) == 0 ||
  length(mouse_mphi_state_col) == 0
) stop("Mouse Mphi metadata columns not found.")

mouse_samples <- c("Sham1", "Sham20", "Tx17", "Tx5")

mouse_mphi_md <- mouse_mphi_md0 %>%
  transmute(
    cell,
    sample = as.character(.data[[mouse_mphi_sample_col]]),
    sender_state = as.character(.data[[mouse_mphi_state_col]])
  ) %>%
  filter(sample %in% mouse_samples)

mouse_hsc_md0 <- mouse_hsc[[]] %>%
  rownames_to_column("cell")

if (
  !"sample_hsc_v610" %in% colnames(mouse_hsc_md0) ||
  !"HSC_state3_v612" %in% colnames(mouse_hsc_md0)
) stop("Mouse HSC sample/state columns not found.")

mouse_hsc_md <- mouse_hsc_md0 %>%
  transmute(
    cell,
    sample = as.character(sample_hsc_v610),
    receiver_state = as.character(HSC_state3_v612)
  ) %>%
  filter(sample %in% mouse_samples)

mouse_mphi_counts <- get_feature_counts_v5(
  mouse_mphi,
  MOUSE_LIGANDS,
  mouse_mphi_md$cell,
  object_label = "mouse_mphi"
)$matrix

mouse_hsc_counts <- get_feature_counts_v5(
  mouse_hsc,
  unique(c(MOUSE_RECEPTORS, MOUSE_CORE, MOUSE_REMODEL, MOUSE_CYCLING)),
  mouse_hsc_md$cell,
  object_label = "mouse_hsc"
)$matrix

mouse_mphi_cp10k <- make_cp10k(
  mouse_mphi_counts,
  resolve_ncount(mouse_mphi, mouse_mphi_md$cell)
)

mouse_hsc_cp10k <- make_cp10k(
  mouse_hsc_counts,
  resolve_ncount(mouse_hsc, mouse_hsc_md$cell)
)

mouse_ligand_state <- summarise_state_ligand(
  mouse_mphi_md, mouse_mphi_cp10k, MOUSE_LIGANDS
)

mouse_hsc_sum <- summarise_hsc(
  mouse_hsc_md, mouse_hsc_cp10k,
  MOUSE_RECEPTORS, MOUSE_CORE, MOUSE_REMODEL, MOUSE_CYCLING
)

write_tsv(
  mouse_ligand_state,
  file.path(DIR_TABLE, "08_Mouse_Mphi_state_PDGF_ligand_burden.tsv")
)

write_tsv(
  mouse_hsc_sum$receptors,
  file.path(DIR_TABLE, "09_Mouse_HSC_PDGF_receptors.tsv")
)

write_tsv(
  mouse_hsc_sum$state_outputs,
  file.path(DIR_TABLE, "10_Mouse_HSC_state_outputs.tsv")
)

write_tsv(
  mouse_hsc_sum$state_fractions,
  file.path(DIR_TABLE, "11_Mouse_HSC_state_fractions.tsv")
)

mouse_total_ligand <- mouse_ligand_state %>%
  group_by(sample, gene) %>%
  summarise(
    total_myeloid_ligand_burden = sum(ligand_burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = gene,
    values_from = total_myeloid_ligand_burden,
    names_prefix = "total_"
  )

mouse_trem2_ligand <- mouse_ligand_state %>%
  filter(sender_state == "Lipid-associated/TREM2-Mphi") %>%
  select(sample, gene, TREM2_ligand_burden = ligand_burden) %>%
  pivot_wider(
    names_from = gene,
    values_from = TREM2_ligand_burden,
    names_prefix = "TREM2_"
  )

mouse_receptor_whole <- mouse_hsc_sum$receptors %>%
  group_by(sample, receptor) %>%
  summarise(
    mean_CP10k = weighted.mean(mean_CP10k, w = cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(sample, receptor, mean_CP10k) %>%
  pivot_wider(
    names_from = receptor,
    values_from = mean_CP10k,
    names_prefix = "HSC_"
  )

mouse_ecm_fraction <- mouse_hsc_sum$state_fractions %>%
  filter(receiver_state == "ECM-activated HSC") %>%
  select(sample, ECM_activated_HSC_fraction = receiver_fraction)

mouse_ecm_outputs <- mouse_hsc_sum$state_outputs %>%
  filter(receiver_state == "ECM-activated HSC") %>%
  select(
    sample,
    ECM_activated_core_collagen = core_collagen,
    ECM_activated_ECM_remodeling = ECM_remodeling,
    ECM_activated_cycling = cycling
  )

mouse_sample_axis <- mouse_total_ligand %>%
  full_join(mouse_trem2_ligand, by = "sample") %>%
  left_join(mouse_receptor_whole, by = "sample") %>%
  left_join(mouse_ecm_fraction, by = "sample") %>%
  left_join(mouse_ecm_outputs, by = "sample") %>%
  left_join(mouse_hsc_sum$whole_hsc, by = "sample") %>%
  mutate(
    group = ifelse(grepl("^Tx", sample), "Tx", "Sham"),
    total_Pdgfa_Pdgfra_product = total_Pdgfa * HSC_Pdgfra,
    total_Pdgfb_Pdgfrb_product = total_Pdgfb * HSC_Pdgfrb,
    TREM2_Pdgfa_Pdgfra_product = TREM2_Pdgfa * HSC_Pdgfra,
    TREM2_Pdgfb_Pdgfrb_product = TREM2_Pdgfb * HSC_Pdgfrb
  )

write_tsv(
  mouse_sample_axis,
  file.path(DIR_TABLE, "12_Mouse_sample_PDGF_axis_summary.tsv")
)

mouse_group_effects <- mouse_sample_axis %>%
  select(
    sample, group,
    total_Pdgfa, total_Pdgfb,
    TREM2_Pdgfa, TREM2_Pdgfb,
    HSC_Pdgfra, HSC_Pdgfrb,
    total_Pdgfa_Pdgfra_product,
    total_Pdgfb_Pdgfrb_product,
    TREM2_Pdgfa_Pdgfra_product,
    TREM2_Pdgfb_Pdgfrb_product,
    ECM_activated_HSC_fraction,
    ECM_activated_core_collagen,
    ECM_activated_ECM_remodeling,
    ECM_activated_cycling,
    wholeHSC_core_collagen,
    wholeHSC_ECM_remodeling,
    wholeHSC_cycling
  ) %>%
  pivot_longer(
    cols = -c(sample, group),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(metric, group) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(mean_value, n_samples)
  ) %>%
  mutate(
    Tx_minus_Sham = mean_value_Tx - mean_value_Sham,
    Tx_over_Sham = ifelse(
      mean_value_Sham == 0,
      NA_real_,
      mean_value_Tx / mean_value_Sham
    ),
    inference = "descriptive_only_n2_per_group"
  )

write_tsv(
  mouse_group_effects,
  file.path(DIR_TABLE, "13_Mouse_PDGF_Sham_vs_Tx_effects.tsv")
)

mouse_predictors <- c(
  "total_Pdgfa",
  "total_Pdgfb",
  "TREM2_Pdgfa",
  "TREM2_Pdgfb",
  "total_Pdgfa_Pdgfra_product",
  "total_Pdgfb_Pdgfrb_product",
  "TREM2_Pdgfa_Pdgfra_product",
  "TREM2_Pdgfb_Pdgfrb_product"
)

mouse_endpoints <- c(
  "ECM_activated_HSC_fraction",
  "ECM_activated_core_collagen",
  "ECM_activated_ECM_remodeling",
  "ECM_activated_cycling",
  "wholeHSC_core_collagen",
  "wholeHSC_ECM_remodeling",
  "wholeHSC_cycling"
)

mouse_cor <- bind_rows(lapply(mouse_predictors, function(pred) {
  bind_rows(lapply(mouse_endpoints, function(endp) {
    res <- safe_spearman(
      mouse_sample_axis[[pred]],
      mouse_sample_axis[[endp]],
      human = FALSE
    )
    tibble(
      species = "Mouse",
      predictor = pred,
      endpoint = endp,
      n = res$n,
      rho = res$rho,
      p_nominal = NA_real_,
      inference = "descriptive_only_n4"
    )
  }))
}))

write_tsv(
  mouse_cor,
  file.path(DIR_TABLE, "14_Mouse_PDGF_HSC_correlations.tsv")
)

human_primary <- human_cor %>%
  filter(
    predictor %in% c(
      "total_PDGFA_PDGFRA_product",
      "TREM2_PDGFA_PDGFRA_product",
      "total_PDGFB_PDGFRB_product",
      "TREM2_PDGFB_PDGFRB_product"
    )
  ) %>%
  mutate(
    axis = ifelse(grepl("PDGFA", predictor), "PDGFA-PDGFRA", "PDGFB-PDGFRB"),
    sender_scope = ifelse(grepl("^TREM2", predictor), "TREM2-associated", "total-myeloid")
  )

mouse_primary <- mouse_cor %>%
  filter(
    predictor %in% c(
      "total_Pdgfa_Pdgfra_product",
      "TREM2_Pdgfa_Pdgfra_product",
      "total_Pdgfb_Pdgfrb_product",
      "TREM2_Pdgfb_Pdgfrb_product"
    )
  ) %>%
  mutate(
    axis = ifelse(grepl("Pdgfa", predictor), "PDGFA-PDGFRA", "PDGFB-PDGFRB"),
    sender_scope = ifelse(grepl("^TREM2", predictor), "TREM2-associated", "total-myeloid")
  )

cross_species <- bind_rows(human_primary, mouse_primary) %>%
  select(
    species, axis, sender_scope, predictor,
    endpoint, n, rho, p_nominal, inference
  )

write_tsv(
  cross_species,
  file.path(DIR_TABLE, "15_PRIMARY_cross_species_PDGF_HSC_summary.tsv")
)

fig_df <- cross_species %>%
  filter(grepl("cycling|fraction|core_collagen", endpoint))

p <- ggplot(
  fig_df,
  aes(x = endpoint, y = rho)
) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2.8) +
  facet_grid(
    species ~ axis + sender_scope,
    scales = "free_x",
    space = "free_x"
  ) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Human-Mouse PDGF-HSC sample-level association",
    x = NULL,
    y = "Spearman rho"
  )

ggsave(
  file.path(DIR_FIG, "01_cross_species_PDGF_HSC_association.pdf"),
  p,
  width = 13,
  height = 7
)

coverage <- tibble(
  species = c(rep("Human", 4), rep("Mouse", 4)),
  compartment = rep(
    c("myeloid_ligand", "myeloid_ligand", "HSC_receptor", "HSC_receptor"),
    2
  ),
  gene = c(
    HUMAN_LIGANDS, HUMAN_RECEPTORS,
    MOUSE_LIGANDS, MOUSE_RECEPTORS
  ),
  measurable = c(
    HUMAN_LIGANDS %in% rownames(human_mphi_cp10k),
    HUMAN_RECEPTORS %in% rownames(human_hsc_cp10k),
    MOUSE_LIGANDS %in% rownames(mouse_mphi_cp10k),
    MOUSE_RECEPTORS %in% rownames(mouse_hsc_cp10k)
  )
)

write_tsv(
  coverage,
  file.path(DIR_PROV, "gene_coverage_v1.0.2.tsv")
)

integrity <- tibble(
  item = c(
    "Human_Mphi_source",
    "Human_HSC_source",
    "Mouse_Mphi_source",
    "Mouse_HSC_source",
    "Human_samples",
    "Human_samples_n",
    "Human_TREM2_associated_states",
    "Human_HSC_state_mapping",
    "Human_PDGF_correlation_multiple_testing",
    "Mouse_samples",
    "Mouse_samples_n",
    "Human_inference",
    "Mouse_group_inference",
    "Mouse_correlation_inference",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    HUMAN_MPHI_RDS,
    HUMAN_HSC_RDS,
    MOUSE_MPHI_RDS,
    MOUSE_HSC_RDS,
    paste(human_samples, collapse = ";"),
    as.character(length(human_samples)),
    paste(human_trem2_states, collapse = ";"),
    "HSC_quiescent=>qHSC;HSC_fibroblast_mesenchymal=>aHSC_mesenchymal",
    "BH correction across 56 exploratory Human predictor-endpoint tests",
    paste(mouse_samples, collapse = ";"),
    as.character(length(mouse_samples)),
    "exploratory_sample_level_n8",
    "descriptive_only_n2_per_group",
    "descriptive_only_n4",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(DIR_PROV, "integrity_summary_v1.0.2.tsv")
)

saveRDS(
  list(
    human_sample_axis = human_sample_axis,
    human_correlations = human_cor,
    mouse_sample_axis = mouse_sample_axis,
    mouse_group_effects = mouse_group_effects,
    mouse_correlations = mouse_cor,
    cross_species = cross_species,
    coverage = coverage
  ),
  file.path(
    DIR_OBJ,
    "Human_Mouse_MASH_Mphi_HSC_PDGF_CrossSpecies_v1.0.2_results.rds"
  )
)

cat("\n============================================================\n")
cat("GENE COVERAGE\n")
cat("============================================================\n")
print(coverage)

cat("\n============================================================\n")
cat("HUMAN HSC STATE RECONSTRUCTION\n")
cat("============================================================\n")
print(human_hsc_state_audit)

cat("\n============================================================\n")
cat("HUMAN TREM2-ASSOCIATED STATES\n")
cat("============================================================\n")
print(human_trem2_states)

cat("\n============================================================\n")
cat("HUMAN TREM2 PDGF SPARSITY\n")
cat("============================================================\n")
print(human_trem2_sparsity)

cat("\n============================================================\n")
cat("HUMAN TOP PDGF-HSC ASSOCIATIONS\n")
cat("============================================================\n")
print(
  human_cor %>%
    arrange(desc(abs(rho))) %>%
    slice_head(n = 20)
)

cat("\n============================================================\n")
cat("MOUSE SHAM vs TX PDGF EFFECTS\n")
cat("============================================================\n")
print(mouse_group_effects)

cat("\n============================================================\n")
cat("MOUSE TOP PDGF-HSC ASSOCIATIONS\n")
cat("============================================================\n")
print(
  mouse_cor %>%
    arrange(desc(abs(rho))) %>%
    slice_head(n = 20)
)

cat("\n============================================================\n")
cat("CROSS-SPECIES PRIMARY SUMMARY\n")
cat("============================================================\n")
print(cross_species)

cat("\nOutput directory:\n")
cat(OUT, "\n\n")
cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
