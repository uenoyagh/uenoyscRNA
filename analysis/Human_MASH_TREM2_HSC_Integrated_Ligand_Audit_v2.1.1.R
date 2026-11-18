#!/usr/bin/env Rscript

# ============================================================
# Human MASH TREM2-Mphi -> HSC integrated ligand audit v2.1.1
#
# PURPOSE
#   Integrate:
#     1) NicheNet ligand activity ranking (v2.0.1)
#     2) actual TREM2-Mphi abundance-weighted ligand burden
#     3) whole-HSC pseudobulk TGF-beta Hallmark response (v1.9)
#     4) whole-HSC pseudobulk ECM / fibrogenic programs
#     5) sender-state specificity across the five canonical
#        myeloid states
#
# PRIMARY QUESTION
#   Which TREM2-Mphi ligands show convergent support from both
#   NicheNet regulatory potential and sample-level sender/receiver
#   covariation?
#
# IMPORTANT
#   - Biological sample is the unit of inference.
#   - Candidate set is prespecified as PRIMARY NicheNet top 30
#     plus TGFB1.
#   - No causal inference.
#   - No reclustering.
#   - No source RDS overwrite.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(Matrix)
  library(ggplot2)
})

set.seed(20261118)

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

MYELOID_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9",
  "objects",
  "Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

HSC_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_HSC_ECM_Proliferation_v1.2",
  "objects",
  "Human_MASH_HSC_ECM_Proliferation_v1.2.rds"
)

NICHE_TOP <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TREM2_HSC_NicheNet_v2.0.1",
  "tables",
  "12_PRIMARY_top30_ligands.tsv"
)

NICHE_ALL <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TREM2_HSC_NicheNet_v2.0.1",
  "tables",
  "10_NicheNet_ligand_activities_all_genesets.tsv"
)

HSC_HALLMARK <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TGFbeta_Cell_PB_Concordance_v1.9",
  "tables",
  "08_whole_sample_HSC_pseudobulk_Hallmark_score.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TREM2_HSC_Integrated_Ligand_Audit_v2.1.1"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG <- file.path(OUT, "figures")
DIR_OBJ <- file.path(OUT, "objects")
DIR_LOG <- file.path(OUT, "logs")
DIR_PROV <- file.path(OUT, "provenance")

dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OBJ, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PROV, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(
  DIR_LOG,
  "Human_MASH_TREM2_HSC_Integrated_Ligand_Audit_v2.1.1.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH integrated ligand audit v2.1.1\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(file.exists(MYELOID_RDS))
stopifnot(file.exists(HSC_RDS))
stopifnot(file.exists(NICHE_TOP))
stopifnot(file.exists(NICHE_ALL))
stopifnot(file.exists(HSC_HALLMARK))

# ------------------------------------------------------------
# 1. Load
# ------------------------------------------------------------

mphi <- readRDS(MYELOID_RDS)
hsc <- readRDS(HSC_RDS)

DefaultAssay(mphi) <- "RNA"
DefaultAssay(hsc) <- "RNA"

mphi <- NormalizeData(
  mphi,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

hsc <- NormalizeData(
  hsc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

top30 <- read_tsv(
  NICHE_TOP,
  show_col_types = FALSE
)

activity_all <- read_tsv(
  NICHE_ALL,
  show_col_types = FALSE
)

primary_activity <- activity_all %>%
  filter(
    geneset == "PRIMARY_rho0.60_p0.10"
  ) %>%
  transmute(
    ligand = test_ligand,
    NicheNet_rank = rank,
    NicheNet_auroc = auroc,
    NicheNet_aupr = aupr,
    NicheNet_aupr_corrected = aupr_corrected,
    NicheNet_pearson = pearson
  )

hallmark <- read_tsv(
  HSC_HALLMARK,
  show_col_types = FALSE
)

required_mphi <- c(
  "sample_v17",
  "NAS_v17",
  "myeloid_state_FINAL_v19"
)

if (length(setdiff(required_mphi, colnames(mphi[[]]))) > 0) {
  stop("Required myeloid metadata missing.")
}

required_hsc <- c(
  "sample",
  "NAS_score"
)

if (length(setdiff(required_hsc, colnames(hsc[[]]))) > 0) {
  stop("Required HSC metadata missing.")
}

if (!all(c(
  "test_ligand",
  "aupr_corrected",
  "rank"
) %in% colnames(top30))) {
  stop("NicheNet top30 table missing required columns.")
}

# ------------------------------------------------------------
# 2. Candidate ligands = top30 + TGFB1
# ------------------------------------------------------------

candidate_ligands <- unique(
  c(
    top30$test_ligand,
    "TGFB1"
  )
)

candidate_ligands_measurable <- intersect(
  candidate_ligands,
  rownames(mphi)
)

candidate_audit <- tibble(
  ligand = candidate_ligands,
  measurable_in_myeloid =
    ligand %in% rownames(mphi),
  NicheNet_top30 =
    ligand %in% top30$test_ligand
) %>%
  left_join(
    primary_activity,
    by = "ligand"
  )

write_tsv(
  candidate_audit,
  file.path(
    DIR_TABLE,
    "01_candidate_ligand_audit.tsv"
  )
)

if (length(candidate_ligands_measurable) < 10) {
  stop("Too few candidate ligands measurable in myeloid object.")
}

# ------------------------------------------------------------
# 3. Myeloid sender burden for all five states
# ------------------------------------------------------------

CANONICAL_STATES <- c(
  "Resident/Kupffer-like_Mphi",
  "Repair/Resolution_Mphi",
  "Lipid-associated/TREM2_Mphi",
  "Inflammatory/Transitional_Mphi",
  "Inflammatory_Monocyte"
)

mmd <- mphi[[]] %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample = as.character(sample_v17),
    NAS = suppressWarnings(as.numeric(NAS_v17)),
    sender_state = as.character(myeloid_state_FINAL_v19)
  ) %>%
  filter(
    sender_state %in% CANONICAL_STATES
  )

expr <- FetchData(
  mphi,
  vars = candidate_ligands_measurable,
  layer = "data"
) %>%
  rownames_to_column("cell")

cnt <- FetchData(
  mphi,
  vars = candidate_ligands_measurable,
  layer = "counts"
) %>%
  rownames_to_column("cell")

total_n <- mmd %>%
  count(
    sample,
    NAS,
    name = "total_myeloid_n"
  )

state_n <- mmd %>%
  count(
    sample,
    NAS,
    sender_state,
    name = "state_n"
  )

expr_long <- expr %>%
  pivot_longer(
    -cell,
    names_to = "ligand",
    values_to = "expr"
  ) %>%
  left_join(
    mmd,
    by = "cell"
  ) %>%
  filter(
    !is.na(sender_state)
  )

cnt_long <- cnt %>%
  pivot_longer(
    -cell,
    names_to = "ligand",
    values_to = "count"
  ) %>%
  left_join(
    mmd,
    by = "cell"
  ) %>%
  filter(
    !is.na(sender_state)
  )

burden_expr <- expr_long %>%
  group_by(
    sample,
    NAS,
    sender_state,
    ligand
  ) %>%
  summarise(
    ligand_mean =
      mean(
        expr,
        na.rm = TRUE
      ),
    ligand_expr_sum =
      sum(
        expr,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

burden_count <- cnt_long %>%
  group_by(
    sample,
    NAS,
    sender_state,
    ligand
  ) %>%
  summarise(
    ligand_detection =
      mean(
        count > 0,
        na.rm = TRUE
      ),
    ligand_detected_n =
      sum(
        count > 0,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

sender_burden <- burden_expr %>%
  full_join(
    burden_count,
    by = c(
      "sample",
      "NAS",
      "sender_state",
      "ligand"
    )
  ) %>%
  left_join(
    state_n,
    by = c(
      "sample",
      "NAS",
      "sender_state"
    )
  ) %>%
  left_join(
    total_n,
    by = c(
      "sample",
      "NAS"
    )
  ) %>%
  mutate(
    state_fraction =
      state_n /
      total_myeloid_n,

    ligand_burden_expr =
      ligand_expr_sum /
      total_myeloid_n,

    ligand_burden_detection =
      ligand_detected_n /
      total_myeloid_n
  ) %>%
  arrange(
    ligand,
    sample,
    sender_state
  )

write_tsv(
  sender_burden,
  file.path(
    DIR_TABLE,
    "02_candidate_ligand_sender_burden_by_sample_state.tsv"
  )
)

# ------------------------------------------------------------
# 4. Whole-HSC pseudobulk ECM and fibrogenic programs
# ------------------------------------------------------------

ECM_GENES <- c(
  "COL1A1","COL1A2","COL3A1","COL5A1","COL5A2",
  "FN1","SPARC","THBS2","LOXL2","TIMP1",
  "TGFBI","CTHRC1","SERPINE1","POSTN"
)

FIBROGENIC_GENES <- c(
  "COL1A1","COL1A2","COL3A1","FN1",
  "THBS2","POSTN","CTHRC1","LOXL2"
)

hsc_counts <- GetAssayData(
  hsc,
  assay = "RNA",
  layer = "counts"
)

hmd <- hsc[[]]

samples <- sort(
  unique(
    as.character(hmd$sample)
  )
)

pb_list <- lapply(
  samples,
  function(s) {

    idx <- which(
      as.character(hmd$sample) == s
    )

    Matrix::rowSums(
      hsc_counts[
        ,
        idx,
        drop = FALSE
      ]
    )
  }
)

pb_counts <- do.call(
  cbind,
  pb_list
)

rownames(pb_counts) <- rownames(hsc_counts)
colnames(pb_counts) <- samples

pb_lib <- colSums(pb_counts)

pb_logcpm <- log2(
  sweep(
    pb_counts,
    2,
    pb_lib,
    "/"
  ) * 1e6 + 1
)

z_fun <- function(x) {
  sx <- sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / sx
}

score_pb_gene_set <- function(mat, genes) {

  g <- intersect(
    genes,
    rownames(mat)
  )

  if (length(g) < 2) {
    stop("Too few genes available for pseudobulk score.")
  }

  sample_x_gene <- t(
    mat[
      g,
      ,
      drop = FALSE
    ]
  )

  z <- apply(
    sample_x_gene,
    2,
    z_fun
  )

  if (is.vector(z)) {
    z <- matrix(
      z,
      ncol = 1
    )
  }

  rowMeans(
    z,
    na.rm = TRUE
  )
}

ecm_pb_score <- score_pb_gene_set(
  pb_logcpm,
  ECM_GENES
)

fib_pb_score <- score_pb_gene_set(
  pb_logcpm,
  FIBROGENIC_GENES
)

hsc_pb <- tibble(
  sample = samples,
  ECM_pseudobulk_score = ecm_pb_score,
  fibrogenic_pseudobulk_score = fib_pb_score
) %>%
  left_join(
    hmd %>%
      transmute(
        sample = as.character(sample),
        NAS = suppressWarnings(
          as.numeric(NAS_score)
        )
      ) %>%
      distinct(),
    by = "sample"
  ) %>%
  left_join(
    hallmark,
    by = c(
      "sample",
      "NAS"
    )
  )

if (!("Hallmark_whole_HSC_pseudobulk_score" %in% colnames(hsc_pb))) {
  stop("v1.9 whole-HSC Hallmark score column not found.")
}

write_tsv(
  hsc_pb,
  file.path(
    DIR_TABLE,
    "03_HSC_pseudobulk_endpoints_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 5. Statistics helper
# ------------------------------------------------------------

safe_spearman <- function(x, y) {

  keep <- complete.cases(x, y)

  if (
    sum(keep) < 4 ||
    length(unique(x[keep])) < 2 ||
    length(unique(y[keep])) < 2
  ) {
    return(
      tibble(
        n = sum(keep),
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  ct <- suppressWarnings(
    cor.test(
      x[keep],
      y[keep],
      method = "spearman",
      exact = FALSE
    )
  )

  tibble(
    n = sum(keep),
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

# ------------------------------------------------------------
# 6. PRIMARY: TREM2-state ligand burden vs HSC endpoints
# ------------------------------------------------------------

TREM2_STATE <- "Lipid-associated/TREM2_Mphi"

trem2_burden <- sender_burden %>%
  filter(
    sender_state == TREM2_STATE
  ) %>%
  inner_join(
    hsc_pb,
    by = c(
      "sample",
      "NAS"
    )
  )

ENDPOINTS <- c(
  "Hallmark_whole_HSC_pseudobulk_score",
  "ECM_pseudobulk_score",
  "fibrogenic_pseudobulk_score",
  "NAS"
)

primary_assoc <- bind_rows(
  lapply(
    candidate_ligands_measurable,
    function(lig) {

      d <- trem2_burden %>%
        filter(
          ligand == lig
        )

      bind_rows(
        lapply(
          ENDPOINTS,
          function(endpoint) {

            safe_spearman(
              d$ligand_burden_expr,
              d[[endpoint]]
            ) %>%
              mutate(
                ligand = lig,
                sender_state = TREM2_STATE,
                sender_metric =
                  "ligand_burden_expr",
                receiver_endpoint =
                  endpoint,
                .before = 1
              )
          }
        )
      )
    }
  )
) %>%
  group_by(
    receiver_endpoint
  ) %>%
  mutate(
    FDR_within_endpoint =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  ungroup()

write_tsv(
  primary_assoc,
  file.path(
    DIR_TABLE,
    "04_PRIMARY_TREM2_ligand_burden_vs_HSC_endpoints.tsv"
  )
)

# ------------------------------------------------------------
# 7. Integrated ranking against the primary HSC Hallmark endpoint
# ------------------------------------------------------------

hallmark_assoc <- primary_assoc %>%
  filter(
    receiver_endpoint ==
      "Hallmark_whole_HSC_pseudobulk_score"
  ) %>%
  select(
    ligand,
    sample_rho = rho,
    sample_p = p_value,
    sample_FDR = FDR_within_endpoint
  )

integrated <- candidate_audit %>%
  filter(
    measurable_in_myeloid
  ) %>%
  left_join(
    hallmark_assoc,
    by = "ligand"
  ) %>%
  mutate(
    sample_correlation_rank =
      min_rank(
        desc(sample_rho)
      ),

    NicheNet_rank_for_integration =
      ifelse(
        is.na(NicheNet_rank),
        max(
          NicheNet_rank,
          na.rm = TRUE
        ) + 1,
        NicheNet_rank
      ),

    integrated_mean_rank =
      (
        NicheNet_rank_for_integration +
          sample_correlation_rank
      ) / 2
  ) %>%
  arrange(
    integrated_mean_rank,
    desc(sample_rho)
  ) %>%
  mutate(
    integrated_rank =
      row_number()
  )

write_tsv(
  integrated,
  file.path(
    DIR_TABLE,
    "05_INTEGRATED_candidate_ligand_ranking.tsv"
  )
)

# ------------------------------------------------------------
# 8. Top candidates: state specificity
# ------------------------------------------------------------

TOP_N <- 10

top_integrated <- integrated %>%
  slice_head(
    n = min(
      TOP_N,
      nrow(integrated)
    )
  ) %>%
  pull(ligand)

top_for_specificity <- unique(
  c(
    top_integrated,
    "TGFB1"
  )
)

state_specificity <- bind_rows(
  lapply(
    top_for_specificity,
    function(lig) {

      d_lig <- sender_burden %>%
        filter(
          ligand == lig
        ) %>%
        inner_join(
          hsc_pb,
          by = c(
            "sample",
            "NAS"
          )
        )

      d_lig %>%
        group_by(
          sender_state
        ) %>%
        group_modify(
          ~ safe_spearman(
            .x$ligand_burden_expr,
            .x$Hallmark_whole_HSC_pseudobulk_score
          )
        ) %>%
        ungroup() %>%
        mutate(
          ligand = lig,
          .before = 1
        )
    }
  )
) %>%
  group_by(
    ligand
  ) %>%
  mutate(
    FDR_across_states =
      p.adjust(
        p_value,
        method = "BH"
      ),

    state_rho_rank =
      min_rank(
        desc(rho)
      )
  ) %>%
  ungroup()

write_tsv(
  state_specificity,
  file.path(
    DIR_TABLE,
    "06_TOP_ligand_sender_state_specificity.tsv"
  )
)

trem2_specificity_summary <- state_specificity %>%
  filter(
    sender_state == TREM2_STATE
  ) %>%
  select(
    ligand,
    TREM2_rho = rho,
    TREM2_p = p_value,
    TREM2_FDR_across_states = FDR_across_states,
    TREM2_state_rho_rank = state_rho_rank
  ) %>%
  arrange(
    TREM2_state_rho_rank,
    desc(TREM2_rho)
  )

write_tsv(
  trem2_specificity_summary,
  file.path(
    DIR_TABLE,
    "07_TOP_ligand_TREM2_specificity_summary.tsv"
  )
)

# ------------------------------------------------------------
# 9. Leave-one-out for top integrated candidates + TGFB1
# ------------------------------------------------------------

loo_ligands <- unique(
  c(
    integrated %>%
      slice_head(n = 5) %>%
      pull(ligand),
    "TGFB1"
  )
)

loo <- bind_rows(
  lapply(
    loo_ligands,
    function(lig) {

      d <- trem2_burden %>%
        filter(
          ligand == lig
        ) %>%
        arrange(sample)

      bind_rows(
        lapply(
          seq_len(nrow(d)),
          function(i) {

            sub <- d[-i, ]

            safe_spearman(
              sub$ligand_burden_expr,
              sub$Hallmark_whole_HSC_pseudobulk_score
            ) %>%
              mutate(
                ligand = lig,
                omitted_sample =
                  d$sample[i],
                omitted_NAS =
                  d$NAS[i],
                .before = 1
              )
          }
        )
      )
    }
  )
)

write_tsv(
  loo,
  file.path(
    DIR_TABLE,
    "08_TOP_ligand_leave_one_out.tsv"
  )
)

loo_summary <- loo %>%
  group_by(ligand) %>%
  summarise(
    tests = n(),
    all_positive =
      all(
        rho > 0,
        na.rm = TRUE
      ),
    min_rho =
      min(
        rho,
        na.rm = TRUE
      ),
    median_rho =
      median(
        rho,
        na.rm = TRUE
      ),
    max_rho =
      max(
        rho,
        na.rm = TRUE
      ),
    p_lt_0_05_n =
      sum(
        p_value < 0.05,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(median_rho)
  )

write_tsv(
  loo_summary,
  file.path(
    DIR_TABLE,
    "09_TOP_ligand_leave_one_out_summary.tsv"
  )
)

# ------------------------------------------------------------
# 10. Focus comparison: top integrated candidates + TGFB1
# ------------------------------------------------------------

focus <- integrated %>%
  filter(
    ligand %in% unique(
      c(
        integrated %>%
          slice_head(n = 10) %>%
          pull(ligand),
        "TGFB1"
      )
    )
  ) %>%
  left_join(
    trem2_specificity_summary,
    by = "ligand"
  ) %>%
  left_join(
    loo_summary %>%
      select(
        ligand,
        LOO_all_positive = all_positive,
        LOO_min_rho = min_rho,
        LOO_median_rho = median_rho,
        LOO_p_lt_0_05_n = p_lt_0_05_n
      ),
    by = "ligand"
  ) %>%
  arrange(
    integrated_rank
  )

write_tsv(
  focus,
  file.path(
    DIR_TABLE,
    "10_FOCUS_integrated_ligand_comparison.tsv"
  )
)

# ------------------------------------------------------------
# 11. Plot
# ------------------------------------------------------------

plot_df <- integrated %>%
  slice_head(
    n = min(
      20,
      nrow(integrated)
    )
  ) %>%
  mutate(
    ligand = factor(
      ligand,
      levels = rev(ligand)
    )
  )

p <- ggplot(
  plot_df,
  aes(
    x = ligand,
    y = sample_rho
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(
    title = "Integrated TREM2-Mphi ligand audit",
    subtitle = "Sample-level ligand burden vs whole-HSC pseudobulk TGF-beta response",
    x = NULL,
    y = "Spearman rho"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_integrated_top20_sample_rho.pdf"
  ),
  p,
  width = 8,
  height = 7
)

# ------------------------------------------------------------
# 12. Save compact object
# ------------------------------------------------------------

saveRDS(
  list(
    candidate_audit = candidate_audit,
    sender_burden = sender_burden,
    hsc_pb = hsc_pb,
    primary_assoc = primary_assoc,
    integrated = integrated,
    state_specificity = state_specificity,
    loo = loo,
    loo_summary = loo_summary,
    focus = focus
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_TREM2_HSC_Integrated_Ligand_Audit_v2.1.1_results.rds"
  )
)

# ------------------------------------------------------------
# 13. Provenance / integrity
# ------------------------------------------------------------

provenance <- tibble(
  item = c(
    "analysis_version",
    "myeloid_source_RDS",
    "HSC_source_RDS",
    "NicheNet_top30_source",
    "NicheNet_all_activity_source",
    "HSC_Hallmark_source",
    "candidate_definition",
    "primary_sender_state",
    "primary_sender_metric",
    "primary_receiver_endpoint",
    "ECM_gene_n",
    "fibrogenic_gene_n",
    "association_unit",
    "source_overwritten",
    "reclustering"
  ),
  value = c(
    "v2.1.1",
    MYELOID_RDS,
    HSC_RDS,
    NICHE_TOP,
    NICHE_ALL,
    HSC_HALLMARK,
    "NicheNet primary top30 plus TGFB1",
    TREM2_STATE,
    "ligand_expr_sum_in_state / total_myeloid_n",
    "whole-HSC pseudobulk TGF-beta Hallmark score from v1.9",
    as.character(
      length(
        intersect(
          ECM_GENES,
          rownames(hsc)
        )
      )
    ),
    as.character(
      length(
        intersect(
          FIBROGENIC_GENES,
          rownames(hsc)
        )
      )
    ),
    "biological_sample",
    "FALSE",
    "none"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "integrated_ligand_provenance_v2.1.1.tsv"
  )
)

integrity <- tibble(
  item = c(
    "candidate_ligands_requested",
    "candidate_ligands_measurable",
    "all_8_HSC_samples_present",
    "all_5_sender_states_present",
    "TGFB1_in_candidate_set",
    "top_integrated_state_specificity_completed",
    "LOO_completed",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    as.character(
      length(
        candidate_ligands
      )
    ),
    as.character(
      length(
        candidate_ligands_measurable
      )
    ),
    as.character(
      nrow(hsc_pb) == 8
    ),
    as.character(
      length(
        intersect(
          CANONICAL_STATES,
          unique(
            sender_burden$sender_state
          )
        )
      ) == 5
    ),
    as.character(
      "TGFB1" %in%
        candidate_ligands_measurable
    ),
    as.character(
      length(
        unique(
          state_specificity$ligand
        )
      ) >=
        length(
          top_integrated
        )
    ),
    as.character(
      nrow(loo) > 0
    ),
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v2.1.1.tsv"
  )
)

# ------------------------------------------------------------
# 14. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("INTEGRATED LIGAND RANKING\n")
cat("============================================================\n")
print(
  integrated %>%
    select(
      integrated_rank,
      ligand,
      NicheNet_rank,
      NicheNet_aupr_corrected,
      sample_rho,
      sample_p,
      sample_FDR,
      sample_correlation_rank,
      integrated_mean_rank
    ) %>%
    head(31)
)

cat("\n")
cat("============================================================\n")
cat("FOCUS COMPARISON\n")
cat("============================================================\n")
print(focus)

cat("\n")
cat("============================================================\n")
cat("TREM2 SOURCE SPECIFICITY\n")
cat("============================================================\n")
print(trem2_specificity_summary)

cat("\n")
cat("============================================================\n")
cat("LEAVE-ONE-OUT SUMMARY\n")
cat("============================================================\n")
print(loo_summary)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
