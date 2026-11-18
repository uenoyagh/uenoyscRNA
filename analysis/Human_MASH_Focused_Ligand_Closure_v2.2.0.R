#!/usr/bin/env Rscript

# ============================================================
# Human MASH focused ligand closure audit v2.2.0
#
# PRIMARY FOCUS:
#   SPP1  = strongest integrated mechanistic candidate
#   CDCP1 = strongest sample-level candidate / novel candidate
#   TGFB1 = strongest robust sample-level TGF-related candidate
#
# PURPOSE
#   Compare these three ligands on identical criteria:
#     1. TREM2-Mphi abundance-weighted ligand burden
#     2. exact 8! permutation association with HSC endpoints
#     3. partial Spearman controlling NAS
#     4. leave-one-sample-out robustness
#     5. burden decomposition
#     6. sender-state specificity
#     7. HSC receptor support from NicheNet LR network
#     8. NicheNet predicted target support in whole-HSC pseudobulk
#
# IMPORTANT
#   Biological sample is the unit of inference.
#   No reclustering.
#   Source RDS files are never overwritten.
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

if (!requireNamespace("nichenetr", quietly = TRUE)) {
  stop("Package 'nichenetr' is required. It was used successfully in v2.0.1.")
}

set.seed(20261118)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

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

HSC_ENDPOINTS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TREM2_HSC_Integrated_Ligand_Audit_v2.1.1",
  "tables",
  "03_HSC_pseudobulk_endpoints_by_sample.tsv"
)

NICHE_RECEIVER_ASSOC <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TREM2_HSC_NicheNet_v2.0.1",
  "tables",
  "01_wholeHSC_gene_association_with_TREM2_TGFB1.tsv"
)

NICHE_PRIMARY_GENESET <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TREM2_HSC_NicheNet_v2.0.1",
  "tables",
  "03_PRIMARY_receiver_geneset.tsv"
)

NICHE_ACTIVITY_ALL <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TREM2_HSC_NicheNet_v2.0.1",
  "tables",
  "10_NicheNet_ligand_activities_all_genesets.tsv"
)

LR_FILE <- file.path(
  ROOT,
  "resources",
  "NicheNet_v2_human",
  "lr_network_human_21122021.rds"
)

LT_FILE <- file.path(
  ROOT,
  "resources",
  "NicheNet_v2_human",
  "ligand_target_matrix_nsga2r_final.rds"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_Focused_Ligand_Closure_v2.2"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG   <- file.path(OUT, "figures")
DIR_OBJ   <- file.path(OUT, "objects")
DIR_LOG   <- file.path(OUT, "logs")
DIR_PROV  <- file.path(OUT, "provenance")

dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OBJ, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PROV, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(
  DIR_LOG,
  "Human_MASH_Focused_Ligand_Closure_v2.2.0.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH focused ligand closure audit v2.2.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

for (f in c(
  MYELOID_RDS,
  HSC_RDS,
  HSC_ENDPOINTS,
  NICHE_PRIMARY_GENESET,
  NICHE_ACTIVITY_ALL,
  LR_FILE,
  LT_FILE
)) {
  stopifnot(file.exists(f))
}

# ------------------------------------------------------------
# 1. Constants
# ------------------------------------------------------------

FOCUS_LIGANDS <- c(
  "SPP1",
  "CDCP1",
  "TGFB1"
)

TREM2_STATE <- "Lipid-associated/TREM2_Mphi"

CANONICAL_STATES <- c(
  "Resident/Kupffer-like_Mphi",
  "Repair/Resolution_Mphi",
  "Lipid-associated/TREM2_Mphi",
  "Inflammatory/Transitional_Mphi",
  "Inflammatory_Monocyte"
)

ENDPOINTS <- c(
  "Hallmark_whole_HSC_pseudobulk_score",
  "ECM_pseudobulk_score",
  "fibrogenic_pseudobulk_score",
  "NAS"
)

# ------------------------------------------------------------
# 2. Load objects/resources
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

hsc_endpoints <- read_tsv(
  HSC_ENDPOINTS,
  show_col_types = FALSE
)

primary_geneset <- read_tsv(
  NICHE_PRIMARY_GENESET,
  show_col_types = FALSE
)$gene

activity_all <- read_tsv(
  NICHE_ACTIVITY_ALL,
  show_col_types = FALSE
)

lr_network <- readRDS(LR_FILE) %>%
  distinct(from, to)

ligand_target_matrix <- readRDS(LT_FILE)

if (!all(FOCUS_LIGANDS %in% rownames(mphi))) {
  stop(
    "Missing focus ligands in myeloid object: ",
    paste(
      setdiff(FOCUS_LIGANDS, rownames(mphi)),
      collapse = ";"
    )
  )
}

if (ncol(hsc) != 7640) {
  warning("Canonical HSC cell count differs from expected 7640.")
}

# ------------------------------------------------------------
# 3. Metadata
# ------------------------------------------------------------

mmd <- mphi[[]] %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample = as.character(sample_v17),
    NAS = suppressWarnings(as.numeric(NAS_v17)),
    sender_state = as.character(myeloid_state_FINAL_v19)
  ) %>%
  filter(sender_state %in% CANONICAL_STATES)

hmd <- hsc[[]] %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample = as.character(sample),
    NAS = suppressWarnings(as.numeric(NAS_score))
  )

# ------------------------------------------------------------
# 4. Focus-ligand sender burden
# ------------------------------------------------------------

mphi_expr <- FetchData(
  mphi,
  vars = FOCUS_LIGANDS,
  layer = "data"
) %>%
  rownames_to_column("cell")

mphi_count <- FetchData(
  mphi,
  vars = FOCUS_LIGANDS,
  layer = "counts"
) %>%
  rownames_to_column("cell")

expr_long <- mphi_expr %>%
  pivot_longer(
    -cell,
    names_to = "ligand",
    values_to = "expr"
  ) %>%
  left_join(mmd, by = "cell") %>%
  filter(!is.na(sender_state))

count_long <- mphi_count %>%
  pivot_longer(
    -cell,
    names_to = "ligand",
    values_to = "count"
  ) %>%
  left_join(mmd, by = "cell") %>%
  filter(!is.na(sender_state))

total_n <- mmd %>%
  count(sample, NAS, name = "total_myeloid_n")

state_n <- mmd %>%
  count(sample, NAS, sender_state, name = "state_n")

expr_sum <- expr_long %>%
  group_by(sample, NAS, sender_state, ligand) %>%
  summarise(
    ligand_mean = mean(expr, na.rm = TRUE),
    ligand_expr_sum = sum(expr, na.rm = TRUE),
    .groups = "drop"
  )

det_sum <- count_long %>%
  group_by(sample, NAS, sender_state, ligand) %>%
  summarise(
    ligand_detection = mean(count > 0, na.rm = TRUE),
    ligand_detected_n = sum(count > 0, na.rm = TRUE),
    .groups = "drop"
  )

burden <- expr_sum %>%
  full_join(
    det_sum,
    by = c("sample","NAS","sender_state","ligand")
  ) %>%
  left_join(
    state_n,
    by = c("sample","NAS","sender_state")
  ) %>%
  left_join(
    total_n,
    by = c("sample","NAS")
  ) %>%
  group_by(sample, NAS, ligand) %>%
  mutate(
    total_ligand_expr_sum =
      sum(ligand_expr_sum, na.rm = TRUE),
    total_ligand_detected_n =
      sum(ligand_detected_n, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    state_fraction = state_n / total_myeloid_n,
    ligand_burden_expr =
      ligand_expr_sum / total_myeloid_n,
    ligand_burden_detection =
      ligand_detected_n / total_myeloid_n,
    ligand_expression_share =
      ifelse(
        total_ligand_expr_sum > 0,
        ligand_expr_sum / total_ligand_expr_sum,
        0
      ),
    ligand_positive_cell_share =
      ifelse(
        total_ligand_detected_n > 0,
        ligand_detected_n / total_ligand_detected_n,
        0
      )
  ) %>%
  arrange(ligand, sample, sender_state)

write_tsv(
  burden,
  file.path(
    DIR_TABLE,
    "01_focus_ligand_burden_by_sample_state.tsv"
  )
)

# ------------------------------------------------------------
# 5. Statistical helpers
# ------------------------------------------------------------

safe_spearman <- function(x, y) {

  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]

  if (
    length(x) < 4 ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(
      tibble(
        n = length(x),
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  ct <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )

  tibble(
    n = length(x),
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

all_permutations <- function(x) {
  if (length(x) == 1) {
    return(matrix(x, nrow = 1))
  }

  out <- vector("list", length(x))

  for (i in seq_along(x)) {
    out[[i]] <- cbind(
      x[i],
      all_permutations(x[-i])
    )
  }

  do.call(rbind, out)
}

exact_perm <- function(x, y) {

  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]

  n <- length(x)

  if (n > 8) {
    stop("Exact permutation function is restricted to n <= 8.")
  }

  rho_obs <- suppressWarnings(
    cor(x, y, method = "spearman")
  )

  perms <- all_permutations(seq_len(n))

  rhos <- apply(
    perms,
    1,
    function(idx) {
      suppressWarnings(
        cor(
          x,
          y[idx],
          method = "spearman"
        )
      )
    }
  )

  tibble(
    n = n,
    rho = rho_obs,
    permutations = length(rhos),
    p_exact =
      mean(
        abs(rhos) >= abs(rho_obs) - 1e-12
      )
  )
}

partial_spearman <- function(x, y, z) {

  keep <- complete.cases(x, y, z)
  x <- x[keep]
  y <- y[keep]
  z <- z[keep]

  n <- length(x)

  if (n < 5) {
    return(
      tibble(
        n = n,
        partial_rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  r_xy <- suppressWarnings(cor(x, y, method = "spearman"))
  r_xz <- suppressWarnings(cor(x, z, method = "spearman"))
  r_yz <- suppressWarnings(cor(y, z, method = "spearman"))

  denom <- sqrt(
    (1 - r_xz^2) *
      (1 - r_yz^2)
  )

  if (!is.finite(denom) || denom <= 0) {
    return(
      tibble(
        n = n,
        partial_rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  rp <- (
    r_xy - r_xz * r_yz
  ) / denom

  rp <- max(
    min(rp, 0.999999),
    -0.999999
  )

  df <- n - 3

  tval <- rp *
    sqrt(
      df /
        (1 - rp^2)
    )

  tibble(
    n = n,
    partial_rho = rp,
    p_value =
      2 * pt(
        -abs(tval),
        df = df
      )
  )
}

# ------------------------------------------------------------
# 6. Primary TREM2 associations with exact permutation
# ------------------------------------------------------------

trem2 <- burden %>%
  filter(sender_state == TREM2_STATE) %>%
  inner_join(
    hsc_endpoints,
    by = c("sample","NAS")
  )

primary <- bind_rows(
  lapply(
    FOCUS_LIGANDS,
    function(lig) {

      d <- trem2 %>%
        filter(ligand == lig) %>%
        arrange(sample)

      bind_rows(
        lapply(
          ENDPOINTS,
          function(ep) {

            asym <- safe_spearman(
              d$ligand_burden_expr,
              d[[ep]]
            )

            exact <- exact_perm(
              d$ligand_burden_expr,
              d[[ep]]
            )

            if (ep != "NAS") {
              part <- partial_spearman(
                d$ligand_burden_expr,
                d[[ep]],
                d$NAS
              )
            } else {
              part <- tibble(
                n = nrow(d),
                partial_rho = NA_real_,
                p_value = NA_real_
              )
            }

            tibble(
              ligand = lig,
              receiver_endpoint = ep,
              n = asym$n,
              rho = asym$rho,
              p_value = asym$p_value,
              p_exact = exact$p_exact,
              partial_rho_NAS = part$partial_rho,
              partial_p_NAS = part$p_value
            )
          }
        )
      )
    }
  )
) %>%
  group_by(receiver_endpoint) %>%
  mutate(
    FDR_exact_within_endpoint =
      p.adjust(
        p_exact,
        method = "BH"
      ),
    FDR_partial_within_endpoint =
      p.adjust(
        partial_p_NAS,
        method = "BH"
      )
  ) %>%
  ungroup()

write_tsv(
  primary,
  file.path(
    DIR_TABLE,
    "02_PRIMARY_focus_ligand_vs_HSC_endpoints.tsv"
  )
)

# ------------------------------------------------------------
# 7. Burden decomposition vs primary Hallmark endpoint
# ------------------------------------------------------------

BURDEN_METRICS <- c(
  "ligand_expression_share",
  "ligand_burden_expr",
  "ligand_burden_detection",
  "ligand_mean",
  "ligand_detection",
  "state_fraction"
)

decomposition <- bind_rows(
  lapply(
    FOCUS_LIGANDS,
    function(lig) {

      d <- trem2 %>%
        filter(ligand == lig)

      bind_rows(
        lapply(
          BURDEN_METRICS,
          function(metric) {

            a <- safe_spearman(
              d[[metric]],
              d$Hallmark_whole_HSC_pseudobulk_score
            )

            e <- exact_perm(
              d[[metric]],
              d$Hallmark_whole_HSC_pseudobulk_score
            )

            tibble(
              ligand = lig,
              sender_metric = metric,
              n = a$n,
              rho = a$rho,
              p_value = a$p_value,
              p_exact = e$p_exact
            )
          }
        )
      )
    }
  )
) %>%
  group_by(ligand) %>%
  mutate(
    FDR_exact_within_ligand =
      p.adjust(
        p_exact,
        method = "BH"
      )
  ) %>%
  ungroup() %>%
  arrange(ligand, FDR_exact_within_ligand, desc(rho))

write_tsv(
  decomposition,
  file.path(
    DIR_TABLE,
    "03_focus_ligand_burden_decomposition.tsv"
  )
)

# ------------------------------------------------------------
# 8. Leave-one-out for all three focus ligands
# ------------------------------------------------------------

loo <- bind_rows(
  lapply(
    FOCUS_LIGANDS,
    function(lig) {

      d <- trem2 %>%
        filter(ligand == lig) %>%
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
                omitted_sample = d$sample[i],
                omitted_NAS = d$NAS[i],
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
    "04_focus_ligand_leave_one_out.tsv"
  )
)

loo_summary <- loo %>%
  group_by(ligand) %>%
  summarise(
    tests = n(),
    all_positive = all(rho > 0, na.rm = TRUE),
    min_rho = min(rho, na.rm = TRUE),
    median_rho = median(rho, na.rm = TRUE),
    max_rho = max(rho, na.rm = TRUE),
    p_lt_0_05_n = sum(p_value < 0.05, na.rm = TRUE),
    p_lt_0_10_n = sum(p_value < 0.10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_rho))

write_tsv(
  loo_summary,
  file.path(
    DIR_TABLE,
    "05_focus_ligand_leave_one_out_summary.tsv"
  )
)

# ------------------------------------------------------------
# 9. Sender-state specificity
# ------------------------------------------------------------

state_assoc <- bind_rows(
  lapply(
    FOCUS_LIGANDS,
    function(lig) {

      d <- burden %>%
        filter(ligand == lig) %>%
        inner_join(
          hsc_endpoints,
          by = c("sample","NAS")
        )

      bind_rows(
        lapply(
          CANONICAL_STATES,
          function(st) {

            ds <- d %>%
              filter(sender_state == st) %>%
              arrange(sample)

            a <- safe_spearman(
              ds$ligand_burden_expr,
              ds$Hallmark_whole_HSC_pseudobulk_score
            )

            e <- exact_perm(
              ds$ligand_burden_expr,
              ds$Hallmark_whole_HSC_pseudobulk_score
            )

            tibble(
              ligand = lig,
              sender_state = st,
              n = a$n,
              rho = a$rho,
              p_value = a$p_value,
              p_exact = e$p_exact
            )
          }
        )
      )
    }
  )
) %>%
  group_by(ligand) %>%
  mutate(
    FDR_exact_across_states =
      p.adjust(
        p_exact,
        method = "BH"
      ),
    rho_rank =
      min_rank(
        desc(rho)
      )
  ) %>%
  ungroup() %>%
  arrange(ligand, rho_rank)

write_tsv(
  state_assoc,
  file.path(
    DIR_TABLE,
    "06_focus_ligand_sender_state_specificity.tsv"
  )
)

# ------------------------------------------------------------
# 10. HSC receptor expression support
# ------------------------------------------------------------

focus_lr <- lr_network %>%
  filter(from %in% FOCUS_LIGANDS) %>%
  distinct(from, to)

focus_receptors <- intersect(
  unique(focus_lr$to),
  rownames(hsc)
)

hsc_receptor_counts <- FetchData(
  hsc,
  vars = focus_receptors,
  layer = "counts"
) %>%
  rownames_to_column("cell") %>%
  pivot_longer(
    -cell,
    names_to = "receptor",
    values_to = "count"
  ) %>%
  left_join(
    hmd,
    by = "cell"
  )

receptor_support <- hsc_receptor_counts %>%
  group_by(
    sample,
    receptor
  ) %>%
  summarise(
    HSC_cells = n(),
    detection_fraction =
      mean(count > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(receptor) %>%
  summarise(
    samples_present = n(),
    samples_detection_ge_0.05 =
      sum(detection_fraction >= 0.05),
    mean_detection_fraction =
      mean(detection_fraction, na.rm = TRUE),
    max_detection_fraction =
      max(detection_fraction, na.rm = TRUE),
    .groups = "drop"
  )

focus_lr_support <- focus_lr %>%
  left_join(
    receptor_support,
    by = c("to" = "receptor")
  ) %>%
  mutate(
    receptor_measurable = to %in% rownames(hsc),
    receptor_primary_expressed =
      !is.na(samples_detection_ge_0.05) &
      samples_detection_ge_0.05 >= 4
  ) %>%
  arrange(
    from,
    desc(receptor_primary_expressed),
    desc(mean_detection_fraction)
  )

write_tsv(
  focus_lr_support,
  file.path(
    DIR_TABLE,
    "07_focus_ligand_HSC_receptor_support.tsv"
  )
)

# ------------------------------------------------------------
# 11. Whole-HSC pseudobulk expression matrix
# ------------------------------------------------------------

hsc_counts <- GetAssayData(
  hsc,
  assay = "RNA",
  layer = "counts"
)

samples <- sort(
  unique(hmd$sample)
)

pb_list <- lapply(
  samples,
  function(s) {

    idx <- which(
      hmd$sample == s
    )

    Matrix::rowSums(
      hsc_counts[, idx, drop = FALSE]
    )
  }
)

pb_counts <- do.call(cbind, pb_list)
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

sample_meta <- hmd %>%
  distinct(sample, NAS) %>%
  arrange(sample)

# ------------------------------------------------------------
# 12. NicheNet targets + actual sample-level target association
# ------------------------------------------------------------

target_results <- bind_rows(
  lapply(
    FOCUS_LIGANDS,
    function(lig) {

      links <- nichenetr::get_weighted_ligand_target_links(
        ligand = lig,
        geneset = primary_geneset,
        ligand_target_matrix = ligand_target_matrix,
        n = 250
      ) %>%
        arrange(desc(weight)) %>%
        slice_head(n = 100)

      d_sender <- trem2 %>%
        filter(ligand == lig) %>%
        select(
          sample,
          NAS,
          ligand_burden_expr
        ) %>%
        arrange(sample)

      bind_rows(
        lapply(
          seq_len(nrow(links)),
          function(i) {

            target <- links$target[i]

            if (!(target %in% rownames(pb_logcpm))) {
              return(
                tibble(
                  ligand = lig,
                  target = target,
                  weight = links$weight[i],
                  n = 0,
                  rho = NA_real_,
                  p_value = NA_real_,
                  p_exact = NA_real_,
                  partial_rho_NAS = NA_real_,
                  partial_p_NAS = NA_real_
                )
              )
            }

            dat <- d_sender %>%
              mutate(
                target_logCPM =
                  pb_logcpm[
                    target,
                    sample
                  ]
              )

            a <- safe_spearman(
              dat$ligand_burden_expr,
              dat$target_logCPM
            )

            e <- exact_perm(
              dat$ligand_burden_expr,
              dat$target_logCPM
            )

            p <- partial_spearman(
              dat$ligand_burden_expr,
              dat$target_logCPM,
              dat$NAS
            )

            tibble(
              ligand = lig,
              target = target,
              weight = links$weight[i],
              n = a$n,
              rho = a$rho,
              p_value = a$p_value,
              p_exact = e$p_exact,
              partial_rho_NAS = p$partial_rho,
              partial_p_NAS = p$p_value
            )
          }
        )
      )
    }
  )
) %>%
  group_by(ligand) %>%
  mutate(
    FDR_exact_within_ligand =
      p.adjust(
        p_exact,
        method = "BH"
      ),
    FDR_partial_within_ligand =
      p.adjust(
        partial_p_NAS,
        method = "BH"
      )
  ) %>%
  ungroup() %>%
  arrange(
    ligand,
    FDR_exact_within_ligand,
    desc(rho),
    desc(weight)
  )

write_tsv(
  target_results,
  file.path(
    DIR_TABLE,
    "08_focus_ligand_NicheNet_target_support.tsv"
  )
)

target_summary <- target_results %>%
  group_by(ligand) %>%
  summarise(
    NicheNet_targets_tested = n(),
    targets_positive_rho =
      sum(rho > 0, na.rm = TRUE),
    targets_exact_p_lt_0.05 =
      sum(p_exact < 0.05, na.rm = TRUE),
    targets_exact_FDR_lt_0.05 =
      sum(
        FDR_exact_within_ligand < 0.05,
        na.rm = TRUE
      ),
    targets_partial_p_lt_0.05 =
      sum(
        partial_p_NAS < 0.05,
        na.rm = TRUE
      ),
    best_target =
      target[
        which.min(
          ifelse(
            is.na(p_exact),
            Inf,
            p_exact
          )
        )
      ],
    best_target_exact_p =
      min(p_exact, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  target_summary,
  file.path(
    DIR_TABLE,
    "09_focus_ligand_target_support_summary.tsv"
  )
)

# ------------------------------------------------------------
# 13. NicheNet activity summary
# ------------------------------------------------------------

activity_focus <- activity_all %>%
  filter(
    test_ligand %in% FOCUS_LIGANDS
  ) %>%
  arrange(
    geneset,
    rank
  )

write_tsv(
  activity_focus,
  file.path(
    DIR_TABLE,
    "10_focus_ligand_NicheNet_activity.tsv"
  )
)

# ------------------------------------------------------------
# 14. Final compact comparison
# ------------------------------------------------------------

primary_hallmark <- primary %>%
  filter(
    receiver_endpoint ==
      "Hallmark_whole_HSC_pseudobulk_score"
  ) %>%
  select(
    ligand,
    rho,
    p_exact,
    FDR_exact_within_endpoint,
    partial_rho_NAS,
    partial_p_NAS
  )

trem2_specificity <- state_assoc %>%
  filter(
    sender_state == TREM2_STATE
  ) %>%
  select(
    ligand,
    TREM2_state_rho_rank = rho_rank,
    TREM2_state_exact_p = p_exact,
    TREM2_state_FDR =
      FDR_exact_across_states
  )

decomp_best <- decomposition %>%
  group_by(ligand) %>%
  arrange(
    FDR_exact_within_ligand,
    desc(rho)
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(
    ligand,
    strongest_burden_component =
      sender_metric,
    strongest_component_rho =
      rho,
    strongest_component_exact_p =
      p_exact
  )

nichenet_primary <- activity_focus %>%
  filter(
    geneset ==
      "PRIMARY_rho0.60_p0.10"
  ) %>%
  select(
    ligand = test_ligand,
    NicheNet_rank = rank,
    NicheNet_aupr_corrected =
      aupr_corrected
  )

final <- primary_hallmark %>%
  left_join(
    nichenet_primary,
    by = "ligand"
  ) %>%
  left_join(
    loo_summary,
    by = "ligand"
  ) %>%
  left_join(
    trem2_specificity,
    by = "ligand"
  ) %>%
  left_join(
    decomp_best,
    by = "ligand"
  ) %>%
  left_join(
    target_summary,
    by = "ligand"
  ) %>%
  arrange(
    NicheNet_rank
  )

write_tsv(
  final,
  file.path(
    DIR_TABLE,
    "11_FINAL_focus_ligand_comparison.tsv"
  )
)

# ------------------------------------------------------------
# 15. Figure
# ------------------------------------------------------------

p <- ggplot(
  trem2 %>%
    filter(ligand %in% FOCUS_LIGANDS),
  aes(
    x = ligand_burden_expr,
    y = Hallmark_whole_HSC_pseudobulk_score,
    label = sample
  )
) +
  geom_point(size = 2.7) +
  geom_text(
    nudge_y = 0.02,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~ ligand,
    scales = "free_x"
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Focused TREM2-Mphi ligand burden vs whole-HSC response",
    x = "Abundance-weighted ligand burden",
    y = "Whole-HSC pseudobulk TGF-beta Hallmark score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_focus_ligands_vs_HSC_Hallmark.pdf"
  ),
  p,
  width = 10,
  height = 4.5
)

# ------------------------------------------------------------
# 16. Save object
# ------------------------------------------------------------

saveRDS(
  list(
    burden = burden,
    primary = primary,
    decomposition = decomposition,
    loo = loo,
    loo_summary = loo_summary,
    state_assoc = state_assoc,
    focus_lr_support = focus_lr_support,
    target_results = target_results,
    target_summary = target_summary,
    activity_focus = activity_focus,
    final = final
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_Focused_Ligand_Closure_v2.2_results.rds"
  )
)

# ------------------------------------------------------------
# 17. Provenance / integrity
# ------------------------------------------------------------

provenance <- tibble(
  item = c(
    "analysis_version",
    "focus_ligands",
    "primary_sender_state",
    "primary_sender_metric",
    "primary_receiver_endpoint",
    "association_unit",
    "exact_permutations_n8",
    "NicheNet_target_geneset",
    "source_overwritten",
    "reclustering"
  ),
  value = c(
    "v2.2.0",
    paste(FOCUS_LIGANDS, collapse = ";"),
    TREM2_STATE,
    "ligand_expr_sum_in_state / total_myeloid_n",
    "whole-HSC pseudobulk TGF-beta Hallmark",
    "biological_sample",
    "40320",
    "v2.0.1 PRIMARY receiver geneset",
    "FALSE",
    "none"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "focused_ligand_provenance_v2.2.tsv"
  )
)

integrity <- tibble(
  item = c(
    "focus_ligands_n",
    "all_focus_ligands_measurable",
    "all_8_samples_present",
    "all_5_sender_states_present",
    "HSC_cells",
    "expected_HSC_cells",
    "exact_permutation_completed",
    "receptor_audit_completed",
    "NicheNet_target_audit_completed",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    as.character(length(FOCUS_LIGANDS)),
    as.character(
      all(
        FOCUS_LIGANDS %in%
          rownames(mphi)
      )
    ),
    as.character(
      length(unique(trem2$sample)) == 8
    ),
    as.character(
      length(
        intersect(
          CANONICAL_STATES,
          unique(burden$sender_state)
        )
      ) == 5
    ),
    as.character(ncol(hsc)),
    "7640",
    as.character(
      all(
        primary$n == 8
      )
    ),
    as.character(
      nrow(focus_lr_support) > 0
    ),
    as.character(
      nrow(target_results) > 0
    ),
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v2.2.tsv"
  )
)

# ------------------------------------------------------------
# 18. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("PRIMARY FOCUS-LIGAND ASSOCIATIONS\n")
cat("============================================================\n")
print(primary)

cat("\n")
cat("============================================================\n")
cat("BURDEN DECOMPOSITION\n")
cat("============================================================\n")
print(decomposition)

cat("\n")
cat("============================================================\n")
cat("LEAVE-ONE-OUT\n")
cat("============================================================\n")
print(loo_summary)

cat("\n")
cat("============================================================\n")
cat("SENDER-STATE SPECIFICITY\n")
cat("============================================================\n")
print(state_assoc)

cat("\n")
cat("============================================================\n")
cat("HSC RECEPTOR SUPPORT\n")
cat("============================================================\n")
print(focus_lr_support)

cat("\n")
cat("============================================================\n")
cat("TARGET SUPPORT SUMMARY\n")
cat("============================================================\n")
print(target_summary)

cat("\n")
cat("============================================================\n")
cat("FINAL COMPARISON\n")
cat("============================================================\n")
print(final)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
