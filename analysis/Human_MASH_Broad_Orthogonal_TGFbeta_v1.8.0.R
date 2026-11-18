#!/usr/bin/env Rscript

# ============================================================
# Human MASH broad orthogonal TGF-beta pathway audit v1.8.0
#
# PURPOSE
#   Broaden the v1.7 orthogonal TGF-beta response validation
#   beyond the 3-gene score (SMAD7/PMEPA1/THBS1).
#
# PRIMARY gene set:
#   MSigDB HALLMARK_TGF_BETA_SIGNALING (54 genes),
#   with the following removed to reduce circularity:
#     - TGFB1   : sender ligand under test
#     - TGFBR1  : receptor used in the LR hypothesis
#     - any gene overlapping the ECM-core score used to define
#       ECM-high HSC (currently SERPINE1 is the overlap)
#
# ANALYSES
#   1. Audit Flex measurability of the Hallmark set
#   2. Construct broad orthogonal Hallmark score
#   3. Test TREM2-Mphi TGFB1 burden vs broad Hallmark score
#   4. Exact permutation test (8! = 40,320)
#   5. Partial Spearman controlling for NAS
#   6. Leave-one-sample-out robustness
#   7. Compare all five myeloid sender states
#   8. ECM-high vs ECM-nonhigh paired pathway comparison
#   9. Paired pseudobulk gene-level audit of Hallmark genes
#
# STATISTICAL PRINCIPLE
#   Biological sample is the unit of association.
#   Cell-level values are summarized within biological samples.
#   No pseudoreplication is used for inference.
#
# GENE SET SOURCE / ATTRIBUTION
#   MSigDB HALLMARK_TGF_BETA_SIGNALING, M5896
#   Broad Institute / MIT / Regents of University of California
#   CC-BY-4.0
#   https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/
#       HALLMARK_TGF_BETA_SIGNALING.html
#
# IMPORTANT
#   - No reclustering
#   - Source RDS is never overwritten
#   - Canonical HSC receiver remains 7,640 cells
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(Matrix)
})

set.seed(20261118)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

HSC_RDS <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_HSC_ECM_Proliferation_v1.2",
  "objects",
  "Human_MASH_HSC_ECM_Proliferation_v1.2.rds"
)

SENDER_TABLE <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_TGFB1_Sender_Burden_Audit_v1.5",
  "tables",
  "02_TGFB1_sender_state_burden_by_sample.tsv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_Broad_Orthogonal_TGFbeta_v1.8"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG   <- file.path(OUT, "figures")
DIR_OBJ   <- file.path(OUT, "objects")
DIR_LOG   <- file.path(OUT, "logs")
DIR_PROV  <- file.path(OUT, "provenance")

dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OBJ,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOG,   recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PROV,  recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(
  DIR_LOG,
  "Human_MASH_Broad_Orthogonal_TGFbeta_v1.8.0.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH broad orthogonal TGF-beta audit v1.8.0\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")

stopifnot(file.exists(HSC_RDS))
stopifnot(file.exists(SENDER_TABLE))

# ------------------------------------------------------------
# 1. Load
# ------------------------------------------------------------

hsc <- readRDS(HSC_RDS)
DefaultAssay(hsc) <- "RNA"

sender <- read_tsv(
  SENDER_TABLE,
  show_col_types = FALSE
)

cat("HSC cells:", ncol(hsc), "\n")
cat("HSC features:", nrow(hsc), "\n\n")

if (ncol(hsc) != 7640) {
  warning(
    "Expected 7640 canonical HSC cells; observed ",
    ncol(hsc)
  )
}

required_hsc <- c(
  "sample",
  "NAS_score",
  "HSC_state_v1.2",
  "ECM_high_primary_v1.2",
  "ECM_total_score_v1.2"
)

required_sender <- c(
  "sample",
  "NAS",
  "sender_state",
  "TGFB1_burden_expr",
  "TGFB1_expression_share"
)

missing_hsc <- setdiff(
  required_hsc,
  colnames(hsc[[]])
)

missing_sender <- setdiff(
  required_sender,
  colnames(sender)
)

if (length(missing_hsc) > 0) {
  stop(
    "Missing HSC metadata: ",
    paste(missing_hsc, collapse = ";")
  )
}

if (length(missing_sender) > 0) {
  stop(
    "Missing sender columns: ",
    paste(missing_sender, collapse = ";")
  )
}

# ------------------------------------------------------------
# 2. Normalize RNA
# ------------------------------------------------------------

hsc <- NormalizeData(
  hsc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ------------------------------------------------------------
# 3. MSigDB HALLMARK_TGF_BETA_SIGNALING
#    Current human Hallmark set: 54 genes
# ------------------------------------------------------------

HALLMARK_TGF_BETA_SIGNALING <- c(
  "ACVR1",
  "APC",
  "ARID4B",
  "BCAR3",
  "BMP2",
  "BMPR1A",
  "BMPR2",
  "CDH1",
  "CDK9",
  "CDKN1C",
  "CTNNB1",
  "ENG",
  "FKBP1A",
  "FNTA",
  "FURIN",
  "HDAC1",
  "HIPK2",
  "ID1",
  "ID2",
  "ID3",
  "IFNGR2",
  "JUNB",
  "KLF10",
  "LEFTY2",
  "LTBP2",
  "MAP3K7",
  "NCOR2",
  "NOG",
  "PMEPA1",
  "PPM1A",
  "PPP1CA",
  "PPP1R15A",
  "RAB31",
  "RHOA",
  "SERPINE1",
  "SKI",
  "SKIL",
  "SLC20A1",
  "SMAD1",
  "SMAD3",
  "SMAD6",
  "SMAD7",
  "SMURF1",
  "SMURF2",
  "SPTBN1",
  "TGFB1",
  "TGFBR1",
  "TGIF1",
  "THBS1",
  "TJP1",
  "TRIM33",
  "UBE2D3",
  "WWTR1",
  "XIAP"
)

ECM_CORE <- c(
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "COL5A1",
  "COL5A2",
  "FN1",
  "SPARC",
  "THBS2",
  "LOXL2",
  "TIMP1",
  "TGFBI",
  "CTHRC1",
  "SERPINE1",
  "POSTN"
)

CIRCULARITY_EXCLUSIONS <- unique(
  c(
    "TGFB1",
    "TGFBR1",
    intersect(
      HALLMARK_TGF_BETA_SIGNALING,
      ECM_CORE
    )
  )
)

HALLMARK_ORTHOGONAL <- setdiff(
  HALLMARK_TGF_BETA_SIGNALING,
  CIRCULARITY_EXCLUSIONS
)

available_full <- intersect(
  HALLMARK_TGF_BETA_SIGNALING,
  rownames(hsc)
)

available_orthogonal <- intersect(
  HALLMARK_ORTHOGONAL,
  rownames(hsc)
)

missing_orthogonal <- setdiff(
  HALLMARK_ORTHOGONAL,
  rownames(hsc)
)

panel_audit <- tibble(
  gene = HALLMARK_TGF_BETA_SIGNALING,
  in_Hallmark = TRUE,
  overlaps_ECM_core =
    gene %in% ECM_CORE,
  excluded_for_circularity =
    gene %in% CIRCULARITY_EXCLUSIONS,
  available_in_HSC =
    gene %in% rownames(hsc),
  used_in_primary_orthogonal =
    gene %in% available_orthogonal
)

write_tsv(
  panel_audit,
  file.path(
    DIR_TABLE,
    "01_Hallmark_TGFbeta_feature_audit.tsv"
  )
)

cat("Hallmark genes requested:", length(HALLMARK_TGF_BETA_SIGNALING), "\n")
cat("Hallmark genes available:", length(available_full), "\n")
cat("Primary orthogonal genes available:", length(available_orthogonal), "\n")
cat("Circularity exclusions:", paste(CIRCULARITY_EXCLUSIONS, collapse = ";"), "\n\n")

if (length(available_orthogonal) < 10) {
  warning(
    "Fewer than 10 orthogonal Hallmark genes are measurable; ",
    "interpret the broad score cautiously."
  )
}

if (length(available_orthogonal) < 5) {
  stop(
    "Too few orthogonal Hallmark genes are measurable."
  )
}

# ------------------------------------------------------------
# 4. Scoring helper
# ------------------------------------------------------------

z_fun <- function(x) {

  sx <- sd(
    x,
    na.rm = TRUE
  )

  if (
    is.na(sx) ||
    sx == 0
  ) {
    return(
      rep(0, length(x))
    )
  }

  (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) / sx
}

calc_equal_weight_zscore <- function(
  object,
  genes
) {

  x <- FetchData(
    object,
    vars = genes,
    layer = "data"
  )

  xz <- as.data.frame(
    lapply(
      x,
      z_fun
    )
  )

  rowMeans(
    xz,
    na.rm = TRUE
  )
}

# ------------------------------------------------------------
# 5. Broad orthogonal Hallmark score
# ------------------------------------------------------------

hsc$Hallmark_TGFbeta_orthogonal_score_v1.8 <-
  calc_equal_weight_zscore(
    hsc,
    available_orthogonal
  )

hsc$Hallmark_TGFbeta_full_score_v1.8 <-
  calc_equal_weight_zscore(
    hsc,
    available_full
  )

hmd <- hsc[[]]

cell_df <- hmd %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample = as.character(sample),
    NAS = suppressWarnings(
      as.numeric(NAS_score)
    ),
    HSC_state = as.character(
      HSC_state_v1.2
    ),
    ECM_high = as.logical(
      ECM_high_primary_v1.2
    ),
    ECM_total_score = as.numeric(
      ECM_total_score_v1.2
    ),
    Hallmark_orthogonal_score =
      as.numeric(
        Hallmark_TGFbeta_orthogonal_score_v1.8
      ),
    Hallmark_full_score =
      as.numeric(
        Hallmark_TGFbeta_full_score_v1.8
      )
  )

# ------------------------------------------------------------
# 6. Sample-level HSC summary
# ------------------------------------------------------------

sample_hsc <- cell_df %>%
  group_by(
    sample,
    NAS
  ) %>%
  summarise(
    HSC_n = n(),

    ECM_high_fraction =
      mean(
        ECM_high,
        na.rm = TRUE
      ),

    ECM_total_mean =
      mean(
        ECM_total_score,
        na.rm = TRUE
      ),

    Hallmark_orthogonal_mean =
      mean(
        Hallmark_orthogonal_score,
        na.rm = TRUE
      ),

    Hallmark_full_mean =
      mean(
        Hallmark_full_score,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  arrange(sample)

write_tsv(
  sample_hsc,
  file.path(
    DIR_TABLE,
    "02_Hallmark_TGFbeta_by_sample.tsv"
  )
)

# ------------------------------------------------------------
# 7. ECM-high vs ECM-nonhigh sample summaries
# ------------------------------------------------------------

ecm_group <- cell_df %>%
  mutate(
    ECM_group = ifelse(
      ECM_high,
      "ECM_high",
      "ECM_nonhigh"
    )
  ) %>%
  group_by(
    sample,
    NAS,
    ECM_group
  ) %>%
  summarise(
    cells = n(),

    Hallmark_orthogonal_mean =
      mean(
        Hallmark_orthogonal_score,
        na.rm = TRUE
      ),

    Hallmark_full_mean =
      mean(
        Hallmark_full_score,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  ecm_group,
  file.path(
    DIR_TABLE,
    "03_Hallmark_TGFbeta_by_ECM_group.tsv"
  )
)

ecm_wide <- ecm_group %>%
  pivot_wider(
    names_from = ECM_group,
    values_from = c(
      cells,
      Hallmark_orthogonal_mean,
      Hallmark_full_mean
    )
  )

paired_pathway <- bind_rows(
  {
    x <- ecm_wide$Hallmark_orthogonal_mean_ECM_high
    y <- ecm_wide$Hallmark_orthogonal_mean_ECM_nonhigh

    wt <- suppressWarnings(
      wilcox.test(
        x,
        y,
        paired = TRUE,
        exact = FALSE
      )
    )

    tibble(
      metric = "Hallmark_orthogonal_mean",
      n_samples = sum(
        complete.cases(x, y)
      ),
      median_delta =
        median(
          x - y,
          na.rm = TRUE
        ),
      positive_delta_samples =
        sum(
          x > y,
          na.rm = TRUE
        ),
      p_value =
        wt$p.value
    )
  },

  {
    x <- ecm_wide$Hallmark_full_mean_ECM_high
    y <- ecm_wide$Hallmark_full_mean_ECM_nonhigh

    wt <- suppressWarnings(
      wilcox.test(
        x,
        y,
        paired = TRUE,
        exact = FALSE
      )
    )

    tibble(
      metric = "Hallmark_full_mean",
      n_samples = sum(
        complete.cases(x, y)
      ),
      median_delta =
        median(
          x - y,
          na.rm = TRUE
        ),
      positive_delta_samples =
        sum(
          x > y,
          na.rm = TRUE
        ),
      p_value =
        wt$p.value
    )
  }
) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  )

write_tsv(
  paired_pathway,
  file.path(
    DIR_TABLE,
    "04_ECMhigh_vs_nonhigh_Hallmark_paired_test.tsv"
  )
)

# ------------------------------------------------------------
# 8. Statistical helpers
# ------------------------------------------------------------

safe_spearman <- function(
  x,
  y
) {

  keep <- complete.cases(
    x,
    y
  )

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
    rho = unname(
      ct$estimate
    ),
    p_value = ct$p.value
  )
}

all_permutations <- function(x) {

  if (length(x) == 1) {
    return(
      matrix(
        x,
        nrow = 1
      )
    )
  }

  out <- vector(
    "list",
    length(x)
  )

  for (i in seq_along(x)) {

    p_rest <- all_permutations(
      x[-i]
    )

    out[[i]] <- cbind(
      x[i],
      p_rest
    )
  }

  do.call(
    rbind,
    out
  )
}

exact_perm <- function(
  x,
  y
) {

  keep <- complete.cases(
    x,
    y
  )

  x <- x[keep]
  y <- y[keep]

  n <- length(x)

  if (n > 8) {
    stop(
      "Exact permutation implemented only for n <= 8."
    )
  }

  rho_obs <- suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )

  perms <- all_permutations(
    seq_len(n)
  )

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
        abs(rhos) >=
          abs(rho_obs) -
          1e-12
      )
  )
}

partial_spearman <- function(
  x,
  y,
  z
) {

  keep <- complete.cases(
    x,
    y,
    z
  )

  x <- x[keep]
  y <- y[keep]
  z <- z[keep]

  n <- length(x)

  if (
    n < 5 ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2 ||
    length(unique(z)) < 2
  ) {
    return(
      tibble(
        n = n,
        partial_rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  r_xy <- suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )

  r_xz <- suppressWarnings(
    cor(
      x,
      z,
      method = "spearman"
    )
  )

  r_yz <- suppressWarnings(
    cor(
      y,
      z,
      method = "spearman"
    )
  )

  denom <- sqrt(
    (1 - r_xz^2) *
      (1 - r_yz^2)
  )

  if (
    !is.finite(denom) ||
    denom <= 0
  ) {
    return(
      tibble(
        n = n,
        partial_rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  rp <- (
    r_xy -
      r_xz *
      r_yz
  ) / denom

  rp <- max(
    min(
      rp,
      0.999999
    ),
    -0.999999
  )

  df <- n - 3

  t_value <- rp *
    sqrt(
      df /
        (1 - rp^2)
    )

  p_value <- 2 *
    pt(
      -abs(t_value),
      df = df
    )

  tibble(
    n = n,
    partial_rho = rp,
    p_value = p_value
  )
}

# ------------------------------------------------------------
# 9. TREM2-focused association
# ------------------------------------------------------------

TREM2_STATE <- "Lipid-associated/TREM2_Mphi"

trem2 <- sender %>%
  filter(
    sender_state == TREM2_STATE
  ) %>%
  select(
    sample,
    NAS,
    TGFB1_burden_expr,
    TGFB1_expression_share
  ) %>%
  inner_join(
    sample_hsc,
    by = c(
      "sample",
      "NAS"
    )
  ) %>%
  arrange(sample)

if (nrow(trem2) != 8) {
  stop(
    "Expected 8 TREM2 sample rows; observed ",
    nrow(trem2)
  )
}

primary_spearman <- safe_spearman(
  trem2$TGFB1_burden_expr,
  trem2$Hallmark_orthogonal_mean
) %>%
  mutate(
    sender_metric =
      "TREM2_TGFB1_burden_expr",
    receiver_metric =
      "Hallmark_TGFbeta_orthogonal_mean",
    .before = 1
  )

primary_exact <- exact_perm(
  trem2$TGFB1_burden_expr,
  trem2$Hallmark_orthogonal_mean
) %>%
  mutate(
    sender_metric =
      "TREM2_TGFB1_burden_expr",
    receiver_metric =
      "Hallmark_TGFbeta_orthogonal_mean",
    .before = 1
  )

primary_partial <- partial_spearman(
  trem2$TGFB1_burden_expr,
  trem2$Hallmark_orthogonal_mean,
  trem2$NAS
) %>%
  mutate(
    sender_metric =
      "TREM2_TGFB1_burden_expr",
    receiver_metric =
      "Hallmark_TGFbeta_orthogonal_mean",
    covariate = "NAS",
    .before = 1
  )

write_tsv(
  primary_spearman,
  file.path(
    DIR_TABLE,
    "05_PRIMARY_TREM2_TGFB1_vs_Hallmark_Spearman.tsv"
  )
)

write_tsv(
  primary_exact,
  file.path(
    DIR_TABLE,
    "06_PRIMARY_TREM2_TGFB1_vs_Hallmark_exact_permutation.tsv"
  )
)

write_tsv(
  primary_partial,
  file.path(
    DIR_TABLE,
    "07_PRIMARY_TREM2_TGFB1_vs_Hallmark_partial_NAS.tsv"
  )
)

# ------------------------------------------------------------
# 10. Leave-one-out
# ------------------------------------------------------------

loo <- bind_rows(
  lapply(
    seq_len(nrow(trem2)),
    function(i) {

      sub <- trem2[-i, ]

      safe_spearman(
        sub$TGFB1_burden_expr,
        sub$Hallmark_orthogonal_mean
      ) %>%
        mutate(
          omitted_sample =
            trem2$sample[i],
          omitted_NAS =
            trem2$NAS[i],
          .before = 1
        )
    }
  )
)

write_tsv(
  loo,
  file.path(
    DIR_TABLE,
    "08_leave_one_sample_out_Hallmark.tsv"
  )
)

loo_summary <- tibble(
  tests = nrow(loo),

  all_positive =
    all(
      loo$rho > 0,
      na.rm = TRUE
    ),

  min_rho =
    min(
      loo$rho,
      na.rm = TRUE
    ),

  median_rho =
    median(
      loo$rho,
      na.rm = TRUE
    ),

  max_rho =
    max(
      loo$rho,
      na.rm = TRUE
    ),

  p_lt_0_05_n =
    sum(
      loo$p_value < 0.05,
      na.rm = TRUE
    ),

  p_lt_0_10_n =
    sum(
      loo$p_value < 0.10,
      na.rm = TRUE
    )
)

write_tsv(
  loo_summary,
  file.path(
    DIR_TABLE,
    "09_leave_one_sample_out_Hallmark_summary.tsv"
  )
)

# ------------------------------------------------------------
# 11. Compare all sender states
# ------------------------------------------------------------

all_sender <- sender %>%
  inner_join(
    sample_hsc %>%
      select(
        sample,
        NAS,
        Hallmark_orthogonal_mean
      ),
    by = c(
      "sample",
      "NAS"
    )
  ) %>%
  group_by(
    sender_state
  ) %>%
  group_modify(
    ~ safe_spearman(
      .x$TGFB1_burden_expr,
      .x$Hallmark_orthogonal_mean
    )
  ) %>%
  ungroup() %>%
  mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  arrange(
    desc(rho)
  )

write_tsv(
  all_sender,
  file.path(
    DIR_TABLE,
    "10_all_sender_states_vs_Hallmark_TGFbeta.tsv"
  )
)

# ------------------------------------------------------------
# 12. Pseudobulk counts by sample x ECM group
# ------------------------------------------------------------

counts <- GetAssayData(
  hsc,
  assay = "RNA",
  layer = "counts"
)

meta <- hsc[[]]

group_id <- paste(
  as.character(
    meta$sample
  ),
  ifelse(
    as.logical(
      meta$ECM_high_primary_v1.2
    ),
    "ECM_high",
    "ECM_nonhigh"
  ),
  sep = "__"
)

groups <- unique(group_id)

pb_list <- lapply(
  groups,
  function(g) {

    cells <- which(
      group_id == g
    )

    Matrix::rowSums(
      counts[, cells, drop = FALSE]
    )
  }
)

pb_counts <- do.call(
  cbind,
  pb_list
)

rownames(pb_counts) <- rownames(counts)
colnames(pb_counts) <- groups

pb_libsize <- colSums(
  pb_counts
)

pb_cpm <- sweep(
  pb_counts,
  2,
  pb_libsize,
  "/"
) * 1e6

pb_logcpm <- log2(
  pb_cpm + 1
)

pb_meta <- tibble(
  pseudobulk = groups
) %>%
  separate(
    pseudobulk,
    into = c(
      "sample",
      "ECM_group"
    ),
    sep = "__",
    remove = FALSE
  ) %>%
  left_join(
    sample_hsc %>%
      select(
        sample,
        NAS
      ),
    by = "sample"
  ) %>%
  mutate(
    library_size =
      pb_libsize[pseudobulk]
  )

write_tsv(
  pb_meta,
  file.path(
    DIR_TABLE,
    "11_pseudobulk_metadata.tsv"
  )
)

# ------------------------------------------------------------
# 13. Paired pseudobulk Hallmark gene effects
# ------------------------------------------------------------

pb_gene_long <- as_tibble(
  t(
    pb_logcpm[
      available_orthogonal,
      ,
      drop = FALSE
    ]
  ),
  rownames = "pseudobulk"
) %>%
  left_join(
    pb_meta,
    by = "pseudobulk"
  ) %>%
  pivot_longer(
    cols = all_of(
      available_orthogonal
    ),
    names_to = "gene",
    values_to = "logCPM"
  )

pb_gene_wide <- pb_gene_long %>%
  select(
    sample,
    NAS,
    ECM_group,
    gene,
    logCPM
  ) %>%
  pivot_wider(
    names_from = ECM_group,
    values_from = logCPM
  ) %>%
  mutate(
    delta_ECMhigh_minus_nonhigh =
      ECM_high -
      ECM_nonhigh
  )

write_tsv(
  pb_gene_wide,
  file.path(
    DIR_TABLE,
    "12_pseudobulk_Hallmark_gene_deltas_by_sample.tsv"
  )
)

pb_gene_tests <- pb_gene_wide %>%
  group_by(gene) %>%
  summarise(
    n_samples =
      sum(
        complete.cases(
          ECM_high,
          ECM_nonhigh
        )
      ),

    mean_delta =
      mean(
        delta_ECMhigh_minus_nonhigh,
        na.rm = TRUE
      ),

    median_delta =
      median(
        delta_ECMhigh_minus_nonhigh,
        na.rm = TRUE
      ),

    positive_delta_samples =
      sum(
        delta_ECMhigh_minus_nonhigh > 0,
        na.rm = TRUE
      ),

    p_value = {
      d <- delta_ECMhigh_minus_nonhigh[
        is.finite(
          delta_ECMhigh_minus_nonhigh
        )
      ]

      if (
        length(d) < 4 ||
        all(d == 0)
      ) {
        NA_real_
      } else {
        suppressWarnings(
          wilcox.test(
            d,
            mu = 0,
            exact = FALSE
          )$p.value
        )
      }
    },

    .groups = "drop"
  ) %>%
  mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  arrange(
    FDR,
    desc(mean_delta)
  )

write_tsv(
  pb_gene_tests,
  file.path(
    DIR_TABLE,
    "13_pseudobulk_Hallmark_gene_paired_tests.tsv"
  )
)

# ------------------------------------------------------------
# 14. Pseudobulk pathway-level score
#
# Gene-wise z-score across the 16 pseudobulk profiles,
# then equal-weight mean across orthogonal Hallmark genes.
# ------------------------------------------------------------

pb_hallmark_matrix <- t(
  pb_logcpm[
    available_orthogonal,
    ,
    drop = FALSE
  ]
)

pb_hallmark_z <- apply(
  pb_hallmark_matrix,
  2,
  z_fun
)

if (is.vector(pb_hallmark_z)) {
  pb_hallmark_z <- matrix(
    pb_hallmark_z,
    ncol = 1
  )
}

pb_pathway_score <- rowMeans(
  pb_hallmark_z,
  na.rm = TRUE
)

pb_pathway <- pb_meta %>%
  mutate(
    Hallmark_orthogonal_pb_score =
      pb_pathway_score[
        pseudobulk
      ]
  )

write_tsv(
  pb_pathway,
  file.path(
    DIR_TABLE,
    "14_pseudobulk_Hallmark_pathway_score.tsv"
  )
)

pb_pathway_wide <- pb_pathway %>%
  select(
    sample,
    NAS,
    ECM_group,
    Hallmark_orthogonal_pb_score
  ) %>%
  pivot_wider(
    names_from = ECM_group,
    values_from =
      Hallmark_orthogonal_pb_score
  ) %>%
  mutate(
    delta_ECMhigh_minus_nonhigh =
      ECM_high -
      ECM_nonhigh
  )

pb_pathway_test <- {
  d <- pb_pathway_wide$
    delta_ECMhigh_minus_nonhigh

  wt <- suppressWarnings(
    wilcox.test(
      d,
      mu = 0,
      exact = FALSE
    )
  )

  tibble(
    n_samples =
      sum(
        is.finite(d)
      ),
    median_delta =
      median(
        d,
        na.rm = TRUE
      ),
    positive_delta_samples =
      sum(
        d > 0,
        na.rm = TRUE
      ),
    p_value =
      wt$p.value
  )
}

write_tsv(
  pb_pathway_test,
  file.path(
    DIR_TABLE,
    "15_PRIMARY_pseudobulk_Hallmark_pathway_paired_test.tsv"
  )
)

# ------------------------------------------------------------
# 15. NAS summary
# ------------------------------------------------------------

NAS_summary <- trem2 %>%
  group_by(NAS) %>%
  summarise(
    samples = n(),

    TREM2_TGFB1_burden =
      mean(
        TGFB1_burden_expr,
        na.rm = TRUE
      ),

    TREM2_TGFB1_share =
      mean(
        TGFB1_expression_share,
        na.rm = TRUE
      ),

    Hallmark_orthogonal_mean =
      mean(
        Hallmark_orthogonal_mean,
        na.rm = TRUE
      ),

    ECM_high_fraction =
      mean(
        ECM_high_fraction,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  NAS_summary,
  file.path(
    DIR_TABLE,
    "16_NAS_mean_Hallmark_TGFbeta_summary.tsv"
  )
)

# ------------------------------------------------------------
# 16. Figures
# ------------------------------------------------------------

p1 <- ggplot(
  trem2,
  aes(
    x = TGFB1_burden_expr,
    y = Hallmark_orthogonal_mean,
    label = sample
  )
) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.01,
    check_overlap = TRUE
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "TREM2-Mphi TGFB1 burden vs broad orthogonal TGF-beta Hallmark",
    x = "TREM2-Mphi TGFB1 burden",
    y = "HSC orthogonal Hallmark TGF-beta score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_TREM2_TGFB1_vs_Hallmark_TGFbeta.pdf"
  ),
  p1,
  width = 7,
  height = 5.5
)

p2 <- ggplot(
  ecm_group,
  aes(
    x = ECM_group,
    y = Hallmark_orthogonal_mean,
    group = sample
  )
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "Broad orthogonal TGF-beta Hallmark in ECM-high HSC",
    x = NULL,
    y = "Hallmark TGF-beta score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_Hallmark_TGFbeta_ECMhigh_vs_nonhigh.pdf"
  ),
  p2,
  width = 6.5,
  height = 5.5
)

p3 <- ggplot(
  all_sender,
  aes(
    x = reorder(
      sender_state,
      rho
    ),
    y = rho
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(
    title = "TGFB1 burden vs broad HSC TGF-beta Hallmark by sender state",
    x = NULL,
    y = "Spearman rho"
  )

ggsave(
  file.path(
    DIR_FIG,
    "03_sender_state_comparison_Hallmark_TGFbeta.pdf"
  ),
  p3,
  width = 8,
  height = 5.5
)

# ------------------------------------------------------------
# 17. Save results
# ------------------------------------------------------------

saveRDS(
  list(
    panel_audit = panel_audit,
    sample_hsc = sample_hsc,
    ecm_group = ecm_group,
    paired_pathway = paired_pathway,
    trem2 = trem2,
    primary_spearman = primary_spearman,
    primary_exact = primary_exact,
    primary_partial = primary_partial,
    loo = loo,
    all_sender = all_sender,
    pb_meta = pb_meta,
    pb_gene_tests = pb_gene_tests,
    pb_pathway = pb_pathway,
    pb_pathway_test = pb_pathway_test,
    NAS_summary = NAS_summary
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_Broad_Orthogonal_TGFbeta_v1.8_results.rds"
  )
)

# ------------------------------------------------------------
# 18. Provenance
# ------------------------------------------------------------

provenance <- tibble(
  item = c(
    "analysis_version",
    "HSC_source_RDS",
    "sender_source_table",
    "HSC_cells",
    "gene_set_source",
    "gene_set_name",
    "gene_set_original_n",
    "circularity_exclusions",
    "primary_orthogonal_available_n",
    "primary_sender_state",
    "primary_sender_metric",
    "association_unit",
    "pseudobulk_definition",
    "exact_permutations",
    "reclustering",
    "source_overwritten"
  ),
  value = c(
    "v1.8.0",
    HSC_RDS,
    SENDER_TABLE,
    as.character(
      ncol(hsc)
    ),
    "MSigDB_CC-BY-4.0",
    "HALLMARK_TGF_BETA_SIGNALING_M5896",
    as.character(
      length(
        HALLMARK_TGF_BETA_SIGNALING
      )
    ),
    paste(
      CIRCULARITY_EXCLUSIONS,
      collapse = ";"
    ),
    as.character(
      length(
        available_orthogonal
      )
    ),
    TREM2_STATE,
    "TGFB1_burden_expr",
    "biological_sample",
    "sum_RNA_counts_by_sample_x_ECM_group",
    "40320_for_n8",
    "none",
    "FALSE"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "broad_orthogonal_TGFbeta_provenance_v1.8.tsv"
  )
)

# ------------------------------------------------------------
# 19. Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "HSC_input_cells",
    "HSC_expected_cells",
    "HSC_count_matches",
    "Hallmark_original_genes_n",
    "Hallmark_available_genes_n",
    "primary_orthogonal_available_genes_n",
    "primary_orthogonal_overlap_ECM_core",
    "TGFB1_in_primary_score",
    "TGFBR1_in_primary_score",
    "all_8_samples_present",
    "exact_permutation_completed",
    "LOO_tests",
    "pseudobulk_profiles_expected",
    "pseudobulk_profiles_observed",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    as.character(
      ncol(hsc)
    ),
    "7640",
    as.character(
      ncol(hsc) == 7640
    ),
    as.character(
      length(
        HALLMARK_TGF_BETA_SIGNALING
      )
    ),
    as.character(
      length(
        available_full
      )
    ),
    as.character(
      length(
        available_orthogonal
      )
    ),
    as.character(
      length(
        intersect(
          available_orthogonal,
          ECM_CORE
        )
      ) > 0
    ),
    as.character(
      "TGFB1" %in%
        available_orthogonal
    ),
    as.character(
      "TGFBR1" %in%
        available_orthogonal
    ),
    as.character(
      nrow(trem2) == 8
    ),
    as.character(
      primary_exact$permutations ==
        40320
    ),
    as.character(
      nrow(loo)
    ),
    "16",
    as.character(
      nrow(pb_meta)
    ),
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.8.tsv"
  )
)

# ------------------------------------------------------------
# 20. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("HALLMARK PANEL AUDIT\n")
cat("============================================================\n")
print(panel_audit)

cat("\n")
cat("============================================================\n")
cat("PRIMARY TREM2 TGFB1 -> HALLMARK ASSOCIATION\n")
cat("============================================================\n")
print(primary_spearman)
print(primary_exact)
print(primary_partial)

cat("\n")
cat("============================================================\n")
cat("LEAVE-ONE-OUT\n")
cat("============================================================\n")
print(loo)
print(loo_summary)

cat("\n")
cat("============================================================\n")
cat("ECM-HIGH vs ECM-NONHIGH PATHWAY SCORE\n")
cat("============================================================\n")
print(paired_pathway)

cat("\n")
cat("============================================================\n")
cat("ALL SENDER STATES\n")
cat("============================================================\n")
print(all_sender)

cat("\n")
cat("============================================================\n")
cat("PSEUDOBULK PATHWAY TEST\n")
cat("============================================================\n")
print(pb_pathway_test)

cat("\n")
cat("============================================================\n")
cat("TOP PSEUDOBULK HALLMARK GENES\n")
cat("============================================================\n")
print(
  pb_gene_tests %>%
    head(30)
)

cat("\n")
cat("============================================================\n")
cat("NAS SUMMARY\n")
cat("============================================================\n")
print(NAS_summary)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
