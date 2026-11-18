#!/usr/bin/env Rscript

# ============================================================
# Human MASH TGF-beta cell-level vs pseudobulk concordance audit
# v1.9.0
#
# PURPOSE
#   Resolve the discordance observed in v1.8:
#
#     cell-level orthogonal Hallmark score:
#       ECM-high > ECM-nonhigh in 8/8 samples
#
#     paired pseudobulk Hallmark score:
#       median delta < 0, positive in only 2/8 samples
#
# ANALYSES
#   1. Audit RNA-depth / detected-feature differences between
#      ECM-high and ECM-nonhigh HSC.
#   2. Compare three cell-level pathway scores:
#        A. global gene-z score (v1.8-like)
#        B. within-sample gene-z score
#        C. depth/sample-residualized gene score
#   3. Recompute sample x ECM-group pseudobulk logCPM.
#   4. Quantify paired gene-wise pseudobulk deltas.
#   5. Construct whole-sample HSC pseudobulk Hallmark score
#      across the 8 biological samples.
#   6. Test TREM2-Mphi TGFB1 burden against the whole-sample
#      pseudobulk Hallmark score:
#        - Spearman
#        - exact permutation (8! = 40,320)
#        - partial Spearman controlling for NAS
#        - leave-one-sample-out
#   7. Compare cell-level ECM-high/nonhigh deltas with
#      pseudobulk ECM-high/nonhigh deltas sample-by-sample.
#
# STATISTICAL PRINCIPLE
#   Biological sample is the unit of inference.
#   Cell-level values are summarized per biological sample.
#
# GENE SET
#   MSigDB HALLMARK_TGF_BETA_SIGNALING (M5896), with
#   TGFB1, TGFBR1, and ECM-core overlaps removed.
#
# IMPORTANT
#   - No reclustering
#   - Canonical HSC receiver remains fixed
#   - Source RDS is never overwritten
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
  "Human_MASH_TGFbeta_Cell_PB_Concordance_v1.9"
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
  "Human_MASH_TGFbeta_Cell_PB_Concordance_v1.9.0.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Human MASH TGF-beta cell/pseudobulk concordance v1.9.0\n")
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

missing_hsc <- setdiff(
  required_hsc,
  colnames(hsc[[]])
)

if (length(missing_hsc) > 0) {
  stop(
    "Missing HSC metadata: ",
    paste(missing_hsc, collapse = ";")
  )
}

# ------------------------------------------------------------
# 2. Normalize and retrieve counts
# ------------------------------------------------------------

hsc <- NormalizeData(
  hsc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

counts <- GetAssayData(
  hsc,
  assay = "RNA",
  layer = "counts"
)

# ------------------------------------------------------------
# 3. Orthogonal Hallmark TGF-beta set
# ------------------------------------------------------------

HALLMARK_TGF_BETA_SIGNALING <- c(
  "ACVR1","APC","ARID4B","BCAR3","BMP2","BMPR1A","BMPR2","CDH1",
  "CDK9","CDKN1C","CTNNB1","ENG","FKBP1A","FNTA","FURIN","HDAC1",
  "HIPK2","ID1","ID2","ID3","IFNGR2","JUNB","KLF10","LEFTY2",
  "LTBP2","MAP3K7","NCOR2","NOG","PMEPA1","PPM1A","PPP1CA",
  "PPP1R15A","RAB31","RHOA","SERPINE1","SKI","SKIL","SLC20A1",
  "SMAD1","SMAD3","SMAD6","SMAD7","SMURF1","SMURF2","SPTBN1",
  "TGFB1","TGFBR1","TGIF1","THBS1","TJP1","TRIM33","UBE2D3",
  "WWTR1","XIAP"
)

ECM_CORE <- c(
  "COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","FN1","SPARC",
  "THBS2","LOXL2","TIMP1","TGFBI","CTHRC1","SERPINE1","POSTN"
)

EXCLUSIONS <- unique(
  c(
    "TGFB1",
    "TGFBR1",
    intersect(
      HALLMARK_TGF_BETA_SIGNALING,
      ECM_CORE
    )
  )
)

ORTHOGONAL_SET <- setdiff(
  HALLMARK_TGF_BETA_SIGNALING,
  EXCLUSIONS
)

genes <- intersect(
  ORTHOGONAL_SET,
  rownames(hsc)
)

if (length(genes) < 10) {
  warning(
    "Only ",
    length(genes),
    " orthogonal Hallmark genes are measurable."
  )
}

if (length(genes) < 5) {
  stop("Too few orthogonal Hallmark genes.")
}

cat("Orthogonal Hallmark genes available:", length(genes), "\n")
cat("Excluded:", paste(EXCLUSIONS, collapse = ";"), "\n\n")

# ------------------------------------------------------------
# 4. Technical covariates
# ------------------------------------------------------------

md <- hsc[[]]

nCount <- if ("nCount_RNA" %in% colnames(md)) {
  as.numeric(md$nCount_RNA)
} else {
  as.numeric(
    Matrix::colSums(counts)
  )
}

nFeature <- if ("nFeature_RNA" %in% colnames(md)) {
  as.numeric(md$nFeature_RNA)
} else {
  as.numeric(
    Matrix::colSums(counts > 0)
  )
}

expr <- FetchData(
  hsc,
  vars = genes,
  layer = "data"
)

count_gene <- t(
  as.matrix(
    counts[
      genes,
      ,
      drop = FALSE
    ]
  )
)

detected_hallmark_n <- rowSums(
  count_gene > 0,
  na.rm = TRUE
)

cell_meta <- tibble(
  cell = colnames(hsc),
  sample = as.character(md$sample),
  NAS = suppressWarnings(
    as.numeric(md$NAS_score)
  ),
  ECM_high = as.logical(
    md$ECM_high_primary_v1.2
  ),
  ECM_total_score = as.numeric(
    md$ECM_total_score_v1.2
  ),
  nCount_RNA = nCount,
  nFeature_RNA = nFeature,
  Hallmark_detected_n = detected_hallmark_n
)

# ------------------------------------------------------------
# 5. Cell-level score A: global per-gene z-score
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

expr_global_z <- as.data.frame(
  lapply(
    expr,
    z_fun
  )
)

score_global <- rowMeans(
  expr_global_z,
  na.rm = TRUE
)

# ------------------------------------------------------------
# 6. Cell-level score B: within-sample per-gene z-score
# ------------------------------------------------------------

expr_within_z <- matrix(
  NA_real_,
  nrow = nrow(expr),
  ncol = ncol(expr),
  dimnames = list(
    rownames(expr),
    colnames(expr)
  )
)

samples <- unique(
  cell_meta$sample
)

for (s in samples) {

  idx <- which(
    cell_meta$sample == s
  )

  for (g in colnames(expr)) {

    expr_within_z[
      idx,
      g
    ] <- z_fun(
      expr[
        idx,
        g
      ]
    )
  }
}

score_within_sample <- rowMeans(
  expr_within_z,
  na.rm = TRUE
)

# ------------------------------------------------------------
# 7. Cell-level score C:
#    residualize each gene against sample + depth + feature count
# ------------------------------------------------------------

sample_factor <- factor(
  cell_meta$sample
)

log_nCount <- log1p(
  cell_meta$nCount_RNA
)

log_nFeature <- log1p(
  cell_meta$nFeature_RNA
)

resid_mat <- matrix(
  NA_real_,
  nrow = nrow(expr),
  ncol = ncol(expr),
  dimnames = list(
    rownames(expr),
    colnames(expr)
  )
)

for (g in colnames(expr)) {

  dat <- data.frame(
    y = expr[[g]],
    sample_factor = sample_factor,
    log_nCount = log_nCount,
    log_nFeature = log_nFeature
  )

  fit <- tryCatch(
    lm(
      y ~ sample_factor + log_nCount + log_nFeature,
      data = dat
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    resid_mat[, g] <- 0
  } else {
    r <- residuals(fit)
    resid_mat[, g] <- z_fun(r)
  }
}

score_residualized <- rowMeans(
  resid_mat,
  na.rm = TRUE
)

cell_meta$Hallmark_global_z <- score_global
cell_meta$Hallmark_within_sample_z <- score_within_sample
cell_meta$Hallmark_residualized <- score_residualized
cell_meta$Hallmark_mean_logexpr <- rowMeans(
  as.matrix(expr),
  na.rm = TRUE
)

write_tsv(
  cell_meta,
  file.path(
    DIR_TABLE,
    "01_cell_level_scores_and_technical_covariates.tsv.gz"
  )
)

# ------------------------------------------------------------
# 8. Technical audit by sample x ECM group
# ------------------------------------------------------------

technical_by_group <- cell_meta %>%
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

    mean_nCount_RNA =
      mean(
        nCount_RNA,
        na.rm = TRUE
      ),

    median_nCount_RNA =
      median(
        nCount_RNA,
        na.rm = TRUE
      ),

    mean_nFeature_RNA =
      mean(
        nFeature_RNA,
        na.rm = TRUE
      ),

    median_nFeature_RNA =
      median(
        nFeature_RNA,
        na.rm = TRUE
      ),

    mean_Hallmark_detected_n =
      mean(
        Hallmark_detected_n,
        na.rm = TRUE
      ),

    Hallmark_global_z =
      mean(
        Hallmark_global_z,
        na.rm = TRUE
      ),

    Hallmark_within_sample_z =
      mean(
        Hallmark_within_sample_z,
        na.rm = TRUE
      ),

    Hallmark_residualized =
      mean(
        Hallmark_residualized,
        na.rm = TRUE
      ),

    Hallmark_mean_logexpr =
      mean(
        Hallmark_mean_logexpr,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  technical_by_group,
  file.path(
    DIR_TABLE,
    "02_sample_ECMgroup_technical_and_scores.tsv"
  )
)

# ------------------------------------------------------------
# 9. Paired tests of technical covariates and cell-level scores
# ------------------------------------------------------------

wide <- technical_by_group %>%
  pivot_wider(
    names_from = ECM_group,
    values_from = c(
      cells,
      mean_nCount_RNA,
      median_nCount_RNA,
      mean_nFeature_RNA,
      median_nFeature_RNA,
      mean_Hallmark_detected_n,
      Hallmark_global_z,
      Hallmark_within_sample_z,
      Hallmark_residualized,
      Hallmark_mean_logexpr
    )
  )

paired_metrics <- c(
  "mean_nCount_RNA",
  "median_nCount_RNA",
  "mean_nFeature_RNA",
  "median_nFeature_RNA",
  "mean_Hallmark_detected_n",
  "Hallmark_global_z",
  "Hallmark_within_sample_z",
  "Hallmark_residualized",
  "Hallmark_mean_logexpr"
)

paired_test <- function(metric) {

  x <- wide[[paste0(
    metric,
    "_ECM_high"
  )]]

  y <- wide[[paste0(
    metric,
    "_ECM_nonhigh"
  )]]

  keep <- complete.cases(
    x,
    y
  )

  d <- x[keep] - y[keep]

  p <- if (
    length(d) >= 4 &&
    !all(d == 0)
  ) {
    suppressWarnings(
      wilcox.test(
        d,
        mu = 0,
        exact = FALSE
      )$p.value
    )
  } else {
    NA_real_
  }

  tibble(
    metric = metric,
    n_samples = length(d),
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
    p_value = p
  )
}

paired_results <- bind_rows(
  lapply(
    paired_metrics,
    paired_test
  )
) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  )

write_tsv(
  paired_results,
  file.path(
    DIR_TABLE,
    "03_PRIMARY_cell_level_and_technical_paired_tests.tsv"
  )
)

# ------------------------------------------------------------
# 10. Pseudobulk sample x ECM group
# ------------------------------------------------------------

group_id <- paste(
  cell_meta$sample,
  ifelse(
    cell_meta$ECM_high,
    "ECM_high",
    "ECM_nonhigh"
  ),
  sep = "__"
)

groups <- unique(
  group_id
)

pb_list <- lapply(
  groups,
  function(grp) {

    idx <- which(
      group_id == grp
    )

    Matrix::rowSums(
      counts[
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

rownames(pb_counts) <- rownames(counts)
colnames(pb_counts) <- groups

pb_lib <- colSums(
  pb_counts
)

pb_logcpm <- log2(
  sweep(
    pb_counts,
    2,
    pb_lib,
    "/"
  ) * 1e6 + 1
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
    cell_meta %>%
      distinct(
        sample,
        NAS
      ),
    by = "sample"
  ) %>%
  mutate(
    library_size =
      pb_lib[pseudobulk]
  )

write_tsv(
  pb_meta,
  file.path(
    DIR_TABLE,
    "04_ECMgroup_pseudobulk_metadata.tsv"
  )
)

# ------------------------------------------------------------
# 11. Gene-wise pseudobulk paired deltas
# ------------------------------------------------------------

pb_gene_long <- as_tibble(
  t(
    pb_logcpm[
      genes,
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
    cols = all_of(genes),
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
    delta =
      ECM_high -
      ECM_nonhigh
  )

write_tsv(
  pb_gene_wide,
  file.path(
    DIR_TABLE,
    "05_pseudobulk_gene_deltas_by_sample.tsv"
  )
)

pb_sample_delta <- pb_gene_wide %>%
  group_by(
    sample,
    NAS
  ) %>%
  summarise(
    mean_gene_delta =
      mean(
        delta,
        na.rm = TRUE
      ),

    median_gene_delta =
      median(
        delta,
        na.rm = TRUE
      ),

    positive_gene_fraction =
      mean(
        delta > 0,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_tsv(
  pb_sample_delta,
  file.path(
    DIR_TABLE,
    "06_PRIMARY_pseudobulk_pathway_delta_by_sample.tsv"
  )
)

pb_delta_test <- bind_rows(
  {
    d <- pb_sample_delta$mean_gene_delta

    tibble(
      metric = "mean_gene_delta",
      n_samples = sum(is.finite(d)),
      median_across_samples =
        median(
          d,
          na.rm = TRUE
        ),
      positive_samples =
        sum(
          d > 0,
          na.rm = TRUE
        ),
      p_value =
        suppressWarnings(
          wilcox.test(
            d,
            mu = 0,
            exact = FALSE
          )$p.value
        )
    )
  },

  {
    d <- pb_sample_delta$median_gene_delta

    tibble(
      metric = "median_gene_delta",
      n_samples = sum(is.finite(d)),
      median_across_samples =
        median(
          d,
          na.rm = TRUE
        ),
      positive_samples =
        sum(
          d > 0,
          na.rm = TRUE
        ),
      p_value =
        suppressWarnings(
          wilcox.test(
            d,
            mu = 0,
            exact = FALSE
          )$p.value
        )
    )
  }
)

write_tsv(
  pb_delta_test,
  file.path(
    DIR_TABLE,
    "07_PRIMARY_pseudobulk_pathway_delta_paired_tests.tsv"
  )
)

# ------------------------------------------------------------
# 12. Whole-sample HSC pseudobulk
# ------------------------------------------------------------

sample_ids <- unique(
  cell_meta$sample
)

sample_pb_list <- lapply(
  sample_ids,
  function(s) {

    idx <- which(
      cell_meta$sample == s
    )

    Matrix::rowSums(
      counts[
        ,
        idx,
        drop = FALSE
      ]
    )
  }
)

sample_pb_counts <- do.call(
  cbind,
  sample_pb_list
)

rownames(sample_pb_counts) <- rownames(counts)
colnames(sample_pb_counts) <- sample_ids

sample_pb_lib <- colSums(
  sample_pb_counts
)

sample_pb_logcpm <- log2(
  sweep(
    sample_pb_counts,
    2,
    sample_pb_lib,
    "/"
  ) * 1e6 + 1
)

sample_gene_mat <- t(
  sample_pb_logcpm[
    genes,
    ,
    drop = FALSE
  ]
)

sample_gene_z <- apply(
  sample_gene_mat,
  2,
  z_fun
)

if (is.vector(sample_gene_z)) {
  sample_gene_z <- matrix(
    sample_gene_z,
    ncol = 1
  )
}

sample_pb_score <- rowMeans(
  sample_gene_z,
  na.rm = TRUE
)

sample_pb <- tibble(
  sample = rownames(
    sample_gene_mat
  ),
  Hallmark_whole_HSC_pseudobulk_score =
    sample_pb_score
) %>%
  left_join(
    cell_meta %>%
      distinct(
        sample,
        NAS
      ),
    by = "sample"
  ) %>%
  arrange(sample)

write_tsv(
  sample_pb,
  file.path(
    DIR_TABLE,
    "08_whole_sample_HSC_pseudobulk_Hallmark_score.tsv"
  )
)

# ------------------------------------------------------------
# 13. Statistical helpers
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

    out[[i]] <- cbind(
      x[i],
      all_permutations(
        x[-i]
      )
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

  r_xy <- cor(
    x,
    y,
    method = "spearman"
  )

  r_xz <- cor(
    x,
    z,
    method = "spearman"
  )

  r_yz <- cor(
    y,
    z,
    method = "spearman"
  )

  denom <- sqrt(
    (1 - r_xz^2) *
      (1 - r_yz^2)
  )

  if (
    !is.finite(denom) ||
    denom <= 0 ||
    n < 5
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

  tval <- rp *
    sqrt(
      df /
        (1 - rp^2)
    )

  tibble(
    n = n,
    partial_rho = rp,
    p_value =
      2 *
      pt(
        -abs(tval),
        df = df
      )
  )
}

# ------------------------------------------------------------
# 14. TREM2 TGFB1 vs whole-sample HSC pseudobulk score
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
    sample_pb,
    by = c(
      "sample",
      "NAS"
    )
  ) %>%
  arrange(sample)

primary_pb_spearman <- safe_spearman(
  trem2$TGFB1_burden_expr,
  trem2$Hallmark_whole_HSC_pseudobulk_score
)

primary_pb_exact <- exact_perm(
  trem2$TGFB1_burden_expr,
  trem2$Hallmark_whole_HSC_pseudobulk_score
)

primary_pb_partial <- partial_spearman(
  trem2$TGFB1_burden_expr,
  trem2$Hallmark_whole_HSC_pseudobulk_score,
  trem2$NAS
)

write_tsv(
  primary_pb_spearman,
  file.path(
    DIR_TABLE,
    "09_PRIMARY_TREM2_vs_wholeHSC_pseudobulk_Spearman.tsv"
  )
)

write_tsv(
  primary_pb_exact,
  file.path(
    DIR_TABLE,
    "10_PRIMARY_TREM2_vs_wholeHSC_pseudobulk_exact.tsv"
  )
)

write_tsv(
  primary_pb_partial,
  file.path(
    DIR_TABLE,
    "11_PRIMARY_TREM2_vs_wholeHSC_pseudobulk_partial_NAS.tsv"
  )
)

# ------------------------------------------------------------
# 15. Leave-one-out whole-sample pseudobulk association
# ------------------------------------------------------------

loo <- bind_rows(
  lapply(
    seq_len(nrow(trem2)),
    function(i) {

      sub <- trem2[-i, ]

      safe_spearman(
        sub$TGFB1_burden_expr,
        sub$Hallmark_whole_HSC_pseudobulk_score
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
    "12_leave_one_out_wholeHSC_pseudobulk.tsv"
  )
)

# ------------------------------------------------------------
# 16. Cell-level vs pseudobulk delta concordance
# ------------------------------------------------------------

cell_delta <- wide %>%
  transmute(
    sample,
    NAS,

    delta_global =
      Hallmark_global_z_ECM_high -
      Hallmark_global_z_ECM_nonhigh,

    delta_within_sample =
      Hallmark_within_sample_z_ECM_high -
      Hallmark_within_sample_z_ECM_nonhigh,

    delta_residualized =
      Hallmark_residualized_ECM_high -
      Hallmark_residualized_ECM_nonhigh,

    delta_mean_logexpr =
      Hallmark_mean_logexpr_ECM_high -
      Hallmark_mean_logexpr_ECM_nonhigh
  ) %>%
  left_join(
    pb_sample_delta,
    by = c(
      "sample",
      "NAS"
    )
  )

write_tsv(
  cell_delta,
  file.path(
    DIR_TABLE,
    "13_cell_vs_pseudobulk_delta_by_sample.tsv"
  )
)

concordance <- bind_rows(
  safe_spearman(
    cell_delta$delta_global,
    cell_delta$mean_gene_delta
  ) %>%
    mutate(
      cell_metric = "delta_global",
      pb_metric = "mean_gene_delta",
      .before = 1
    ),

  safe_spearman(
    cell_delta$delta_within_sample,
    cell_delta$mean_gene_delta
  ) %>%
    mutate(
      cell_metric = "delta_within_sample",
      pb_metric = "mean_gene_delta",
      .before = 1
    ),

  safe_spearman(
    cell_delta$delta_residualized,
    cell_delta$mean_gene_delta
  ) %>%
    mutate(
      cell_metric = "delta_residualized",
      pb_metric = "mean_gene_delta",
      .before = 1
    ),

  safe_spearman(
    cell_delta$delta_mean_logexpr,
    cell_delta$mean_gene_delta
  ) %>%
    mutate(
      cell_metric = "delta_mean_logexpr",
      pb_metric = "mean_gene_delta",
      .before = 1
    )
)

write_tsv(
  concordance,
  file.path(
    DIR_TABLE,
    "14_cell_vs_pseudobulk_delta_concordance.tsv"
  )
)

# ------------------------------------------------------------
# 17. Gene-level pseudobulk summary
# ------------------------------------------------------------

pb_gene_summary <- pb_gene_wide %>%
  group_by(gene) %>%
  summarise(
    n_samples =
      sum(
        is.finite(delta)
      ),

    mean_delta =
      mean(
        delta,
        na.rm = TRUE
      ),

    median_delta =
      median(
        delta,
        na.rm = TRUE
      ),

    positive_samples =
      sum(
        delta > 0,
        na.rm = TRUE
      ),

    p_value = {
      d <- delta[
        is.finite(delta)
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
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(
    FDR,
    desc(mean_delta)
  )

write_tsv(
  pb_gene_summary,
  file.path(
    DIR_TABLE,
    "15_pseudobulk_gene_summary.tsv"
  )
)

# ------------------------------------------------------------
# 18. Figures
# ------------------------------------------------------------

p1 <- ggplot(
  technical_by_group,
  aes(
    x = ECM_group,
    y = mean_nFeature_RNA,
    group = sample
  )
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "RNA feature complexity by ECM functional state",
    x = NULL,
    y = "Mean nFeature_RNA"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_nFeature_ECMhigh_vs_nonhigh.pdf"
  ),
  p1,
  width = 6.5,
  height = 5.5
)

p2 <- ggplot(
  trem2,
  aes(
    x = TGFB1_burden_expr,
    y = Hallmark_whole_HSC_pseudobulk_score,
    label = sample
  )
) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.03,
    check_overlap = TRUE
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "TREM2-Mphi TGFB1 vs whole-HSC pseudobulk TGF-beta score",
    x = "TREM2-Mphi TGFB1 burden",
    y = "Whole-HSC pseudobulk Hallmark score"
  )

ggsave(
  file.path(
    DIR_FIG,
    "02_TREM2_vs_wholeHSC_pseudobulk_TGFbeta.pdf"
  ),
  p2,
  width = 7,
  height = 5.5
)

p3 <- ggplot(
  cell_delta,
  aes(
    x = delta_residualized,
    y = mean_gene_delta,
    label = sample
  )
) +
  geom_point(size = 3) +
  geom_text(
    nudge_y = 0.03,
    check_overlap = TRUE
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Cell-level vs pseudobulk TGF-beta delta",
    x = "Depth/sample-residualized cell-level delta",
    y = "Pseudobulk mean gene delta"
  )

ggsave(
  file.path(
    DIR_FIG,
    "03_cell_vs_pseudobulk_delta.pdf"
  ),
  p3,
  width = 7,
  height = 5.5
)

# ------------------------------------------------------------
# 19. Save results
# ------------------------------------------------------------

saveRDS(
  list(
    genes = genes,
    technical_by_group = technical_by_group,
    paired_results = paired_results,
    pb_meta = pb_meta,
    pb_sample_delta = pb_sample_delta,
    pb_delta_test = pb_delta_test,
    sample_pb = sample_pb,
    trem2 = trem2,
    primary_pb_spearman = primary_pb_spearman,
    primary_pb_exact = primary_pb_exact,
    primary_pb_partial = primary_pb_partial,
    loo = loo,
    cell_delta = cell_delta,
    concordance = concordance,
    pb_gene_summary = pb_gene_summary
  ),
  file.path(
    DIR_OBJ,
    "Human_MASH_TGFbeta_Cell_PB_Concordance_v1.9_results.rds"
  )
)

# ------------------------------------------------------------
# 20. Provenance
# ------------------------------------------------------------

provenance <- tibble(
  item = c(
    "analysis_version",
    "HSC_source_RDS",
    "sender_source_table",
    "HSC_cells",
    "orthogonal_Hallmark_genes_n",
    "circularity_exclusions",
    "cell_score_A",
    "cell_score_B",
    "cell_score_C",
    "pseudobulk_group_definition",
    "whole_sample_pseudobulk_definition",
    "association_unit",
    "source_overwritten",
    "reclustering"
  ),
  value = c(
    "v1.9.0",
    HSC_RDS,
    SENDER_TABLE,
    as.character(
      ncol(hsc)
    ),
    as.character(
      length(genes)
    ),
    paste(
      EXCLUSIONS,
      collapse = ";"
    ),
    "global_per_gene_z_on_log_normalized_expression",
    "within_sample_per_gene_z",
    "per_gene_residuals_from_sample_plus_log_nCount_plus_log_nFeature_then_z",
    "sample_x_ECM_high_status_sum_counts_logCPM",
    "sample_sum_counts_logCPM_gene_z_across_8_samples",
    "biological_sample",
    "FALSE",
    "none"
  )
)

write_tsv(
  provenance,
  file.path(
    DIR_PROV,
    "concordance_provenance_v1.9.tsv"
  )
)

# ------------------------------------------------------------
# 21. Integrity
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "HSC_input_cells",
    "HSC_expected_cells",
    "HSC_count_matches",
    "orthogonal_genes_n",
    "all_8_samples_present",
    "ECMgroup_pseudobulk_profiles",
    "expected_ECMgroup_pseudobulk_profiles",
    "whole_sample_pseudobulk_profiles",
    "expected_whole_sample_profiles",
    "exact_permutation_completed",
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
      length(genes)
    ),
    as.character(
      length(
        unique(
          cell_meta$sample
        )
      ) == 8
    ),
    as.character(
      nrow(pb_meta)
    ),
    "16",
    as.character(
      nrow(sample_pb)
    ),
    "8",
    as.character(
      primary_pb_exact$permutations ==
        40320
    ),
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v1.9.tsv"
  )
)

# ------------------------------------------------------------
# 22. Console summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CELL-LEVEL / TECHNICAL PAIRED TESTS\n")
cat("============================================================\n")
print(paired_results)

cat("\n")
cat("============================================================\n")
cat("PSEUDOBULK PATHWAY DELTA BY SAMPLE\n")
cat("============================================================\n")
print(pb_sample_delta)
print(pb_delta_test)

cat("\n")
cat("============================================================\n")
cat("TREM2 TGFB1 -> WHOLE-HSC PSEUDOBULK HALLMARK\n")
cat("============================================================\n")
print(primary_pb_spearman)
print(primary_pb_exact)
print(primary_pb_partial)

cat("\n")
cat("============================================================\n")
cat("LEAVE-ONE-OUT WHOLE-HSC PSEUDOBULK\n")
cat("============================================================\n")
print(loo)

cat("\n")
cat("============================================================\n")
cat("CELL vs PSEUDOBULK DELTA CONCORDANCE\n")
cat("============================================================\n")
print(concordance)

cat("\n")
cat("============================================================\n")
cat("TOP PSEUDOBULK GENES\n")
cat("============================================================\n")
print(
  pb_gene_summary %>%
    head(30)
)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
