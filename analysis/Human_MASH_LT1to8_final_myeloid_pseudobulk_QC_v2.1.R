#!/usr/bin/env Rscript

# Human_MASH_LT1to8_final_myeloid_pseudobulk_QC_v2.1.R
#
# PURPOSE
# -------
# Build sample-level pseudobulk count matrices for each FINAL Human MΦ/Monocyte
# state and audit whether state-specific pseudobulk DE is technically feasible.
#
# INPUT:
#   Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds
#
# THIS SCRIPT DOES:
#   - aggregate RAW RNA counts by biological sample within each final state
#   - save pseudobulk count matrices as a single RDS list
#   - report cells/library size/detected genes for every state x sample
#   - calculate descriptive logCPM sample correlations
#   - calculate descriptive PCA coordinates within each state
#   - flag small state x sample groups
#
# THIS SCRIPT DOES NOT:
#   - perform DE
#   - perform pathway analysis
#   - modify/save a Seurat object
#   - use cells as statistical replicates
#
# Statistical unit:
#   biological sample (LT1-LT8), never individual cells.
#
# Seurat v5:
#   GetAssayData(..., layer="counts"), never slot=.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v2.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9/",
  "objects/Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloid_PseudobulkQC_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

PB_RDS <- file.path(
  OBJDIR,
  paste0("Human_MASH_LT1to8_FinalMyeloid_pseudobulk_counts_", VERSION, ".rds")
)

if (!file.exists(INPUT_RDS)) stop("Input v1.9 RDS not found:\n", INPUT_RDS)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

message("============================================================")
message("Human MASH LT1-LT8 final myeloid pseudobulk QC ", VERSION)
message("============================================================")
message("NO DE; biological-sample aggregation only.")

obj <- readRDS(INPUT_RDS)

required_meta <- c(
  "myeloid_state_FINAL_v19",
  "sample_v17",
  "NAS_v17"
)

miss <- setdiff(required_meta, colnames(obj[[]]))
if (length(miss)) stop("Missing required metadata: ", paste(miss, collapse = ", "))

if (!"RNA" %in% Assays(obj)) stop("RNA assay absent.")
if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA counts layer absent.")

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")

sample_id <- as.character(obj$sample_v17)
nas <- as.numeric(obj$NAS_v17)
state <- as.character(obj$myeloid_state_FINAL_v19)

SAMPLE_ORDER <- c("LT7", "LT8", "LT5", "LT6", "LT3", "LT4", "LT1", "LT2")
EXPECTED_NAS <- c(
  LT7 = 0, LT8 = 0,
  LT5 = 2, LT6 = 2,
  LT3 = 4, LT4 = 4,
  LT1 = 6, LT2 = 6
)

STATE_ORDER <- c(
  "Resident/Kupffer-like_Mphi",
  "Repair/Resolution_Mphi",
  "Lipid-associated/TREM2_Mphi",
  "Inflammatory/Transitional_Mphi",
  "Inflammatory_Monocyte"
)

if (!setequal(unique(sample_id), SAMPLE_ORDER)) {
  stop("Unexpected sample universe.")
}
if (!setequal(unique(state), STATE_ORDER)) {
  stop("Unexpected final-state universe.")
}

# ------------------------------------------------------------------
# Aggregate RAW counts by state x biological sample.
# ------------------------------------------------------------------

pb_list <- list()
qc_rows <- list()
count_manifest <- list()

for (st in STATE_ORDER) {
  message("Aggregating state: ", st)

  mat_state <- counts[, state == st, drop = FALSE]
  sample_state <- sample_id[state == st]

  pb <- matrix(
    0,
    nrow = nrow(mat_state),
    ncol = length(SAMPLE_ORDER),
    dimnames = list(rownames(mat_state), SAMPLE_ORDER)
  )

  for (s in SAMPLE_ORDER) {
    idx <- which(sample_state == s)

    if (length(idx)) {
      pb[, s] <- Matrix::rowSums(mat_state[, idx, drop = FALSE])
    }

    libsize <- sum(pb[, s])
    genes_detected <- sum(pb[, s] > 0)

    qc_rows[[length(qc_rows) + 1L]] <- data.frame(
      final_state = st,
      sample = s,
      NAS = unname(EXPECTED_NAS[[s]]),
      n_cells = length(idx),
      pseudobulk_library_size = libsize,
      genes_detected = genes_detected,
      flag_cells_lt20 = length(idx) < 20,
      flag_cells_lt30 = length(idx) < 30,
      stringsAsFactors = FALSE
    )
  }

  pb_list[[st]] <- pb

  count_manifest[[st]] <- data.frame(
    final_state = st,
    n_genes = nrow(pb),
    n_samples = ncol(pb),
    total_counts = sum(pb),
    stringsAsFactors = FALSE
  )
}

saveRDS(pb_list, PB_RDS, compress = FALSE)

qc <- do.call(rbind, qc_rows)
rownames(qc) <- NULL

write.csv(
  qc,
  file.path(TABDIR, "01_pseudobulk_library_QC_v2.1.csv"),
  row.names = FALSE
)

manifest <- do.call(rbind, count_manifest)
rownames(manifest) <- NULL

write.csv(
  manifest,
  file.path(TABDIR, "02_pseudobulk_count_manifest_v2.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Descriptive sample correlations and PCA within each state.
# No statistical inference is performed here.
# ------------------------------------------------------------------

cor_rows <- list()
pca_rows <- list()

for (st in STATE_ORDER) {
  pb <- pb_list[[st]]

  libsize <- colSums(pb)

  # logCPM, with no between-sample normalization factor.
  # This is QC only, not a DE input transformation.
  logcpm <- log2(
    sweep(pb, 2, libsize, "/") * 1e6 + 1
  )

  # Remove genes with no variation across the eight pseudobulks.
  gene_var <- apply(logcpm, 1, var)
  gene_var[!is.finite(gene_var)] <- 0

  keep <- gene_var > 0

  if (sum(keep) < 50) {
    warning("Very few variable genes for state: ", st)
    next
  }

  # Sample correlation using all variable genes.
  cr <- cor(logcpm[keep, , drop = FALSE], method = "pearson")

  for (i in seq_along(SAMPLE_ORDER)) {
    for (j in seq_along(SAMPLE_ORDER)) {
      cor_rows[[length(cor_rows) + 1L]] <- data.frame(
        final_state = st,
        sample1 = SAMPLE_ORDER[i],
        sample2 = SAMPLE_ORDER[j],
        correlation = cr[SAMPLE_ORDER[i], SAMPLE_ORDER[j]],
        stringsAsFactors = FALSE
      )
    }
  }

  # PCA on top 2,000 variable genes or all available if fewer.
  top_n <- min(2000L, sum(keep))
  top_genes <- names(
    sort(gene_var[keep], decreasing = TRUE)
  )[seq_len(top_n)]

  x <- t(logcpm[top_genes, , drop = FALSE])
  pc <- prcomp(x, center = TRUE, scale. = FALSE)

  var_exp <- (pc$sdev^2) / sum(pc$sdev^2)

  pca_rows[[st]] <- data.frame(
    final_state = st,
    sample = rownames(pc$x),
    NAS = unname(EXPECTED_NAS[rownames(pc$x)]),
    PC1 = pc$x[, 1],
    PC2 = pc$x[, 2],
    PC1_variance = var_exp[1],
    PC2_variance = var_exp[2],
    stringsAsFactors = FALSE
  )
}

cor_df <- do.call(rbind, cor_rows)
rownames(cor_df) <- NULL

write.csv(
  cor_df,
  file.path(TABDIR, "03_pseudobulk_sample_correlations_v2.1.csv"),
  row.names = FALSE
)

pca_df <- do.call(rbind, pca_rows)
rownames(pca_df) <- NULL

write.csv(
  pca_df,
  file.path(TABDIR, "04_pseudobulk_PCA_coordinates_v2.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------

# Faceted PCA by final state.
pca_df$final_state <- factor(pca_df$final_state, levels = STATE_ORDER)
pca_df$NAS_factor <- factor(pca_df$NAS, levels = c(0, 2, 4, 6))

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = NAS_factor,
    label = sample
  )
) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  facet_wrap(~ final_state, scales = "free", ncol = 2) +
  theme_classic(base_size = 10) +
  theme(legend.title = element_blank()) +
  labs(
    x = "PC1",
    y = "PC2"
  ) +
  ggtitle(
    "Human MASH final myeloid-state pseudobulk PCA",
    "QC only; 8 biological samples per state"
  )

ggsave(
  file.path(FIGDIR, "01_final_state_pseudobulk_PCA_v2.1.pdf"),
  p_pca,
  width = 10.5,
  height = 10
)

# Library-size QC.
qc$final_state <- factor(qc$final_state, levels = STATE_ORDER)
qc$sample <- factor(qc$sample, levels = SAMPLE_ORDER)

p_lib <- ggplot(
  qc,
  aes(
    x = sample,
    y = pseudobulk_library_size,
    fill = final_state
  )
) +
  geom_col() +
  facet_wrap(~ final_state, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 10) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "Pseudobulk library size") +
  ggtitle(
    "Human MASH final myeloid-state pseudobulk library QC",
    "Biological samples ordered NAS0 -> NAS2 -> NAS4 -> NAS6"
  )

ggsave(
  file.path(FIGDIR, "02_pseudobulk_library_size_QC_v2.1.pdf"),
  p_lib,
  width = 10.5,
  height = 9
)

# ------------------------------------------------------------------
# Feasibility summary
# ------------------------------------------------------------------

state_feas <- do.call(
  rbind,
  lapply(STATE_ORDER, function(st) {
    z <- qc[as.character(qc$final_state) == st, , drop = FALSE]

    data.frame(
      final_state = st,
      min_cells_per_sample = min(z$n_cells),
      max_cells_per_sample = max(z$n_cells),
      min_library_size = min(z$pseudobulk_library_size),
      min_genes_detected = min(z$genes_detected),
      n_samples_lt20_cells = sum(z$n_cells < 20),
      n_samples_lt30_cells = sum(z$n_cells < 30),
      all_8_samples_present = all(z$n_cells > 0),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  state_feas,
  file.path(TABDIR, "05_state_pseudobulk_feasibility_v2.1.csv"),
  row.names = FALSE
)

summary_lines <- c(
  paste0("Human MASH LT1-LT8 final myeloid pseudobulk QC ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "STATUS: SUCCESS",
  "",
  paste0("Input v1.9 RDS: ", INPUT_RDS),
  paste0("Pseudobulk counts RDS: ", PB_RDS),
  "",
  "STATISTICAL UNIT:",
  "Biological sample (LT1-LT8), never cells.",
  "",
  "STATE FEASIBILITY:",
  paste(
    capture.output(print(state_feas, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "IMPORTANT:",
  "- v2.1 performs aggregation/QC only; NO DE was performed.",
  "- All final states are represented in all 8 biological samples.",
  "- Small state x sample groups are explicitly flagged.",
  "- logCPM/PCA/correlations are descriptive QC only.",
  "",
  "NEXT DECISION GATE:",
  "Review library sizes, sample PCA, and minimum cells per state/sample.",
  "Then choose which final states are robust enough for sample-level pseudobulk DE."
)

writeLines(
  summary_lines,
  file.path(TABDIR, "06_pseudobulk_QC_summary_v2.1.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v2.1 COMPLETE")
message("Pseudobulk matrices saved: ", PB_RDS)
message("NO DE performed.")
message("============================================================")
