#!/usr/bin/env Rscript

# ==============================================================================
# Human_MASH_v1.4.4A.2_ProbeAware_LigandTransitionAudit.R
#
# Purpose
#   Probe/molecule-aware descriptive audit of candidate human graft ligands:
#     SfLB1/SfLB2 = pre-transplant, non-engineered cultured organoids
#     Tx17/Tx5    = engrafted human cells derived from IL10-engineered organoids
#
# Primary quantitative evidence
#   v1.4.3B.2C.3 molecule reconstruction
#   - strict shared-TRUE probe sensitivity for cross-sample comparison
#   - standard MEX molecule/cell evidence for IL10 anchor resolution
#
# Causal boundary
#   SfLB -> Tx is NOT an IL10-only perturbation. It combines:
#     transgene status + transplantation/engraftment + host microenvironment
#     + post-engraftment selection/state transition
#     + possible residual innate xenogeneic immune pressure.
#
#   v1.4.4A.1.1 showed major lineage-composition shifts. Therefore,
#   whole-graft ligand changes are descriptive and MUST NOT be interpreted
#   as same-lineage transcriptional induction. Source-lineage decomposition
#   is deferred to v1.4.4A.3.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

B2C3_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.3B.2C.3/",
  "analysis/v1.4.3B.2C.3_prioritized_ligand_molecule_reconstruction/tables"
)

STRICT_TSV <- file.path(
  B2C3_DIR,
  "v1.4.3B.2C.3_crosssample_strict_shared_TRUE_summary.tsv"
)

GENE_TSV <- file.path(
  B2C3_DIR,
  "v1.4.3B.2C.3_gene_reconstruction_by_sample.tsv"
)

IL10_TSV <- file.path(
  B2C3_DIR,
  "v1.4.3B.2C.3_IL10_resolution.tsv"
)

OUT_ROOT <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.2/",
  "analysis/v1.4.4A.2_probe_aware_ligand_transition_audit"
)

TAB_DIR <- file.path(OUT_ROOT, "tables")
FIG_DIR <- file.path(OUT_ROOT, "figures")
LOG_DIR <- file.path(OUT_ROOT, "logs")

dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, "v1.4.4A.2_run_log.txt")
sink(LOG_FILE, split = TRUE)
on.exit(sink(), add = TRUE)

cat("============================================================\n")
cat("Human MASH v1.4.4A.2\n")
cat("Probe-aware ligand transition audit\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("B.2C.3:", B2C3_DIR, "\n")
cat("Output:", OUT_ROOT, "\n\n")

for (f in c(STRICT_TSV, GENE_TSV, IL10_TSV)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

# ------------------------------------------------------------------------------
# 1. Robust TSV readers
#    Only use the validated leading fields needed for this analysis.
#    This avoids dependence on any malformed/concatenated trailing header fields.
# ------------------------------------------------------------------------------

read_first_n_fields <- function(path, n_fields, col_names) {
  x <- readLines(path, warn = FALSE)
  if (length(x) < 2) stop("No data rows in ", path)

  rows <- strsplit(x[-1], "\t", fixed = TRUE)

  bad <- which(lengths(rows) < n_fields)
  if (length(bad) > 0) {
    stop(
      "Rows with fewer than ", n_fields, " fields in ", basename(path),
      ": ", paste(head(bad + 1L, 20), collapse = ", ")
    )
  }

  mat <- do.call(
    rbind,
    lapply(rows, function(z) z[seq_len(n_fields)])
  )

  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(out) <- col_names
  out
}

strict <- read_first_n_fields(
  STRICT_TSV,
  12,
  c(
    "gene",
    "SfLB1_MEX_UMI",
    "SfLB2_MEX_UMI",
    "Tx17_MEX_UMI",
    "Tx5_MEX_UMI",
    "SfLB1_shared_TRUE_UMI",
    "SfLB2_shared_TRUE_UMI",
    "Tx17_shared_TRUE_UMI",
    "Tx5_shared_TRUE_UMI",
    "strict_Tx_mean_UMI_per_cell",
    "strict_SfLB_mean_UMI_per_cell",
    "strict_Tx_vs_SfLB_log2FC_per_cell_descriptive"
  )
)

gene_rec <- read_first_n_fields(
  GENE_TSV,
  5,
  c(
    "sample",
    "gene",
    "canonical_cells",
    "MEX_UMI",
    "MEX_positive_cells"
  )
)

il10_res <- read_first_n_fields(
  IL10_TSV,
  5,
  c(
    "sample",
    "gene",
    "canonical_cells",
    "MEX_UMI",
    "MEX_positive_cells"
  )
)

num_cols_strict <- setdiff(names(strict), "gene")
for (cc in num_cols_strict) strict[[cc]] <- suppressWarnings(as.numeric(strict[[cc]]))

for (cc in c("canonical_cells", "MEX_UMI", "MEX_positive_cells")) {
  gene_rec[[cc]] <- suppressWarnings(as.numeric(gene_rec[[cc]]))
  il10_res[[cc]] <- suppressWarnings(as.numeric(il10_res[[cc]]))
}

# ------------------------------------------------------------------------------
# 2. Candidate and sample validation
# ------------------------------------------------------------------------------

candidate_order <- c(
  "IL10",
  "FGF1",
  "TNFSF12",
  "VEGFA",
  "CXCL12",
  "TGFB3",
  "BMP5",
  "BMP2",
  "HBEGF",
  "BTC",
  "FST",
  "IGFBP2",
  "IGFBP4",
  "INHBA",
  "VEGFB"
)

sample_order <- c("SfLB1", "SfLB2", "Tx17", "Tx5")

missing_strict <- setdiff(candidate_order, strict$gene)
if (length(missing_strict) > 0) {
  stop("Missing candidate genes from strict summary: ",
       paste(missing_strict, collapse = ", "))
}

gene_rec <- gene_rec[
  gene_rec$gene %in% candidate_order &
    gene_rec$sample %in% sample_order,
  ,
  drop = FALSE
]

if (!setequal(unique(gene_rec$sample), sample_order)) {
  stop("Unexpected sample set in gene reconstruction.")
}

if (!setequal(unique(gene_rec$gene), candidate_order)) {
  stop("Candidate gene mismatch in gene reconstruction.")
}

# Each candidate x sample should be unique.
dup <- duplicated(gene_rec[, c("sample", "gene")])
if (any(dup)) stop("Duplicated sample-gene rows in gene reconstruction.")

expected_n <- length(candidate_order) * length(sample_order)
if (nrow(gene_rec) != expected_n) {
  stop("Unexpected sample-gene row count: ", nrow(gene_rec),
       "; expected ", expected_n)
}

# IL10 resolution must exactly agree with gene reconstruction for first 5 fields.
il10_gene <- gene_rec[gene_rec$gene == "IL10", ]
il10_gene <- il10_gene[match(sample_order, il10_gene$sample), ]
il10_res2 <- il10_res[il10_res$gene == "IL10" & il10_res$sample %in% sample_order, ]
il10_res2 <- il10_res2[match(sample_order, il10_res2$sample), ]

if (nrow(il10_res2) != 4) stop("IL10 resolution does not contain all four samples.")

same_il10 <- all(
  il10_gene$canonical_cells == il10_res2$canonical_cells &
    il10_gene$MEX_UMI == il10_res2$MEX_UMI &
    il10_gene$MEX_positive_cells == il10_res2$MEX_positive_cells
)

if (!same_il10) {
  stop("IL10 resolution and gene reconstruction disagree in validated leading fields.")
}

cat("Candidate genes:", length(candidate_order), "/ 15 validated\n")
cat("Sample-gene rows:", nrow(gene_rec), "/ 60 validated\n")
cat("IL10 cross-file calibration: PASS\n\n")

# ------------------------------------------------------------------------------
# 3. Build sample-level metrics
# ------------------------------------------------------------------------------

strict15 <- strict[match(candidate_order, strict$gene), ]

# Per-sample canonical cell counts from molecule reconstruction.
sample_cells <- aggregate(
  canonical_cells ~ sample,
  gene_rec,
  FUN = function(x) unique(x)[1]
)
sample_cells <- sample_cells[match(sample_order, sample_cells$sample), ]

if (anyNA(sample_cells$canonical_cells)) stop("Missing sample cell counts.")

cell_map <- setNames(sample_cells$canonical_cells, sample_cells$sample)

# Long table from standard MEX evidence
mex_long <- gene_rec
mex_long$sample <- factor(mex_long$sample, levels = sample_order)
mex_long$MEX_UMI_per_1000_cells <- 1000 * mex_long$MEX_UMI / mex_long$canonical_cells
mex_long$MEX_positive_fraction <- mex_long$MEX_positive_cells / mex_long$canonical_cells
mex_long$group <- ifelse(
  mex_long$sample %in% c("SfLB1", "SfLB2"),
  "PreTransplant_SfLB",
  "Engrafted_Tx"
)

# Long table from strict shared-TRUE cross-sample sensitivity
strict_long <- do.call(
  rbind,
  lapply(seq_len(nrow(strict15)), function(i) {
    g <- strict15$gene[i]

    umi <- c(
      SfLB1 = strict15$SfLB1_shared_TRUE_UMI[i],
      SfLB2 = strict15$SfLB2_shared_TRUE_UMI[i],
      Tx17  = strict15$Tx17_shared_TRUE_UMI[i],
      Tx5   = strict15$Tx5_shared_TRUE_UMI[i]
    )

    data.frame(
      gene = g,
      sample = sample_order,
      canonical_cells = as.numeric(cell_map[sample_order]),
      shared_TRUE_UMI = as.numeric(umi[sample_order]),
      shared_TRUE_UMI_per_1000_cells =
        1000 * as.numeric(umi[sample_order]) / as.numeric(cell_map[sample_order]),
      stringsAsFactors = FALSE
    )
  })
)

strict_long$sample <- factor(strict_long$sample, levels = sample_order)
strict_long$group <- ifelse(
  strict_long$sample %in% c("SfLB1", "SfLB2"),
  "PreTransplant_SfLB",
  "Engrafted_Tx"
)

write.csv(
  mex_long,
  file.path(TAB_DIR, "01_standard_MEX_sample_metrics_v1.4.4A.2.csv"),
  row.names = FALSE
)

write.csv(
  strict_long,
  file.path(TAB_DIR, "02_strict_shared_TRUE_sample_metrics_v1.4.4A.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Replicate-aware transition classification
# ------------------------------------------------------------------------------

classify_gene <- function(g) {
  z <- strict_long[strict_long$gene == g, ]

  sr <- setNames(
    z$shared_TRUE_UMI_per_1000_cells,
    as.character(z$sample)
  )

  sflb <- sr[c("SfLB1", "SfLB2")]
  tx   <- sr[c("Tx17", "Tx5")]

  row_s <- strict15[strict15$gene == g, ]

  preexisting_both <- all(sflb > 0)
  tx_detected_both  <- all(tx > 0)

  complete_higher <- min(tx) > max(sflb)
  complete_lower  <- max(tx) < min(sflb)

  log2fc <- row_s$strict_Tx_vs_SfLB_log2FC_per_cell_descriptive

  replicate_relation <- if (complete_higher) {
    "Tx_higher_complete_separation"
  } else if (complete_lower) {
    "Tx_lower_complete_separation"
  } else {
    "Replicate_ranges_overlap"
  }

  if (g == "IL10") {
    transition_call <- paste0(
      "Engineered_IL10_anchor; Tx-associated MEX signal in both Tx samples; ",
      "strict shared-probe sensitivity heterogeneous"
    )
  } else if (preexisting_both && complete_higher) {
    transition_call <- "Pre-existing_in_SfLB; Tx_higher_complete_replicate_range_separation"
  } else if (preexisting_both && complete_lower) {
    transition_call <- "Pre-existing_in_SfLB; Tx_lower_complete_replicate_range_separation"
  } else if (preexisting_both && log2fc > 0) {
    transition_call <- "Pre-existing_in_SfLB; Tx_mean_higher_but_replicate_ranges_overlap"
  } else if (preexisting_both && log2fc < 0) {
    transition_call <- "Pre-existing_in_SfLB; no_clear_Tx_enhancement_SfLB_mean_higher"
  } else {
    transition_call <- "Other_or_low_detection_pattern"
  }

  data.frame(
    gene = g,
    SfLB1_shared_TRUE_UMI_per_1000 = sflb["SfLB1"],
    SfLB2_shared_TRUE_UMI_per_1000 = sflb["SfLB2"],
    Tx17_shared_TRUE_UMI_per_1000 = tx["Tx17"],
    Tx5_shared_TRUE_UMI_per_1000 = tx["Tx5"],
    SfLB_mean_shared_TRUE_UMI_per_cell =
      row_s$strict_SfLB_mean_UMI_per_cell,
    Tx_mean_shared_TRUE_UMI_per_cell =
      row_s$strict_Tx_mean_UMI_per_cell,
    Tx_vs_SfLB_log2FC_per_cell_descriptive = log2fc,
    preexisting_in_both_SfLB_strict = preexisting_both,
    detected_in_both_Tx_strict = tx_detected_both,
    replicate_relation = replicate_relation,
    transition_call = transition_call,
    stringsAsFactors = FALSE
  )
}

classification <- do.call(
  rbind,
  lapply(candidate_order, classify_gene)
)

# Add standard MEX IL10 details.
il10_m <- mex_long[mex_long$gene == "IL10", ]
il10_m <- il10_m[match(sample_order, as.character(il10_m$sample)), ]

classification$IL10_special_boundary <- ""
classification$IL10_special_boundary[classification$gene == "IL10"] <-
  paste0(
    "SfLB1 contains rare baseline MEX signal (2 UMI / 2 cells) and SfLB2 is zero; ",
    "Tx17 has 18 UMI / 13 cells and Tx5 has 4 UMI / 4 cells. ",
    "Strict shared-TRUE sensitivity is 0/0/5/0 UMI for SfLB1/SfLB2/Tx17/Tx5. ",
    "Therefore IL10 is treated as the engineered anchor by study design plus Tx-associated ",
    "molecule evidence, but individual IL10 molecules are not labeled transgene-specific."
  )

write.csv(
  classification,
  file.path(TAB_DIR, "03_ligand_transition_classification_v1.4.4A.2.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Compact biological summary
# ------------------------------------------------------------------------------

non_il10 <- classification[classification$gene != "IL10", ]

preexisting_genes <- non_il10$gene[non_il10$preexisting_in_both_SfLB_strict]

complete_up <- non_il10$gene[
  non_il10$replicate_relation == "Tx_higher_complete_separation"
]

complete_down <- non_il10$gene[
  non_il10$replicate_relation == "Tx_lower_complete_separation"
]

overlap_tx_higher <- non_il10$gene[
  non_il10$replicate_relation == "Replicate_ranges_overlap" &
    non_il10$Tx_vs_SfLB_log2FC_per_cell_descriptive > 0
]

overlap_sflb_higher <- non_il10$gene[
  non_il10$replicate_relation == "Replicate_ranges_overlap" &
    non_il10$Tx_vs_SfLB_log2FC_per_cell_descriptive < 0
]

summary_tbl <- data.frame(
  item = c(
    "Non-IL10 candidate ligands pre-existing in both SfLB samples under strict shared-TRUE sensitivity",
    "Pre-existing ligands with Tx-higher complete replicate-range separation",
    "Pre-existing ligands with Tx-mean higher but overlapping replicate ranges",
    "Pre-existing ligands without clear Tx enhancement (SfLB mean higher)",
    "Pre-existing ligands with Tx-lower complete replicate-range separation",
    "IL10 engineered anchor"
  ),
  n = c(
    length(preexisting_genes),
    length(complete_up),
    length(overlap_tx_higher),
    length(overlap_sflb_higher),
    length(complete_down),
    1
  ),
  genes = c(
    paste(preexisting_genes, collapse = ";"),
    paste(complete_up, collapse = ";"),
    paste(overlap_tx_higher, collapse = ";"),
    paste(overlap_sflb_higher, collapse = ";"),
    paste(complete_down, collapse = ";"),
    "IL10"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_tbl,
  file.path(TAB_DIR, "04_ligand_transition_key_summary_v1.4.4A.2.csv"),
  row.names = FALSE
)

cat("=== KEY SUMMARY ===\n")
print(summary_tbl, row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. IL10 anchor table
# ------------------------------------------------------------------------------

il10_table <- mex_long[mex_long$gene == "IL10", c(
  "sample",
  "canonical_cells",
  "MEX_UMI",
  "MEX_positive_cells",
  "MEX_UMI_per_1000_cells",
  "MEX_positive_fraction"
)]

il10_strict <- strict_long[strict_long$gene == "IL10", c(
  "sample",
  "shared_TRUE_UMI",
  "shared_TRUE_UMI_per_1000_cells"
)]

il10_table <- merge(
  il10_table,
  il10_strict,
  by = "sample",
  all = TRUE,
  sort = FALSE
)
il10_table$sample <- factor(il10_table$sample, levels = sample_order)
il10_table <- il10_table[order(il10_table$sample), ]

write.csv(
  il10_table,
  file.path(TAB_DIR, "05_IL10_engineered_anchor_resolution_v1.4.4A.2.csv"),
  row.names = FALSE
)

cat("\n=== IL10 ENGINEERED ANCHOR ===\n")
print(il10_table, row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. Figures
# ------------------------------------------------------------------------------

# ---- Fig 1. Strict shared-TRUE sample heatmap
heat <- strict_long
gene_order_plot <- classification$gene[
  order(
    classification$gene == "IL10",
    classification$Tx_vs_SfLB_log2FC_per_cell_descriptive,
    decreasing = TRUE
  )
]
# Keep IL10 first.
gene_order_plot <- c(
  "IL10",
  setdiff(
    classification$gene[
      order(classification$Tx_vs_SfLB_log2FC_per_cell_descriptive, decreasing = TRUE)
    ],
    "IL10"
  )
)

heat$gene <- factor(heat$gene, levels = rev(gene_order_plot))
heat$sample <- factor(heat$sample, levels = sample_order)
heat$plot_value <- log10(1 + heat$shared_TRUE_UMI_per_1000_cells)

p1 <- ggplot(
  heat,
  aes(x = sample, y = gene, fill = plot_value)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%.2f", shared_TRUE_UMI_per_1000_cells)),
    size = 2.8
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = median(heat$plot_value, na.rm = TRUE)
  ) +
  labs(
    title = "Probe-aware graft ligand abundance",
    subtitle = "Strict shared-TRUE probe sensitivity; values = UMI per 1,000 cells",
    x = NULL,
    y = NULL,
    fill = "log10(1 + UMI/1000)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG_DIR, "Fig01_strict_shared_TRUE_ligand_heatmap_v1.4.4A.2.pdf"),
  p1, width = 8.6, height = 8.2, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig01_strict_shared_TRUE_ligand_heatmap_v1.4.4A.2.png"),
  p1, width = 8.6, height = 8.2, dpi = 300
)

# ---- Fig 2. Replicate-aware strict shared-TRUE abundance, non-IL10
rp <- strict_long[strict_long$gene != "IL10", ]
rp$display_group <- factor(
  ifelse(rp$sample %in% c("SfLB1", "SfLB2"), "Pre-transplant SfLB", "Engrafted Tx"),
  levels = c("Pre-transplant SfLB", "Engrafted Tx")
)
rp$gene <- factor(
  rp$gene,
  levels = rev(
    classification$gene[
      classification$gene != "IL10"
    ][order(
      classification$Tx_vs_SfLB_log2FC_per_cell_descriptive[
        classification$gene != "IL10"
      ],
      decreasing = TRUE
    )]
  )
)

p2 <- ggplot(
  rp,
  aes(x = display_group, y = shared_TRUE_UMI_per_1000_cells)
) +
  geom_point(
    aes(shape = sample),
    size = 2.5,
    position = position_jitter(width = 0.05, height = 0)
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.38,
    linewidth = 0.5
  ) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  labs(
    title = "Replicate-aware ligand transition",
    subtitle = "Strict shared-TRUE UMI per 1,000 cells; descriptive n=2 vs n=2",
    x = NULL,
    y = "UMI per 1,000 cells",
    shape = "Sample"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG_DIR, "Fig02_replicate_aware_nonIL10_ligand_transition_v1.4.4A.2.pdf"),
  p2, width = 12.5, height = 10.0, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig02_replicate_aware_nonIL10_ligand_transition_v1.4.4A.2.png"),
  p2, width = 12.5, height = 10.0, dpi = 300
)

# ---- Fig 3. Descriptive log2FC
fc <- classification[classification$gene != "IL10", ]
fc$gene <- factor(
  fc$gene,
  levels = fc$gene[order(fc$Tx_vs_SfLB_log2FC_per_cell_descriptive)]
)

p3 <- ggplot(
  fc,
  aes(x = Tx_vs_SfLB_log2FC_per_cell_descriptive, y = gene)
) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_segment(
    aes(x = 0, xend = Tx_vs_SfLB_log2FC_per_cell_descriptive, yend = gene),
    linewidth = 0.7
  ) +
  geom_point(size = 3) +
  labs(
    title = "Descriptive SfLB-to-Tx ligand shift",
    subtitle = "Cross-sample strict shared-TRUE sensitivity; no inferential test",
    x = "Tx vs SfLB log2 fold-change per cell",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  file.path(FIG_DIR, "Fig03_strict_shared_TRUE_descriptive_log2FC_v1.4.4A.2.pdf"),
  p3, width = 8.4, height = 6.8, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig03_strict_shared_TRUE_descriptive_log2FC_v1.4.4A.2.png"),
  p3, width = 8.4, height = 6.8, dpi = 300
)

# ---- Fig 4. IL10 standard-MEX anchor
ilplot <- il10_table
ilplot$sample <- factor(ilplot$sample, levels = sample_order)

p4 <- ggplot(ilplot, aes(x = sample, y = MEX_UMI_per_1000_cells)) +
  geom_col(width = 0.68) +
  geom_text(
    aes(label = paste0(MEX_UMI, " UMI\n", MEX_positive_cells, " cells")),
    vjust = -0.25,
    size = 3.2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "IL10 engineered-anchor resolution",
    subtitle = "Standard filtered gene-expression matrix reproduced at molecule/probe level",
    x = NULL,
    y = "IL10 UMI per 1,000 cells"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG_DIR, "Fig04_IL10_MEX_engineered_anchor_v1.4.4A.2.pdf"),
  p4, width = 7.2, height = 5.6, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig04_IL10_MEX_engineered_anchor_v1.4.4A.2.png"),
  p4, width = 7.2, height = 5.6, dpi = 300
)

# ------------------------------------------------------------------------------
# 8. Interpretation note
# ------------------------------------------------------------------------------

note <- c(
  "Human MASH v1.4.4A.2 interpretation notes",
  "",
  "Biological comparison:",
  "- SfLB1/SfLB2 = pre-transplant, non-engineered cultured organoids.",
  "- Tx17/Tx5 = engrafted human cells derived from IL10-engineered organoids.",
  "",
  "Primary evidence:",
  "- Cross-sample ligand comparison uses B.2C.3 strict shared-TRUE probe sensitivity.",
  "- IL10 is resolved separately with standard MEX molecule/cell evidence because the",
  "  all-four shared-TRUE sensitivity is intentionally stringent and is zero in Tx5.",
  "",
  "Critical causal boundary:",
  "- SfLB -> Tx is NOT an IL10-only perturbation.",
  "- The contrast combines transgene status, transplantation/engraftment, host microenvironment,",
  "  post-engraftment selection/state transition, and possible residual innate xenogeneic immune pressure.",
  "- v1.4.4A.1.1 demonstrated major lineage-composition shifts.",
  "- Therefore whole-graft ligand changes do NOT establish same-lineage induction.",
  "- Source-lineage decomposition is required in v1.4.4A.3.",
  "",
  "IL10 boundary:",
  "- SfLB1 contains rare baseline MEX signal; SfLB2 is zero.",
  "- Both Tx samples contain standard-MEX IL10 signal reproduced at molecule/probe level.",
  "- IL10 is the engineered anchor by study design plus Tx-associated evidence.",
  "- Individual IL10 molecules are not called transgene-specific from gene-level/probe evidence alone.",
  "",
  "Classification rule:",
  "- 'Complete separation higher' means the lower of the two Tx per-cell values exceeds",
  "  the higher of the two SfLB per-cell values.",
  "- 'Complete separation lower' is the converse.",
  "- Overlapping ranges are retained as descriptive and are not over-interpreted.",
  "- No inferential test is performed (n=2 biological samples per group)."
)

writeLines(
  note,
  file.path(OUT_ROOT, "README_interpretation_v1.4.4A.2.txt"),
  useBytes = TRUE
)

cat("\n============================================================\n")
cat("v1.4.4A.2 COMPLETE\n")
cat("============================================================\n")
cat("Output:", OUT_ROOT, "\n")
cat("Tables:", length(list.files(TAB_DIR)), "\n")
cat("Figures:", length(list.files(FIG_DIR)), "\n")
cat("No inferential test performed (n=2 per group).\n")
cat("Whole-graft changes require source-lineage decomposition in v1.4.4A.3.\n")
