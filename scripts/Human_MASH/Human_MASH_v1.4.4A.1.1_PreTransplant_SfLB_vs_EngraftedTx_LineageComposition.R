#!/usr/bin/env Rscript

# ==============================================================================
# Human_MASH_v1.4.4A.1.1_PreTransplant_SfLB_vs_EngraftedTx_LineageComposition.R
#
# Purpose
#   Descriptive, replicate-aware composition audit comparing:
#     SfLB1/SfLB2 = pre-transplant, non-engineered cultured organoids
#     Tx17/Tx5    = engrafted human cells derived from IL10-engineered organoids
#
# Important causal boundary
#   SfLB -> Tx is NOT interpreted as an IL10-only perturbation.
#   The contrast combines transgene status, transplantation/engraftment,
#   host-environment exposure, and post-engraftment cell selection/state change.
#
# Input frozen upstream
#   Human_MASH_iPSC_derived_v1.4.2A.1_hepatocyte_state_interpreted.rds
#
# Outputs
#   Tables + descriptive figures only. No re-clustering and no new annotation.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

IN_RDS <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2A.1/objects/Human_MASH_iPSC_derived_v1.4.2A.1_hepatocyte_state_interpreted.rds"

OUT_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.1.1/analysis/v1.4.4A.1.1_pretransplant_SfLB_vs_engraftedTx_lineage_composition"
TAB_DIR  <- file.path(OUT_ROOT, "tables")
FIG_DIR  <- file.path(OUT_ROOT, "figures")
LOG_DIR  <- file.path(OUT_ROOT, "logs")

dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, "v1.4.4A.1.1_run_log.txt")
sink(LOG_FILE, split = TRUE)
on.exit(sink(), add = TRUE)

cat("============================================================\n")
cat("Human MASH v1.4.4A.1\n")
cat("Pre-transplant SfLB vs engrafted Tx lineage composition audit\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Input:", IN_RDS, "\n")
cat("Output:", OUT_ROOT, "\n\n")

stopifnot(file.exists(IN_RDS))
obj <- readRDS(IN_RDS)

cat("Object class:", paste(class(obj), collapse = ", "), "\n")
cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")
cat("Assays:", paste(names(obj@assays), collapse = ", "), "\n")
cat("Reductions:", paste(names(obj@reductions), collapse = ", "), "\n\n")

required_md <- c(
  "sample_v141",
  "broad_lineage_manual_v142A",
  "manual_state_v142A",
  "annotation_status_v142A"
)

missing_md <- setdiff(required_md, colnames(obj[[]]))
if (length(missing_md) > 0) {
  stop("Missing required metadata columns: ", paste(missing_md, collapse = ", "))
}

md <- obj[[]]
md$cell <- rownames(md)

sample_order <- c("SfLB1", "SfLB2", "Tx17", "Tx5")
sample_values <- sort(unique(as.character(md$sample_v141)))

cat("Observed samples:", paste(sample_values, collapse = ", "), "\n")
if (!setequal(sample_values, sample_order)) {
  stop(
    "Unexpected sample set. Observed: ",
    paste(sample_values, collapse = ", "),
    " ; expected: ",
    paste(sample_order, collapse = ", ")
  )
}

if (nrow(md) != 43415L) {
  stop("Unexpected cell count: ", nrow(md), " (expected 43415)")
}

status_values <- unique(as.character(md$annotation_status_v142A))
cat("Annotation status:", paste(status_values, collapse = " | "), "\n")
if (!all(grepl("HEPATOCYTE_STATE_INTERPRETATION_FROZEN", status_values))) {
  warning("Unexpected annotation_status_v142A value.")
}

md$sample <- factor(as.character(md$sample_v141), levels = sample_order)
md$group <- ifelse(
  as.character(md$sample) %in% c("SfLB1", "SfLB2"),
  "PreTransplant_SfLB",
  "Engrafted_Tx"
)
md$group <- factor(md$group, levels = c("PreTransplant_SfLB", "Engrafted_Tx"))

md$lineage <- as.character(md$broad_lineage_manual_v142A)
md$state   <- as.character(md$manual_state_v142A)

if (anyNA(md$lineage)) stop("NA found in broad_lineage_manual_v142A")
if (anyNA(md$state))   stop("NA found in manual_state_v142A")

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
}

tabulate_by_sample <- function(meta, category_col) {
  samples <- levels(meta$sample)
  all_categories <- sort(unique(as.character(meta[[category_col]])))

  out <- do.call(rbind, lapply(samples, function(s) {
    ii <- meta$sample == s

    # Explicit factor levels are essential: categories absent from a sample
    # must be retained as n=0/fraction=0, rather than silently dropped.
    vv <- factor(as.character(meta[[category_col]][ii]), levels = all_categories)
    tt <- table(vv, useNA = "no")

    data.frame(
      sample = s,
      category = names(tt),
      n_cells = as.integer(tt),
      sample_total = sum(ii),
      fraction = as.integer(tt) / sum(ii),
      stringsAsFactors = FALSE
    )
  }))
  out$sample <- factor(out$sample, levels = sample_order)
  out$group <- ifelse(
    out$sample %in% c("SfLB1", "SfLB2"),
    "PreTransplant_SfLB",
    "Engrafted_Tx"
  )
  out$group <- factor(out$group, levels = c("PreTransplant_SfLB", "Engrafted_Tx"))
  out
}

summarize_group_replicates <- function(tbl) {
  cats <- unique(tbl$category)
  do.call(rbind, lapply(cats, function(cc) {
    do.call(rbind, lapply(levels(tbl$group), function(gg) {
      x <- tbl[tbl$category == cc & tbl$group == gg, "fraction", drop = TRUE]
      data.frame(
        category = cc,
        group = gg,
        n_replicates = length(x),
        mean_fraction = mean(x),
        sd_fraction = if (length(x) > 1) sd(x) else NA_real_,
        min_fraction = min(x),
        max_fraction = max(x),
        replicate_values = paste(sprintf("%.6f", x), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

sample_counts <- as.data.frame(table(md$sample), stringsAsFactors = FALSE)
colnames(sample_counts) <- c("sample", "n_cells")
sample_counts$sample <- factor(sample_counts$sample, levels = sample_order)
sample_counts$group <- ifelse(
  sample_counts$sample %in% c("SfLB1", "SfLB2"),
  "PreTransplant_SfLB",
  "Engrafted_Tx"
)
safe_write_csv(sample_counts, file.path(TAB_DIR, "01_sample_cell_counts_v1.4.4A.1.1.csv"))

cat("\n=== SAMPLE CELL COUNTS ===\n")
print(sample_counts)

lineage_tbl <- tabulate_by_sample(md, "lineage")
expected_lineage_rows <- length(sample_order) * length(unique(md$lineage))
if (nrow(lineage_tbl) != expected_lineage_rows) {
  stop("Broad-lineage zero-completion failed: ", nrow(lineage_tbl),
       " rows; expected ", expected_lineage_rows)
}
lineage_tbl$category <- as.character(lineage_tbl$category)
safe_write_csv(lineage_tbl, file.path(TAB_DIR, "02_broad_lineage_counts_fractions_by_sample_v1.4.4A.1.1.csv"))

lineage_group <- summarize_group_replicates(lineage_tbl)
safe_write_csv(lineage_group, file.path(TAB_DIR, "03_broad_lineage_group_replicate_summary_v1.4.4A.1.1.csv"))

lineages <- sort(unique(lineage_tbl$category))
wide_lineage <- data.frame(lineage = lineages, stringsAsFactors = FALSE)

for (s in sample_order) {
  tmp <- lineage_tbl[lineage_tbl$sample == s, c("category", "fraction")]
  names(tmp)[2] <- s
  wide_lineage <- merge(wide_lineage, tmp, by.x = "lineage", by.y = "category",
                        all.x = TRUE, sort = FALSE)
}
wide_lineage[is.na(wide_lineage)] <- 0
wide_lineage$SfLB_mean <- rowMeans(wide_lineage[, c("SfLB1", "SfLB2"), drop = FALSE])
wide_lineage$Tx_mean <- rowMeans(wide_lineage[, c("Tx17", "Tx5"), drop = FALSE])
wide_lineage$Tx_minus_SfLB <- wide_lineage$Tx_mean - wide_lineage$SfLB_mean

safe_write_csv(wide_lineage, file.path(TAB_DIR, "04_broad_lineage_SfLB_vs_Tx_descriptive_delta_v1.4.4A.1.1.csv"))

cat("\n=== BROAD LINEAGE DESCRIPTIVE SUMMARY ===\n")
print(wide_lineage[order(abs(wide_lineage$Tx_minus_SfLB), decreasing = TRUE), ], row.names = FALSE)

state_tbl <- tabulate_by_sample(md, "state")
expected_state_rows <- length(sample_order) * length(unique(md$state))
if (nrow(state_tbl) != expected_state_rows) {
  stop("State zero-completion failed: ", nrow(state_tbl),
       " rows; expected ", expected_state_rows)
}
safe_write_csv(state_tbl, file.path(TAB_DIR, "05_manual_state_counts_fractions_by_sample_v1.4.4A.1.1.csv"))

state_group <- summarize_group_replicates(state_tbl)
safe_write_csv(state_group, file.path(TAB_DIR, "06_manual_state_group_replicate_summary_v1.4.4A.1.1.csv"))

states <- sort(unique(state_tbl$category))
wide_state <- data.frame(state = states, stringsAsFactors = FALSE)

for (s in sample_order) {
  tmp <- state_tbl[state_tbl$sample == s, c("category", "fraction")]
  names(tmp)[2] <- s
  wide_state <- merge(wide_state, tmp, by.x = "state", by.y = "category",
                      all.x = TRUE, sort = FALSE)
}
wide_state[is.na(wide_state)] <- 0
wide_state$SfLB_mean <- rowMeans(wide_state[, c("SfLB1", "SfLB2"), drop = FALSE])
wide_state$Tx_mean   <- rowMeans(wide_state[, c("Tx17", "Tx5"), drop = FALSE])
wide_state$Tx_minus_SfLB <- wide_state$Tx_mean - wide_state$SfLB_mean

state_bias <- do.call(rbind, lapply(states, function(st) {
  x <- md[md$state == st, , drop = FALSE]
  tt <- table(x$sample)
  frac <- tt / sum(tt)
  data.frame(
    state = st,
    state_n_cells = nrow(x),
    dominant_sample = names(frac)[which.max(frac)],
    dominant_sample_fraction = max(as.numeric(frac)),
    sample_bias_ge70 = max(as.numeric(frac)) >= 0.70,
    stringsAsFactors = FALSE
  )
}))

wide_state <- merge(wide_state, state_bias, by = "state", all.x = TRUE, sort = FALSE)
safe_write_csv(wide_state, file.path(TAB_DIR, "07_manual_state_SfLB_vs_Tx_descriptive_delta_and_bias_v1.4.4A.1.1.csv"))

cat("\n=== MANUAL STATE DESCRIPTIVE SUMMARY ===\n")
print(
  wide_state[order(abs(wide_state$Tx_minus_SfLB), decreasing = TRUE),
             c("state","SfLB1","SfLB2","Tx17","Tx5",
               "SfLB_mean","Tx_mean","Tx_minus_SfLB",
               "dominant_sample","dominant_sample_fraction","sample_bias_ge70")],
  row.names = FALSE
)

lineage_levels <- c(
  "Hepatocyte_like",
  "Hepatobiliary_like",
  "Mesenchymal_like",
  "Endothelial_like",
  "Myeloid_like",
  "Unresolved"
)
lineage_levels <- lineage_levels[lineage_levels %in% unique(lineage_tbl$category)]

lineage_tbl$category <- factor(lineage_tbl$category, levels = lineage_levels)

palette_lineage <- c(
  "Hepatocyte_like"    = "#00A6A6",
  "Hepatobiliary_like" = "#00CC66",
  "Mesenchymal_like"   = "#FF4FA3",
  "Endothelial_like"   = "#00BFFF",
  "Myeloid_like"       = "#E74C3C",
  "Unresolved"         = "#8C8C8C"
)

p1 <- ggplot(lineage_tbl, aes(x = sample, y = fraction, fill = category)) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = palette_lineage, drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), expand = c(0, 0)) +
  labs(
    title = "Broad-lineage composition: pre-transplant SfLB vs engrafted Tx",
    subtitle = "Descriptive cell fractions; SfLB -> Tx is not an IL10-only perturbation",
    x = NULL, y = "Fraction of cells", fill = "Frozen broad lineage"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold"),
        legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "Fig01_broad_lineage_composition_by_sample_v1.4.4A.1.1.pdf"),
       p1, width = 10.5, height = 6.2, useDingbats = FALSE)
ggsave(file.path(FIG_DIR, "Fig01_broad_lineage_composition_by_sample_v1.4.4A.1.1.png"),
       p1, width = 10.5, height = 6.2, dpi = 300)

plot_group <- lineage_tbl
plot_group$group <- factor(
  ifelse(plot_group$sample %in% c("SfLB1","SfLB2"),
         "Pre-transplant SfLB", "Engrafted Tx"),
  levels = c("Pre-transplant SfLB", "Engrafted Tx")
)

p2 <- ggplot(plot_group, aes(x = group, y = fraction)) +
  geom_point(aes(shape = sample), size = 2.7,
             position = position_jitter(width = 0.06, height = 0)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.42,
               linewidth = 0.45, middle.linewidth = 0.85) +
  facet_wrap(~ category, scales = "free_y", ncol = 3) +
  scale_y_continuous(labels = function(x) sprintf("%.1f%%", 100*x)) +
  labs(
    title = "Replicate-aware broad-lineage composition",
    subtitle = "Two biological samples per group; descriptive means only",
    x = NULL, y = "Cell fraction", shape = "Sample"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "Fig02_broad_lineage_replicate_aware_SfLB_vs_Tx_v1.4.4A.1.1.pdf"),
       p2, width = 11.5, height = 8.0, useDingbats = FALSE)
ggsave(file.path(FIG_DIR, "Fig02_broad_lineage_replicate_aware_SfLB_vs_Tx_v1.4.4A.1.1.png"),
       p2, width = 11.5, height = 8.0, dpi = 300)

state_heat <- state_tbl
state_heat$sample <- factor(state_heat$sample, levels = sample_order)
state_order <- wide_state$state[order(wide_state$Tx_minus_SfLB, decreasing = TRUE)]
state_heat$category <- factor(state_heat$category, levels = rev(state_order))

p3 <- ggplot(state_heat, aes(x = sample, y = category, fill = fraction)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f", 100*fraction)), size = 2.8) +
  scale_fill_gradient(low = "white", high = "#FF1A1A") +
  labs(
    title = "Frozen iPSC-state composition",
    subtitle = "Cell fraction within each sample; values are percentages",
    x = NULL, y = NULL, fill = "Fraction"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 8.5),
        plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "Fig03_manual_state_fraction_heatmap_v1.4.4A.1.1.pdf"),
       p3, width = 8.5, height = 9.5, useDingbats = FALSE)
ggsave(file.path(FIG_DIR, "Fig03_manual_state_fraction_heatmap_v1.4.4A.1.1.png"),
       p3, width = 8.5, height = 9.5, dpi = 300)

if ("umap_iPSC_v141" %in% names(obj@reductions)) {
  emb <- as.data.frame(Embeddings(obj, "umap_iPSC_v141"))
  emb$cell <- rownames(emb)
  emb$sample <- md[emb$cell, "sample", drop = TRUE]
  colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")

  p4 <- ggplot(emb, aes(UMAP_1, UMAP_2, color = sample)) +
    geom_point(size = 0.18, alpha = 0.75) +
    coord_equal() +
    labs(
      title = "Existing iPSC-only RNA UMAP",
      subtitle = "No re-clustering; frozen v1.4.1/v1.4.2 embedding",
      x = "UMAP 1", y = "UMAP 2", color = "Sample"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")

  ggsave(file.path(FIG_DIR, "Fig04_existing_iPSC_UMAP_by_sample_v1.4.4A.1.1.pdf"),
         p4, width = 8.2, height = 6.7, useDingbats = FALSE)
  ggsave(file.path(FIG_DIR, "Fig04_existing_iPSC_UMAP_by_sample_v1.4.4A.1.1.png"),
         p4, width = 8.2, height = 6.7, dpi = 300)
}

note <- c(
  "Human MASH v1.4.4A.1.1 analysis notes",
  "",
  "Biological definitions:",
  "- SfLB1/SfLB2 = pre-transplant, non-engineered cultured iPSC liver organoids.",
  "- Tx17/Tx5 = engrafted human graft cells derived from IL10-engineered organoids.",
  "",
  "Interpretation boundary:",
  "- SfLB -> Tx is not an IL10-only perturbation.",
  "- The contrast includes transgene status, transplantation/engraftment, host environment, and post-engraftment selection/state changes.",
  "- Composition results are descriptive; there are two biological samples per group.",
  "- Frozen v1.4.2A.1 annotations are reused without re-clustering.",
  "",
  "Next planned modules:",
  "- v1.4.4A.2: B.2C.3 probe-aware ligand comparison.",
  "- v1.4.4A.3: ligand source-lineage atlas.",
  "- v1.4.4A.4/A.5: measurable hepatocyte-function panel and LT7/LT8 reference comparison."
)

writeLines(note, file.path(OUT_ROOT, "README_interpretation_v1.4.4A.1.1.txt"), useBytes = TRUE)

cat("\n============================================================\n")
cat("v1.4.4A.1.1 COMPLETE\n")
cat("============================================================\n")
cat("Output:", OUT_ROOT, "\n")
cat("Tables:", length(list.files(TAB_DIR)), "\n")
cat("Figures:", length(list.files(FIG_DIR)), "\n")
cat("No re-clustering performed.\n")
cat("Absent sample-category combinations explicitly encoded as zero.\n")
cat("No inferential composition test performed (n=2 per group).\n")
