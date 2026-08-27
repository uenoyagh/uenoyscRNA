#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
set.seed(5900)

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

# Mouse MASH macrophage staining-marker integrated summary
# Version: v5.9.1
#
# No reclustering / reintegration / re-UMAP.
# Whole-cell: frozen parent v5.1.1
# Macrophage: Clean-B FINAL v4.14.5
#
# v5.9.1 additions
#   - Explicit Sham-vs-Tx comparison block
#   - Sham-vs-Tx Violin plots for selected and reference markers
#   - Sham1/Sham20/Tx17/Tx5 focused allMphi_positive_fraction_by_sample
#   - Sample-level positive-fraction and mean-expression heatmaps
#   - Sham-vs-Tx two-column positive-fraction / mean-expression heatmaps
#   - Sham-vs-Tx delta table
#
# Selected priority markers:
# Anti-inflammatory: Cd163 > Clec4f >> Slc40a1
# Inflammatory: Spn >> Ccr2
# ECM-associated inflammatory: F13a1 > Tgfb1 > Vcan >>> Pstpip1 > Rasgrp4 > Thbs1
# Repair/Resolution: Mmp12 > Fabp7 >>> Dabp7 > Ksr2
# Lipid-associated/TREM2: Gpnmb > Trem2 >>> Fabp5 > Atp6v0d2
#
# Dabp7 is retained exactly as requested. If absent, it is logged and skipped.

msg <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(...))
}

save_pdf <- function(p, file, width, height) {
  pdf(file, width = width, height = height, useDingbats = FALSE)
  print(p)
  dev.off()
}

find_first_existing <- function(x) {
  y <- x[file.exists(x)]
  if (!length(y)) stop("Input not found:\n", paste(x, collapse = "\n"))
  y[[1]]
}

resolve_first <- function(x, candidates, what) {
  y <- candidates[candidates %in% x]
  if (!length(y)) stop("Could not resolve ", what, ": ", paste(candidates, collapse = ", "))
  y[[1]]
}

canonical_condition <- function(x) {
  x <- as.character(x)
  case_when(
    grepl("^STD", x, ignore.case = TRUE) ~ "STD",
    grepl("CDHFD|CDAHFD", x, ignore.case = TRUE) ~ "CDAHFD",
    grepl("^Sham", x, ignore.case = TRUE) ~ "Sham",
    grepl("^Tx", x, ignore.case = TRUE) ~ "Tx",
    TRUE ~ x
  )
}

resolve_genes <- function(requested, available) {
  al <- tolower(available)
  rr <- vapply(requested, function(g) {
    i <- match(tolower(g), al)
    if (is.na(i)) NA_character_ else available[[i]]
  }, character(1))
  tibble(requested_gene = requested, resolved_gene = unname(rr), present = !is.na(rr))
}

row_zscore <- function(mat) {
  z <- t(apply(mat, 1, function(x) {
    sx <- sd(x, na.rm = TRUE)
    if (!is.finite(sx) || sx == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / sx
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

feature_pages <- function(obj, genes, reduction, file, title, pt = 0.2, ncol = 3, per_page = 12) {
  if (!length(genes)) return(invisible(NULL))
  pdf(file, width = 13, height = 4.1 * ceiling(per_page / ncol), useDingbats = FALSE)
  chunks <- split(seq_along(genes), ceiling(seq_along(genes) / per_page))
  for (i in seq_along(chunks)) {
    gs <- genes[chunks[[i]]]
    p <- FeaturePlot(
      obj, features = gs, reduction = reduction, ncol = ncol,
      order = TRUE, min.cutoff = "q05", max.cutoff = "q98",
      cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
      raster = FALSE, pt.size = pt
    ) & theme_classic(base_size = 8)
    p <- p + plot_annotation(title = paste0(title, " | page ", i, "/", length(chunks)))
    print(p)
  }
  dev.off()
}

split_feature_pdf <- function(obj, genes, reduction, split_by, file, title, pt = 0.2) {
  if (!length(genes)) return(invisible(NULL))
  pdf(file, width = 11, height = 5.2, useDingbats = FALSE)
  for (g in genes) {
    p <- FeaturePlot(
      obj, features = g, reduction = reduction, split.by = split_by,
      keep.scale = "all", ncol = 2, order = TRUE,
      min.cutoff = "q05", max.cutoff = "q98",
      cols = c("#0033FF", "#FFFFFF", "#FF1A1A"),
      raster = FALSE, pt.size = pt
    ) & theme_classic(base_size = 8)
    p <- p + plot_annotation(title = paste0(title, " | ", g))
    print(p)
  }
  dev.off()
}

summarize_by_group <- function(obj, genes, groups, group_type, assay = "RNA") {
  counts <- GetAssayData(obj, assay = assay, layer = "counts")
  data <- GetAssayData(obj, assay = assay, layer = "data")
  lev <- unique(as.character(groups))
  lev <- lev[!is.na(lev)]
  ans <- list()
  k <- 1L
  for (grp in lev) {
    idx <- which(as.character(groups) == grp)
    for (g in genes) {
      xc <- counts[g, idx, drop = FALSE]
      xd <- data[g, idx, drop = FALSE]
      n <- length(idx)
      np <- Matrix::nnzero(xc)
      ans[[k]] <- tibble(
        gene = g, group = grp, group_type = group_type,
        n_cells = n, n_positive = np,
        positive_fraction = np / n,
        positive_percent = 100 * np / n,
        mean_expression = Matrix::sum(xd) / n
      )
      k <- k + 1L
    }
  }
  bind_rows(ans)
}

heatmap_z <- function(summary_df, group_col, value_col, genes, group_levels, title, file, width, height) {
  w <- summary_df %>%
    select(gene, all_of(group_col), all_of(value_col)) %>%
    pivot_wider(names_from = all_of(group_col), values_from = all_of(value_col))
  mat <- as.matrix(w[, setdiff(colnames(w), "gene"), drop = FALSE])
  rownames(mat) <- w$gene
  z <- row_zscore(mat)
  long <- as.data.frame(z) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "group", values_to = "zscore")
  long$gene <- factor(long$gene, levels = genes)
  long$group <- factor(long$group, levels = group_levels)
  p <- ggplot(long, aes(gene, group, fill = zscore)) +
    geom_tile(linewidth = 0.25) +
    scale_fill_gradient2(low = "#0033FF", mid = "#FFFFFF", high = "#FF1A1A", midpoint = 0) +
    labs(title = title, x = NULL, y = NULL, fill = "z-score") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  save_pdf(p, file, width, height)
  long
}

# Paths
ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

WHOLE_RDS <- file.path(
  ROOT, "Mouse_MASH_RDS", "WholeCell_Layer1_ParentFreeze_v5.1.1", "RDS",
  "Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds"
)

MPHI_RDS <- find_first_existing(c(
  file.path(ROOT, "Mouse_MASH_Mphi_RDS", "Mphi_Res2_CleanB_FINAL_v4.14.5", "RDS",
            "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"),
  file.path(ROOT, "Mouse_MASH_Mphi_RDS", "Mphi_Res2_CleanB_FINAL_v4.14.5", "RDS",
            "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.RDS")
))

OUT <- file.path(ROOT, "Mouse_MASH_Mphi_RDS", "Mphi_StainingMarker_IntegratedSummary_v5.9.1")
FIG <- file.path(OUT, "Figures")
TAB <- file.path(OUT, "Tables")
LOG <- file.path(OUT, "Logs")
for (d in c(OUT, FIG, TAB, LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

ASSAY <- "RNA"
WHOLE_ANNOT <- "wholecell_layer1_FINAL_v511"
MPHI_CLASS <- "macrophage_class_Res2_FINAL_v4145_char"
WHOLE_UMAP <- "umapRPCA"
MPHI_UMAP <- "umapRPCA"

SELECTED <- list(
  "Anti-inflammatory-Mphi" = c("Cd163", "Clec4f", "Slc40a1"),
  "Inflammatory-Mphi" = c("Spn", "Ccr2"),
  "ECM-associated inflammatory-Mphi" = c("F13a1", "Tgfb1", "Vcan", "Pstpip1", "Rasgrp4", "Thbs1"),
  "Repair/Resolution-Mphi" = c("Mmp12", "Fabp7", "Dabp7", "Ksr2"),
  "Lipid-associated/TREM2-Mphi" = c("Gpnmb", "Trem2", "Fabp5", "Atp6v0d2")
)

REFERENCE <- list(
  "Pan-macrophage" = c("Cd68", "Adgre1"),
  "Conventional-M1-M2" = c("Cd80", "Cd86", "Mrc1"),
  "Inflammatory-reference" = c("Il1b", "Tnf", "Nfkbia", "Ccl3", "Ccl4", "Stat1", "Cxcl2", "Ptgs2"),
  "ECM-fibrotic-reference" = c("Fn1", "Spp1", "Mmp14", "Lgals3"),
  "Repair-resolution-reference" = c("Mertk", "Mfge8", "Gas6", "Igf1", "Hmox1", "Il10ra", "Tgm2"),
  "Lipid-TREM2-reference" = c("Cd9", "Lpl", "Apoc1"),
  "Resident-homeostatic-reference" = c("Timd4", "Folr2", "Vsig4", "Marco")
)

SEL_ORDER <- unique(unlist(SELECTED, use.names = FALSE))
REF_ORDER <- unique(unlist(REFERENCE, use.names = FALSE))
ALL_REQ <- unique(c(SEL_ORDER, REF_ORDER))

manifest <- bind_rows(
  lapply(names(SELECTED), function(g) tibble(
    marker_set = "Selected", marker_group = g,
    priority_order = seq_along(SELECTED[[g]]), requested_gene = SELECTED[[g]]
  )),
  lapply(names(REFERENCE), function(g) tibble(
    marker_set = "Reference", marker_group = g,
    priority_order = seq_along(REFERENCE[[g]]), requested_gene = REFERENCE[[g]]
  ))
)
write.csv(manifest, file.path(TAB, "00_marker_manifest_v5.9.1.csv"), row.names = FALSE)

# Load
if (!file.exists(WHOLE_RDS)) stop("Whole-cell RDS not found: ", WHOLE_RDS)
msg("Loading whole-cell RDS...")
whole <- readRDS(WHOLE_RDS)
msg("Loading Clean-B macrophage RDS...")
mphi <- readRDS(MPHI_RDS)

DefaultAssay(whole) <- ASSAY
DefaultAssay(mphi) <- ASSAY

stopifnot(WHOLE_ANNOT %in% colnames(whole@meta.data))
stopifnot(MPHI_CLASS %in% colnames(mphi@meta.data))
stopifnot(WHOLE_UMAP %in% Reductions(whole))
stopifnot(MPHI_UMAP %in% Reductions(mphi))

SAMPLE_COL <- resolve_first(
  colnames(mphi@meta.data),
  c("sample_for_annotation", "sample_FIXED2", "sample"),
  "Mphi sample column"
)
COND_COL <- resolve_first(
  colnames(mphi@meta.data),
  c("condition_FIXED2", "condition", "sample_for_annotation"),
  "Mphi condition column"
)

whole$layer1_v591 <- as.character(whole@meta.data[[WHOLE_ANNOT]])
mphi$class_v591 <- as.character(mphi@meta.data[[MPHI_CLASS]])
mphi$sample_v591 <- as.character(mphi@meta.data[[SAMPLE_COL]])
mphi$condition_v591 <- canonical_condition(mphi@meta.data[[COND_COL]])

CLASS_LEVELS <- c(
  "Anti-inflammatory-Mphi",
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi",
  "Other"
)
mphi$class_v591 <- factor(mphi$class_v591, levels = CLASS_LEVELS)

msg("Whole-cell cells: ", ncol(whole))
msg("Mphi cells: ", ncol(mphi))
print(table(mphi$class_v591, useNA = "ifany"))
print(table(mphi$sample_v591, useNA = "ifany"))

# Gene audit
wa <- resolve_genes(ALL_REQ, rownames(whole)) %>%
  rename(whole_resolved_gene = resolved_gene, whole_present = present)
ma <- resolve_genes(ALL_REQ, rownames(mphi)) %>%
  rename(mphi_resolved_gene = resolved_gene, mphi_present = present)

audit <- manifest %>%
  left_join(wa, by = "requested_gene") %>%
  left_join(ma, by = "requested_gene")

write.csv(audit, file.path(TAB, "01_gene_audit_wholecell_and_Mphi_v5.9.1.csv"), row.names = FALSE)
print(audit)

resolve_ordered <- function(audit_df, requested_order, resolved_col, present_col) {
  audit_df %>%
    filter(requested_gene %in% requested_order, .data[[present_col]]) %>%
    arrange(match(requested_gene, requested_order)) %>%
    pull(all_of(resolved_col)) %>%
    unique()
}

whole_sel <- resolve_ordered(wa, SEL_ORDER, "whole_resolved_gene", "whole_present")
whole_ref <- resolve_ordered(wa, REF_ORDER, "whole_resolved_gene", "whole_present")
mphi_sel <- resolve_ordered(ma, SEL_ORDER, "mphi_resolved_gene", "mphi_present")
mphi_ref <- resolve_ordered(ma, REF_ORDER, "mphi_resolved_gene", "mphi_present")
whole_all <- unique(c(whole_sel, whole_ref))
mphi_all <- unique(c(mphi_sel, mphi_ref))

# 1. Whole-cell UMAP
msg("Whole-cell UMAP...")
feature_pages(whole, whole_sel, WHOLE_UMAP,
              file.path(FIG, "01A_WHOLECELL_UMAP_selected_v5.9.1.pdf"),
              "Whole-cell UMAP | selected", pt = 0.13)
feature_pages(whole, whole_ref, WHOLE_UMAP,
              file.path(FIG, "01B_WHOLECELL_UMAP_reference_v5.9.1.pdf"),
              "Whole-cell UMAP | reference", pt = 0.13)

# 2. Whole-cell DotPlot + quantitative table
whole$layer1_v591 <- factor(whole$layer1_v591, levels = sort(unique(as.character(whole$layer1_v591))))

if (length(whole_sel)) {
  p <- DotPlot(whole, features = whole_sel, group.by = "layer1_v591", assay = ASSAY) +
    scale_color_gradient2(low = "#0033FF", mid = "#FFFFFF", high = "#FF1A1A", midpoint = 0) +
    RotatedAxis() + labs(title = "Whole-cell Layer1 | selected markers", x = NULL, y = NULL) +
    theme_classic(base_size = 8)
  save_pdf(p, file.path(FIG, "02A_WHOLECELL_Layer1_DotPlot_selected_v5.9.1.pdf"), 13, 8)
}

if (length(whole_ref)) {
  p <- DotPlot(whole, features = whole_ref, group.by = "layer1_v591", assay = ASSAY) +
    scale_color_gradient2(low = "#0033FF", mid = "#FFFFFF", high = "#FF1A1A", midpoint = 0) +
    RotatedAxis() + labs(title = "Whole-cell Layer1 | reference markers", x = NULL, y = NULL) +
    theme_classic(base_size = 8)
  save_pdf(p, file.path(FIG, "02B_WHOLECELL_Layer1_DotPlot_reference_v5.9.1.pdf"), 16, 8)
}

whole_sum <- summarize_by_group(whole, whole_all, whole$layer1_v591, "wholecell_Layer1", ASSAY)
write.csv(whole_sum, file.path(TAB, "02_wholecell_expression_positive_fraction_by_Layer1_v5.9.1.csv"), row.names = FALSE)

whole_sel_sum <- whole_sum %>% filter(gene %in% whole_sel)
if (nrow(whole_sel_sum)) {
  p <- ggplot(whole_sel_sum, aes(factor(gene, levels = whole_sel), group, fill = positive_percent)) +
    geom_tile(linewidth = 0.25) +
    scale_fill_gradientn(colours = c("#0033FF", "#FFFFFF", "#FF1A1A")) +
    labs(title = "Whole-cell Layer1 | selected marker RNA-positive fraction",
         x = NULL, y = NULL, fill = "% positive") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  save_pdf(p, file.path(FIG, "03_WHOLECELL_Layer1_positive_fraction_heatmap_selected_v5.9.1.pdf"), 13, 8)
}

# 3. Mphi UMAP
msg("Mphi UMAP...")
feature_pages(mphi, mphi_sel, MPHI_UMAP,
              file.path(FIG, "04A_MPHI_UMAP_selected_v5.9.1.pdf"),
              "Clean-B Mphi UMAP | selected", pt = 0.27)
feature_pages(mphi, mphi_ref, MPHI_UMAP,
              file.path(FIG, "04B_MPHI_UMAP_reference_v5.9.1.pdf"),
              "Clean-B Mphi UMAP | reference", pt = 0.27)

# 4. Mphi Sham vs Tx UMAP
mphi$ShamTx_v591 <- ifelse(mphi$condition_v591 %in% c("Sham", "Tx"), mphi$condition_v591, NA_character_)
mphi$ShamTx_v591 <- factor(mphi$ShamTx_v591, levels = c("Sham", "Tx"))
mphi_st <- subset(mphi, cells = colnames(mphi)[!is.na(mphi$ShamTx_v591)])

msg("Mphi Sham vs Tx UMAP...")
split_feature_pdf(mphi_st, mphi_sel, MPHI_UMAP, "ShamTx_v591",
                  file.path(FIG, "05A_MPHI_UMAP_SHAM_vs_TX_selected_v5.9.1.pdf"),
                  "Clean-B Mphi Sham vs Tx | selected", pt = 0.24)
split_feature_pdf(mphi_st, mphi_ref, MPHI_UMAP, "ShamTx_v591",
                  file.path(FIG, "05B_MPHI_UMAP_SHAM_vs_TX_reference_v5.9.1.pdf"),
                  "Clean-B Mphi Sham vs Tx | reference", pt = 0.24)

# 5. Mphi subtype DotPlot
msg("Mphi subtype DotPlot...")
if (length(mphi_sel)) {
  p <- DotPlot(mphi, features = mphi_sel, group.by = "class_v591", assay = ASSAY) +
    scale_color_gradient2(low = "#0033FF", mid = "#FFFFFF", high = "#FF1A1A", midpoint = 0) +
    RotatedAxis() + labs(title = "Mphi subtype | selected markers", x = NULL, y = NULL) +
    theme_classic(base_size = 8)
  save_pdf(p, file.path(FIG, "06A_MPHI_subtype_DotPlot_selected_v5.9.1.pdf"), 13, 5.5)
}
if (length(mphi_ref)) {
  p <- DotPlot(mphi, features = mphi_ref, group.by = "class_v591", assay = ASSAY) +
    scale_color_gradient2(low = "#0033FF", mid = "#FFFFFF", high = "#FF1A1A", midpoint = 0) +
    RotatedAxis() + labs(title = "Mphi subtype | reference markers", x = NULL, y = NULL) +
    theme_classic(base_size = 8)
  save_pdf(p, file.path(FIG, "06B_MPHI_subtype_DotPlot_reference_v5.9.1.pdf"), 16, 5.5)
}

# 6. Violin by marker group
msg("Violin plots...")
pdf(file.path(FIG, "07A_MPHI_subtype_Violin_selected_by_group_v5.9.1.pdf"), 12, 7, useDingbats = FALSE)
for (grp in names(SELECTED)) {
  req <- SELECTED[[grp]]
  gs <- ma %>% filter(requested_gene %in% req, mphi_present) %>%
    arrange(match(requested_gene, req)) %>% pull(mphi_resolved_gene) %>% unique()
  if (!length(gs)) next
  p <- VlnPlot(mphi, features = gs, group.by = "class_v591", pt.size = 0,
               ncol = min(2, length(gs)), assay = ASSAY) &
    theme_classic(base_size = 8) &
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p <- p + plot_annotation(title = paste0("Selected | ", grp))
  print(p)
}
dev.off()

pdf(file.path(FIG, "07B_MPHI_subtype_Violin_reference_by_group_v5.9.1.pdf"), 12, 7, useDingbats = FALSE)
for (grp in names(REFERENCE)) {
  req <- REFERENCE[[grp]]
  gs <- ma %>% filter(requested_gene %in% req, mphi_present) %>%
    arrange(match(requested_gene, req)) %>% pull(mphi_resolved_gene) %>% unique()
  if (!length(gs)) next
  p <- VlnPlot(mphi, features = gs, group.by = "class_v591", pt.size = 0,
               ncol = min(2, length(gs)), assay = ASSAY) &
    theme_classic(base_size = 8) &
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p <- p + plot_annotation(title = paste0("Reference | ", grp))
  print(p)
}
dev.off()

# 7. Mphi subtype quantitative table
msg("Subtype summaries...")
sub_sum <- summarize_by_group(mphi, mphi_all, mphi$class_v591, "Mphi_subtype", ASSAY)
write.csv(sub_sum, file.path(TAB, "03_Mphi_expression_positive_fraction_by_subtype_v5.9.1.csv"), row.names = FALSE)

# 8. allMphi_positive_fraction_by_sample
msg("allMphi_positive_fraction_by_sample...")
sample_sum <- summarize_by_group(mphi, mphi_all, mphi$sample_v591, "allMphi_sample", ASSAY) %>%
  rename(sample = group) %>%
  left_join(tibble(sample = mphi$sample_v591, condition = mphi$condition_v591) %>% distinct(),
            by = "sample")
write.csv(sample_sum, file.path(TAB, "04_allMphi_positive_fraction_by_sample_v5.9.1.csv"), row.names = FALSE)

preferred_samples <- c("STD_rep1", "STD", "CDHFD_rep1", "CDAHFD_rep1", "CDHFD", "CDAHFD",
                       "Sham1", "Sham20", "Tx17", "Tx5")
sample_levels <- unique(c(
  preferred_samples[preferred_samples %in% unique(mphi$sample_v591)],
  setdiff(sort(unique(mphi$sample_v591)), preferred_samples)
))

ss <- sample_sum %>% filter(gene %in% mphi_sel)
if (nrow(ss)) {
  p <- ggplot(ss, aes(factor(sample, levels = sample_levels), positive_fraction, group = 1)) +
    geom_line(linewidth = 0.45) + geom_point(size = 1.8) +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    labs(title = "allMphi positive fraction by sample | selected markers",
         x = NULL, y = "RNA-positive fraction") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  save_pdf(p, file.path(FIG, "08_allMphi_positive_fraction_by_sample_selected_v5.9.1.pdf"), 14, 11)
}

# 9-11. Heatmaps
mean_z_long <- heatmap_z(
  sub_sum, "group", "mean_expression", mphi_all, CLASS_LEVELS,
  "Mphi subtype mean expression | row z-score",
  file.path(FIG, "09_MPHI_subtype_mean_expression_Zscore_heatmap_v5.9.1.pdf"),
  18, 5.5
)
write.csv(mean_z_long, file.path(TAB, "05_MPHI_subtype_mean_expression_Zscore_long_v5.9.1.csv"), row.names = FALSE)

pos_z_long <- heatmap_z(
  sub_sum, "group", "positive_fraction", mphi_all, CLASS_LEVELS,
  "Mphi subtype RNA-positive fraction | row z-score",
  file.path(FIG, "10_MPHI_subtype_positive_fraction_Zscore_heatmap_v5.9.1.pdf"),
  18, 5.5
)
write.csv(pos_z_long, file.path(TAB, "06_MPHI_subtype_positive_fraction_Zscore_long_v5.9.1.csv"), row.names = FALSE)

sample_for_heat <- sample_sum %>% rename(group = sample)
sample_z_long <- heatmap_z(
  sample_for_heat, "group", "positive_fraction", mphi_all, sample_levels,
  "allMphi positive fraction by sample | row z-score",
  file.path(FIG, "11_allMphi_positive_fraction_by_sample_Zscore_heatmap_v5.9.1.pdf"),
  18, 6
)
write.csv(sample_z_long, file.path(TAB, "07_allMphi_positive_fraction_by_sample_Zscore_long_v5.9.1.csv"), row.names = FALSE)

# 12. Condition summary + Sham/Tx change table
msg("Condition summaries...")
cond_sum <- summarize_by_group(mphi, mphi_all, mphi$condition_v591, "allMphi_condition", ASSAY)
write.csv(cond_sum, file.path(TAB, "08_allMphi_expression_positive_fraction_by_condition_v5.9.1.csv"), row.names = FALSE)

st <- cond_sum %>%
  filter(group %in% c("Sham", "Tx")) %>%
  select(gene, group, positive_fraction, mean_expression) %>%
  pivot_wider(names_from = group, values_from = c(positive_fraction, mean_expression)) %>%
  mutate(
    delta_positive_fraction_Tx_minus_Sham = positive_fraction_Tx - positive_fraction_Sham,
    log2_mean_expression_Tx_vs_Sham =
      log2((mean_expression_Tx + 1e-6) / (mean_expression_Sham + 1e-6))
  )
write.csv(st, file.path(TAB, "09_allMphi_SHAM_vs_TX_marker_summary_v5.9.1.csv"), row.names = FALSE)


# 13. Explicit Sham vs Tx comparison block
# ----------------------------------------------------------------------

msg("Explicit Sham vs Tx comparison block...")

# Focused four biological samples
SHAM_TX_SAMPLE_ORDER <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

sham_tx_sample_levels <- SHAM_TX_SAMPLE_ORDER[
  SHAM_TX_SAMPLE_ORDER %in%
    unique(
      mphi$sample_v591
    )
]

# --------------------------------------------------------------
# 13A. All-Mphi Sham vs Tx violin plots
# --------------------------------------------------------------

pdf(
  file.path(
    FIG,
    "12A_allMphi_Violin_SHAM_vs_TX_selected_v5.9.1.pdf"
  ),
  width = 10,
  height = 5.5,
  useDingbats = FALSE
)

for (g in mphi_sel) {

  p <- VlnPlot(
    mphi_st,
    features = g,
    group.by = "ShamTx_v591",
    pt.size = 0,
    assay = ASSAY
  ) +
    labs(
      title = paste0(
        g,
        " | all Mphi | Sham vs Tx"
      ),
      x = NULL,
      y = "RNA expression"
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  print(p)
}

dev.off()

pdf(
  file.path(
    FIG,
    "12B_allMphi_Violin_SHAM_vs_TX_reference_v5.9.1.pdf"
  ),
  width = 10,
  height = 5.5,
  useDingbats = FALSE
)

for (g in mphi_ref) {

  p <- VlnPlot(
    mphi_st,
    features = g,
    group.by = "ShamTx_v591",
    pt.size = 0,
    assay = ASSAY
  ) +
    labs(
      title = paste0(
        g,
        " | all Mphi | Sham vs Tx | reference"
      ),
      x = NULL,
      y = "RNA expression"
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  print(p)
}

dev.off()

# --------------------------------------------------------------
# 13B. Sham vs Tx within each macrophage subtype
# --------------------------------------------------------------

pdf(
  file.path(
    FIG,
    "13_MPHI_subtype_Violin_SHAM_vs_TX_selected_v5.9.1.pdf"
  ),
  width = 14,
  height = 6.5,
  useDingbats = FALSE
)

for (g in mphi_sel) {

  p <- VlnPlot(
    mphi_st,
    features = g,
    group.by = "class_v591",
    split.by = "ShamTx_v591",
    split.plot = TRUE,
    pt.size = 0,
    assay = ASSAY
  ) +
    labs(
      title = paste0(
        g,
        " | Mphi subtype | Sham vs Tx"
      ),
      x = NULL,
      y = "RNA expression"
    ) +
    theme_classic(
      base_size = 8
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  print(p)
}

dev.off()

# --------------------------------------------------------------
# 13C. Focused allMphi_positive_fraction_by_sample:
#      Sham1 / Sham20 / Tx17 / Tx5
# --------------------------------------------------------------

focused_sample_sum <- sample_sum %>%
  filter(
    sample %in%
      sham_tx_sample_levels
  )

write.csv(
  focused_sample_sum,
  file.path(
    TAB,
    "10_allMphi_positive_fraction_by_sample_SHAM_TX_only_v5.9.1.csv"
  ),
  row.names = FALSE
)

focused_selected <- focused_sample_sum %>%
  filter(
    gene %in%
      mphi_sel
  )

if (
  nrow(
    focused_selected
  )
) {

  p <- ggplot(
    focused_selected,
    aes(
      x = factor(
        sample,
        levels = sham_tx_sample_levels
      ),
      y = positive_fraction,
      group = 1
    )
  ) +
    geom_line(
      linewidth = 0.5
    ) +
    geom_point(
      size = 2
    ) +
    facet_wrap(
      ~ gene,
      scales = "free_y",
      ncol = 4
    ) +
    labs(
      title =
        "allMphi positive fraction | Sham1/Sham20 vs Tx17/Tx5",
      x = NULL,
      y = "RNA-positive fraction"
    ) +
    theme_classic(
      base_size = 8
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  save_pdf(
    p,
    file.path(
      FIG,
      "14_allMphi_positive_fraction_by_sample_SHAM_TX_selected_v5.9.1.pdf"
    ),
    14,
    11
  )
}

# --------------------------------------------------------------
# 13D. Focused positive-fraction heatmap by biological sample
# --------------------------------------------------------------

if (
  nrow(
    focused_sample_sum
  )
) {

  pos_sample <- focused_sample_sum %>%
    filter(
      gene %in%
        mphi_all
    ) %>%
    mutate(
      gene = factor(
        gene,
        levels = mphi_all
      ),
      sample = factor(
        sample,
        levels = sham_tx_sample_levels
      )
    )

  p <- ggplot(
    pos_sample,
    aes(
      x = gene,
      y = sample,
      fill = positive_fraction
    )
  ) +
    geom_tile(
      linewidth = 0.25
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
        "allMphi positive fraction | Sham1/Sham20 vs Tx17/Tx5",
      x = NULL,
      y = NULL,
      fill = "Positive fraction"
    ) +
    theme_classic(
      base_size = 8
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  save_pdf(
    p,
    file.path(
      FIG,
      "15_allMphi_positive_fraction_SHAM_TX_sample_heatmap_v5.9.1.pdf"
    ),
    18,
    5
  )
}

# --------------------------------------------------------------
# 13E. Focused mean-expression heatmap by biological sample
# --------------------------------------------------------------

if (
  nrow(
    focused_sample_sum
  )
) {

  mean_sample <- focused_sample_sum %>%
    filter(
      gene %in%
        mphi_all
    ) %>%
    mutate(
      gene = factor(
        gene,
        levels = mphi_all
      ),
      sample = factor(
        sample,
        levels = sham_tx_sample_levels
      )
    )

  p <- ggplot(
    mean_sample,
    aes(
      x = gene,
      y = sample,
      fill = mean_expression
    )
  ) +
    geom_tile(
      linewidth = 0.25
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
        "allMphi mean expression | Sham1/Sham20 vs Tx17/Tx5",
      x = NULL,
      y = NULL,
      fill = "Mean expression"
    ) +
    theme_classic(
      base_size = 8
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  save_pdf(
    p,
    file.path(
      FIG,
      "16_allMphi_mean_expression_SHAM_TX_sample_heatmap_v5.9.1.pdf"
    ),
    18,
    5
  )
}

# --------------------------------------------------------------
# 13F. Sham vs Tx condition-level two-column heatmaps
# --------------------------------------------------------------

sham_tx_condition_raw <- cond_sum %>%
  filter(
    group %in%
      c(
        "Sham",
        "Tx"
      )
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Sham",
        "Tx"
      )
    ),
    gene = factor(
      gene,
      levels = mphi_all
    )
  )

write.csv(
  sham_tx_condition_raw,
  file.path(
    TAB,
    "11_allMphi_SHAM_vs_TX_condition_raw_v5.9.1.csv"
  ),
  row.names = FALSE
)

if (
  nrow(
    sham_tx_condition_raw
  )
) {

  p <- ggplot(
    sham_tx_condition_raw,
    aes(
      x = gene,
      y = group,
      fill = positive_fraction
    )
  ) +
    geom_tile(
      linewidth = 0.3
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
        "allMphi RNA-positive fraction | Sham vs Tx",
      x = NULL,
      y = NULL,
      fill = "Positive fraction"
    ) +
    theme_classic(
      base_size = 8
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  save_pdf(
    p,
    file.path(
      FIG,
      "17A_allMphi_SHAM_vs_TX_positive_fraction_heatmap_v5.9.1.pdf"
    ),
    18,
    3.3
  )

  p <- ggplot(
    sham_tx_condition_raw,
    aes(
      x = gene,
      y = group,
      fill = mean_expression
    )
  ) +
    geom_tile(
      linewidth = 0.3
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
        "allMphi mean RNA expression | Sham vs Tx",
      x = NULL,
      y = NULL,
      fill = "Mean expression"
    ) +
    theme_classic(
      base_size = 8
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )

  save_pdf(
    p,
    file.path(
      FIG,
      "17B_allMphi_SHAM_vs_TX_mean_expression_heatmap_v5.9.1.pdf"
    ),
    18,
    3.3
  )
}

# --------------------------------------------------------------
# 13G. Sham vs Tx delta summary
# --------------------------------------------------------------

sham_tx_delta <- st %>%
  mutate(
    selected_marker =
      gene %in%
        mphi_sel
  ) %>%
  arrange(
    desc(
      selected_marker
    ),
    desc(
      abs(
        delta_positive_fraction_Tx_minus_Sham
      )
    )
  )

write.csv(
  sham_tx_delta,
  file.path(
    TAB,
    "12_allMphi_SHAM_vs_TX_delta_summary_v5.9.1.csv"
  ),
  row.names = FALSE
)

# Index / metadata
index <- tibble(
  output_type = c(
    rep("Figure", 27),
    rep("Table", 13)
  ),
  file = c(
    "01A_WHOLECELL_UMAP_selected_v5.9.1.pdf",
    "01B_WHOLECELL_UMAP_reference_v5.9.1.pdf",
    "02A_WHOLECELL_Layer1_DotPlot_selected_v5.9.1.pdf",
    "02B_WHOLECELL_Layer1_DotPlot_reference_v5.9.1.pdf",
    "03_WHOLECELL_Layer1_positive_fraction_heatmap_selected_v5.9.1.pdf",
    "04A_MPHI_UMAP_selected_v5.9.1.pdf",
    "04B_MPHI_UMAP_reference_v5.9.1.pdf",
    "05A_MPHI_UMAP_SHAM_vs_TX_selected_v5.9.1.pdf",
    "05B_MPHI_UMAP_SHAM_vs_TX_reference_v5.9.1.pdf",
    "06A_MPHI_subtype_DotPlot_selected_v5.9.1.pdf",
    "06B_MPHI_subtype_DotPlot_reference_v5.9.1.pdf",
    "07A_MPHI_subtype_Violin_selected_by_group_v5.9.1.pdf",
    "07B_MPHI_subtype_Violin_reference_by_group_v5.9.1.pdf",
    "08_allMphi_positive_fraction_by_sample_selected_v5.9.1.pdf",
    "09_MPHI_subtype_mean_expression_Zscore_heatmap_v5.9.1.pdf",
    "10_MPHI_subtype_positive_fraction_Zscore_heatmap_v5.9.1.pdf",
    "11_allMphi_positive_fraction_by_sample_Zscore_heatmap_v5.9.1.pdf",
    "12A_allMphi_Violin_SHAM_vs_TX_selected_v5.9.1.pdf",
    "12B_allMphi_Violin_SHAM_vs_TX_reference_v5.9.1.pdf",
    "13_MPHI_subtype_Violin_SHAM_vs_TX_selected_v5.9.1.pdf",
    "14_allMphi_positive_fraction_by_sample_SHAM_TX_selected_v5.9.1.pdf",
    "15_allMphi_positive_fraction_SHAM_TX_sample_heatmap_v5.9.1.pdf",
    "16_allMphi_mean_expression_SHAM_TX_sample_heatmap_v5.9.1.pdf",
    "17A_allMphi_SHAM_vs_TX_positive_fraction_heatmap_v5.9.1.pdf",
    "17B_allMphi_SHAM_vs_TX_mean_expression_heatmap_v5.9.1.pdf",
    "05A_MPHI_UMAP_SHAM_vs_TX_selected_v5.9.1.pdf",
    "05B_MPHI_UMAP_SHAM_vs_TX_reference_v5.9.1.pdf",

    "00_marker_manifest_v5.9.1.csv",
    "01_gene_audit_wholecell_and_Mphi_v5.9.1.csv",
    "02_wholecell_expression_positive_fraction_by_Layer1_v5.9.1.csv",
    "03_Mphi_expression_positive_fraction_by_subtype_v5.9.1.csv",
    "04_allMphi_positive_fraction_by_sample_v5.9.1.csv",
    "05_MPHI_subtype_mean_expression_Zscore_long_v5.9.1.csv",
    "06_MPHI_subtype_positive_fraction_Zscore_long_v5.9.1.csv",
    "07_allMphi_positive_fraction_by_sample_Zscore_long_v5.9.1.csv",
    "08_allMphi_expression_positive_fraction_by_condition_v5.9.1.csv",
    "09_allMphi_SHAM_vs_TX_marker_summary_v5.9.1.csv",
    "10_allMphi_positive_fraction_by_sample_SHAM_TX_only_v5.9.1.csv",
    "11_allMphi_SHAM_vs_TX_condition_raw_v5.9.1.csv",
    "12_allMphi_SHAM_vs_TX_delta_summary_v5.9.1.csv"
  )
)
write.csv(index, file.path(OUT, "OUTPUT_INDEX_v5.9.1.csv"), row.names = FALSE)

metadata <- tibble(
  parameter = c("version", "whole_RDS", "mphi_RDS", "whole_annotation", "mphi_class",
                "mphi_sample_col", "mphi_condition_col", "whole_UMAP", "mphi_UMAP"),
  value = c("v5.9.1", WHOLE_RDS, MPHI_RDS, WHOLE_ANNOT, MPHI_CLASS,
            SAMPLE_COL, COND_COL, WHOLE_UMAP, MPHI_UMAP)
)
write.csv(metadata, file.path(LOG, "analysis_metadata_v5.9.1.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(LOG, "sessionInfo_v5.9.1.txt"))

msg("DONE.")
msg("Output: ", OUT)
