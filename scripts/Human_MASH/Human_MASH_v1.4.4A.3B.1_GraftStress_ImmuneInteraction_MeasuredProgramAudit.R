#!/usr/bin/env Rscript

# ==============================================================================
# Human_MASH_v1.4.4A.3B.1_GraftStress_ImmuneInteraction_MeasuredProgramAudit.R
#
# Purpose
#   Descriptive, replicate-aware audit of probe-measurable human graft programs:
#     - IFN response
#     - TNF/NF-kB response
#     - cellular / oxidative stress
#     - apoptosis / death / survival
#     - hypoxia / ischemia adaptation
#     - graft-host innate interface (exploratory)
#
#   Antigen-presentation coverage is only 4/16 genes in v1.4.4A.3B.0 and is
#   therefore NOT treated as a primary aggregate module. Its four measured genes
#   are retained as exploratory gene-level readouts only.
#
# Biological comparison
#   SfLB1/SfLB2 = pre-transplant, non-engineered cultured organoids
#   Tx17/Tx5    = engrafted human cells derived from IL10-engineered organoids
#
# Critical causal boundary
#   SfLB -> Tx is NOT an IL10-only perturbation. It combines:
#     transgene status + transplantation/engraftment + host microenvironment
#     + post-engraftment selection/state transition
#     + possible residual innate xenogeneic immune pressure in NOD/scid recipients.
#
# Interpretation boundary
#   This analysis can identify graft stress / inflammatory-adaptation signatures.
#   It CANNOT by itself establish xenogeneic rejection, because key classical
#   antigen-presentation genes (including B2M and multiple HLA genes) are outside
#   the 9,137-feature Flex space.
#
# Statistics
#   n=2 biological samples per group; descriptive only, no inferential test.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_iPSC_derived_v1.4.2A.1/objects/",
  "Human_MASH_iPSC_derived_v1.4.2A.1_hepatocyte_state_interpreted.rds"
)

PREFLIGHT_DIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.3B.0/",
  "analysis/v1.4.4A.3B.0_graft_stress_xenogeneic_pressure_feature_preflight/tables"
)

COVERAGE_CSV <- file.path(
  PREFLIGHT_DIR,
  "01_candidate_gene_feature_coverage_long_v1.4.4A.3B.0.csv"
)

OUT_ROOT <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4A.3B.1/",
  "analysis/v1.4.4A.3B.1_graft_stress_immune_interaction_measured_program_audit"
)

TAB_DIR <- file.path(OUT_ROOT, "tables")
FIG_DIR <- file.path(OUT_ROOT, "figures")
LOG_DIR <- file.path(OUT_ROOT, "logs")

dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, "v1.4.4A.3B.1_run_log.txt")
sink(LOG_FILE, split = TRUE)
on.exit(sink(), add = TRUE)

cat("============================================================\n")
cat("Human MASH v1.4.4A.3B.1\n")
cat("Graft stress / immune-interaction measured-program audit\n")
cat("============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Input:", INPUT_RDS, "\n")
cat("Preflight:", COVERAGE_CSV, "\n")
cat("Output:", OUT_ROOT, "\n\n")

for (f in c(INPUT_RDS, COVERAGE_CSV)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

# ------------------------------------------------------------------------------
# 1. Load frozen object and preflight
# ------------------------------------------------------------------------------

obj <- readRDS(INPUT_RDS)

if (!inherits(obj, "Seurat")) stop("Input object is not a Seurat object.")
if (!"RNA" %in% Assays(obj)) stop("RNA assay missing.")
if (ncol(obj) != 43415L) stop("Unexpected cell count: ", ncol(obj))

rna_features <- rownames(obj[["RNA"]])
if (length(rna_features) != 9137L) {
  stop("Unexpected RNA feature count: ", length(rna_features))
}

required_meta <- c(
  "sample_v141",
  "broad_lineage_manual_v142A",
  "manual_state_v142A",
  "annotation_status_v142A"
)
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing metadata: ", paste(missing_meta, collapse = ", "))
}

coverage <- read.csv(
  COVERAGE_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

needed_cov <- c("panel", "gene", "present_in_Flex_RNA")
if (!all(needed_cov %in% colnames(coverage))) {
  stop("Preflight coverage file schema mismatch.")
}

coverage$present_in_Flex_RNA <- as.logical(coverage$present_in_Flex_RNA)

expected_present_n <- c(
  IFN_response = 11L,
  Antigen_presentation = 4L,
  TNF_NFkB_response = 10L,
  Cellular_oxidative_stress = 12L,
  Apoptosis_death_survival = 8L,
  Hypoxia_ischemia_adaptation = 10L,
  Graft_host_innate_interface = 7L
)

observed_present_n <- tapply(
  coverage$present_in_Flex_RNA,
  coverage$panel,
  sum
)

for (nm in names(expected_present_n)) {
  if (is.na(observed_present_n[nm]) ||
      observed_present_n[nm] != expected_present_n[nm]) {
    stop(
      "Preflight measured-gene count mismatch for ", nm,
      ": observed ", observed_present_n[nm],
      ", expected ", expected_present_n[nm]
    )
  }
}

cat("Preflight measured-gene counts: PASS\n")

# ------------------------------------------------------------------------------
# 2. Freeze measured-gene modules
# ------------------------------------------------------------------------------

primary_modules <- c(
  "IFN_response",
  "TNF_NFkB_response",
  "Cellular_oxidative_stress",
  "Apoptosis_death_survival",
  "Hypoxia_ischemia_adaptation"
)

exploratory_modules <- c(
  "Graft_host_innate_interface"
)

antigen_panel <- "Antigen_presentation"

module_genes <- lapply(
  c(primary_modules, exploratory_modules),
  function(mm) {
    coverage$gene[
      coverage$panel == mm &
        coverage$present_in_Flex_RNA
    ]
  }
)
names(module_genes) <- c(primary_modules, exploratory_modules)

antigen_genes <- coverage$gene[
  coverage$panel == antigen_panel &
    coverage$present_in_Flex_RNA
]

freeze_tbl <- do.call(
  rbind,
  lapply(names(module_genes), function(mm) {
    data.frame(
      program = mm,
      gene = module_genes[[mm]],
      analysis_role = ifelse(
        mm %in% primary_modules,
        "PRIMARY_MEASURED_PROGRAM",
        "EXPLORATORY_MEASURED_PROGRAM"
      ),
      stringsAsFactors = FALSE
    )
  })
)

antigen_freeze <- data.frame(
  program = antigen_panel,
  gene = antigen_genes,
  analysis_role = "GENE_LEVEL_ONLY_INSUFFICIENT_COVERAGE_FOR_PRIMARY_MODULE",
  stringsAsFactors = FALSE
)

freeze_all <- rbind(freeze_tbl, antigen_freeze)

write.csv(
  freeze_all,
  file.path(TAB_DIR, "00_measured_program_definition_freeze_v1.4.4A.3B.1.csv"),
  row.names = FALSE
)

cat("\n=== MEASURED PROGRAM DEFINITIONS ===\n")
for (mm in names(module_genes)) {
  cat(mm, " (n=", length(module_genes[[mm]]), "): ",
      paste(module_genes[[mm]], collapse = ", "), "\n", sep = "")
}
cat(
  antigen_panel, " (gene-level only; n=", length(antigen_genes), "): ",
  paste(antigen_genes, collapse = ", "), "\n", sep = ""
)

# ------------------------------------------------------------------------------
# 3. Extract raw counts and compute transparent cell-level normalized expression
# ------------------------------------------------------------------------------

counts <- tryCatch(
  {
    LayerData(obj, assay = "RNA", layer = "counts")
  },
  error = function(e1) {
    tryCatch(
      {
        GetAssayData(obj, assay = "RNA", layer = "counts")
      },
      error = function(e2) {
        GetAssayData(obj, assay = "RNA", slot = "counts")
      }
    )
  }
)

if (ncol(counts) != ncol(obj) || nrow(counts) != length(rna_features)) {
  stop("RNA counts matrix dimension mismatch.")
}

libsize <- Matrix::colSums(counts)
if (any(!is.finite(libsize)) || any(libsize <= 0)) {
  stop("Invalid RNA library sizes detected.")
}

all_measured_genes <- sort(unique(c(
  unlist(module_genes, use.names = FALSE),
  antigen_genes
)))

missing_from_counts <- setdiff(all_measured_genes, rownames(counts))
if (length(missing_from_counts) > 0) {
  stop(
    "Measured preflight genes missing from counts matrix: ",
    paste(missing_from_counts, collapse = ", ")
  )
}

sel_counts <- counts[all_measured_genes, , drop = FALSE]

# 60 x 43,415 scale is small enough for a dense selected-gene matrix.
sel_dense <- as.matrix(sel_counts)
log_cp10k <- log1p(
  sweep(sel_dense, 2, libsize, "/") * 10000
)

# Equal-weight gene z-score across all iPSC-derived cells.
gene_mu <- rowMeans(log_cp10k)
gene_sd <- apply(log_cp10k, 1, sd)

zero_sd <- names(gene_sd)[!is.finite(gene_sd) | gene_sd == 0]
if (length(zero_sd) > 0) {
  stop(
    "Measured genes with zero/non-finite SD: ",
    paste(zero_sd, collapse = ", ")
  )
}

gene_z <- sweep(log_cp10k, 1, gene_mu, "-")
gene_z <- sweep(gene_z, 1, gene_sd, "/")

# ------------------------------------------------------------------------------
# 4. Cell-level program scores
# ------------------------------------------------------------------------------

program_scores <- do.call(
  cbind,
  lapply(names(module_genes), function(mm) {
    colMeans(gene_z[module_genes[[mm]], , drop = FALSE])
  })
)

colnames(program_scores) <- names(module_genes)
rownames(program_scores) <- colnames(obj)

md <- obj@meta.data
md$cell_id <- rownames(md)

score_df <- cbind(
  md[, c(
    "cell_id",
    "sample_v141",
    "broad_lineage_manual_v142A",
    "manual_state_v142A",
    "annotation_status_v142A"
  )],
  as.data.frame(program_scores)
)

names(score_df)[names(score_df) == "sample_v141"] <- "sample"
names(score_df)[names(score_df) == "broad_lineage_manual_v142A"] <- "lineage"
names(score_df)[names(score_df) == "manual_state_v142A"] <- "state"

sample_order <- c("SfLB1", "SfLB2", "Tx17", "Tx5")
score_df$sample <- factor(score_df$sample, levels = sample_order)
score_df$bio_group <- ifelse(
  score_df$sample %in% c("SfLB1", "SfLB2"),
  "PreTransplant_SfLB",
  "Engrafted_Tx"
)

# ------------------------------------------------------------------------------
# 5. Aggregation helpers
# ------------------------------------------------------------------------------

summarize_scores <- function(df, group_cols) {
  programs <- names(module_genes)

  key <- interaction(
    df[, group_cols, drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )

  idx <- split(seq_len(nrow(df)), key)

  out <- do.call(
    rbind,
    lapply(idx, function(ii) {
      base <- df[ii[1], group_cols, drop = FALSE]

      vals <- do.call(
        cbind,
        lapply(programs, function(mm) {
          x <- df[[mm]][ii]
          c(
            mean = mean(x, na.rm = TRUE),
            median = median(x, na.rm = TRUE)
          )
        })
      )

      # Flatten as program_mean / program_median.
      flat <- c()
      for (j in seq_along(programs)) {
        flat[paste0(programs[j], "_mean")] <- vals["mean", j]
        flat[paste0(programs[j], "_median")] <- vals["median", j]
      }

      cbind(
        base,
        n_cells = length(ii),
        as.data.frame(as.list(flat), check.names = FALSE)
      )
    })
  )

  rownames(out) <- NULL
  out
}

whole_sample <- summarize_scores(score_df, c("sample", "bio_group"))

lineage_sample <- summarize_scores(
  score_df,
  c("sample", "bio_group", "lineage")
)

state_sample <- summarize_scores(
  score_df,
  c("sample", "bio_group", "state")
)

write.csv(
  whole_sample,
  file.path(TAB_DIR, "01_whole_sample_measured_program_scores_v1.4.4A.3B.1.csv"),
  row.names = FALSE
)

write.csv(
  lineage_sample,
  file.path(TAB_DIR, "02_broad_lineage_measured_program_scores_by_sample_v1.4.4A.3B.1.csv"),
  row.names = FALSE
)

write.csv(
  state_sample,
  file.path(TAB_DIR, "03_manual_state_measured_program_scores_by_sample_v1.4.4A.3B.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Replicate-aware descriptive comparison
# ------------------------------------------------------------------------------

classify_2v2 <- function(sf1, sf2, tx17, tx5) {
  s <- c(sf1, sf2)
  t <- c(tx17, tx5)

  if (anyNA(c(s, t))) return("Not_comparable_missing_sample")

  if (min(t) > max(s)) {
    "Tx_higher_complete_replicate_range_separation"
  } else if (max(t) < min(s)) {
    "Tx_lower_complete_replicate_range_separation"
  } else {
    "Replicate_ranges_overlap"
  }
}

compare_from_table <- function(tab, strat_col = NULL) {
  programs <- names(module_genes)

  strata <- if (is.null(strat_col)) {
    "__WHOLE__"
  } else {
    unique(as.character(tab[[strat_col]]))
  }

  out <- list()
  kk <- 1L

  for (ss in strata) {
    z <- if (is.null(strat_col)) {
      tab
    } else {
      tab[as.character(tab[[strat_col]]) == ss, , drop = FALSE]
    }

    if (!all(sample_order %in% as.character(z$sample))) next

    z <- z[match(sample_order, as.character(z$sample)), , drop = FALSE]

    for (mm in programs) {
      coln <- paste0(mm, "_mean")
      vals <- setNames(as.numeric(z[[coln]]), as.character(z$sample))

      out[[kk]] <- data.frame(
        stratum = ss,
        program = mm,
        SfLB1 = vals["SfLB1"],
        SfLB2 = vals["SfLB2"],
        Tx17 = vals["Tx17"],
        Tx5 = vals["Tx5"],
        SfLB_mean = mean(vals[c("SfLB1", "SfLB2")]),
        Tx_mean = mean(vals[c("Tx17", "Tx5")]),
        Tx_minus_SfLB = mean(vals[c("Tx17", "Tx5")]) -
          mean(vals[c("SfLB1", "SfLB2")]),
        replicate_relation = classify_2v2(
          vals["SfLB1"], vals["SfLB2"],
          vals["Tx17"], vals["Tx5"]
        ),
        stringsAsFactors = FALSE
      )
      kk <- kk + 1L
    }
  }

  do.call(rbind, out)
}

whole_compare <- compare_from_table(whole_sample, NULL)
whole_compare$interpretation_boundary <-
  "WHOLE_GRAFT_COMPOSITION_CONFOUNDED_DESCRIPTIVE_ONLY"

write.csv(
  whole_compare,
  file.path(TAB_DIR, "04_whole_graft_program_transition_COMPOSITION_CONFOUNDED_v1.4.4A.3B.1.csv"),
  row.names = FALSE
)

# Primary within-lineage comparison only in broad lineages retained as major
# compartments in both SfLB and Tx (from v1.4.4A.3.1).
comparable_lineages <- c("Hepatocyte_like", "Mesenchymal_like")

lineage_compare_input <- lineage_sample[
  as.character(lineage_sample$lineage) %in% comparable_lineages,
  ,
  drop = FALSE
]

lineage_compare <- compare_from_table(
  lineage_compare_input,
  "lineage"
)

lineage_compare$interpretation_boundary <-
  "WITHIN_BROAD_LINEAGE_ENGRAFTMENT_ASSOCIATED_NOT_IL10_CAUSAL"

write.csv(
  lineage_compare,
  file.path(TAB_DIR, "05_comparable_lineage_program_transition_v1.4.4A.3B.1.csv"),
  row.names = FALSE
)

cat("\n=== PRIMARY WITHIN-LINEAGE PROGRAM TRANSITIONS ===\n")
print(
  lineage_compare[
    lineage_compare$program %in% primary_modules,
    c(
      "stratum", "program",
      "SfLB1", "SfLB2", "Tx17", "Tx5",
      "Tx_minus_SfLB", "replicate_relation"
    )
  ],
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Gene-level audit: sample x lineage
# ------------------------------------------------------------------------------

gene_metric_rows <- list()
kk <- 1L

for (ss in sample_order) {
  for (ll in unique(as.character(score_df$lineage))) {
    ii <- which(
      as.character(score_df$sample) == ss &
        as.character(score_df$lineage) == ll
    )

    if (length(ii) == 0) next

    for (gg in all_measured_genes) {
      gidx <- match(gg, rownames(sel_dense))
      x_counts <- sel_dense[gidx, ii]
      x_log <- log_cp10k[gidx, ii]

      gene_metric_rows[[kk]] <- data.frame(
        sample = ss,
        lineage = ll,
        gene = gg,
        n_cells = length(ii),
        positive_cells = sum(x_counts > 0),
        positive_fraction = mean(x_counts > 0),
        mean_logCP10k = mean(x_log),
        median_logCP10k = median(x_log),
        total_UMI = sum(x_counts),
        UMI_per_1000_cells = 1000 * sum(x_counts) / length(ii),
        stringsAsFactors = FALSE
      )
      kk <- kk + 1L
    }
  }
}

gene_metrics <- do.call(rbind, gene_metric_rows)

write.csv(
  gene_metrics,
  file.path(TAB_DIR, "06_measured_gene_metrics_by_sample_lineage_v1.4.4A.3B.1.csv"),
  row.names = FALSE
)

# Antigen-presentation gene-level exploratory table only.
antigen_gene_metrics <- gene_metrics[
  gene_metrics$gene %in% antigen_genes,
  ,
  drop = FALSE
]

write.csv(
  antigen_gene_metrics,
  file.path(TAB_DIR, "07_antigen_presentation_MEASURED_GENES_ONLY_exploratory_v1.4.4A.3B.1.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Module overlap audit
# ------------------------------------------------------------------------------

module_names <- names(module_genes)
overlap <- matrix(
  0L,
  nrow = length(module_names),
  ncol = length(module_names),
  dimnames = list(module_names, module_names)
)

for (i in seq_along(module_names)) {
  for (j in seq_along(module_names)) {
    overlap[i, j] <- length(intersect(
      module_genes[[module_names[i]]],
      module_genes[[module_names[j]]]
    ))
  }
}

write.csv(
  overlap,
  file.path(TAB_DIR, "08_measured_program_gene_overlap_matrix_v1.4.4A.3B.1.csv"),
  row.names = TRUE
)

# ------------------------------------------------------------------------------
# 9. Figures
# ------------------------------------------------------------------------------

# Long plotting table for sample-level means.
make_long <- function(tab, strat_label_col = NULL) {
  rows <- list()
  k <- 1L

  for (i in seq_len(nrow(tab))) {
    for (mm in names(module_genes)) {
      rows[[k]] <- data.frame(
        sample = as.character(tab$sample[i]),
        bio_group = as.character(tab$bio_group[i]),
        stratum = if (is.null(strat_label_col)) {
          "Whole graft"
        } else {
          as.character(tab[[strat_label_col]][i])
        },
        program = mm,
        score = as.numeric(tab[[paste0(mm, "_mean")]][i]),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  do.call(rbind, rows)
}

whole_long <- make_long(whole_sample)
whole_long$bio_group <- factor(
  whole_long$bio_group,
  levels = c("PreTransplant_SfLB", "Engrafted_Tx")
)
whole_long$sample <- factor(whole_long$sample, levels = sample_order)

p1 <- ggplot(
  whole_long,
  aes(x = bio_group, y = score)
) +
  geom_point(
    aes(shape = sample),
    size = 2.6,
    position = position_jitter(width = 0.05, height = 0)
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.38,
    linewidth = 0.45
  ) +
  facet_wrap(~ program, scales = "free_y", ncol = 3) +
  labs(
    title = "Whole-graft measured stress / immune-interaction programs",
    subtitle = "Composition-confounded descriptive view; n=2 per group",
    x = NULL,
    y = "Mean equal-weight gene z-score",
    shape = "Sample"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG_DIR, "Fig01_whole_graft_programs_COMPOSITION_CONFOUNDED_v1.4.4A.3B.1.pdf"),
  p1, width = 11.5, height = 7.5, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig01_whole_graft_programs_COMPOSITION_CONFOUNDED_v1.4.4A.3B.1.png"),
  p1, width = 11.5, height = 7.5, dpi = 300
)

lineage_long <- make_long(
  lineage_compare_input,
  "lineage"
)
lineage_long$bio_group <- factor(
  lineage_long$bio_group,
  levels = c("PreTransplant_SfLB", "Engrafted_Tx")
)
lineage_long$sample <- factor(lineage_long$sample, levels = sample_order)

p2 <- ggplot(
  lineage_long,
  aes(x = bio_group, y = score)
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
    linewidth = 0.45
  ) +
  facet_grid(program ~ stratum, scales = "free_y") +
  labs(
    title = "Within-lineage graft stress / immune-interaction programs",
    subtitle = "Hepatocyte-like and Mesenchymal-like compartments; descriptive n=2 vs n=2",
    x = NULL,
    y = "Mean equal-weight gene z-score",
    shape = "Sample"
  ) +
  theme_classic(base_size = 9.2) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG_DIR, "Fig02_within_comparable_lineage_programs_v1.4.4A.3B.1.pdf"),
  p2, width = 9.5, height = 13.0, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig02_within_comparable_lineage_programs_v1.4.4A.3B.1.png"),
  p2, width = 9.5, height = 13.0, dpi = 280
)

# Heatmap of sample-lineage module means.
heat <- make_long(lineage_sample, "lineage")
heat$sample_lineage <- paste(heat$sample, heat$stratum, sep = " | ")
heat$program <- factor(heat$program, levels = rev(names(module_genes)))

sample_lineage_order <- unique(
  paste(
    lineage_sample$sample,
    lineage_sample$lineage,
    sep = " | "
  )
)
heat$sample_lineage <- factor(
  heat$sample_lineage,
  levels = sample_lineage_order
)

p3 <- ggplot(
  heat,
  aes(x = sample_lineage, y = program, fill = score)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0
  ) +
  labs(
    title = "Measured graft program atlas by sample and broad lineage",
    subtitle = "Equal-weight gene z-score; blue = lower, red = higher relative to all iPSC-derived cells",
    x = "Sample | broad lineage",
    y = NULL,
    fill = "Program score"
  ) +
  theme_minimal(base_size = 9.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7.2),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG_DIR, "Fig03_sample_lineage_program_heatmap_v1.4.4A.3B.1.pdf"),
  p3, width = 15.0, height = 6.7, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Fig03_sample_lineage_program_heatmap_v1.4.4A.3B.1.png"),
  p3, width = 15.0, height = 6.7, dpi = 300
)

# ------------------------------------------------------------------------------
# 10. Interpretation note
# ------------------------------------------------------------------------------

note <- c(
  "Human MASH v1.4.4A.3B.1 interpretation notes",
  "",
  "Primary measured programs:",
  "- IFN response",
  "- TNF/NF-kB response",
  "- cellular / oxidative stress",
  "- apoptosis / death / survival",
  "- hypoxia / ischemia adaptation",
  "",
  "Exploratory measured program:",
  "- graft-host innate interface",
  "",
  "Antigen presentation:",
  "- Only 4/16 candidate genes are measured in the Flex feature space.",
  "- B2M, classical HLA-A/B/C, multiple HLA-II genes, CIITA, TAP1, PSMB8 and PSMB9 are not measured.",
  "- Therefore antigen presentation is retained as gene-level exploratory evidence only.",
  "- Absence from the Flex feature space is not biological absence.",
  "",
  "Scoring:",
  "- Raw RNA counts are library-normalized to log1p(CP10k).",
  "- Each measured gene is z-standardized across all 43,415 iPSC-derived cells.",
  "- Each program score is the equal-weight mean of its measured gene z-scores.",
  "- No AddModuleScore control-gene background is used because this is a targeted Flex feature space.",
  "",
  "Primary comparison:",
  "- Whole-graft scores are composition-confounded and secondary/descriptive.",
  "- Primary SfLB-vs-Tx comparisons are within Hepatocyte-like and Mesenchymal-like broad lineages,",
  "  the two major lineages retained in both groups after v1.4.4A.3.1.",
  "- Myeloid-like and Unresolved are Tx-emergent and are characterized descriptively rather than",
  "  interpreted as within-lineage increases.",
  "",
  "Causal boundary:",
  "- SfLB -> Tx is NOT an IL10-only perturbation.",
  "- Increased stress/inflammatory programs are called engraftment-associated.",
  "- Human graft transcriptomics alone cannot prove residual xenogeneic rejection.",
  "- Possible residual innate xenogeneic immune pressure must be integrated with mouse-host immune data.",
  "",
  "Statistics:",
  "- n=2 biological samples per group.",
  "- No inferential test is performed.",
  "- Complete replicate-range separation is descriptive evidence, not a p-value."
)

writeLines(
  note,
  file.path(OUT_ROOT, "README_interpretation_v1.4.4A.3B.1.txt"),
  useBytes = TRUE
)

cat("\n============================================================\n")
cat("v1.4.4A.3B.1 COMPLETE\n")
cat("============================================================\n")
cat("Output:", OUT_ROOT, "\n")
cat("Tables:", length(list.files(TAB_DIR)), "\n")
cat("Figures:", length(list.files(FIG_DIR)), "\n")
cat("Primary modules:", length(primary_modules), "\n")
cat("Exploratory aggregate modules:", length(exploratory_modules), "\n")
cat("Antigen presentation: gene-level only (4/16 measured).\n")
cat("No inferential test performed (n=2 per group).\n")
cat("Whole-graft scores are composition-confounded; comparable-lineage results are primary.\n")
