#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ Res2 contamination audit
# v4.11.6
#
# PURPOSE
#   Audit non-MΦ candidate cells identified in v4.11.5 before any exclusion.
#
#   Integrates:
#     - fixed Res2.0 cluster
#     - B / T / NK / neutrophil conservative candidate flags
#     - macrophage identity genes / macrophage identity score
#     - biological sample
#     - fixed mphi.umap.rpca coordinates
#
# IMPORTANT
#   - NO cells are removed.
#   - NO clustering is recalculated.
#   - NO UMAP is recalculated.
#   - v4.8.4 annotation remains unchanged.
#
# MAIN QUESTIONS
#   1) Which Res2 clusters contain candidate non-MΦ cells?
#   2) Do candidate cells form cluster-level islands or scattered cells?
#   3) Are B/T/NK/neutrophil flags overlapping in the same cells?
#   4) Do flagged cells retain a macrophage identity program?
#   5) Are flags sample-specific?
#   6) Which candidates should be considered for v4.11.7 clean-MΦ?
#
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4116)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

PROJECT_DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
    PROJECT_DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_contamination_QC_v4.11.5",
    "Mouse_Mphi_Res2_contamination_QC_annotated_v4.11.5.rds"
)

FALLBACK_RDS <- file.path(
    PROJECT_DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_manual_annotation_v4.8.4",
    "RDS",
    "Mouse_Mphi_Res2_manual_class_annotated_v4.8.4.rds"
)

OUTPUT_DIR <- file.path(
    PROJECT_DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_contamination_audit_v4.11.6"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "dplyr",
    "tidyr",
    "tibble",
    "ggplot2",
    "patchwork",
    "scales"
)

missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(missing_packages) > 0L) {
    stop(
        "Missing package(s): ",
        paste(missing_packages, collapse = ", ")
    )
}

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(patchwork)
    library(scales)
})

# ------------------------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------------------------

msg <- function(...) {
    message(
        "[",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        "] ",
        paste0(...)
    )
}

first_existing <- function(x, candidates) {
    hit <- candidates[candidates %in% x]
    if (length(hit) == 0L) {
        return(NA_character_)
    }
    hit[[1]]
}

get_data_layer <- function(object, assay = "RNA") {
    x <- tryCatch(
        LayerData(
            object,
            assay = assay,
            layer = "data"
        ),
        error = function(e) NULL
    )

    if (is.null(x)) {
        x <- tryCatch(
            GetAssayData(
                object,
                assay = assay,
                slot = "data"
            ),
            error = function(e) NULL
        )
    }

    x
}

mean_gene_score <- function(mat, genes) {
    genes <- intersect(
        genes,
        rownames(mat)
    )

    if (length(genes) == 0L) {
        return(
            rep(
                NA_real_,
                ncol(mat)
            )
        )
    }

    as.numeric(
        Matrix::colMeans(
            mat[
                genes,
                ,
                drop = FALSE
            ]
        )
    )
}

safe_pct <- function(n, d) {
    ifelse(
        d > 0,
        100 * n / d,
        NA_real_
    )
}

save_pdf <- function(filename, plot, width, height) {
    ggsave(
        filename = file.path(FIG_DIR, filename),
        plot = plot,
        device = cairo_pdf,
        width = width,
        height = height,
        units = "in",
        limitsize = FALSE
    )
}

# ------------------------------------------------------------------------------
# 3. Load QC-annotated object
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
    stop(
        "v4.11.5 QC-annotated RDS was not found:\n",
        INPUT_RDS,
        "\nRun v4.11.5 first."
    )
}

msg("Loading v4.11.5 QC object: ", INPUT_RDS)

mphi <- readRDS(INPUT_RDS)

if (!inherits(mphi, "Seurat")) {
    stop("Input object is not a Seurat object.")
}

DefaultAssay(mphi) <- "RNA"

# ------------------------------------------------------------------------------
# 4. Detect metadata columns
# ------------------------------------------------------------------------------

meta_cols <- colnames(mphi@meta.data)

SAMPLE_COL <- first_existing(
    meta_cols,
    c(
        "sample_4group",
        "sample",
        "sample_id",
        "orig.ident"
    )
)

CLASS_COL <- first_existing(
    meta_cols,
    c(
        "macrophage_class_Res2_v484",
        "macrophage_class_v484",
        "manual_class_v484",
        "macrophage_class"
    )
)

CLUSTER_COL <- first_existing(
    meta_cols,
    c(
        "cluster_res2",
        "mphi_rpca_res_2",
        "mphi_rpca_res_2.0",
        "mphi_rpca_res_2.0_cluster",
        "integratedRPCA_snn_res.2",
        "integratedRPCA_snn_res.2.0",
        "seurat_clusters"
    )
)

required_flag_cols <- c(
    "QC_B_candidate",
    "QC_T_candidate",
    "QC_NK_candidate",
    "QC_Neutrophil_candidate",
    "QC_any_nonMphi_candidate"
)

missing_flags <- setdiff(
    required_flag_cols,
    meta_cols
)

if (is.na(SAMPLE_COL)) {
    stop("Biological sample column not found.")
}

if (is.na(CLASS_COL)) {
    stop("v4.8.4 macrophage class column not found.")
}

if (is.na(CLUSTER_COL)) {
    stop(
        "Could not identify Res2.0 cluster column.\nAvailable metadata columns include:\n",
        paste(meta_cols, collapse = "\n")
    )
}

if (length(missing_flags) > 0L) {
    stop(
        "v4.11.5 QC flag column(s) missing: ",
        paste(missing_flags, collapse = ", ")
    )
}

msg("SAMPLE_COL = ", SAMPLE_COL)
msg("CLASS_COL = ", CLASS_COL)
msg("CLUSTER_COL = ", CLUSTER_COL)

# ------------------------------------------------------------------------------
# 5. Fixed RPCA UMAP
# ------------------------------------------------------------------------------

reduction_names <- Reductions(mphi)

candidate_reductions <- c(
    "mphi.umap.rpca",
    "MphiRPCAUMAP",
    "mphiRPCAUMAP",
    "mphi_rpca_umap",
    "umapRPCA",
    "UMAPRPCA",
    "umap.rpca",
    "umap"
)

reduction_hit <- candidate_reductions[
    candidate_reductions %in% reduction_names
]

if (length(reduction_hit) == 0L) {
    stop(
        "Fixed RPCA UMAP reduction not found.\nAvailable reductions: ",
        paste(reduction_names, collapse = ", ")
    )
}

UMAP_REDUCTION <- reduction_hit[[1]]

emb <- Embeddings(
    mphi,
    reduction = UMAP_REDUCTION
)

if (!identical(
    rownames(emb),
    colnames(mphi)
)) {
    stop("UMAP cell IDs do not match Seurat cell order.")
}

msg("Using fixed RPCA UMAP: ", UMAP_REDUCTION)

# ------------------------------------------------------------------------------
# 6. RNA normalized data and MΦ identity genes
# ------------------------------------------------------------------------------

rna_data <- get_data_layer(
    mphi,
    assay = "RNA"
)

if (is.null(rna_data)) {
    stop("RNA normalized data layer not available.")
}

if (!identical(
    colnames(rna_data),
    colnames(mphi)
)) {
    stop("RNA data cell order does not match Seurat object.")
}

# Core macrophage identity.
# Chosen to avoid depending only on inflammatory / lipid / M2 state.
MPHI_IDENTITY_GENES <- c(
    "Adgre1",
    "Csf1r",
    "Fcgr1",
    "Cd68",
    "Lyz2",
    "C1qa",
    "C1qb",
    "C1qc",
    "Aif1",
    "Tyrobp",
    "Ctss",
    "Fcer1g"
)

MPHI_IDENTITY_USE <- intersect(
    MPHI_IDENTITY_GENES,
    rownames(rna_data)
)

if (length(MPHI_IDENTITY_USE) < 5L) {
    stop(
        "Too few macrophage identity genes detected: ",
        paste(MPHI_IDENTITY_USE, collapse = ", ")
    )
}

msg(
    "Macrophage identity genes detected: ",
    paste(MPHI_IDENTITY_USE, collapse = ", ")
)

mphi$QC_Mphi_identity_score <- mean_gene_score(
    rna_data,
    MPHI_IDENTITY_USE
)

# Number of macrophage identity genes expressed per cell.
mphi_identity_mat <- rna_data[
    MPHI_IDENTITY_USE,
    ,
    drop = FALSE
]

mphi$QC_Mphi_identity_n <- as.integer(
    Matrix::colSums(
        mphi_identity_mat > 0
    )
)

# ------------------------------------------------------------------------------
# 7. Canonical audit metadata
# ------------------------------------------------------------------------------

audit_df <- mphi@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = as.character(.data[[SAMPLE_COL]]),
        macrophage_class = as.character(.data[[CLASS_COL]]),
        cluster = as.character(.data[[CLUSTER_COL]]),
        B_candidate = QC_B_candidate,
        T_candidate = QC_T_candidate,
        NK_candidate = QC_NK_candidate,
        Neutrophil_candidate = QC_Neutrophil_candidate,
        any_nonMphi_candidate = QC_any_nonMphi_candidate,
        Mphi_identity_score = QC_Mphi_identity_score,
        Mphi_identity_n = QC_Mphi_identity_n
    )

# Numeric-like cluster sort.
cluster_levels <- unique(audit_df$cluster)

suppressWarnings({
    cluster_num <- as.numeric(cluster_levels)
})

if (all(!is.na(cluster_num))) {
    cluster_levels <- cluster_levels[
        order(cluster_num)
    ]
} else {
    cluster_levels <- sort(cluster_levels)
}

audit_df$cluster <- factor(
    audit_df$cluster,
    levels = cluster_levels
)

# ------------------------------------------------------------------------------
# 8. Flag overlap / doublet-like pattern
# ------------------------------------------------------------------------------

audit_df <- audit_df %>%
    mutate(
        n_lineage_flags =
            as.integer(B_candidate) +
            as.integer(T_candidate) +
            as.integer(NK_candidate) +
            as.integer(Neutrophil_candidate),

        flag_pattern = case_when(
            n_lineage_flags == 0L ~ "None",
            n_lineage_flags == 1L & B_candidate ~ "B only",
            n_lineage_flags == 1L & T_candidate ~ "T only",
            n_lineage_flags == 1L & NK_candidate ~ "NK only",
            n_lineage_flags == 1L & Neutrophil_candidate ~ "Neutrophil only",
            n_lineage_flags >= 2L ~ "Multi-lineage",
            TRUE ~ "Other"
        )
    )

# Push back to Seurat metadata.
mphi$QC_n_lineage_flags_v4116 <- audit_df$n_lineage_flags
mphi$QC_flag_pattern_v4116 <- audit_df$flag_pattern

write.csv(
    audit_df,
    file.path(
        TAB_DIR,
        "01_cell_level_contamination_audit_v4.11.6.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Cluster-level contamination burden
# ------------------------------------------------------------------------------

cluster_summary <- audit_df %>%
    group_by(cluster) %>%
    summarise(
        total_cells = n(),

        B_n = sum(B_candidate),
        T_n = sum(T_candidate),
        NK_n = sum(NK_candidate),
        Neutrophil_n = sum(Neutrophil_candidate),
        Any_n = sum(any_nonMphi_candidate),
        Multi_lineage_n = sum(n_lineage_flags >= 2L),

        B_percent = safe_pct(B_n, total_cells),
        T_percent = safe_pct(T_n, total_cells),
        NK_percent = safe_pct(NK_n, total_cells),
        Neutrophil_percent = safe_pct(Neutrophil_n, total_cells),
        Any_percent = safe_pct(Any_n, total_cells),
        Multi_lineage_percent = safe_pct(
            Multi_lineage_n,
            total_cells
        ),

        median_Mphi_identity =
            median(
                Mphi_identity_score,
                na.rm = TRUE
            ),

        median_Mphi_identity_n =
            median(
                Mphi_identity_n,
                na.rm = TRUE
            ),

        .groups = "drop"
    )

write.csv(
    cluster_summary,
    file.path(
        TAB_DIR,
        "02_contamination_burden_by_Res2_cluster_v4.11.6.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Cluster x lineage long table
# ------------------------------------------------------------------------------

cluster_lineage_long <- cluster_summary %>%
    select(
        cluster,
        total_cells,
        B_percent,
        T_percent,
        NK_percent,
        Neutrophil_percent
    ) %>%
    pivot_longer(
        cols = c(
            B_percent,
            T_percent,
            NK_percent,
            Neutrophil_percent
        ),
        names_to = "lineage",
        values_to = "candidate_percent"
    ) %>%
    mutate(
        lineage = recode(
            lineage,
            B_percent = "B",
            T_percent = "T",
            NK_percent = "NK",
            Neutrophil_percent = "Neutrophil"
        )
    )

# ------------------------------------------------------------------------------
# 11. Cluster-level contamination heatmap-style plot
# ------------------------------------------------------------------------------

p_cluster_lineage <- ggplot(
    cluster_lineage_long,
    aes(
        x = cluster,
        y = lineage,
        fill = candidate_percent
    )
) +
    geom_tile() +
    geom_text(
        aes(
            label = sprintf(
                "%.1f",
                candidate_percent
            )
        ),
        size = 3
    ) +
    scale_fill_gradient(
        low = "white",
        high = "black",
        name = "% candidate"
    ) +
    labs(
        title = "Lineage-candidate burden by Res2.0 cluster",
        subtitle = "% of cells within each cluster flagged by conservative v4.11.5 criteria",
        x = "Res2 cluster",
        y = NULL
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        axis.text.x = element_text(
            angle = 0
        )
    )

save_pdf(
    "01_candidate_burden_by_Res2_cluster_v4.11.6.pdf",
    p_cluster_lineage,
    14,
    5.5
)

# ------------------------------------------------------------------------------
# 12. Any-candidate + macrophage identity by cluster
# ------------------------------------------------------------------------------

cluster_identity_long <- audit_df %>%
    group_by(
        cluster,
        any_nonMphi_candidate
    ) %>%
    summarise(
        n_cells = n(),
        mean_Mphi_identity =
            mean(
                Mphi_identity_score,
                na.rm = TRUE
            ),
        median_Mphi_identity =
            median(
                Mphi_identity_score,
                na.rm = TRUE
            ),
        mean_Mphi_identity_n =
            mean(
                Mphi_identity_n,
                na.rm = TRUE
            ),
        .groups = "drop"
    )

write.csv(
    cluster_identity_long,
    file.path(
        TAB_DIR,
        "03_Mphi_identity_candidate_vs_non_candidate_by_cluster_v4.11.6.csv"
    ),
    row.names = FALSE
)

p_identity_cluster <- ggplot(
    audit_df,
    aes(
        x = cluster,
        y = Mphi_identity_score,
        fill = any_nonMphi_candidate
    )
) +
    geom_boxplot(
        outlier.shape = NA,
        width = 0.72
    ) +
    scale_fill_manual(
        values = c(
            "FALSE" = "grey85",
            "TRUE" = "grey35"
        ),
        labels = c(
            "FALSE" = "Not flagged",
            "TRUE" = "Candidate non-MΦ"
        ),
        name = NULL
    ) +
    labs(
        title = "Macrophage identity retained in candidate cells",
        subtitle = "Mean normalized expression of core MΦ identity genes",
        x = "Res2 cluster",
        y = "MΦ identity score"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        legend.position = "top"
    )

save_pdf(
    "02_Mphi_identity_by_cluster_candidate_status_v4.11.6.pdf",
    p_identity_cluster,
    14,
    7
)

# ------------------------------------------------------------------------------
# 13. UMAP: cluster number + candidate overlay
# ------------------------------------------------------------------------------

umap_df <- tibble(
    cell = rownames(emb),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2]
) %>%
    left_join(
        audit_df,
        by = "cell"
    )

cluster_centers <- umap_df %>%
    group_by(cluster) %>%
    summarise(
        UMAP_1 = median(
            UMAP_1,
            na.rm = TRUE
        ),
        UMAP_2 = median(
            UMAP_2,
            na.rm = TRUE
        ),
        .groups = "drop"
    )

p_cluster_candidate <- ggplot(
    umap_df,
    aes(
        UMAP_1,
        UMAP_2
    )
) +
    geom_point(
        color = "grey87",
        size = 0.28,
        alpha = 0.55
    ) +
    geom_point(
        data = umap_df %>%
            filter(
                any_nonMphi_candidate
            ),
        color = "#FF1A1A",
        size = 0.85,
        alpha = 0.95
    ) +
    geom_label(
        data = cluster_centers,
        aes(
            label = cluster
        ),
        size = 3.7,
        fontface = "bold",
        label.size = 0.20,
        fill = "white"
    ) +
    coord_equal() +
    labs(
        title = "Candidate non-MΦ cells by fixed Res2.0 cluster",
        subtitle = "Red = conservative v4.11.5 candidate; labels = Res2 cluster",
        x = paste0(
            UMAP_REDUCTION,
            "_1"
        ),
        y = paste0(
            UMAP_REDUCTION,
            "_2"
        )
    ) +
    theme_classic(
        base_size = 11
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        )
    )

save_pdf(
    "03_candidate_cells_with_Res2_cluster_labels_v4.11.6.pdf",
    p_cluster_candidate,
    9,
    7.5
)

# ------------------------------------------------------------------------------
# 14. UMAP: flag pattern
# ------------------------------------------------------------------------------

pattern_order <- c(
    "None",
    "B only",
    "T only",
    "NK only",
    "Neutrophil only",
    "Multi-lineage",
    "Other"
)

umap_df$flag_pattern <- factor(
    umap_df$flag_pattern,
    levels = pattern_order
)

p_flag_pattern <- ggplot(
    umap_df,
    aes(
        UMAP_1,
        UMAP_2
    )
) +
    geom_point(
        data = umap_df %>%
            filter(
                flag_pattern == "None"
            ),
        color = "grey88",
        size = 0.26,
        alpha = 0.40
    ) +
    geom_point(
        data = umap_df %>%
            filter(
                flag_pattern != "None"
            ),
        aes(
            shape = flag_pattern
        ),
        size = 1.05,
        alpha = 0.95
    ) +
    coord_equal() +
    labs(
        title = "Lineage flag pattern on fixed RPCA UMAP",
        subtitle = "Multi-lineage cells are especially suspicious for doublet/multiplet contamination",
        x = NULL,
        y = NULL,
        shape = "QC pattern"
    ) +
    theme_classic(
        base_size = 11
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        )
    )

save_pdf(
    "04_lineage_flag_pattern_fixed_RPCA_UMAP_v4.11.6.pdf",
    p_flag_pattern,
    9,
    7.5
)

# ------------------------------------------------------------------------------
# 15. Flag-overlap table
# ------------------------------------------------------------------------------

flag_overlap <- audit_df %>%
    count(
        B_candidate,
        T_candidate,
        NK_candidate,
        Neutrophil_candidate,
        flag_pattern,
        name = "n_cells"
    ) %>%
    mutate(
        percent_total =
            100 *
            n_cells /
            nrow(audit_df)
    ) %>%
    arrange(
        desc(n_cells)
    )

write.csv(
    flag_overlap,
    file.path(
        TAB_DIR,
        "04_lineage_flag_overlap_v4.11.6.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 16. Candidate burden by sample
# ------------------------------------------------------------------------------

sample_summary <- audit_df %>%
    group_by(sample) %>%
    summarise(
        total_cells = n(),
        Any_n = sum(any_nonMphi_candidate),
        Multi_lineage_n = sum(n_lineage_flags >= 2L),
        B_n = sum(B_candidate),
        T_n = sum(T_candidate),
        NK_n = sum(NK_candidate),
        Neutrophil_n = sum(Neutrophil_candidate),
        Any_percent = safe_pct(
            Any_n,
            total_cells
        ),
        Multi_lineage_percent = safe_pct(
            Multi_lineage_n,
            total_cells
        ),
        .groups = "drop"
    )

write.csv(
    sample_summary,
    file.path(
        TAB_DIR,
        "05_contamination_burden_by_sample_v4.11.6.csv"
    ),
    row.names = FALSE
)

sample_long <- sample_summary %>%
    select(
        sample,
        B_n,
        T_n,
        NK_n,
        Neutrophil_n
    ) %>%
    pivot_longer(
        cols = c(
            B_n,
            T_n,
            NK_n,
            Neutrophil_n
        ),
        names_to = "lineage",
        values_to = "n_candidate"
    ) %>%
    mutate(
        lineage = recode(
            lineage,
            B_n = "B",
            T_n = "T",
            NK_n = "NK",
            Neutrophil_n = "Neutrophil"
        )
    )

p_sample <- ggplot(
    sample_long,
    aes(
        x = sample,
        y = n_candidate,
        fill = lineage
    )
) +
    geom_col(
        position = "dodge"
    ) +
    labs(
        title = "Candidate non-MΦ cells by biological sample",
        x = NULL,
        y = "Candidate cell count",
        fill = "Lineage"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        axis.text.x = element_text(
            angle = 35,
            hjust = 1
        )
    )

save_pdf(
    "05_candidate_counts_by_sample_v4.11.6.pdf",
    p_sample,
    10,
    6
)

# ------------------------------------------------------------------------------
# 17. Sample x cluster candidate burden
# ------------------------------------------------------------------------------

sample_cluster_summary <- audit_df %>%
    group_by(
        sample,
        cluster
    ) %>%
    summarise(
        total_cells = n(),
        Any_n = sum(any_nonMphi_candidate),
        B_n = sum(B_candidate),
        T_n = sum(T_candidate),
        NK_n = sum(NK_candidate),
        Neutrophil_n = sum(Neutrophil_candidate),
        Any_percent = safe_pct(
            Any_n,
            total_cells
        ),
        .groups = "drop"
    )

write.csv(
    sample_cluster_summary,
    file.path(
        TAB_DIR,
        "06_contamination_burden_sample_x_cluster_v4.11.6.csv"
    ),
    row.names = FALSE
)

p_sample_cluster <- ggplot(
    sample_cluster_summary,
    aes(
        x = cluster,
        y = sample,
        fill = Any_percent
    )
) +
    geom_tile() +
    geom_text(
        aes(
            label = sprintf(
                "%.0f",
                Any_percent
            )
        ),
        size = 2.8
    ) +
    scale_fill_gradient(
        low = "white",
        high = "black",
        name = "% candidate"
    ) +
    labs(
        title = "Candidate burden by sample × Res2 cluster",
        x = "Res2 cluster",
        y = NULL
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        )
    )

save_pdf(
    "06_candidate_percent_sample_x_cluster_v4.11.6.pdf",
    p_sample_cluster,
    14,
    5.8
)

# ------------------------------------------------------------------------------
# 18. Macrophage identity by flag pattern
# ------------------------------------------------------------------------------

p_identity_pattern <- ggplot(
    audit_df,
    aes(
        x = flag_pattern,
        y = Mphi_identity_score
    )
) +
    geom_violin(
        scale = "width",
        trim = TRUE
    ) +
    geom_boxplot(
        width = 0.12,
        outlier.shape = NA
    ) +
    labs(
        title = "Macrophage identity by contamination flag pattern",
        subtitle = "Low MΦ identity + multi-lineage signal strengthens a non-MΦ/doublet interpretation",
        x = NULL,
        y = "MΦ identity score"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        axis.text.x = element_text(
            angle = 35,
            hjust = 1
        )
    )

save_pdf(
    "07_Mphi_identity_by_flag_pattern_v4.11.6.pdf",
    p_identity_pattern,
    10,
    6
)

# ------------------------------------------------------------------------------
# 19. Candidate-cell detailed audit table
# ------------------------------------------------------------------------------

candidate_detailed <- audit_df %>%
    filter(
        any_nonMphi_candidate
    ) %>%
    arrange(
        desc(n_lineage_flags),
        Mphi_identity_score,
        cluster,
        sample
    )

write.csv(
    candidate_detailed,
    file.path(
        TAB_DIR,
        "07_candidate_cell_detailed_audit_v4.11.6.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 20. Preliminary exclusion-risk classification
#
# IMPORTANT:
# This is an AUDIT category, not automatic filtering.
#
# Tier A:
#   Multi-lineage flag + low macrophage identity
#
# Tier B:
#   Single-lineage flag + low macrophage identity
#
# Tier C:
#   Candidate flag but macrophage identity retained
#
# Low MΦ identity is defined relative to candidate/noncandidate distribution
# using the 10th percentile of all MΦ identity scores.
# ------------------------------------------------------------------------------

MPHI_IDENTITY_LOW_CUTOFF <- quantile(
    audit_df$Mphi_identity_score,
    probs = 0.10,
    na.rm = TRUE
)

audit_df <- audit_df %>%
    mutate(
        Mphi_identity_low =
            Mphi_identity_score <=
            MPHI_IDENTITY_LOW_CUTOFF,

        exclusion_audit_tier = case_when(
            !any_nonMphi_candidate ~ "Not flagged",

            n_lineage_flags >= 2L &
                Mphi_identity_low ~
                "Tier A: high suspicion",

            n_lineage_flags == 1L &
                Mphi_identity_low ~
                "Tier B: moderate suspicion",

            any_nonMphi_candidate &
                !Mphi_identity_low ~
                "Tier C: MΦ identity retained",

            TRUE ~
                "Unclassified"
        )
    )

tier_summary <- audit_df %>%
    count(
        exclusion_audit_tier,
        name = "n_cells"
    ) %>%
    mutate(
        percent_total =
            100 *
            n_cells /
            nrow(audit_df)
    )

write.csv(
    tier_summary,
    file.path(
        TAB_DIR,
        "08_preliminary_exclusion_audit_tier_summary_v4.11.6.csv"
    ),
    row.names = FALSE
)

write.csv(
    audit_df,
    file.path(
        TAB_DIR,
        "09_all_cells_with_preliminary_exclusion_audit_tier_v4.11.6.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 21. UMAP audit tier
# ------------------------------------------------------------------------------

umap_tier_df <- tibble(
    cell = rownames(emb),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2]
) %>%
    left_join(
        audit_df %>%
            select(
                cell,
                exclusion_audit_tier
            ),
        by = "cell"
    )

tier_order <- c(
    "Not flagged",
    "Tier C: MΦ identity retained",
    "Tier B: moderate suspicion",
    "Tier A: high suspicion"
)

umap_tier_df$exclusion_audit_tier <- factor(
    umap_tier_df$exclusion_audit_tier,
    levels = tier_order
)

p_tier <- ggplot(
    umap_tier_df,
    aes(
        UMAP_1,
        UMAP_2
    )
) +
    geom_point(
        data = umap_tier_df %>%
            filter(
                exclusion_audit_tier ==
                    "Not flagged"
            ),
        color = "grey88",
        size = 0.25,
        alpha = 0.40
    ) +
    geom_point(
        data = umap_tier_df %>%
            filter(
                exclusion_audit_tier !=
                    "Not flagged"
            ),
        aes(
            shape = exclusion_audit_tier
        ),
        size = 1.05,
        alpha = 0.95
    ) +
    coord_equal() +
    labs(
        title = "Preliminary contamination audit tiers",
        subtitle = "Audit only — no cells removed",
        x = NULL,
        y = NULL,
        shape = "Audit tier"
    ) +
    theme_classic(
        base_size = 11
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        )
    )

save_pdf(
    "08_preliminary_exclusion_audit_tiers_fixed_RPCA_UMAP_v4.11.6.pdf",
    p_tier,
    9,
    7.5
)

# ------------------------------------------------------------------------------
# 22. Cluster-level audit summary with tier burden
# ------------------------------------------------------------------------------

cluster_tier_summary <- audit_df %>%
    group_by(
        cluster
    ) %>%
    summarise(
        total_cells = n(),

        TierA_n =
            sum(
                exclusion_audit_tier ==
                    "Tier A: high suspicion"
            ),

        TierB_n =
            sum(
                exclusion_audit_tier ==
                    "Tier B: moderate suspicion"
            ),

        TierC_n =
            sum(
                exclusion_audit_tier ==
                    "Tier C: MΦ identity retained"
            ),

        TierA_percent =
            safe_pct(
                TierA_n,
                total_cells
            ),

        TierB_percent =
            safe_pct(
                TierB_n,
                total_cells
            ),

        TierC_percent =
            safe_pct(
                TierC_n,
                total_cells
            ),

        .groups = "drop"
    )

write.csv(
    cluster_tier_summary,
    file.path(
        TAB_DIR,
        "10_audit_tier_burden_by_Res2_cluster_v4.11.6.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 23. Summary figure
# ------------------------------------------------------------------------------

p_summary <- (
    p_cluster_candidate /
    p_cluster_lineage /
    p_identity_cluster /
    p_sample_cluster /
    p_identity_pattern /
    p_tier
) +
    patchwork::plot_layout(
        heights = c(
            1.15,
            0.75,
            0.90,
            0.75,
            0.80,
            1.05
        )
    ) +
    patchwork::plot_annotation(
        title = "Mouse MASH MΦ contamination audit v4.11.6",
        subtitle = "Res2 cluster × lineage candidate × MΦ identity × biological sample",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 18
            )
        )
    )

save_pdf(
    "09_contamination_audit_summary_v4.11.6.pdf",
    p_summary,
    16,
    36
)

# ------------------------------------------------------------------------------
# 24. Save audit-annotated object
# ------------------------------------------------------------------------------

# Add tier back to Seurat metadata in original cell order.
tier_map <- audit_df$exclusion_audit_tier
names(tier_map) <- audit_df$cell

mphi$QC_exclusion_audit_tier_v4116 <- tier_map[
    colnames(mphi)
]

mphi$QC_Mphi_identity_score_v4116 <- mphi$QC_Mphi_identity_score
mphi$QC_Mphi_identity_n_v4116 <- mphi$QC_Mphi_identity_n

AUDIT_RDS <- file.path(
    OUTPUT_DIR,
    "Mouse_Mphi_Res2_contamination_audit_annotated_v4.11.6.rds"
)

saveRDS(
    mphi,
    AUDIT_RDS
)

# ------------------------------------------------------------------------------
# 25. README
# ------------------------------------------------------------------------------

readme <- c(
    "Mouse MASH MΦ Res2 contamination audit v4.11.6",
    "",
    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",
    "No cells were removed.",
    "No clustering or UMAP was recalculated.",
    "",
    paste0(
        "Fixed RPCA UMAP: ",
        UMAP_REDUCTION
    ),
    paste0(
        "Res2 cluster metadata: ",
        CLUSTER_COL
    ),
    "",
    "Audit dimensions:",
    "  B / T / NK / neutrophil conservative candidate flags",
    "  lineage-flag overlap",
    "  macrophage identity score",
    "  Res2 cluster",
    "  biological sample",
    "",
    "Macrophage identity genes:",
    paste0(
        "  ",
        paste(
            MPHI_IDENTITY_USE,
            collapse = ", "
        )
    ),
    "",
    paste0(
        "Preliminary low-MΦ-identity cutoff (10th percentile): ",
        signif(
            MPHI_IDENTITY_LOW_CUTOFF,
            5
        )
    ),
    "",
    "Preliminary audit tiers:",
    "  Tier A = multi-lineage candidate + low MΦ identity",
    "  Tier B = single-lineage candidate + low MΦ identity",
    "  Tier C = candidate flag but MΦ identity retained",
    "",
    "These tiers are NOT automatic exclusion criteria.",
    "Use them to define the clean-MΦ rule in the next version.",
    "",
    "Primary outputs:",
    "  01 candidate burden by Res2 cluster",
    "  02 MΦ identity by cluster and candidate status",
    "  03 candidate cells with Res2 labels",
    "  04 lineage flag pattern UMAP",
    "  05 candidate counts by sample",
    "  06 sample x cluster candidate burden",
    "  07 MΦ identity by flag pattern",
    "  08 preliminary audit tier UMAP",
    "  09 summary",
    "",
    "Next step:",
    "  Review Tier A/B/C distribution and decide v4.11.7 clean-MΦ rule."
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.11.6.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.11.6.txt"
    )
)

# ------------------------------------------------------------------------------
# 26. Final
# ------------------------------------------------------------------------------

msg("DONE.")
msg("Output: ", OUTPUT_DIR)

print(
    tier_summary
)

print(
    cluster_summary %>%
        arrange(
            desc(Any_percent)
        )
)

