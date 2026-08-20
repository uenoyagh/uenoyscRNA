#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ Res2 contamination QC
# v4.11.5
#
# Purpose:
#   QC ONLY.
#   Examine possible B / T / NK / neutrophil contamination or doublets
#   in the fixed v4.8.4 MΦ-only object.
#
# IMPORTANT:
#   - Res2.0 clustering is NOT changed.
#   - v4.8.4 macrophage annotation is NOT changed.
#   - No cells are removed.
#   - No PCA / RPCA / UMAP is recalculated.
#   - Existing MΦ-only RPCA UMAP coordinates are used.
#
# v4.11.5 FIX:
#   - AddModuleScore() is NOT used for contamination scores.
#   - Each contamination score is calculated directly as the mean normalized
#     RNA expression of the corresponding lineage-marker set in each cell.
#
# Output:
#   1. contamination score UMAPs
#   2. individual marker UMAPs
#   3. contamination score by MΦ subtype
#   4. contamination score by biological sample
#   5. candidate contaminating-cell table
#   6. contamination burden summary
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

INPUT_RDS <- paste0(
    "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
    "Mouse_MASH_Mphi_RDS/",
    "Mphi_Res2_manual_annotation_v4.8.4/",
    "RDS/",
    "Mouse_Mphi_Res2_manual_class_annotated_v4.8.4.rds"
)

OUTPUT_DIR <- paste0(
    "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
    "Mouse_MASH_Mphi_RDS/",
    "Mphi_Res2_contamination_QC_v4.11.5"
)

FIG_DIR <- file.path(
    OUTPUT_DIR,
    "Figures"
)

TAB_DIR <- file.path(
    OUTPUT_DIR,
    "Tables"
)

LOG_DIR <- file.path(
    OUTPUT_DIR,
    "Logs"
)

dir.create(
    FIG_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    TAB_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    LOG_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


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
        paste(
            missing_packages,
            collapse = ", "
        )
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
# 2. Utility
# ------------------------------------------------------------------------------

msg <- function(...) {

    message(
        "[",
        format(
            Sys.time(),
            "%Y-%m-%d %H:%M:%S"
        ),
        "] ",
        paste0(...)
    )
}


# ------------------------------------------------------------------------------
# 3. Load fixed v4.8.4 object
# ------------------------------------------------------------------------------

msg(
    "Loading fixed v4.8.4 RDS..."
)

if (!file.exists(INPUT_RDS)) {
    stop(
        "Input RDS not found: ",
        INPUT_RDS
    )
}

mphi <- readRDS(
    INPUT_RDS
)

if (!inherits(mphi, "Seurat")) {
    stop(
        "Input object is not a Seurat object."
    )
}

msg(
    "Cells: ",
    ncol(mphi)
)

msg(
    "Genes: ",
    nrow(mphi)
)

if (!"RNA" %in% Assays(mphi)) {

    stop(
        "RNA assay not found."
    )
}

DefaultAssay(mphi) <- "RNA"


# ------------------------------------------------------------------------------
# 4. Verify fixed metadata
# ------------------------------------------------------------------------------

required_meta <- c(
    "sample_4group",
    "macrophage_class_Res2_v484"
)

missing_meta <- setdiff(
    required_meta,
    colnames(mphi@meta.data)
)

if (length(missing_meta) > 0L) {

    stop(
        "Required metadata missing: ",
        paste(
            missing_meta,
            collapse = ", "
        )
    )
}

msg(
    "Biological samples: ",
    paste(
        unique(mphi$sample_4group),
        collapse = ", "
    )
)

msg(
    "Fixed MΦ classes: ",
    paste(
        unique(mphi$macrophage_class_Res2_v484),
        collapse = ", "
    )
)


# ------------------------------------------------------------------------------
# 5. Detect fixed RPCA UMAP reduction
#
# Priority:
#   mphi.umap.rpca = v4.8.x reference coordinate system
# ------------------------------------------------------------------------------

reduction_names <- Reductions(
    mphi
)

msg(
    "Available reductions: ",
    paste(
        reduction_names,
        collapse = ", "
    )
)

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
    candidate_reductions %in%
        reduction_names
]

if (length(reduction_hit) == 0L) {

    stop(
        paste0(
            "Could not identify the fixed RPCA UMAP reduction.\n",
            "Available reductions: ",
            paste(
                reduction_names,
                collapse = ", "
            )
        )
    )
}

UMAP_REDUCTION <- reduction_hit[[1]]

msg(
    "Using existing UMAP reduction: ",
    UMAP_REDUCTION
)


# ------------------------------------------------------------------------------
# 6. Verify existing UMAP coordinates
# ------------------------------------------------------------------------------

umap_embeddings <- Embeddings(
    mphi,
    reduction = UMAP_REDUCTION
)

if (ncol(umap_embeddings) < 2L) {

    stop(
        "Selected reduction has fewer than two dimensions."
    )
}

if (!identical(
    rownames(umap_embeddings),
    colnames(mphi)
)) {

    stop(
        "UMAP cell IDs do not match Seurat cell IDs."
    )
}

msg(
    "Existing RPCA UMAP coordinates verified."
)

write.csv(
    data.frame(
        cell = rownames(umap_embeddings),
        UMAP_1 = umap_embeddings[, 1],
        UMAP_2 = umap_embeddings[, 2],
        check.names = FALSE
    ),
    file.path(
        TAB_DIR,
        "00_fixed_RPCA_UMAP_coordinates_v4.11.5.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7. Define contamination marker programs
# ------------------------------------------------------------------------------

CONTAMINATION_PROGRAMS <- list(

    B_cell = c(
        "Cd19",
        "Ms4a1",
        "Cd79a",
        "Cd79b",
        "Cd37",
        "Cd22",
        "Cd74",
        "H2-Aa",
        "H2-Ab1",
        "Ebf1",
        "Pax5",
        "Pou2af1",
        "Bank1",
        "Fcrl1",
        "Fcrl5",
        "Igkc",
        "Ighm"
    ),

    T_cell = c(
        "Cd3d",
        "Cd3e",
        "Cd3g",
        "Trac",
        "Trbc1",
        "Trbc2",
        "Cd247",
        "Lck",
        "Lat",
        "Cd2"
    ),

    NK_cell = c(
        "Nkg7",
        "Klrd1",
        "Klrk1",
        "Klrb1c",
        "Ncr1",
        "Prf1",
        "Gzma",
        "Gzmb",
        "Ccl5",
        "Xcl1",
        "Xcl2"
    ),

    Neutrophil = c(
        "S100a8",
        "S100a9",
        "Retnlg",
        "Ltf",
        "Mpo",
        "Elane",
        "Camp",
        "Ngp",
        "Ly6g",
        "Csf3r",
        "Mmp8",
        "Mmp9"
    )
)


# ------------------------------------------------------------------------------
# 8. Detect available genes
# ------------------------------------------------------------------------------

available_genes <- rownames(
    mphi[["RNA"]]
)

PROGRAM_USE <- lapply(
    CONTAMINATION_PROGRAMS,
    function(x) {

        intersect(
            x,
            available_genes
        )
    }
)

program_detection <- bind_rows(
    lapply(
        names(CONTAMINATION_PROGRAMS),
        function(program_name) {

            expected <- CONTAMINATION_PROGRAMS[[program_name]]
            detected <- PROGRAM_USE[[program_name]]

            tibble(
                program = program_name,
                n_expected = length(expected),
                n_detected = length(detected),
                detected_genes = paste(
                    detected,
                    collapse = ";"
                ),
                missing_genes = paste(
                    setdiff(
                        expected,
                        detected
                    ),
                    collapse = ";"
                )
            )
        }
    )
)

write.csv(
    program_detection,
    file.path(
        TAB_DIR,
        "01_contamination_marker_detection_v4.11.5.csv"
    ),
    row.names = FALSE
)

print(
    program_detection
)


# ------------------------------------------------------------------------------
# 9. Calculate contamination scores directly from normalized RNA data
#
# v4.11.5 FIX:
#   AddModuleScore() is NOT used.
#
# For each cell:
#   contamination score =
#       mean normalized expression of detected lineage-specific marker genes
#
# This is used ONLY as a QC visualization metric.
# It does NOT alter clustering or annotation.
# ------------------------------------------------------------------------------

rna_data <- tryCatch(
    LayerData(
        mphi,
        assay = "RNA",
        layer = "data"
    ),
    error = function(e) NULL
)

if (
    is.null(rna_data) ||
    nrow(rna_data) == 0L
) {

    rna_data <- tryCatch(
        GetAssayData(
            mphi,
            assay = "RNA",
            slot = "data"
        ),
        error = function(e) NULL
    )
}

if (
    is.null(rna_data) ||
    nrow(rna_data) == 0L
) {

    stop(
        "Normalized RNA data layer could not be retrieved."
    )
}

if (!identical(
    colnames(rna_data),
    colnames(mphi)
)) {

    stop(
        "Normalized RNA data cell order does not match Seurat object."
    )
}


# ------------------------------------------------------------------------------
# 9A. Helper: mean normalized expression per cell
# ------------------------------------------------------------------------------

calculate_lineage_score <- function(
    marker_genes,
    score_name
) {

    marker_genes <- intersect(
        marker_genes,
        rownames(rna_data)
    )

    if (length(marker_genes) < 2L) {

        stop(
            score_name,
            ": fewer than 2 marker genes detected. Detected genes = ",
            paste(
                marker_genes,
                collapse = ", "
            )
        )
    }

    message(
        score_name,
        ": using ",
        length(marker_genes),
        " genes: ",
        paste(
            marker_genes,
            collapse = ", "
        )
    )

    score <- Matrix::colMeans(
        rna_data[
            marker_genes,
            ,
            drop = FALSE
        ]
    )

    as.numeric(score)
}


# ------------------------------------------------------------------------------
# 9B. Calculate four contamination scores
# ------------------------------------------------------------------------------

mphi$QC_B_cell_score <- calculate_lineage_score(
    PROGRAM_USE[["B_cell"]],
    "B-cell QC"
)

mphi$QC_T_cell_score <- calculate_lineage_score(
    PROGRAM_USE[["T_cell"]],
    "T-cell QC"
)

mphi$QC_NK_cell_score <- calculate_lineage_score(
    PROGRAM_USE[["NK_cell"]],
    "NK-cell QC"
)

mphi$QC_Neutrophil_score <- calculate_lineage_score(
    PROGRAM_USE[["Neutrophil"]],
    "Neutrophil QC"
)


# ------------------------------------------------------------------------------
# 10. Score columns and sanity checks
# ------------------------------------------------------------------------------

SCORE_COLS <- c(
    B_cell = "QC_B_cell_score",
    T_cell = "QC_T_cell_score",
    NK_cell = "QC_NK_cell_score",
    Neutrophil = "QC_Neutrophil_score"
)

missing_scores <- setdiff(
    unname(SCORE_COLS),
    colnames(mphi@meta.data)
)

if (length(missing_scores) > 0L) {

    stop(
        "Missing contamination score column(s): ",
        paste(
            missing_scores,
            collapse = ", "
        )
    )
}

score_summary <- mphi@meta.data %>%
    summarise(
        B_min = min(
            QC_B_cell_score,
            na.rm = TRUE
        ),
        B_median = median(
            QC_B_cell_score,
            na.rm = TRUE
        ),
        B_max = max(
            QC_B_cell_score,
            na.rm = TRUE
        ),

        T_min = min(
            QC_T_cell_score,
            na.rm = TRUE
        ),
        T_median = median(
            QC_T_cell_score,
            na.rm = TRUE
        ),
        T_max = max(
            QC_T_cell_score,
            na.rm = TRUE
        ),

        NK_min = min(
            QC_NK_cell_score,
            na.rm = TRUE
        ),
        NK_median = median(
            QC_NK_cell_score,
            na.rm = TRUE
        ),
        NK_max = max(
            QC_NK_cell_score,
            na.rm = TRUE
        ),

        Neutrophil_min = min(
            QC_Neutrophil_score,
            na.rm = TRUE
        ),
        Neutrophil_median = median(
            QC_Neutrophil_score,
            na.rm = TRUE
        ),
        Neutrophil_max = max(
            QC_Neutrophil_score,
            na.rm = TRUE
        )
    )

print(
    score_summary
)

write.csv(
    score_summary,
    file.path(
        TAB_DIR,
        "01b_contamination_score_summary_v4.11.5.csv"
    ),
    row.names = FALSE
)

message(
    "Contamination score columns successfully generated:"
)

print(
    SCORE_COLS
)


# ------------------------------------------------------------------------------
# 11. Fixed UMAP plotting function
# ------------------------------------------------------------------------------

plot_score_umap <- function(
    object,
    score_col,
    title_text
) {

    emb <- Embeddings(
        object,
        reduction = UMAP_REDUCTION
    )

    plot_df <- tibble(
        cell = rownames(emb),
        UMAP_1 = emb[, 1],
        UMAP_2 = emb[, 2],
        score = object@meta.data[
            rownames(emb),
            score_col
        ]
    )

    score_limits <- quantile(
        plot_df$score,
        probs = c(
            0.02,
            0.98
        ),
        na.rm = TRUE
    )

    if (
        !all(is.finite(score_limits)) ||
        score_limits[[1]] == score_limits[[2]]
    ) {

        score_limits <- range(
            plot_df$score,
            na.rm = TRUE
        )
    }

    ggplot(
        plot_df,
        aes(
            x = UMAP_1,
            y = UMAP_2,
            color = score
        )
    ) +
        geom_point(
            size = 0.45,
            alpha = 0.90
        ) +
        scale_color_gradient2(
            low = "#0033FF",
            mid = "#FFFFFF",
            high = "#FF1A1A",
            midpoint = median(
                plot_df$score,
                na.rm = TRUE
            ),
            limits = score_limits,
            oob = scales::squish
        ) +
        coord_equal() +
        labs(
            title = title_text,
            subtitle = paste0(
                "Fixed MΦ-only RPCA UMAP | ",
                "blue = low, white = midpoint, red = high"
            ),
            x = paste0(
                UMAP_REDUCTION,
                "_1"
            ),
            y = paste0(
                UMAP_REDUCTION,
                "_2"
            ),
            color = "QC score"
        ) +
        theme_classic(
            base_size = 11
        ) +
        theme(
            plot.title = element_text(
                face = "bold",
                size = 14
            ),
            plot.subtitle = element_text(
                size = 10
            )
        )
}


# ------------------------------------------------------------------------------
# 12. Contamination-score UMAPs
# ------------------------------------------------------------------------------

p_B <- plot_score_umap(
    mphi,
    SCORE_COLS[["B_cell"]],
    "B-cell contamination score"
)

p_T <- plot_score_umap(
    mphi,
    SCORE_COLS[["T_cell"]],
    "T-cell contamination score"
)

p_NK <- plot_score_umap(
    mphi,
    SCORE_COLS[["NK_cell"]],
    "NK-cell contamination score"
)

p_Neu <- plot_score_umap(
    mphi,
    SCORE_COLS[["Neutrophil"]],
    "Neutrophil contamination score"
)

p_contam_all <- (
    p_B |
    p_T
) / (
    p_NK |
    p_Neu
) +
    plot_annotation(
        title = "Potential non-MΦ contamination on fixed MΦ-only RPCA UMAP",
        subtitle = "QC only — no cells removed, no clustering or UMAP recalculated"
    )

ggsave(
    file.path(
        FIG_DIR,
        "01_contamination_scores_fixed_RPCA_UMAP_v4.11.5.pdf"
    ),
    p_contam_all,
    width = 13,
    height = 11,
    units = "in",
    device = cairo_pdf
)

ggsave(
    file.path(
        FIG_DIR,
        "01_contamination_scores_fixed_RPCA_UMAP_v4.11.5.jpg"
    ),
    p_contam_all,
    width = 13,
    height = 11,
    units = "in",
    dpi = 350,
    quality = 95
)


# ------------------------------------------------------------------------------
# 13. Individual high-specificity marker genes
# ------------------------------------------------------------------------------

MARKERS_TO_PLOT <- c(

    # B
    "Cd19",
    "Cd79a",
    "Pax5",
    "Ebf1",

    # T
    "Cd3d",
    "Cd3e",

    # NK
    "Nkg7",
    "Klrd1",
    "Ncr1",

    # Neutrophil
    "S100a8",
    "S100a9",
    "Retnlg",
    "Ly6g",
    "Mpo"
)

MARKERS_TO_PLOT <- intersect(
    MARKERS_TO_PLOT,
    available_genes
)

msg(
    "Individual contamination markers: ",
    paste(
        MARKERS_TO_PLOT,
        collapse = ", "
    )
)


# ------------------------------------------------------------------------------
# 14. Custom individual-gene UMAP
# ------------------------------------------------------------------------------

plot_gene_umap <- function(
    gene_name
) {

    emb <- Embeddings(
        mphi,
        reduction = UMAP_REDUCTION
    )

    expression_now <- as.numeric(
        rna_data[
            gene_name,
            rownames(emb)
        ]
    )

    plot_df <- tibble(
        cell = rownames(emb),
        UMAP_1 = emb[, 1],
        UMAP_2 = emb[, 2],
        expression = expression_now
    )

    upper_limit <- quantile(
        plot_df$expression,
        0.99,
        na.rm = TRUE
    )

    if (
        !is.finite(upper_limit) ||
        upper_limit <= 0
    ) {

        upper_limit <- max(
            plot_df$expression,
            na.rm = TRUE
        )
    }

    if (
        !is.finite(upper_limit) ||
        upper_limit <= 0
    ) {

        upper_limit <- 1
    }

    ggplot(
        plot_df,
        aes(
            UMAP_1,
            UMAP_2,
            color = expression
        )
    ) +
        geom_point(
            size = 0.35,
            alpha = 0.85
        ) +
        scale_color_gradientn(
            colours = c(
                "#0033FF",
                "#FFFFFF",
                "#FF1A1A"
            ),
            limits = c(
                0,
                upper_limit
            ),
            oob = scales::squish
        ) +
        coord_equal() +
        labs(
            title = gene_name,
            subtitle = "Normalized RNA expression | fixed RPCA UMAP",
            x = NULL,
            y = NULL,
            color = "Expression"
        ) +
        theme_classic(
            base_size = 9
        ) +
        theme(
            plot.title = element_text(
                face = "bold"
            )
        )
}


gene_plots <- lapply(
    MARKERS_TO_PLOT,
    plot_gene_umap
)

p_gene_grid <- wrap_plots(
    gene_plots,
    ncol = 4
) +
    plot_annotation(
        title = "Individual non-MΦ lineage markers on fixed RPCA UMAP"
    )

ggsave(
    file.path(
        FIG_DIR,
        "02_individual_contamination_markers_fixed_RPCA_UMAP_v4.11.5.pdf"
    ),
    p_gene_grid,
    width = 14,
    height = max(
        8,
        ceiling(
            length(gene_plots) / 4
        ) * 3.5
    ),
    units = "in",
    device = cairo_pdf
)


# ------------------------------------------------------------------------------
# 15. Long-form score data
# ------------------------------------------------------------------------------

score_df <- mphi@meta.data %>%
    rownames_to_column(
        "cell"
    ) %>%
    select(
        cell,
        sample_4group,
        macrophage_class_Res2_v484,
        all_of(
            unname(SCORE_COLS)
        )
    ) %>%
    pivot_longer(
        cols = all_of(
            unname(SCORE_COLS)
        ),
        names_to = "score_column",
        values_to = "score"
    ) %>%
    mutate(
        lineage = case_when(
            score_column == SCORE_COLS[["B_cell"]] ~ "B cell",
            score_column == SCORE_COLS[["T_cell"]] ~ "T cell",
            score_column == SCORE_COLS[["NK_cell"]] ~ "NK cell",
            score_column == SCORE_COLS[["Neutrophil"]] ~ "Neutrophil",
            TRUE ~ score_column
        )
    )

write.csv(
    score_df,
    file.path(
        TAB_DIR,
        "02_cell_level_contamination_scores_long_v4.11.5.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 16. Contamination score by fixed MΦ subtype
# ------------------------------------------------------------------------------

p_subtype <- ggplot(
    score_df,
    aes(
        x = macrophage_class_Res2_v484,
        y = score
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
    facet_wrap(
        ~ lineage,
        scales = "free_y",
        ncol = 2
    ) +
    labs(
        title = "Potential contamination scores by fixed MΦ subtype",
        x = NULL,
        y = "Mean normalized lineage-marker expression"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),
        plot.title = element_text(
            face = "bold"
        )
    )

ggsave(
    file.path(
        FIG_DIR,
        "03_contamination_scores_by_Mphi_subtype_v4.11.5.pdf"
    ),
    p_subtype,
    width = 12,
    height = 8,
    units = "in",
    device = cairo_pdf
)


# ------------------------------------------------------------------------------
# 17. Contamination score by biological sample
# ------------------------------------------------------------------------------

p_sample <- ggplot(
    score_df,
    aes(
        x = sample_4group,
        y = score
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
    facet_wrap(
        ~ lineage,
        scales = "free_y",
        ncol = 2
    ) +
    labs(
        title = "Potential contamination scores by biological sample",
        x = NULL,
        y = "Mean normalized lineage-marker expression"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),
        plot.title = element_text(
            face = "bold"
        )
    )

ggsave(
    file.path(
        FIG_DIR,
        "04_contamination_scores_by_sample_v4.11.5.pdf"
    ),
    p_sample,
    width = 12,
    height = 8,
    units = "in",
    device = cairo_pdf
)


# ------------------------------------------------------------------------------
# 18. Cell-level expression of highly specific markers
#
# Module / mean-expression scores alone should NOT determine cell removal.
# ------------------------------------------------------------------------------

SPECIFIC_MARKERS <- list(

    B_cell = intersect(
        c(
            "Cd19",
            "Cd79a",
            "Cd79b",
            "Ms4a1",
            "Pax5",
            "Ebf1"
        ),
        available_genes
    ),

    T_cell = intersect(
        c(
            "Cd3d",
            "Cd3e",
            "Cd3g",
            "Trac"
        ),
        available_genes
    ),

    NK_cell = intersect(
        c(
            "Ncr1",
            "Klrd1",
            "Klrk1"
        ),
        available_genes
    ),

    Neutrophil = intersect(
        c(
            "Ly6g",
            "Mpo",
            "Elane",
            "Csf3r",
            "Retnlg"
        ),
        available_genes
    )
)

specific_marker_audit <- bind_rows(
    lapply(
        names(SPECIFIC_MARKERS),
        function(lineage_name) {

            tibble(
                lineage = lineage_name,
                gene = SPECIFIC_MARKERS[[lineage_name]]
            )
        }
    )
)

write.csv(
    specific_marker_audit,
    file.path(
        TAB_DIR,
        "03_specific_marker_panel_v4.11.5.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 19. Number of specific markers detected per cell
# ------------------------------------------------------------------------------

count_positive_markers <- function(
    marker_genes
) {

    if (length(marker_genes) == 0L) {

        return(
            rep(
                0L,
                ncol(mphi)
            )
        )
    }

    mat <- rna_data[
        marker_genes,
        ,
        drop = FALSE
    ]

    as.integer(
        Matrix::colSums(
            mat > 0
        )
    )
}

mphi$QC_B_specific_n <- count_positive_markers(
    SPECIFIC_MARKERS[["B_cell"]]
)

mphi$QC_T_specific_n <- count_positive_markers(
    SPECIFIC_MARKERS[["T_cell"]]
)

mphi$QC_NK_specific_n <- count_positive_markers(
    SPECIFIC_MARKERS[["NK_cell"]]
)

mphi$QC_Neutrophil_specific_n <- count_positive_markers(
    SPECIFIC_MARKERS[["Neutrophil"]]
)


# ------------------------------------------------------------------------------
# 20. Conservative candidate flags
#
# IMPORTANT:
# These are QC flags, NOT automatic exclusion criteria.
#
# Requiring >=2 specific lineage markers reduces the chance that a single
# ambient transcript alone labels a macrophage as contaminating.
# ------------------------------------------------------------------------------

mphi$QC_B_candidate <- mphi$QC_B_specific_n >= 2L

mphi$QC_T_candidate <- mphi$QC_T_specific_n >= 2L

mphi$QC_NK_candidate <- mphi$QC_NK_specific_n >= 2L

mphi$QC_Neutrophil_candidate <- mphi$QC_Neutrophil_specific_n >= 2L

mphi$QC_any_nonMphi_candidate <- (
    mphi$QC_B_candidate |
    mphi$QC_T_candidate |
    mphi$QC_NK_candidate |
    mphi$QC_Neutrophil_candidate
)


# ------------------------------------------------------------------------------
# 21. Candidate-cell table
# ------------------------------------------------------------------------------

candidate_table <- mphi@meta.data %>%
    rownames_to_column(
        "cell"
    ) %>%
    filter(
        QC_any_nonMphi_candidate
    ) %>%
    select(
        cell,
        sample_4group,
        macrophage_class_Res2_v484,
        QC_B_cell_score,
        QC_T_cell_score,
        QC_NK_cell_score,
        QC_Neutrophil_score,
        QC_B_specific_n,
        QC_T_specific_n,
        QC_NK_specific_n,
        QC_Neutrophil_specific_n,
        QC_B_candidate,
        QC_T_candidate,
        QC_NK_candidate,
        QC_Neutrophil_candidate
    )

write.csv(
    candidate_table,
    file.path(
        TAB_DIR,
        "04_candidate_nonMphi_cells_v4.11.5.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 22. Overall contamination burden
# ------------------------------------------------------------------------------

overall_summary <- tibble(

    category = c(
        "B_cell",
        "T_cell",
        "NK_cell",
        "Neutrophil",
        "Any_nonMphi"
    ),

    n_cells = c(
        sum(
            mphi$QC_B_candidate
        ),
        sum(
            mphi$QC_T_candidate
        ),
        sum(
            mphi$QC_NK_candidate
        ),
        sum(
            mphi$QC_Neutrophil_candidate
        ),
        sum(
            mphi$QC_any_nonMphi_candidate
        )
    )
) %>%
    mutate(
        total_Mphi_cells = ncol(mphi),
        percent = n_cells /
            total_Mphi_cells *
            100
    )

write.csv(
    overall_summary,
    file.path(
        TAB_DIR,
        "05_overall_contamination_burden_v4.11.5.csv"
    ),
    row.names = FALSE
)

print(
    overall_summary
)


# ------------------------------------------------------------------------------
# 23. Contamination burden by sample
# ------------------------------------------------------------------------------

sample_summary <- mphi@meta.data %>%
    rownames_to_column(
        "cell"
    ) %>%
    group_by(
        sample_4group
    ) %>%
    summarise(
        total_cells = n(),

        B_candidate = sum(
            QC_B_candidate
        ),

        T_candidate = sum(
            QC_T_candidate
        ),

        NK_candidate = sum(
            QC_NK_candidate
        ),

        Neutrophil_candidate = sum(
            QC_Neutrophil_candidate
        ),

        Any_candidate = sum(
            QC_any_nonMphi_candidate
        ),

        Any_percent = 100 *
            Any_candidate /
            total_cells,

        .groups = "drop"
    )

write.csv(
    sample_summary,
    file.path(
        TAB_DIR,
        "06_contamination_burden_by_sample_v4.11.5.csv"
    ),
    row.names = FALSE
)

print(
    sample_summary
)


# ------------------------------------------------------------------------------
# 24. Contamination burden by MΦ subtype
# ------------------------------------------------------------------------------

subtype_summary <- mphi@meta.data %>%
    rownames_to_column(
        "cell"
    ) %>%
    group_by(
        macrophage_class_Res2_v484
    ) %>%
    summarise(
        total_cells = n(),

        B_candidate = sum(
            QC_B_candidate
        ),

        T_candidate = sum(
            QC_T_candidate
        ),

        NK_candidate = sum(
            QC_NK_candidate
        ),

        Neutrophil_candidate = sum(
            QC_Neutrophil_candidate
        ),

        Any_candidate = sum(
            QC_any_nonMphi_candidate
        ),

        Any_percent = 100 *
            Any_candidate /
            total_cells,

        .groups = "drop"
    )

write.csv(
    subtype_summary,
    file.path(
        TAB_DIR,
        "07_contamination_burden_by_Mphi_subtype_v4.11.5.csv"
    ),
    row.names = FALSE
)

print(
    subtype_summary
)


# ------------------------------------------------------------------------------
# 25. Sample × subtype contamination burden
# ------------------------------------------------------------------------------

sample_subtype_summary <- mphi@meta.data %>%
    rownames_to_column(
        "cell"
    ) %>%
    group_by(
        sample_4group,
        macrophage_class_Res2_v484
    ) %>%
    summarise(
        total_cells = n(),

        B_candidate = sum(
            QC_B_candidate
        ),

        T_candidate = sum(
            QC_T_candidate
        ),

        NK_candidate = sum(
            QC_NK_candidate
        ),

        Neutrophil_candidate = sum(
            QC_Neutrophil_candidate
        ),

        Any_candidate = sum(
            QC_any_nonMphi_candidate
        ),

        Any_percent = 100 *
            Any_candidate /
            total_cells,

        .groups = "drop"
    )

write.csv(
    sample_subtype_summary,
    file.path(
        TAB_DIR,
        "08_contamination_burden_sample_x_subtype_v4.11.5.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 26. UMAP: conservative candidate cells
# ------------------------------------------------------------------------------

emb <- Embeddings(
    mphi,
    reduction = UMAP_REDUCTION
)

candidate_umap_df <- tibble(
    cell = rownames(emb),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    candidate = ifelse(
        mphi@meta.data[
            rownames(emb),
            "QC_any_nonMphi_candidate"
        ],
        "Candidate non-MΦ",
        "Other MΦ"
    )
)

p_candidate <- ggplot(
    candidate_umap_df,
    aes(
        UMAP_1,
        UMAP_2
    )
) +
    geom_point(
        data = candidate_umap_df %>%
            filter(
                candidate ==
                    "Other MΦ"
            ),
        color = "grey85",
        size = 0.35,
        alpha = 0.5
    ) +
    geom_point(
        data = candidate_umap_df %>%
            filter(
                candidate ==
                    "Candidate non-MΦ"
            ),
        color = "#FF1A1A",
        size = 0.8,
        alpha = 0.9
    ) +
    coord_equal() +
    labs(
        title = "Conservative candidate non-MΦ cells",
        subtitle = "Red = ≥2 lineage-specific markers from B/T/NK/neutrophil panel",
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

ggsave(
    file.path(
        FIG_DIR,
        "05_candidate_nonMphi_cells_fixed_RPCA_UMAP_v4.11.5.pdf"
    ),
    p_candidate,
    width = 8,
    height = 7,
    units = "in",
    device = cairo_pdf
)


# ------------------------------------------------------------------------------
# 27. UMAP: candidate lineage separately
# ------------------------------------------------------------------------------

candidate_lineage_df <- tibble(
    cell = rownames(emb),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    B = mphi@meta.data[
        rownames(emb),
        "QC_B_candidate"
    ],
    T = mphi@meta.data[
        rownames(emb),
        "QC_T_candidate"
    ],
    NK = mphi@meta.data[
        rownames(emb),
        "QC_NK_candidate"
    ],
    Neutrophil = mphi@meta.data[
        rownames(emb),
        "QC_Neutrophil_candidate"
    ]
) %>%
    pivot_longer(
        cols = c(
            B,
            T,
            NK,
            Neutrophil
        ),
        names_to = "lineage",
        values_to = "candidate"
    )

p_candidate_by_lineage <- ggplot(
    candidate_lineage_df,
    aes(
        UMAP_1,
        UMAP_2
    )
) +
    geom_point(
        color = "grey88",
        size = 0.25,
        alpha = 0.35
    ) +
    geom_point(
        data = candidate_lineage_df %>%
            filter(
                candidate
            ),
        color = "#FF1A1A",
        size = 0.75,
        alpha = 0.9
    ) +
    facet_wrap(
        ~ lineage,
        ncol = 2
    ) +
    coord_equal() +
    labs(
        title = "Conservative lineage-specific contamination candidates",
        subtitle = "Red = ≥2 specific markers expressed in the same cell",
        x = NULL,
        y = NULL
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        strip.text = element_text(
            face = "bold"
        )
    )

ggsave(
    file.path(
        FIG_DIR,
        "06_candidate_nonMphi_cells_by_lineage_fixed_RPCA_UMAP_v4.11.5.pdf"
    ),
    p_candidate_by_lineage,
    width = 11,
    height = 9,
    units = "in",
    device = cairo_pdf
)


# ------------------------------------------------------------------------------
# 28. Save QC-annotated object
#
# This is a NEW QC object.
# Original v4.8.4 RDS is never overwritten.
# ------------------------------------------------------------------------------

QC_RDS <- file.path(
    OUTPUT_DIR,
    "Mouse_Mphi_Res2_contamination_QC_annotated_v4.11.5.rds"
)

saveRDS(
    mphi,
    QC_RDS
)


# ------------------------------------------------------------------------------
# 29. README
# ------------------------------------------------------------------------------

readme <- c(

    "Mouse MASH MΦ Res2 contamination QC v4.11.5",
    "",

    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",

    "Purpose:",
    "  Identify possible B/T/NK/neutrophil contamination or doublets.",
    "",

    "Fixed framework:",
    "  Res2.0 clustering unchanged.",
    "  v4.8.4 macrophage annotation unchanged.",
    paste0(
        "  Existing MΦ-only RPCA UMAP used: ",
        UMAP_REDUCTION
    ),
    "  PCA/RPCA/UMAP NOT recalculated.",
    "  No cells removed.",
    "",

    "QC scores:",
    "  AddModuleScore is NOT used.",
    "  Score = mean normalized expression of lineage-marker genes per cell.",
    "",

    "QC programs:",
    "  B cell",
    "  T cell",
    "  NK cell",
    "  Neutrophil",
    "",

    "Conservative candidate definition:",
    "  >= 2 expressed lineage-specific markers in the same cell.",
    "  This is a QC flag, NOT an automatic exclusion criterion.",
    "",

    "Important:",
    "  Shared macrophage/immune genes can produce lineage-score signal.",
    "  Therefore candidate removal should only be considered after",
    "  inspection of UMAP localization, specific marker co-expression,",
    "  sample dependence and macrophage identity.",
    "",

    "Primary figures:",
    "  01 contamination scores on fixed RPCA UMAP",
    "  02 individual marker UMAP",
    "  03 contamination score by MΦ subtype",
    "  04 contamination score by sample",
    "  05 any conservative candidate-cell UMAP",
    "  06 lineage-specific candidate-cell UMAP",
    "",

    "No biological cells were excluded in v4.11.5."
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.11.5.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.11.5.txt"
    )
)


# ------------------------------------------------------------------------------
# 30. Final
# ------------------------------------------------------------------------------

msg(
    "DONE."
)

msg(
    "Output: ",
    OUTPUT_DIR
)

msg(
    "Candidate non-MΦ cells: ",
    sum(
        mphi$QC_any_nonMphi_candidate
    ),
    " / ",
    ncol(mphi),
    " (",
    round(
        100 *
            mean(
                mphi$QC_any_nonMphi_candidate
            ),
        3
    ),
    "%)"
)
