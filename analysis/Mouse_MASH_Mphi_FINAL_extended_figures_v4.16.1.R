#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# FINAL extended figure regeneration
# v4.16.1
#
# PURPOSE
#   Regenerate and consolidate the major macrophage figures used throughout
#   v4.12-v4.15 into one FINAL output tree, using the finalized v4.14.5
#   macrophage annotation and the v4.15.5 publication UMAP style.
#
# FINAL PARENT FRAMEWORK
#   Res2.0
#
# FINAL ANNOTATION COLUMN
#   macrophage_class_Res2_FINAL_v4145
#
# FINAL CLASSES
#   1) Inflammatory-MΦ
#   2) Anti-inflammatory-MΦ
#   3) ECM-associated inflammatory-MΦ
#   4) Repair/Resolution-MΦ
#   5) Lipid-associated/TREM2-MΦ
#   + Other
#
# OUTPUT SECTIONS
#   01_UMAP
#   02_Composition
#   03_M1_M2_programs
#   04_AntiInflammatory_heterogeneity
#   05_IL10_response
#   06_Pseudobulk_DE
#   07_Pathway_enrichment
#   08_Robustness
#   Tables
#   Logs
#
# IMPORTANT
#   - No reclustering of the parent MΦ object.
#   - No change in final parent MΦ annotation.
#   - STD vs CDAHFD is descriptive n=1 vs n=1.
#   - Sham vs Tx uses biological sample-level summaries (n=2 vs n=2).
#   - Optional legacy/secondary inputs are skipped if not found.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4160)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

FINAL_RDS <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

ANTI_RDS <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "AntiInflammatory_functional_state_audit_CleanB_v4.14.0.1",
    "RDS",
    "Mouse_Mphi_AntiInflammatory_Res1.2_functional_state_annotated_v4.14.0.1.rds"
)

PATHWAY_TABLE <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_subtype_pathway_enrichment_CleanB_v4.14.2",
    "Tables",
    "05_fgsea_Tx_vs_Sham_all_subtypes_v4.14.2.csv"
)

ROBUSTNESS_DIR <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_final_classification_robustness_v4.14.3",
    "Tables"
)

OUTPUT_DIR <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "FINAL_Extended_Figures_v4.16.1"
)

DIR_01 <- file.path(OUTPUT_DIR, "01_UMAP")
DIR_02 <- file.path(OUTPUT_DIR, "02_Composition")
DIR_03 <- file.path(OUTPUT_DIR, "03_M1_M2_programs")
DIR_04 <- file.path(OUTPUT_DIR, "04_AntiInflammatory_heterogeneity")
DIR_05 <- file.path(OUTPUT_DIR, "05_IL10_response")
DIR_06 <- file.path(OUTPUT_DIR, "06_Pseudobulk_DE")
DIR_07 <- file.path(OUTPUT_DIR, "07_Pathway_enrichment")
DIR_08 <- file.path(OUTPUT_DIR, "08_Robustness")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

for (d in c(
    DIR_01, DIR_02, DIR_03, DIR_04,
    DIR_05, DIR_06, DIR_07, DIR_08,
    TAB_DIR, LOG_DIR
)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

ASSAY_USE <- "RNA"
FINAL_CLASS_COL <- "macrophage_class_Res2_FINAL_v4145"

# ------------------------------------------------------------------------------
# 1. Ordering / colors
# ------------------------------------------------------------------------------

CLASS_ORDER <- c(
    "Inflammatory-Mphi",
    "Anti-inflammatory-Mphi",
    "ECM-associated inflammatory-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi",
    "Other"
)

CLASS_LABELS <- c(
    "Inflammatory-Mphi" = "Inflammatory-MΦ",
    "Anti-inflammatory-Mphi" = "Anti-inflammatory-MΦ",
    "ECM-associated inflammatory-Mphi" = "ECM-associated inflammatory-MΦ",
    "Repair/Resolution-Mphi" = "Repair/Resolution-MΦ",
    "Lipid-associated/TREM2-Mphi" = "Lipid-associated/TREM2-MΦ",
    "Other" = "Other"
)

CLASS_COLORS <- c(
    "Inflammatory-Mphi" = "#FF2D2D",
    "Anti-inflammatory-Mphi" = "#00AEEF",
    "ECM-associated inflammatory-Mphi" = "#FF8C00",
    "Repair/Resolution-Mphi" = "#00C853",
    "Lipid-associated/TREM2-Mphi" = "#7B2CBF",
    "Other" = "#808080"
)

CONDITION_ORDER <- c("STD", "CDAHFD", "Sham", "Tx")

CONDITION_COLORS <- c(
    "STD" = "#0066FF",
    "CDAHFD" = "#FF1A1A",
    "Sham" = "#FF1A1A",
    "Tx" = "#0066FF"
)

SAMPLE_ORDER <- c(
    "STD_rep1",
    "CDHFD_rep1",
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
)

# v4.15.5 UMAP publication window / density
UMAP_XLIM <- c(-8.8, -2.0)
UMAP_YLIM <- c(-7.1, -0.5)
UMAP_PT_SIZE_OTHER <- 0.60
UMAP_PT_SIZE_MAJOR <- 1.20
UMAP_ALPHA_OTHER <- 0.42
UMAP_ALPHA_MAJOR <- 0.95
UMAP_LEGEND_POINT_SIZE <- 4.4

# ------------------------------------------------------------------------------
# 2. Functional programs
# ------------------------------------------------------------------------------

PROGRAMS <- list(

    M1 = c(
        "Il1b","Tnf","Ccl2","Cxcl10",
        "Stat1","Cd80","Cd86","Nos2"
    ),

    M2 = c(
        "Mrc1","Cd163","Il1rn","Mertk",
        "Igf1","Hmox1","Klf4","Maf"
    ),

    Anti_inflammatory = c(
        "Mrc1","Cd163","Il1rn","Mertk",
        "Igf1","Hmox1","Klf4","Maf"
    ),

    IL10_STAT3 = c(
        "Il10ra","Il10rb","Jak1","Tyk2",
        "Stat3","Socs3","Bcl3","Il1rn"
    ),

    Repair_Resolution = c(
        "Mertk","Axl","Mfge8","Gas6",
        "Igf1","Hmox1","Mmp13","Mmp14","Plau"
    ),

    Inflammatory = c(
        "Il1b","Tnf","Ccl2","Cxcl10",
        "Stat1","Cd80","Cd86"
    ),

    ECM_associated_inflammatory = c(
        "Thbs1","Fn1","Tgfb1",
        "Col1a1","Col1a2","Col3a1",
        "Il1b","Stat1","Cxcl10"
    ),

    Lipid_TREM2 = c(
        "Trem2","Gpnmb","Cd9","Lpl",
        "Apoe","Fabp5","Abca1","Plin2"
    ),

    Efferocytosis = c(
        "Mertk","Axl","Mfge8","Gas6",
        "C1qa","C1qb","C1qc","Lgals3"
    )
)

ANTI_DOTPLOT_GENES <- c(
    "Apoe","C1qa","Cd163","Cd9","Cxcl10",
    "Gas6","Gpnmb","Hmox1","Igf1","Il1rn",
    "Il1b","Il10ra","Mertk","Mfge8","Mrc1",
    "Mmp12","Mmp13","Mmp14","Plau","Socs3",
    "Spp1","Tgfb1","Thbs1","Tnf","Trem2"
)

# ------------------------------------------------------------------------------
# 3. Packages
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
    "scales",
    "pheatmap"
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
    library(pheatmap)
})

# ------------------------------------------------------------------------------
# 4. Helpers
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
    if (length(hit) == 0L) return(NA_character_)
    hit[[1]]
}

canonical_condition <- function(sample_name) {
    x <- as.character(sample_name)
    out <- rep(NA_character_, length(x))
    out[grepl("^STD", x, ignore.case = TRUE)] <- "STD"
    out[grepl("CDAHFD|CDHFD", x, ignore.case = TRUE)] <- "CDAHFD"
    out[grepl("^Sham", x, ignore.case = TRUE)] <- "Sham"
    out[grepl("^Tx", x, ignore.case = TRUE)] <- "Tx"
    factor(out, levels = CONDITION_ORDER)
}

get_layer_safe <- function(object, assay, layer) {
    x <- tryCatch(
        LayerData(
            object,
            assay = assay,
            layer = layer
        ),
        error = function(e) NULL
    )

    if (is.null(x)) {
        x <- tryCatch(
            GetAssayData(
                object,
                assay = assay,
                slot = layer
            ),
            error = function(e) NULL
        )
    }

    x
}

save_pdf <- function(path, plot, width, height) {
    ggsave(
        filename = path,
        plot = plot,
        device = cairo_pdf,
        width = width,
        height = height,
        units = "in",
        limitsize = FALSE
    )
}

mean_expression_score <- function(mat, genes) {
    genes <- intersect(genes, rownames(mat))
    if (length(genes) == 0L) {
        return(rep(NA_real_, ncol(mat)))
    }

    as.numeric(
        Matrix::colMeans(
            mat[genes, , drop = FALSE]
        )
    )
}

safe_z <- function(x) {
    z <- as.numeric(scale(x))
    z[!is.finite(z)] <- 0
    z
}

safe_log2fc <- function(a, b, pseudo = 1) {
    log2((a + pseudo) / (b + pseudo))
}

# ------------------------------------------------------------------------------
# 5. Load FINAL MΦ object
# ------------------------------------------------------------------------------

if (!file.exists(FINAL_RDS)) {
    stop("FINAL RDS not found:\n", FINAL_RDS)
}

msg("Loading FINAL MΦ RDS: ", FINAL_RDS)

mphi <- readRDS(FINAL_RDS)

if (!inherits(mphi, "Seurat")) {
    stop("FINAL input is not a Seurat object.")
}

if (!FINAL_CLASS_COL %in% colnames(mphi@meta.data)) {
    stop("FINAL annotation column missing: ", FINAL_CLASS_COL)
}

DefaultAssay(mphi) <- ASSAY_USE

mphi <- JoinLayers(
    mphi,
    assay = ASSAY_USE
)

if (!"data" %in% Layers(mphi[[ASSAY_USE]])) {
    mphi <- NormalizeData(
        mphi,
        assay = ASSAY_USE,
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
    )
}

rna_data <- get_layer_safe(
    mphi,
    ASSAY_USE,
    "data"
)

counts <- get_layer_safe(
    mphi,
    ASSAY_USE,
    "counts"
)

if (is.null(rna_data)) stop("RNA data layer not found.")
if (is.null(counts)) stop("RNA counts layer not found.")

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

if (is.na(SAMPLE_COL)) {
    stop("Sample column not found.")
}

available_reductions <- Reductions(mphi)

REDUCTION_USE <- first_existing(
    available_reductions,
    c(
        "umapRPCA",
        "mphi.umap.rpca",
        "rpca.umap",
        "umap.rpca",
        "umap"
    )
)

if (is.na(REDUCTION_USE)) {
    umap_like <- available_reductions[
        grepl(
            "umap",
            available_reductions,
            ignore.case = TRUE
        )
    ]

    if (length(umap_like) == 0L) {
        stop("No UMAP reduction found.")
    }

    REDUCTION_USE <- umap_like[[1]]
}

mphi$sample_FINAL_v4160 <- as.character(
    mphi@meta.data[[SAMPLE_COL]]
)

mphi$condition_FINAL_v4160 <- canonical_condition(
    mphi$sample_FINAL_v4160
)

mphi$class_FINAL_v4160 <- factor(
    as.character(
        mphi@meta.data[[FINAL_CLASS_COL]]
    ),
    levels = CLASS_ORDER
)

meta <- mphi@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = sample_FINAL_v4160,
        condition = as.character(condition_FINAL_v4160),
        macrophage_class = as.character(class_FINAL_v4160)
    )

# ------------------------------------------------------------------------------
# 6. 01_UMAP
# ------------------------------------------------------------------------------

msg("01_UMAP")

emb <- as.data.frame(
    Embeddings(
        mphi,
        reduction = REDUCTION_USE
    )[, 1:2, drop = FALSE]
)

colnames(emb) <- c("UMAP_1", "UMAP_2")

umap_df <- emb %>%
    rownames_to_column("cell") %>%
    left_join(
        meta,
        by = "cell"
    ) %>%
    mutate(
        macrophage_class = factor(
            macrophage_class,
            levels = CLASS_ORDER
        ),
        condition = factor(
            condition,
            levels = CONDITION_ORDER
        )
    )

umap_other <- umap_df %>%
    filter(macrophage_class == "Other")

umap_major <- umap_df %>%
    filter(macrophage_class != "Other")

p_umap <- ggplot() +
    geom_point(
        data = umap_other,
        aes(
            x = UMAP_1,
            y = UMAP_2,
            color = macrophage_class
        ),
        size = UMAP_PT_SIZE_OTHER,
        alpha = UMAP_ALPHA_OTHER,
        stroke = 0
    ) +
    geom_point(
        data = umap_major,
        aes(
            x = UMAP_1,
            y = UMAP_2,
            color = macrophage_class
        ),
        size = UMAP_PT_SIZE_MAJOR,
        alpha = UMAP_ALPHA_MAJOR,
        stroke = 0
    ) +
    scale_color_manual(
        values = CLASS_COLORS,
        breaks = CLASS_ORDER,
        labels = CLASS_LABELS[CLASS_ORDER],
        drop = FALSE
    ) +
    guides(
        color = guide_legend(
            override.aes = list(
                size = UMAP_LEGEND_POINT_SIZE,
                alpha = 1
            )
        )
    ) +
    coord_cartesian(
        xlim = UMAP_XLIM,
        ylim = UMAP_YLIM,
        expand = FALSE
    ) +
    labs(
        title = "Mouse MASH macrophages",
        subtitle = "Clean-B FINAL annotation | Res2.0",
        x = "UMAP 1",
        y = "UMAP 2",
        color = "MΦ subtype"
    ) +
    theme_classic(base_size = 11) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 14
        ),
        legend.title = element_text(
            face = "bold"
        )
    )

save_pdf(
    file.path(
        DIR_01,
        "01_FINAL_Mphi_subtype_UMAP_v4.16.1.pdf"
    ),
    p_umap,
    9.5,
    7
)

p_umap_condition <- ggplot() +
    geom_point(
        data = umap_other,
        aes(
            x = UMAP_1,
            y = UMAP_2,
            color = macrophage_class
        ),
        size = UMAP_PT_SIZE_OTHER,
        alpha = UMAP_ALPHA_OTHER,
        stroke = 0
    ) +
    geom_point(
        data = umap_major,
        aes(
            x = UMAP_1,
            y = UMAP_2,
            color = macrophage_class
        ),
        size = UMAP_PT_SIZE_MAJOR,
        alpha = UMAP_ALPHA_MAJOR,
        stroke = 0
    ) +
    facet_wrap(
        ~ condition,
        ncol = 2,
        drop = FALSE
    ) +
    scale_color_manual(
        values = CLASS_COLORS,
        breaks = CLASS_ORDER,
        labels = CLASS_LABELS[CLASS_ORDER],
        drop = FALSE
    ) +
    guides(
        color = guide_legend(
            override.aes = list(
                size = UMAP_LEGEND_POINT_SIZE,
                alpha = 1
            )
        )
    ) +
    coord_cartesian(
        xlim = UMAP_XLIM,
        ylim = UMAP_YLIM,
        expand = FALSE
    ) +
    labs(
        title = "FINAL MΦ subtype distribution by condition",
        x = "UMAP 1",
        y = "UMAP 2",
        color = "MΦ subtype"
    ) +
    theme_classic(base_size = 9) +
    theme(
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")
    )

save_pdf(
    file.path(
        DIR_01,
        "02_FINAL_Mphi_subtype_UMAP_by_condition_v4.16.1.pdf"
    ),
    p_umap_condition,
    12,
    9
)

# ------------------------------------------------------------------------------
# 7. 02_Composition
# ------------------------------------------------------------------------------

msg("02_Composition")

sample_totals <- meta %>%
    count(
        sample,
        name = "total_Mphi"
    )

composition <- meta %>%
    count(
        sample,
        condition,
        macrophage_class,
        name = "n_cells"
    ) %>%
    complete(
        sample = SAMPLE_ORDER,
        macrophage_class = CLASS_ORDER,
        fill = list(n_cells = 0L)
    ) %>%
    mutate(
        condition = as.character(
            canonical_condition(sample)
        )
    ) %>%
    left_join(
        sample_totals,
        by = "sample"
    ) %>%
    mutate(
        fraction = n_cells / total_Mphi,
        percent = 100 * fraction,
        sample = factor(
            sample,
            levels = SAMPLE_ORDER
        ),
        condition = factor(
            condition,
            levels = CONDITION_ORDER
        ),
        macrophage_class = factor(
            macrophage_class,
            levels = CLASS_ORDER
        )
    )

write.csv(
    composition,
    file.path(
        TAB_DIR,
        "01_FINAL_subtype_composition_by_sample_v4.16.1.csv"
    ),
    row.names = FALSE
)

p_comp <- ggplot(
    composition,
    aes(
        x = sample,
        y = percent,
        fill = macrophage_class
    )
) +
    geom_col(width = 0.72) +
    scale_fill_manual(
        values = CLASS_COLORS,
        breaks = CLASS_ORDER,
        labels = CLASS_LABELS[CLASS_ORDER],
        drop = FALSE
    ) +
    labs(
        title = "FINAL MΦ subtype composition",
        subtitle = "Biological-sample level",
        x = NULL,
        y = "% of Clean-B MΦ",
        fill = "MΦ subtype"
    ) +
    theme_classic(base_size = 10) +
    theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(
            angle = 35,
            hjust = 1
        )
    )

save_pdf(
    file.path(
        DIR_02,
        "01_FINAL_subtype_composition_by_sample_v4.16.1.pdf"
    ),
    p_comp,
    10,
    6
)

metric_counts <- meta %>%
    count(
        sample,
        condition,
        macrophage_class,
        name = "n_cells"
    ) %>%
    complete(
        sample = SAMPLE_ORDER,
        macrophage_class = CLASS_ORDER,
        fill = list(n_cells = 0L)
    ) %>%
    mutate(
        condition = as.character(
            canonical_condition(sample)
        )
    ) %>%
    pivot_wider(
        names_from = macrophage_class,
        values_from = n_cells,
        values_fill = 0
    ) %>%
    left_join(
        sample_totals,
        by = "sample"
    ) %>%
    mutate(
        M1_n = .data[["Inflammatory-Mphi"]],
        M2_n = .data[["Anti-inflammatory-Mphi"]],

        AntiInflammatory_percent_total_Mphi =
            100 * M2_n / total_Mphi,

        M2_percent_within_M1_M2 =
            100 * M2_n / (M1_n + M2_n),

        M2_to_M1_ratio =
            (M2_n + 0.5) / (M1_n + 0.5),

        sample = factor(
            sample,
            levels = SAMPLE_ORDER
        ),

        condition = factor(
            condition,
            levels = CONDITION_ORDER
        )
    )

write.csv(
    metric_counts,
    file.path(
        TAB_DIR,
        "02_FINAL_M1_M2_abundance_metrics_v4.16.1.csv"
    ),
    row.names = FALSE
)

make_two_panel_metric <- function(
    df,
    y_col,
    y_label,
    panel_letter,
    panel_title,
    file_name
) {

    disease <- df %>%
        filter(
            condition %in%
                c("STD", "CDAHFD")
        )

    treatment <- df %>%
        filter(
            condition %in%
                c("Sham", "Tx")
        )

    p1 <- ggplot(
        disease,
        aes(
            x = condition,
            y = .data[[y_col]],
            group = 1
        )
    ) +
        geom_line(
            linewidth = 0.6
        ) +
        geom_point(
            aes(color = condition),
            size = 3
        ) +
        scale_color_manual(
            values = CONDITION_COLORS,
            guide = "none"
        ) +
        labs(
            title = paste0(
                panel_letter,
                ". ",
                panel_title
            ),
            subtitle = "Disease comparison: STD → CDAHFD",
            x = NULL,
            y = y_label
        ) +
        theme_classic(base_size = 10) +
        theme(
            plot.title = element_text(
                face = "bold"
            )
        )

    means <- treatment %>%
        group_by(condition) %>%
        summarise(
            mean_value = mean(
                .data[[y_col]],
                na.rm = TRUE
            ),
            .groups = "drop"
        )

    p2 <- ggplot(
        treatment,
        aes(
            x = condition,
            y = .data[[y_col]]
        )
    ) +
        geom_line(
            data = means,
            aes(
                x = condition,
                y = mean_value,
                group = 1
            ),
            inherit.aes = FALSE,
            linewidth = 0.6
        ) +
        geom_point(
            data = means,
            aes(
                x = condition,
                y = mean_value,
                color = condition
            ),
            inherit.aes = FALSE,
            size = 3
        ) +
        geom_point(
            aes(shape = sample),
            size = 2.6,
            color = "black"
        ) +
        scale_color_manual(
            values = CONDITION_COLORS,
            guide = "none"
        ) +
        labs(
            title = paste0(
                panel_letter,
                ". ",
                panel_title
            ),
            subtitle = "Treatment comparison: Sham → Tx | points = biological samples",
            x = NULL,
            y = y_label,
            shape = "Sample"
        ) +
        theme_classic(base_size = 10) +
        theme(
            plot.title = element_text(
                face = "bold"
            )
        )

    p <- p1 + p2 +
        patchwork::plot_layout(
            widths = c(1, 1.15)
        )

    save_pdf(
        file.path(
            DIR_02,
            file_name
        ),
        p,
        12,
        5.8
    )
}

make_two_panel_metric(
    metric_counts,
    "AntiInflammatory_percent_total_Mphi",
    "Anti-inflammatory-MΦ / all MΦ (%)",
    "A",
    "Anti-inflammatory-MΦ fraction",
    "02_A_AntiInflammatory_fraction_v4.16.1.pdf"
)

make_two_panel_metric(
    metric_counts,
    "M2_percent_within_M1_M2",
    "M2 / (M1 + M2) (%)",
    "B",
    "M2 polarization fraction",
    "03_B_M2_polarization_fraction_v4.16.1.pdf"
)

make_two_panel_metric(
    metric_counts,
    "M2_to_M1_ratio",
    "(M2 + 0.5) / (M1 + 0.5)",
    "C",
    "M2 / M1 ratio",
    "04_C_M2_to_M1_ratio_v4.16.1.pdf"
)

# ------------------------------------------------------------------------------
# 8. 03_M1_M2_programs
# ------------------------------------------------------------------------------

msg("03_M1_M2_programs")

mphi$M1_score_FINAL_v4160 <- mean_expression_score(
    rna_data,
    PROGRAMS$M1
)

mphi$M2_score_FINAL_v4160 <- mean_expression_score(
    rna_data,
    PROGRAMS$M2
)

m1m2_df <- mphi@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = sample_FINAL_v4160,
        condition = as.character(
            condition_FINAL_v4160
        ),
        M1_score = M1_score_FINAL_v4160,
        M2_score = M2_score_FINAL_v4160
    ) %>%
    filter(
        condition %in%
            c("STD", "CDAHFD")
    ) %>%
    pivot_longer(
        cols = c(
            M1_score,
            M2_score
        ),
        names_to = "program",
        values_to = "score"
    ) %>%
    mutate(
        program = recode(
            program,
            M1_score = "M1 score",
            M2_score = "M2 score"
        ),
        condition = factor(
            condition,
            levels = c(
                "STD",
                "CDAHFD"
            )
        )
    )

write.csv(
    m1m2_df,
    file.path(
        TAB_DIR,
        "03_M1_M2_program_cell_scores_STD_CDAHFD_v4.16.1.csv"
    ),
    row.names = FALSE
)

p_m1m2 <- ggplot(
    m1m2_df,
    aes(
        x = condition,
        y = score,
        fill = condition
    )
) +
    geom_violin(
        trim = TRUE,
        scale = "width",
        linewidth = 0.25
    ) +
    geom_boxplot(
        width = 0.12,
        outlier.shape = NA,
        alpha = 0.55
    ) +
    facet_wrap(
        ~ program,
        nrow = 1,
        scales = "free_y"
    ) +
    scale_fill_manual(
        values = c(
            STD = "#3366FF",
            CDAHFD = "#FF3B30"
        )
    ) +
    labs(
        title = "M1 / M2 marker programs: STD vs CDAHFD",
        x = NULL,
        y = "Mean normalized marker expression"
    ) +
    theme_classic(base_size = 10) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        strip.text = element_text(
            face = "bold"
        ),
        legend.position = "none"
    )

save_pdf(
    file.path(
        DIR_03,
        "01_M1_M2_marker_programs_STD_vs_CDAHFD_v4.16.1.pdf"
    ),
    p_m1m2,
    8.5,
    5.2
)


# ------------------------------------------------------------------------------
# 8b. M1 / M2 marker programs: Sham vs Tx
# ------------------------------------------------------------------------------

msg("03b_M1_M2_programs_Sham_vs_Tx")

m1m2_sham_tx_df <- mphi@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = sample_FINAL_v4160,
        condition = as.character(
            condition_FINAL_v4160
        ),
        M1_score = M1_score_FINAL_v4160,
        M2_score = M2_score_FINAL_v4160
    ) %>%
    filter(
        condition %in%
            c(
                "Sham",
                "Tx"
            )
    ) %>%
    pivot_longer(
        cols = c(
            M1_score,
            M2_score
        ),
        names_to = "program",
        values_to = "score"
    ) %>%
    mutate(
        program = recode(
            program,
            M1_score = "M1 score",
            M2_score = "M2 score"
        ),
        condition = factor(
            condition,
            levels = c(
                "Sham",
                "Tx"
            )
        ),
        sample = factor(
            sample,
            levels = c(
                "Sham1",
                "Sham20",
                "Tx17",
                "Tx5"
            )
        )
    )

write.csv(
    m1m2_sham_tx_df,
    file.path(
        TAB_DIR,
        "03b_M1_M2_program_cell_scores_Sham_Tx_v4.16.1.csv"
    ),
    row.names = FALSE
)

p_m1m2_sham_tx <- ggplot(
    m1m2_sham_tx_df,
    aes(
        x = condition,
        y = score,
        fill = condition
    )
) +
    geom_violin(
        trim = TRUE,
        scale = "width",
        linewidth = 0.25
    ) +
    geom_boxplot(
        width = 0.12,
        outlier.shape = NA,
        alpha = 0.55
    ) +
    facet_wrap(
        ~ program,
        nrow = 1,
        scales = "free_y"
    ) +
    scale_fill_manual(
        values = c(
            Sham = "#FF3B30",
            Tx = "#3366FF"
        )
    ) +
    labs(
        title = "M1 / M2 marker programs: Sham vs Tx",
        x = NULL,
        y = "Mean normalized marker expression"
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
        ),
        legend.position = "none"
    )

save_pdf(
    file.path(
        DIR_03,
        "02_M1_M2_marker_programs_Sham_vs_Tx_v4.16.1.pdf"
    ),
    p_m1m2_sham_tx,
    8.5,
    5.2
)

# ------------------------------------------------------------------------------
# 8c. Biological sample-level M1 / M2 scores: Sham vs Tx
#
# Each point is one biological sample.
# No cell-level inferential statistics are performed here.
# ------------------------------------------------------------------------------

msg("03c_M1_M2_sample_level_Sham_vs_Tx")

m1m2_sample_level <- m1m2_sham_tx_df %>%
    group_by(
        sample,
        condition,
        program
    ) %>%
    summarise(
        n_cells = n(),
        mean_score = mean(
            score,
            na.rm = TRUE
        ),
        median_score = median(
            score,
            na.rm = TRUE
        ),
        .groups = "drop"
    )

write.csv(
    m1m2_sample_level,
    file.path(
        TAB_DIR,
        "03c_M1_M2_program_sample_level_Sham_Tx_v4.16.1.csv"
    ),
    row.names = FALSE
)

m1m2_condition_means <- m1m2_sample_level %>%
    group_by(
        condition,
        program
    ) %>%
    summarise(
        mean_of_samples = mean(
            mean_score,
            na.rm = TRUE
        ),
        .groups = "drop"
    )

p_m1m2_sample <- ggplot(
    m1m2_sample_level,
    aes(
        x = condition,
        y = mean_score
    )
) +
    geom_line(
        data = m1m2_condition_means,
        aes(
            x = condition,
            y = mean_of_samples,
            group = 1
        ),
        inherit.aes = FALSE,
        linewidth = 0.65,
        color = "black"
    ) +
    geom_point(
        data = m1m2_condition_means,
        aes(
            x = condition,
            y = mean_of_samples,
            fill = condition
        ),
        inherit.aes = FALSE,
        shape = 21,
        size = 3.6,
        stroke = 0.6,
        color = "black"
    ) +
    geom_point(
        aes(
            shape = sample
        ),
        size = 3.0,
        color = "black",
        position = position_jitter(
            width = 0.035,
            height = 0
        )
    ) +
    facet_wrap(
        ~ program,
        nrow = 1,
        scales = "free_y"
    ) +
    scale_fill_manual(
        values = c(
            Sham = "#FF3B30",
            Tx = "#3366FF"
        ),
        guide = "none"
    ) +
    scale_shape_manual(
        values = c(
            "Sham1" = 16,
            "Sham20" = 17,
            "Tx17" = 15,
            "Tx5" = 3
        ),
        drop = FALSE
    ) +
    labs(
        title = "M1 / M2 marker programs: biological sample-level Sham vs Tx",
        subtitle = "Black symbols = individual biological samples | colored circles = group means",
        x = NULL,
        y = "Sample mean marker-program score",
        shape = "Sample"
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

save_pdf(
    file.path(
        DIR_03,
        "03_M1_M2_marker_programs_sample_level_Sham_vs_Tx_v4.16.1.pdf"
    ),
    p_m1m2_sample,
    9.5,
    5.5
)

# ------------------------------------------------------------------------------
# 9. 04_AntiInflammatory_heterogeneity
# ------------------------------------------------------------------------------

msg("04_AntiInflammatory_heterogeneity")

if (file.exists(ANTI_RDS)) {

    msg("Loading Anti-inflammatory secondary-state RDS: ", ANTI_RDS)

    anti <- readRDS(ANTI_RDS)
    DefaultAssay(anti) <- ASSAY_USE

    anti <- JoinLayers(
        anti,
        assay = ASSAY_USE
    )

    if (!"data" %in%
        Layers(anti[[ASSAY_USE]])) {

        anti <- NormalizeData(
            anti,
            assay = ASSAY_USE,
            normalization.method = "LogNormalize",
            scale.factor = 10000,
            verbose = FALSE
        )
    }

    anti_data <- get_layer_safe(
        anti,
        ASSAY_USE,
        "data"
    )

    anti_meta_cols <- colnames(
        anti@meta.data
    )

    ANTI_SAMPLE_COL <- first_existing(
        anti_meta_cols,
        c(
            "sample_v4140",
            "sample_4group",
            "sample",
            "orig.ident"
        )
    )

    ANTI_CLUSTER_COL <- first_existing(
        anti_meta_cols,
        c(
            "anti_subcluster_v4140",
            "anti_subcluster",
            "seurat_clusters"
        )
    )

    ANTI_STATE_COL <- first_existing(
        anti_meta_cols,
        c(
            "anti_functional_state_v41401",
            "anti_functional_state"
        )
    )

    anti_reductions <- Reductions(anti)

    ANTI_UMAP <- first_existing(
        anti_reductions,
        c(
            "anti.umap.rpca",
            "anti.umap",
            "umap"
        )
    )

    if (is.na(ANTI_UMAP)) {
        anti_umap_like <- anti_reductions[
            grepl(
                "umap",
                anti_reductions,
                ignore.case = TRUE
            )
        ]

        if (length(anti_umap_like) > 0L) {
            ANTI_UMAP <- anti_umap_like[[1]]
        }
    }

    if (!is.na(ANTI_CLUSTER_COL) &&
        !is.null(anti_data)) {

        anti_cluster <- as.character(
            anti@meta.data[[ANTI_CLUSTER_COL]]
        )

        anti_score_df <- anti@meta.data %>%
            rownames_to_column("cell") %>%
            transmute(
                cell = cell,
                subcluster = anti_cluster
            )

        anti_program_names <- c(
            "Anti_inflammatory",
            "IL10_STAT3",
            "Repair_Resolution",
            "Inflammatory",
            "ECM_associated_inflammatory",
            "Lipid_TREM2",
            "Efferocytosis"
        )

        for (program_now in anti_program_names) {
            anti_score_df[[program_now]] <-
                mean_expression_score(
                    anti_data,
                    PROGRAMS[[program_now]]
                )
        }

        cluster_scores <- anti_score_df %>%
            pivot_longer(
                cols = all_of(
                    anti_program_names
                ),
                names_to = "program",
                values_to = "score"
            ) %>%
            group_by(
                subcluster,
                program
            ) %>%
            summarise(
                mean_score = mean(
                    score,
                    na.rm = TRUE
                ),
                .groups = "drop"
            )

        write.csv(
            cluster_scores,
            file.path(
                TAB_DIR,
                "04_AntiInflammatory_subcluster_program_scores_v4.16.1.csv"
            ),
            row.names = FALSE
        )

        heat_wide <- cluster_scores %>%
            pivot_wider(
                names_from = subcluster,
                values_from = mean_score
            )

        heat_mat <- as.matrix(
            heat_wide[
                ,
                -1,
                drop = FALSE
            ]
        )

        rownames(
            heat_mat
        ) <- heat_wide$program

        heat_z <- t(
            scale(
                t(
                    heat_mat
                )
            )
        )

        heat_z[
            !is.finite(
                heat_z
            )
        ] <- 0

        heat_plot <- pmax(
            pmin(
                heat_z,
                2
            ),
            -2
        )

        grDevices::cairo_pdf(
            file.path(
                DIR_04,
                "01_AntiInflammatory_internal_state_heatmap_v4.16.1.pdf"
            ),
            width = 8,
            height = 6
        )

        pheatmap::pheatmap(
            heat_plot,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            color = grDevices::colorRampPalette(
                c(
                    "#0033FF",
                    "#FFFFFF",
                    "#FF1A1A"
                )
            )(101),
            breaks = seq(
                -2,
                2,
                length.out = 102
            ),
            border_color = "white",
            fontsize_row = 8,
            fontsize_col = 7,
            angle_col = 45,
            main = paste0(
                "Anti-inflammatory-MΦ internal functional states\n",
                "subcluster mean score | row z-score"
            )
        )

        grDevices::dev.off()

        # Feature-program UMAP panel
        if (!is.na(ANTI_UMAP)) {

            anti_emb <- as.data.frame(
                Embeddings(
                    anti,
                    reduction = ANTI_UMAP
                )[, 1:2, drop = FALSE]
            )

            colnames(
                anti_emb
            ) <- c(
                "UMAP_1",
                "UMAP_2"
            )

            anti_plot_df <- anti_emb %>%
                rownames_to_column(
                    "cell"
                ) %>%
                left_join(
                    anti_score_df,
                    by = "cell"
                )

            program_plots <- list()

            for (program_now in anti_program_names) {

                p <- ggplot(
                    anti_plot_df,
                    aes(
                        x = UMAP_1,
                        y = UMAP_2,
                        color = .data[[program_now]]
                    )
                ) +
                    geom_point(
                        size = 0.42,
                        alpha = 0.85
                    ) +
                    scale_color_gradient2(
                        low = "#0033FF",
                        mid = "#FFFFFF",
                        high = "#FF1A1A",
                        midpoint = median(
                            anti_plot_df[[program_now]],
                            na.rm = TRUE
                        ),
                        name = "Score"
                    ) +
                    labs(
                        title = gsub(
                            "_",
                            " ",
                            program_now
                        ),
                        x = NULL,
                        y = NULL
                    ) +
                    theme_classic(
                        base_size = 8
                    ) +
                    theme(
                        plot.title = element_text(
                            face = "bold",
                            hjust = 0.5
                        )
                    )

                program_plots[[program_now]] <- p
            }

            p_program_umap <- wrap_plots(
                program_plots,
                ncol = 3
            ) +
                plot_annotation(
                    title = "Anti-inflammatory-MΦ functional programs on RPCA UMAP",
                    theme = theme(
                        plot.title = element_text(
                            face = "bold"
                        )
                    )
                )

            save_pdf(
                file.path(
                    DIR_04,
                    "02_AntiInflammatory_functional_program_UMAPs_v4.16.1.pdf"
                ),
                p_program_umap,
                14,
                13
            )
        }

        # Gene-level validation DotPlot
        genes_use <- intersect(
            ANTI_DOTPLOT_GENES,
            rownames(
                anti_data
            )
        )

        dot_rows <- list()

        for (cl in unique(
            anti_cluster
        )) {

            cells_now <- colnames(
                anti
            )[
                anti_cluster == cl
            ]

            if (length(cells_now) == 0L) next

            mat <- anti_data[
                genes_use,
                cells_now,
                drop = FALSE
            ]

            dot_rows[[cl]] <- tibble(
                subcluster = cl,
                gene = genes_use,
                avg_expr = as.numeric(
                    Matrix::rowMeans(
                        mat
                    )
                ),
                pct_expr = 100 * as.numeric(
                    Matrix::rowMeans(
                        mat > 0
                    )
                )
            )
        }

        dot_df <- bind_rows(
            dot_rows
        ) %>%
            group_by(
                gene
            ) %>%
            mutate(
                avg_expr_z = safe_z(
                    avg_expr
                )
            ) %>%
            ungroup()

        write.csv(
            dot_df,
            file.path(
                TAB_DIR,
                "05_AntiInflammatory_gene_validation_DotPlot_numeric_v4.16.1.csv"
            ),
            row.names = FALSE
        )

        p_dot <- ggplot(
            dot_df,
            aes(
                x = gene,
                y = subcluster
            )
        ) +
            geom_point(
                aes(
                    size = pct_expr,
                    color = avg_expr_z
                )
            ) +
            scale_size_continuous(
                range = c(
                    0.2,
                    7
                ),
                limits = c(
                    0,
                    100
                ),
                name = "% expressed"
            ) +
            scale_color_gradient2(
                low = "#0033FF",
                mid = "#FFFFFF",
                high = "#FF1A1A",
                midpoint = 0,
                name = "Average\nexpression\nz-score"
            ) +
            labs(
                title = "Anti-inflammatory-MΦ internal-state gene validation",
                x = NULL,
                y = "Subcluster"
            ) +
            theme_classic(
                base_size = 8
            ) +
            theme(
                plot.title = element_text(
                    face = "bold"
                ),
                axis.text.x = element_text(
                    angle = 60,
                    hjust = 1
                )
            )

        save_pdf(
            file.path(
                DIR_04,
                "03_AntiInflammatory_internal_state_gene_DotPlot_v4.16.1.pdf"
            ),
            p_dot,
            15,
            6
        )
    }

} else {
    msg("SKIP 04 secondary-state analyses: Anti RDS not found.")
}

# ------------------------------------------------------------------------------
# 10. 05_IL10_response
# ------------------------------------------------------------------------------

msg("05_IL10_response")

anti_cells_final <- meta$cell[
    meta$macrophage_class ==
        "Anti-inflammatory-Mphi"
]

if (length(anti_cells_final) > 0L) {

    anti_final_data <- rna_data[
        ,
        anti_cells_final,
        drop = FALSE
    ]

    il10_score <- mean_expression_score(
        anti_final_data,
        PROGRAMS$IL10_STAT3
    )

    il10_df <- meta %>%
        filter(
            cell %in%
                anti_cells_final
        ) %>%
        mutate(
            IL10_response = il10_score
        )

    # sample-wise median split
    il10_df <- il10_df %>%
        group_by(
            sample
        ) %>%
        mutate(
            IL10_group = ifelse(
                IL10_response >
                    median(
                        IL10_response,
                        na.rm = TRUE
                    ),
                "IL10-response-high",
                "IL10-response-low"
            )
        ) %>%
        ungroup()

    program_for_difference <- c(
        "Repair_Resolution",
        "Lipid_TREM2",
        "Inflammatory",
        "IL10_STAT3",
        "ECM_associated_inflammatory",
        "Anti_inflammatory"
    )

    score_matrix <- tibble(
        cell = anti_cells_final
    )

    for (program_now in program_for_difference) {
        score_matrix[[program_now]] <-
            mean_expression_score(
                anti_final_data,
                PROGRAMS[[program_now]]
            )
    }

    score_matrix <- score_matrix %>%
        left_join(
            il10_df %>%
                select(
                    cell,
                    sample,
                    condition,
                    IL10_group
                ),
            by = "cell"
        )

    diff_df <- score_matrix %>%
        pivot_longer(
            cols = all_of(
                program_for_difference
            ),
            names_to = "program",
            values_to = "score"
        ) %>%
        group_by(
            sample,
            condition,
            program,
            IL10_group
        ) %>%
        summarise(
            mean_score = mean(
                score,
                na.rm = TRUE
            ),
            .groups = "drop"
        ) %>%
        pivot_wider(
            names_from = IL10_group,
            values_from = mean_score
        ) %>%
        mutate(
            score_difference =
                `IL10-response-high` -
                `IL10-response-low`
        )

    write.csv(
        diff_df,
        file.path(
            TAB_DIR,
            "06_IL10_high_minus_low_program_differences_v4.16.1.csv"
        ),
        row.names = FALSE
    )

    p_il10diff <- ggplot(
        diff_df,
        aes(
            x = score_difference,
            y = program,
            shape = sample
        )
    ) +
        geom_vline(
            xintercept = 0,
            linetype = 2,
            linewidth = 0.4
        ) +
        geom_point(
            size = 2.6
        ) +
        labs(
            title = "IL10-response-high M2 minus IL10-response-low M2",
            subtitle = "Individual biological samples shown; positive = higher in IL10-high M2",
            x = "Score difference: IL10-response-high M2 minus low M2",
            y = NULL,
            shape = "Sample"
        ) +
        theme_classic(
            base_size = 9
        ) +
        theme(
            plot.title = element_text(
                face = "bold"
            )
        )

    save_pdf(
        file.path(
            DIR_05,
            "01_IL10_response_high_minus_low_M2_v4.16.1.pdf"
        ),
        p_il10diff,
        10,
        6
    )
}

# ------------------------------------------------------------------------------
# 11. 06_Pseudobulk_DE
# ------------------------------------------------------------------------------

msg("06_Pseudobulk_DE")

pb_rows <- list()

for (class_now in CLASS_ORDER[
    CLASS_ORDER != "Other"
]) {

    for (sample_now in SAMPLE_ORDER) {

        cells_now <- meta$cell[
            meta$macrophage_class ==
                class_now &
            meta$sample ==
                sample_now
        ]

        if (length(cells_now) == 0L) next

        pb <- Matrix::rowSums(
            counts[
                ,
                cells_now,
                drop = FALSE
            ]
        )

        lib <- sum(pb)

        pb_rows[[paste(
            class_now,
            sample_now,
            sep = "__"
        )]] <- tibble(
            macrophage_class = class_now,
            sample = sample_now,
            condition = as.character(
                canonical_condition(
                    sample_now
                )
            ),
            gene = rownames(
                counts
            ),
            raw_count = as.numeric(
                pb
            ),
            CPM = as.numeric(
                pb
            ) / lib * 1e6,
            logCPM = log2(
                as.numeric(
                    pb
                ) / lib * 1e6 + 1
            )
        )
    }
}

pb_long <- bind_rows(
    pb_rows
)

write.csv(
    pb_long,
    file.path(
        TAB_DIR,
        "07_FINAL_sample_level_pseudobulk_all_genes_v4.16.1.csv"
    ),
    row.names = FALSE
)

effect_rows <- list()

for (class_now in CLASS_ORDER[
    CLASS_ORDER != "Other"
]) {

    dat <- pb_long %>%
        filter(
            macrophage_class ==
                class_now,
            sample %in%
                c(
                    "Sham1",
                    "Sham20",
                    "Tx17",
                    "Tx5"
                )
        ) %>%
        select(
            gene,
            sample,
            logCPM
        ) %>%
        pivot_wider(
            names_from = sample,
            values_from = logCPM
        )

    needed <- c(
        "Sham1",
        "Sham20",
        "Tx17",
        "Tx5"
    )

    if (!all(
        needed %in%
            colnames(
                dat
            )
    )) next

    effect_rows[[class_now]] <- dat %>%
        mutate(
            Sham_mean =
                rowMeans(
                    cbind(
                        Sham1,
                        Sham20
                    ),
                    na.rm = TRUE
                ),

            Tx_mean =
                rowMeans(
                    cbind(
                        Tx17,
                        Tx5
                    ),
                    na.rm = TRUE
                ),

            log2FC_Tx_vs_Sham =
                Tx_mean -
                Sham_mean,

            log2FC_Tx17_vs_Sham =
                Tx17 -
                Sham_mean,

            log2FC_Tx5_vs_Sham =
                Tx5 -
                Sham_mean,

            concordant =
                sign(
                    log2FC_Tx17_vs_Sham
                ) ==
                sign(
                    log2FC_Tx5_vs_Sham
                ),

            macrophage_class =
                class_now
        )
}

effect_all <- bind_rows(
    effect_rows
)

write.csv(
    effect_all,
    file.path(
        TAB_DIR,
        "08_FINAL_pseudobulk_effects_Tx_vs_Sham_v4.16.1.csv"
    ),
    row.names = FALSE
)

top_effect <- effect_all %>%
    filter(
        concordant
    ) %>%
    group_by(
        macrophage_class
    ) %>%
    arrange(
        desc(
            abs(
                log2FC_Tx_vs_Sham
            )
        ),
        .by_group = TRUE
    ) %>%
    slice_head(
        n = 20
    ) %>%
    ungroup()

write.csv(
    top_effect,
    file.path(
        TAB_DIR,
        "09_FINAL_top20_concordant_pseudobulk_effects_v4.16.1.csv"
    ),
    row.names = FALSE
)

p_effect <- ggplot(
    top_effect,
    aes(
        x = log2FC_Tx_vs_Sham,
        y = reorder(
            gene,
            log2FC_Tx_vs_Sham
        ),
        shape = log2FC_Tx_vs_Sham > 0
    )
) +
    geom_vline(
        xintercept = 0,
        linetype = 2,
        linewidth = 0.4
    ) +
    geom_segment(
        aes(
            x = 0,
            xend = log2FC_Tx_vs_Sham,
            yend = gene
        ),
        linewidth = 0.45
    ) +
    geom_point(
        size = 2.1
    ) +
    facet_wrap(
        ~ macrophage_class,
        scales = "free_y",
        ncol = 3,
        labeller = as_labeller(
            CLASS_LABELS
        )
    ) +
    labs(
        title = "Clean-B FINAL MΦ subtype pseudobulk: Tx vs Sham",
        subtitle = "Top effect sizes | Tx17/Tx5 directionally concordant",
        x = "log2FC (Tx / Sham)",
        y = NULL,
        shape = "Direction"
    ) +
    theme_classic(
        base_size = 8
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        strip.text = element_text(
            face = "bold"
        )
    )

save_pdf(
    file.path(
        DIR_06,
        "01_FINAL_subtype_pseudobulk_Tx_vs_Sham_top_effects_v4.16.1.pdf"
    ),
    p_effect,
    12,
    12
)

p_rep <- ggplot(
    effect_all,
    aes(
        x = log2FC_Tx17_vs_Sham,
        y = log2FC_Tx5_vs_Sham
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = 2,
        linewidth = 0.35
    ) +
    geom_vline(
        xintercept = 0,
        linetype = 2,
        linewidth = 0.35
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = 3,
        linewidth = 0.35
    ) +
    geom_point(
        size = 0.65,
        alpha = 0.45
    ) +
    facet_wrap(
        ~ macrophage_class,
        scales = "free",
        ncol = 3,
        labeller = as_labeller(
            CLASS_LABELS
        )
    ) +
    labs(
        title = "Tx replicate consistency by MΦ subtype",
        subtitle = "Each gene: Tx sample minus mean Sham log2(CPM+1)",
        x = "Tx17 − mean Sham",
        y = "Tx5 − mean Sham"
    ) +
    theme_classic(
        base_size = 8
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        strip.text = element_text(
            face = "bold"
        )
    )

save_pdf(
    file.path(
        DIR_06,
        "02_FINAL_Tx_replicate_consistency_by_subtype_v4.16.1.pdf"
    ),
    p_rep,
    12,
    9
)

# ------------------------------------------------------------------------------
# 12. 07_Pathway_enrichment
# ------------------------------------------------------------------------------

msg("07_Pathway_enrichment")

if (file.exists(PATHWAY_TABLE)) {

    pathway <- read.csv(
        PATHWAY_TABLE,
        check.names = FALSE
    )

    if (all(
        c(
            "macrophage_class",
            "pathway",
            "NES"
        ) %in%
            colnames(
                pathway
            )
    )) {

        pathway$macrophage_class <- recode(
            pathway$macrophage_class,
            "Fibrogenic-Mphi" =
                "ECM-associated inflammatory-Mphi"
        )

        pathway$pathway_label <- if (
            "pathway_label" %in%
                colnames(
                    pathway
                )
        ) {
            pathway$pathway_label
        } else {
            gsub(
                "_",
                " ",
                pathway$pathway
            )
        }

        mechanism_regex <- paste(
            c(
                "TNFA",
                "NFKB",
                "INFLAMMATORY",
                "IL6",
                "JAK",
                "STAT3",
                "TGF",
                "FATTY_ACID",
                "CHOLESTEROL",
                "BILE_ACID",
                "PPAR",
                "OXIDATIVE",
                "PEROXISOME",
                "HYPOXIA",
                "ANGIOGEN",
                "APOPTOSIS",
                "INTERFERON",
                "EXTRACELLULAR_MATRIX",
                "ECM"
            ),
            collapse = "|"
        )

        pathway_focus <- pathway %>%
            filter(
                grepl(
                    mechanism_regex,
                    pathway,
                    ignore.case = TRUE
                ),
                macrophage_class %in%
                    CLASS_ORDER[
                        CLASS_ORDER != "Other"
                    ]
            ) %>%
            group_by(
                macrophage_class
            ) %>%
            arrange(
                if ("padj" %in% colnames(.)) {
                    padj
                } else {
                    desc(
                        abs(
                            NES
                        )
                    )
                },
                .by_group = TRUE
            ) %>%
            slice_head(
                n = 18
            ) %>%
            ungroup()

        write.csv(
            pathway_focus,
            file.path(
                TAB_DIR,
                "10_FINAL_mechanism_focused_pathways_v4.16.1.csv"
            ),
            row.names = FALSE
        )

        selected_pathways <- unique(
            pathway_focus$pathway
        )

        heat <- pathway %>%
            filter(
                pathway %in%
                    selected_pathways,
                macrophage_class %in%
                    CLASS_ORDER[
                        CLASS_ORDER != "Other"
                    ]
            ) %>%
            select(
                pathway,
                pathway_label,
                macrophage_class,
                NES
            ) %>%
            distinct() %>%
            pivot_wider(
                names_from = macrophage_class,
                values_from = NES
            )

        if (nrow(heat) > 0L) {

            cols_use <- intersect(
                CLASS_ORDER[
                    CLASS_ORDER != "Other"
                ],
                colnames(
                    heat
                )
            )

            mat <- as.matrix(
                heat[
                    ,
                    cols_use,
                    drop = FALSE
                ]
            )

            rownames(
                mat
            ) <- make.unique(
                heat$pathway_label
            )

            mat[
                !is.finite(
                    mat
                )
            ] <- 0

            lim <- max(
                2,
                quantile(
                    abs(
                        mat
                    ),
                    0.95,
                    na.rm = TRUE
                )
            )

            mat_plot <- pmax(
                pmin(
                    mat,
                    lim
                ),
                -lim
            )

            grDevices::cairo_pdf(
                file.path(
                    DIR_07,
                    "01_FINAL_mechanism_focused_pathway_NES_heatmap_v4.16.1.pdf"
                ),
                width = 9,
                height = 12
            )

            pheatmap::pheatmap(
                mat_plot,
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                color = grDevices::colorRampPalette(
                    c(
                        "#0033FF",
                        "#FFFFFF",
                        "#FF1A1A"
                    )
                )(101),
                breaks = seq(
                    -lim,
                    lim,
                    length.out = 102
                ),
                border_color = "white",
                fontsize_row = 6.5,
                fontsize_col = 8,
                angle_col = 45,
                main = paste0(
                    "MΦ subtype pathway remodeling: Sham → Tx\n",
                    "NES > 0 = enriched in Tx"
                )
            )

            grDevices::dev.off()
        }
    }

} else {
    msg("SKIP 07 pathway enrichment: pathway table not found.")
}

# ------------------------------------------------------------------------------
# 13. 08_Robustness
# ------------------------------------------------------------------------------

msg("08_Robustness")

robustness_files <- c(
    composition =
        file.path(
            ROBUSTNESS_DIR,
            "01_composition_sensitivity_all_cleaning_conditions_v4.14.3.csv"
        ),

    contamination =
        file.path(
            ROBUSTNESS_DIR,
            "03_contamination_candidate_burden_by_subtype_v4.14.3.csv"
        ),

    purity =
        file.path(
            ROBUSTNESS_DIR,
            "04_Res2_cluster_class_purity_entropy_CleanB_v4.14.3.csv"
        ),

    specificity =
        file.path(
            ROBUSTNESS_DIR,
            "07_expected_program_specificity_margin_v4.14.3.csv"
        )
)

robustness_present <- file.exists(
    robustness_files
)

write.csv(
    tibble(
        item = names(
            robustness_files
        ),
        file = unname(
            robustness_files
        ),
        present = robustness_present
    ),
    file.path(
        TAB_DIR,
        "11_robustness_input_audit_v4.16.1.csv"
    ),
    row.names = FALSE
)

# Composition sensitivity plot
if (file.exists(
    robustness_files[["composition"]]
)) {

    comp_sens <- read.csv(
        robustness_files[["composition"]],
        check.names = FALSE
    )

    comp_sens$macrophage_class <- recode(
        comp_sens$macrophage_class,
        "Fibrogenic-Mphi" =
            "ECM-associated inflammatory-Mphi"
    )

    p_rob_comp <- ggplot(
        comp_sens,
        aes(
            x = analysis,
            y = percent,
            group = sample,
            shape = sample
        )
    ) +
        geom_line(
            linewidth = 0.45,
            alpha = 0.65
        ) +
        geom_point(
            size = 2.2
        ) +
        facet_wrap(
            ~ macrophage_class,
            scales = "free_y",
            ncol = 3,
            labeller = as_labeller(
                CLASS_LABELS
            )
        ) +
        labs(
            title = "MΦ subtype abundance robustness to cleaning",
            subtitle = "Original / Clean-A / Clean-B / Clean-C",
            x = NULL,
            y = "% retained MΦ",
            shape = "Sample"
        ) +
        theme_classic(
            base_size = 8
        ) +
        theme(
            plot.title = element_text(
                face = "bold"
            ),
            strip.text = element_text(
                face = "bold"
            )
        )

    save_pdf(
        file.path(
            DIR_08,
            "01_FINAL_composition_cleaning_sensitivity_v4.16.1.pdf"
        ),
        p_rob_comp,
        11,
        8
    )
}

# contamination
if (file.exists(
    robustness_files[["contamination"]]
)) {

    contam <- read.csv(
        robustness_files[["contamination"]],
        check.names = FALSE
    )

    contam$macrophage_class <- recode(
        contam$macrophage_class,
        "Fibrogenic-Mphi" =
            "ECM-associated inflammatory-Mphi"
    )

    contam_long <- contam %>%
        select(
            macrophage_class,
            ends_with("_pct")
        ) %>%
        pivot_longer(
            cols = -macrophage_class,
            names_to = "QC_type",
            values_to = "percent"
        )

    p_contam <- ggplot(
        contam_long,
        aes(
            x = macrophage_class,
            y = percent,
            shape = QC_type
        )
    ) +
        geom_point(
            size = 2.4,
            position = position_dodge(
                width = 0.35
            )
        ) +
        labs(
            title = "Contamination-candidate burden by FINAL MΦ subtype",
            x = NULL,
            y = "% cells flagged",
            shape = "QC flag"
        ) +
        theme_classic(
            base_size = 8
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
        file.path(
            DIR_08,
            "02_FINAL_contamination_candidate_burden_v4.16.1.pdf"
        ),
        p_contam,
        9,
        5.5
    )
}

# ------------------------------------------------------------------------------
# 14. Figure index
# ------------------------------------------------------------------------------

figure_files <- list.files(
    OUTPUT_DIR,
    pattern = "\\.pdf$",
    recursive = TRUE,
    full.names = TRUE
)

figure_index <- tibble(
    file = figure_files,
    relative_path = sub(
        paste0(
            "^",
            gsub(
                "([.()+*?^$|{}\\[\\]\\\\])",
                "\\\\\\1",
                OUTPUT_DIR
            ),
            "/?"
        ),
        "",
        figure_files
    )
) %>%
    arrange(
        relative_path
    )

write.csv(
    figure_index,
    file.path(
        OUTPUT_DIR,
        "FINAL_EXTENDED_FIGURE_INDEX_v4.16.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 15. README / session
# ------------------------------------------------------------------------------

readme <- c(
    "Mouse MASH MΦ FINAL extended figures v4.16.1",
    "",
    paste0(
        "FINAL RDS: ",
        FINAL_RDS
    ),
    "",
    paste0(
        "FINAL annotation: ",
        FINAL_CLASS_COL
    ),
    "",
    "Final parent classification:",
    "  Inflammatory-MΦ",
    "  Anti-inflammatory-MΦ",
    "  ECM-associated inflammatory-MΦ",
    "  Repair/Resolution-MΦ",
    "  Lipid-associated/TREM2-MΦ",
    "  Other",
    "",
    "Final UMAP style:",
    "  fixed window x=-8.8..-2.0, y=-7.1..-0.5",
    "  major dots=1.20, Other=0.60",
    "  Inflammatory=red",
    "  Anti-inflammatory=cyan",
    "  ECM-associated inflammatory=orange",
    "  Repair/Resolution=green",
    "  Lipid/TREM2=purple",
    "  Other=gray",
    "",
    "Sections:",
    "  01 UMAP",
    "  02 Composition / M2 metrics",
    "  03 M1/M2 programs (STD vs CDAHFD; Sham vs Tx; sample-level Sham vs Tx)",
    "  04 Anti-inflammatory internal heterogeneity",
    "  05 IL10-response-high vs low M2",
    "  06 Sample-level pseudobulk and Tx replicate consistency",
    "  07 Pathway enrichment",
    "  08 Classification robustness",
    "",
    "Interpretation:",
    "  STD vs CDAHFD is descriptive (n=1 vs n=1).",
    "  Sham vs Tx summaries use biological samples (n=2 vs n=2).",
    "",
    "Optional inputs are skipped if absent."
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_FINAL_EXTENDED_v4.16.1.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.16.1.txt"
    )
)

# ------------------------------------------------------------------------------
# 16. Final
# ------------------------------------------------------------------------------

msg("DONE.")
msg("Output: ", OUTPUT_DIR)
msg(
    "Figure index: ",
    file.path(
        OUTPUT_DIR,
        "FINAL_EXTENDED_FIGURE_INDEX_v4.16.1.csv"
    )
)

print(
    figure_index
)
