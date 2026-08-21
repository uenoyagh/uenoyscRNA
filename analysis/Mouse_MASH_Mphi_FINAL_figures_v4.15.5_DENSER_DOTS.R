#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# FINAL figure regeneration from Clean-B FINAL annotation
# v4.15.5
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4150)

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

OUTPUT_DIR <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "Final_Figures_v4.15.5"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"
FINAL_CLASS_COL <- "macrophage_class_Res2_FINAL_v4145"

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

CONDITION_ORDER <- c("STD", "CDAHFD", "Sham", "Tx")

SAMPLE_ORDER <- c(
    "STD_rep1",
    "CDHFD_rep1",
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
)

CLASS_COLORS <- c(
    "Inflammatory-Mphi" = "#FF2D2D",
    "Anti-inflammatory-Mphi" = "#00AEEF",
    "ECM-associated inflammatory-Mphi" = "#FF8C00",
    "Repair/Resolution-Mphi" = "#00C853",
    "Lipid-associated/TREM2-Mphi" = "#7B2CBF",
    "Other" = "#808080"
)

CONDITION_COLORS <- c(
    "STD" = "#0066FF",
    "CDAHFD" = "#FF1A1A",
    "Sham" = "#FF1A1A",
    "Tx" = "#0066FF"
)


# ------------------------------------------------------------------------------
# UMAP display parameters
#
# IMPORTANT:
#   These settings affect DISPLAY ONLY.
#   No cells are removed from the Seurat object or downstream quantification.
#
#   Fixed plotting window intentionally excludes the sparse upper-right
#   "Other" tail from the publication UMAP while preserving all cells in the
#   Seurat object and all quantitative analyses.
# ------------------------------------------------------------------------------

UMAP_XLIM <- c(
    -8.8,
    -2.0
)

UMAP_YLIM <- c(
    -7.1,
    -0.5
)

UMAP_PT_SIZE_OTHER <- 0.60
UMAP_PT_SIZE_MAJOR <- 1.20

UMAP_ALPHA_OTHER <- 0.42
UMAP_ALPHA_MAJOR <- 0.95

UMAP_LEGEND_POINT_SIZE <- 4.4

IL10_RESPONSE_GENES <- c(
    "Il10ra","Il10rb","Jak1","Tyk2","Stat3","Socs3","Bcl3","Il1rn"
)

GENE_BLOCKS <- list(
    Inflammatory = c("Il1b","Tnf","Ccl2","Cxcl10"),
    Anti_inflammatory = c("Mrc1","Cd163","Il1rn","Mertk","Igf1","Hmox1"),
    IL10_STAT3 = c("Stat3","Socs3"),
    ECM_inflammatory = c("Thbs1","Fn1","Tgfb1","Col1a1","Col1a2","Col3a1"),
    Repair_Resolution = c("Mfge8","Gas6","Mmp13","Mmp14","Plau"),
    Lipid_TREM2 = c("Trem2","Gpnmb","Cd9","Lpl","Apoe")
)

GENE_ORDER <- unique(unlist(GENE_BLOCKS, use.names = FALSE))

required_packages <- c(
    "Seurat","SeuratObject","Matrix","dplyr","tidyr","tibble",
    "ggplot2","patchwork","scales","pheatmap"
)

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
    stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
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

msg <- function(...) {
    message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(...))
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
        LayerData(object, assay = assay, layer = layer),
        error = function(e) NULL
    )
    if (is.null(x)) {
        x <- tryCatch(
            GetAssayData(object, assay = assay, slot = layer),
            error = function(e) NULL
        )
    }
    x
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

mean_expression_score <- function(mat, genes) {
    genes <- intersect(genes, rownames(mat))
    if (length(genes) == 0L) return(rep(NA_real_, ncol(mat)))
    as.numeric(Matrix::colMeans(mat[genes, , drop = FALSE]))
}


# ------------------------------------------------------------------------------
# 5. Load FINAL Clean-B MΦ object
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
    stop(
        "FINAL RDS not found:\n",
        INPUT_RDS
    )
}

msg(
    "Loading FINAL RDS: ",
    INPUT_RDS
)

mphi <- readRDS(
    INPUT_RDS
)

if (!inherits(
    mphi,
    "Seurat"
)) {
    stop(
        "Input object is not a Seurat object."
    )
}

if (!FINAL_CLASS_COL %in%
    colnames(
        mphi@meta.data
    )) {

    stop(
        "FINAL annotation column missing: ",
        FINAL_CLASS_COL
    )
}

DefaultAssay(
    mphi
) <- ASSAY_USE


# ------------------------------------------------------------------------------
# 6. Resolve sample column and UMAP reduction
# ------------------------------------------------------------------------------

meta_cols <- colnames(
    mphi@meta.data
)

SAMPLE_COL <- first_existing(
    meta_cols,
    c(
        "sample_4group",
        "sample",
        "sample_id",
        "orig.ident"
    )
)

if (is.na(
    SAMPLE_COL
)) {
    stop(
        "Sample metadata column not found."
    )
}

available_reductions <- Reductions(
    mphi
)

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

if (is.na(
    REDUCTION_USE
)) {

    umap_like <- available_reductions[
        grepl(
            "umap",
            available_reductions,
            ignore.case = TRUE
        )
    ]

    if (length(
        umap_like
    ) == 0L) {
        stop(
            "No UMAP-like reduction found."
        )
    }

    REDUCTION_USE <- umap_like[[1]]
}

msg(
    "Sample column: ",
    SAMPLE_COL
)

msg(
    "UMAP reduction: ",
    REDUCTION_USE
)


# ------------------------------------------------------------------------------
# 7. FINAL metadata used by all figures
# ------------------------------------------------------------------------------

mphi$sample_FINAL_v4150 <- as.character(
    mphi@meta.data[[SAMPLE_COL]]
)

mphi$condition_FINAL_v4150 <- canonical_condition(
    mphi$sample_FINAL_v4150
)

mphi$class_FINAL_v4150 <- factor(
    as.character(
        mphi@meta.data[[FINAL_CLASS_COL]]
    ),
    levels = CLASS_ORDER
)

if (any(
    is.na(
        mphi$class_FINAL_v4150
    )
)) {
    stop(
        "NA generated while constructing FINAL macrophage class factor."
    )
}

meta <- mphi@meta.data %>%
    rownames_to_column(
        "cell"
    ) %>%
    transmute(
        cell = cell,
        sample = sample_FINAL_v4150,
        condition = as.character(
            condition_FINAL_v4150
        ),
        macrophage_class = as.character(
            class_FINAL_v4150
        )
    )


# ------------------------------------------------------------------------------
# 8. FINAL class-count audit
# ------------------------------------------------------------------------------

final_count_table <- meta %>%
    count(
        macrophage_class,
        name = "n_cells"
    ) %>%
    mutate(
        fraction =
            n_cells /
            sum(
                n_cells
            ),

        percent =
            100 *
            fraction
    )

write.csv(
    final_count_table,
    file.path(
        TAB_DIR,
        "00_FINAL_Mphi_class_counts_v4.15.5.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 11. UMAP display limits
#
# Fixed publication window.
# DISPLAY ONLY: all cells remain in the Seurat object and downstream analyses.
# ------------------------------------------------------------------------------

msg(
    "UMAP fixed display x-range: ",
    paste(
        UMAP_XLIM,
        collapse = " to "
    )
)

msg(
    "UMAP fixed display y-range: ",
    paste(
        UMAP_YLIM,
        collapse = " to "
    )
)

write.csv(
    tibble(
        axis = c(
            "x",
            "y"
        ),
        min = c(
            UMAP_XLIM[1],
            UMAP_YLIM[1]
        ),
        max = c(
            UMAP_XLIM[2],
            UMAP_YLIM[2]
        )
    ),
    file.path(
        TAB_DIR,
        "00b_FINAL_UMAP_display_limits_v4.15.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. Prepare UMAP data
#
# Other is drawn first as a background layer.
# Major five MΦ classes are drawn afterward.
# ------------------------------------------------------------------------------

umap_emb <- as.data.frame(
    Embeddings(
        mphi,
        reduction = REDUCTION_USE
    )[
        ,
        1:2,
        drop = FALSE
    ]
)

colnames(
    umap_emb
) <- c(
    "UMAP_1",
    "UMAP_2"
)

umap_df <- umap_emb %>%
    rownames_to_column(
        "cell"
    ) %>%
    left_join(
        meta %>%
            select(
                cell,
                sample,
                condition,
                macrophage_class
            ),
        by = "cell"
    ) %>%
    mutate(
        macrophage_class =
            factor(
                macrophage_class,
                levels =
                    CLASS_ORDER
            ),

        condition =
            factor(
                condition,
                levels =
                    CONDITION_ORDER
            )
    )

umap_other <- umap_df %>%
    filter(
        macrophage_class ==
            "Other"
    )

umap_major <- umap_df %>%
    filter(
        macrophage_class !=
            "Other"
    )

# ------------------------------------------------------------------------------
# 13. Figure 01: FINAL MΦ subtype UMAP
# ------------------------------------------------------------------------------

p_umap <- ggplot() +

    geom_point(
        data =
            umap_other,
        aes(
            x =
                UMAP_1,
            y =
                UMAP_2,
            color =
                macrophage_class
        ),
        size =
            UMAP_PT_SIZE_OTHER,
        alpha =
            UMAP_ALPHA_OTHER,
        stroke =
            0
    ) +

    geom_point(
        data =
            umap_major,
        aes(
            x =
                UMAP_1,
            y =
                UMAP_2,
            color =
                macrophage_class
        ),
        size =
            UMAP_PT_SIZE_MAJOR,
        alpha =
            UMAP_ALPHA_MAJOR,
        stroke =
            0
    ) +

    scale_color_manual(
        values =
            CLASS_COLORS,
        breaks =
            CLASS_ORDER,
        labels =
            CLASS_LABELS[
                CLASS_ORDER
            ],
        drop =
            FALSE
    ) +
    guides(
        color =
            guide_legend(
                override.aes = list(
                    size =
                        UMAP_LEGEND_POINT_SIZE,
                    alpha =
                        1
                )
            )
    ) +

    coord_cartesian(
        xlim =
            UMAP_XLIM,
        ylim =
            UMAP_YLIM,
        expand =
            FALSE
    ) +

    labs(
        title =
            "Mouse MASH macrophages",
        subtitle =
            "Clean-B FINAL annotation | Res2.0",
        x =
            "UMAP 1",
        y =
            "UMAP 2",
        color =
            "MΦ subtype"
    ) +

    theme_classic(
        base_size =
            11
    ) +

    theme(
        plot.title =
            element_text(
                face =
                    "bold",
                size =
                    14
            ),

        axis.title =
            element_text(
                size =
                    11
            ),

        axis.text =
            element_text(
                size =
                    9
            ),

        legend.title =
            element_text(
                face =
                    "bold",
                size =
                    10
            ),

        legend.text =
            element_text(
                size =
                    9.5
            ),

        legend.key.height =
            grid::unit(
                0.55,
                "cm"
            )
    )

save_pdf(
    "01_FINAL_Mphi_subtype_UMAP_v4.15.5.pdf",
    p_umap,
    9.5,
    7
)

# ------------------------------------------------------------------------------
# 14. Figure 02: condition-split UMAP
#
# Same x/y limits as Figure 01.
# ------------------------------------------------------------------------------

p_condition_umap <- ggplot() +

    geom_point(
        data =
            umap_other,
        aes(
            x =
                UMAP_1,
            y =
                UMAP_2,
            color =
                macrophage_class
        ),
        size =
            UMAP_PT_SIZE_OTHER,
        alpha =
            UMAP_ALPHA_OTHER,
        stroke =
            0
    ) +

    geom_point(
        data =
            umap_major,
        aes(
            x =
                UMAP_1,
            y =
                UMAP_2,
            color =
                macrophage_class
        ),
        size =
            UMAP_PT_SIZE_MAJOR,
        alpha =
            UMAP_ALPHA_MAJOR,
        stroke =
            0
    ) +

    facet_wrap(
        ~ condition,
        ncol =
            2,
        drop =
            FALSE
    ) +

    scale_color_manual(
        values =
            CLASS_COLORS,
        breaks =
            CLASS_ORDER,
        labels =
            CLASS_LABELS[
                CLASS_ORDER
            ],
        drop =
            FALSE
    ) +
    guides(
        color =
            guide_legend(
                override.aes = list(
                    size =
                        UMAP_LEGEND_POINT_SIZE,
                    alpha =
                        1
                )
            )
    ) +

    coord_cartesian(
        xlim =
            UMAP_XLIM,
        ylim =
            UMAP_YLIM,
        expand =
            FALSE
    ) +

    labs(
        title =
            "FINAL MΦ subtype distribution by condition",
        subtitle =
            "Identical UMAP coordinates and display limits across conditions",
        x =
            "UMAP 1",
        y =
            "UMAP 2",
        color =
            "MΦ subtype"
    ) +

    theme_classic(
        base_size =
            9
    ) +

    theme(
        plot.title =
            element_text(
                face =
                    "bold"
            ),

        strip.text =
            element_text(
                face =
                    "bold",
                size =
                    10
            ),

        legend.title =
            element_text(
                face =
                    "bold"
            ),

        legend.text =
            element_text(
                size =
                    8.5
            )
    )

save_pdf(
    "02_FINAL_Mphi_subtype_UMAP_by_condition_v4.15.5.pdf",
    p_condition_umap,
    12,
    9
)

# Composition
sample_totals <- meta %>%
    count(sample, name = "total_Mphi")

composition <- meta %>%
    count(sample, condition, macrophage_class, name = "n_cells") %>%
    complete(
        sample = SAMPLE_ORDER,
        macrophage_class = CLASS_ORDER,
        fill = list(n_cells = 0L)
    ) %>%
    mutate(
        condition = as.character(canonical_condition(sample))
    ) %>%
    left_join(sample_totals, by = "sample") %>%
    mutate(
        fraction = n_cells / total_Mphi,
        percent = 100 * fraction,
        sample = factor(sample, levels = SAMPLE_ORDER),
        condition = factor(condition, levels = CONDITION_ORDER),
        macrophage_class = factor(macrophage_class, levels = CLASS_ORDER)
    )

write.csv(
    composition,
    file.path(TAB_DIR, "01_FINAL_subtype_composition_by_sample_v4.15.0.csv"),
    row.names = FALSE
)

# 03 sample composition
p_comp_sample <- ggplot(
    composition,
    aes(x = sample, y = percent, fill = macrophage_class)
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
        axis.text.x = element_text(angle = 35, hjust = 1)
    )

save_pdf(
    "03_FINAL_Mphi_subtype_composition_by_sample_v4.15.0.pdf",
    p_comp_sample,
    10,
    6
)

make_pair_change_plot <- function(comp, condition_a, condition_b, title_text, filename) {
    df <- comp %>%
        filter(condition %in% c(condition_a, condition_b)) %>%
        mutate(condition = factor(condition, levels = c(condition_a, condition_b)))

    p <- ggplot(
        df,
        aes(x = condition, y = percent, group = sample)
    ) +
        geom_line(linewidth = 0.55, alpha = 0.65) +
        geom_point(
            aes(color = macrophage_class),
            size = 2.8
        ) +
        facet_wrap(
            ~ macrophage_class,
            scales = "free_y",
            ncol = 3,
            labeller = as_labeller(CLASS_LABELS)
        ) +
        scale_color_manual(values = CLASS_COLORS, guide = "none") +
        labs(
            title = title_text,
            subtitle = "Each point = biological sample",
            x = NULL,
            y = "% of Clean-B MΦ"
        ) +
        theme_classic(base_size = 9) +
        theme(
            plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold", size = 8)
        )

    save_pdf(filename, p, 11, 8)
    invisible(p)
}

# 04 disease
make_pair_change_plot(
    composition,
    "STD",
    "CDAHFD",
    "MΦ subtype remodeling: STD → CDAHFD",
    "04_FINAL_Mphi_subtype_change_STD_to_CDAHFD_v4.15.0.pdf"
)

# 05 treatment
make_pair_change_plot(
    composition,
    "Sham",
    "Tx",
    "MΦ subtype remodeling: Sham → Tx",
    "05_FINAL_Mphi_subtype_change_Sham_to_Tx_v4.15.0.pdf"
)

# M1/M2 metrics
metric_counts <- meta %>%
    count(sample, condition, macrophage_class, name = "n_cells") %>%
    complete(
        sample = SAMPLE_ORDER,
        macrophage_class = CLASS_ORDER,
        fill = list(n_cells = 0L)
    ) %>%
    mutate(condition = as.character(canonical_condition(sample))) %>%
    pivot_wider(
        names_from = macrophage_class,
        values_from = n_cells,
        values_fill = 0
    ) %>%
    left_join(sample_totals, by = "sample") %>%
    mutate(
        M1_n = .data[["Inflammatory-Mphi"]],
        M2_n = .data[["Anti-inflammatory-Mphi"]],
        AntiInflammatory_fraction_total_Mphi = M2_n / total_Mphi,
        AntiInflammatory_percent_total_Mphi = 100 * AntiInflammatory_fraction_total_Mphi,
        M2_fraction_within_M1_M2 = M2_n / (M1_n + M2_n),
        M2_percent_within_M1_M2 = 100 * M2_fraction_within_M1_M2,
        M2_to_M1_ratio = (M2_n + 0.5) / (M1_n + 0.5),
        log2_M2_to_M1_ratio = log2(M2_to_M1_ratio),
        sample = factor(sample, levels = SAMPLE_ORDER),
        condition = factor(condition, levels = CONDITION_ORDER)
    )

write.csv(
    metric_counts,
    file.path(TAB_DIR, "02_FINAL_M2_metrics_by_sample_v4.15.0.csv"),
    row.names = FALSE
)

make_metric_plot <- function(df, y_col, y_label, title_text, filename) {
    p <- ggplot(
        df,
        aes(x = condition, y = .data[[y_col]])
    ) +
        geom_point(
            aes(shape = sample, color = condition),
            size = 3.2
        ) +
        scale_color_manual(values = CONDITION_COLORS, drop = FALSE) +
        labs(
            title = title_text,
            subtitle = "Each point = biological sample",
            x = NULL,
            y = y_label,
            shape = "Sample",
            color = "Condition"
        ) +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(face = "bold"))

    save_pdf(filename, p, 8, 5.5)
    invisible(p)
}

# 06 Anti-inflammatory / total MΦ
make_metric_plot(
    metric_counts,
    "AntiInflammatory_percent_total_Mphi",
    "Anti-inflammatory-MΦ / total MΦ (%)",
    "FINAL Anti-inflammatory-MΦ fraction",
    "06_FINAL_AntiInflammatory_fraction_total_Mphi_v4.15.0.pdf"
)

# 07 M2 fraction within strict M1+M2
make_metric_plot(
    metric_counts,
    "M2_percent_within_M1_M2",
    "Anti-inflammatory-MΦ / (Inflammatory + Anti-inflammatory) (%)",
    "FINAL M2 fraction within M1+M2",
    "07_FINAL_M2_fraction_within_M1_M2_v4.15.0.pdf"
)

# 08 M2/M1 ratio
make_metric_plot(
    metric_counts,
    "M2_to_M1_ratio",
    "Anti-inflammatory-MΦ / Inflammatory-MΦ ratio",
    "FINAL M2/M1 ratio",
    "08_FINAL_M2_to_M1_ratio_v4.15.0.pdf"
)

# RNA layer / IL10 response
mphi <- JoinLayers(mphi, assay = ASSAY_USE)

if (!"data" %in% Layers(mphi[[ASSAY_USE]])) {
    mphi <- NormalizeData(
        mphi,
        assay = ASSAY_USE,
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
    )
}

rna_data <- get_layer_safe(mphi, ASSAY_USE, "data")
if (is.null(rna_data)) stop("Normalized RNA data layer not available.")

IL10_USE <- intersect(IL10_RESPONSE_GENES, rownames(rna_data))

write.csv(
    tibble(
        gene = IL10_RESPONSE_GENES,
        detected = IL10_RESPONSE_GENES %in% rownames(rna_data)
    ),
    file.path(TAB_DIR, "03_FINAL_IL10_response_gene_audit_v4.15.0.csv"),
    row.names = FALSE
)

mphi$IL10_response_FINAL_v4150 <- mean_expression_score(rna_data, IL10_USE)

il10_sample <- mphi@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = sample_FINAL_v4150,
        condition = as.character(condition_FINAL_v4150),
        macrophage_class = as.character(class_FINAL_v4150),
        IL10_response = IL10_response_FINAL_v4150
    ) %>%
    group_by(sample, condition, macrophage_class) %>%
    summarise(
        n_cells = n(),
        mean_IL10_response = mean(IL10_response, na.rm = TRUE),
        median_IL10_response = median(IL10_response, na.rm = TRUE),
        .groups = "drop"
    )

write.csv(
    il10_sample,
    file.path(TAB_DIR, "04_FINAL_IL10_response_by_sample_subtype_v4.15.0.csv"),
    row.names = FALSE
)

# 09 IL10 response
p_il10 <- ggplot(
    il10_sample,
    aes(
        x = factor(condition, levels = CONDITION_ORDER),
        y = mean_IL10_response
    )
) +
    geom_point(
        aes(
            shape = sample,
            color = factor(macrophage_class, levels = CLASS_ORDER)
        ),
        size = 2.7
    ) +
    facet_wrap(
        ~ macrophage_class,
        scales = "free_y",
        ncol = 3,
        labeller = as_labeller(CLASS_LABELS)
    ) +
    scale_color_manual(values = CLASS_COLORS, guide = "none") +
    labs(
        title = "FINAL IL10-response program across MΦ subtypes",
        subtitle = paste0(
            "Mean normalized expression of: ",
            paste(IL10_USE, collapse = ", ")
        ),
        x = NULL,
        y = "IL10-response score",
        shape = "Sample"
    ) +
    theme_classic(base_size = 9) +
    theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        strip.text = element_text(face = "bold", size = 8)
    )

save_pdf(
    "09_FINAL_IL10_response_by_subtype_condition_v4.15.0.pdf",
    p_il10,
    11,
    8
)

# 10 fixed-order heatmap
GENE_USE <- intersect(GENE_ORDER, rownames(rna_data))

sample_mean_rows <- list()

for (class_now in CLASS_ORDER) {
    for (sample_now in SAMPLE_ORDER) {

        cells_now <- meta$cell[
            meta$macrophage_class == class_now &
            meta$sample == sample_now
        ]

        if (length(cells_now) == 0L) next

        mat <- rna_data[
            GENE_USE,
            cells_now,
            drop = FALSE
        ]

        sample_mean_rows[[paste(class_now, sample_now, sep = "__")]] <- tibble(
            macrophage_class = class_now,
            sample = sample_now,
            gene = GENE_USE,
            mean_expr = as.numeric(Matrix::rowMeans(mat))
        )
    }
}

sample_mean <- bind_rows(sample_mean_rows)

write.csv(
    sample_mean,
    file.path(TAB_DIR, "05_FINAL_representative_gene_sample_means_v4.15.0.csv"),
    row.names = FALSE
)

sample_mean <- sample_mean %>%
    mutate(column_id = paste(macrophage_class, sample, sep = " | "))

COLUMN_ORDER <- as.vector(
    unlist(
        lapply(
            CLASS_ORDER,
            function(class_now) {
                paste(class_now, SAMPLE_ORDER, sep = " | ")
            }
        ),
        use.names = FALSE
    )
)

heat_wide <- sample_mean %>%
    select(gene, column_id, mean_expr) %>%
    pivot_wider(
        names_from = column_id,
        values_from = mean_expr
    )

heat_mat <- as.matrix(
    heat_wide[, -1, drop = FALSE]
)

rownames(heat_mat) <- heat_wide$gene

heat_mat <- heat_mat[
    intersect(GENE_ORDER, rownames(heat_mat)),
    intersect(COLUMN_ORDER, colnames(heat_mat)),
    drop = FALSE
]

heat_z <- t(scale(t(heat_mat)))
heat_z[!is.finite(heat_z)] <- 0

HEAT_LIMIT <- 2

heat_plot <- pmax(
    pmin(heat_z, HEAT_LIMIT),
    -HEAT_LIMIT
)

write.csv(
    data.frame(
        gene = rownames(heat_z),
        heat_z,
        check.names = FALSE
    ),
    file.path(TAB_DIR, "06_FINAL_fixed_order_gene_zscore_matrix_v4.15.0.csv"),
    row.names = FALSE
)

detected_by_block <- lapply(
    GENE_BLOCKS,
    function(genes) intersect(genes, rownames(heat_plot))
)

row_gaps <- cumsum(lengths(detected_by_block))
row_gaps <- row_gaps[row_gaps < nrow(heat_plot)]

n_sample_per_class <- length(SAMPLE_ORDER)

col_gaps <- seq(
    n_sample_per_class,
    ncol(heat_plot) - n_sample_per_class,
    by = n_sample_per_class
)

annotation_col <- tibble(
    column_id = colnames(heat_plot)
) %>%
    separate(
        column_id,
        into = c("macrophage_class", "sample"),
        sep = " \\| ",
        remove = FALSE
    ) %>%
    mutate(
        condition = as.character(canonical_condition(sample))
    ) %>%
    select(macrophage_class, condition) %>%
    as.data.frame()

rownames(annotation_col) <- colnames(heat_plot)

annotation_colors <- list(
    macrophage_class = CLASS_COLORS,
    condition = CONDITION_COLORS
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "10_FINAL_fixed_order_representative_gene_heatmap_v4.15.0.pdf"
    ),
    width = 20,
    height = 9
)

pheatmap::pheatmap(
    heat_plot,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = row_gaps,
    gaps_col = col_gaps,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    color = grDevices::colorRampPalette(
        c("#0033FF", "#FFFFFF", "#FF1A1A")
    )(101),
    breaks = seq(
        -HEAT_LIMIT,
        HEAT_LIMIT,
        length.out = 102
    ),
    border_color = "white",
    fontsize_row = 8,
    fontsize_col = 5.8,
    angle_col = 45,
    main = paste0(
        "FINAL MΦ representative gene programs\n",
        "Clean-B Res2.0 | fixed class/sample order | row z-score"
    )
)

grDevices::dev.off()

figure_index <- tibble(
    figure = sprintf("%02d", 1:10),
    file = c(
        "01_FINAL_Mphi_subtype_UMAP_v4.15.5.pdf",
        "02_FINAL_Mphi_subtype_UMAP_by_condition_v4.15.5.pdf",
        "03_FINAL_Mphi_subtype_composition_by_sample_v4.15.0.pdf",
        "04_FINAL_Mphi_subtype_change_STD_to_CDAHFD_v4.15.0.pdf",
        "05_FINAL_Mphi_subtype_change_Sham_to_Tx_v4.15.0.pdf",
        "06_FINAL_AntiInflammatory_fraction_total_Mphi_v4.15.0.pdf",
        "07_FINAL_M2_fraction_within_M1_M2_v4.15.0.pdf",
        "08_FINAL_M2_to_M1_ratio_v4.15.0.pdf",
        "09_FINAL_IL10_response_by_subtype_condition_v4.15.0.pdf",
        "10_FINAL_fixed_order_representative_gene_heatmap_v4.15.0.pdf"
    ),
    purpose = c(
        "FINAL subtype UMAP",
        "Condition-split FINAL subtype UMAP",
        "Sample-level subtype composition",
        "Disease-associated subtype change",
        "Treatment-associated subtype change",
        "Anti-inflammatory-MΦ / total MΦ",
        "Anti-inflammatory / (Inflammatory + Anti-inflammatory)",
        "Anti-inflammatory / Inflammatory ratio",
        "IL10-response program",
        "FINAL fixed-order representative-gene heatmap"
    )
)

write.csv(
    figure_index,
    file.path(OUTPUT_DIR, "FINAL_FIGURE_INDEX_v4.15.5.csv"),
    row.names = FALSE
)

readme <- c(
    "Mouse MASH MΦ FINAL figures v4.15.5",
    "",
    paste0("Input FINAL RDS: ", INPUT_RDS),
    "",
    paste0("FINAL annotation column: ", FINAL_CLASS_COL),
    "",
    paste0("UMAP reduction used: ", REDUCTION_USE),
    "",
    "FINAL parent clustering:",
    "  Res2.0",
    "",
    "FINAL classification:",
    "  Inflammatory-MΦ",
    "  Anti-inflammatory-MΦ",
    "  ECM-associated inflammatory-MΦ",
    "  Repair/Resolution-MΦ",
    "  Lipid-associated/TREM2-MΦ",
    "  Other",
    "",
    "Definitions:",
    "  M1 = Inflammatory-MΦ",
    "  M2 = Anti-inflammatory-MΦ",
    "  M2 fraction within M1+M2 = M2 / (M1 + M2)",
    "  M2/M1 ratio = (M2 + 0.5) / (M1 + 0.5)",
    "",
    "Statistical interpretation:",
    "  STD vs CDAHFD: descriptive n=1 vs n=1",
    "  Sham vs Tx: biological sample-level n=2 vs n=2",
    "",
    "UMAP display policy:",
    "  Figure 01/02 use a fixed publication window:",
    "    UMAP 1 = -8.8 to -2.0",
    "    UMAP 2 = -7.1 to -0.5",
    "  This excludes the sparse upper-right Other tail from display only.",
    "  No cells are removed from the Seurat object or quantitative analyses.",
    "  Other is drawn first as background; major five MΦ classes are overlaid.",
    "  Major-class dots are enlarged to 1.20 for denser publication rendering.",
    "  Legend point symbols are enlarged independently to 4.4 for readability.",
    "  Figure 01 and Figure 02 use identical x/y display limits.",
    "",
    "Important:",
    "  No reclustering and no reannotation are performed here.",
    "  This script only regenerates FINAL figures from v4.14.5."
)

writeLines(
    readme,
    file.path(OUTPUT_DIR, "README_FINAL_FIGURES_v4.15.5.txt")
)

capture.output(
    sessionInfo(),
    file = file.path(LOG_DIR, "sessionInfo_v4.15.5.txt")
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)
msg("FINAL figure index: ", file.path(OUTPUT_DIR, "FINAL_FIGURE_INDEX_v4.15.5.csv"))

print(figure_index)
print(final_count_table)

print(
    metric_counts %>%
        select(
            sample,
            condition,
            M1_n,
            M2_n,
            AntiInflammatory_percent_total_Mphi,
            M2_percent_within_M1_M2,
            M2_to_M1_ratio
        )
)
