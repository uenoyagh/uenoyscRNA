#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Manuscript visual refinement
# v4.16.5
#
# PURPOSE
#   1) Improve presentation of compositional CLR results for publication.
#   2) Improve FINAL subtype gene-level DotPlot readability and density.
#
# CLR OUTPUTS
#   01 Mean ΔCLR effect-size plot (Tx - Sham)
#   02 Individual samples + group mean + observed min-max range
#   03 Manuscript 2-panel: ΔCLR effect-size + Aitchison PCA
#
# DOTPLOT OUTPUT
#   04 FINAL subtype marker DotPlot, denser layout and larger subtype labels
#
# IMPORTANT
#   - Final annotation remains v4.14.5 / Res2.0.
#   - n=2 Sham and n=2 Tx.
#   - No SEM / 95% CI is shown because n=2 does not support a stable estimate
#     of sampling uncertainty.
#   - Range bars represent observed biological min-max only.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4164)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

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
    "FINAL_Manuscript_Visual_Refinement_v4.16.5"
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

ASSAY_USE <- "RNA"
FINAL_CLASS_COL <- "macrophage_class_Res2_FINAL_v4145"

# ------------------------------------------------------------------------------
# 1. Orders / labels / colors
# ------------------------------------------------------------------------------

CLASS_ORDER <- c(
    "Inflammatory-Mphi",
    "Anti-inflammatory-Mphi",
    "ECM-associated inflammatory-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi",
    "Other"
)

CLASS_ORDER_PLOT <- c(
    "Anti-inflammatory-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi",
    "Other",
    "Inflammatory-Mphi",
    "ECM-associated inflammatory-Mphi"
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

CONDITION_COLORS <- c(
    "Sham" = "#FF3B30",
    "Tx" = "#3366FF"
)

SAMPLE_SHAPES <- c(
    "Sham1" = 16,
    "Sham20" = 17,
    "Tx17" = 15,
    "Tx5" = 3
)

SAMPLE_ORDER_ALL <- c(
    "STD_rep1",
    "CDHFD_rep1",
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
)

# ------------------------------------------------------------------------------
# 2. Marker genes for improved FINAL DotPlot
# ------------------------------------------------------------------------------

MARKER_BLOCKS <- list(

    Inflammatory = c(
        "Il1b",
        "Tnf",
        "Cxcl10",
        "Ccl2",
        "Stat1"
    ),

    Anti_inflammatory = c(
        "Cd163",
        "Mrc1",
        "Il1rn",
        "Mertk",
        "Igf1",
        "Hmox1"
    ),

    ECM_associated = c(
        "Thbs1",
        "Fn1",
        "Tgfb1",
        "Col1a1",
        "Col1a2",
        "Col3a1"
    ),

    Repair_Resolution = c(
        "Mfge8",
        "Gas6",
        "Mmp13",
        "Mmp14",
        "Plau"
    ),

    Lipid_TREM2 = c(
        "Trem2",
        "Gpnmb",
        "Cd9",
        "Lpl",
        "Apoe"
    )
)

MARKER_ORDER <- unique(
    unlist(
        MARKER_BLOCKS,
        use.names = FALSE
    )
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

if (length(
    missing_packages
) > 0L) {

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
# 4. Helpers
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


first_existing <- function(
    x,
    candidates
) {

    hit <- candidates[
        candidates %in%
            x
    ]

    if (length(
        hit
    ) == 0L) {

        return(
            NA_character_
        )
    }

    hit[[1]]
}


get_layer_safe <- function(
    object,
    assay,
    layer
) {

    x <- tryCatch(
        LayerData(
            object,
            assay = assay,
            layer = layer
        ),
        error = function(e) NULL
    )

    if (is.null(
        x
    )) {

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


canonical_condition <- function(
    sample_name
) {

    x <- as.character(
        sample_name
    )

    out <- rep(
        NA_character_,
        length(
            x
        )
    )

    out[
        grepl(
            "^STD",
            x,
            ignore.case = TRUE
        )
    ] <- "STD"

    out[
        grepl(
            "CDHFD|CDAHFD",
            x,
            ignore.case = TRUE
        )
    ] <- "CDAHFD"

    out[
        grepl(
            "^Sham",
            x,
            ignore.case = TRUE
        )
    ] <- "Sham"

    out[
        grepl(
            "^Tx",
            x,
            ignore.case = TRUE
        )
    ] <- "Tx"

    out
}


save_pdf <- function(
    filename,
    plot,
    width,
    height
) {

    ggsave(
        filename = file.path(
            FIG_DIR,
            filename
        ),
        plot = plot,
        device = cairo_pdf,
        width = width,
        height = height,
        units = "in",
        limitsize = FALSE
    )
}

# ------------------------------------------------------------------------------
# 5. Load FINAL object
# ------------------------------------------------------------------------------

if (!file.exists(
    INPUT_RDS
)) {

    stop(
        "FINAL RDS not found:\n",
        INPUT_RDS
    )
}

msg(
    "Loading FINAL MΦ RDS: ",
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
        "Input is not a Seurat object."
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

mphi <- JoinLayers(
    mphi,
    assay = ASSAY_USE
)

if (!"data" %in%
    Layers(
        mphi[[ASSAY_USE]]
    )) {

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

if (is.null(
    rna_data
)) {

    stop(
        "RNA normalized data layer unavailable."
    )
}

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

meta <- mphi@meta.data %>%
    rownames_to_column(
        "cell"
    ) %>%
    transmute(
        cell = cell,

        sample =
            as.character(
                .data[[SAMPLE_COL]]
            ),

        condition =
            canonical_condition(
                .data[[SAMPLE_COL]]
            ),

        macrophage_class =
            as.character(
                .data[[FINAL_CLASS_COL]]
            )
    )

# ==============================================================================
# PART A. Compositional CLR refinement
# ==============================================================================

msg(
    "PART A: compositional CLR refinement"
)

# ------------------------------------------------------------------------------
# 6. Sample x subtype counts
# ------------------------------------------------------------------------------

comp_counts <- meta %>%
    count(
        sample,
        macrophage_class,
        name = "n_cells"
    ) %>%
    complete(
        sample = SAMPLE_ORDER_ALL,
        macrophage_class = CLASS_ORDER,
        fill = list(
            n_cells = 0L
        )
    ) %>%
    mutate(
        condition =
            canonical_condition(
                sample
            )
    )

write.csv(
    comp_counts,
    file.path(
        TAB_DIR,
        "01_sample_subtype_counts_v4.16.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. CLR transform
# ------------------------------------------------------------------------------

PSEUDOCOUNT <- 0.5

clr_df <- comp_counts %>%
    group_by(
        sample
    ) %>%
    mutate(
        log_count =
            log(
                n_cells +
                    PSEUDOCOUNT
            ),

        CLR =
            log_count -
            mean(
                log_count
            )
    ) %>%
    ungroup() %>%
    mutate(
        macrophage_class =
            factor(
                macrophage_class,
                levels =
                    CLASS_ORDER
            ),

        sample =
            factor(
                sample,
                levels =
                    SAMPLE_ORDER_ALL
            )
    )

write.csv(
    clr_df,
    file.path(
        TAB_DIR,
        "02_sample_subtype_CLR_v4.16.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Sham / Tx summaries
# ------------------------------------------------------------------------------

clr_sham_tx <- clr_df %>%
    filter(
        condition %in%
            c(
                "Sham",
                "Tx"
            )
    ) %>%
    mutate(
        condition =
            factor(
                condition,
                levels =
                    c(
                        "Sham",
                        "Tx"
                    )
            )
    )

clr_summary <- clr_sham_tx %>%
    group_by(
        macrophage_class,
        condition
    ) %>%
    summarise(
        n_samples =
            n(),

        mean_CLR =
            mean(
                CLR
            ),

        min_CLR =
            min(
                CLR
            ),

        max_CLR =
            max(
                CLR
            ),

        .groups =
            "drop"
    )

write.csv(
    clr_summary,
    file.path(
        TAB_DIR,
        "03_CLR_group_summary_mean_range_v4.16.5.csv"
    ),
    row.names = FALSE
)

effect_df <- clr_summary %>%
    select(
        macrophage_class,
        condition,
        mean_CLR
    ) %>%
    pivot_wider(
        names_from =
            condition,

        values_from =
            mean_CLR
    ) %>%
    mutate(
        delta_CLR =
            Tx -
            Sham,

        direction =
            case_when(
                delta_CLR >
                    0 ~
                    "Higher in Tx",

                delta_CLR <
                    0 ~
                    "Lower in Tx",

                TRUE ~
                    "No change"
            )
    ) %>%
    mutate(
        macrophage_class =
            factor(
                macrophage_class,
                levels =
                    rev(
                        CLASS_ORDER_PLOT
                    )
            )
    )

write.csv(
    effect_df,
    file.path(
        TAB_DIR,
        "04_CLR_delta_effect_Tx_minus_Sham_v4.16.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. FIGURE 01
# Mean ΔCLR horizontal effect-size plot
# ------------------------------------------------------------------------------

p_delta <- ggplot(
    effect_df,
    aes(
        x = delta_CLR,
        y = macrophage_class
    )
) +
    geom_vline(
        xintercept = 0,
        linetype = 2,
        linewidth = 0.45,
        color = "grey40"
    ) +
    geom_segment(
        aes(
            x = 0,
            xend =
                delta_CLR,
            yend =
                macrophage_class
        ),
        linewidth = 0.8,
        color = "grey40"
    ) +
    geom_point(
        aes(
            fill =
                direction
        ),
        shape = 21,
        size = 4.5,
        stroke = 0.7,
        color = "black"
    ) +
    scale_fill_manual(
        values = c(
            "Higher in Tx" =
                "#3366FF",

            "Lower in Tx" =
                "#FF3B30",

            "No change" =
                "grey70"
        )
    ) +
    scale_y_discrete(
        labels =
            function(x) {
                CLASS_LABELS[
                    x
                ]
            }
    ) +
    labs(
        title =
            "MΦ compositional effect: Sham → Tx",
        subtitle =
            "Mean ΔCLR = mean(Tx) − mean(Sham); positive values indicate relative enrichment in Tx",
        x =
            "Mean ΔCLR (Tx − Sham)",
        y =
            NULL,
        fill =
            NULL
    ) +
    theme_classic(
        base_size = 11
    ) +
    theme(
        plot.title =
            element_text(
                face =
                    "bold",
                size =
                    13
            ),

        axis.text.y =
            element_text(
                face =
                    "bold",
                size =
                    10.5
            ),

        axis.text.x =
            element_text(
                size =
                    9.5
            ),

        legend.position =
            "bottom"
    )

save_pdf(
    "01_FINAL_CLR_delta_effect_size_v4.16.5.pdf",
    p_delta,
    8.5,
    6
)

# ------------------------------------------------------------------------------
# 10. FIGURE 02
# Individual points + group mean + observed min-max range
# ------------------------------------------------------------------------------

clr_summary_plot <- clr_summary %>%
    mutate(
        macrophage_class =
            factor(
                macrophage_class,
                levels =
                    rev(
                        CLASS_ORDER_PLOT
                    )
            )
    )

clr_individual_plot <- clr_sham_tx %>%
    mutate(
        macrophage_class =
            factor(
                macrophage_class,
                levels =
                    rev(
                        CLASS_ORDER_PLOT
                    )
            )
    )

position_map <- c(
    "Sham" = -0.13,
    "Tx" = 0.13
)

base_y <- as.numeric(
    clr_summary_plot$macrophage_class
)

clr_summary_plot$y_numeric <-
    base_y +
    position_map[
        as.character(
            clr_summary_plot$condition
        )
    ]

clr_individual_plot$y_numeric <-
    as.numeric(
        clr_individual_plot$macrophage_class
    ) +
    position_map[
        as.character(
            clr_individual_plot$condition
        )
    ]

p_range <- ggplot() +
    geom_vline(
        xintercept = 0,
        linetype = 3,
        linewidth = 0.35,
        color = "grey70"
    ) +
    geom_errorbarh(
        data =
            clr_summary_plot,
        aes(
            y =
                y_numeric,
            xmin =
                min_CLR,
            xmax =
                max_CLR,
            color =
                condition
        ),
        height = 0,
        linewidth = 0.9
    ) +
    geom_point(
        data =
            clr_summary_plot,
        aes(
            x =
                mean_CLR,
            y =
                y_numeric,
            fill =
                condition
        ),
        shape = 21,
        size = 4.2,
        stroke = 0.7,
        color = "black"
    ) +
    geom_point(
        data =
            clr_individual_plot,
        aes(
            x =
                CLR,
            y =
                y_numeric,
            shape =
                sample,
            color =
                condition
        ),
        size = 2.7,
        stroke = 0.75
    ) +
    scale_color_manual(
        values =
            CONDITION_COLORS
    ) +
    scale_fill_manual(
        values =
            CONDITION_COLORS
    ) +
    scale_shape_manual(
        values =
            SAMPLE_SHAPES
    ) +
    scale_y_continuous(
        breaks =
            seq_along(
                rev(
                    CLASS_ORDER_PLOT
                )
            ),
        labels =
            CLASS_LABELS[
                rev(
                    CLASS_ORDER_PLOT
                )
            ],
        expand =
            expansion(
                add =
                    0.6
            )
    ) +
    labs(
        title =
            "FINAL MΦ compositional remodeling: Sham vs Tx",
        subtitle =
            "Large circles = group mean | horizontal bars = observed min–max | symbols = biological samples",
        x =
            "CLR abundance",
        y =
            NULL,
        color =
            "Condition",
        fill =
            "Condition",
        shape =
            "Sample"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title =
            element_text(
                face =
                    "bold"
            ),

        axis.text.y =
            element_text(
                face =
                    "bold",
                size =
                    9.5
            ),

        legend.position =
            "right"
    )

save_pdf(
    "02_FINAL_CLR_individual_mean_observed_range_v4.16.5.pdf",
    p_range,
    10,
    6.5
)

# ------------------------------------------------------------------------------
# 11. Aitchison PCA
# ------------------------------------------------------------------------------

clr_mat <- clr_df %>%
    select(
        sample,
        macrophage_class,
        CLR
    ) %>%
    pivot_wider(
        names_from =
            macrophage_class,
        values_from =
            CLR
    ) %>%
    arrange(
        factor(
            sample,
            levels =
                SAMPLE_ORDER_ALL
        )
    )

clr_numeric <- as.matrix(
    clr_mat[
        ,
        -1,
        drop = FALSE
    ]
)

rownames(
    clr_numeric
) <- clr_mat$sample

pca <- prcomp(
    clr_numeric,
    center = FALSE,
    scale. = FALSE
)

pca_df <- as.data.frame(
    pca$x[
        ,
        1:2,
        drop = FALSE
    ]
) %>%
    rownames_to_column(
        "sample"
    ) %>%
    mutate(
        condition =
            canonical_condition(
                sample
            )
    )

var_exp <- (
    pca$sdev^2 /
        sum(
            pca$sdev^2
        )
) * 100

PCA_COLORS <- c(
    "STD" = "#3366FF",
    "CDAHFD" = "#FF3B30",
    "Sham" = "#FF3B30",
    "Tx" = "#3366FF"
)

PCA_SHAPES <- c(
    "STD_rep1" = 16,
    "CDHFD_rep1" = 17,
    "Sham1" = 16,
    "Sham20" = 17,
    "Tx17" = 15,
    "Tx5" = 3
)

p_pca <- ggplot(
    pca_df,
    aes(
        x = PC1,
        y = PC2,
        color = condition,
        shape = sample
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = 3,
        linewidth = 0.35,
        color = "grey70"
    ) +
    geom_vline(
        xintercept = 0,
        linetype = 3,
        linewidth = 0.35,
        color = "grey70"
    ) +
    geom_point(
        size = 3.6,
        stroke = 0.8
    ) +
    geom_text(
        aes(
            label = sample
        ),
        nudge_y = 0.10,
        size = 3.0,
        show.legend = FALSE
    ) +
    scale_color_manual(
        values =
            PCA_COLORS
    ) +
    scale_shape_manual(
        values =
            PCA_SHAPES
    ) +
    labs(
        title =
            "Aitchison PCA of FINAL MΦ composition",
        subtitle =
            "Each point = biological sample",
        x =
            paste0(
                "PC1 (",
                round(
                    var_exp[1],
                    1
                ),
                "%)"
            ),
        y =
            paste0(
                "PC2 (",
                round(
                    var_exp[2],
                    1
                ),
                "%)"
            ),
        color =
            "Condition",
        shape =
            "Sample"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title =
            element_text(
                face =
                    "bold"
            )
    )

# ------------------------------------------------------------------------------
# 12. FIGURE 03
# Manuscript 2-panel: ΔCLR effect + Aitchison PCA
# ------------------------------------------------------------------------------

p_manuscript <- (
    p_delta +
        theme(
            legend.position =
                "bottom"
        )
) +
    (
        p_pca +
            theme(
                legend.position =
                    "right"
            )
    ) +
    patchwork::plot_layout(
        widths =
            c(
                1,
                1.05
            )
    ) +
    patchwork::plot_annotation(
        title =
            "Macrophage compositional remodeling after transplantation",
        theme =
            theme(
                plot.title =
                    element_text(
                        face =
                            "bold",
                        size =
                            14
                    )
            )
    )

save_pdf(
    "03_FINAL_CLR_effect_plus_Aitchison_PCA_manuscript_v4.16.5.pdf",
    p_manuscript,
    14,
    6.5
)

# ==============================================================================
# PART B. Improved FINAL subtype DotPlot
# ==============================================================================

msg(
    "PART B: improved FINAL subtype DotPlot"
)

genes_use <- intersect(
    MARKER_ORDER,
    rownames(
        rna_data
    )
)

dot_rows <- list()

for (class_now in CLASS_ORDER) {

    cells_now <- meta$cell[
        meta$macrophage_class ==
            class_now
    ]

    if (length(
        cells_now
    ) == 0L) {

        next
    }

    mat <- rna_data[
        genes_use,
        cells_now,
        drop = FALSE
    ]

    dot_rows[[class_now]] <- tibble(
        macrophage_class =
            class_now,

        gene =
            genes_use,

        avg_expr =
            as.numeric(
                Matrix::rowMeans(
                    mat
                )
            ),

        pct_expr =
            100 *
            as.numeric(
                Matrix::rowMeans(
                    mat >
                        0
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
        avg_expr_z =
            as.numeric(
                scale(
                    avg_expr
                )
            )
    ) %>%
    ungroup() %>%
    mutate(
        avg_expr_z =
            ifelse(
                is.finite(
                    avg_expr_z
                ),
                avg_expr_z,
                0
            ),

        macrophage_class =
            factor(
                macrophage_class,
                levels =
                    rev(
                        CLASS_ORDER
                    )
            ),

        gene =
            factor(
                gene,
                levels =
                    MARKER_ORDER
            )
    )

write.csv(
    dot_df,
    file.path(
        TAB_DIR,
        "05_FINAL_subtype_DotPlot_numeric_v4.16.5.csv"
    ),
    row.names = FALSE
)

# gene-block separators
block_lengths <- lengths(
    lapply(
        MARKER_BLOCKS,
        function(x) {
            intersect(
                x,
                genes_use
            )
        }
    )
)

block_breaks <- cumsum(
    block_lengths
)

block_breaks <- block_breaks[
    block_breaks <
        length(
            genes_use
        )
]

p_dot <- ggplot(
    dot_df,
    aes(
        x =
            gene,
        y =
            macrophage_class
    )
) +
    geom_point(
        aes(
            size =
                pct_expr,
            color =
                avg_expr_z
        ),
        alpha = 0.95
    ) +
    geom_vline(
        xintercept =
            block_breaks +
            0.5,
        linewidth =
            0.35,
        color =
            "grey85"
    ) +
    scale_size_continuous(
        range =
            c(
                1.0,
                10.0
            ),
        limits =
            c(
                0,
                100
            ),
        breaks =
            c(
                25,
                50,
                75,
                100
            ),
        name =
            "% expressed"
    ) +
    scale_color_gradient2(
        low =
            "#0033FF",
        mid =
            "#FFFFFF",
        high =
            "#FF1A1A",
        midpoint =
            0,
        limits =
            c(
                -2.5,
                2.5
            ),
        oob =
            scales::squish,
        name =
            "Average\nexpression\nz-score"
    ) +
    scale_y_discrete(
        labels =
            function(x) {
                CLASS_LABELS[
                    x
                ]
            },
        expand =
            expansion(
                add =
                    c(
                        0.55,
                        0.55
                    )
            )
    ) +
    scale_x_discrete(
        expand =
            expansion(
                add =
                    c(
                        0.70,
                        0.70
                    )
            )
    ) +
    labs(
        title =
            "FINAL MΦ subtype gene-level validation",
        subtitle =
            "Clean-B Res2.0 final annotation",
        x =
            NULL,
        y =
            NULL
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

        plot.subtitle =
            element_text(
                size =
                    10.5
            ),

        axis.text.y =
            element_text(
                face =
                    "plain",
                size =
                    12.5,
                margin =
                    margin(
                        r =
                            5
                    )
            ),

        axis.text.x =
            element_text(
                angle =
                    60,
                hjust =
                    1,
                vjust =
                    1,
                size =
                    9.5,
                face =
                    "plain"
            ),

        axis.ticks =
            element_blank(),

        legend.title =
            element_text(
                size =
                    10
            ),

        legend.text =
            element_text(
                size =
                    9.5
            ),

        panel.grid =
            element_blank(),

        plot.margin =
            margin(
                10,
                12,
                10,
                12
            )
    )

save_pdf(
    "04_FINAL_subtype_marker_DotPlot_dense_unclipped_v4.16.5.pdf",
    p_dot,
    13.5,
    5.8
)

# ------------------------------------------------------------------------------
# 13. README
# ------------------------------------------------------------------------------

readme <- c(
    "Mouse MASH MΦ manuscript visual refinement v4.16.5",
    "",
    "CLR figures:",
    "  01 mean ΔCLR effect-size plot",
    "  02 individual samples + group mean + observed min-max",
    "  03 ΔCLR effect-size + Aitchison PCA manuscript panel",
    "",
    "DotPlot:",
    "  04 denser FINAL subtype marker DotPlot",
    "  subtype label size retained at 12.5, regular weight",
    "  gene label size retained at 9.5, regular weight",
    "  dots remain enlarged",
    "  edge padding increased to prevent clipping of large dots",
    "  gene-block separators retained",
    "",
    "Important:",
    "  n=2 Sham and n=2 Tx.",
    "  No SEM or 95% CI is plotted.",
    "  Range bars represent observed biological min-max only."
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_visual_refinement_v4.16.5.txt"
    )
)

capture.output(
    sessionInfo(),
    file =
        file.path(
            LOG_DIR,
            "sessionInfo_v4.16.5.txt"
        )
)

msg(
    "DONE."
)

msg(
    "Output: ",
    OUTPUT_DIR
)
