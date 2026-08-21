#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Manuscript-gap completion analysis
# FINAL Clean-B / Res2.0 / v4.16.2
#
# PURPOSE
#   Add the few analyses that remain most useful for a manuscript-level
#   description of macrophage-compartment remodeling after transplantation.
#
# ADDITIONAL OUTPUTS
#   1) Formal sample-level compositional effect (CLR / log-ratio)
#   2) FINAL subtype gene-level DotPlot
#   3) FINAL subtype x biological-sample marker heatmap
#   4) Key mechanistic gene sample-level pseudobulk consistency: Sham -> Tx
#
# IMPORTANT
#   - Uses FINAL v4.14.5 annotation.
#   - Parent clustering remains Res2.0.
#   - No cell-level p-value is used as biological-replicate evidence.
#   - Sham vs Tx has n=2 vs n=2; results are effect-size / consistency focused.
#   - STD vs CDAHFD has n=1 vs n=1 and is descriptive only.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4162)

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
    "FINAL_Manuscript_Gap_Analysis_v4.16.2"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"
FINAL_CLASS_COL <- "macrophage_class_Res2_FINAL_v4145"

# ------------------------------------------------------------------------------
# 1. Class / sample definitions
# ------------------------------------------------------------------------------

CLASS_ORDER <- c(
    "Inflammatory-Mphi",
    "Anti-inflammatory-Mphi",
    "ECM-associated inflammatory-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi",
    "Other"
)

BIOLOGICAL_CLASS_ORDER <- CLASS_ORDER[
    CLASS_ORDER != "Other"
]

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

SAMPLE_ORDER <- c(
    "STD_rep1",
    "CDHFD_rep1",
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
)

CONDITION_ORDER <- c(
    "STD",
    "CDAHFD",
    "Sham",
    "Tx"
)

CONDITION_COLORS <- c(
    "STD" = "#3366FF",
    "CDAHFD" = "#FF3B30",
    "Sham" = "#FF3B30",
    "Tx" = "#3366FF"
)

# ------------------------------------------------------------------------------
# 2. Gene panels
# ------------------------------------------------------------------------------

FINAL_MARKER_BLOCKS <- list(

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

    ECM_associated_inflammatory = c(
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

FINAL_MARKER_ORDER <- unique(
    unlist(
        FINAL_MARKER_BLOCKS,
        use.names = FALSE
    )
)

MECHANISTIC_GENES <- unique(c(
    "Cd163",
    "Mrc1",
    "Il1rn",
    "Mertk",
    "Igf1",
    "Hmox1",
    "Il10ra",
    "Stat3",
    "Socs3",
    "Thbs1",
    "Fn1",
    "Tgfb1",
    "Il1b",
    "Stat1",
    "Cxcl10",
    "Trem2",
    "Gpnmb",
    "Cd9",
    "Lpl",
    "Apoe",
    "Mfge8",
    "Gas6",
    "Mmp13",
    "Mmp14",
    "Plau"
))

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

    hit <- candidates[
        candidates %in% x
    ]

    if (length(hit) == 0L) {
        return(NA_character_)
    }

    hit[[1]]
}

canonical_condition <- function(x) {

    x <- as.character(x)

    out <- rep(
        NA_character_,
        length(x)
    )

    out[
        grepl("^STD", x, ignore.case = TRUE)
    ] <- "STD"

    out[
        grepl("CDHFD|CDAHFD", x, ignore.case = TRUE)
    ] <- "CDAHFD"

    out[
        grepl("^Sham", x, ignore.case = TRUE)
    ] <- "Sham"

    out[
        grepl("^Tx", x, ignore.case = TRUE)
    ] <- "Tx"

    factor(
        out,
        levels = CONDITION_ORDER
    )
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

row_zscore <- function(mat) {

    z <- t(
        scale(
            t(mat)
        )
    )

    z[
        !is.finite(z)
    ] <- 0

    z
}

# ------------------------------------------------------------------------------
# 5. Load FINAL object
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
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

counts <- get_layer_safe(
    mphi,
    ASSAY_USE,
    "counts"
)

if (is.null(rna_data)) {
    stop("RNA data layer unavailable.")
}

if (is.null(counts)) {
    stop("RNA counts layer unavailable.")
}

# ------------------------------------------------------------------------------
# 6. Metadata
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

if (is.na(SAMPLE_COL)) {
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
        sample = as.character(
            .data[[SAMPLE_COL]]
        ),
        condition = as.character(
            canonical_condition(
                .data[[SAMPLE_COL]]
            )
        ),
        macrophage_class = as.character(
            .data[[FINAL_CLASS_COL]]
        )
    )

# ==============================================================================
# PART 1
# Formal sample-level compositional effect
# ==============================================================================

msg(
    "PART 1: formal compositional effect"
)

# ------------------------------------------------------------------------------
# 7. Sample x subtype count matrix
# ------------------------------------------------------------------------------

comp_counts <- meta %>%
    count(
        sample,
        macrophage_class,
        name = "n_cells"
    ) %>%
    complete(
        sample = SAMPLE_ORDER,
        macrophage_class = CLASS_ORDER,
        fill = list(
            n_cells = 0L
        )
    ) %>%
    mutate(
        condition = as.character(
            canonical_condition(
                sample
            )
        )
    )

write.csv(
    comp_counts,
    file.path(
        TAB_DIR,
        "01_sample_subtype_counts_v4.16.2.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. CLR transform
#
# pseudo-count = 0.5 cells.
# We retain Other in the compositional denominator because it is a real
# retained Clean-B compartment. Biological interpretation focuses on 5 classes.
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

        geometric_log_mean =
            mean(
                log_count
            ),

        CLR =
            log_count -
            geometric_log_mean
    ) %>%
    ungroup() %>%
    mutate(
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
    clr_df,
    file.path(
        TAB_DIR,
        "02_sample_subtype_CLR_v4.16.2.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Sham -> Tx CLR effect
# ------------------------------------------------------------------------------

clr_effect <- clr_df %>%
    filter(
        condition %in%
            c(
                "Sham",
                "Tx"
            )
    ) %>%
    group_by(
        macrophage_class,
        condition
    ) %>%
    summarise(
        n_samples = sum(
            is.finite(
                CLR
            )
        ),
        mean_CLR = mean(
            CLR,
            na.rm = TRUE
        ),
        min_CLR = min(
            CLR,
            na.rm = TRUE
        ),
        max_CLR = max(
            CLR,
            na.rm = TRUE
        ),
        .groups = "drop"
    ) %>%
    pivot_wider(
        names_from = condition,
        values_from = c(
            n_samples,
            mean_CLR,
            min_CLR,
            max_CLR
        ),
        names_sep = "_"
    ) %>%
    mutate(
        delta_CLR_Tx_minus_Sham =
            mean_CLR_Tx -
            mean_CLR_Sham
    ) %>%
    arrange(
        desc(
            delta_CLR_Tx_minus_Sham
        )
    )

write.csv(
    clr_effect,
    file.path(
        TAB_DIR,
        "03_CLR_effect_Tx_minus_Sham_v4.16.2.csv"
    ),
    row.names = FALSE
)

p_clr <- ggplot(
    clr_df %>%
        filter(
            condition %in%
                c(
                    "Sham",
                    "Tx"
                )
        ),
    aes(
        x = condition,
        y = CLR
    )
) +
    geom_line(
        aes(
            group = sample
        ),
        linewidth = 0.4,
        alpha = 0.35
    ) +
    geom_point(
        aes(
            shape = sample,
            color = condition
        ),
        size = 2.8
    ) +
    facet_wrap(
        ~ macrophage_class,
        ncol = 3,
        scales = "free_y",
        labeller = as_labeller(
            CLASS_LABELS
        )
    ) +
    scale_color_manual(
        values = CONDITION_COLORS[
            c(
                "Sham",
                "Tx"
            )
        ]
    ) +
    labs(
        title =
            "FINAL MΦ compositional remodeling: Sham vs Tx",
        subtitle =
            "Centered log-ratio (CLR); each point = biological sample",
        x = NULL,
        y = "CLR abundance",
        shape = "Sample",
        color = "Condition"
    ) +
    theme_classic(
        base_size = 9
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        strip.text = element_text(
            face = "bold",
            size = 8
        )
    )

save_pdf(
    "01_FINAL_compositional_CLR_Sham_vs_Tx_v4.16.2.pdf",
    p_clr,
    11,
    8
)

# ------------------------------------------------------------------------------
# 10. Aitchison PCA
# ------------------------------------------------------------------------------

clr_mat <- clr_df %>%
    select(
        sample,
        macrophage_class,
        CLR
    ) %>%
    pivot_wider(
        names_from = macrophage_class,
        values_from = CLR
    ) %>%
    arrange(
        factor(
            sample,
            levels = SAMPLE_ORDER
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
        condition = canonical_condition(
            sample
        )
    )

write.csv(
    pca_df,
    file.path(
        TAB_DIR,
        "04_Aitchison_CLR_PCA_coordinates_v4.16.2.csv"
    ),
    row.names = FALSE
)

var_exp <- (
    pca$sdev^2 /
        sum(
            pca$sdev^2
        )
) * 100

p_pca <- ggplot(
    pca_df,
    aes(
        x = PC1,
        y = PC2,
        color = condition,
        label = sample
    )
) +
    geom_point(
        size = 3.2
    ) +
    geom_text(
        nudge_y = 0.08,
        size = 3
    ) +
    scale_color_manual(
        values = CONDITION_COLORS,
        drop = FALSE
    ) +
    labs(
        title =
            "Aitchison geometry of FINAL MΦ composition",
        subtitle =
            "PCA of sample-level CLR-transformed subtype abundances",
        x = paste0(
            "PC1 (",
            round(
                var_exp[1],
                1
            ),
            "%)"
        ),
        y = paste0(
            "PC2 (",
            round(
                var_exp[2],
                1
            ),
            "%)"
        ),
        color = "Condition"
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
    "02_FINAL_Aitchison_CLR_PCA_v4.16.2.pdf",
    p_pca,
    7.5,
    6
)

# ==============================================================================
# PART 2
# FINAL subtype gene-level validation DotPlot
# ==============================================================================

msg(
    "PART 2: FINAL subtype gene-level DotPlot"
)

genes_use <- intersect(
    FINAL_MARKER_ORDER,
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
        avg_expr_z =
            as.numeric(
                scale(
                    avg_expr
                )
            )
    ) %>%
    ungroup() %>%
    mutate(
        avg_expr_z = ifelse(
            is.finite(
                avg_expr_z
            ),
            avg_expr_z,
            0
        ),

        macrophage_class = factor(
            macrophage_class,
            levels = rev(
                CLASS_ORDER
            )
        ),

        gene = factor(
            gene,
            levels = FINAL_MARKER_ORDER
        )
    )

write.csv(
    dot_df,
    file.path(
        TAB_DIR,
        "05_FINAL_subtype_marker_DotPlot_numeric_v4.16.2.csv"
    ),
    row.names = FALSE
)

p_dot <- ggplot(
    dot_df,
    aes(
        x = gene,
        y = macrophage_class
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
            0.3,
            7.5
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
        limits = c(
            -2.5,
            2.5
        ),
        oob = scales::squish,
        name = "Average\nexpression\nz-score"
    ) +
    scale_y_discrete(
        labels = function(x) {
            CLASS_LABELS[
                x
            ]
        }
    ) +
    labs(
        title =
            "FINAL MΦ subtype gene-level validation",
        subtitle =
            "Clean-B Res2.0 final annotation",
        x = NULL,
        y = NULL
    ) +
    theme_classic(
        base_size = 9
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
    "03_FINAL_subtype_marker_gene_DotPlot_v4.16.2.pdf",
    p_dot,
    15,
    5.5
)

# ==============================================================================
# PART 3
# FINAL subtype x biological-sample marker heatmap
# ==============================================================================

msg(
    "PART 3: sample-level marker heatmap"
)

sample_mean_rows <- list()

for (class_now in BIOLOGICAL_CLASS_ORDER) {

    for (sample_now in SAMPLE_ORDER) {

        cells_now <- meta$cell[
            meta$macrophage_class ==
                class_now &
            meta$sample ==
                sample_now
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

        sample_mean_rows[[
            paste(
                class_now,
                sample_now,
                sep = "__"
            )
        ]] <- tibble(
            macrophage_class =
                class_now,

            sample =
                sample_now,

            condition =
                as.character(
                    canonical_condition(
                        sample_now
                    )
                ),

            gene =
                genes_use,

            mean_expr =
                as.numeric(
                    Matrix::rowMeans(
                        mat
                    )
                )
        )
    }
}

sample_mean <- bind_rows(
    sample_mean_rows
)

write.csv(
    sample_mean,
    file.path(
        TAB_DIR,
        "06_FINAL_subtype_sample_marker_means_v4.16.2.csv"
    ),
    row.names = FALSE
)

# Heatmap by biological sample, averaging over final class-specific cells.
# Columns: class | sample
sample_mean <- sample_mean %>%
    mutate(
        column_id = paste(
            macrophage_class,
            sample,
            sep = " | "
        )
    )

column_order <- as.vector(
    unlist(
        lapply(
            BIOLOGICAL_CLASS_ORDER,
            function(class_now) {
                paste(
                    class_now,
                    SAMPLE_ORDER,
                    sep = " | "
                )
            }
        ),
        use.names = FALSE
    )
)

heat_wide <- sample_mean %>%
    select(
        gene,
        column_id,
        mean_expr
    ) %>%
    pivot_wider(
        names_from = column_id,
        values_from = mean_expr
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
) <- heat_wide$gene

heat_mat <- heat_mat[
    intersect(
        FINAL_MARKER_ORDER,
        rownames(
            heat_mat
        )
    ),
    intersect(
        column_order,
        colnames(
            heat_mat
        )
    ),
    drop = FALSE
]

heat_z <- row_zscore(
    heat_mat
)

HEAT_LIMIT <- 2

heat_plot <- pmax(
    pmin(
        heat_z,
        HEAT_LIMIT
    ),
    -HEAT_LIMIT
)

# Gene block gaps
detected_by_block <- lapply(
    FINAL_MARKER_BLOCKS,
    function(x) {
        intersect(
            x,
            rownames(
                heat_plot
            )
        )
    }
)

row_gaps <- cumsum(
    lengths(
        detected_by_block
    )
)

row_gaps <- row_gaps[
    row_gaps <
        nrow(
            heat_plot
        )
]

col_gaps <- seq(
    length(
        SAMPLE_ORDER
    ),
    ncol(
        heat_plot
    ) -
        length(
            SAMPLE_ORDER
        ),
    by =
        length(
            SAMPLE_ORDER
        )
)

annotation_col <- tibble(
    column_id =
        colnames(
            heat_plot
        )
) %>%
    separate(
        column_id,
        into = c(
            "macrophage_class",
            "sample"
        ),
        sep = " \\| ",
        remove = FALSE
    ) %>%
    mutate(
        condition = as.character(
            canonical_condition(
                sample
            )
        )
    ) %>%
    select(
        macrophage_class,
        condition
    ) %>%
    as.data.frame()

rownames(
    annotation_col
) <- colnames(
    heat_plot
)

annotation_colors <- list(
    macrophage_class =
        CLASS_COLORS,

    condition =
        CONDITION_COLORS
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "04_FINAL_subtype_sample_marker_heatmap_v4.16.2.pdf"
    ),
    width = 18,
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
        c(
            "#0033FF",
            "#FFFFFF",
            "#FF1A1A"
        )
    )(101),
    breaks = seq(
        -HEAT_LIMIT,
        HEAT_LIMIT,
        length.out = 102
    ),
    border_color = "white",
    fontsize_row = 7.5,
    fontsize_col = 5.5,
    angle_col = 45,
    main = paste0(
        "FINAL MΦ subtype marker programs by biological sample\n",
        "fixed class/sample order | row z-score"
    )
)

grDevices::dev.off()

# ==============================================================================
# PART 4
# Key mechanistic gene sample-level pseudobulk consistency
# ==============================================================================

msg(
    "PART 4: mechanistic-gene sample-level consistency"
)

mech_use <- intersect(
    MECHANISTIC_GENES,
    rownames(
        counts
    )
)

pb_rows <- list()

for (class_now in BIOLOGICAL_CLASS_ORDER) {

    for (sample_now in c(
        "Sham1",
        "Sham20",
        "Tx17",
        "Tx5"
    )) {

        cells_now <- meta$cell[
            meta$macrophage_class ==
                class_now &
            meta$sample ==
                sample_now
        ]

        if (length(
            cells_now
        ) == 0L) {
            next
        }

        pb <- Matrix::rowSums(
            counts[
                mech_use,
                cells_now,
                drop = FALSE
            ]
        )

        lib <- sum(
            Matrix::rowSums(
                counts[
                    ,
                    cells_now,
                    drop = FALSE
                ]
            )
        )

        pb_rows[[
            paste(
                class_now,
                sample_now,
                sep = "__"
            )
        ]] <- tibble(
            macrophage_class =
                class_now,

            sample =
                sample_now,

            condition =
                as.character(
                    canonical_condition(
                        sample_now
                    )
                ),

            gene =
                mech_use,

            CPM =
                as.numeric(
                    pb
                ) /
                lib *
                1e6,

            logCPM =
                log2(
                    as.numeric(
                        pb
                    ) /
                    lib *
                    1e6 +
                    1
                )
        )
    }
}

pb_mech <- bind_rows(
    pb_rows
)

write.csv(
    pb_mech,
    file.path(
        TAB_DIR,
        "07_mechanistic_gene_pseudobulk_by_sample_v4.16.2.csv"
    ),
    row.names = FALSE
)

effect_mech <- pb_mech %>%
    select(
        macrophage_class,
        gene,
        sample,
        logCPM
    ) %>%
    pivot_wider(
        names_from = sample,
        values_from = logCPM
    ) %>%
    filter(
        !is.na(
            Sham1
        ),
        !is.na(
            Sham20
        ),
        !is.na(
            Tx17
        ),
        !is.na(
            Tx5
        )
    ) %>%
    mutate(
        Sham_mean =
            (
                Sham1 +
                Sham20
            ) /
            2,

        Tx_mean =
            (
                Tx17 +
                Tx5
            ) /
            2,

        Tx17_minus_Sham =
            Tx17 -
            Sham_mean,

        Tx5_minus_Sham =
            Tx5 -
            Sham_mean,

        mean_effect_Tx_minus_Sham =
            Tx_mean -
            Sham_mean,

        replicate_concordant =
            sign(
                Tx17_minus_Sham
            ) ==
            sign(
                Tx5_minus_Sham
            )
    )

write.csv(
    effect_mech,
    file.path(
        TAB_DIR,
        "08_mechanistic_gene_Tx_vs_Sham_effects_v4.16.2.csv"
    ),
    row.names = FALSE
)

p_mech <- effect_mech %>%
    filter(
        replicate_concordant
    ) %>%
    mutate(
        macrophage_class = factor(
            macrophage_class,
            levels =
                BIOLOGICAL_CLASS_ORDER
        )
    ) %>%
    ggplot(
        aes(
            x = mean_effect_Tx_minus_Sham,
            y = gene
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
            xend = mean_effect_Tx_minus_Sham,
            yend = gene
        ),
        linewidth = 0.45
    ) +
    geom_point(
        aes(
            shape =
                mean_effect_Tx_minus_Sham >
                0
        ),
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
        title =
            "Mechanistically relevant MΦ genes: Sham → Tx",
        subtitle =
            "Only genes with concordant direction in Tx17 and Tx5",
        x =
            "Mean log2(CPM+1) difference: Tx − Sham",
        y = NULL,
        shape =
            "Higher in Tx"
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
    "05_FINAL_mechanistic_gene_Tx_vs_Sham_concordance_v4.16.2.pdf",
    p_mech,
    12,
    10
)

# ------------------------------------------------------------------------------
# 11. Manuscript-readiness summary table
# ------------------------------------------------------------------------------

summary_table <- tibble(
    item = c(
        "FINAL subtype UMAP",
        "Sample-level subtype composition",
        "Anti-inflammatory-MΦ fraction",
        "M2/(M1+M2)",
        "M2/M1 ratio",
        "M1/M2 transcriptional programs",
        "Anti-inflammatory-MΦ heterogeneity",
        "IL10-response-high vs low M2",
        "Subtype pseudobulk Tx vs Sham",
        "Tx replicate consistency",
        "Subtype pathway enrichment",
        "Formal sample-level compositional effect",
        "FINAL subtype gene-level validation",
        "Mechanistic gene sample-level consistency",
        "Total MΦ abundance among whole liver"
    ),

    status_after_v4162 = c(
        rep(
            "Available",
            11
        ),
        "Added in v4.16.2",
        "Added in v4.16.2",
        "Added in v4.16.2",
        "Requires Whole-Liver object"
    ),

    manuscript_role = c(
        "Core phenotype",
        "Core phenotype",
        "Core phenotype",
        "Core phenotype",
        "Core phenotype",
        "Functional support",
        "Mechanistic / supplementary",
        "Mechanistic support",
        "Subtype response",
        "Replicate QC",
        "Mechanistic pathway support",
        "Composition-aware support",
        "Annotation validation",
        "Mechanistic consistency",
        "Whole-liver context"
    )
)

write.csv(
    summary_table,
    file.path(
        OUTPUT_DIR,
        "MANUSCRIPT_READINESS_SUMMARY_v4.16.2.csv"
    ),
    row.names = FALSE
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.16.2.txt"
    )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)

print(
    summary_table
)
