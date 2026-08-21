#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Gene-level DotPlot / heatmap
# v4.13.2
#
# PRIMARY DATASET
#   Clean-B
#
# PURPOSE
#   Validate functional classes at gene level and show which genes drive
#   each module/signature.
#
# OUTPUT
#   1) subtype × condition custom DotPlot
#   2) subtype × biological-sample custom DotPlot
#   3) sample-level pseudobulk heatmap
#   4) condition-level heatmap
#   5) compact representative-gene heatmap
#
# IMPORTANT
#   - Seurat DotPlot() is NOT used.
#   - DotPlot values are calculated directly from normalized RNA.
#   - Heatmaps are biological-sample level where applicable.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4132)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_cleaning_sensitivity_v4.11.7",
    "RDS",
    "Mouse_Mphi_Res2_Clean_B_v4.11.7.rds"
)

OUTPUT_DIR <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_gene_level_DotPlot_heatmap_CleanB_v4.13.2"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"

CONDITION_ORDER <- c(
    "STD",
    "CDAHFD",
    "Sham",
    "Tx"
)

SAMPLE_ORDER <- c(
    "STD_rep1",
    "CDHFD_rep1",
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
)

SUBTYPE_ORDER <- c(
    "Inflammatory-Mphi",
    "Anti-inflammatory-Mphi",
    "Fibrogenic-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi"
)

SUBTYPE_LABELS <- c(
    "Inflammatory-Mphi"           = "Inflammatory-MΦ",
    "Anti-inflammatory-Mphi"      = "Anti-inflammatory-MΦ",
    "Fibrogenic-Mphi"             = "Fibrogenic-MΦ",
    "Repair/Resolution-Mphi"      = "Repair/Resolution-MΦ",
    "Lipid-associated/TREM2-Mphi" = "Lipid-associated/TREM2-MΦ"
)

# ------------------------------------------------------------------------------
# 1. Representative gene sets
# ------------------------------------------------------------------------------

GENE_SETS <- list(

    Inflammatory = c(
        "Il1b",
        "Tnf",
        "Ccl2",
        "Cxcl10",
        "Nos2",
        "Cd80",
        "Cd86",
        "Stat1"
    ),

    Anti_inflammatory = c(
        "Mrc1",
        "Cd163",
        "Il1rn",
        "Arg1",
        "Mertk",
        "Igf1",
        "Hmox1",
        "Klf4",
        "Maf"
    ),

    Fibrogenic = c(
        "Spp1",
        "Tgfb1",
        "Pdgfb",
        "Thbs1",
        "Lgals3",
        "Gpnmb",
        "Mmp12"
    ),

    Repair_Resolution = c(
        "Mertk",
        "Axl",
        "Mfge8",
        "Gas6",
        "Igf1",
        "Hmox1",
        "Mmp12",
        "Mmp13",
        "Mmp14",
        "Plau"
    ),

    Lipid_TREM2 = c(
        "Trem2",
        "Gpnmb",
        "Cd9",
        "Lpl",
        "Apoe",
        "Fabp5",
        "Abca1",
        "Plin2"
    ),

    IL10_STAT3 = c(
        "Il10ra",
        "Il10rb",
        "Jak1",
        "Tyk2",
        "Stat3",
        "Socs3",
        "Bcl3",
        "Il1rn"
    ),

    Efferocytosis = c(
        "Mertk",
        "Axl",
        "Mfge8",
        "Gas6",
        "Marco",
        "Cd36",
        "Lrp1",
        "C1qa",
        "C1qb",
        "C1qc"
    )
)

COMPACT_GENES <- c(
    "Il1b",
    "Tnf",
    "Ccl2",
    "Cxcl10",
    "Mrc1",
    "Cd163",
    "Il1rn",
    "Mertk",
    "Igf1",
    "Hmox1",
    "Spp1",
    "Tgfb1",
    "Pdgfb",
    "Thbs1",
    "Mfge8",
    "Gas6",
    "Mmp12",
    "Mmp13",
    "Trem2",
    "Gpnmb",
    "Cd9",
    "Lpl",
    "Apoe",
    "Stat3",
    "Socs3"
)

# ------------------------------------------------------------------------------
# 2. Packages
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
# 3. Helpers
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
# 4. Load Clean-B
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
    stop(
        "Clean-B RDS not found:\n",
        INPUT_RDS
    )
}

msg("Loading: ", INPUT_RDS)

mphi <- readRDS(INPUT_RDS)

DefaultAssay(mphi) <- ASSAY_USE

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

if (is.na(SAMPLE_COL)) stop("Sample column not found.")
if (is.na(CLASS_COL)) stop("MΦ class column not found.")

rna_data <- get_layer_safe(
    mphi,
    ASSAY_USE,
    "data"
)

if (is.null(rna_data)) {
    stop("Normalized RNA data layer not found.")
}

if (!identical(
    colnames(rna_data),
    colnames(mphi)
)) {
    stop("RNA data cell order mismatch.")
}

# ------------------------------------------------------------------------------
# 5. Canonical metadata
# ------------------------------------------------------------------------------

meta <- mphi@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = as.character(.data[[SAMPLE_COL]]),
        condition = canonical_condition(
            as.character(.data[[SAMPLE_COL]])
        ),
        macrophage_class = as.character(.data[[CLASS_COL]])
    ) %>%
    filter(
        macrophage_class %in%
            SUBTYPE_ORDER
    )

# ------------------------------------------------------------------------------
# 6. Detect marker genes
# ------------------------------------------------------------------------------

GENE_SETS_USE <- lapply(
    GENE_SETS,
    intersect,
    y = rownames(rna_data)
)

gene_membership <- bind_rows(
    lapply(
        names(GENE_SETS_USE),
        function(set_name) {
            tibble(
                module = set_name,
                gene = GENE_SETS_USE[[set_name]]
            )
        }
    )
)

write.csv(
    gene_membership,
    file.path(
        TAB_DIR,
        "01_gene_set_membership_v4.13.2.csv"
    ),
    row.names = FALSE
)

ALL_GENES <- unique(
    gene_membership$gene
)

COMPACT_USE <- intersect(
    COMPACT_GENES,
    ALL_GENES
)

# ------------------------------------------------------------------------------
# 7. Custom DotPlot: subtype × condition
# ------------------------------------------------------------------------------

group_condition_levels <- as.vector(
    outer(
        SUBTYPE_ORDER,
        CONDITION_ORDER,
        function(a, b) {
            paste(
                a,
                b,
                sep = " | "
            )
        }
    )
)

dot_condition_rows <- list()

for (group_now in group_condition_levels) {

    parts <- strsplit(
        group_now,
        " \\| "
    )[[1]]

    subtype_now <- parts[[1]]
    condition_now <- parts[[2]]

    cells <- meta$cell[
        meta$macrophage_class ==
            subtype_now &
        as.character(meta$condition) ==
            condition_now
    ]

    if (length(cells) == 0L) next

    mat <- rna_data[
        ALL_GENES,
        cells,
        drop = FALSE
    ]

    dot_condition_rows[[group_now]] <- tibble(
        group = group_now,
        gene = ALL_GENES,
        avg_expr = as.numeric(
            Matrix::rowMeans(mat)
        ),
        pct_expr = 100 *
            as.numeric(
                Matrix::rowMeans(
                    mat > 0
                )
            )
    )
}

dot_condition <- bind_rows(
    dot_condition_rows
) %>%
    group_by(gene) %>%
    mutate(
        avg_expr_z = {
            z <- as.numeric(scale(avg_expr))
            z[!is.finite(z)] <- 0
            pmax(
                pmin(z, 2.5),
                -2.5
            )
        }
    ) %>%
    ungroup() %>%
    mutate(
        group = factor(
            group,
            levels = group_condition_levels
        ),
        gene = factor(
            gene,
            levels = ALL_GENES
        )
    )

write.csv(
    dot_condition,
    file.path(
        TAB_DIR,
        "02_DotPlot_subtype_x_condition_numeric_v4.13.2.csv"
    ),
    row.names = FALSE
)

p_dot_condition <- ggplot(
    dot_condition,
    aes(
        x = gene,
        y = group
    )
) +
    geom_point(
        aes(
            size = pct_expr,
            color = avg_expr_z
        )
    ) +
    scale_size_continuous(
        range = c(0.3, 7),
        limits = c(0, 100),
        name = "% expressed"
    ) +
    scale_color_gradient2(
        low = "#0033FF",
        mid = "#FFFFFF",
        high = "#FF1A1A",
        midpoint = 0,
        limits = c(-2.5, 2.5),
        oob = scales::squish,
        name = "Average\nexpression\nz-score"
    ) +
    labs(
        title = "Clean-B MΦ gene-level functional validation",
        subtitle = "Subtype × condition | custom DotPlot",
        x = NULL,
        y = NULL
    ) +
    theme_classic(
        base_size = 8.5
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        axis.text.x = element_text(
            angle = 60,
            hjust = 1
        ),
        axis.text.y = element_text(
            size = 6.5
        )
    )

save_pdf(
    "01_DotPlot_subtype_x_condition_v4.13.2.pdf",
    p_dot_condition,
    22,
    10.5
)

# ------------------------------------------------------------------------------
# 8. Custom DotPlot: subtype × biological sample
# ------------------------------------------------------------------------------

group_sample_levels <- as.vector(
    outer(
        SUBTYPE_ORDER,
        SAMPLE_ORDER,
        function(a, b) {
            paste(
                a,
                b,
                sep = " | "
            )
        }
    )
)

dot_sample_rows <- list()

for (group_now in group_sample_levels) {

    parts <- strsplit(
        group_now,
        " \\| "
    )[[1]]

    subtype_now <- parts[[1]]
    sample_now <- parts[[2]]

    cells <- meta$cell[
        meta$macrophage_class ==
            subtype_now &
        meta$sample ==
            sample_now
    ]

    if (length(cells) == 0L) next

    mat <- rna_data[
        ALL_GENES,
        cells,
        drop = FALSE
    ]

    dot_sample_rows[[group_now]] <- tibble(
        group = group_now,
        gene = ALL_GENES,
        avg_expr = as.numeric(
            Matrix::rowMeans(mat)
        ),
        pct_expr = 100 *
            as.numeric(
                Matrix::rowMeans(
                    mat > 0
                )
            )
    )
}

dot_sample <- bind_rows(
    dot_sample_rows
) %>%
    group_by(gene) %>%
    mutate(
        avg_expr_z = {
            z <- as.numeric(scale(avg_expr))
            z[!is.finite(z)] <- 0
            pmax(
                pmin(z, 2.5),
                -2.5
            )
        }
    ) %>%
    ungroup() %>%
    mutate(
        group = factor(
            group,
            levels = group_sample_levels
        ),
        gene = factor(
            gene,
            levels = ALL_GENES
        )
    )

write.csv(
    dot_sample,
    file.path(
        TAB_DIR,
        "03_DotPlot_subtype_x_sample_numeric_v4.13.2.csv"
    ),
    row.names = FALSE
)

p_dot_sample <- ggplot(
    dot_sample,
    aes(
        x = gene,
        y = group
    )
) +
    geom_point(
        aes(
            size = pct_expr,
            color = avg_expr_z
        )
    ) +
    scale_size_continuous(
        range = c(0.2, 6),
        limits = c(0, 100),
        name = "% expressed"
    ) +
    scale_color_gradient2(
        low = "#0033FF",
        mid = "#FFFFFF",
        high = "#FF1A1A",
        midpoint = 0,
        limits = c(-2.5, 2.5),
        oob = scales::squish,
        name = "Average\nexpression\nz-score"
    ) +
    labs(
        title = "Clean-B MΦ gene-level validation by biological sample",
        subtitle = "Subtype × sample | custom DotPlot",
        x = NULL,
        y = NULL
    ) +
    theme_classic(
        base_size = 7.5
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        axis.text.x = element_text(
            angle = 60,
            hjust = 1
        ),
        axis.text.y = element_text(
            size = 5.8
        )
    )

save_pdf(
    "02_DotPlot_subtype_x_sample_v4.13.2.pdf",
    p_dot_sample,
    23,
    14
)

# ------------------------------------------------------------------------------
# 9. Sample-level mean-expression heatmap
# ------------------------------------------------------------------------------

sample_mean_rows <- list()

for (subtype_now in SUBTYPE_ORDER) {

    for (sample_now in SAMPLE_ORDER) {

        cells <- meta$cell[
            meta$macrophage_class ==
                subtype_now &
            meta$sample ==
                sample_now
        ]

        if (length(cells) == 0L) next

        mat <- rna_data[
            ALL_GENES,
            cells,
            drop = FALSE
        ]

        sample_mean_rows[[paste(
            subtype_now,
            sample_now,
            sep = "__"
        )]] <- tibble(
            macrophage_class = subtype_now,
            sample = sample_now,
            gene = ALL_GENES,
            mean_expr = as.numeric(
                Matrix::rowMeans(mat)
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
        "04_sample_level_mean_expression_v4.13.2.csv"
    ),
    row.names = FALSE
)

sample_heat <- sample_mean %>%
    mutate(
        column_id = paste(
            macrophage_class,
            sample,
            sep = " | "
        )
    ) %>%
    select(
        gene,
        column_id,
        mean_expr
    ) %>%
    pivot_wider(
        names_from = column_id,
        values_from = mean_expr
    )

sample_heat_mat <- as.matrix(
    sample_heat[
        ,
        -1,
        drop = FALSE
    ]
)

rownames(sample_heat_mat) <- sample_heat$gene

sample_heat_z <- t(
    scale(
        t(
            sample_heat_mat
        )
    )
)

sample_heat_z[
    !is.finite(
        sample_heat_z
    )
] <- 0

SAMPLE_HEAT_LIMIT <- 2.0

sample_heat_plot <- pmax(
    pmin(
        sample_heat_z,
        SAMPLE_HEAT_LIMIT
    ),
    -SAMPLE_HEAT_LIMIT
)

write.csv(
    data.frame(
        gene = rownames(sample_heat_z),
        sample_heat_z,
        check.names = FALSE
    ),
    file.path(
        TAB_DIR,
        "05_sample_level_gene_zscore_matrix_v4.13.2.csv"
    ),
    row.names = FALSE
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "03_sample_level_gene_heatmap_v4.13.2.pdf"
    ),
    width = 18,
    height = 12
)

pheatmap::pheatmap(
    sample_heat_plot,
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
        -SAMPLE_HEAT_LIMIT,
        SAMPLE_HEAT_LIMIT,
        length.out = 102
    ),
    border_color = "white",
    fontsize_row = 7,
    fontsize_col = 6,
    angle_col = 45,
    main = paste0(
        "Clean-B MΦ gene-level heatmap by biological sample\n",
        "row z-score | display clipped at +/-2"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 10. Condition-level heatmap
# ------------------------------------------------------------------------------

condition_mean <- sample_mean %>%
    mutate(
        condition = canonical_condition(
            sample
        )
    ) %>%
    group_by(
        macrophage_class,
        condition,
        gene
    ) %>%
    summarise(
        mean_expr = mean(
            mean_expr,
            na.rm = TRUE
        ),
        .groups = "drop"
    ) %>%
    mutate(
        column_id = paste(
            macrophage_class,
            condition,
            sep = " | "
        )
    ) %>%
    select(
        gene,
        column_id,
        mean_expr
    ) %>%
    pivot_wider(
        names_from = column_id,
        values_from = mean_expr
    )

condition_heat_mat <- as.matrix(
    condition_mean[
        ,
        -1,
        drop = FALSE
    ]
)

rownames(condition_heat_mat) <- condition_mean$gene

condition_heat_z <- t(
    scale(
        t(
            condition_heat_mat
        )
    )
)

condition_heat_z[
    !is.finite(
        condition_heat_z
    )
] <- 0

condition_heat_plot <- pmax(
    pmin(
        condition_heat_z,
        2
    ),
    -2
)

write.csv(
    data.frame(
        gene = rownames(condition_heat_z),
        condition_heat_z,
        check.names = FALSE
    ),
    file.path(
        TAB_DIR,
        "06_condition_level_gene_zscore_matrix_v4.13.2.csv"
    ),
    row.names = FALSE
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "04_condition_level_gene_heatmap_v4.13.2.pdf"
    ),
    width = 15,
    height = 11
)

pheatmap::pheatmap(
    condition_heat_plot,
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
    fontsize_row = 7,
    fontsize_col = 7,
    angle_col = 45,
    main = paste0(
        "Clean-B MΦ gene-level heatmap by subtype x condition\n",
        "row z-score | display clipped at +/-2"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 11. Compact representative-gene heatmap
# ------------------------------------------------------------------------------

compact_mat <- sample_heat_z[
    intersect(
        COMPACT_USE,
        rownames(sample_heat_z)
    ),
    ,
    drop = FALSE
]

compact_plot <- pmax(
    pmin(
        compact_mat,
        2
    ),
    -2
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "05_compact_representative_gene_heatmap_v4.13.2.pdf"
    ),
    width = 17,
    height = 8.5
)

pheatmap::pheatmap(
    compact_plot,
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
    fontsize_col = 6,
    angle_col = 45,
    main = paste0(
        "Representative functional genes across Clean-B MΦ subtypes\n",
        "biological sample-level mean expression | row z-score"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 12. README
# ------------------------------------------------------------------------------

readme <- c(
    "Mouse MASH MΦ gene-level DotPlot / heatmap v4.13.2",
    "",
    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",
    "Primary dataset:",
    "  Clean-B",
    "",
    "Purpose:",
    "  gene-level validation of MΦ functional programs",
    "",
    "DotPlot:",
    "  custom implementation",
    "  Seurat DotPlot() is NOT used",
    "  dot size = percent expressed",
    "  color = gene-wise scaled average normalized expression",
    "",
    "Heatmap:",
    "  biological-sample mean normalized expression",
    "  row z-score",
    "  display clipped at +/-2",
    "",
    "Modules:",
    "  Inflammatory",
    "  Anti-inflammatory",
    "  Fibrogenic",
    "  Repair/Resolution",
    "  Lipid/TREM2",
    "  IL10/STAT3",
    "  Efferocytosis",
    "",
    "Primary outputs:",
    "  01 subtype x condition DotPlot",
    "  02 subtype x sample DotPlot",
    "  03 sample-level heatmap",
    "  04 condition-level heatmap",
    "  05 compact representative-gene heatmap"
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.13.2.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.13.2.txt"
    )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)
