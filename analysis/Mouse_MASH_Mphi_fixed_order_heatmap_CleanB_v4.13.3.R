#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Publication-oriented fixed-order representative-gene heatmap
# v4.13.3
#
# PURPOSE
#   Rebuild the compact representative-gene heatmap from v4.13.2 for
#   publication use.
#
# DESIGN
#   - Clean-B primary dataset
#   - NO row clustering
#   - NO column clustering
#   - fixed biological gene order
#   - fixed MΦ subtype order
#   - fixed biological-sample order
#   - row z-score
#   - display clipped at +/-2
#
# OUTPUT
#   01 all-sample fixed-order heatmap
#   02 Sham/Tx treatment-focused fixed-order heatmap
#   03 condition-mean fixed-order heatmap
#   04 numeric z-score matrix
#
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4133)

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
    "Mphi_fixed_order_heatmap_CleanB_v4.13.3"
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

# ------------------------------------------------------------------------------
# 1. Fixed biological order
# ------------------------------------------------------------------------------

SUBTYPE_ORDER <- c(
    "Inflammatory-Mphi",
    "Anti-inflammatory-Mphi",
    "Fibrogenic-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi"
)

SUBTYPE_LABELS <- c(
    "Inflammatory-Mphi" =
        "Inflammatory-MΦ",
    "Anti-inflammatory-Mphi" =
        "Anti-inflammatory-MΦ",
    "Fibrogenic-Mphi" =
        "Fibrogenic-MΦ",
    "Repair/Resolution-Mphi" =
        "Repair/Resolution-MΦ",
    "Lipid-associated/TREM2-Mphi" =
        "Lipid-associated/TREM2-MΦ"
)

SAMPLE_ORDER_ALL <- c(
    "STD_rep1",
    "CDHFD_rep1",
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
)

SAMPLE_ORDER_TREATMENT <- c(
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

# ------------------------------------------------------------------------------
# 2. Publication representative genes
#
# The same genes are kept in a biologically interpretable fixed order.
# Duplicated genes are assigned once only, to the first relevant program.
# ------------------------------------------------------------------------------

GENE_BLOCKS <- list(

    Inflammatory = c(
        "Il1b",
        "Tnf",
        "Ccl2",
        "Cxcl10"
    ),

    Anti_inflammatory = c(
        "Mrc1",
        "Cd163",
        "Il1rn",
        "Mertk",
        "Igf1",
        "Hmox1"
    ),

    IL10_STAT3 = c(
        "Stat3",
        "Socs3"
    ),

    Repair_Resolution = c(
        "Mfge8",
        "Gas6",
        "Mmp13"
    ),

    Fibrogenic = c(
        "Spp1",
        "Tgfb1",
        "Pdgfb",
        "Thbs1",
        "Mmp12"
    ),

    Lipid_TREM2 = c(
        "Trem2",
        "Gpnmb",
        "Cd9",
        "Lpl",
        "Apoe"
    )
)

# Flatten while preserving first occurrence.
GENE_ORDER <- unique(
    unlist(
        GENE_BLOCKS,
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
    library(pheatmap)
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

    if (length(hit) == 0L) {

        return(
            NA_character_
        )
    }

    hit[[1]]
}


canonical_condition <- function(
    sample_name
) {

    x <- as.character(
        sample_name
    )

    out <- rep(
        NA_character_,
        length(x)
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
            "CDAHFD|CDHFD",
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


# ------------------------------------------------------------------------------
# 5. Load Clean-B
# ------------------------------------------------------------------------------

if (!file.exists(
    INPUT_RDS
)) {

    stop(
        "Clean-B RDS not found:\n",
        INPUT_RDS
    )
}

msg(
    "Loading Clean-B: ",
    INPUT_RDS
)

mphi <- readRDS(
    INPUT_RDS
)

DefaultAssay(
    mphi
) <- ASSAY_USE

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

CLASS_COL <- first_existing(
    meta_cols,
    c(
        "macrophage_class_Res2_v484",
        "macrophage_class_v484",
        "manual_class_v484",
        "macrophage_class"
    )
)

if (is.na(
    SAMPLE_COL
)) {

    stop(
        "Sample metadata column not found."
    )
}

if (is.na(
    CLASS_COL
)) {

    stop(
        "MΦ class metadata column not found."
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
        "Normalized RNA data layer not found."
    )
}

if (!identical(
    colnames(rna_data),
    colnames(mphi)
)) {

    stop(
        "RNA data cell order does not match Seurat object."
    )
}

# ------------------------------------------------------------------------------
# 6. Metadata
# ------------------------------------------------------------------------------

meta <- mphi@meta.data %>%
    rownames_to_column(
        "cell"
    ) %>%
    transmute(
        cell = cell,

        sample = as.character(
            .data[[SAMPLE_COL]]
        ),

        condition = canonical_condition(
            as.character(
                .data[[SAMPLE_COL]]
            )
        ),

        macrophage_class = as.character(
            .data[[CLASS_COL]]
        )
    ) %>%
    filter(
        macrophage_class %in%
            SUBTYPE_ORDER
    )

# ------------------------------------------------------------------------------
# 7. Gene availability
# ------------------------------------------------------------------------------

GENE_USE <- intersect(
    GENE_ORDER,
    rownames(
        rna_data
    )
)

missing_genes <- setdiff(
    GENE_ORDER,
    GENE_USE
)

gene_audit <- bind_rows(
    lapply(
        names(
            GENE_BLOCKS
        ),
        function(
            block_name
        ) {

            tibble(
                block =
                    block_name,

                gene =
                    GENE_BLOCKS[[block_name]],

                detected =
                    GENE_BLOCKS[[block_name]] %in%
                    rownames(
                        rna_data
                    )
            )
        }
    )
)

write.csv(
    gene_audit,
    file.path(
        TAB_DIR,
        "01_gene_availability_v4.13.3.csv"
    ),
    row.names = FALSE
)

msg(
    "Representative genes detected: ",
    length(
        GENE_USE
    ),
    " / ",
    length(
        GENE_ORDER
    )
)

if (length(
    missing_genes
) > 0L) {

    msg(
        "Missing genes: ",
        paste(
            missing_genes,
            collapse = ", "
        )
    )
}

# ------------------------------------------------------------------------------
# 8. Sample-level mean expression
# ------------------------------------------------------------------------------

sample_mean_rows <- list()

for (
    subtype_now in SUBTYPE_ORDER
) {

    for (
        sample_now in SAMPLE_ORDER_ALL
    ) {

        cells <- meta$cell[
            meta$macrophage_class ==
                subtype_now &
            meta$sample ==
                sample_now
        ]

        if (length(
            cells
        ) == 0L) {

            next
        }

        mat <- rna_data[
            GENE_USE,
            cells,
            drop = FALSE
        ]

        sample_mean_rows[[paste(
            subtype_now,
            sample_now,
            sep = "__"
        )]] <- tibble(
            macrophage_class =
                subtype_now,

            sample =
                sample_now,

            condition =
                as.character(
                    canonical_condition(
                        sample_now
                    )
                ),

            gene =
                GENE_USE,

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
        "02_sample_level_mean_expression_v4.13.3.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Build fixed-order matrix
# ------------------------------------------------------------------------------

sample_mean <- sample_mean %>%
    mutate(
        column_id = paste(
            macrophage_class,
            sample,
            sep = " | "
        )
    )

COLUMN_ORDER_ALL <- as.vector(
    unlist(
        lapply(
            SUBTYPE_ORDER,
            function(
                subtype_now
            ) {

                paste(
                    subtype_now,
                    SAMPLE_ORDER_ALL,
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
        names_from =
            column_id,
        values_from =
            mean_expr
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

available_columns <- intersect(
    COLUMN_ORDER_ALL,
    colnames(
        heat_mat
    )
)

heat_mat <- heat_mat[
    intersect(
        GENE_ORDER,
        rownames(
            heat_mat
        )
    ),
    available_columns,
    drop = FALSE
]

# ------------------------------------------------------------------------------
# 10. Gene-wise z-score
# ------------------------------------------------------------------------------

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

HEAT_LIMIT <- 2.0

heat_plot <- pmax(
    pmin(
        heat_z,
        HEAT_LIMIT
    ),
    -HEAT_LIMIT
)

write.csv(
    data.frame(
        gene =
            rownames(
                heat_z
            ),
        heat_z,
        check.names =
            FALSE
    ),
    file.path(
        TAB_DIR,
        "03_fixed_order_gene_zscore_matrix_v4.13.3.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. Column annotation
# ------------------------------------------------------------------------------

column_annotation <- tibble(
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
        condition =
            as.character(
                canonical_condition(
                    sample
                )
            )
    )

annotation_col <- column_annotation %>%
    select(
        macrophage_class,
        condition
    ) %>%
    as.data.frame()

rownames(
    annotation_col
) <- column_annotation$column_id

annotation_colors <- list(
    macrophage_class = c(
        "Inflammatory-Mphi" =
            "#E41A1C",
        "Anti-inflammatory-Mphi" =
            "#00AEEF",
        "Fibrogenic-Mphi" =
            "#FF1493",
        "Repair/Resolution-Mphi" =
            "#00C853",
        "Lipid-associated/TREM2-Mphi" =
            "#7B2CBF"
    ),

    condition = c(
        "STD" =
            "#0066FF",
        "CDAHFD" =
            "#FF1A1A",
        "Sham" =
            "#FF1A1A",
        "Tx" =
            "#0066FF"
    )
)

# ------------------------------------------------------------------------------
# 12. Row gaps based on biological programs
# ------------------------------------------------------------------------------

detected_by_block <- lapply(
    GENE_BLOCKS,
    function(
        genes
    ) {

        intersect(
            genes,
            rownames(
                heat_plot
            )
        )
    }
)

block_lengths <- lengths(
    detected_by_block
)

row_gaps <- cumsum(
    block_lengths
)

row_gaps <- row_gaps[
    row_gaps <
        nrow(
            heat_plot
        )
]

# Column gaps after each subtype block.
n_sample_all <- length(
    SAMPLE_ORDER_ALL
)

col_gaps_all <- seq(
    n_sample_all,
    ncol(
        heat_plot
    ) -
        n_sample_all,
    by =
        n_sample_all
)

# ------------------------------------------------------------------------------
# 13. Figure 01: all samples, fixed order
# ------------------------------------------------------------------------------

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "01_fixed_order_representative_gene_heatmap_all_samples_v4.13.3.pdf"
    ),
    width = 18,
    height = 8.5
)

pheatmap::pheatmap(
    heat_plot,

    cluster_rows =
        FALSE,

    cluster_cols =
        FALSE,

    gaps_row =
        row_gaps,

    gaps_col =
        col_gaps_all,

    annotation_col =
        annotation_col[
            colnames(
                heat_plot
            ),
            ,
            drop = FALSE
        ],

    annotation_colors =
        annotation_colors,

    color =
        grDevices::colorRampPalette(
            c(
                "#0033FF",
                "#FFFFFF",
                "#FF1A1A"
            )
        )(101),

    breaks =
        seq(
            -HEAT_LIMIT,
            HEAT_LIMIT,
            length.out = 102
        ),

    border_color =
        "white",

    fontsize_row =
        8.5,

    fontsize_col =
        6.2,

    angle_col =
        45,

    main =
        paste0(
            "Representative functional genes across Clean-B MΦ subtypes\n",
            "fixed biological order | sample-level mean expression | row z-score"
        )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 14. Figure 02: Sham / Tx only
# ------------------------------------------------------------------------------

COLUMN_ORDER_TREATMENT <- as.vector(
    unlist(
        lapply(
            SUBTYPE_ORDER,
            function(
                subtype_now
            ) {

                paste(
                    subtype_now,
                    SAMPLE_ORDER_TREATMENT,
                    sep = " | "
                )
            }
        ),
        use.names = FALSE
    )
)

treatment_columns <- intersect(
    COLUMN_ORDER_TREATMENT,
    colnames(
        heat_plot
    )
)

heat_treatment <- heat_plot[
    ,
    treatment_columns,
    drop = FALSE
]

annotation_treatment <- annotation_col[
    treatment_columns,
    ,
    drop = FALSE
]

n_sample_treatment <- length(
    SAMPLE_ORDER_TREATMENT
)

col_gaps_treatment <- seq(
    n_sample_treatment,
    ncol(
        heat_treatment
    ) -
        n_sample_treatment,
    by =
        n_sample_treatment
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "02_fixed_order_representative_gene_heatmap_Sham_Tx_v4.13.3.pdf"
    ),
    width = 14,
    height = 8.5
)

pheatmap::pheatmap(
    heat_treatment,

    cluster_rows =
        FALSE,

    cluster_cols =
        FALSE,

    gaps_row =
        row_gaps,

    gaps_col =
        col_gaps_treatment,

    annotation_col =
        annotation_treatment,

    annotation_colors =
        annotation_colors,

    color =
        grDevices::colorRampPalette(
            c(
                "#0033FF",
                "#FFFFFF",
                "#FF1A1A"
            )
        )(101),

    breaks =
        seq(
            -HEAT_LIMIT,
            HEAT_LIMIT,
            length.out = 102
        ),

    border_color =
        "white",

    fontsize_row =
        8.5,

    fontsize_col =
        7,

    angle_col =
        45,

    main =
        paste0(
            "Treatment-associated MΦ functional remodeling\n",
            "Sham1 / Sham20 → Tx17 / Tx5 | fixed subtype and gene order"
        )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 15. Figure 03: condition means
# ------------------------------------------------------------------------------

condition_mean <- sample_mean %>%
    mutate(
        condition = factor(
            condition,
            levels = CONDITION_ORDER
        )
    ) %>%
    group_by(
        macrophage_class,
        condition,
        gene
    ) %>%
    summarise(
        mean_expr =
            mean(
                mean_expr,
                na.rm = TRUE
            ),
        .groups =
            "drop"
    ) %>%
    mutate(
        column_id =
            paste(
                macrophage_class,
                condition,
                sep = " | "
            )
    )

CONDITION_COLUMN_ORDER <- as.vector(
    unlist(
        lapply(
            SUBTYPE_ORDER,
            function(
                subtype_now
            ) {

                paste(
                    subtype_now,
                    CONDITION_ORDER,
                    sep = " | "
                )
            }
        ),
        use.names = FALSE
    )
)

condition_wide <- condition_mean %>%
    select(
        gene,
        column_id,
        mean_expr
    ) %>%
    pivot_wider(
        names_from =
            column_id,
        values_from =
            mean_expr
    )

condition_mat <- as.matrix(
    condition_wide[
        ,
        -1,
        drop = FALSE
    ]
)

rownames(
    condition_mat
) <- condition_wide$gene

condition_mat <- condition_mat[
    intersect(
        GENE_ORDER,
        rownames(
            condition_mat
        )
    ),
    intersect(
        CONDITION_COLUMN_ORDER,
        colnames(
            condition_mat
        )
    ),
    drop = FALSE
]

condition_z <- t(
    scale(
        t(
            condition_mat
        )
    )
)

condition_z[
    !is.finite(
        condition_z
    )
] <- 0

condition_plot <- pmax(
    pmin(
        condition_z,
        HEAT_LIMIT
    ),
    -HEAT_LIMIT
)

condition_annotation <- tibble(
    column_id =
        colnames(
            condition_plot
        )
) %>%
    separate(
        column_id,
        into = c(
            "macrophage_class",
            "condition"
        ),
        sep = " \\| ",
        remove = FALSE
    )

annotation_condition <- condition_annotation %>%
    select(
        macrophage_class,
        condition
    ) %>%
    as.data.frame()

rownames(
    annotation_condition
) <- condition_annotation$column_id

n_condition <- length(
    CONDITION_ORDER
)

col_gaps_condition <- seq(
    n_condition,
    ncol(
        condition_plot
    ) -
        n_condition,
    by =
        n_condition
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "03_fixed_order_representative_gene_heatmap_condition_means_v4.13.3.pdf"
    ),
    width = 14,
    height = 8.5
)

pheatmap::pheatmap(
    condition_plot,

    cluster_rows =
        FALSE,

    cluster_cols =
        FALSE,

    gaps_row =
        row_gaps,

    gaps_col =
        col_gaps_condition,

    annotation_col =
        annotation_condition,

    annotation_colors =
        annotation_colors,

    color =
        grDevices::colorRampPalette(
            c(
                "#0033FF",
                "#FFFFFF",
                "#FF1A1A"
            )
        )(101),

    breaks =
        seq(
            -HEAT_LIMIT,
            HEAT_LIMIT,
            length.out = 102
        ),

    border_color =
        "white",

    fontsize_row =
        8.5,

    fontsize_col =
        7,

    angle_col =
        45,

    main =
        paste0(
            "Representative functional genes across Clean-B MΦ subtypes\n",
            "fixed subtype / condition / gene order"
        )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 16. README
# ------------------------------------------------------------------------------

readme <- c(

    "Mouse MASH MΦ fixed-order representative-gene heatmap v4.13.3",
    "",

    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",

    "Purpose:",
    "  publication-oriented replacement of v4.13.2 compact heatmap",
    "",

    "Design:",
    "  no row clustering",
    "  no column clustering",
    "  fixed biological gene order",
    "  fixed MΦ subtype order",
    "  fixed biological sample order",
    "  row z-score",
    "  display clipped at +/-2",
    "",

    "Gene-block order:",
    "  Inflammatory",
    "  Anti-inflammatory",
    "  IL10/STAT3",
    "  Repair/Resolution",
    "  Fibrogenic",
    "  Lipid/TREM2",
    "",

    "Subtype order:",
    "  Inflammatory-MΦ",
    "  Anti-inflammatory-MΦ",
    "  Fibrogenic-MΦ",
    "  Repair/Resolution-MΦ",
    "  Lipid-associated/TREM2-MΦ",
    "",

    "Primary publication candidate:",
    "  02_fixed_order_representative_gene_heatmap_Sham_Tx_v4.13.3.pdf",
    "",

    "Other outputs:",
    "  01 all biological samples",
    "  03 condition means",
    "",

    "Important:",
    "  sample-level Sham/Tx version retains biological replicate visibility.",
    "  condition-mean version is for visual summary only."
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.13.3.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.13.3.txt"
    )
)

msg(
    "DONE."
)

msg(
    "Output: ",
    OUTPUT_DIR
)

msg(
    "Primary publication candidate: ",
    file.path(
        FIG_DIR,
        "02_fixed_order_representative_gene_heatmap_Sham_Tx_v4.13.3.pdf"
    )
)
