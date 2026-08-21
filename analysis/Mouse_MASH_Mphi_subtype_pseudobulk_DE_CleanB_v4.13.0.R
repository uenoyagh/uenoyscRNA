#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Clean-B subtype-specific sample-level pseudobulk DE
# v4.13.0
#
# PRIMARY DATASET
#   Clean-B from v4.11.7
#   - Res2 cluster 24 removed
#   - Res2 cluster 27 removed
#
# FIXED FRAMEWORK
#   - MΦ-only RPCA Res2.0
#   - v4.8.4 manual MΦ annotation unchanged
#   - NO reclustering
#   - biological sample is the pseudobulk replicate
#
# ANALYSES
#   A. STD vs CDAHFD
#      - n=1 vs n=1
#      - descriptive CPM / log2FC only
#      - NO p-values
#
#   B. Sham vs Tx
#      - Sham1, Sham20 vs Tx17, Tx5
#      - descriptive CPM / log2FC
#      - exploratory edgeR if available
#      - replicate consistency
#
# SUBTYPES
#   1. Inflammatory-Mphi
#   2. Anti-inflammatory-Mphi
#   3. Fibrogenic-Mphi
#   4. Repair/Resolution-Mphi
#   5. Lipid-associated/TREM2-Mphi
#
# IMPORTANT
#   - Pseudobulk is built explicitly using Matrix::rowSums(raw counts).
#   - No Seurat aggregation helper is used.
#   - Effect size + replicate consistency are primary.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4130)

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
    "Mphi_subtype_pseudobulk_DE_CleanB_v4.13.0"
)

FIG_DIR <- file.path(
    OUTPUT_DIR,
    "Figures"
)

TAB_DIR <- file.path(
    OUTPUT_DIR,
    "Tables"
)

PB_DIR <- file.path(
    OUTPUT_DIR,
    "Pseudobulk"
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
    PB_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    LOG_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

ASSAY_USE <- "RNA"

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

SAMPLE_ORDER <- c(
    "STD_rep1",
    "CDHFD_rep1",
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
)

DISEASE_SAMPLES <- c(
    "STD_rep1",
    "CDHFD_rep1"
)

SHAM_SAMPLES <- c(
    "Sham1",
    "Sham20"
)

TX_SAMPLES <- c(
    "Tx17",
    "Tx5"
)

CPM_PSEUDOCOUNT <- 1

TOP_N <- 20L

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
        "Missing required package(s): ",
        paste(
            missing_packages,
            collapse = ", "
        )
    )
}

HAS_EDGER <- requireNamespace(
    "edgeR",
    quietly = TRUE
)

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

if (HAS_EDGER) {

    suppressPackageStartupMessages(
        library(edgeR)
    )
}

# ------------------------------------------------------------------------------
# 2. Helpers
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

    out
}


calc_cpm <- function(
    count_matrix
) {

    lib <- colSums(
        count_matrix
    )

    if (any(
        lib <= 0
    )) {

        stop(
            "At least one pseudobulk library has zero total counts."
        )
    }

    sweep(
        count_matrix,
        2,
        lib / 1e6,
        "/"
    )
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


safe_filename <- function(
    x
) {

    gsub(
        "[^A-Za-z0-9]+",
        "_",
        x
    )
}


# ------------------------------------------------------------------------------
# 3. Load Clean-B
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

if (!inherits(
    mphi,
    "Seurat"
)) {

    stop(
        "Input is not a Seurat object."
    )
}

DefaultAssay(
    mphi
) <- ASSAY_USE


# ------------------------------------------------------------------------------
# 4. Detect metadata
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
        "v4.8.4 macrophage class metadata column not found."
    )
}

msg(
    "SAMPLE_COL = ",
    SAMPLE_COL
)

msg(
    "CLASS_COL = ",
    CLASS_COL
)


# ------------------------------------------------------------------------------
# 5. Raw counts
# ------------------------------------------------------------------------------

counts <- get_layer_safe(
    mphi,
    ASSAY_USE,
    "counts"
)

if (is.null(
    counts
)) {

    stop(
        "RNA counts layer not found."
    )
}

if (!inherits(
    counts,
    "Matrix"
)) {

    counts <- Matrix::Matrix(
        counts,
        sparse = TRUE
    )
}

if (!identical(
    colnames(counts),
    colnames(mphi)
)) {

    stop(
        "Counts cell order does not match Seurat object."
    )
}

msg(
    "Raw counts matrix: ",
    nrow(counts),
    " genes x ",
    ncol(counts),
    " cells"
)


# ------------------------------------------------------------------------------
# 6. Canonical metadata
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
    )

if (anyNA(
    meta$condition
)) {

    stop(
        "Could not assign condition for sample(s): ",
        paste(
            unique(
                meta$sample[
                    is.na(
                        meta$condition
                    )
                ]
            ),
            collapse = ", "
        )
    )
}

missing_samples <- setdiff(
    SAMPLE_ORDER,
    unique(
        meta$sample
    )
)

if (length(
    missing_samples
) > 0L) {

    stop(
        "Expected sample(s) missing: ",
        paste(
            missing_samples,
            collapse = ", "
        )
    )
}

missing_subtypes <- setdiff(
    SUBTYPE_ORDER,
    unique(
        meta$macrophage_class
    )
)

if (length(
    missing_subtypes
) > 0L) {

    stop(
        "Expected MΦ subtype(s) missing: ",
        paste(
            missing_subtypes,
            collapse = ", "
        )
    )
}


# ------------------------------------------------------------------------------
# 7. Cell-count QC
# ------------------------------------------------------------------------------

cell_count_qc <- meta %>%
    filter(
        macrophage_class %in%
            SUBTYPE_ORDER
    ) %>%
    count(
        sample,
        condition,
        macrophage_class,
        name = "n_cells"
    ) %>%
    arrange(
        factor(
            macrophage_class,
            levels =
                SUBTYPE_ORDER
        ),
        factor(
            sample,
            levels =
                SAMPLE_ORDER
        )
    )

write.csv(
    cell_count_qc,
    file.path(
        TAB_DIR,
        "00_cell_count_QC_sample_x_subtype_v4.13.0.csv"
    ),
    row.names = FALSE
)

print(
    cell_count_qc
)


# ------------------------------------------------------------------------------
# 8. SAFE pseudobulk builder
#
# Explicit raw-count aggregation:
#
#   pb <- Matrix::rowSums(
#       counts[, cells, drop = FALSE]
#   )
#
# ------------------------------------------------------------------------------

build_subtype_pseudobulk <- function(
    subtype_name
) {

    pb_list <- vector(
        "list",
        length(
            SAMPLE_ORDER
        )
    )

    names(
        pb_list
    ) <- SAMPLE_ORDER

    qc_list <- vector(
        "list",
        length(
            SAMPLE_ORDER
        )
    )

    for (
        i in seq_along(
            SAMPLE_ORDER
        )
    ) {

        smp <- SAMPLE_ORDER[[i]]

        cells <- meta$cell[
            meta$sample ==
                smp &
            meta$macrophage_class ==
                subtype_name
        ]

        cells <- intersect(
            cells,
            colnames(
                counts
            )
        )

        if (length(
            cells
        ) == 0L) {

            stop(
                "No cells found for subtype/sample: ",
                subtype_name,
                " / ",
                smp
            )
        }

        pb <- Matrix::rowSums(
            counts[
                ,
                cells,
                drop = FALSE
            ]
        )

        pb_list[[smp]] <- pb

        qc_list[[i]] <- tibble(
            macrophage_class =
                subtype_name,

            sample =
                smp,

            condition =
                canonical_condition(
                    smp
                ),

            n_cells =
                length(
                    cells
                ),

            library_size =
                sum(
                    pb
                ),

            mean_counts_per_cell =
                sum(
                    pb
                ) /
                length(
                    cells
                )
        )
    }

    pb_mat <- do.call(
        cbind,
        pb_list
    )

    rownames(
        pb_mat
    ) <- rownames(
        counts
    )

    colnames(
        pb_mat
    ) <- SAMPLE_ORDER

    # Safety: no two libraries should be identical.
    for (
        i in seq_len(
            ncol(
                pb_mat
            ) - 1L
        )
    ) {

        for (
            j in (
                i + 1L
            ):ncol(
                pb_mat
            )
        ) {

            if (identical(
                pb_mat[, i],
                pb_mat[, j]
            )) {

                stop(
                    "Identical pseudobulk libraries detected: ",
                    subtype_name,
                    " / ",
                    colnames(pb_mat)[[i]],
                    " vs ",
                    colnames(pb_mat)[[j]]
                )
            }
        }
    }

    list(
        counts = pb_mat,
        qc = bind_rows(
            qc_list
        )
    )
}


# ------------------------------------------------------------------------------
# 9. Build all subtype pseudobulk matrices
# ------------------------------------------------------------------------------

PB_OBJECTS <- list()

PB_QC_ALL <- list()

for (
    subtype_name in SUBTYPE_ORDER
) {

    msg(
        "Building pseudobulk: ",
        subtype_name
    )

    pb_now <- build_subtype_pseudobulk(
        subtype_name
    )

    PB_OBJECTS[[subtype_name]] <-
        pb_now$counts

    PB_QC_ALL[[subtype_name]] <-
        pb_now$qc

    subtype_safe <- safe_filename(
        subtype_name
    )

    write.csv(
        data.frame(
            gene = rownames(
                pb_now$counts
            ),
            pb_now$counts,
            check.names = FALSE
        ),
        file.path(
            PB_DIR,
            paste0(
                "PB_raw_counts_",
                subtype_safe,
                "_v4.13.0.csv"
            )
        ),
        row.names = FALSE
    )
}

pb_qc_all <- bind_rows(
    PB_QC_ALL
)

write.csv(
    pb_qc_all,
    file.path(
        PB_DIR,
        "01_pseudobulk_QC_all_subtypes_v4.13.0.csv"
    ),
    row.names = FALSE
)

msg(
    "All subtype pseudobulk matrices passed identity checks."
)


# ------------------------------------------------------------------------------
# 10. Descriptive DE function
# ------------------------------------------------------------------------------

run_descriptive_de <- function(
    pb_counts,
    subtype_name
) {

    pb_cpm <- calc_cpm(
        pb_counts
    )

    # Disease: one sample each.
    std_cpm <- pb_cpm[
        ,
        "STD_rep1"
    ]

    cdahfd_cpm <- pb_cpm[
        ,
        "CDHFD_rep1"
    ]

    disease_df <- tibble(
        macrophage_class =
            subtype_name,

        gene =
            rownames(
                pb_cpm
            ),

        CPM_STD =
            as.numeric(
                std_cpm
            ),

        CPM_CDAHFD =
            as.numeric(
                cdahfd_cpm
            ),

        log2FC_CDAHFD_vs_STD =
            log2(
                (
                    cdahfd_cpm +
                        CPM_PSEUDOCOUNT
                ) /
                (
                    std_cpm +
                        CPM_PSEUDOCOUNT
                )
            )
    ) %>%
        arrange(
            desc(
                abs(
                    log2FC_CDAHFD_vs_STD
                )
            )
        )

    # Treatment: mean of two biological samples.
    sham_mean <- rowMeans(
        pb_cpm[
            ,
            SHAM_SAMPLES,
            drop = FALSE
        ]
    )

    tx_mean <- rowMeans(
        pb_cpm[
            ,
            TX_SAMPLES,
            drop = FALSE
        ]
    )

    treatment_df <- tibble(
        macrophage_class =
            subtype_name,

        gene =
            rownames(
                pb_cpm
            ),

        CPM_Sham1 =
            as.numeric(
                pb_cpm[
                    ,
                    "Sham1"
                ]
            ),

        CPM_Sham20 =
            as.numeric(
                pb_cpm[
                    ,
                    "Sham20"
                ]
            ),

        CPM_Tx17 =
            as.numeric(
                pb_cpm[
                    ,
                    "Tx17"
                ]
            ),

        CPM_Tx5 =
            as.numeric(
                pb_cpm[
                    ,
                    "Tx5"
                ]
            ),

        mean_CPM_Sham =
            as.numeric(
                sham_mean
            ),

        mean_CPM_Tx =
            as.numeric(
                tx_mean
            ),

        log2FC_Tx_vs_Sham =
            log2(
                (
                    tx_mean +
                        CPM_PSEUDOCOUNT
                ) /
                (
                    sham_mean +
                        CPM_PSEUDOCOUNT
                )
            )
    ) %>%
        arrange(
            desc(
                abs(
                    log2FC_Tx_vs_Sham
                )
            )
        )

    list(
        cpm = pb_cpm,
        disease = disease_df,
        treatment = treatment_df
    )
}


# ------------------------------------------------------------------------------
# 11. Run descriptive DE for all subtypes
# ------------------------------------------------------------------------------

DESC_RESULTS <- list()

DISEASE_ALL <- list()

TREATMENT_ALL <- list()

for (
    subtype_name in SUBTYPE_ORDER
) {

    msg(
        "Descriptive DE: ",
        subtype_name
    )

    desc_now <- run_descriptive_de(
        PB_OBJECTS[[subtype_name]],
        subtype_name
    )

    DESC_RESULTS[[subtype_name]] <-
        desc_now

    DISEASE_ALL[[subtype_name]] <-
        desc_now$disease

    TREATMENT_ALL[[subtype_name]] <-
        desc_now$treatment

    subtype_safe <- safe_filename(
        subtype_name
    )

    write.csv(
        data.frame(
            gene =
                rownames(
                    desc_now$cpm
                ),
            desc_now$cpm,
            check.names =
                FALSE
        ),
        file.path(
            PB_DIR,
            paste0(
                "PB_CPM_",
                subtype_safe,
                "_v4.13.0.csv"
            )
        ),
        row.names = FALSE
    )

    write.csv(
        desc_now$disease,
        file.path(
            TAB_DIR,
            paste0(
                "DE_descriptive_CDAHFD_vs_STD_",
                subtype_safe,
                "_v4.13.0.csv"
            )
        ),
        row.names = FALSE
    )

    write.csv(
        desc_now$treatment,
        file.path(
            TAB_DIR,
            paste0(
                "DE_descriptive_Tx_vs_Sham_",
                subtype_safe,
                "_v4.13.0.csv"
            )
        ),
        row.names = FALSE
    )
}

disease_all <- bind_rows(
    DISEASE_ALL
)

treatment_all <- bind_rows(
    TREATMENT_ALL
)

write.csv(
    disease_all,
    file.path(
        TAB_DIR,
        "02_DE_descriptive_CDAHFD_vs_STD_all_subtypes_v4.13.0.csv"
    ),
    row.names = FALSE
)

write.csv(
    treatment_all,
    file.path(
        TAB_DIR,
        "03_DE_descriptive_Tx_vs_Sham_all_subtypes_v4.13.0.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 12. Treatment replicate consistency
#
# Compare each Tx sample against the mean of the two Sham samples.
# ------------------------------------------------------------------------------

CONSISTENCY_ALL <- list()

for (
    subtype_name in SUBTYPE_ORDER
) {

    pb_cpm <- DESC_RESULTS[[subtype_name]]$cpm

    log2cpm <- log2(
        pb_cpm +
            CPM_PSEUDOCOUNT
    )

    sham_mean_log <- rowMeans(
        log2cpm[
            ,
            SHAM_SAMPLES,
            drop = FALSE
        ]
    )

    tx17_delta <- log2cpm[
        ,
        "Tx17"
    ] -
        sham_mean_log

    tx5_delta <- log2cpm[
        ,
        "Tx5"
    ] -
        sham_mean_log

    consistency_now <- tibble(
        macrophage_class =
            subtype_name,

        gene =
            rownames(
                log2cpm
            ),

        Tx17_minus_ShamMean =
            as.numeric(
                tx17_delta
            ),

        Tx5_minus_ShamMean =
            as.numeric(
                tx5_delta
            ),

        direction_consistent =
            sign(
                Tx17_minus_ShamMean
            ) ==
            sign(
                Tx5_minus_ShamMean
            ),

        consistent_direction =
            case_when(

                direction_consistent &
                    Tx17_minus_ShamMean >
                        0 ~
                    "Up in both Tx",

                direction_consistent &
                    Tx17_minus_ShamMean <
                        0 ~
                    "Down in both Tx",

                TRUE ~
                    "Discordant"
            )
    ) %>%
        left_join(
            TREATMENT_ALL[[subtype_name]] %>%
                select(
                    gene,
                    log2FC_Tx_vs_Sham
                ),
            by = "gene"
        )

    CONSISTENCY_ALL[[subtype_name]] <-
        consistency_now

    write.csv(
        consistency_now,
        file.path(
            TAB_DIR,
            paste0(
                "Replicate_consistency_Tx_vs_Sham_",
                safe_filename(
                    subtype_name
                ),
                "_v4.13.0.csv"
            )
        ),
        row.names = FALSE
    )
}

consistency_all <- bind_rows(
    CONSISTENCY_ALL
)

write.csv(
    consistency_all,
    file.path(
        TAB_DIR,
        "04_Tx_vs_Sham_replicate_consistency_all_subtypes_v4.13.0.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 13. Exploratory edgeR for Sham vs Tx
# ------------------------------------------------------------------------------

EDGER_ALL <- list()

if (HAS_EDGER) {

    msg(
        "edgeR available: running exploratory Sham vs Tx DE."
    )

    group <- factor(
        c(
            "Sham",
            "Sham",
            "Tx",
            "Tx"
        ),
        levels = c(
            "Sham",
            "Tx"
        )
    )

    design <- model.matrix(
        ~ group
    )

    for (
        subtype_name in SUBTYPE_ORDER
    ) {

        pb_counts <- PB_OBJECTS[[subtype_name]][
            ,
            c(
                SHAM_SAMPLES,
                TX_SAMPLES
            ),
            drop = FALSE
        ]

        y <- edgeR::DGEList(
            counts = round(
                as.matrix(
                    pb_counts
                )
            ),
            group = group
        )

        keep <- edgeR::filterByExpr(
            y,
            group = group
        )

        y <- y[
            keep,
            ,
            keep.lib.sizes = FALSE
        ]

        y <- edgeR::calcNormFactors(
            y
        )

        y <- edgeR::estimateDisp(
            y,
            design,
            robust = TRUE
        )

        fit <- edgeR::glmQLFit(
            y,
            design,
            robust = TRUE
        )

        qlf <- edgeR::glmQLFTest(
            fit,
            coef = 2
        )

        edger_now <- edgeR::topTags(
            qlf,
            n = Inf,
            sort.by = "PValue"
        )$table %>%
            rownames_to_column(
                "gene"
            ) %>%
            rename(
                edgeR_log2FC =
                    logFC,

                edgeR_logCPM =
                    logCPM,

                edgeR_F =
                    F,

                edgeR_PValue =
                    PValue,

                edgeR_FDR =
                    FDR
            ) %>%
            mutate(
                macrophage_class =
                    subtype_name,
                .before =
                    gene
            ) %>%
            left_join(
                TREATMENT_ALL[[subtype_name]] %>%
                    select(
                        gene,
                        mean_CPM_Sham,
                        mean_CPM_Tx,
                        log2FC_Tx_vs_Sham
                    ),
                by = "gene"
            ) %>%
            left_join(
                CONSISTENCY_ALL[[subtype_name]] %>%
                    select(
                        gene,
                        Tx17_minus_ShamMean,
                        Tx5_minus_ShamMean,
                        direction_consistent,
                        consistent_direction
                    ),
                by = "gene"
            )

        EDGER_ALL[[subtype_name]] <-
            edger_now

        write.csv(
            edger_now,
            file.path(
                TAB_DIR,
                paste0(
                    "edgeR_exploratory_Tx_vs_Sham_",
                    safe_filename(
                        subtype_name
                    ),
                    "_v4.13.0.csv"
                )
            ),
            row.names = FALSE
        )
    }

    edger_all <- bind_rows(
        EDGER_ALL
    )

    write.csv(
        edger_all,
        file.path(
            TAB_DIR,
            "05_edgeR_exploratory_Tx_vs_Sham_all_subtypes_v4.13.0.csv"
        ),
        row.names = FALSE
    )

} else {

    msg(
        paste0(
            "edgeR not installed. ",
            "Descriptive DE and replicate-consistency analyses will still be produced."
        )
    )
}


# ------------------------------------------------------------------------------
# 14. Top-effect genes: disease
# ------------------------------------------------------------------------------

disease_top <- disease_all %>%
    group_by(
        macrophage_class
    ) %>%
    arrange(
        desc(
            abs(
                log2FC_CDAHFD_vs_STD
            )
        ),
        .by_group = TRUE
    ) %>%
    slice_head(
        n = TOP_N
    ) %>%
    ungroup()

write.csv(
    disease_top,
    file.path(
        TAB_DIR,
        "06_top_effect_genes_CDAHFD_vs_STD_v4.13.0.csv"
    ),
    row.names = FALSE
)

p_disease_top <- ggplot(
    disease_top,
    aes(
        x = log2FC_CDAHFD_vs_STD,
        y = reorder(
            gene,
            log2FC_CDAHFD_vs_STD
        )
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
            xend =
                log2FC_CDAHFD_vs_STD,
            yend =
                reorder(
                    gene,
                    log2FC_CDAHFD_vs_STD
                )
        ),
        linewidth = 0.5
    ) +
    geom_point(
        size = 2
    ) +
    facet_wrap(
        ~ macrophage_class,
        scales = "free_y",
        ncol = 3,
        labeller = as_labeller(
            SUBTYPE_LABELS
        )
    ) +
    labs(
        title =
            "Clean-B MΦ subtype pseudobulk: CDAHFD vs STD",
        subtitle =
            "Top absolute effect sizes | descriptive n=1 vs n=1",
        x =
            "log2FC (CDAHFD / STD)",
        y =
            NULL
    ) +
    theme_classic(
        base_size = 9
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            ),
        strip.text =
            element_text(
                face = "bold"
            )
    )

save_pdf(
    "01_top_effect_genes_CDAHFD_vs_STD_v4.13.0.pdf",
    p_disease_top,
    13,
    12
)


# ------------------------------------------------------------------------------
# 15. Top-effect genes: treatment
#
# Prefer genes consistent in BOTH Tx samples.
# If fewer than TOP_N are available, include strongest discordant genes.
# ------------------------------------------------------------------------------

treatment_consistency_merged <- treatment_all %>%
    left_join(
        consistency_all %>%
            select(
                macrophage_class,
                gene,
                Tx17_minus_ShamMean,
                Tx5_minus_ShamMean,
                direction_consistent,
                consistent_direction
            ),
        by = c(
            "macrophage_class",
            "gene"
        )
    )

treatment_top <- bind_rows(
    lapply(
        SUBTYPE_ORDER,
        function(
            subtype_name
        ) {

            df <- treatment_consistency_merged %>%
                filter(
                    macrophage_class ==
                        subtype_name
                )

            df_consistent <- df %>%
                filter(
                    direction_consistent
                ) %>%
                arrange(
                    desc(
                        abs(
                            log2FC_Tx_vs_Sham
                        )
                    )
                )

            if (nrow(
                df_consistent
            ) >= TOP_N) {

                return(
                    df_consistent %>%
                        slice_head(
                            n = TOP_N
                        )
                )
            }

            remaining_n <- TOP_N -
                nrow(
                    df_consistent
                )

            df_other <- df %>%
                filter(
                    !direction_consistent
                ) %>%
                arrange(
                    desc(
                        abs(
                            log2FC_Tx_vs_Sham
                        )
                    )
                ) %>%
                slice_head(
                    n = remaining_n
                )

            bind_rows(
                df_consistent,
                df_other
            )
        }
    )
)

write.csv(
    treatment_top,
    file.path(
        TAB_DIR,
        "07_top_effect_genes_Tx_vs_Sham_v4.13.0.csv"
    ),
    row.names = FALSE
)

p_treatment_top <- ggplot(
    treatment_top,
    aes(
        x = log2FC_Tx_vs_Sham,
        y = reorder(
            gene,
            log2FC_Tx_vs_Sham
        ),
        shape =
            consistent_direction
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
            xend =
                log2FC_Tx_vs_Sham,
            yend =
                reorder(
                    gene,
                    log2FC_Tx_vs_Sham
                )
        ),
        linewidth = 0.5
    ) +
    geom_point(
        size = 2
    ) +
    facet_wrap(
        ~ macrophage_class,
        scales = "free_y",
        ncol = 3,
        labeller = as_labeller(
            SUBTYPE_LABELS
        )
    ) +
    labs(
        title =
            "Clean-B MΦ subtype pseudobulk: Tx vs Sham",
        subtitle =
            "Top effect sizes | shape indicates biological-replicate consistency",
        x =
            "log2FC (Tx / Sham)",
        y =
            NULL,
        shape =
            "Tx replicate\nconsistency"
    ) +
    theme_classic(
        base_size = 9
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            ),
        strip.text =
            element_text(
                face = "bold"
            )
    )

save_pdf(
    "02_top_effect_genes_Tx_vs_Sham_v4.13.0.pdf",
    p_treatment_top,
    13,
    12
)


# ------------------------------------------------------------------------------
# 16. Replicate-consistency scatterplots
# ------------------------------------------------------------------------------

# Limit to genes with reasonable expression in Sham or Tx to reduce
# low-count visual noise.
consistency_plot_df <- treatment_consistency_merged %>%
    filter(
        mean_CPM_Sham >= 1 |
            mean_CPM_Tx >= 1
    )

p_consistency <- ggplot(
    consistency_plot_df,
    aes(
        x = Tx17_minus_ShamMean,
        y = Tx5_minus_ShamMean
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
        alpha = 0.35,
        size = 0.9
    ) +
    facet_wrap(
        ~ macrophage_class,
        scales = "free",
        ncol = 3,
        labeller = as_labeller(
            SUBTYPE_LABELS
        )
    ) +
    labs(
        title =
            "Tx replicate consistency by MΦ subtype",
        subtitle =
            "Each gene: Tx sample minus mean Sham log2(CPM+1)",
        x =
            "Tx17 − mean Sham",
        y =
            "Tx5 − mean Sham"
    ) +
    theme_classic(
        base_size = 9
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            ),
        strip.text =
            element_text(
                face = "bold"
            )
    )

save_pdf(
    "03_Tx_replicate_consistency_all_subtypes_v4.13.0.pdf",
    p_consistency,
    12,
    10
)


# ------------------------------------------------------------------------------
# 17. Exploratory edgeR volcano-style plots
# ------------------------------------------------------------------------------

if (HAS_EDGER) {

    edger_plot_df <- edger_all %>%
        mutate(
            minus_log10_P =
                -log10(
                    pmax(
                        edgeR_PValue,
                        1e-300
                    )
                )
        )

    p_edger <- ggplot(
        edger_plot_df,
        aes(
            x = edgeR_log2FC,
            y = minus_log10_P
        )
    ) +
        geom_vline(
            xintercept = 0,
            linetype = 2,
            linewidth = 0.35
        ) +
        geom_point(
            alpha = 0.40,
            size = 0.9
        ) +
        facet_wrap(
            ~ macrophage_class,
            scales = "free",
            ncol = 3,
            labeller = as_labeller(
                SUBTYPE_LABELS
            )
        ) +
        labs(
            title =
                "Exploratory edgeR pseudobulk DE: Tx vs Sham",
            subtitle =
                "n=2 Sham vs n=2 Tx | interpret p/FDR cautiously",
            x =
                "edgeR log2FC (Tx / Sham)",
            y =
                "-log10(P)"
        ) +
        theme_classic(
            base_size = 9
        ) +
        theme(
            plot.title =
                element_text(
                    face = "bold"
                ),
            strip.text =
                element_text(
                    face = "bold"
                )
        )

    save_pdf(
        "04_edgeR_exploratory_Tx_vs_Sham_all_subtypes_v4.13.0.pdf",
        p_edger,
        12,
        10
    )
}


# ------------------------------------------------------------------------------
# 18. Sample correlation heatmaps
# ------------------------------------------------------------------------------

correlation_tables <- list()

correlation_plots <- list()

for (
    subtype_name in SUBTYPE_ORDER
) {

    pb_cpm <- DESC_RESULTS[[subtype_name]]$cpm

    log2cpm <- log2(
        pb_cpm +
            CPM_PSEUDOCOUNT
    )

    variances <- apply(
        log2cpm,
        1,
        var
    )

    variable_genes <- names(
        variances
    )[
        variances > 0
    ]

    if (length(
        variable_genes
    ) < 10L) {

        variable_genes <- rownames(
            log2cpm
        )
    }

    cor_mat <- cor(
        log2cpm[
            variable_genes,
            ,
            drop = FALSE
        ],
        method = "pearson"
    )

    cor_df <- as.data.frame(
        as.table(
            cor_mat
        )
    )

    colnames(
        cor_df
    ) <- c(
        "Sample1",
        "Sample2",
        "Correlation"
    )

    cor_df$macrophage_class <-
        subtype_name

    correlation_tables[[subtype_name]] <-
        cor_df

    correlation_plots[[subtype_name]] <-
        ggplot(
            cor_df,
            aes(
                x = Sample1,
                y = Sample2,
                fill = Correlation
            )
        ) +
        geom_tile() +
        geom_text(
            aes(
                label =
                    sprintf(
                        "%.3f",
                        Correlation
                    )
            ),
            size = 3.2
        ) +
        scale_fill_gradient(
            low = "white",
            high = "black",
            limits = c(
                min(
                    cor_mat,
                    na.rm = TRUE
                ),
                1
            )
        ) +
        labs(
            title =
                SUBTYPE_LABELS[[subtype_name]],
            x = NULL,
            y = NULL
        ) +
        theme_classic(
            base_size = 8.5
        ) +
        theme(
            plot.title =
                element_text(
                    face = "bold",
                    hjust = 0.5
                ),
            axis.text.x =
                element_text(
                    angle = 45,
                    hjust = 1
                )
        )
}

correlation_all <- bind_rows(
    correlation_tables
)

write.csv(
    correlation_all,
    file.path(
        TAB_DIR,
        "08_sample_correlation_all_subtypes_v4.13.0.csv"
    ),
    row.names = FALSE
)

p_cor_all <- wrap_plots(
    correlation_plots,
    ncol = 2
) +
    plot_annotation(
        title =
            "Sample-level pseudobulk correlation by MΦ subtype"
    )

save_pdf(
    "05_sample_correlation_all_subtypes_v4.13.0.pdf",
    p_cor_all,
    12,
    15
)


# ------------------------------------------------------------------------------
# 19. Conservative high-priority treatment DE table
#
# Primary prioritization:
#   - mean CPM >= 1 in Sham or Tx
#   - |log2FC| >= 0.5
#   - both Tx samples change in same direction
#
# edgeR FDR is added if available, but is NOT required.
# ------------------------------------------------------------------------------

priority_treatment <- treatment_consistency_merged %>%
    filter(
        mean_CPM_Sham >= 1 |
            mean_CPM_Tx >= 1
    ) %>%
    filter(
        abs(
            log2FC_Tx_vs_Sham
        ) >= 0.5
    ) %>%
    filter(
        direction_consistent
    ) %>%
    arrange(
        macrophage_class,
        desc(
            abs(
                log2FC_Tx_vs_Sham
            )
        )
    )

if (HAS_EDGER) {

    priority_treatment <- priority_treatment %>%
        left_join(
            edger_all %>%
                select(
                    macrophage_class,
                    gene,
                    edgeR_log2FC,
                    edgeR_PValue,
                    edgeR_FDR
                ),
            by = c(
                "macrophage_class",
                "gene"
            )
        )
}

write.csv(
    priority_treatment,
    file.path(
        TAB_DIR,
        "09_high_priority_Tx_vs_Sham_genes_v4.13.0.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 20. Conservative disease-effect table
#
# n=1 vs n=1:
#   no p-value
#   filter only on expression and effect size
# ------------------------------------------------------------------------------

priority_disease <- disease_all %>%
    filter(
        CPM_STD >= 1 |
            CPM_CDAHFD >= 1
    ) %>%
    filter(
        abs(
            log2FC_CDAHFD_vs_STD
        ) >= 0.5
    ) %>%
    arrange(
        macrophage_class,
        desc(
            abs(
                log2FC_CDAHFD_vs_STD
            )
        )
    )

write.csv(
    priority_disease,
    file.path(
        TAB_DIR,
        "10_high_priority_CDAHFD_vs_STD_genes_v4.13.0.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 21. Summary counts
# ------------------------------------------------------------------------------

summary_treatment <- priority_treatment %>%
    group_by(
        macrophage_class
    ) %>%
    summarise(
        n_priority_genes =
            n(),

        n_up =
            sum(
                log2FC_Tx_vs_Sham >
                    0
            ),

        n_down =
            sum(
                log2FC_Tx_vs_Sham <
                    0
            ),

        median_abs_log2FC =
            median(
                abs(
                    log2FC_Tx_vs_Sham
                ),
                na.rm = TRUE
            ),

        .groups =
            "drop"
    )

summary_disease <- priority_disease %>%
    group_by(
        macrophage_class
    ) %>%
    summarise(
        n_priority_genes =
            n(),

        n_up =
            sum(
                log2FC_CDAHFD_vs_STD >
                    0
            ),

        n_down =
            sum(
                log2FC_CDAHFD_vs_STD <
                    0
            ),

        median_abs_log2FC =
            median(
                abs(
                    log2FC_CDAHFD_vs_STD
                ),
                na.rm = TRUE
            ),

        .groups =
            "drop"
    )

write.csv(
    summary_treatment,
    file.path(
        TAB_DIR,
        "11_summary_Tx_vs_Sham_priority_genes_v4.13.0.csv"
    ),
    row.names = FALSE
)

write.csv(
    summary_disease,
    file.path(
        TAB_DIR,
        "12_summary_CDAHFD_vs_STD_priority_genes_v4.13.0.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------------------------
# 22. Master figure
# ------------------------------------------------------------------------------

master_plots <- list(
    p_disease_top,
    p_treatment_top,
    p_consistency,
    p_cor_all
)

if (HAS_EDGER) {

    master_plots <- append(
        master_plots,
        list(
            p_edger
        )
    )
}

p_master <- wrap_plots(
    master_plots,
    ncol = 1
) +
    plot_annotation(
        title =
            "Clean-B MΦ subtype-specific sample-level pseudobulk DE v4.13.0",
        subtitle =
            paste0(
                "Disease: STD vs CDAHFD descriptive | ",
                "Treatment: Sham1/Sham20 vs Tx17/Tx5"
            ),
        theme =
            theme(
                plot.title =
                    element_text(
                        face = "bold",
                        size = 18
                    )
            )
    )

save_pdf(
    "06_Mphi_subtype_pseudobulk_DE_master_v4.13.0.pdf",
    p_master,
    15,
    ifelse(
        HAS_EDGER,
        46,
        36
    )
)


# ------------------------------------------------------------------------------
# 23. README
# ------------------------------------------------------------------------------

readme <- c(

    "Mouse MASH MΦ subtype-specific sample-level pseudobulk DE v4.13.0",
    "",

    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",

    "Primary dataset:",
    "  Clean-B from v4.11.7",
    "  Res2 cluster 24 removed",
    "  Res2 cluster 27 removed",
    "  v4.8.4 macrophage classification unchanged",
    "",

    "Subtypes:",
    "  Inflammatory-MΦ",
    "  Anti-inflammatory-MΦ",
    "  Fibrogenic-MΦ",
    "  Repair/Resolution-MΦ",
    "  Lipid-associated/TREM2-MΦ",
    "",

    "Pseudobulk:",
    "  biological sample = replicate",
    "  raw counts explicitly summed with Matrix::rowSums",
    "  no Seurat aggregation helper used",
    "",

    "Disease analysis:",
    "  STD_rep1 vs CDHFD_rep1",
    "  descriptive CPM and log2FC only",
    "  no p-values because n=1 vs n=1",
    "",

    "Treatment analysis:",
    "  Sham1/Sham20 vs Tx17/Tx5",
    "  descriptive CPM and log2FC",
    "  replicate consistency",
    paste0(
        "  exploratory edgeR: ",
        ifelse(
            HAS_EDGER,
            "ON",
            "OFF"
        )
    ),
    "",

    "Primary interpretation:",
    "  prioritize effect size + biological replicate consistency",
    "  edgeR P/FDR are exploratory because n=2 vs n=2",
    "",

    "High-priority Tx/Sham gene rule:",
    "  mean CPM >= 1 in Sham or Tx",
    "  abs(log2FC) >= 0.5",
    "  Tx17 and Tx5 move in the same direction relative to mean Sham",
    "",

    "Primary outputs:",
    "  01 top disease effect genes",
    "  02 top treatment effect genes",
    "  03 Tx replicate consistency",
    "  04 exploratory edgeR",
    "  05 pseudobulk sample correlations",
    "  06 master figure",
    "",

    "Key tables:",
    "  02 disease DE all subtypes",
    "  03 treatment DE all subtypes",
    "  04 replicate consistency",
    "  05 exploratory edgeR all subtypes",
    "  09 high-priority Tx/Sham genes",
    "  10 high-priority CDAHFD/STD genes",
    "",

    "Next planned analysis:",
    "  formal subtype compositional analysis v4.13.1"
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.13.0.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.13.0.txt"
    )
)


# ------------------------------------------------------------------------------
# 24. Final
# ------------------------------------------------------------------------------

msg(
    "DONE."
)

msg(
    "Output: ",
    OUTPUT_DIR
)

msg(
    "edgeR exploratory treatment analysis: ",
    ifelse(
        HAS_EDGER,
        "ON",
        "OFF"
    )
)

print(
    summary_treatment
)

print(
    summary_disease
)
