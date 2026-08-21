#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Final annotation update after Fibrogenic-focused audit
# Clean-B
# v4.14.5
#
# PURPOSE
#   Update the biological annotation name:
#
#     Fibrogenic-MΦ
#       ->
#     ECM-associated inflammatory-MΦ
#
#   based on v4.14.4 focused audit.
#
# IMPORTANT
#   - Res2.0 clustering is NOT changed.
#   - Cell membership is NOT changed.
#   - Historical v4.8.4 / v4.10-v4.14 annotation columns are preserved.
#   - A new FINAL annotation column is added.
#   - Historical RDS files are NOT overwritten.
#
# FINAL FIVE-CLASS FRAMEWORK
#   1) Inflammatory-MΦ
#   2) Anti-inflammatory-MΦ
#   3) ECM-associated inflammatory-MΦ
#   4) Repair/Resolution-MΦ
#   5) Lipid-associated/TREM2-MΦ
#   + Other
#
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4145)

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
    "Mphi_Res2_CleanB_FINAL_v4.14.5"
)

RDS_DIR <- file.path(
    OUTPUT_DIR,
    "RDS"
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
    RDS_DIR,
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

OUTPUT_RDS <- file.path(
    RDS_DIR,
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

# ------------------------------------------------------------------------------
# 1. Final naming policy
# ------------------------------------------------------------------------------

OLD_CLASS_NAME <- "Fibrogenic-Mphi"

NEW_CLASS_NAME <- "ECM-associated inflammatory-Mphi"

FINAL_CLASS_COL <- "macrophage_class_Res2_FINAL_v4145"

FINAL_CLASS_ORDER <- c(
    "Inflammatory-Mphi",
    "Anti-inflammatory-Mphi",
    "ECM-associated inflammatory-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi",
    "Other"
)

# ------------------------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
    "Seurat",
    "SeuratObject",
    "dplyr",
    "tibble"
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
    library(dplyr)
    library(tibble)
})

# ------------------------------------------------------------------------------
# 3. Helpers
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

# ------------------------------------------------------------------------------
# 4. Load Clean-B
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
    "Loading Clean-B RDS: ",
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

meta_cols <- colnames(
    mphi@meta.data
)

# ------------------------------------------------------------------------------
# 5. Detect source annotation column
# ------------------------------------------------------------------------------

SOURCE_CLASS_COL <- first_existing(
    meta_cols,
    c(
        "macrophage_class_Res2_v484",
        "macrophage_class_v484",
        "manual_class_v484",
        "macrophage_class"
    )
)

if (is.na(
    SOURCE_CLASS_COL
)) {

    stop(
        "Could not identify historical macrophage-class column."
    )
}

msg(
    "Historical source annotation column: ",
    SOURCE_CLASS_COL
)

# ------------------------------------------------------------------------------
# 6. Detect Res2 cluster column
# ------------------------------------------------------------------------------

CLUSTER_COL <- first_existing(
    meta_cols,
    c(
        "cluster_res2",
        "mphi_rpca_res_2",
        "mphi_rpca_res_2.0",
        "integratedRPCA_snn_res.2",
        "integratedRPCA_snn_res.2.0",
        "seurat_clusters"
    )
)

if (is.na(
    CLUSTER_COL
)) {

    stop(
        "Could not identify Res2 cluster column."
    )
}

msg(
    "Res2 cluster column: ",
    CLUSTER_COL
)

# ------------------------------------------------------------------------------
# 7. Pre-update audit
# ------------------------------------------------------------------------------

old_annotation <- as.character(
    mphi@meta.data[[SOURCE_CLASS_COL]]
)

old_counts <- as.data.frame(
    table(
        old_annotation,
        useNA = "ifany"
    )
)

colnames(
    old_counts
) <- c(
    "annotation",
    "n_cells"
)

write.csv(
    old_counts,
    file.path(
        TAB_DIR,
        "01_annotation_counts_BEFORE_v4.14.5.csv"
    ),
    row.names = FALSE
)

n_old_fib <- sum(
    old_annotation ==
        OLD_CLASS_NAME,
    na.rm = TRUE
)

msg(
    "Cells previously labeled ",
    OLD_CLASS_NAME,
    ": ",
    n_old_fib
)

if (n_old_fib == 0L) {

    stop(
        "No cells with historical label '",
        OLD_CLASS_NAME,
        "' were found."
    )
}

# ------------------------------------------------------------------------------
# 8. Create FINAL annotation
# ------------------------------------------------------------------------------

final_annotation <- old_annotation

final_annotation[
    final_annotation ==
        OLD_CLASS_NAME
] <- NEW_CLASS_NAME

mphi@meta.data[[FINAL_CLASS_COL]] <- factor(
    final_annotation,
    levels = FINAL_CLASS_ORDER
)

# Also add a character copy for safer export/interoperability.
mphi$macrophage_class_Res2_FINAL_v4145_char <- as.character(
    mphi@meta.data[[FINAL_CLASS_COL]]
)

# ------------------------------------------------------------------------------
# 9. Verify that ONLY the label changed
# ------------------------------------------------------------------------------

n_changed <- sum(
    old_annotation !=
        as.character(
            mphi@meta.data[[FINAL_CLASS_COL]]
        ),
    na.rm = TRUE
)

if (n_changed != n_old_fib) {

    stop(
        "Unexpected number of changed cells.\n",
        "Expected: ",
        n_old_fib,
        "\nObserved: ",
        n_changed
    )
}

unchanged_nonfib <- all(
    old_annotation[
        old_annotation !=
            OLD_CLASS_NAME
    ] ==
        as.character(
            mphi@meta.data[[FINAL_CLASS_COL]]
        )[
            old_annotation !=
                OLD_CLASS_NAME
        ]
)

if (!unchanged_nonfib) {

    stop(
        "At least one non-Fibrogenic cell changed annotation unexpectedly."
    )
}

msg(
    "Verified: only historical Fibrogenic-Mphi cells were renamed."
)

# ------------------------------------------------------------------------------
# 10. Res2 cluster audit
# ------------------------------------------------------------------------------

cluster_annotation_audit <- mphi@meta.data %>%
    rownames_to_column(
        "cell"
    ) %>%
    transmute(
        cell = cell,
        Res2_cluster = as.character(
            .data[[CLUSTER_COL]]
        ),
        historical_class = as.character(
            .data[[SOURCE_CLASS_COL]]
        ),
        final_class = as.character(
            .data[[FINAL_CLASS_COL]]
        )
    ) %>%
    count(
        Res2_cluster,
        historical_class,
        final_class,
        name = "n_cells"
    ) %>%
    arrange(
        suppressWarnings(
            as.numeric(
                Res2_cluster
            )
        )
    )

write.csv(
    cluster_annotation_audit,
    file.path(
        TAB_DIR,
        "02_Res2_cluster_annotation_mapping_v4.14.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. Final annotation counts
# ------------------------------------------------------------------------------

final_counts <- mphi@meta.data %>%
    rownames_to_column(
        "cell"
    ) %>%
    count(
        .data[[FINAL_CLASS_COL]],
        name = "n_cells"
    )

colnames(
    final_counts
)[1] <- "final_annotation"

write.csv(
    final_counts,
    file.path(
        TAB_DIR,
        "03_annotation_counts_FINAL_v4.14.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. Explicit old -> new mapping table
# ------------------------------------------------------------------------------

annotation_mapping <- tibble(
    historical_annotation = c(
        "Inflammatory-Mphi",
        "Anti-inflammatory-Mphi",
        "Fibrogenic-Mphi",
        "Repair/Resolution-Mphi",
        "Lipid-associated/TREM2-Mphi",
        "Other"
    ),

    final_annotation = c(
        "Inflammatory-Mphi",
        "Anti-inflammatory-Mphi",
        "ECM-associated inflammatory-Mphi",
        "Repair/Resolution-Mphi",
        "Lipid-associated/TREM2-Mphi",
        "Other"
    ),

    action = c(
        "unchanged",
        "unchanged",
        "renamed after v4.14.4 focused audit",
        "unchanged",
        "unchanged",
        "unchanged"
    )
)

write.csv(
    annotation_mapping,
    file.path(
        TAB_DIR,
        "04_historical_to_FINAL_annotation_mapping_v4.14.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 13. Final framework audit
# ------------------------------------------------------------------------------

unexpected_final <- setdiff(
    unique(
        as.character(
            mphi@meta.data[[FINAL_CLASS_COL]]
        )
    ),
    FINAL_CLASS_ORDER
)

unexpected_final <- unexpected_final[
    !is.na(
        unexpected_final
    )
]

if (length(
    unexpected_final
) > 0L) {

    stop(
        "Unexpected FINAL annotation(s): ",
        paste(
            unexpected_final,
            collapse = ", "
        )
    )
}

if (any(
    is.na(
        mphi@meta.data[[FINAL_CLASS_COL]]
    )
)) {

    stop(
        "NA values were generated in FINAL annotation."
    )
}

# ------------------------------------------------------------------------------
# 14. Preserve analysis provenance
# ------------------------------------------------------------------------------

mphi$annotation_version_FINAL <- "v4.14.5"

mphi$annotation_parent_resolution_FINAL <- "Res2.0"

mphi$annotation_policy_FINAL <- paste0(
    "v4.8.4 five-class framework retained; ",
    "Fibrogenic-Mphi renamed to ECM-associated inflammatory-Mphi ",
    "after v4.14.4 focused audit; cluster membership unchanged."
)

# ------------------------------------------------------------------------------
# 15. Save FINAL Clean-B RDS
# ------------------------------------------------------------------------------

msg(
    "Saving FINAL Clean-B RDS..."
)

saveRDS(
    mphi,
    OUTPUT_RDS
)

if (!file.exists(
    OUTPUT_RDS
)) {

    stop(
        "FINAL RDS was not created."
    )
}

# ------------------------------------------------------------------------------
# 16. Reload validation
# ------------------------------------------------------------------------------

msg(
    "Reloading FINAL RDS for validation..."
)

mphi_check <- readRDS(
    OUTPUT_RDS
)

if (!FINAL_CLASS_COL %in%
    colnames(
        mphi_check@meta.data
    )) {

    stop(
        "FINAL class column missing after reload."
    )
}

if (!identical(
    colnames(
        mphi
    ),
    colnames(
        mphi_check
    )
)) {

    stop(
        "Cell IDs changed after saving/reloading FINAL RDS."
    )
}

if (!identical(
    rownames(
        mphi
    ),
    rownames(
        mphi_check
    )
)) {

    stop(
        "Gene IDs changed after saving/reloading FINAL RDS."
    )
}

if (!identical(
    as.character(
        mphi@meta.data[[FINAL_CLASS_COL]]
    ),
    as.character(
        mphi_check@meta.data[[FINAL_CLASS_COL]]
    )
)) {

    stop(
        "FINAL annotation changed after save/reload."
    )
}

msg(
    "Reload validation passed."
)

# ------------------------------------------------------------------------------
# 17. Final audit summary
# ------------------------------------------------------------------------------

final_summary <- tibble(
    item = c(
        "Input RDS",
        "Output RDS",
        "Source annotation column",
        "Final annotation column",
        "Parent resolution",
        "Cells renamed",
        "Old name",
        "New name",
        "Cluster membership changed",
        "Cell IDs preserved",
        "Gene IDs preserved"
    ),

    value = c(
        INPUT_RDS,
        OUTPUT_RDS,
        SOURCE_CLASS_COL,
        FINAL_CLASS_COL,
        "Res2.0",
        as.character(
            n_changed
        ),
        OLD_CLASS_NAME,
        NEW_CLASS_NAME,
        "No",
        "Yes",
        "Yes"
    )
)

write.csv(
    final_summary,
    file.path(
        TAB_DIR,
        "05_FINAL_annotation_update_summary_v4.14.5.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 18. README
# ------------------------------------------------------------------------------

readme <- c(

    "Mouse MASH MΦ Clean-B FINAL annotation v4.14.5",
    "",

    paste0(
        "Input RDS: ",
        INPUT_RDS
    ),

    paste0(
        "Output RDS: ",
        OUTPUT_RDS
    ),

    "",

    "Final parent clustering:",
    "  Res2.0",
    "",

    "Final parent MΦ annotation:",
    "  1. Inflammatory-MΦ",
    "  2. Anti-inflammatory-MΦ",
    "  3. ECM-associated inflammatory-MΦ",
    "  4. Repair/Resolution-MΦ",
    "  5. Lipid-associated/TREM2-MΦ",
    "  + Other",
    "",

    "Annotation change:",
    "  Fibrogenic-MΦ",
    "    ->",
    "  ECM-associated inflammatory-MΦ",
    "",

    "Rationale:",
    "  v4.14.4 focused audit showed that Res2 clusters 3/21",
    "  retain ECM/fibrosis-associated genes including Thbs1, Fn1,",
    "  Col1a1/Col1a2/Col3a1 and Tgfb1, while also showing",
    "  inflammatory features such as Il1b, Stat1 and Cxcl10.",
    "  The population is therefore better represented as an",
    "  ECM-associated inflammatory macrophage state than as a",
    "  purely Fibrogenic-MΦ population.",
    "",

    "Important:",
    "  No cell was reassigned to another Res2 cluster.",
    "  No parent MΦ class membership was changed except the class name.",
    "  Historical annotation columns are preserved for backward compatibility.",
    "",

    paste0(
        "Final analysis column: ",
        FINAL_CLASS_COL
    ),
    "",

    "Use this FINAL annotation column for all subsequent analyses."
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_FINAL_v4.14.5.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.14.5.txt"
    )
)

# ------------------------------------------------------------------------------
# 19. Final messages
# ------------------------------------------------------------------------------

msg(
    "DONE."
)

msg(
    "FINAL RDS: ",
    OUTPUT_RDS
)

msg(
    "FINAL annotation column: ",
    FINAL_CLASS_COL
)

msg(
    "Renamed cells: ",
    n_changed
)

print(
    annotation_mapping
)

print(
    final_counts
)
