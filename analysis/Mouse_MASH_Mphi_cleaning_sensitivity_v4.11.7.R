#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ cleaning sensitivity analysis
# v4.11.7
#
# PURPOSE
#   Compare macrophage conclusions across four contamination-cleaning levels:
#
#     Original : no cells removed
#     Clean-A  : remove Res2 cluster 24 only
#     Clean-B  : remove Res2 clusters 24 + 27
#     Clean-C  : Clean-B + strict cell-level contamination candidates
#
#   The goal is NOT to choose the most aggressive cleaning condition.
#   The goal is to test whether major biological conclusions are robust.
#
# PRIMARY ROBUSTNESS QUESTIONS
#   1) Does Anti-inflammatory-MΦ increase from Sham -> Tx remain?
#   2) Does M2/(M1+M2) increase from Sham -> Tx remain?
#   3) Does log2(M2/M1) increase from Sham -> Tx remain?
#   4) How much does Fibrogenic-MΦ change when cluster 27 is removed?
#   5) Are major functional-module conclusions retained?
#   6) Do pseudobulk top genes become cleaner after contamination removal?
#
# IMPORTANT
#   - Res2.0 clustering is NOT recalculated.
#   - v4.8.4 macrophage annotation is NOT recalculated.
#   - fixed mphi.umap.rpca coordinates are retained.
#   - raw-count pseudobulk uses explicit Matrix::rowSums(counts[, cells]).
#   - STD vs CDAHFD is descriptive only if n=1 per condition.
#
# CLEAN-C STRICT CELL RULE
#   Outside clusters 24 and 27, remove a cell only if:
#
#       A) >= 2 lineage candidate FLAGS
#
#          OR
#
#       B) B-cell candidate AND low MΦ identity
#
#          OR
#
#       C) T-cell candidate AND low MΦ identity
#
#   NK-only and neutrophil-only cells are NOT automatically removed because
#   NK/inflammatory/neutrophil-related genes can overlap with activated myeloid
#   biology. This makes Clean-C intentionally conservative.
#
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4117)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

PROJECT_DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
    PROJECT_DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_contamination_audit_v4.11.6",
    "Mouse_Mphi_Res2_contamination_audit_annotated_v4.11.6.rds"
)

OUTPUT_DIR <- file.path(
    PROJECT_DATA_ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_cleaning_sensitivity_v4.11.7"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
PB_DIR  <- file.path(OUTPUT_DIR, "Pseudobulk")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PB_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"

CONDITION_ORDER <- c(
    "STD",
    "CDAHFD",
    "Sham",
    "Tx"
)

ANALYSIS_ORDER <- c(
    "Original",
    "Clean-A",
    "Clean-B",
    "Clean-C"
)

MAJOR_CLASS_ORDER <- c(
    "Inflammatory-Mphi",
    "Anti-inflammatory-Mphi",
    "Fibrogenic-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi"
)

DISPLAY_CLASS <- c(
    "Inflammatory-Mphi"           = "Inflammatory-MΦ",
    "Anti-inflammatory-Mphi"      = "Anti-inflammatory-MΦ",
    "Fibrogenic-Mphi"             = "Fibrogenic-MΦ",
    "Repair/Resolution-Mphi"      = "Repair/Resolution-MΦ",
    "Lipid-associated/TREM2-Mphi" = "Lipid-associated/TREM2-MΦ",
    "Other"                       = "Other"
)

# ------------------------------------------------------------------------------
# 1. Functional marker sets
# ------------------------------------------------------------------------------

FUNCTIONAL_PROGRAMS <- list(

    Inflammatory = c(
        "Il1b","Tnf","Ccl2","Ccl3","Ccl4","Cxcl10",
        "Nos2","Cd80","Cd86","Stat1"
    ),

    Anti_inflammatory = c(
        "Mrc1","Cd163","Il1rn","Retnla","Chil3","Arg1",
        "Mertk","Igf1","Hmox1","Klf4","Maf"
    ),

    Fibrogenic = c(
        "Spp1","Tgfb1","Pdgfb","Thbs1","Lgals3",
        "Gpnmb","Mmp12","Mmp14","Ctsb"
    ),

    Repair_Resolution = c(
        "Mertk","Axl","Mfge8","Gas6","Igf1","Hmox1",
        "Mmp12","Mmp13","Mmp14","Plau"
    ),

    Lipid_TREM2 = c(
        "Trem2","Gpnmb","Cd9","Lpl","Apoe",
        "Fabp5","Abca1","Plin2","Ctsd"
    ),

    IL10_response = c(
        "Il10ra","Il10rb","Jak1","Tyk2","Stat3",
        "Socs3","Bcl3","Il1rn"
    )
)

QC_GENES <- c(
    "Adgre1","Lyz2","Csf1r","Il1b","Tnf",
    "Mrc1","Cd163","Spp1","Trem2","Apoe"
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
    if (length(hit) == 0L) {
        return(NA_character_)
    }
    hit[[1]]
}

canonical_condition <- function(sample) {
    smp <- as.character(sample)
    out <- rep(NA_character_, length(smp))

    out[grepl("^STD", smp, ignore.case = TRUE)] <- "STD"
    out[grepl("CDAHFD|CDHFD", smp, ignore.case = TRUE)] <- "CDAHFD"
    out[grepl("^Sham", smp, ignore.case = TRUE)] <- "Sham"
    out[grepl("^Tx", smp, ignore.case = TRUE)] <- "Tx"

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

safe_pct <- function(n, d) {
    ifelse(
        d > 0,
        100 * n / d,
        NA_real_
    )
}

safe_log2_ratio <- function(a, b, pseudo = 0.5) {
    log2(
        (a + pseudo) /
        (b + pseudo)
    )
}

safe_filename <- function(x) {
    gsub(
        "[^A-Za-z0-9]+",
        "_",
        x
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

mean_expression_score <- function(
    normalized_matrix,
    genes
) {
    genes <- intersect(
        genes,
        rownames(normalized_matrix)
    )

    if (length(genes) == 0L) {
        return(
            rep(
                NA_real_,
                ncol(normalized_matrix)
            )
        )
    }

    as.numeric(
        Matrix::colMeans(
            normalized_matrix[
                genes,
                ,
                drop = FALSE
            ]
        )
    )
}

calc_cpm <- function(pb_counts) {
    libs <- colSums(pb_counts)

    if (any(libs <= 0)) {
        stop(
            "At least one pseudobulk library has zero total counts."
        )
    }

    sweep(
        pb_counts,
        2,
        libs / 1e6,
        "/"
    )
}

# ------------------------------------------------------------------------------
# 4. Load v4.11.6 audit object
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
    stop(
        "Input v4.11.6 audit RDS not found:\n",
        INPUT_RDS
    )
}

msg("Loading: ", INPUT_RDS)

mphi <- readRDS(
    INPUT_RDS
)

if (!inherits(mphi, "Seurat")) {
    stop("Input is not a Seurat object.")
}

DefaultAssay(mphi) <- ASSAY_USE

# ------------------------------------------------------------------------------
# 5. Detect metadata
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

if (is.na(SAMPLE_COL)) {
    stop("Sample column not found.")
}

if (is.na(CLASS_COL)) {
    stop("Fixed v4.8.4 macrophage class column not found.")
}

if (is.na(CLUSTER_COL)) {
    stop("Res2.0 cluster column not found.")
}

required_qc_cols <- c(
    "QC_B_candidate",
    "QC_T_candidate",
    "QC_NK_candidate",
    "QC_Neutrophil_candidate",
    "QC_any_nonMphi_candidate",
    "QC_n_lineage_flags_v4116",
    "QC_Mphi_identity_score_v4116"
)

missing_qc <- setdiff(
    required_qc_cols,
    meta_cols
)

if (length(missing_qc) > 0L) {
    stop(
        "Required v4.11.6 QC column(s) missing: ",
        paste(missing_qc, collapse = ", ")
    )
}

msg("SAMPLE_COL = ", SAMPLE_COL)
msg("CLASS_COL = ", CLASS_COL)
msg("CLUSTER_COL = ", CLUSTER_COL)

# ------------------------------------------------------------------------------
# 6. Validate RNA matrices
# ------------------------------------------------------------------------------

counts <- get_layer_safe(
    mphi,
    ASSAY_USE,
    "counts"
)

data_norm <- get_layer_safe(
    mphi,
    ASSAY_USE,
    "data"
)

if (is.null(counts)) {
    stop("RNA counts layer not found.")
}

if (is.null(data_norm)) {
    stop("RNA normalized data layer not found.")
}

if (!identical(
    colnames(counts),
    colnames(mphi)
)) {
    stop(
        "RNA counts cell order does not match Seurat object."
    )
}

if (!identical(
    colnames(data_norm),
    colnames(mphi)
)) {
    stop(
        "RNA normalized data cell order does not match Seurat object."
    )
}

# ------------------------------------------------------------------------------
# 7. Canonical metadata and strict Clean-C rule
# ------------------------------------------------------------------------------

meta <- mphi@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = as.character(.data[[SAMPLE_COL]]),
        condition = canonical_condition(
            as.character(.data[[SAMPLE_COL]])
        ),
        macrophage_class = as.character(.data[[CLASS_COL]]),
        cluster = as.character(.data[[CLUSTER_COL]]),

        B_candidate = QC_B_candidate,
        T_candidate = QC_T_candidate,
        NK_candidate = QC_NK_candidate,
        Neutrophil_candidate = QC_Neutrophil_candidate,
        any_nonMphi_candidate = QC_any_nonMphi_candidate,

        n_lineage_flags = QC_n_lineage_flags_v4116,
        Mphi_identity_score =
            QC_Mphi_identity_score_v4116
    )

if (anyNA(meta$condition)) {
    stop(
        "Could not assign condition for sample(s): ",
        paste(
            unique(
                meta$sample[
                    is.na(meta$condition)
                ]
            ),
            collapse = ", "
        )
    )
}

# Low-MΦ identity cutoff is defined from the full original dataset.
MPHI_IDENTITY_LOW_CUTOFF <- quantile(
    meta$Mphi_identity_score,
    probs = 0.10,
    na.rm = TRUE
)

meta <- meta %>%
    mutate(
        Mphi_identity_low =
            Mphi_identity_score <=
            MPHI_IDENTITY_LOW_CUTOFF,

        strict_cell_candidate = (
            n_lineage_flags >= 2L
        ) | (
            B_candidate &
                Mphi_identity_low
        ) | (
            T_candidate &
                Mphi_identity_low
        )
    )

# ------------------------------------------------------------------------------
# 8. Define the four sensitivity-analysis masks
# ------------------------------------------------------------------------------

meta <- meta %>%
    mutate(
        keep_Original = TRUE,

        keep_Clean_A =
            cluster != "24",

        keep_Clean_B =
            !(cluster %in% c("24", "27")),

        keep_Clean_C =
            !(cluster %in% c("24", "27")) &
            !strict_cell_candidate
    )

analysis_keep_cols <- c(
    "Original" = "keep_Original",
    "Clean-A"  = "keep_Clean_A",
    "Clean-B"  = "keep_Clean_B",
    "Clean-C"  = "keep_Clean_C"
)

# ------------------------------------------------------------------------------
# 9. Cleaning rule audit
# ------------------------------------------------------------------------------

cleaning_rule_audit <- tibble(
    analysis = ANALYSIS_ORDER,
    rule = c(
        "No cells removed",
        "Remove cluster 24",
        "Remove clusters 24 and 27",
        paste0(
            "Remove clusters 24 and 27 + strict cell candidates: ",
            "multi-lineage OR B/T candidate with low MΦ identity"
        )
    )
)

write.csv(
    cleaning_rule_audit,
    file.path(
        TAB_DIR,
        "00_cleaning_rules_v4.11.7.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Cell retention summary
# ------------------------------------------------------------------------------

retention_summary <- bind_rows(
    lapply(
        ANALYSIS_ORDER,
        function(analysis_name) {

            keep_col <- analysis_keep_cols[[analysis_name]]
            keep_vec <- meta[[keep_col]]

            tibble(
                analysis = analysis_name,
                total_original_cells = nrow(meta),
                retained_cells = sum(keep_vec),
                removed_cells = sum(!keep_vec),
                retained_percent =
                    100 *
                    mean(keep_vec),
                removed_percent =
                    100 *
                    mean(!keep_vec)
            )
        }
    )
)

write.csv(
    retention_summary,
    file.path(
        TAB_DIR,
        "01_cell_retention_summary_v4.11.7.csv"
    ),
    row.names = FALSE
)

print(retention_summary)

# ------------------------------------------------------------------------------
# 11. Removed cells by sample / subtype / cluster
# ------------------------------------------------------------------------------

removal_detail <- bind_rows(
    lapply(
        ANALYSIS_ORDER[-1],
        function(analysis_name) {

            keep_col <- analysis_keep_cols[[analysis_name]]

            meta %>%
                filter(
                    !.data[[keep_col]]
                ) %>%
                count(
                    sample,
                    condition,
                    macrophage_class,
                    cluster,
                    name = "n_removed"
                ) %>%
                mutate(
                    analysis = analysis_name,
                    .before = sample
                )
        }
    )
)

write.csv(
    removal_detail,
    file.path(
        TAB_DIR,
        "02_removed_cells_by_sample_subtype_cluster_v4.11.7.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. Build long retained-cell table
# ------------------------------------------------------------------------------

retained_long <- bind_rows(
    lapply(
        ANALYSIS_ORDER,
        function(analysis_name) {

            keep_col <- analysis_keep_cols[[analysis_name]]

            meta %>%
                filter(
                    .data[[keep_col]]
                ) %>%
                mutate(
                    analysis = analysis_name,
                    .before = cell
                )
        }
    )
)

retained_long$analysis <- factor(
    retained_long$analysis,
    levels = ANALYSIS_ORDER
)

# ------------------------------------------------------------------------------
# 13. Sample-level subtype counts and fractions
# ------------------------------------------------------------------------------

sample_total <- retained_long %>%
    group_by(
        analysis,
        sample,
        condition
    ) %>%
    summarise(
        total_Mphi = n(),
        .groups = "drop"
    )

subtype_counts <- retained_long %>%
    count(
        analysis,
        sample,
        condition,
        macrophage_class,
        name = "n_cells"
    ) %>%
    left_join(
        sample_total,
        by = c(
            "analysis",
            "sample",
            "condition"
        )
    ) %>%
    mutate(
        fraction_all_Mphi =
            n_cells /
            total_Mphi,
        percent_all_Mphi =
            100 *
            fraction_all_Mphi
    )

write.csv(
    subtype_counts,
    file.path(
        TAB_DIR,
        "03_subtype_counts_and_fractions_by_sample_v4.11.7.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 14. M1 / M2 balance
# ------------------------------------------------------------------------------

m1m2_table <- subtype_counts %>%
    filter(
        macrophage_class %in%
            c(
                "Inflammatory-Mphi",
                "Anti-inflammatory-Mphi"
            )
    ) %>%
    select(
        analysis,
        sample,
        condition,
        macrophage_class,
        n_cells
    ) %>%
    pivot_wider(
        names_from = macrophage_class,
        values_from = n_cells,
        values_fill = 0
    ) %>%
    rename(
        M1 = `Inflammatory-Mphi`,
        M2 = `Anti-inflammatory-Mphi`
    ) %>%
    mutate(
        M2_fraction_M1M2 =
            M2 /
            pmax(
                M1 + M2,
                1
            ),

        M2_percent_M1M2 =
            100 *
            M2_fraction_M1M2,

        M2_M1_ratio =
            (M2 + 0.5) /
            (M1 + 0.5),

        log2_M2_M1 =
            log2(
                M2_M1_ratio
            )
    )

write.csv(
    m1m2_table,
    file.path(
        TAB_DIR,
        "04_M1_M2_balance_by_sample_v4.11.7.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 15. Condition-level subtype summary
# ------------------------------------------------------------------------------

condition_subtype_summary <- subtype_counts %>%
    filter(
        macrophage_class %in%
            MAJOR_CLASS_ORDER
    ) %>%
    group_by(
        analysis,
        condition,
        macrophage_class
    ) %>%
    summarise(
        n_biological_samples =
            n_distinct(sample),

        mean_n_cells =
            mean(n_cells),

        mean_percent_all_Mphi =
            mean(percent_all_Mphi),

        min_percent_all_Mphi =
            min(percent_all_Mphi),

        max_percent_all_Mphi =
            max(percent_all_Mphi),

        .groups = "drop"
    )

write.csv(
    condition_subtype_summary,
    file.path(
        TAB_DIR,
        "05_condition_subtype_summary_v4.11.7.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 16. Figure: subtype abundance sensitivity
# ------------------------------------------------------------------------------

p_subtype_sensitivity <- ggplot(
    subtype_counts %>%
        filter(
            macrophage_class %in%
                MAJOR_CLASS_ORDER
        ),
    aes(
        x = condition,
        y = percent_all_Mphi,
        group = sample
    )
) +
    geom_point(
        size = 2.1
    ) +
    facet_grid(
        macrophage_class ~ analysis,
        scales = "free_y",
        labeller = labeller(
            macrophage_class =
                as_labeller(
                    DISPLAY_CLASS
                )
        )
    ) +
    labs(
        title =
            "MΦ subtype abundance sensitivity to contamination cleaning",
        subtitle =
            "Each point = biological sample | y = subtype fraction of all retained MΦ",
        x = NULL,
        y = "% of retained MΦ"
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
                face = "bold",
                size = 8
            ),
        axis.text.x =
            element_text(
                angle = 35,
                hjust = 1
            )
    )

save_pdf(
    "01_subtype_abundance_cleaning_sensitivity_v4.11.7.pdf",
    p_subtype_sensitivity,
    15,
    15
)

# ------------------------------------------------------------------------------
# 17. Figure: Anti-inflammatory MΦ abundance
# ------------------------------------------------------------------------------

anti_df <- subtype_counts %>%
    filter(
        macrophage_class ==
            "Anti-inflammatory-Mphi"
    )

p_anti <- ggplot(
    anti_df,
    aes(
        x = condition,
        y = percent_all_Mphi,
        group = analysis
    )
) +
    geom_point(
        aes(
            shape = sample
        ),
        size = 2.5
    ) +
    facet_wrap(
        ~ analysis,
        nrow = 1
    ) +
    labs(
        title =
            "Anti-inflammatory-MΦ abundance across cleaning conditions",
        subtitle =
            "Primary robustness test for Sham → Tx M2 recovery",
        x = NULL,
        y = "Anti-inflammatory-MΦ / all retained MΦ (%)",
        shape = "Sample"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            ),
        axis.text.x =
            element_text(
                angle = 35,
                hjust = 1
            )
    )

save_pdf(
    "02_Anti_inflammatory_Mphi_abundance_sensitivity_v4.11.7.pdf",
    p_anti,
    13,
    5
)

# ------------------------------------------------------------------------------
# 18. Figure: M2/(M1+M2)
# ------------------------------------------------------------------------------

p_m2_fraction <- ggplot(
    m1m2_table,
    aes(
        x = condition,
        y = M2_percent_M1M2
    )
) +
    geom_point(
        aes(
            shape = sample
        ),
        size = 2.5
    ) +
    facet_wrap(
        ~ analysis,
        nrow = 1
    ) +
    labs(
        title =
            "M2 / (M1 + M2) sensitivity analysis",
        x = NULL,
        y = "M2 / (M1 + M2) (%)",
        shape = "Sample"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            ),
        axis.text.x =
            element_text(
                angle = 35,
                hjust = 1
            )
    )

save_pdf(
    "03_M2_fraction_M1M2_sensitivity_v4.11.7.pdf",
    p_m2_fraction,
    13,
    5
)

# ------------------------------------------------------------------------------
# 19. Figure: log2(M2/M1)
# ------------------------------------------------------------------------------

p_m2m1 <- ggplot(
    m1m2_table,
    aes(
        x = condition,
        y = log2_M2_M1
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = 2,
        linewidth = 0.4
    ) +
    geom_point(
        aes(
            shape = sample
        ),
        size = 2.5
    ) +
    facet_wrap(
        ~ analysis,
        nrow = 1
    ) +
    labs(
        title =
            "log2(M2/M1) sensitivity analysis",
        x = NULL,
        y = "log2[(M2 + 0.5)/(M1 + 0.5)]",
        shape = "Sample"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            ),
        axis.text.x =
            element_text(
                angle = 35,
                hjust = 1
            )
    )

save_pdf(
    "04_log2_M2_M1_sensitivity_v4.11.7.pdf",
    p_m2m1,
    13,
    5
)

# ------------------------------------------------------------------------------
# 20. Figure: Fibrogenic-MΦ sensitivity
# ------------------------------------------------------------------------------

fib_df <- subtype_counts %>%
    filter(
        macrophage_class ==
            "Fibrogenic-Mphi"
    )

p_fib <- ggplot(
    fib_df,
    aes(
        x = condition,
        y = percent_all_Mphi
    )
) +
    geom_point(
        aes(
            shape = sample
        ),
        size = 2.5
    ) +
    facet_wrap(
        ~ analysis,
        nrow = 1
    ) +
    labs(
        title =
            "Fibrogenic-MΦ abundance after progressive contamination cleaning",
        subtitle =
            "Particularly important because Res2 cluster 27 was classified as Fibrogenic-MΦ",
        x = NULL,
        y = "Fibrogenic-MΦ / all retained MΦ (%)",
        shape = "Sample"
    ) +
    theme_classic(
        base_size = 10
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            ),
        axis.text.x =
            element_text(
                angle = 35,
                hjust = 1
            )
    )

save_pdf(
    "05_Fibrogenic_Mphi_abundance_sensitivity_v4.11.7.pdf",
    p_fib,
    13,
    5
)

# ------------------------------------------------------------------------------
# 21. Functional program scores from normalized RNA
# ------------------------------------------------------------------------------

PROGRAM_USE <- lapply(
    FUNCTIONAL_PROGRAMS,
    intersect,
    y = rownames(data_norm)
)

program_detection <- bind_rows(
    lapply(
        names(FUNCTIONAL_PROGRAMS),
        function(program_name) {

            tibble(
                program = program_name,
                detected_genes =
                    paste(
                        PROGRAM_USE[[program_name]],
                        collapse = ";"
                    ),
                n_detected =
                    length(
                        PROGRAM_USE[[program_name]]
                    )
            )
        }
    )
)

write.csv(
    program_detection,
    file.path(
        TAB_DIR,
        "06_functional_program_gene_detection_v4.11.7.csv"
    ),
    row.names = FALSE
)

program_scores <- tibble(
    cell = colnames(mphi)
)

for (program_name in names(PROGRAM_USE)) {

    score_now <- mean_expression_score(
        data_norm,
        PROGRAM_USE[[program_name]]
    )

    program_scores[[program_name]] <- score_now
}

program_scores <- program_scores %>%
    left_join(
        meta %>%
            select(
                cell,
                sample,
                condition,
                macrophage_class
            ),
        by = "cell"
    )

# ------------------------------------------------------------------------------
# 22. Functional score summary under each cleaning condition
# ------------------------------------------------------------------------------

program_long <- program_scores %>%
    pivot_longer(
        cols = all_of(
            names(PROGRAM_USE)
        ),
        names_to = "program",
        values_to = "score"
    )

functional_sensitivity <- bind_rows(
    lapply(
        ANALYSIS_ORDER,
        function(analysis_name) {

            keep_col <- analysis_keep_cols[[analysis_name]]

            keep_cells <- meta$cell[
                meta[[keep_col]]
            ]

            program_long %>%
                filter(
                    cell %in% keep_cells,
                    macrophage_class %in%
                        MAJOR_CLASS_ORDER
                ) %>%
                group_by(
                    sample,
                    condition,
                    macrophage_class,
                    program
                ) %>%
                summarise(
                    mean_score =
                        mean(
                            score,
                            na.rm = TRUE
                        ),
                    .groups = "drop"
                ) %>%
                mutate(
                    analysis = analysis_name,
                    .before = sample
                )
        }
    )
)

write.csv(
    functional_sensitivity,
    file.path(
        TAB_DIR,
        "07_functional_program_scores_by_sample_subtype_v4.11.7.csv"
    ),
    row.names = FALSE
)

# Focus on Anti-inflammatory MΦ program behavior.
p_functional_anti <- ggplot(
    functional_sensitivity %>%
        filter(
            macrophage_class ==
                "Anti-inflammatory-Mphi"
        ),
    aes(
        x = condition,
        y = mean_score
    )
) +
    geom_point(
        aes(
            shape = sample
        ),
        size = 2.2
    ) +
    facet_grid(
        program ~ analysis,
        scales = "free_y"
    ) +
    labs(
        title =
            "Anti-inflammatory-MΦ functional programs across cleaning conditions",
        x = NULL,
        y = "Mean normalized program expression",
        shape = "Sample"
    ) +
    theme_classic(
        base_size = 8.5
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            ),
        strip.text =
            element_text(
                face = "bold",
                size = 7.5
            ),
        axis.text.x =
            element_text(
                angle = 35,
                hjust = 1
            )
    )

save_pdf(
    "06_Anti_inflammatory_Mphi_functional_program_sensitivity_v4.11.7.pdf",
    p_functional_anti,
    15,
    14
)

# ------------------------------------------------------------------------------
# 23. SAFE pseudobulk builder
# ------------------------------------------------------------------------------

build_safe_pseudobulk <- function(
    analysis_name,
    keep_cells
) {

    meta_keep <- meta %>%
        filter(
            cell %in%
                keep_cells,
            macrophage_class %in%
                MAJOR_CLASS_ORDER
        )

    pb_meta <- meta_keep %>%
        distinct(
            sample,
            condition,
            macrophage_class
        ) %>%
        arrange(
            factor(
                macrophage_class,
                levels =
                    MAJOR_CLASS_ORDER
            ),
            condition,
            sample
        ) %>%
        mutate(
            pb_id =
                paste(
                    analysis_name,
                    macrophage_class,
                    sample,
                    sep = "__"
                )
        )

    pb_list <- vector(
        "list",
        nrow(pb_meta)
    )

    names(pb_list) <- pb_meta$pb_id

    pb_qc <- vector(
        "list",
        nrow(pb_meta)
    )

    for (i in seq_len(nrow(pb_meta))) {

        smp <- pb_meta$sample[[i]]
        subtype <- pb_meta$macrophage_class[[i]]
        pb_id <- pb_meta$pb_id[[i]]

        cells <- meta_keep$cell[
            meta_keep$sample == smp &
                meta_keep$macrophage_class ==
                    subtype
        ]

        cells <- intersect(
            cells,
            colnames(counts)
        )

        if (length(cells) == 0L) {
            next
        }

        pb <- Matrix::rowSums(
            counts[
                ,
                cells,
                drop = FALSE
            ]
        )

        pb_list[[pb_id]] <- pb

        pb_qc[[i]] <- tibble(
            analysis = analysis_name,
            pb_id = pb_id,
            sample = smp,
            condition = as.character(
                pb_meta$condition[[i]]
            ),
            macrophage_class = subtype,
            n_cells = length(cells),
            library_size = sum(pb)
        )
    }

    valid <- !vapply(
        pb_list,
        is.null,
        logical(1)
    )

    pb_list <- pb_list[valid]

    pb_mat <- do.call(
        cbind,
        pb_list
    )

    rownames(pb_mat) <- rownames(counts)
    colnames(pb_mat) <- names(pb_list)

    list(
        counts = pb_mat,
        qc = bind_rows(pb_qc),
        meta = pb_meta %>%
            filter(
                pb_id %in%
                    colnames(pb_mat)
            )
    )
}

# ------------------------------------------------------------------------------
# 24. Build pseudobulk for all four cleaning conditions
# ------------------------------------------------------------------------------

pb_results <- list()

for (analysis_name in ANALYSIS_ORDER) {

    keep_col <- analysis_keep_cols[[analysis_name]]

    keep_cells <- meta$cell[
        meta[[keep_col]]
    ]

    msg(
        "Building safe pseudobulk: ",
        analysis_name
    )

    pb_results[[analysis_name]] <-
        build_safe_pseudobulk(
            analysis_name,
            keep_cells
        )

    write.csv(
        pb_results[[analysis_name]]$qc,
        file.path(
            PB_DIR,
            paste0(
                "PB_QC_",
                safe_filename(analysis_name),
                "_v4.11.7.csv"
            )
        ),
        row.names = FALSE
    )
}

# ------------------------------------------------------------------------------
# 25. Pseudobulk effect-size comparison
# ------------------------------------------------------------------------------

pb_effect_all <- list()

for (analysis_name in ANALYSIS_ORDER) {

    pb <- pb_results[[analysis_name]]

    pb_cpm <- calc_cpm(
        pb$counts
    )

    pb_meta <- pb$meta

    for (subtype in MAJOR_CLASS_ORDER) {

        get_mean_cpm <- function(cond) {

            cols <- pb_meta$pb_id[
                pb_meta$macrophage_class ==
                    subtype &
                    as.character(
                        pb_meta$condition
                    ) ==
                    cond
            ]

            cols <- intersect(
                cols,
                colnames(pb_cpm)
            )

            if (length(cols) == 0L) {
                return(
                    rep(
                        NA_real_,
                        nrow(pb_cpm)
                    )
                )
            }

            if (length(cols) == 1L) {
                return(
                    as.numeric(
                        pb_cpm[, cols]
                    )
                )
            }

            rowMeans(
                pb_cpm[
                    ,
                    cols,
                    drop = FALSE
                ]
            )
        }

        std <- get_mean_cpm("STD")
        cdahfd <- get_mean_cpm("CDAHFD")
        sham <- get_mean_cpm("Sham")
        tx <- get_mean_cpm("Tx")

        pb_effect_all[[paste(
            analysis_name,
            subtype,
            sep = "__"
        )]] <- tibble(
            analysis = analysis_name,
            macrophage_class = subtype,
            gene = rownames(pb_cpm),

            log2FC_CDAHFD_vs_STD =
                log2(
                    (cdahfd + 1) /
                        (std + 1)
                ),

            log2FC_Tx_vs_Sham =
                log2(
                    (tx + 1) /
                        (sham + 1)
                )
        )
    }
}

pb_effect_table <- bind_rows(
    pb_effect_all
)

write.csv(
    pb_effect_table,
    file.path(
        PB_DIR,
        "08_pseudobulk_effect_sizes_all_cleaning_conditions_v4.11.7.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 26. QC-gene pseudobulk robustness
# ------------------------------------------------------------------------------

qc_genes_use <- intersect(
    QC_GENES,
    rownames(counts)
)

qc_pb <- pb_effect_table %>%
    filter(
        gene %in%
            qc_genes_use
    )

write.csv(
    qc_pb,
    file.path(
        PB_DIR,
        "09_QC_gene_pseudobulk_effect_sizes_v4.11.7.csv"
    ),
    row.names = FALSE
)

p_qc_pb <- ggplot(
    qc_pb,
    aes(
        x = gene,
        y = log2FC_Tx_vs_Sham
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = 2,
        linewidth = 0.4
    ) +
    geom_point(
        size = 2
    ) +
    facet_grid(
        macrophage_class ~ analysis,
        scales = "free_y",
        labeller = labeller(
            macrophage_class =
                as_labeller(
                    DISPLAY_CLASS
                )
        )
    ) +
    labs(
        title =
            "Key MΦ-gene pseudobulk robustness: Sham → Tx",
        x = NULL,
        y = "log2FC (Tx / Sham)"
    ) +
    theme_classic(
        base_size = 8.5
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            ),
        axis.text.x =
            element_text(
                angle = 60,
                hjust = 1
            ),
        strip.text =
            element_text(
                face = "bold",
                size = 7.5
            )
    )

save_pdf(
    "07_QC_gene_pseudobulk_Tx_vs_Sham_sensitivity_v4.11.7.pdf",
    p_qc_pb,
    15,
    12
)

# ------------------------------------------------------------------------------
# 27. Key robustness summary table
# ------------------------------------------------------------------------------

anti_summary <- anti_df %>%
    select(
        analysis,
        sample,
        condition,
        anti_percent = percent_all_Mphi
    )

fib_summary <- fib_df %>%
    select(
        analysis,
        sample,
        condition,
        fib_percent = percent_all_Mphi
    )

key_summary <- m1m2_table %>%
    select(
        analysis,
        sample,
        condition,
        M1,
        M2,
        M2_percent_M1M2,
        log2_M2_M1
    ) %>%
    left_join(
        anti_summary,
        by = c(
            "analysis",
            "sample",
            "condition"
        )
    ) %>%
    left_join(
        fib_summary,
        by = c(
            "analysis",
            "sample",
            "condition"
        )
    )

write.csv(
    key_summary,
    file.path(
        TAB_DIR,
        "08_key_robustness_metrics_by_sample_v4.11.7.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 28. Change summary: Sham -> Tx
# ------------------------------------------------------------------------------

tx_summary <- key_summary %>%
    filter(
        condition %in%
            c(
                "Sham",
                "Tx"
            )
    ) %>%
    group_by(
        analysis,
        condition
    ) %>%
    summarise(
        n_samples =
            n_distinct(sample),

        mean_anti_percent =
            mean(
                anti_percent,
                na.rm = TRUE
            ),

        mean_M2_percent_M1M2 =
            mean(
                M2_percent_M1M2,
                na.rm = TRUE
            ),

        mean_log2_M2_M1 =
            mean(
                log2_M2_M1,
                na.rm = TRUE
            ),

        mean_fib_percent =
            mean(
                fib_percent,
                na.rm = TRUE
            ),

        .groups = "drop"
    ) %>%
    pivot_wider(
        names_from =
            condition,
        values_from =
            c(
                mean_anti_percent,
                mean_M2_percent_M1M2,
                mean_log2_M2_M1,
                mean_fib_percent
            )
    ) %>%
    mutate(
        delta_Anti_Mphi =
            mean_anti_percent_Tx -
            mean_anti_percent_Sham,

        delta_M2_fraction =
            mean_M2_percent_M1M2_Tx -
            mean_M2_percent_M1M2_Sham,

        delta_log2_M2_M1 =
            mean_log2_M2_M1_Tx -
            mean_log2_M2_M1_Sham,

        delta_Fibrogenic_Mphi =
            mean_fib_percent_Tx -
            mean_fib_percent_Sham
    )

write.csv(
    tx_summary,
    file.path(
        TAB_DIR,
        "09_Sham_to_Tx_robustness_delta_summary_v4.11.7.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 29. Robustness delta figure
# ------------------------------------------------------------------------------

delta_long <- tx_summary %>%
    select(
        analysis,
        delta_Anti_Mphi,
        delta_M2_fraction,
        delta_log2_M2_M1,
        delta_Fibrogenic_Mphi
    ) %>%
    pivot_longer(
        cols = -analysis,
        names_to = "metric",
        values_to = "delta"
    )

p_delta <- ggplot(
    delta_long,
    aes(
        x = analysis,
        y = delta,
        group = metric
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = 2,
        linewidth = 0.4
    ) +
    geom_line(
        linewidth = 0.65
    ) +
    geom_point(
        size = 2.4
    ) +
    facet_wrap(
        ~ metric,
        scales = "free_y",
        ncol = 2
    ) +
    labs(
        title =
            "Robustness of Sham → Tx conclusions to contamination cleaning",
        subtitle =
            "Positive delta = increase in Tx relative to Sham",
        x = NULL,
        y = "Tx − Sham"
    ) +
    theme_classic(
        base_size = 10
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
    "08_Sham_to_Tx_delta_robustness_summary_v4.11.7.pdf",
    p_delta,
    10,
    8
)

# ------------------------------------------------------------------------------
# 30. Save retained cell lists / RDS objects
# ------------------------------------------------------------------------------

for (analysis_name in ANALYSIS_ORDER) {

    keep_col <- analysis_keep_cols[[analysis_name]]

    keep_cells <- meta$cell[
        meta[[keep_col]]
    ]

    writeLines(
        keep_cells,
        file.path(
            TAB_DIR,
            paste0(
                "retained_cells_",
                safe_filename(analysis_name),
                "_v4.11.7.txt"
            )
        )
    )

    obj_clean <- subset(
        mphi,
        cells = keep_cells
    )

    saveRDS(
        obj_clean,
        file.path(
            RDS_DIR,
            paste0(
                "Mouse_Mphi_Res2_",
                safe_filename(analysis_name),
                "_v4.11.7.rds"
            )
        )
    )

    rm(obj_clean)
    gc()
}

# ------------------------------------------------------------------------------
# 31. Master summary figure
# ------------------------------------------------------------------------------

p_summary <- (
    p_anti /
    p_m2_fraction /
    p_m2m1 /
    p_fib /
    p_delta
) +
    patchwork::plot_layout(
        heights = c(
            0.75,
            0.75,
            0.75,
            0.75,
            1.0
        )
    ) +
    patchwork::plot_annotation(
        title =
            "Mouse MASH MΦ cleaning sensitivity analysis v4.11.7",
        subtitle =
            "Original vs Clean-A vs Clean-B vs Clean-C",
        theme = theme(
            plot.title =
                element_text(
                    face = "bold",
                    size = 18
                )
        )
    )

save_pdf(
    "09_cleaning_sensitivity_master_summary_v4.11.7.pdf",
    p_summary,
    14,
    25
)

# ------------------------------------------------------------------------------
# 32. README
# ------------------------------------------------------------------------------

readme <- c(

    "Mouse MASH MΦ cleaning sensitivity analysis v4.11.7",
    "",

    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",

    "Cleaning conditions:",
    "  Original: no cells removed",
    "  Clean-A : remove Res2 cluster 24",
    "  Clean-B : remove Res2 clusters 24 + 27",
    paste0(
        "  Clean-C : Clean-B + strict cell-level candidates"
    ),
    "",

    "Clean-C strict cell rule:",
    "  remove if >=2 lineage candidate flags",
    "  OR B candidate + low MΦ identity",
    "  OR T candidate + low MΦ identity",
    "  NK-only and neutrophil-only are retained by default",
    "",

    paste0(
        "Low MΦ identity cutoff = 10th percentile of original dataset: ",
        signif(
            MPHI_IDENTITY_LOW_CUTOFF,
            5
        )
    ),
    "",

    "Fixed framework:",
    "  Res2.0 clustering unchanged",
    "  v4.8.4 MΦ annotation unchanged",
    "  no UMAP recalculation",
    "  safe explicit rowSums pseudobulk",
    "",

    "Primary robustness endpoints:",
    "  Anti-inflammatory-MΦ / all MΦ",
    "  M2 / (M1 + M2)",
    "  log2(M2/M1)",
    "  Fibrogenic-MΦ / all MΦ",
    "  Anti-inflammatory-MΦ functional programs",
    "  key-gene pseudobulk Sham vs Tx",
    "",

    "Primary figures:",
    "  01 subtype abundance sensitivity",
    "  02 Anti-inflammatory-MΦ abundance sensitivity",
    "  03 M2/(M1+M2) sensitivity",
    "  04 log2(M2/M1) sensitivity",
    "  05 Fibrogenic-MΦ sensitivity",
    "  06 Anti-inflammatory-MΦ functional-program sensitivity",
    "  07 key-gene pseudobulk sensitivity",
    "  08 Sham->Tx delta robustness",
    "  09 master summary",
    "",

    "Interpretation principle:",
    "  A biological conclusion is considered robust when its direction",
    "  is retained across Original, Clean-A, Clean-B and Clean-C.",
    "",
    "No clustering or manual MΦ-class assignment is changed in v4.11.7."
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.11.7.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.11.7.txt"
    )
)

# ------------------------------------------------------------------------------
# 33. Final
# ------------------------------------------------------------------------------

msg("DONE.")
msg("Output: ", OUTPUT_DIR)

print(retention_summary)
print(tx_summary)
