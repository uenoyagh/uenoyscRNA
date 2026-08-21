#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Formal compositional analysis
# v4.13.1
#
# PRIMARY GOAL
#   Quantify MΦ subtype-composition changes using BIOLOGICAL SAMPLES as units.
#
# PRIMARY DATASET
#   Clean-B from v4.11.7
#
# SENSITIVITY
#   Original / Clean-A / Clean-B / Clean-C
#
# METHODS
#   1) sample-level subtype counts and proportions
#   2) centered log-ratio (CLR) transformation
#   3) Aitchison distance
#   4) composition PCA in CLR space
#   5) subtype-specific Tx-Sham CLR effects
#   6) biologically interpretable balances:
#        M2/M1
#        M2/(M1+M2)
#        Anti-inflammatory vs Fibrogenic
#        Repair/Resolution vs Fibrogenic
#        Anti-inflammatory + Repair vs Inflammatory + Fibrogenic
#   7) Original/Clean-A/B/C robustness
#
# IMPORTANT
#   - Biological sample is the unit of inference.
#   - Sham vs Tx = n=2 vs n=2.
#   - STD vs CDAHFD = n=1 vs n=1 and is descriptive only.
#   - No cell-level p-values are produced.
#
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4131)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_contamination_audit_v4.11.6",
    "Mouse_Mphi_Res2_contamination_audit_annotated_v4.11.6.rds"
)

OUTPUT_DIR <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_formal_compositional_analysis_v4.13.1"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

PSEUDOCOUNT <- 0.5

ANALYSIS_ORDER <- c(
    "Original",
    "Clean-A",
    "Clean-B",
    "Clean-C"
)

CONDITION_ORDER <- c(
    "STD",
    "CDAHFD",
    "Sham",
    "Tx"
)

SUBTYPE_ORDER <- c(
    "Inflammatory-Mphi",
    "Anti-inflammatory-Mphi",
    "Fibrogenic-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi",
    "Other"
)

SUBTYPE_LABELS <- c(
    "Inflammatory-Mphi"           = "Inflammatory-MΦ",
    "Anti-inflammatory-Mphi"      = "Anti-inflammatory-MΦ",
    "Fibrogenic-Mphi"             = "Fibrogenic-MΦ",
    "Repair/Resolution-Mphi"      = "Repair/Resolution-MΦ",
    "Lipid-associated/TREM2-Mphi" = "Lipid-associated/TREM2-MΦ",
    "Other"                       = "Other"
)

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
    "Seurat",
    "SeuratObject",
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
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(patchwork)
    library(scales)
})

# ------------------------------------------------------------------------------
# 2. Helpers
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

clr_transform <- function(counts, pseudocount = 0.5) {
    x <- counts + pseudocount
    p <- x / sum(x)
    lp <- log(p)
    lp - mean(lp)
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

safe_pct <- function(n, d) {
    ifelse(d > 0, 100 * n / d, NA_real_)
}

# ------------------------------------------------------------------------------
# 3. Load audit object
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
    stop(
        "Input v4.11.6 audit RDS not found:\n",
        INPUT_RDS
    )
}

msg("Loading: ", INPUT_RDS)

mphi <- readRDS(INPUT_RDS)

if (!inherits(mphi, "Seurat")) {
    stop("Input is not a Seurat object.")
}

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

if (is.na(SAMPLE_COL)) stop("Sample column not found.")
if (is.na(CLASS_COL)) stop("MΦ class column not found.")
if (is.na(CLUSTER_COL)) stop("Res2 cluster column not found.")

required_qc <- c(
    "QC_B_candidate",
    "QC_T_candidate",
    "QC_NK_candidate",
    "QC_Neutrophil_candidate",
    "QC_n_lineage_flags_v4116",
    "QC_Mphi_identity_score_v4116"
)

missing_qc <- setdiff(
    required_qc,
    meta_cols
)

if (length(missing_qc) > 0L) {
    stop(
        "Required v4.11.6 QC columns missing: ",
        paste(missing_qc, collapse = ", ")
    )
}

# ------------------------------------------------------------------------------
# 4. Canonical metadata + cleaning masks
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
        n_lineage_flags = QC_n_lineage_flags_v4116,
        Mphi_identity_score = QC_Mphi_identity_score_v4116
    )

if (anyNA(meta$condition)) {
    stop("Unresolved biological sample condition.")
}

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
        ),

        keep_Original = TRUE,

        keep_Clean_A =
            cluster != "24",

        keep_Clean_B =
            !(cluster %in% c("24", "27")),

        keep_Clean_C =
            !(cluster %in% c("24", "27")) &
            !strict_cell_candidate
    )

KEEP_COLS <- c(
    "Original" = "keep_Original",
    "Clean-A"  = "keep_Clean_A",
    "Clean-B"  = "keep_Clean_B",
    "Clean-C"  = "keep_Clean_C"
)

# ------------------------------------------------------------------------------
# 5. Build complete sample × subtype tables
#
# v4.13.1 FIX:
#   Do NOT carry a factor-valued condition through complete() + ifelse().
#   condition is reconstructed explicitly from sample AFTER completion.
# ------------------------------------------------------------------------------

composition_tables <- list()

for (analysis_name in ANALYSIS_ORDER) {

    keep_col <- KEEP_COLS[[analysis_name]]

    df_keep <- meta %>%
        filter(
            .data[[keep_col]]
        )

    sample_totals <- df_keep %>%
        count(
            sample,
            name = "total_Mphi"
        )

    comp <- df_keep %>%
        count(
            sample,
            macrophage_class,
            name = "n_cells"
        ) %>%
        complete(
            sample,
            macrophage_class = SUBTYPE_ORDER,
            fill = list(
                n_cells = 0L
            )
        ) %>%
        mutate(
            # Reconstruct from sample as CHARACTER first.
            condition = as.character(
                canonical_condition(sample)
            )
        ) %>%
        left_join(
            sample_totals,
            by = "sample"
        ) %>%
        mutate(
            analysis = analysis_name,
            fraction =
                n_cells /
                total_Mphi,
            percent =
                100 *
                fraction
        ) %>%
        select(
            analysis,
            sample,
            condition,
            macrophage_class,
            n_cells,
            total_Mphi,
            fraction,
            percent
        )

    if (anyNA(comp$condition)) {
        stop(
            "Condition assignment failed in ",
            analysis_name,
            " for sample(s): ",
            paste(
                unique(
                    comp$sample[
                        is.na(comp$condition)
                    ]
                ),
                collapse = ", "
            )
        )
    }

    composition_tables[[analysis_name]] <- comp
}

composition_long <- bind_rows(
    composition_tables
)

composition_long$analysis <- factor(
    composition_long$analysis,
    levels = ANALYSIS_ORDER
)

composition_long$condition <- factor(
    composition_long$condition,
    levels = CONDITION_ORDER
)

composition_long$macrophage_class <- factor(
    composition_long$macrophage_class,
    levels = SUBTYPE_ORDER
)

write.csv(
    composition_long,
    file.path(
        TAB_DIR,
        "01_sample_level_subtype_composition_all_cleaning_conditions_v4.13.1.csv"
    ),
    row.names = FALSE
)

# Audit condition assignment before CLR.
condition_audit <- composition_long %>%
    distinct(
        analysis,
        sample,
        condition
    ) %>%
    arrange(
        analysis,
        factor(
            condition,
            levels = CONDITION_ORDER
        ),
        sample
    )

write.csv(
    condition_audit,
    file.path(
        TAB_DIR,
        "01b_condition_assignment_audit_v4.13.1.csv"
    ),
    row.names = FALSE
)

print(condition_audit)

# ------------------------------------------------------------------------------
# 6. CLR transformation
# ------------------------------------------------------------------------------

clr_rows <- list()

for (analysis_name in ANALYSIS_ORDER) {

    comp_now <- composition_long %>%
        filter(
            as.character(analysis) ==
                analysis_name
        )

    for (smp in unique(comp_now$sample)) {

        sub <- comp_now %>%
            filter(
                sample == smp
            ) %>%
            arrange(
                factor(
                    macrophage_class,
                    levels = SUBTYPE_ORDER
                )
            )

        count_vec <- sub$n_cells

        names(count_vec) <- as.character(
            sub$macrophage_class
        )

        clr_vec <- clr_transform(
            count_vec,
            pseudocount = PSEUDOCOUNT
        )

        condition_now <- unique(
            as.character(
                sub$condition
            )
        )

        condition_now <- condition_now[
            !is.na(condition_now)
        ]

        if (length(condition_now) != 1L) {
            stop(
                "Expected exactly one condition for ",
                analysis_name,
                " / ",
                smp,
                " but found: ",
                paste(
                    condition_now,
                    collapse = ", "
                )
            )
        }

        clr_rows[[paste(
            analysis_name,
            smp,
            sep = "__"
        )]] <- tibble(
            analysis = analysis_name,
            sample = smp,
            condition = condition_now[[1]],
            macrophage_class = names(clr_vec),
            clr = as.numeric(clr_vec)
        )
    }
}

clr_long <- bind_rows(
    clr_rows
)

clr_long$analysis <- factor(
    clr_long$analysis,
    levels = ANALYSIS_ORDER
)

clr_long$condition <- factor(
    clr_long$condition,
    levels = CONDITION_ORDER
)

clr_long$macrophage_class <- factor(
    clr_long$macrophage_class,
    levels = SUBTYPE_ORDER
)

write.csv(
    clr_long,
    file.path(
        TAB_DIR,
        "02_CLR_transformed_subtype_composition_v4.13.1.csv"
    ),
    row.names = FALSE
)

# Critical audit before downstream analysis.
clr_condition_audit <- clr_long %>%
    count(
        analysis,
        condition,
        name = "n_subtype_rows"
    )

write.csv(
    clr_condition_audit,
    file.path(
        TAB_DIR,
        "02b_CLR_condition_audit_v4.13.1.csv"
    ),
    row.names = FALSE
)

print(clr_condition_audit)

# ------------------------------------------------------------------------------
# 7. Primary Clean-B CLR effect
# ------------------------------------------------------------------------------

cleanB_clr <- clr_long %>%
    filter(
        as.character(analysis) ==
            "Clean-B"
    )

# Explicitly verify Sham / Tx are present.
cleanB_conditions_present <- sort(
    unique(
        as.character(
            cleanB_clr$condition
        )
    )
)

cleanB_conditions_present <- cleanB_conditions_present[
    !is.na(cleanB_conditions_present)
]

if (!all(
    c(
        "Sham",
        "Tx"
    ) %in%
        cleanB_conditions_present
)) {

    stop(
        paste0(
            "Clean-B CLR table does not contain both Sham and Tx.\n",
            "Conditions present: ",
            paste(
                cleanB_conditions_present,
                collapse = ", "
            ),
            "\nCheck 01b_condition_assignment_audit_v4.13.1.csv and ",
            "02b_CLR_condition_audit_v4.13.1.csv."
        )
    )
}

cleanB_effect_long <- cleanB_clr %>%
    filter(
        as.character(condition) %in%
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
        n_samples = n_distinct(sample),

        mean_CLR = mean(
            clr,
            na.rm = TRUE
        ),

        min_CLR = min(
            clr,
            na.rm = TRUE
        ),

        max_CLR = max(
            clr,
            na.rm = TRUE
        ),

        .groups = "drop"
    )

write.csv(
    cleanB_effect_long,
    file.path(
        TAB_DIR,
        "03a_CleanB_subtype_CLR_by_condition_v4.13.1.csv"
    ),
    row.names = FALSE
)

cleanB_effect <- cleanB_effect_long %>%
    pivot_wider(
        names_from = condition,
        values_from = c(
            n_samples,
            mean_CLR,
            min_CLR,
            max_CLR
        ),
        names_sep = "_"
    )

required_cleanB_cols <- c(
    "mean_CLR_Sham",
    "mean_CLR_Tx"
)

missing_cleanB_cols <- setdiff(
    required_cleanB_cols,
    colnames(cleanB_effect)
)

if (length(missing_cleanB_cols) > 0L) {

    stop(
        "Expected CLR summary columns were not created: ",
        paste(
            missing_cleanB_cols,
            collapse = ", "
        ),
        "\nActual columns: ",
        paste(
            colnames(cleanB_effect),
            collapse = ", "
        )
    )
}

cleanB_effect <- cleanB_effect %>%
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
    cleanB_effect,
    file.path(
        TAB_DIR,
        "03_CleanB_subtype_CLR_effect_Tx_vs_Sham_v4.13.1.csv"
    ),
    row.names = FALSE
)

print(cleanB_effect)

# ------------------------------------------------------------------------------
# 8. CLR effect sensitivity across Original/A/B/C
# ------------------------------------------------------------------------------

clr_effect_sensitivity <- clr_long %>%
    filter(
        condition %in%
            c(
                "Sham",
                "Tx"
            )
    ) %>%
    group_by(
        analysis,
        macrophage_class,
        condition
    ) %>%
    summarise(
        mean_CLR = mean(clr),
        .groups = "drop"
    ) %>%
    pivot_wider(
        names_from = condition,
        values_from = mean_CLR
    ) %>%
    mutate(
        delta_CLR_Tx_minus_Sham =
            Tx -
            Sham
    )

write.csv(
    clr_effect_sensitivity,
    file.path(
        TAB_DIR,
        "04_CLR_effect_sensitivity_all_cleaning_conditions_v4.13.1.csv"
    ),
    row.names = FALSE
)

p_clr_effect <- ggplot(
    clr_effect_sensitivity,
    aes(
        x = analysis,
        y = delta_CLR_Tx_minus_Sham,
        group = macrophage_class
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
        size = 2.5
    ) +
    facet_wrap(
        ~ macrophage_class,
        ncol = 3,
        labeller = as_labeller(
            SUBTYPE_LABELS
        )
    ) +
    labs(
        title = "MΦ subtype compositional effect: Tx − Sham",
        subtitle = "Centered log-ratio (CLR); biological samples are the units",
        x = NULL,
        y = "ΔCLR (Tx − Sham)"
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
    "01_CLR_effect_sensitivity_v4.13.1.pdf",
    p_clr_effect,
    11,
    8
)

# ------------------------------------------------------------------------------
# 9. Aitchison distance + composition PCA
# ------------------------------------------------------------------------------

cleanB_clr_wide <- cleanB_clr %>%
    select(
        sample,
        condition,
        macrophage_class,
        clr
    ) %>%
    pivot_wider(
        names_from = macrophage_class,
        values_from = clr
    )

clr_matrix <- as.matrix(
    cleanB_clr_wide[
        ,
        setdiff(
            colnames(cleanB_clr_wide),
            c(
                "sample",
                "condition"
            )
        ),
        drop = FALSE
    ]
)

rownames(clr_matrix) <- cleanB_clr_wide$sample

aitchison_distance <- as.matrix(
    dist(
        clr_matrix,
        method = "euclidean"
    )
)

write.csv(
    aitchison_distance,
    file.path(
        TAB_DIR,
        "05_CleanB_Aitchison_distance_matrix_v4.13.1.csv"
    )
)

pca <- prcomp(
    clr_matrix,
    center = TRUE,
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
    left_join(
        cleanB_clr_wide %>%
            select(
                sample,
                condition
            ),
        by = "sample"
    )

var_exp <- 100 *
    (
        pca$sdev^2 /
        sum(
            pca$sdev^2
        )
    )

write.csv(
    pca_df,
    file.path(
        TAB_DIR,
        "06_CleanB_CLR_PCA_coordinates_v4.13.1.csv"
    ),
    row.names = FALSE
)

p_pca <- ggplot(
    pca_df,
    aes(
        x = PC1,
        y = PC2,
        shape = sample
    )
) +
    geom_point(
        size = 4
    ) +
    geom_text(
        aes(
            label = sample
        ),
        nudge_y = 0.12,
        size = 3.2
    ) +
    labs(
        title = "Clean-B MΦ composition PCA",
        subtitle = "PCA in CLR space",
        x = paste0(
            "PC1 (",
            round(
                var_exp[[1]],
                1
            ),
            "%)"
        ),
        y = paste0(
            "PC2 (",
            round(
                var_exp[[2]],
                1
            ),
            "%)"
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

save_pdf(
    "02_CleanB_CLR_PCA_v4.13.1.pdf",
    p_pca,
    7,
    6
)

# ------------------------------------------------------------------------------
# 10. Aitchison distance heatmap
# ------------------------------------------------------------------------------

aitch_df <- as.data.frame(
    as.table(
        aitchison_distance
    )
)

colnames(aitch_df) <- c(
    "Sample1",
    "Sample2",
    "Aitchison_distance"
)

p_aitch <- ggplot(
    aitch_df,
    aes(
        x = Sample1,
        y = Sample2,
        fill = Aitchison_distance
    )
) +
    geom_tile() +
    geom_text(
        aes(
            label = sprintf(
                "%.2f",
                Aitchison_distance
            )
        ),
        size = 3.3
    ) +
    scale_fill_gradient(
        low = "white",
        high = "black"
    ) +
    labs(
        title = "Clean-B Aitchison distance",
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
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        )
    )

save_pdf(
    "03_CleanB_Aitchison_distance_heatmap_v4.13.1.pdf",
    p_aitch,
    7,
    6
)

# ------------------------------------------------------------------------------
# 11. Biologically interpretable balance metrics
# ------------------------------------------------------------------------------

build_balance_table <- function(comp_df) {

    wide <- comp_df %>%
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
        )

    wide %>%
        mutate(
            M2_M1_logratio =
                log(
                    (
                        `Anti-inflammatory-Mphi` +
                        PSEUDOCOUNT
                    ) /
                    (
                        `Inflammatory-Mphi` +
                        PSEUDOCOUNT
                    )
                ),

            M2_fraction_M1M2 =
                (
                    `Anti-inflammatory-Mphi` +
                    PSEUDOCOUNT
                ) /
                (
                    `Anti-inflammatory-Mphi` +
                    `Inflammatory-Mphi` +
                    2 * PSEUDOCOUNT
                ),

            Anti_vs_Fib_logratio =
                log(
                    (
                        `Anti-inflammatory-Mphi` +
                        PSEUDOCOUNT
                    ) /
                    (
                        `Fibrogenic-Mphi` +
                        PSEUDOCOUNT
                    )
                ),

            Repair_vs_Fib_logratio =
                log(
                    (
                        `Repair/Resolution-Mphi` +
                        PSEUDOCOUNT
                    ) /
                    (
                        `Fibrogenic-Mphi` +
                        PSEUDOCOUNT
                    )
                ),

            Protective_vs_InflamFib_logratio =
                log(
                    (
                        `Anti-inflammatory-Mphi` +
                        `Repair/Resolution-Mphi` +
                        PSEUDOCOUNT
                    ) /
                    (
                        `Inflammatory-Mphi` +
                        `Fibrogenic-Mphi` +
                        PSEUDOCOUNT
                    )
                )
        )
}

balance_table <- build_balance_table(
    composition_long
)

write.csv(
    balance_table,
    file.path(
        TAB_DIR,
        "07_sample_level_balance_metrics_all_cleaning_conditions_v4.13.1.csv"
    ),
    row.names = FALSE
)

balance_long <- balance_table %>%
    select(
        analysis,
        sample,
        condition,
        M2_M1_logratio,
        M2_fraction_M1M2,
        Anti_vs_Fib_logratio,
        Repair_vs_Fib_logratio,
        Protective_vs_InflamFib_logratio
    ) %>%
    pivot_longer(
        cols = c(
            M2_M1_logratio,
            M2_fraction_M1M2,
            Anti_vs_Fib_logratio,
            Repair_vs_Fib_logratio,
            Protective_vs_InflamFib_logratio
        ),
        names_to = "metric",
        values_to = "value"
    )

# ------------------------------------------------------------------------------
# 12. Clean-B balance plots
# ------------------------------------------------------------------------------

p_balances <- ggplot(
    balance_long %>%
        filter(
            analysis == "Clean-B"
        ),
    aes(
        x = factor(
            condition,
            levels = CONDITION_ORDER
        ),
        y = value
    )
) +
    geom_point(
        aes(
            shape = sample
        ),
        size = 2.8
    ) +
    facet_wrap(
        ~ metric,
        scales = "free_y",
        ncol = 2
    ) +
    labs(
        title = "Clean-B MΦ compositional balances",
        subtitle = "Each point = biological sample",
        x = NULL,
        y = NULL,
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
            face = "bold",
            size = 8.5
        )
    )

save_pdf(
    "04_CleanB_compositional_balances_v4.13.1.pdf",
    p_balances,
    10,
    9
)

# ------------------------------------------------------------------------------
# 13. Tx-Sham balance effect sensitivity
# ------------------------------------------------------------------------------

balance_effect <- balance_long %>%
    filter(
        condition %in%
            c(
                "Sham",
                "Tx"
            )
    ) %>%
    group_by(
        analysis,
        metric,
        condition
    ) %>%
    summarise(
        mean_value = mean(value),
        .groups = "drop"
    ) %>%
    pivot_wider(
        names_from = condition,
        values_from = mean_value
    ) %>%
    mutate(
        delta_Tx_minus_Sham =
            Tx -
            Sham
    )

write.csv(
    balance_effect,
    file.path(
        TAB_DIR,
        "08_Tx_minus_Sham_balance_effect_sensitivity_v4.13.1.csv"
    ),
    row.names = FALSE
)

p_balance_sens <- ggplot(
    balance_effect,
    aes(
        x = analysis,
        y = delta_Tx_minus_Sham,
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
        size = 2.5
    ) +
    facet_wrap(
        ~ metric,
        scales = "free_y",
        ncol = 2
    ) +
    labs(
        title = "Robustness of compositional balances to cleaning",
        subtitle = "Tx − Sham effect across Original / Clean-A / Clean-B / Clean-C",
        x = NULL,
        y = "Δ balance (Tx − Sham)"
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
    "05_balance_effect_sensitivity_v4.13.1.pdf",
    p_balance_sens,
    10,
    8
)

# ------------------------------------------------------------------------------
# 14. Stacked composition plot
# ------------------------------------------------------------------------------

p_stack <- ggplot(
    composition_long %>%
        filter(
            analysis == "Clean-B"
        ),
    aes(
        x = sample,
        y = percent,
        fill = macrophage_class
    )
) +
    geom_col(
        width = 0.75
    ) +
    scale_fill_manual(
        values = c(
            "Inflammatory-Mphi" = "#E41A1C",
            "Anti-inflammatory-Mphi" = "#00AEEF",
            "Fibrogenic-Mphi" = "#FF1493",
            "Repair/Resolution-Mphi" = "#00C853",
            "Lipid-associated/TREM2-Mphi" = "#7B2CBF",
            "Other" = "#8C8C8C"
        ),
        labels = SUBTYPE_LABELS
    ) +
    labs(
        title = "Clean-B MΦ sample-level composition",
        subtitle = "Compositional view; each bar = one biological sample",
        x = NULL,
        y = "% of retained MΦ",
        fill = "MΦ subtype"
    ) +
    theme_classic(
        base_size = 10
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
    "06_CleanB_stacked_sample_composition_v4.13.1.pdf",
    p_stack,
    9,
    6
)

# ------------------------------------------------------------------------------
# 15. Master summary
# ------------------------------------------------------------------------------

p_master <- (
    p_stack /
    p_pca /
    p_clr_effect /
    p_balances /
    p_balance_sens
) +
    patchwork::plot_layout(
        heights = c(
            0.8,
            0.8,
            1.0,
            1.2,
            1.0
        )
    ) +
    patchwork::plot_annotation(
        title = "Mouse MASH MΦ formal compositional analysis v4.13.1",
        subtitle = "Biological sample-level CLR / Aitchison / log-ratio analysis",
        theme = theme(
            plot.title = element_text(
                face = "bold",
                size = 18
            )
        )
    )

save_pdf(
    "07_formal_compositional_analysis_master_v4.13.1.pdf",
    p_master,
    13,
    28
)

# ------------------------------------------------------------------------------
# 16. README
# ------------------------------------------------------------------------------

readme <- c(
    "Mouse MASH MΦ formal compositional analysis v4.13.1",
    "",
    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",
    "Primary dataset:",
    "  Clean-B",
    "",
    "Sensitivity datasets:",
    "  Original",
    "  Clean-A",
    "  Clean-B",
    "  Clean-C",
    "",
    "Statistical unit:",
    "  biological sample",
    "",
    "Methods:",
    "  sample-level subtype composition",
    paste0(
        "  CLR with pseudocount = ",
        PSEUDOCOUNT
    ),
    "  Aitchison distance",
    "  CLR-space PCA",
    "  biologically interpretable log-ratio balances",
    "",
    "Important limitation:",
    "  Sham vs Tx = n=2 vs n=2.",
    "  STD vs CDAHFD = n=1 vs n=1 and is descriptive only.",
    "  No cell-level inferential p-values are generated.",
    "",
    "Primary outputs:",
    "  01 CLR effect sensitivity",
    "  02 Clean-B CLR PCA",
    "  03 Aitchison distance",
    "  04 Clean-B compositional balances",
    "  05 balance sensitivity",
    "  06 stacked composition",
    "  07 master summary"
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.13.1.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.13.1.txt"
    )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)

print(cleanB_effect)
print(balance_effect)
