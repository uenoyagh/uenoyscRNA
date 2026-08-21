#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Final classification robustness audit
# v4.14.3
#
# PURPOSE
#   Integrate the major QC/robustness checks supporting the fixed v4.8.4
#   five-subtype macrophage classification.
#
# FIXED PARENT CLASSIFICATION
#   1) Inflammatory-MΦ
#   2) Anti-inflammatory-MΦ
#   3) Fibrogenic-MΦ
#   4) Repair/Resolution-MΦ
#   5) Lipid-associated/TREM2-MΦ
#   + Other
#
# AUDIT DOMAINS
#   A. Original / Clean-A / Clean-B / Clean-C abundance sensitivity
#   B. contamination-candidate burden by subtype
#   C. cluster-to-class purity and entropy
#   D. marker/module program specificity
#   E. sample-level reproducibility
#   F. Anti-inflammatory-MΦ internal heterogeneity is treated as a secondary
#      state layer and does NOT redefine the parent class
#
# IMPORTANT
#   This script is a descriptive robustness audit, not a formal classifier
#   benchmark against an external ground truth.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4143)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_AUDIT_RDS <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_contamination_audit_v4.11.6",
    "Mouse_Mphi_Res2_contamination_audit_annotated_v4.11.6.rds"
)

ANTI_STATE_RDS <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "AntiInflammatory_functional_state_audit_CleanB_v4.14.0.1",
    "RDS",
    "Mouse_Mphi_AntiInflammatory_Res1.2_functional_state_annotated_v4.14.0.1.rds"
)

OUTPUT_DIR <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_final_classification_robustness_v4.14.3"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"

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

PROGRAMS <- list(

    Inflammatory = c(
        "Il1b","Tnf","Ccl2","Cxcl10",
        "Nos2","Cd80","Cd86","Stat1"
    ),

    Anti_inflammatory = c(
        "Mrc1","Cd163","Il1rn","Arg1",
        "Mertk","Igf1","Hmox1","Klf4","Maf"
    ),

    Fibrogenic = c(
        "Spp1","Tgfb1","Pdgfb","Thbs1",
        "Lgals3","Gpnmb","Mmp12"
    ),

    Repair_Resolution = c(
        "Mertk","Axl","Mfge8","Gas6",
        "Igf1","Hmox1","Mmp13","Mmp14","Plau"
    ),

    Lipid_TREM2 = c(
        "Trem2","Gpnmb","Cd9","Lpl",
        "Apoe","Fabp5","Abca1","Plin2"
    )
)

CLASS_TO_PROGRAM <- c(
    "Inflammatory-Mphi" =
        "Inflammatory",
    "Anti-inflammatory-Mphi" =
        "Anti_inflammatory",
    "Fibrogenic-Mphi" =
        "Fibrogenic",
    "Repair/Resolution-Mphi" =
        "Repair_Resolution",
    "Lipid-associated/TREM2-Mphi" =
        "Lipid_TREM2"
)

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

    out
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

mean_expression_score <- function(mat, genes) {

    genes <- intersect(
        genes,
        rownames(mat)
    )

    if (length(genes) == 0L) {
        return(
            rep(
                NA_real_,
                ncol(mat)
            )
        )
    }

    as.numeric(
        Matrix::colMeans(
            mat[
                genes,
                ,
                drop = FALSE
            ]
        )
    )
}

entropy_norm <- function(x) {

    p <- x / sum(x)
    p <- p[
        p > 0
    ]

    if (length(p) <= 1L) return(0)

    -sum(
        p * log(p)
    ) /
        log(
            length(
                SUBTYPE_ORDER
            )
        )
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
# 3. Load v4.11.6 audit object
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_AUDIT_RDS)) {
    stop(
        "v4.11.6 audit RDS not found:\n",
        INPUT_AUDIT_RDS
    )
}

msg("Loading audit object: ", INPUT_AUDIT_RDS)

mphi <- readRDS(
    INPUT_AUDIT_RDS
)

DefaultAssay(mphi) <- ASSAY_USE

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

if (is.na(SAMPLE_COL)) stop("Sample column not found.")
if (is.na(CLASS_COL)) stop("Class column not found.")
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
        "Missing QC metadata: ",
        paste(
            missing_qc,
            collapse = ", "
        )
    )
}

# ------------------------------------------------------------------------------
# 4. Canonical metadata + cleaning masks
# ------------------------------------------------------------------------------

meta <- mphi@meta.data %>%
    rownames_to_column("cell") %>%
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
        ),
        cluster = as.character(
            .data[[CLUSTER_COL]]
        ),
        B_candidate = QC_B_candidate,
        T_candidate = QC_T_candidate,
        NK_candidate = QC_NK_candidate,
        Neutrophil_candidate = QC_Neutrophil_candidate,
        n_lineage_flags = QC_n_lineage_flags_v4116,
        Mphi_identity_score = QC_Mphi_identity_score_v4116
    )

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
            !(cluster %in%
                c(
                    "24",
                    "27"
                )),

        keep_Clean_C =
            !(cluster %in%
                c(
                    "24",
                    "27"
                )) &
            !strict_cell_candidate
    )

KEEP_COLS <- c(
    "Original" = "keep_Original",
    "Clean-A"  = "keep_Clean_A",
    "Clean-B"  = "keep_Clean_B",
    "Clean-C"  = "keep_Clean_C"
)

# ------------------------------------------------------------------------------
# 5. Abundance sensitivity
# ------------------------------------------------------------------------------

composition_list <- list()

for (analysis_name in ANALYSIS_ORDER) {

    keep_col <- KEEP_COLS[[analysis_name]]

    df <- meta %>%
        filter(
            .data[[keep_col]]
        )

    totals <- df %>%
        count(
            sample,
            name = "total_Mphi"
        )

    comp <- df %>%
        count(
            sample,
            condition,
            macrophage_class,
            name = "n_cells"
        ) %>%
        complete(
            sample,
            macrophage_class =
                SUBTYPE_ORDER,
            fill = list(
                n_cells = 0L
            )
        ) %>%
        mutate(
            condition =
                canonical_condition(
                    sample
                )
        ) %>%
        left_join(
            totals,
            by = "sample"
        ) %>%
        mutate(
            analysis =
                analysis_name,
            percent =
                100 *
                n_cells /
                total_Mphi
        )

    composition_list[[analysis_name]] <- comp
}

composition <- bind_rows(
    composition_list
)

write.csv(
    composition,
    file.path(
        TAB_DIR,
        "01_composition_sensitivity_all_cleaning_conditions_v4.14.3.csv"
    ),
    row.names = FALSE
)

composition_wide <- composition %>%
    select(
        analysis,
        sample,
        macrophage_class,
        percent
    ) %>%
    pivot_wider(
        names_from = analysis,
        values_from = percent
    ) %>%
    mutate(
        max_abs_delta_vs_CleanB =
            pmax(
                abs(
                    Original -
                        `Clean-B`
                ),
                abs(
                    `Clean-A` -
                        `Clean-B`
                ),
                abs(
                    `Clean-C` -
                        `Clean-B`
                ),
                na.rm = TRUE
            )
    )

write.csv(
    composition_wide,
    file.path(
        TAB_DIR,
        "02_composition_delta_vs_CleanB_v4.14.3.csv"
    ),
    row.names = FALSE
)

p_comp <- ggplot(
    composition,
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
        size = 2.3
    ) +
    facet_wrap(
        ~ macrophage_class,
        scales = "free_y",
        ncol = 3,
        labeller =
            as_labeller(
                SUBTYPE_LABELS
            )
    ) +
    labs(
        title =
            "MΦ subtype abundance robustness to cleaning",
        subtitle =
            "Each line = biological sample",
        x =
            NULL,
        y =
            "% of retained MΦ",
        shape =
            "Sample"
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
    "01_composition_cleaning_sensitivity_v4.14.3.pdf",
    p_comp,
    11,
    8
)

# ------------------------------------------------------------------------------
# 6. Contamination burden by parent subtype
# ------------------------------------------------------------------------------

contam_by_class <- meta %>%
    group_by(
        macrophage_class
    ) %>%
    summarise(
        n_cells = n(),
        B_pct = 100 * mean(
            B_candidate,
            na.rm = TRUE
        ),
        T_pct = 100 * mean(
            T_candidate,
            na.rm = TRUE
        ),
        NK_pct = 100 * mean(
            NK_candidate,
            na.rm = TRUE
        ),
        Neutrophil_pct = 100 * mean(
            Neutrophil_candidate,
            na.rm = TRUE
        ),
        strict_candidate_pct =
            100 *
            mean(
                strict_cell_candidate,
                na.rm = TRUE
            ),
        .groups = "drop"
    )

write.csv(
    contam_by_class,
    file.path(
        TAB_DIR,
        "03_contamination_candidate_burden_by_subtype_v4.14.3.csv"
    ),
    row.names = FALSE
)

contam_long <- contam_by_class %>%
    select(
        macrophage_class,
        B_pct,
        T_pct,
        NK_pct,
        Neutrophil_pct,
        strict_candidate_pct
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
        shape = QC_type,
        group = QC_type
    )
) +
    geom_point(
        size = 2.5,
        position =
            position_dodge(
                width = 0.35
            )
    ) +
    labs(
        title =
            "Contamination-candidate burden by MΦ subtype",
        x =
            NULL,
        y =
            "% cells flagged",
        shape =
            "QC flag"
    ) +
    theme_classic(
        base_size = 9
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
    "02_contamination_burden_by_subtype_v4.14.3.pdf",
    p_contam,
    9,
    5.5
)

# ------------------------------------------------------------------------------
# 7. Cluster -> class purity / entropy in Clean-B
# ------------------------------------------------------------------------------

cleanB_meta <- meta %>%
    filter(
        keep_Clean_B
    )

cluster_class_counts <- cleanB_meta %>%
    count(
        cluster,
        macrophage_class,
        name = "n_cells"
    )

cluster_audit <- cluster_class_counts %>%
    group_by(
        cluster
    ) %>%
    summarise(
        total_cells = sum(
            n_cells
        ),
        dominant_class =
            macrophage_class[
                which.max(
                    n_cells
                )
            ],
        dominant_cells =
            max(
                n_cells
            ),
        purity =
            dominant_cells /
            total_cells,
        normalized_entropy =
            entropy_norm(
                n_cells
            ),
        n_classes_present =
            sum(
                n_cells > 0
            ),
        .groups = "drop"
    ) %>%
    arrange(
        as.numeric(
            cluster
        )
    )

write.csv(
    cluster_audit,
    file.path(
        TAB_DIR,
        "04_Res2_cluster_class_purity_entropy_CleanB_v4.14.3.csv"
    ),
    row.names = FALSE
)

p_purity <- ggplot(
    cluster_audit,
    aes(
        x = reorder(
            cluster,
            purity
        ),
        y = purity
    )
) +
    geom_col(
        width = 0.75
    ) +
    geom_hline(
        yintercept = 0.8,
        linetype = 2,
        linewidth = 0.4
    ) +
    coord_flip() +
    labs(
        title =
            "Res2 cluster-to-class purity in Clean-B",
        subtitle =
            "Dashed line = 0.80 descriptive reference",
        x =
            "Res2 cluster",
        y =
            "Dominant-class fraction"
    ) +
    theme_classic(
        base_size = 9
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            )
    )

save_pdf(
    "03_Res2_cluster_class_purity_v4.14.3.pdf",
    p_purity,
    7,
    8
)

# ------------------------------------------------------------------------------
# 8. Functional program specificity in Clean-B
# ------------------------------------------------------------------------------

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

if (is.null(rna_data)) {
    stop(
        "RNA normalized data layer unavailable."
    )
}

cleanB_cells <- cleanB_meta$cell

rna_cleanB <- rna_data[
    ,
    cleanB_cells,
    drop = FALSE
]

program_scores <- tibble(
    cell = cleanB_cells
)

for (program_name in names(PROGRAMS)) {

    program_scores[[program_name]] <-
        mean_expression_score(
            rna_cleanB,
            PROGRAMS[[program_name]]
        )
}

program_scores <- program_scores %>%
    left_join(
        cleanB_meta %>%
            select(
                cell,
                sample,
                condition,
                macrophage_class,
                cluster
            ),
        by = "cell"
    )

write.csv(
    program_scores,
    file.path(
        TAB_DIR,
        "05_cell_level_program_scores_CleanB_v4.14.3.csv"
    ),
    row.names = FALSE
)

class_program_means <- program_scores %>%
    filter(
        macrophage_class %in%
            names(
                CLASS_TO_PROGRAM
            )
    ) %>%
    pivot_longer(
        cols = all_of(
            names(PROGRAMS)
        ),
        names_to = "program",
        values_to = "score"
    ) %>%
    group_by(
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
    )

write.csv(
    class_program_means,
    file.path(
        TAB_DIR,
        "06_class_program_mean_scores_CleanB_v4.14.3.csv"
    ),
    row.names = FALSE
)

class_program_wide <- class_program_means %>%
    pivot_wider(
        names_from = macrophage_class,
        values_from = mean_score
    )

prog_mat <- as.matrix(
    class_program_wide[
        ,
        intersect(
            names(
                CLASS_TO_PROGRAM
            ),
            colnames(
                class_program_wide
            )
        ),
        drop = FALSE
    ]
)

rownames(
    prog_mat
) <- class_program_wide$program

prog_z <- t(
    scale(
        t(
            prog_mat
        )
    )
)

prog_z[
    !is.finite(
        prog_z
    )
] <- 0

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "04_classification_program_specificity_heatmap_v4.14.3.pdf"
    ),
    width = 8,
    height = 6
)

pheatmap::pheatmap(
    prog_z,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
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
    fontsize_row = 9,
    fontsize_col = 8,
    angle_col = 45,
    main = paste0(
        "Functional-program specificity of fixed MΦ classes\n",
        "Clean-B | row z-score"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 9. Own-program specificity margin
# ------------------------------------------------------------------------------

specificity_rows <- list()

for (class_now in names(
    CLASS_TO_PROGRAM
)) {

    program_now <-
        CLASS_TO_PROGRAM[[class_now]]

    own_score <- class_program_means %>%
        filter(
            macrophage_class ==
                class_now,
            program ==
                program_now
        ) %>%
        pull(
            mean_score
        )

    other_score <- class_program_means %>%
        filter(
            macrophage_class !=
                class_now,
            program ==
                program_now
        ) %>%
        summarise(
            max_other =
                max(
                    mean_score,
                    na.rm = TRUE
                )
        ) %>%
        pull(
            max_other
        )

    specificity_rows[[class_now]] <- tibble(
        macrophage_class =
            class_now,
        expected_program =
            program_now,
        own_mean_score =
            own_score,
        max_other_class_score =
            other_score,
        specificity_margin =
            own_score -
            other_score
    )
}

specificity <- bind_rows(
    specificity_rows
)

write.csv(
    specificity,
    file.path(
        TAB_DIR,
        "07_expected_program_specificity_margin_v4.14.3.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Sample-level program reproducibility
# ------------------------------------------------------------------------------

sample_program <- program_scores %>%
    filter(
        macrophage_class %in%
            names(
                CLASS_TO_PROGRAM
            )
    ) %>%
    pivot_longer(
        cols = all_of(
            names(PROGRAMS)
        ),
        names_to = "program",
        values_to = "score"
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
    )

write.csv(
    sample_program,
    file.path(
        TAB_DIR,
        "08_sample_level_class_program_means_v4.14.3.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. Anti-inflammatory secondary-state audit
# ------------------------------------------------------------------------------

anti_summary <- tibble(
    anti_state_RDS_present =
        file.exists(
            ANTI_STATE_RDS
        ),
    interpretation =
        "Secondary Res1.2 heterogeneity layer; parent Anti-inflammatory-MΦ class unchanged"
)

if (file.exists(
    ANTI_STATE_RDS
)) {

    anti <- readRDS(
        ANTI_STATE_RDS
    )

    if ("anti_functional_state_v41401" %in%
        colnames(
            anti@meta.data
        )) {

        anti_state_counts <- anti@meta.data %>%
            rownames_to_column("cell") %>%
            count(
                sample_v4140,
                anti_functional_state_v41401,
                name = "n_cells"
            )

        write.csv(
            anti_state_counts,
            file.path(
                TAB_DIR,
                "09_AntiInflammatory_secondary_state_counts_v4.14.3.csv"
            ),
            row.names = FALSE
        )
    }
}

write.csv(
    anti_summary,
    file.path(
        TAB_DIR,
        "10_AntiInflammatory_secondary_state_policy_v4.14.3.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. Integrated descriptive robustness summary
# ------------------------------------------------------------------------------

composition_summary <- composition_wide %>%
    group_by(
        macrophage_class
    ) %>%
    summarise(
        median_max_abs_delta_pct =
            median(
                max_abs_delta_vs_CleanB,
                na.rm = TRUE
            ),
        max_abs_delta_pct =
            max(
                max_abs_delta_vs_CleanB,
                na.rm = TRUE
            ),
        .groups = "drop"
    )

purity_summary <- cluster_audit %>%
    group_by(
        dominant_class
    ) %>%
    summarise(
        n_clusters =
            n(),
        median_purity =
            median(
                purity,
                na.rm = TRUE
            ),
        min_purity =
            min(
                purity,
                na.rm = TRUE
            ),
        .groups = "drop"
    )

robustness_summary <- tibble(
    domain = c(
        "Cleaning sensitivity",
        "Contamination burden",
        "Cluster purity",
        "Functional program specificity",
        "Secondary heterogeneity policy"
    ),
    interpretation = c(
        "Review subtype abundance stability across Original/Clean-A/B/C.",
        "Review lineage-candidate flags by subtype; Clean-B/C sensitivity determines impact.",
        "Review dominant-class purity and entropy for each Res2 cluster.",
        "Expected class programs should show positive class-specificity margins where biologically separable.",
        "Anti-inflammatory Res1.2 metastates are nested states and do not redefine the parent five-class framework."
    )
)

write.csv(
    robustness_summary,
    file.path(
        TAB_DIR,
        "11_final_robustness_interpretation_framework_v4.14.3.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 13. Master figure
# ------------------------------------------------------------------------------

p_master <- (
    p_comp /
    p_contam /
    p_purity
) +
    patchwork::plot_layout(
        heights = c(
            1.25,
            0.8,
            1.0
        )
    ) +
    patchwork::plot_annotation(
        title =
            "Mouse MASH MΦ final classification robustness audit",
        subtitle =
            "Cleaning sensitivity | contamination burden | Res2 cluster purity",
        theme = theme(
            plot.title =
                element_text(
                    face = "bold",
                    size = 16
                )
        )
    )

save_pdf(
    "05_final_classification_robustness_master_v4.14.3.pdf",
    p_master,
    12,
    22
)

# ------------------------------------------------------------------------------
# 14. README
# ------------------------------------------------------------------------------

readme <- c(
    "Mouse MASH MΦ final classification robustness audit v4.14.3",
    "",
    paste0(
        "Input audit RDS: ",
        INPUT_AUDIT_RDS
    ),
    "",
    "Fixed parent framework:",
    "  Res2.0",
    "  Inflammatory-MΦ",
    "  Anti-inflammatory-MΦ",
    "  Fibrogenic-MΦ",
    "  Repair/Resolution-MΦ",
    "  Lipid-associated/TREM2-MΦ",
    "  Other",
    "",
    "Audit domains:",
    "  Original/Clean-A/B/C abundance sensitivity",
    "  B/T/NK/neutrophil candidate burden",
    "  Res2 cluster-to-class purity and entropy",
    "  expected functional-program specificity",
    "  biological-sample program reproducibility",
    "  nested Anti-inflammatory-MΦ Res1.2 states",
    "",
    "Important:",
    "  This is a descriptive robustness audit.",
    "  It is not an external classifier validation.",
    "  Anti-inflammatory Res1.2 metastates remain nested secondary states.",
    "",
    "Primary outputs:",
    "  01 composition cleaning sensitivity",
    "  02 contamination burden",
    "  03 cluster purity",
    "  04 program specificity heatmap",
    "  05 master audit"
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.14.3.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.14.3.txt"
    )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)

print(
    composition_summary
)

print(
    specificity
)

print(
    purity_summary
)
