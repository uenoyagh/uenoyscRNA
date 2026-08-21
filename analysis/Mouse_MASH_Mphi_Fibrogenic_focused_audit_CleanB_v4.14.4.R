#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Fibrogenic-MΦ focused classification audit
# Clean-B / Res2.0
# v4.14.4
#
# PURPOSE
#   Resolve the only remaining robustness concern from v4.14.3:
#   Fibrogenic-MΦ showed a negative program-specificity margin.
#
# FOCUS
#   Fibrogenic-MΦ Res2 clusters: 3 and 21
#
# COMPARATORS
#   Repair/Resolution-MΦ
#   Lipid-associated/TREM2-MΦ
#
# ANALYSES
#   1) cluster/sample cell counts
#   2) representative-gene DotPlot
#   3) cluster-level mean-expression heatmap
#   4) sample-level pseudobulk CPM heatmap
#   5) targeted Fibrogenic-vs-comparator gene effect tables
#   6) program score comparison
#   7) final audit summary
#
# IMPORTANT
#   This is a focused annotation audit.
#   Parent Res2.0 cluster identities are not modified.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4144)

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
    "Mphi_Fibrogenic_focused_audit_CleanB_v4.14.4"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"

TARGET_CLASSES <- c(
    "Fibrogenic-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi"
)

FIBROGENIC_CLUSTERS <- c(
    "3",
    "21"
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

# ------------------------------------------------------------------------------
# 1. Gene panels
# ------------------------------------------------------------------------------

GENE_PANELS <- list(

    Fibrogenic_ECM = c(
        "Spp1",
        "Tgfb1",
        "Pdgfb",
        "Pdgfa",
        "Pdgfc",
        "Thbs1",
        "Timp1",
        "Fn1",
        "Col1a1",
        "Col1a2",
        "Col3a1",
        "Mmp9",
        "Mmp12",
        "Mmp14",
        "Lgals3"
    ),

    Lipid_TREM2 = c(
        "Trem2",
        "Gpnmb",
        "Cd9",
        "Lpl",
        "Apoe",
        "Fabp5",
        "Abca1",
        "Plin2",
        "Ctsd"
    ),

    Repair_Resolution = c(
        "Mertk",
        "Axl",
        "Mfge8",
        "Gas6",
        "Igf1",
        "Hmox1",
        "Mmp13",
        "Mmp14",
        "Plau"
    ),

    Inflammatory = c(
        "Il1b",
        "Tnf",
        "Ccl2",
        "Cxcl10",
        "Stat1"
    ),

    Macrophage_identity = c(
        "Adgre1",
        "Lyz2",
        "Csf1r",
        "C1qa",
        "C1qb",
        "C1qc"
    )
)

REPRESENTATIVE_GENES <- unique(
    unlist(
        GENE_PANELS,
        use.names = FALSE
    )
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

    factor(
        out,
        levels = CONDITION_ORDER
    )
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

mphi <- JoinLayers(
    mphi,
    assay = ASSAY_USE
)

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
if (is.na(CLASS_COL)) stop("Class column not found.")
if (is.na(CLUSTER_COL)) stop("Res2 cluster column not found.")

# ------------------------------------------------------------------------------
# 5. Metadata audit
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
        )
    )

fib_cluster_classes <- meta %>%
    filter(
        cluster %in%
            FIBROGENIC_CLUSTERS
    ) %>%
    count(
        cluster,
        macrophage_class,
        name = "n_cells"
    )

write.csv(
    fib_cluster_classes,
    file.path(
        TAB_DIR,
        "00_fibrogenic_cluster_class_audit_v4.14.4.csv"
    ),
    row.names = FALSE
)

print(fib_cluster_classes)

if (!all(
    fib_cluster_classes$macrophage_class ==
        "Fibrogenic-Mphi"
)) {
    warning(
        "One or more target clusters are not exclusively labeled Fibrogenic-Mphi."
    )
}

# ------------------------------------------------------------------------------
# 6. Subset focused comparison
# ------------------------------------------------------------------------------

focus_cells <- meta$cell[
    meta$macrophage_class %in%
        TARGET_CLASSES
]

focus <- subset(
    mphi,
    cells = focus_cells
)

DefaultAssay(focus) <- ASSAY_USE

meta_focus <- meta %>%
    filter(
        cell %in%
            focus_cells
    )

# ------------------------------------------------------------------------------
# 7. Cell counts
# ------------------------------------------------------------------------------

count_table <- meta_focus %>%
    count(
        macrophage_class,
        cluster,
        sample,
        condition,
        name = "n_cells"
    ) %>%
    arrange(
        macrophage_class,
        as.numeric(
            cluster
        ),
        factor(
            sample,
            levels = SAMPLE_ORDER
        )
    )

write.csv(
    count_table,
    file.path(
        TAB_DIR,
        "01_cell_counts_class_cluster_sample_v4.14.4.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Normalized RNA data
# ------------------------------------------------------------------------------

if (!"data" %in%
    Layers(
        focus[[ASSAY_USE]]
    )) {

    focus <- NormalizeData(
        focus,
        assay = ASSAY_USE,
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
    )
}

rna_data <- get_layer_safe(
    focus,
    ASSAY_USE,
    "data"
)

counts <- get_layer_safe(
    focus,
    ASSAY_USE,
    "counts"
)

if (is.null(rna_data)) stop("RNA data layer missing.")
if (is.null(counts)) stop("RNA counts layer missing.")

GENE_USE <- intersect(
    REPRESENTATIVE_GENES,
    rownames(rna_data)
)

# ------------------------------------------------------------------------------
# 9. Representative-gene custom DotPlot by class
# ------------------------------------------------------------------------------

dot_rows <- list()

for (class_now in TARGET_CLASSES) {

    cells_now <- meta_focus$cell[
        meta_focus$macrophage_class ==
            class_now
    ]

    mat <- rna_data[
        GENE_USE,
        cells_now,
        drop = FALSE
    ]

    dot_rows[[class_now]] <- tibble(
        macrophage_class = class_now,
        gene = GENE_USE,
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
        avg_expr_z = {
            z <- as.numeric(
                scale(
                    avg_expr
                )
            )
            z[!is.finite(z)] <- 0
            pmax(
                pmin(
                    z,
                    2.5
                ),
                -2.5
            )
        }
    ) %>%
    ungroup()

write.csv(
    dot_df,
    file.path(
        TAB_DIR,
        "02_gene_DotPlot_numeric_by_class_v4.14.4.csv"
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
            7
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
    labs(
        title =
            "Fibrogenic-MΦ focused gene-level audit",
        subtitle =
            "Fibrogenic vs Repair/Resolution vs Lipid-associated/TREM2",
        x = NULL,
        y = NULL
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
                angle = 60,
                hjust = 1
            )
    )

save_pdf(
    "01_Fibrogenic_focused_gene_DotPlot_v4.14.4.pdf",
    p_dot,
    15,
    5.5
)

# ------------------------------------------------------------------------------
# 10. Cluster-level mean-expression heatmap
# ------------------------------------------------------------------------------

cluster_expr_rows <- list()

clusters_focus <- sort(
    unique(
        meta_focus$cluster
    )
)

for (cl in clusters_focus) {

    cells_now <- meta_focus$cell[
        meta_focus$cluster ==
            cl
    ]

    mat <- rna_data[
        GENE_USE,
        cells_now,
        drop = FALSE
    ]

    cluster_expr_rows[[cl]] <- tibble(
        cluster = cl,
        gene = GENE_USE,
        mean_expr =
            as.numeric(
                Matrix::rowMeans(
                    mat
                )
            )
    )
}

cluster_expr <- bind_rows(
    cluster_expr_rows
)

cluster_expr_wide <- cluster_expr %>%
    pivot_wider(
        names_from = cluster,
        values_from = mean_expr
    )

cluster_mat <- as.matrix(
    cluster_expr_wide[
        ,
        -1,
        drop = FALSE
    ]
)

rownames(
    cluster_mat
) <- cluster_expr_wide$gene

cluster_z <- t(
    scale(
        t(
            cluster_mat
        )
    )
)

cluster_z[
    !is.finite(
        cluster_z
    )
] <- 0

cluster_z_plot <- pmax(
    pmin(
        cluster_z,
        2
    ),
    -2
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "02_cluster_level_gene_heatmap_v4.14.4.pdf"
    ),
    width = 12,
    height = 10
)

pheatmap::pheatmap(
    cluster_z_plot,
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
    fontsize_row = 7.5,
    fontsize_col = 8,
    angle_col = 45,
    main = paste0(
        "Fibrogenic-focused Res2 cluster gene profiles\n",
        "row z-score"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 11. Sample-level pseudobulk CPM
# ------------------------------------------------------------------------------

pb_rows <- list()

for (class_now in TARGET_CLASSES) {

    for (sample_now in SAMPLE_ORDER) {

        cells_now <- meta_focus$cell[
            meta_focus$macrophage_class ==
                class_now &
            meta_focus$sample ==
                sample_now
        ]

        if (length(cells_now) == 0L) {
            next
        }

        pb <- Matrix::rowSums(
            counts[
                ,
                cells_now,
                drop = FALSE
            ]
        )

        lib_size <- sum(pb)

        pb_rows[[paste(
            class_now,
            sample_now,
            sep = "__"
        )]] <- tibble(
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
                rownames(counts),
            raw_count =
                as.numeric(pb),
            library_size =
                lib_size,
            CPM =
                as.numeric(pb) /
                lib_size *
                1e6
        )
    }
}

pb_long <- bind_rows(
    pb_rows
)

write.csv(
    pb_long %>%
        filter(
            gene %in%
                GENE_USE
        ),
    file.path(
        TAB_DIR,
        "03_target_gene_pseudobulk_CPM_by_class_sample_v4.14.4.csv"
    ),
    row.names = FALSE
)

pb_heat <- pb_long %>%
    filter(
        gene %in%
            GENE_USE
    ) %>%
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
        CPM
    ) %>%
    pivot_wider(
        names_from = column_id,
        values_from = CPM
    )

pb_mat <- as.matrix(
    pb_heat[
        ,
        -1,
        drop = FALSE
    ]
)

rownames(
    pb_mat
) <- pb_heat$gene

pb_log <- log2(
    pb_mat + 1
)

pb_z <- t(
    scale(
        t(
            pb_log
        )
    )
)

pb_z[
    !is.finite(
        pb_z
    )
] <- 0

pb_z_plot <- pmax(
    pmin(
        pb_z,
        2
    ),
    -2
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "03_sample_level_pseudobulk_gene_heatmap_v4.14.4.pdf"
    ),
    width = 16,
    height = 10
)

pheatmap::pheatmap(
    pb_z_plot,
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
    fontsize_row = 7.5,
    fontsize_col = 6,
    angle_col = 45,
    main = paste0(
        "Fibrogenic-focused sample-level pseudobulk expression\n",
        "log2(CPM+1) | row z-score"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 12. Targeted class contrasts
# ------------------------------------------------------------------------------

class_mean <- pb_long %>%
    filter(
        gene %in%
            GENE_USE,
        sample %in%
            c(
                "Sham1",
                "Sham20",
                "Tx17",
                "Tx5"
            )
    ) %>%
    group_by(
        macrophage_class,
        gene
    ) %>%
    summarise(
        mean_CPM =
            mean(
                CPM,
                na.rm = TRUE
            ),
        .groups = "drop"
    )

class_wide <- class_mean %>%
    pivot_wider(
        names_from = macrophage_class,
        values_from = mean_CPM
    )

required_cols <- TARGET_CLASSES

if (!all(
    required_cols %in%
        colnames(
            class_wide
        )
)) {
    stop(
        "Missing one or more target classes in class-level contrast."
    )
}

contrast_table <- class_wide %>%
    mutate(
        log2FC_Fibrogenic_vs_Repair =
            log2(
                (
                    .data[["Fibrogenic-Mphi"]] + 1
                ) /
                (
                    .data[["Repair/Resolution-Mphi"]] + 1
                )
            ),

        log2FC_Fibrogenic_vs_LipidTREM2 =
            log2(
                (
                    .data[["Fibrogenic-Mphi"]] + 1
                ) /
                (
                    .data[["Lipid-associated/TREM2-Mphi"]] + 1
                )
            ),

        min_Fibrogenic_specificity =
            pmin(
                log2FC_Fibrogenic_vs_Repair,
                log2FC_Fibrogenic_vs_LipidTREM2
            )
    ) %>%
    arrange(
        desc(
            min_Fibrogenic_specificity
        )
    )

write.csv(
    contrast_table,
    file.path(
        TAB_DIR,
        "04_Fibrogenic_vs_comparators_target_gene_effects_v4.14.4.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 13. Program-score audit
# ------------------------------------------------------------------------------

score_df <- tibble(
    cell =
        colnames(
            focus
        )
)

for (program_name in names(
    GENE_PANELS
)) {

    score_df[[program_name]] <-
        mean_expression_score(
            rna_data,
            GENE_PANELS[[program_name]]
        )
}

score_df <- score_df %>%
    left_join(
        meta_focus %>%
            select(
                cell,
                sample,
                condition,
                macrophage_class,
                cluster
            ),
        by = "cell"
    )

program_summary <- score_df %>%
    pivot_longer(
        cols = all_of(
            names(
                GENE_PANELS
            )
        ),
        names_to =
            "program",
        values_to =
            "score"
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
    program_summary,
    file.path(
        TAB_DIR,
        "05_program_score_summary_by_class_v4.14.4.csv"
    ),
    row.names = FALSE
)

program_wide <- program_summary %>%
    pivot_wider(
        names_from =
            macrophage_class,
        values_from =
            mean_score
    )

program_mat <- as.matrix(
    program_wide[
        ,
        TARGET_CLASSES,
        drop = FALSE
    ]
)

rownames(
    program_mat
) <- program_wide$program

program_z <- t(
    scale(
        t(
            program_mat
        )
    )
)

program_z[
    !is.finite(
        program_z
    )
] <- 0

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "04_program_specificity_focused_heatmap_v4.14.4.pdf"
    ),
    width = 7.5,
    height = 5.5
)

pheatmap::pheatmap(
    program_z,
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
        "Focused functional-program specificity audit\n",
        "Fibrogenic vs Repair/Resolution vs Lipid/TREM2"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 14. Fibrogenic-panel specificity score
# ------------------------------------------------------------------------------

fib_program <- program_summary %>%
    filter(
        program ==
            "Fibrogenic_ECM"
    )

fib_specificity <- fib_program %>%
    summarise(
        Fibrogenic_score =
            mean_score[
                macrophage_class ==
                    "Fibrogenic-Mphi"
            ],
        Repair_score =
            mean_score[
                macrophage_class ==
                    "Repair/Resolution-Mphi"
            ],
        LipidTREM2_score =
            mean_score[
                macrophage_class ==
                    "Lipid-associated/TREM2-Mphi"
            ],
        margin_vs_best_comparator =
            Fibrogenic_score -
            max(
                Repair_score,
                LipidTREM2_score
            )
    )

write.csv(
    fib_specificity,
    file.path(
        TAB_DIR,
        "06_Fibrogenic_ECM_program_specificity_margin_v4.14.4.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 15. Final audit summary
# ------------------------------------------------------------------------------

top_specific_genes <- contrast_table %>%
    filter(
        min_Fibrogenic_specificity > 0
    ) %>%
    slice_head(
        n = 20
    )

write.csv(
    top_specific_genes,
    file.path(
        TAB_DIR,
        "07_top_Fibrogenic_specific_candidate_genes_v4.14.4.csv"
    ),
    row.names = FALSE
)

audit_summary <- tibble(
    metric = c(
        "Fibrogenic Res2 clusters",
        "Fibrogenic ECM program margin",
        "Number of target genes higher than both comparators",
        "Interpretation rule"
    ),
    value = c(
        paste(
            FIBROGENIC_CLUSTERS,
            collapse = ", "
        ),
        as.character(
            round(
                fib_specificity$margin_vs_best_comparator,
                4
            )
        ),
        as.character(
            nrow(
                contrast_table %>%
                    filter(
                        min_Fibrogenic_specificity > 0
                    )
            )
        ),
        "Retain Fibrogenic-Mphi if gene-level ECM/fibrogenic identity is coherent despite overlap with remodeling/lipid programs."
    )
)

write.csv(
    audit_summary,
    file.path(
        TAB_DIR,
        "08_final_Fibrogenic_audit_summary_v4.14.4.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 16. README
# ------------------------------------------------------------------------------

readme <- c(

    "Mouse MASH MΦ Fibrogenic-focused audit v4.14.4",
    "",

    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",

    "Purpose:",
    "  Resolve negative Fibrogenic program-specificity margin seen in v4.14.3.",
    "",

    "Focused classes:",
    "  Fibrogenic-MΦ",
    "  Repair/Resolution-MΦ",
    "  Lipid-associated/TREM2-MΦ",
    "",

    "Fibrogenic Res2 clusters:",
    paste0(
        "  ",
        paste(
            FIBROGENIC_CLUSTERS,
            collapse = ", "
        )
    ),
    "",

    "Primary outputs:",
    "  01 gene-level DotPlot",
    "  02 cluster-level heatmap",
    "  03 sample-level pseudobulk heatmap",
    "  04 focused program-specificity heatmap",
    "",

    "Key tables:",
    "  04 Fibrogenic-vs-comparator gene effects",
    "  06 Fibrogenic ECM specificity margin",
    "  07 top Fibrogenic-specific candidate genes",
    "  08 final audit summary",
    "",

    "Interpretation:",
    "  Shared fibrogenic/remodeling/lipid programs are expected.",
    "  The parent Fibrogenic-MΦ label should be retained only if a coherent",
    "  gene-level fibrogenic/ECM-regulatory signature remains distinguishable."
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.14.4.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.14.4.txt"
    )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)

print(
    fib_specificity
)

print(
    top_specific_genes
)
