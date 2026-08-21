#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Anti-inflammatory-MΦ functional-state annotation audit
# Clean-B / Res1.2
# v4.14.0.1
#
# PURPOSE
#   Collapse the 14 Res1.2 Anti-inflammatory-MΦ subclusters into a smaller
#   number of biologically interpretable functional states.
#
# PRINCIPLE
#   - Parent MΦ classification remains fixed at Res2.0.
#   - Anti-inflammatory-MΦ internal clustering remains fixed at Res1.2.
#   - Res1.2 cluster IDs are preserved.
#   - Functional-state annotation is ADDED as a higher-level layer.
#
# DATA-DRIVEN STRATEGY
#   1) compute seven functional-program scores at cell level
#   2) summarize by Res1.2 subcluster
#   3) row-standardize program profiles across subclusters
#   4) hierarchical clustering of subclusters in program space
#   5) cut tree into K = 5 functional metastates
#   6) assign provisional biological names from metastate program profiles
#   7) aggregate Sham/Tx abundance at metastate level
#
# IMPORTANT
#   Provisional metastate names should be reviewed before publication.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(41401)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "AntiInflammatory_heterogeneity_CleanB_v4.14.0",
    "RDS",
    "Mouse_Mphi_AntiInflammatory_heterogeneity_CleanB_v4.14.0.rds"
)

OUTPUT_DIR <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "AntiInflammatory_functional_state_audit_CleanB_v4.14.0.1"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"
REDUCTION_USE <- "anti.umap.rpca"

N_METASTATES <- 5L

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

# ------------------------------------------------------------------------------
# 1. Functional programs
# ------------------------------------------------------------------------------

PROGRAMS <- list(

    Anti_inflammatory = c(
        "Mrc1","Cd163","Il1rn","Arg1","Mertk",
        "Igf1","Hmox1","Klf4","Maf"
    ),

    IL10_STAT3 = c(
        "Il10ra","Il10rb","Jak1","Tyk2",
        "Stat3","Socs3","Bcl3","Il1rn"
    ),

    Repair_Resolution = c(
        "Mertk","Axl","Mfge8","Gas6","Igf1",
        "Hmox1","Mmp12","Mmp13","Mmp14","Plau"
    ),

    Inflammatory = c(
        "Il1b","Tnf","Ccl2","Ccl3","Ccl4",
        "Cxcl10","Nos2","Cd80","Cd86","Stat1"
    ),

    Fibrogenic = c(
        "Spp1","Tgfb1","Pdgfb","Thbs1",
        "Lgals3","Gpnmb","Mmp12","Mmp14","Ctsb"
    ),

    Lipid_TREM2 = c(
        "Trem2","Gpnmb","Cd9","Lpl","Apoe",
        "Fabp5","Abca1","Plin2","Ctsd"
    ),

    Efferocytosis = c(
        "Mertk","Axl","Mfge8","Gas6","Marco",
        "Cd36","Lrp1","C1qa","C1qb","C1qc"
    )
)

PROGRAM_ORDER <- names(PROGRAMS)

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

safe_scale <- function(x) {
    z <- as.numeric(scale(x))
    z[!is.finite(z)] <- 0
    z
}

# ------------------------------------------------------------------------------
# 4. Load final Res1.2 object
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
    stop(
        "Final Res1.2 RDS not found:\n",
        INPUT_RDS
    )
}

msg("Loading final Res1.2 object: ", INPUT_RDS)

anti <- readRDS(INPUT_RDS)

if (!inherits(anti, "Seurat")) {
    stop("Input is not a Seurat object.")
}

if (!"anti_subcluster_v4140" %in%
    colnames(anti@meta.data)) {

    stop(
        "anti_subcluster_v4140 is missing."
    )
}

if (!REDUCTION_USE %in%
    Reductions(anti)) {

    stop(
        "Required reduction missing: ",
        REDUCTION_USE
    )
}

# ------------------------------------------------------------------------------
# 5. Join RNA layers / normalized data
# ------------------------------------------------------------------------------

DefaultAssay(anti) <- ASSAY_USE

anti <- JoinLayers(
    anti,
    assay = ASSAY_USE
)

if (!"data" %in%
    Layers(
        anti[[ASSAY_USE]]
    )) {

    anti <- NormalizeData(
        anti,
        assay = ASSAY_USE,
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
    )
}

rna_data <- get_layer_safe(
    anti,
    ASSAY_USE,
    "data"
)

if (is.null(rna_data)) {
    stop("Normalized RNA data layer not found.")
}

# ------------------------------------------------------------------------------
# 6. Program gene audit
# ------------------------------------------------------------------------------

PROGRAM_USE <- lapply(
    PROGRAMS,
    intersect,
    y = rownames(rna_data)
)

program_gene_audit <- bind_rows(
    lapply(
        names(PROGRAM_USE),
        function(program_name) {

            tibble(
                program = program_name,
                n_requested =
                    length(
                        PROGRAMS[[program_name]]
                    ),
                n_detected =
                    length(
                        PROGRAM_USE[[program_name]]
                    ),
                detected =
                    paste(
                        PROGRAM_USE[[program_name]],
                        collapse = ";"
                    ),
                missing =
                    paste(
                        setdiff(
                            PROGRAMS[[program_name]],
                            PROGRAM_USE[[program_name]]
                        ),
                        collapse = ";"
                    )
            )
        }
    )
)

write.csv(
    program_gene_audit,
    file.path(
        TAB_DIR,
        "01_program_gene_audit_v4.14.0.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Cell-level program scores
# ------------------------------------------------------------------------------

score_df <- tibble(
    cell = colnames(anti)
)

for (program_name in names(PROGRAM_USE)) {

    score_df[[program_name]] <-
        mean_expression_score(
            rna_data,
            PROGRAM_USE[[program_name]]
        )
}

meta_small <- anti@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = as.character(
            sample_v4140
        ),
        condition = as.character(
            condition_v4140
        ),
        subcluster = as.character(
            anti_subcluster_v4140
        )
    )

score_df <- score_df %>%
    left_join(
        meta_small,
        by = "cell"
    )

write.csv(
    score_df,
    file.path(
        TAB_DIR,
        "02_cell_level_program_scores_v4.14.0.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Subcluster-level mean program profiles
# ------------------------------------------------------------------------------

subcluster_program_long <- score_df %>%
    pivot_longer(
        cols = all_of(
            PROGRAM_ORDER
        ),
        names_to = "program",
        values_to = "score"
    ) %>%
    group_by(
        subcluster,
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
    subcluster_program_long,
    file.path(
        TAB_DIR,
        "03_subcluster_program_mean_scores_v4.14.0.1.csv"
    ),
    row.names = FALSE
)

subcluster_program_wide <- subcluster_program_long %>%
    pivot_wider(
        names_from = program,
        values_from = mean_score
    ) %>%
    arrange(
        as.numeric(
            subcluster
        )
    )

program_mat <- as.matrix(
    subcluster_program_wide[
        ,
        PROGRAM_ORDER,
        drop = FALSE
    ]
)

rownames(program_mat) <-
    subcluster_program_wide$subcluster

# ------------------------------------------------------------------------------
# 9. Program-wise z-score across subclusters
# ------------------------------------------------------------------------------

program_z <- scale(
    program_mat,
    center = TRUE,
    scale = TRUE
)

program_z <- as.matrix(
    program_z
)

program_z[
    !is.finite(
        program_z
    )
] <- 0

write.csv(
    data.frame(
        subcluster = rownames(program_z),
        program_z,
        check.names = FALSE
    ),
    file.path(
        TAB_DIR,
        "04_subcluster_program_zscores_v4.14.0.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Hierarchical clustering into functional metastates
# ------------------------------------------------------------------------------

subcluster_dist <- dist(
    program_z,
    method = "euclidean"
)

subcluster_hc <- hclust(
    subcluster_dist,
    method = "ward.D2"
)

metastate_id <- cutree(
    subcluster_hc,
    k = N_METASTATES
)

metastate_map <- tibble(
    subcluster = names(
        metastate_id
    ),
    metastate_id = as.integer(
        metastate_id
    )
) %>%
    arrange(
        metastate_id,
        as.numeric(
            subcluster
        )
    )

# ------------------------------------------------------------------------------
# 11. Metastate program profiles
# ------------------------------------------------------------------------------

metastate_profile <- subcluster_program_long %>%
    left_join(
        metastate_map,
        by = "subcluster"
    ) %>%
    group_by(
        metastate_id,
        program
    ) %>%
    summarise(
        mean_score =
            mean(
                mean_score,
                na.rm = TRUE
            ),
        .groups = "drop"
    ) %>%
    group_by(
        program
    ) %>%
    mutate(
        z =
            safe_scale(
                mean_score
            )
    ) %>%
    ungroup()

write.csv(
    metastate_profile,
    file.path(
        TAB_DIR,
        "05_metastate_program_profiles_v4.14.0.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. Data-driven provisional metastate names
#
# Naming logic:
#   - Fibrogenic + Lipid_TREM2 high -> Fibrogenic/Lipid-TREM2-like
#   - IL10_STAT3 + Inflammatory high -> IL10/STAT3-responsive activated
#   - Anti-inflammatory + Efferocytosis high -> Anti-inflammatory/Efferocytosis
#   - Repair/Resolution dominant -> Repair/Resolution
#   - otherwise -> Transitional/Mixed
#
# These are PROVISIONAL and should be reviewed.
# ------------------------------------------------------------------------------

meta_profile_wide <- metastate_profile %>%
    select(
        metastate_id,
        program,
        z
    ) %>%
    pivot_wider(
        names_from = program,
        values_from = z
    )

assign_metastate_name <- function(df_row) {

    fib_lipid <- mean(
        c(
            df_row[["Fibrogenic"]],
            df_row[["Lipid_TREM2"]]
        ),
        na.rm = TRUE
    )

    il10_inflam <- mean(
        c(
            df_row[["IL10_STAT3"]],
            df_row[["Inflammatory"]]
        ),
        na.rm = TRUE
    )

    anti_eff <- mean(
        c(
            df_row[["Anti_inflammatory"]],
            df_row[["Efferocytosis"]]
        ),
        na.rm = TRUE
    )

    repair <- df_row[["Repair_Resolution"]]

    candidate_scores <- c(
        Fibrogenic_Lipid_TREM2 = fib_lipid,
        IL10_STAT3_Activated = il10_inflam,
        Anti_inflammatory_Efferocytosis = anti_eff,
        Repair_Resolution = repair
    )

    best <- names(
        which.max(
            candidate_scores
        )
    )

    best_value <- max(
        candidate_scores,
        na.rm = TRUE
    )

    if (!is.finite(best_value) ||
        best_value < 0.25) {

        return(
            "Transitional/Mixed"
        )
    }

    switch(
        best,
        Fibrogenic_Lipid_TREM2 =
            "Fibrogenic/Lipid-TREM2-like",

        IL10_STAT3_Activated =
            "IL10/STAT3-responsive activated",

        Anti_inflammatory_Efferocytosis =
            "Anti-inflammatory/Efferocytosis-high",

        Repair_Resolution =
            "Repair/Resolution-high",

        "Transitional/Mixed"
    )
}

metastate_name_table <- meta_profile_wide

metastate_name_table$provisional_state <- vapply(
    seq_len(
        nrow(
            metastate_name_table
        )
    ),
    function(i) {

        assign_metastate_name(
            as.list(
                metastate_name_table[
                    i,
                    ,
                    drop = FALSE
                ]
            )
        )
    },
    character(1)
)

# If duplicate provisional names occur, append metastate id so states remain
# unambiguous during audit.
duplicate_name_flag <- duplicated(
    metastate_name_table$provisional_state
) |
    duplicated(
        metastate_name_table$provisional_state,
        fromLast = TRUE
    )

metastate_name_table$functional_state <- metastate_name_table$provisional_state

metastate_name_table$functional_state[
    duplicate_name_flag
] <- paste0(
    metastate_name_table$provisional_state[
        duplicate_name_flag
    ],
    " [M",
    metastate_name_table$metastate_id[
        duplicate_name_flag
    ],
    "]"
)

write.csv(
    metastate_name_table,
    file.path(
        TAB_DIR,
        "06_metastate_provisional_names_v4.14.0.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 13. Final audit mapping: Res1.2 cluster -> functional metastate
# ------------------------------------------------------------------------------

cluster_state_map <- metastate_map %>%
    left_join(
        metastate_name_table %>%
            select(
                metastate_id,
                functional_state
            ),
        by = "metastate_id"
    ) %>%
    arrange(
        as.numeric(
            subcluster
        )
    )

write.csv(
    cluster_state_map,
    file.path(
        TAB_DIR,
        "07_Res1.2_cluster_to_functional_state_mapping_v4.14.0.1.csv"
    ),
    row.names = FALSE
)

print(
    cluster_state_map
)

# ------------------------------------------------------------------------------
# 14. Add state annotation to Seurat object
# ------------------------------------------------------------------------------

state_lookup <- setNames(
    cluster_state_map$functional_state,
    cluster_state_map$subcluster
)

anti$anti_functional_state_v41401 <- unname(
    state_lookup[
        as.character(
            anti$anti_subcluster_v4140
        )
    ]
)

anti$anti_functional_state_v41401 <- factor(
    anti$anti_functional_state_v41401,
    levels = unique(
        cluster_state_map$functional_state
    )
)

# ------------------------------------------------------------------------------
# 15. UMAP: Res1.2 cluster and functional metastate
# ------------------------------------------------------------------------------

p_cluster <- DimPlot(
    anti,
    reduction = REDUCTION_USE,
    group.by = "anti_subcluster_v4140",
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    pt.size = 0.6
) +
    NoLegend() +
    labs(
        title =
            "Anti-inflammatory-MΦ Res1.2 subclusters"
    ) +
    theme_classic(
        base_size = 10
    )

p_state <- DimPlot(
    anti,
    reduction = REDUCTION_USE,
    group.by = "anti_functional_state_v41401",
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    pt.size = 0.6
) +
    labs(
        title =
            "Anti-inflammatory-MΦ functional metastates",
        subtitle =
            "Res1.2 subclusters collapsed by functional-program similarity",
        color =
            "Functional state"
    ) +
    theme_classic(
        base_size = 10
    )

p_umap_pair <- p_cluster + p_state

save_pdf(
    "01_Res1.2_cluster_vs_functional_metastate_UMAP_v4.14.0.1.pdf",
    p_umap_pair,
    14,
    6.5
)

# ------------------------------------------------------------------------------
# 16. Program heatmap with metastate annotation
# ------------------------------------------------------------------------------

annotation_row <- cluster_state_map %>%
    arrange(
        match(
            subcluster,
            rownames(
                program_z
            )
        )
    ) %>%
    select(
        functional_state
    ) %>%
    as.data.frame()

rownames(
    annotation_row
) <- cluster_state_map$subcluster

program_z_plot <- pmax(
    pmin(
        program_z,
        2
    ),
    -2
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "02_Res1.2_subcluster_program_heatmap_with_metastate_v4.14.0.1.pdf"
    ),
    width = 9,
    height = 7
)

pheatmap::pheatmap(
    program_z_plot,
    cluster_rows = subcluster_hc,
    cluster_cols = TRUE,
    annotation_row = annotation_row[
        rownames(
            program_z_plot
        ),
        ,
        drop = FALSE
    ],
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
    fontsize_row = 8.5,
    fontsize_col = 8,
    angle_col = 45,
    main = paste0(
        "Res1.2 Anti-inflammatory-MΦ functional profiles\n",
        "subclusters grouped into 5 metastates"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 17. Metastate-level program heatmap
# ------------------------------------------------------------------------------

metastate_heat <- metastate_profile %>%
    select(
        metastate_id,
        program,
        z
    ) %>%
    pivot_wider(
        names_from = metastate_id,
        values_from = z
    )

metastate_heat_mat <- as.matrix(
    metastate_heat[
        ,
        -1,
        drop = FALSE
    ]
)

rownames(
    metastate_heat_mat
) <- metastate_heat$program

metastate_labels <- metastate_name_table %>%
    arrange(
        metastate_id
    )

colnames(
    metastate_heat_mat
) <- metastate_labels$functional_state

metastate_heat_plot <- pmax(
    pmin(
        metastate_heat_mat,
        2
    ),
    -2
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "03_functional_metastate_program_heatmap_v4.14.0.1.pdf"
    ),
    width = 9,
    height = 6.5
)

pheatmap::pheatmap(
    metastate_heat_plot,
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
    main =
        "Anti-inflammatory-MΦ functional metastates"
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 18. Sample-level metastate abundance
# ------------------------------------------------------------------------------

abundance <- anti@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = as.character(
            sample_v4140
        ),
        condition = as.character(
            condition_v4140
        ),
        functional_state = as.character(
            anti_functional_state_v41401
        )
    ) %>%
    count(
        sample,
        condition,
        functional_state,
        name = "n_cells"
    )

sample_totals <- abundance %>%
    group_by(
        sample,
        condition
    ) %>%
    summarise(
        total_Anti_inflammatory_Mphi =
            sum(
                n_cells
            ),
        .groups = "drop"
    )

abundance <- abundance %>%
    left_join(
        sample_totals,
        by = c(
            "sample",
            "condition"
        )
    ) %>%
    mutate(
        fraction =
            n_cells /
            total_Anti_inflammatory_Mphi,

        percent =
            100 *
            fraction
    )

write.csv(
    abundance,
    file.path(
        TAB_DIR,
        "08_functional_metastate_abundance_by_sample_v4.14.0.1.csv"
    ),
    row.names = FALSE
)

p_abundance <- ggplot(
    abundance,
    aes(
        x = factor(
            condition,
            levels = CONDITION_ORDER
        ),
        y = percent
    )
) +
    geom_point(
        aes(
            shape = sample
        ),
        size = 3
    ) +
    facet_wrap(
        ~ functional_state,
        scales = "free_y",
        ncol = 3
    ) +
    labs(
        title =
            "Anti-inflammatory-MΦ functional-state abundance",
        subtitle =
            "Each point = biological sample",
        x =
            NULL,
        y =
            "Functional state / Anti-inflammatory-MΦ (%)",
        shape =
            "Sample"
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
        ),
        strip.text = element_text(
            face = "bold",
            size = 8
        )
    )

save_pdf(
    "04_functional_metastate_abundance_by_condition_v4.14.0.1.pdf",
    p_abundance,
    11,
    8
)

# ------------------------------------------------------------------------------
# 19. Sham -> Tx metastate effect
# ------------------------------------------------------------------------------

treatment_effect <- abundance %>%
    filter(
        condition %in%
            c(
                "Sham",
                "Tx"
            )
    ) %>%
    group_by(
        functional_state,
        condition
    ) %>%
    summarise(
        mean_percent =
            mean(
                percent,
                na.rm = TRUE
            ),
        min_percent =
            min(
                percent,
                na.rm = TRUE
            ),
        max_percent =
            max(
                percent,
                na.rm = TRUE
            ),
        .groups = "drop"
    ) %>%
    pivot_wider(
        names_from = condition,
        values_from = c(
            mean_percent,
            min_percent,
            max_percent
        ),
        names_sep = "_"
    ) %>%
    mutate(
        delta_percent_Tx_minus_Sham =
            mean_percent_Tx -
            mean_percent_Sham
    ) %>%
    arrange(
        desc(
            abs(
                delta_percent_Tx_minus_Sham
            )
        )
    )

write.csv(
    treatment_effect,
    file.path(
        TAB_DIR,
        "09_functional_metastate_Tx_vs_Sham_effect_v4.14.0.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 20. Stacked sample composition
# ------------------------------------------------------------------------------

p_stack <- ggplot(
    abundance,
    aes(
        x = factor(
            sample,
            levels = SAMPLE_ORDER
        ),
        y = percent,
        fill = functional_state
    )
) +
    geom_col(
        width = 0.75
    ) +
    labs(
        title =
            "Anti-inflammatory-MΦ functional-state composition",
        subtitle =
            "Res1.2 subclusters collapsed into 5 metastates",
        x =
            NULL,
        y =
            "% of Anti-inflammatory-MΦ",
        fill =
            "Functional state"
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
    "05_functional_metastate_stacked_composition_v4.14.0.1.pdf",
    p_stack,
    9,
    6
)

# ------------------------------------------------------------------------------
# 21. Save annotated object
# ------------------------------------------------------------------------------

OUTPUT_RDS <- file.path(
    RDS_DIR,
    "Mouse_Mphi_AntiInflammatory_Res1.2_functional_state_annotated_v4.14.0.1.rds"
)

saveRDS(
    anti,
    OUTPUT_RDS
)

# ------------------------------------------------------------------------------
# 22. README
# ------------------------------------------------------------------------------

readme <- c(

    "Mouse MASH Anti-inflammatory-MΦ functional-state audit v4.14.0.1",
    "",

    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",

    "Fixed framework:",
    "  Parent MΦ classification = Res2.0",
    "  Anti-inflammatory-MΦ internal clustering = Res1.2",
    "  Res1.2 cluster IDs are preserved",
    "",

    paste0(
        "Functional metastates = ",
        N_METASTATES
    ),
    "",

    "Method:",
    "  seven functional-program profiles",
    "  subcluster-level mean scores",
    "  program-wise z-score",
    "  Ward.D2 hierarchical clustering",
    "  cutree into five metastates",
    "  provisional biological naming",
    "",

    "Important:",
    "  metastate names are provisional and should be reviewed.",
    "  Publication conclusions should use biological-sample-level abundance.",
    "",

    "Primary outputs:",
    "  01 Res1.2 cluster vs metastate UMAP",
    "  02 subcluster functional heatmap + metastate annotation",
    "  03 metastate functional heatmap",
    "  04 metastate abundance by condition",
    "  05 stacked metastate composition",
    "",

    "Key mapping table:",
    "  07_Res1.2_cluster_to_functional_state_mapping_v4.14.0.1.csv",
    "",

    "Next:",
    "  review metastate names, then proceed to output-ligand analysis v4.14.1"
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.14.0.1.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.14.0.1.txt"
    )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)
msg("Annotated RDS: ", OUTPUT_RDS)

print(
    cluster_state_map
)

print(
    treatment_effect
)
