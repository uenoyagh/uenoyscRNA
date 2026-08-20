#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Clean-B Anti-inflammatory-MΦ mechanistic decomposition
# Sham vs Tx
# v4.12.1
#
# PURPOSE
#   Decompose the v4.12.0 result:
#     - Anti-inflammatory-MΦ abundance increases after Tx
#     - but IL10/STAT3 program does not uniformly increase
#     - inflammatory + repair/resolution programs partly increase
#
# ANALYSIS
#   1) SAFE sample-level pseudobulk
#   2) gene-by-gene decomposition of:
#        Inflammatory
#        Anti-inflammatory
#        IL10/STAT3
#        Efferocytosis
#        Repair/Resolution
#        Fibrogenic
#        Lipid/TREM2
#   3) biological-replicate direction consistency
#   4) HSC-relevant macrophage ligand screen
#   5) compact mechanistic summary
#
# IMPORTANT
#   - Clean-B is primary dataset
#   - Res2.0 / v4.8.4 annotation unchanged
#   - n=2 Sham vs n=2 Tx
#   - effect size + replicate consistency are primary
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4121)

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
    "CleanB_Anti_inflammatory_Mphi_mechanism_v4.12.1"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
PB_DIR  <- file.path(OUTPUT_DIR, "Pseudobulk")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"
TARGET_CLASS <- "Anti-inflammatory-Mphi"

SHAM_SAMPLES <- c("Sham1", "Sham20")
TX_SAMPLES   <- c("Tx17", "Tx5")
SAMPLE_ORDER <- c(SHAM_SAMPLES, TX_SAMPLES)

# ------------------------------------------------------------------------------
# 1. Gene programs
# ------------------------------------------------------------------------------

PROGRAMS <- list(

    Inflammatory = c(
        "Il1b","Tnf","Ccl2","Ccl3","Ccl4","Cxcl10",
        "Nos2","Cd80","Cd86","Stat1"
    ),

    Anti_inflammatory = c(
        "Mrc1","Cd163","Il1rn","Retnla","Chil3","Arg1",
        "Mertk","Igf1","Hmox1","Klf4","Maf"
    ),

    IL10_STAT3 = c(
        "Il10","Il10ra","Il10rb","Jak1","Tyk2",
        "Stat3","Socs3","Bcl3","Il1rn"
    ),

    Efferocytosis = c(
        "Mertk","Axl","Mfge8","Gas6","Marco",
        "Cd36","Lrp1","C1qa","C1qb","C1qc"
    ),

    Repair_Resolution = c(
        "Mertk","Axl","Mfge8","Gas6","Igf1","Hmox1",
        "Mmp12","Mmp13","Mmp14","Plau"
    ),

    Fibrogenic = c(
        "Spp1","Tgfb1","Pdgfb","Thbs1","Lgals3",
        "Gpnmb","Mmp12","Mmp14","Ctsb"
    ),

    Lipid_TREM2 = c(
        "Trem2","Gpnmb","Cd9","Lpl","Apoe",
        "Fabp5","Abca1","Plin2","Ctsd"
    )
)

HSC_LIGANDS <- list(

    Profibrotic_HSC = c(
        "Tgfb1",
        "Pdgfb",
        "Spp1",
        "Thbs1",
        "Il1b",
        "Tnf"
    ),

    Resolution_HSC = c(
        "Igf1",
        "Gas6",
        "Mfge8",
        "Il10",
        "Areg"
    ),

    Chemotactic_HSC = c(
        "Ccl2",
        "Ccl3",
        "Ccl4",
        "Cxcl10"
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
    if (length(hit) == 0L) return(NA_character_)
    hit[[1]]
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

calc_cpm <- function(count_matrix) {
    libs <- colSums(count_matrix)

    if (any(libs <= 0)) {
        stop("At least one pseudobulk library has zero counts.")
    }

    sweep(
        count_matrix,
        2,
        libs / 1e6,
        "/"
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
# 4. Load Clean-B
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
    stop(
        "Clean-B RDS not found: ",
        INPUT_RDS
    )
}

msg("Loading Clean-B...")

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

if (is.na(SAMPLE_COL)) {
    stop("Sample metadata column not found.")
}

if (is.na(CLASS_COL)) {
    stop("Macrophage class metadata column not found.")
}

counts <- get_layer_safe(
    mphi,
    ASSAY_USE,
    "counts"
)

if (is.null(counts)) {
    stop("RNA counts layer not found.")
}

if (!identical(
    colnames(counts),
    colnames(mphi)
)) {
    stop("Counts cell order does not match Seurat object.")
}

# ------------------------------------------------------------------------------
# 5. Target cells
# ------------------------------------------------------------------------------

meta <- mphi@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = as.character(.data[[SAMPLE_COL]]),
        macrophage_class = as.character(.data[[CLASS_COL]])
    ) %>%
    mutate(
        condition = case_when(
            sample %in% SHAM_SAMPLES ~ "Sham",
            sample %in% TX_SAMPLES ~ "Tx",
            TRUE ~ NA_character_
        )
    )

target_meta <- meta %>%
    filter(
        macrophage_class == TARGET_CLASS,
        !is.na(condition)
    )

target_counts <- target_meta %>%
    count(
        sample,
        condition,
        name = "n_cells"
    )

write.csv(
    target_counts,
    file.path(
        TAB_DIR,
        "00_Anti_inflammatory_Mphi_cell_counts_v4.12.1.csv"
    ),
    row.names = FALSE
)

print(target_counts)

# ------------------------------------------------------------------------------
# 6. SAFE pseudobulk
# ------------------------------------------------------------------------------

pb_list <- list()

for (smp in SAMPLE_ORDER) {

    cells <- target_meta$cell[
        target_meta$sample == smp
    ]

    cells <- intersect(
        cells,
        colnames(counts)
    )

    if (length(cells) == 0L) {
        stop(
            "No Anti-inflammatory-MΦ cells in ",
            smp
        )
    }

    pb_list[[smp]] <- Matrix::rowSums(
        counts[
            ,
            cells,
            drop = FALSE
        ]
    )
}

pb_counts <- do.call(
    cbind,
    pb_list
)

rownames(pb_counts) <- rownames(counts)
colnames(pb_counts) <- SAMPLE_ORDER

for (i in seq_len(ncol(pb_counts) - 1L)) {
    for (j in (i + 1L):ncol(pb_counts)) {
        if (identical(
            pb_counts[, i],
            pb_counts[, j]
        )) {
            stop(
                "Identical pseudobulk columns detected: ",
                colnames(pb_counts)[[i]],
                " / ",
                colnames(pb_counts)[[j]]
            )
        }
    }
}

pb_cpm <- calc_cpm(pb_counts)
pb_log2cpm <- log2(pb_cpm + 1)

write.csv(
    data.frame(
        gene = rownames(pb_counts),
        pb_counts,
        check.names = FALSE
    ),
    file.path(
        PB_DIR,
        "01_pseudobulk_raw_counts_v4.12.1.csv"
    ),
    row.names = FALSE
)

write.csv(
    data.frame(
        gene = rownames(pb_cpm),
        pb_cpm,
        check.names = FALSE
    ),
    file.path(
        PB_DIR,
        "02_pseudobulk_CPM_v4.12.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Gene-level Sham vs Tx effect size
# ------------------------------------------------------------------------------

mean_sham <- rowMeans(
    pb_cpm[
        ,
        SHAM_SAMPLES,
        drop = FALSE
    ]
)

mean_tx <- rowMeans(
    pb_cpm[
        ,
        TX_SAMPLES,
        drop = FALSE
    ]
)

gene_effect <- tibble(
    gene = rownames(pb_cpm),
    mean_CPM_Sham = as.numeric(mean_sham),
    mean_CPM_Tx = as.numeric(mean_tx),
    log2FC_Tx_vs_Sham = log2(
        (mean_tx + 1) /
        (mean_sham + 1)
    )
)

write.csv(
    gene_effect,
    file.path(
        PB_DIR,
        "03_gene_effect_size_Tx_vs_Sham_v4.12.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Program gene table
# ------------------------------------------------------------------------------

PROGRAM_USE <- lapply(
    PROGRAMS,
    intersect,
    y = rownames(pb_cpm)
)

program_gene_table <- bind_rows(
    lapply(
        names(PROGRAM_USE),
        function(program_name) {

            tibble(
                program = program_name,
                gene = PROGRAM_USE[[program_name]]
            )
        }
    )
)

program_gene_effect <- program_gene_table %>%
    left_join(
        gene_effect,
        by = "gene"
    )

write.csv(
    program_gene_effect,
    file.path(
        TAB_DIR,
        "04_program_gene_effect_sizes_v4.12.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Program-level gene-by-gene effect plot
# ------------------------------------------------------------------------------

p_program_genes <- ggplot(
    program_gene_effect,
    aes(
        x = log2FC_Tx_vs_Sham,
        y = reorder(
            gene,
            log2FC_Tx_vs_Sham
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
            xend = log2FC_Tx_vs_Sham,
            yend = reorder(
                gene,
                log2FC_Tx_vs_Sham
            )
        ),
        linewidth = 0.5
    ) +
    geom_point(
        size = 2.1
    ) +
    facet_wrap(
        ~ program,
        scales = "free_y",
        ncol = 3
    ) +
    labs(
        title =
            "Anti-inflammatory-MΦ program decomposition: Sham → Tx",
        subtitle =
            "Gene-level SAFE pseudobulk effect sizes",
        x =
            "log2FC (Tx / Sham)",
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
    "01_program_gene_effect_sizes_v4.12.1.pdf",
    p_program_genes,
    13,
    13
)

# ------------------------------------------------------------------------------
# 10. Four-sample gene-level table
# ------------------------------------------------------------------------------

all_program_genes <- unique(
    program_gene_table$gene
)

sample_gene_long <- as.data.frame(
    pb_log2cpm[
        all_program_genes,
        ,
        drop = FALSE
    ]
) %>%
    rownames_to_column("gene") %>%
    pivot_longer(
        cols = -gene,
        names_to = "sample",
        values_to = "log2CPM"
    ) %>%
    mutate(
        condition = ifelse(
            sample %in% SHAM_SAMPLES,
            "Sham",
            "Tx"
        )
    ) %>%
    left_join(
        program_gene_table,
        by = "gene"
    )

write.csv(
    sample_gene_long,
    file.path(
        TAB_DIR,
        "05_program_gene_sample_level_log2CPM_v4.12.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. Replicate consistency analysis
#
# We compare:
#   Tx17 vs mean Sham
#   Tx5  vs mean Sham
#
# A gene is "direction-consistent" if both Tx samples move in the same
# direction relative to the Sham mean.
# ------------------------------------------------------------------------------

sham_mean_log <- rowMeans(
    pb_log2cpm[
        ,
        SHAM_SAMPLES,
        drop = FALSE
    ]
)

tx17_delta <- pb_log2cpm[, "Tx17"] - sham_mean_log
tx5_delta  <- pb_log2cpm[, "Tx5"]  - sham_mean_log

consistency <- tibble(
    gene = rownames(pb_log2cpm),
    Tx17_minus_ShamMean =
        as.numeric(tx17_delta),
    Tx5_minus_ShamMean =
        as.numeric(tx5_delta),

    direction_consistent =
        sign(Tx17_minus_ShamMean) ==
        sign(Tx5_minus_ShamMean),

    consistent_direction = case_when(
        direction_consistent &
            Tx17_minus_ShamMean > 0 ~
            "Up in both Tx",

        direction_consistent &
            Tx17_minus_ShamMean < 0 ~
            "Down in both Tx",

        TRUE ~
            "Discordant"
    )
)

program_consistency <- program_gene_table %>%
    left_join(
        consistency,
        by = "gene"
    ) %>%
    left_join(
        gene_effect %>%
            select(
                gene,
                log2FC_Tx_vs_Sham
            ),
        by = "gene"
    )

write.csv(
    program_consistency,
    file.path(
        TAB_DIR,
        "06_program_gene_replicate_consistency_v4.12.1.csv"
    ),
    row.names = FALSE
)

consistency_summary <- program_consistency %>%
    group_by(program) %>%
    summarise(
        n_genes = n(),
        n_up_both =
            sum(
                consistent_direction ==
                    "Up in both Tx",
                na.rm = TRUE
            ),
        n_down_both =
            sum(
                consistent_direction ==
                    "Down in both Tx",
                na.rm = TRUE
            ),
        n_discordant =
            sum(
                consistent_direction ==
                    "Discordant",
                na.rm = TRUE
            ),
        pct_consistent =
            100 *
            (
                n_up_both +
                n_down_both
            ) /
            n_genes,
        .groups = "drop"
    )

write.csv(
    consistency_summary,
    file.path(
        TAB_DIR,
        "07_program_replicate_consistency_summary_v4.12.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. Consistency plot
# ------------------------------------------------------------------------------

p_consistency <- ggplot(
    program_consistency,
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
        size = 1.9
    ) +
    facet_wrap(
        ~ program,
        scales = "free",
        ncol = 3
    ) +
    labs(
        title =
            "Biological-replicate consistency of Anti-inflammatory-MΦ gene changes",
        subtitle =
            "Each axis = Tx sample minus mean Sham expression",
        x =
            "Tx17 − mean Sham log2CPM",
        y =
            "Tx5 − mean Sham log2CPM"
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
    "02_program_gene_replicate_consistency_v4.12.1.pdf",
    p_consistency,
    12,
    11
)

# ------------------------------------------------------------------------------
# 13. HSC-relevant ligand screen
# ------------------------------------------------------------------------------

HSC_LIGAND_USE <- lapply(
    HSC_LIGANDS,
    intersect,
    y = rownames(pb_cpm)
)

hsc_ligand_membership <- bind_rows(
    lapply(
        names(HSC_LIGAND_USE),
        function(group_name) {

            tibble(
                ligand_group = group_name,
                gene = HSC_LIGAND_USE[[group_name]]
            )
        }
    )
)

hsc_ligand_effect <- hsc_ligand_membership %>%
    left_join(
        gene_effect,
        by = "gene"
    ) %>%
    left_join(
        consistency,
        by = "gene"
    )

write.csv(
    hsc_ligand_effect,
    file.path(
        TAB_DIR,
        "08_HSC_relevant_Mphi_ligand_effects_v4.12.1.csv"
    ),
    row.names = FALSE
)

p_hsc_ligand <- ggplot(
    hsc_ligand_effect,
    aes(
        x = log2FC_Tx_vs_Sham,
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
            xend = log2FC_Tx_vs_Sham,
            yend = gene
        ),
        linewidth = 0.55
    ) +
    geom_point(
        aes(
            shape = consistent_direction
        ),
        size = 2.8
    ) +
    facet_wrap(
        ~ ligand_group,
        scales = "free_y",
        ncol = 1
    ) +
    labs(
        title =
            "Anti-inflammatory-MΦ → HSC candidate ligand shift",
        subtitle =
            "Clean-B SAFE pseudobulk | shape = replicate consistency",
        x =
            "log2FC (Tx / Sham)",
        y =
            NULL,
        shape =
            "Tx replicate\nconsistency"
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
    "03_HSC_relevant_Mphi_ligand_shift_v4.12.1.pdf",
    p_hsc_ligand,
    9,
    9
)

# ------------------------------------------------------------------------------
# 14. Key mechanistic gene matrix
# ------------------------------------------------------------------------------

MECHANISTIC_GENES <- unique(
    c(
        "Il1b",
        "Tnf",
        "Ccl2",
        "Ccl3",
        "Ccl4",
        "Cxcl10",
        "Mrc1",
        "Cd163",
        "Il1rn",
        "Mertk",
        "Igf1",
        "Hmox1",
        "Maf",
        "Il10ra",
        "Il10rb",
        "Stat3",
        "Socs3",
        "Mfge8",
        "Gas6",
        "Spp1",
        "Tgfb1",
        "Pdgfb",
        "Trem2",
        "Gpnmb",
        "Apoe"
    )
)

MECHANISTIC_USE <- intersect(
    MECHANISTIC_GENES,
    rownames(pb_log2cpm)
)

mech_z <- pb_log2cpm[
    MECHANISTIC_USE,
    SAMPLE_ORDER,
    drop = FALSE
]

mech_z <- t(
    scale(
        t(
            mech_z
        )
    )
)

mech_z[
    !is.finite(
        mech_z
    )
] <- 0

mech_df <- as.data.frame(
    mech_z
) %>%
    rownames_to_column(
        "gene"
    ) %>%
    pivot_longer(
        cols = -gene,
        names_to = "sample",
        values_to = "z"
    ) %>%
    mutate(
        sample = factor(
            sample,
            levels = SAMPLE_ORDER
        ),
        gene = factor(
            gene,
            levels = rev(
                MECHANISTIC_USE
            )
        )
    )

p_mech_heat <- ggplot(
    mech_df,
    aes(
        x = sample,
        y = gene,
        fill = z
    )
) +
    geom_tile() +
    scale_fill_gradient2(
        low = "#0033FF",
        mid = "#FFFFFF",
        high = "#FF1A1A",
        midpoint = 0,
        limits = c(
            -2,
            2
        ),
        oob = scales::squish,
        name = "Gene-wise\nz-score"
    ) +
    labs(
        title =
            "Mechanistic gene matrix: Anti-inflammatory-MΦ",
        subtitle =
            "Clean-B sample-level pseudobulk",
        x = NULL,
        y = NULL
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
    "04_mechanistic_gene_heatmap_v4.12.1.pdf",
    p_mech_heat,
    7,
    10
)

# ------------------------------------------------------------------------------
# 15. Mechanistic summary table
# ------------------------------------------------------------------------------

mechanistic_summary <- program_consistency %>%
    group_by(
        program
    ) %>%
    summarise(
        mean_log2FC =
            mean(
                log2FC_Tx_vs_Sham,
                na.rm = TRUE
            ),

        median_log2FC =
            median(
                log2FC_Tx_vs_Sham,
                na.rm = TRUE
            ),

        n_genes =
            n(),

        n_consistent_up =
            sum(
                consistent_direction ==
                    "Up in both Tx",
                na.rm = TRUE
            ),

        n_consistent_down =
            sum(
                consistent_direction ==
                    "Down in both Tx",
                na.rm = TRUE
            ),

        pct_direction_consistent =
            100 *
            mean(
                direction_consistent,
                na.rm = TRUE
            ),

        .groups =
            "drop"
    ) %>%
    arrange(
        desc(
            mean_log2FC
        )
    )

write.csv(
    mechanistic_summary,
    file.path(
        TAB_DIR,
        "09_mechanistic_program_summary_v4.12.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 16. Master figure
# ------------------------------------------------------------------------------

p_summary <- (
    p_program_genes /
    p_consistency /
    p_hsc_ligand /
    p_mech_heat
) +
    patchwork::plot_layout(
        heights = c(
            1.1,
            1.0,
            0.8,
            0.9
        )
    ) +
    patchwork::plot_annotation(
        title =
            "Clean-B Anti-inflammatory-MΦ mechanistic decomposition v4.12.1",
        subtitle =
            "Sham vs Tx | gene-level effects, replicate consistency, HSC-facing ligands",
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
    "05_mechanistic_decomposition_master_v4.12.1.pdf",
    p_summary,
    14,
    31
)

# ------------------------------------------------------------------------------
# 17. README
# ------------------------------------------------------------------------------

readme <- c(

    "Clean-B Anti-inflammatory-MΦ mechanistic decomposition v4.12.1",
    "",

    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",

    "Goal:",
    "  Determine what the Tx-expanded Anti-inflammatory-MΦ population",
    "  is functionally doing before MΦ-HSC interaction analysis.",
    "",

    "Primary analytical principles:",
    "  Biological sample is the replicate.",
    "  Sham1/Sham20 vs Tx17/Tx5.",
    "  SAFE pseudobulk = explicit Matrix::rowSums(raw counts).",
    "  Effect size + replicate consistency are primary.",
    "",

    "Programs:",
    "  Inflammatory",
    "  Anti-inflammatory",
    "  IL10/STAT3",
    "  Efferocytosis",
    "  Repair/Resolution",
    "  Fibrogenic",
    "  Lipid/TREM2",
    "",

    "HSC-facing ligand groups:",
    "  Profibrotic_HSC",
    "  Resolution_HSC",
    "  Chemotactic_HSC",
    "",

    "Primary outputs:",
    "  01 program gene effect sizes",
    "  02 biological replicate consistency",
    "  03 HSC-relevant macrophage ligand shift",
    "  04 mechanistic gene heatmap",
    "  05 master summary",
    "",

    "Next decision:",
    "  Use the HSC-facing ligand results to choose ligand-receptor pairs",
    "  for direct MΦ-HSC interaction analysis in the whole-liver dataset."
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.12.1.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.12.1.txt"
    )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)

print(mechanistic_summary)
