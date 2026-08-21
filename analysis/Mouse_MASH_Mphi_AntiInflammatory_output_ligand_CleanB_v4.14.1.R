#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Anti-inflammatory-MΦ functional-state output ligand analysis
# Clean-B / Res1.2 metastates
# v4.14.1
#
# PURPOSE
#   Identify candidate output ligands from Anti-inflammatory-MΦ functional states
#   and quantify their sample-level expression changes.
#
# FIXED FRAMEWORK
#   Parent MΦ classification = Res2.0
#   Anti-inflammatory-MΦ internal clustering = Res1.2
#   Functional metastates = v4.14.0.1
#
# PRIMARY CONTRAST
#   Sham vs Tx (n=2 vs n=2)
#
# SECONDARY DESCRIPTIVE CONTRAST
#   STD vs CDAHFD (n=1 vs n=1)
#
# OUTPUT
#   1) sample-level pseudobulk ligand CPM
#   2) Tx vs Sham effect-size tables
#   3) replicate-consistency tables
#   4) ligand heatmaps
#   5) candidate output ligands to HSC / LSEC / Hepatocyte
#
# IMPORTANT
#   This is expression-based ligand screening.
#   It is NOT yet receptor-aware cell-cell interaction inference.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4141)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "AntiInflammatory_functional_state_audit_CleanB_v4.14.0.1",
    "RDS",
    "Mouse_Mphi_AntiInflammatory_Res1.2_functional_state_annotated_v4.14.0.1.rds"
)

OUTPUT_DIR <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "AntiInflammatory_output_ligand_CleanB_v4.14.1"
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

MIN_CPM <- 1
MIN_ABS_LOG2FC <- 0.5

# ------------------------------------------------------------------------------
# 1. Candidate ligand sets
#
# These are expression-screening candidates. Receptor-aware validation comes
# later when HSC/LSEC/Hepatocyte whole-liver objects are finalized.
# ------------------------------------------------------------------------------

LIGAND_SETS <- list(

    HSC_candidate = c(
        "Tgfb1",
        "Pdgfa",
        "Pdgfb",
        "Pdgfc",
        "Igf1",
        "Hbegf",
        "Areg",
        "Ereg",
        "Tnf",
        "Il1b",
        "Il6",
        "Osm",
        "Ccl2",
        "Ccl3",
        "Ccl4",
        "Ccl7",
        "Ccl8",
        "Ccl12",
        "Cxcl2",
        "Cxcl10",
        "Spp1",
        "Gpnmb",
        "Thbs1",
        "Gas6",
        "Mfge8"
    ),

    LSEC_candidate = c(
        "Vegfa",
        "Vegfb",
        "Pgf",
        "Angpt1",
        "Angpt2",
        "Tnf",
        "Il1b",
        "Il6",
        "Osm",
        "Ccl2",
        "Ccl3",
        "Ccl4",
        "Cxcl2",
        "Cxcl10",
        "Igf1",
        "Hbegf",
        "Areg",
        "Gas6",
        "Spp1",
        "Tgfb1"
    ),

    Hepatocyte_candidate = c(
        "Il6",
        "Osm",
        "Tnf",
        "Il1b",
        "Tgfb1",
        "Igf1",
        "Hbegf",
        "Areg",
        "Ereg",
        "Gas6",
        "Mfge8",
        "Spp1",
        "Ccl2",
        "Cxcl10"
    ),

    Repair_resolution_output = c(
        "Igf1",
        "Gas6",
        "Mfge8",
        "Hbegf",
        "Areg",
        "Ereg",
        "Il1rn",
        "Tgfb1",
        "Mmp12",
        "Mmp13",
        "Mmp14",
        "Plau"
    ),

    Fibrogenic_output = c(
        "Tgfb1",
        "Pdgfa",
        "Pdgfb",
        "Pdgfc",
        "Spp1",
        "Thbs1",
        "Lgals3",
        "Gpnmb",
        "Tnf",
        "Il1b",
        "Ccl2",
        "Ccl7",
        "Ccl8",
        "Ccl12"
    ),

    Inflammatory_output = c(
        "Tnf",
        "Il1b",
        "Il6",
        "Osm",
        "Ccl2",
        "Ccl3",
        "Ccl4",
        "Ccl7",
        "Ccl8",
        "Ccl12",
        "Cxcl2",
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

safe_z <- function(x) {
    z <- as.numeric(scale(x))
    z[!is.finite(z)] <- 0
    z
}

# ------------------------------------------------------------------------------
# 4. Load annotated Res1.2 metastate object
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
    stop(
        "Input RDS not found:\n",
        INPUT_RDS
    )
}

msg("Loading: ", INPUT_RDS)

anti <- readRDS(INPUT_RDS)

if (!inherits(anti, "Seurat")) {
    stop("Input is not a Seurat object.")
}

required_meta <- c(
    "sample_v4140",
    "condition_v4140",
    "anti_subcluster_v4140",
    "anti_functional_state_v41401"
)

missing_meta <- setdiff(
    required_meta,
    colnames(anti@meta.data)
)

if (length(missing_meta) > 0L) {
    stop(
        "Missing metadata column(s): ",
        paste(missing_meta, collapse = ", ")
    )
}

# ------------------------------------------------------------------------------
# 5. Join RNA layers
# ------------------------------------------------------------------------------

DefaultAssay(anti) <- ASSAY_USE

anti <- JoinLayers(
    anti,
    assay = ASSAY_USE
)

counts <- get_layer_safe(
    anti,
    ASSAY_USE,
    "counts"
)

if (is.null(counts)) {
    stop("RNA counts layer not found.")
}

if (!identical(
    colnames(counts),
    colnames(anti)
)) {
    stop("RNA count matrix cell order mismatch.")
}

# ------------------------------------------------------------------------------
# 6. Metadata
# ------------------------------------------------------------------------------

meta <- anti@meta.data %>%
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
        ),
        functional_state = as.character(
            anti_functional_state_v41401
        )
    )

STATE_ORDER <- unique(
    as.character(
        anti$anti_functional_state_v41401
    )
)

STATE_ORDER <- STATE_ORDER[
    !is.na(STATE_ORDER)
]

# ------------------------------------------------------------------------------
# 7. Ligand-gene audit
# ------------------------------------------------------------------------------

LIGAND_USE <- lapply(
    LIGAND_SETS,
    intersect,
    y = rownames(counts)
)

ligand_audit <- bind_rows(
    lapply(
        names(LIGAND_USE),
        function(set_name) {

            tibble(
                ligand_set = set_name,
                n_requested =
                    length(
                        LIGAND_SETS[[set_name]]
                    ),
                n_detected =
                    length(
                        LIGAND_USE[[set_name]]
                    ),
                detected =
                    paste(
                        LIGAND_USE[[set_name]],
                        collapse = ";"
                    ),
                missing =
                    paste(
                        setdiff(
                            LIGAND_SETS[[set_name]],
                            LIGAND_USE[[set_name]]
                        ),
                        collapse = ";"
                    )
            )
        }
    )
)

write.csv(
    ligand_audit,
    file.path(
        TAB_DIR,
        "01_ligand_gene_audit_v4.14.1.csv"
    ),
    row.names = FALSE
)

ALL_LIGANDS <- unique(
    unlist(
        LIGAND_USE,
        use.names = FALSE
    )
)

# ------------------------------------------------------------------------------
# 8. Sample-level pseudobulk counts by functional state
# ------------------------------------------------------------------------------

pb_rows <- list()

for (state_now in STATE_ORDER) {

    for (sample_now in SAMPLE_ORDER) {

        cells_now <- meta$cell[
            meta$functional_state ==
                state_now &
            meta$sample ==
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
            state_now,
            sample_now,
            sep = "__"
        )]] <- tibble(
            functional_state = state_now,
            sample = sample_now,
            condition = as.character(
                canonical_condition(
                    sample_now
                )
            ),
            gene = rownames(counts),
            raw_count = as.numeric(pb),
            library_size = lib_size,
            CPM =
                as.numeric(pb) /
                lib_size *
                1e6,
            n_cells =
                length(
                    cells_now
                )
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
                ALL_LIGANDS
        ),
    file.path(
        TAB_DIR,
        "02_ligand_pseudobulk_CPM_by_state_sample_v4.14.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Sham vs Tx sample-level effect sizes
# ------------------------------------------------------------------------------

tx_effect_rows <- list()

for (state_now in STATE_ORDER) {

    dat_state <- pb_long %>%
        filter(
            functional_state ==
                state_now,
            gene %in%
                ALL_LIGANDS,
            sample %in%
                c(
                    "Sham1",
                    "Sham20",
                    "Tx17",
                    "Tx5"
                )
        )

    if (nrow(dat_state) == 0L) {
        next
    }

    wide <- dat_state %>%
        select(
            gene,
            sample,
            CPM
        ) %>%
        pivot_wider(
            names_from = sample,
            values_from = CPM
        )

    required_cols <- c(
        "Sham1",
        "Sham20",
        "Tx17",
        "Tx5"
    )

    if (!all(
        required_cols %in%
            colnames(wide)
    )) {
        next
    }

    res <- wide %>%
        mutate(
            Sham_mean =
                rowMeans(
                    cbind(
                        Sham1,
                        Sham20
                    ),
                    na.rm = TRUE
                ),

            Tx_mean =
                rowMeans(
                    cbind(
                        Tx17,
                        Tx5
                    ),
                    na.rm = TRUE
                ),

            log2FC_Tx_vs_Sham =
                log2(
                    (
                        Tx_mean + 1
                    ) /
                    (
                        Sham_mean + 1
                    )
                ),

            log2FC_Tx17_vs_Sham =
                log2(
                    (
                        Tx17 + 1
                    ) /
                    (
                        Sham_mean + 1
                    )
                ),

            log2FC_Tx5_vs_Sham =
                log2(
                    (
                        Tx5 + 1
                    ) /
                    (
                        Sham_mean + 1
                    )
                ),

            replicate_concordant =
                sign(
                    log2FC_Tx17_vs_Sham
                ) ==
                sign(
                    log2FC_Tx5_vs_Sham
                ),

            functional_state =
                state_now
        )

    tx_effect_rows[[state_now]] <- res
}

tx_effect <- bind_rows(
    tx_effect_rows
)

write.csv(
    tx_effect,
    file.path(
        TAB_DIR,
        "03_ligand_effect_Tx_vs_Sham_by_functional_state_v4.14.1.csv"
    ),
    row.names = FALSE
)

high_priority_tx <- tx_effect %>%
    filter(
        pmax(
            Sham_mean,
            Tx_mean,
            na.rm = TRUE
        ) >=
            MIN_CPM,
        abs(
            log2FC_Tx_vs_Sham
        ) >=
            MIN_ABS_LOG2FC,
        replicate_concordant
    ) %>%
    arrange(
        functional_state,
        desc(
            abs(
                log2FC_Tx_vs_Sham
            )
        )
    )

write.csv(
    high_priority_tx,
    file.path(
        TAB_DIR,
        "04_high_priority_concordant_ligands_Tx_vs_Sham_v4.14.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. STD vs CDAHFD descriptive effect
# ------------------------------------------------------------------------------

disease_effect_rows <- list()

for (state_now in STATE_ORDER) {

    dat_state <- pb_long %>%
        filter(
            functional_state ==
                state_now,
            gene %in%
                ALL_LIGANDS,
            sample %in%
                c(
                    "STD_rep1",
                    "CDHFD_rep1"
                )
        )

    wide <- dat_state %>%
        select(
            gene,
            sample,
            CPM
        ) %>%
        pivot_wider(
            names_from = sample,
            values_from = CPM
        )

    if (!all(
        c(
            "STD_rep1",
            "CDHFD_rep1"
        ) %in%
            colnames(wide)
    )) {
        next
    }

    disease_effect_rows[[state_now]] <- wide %>%
        mutate(
            log2FC_CDAHFD_vs_STD =
                log2(
                    (
                        CDHFD_rep1 + 1
                    ) /
                    (
                        STD_rep1 + 1
                    )
                ),
            functional_state =
                state_now
        )
}

disease_effect <- bind_rows(
    disease_effect_rows
)

write.csv(
    disease_effect,
    file.path(
        TAB_DIR,
        "05_ligand_effect_CDAHFD_vs_STD_descriptive_v4.14.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. Ligand membership table
# ------------------------------------------------------------------------------

ligand_membership <- bind_rows(
    lapply(
        names(LIGAND_USE),
        function(set_name) {

            tibble(
                ligand_set = set_name,
                gene = LIGAND_USE[[set_name]]
            )
        }
    )
)

write.csv(
    ligand_membership,
    file.path(
        TAB_DIR,
        "06_ligand_membership_v4.14.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 12. Tx vs Sham effect-size plot
# ------------------------------------------------------------------------------

plot_tx <- tx_effect %>%
    filter(
        pmax(
            Sham_mean,
            Tx_mean,
            na.rm = TRUE
        ) >=
            MIN_CPM
    ) %>%
    mutate(
        direction =
            ifelse(
                log2FC_Tx_vs_Sham > 0,
                "Up in Tx",
                "Down in Tx"
            )
    )

p_effect <- ggplot(
    plot_tx,
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
    geom_point(
        aes(
            shape = replicate_concordant
        ),
        size = 2
    ) +
    facet_wrap(
        ~ functional_state,
        scales = "free_y",
        ncol = 2
    ) +
    labs(
        title =
            "Anti-inflammatory-MΦ output ligand remodeling",
        subtitle =
            "Sham → Tx | sample-level pseudobulk CPM",
        x =
            "log2FC (Tx mean / Sham mean)",
        y =
            NULL,
        shape =
            "Tx17/Tx5\nconcordant"
    ) +
    theme_classic(
        base_size = 9
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        axis.text.y = element_text(
            size = 6.5
        ),
        strip.text = element_text(
            face = "bold",
            size = 7.5
        )
    )

save_pdf(
    "01_output_ligand_effect_Tx_vs_Sham_v4.14.1.pdf",
    p_effect,
    11,
    13
)

# ------------------------------------------------------------------------------
# 13. Replicate-consistency plot
# ------------------------------------------------------------------------------

p_consistency <- ggplot(
    tx_effect %>%
        filter(
            pmax(
                Sham_mean,
                Tx_mean,
                na.rm = TRUE
            ) >=
                MIN_CPM
        ),
    aes(
        x = log2FC_Tx17_vs_Sham,
        y = log2FC_Tx5_vs_Sham
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
        size = 1.8
    ) +
    facet_wrap(
        ~ functional_state,
        ncol = 3
    ) +
    labs(
        title =
            "Tx replicate consistency of output-ligand changes",
        x =
            "log2FC Tx17 vs mean Sham",
        y =
            "log2FC Tx5 vs mean Sham"
    ) +
    theme_classic(
        base_size = 9
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        strip.text = element_text(
            face = "bold",
            size = 7
        )
    )

save_pdf(
    "02_output_ligand_Tx_replicate_consistency_v4.14.1.pdf",
    p_consistency,
    11,
    8
)

# ------------------------------------------------------------------------------
# 14. Sample-level ligand heatmap
# ------------------------------------------------------------------------------

heat_source <- pb_long %>%
    filter(
        gene %in%
            ALL_LIGANDS
    ) %>%
    mutate(
        column_id = paste(
            functional_state,
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

heat_mat <- as.matrix(
    heat_source[
        ,
        -1,
        drop = FALSE
    ]
)

rownames(
    heat_mat
) <- heat_source$gene

heat_log <- log2(
    heat_mat + 1
)

heat_z <- t(
    scale(
        t(
            heat_log
        )
    )
)

heat_z[
    !is.finite(
        heat_z
    )
] <- 0

heat_plot <- pmax(
    pmin(
        heat_z,
        2
    ),
    -2
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
        "07_ligand_sample_level_zscore_matrix_v4.14.1.csv"
    ),
    row.names = FALSE
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "03_output_ligand_sample_level_heatmap_v4.14.1.pdf"
    ),
    width = 18,
    height = 11
)

pheatmap::pheatmap(
    heat_plot,
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
    fontsize_col = 5.5,
    angle_col = 45,
    main = paste0(
        "Anti-inflammatory-MΦ output ligand expression\n",
        "sample-level pseudobulk log2(CPM+1) | row z-score"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 15. Candidate-target heatmaps
# ------------------------------------------------------------------------------

for (set_name in c(
    "HSC_candidate",
    "LSEC_candidate",
    "Hepatocyte_candidate"
)) {

    genes_now <- intersect(
        LIGAND_USE[[set_name]],
        rownames(
            heat_plot
        )
    )

    if (length(genes_now) == 0L) {
        next
    }

    mat_now <- heat_plot[
        genes_now,
        ,
        drop = FALSE
    ]

    outfile <- paste0(
        "04_",
        set_name,
        "_ligand_heatmap_v4.14.1.pdf"
    )

    grDevices::cairo_pdf(
        file.path(
            FIG_DIR,
            outfile
        ),
        width = 16,
        height = 8
    )

    pheatmap::pheatmap(
        mat_now,
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
        fontsize_col = 5.5,
        angle_col = 45,
        main = paste0(
            set_name,
            " | Anti-inflammatory-MΦ output ligands"
        )
    )

    grDevices::dev.off()
}

# ------------------------------------------------------------------------------
# 16. High-priority target-specific tables
# ------------------------------------------------------------------------------

high_priority_target <- high_priority_tx %>%
    left_join(
        ligand_membership,
        by = "gene"
    ) %>%
    arrange(
        ligand_set,
        functional_state,
        desc(
            abs(
                log2FC_Tx_vs_Sham
            )
        )
    )

write.csv(
    high_priority_target,
    file.path(
        TAB_DIR,
        "08_high_priority_target_specific_ligands_v4.14.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 17. README
# ------------------------------------------------------------------------------

readme <- c(

    "Mouse MASH Anti-inflammatory-MΦ output ligand analysis v4.14.1",
    "",

    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",

    "Fixed framework:",
    "  Parent MΦ classification = Res2.0",
    "  Anti-inflammatory-MΦ internal clustering = Res1.2",
    "  Functional metastates = v4.14.0.1",
    "",

    "Primary comparison:",
    "  Sham vs Tx",
    "  biological samples: Sham1, Sham20, Tx17, Tx5",
    "",

    "Secondary descriptive comparison:",
    "  STD vs CDAHFD",
    "  n=1 vs n=1",
    "",

    "Primary analysis:",
    "  sample-level pseudobulk raw counts",
    "  CPM normalization",
    "  Tx mean vs Sham mean log2FC",
    "  Tx17 / Tx5 direction concordance",
    "",

    paste0(
        "High-priority criteria: CPM >= ",
        MIN_CPM,
        ", |log2FC| >= ",
        MIN_ABS_LOG2FC,
        ", Tx17/Tx5 concordant"
    ),
    "",

    "Candidate target screens:",
    "  HSC",
    "  LSEC",
    "  Hepatocyte",
    "",

    "Important limitation:",
    "  This is ligand-expression screening.",
    "  Receptor-aware cell-cell interaction inference is NOT performed here.",
    "",

    "Primary outputs:",
    "  01 ligand effect-size plot",
    "  02 Tx replicate consistency",
    "  03 all-ligand sample-level heatmap",
    "  04 target-specific ligand heatmaps",
    "",

    "Key tables:",
    "  03 ligand effect Tx vs Sham",
    "  04 high-priority concordant ligands",
    "  08 target-specific high-priority ligands",
    "",

    "Next:",
    "  review candidate ligands and proceed to receptor-aware HSC/LSEC/Hepatocyte interaction analysis after whole-liver Res2.0 annotation is finalized"
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.14.1.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.14.1.txt"
    )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)

print(
    high_priority_tx
)
