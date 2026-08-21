#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Subtype-specific pathway enrichment
# Clean-B
# v4.14.2
#
# PRIMARY CONTRAST
#   Sham vs Tx, biological samples:
#     Sham1, Sham20, Tx17, Tx5
#
# SECONDARY DESCRIPTIVE CONTRAST
#   STD_rep1 vs CDHFD_rep1
#
# RANKING
#   sample-level pseudobulk log2(CPM + 1)
#   rank = mean(Tx) - mean(Sham)
#
# ENRICHMENT
#   Hallmark
#   Reactome
#   KEGG
#
# IMPORTANT
#   - Biological sample is the unit used to build ranks.
#   - With n=2 vs n=2, pathway results are exploratory/descriptive.
#   - No cell-level p-values are used for treatment inference.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4142)

# ------------------------------------------------------------------------------
# 0. Paths / parameters
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
    "Mphi_subtype_pathway_enrichment_CleanB_v4.14.2"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"

SUBTYPE_ORDER <- c(
    "Inflammatory-Mphi",
    "Anti-inflammatory-Mphi",
    "Fibrogenic-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi"
)

SUBTYPE_LABELS <- c(
    "Inflammatory-Mphi"           = "Inflammatory-MΦ",
    "Anti-inflammatory-Mphi"      = "Anti-inflammatory-MΦ",
    "Fibrogenic-Mphi"             = "Fibrogenic-MΦ",
    "Repair/Resolution-Mphi"      = "Repair/Resolution-MΦ",
    "Lipid-associated/TREM2-Mphi" = "Lipid-associated/TREM2-MΦ"
)

SAMPLE_ORDER <- c(
    "STD_rep1",
    "CDHFD_rep1",
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
)

MIN_CPM_FOR_RANK <- 0.5
FGSEA_MIN_SIZE <- 10L
FGSEA_MAX_SIZE <- 500L
FGSEA_EPS <- 0

TOP_N_PATHWAYS_PER_SUBTYPE <- 10L

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
    "fgsea",
    "msigdbr"
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
        paste(missing_packages, collapse = ", "),
        "\nInstall before running v4.14.2.",
        "\nBioconductor package fgsea may require BiocManager::install('fgsea')."
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
    library(fgsea)
    library(msigdbr)
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

    out
}

safe_pathway_name <- function(x) {
    x <- gsub("^HALLMARK_", "", x)
    x <- gsub("^REACTOME_", "", x)
    x <- gsub("^KEGG_", "", x)
    x <- gsub("_", " ", x)
    tools::toTitleCase(tolower(x))
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


flatten_fgsea_for_csv <- function(df) {

    out <- df

    if ("leadingEdge" %in%
        colnames(out)) {

        out$leadingEdge <- vapply(
            out$leadingEdge,
            function(x) {

                if (length(x) == 0L) {
                    return("")
                }

                paste(
                    x,
                    collapse = ";"
                )
            },
            character(1)
        )
    }

    out
}

# ------------------------------------------------------------------------------
# 3. Load Clean-B object
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

if (is.na(SAMPLE_COL)) stop("Sample metadata column not found.")
if (is.na(CLASS_COL)) stop("MΦ class metadata column not found.")

mphi <- JoinLayers(
    mphi,
    assay = ASSAY_USE
)

counts <- get_layer_safe(
    mphi,
    ASSAY_USE,
    "counts"
)

if (is.null(counts)) {
    stop("RNA counts layer not found.")
}

# ------------------------------------------------------------------------------
# 4. Canonical metadata
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
        )
    ) %>%
    filter(
        macrophage_class %in%
            SUBTYPE_ORDER
    )

# ------------------------------------------------------------------------------
# 5. Build sample-level pseudobulk matrix
# ------------------------------------------------------------------------------

pb_rows <- list()

for (subtype_now in SUBTYPE_ORDER) {

    for (sample_now in SAMPLE_ORDER) {

        cells_now <- meta$cell[
            meta$macrophage_class ==
                subtype_now &
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
            subtype_now,
            sample_now,
            sep = "__"
        )]] <- tibble(
            macrophage_class = subtype_now,
            sample = sample_now,
            condition = canonical_condition(sample_now),
            gene = rownames(counts),
            raw_count = as.numeric(pb),
            library_size = lib_size,
            CPM = as.numeric(pb) / lib_size * 1e6,
            logCPM = log2(
                as.numeric(pb) / lib_size * 1e6 + 1
            ),
            n_cells = length(cells_now)
        )
    }
}

pb_long <- bind_rows(pb_rows)

write.csv(
    pb_long,
    file.path(
        TAB_DIR,
        "01_sample_level_pseudobulk_all_genes_v4.14.2.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Build Sham vs Tx ranks
# ------------------------------------------------------------------------------

rank_tables <- list()

for (subtype_now in SUBTYPE_ORDER) {

    dat <- pb_long %>%
        filter(
            macrophage_class == subtype_now,
            sample %in%
                c(
                    "Sham1",
                    "Sham20",
                    "Tx17",
                    "Tx5"
                )
        ) %>%
        select(
            gene,
            sample,
            CPM,
            logCPM
        )

    wide_cpm <- dat %>%
        select(
            gene,
            sample,
            CPM
        ) %>%
        pivot_wider(
            names_from = sample,
            values_from = CPM
        )

    wide_log <- dat %>%
        select(
            gene,
            sample,
            logCPM
        ) %>%
        pivot_wider(
            names_from = sample,
            values_from = logCPM
        )

    required_samples <- c(
        "Sham1",
        "Sham20",
        "Tx17",
        "Tx5"
    )

    if (!all(required_samples %in%
        colnames(wide_log))) {
        stop(
            "Missing biological sample(s) in subtype ",
            subtype_now
        )
    }

    rank_df <- wide_log %>%
        mutate(
            Sham_mean_logCPM =
                rowMeans(
                    cbind(
                        Sham1,
                        Sham20
                    ),
                    na.rm = TRUE
                ),

            Tx_mean_logCPM =
                rowMeans(
                    cbind(
                        Tx17,
                        Tx5
                    ),
                    na.rm = TRUE
                ),

            rank_Tx_vs_Sham =
                Tx_mean_logCPM -
                Sham_mean_logCPM,

            rank_Tx17_vs_Sham =
                Tx17 -
                Sham_mean_logCPM,

            rank_Tx5_vs_Sham =
                Tx5 -
                Sham_mean_logCPM
        ) %>%
        left_join(
            wide_cpm %>%
                mutate(
                    max_CPM =
                        pmax(
                            Sham1,
                            Sham20,
                            Tx17,
                            Tx5,
                            na.rm = TRUE
                        )
                ) %>%
                select(
                    gene,
                    max_CPM
                ),
            by = "gene"
        ) %>%
        mutate(
            macrophage_class =
                subtype_now,
            replicate_concordant =
                sign(
                    rank_Tx17_vs_Sham
                ) ==
                sign(
                    rank_Tx5_vs_Sham
                )
        ) %>%
        arrange(
            desc(
                rank_Tx_vs_Sham
            )
        )

    rank_tables[[subtype_now]] <- rank_df
}

rank_all <- bind_rows(rank_tables)

write.csv(
    rank_all,
    file.path(
        TAB_DIR,
        "02_gene_rank_Tx_vs_Sham_by_subtype_v4.14.2.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Build disease descriptive ranks
# ------------------------------------------------------------------------------

disease_rank_tables <- list()

for (subtype_now in SUBTYPE_ORDER) {

    dat <- pb_long %>%
        filter(
            macrophage_class ==
                subtype_now,
            sample %in%
                c(
                    "STD_rep1",
                    "CDHFD_rep1"
                )
        ) %>%
        select(
            gene,
            sample,
            logCPM
        ) %>%
        pivot_wider(
            names_from = sample,
            values_from = logCPM
        )

    if (!all(
        c(
            "STD_rep1",
            "CDHFD_rep1"
        ) %in%
            colnames(dat)
    )) {
        next
    }

    disease_rank_tables[[subtype_now]] <- dat %>%
        mutate(
            rank_CDAHFD_vs_STD =
                CDHFD_rep1 -
                STD_rep1,
            macrophage_class =
                subtype_now
        )
}

disease_rank_all <- bind_rows(
    disease_rank_tables
)

write.csv(
    disease_rank_all,
    file.path(
        TAB_DIR,
        "03_gene_rank_CDAHFD_vs_STD_descriptive_v4.14.2.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. MSigDB gene sets
# ------------------------------------------------------------------------------

msg("Loading MSigDB mouse gene sets...")

msig_all <- msigdbr::msigdbr(
    species = "Mus musculus"
)

required_msig_cols <- c(
    "gs_name",
    "gene_symbol"
)

if (!all(
    required_msig_cols %in%
        colnames(msig_all)
)) {
    stop(
        "Unexpected msigdbr schema. Required columns absent."
    )
}

collection_col <- if (
    "gs_collection" %in%
        colnames(msig_all)
) {
    "gs_collection"
} else if (
    "gs_cat" %in%
        colnames(msig_all)
) {
    "gs_cat"
} else {
    NA_character_
}

subcollection_col <- if (
    "gs_subcollection" %in%
        colnames(msig_all)
) {
    "gs_subcollection"
} else if (
    "gs_subcat" %in%
        colnames(msig_all)
) {
    "gs_subcat"
} else {
    NA_character_
}

if (is.na(collection_col)) {
    stop(
        "Could not identify MSigDB collection column."
    )
}

hallmark_df <- msig_all %>%
    filter(
        .data[[collection_col]] ==
            "H"
    )

reactome_df <- msig_all %>%
    filter(
        grepl(
            "REACTOME",
            gs_name,
            ignore.case = TRUE
        ) |
            (
                !is.na(subcollection_col) &
                grepl(
                    "REACTOME",
                    .data[[subcollection_col]],
                    ignore.case = TRUE
                )
            )
    )

kegg_df <- msig_all %>%
    filter(
        grepl(
            "KEGG",
            gs_name,
            ignore.case = TRUE
        ) |
            (
                !is.na(subcollection_col) &
                grepl(
                    "KEGG",
                    .data[[subcollection_col]],
                    ignore.case = TRUE
                )
            )
    )

make_pathways <- function(df) {
    split(
        df$gene_symbol,
        df$gs_name
    )
}

PATHWAY_COLLECTIONS <- list(
    Hallmark = make_pathways(
        hallmark_df
    ),
    Reactome = make_pathways(
        reactome_df
    ),
    KEGG = make_pathways(
        kegg_df
    )
)

pathway_audit <- tibble(
    collection = names(
        PATHWAY_COLLECTIONS
    ),
    n_pathways = vapply(
        PATHWAY_COLLECTIONS,
        length,
        integer(1)
    )
)

write.csv(
    pathway_audit,
    file.path(
        TAB_DIR,
        "04_pathway_collection_audit_v4.14.2.csv"
    ),
    row.names = FALSE
)

print(pathway_audit)

# ------------------------------------------------------------------------------
# 9. Run fgsea for Sham vs Tx
# ------------------------------------------------------------------------------

fgsea_results <- list()

for (subtype_now in SUBTYPE_ORDER) {

    rank_df <- rank_tables[[subtype_now]] %>%
        filter(
            max_CPM >=
                MIN_CPM_FOR_RANK
        )

    stats <- rank_df$rank_Tx_vs_Sham
    names(stats) <- rank_df$gene

    stats <- stats[
        !is.na(stats)
    ]

    stats <- sort(
        stats,
        decreasing = TRUE
    )

    stats <- stats[
        !duplicated(
            names(stats)
        )
    ]

    for (collection_now in names(
        PATHWAY_COLLECTIONS
    )) {

        pathways_now <-
            PATHWAY_COLLECTIONS[[collection_now]]

        if (length(pathways_now) == 0L) {
            next
        }

        res <- suppressWarnings(
            fgsea::fgseaMultilevel(
                pathways = pathways_now,
                stats = stats,
                minSize = FGSEA_MIN_SIZE,
                maxSize = FGSEA_MAX_SIZE,
                eps = FGSEA_EPS
            )
        )

        res <- as_tibble(res) %>%
            mutate(
                macrophage_class =
                    subtype_now,
                collection =
                    collection_now,
                pathway_label =
                    safe_pathway_name(
                        pathway
                    )
            )

        fgsea_results[[paste(
            subtype_now,
            collection_now,
            sep = "__"
        )]] <- res
    }
}

fgsea_all <- bind_rows(
    fgsea_results
)

fgsea_all_export <- flatten_fgsea_for_csv(
    fgsea_all
)

write.csv(
    fgsea_all_export,
    file.path(
        TAB_DIR,
        "05_fgsea_Tx_vs_Sham_all_subtypes_v4.14.2.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Top pathways per subtype
# ------------------------------------------------------------------------------

top_pathways <- fgsea_all %>%
    filter(
        is.finite(
            NES
        )
    ) %>%
    group_by(
        macrophage_class
    ) %>%
    arrange(
        padj,
        desc(
            abs(
                NES
            )
        ),
        .by_group = TRUE
    ) %>%
    slice_head(
        n = TOP_N_PATHWAYS_PER_SUBTYPE
    ) %>%
    ungroup()

write.csv(
    flatten_fgsea_for_csv(
        top_pathways
    ),
    file.path(
        TAB_DIR,
        "06_top_pathways_Tx_vs_Sham_v4.14.2.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 11. NES heatmap
# ------------------------------------------------------------------------------

heat_source <- top_pathways %>%
    distinct(
        collection,
        pathway,
        pathway_label
    )

heat_pathways <- unique(
    top_pathways$pathway
)

nes_df <- fgsea_all %>%
    filter(
        pathway %in%
            heat_pathways
    ) %>%
    select(
        pathway,
        pathway_label,
        macrophage_class,
        NES
    ) %>%
    distinct() %>%
    pivot_wider(
        names_from =
            macrophage_class,
        values_from =
            NES
    )

if (nrow(nes_df) > 0L) {

    nes_mat <- as.matrix(
        nes_df[
            ,
            intersect(
                SUBTYPE_ORDER,
                colnames(
                    nes_df
                )
            ),
            drop = FALSE
        ]
    )

    rownames(
        nes_mat
    ) <- make.unique(
        nes_df$pathway_label
    )

    nes_mat[
        !is.finite(
            nes_mat
        )
    ] <- 0

    NES_LIMIT <- max(
        2,
        quantile(
            abs(
                nes_mat
            ),
            probs = 0.95,
            na.rm = TRUE
        )
    )

    nes_plot <- pmax(
        pmin(
            nes_mat,
            NES_LIMIT
        ),
        -NES_LIMIT
    )

    grDevices::cairo_pdf(
        file.path(
            FIG_DIR,
            "01_pathway_NES_heatmap_Tx_vs_Sham_v4.14.2.pdf"
        ),
        width = 9,
        height = 14
    )

    pheatmap::pheatmap(
        nes_plot,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = grDevices::colorRampPalette(
            c(
                "#0033FF",
                "#FFFFFF",
                "#FF1A1A"
            )
        )(101),
        breaks = seq(
            -NES_LIMIT,
            NES_LIMIT,
            length.out = 102
        ),
        border_color = "white",
        fontsize_row = 6.5,
        fontsize_col = 8,
        angle_col = 45,
        main = paste0(
            "MΦ subtype pathway remodeling: Sham → Tx\n",
            "NES > 0 = enriched in Tx | sample-level pseudobulk rank"
        )
    )

    grDevices::dev.off()
}

# ------------------------------------------------------------------------------
# 12. Hallmark overview plot
# ------------------------------------------------------------------------------

hallmark_plot_df <- fgsea_all %>%
    filter(
        collection == "Hallmark",
        is.finite(
            NES
        )
    ) %>%
    mutate(
        macrophage_class =
            factor(
                macrophage_class,
                levels =
                    SUBTYPE_ORDER
            ),
        significant =
            !is.na(padj) &
            padj < 0.10
    )

p_hallmark <- ggplot(
    hallmark_plot_df,
    aes(
        x = NES,
        y = reorder(
            pathway_label,
            NES
        )
    )
) +
    geom_vline(
        xintercept = 0,
        linetype = 2,
        linewidth = 0.35
    ) +
    geom_point(
        aes(
            size = -log10(
                pmax(
                    padj,
                    1e-300
                )
            ),
            shape = significant
        ),
        alpha = 0.85
    ) +
    facet_wrap(
        ~ macrophage_class,
        scales = "free_y",
        ncol = 2,
        labeller =
            as_labeller(
                SUBTYPE_LABELS
            )
    ) +
    labs(
        title =
            "Hallmark pathway remodeling by MΦ subtype",
        subtitle =
            "Sham → Tx | NES > 0 indicates enrichment in Tx",
        x =
            "Normalized enrichment score (NES)",
        y =
            NULL,
        size =
            "-log10(adj. P)",
        shape =
            "FDR < 0.10"
    ) +
    theme_classic(
        base_size = 9
    ) +
    theme(
        plot.title =
            element_text(
                face = "bold"
            ),
        axis.text.y =
            element_text(
                size = 6
            ),
        strip.text =
            element_text(
                face = "bold"
            )
    )

save_pdf(
    "02_Hallmark_pathway_overview_Tx_vs_Sham_v4.14.2.pdf",
    p_hallmark,
    12,
    16
)

# ------------------------------------------------------------------------------
# 13. Curated mechanism screen
# ------------------------------------------------------------------------------

mechanism_terms <- c(
    "IL6",
    "JAK",
    "STAT3",
    "TNFA",
    "NFKB",
    "INFLAMMATORY",
    "INTERFERON",
    "TGF",
    "COLLAGEN",
    "EXTRACELLULAR_MATRIX",
    "ECM",
    "FATTY_ACID",
    "CHOLESTEROL",
    "OXIDATIVE",
    "PEROXISOME",
    "APOPTOSIS",
    "PHAGOCYT",
    "ENDOCYT",
    "LYSOSOM",
    "ANGIOGEN",
    "WOUND",
    "REPAIR"
)

mechanism_regex <- paste(
    mechanism_terms,
    collapse = "|"
)

mechanism_df <- fgsea_all %>%
    filter(
        grepl(
            mechanism_regex,
            pathway,
            ignore.case = TRUE
        )
    ) %>%
    arrange(
        macrophage_class,
        padj,
        desc(
            abs(
                NES
            )
        )
    )

write.csv(
    flatten_fgsea_for_csv(
        mechanism_df
    ),
    file.path(
        TAB_DIR,
        "07_mechanism_focused_pathways_v4.14.2.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 14. README
# ------------------------------------------------------------------------------

readme <- c(
    "Mouse MASH MΦ subtype-specific pathway enrichment v4.14.2",
    "",
    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",
    "Primary contrast:",
    "  Sham1 + Sham20 vs Tx17 + Tx5",
    "",
    "Ranking:",
    "  biological-sample pseudobulk log2(CPM+1)",
    "  rank = mean(Tx) - mean(Sham)",
    "",
    "Collections:",
    "  Hallmark",
    "  Reactome",
    "  KEGG",
    "",
    "Interpretation:",
    "  NES > 0 = pathway enriched in Tx",
    "  NES < 0 = pathway enriched in Sham",
    "",
    "Important limitation:",
    "  n=2 vs n=2.",
    "  GSEA is exploratory/descriptive.",
    "  Do not interpret pathway FDR as definitive treatment validation.",
    "",
    "Primary outputs:",
    "  01 pathway NES heatmap",
    "  02 Hallmark overview",
    "",
    "Key tables:",
    "  05 all fgsea results",
    "  06 top pathways",
    "  07 mechanism-focused pathways"
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.14.2.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.14.2.txt"
    )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)

print(
    top_pathways
)
