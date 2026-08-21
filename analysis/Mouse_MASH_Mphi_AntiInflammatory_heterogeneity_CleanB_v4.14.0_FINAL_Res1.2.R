#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Anti-inflammatory-MΦ heterogeneity analysis
# Clean-B
# v4.14.0
#
# PURPOSE
#   Resolve internal states within the fixed v4.8.4 Anti-inflammatory-MΦ class.
#
# IMPORTANT
#   - Parent 5-subtype classification is NOT changed.
#   - This is a SECONDARY subclustering analysis within Anti-inflammatory-MΦ.
#   - Clean-B is the primary dataset.
#   - Biological sample remains the unit for abundance interpretation.
#
# WORKFLOW
#   1) Load Clean-B
#   2) Extract Anti-inflammatory-MΦ only
#   3) Re-normalize / variable features / scale / PCA
#   4) RPCA-style reintegration across biological samples
#   5) UMAP + clustering resolution sweep
#   6) choose default resolution = 0.8
#   7) subcluster markers
#   8) sample-level subcluster abundance
#   9) functional program mapping
#  10) representative-gene DotPlot
#  11) sample-level heatmap
#  12) save annotated object
#
# NOTE
#   Resolution sweep extends through Res2.5.
#   Final internal-state resolution is fixed at Res1.2.
#   Parent MΦ classification remains fixed at Res2.0.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4140)

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
    "AntiInflammatory_heterogeneity_CleanB_v4.14.0"
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
TARGET_CLASS <- "Anti-inflammatory-Mphi"

NFEATURES <- 3000L
NPCS <- 30L
DIMS_USE <- 1:20

# RPCA integration weight.
# Anti-inflammatory-MΦ subset contains fewer anchors than the full MΦ object,
# so use a k.weight smaller than the available anchor-cell count.
K_WEIGHT <- 50L

RESOLUTIONS <- c(
    0.4,
    0.6,
    0.8,
    1.0,
    1.2,
    1.5,
    2.0,
    2.5
)

# Final internal-state resolution.
# Parent MΦ classification remains fixed at Res2.0.
# Anti-inflammatory-MΦ internal heterogeneity is fixed here at Res1.2.
DEFAULT_RESOLUTION <- 1.2

MIN_PCT_MARKER <- 0.15
LOGFC_THRESHOLD_MARKER <- 0.25

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
# 1. Functional gene programs
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

REPRESENTATIVE_GENES <- unique(
    c(
        "Mrc1","Cd163","Il1rn","Mertk","Igf1","Hmox1",
        "Il10ra","Il10rb","Stat3","Socs3",
        "Mfge8","Gas6","Mmp12","Mmp13","Mmp14",
        "Il1b","Tnf","Ccl2","Cxcl10",
        "Spp1","Tgfb1","Pdgfb","Thbs1",
        "Trem2","Gpnmb","Cd9","Lpl","Apoe"
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

safe_filename <- function(x) {

    gsub(
        "[^A-Za-z0-9]+",
        "_",
        x
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

msg("Loading Clean-B: ", INPUT_RDS)

mphi <- readRDS(INPUT_RDS)

if (!inherits(mphi, "Seurat")) {
    stop("Input is not a Seurat object.")
}

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

msg("SAMPLE_COL = ", SAMPLE_COL)
msg("CLASS_COL = ", CLASS_COL)

# ------------------------------------------------------------------------------
# 5. Extract Anti-inflammatory-MΦ
# ------------------------------------------------------------------------------

target_cells <- rownames(mphi@meta.data)[
    as.character(
        mphi@meta.data[[CLASS_COL]]
    ) ==
        TARGET_CLASS
]

if (length(target_cells) < 100L) {
    stop(
        "Too few Anti-inflammatory-MΦ cells: ",
        length(target_cells)
    )
}

anti <- subset(
    mphi,
    cells = target_cells
)

DefaultAssay(anti) <- ASSAY_USE

anti$sample_v4140 <- as.character(
    anti@meta.data[[SAMPLE_COL]]
)

anti$condition_v4140 <- canonical_condition(
    anti$sample_v4140
)

msg(
    "Anti-inflammatory-MΦ cells: ",
    ncol(anti)
)

cell_count_sample <- anti@meta.data %>%
    rownames_to_column("cell") %>%
    count(
        sample_v4140,
        condition_v4140,
        name = "n_cells"
    )

write.csv(
    cell_count_sample,
    file.path(
        TAB_DIR,
        "00_Anti_inflammatory_Mphi_cell_counts_by_sample_v4.14.0.csv"
    ),
    row.names = FALSE
)

print(cell_count_sample)

# ------------------------------------------------------------------------------
# 6. Rebuild expression analysis from RNA
# ------------------------------------------------------------------------------

# Remove previous integration/clustering reductions where possible by
# creating a fresh Seurat object from raw counts + metadata.
counts <- get_layer_safe(
    anti,
    ASSAY_USE,
    "counts"
)

if (is.null(counts)) {
    stop("RNA counts layer not found.")
}

meta_anti <- anti@meta.data

anti_fresh <- CreateSeuratObject(
    counts = counts,
    meta.data = meta_anti,
    assay = "RNA",
    project = "AntiInflammatory_Mphi_CleanB_v4140"
)

anti_fresh$sample_v4140 <- as.character(
    anti_fresh@meta.data[[SAMPLE_COL]]
)

anti_fresh$condition_v4140 <- canonical_condition(
    anti_fresh$sample_v4140
)

# ------------------------------------------------------------------------------
# 7. Split by biological sample and preprocess
# ------------------------------------------------------------------------------

anti_list <- SplitObject(
    anti_fresh,
    split.by = "sample_v4140"
)

anti_list <- lapply(
    anti_list,
    function(obj) {

        obj <- NormalizeData(
            obj,
            normalization.method = "LogNormalize",
            scale.factor = 10000,
            verbose = FALSE
        )

        obj <- FindVariableFeatures(
            obj,
            selection.method = "vst",
            nfeatures = NFEATURES,
            verbose = FALSE
        )

        obj
    }
)

# ------------------------------------------------------------------------------
# 8. RPCA integration
# ------------------------------------------------------------------------------

msg("Finding RPCA integration features...")

integration_features <- SelectIntegrationFeatures(
    object.list = anti_list,
    nfeatures = NFEATURES
)

anti_list <- lapply(
    anti_list,
    function(obj) {

        obj <- ScaleData(
            obj,
            features = integration_features,
            verbose = FALSE
        )

        obj <- RunPCA(
            obj,
            features = integration_features,
            npcs = NPCS,
            verbose = FALSE
        )

        obj
    }
)

msg("Finding RPCA anchors...")

anchors <- FindIntegrationAnchors(
    object.list = anti_list,
    anchor.features = integration_features,
    reduction = "rpca",
    dims = DIMS_USE,
    verbose = FALSE
)

msg("Integrating Anti-inflammatory-MΦ...")

anti_int <- IntegrateData(
    anchorset = anchors,
    dims = DIMS_USE,
    new.assay.name = "integratedRPCA",
    k.weight = K_WEIGHT,
    verbose = FALSE
)

DefaultAssay(anti_int) <- "integratedRPCA"

anti_int <- ScaleData(
    anti_int,
    verbose = FALSE
)

anti_int <- RunPCA(
    anti_int,
    npcs = NPCS,
    verbose = FALSE
)

anti_int <- FindNeighbors(
    anti_int,
    reduction = "pca",
    dims = DIMS_USE,
    verbose = FALSE
)

anti_int <- RunUMAP(
    anti_int,
    reduction = "pca",
    dims = DIMS_USE,
    reduction.name = "anti.umap.rpca",
    reduction.key = "AntiRPCAUMAP_",
    umap.method = "uwot",
    metric = "cosine",
    seed.use = 4140,
    verbose = FALSE
)

# ------------------------------------------------------------------------------
# 9. Resolution sweep
# ------------------------------------------------------------------------------

for (res in RESOLUTIONS) {

    anti_int <- FindClusters(
        anti_int,
        resolution = res,
        algorithm = 1,
        verbose = FALSE
    )

    source_col <- paste0(
        "integratedRPCA_snn_res.",
        res
    )

    target_col <- paste0(
        "anti_res_",
        gsub(
            "\\.",
            "_",
            as.character(res)
        ),
        "_v4140"
    )

    if (source_col %in% colnames(anti_int@meta.data)) {

        anti_int@meta.data[[target_col]] <-
            as.character(
                anti_int@meta.data[[source_col]]
            )
    }
}

# ------------------------------------------------------------------------------
# 10. Resolution comparison UMAPs
# ------------------------------------------------------------------------------

resolution_plots <- list()
resolution_qc <- list()

for (res in RESOLUTIONS) {

    col_now <- paste0(
        "anti_res_",
        gsub(
            "\\.",
            "_",
            as.character(res)
        ),
        "_v4140"
    )

    if (!col_now %in% colnames(anti_int@meta.data)) {
        next
    }

    ids_now <- factor(
        anti_int@meta.data[[col_now]]
    )

    Idents(anti_int) <- ids_now

    n_clusters <- nlevels(ids_now)

    cluster_sizes <- table(ids_now)

    resolution_qc[[as.character(res)]] <- tibble(
        resolution = res,
        n_clusters = n_clusters,
        smallest_cluster = min(cluster_sizes),
        median_cluster = median(cluster_sizes),
        largest_cluster = max(cluster_sizes)
    )

    p <- DimPlot(
        anti_int,
        reduction = "anti.umap.rpca",
        group.by = col_now,
        label = TRUE,
        repel = TRUE,
        raster = FALSE,
        pt.size = 0.55
    ) +
        NoLegend() +
        labs(
            title = paste0(
                "Anti-inflammatory-MΦ RPCA | resolution ",
                res
            )
        ) +
        theme_classic(
            base_size = 10
        ) +
        theme(
            plot.title = element_text(
                face = "bold"
            )
        )

    resolution_plots[[as.character(res)]] <- p
}

resolution_qc_df <- bind_rows(
    resolution_qc
)

write.csv(
    resolution_qc_df,
    file.path(
        TAB_DIR,
        "01_resolution_sweep_QC_v4.14.0.csv"
    ),
    row.names = FALSE
)

p_resolution_grid <- wrap_plots(
    resolution_plots,
    ncol = 3
) +
    plot_annotation(
        title =
            "Anti-inflammatory-MΦ internal subclustering resolution sweep",
        subtitle =
            "Clean-B | sample-level RPCA integration"
    )

save_pdf(
    "01_Anti_inflammatory_Mphi_resolution_sweep_v4.14.0.pdf",
    p_resolution_grid,
    15,
    10
)

# ------------------------------------------------------------------------------
# 11. Choose default resolution
# ------------------------------------------------------------------------------

DEFAULT_COL <- paste0(
    "anti_res_",
    gsub(
        "\\.",
        "_",
        as.character(
            DEFAULT_RESOLUTION
        )
    ),
    "_v4140"
)

if (!DEFAULT_COL %in% colnames(anti_int@meta.data)) {
    stop(
        "Default resolution column missing: ",
        DEFAULT_COL
    )
}

anti_int$anti_subcluster_v4140 <- factor(
    anti_int@meta.data[[DEFAULT_COL]]
)

Idents(anti_int) <- anti_int$anti_subcluster_v4140

# ------------------------------------------------------------------------------
# 12. Default UMAP
# ------------------------------------------------------------------------------

p_default <- DimPlot(
    anti_int,
    reduction = "anti.umap.rpca",
    group.by = "anti_subcluster_v4140",
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    pt.size = 0.65
) +
    labs(
        title =
            "Anti-inflammatory-MΦ internal states",
        subtitle =
            paste0(
                "Clean-B | RPCA | final internal-state resolution ",
                DEFAULT_RESOLUTION
            ),
        color =
            "Subcluster"
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
    "02_Anti_inflammatory_Mphi_final_Res1.2_subcluster_UMAP_v4.14.0.pdf",
    p_default,
    9,
    7
)

# ------------------------------------------------------------------------------
# 13. Condition UMAP
# ------------------------------------------------------------------------------

p_condition <- DimPlot(
    anti_int,
    reduction = "anti.umap.rpca",
    group.by = "condition_v4140",
    split.by = "condition_v4140",
    raster = FALSE,
    pt.size = 0.55,
    ncol = 2
) +
    labs(
        title =
            "Anti-inflammatory-MΦ distribution by condition"
    ) +
    theme_classic(
        base_size = 10
    )

save_pdf(
    "03_Anti_inflammatory_Mphi_condition_UMAP_v4.14.0.pdf",
    p_condition,
    12,
    9
)

# ------------------------------------------------------------------------------
# 13b. Final Res1.2 subcluster audit
# ------------------------------------------------------------------------------

final_cluster_counts <- anti_int@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = sample_v4140,
        condition = as.character(
            condition_v4140
        ),
        subcluster = as.character(
            anti_subcluster_v4140
        )
    ) %>%
    count(
        sample,
        condition,
        subcluster,
        name = "n_cells"
    ) %>%
    arrange(
        as.numeric(
            subcluster
        ),
        factor(
            sample,
            levels = SAMPLE_ORDER
        )
    )

write.csv(
    final_cluster_counts,
    file.path(
        TAB_DIR,
        "01b_final_Res1.2_subcluster_cell_counts_by_sample_v4.14.0.csv"
    ),
    row.names = FALSE
)

final_cluster_totals <- final_cluster_counts %>%
    group_by(
        subcluster
    ) %>%
    summarise(
        total_cells = sum(
            n_cells
        ),
        n_samples_present = n_distinct(
            sample[
                n_cells > 0
            ]
        ),
        .groups = "drop"
    ) %>%
    arrange(
        as.numeric(
            subcluster
        )
    )

write.csv(
    final_cluster_totals,
    file.path(
        TAB_DIR,
        "01c_final_Res1.2_subcluster_total_cells_v4.14.0.csv"
    ),
    row.names = FALSE
)

print(
    final_cluster_totals
)

# ------------------------------------------------------------------------------
# 14. Subcluster markers
#
# Seurat v5:
# RNA layers can remain split after sample-level integration.
# JoinLayers() is required before FindAllMarkers().
# ------------------------------------------------------------------------------

DefaultAssay(
    anti_int
) <- "RNA"

msg(
    "Joining RNA layers before marker analysis..."
)

anti_int <- JoinLayers(
    anti_int,
    assay = "RNA"
)

msg(
    "RNA layers after JoinLayers: ",
    paste(
        Layers(
            anti_int[["RNA"]]
        ),
        collapse = ", "
    )
)

if (!"data" %in%
    Layers(
        anti_int[["RNA"]]
    )) {

    msg(
        "Normalized RNA data layer not found; running NormalizeData()."
    )

    anti_int <- NormalizeData(
        anti_int,
        assay = "RNA",
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
    )
}

Idents(
    anti_int
) <- anti_int$anti_subcluster_v4140

print(
    table(
        Idents(
            anti_int
        )
    )
)

if (
    length(
        unique(
            Idents(
                anti_int
            )
        )
    ) < 2L
) {

    stop(
        "Fewer than 2 Anti-inflammatory-MΦ subclusters are present."
    )
}

msg(
    "Finding Anti-inflammatory-MΦ subcluster markers..."
)

markers <- FindAllMarkers(
    object = anti_int,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = MIN_PCT_MARKER,
    logfc.threshold = LOGFC_THRESHOLD_MARKER,
    test.use = "wilcox",
    verbose = TRUE
)

if (
    is.null(
        markers
    ) ||
    nrow(
        markers
    ) == 0L
) {

    stop(
        paste0(
            "FindAllMarkers returned no markers after JoinLayers().\n",
            "Check RNA layers, identity assignments, and marker thresholds."
        )
    )
}

write.csv(
    markers,
    file.path(
        TAB_DIR,
        "02_Anti_inflammatory_subcluster_markers_all_v4.14.0.csv"
    ),
    row.names = FALSE
)

top_markers <- markers %>%
    group_by(
        cluster
    ) %>%
    arrange(
        desc(
            avg_log2FC
        ),
        .by_group = TRUE
    ) %>%
    slice_head(
        n = 15
    ) %>%
    ungroup()

write.csv(
    top_markers,
    file.path(
        TAB_DIR,
        "03_Anti_inflammatory_subcluster_top15_markers_v4.14.0.csv"
    ),
    row.names = FALSE
)

msg(
    "Marker analysis completed: ",
    nrow(
        markers
    ),
    " marker rows."
)

# ------------------------------------------------------------------------------
# 15. Sample-level subcluster abundance
# ------------------------------------------------------------------------------

abundance <- anti_int@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = sample_v4140,
        condition = as.character(
            condition_v4140
        ),
        subcluster = as.character(
            anti_subcluster_v4140
        )
    ) %>%
    count(
        sample,
        condition,
        subcluster,
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
        fraction_within_Anti_inflammatory =
            n_cells /
            total_Anti_inflammatory_Mphi,

        percent_within_Anti_inflammatory =
            100 *
            fraction_within_Anti_inflammatory
    )

write.csv(
    abundance,
    file.path(
        TAB_DIR,
        "04_subcluster_abundance_by_sample_v4.14.0.csv"
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
        y = percent_within_Anti_inflammatory
    )
) +
    geom_point(
        aes(
            shape = sample
        ),
        size = 2.7
    ) +
    facet_wrap(
        ~ subcluster,
        scales = "free_y",
        ncol = 3
    ) +
    labs(
        title =
            "Anti-inflammatory-MΦ internal-state abundance",
        subtitle =
            "Each point = biological sample",
        x =
            NULL,
        y =
            "Subcluster / Anti-inflammatory-MΦ (%)",
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
            face = "bold"
        )
    )

save_pdf(
    "04_Anti_inflammatory_subcluster_abundance_v4.14.0.pdf",
    p_abundance,
    12,
    9
)

# ------------------------------------------------------------------------------
# 16. Functional program scores
# ------------------------------------------------------------------------------

rna_data <- get_layer_safe(
    anti_int,
    "RNA",
    "data"
)

if (is.null(rna_data)) {
    stop("Normalized RNA data layer not available after integration.")
}

PROGRAM_USE <- lapply(
    PROGRAMS,
    intersect,
    y = rownames(
        rna_data
    )
)

program_detection <- bind_rows(
    lapply(
        names(PROGRAM_USE),
        function(program_name) {

            tibble(
                program =
                    program_name,

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
    program_detection,
    file.path(
        TAB_DIR,
        "05_program_gene_detection_v4.14.0.csv"
    ),
    row.names = FALSE
)

score_df <- tibble(
    cell = colnames(
        anti_int
    )
)

for (program_name in names(PROGRAM_USE)) {

    score_df[[program_name]] <-
        mean_expression_score(
            rna_data,
            PROGRAM_USE[[program_name]]
        )
}

score_df <- score_df %>%
    left_join(
        anti_int@meta.data %>%
            rownames_to_column("cell") %>%
            transmute(
                cell = cell,
                sample = sample_v4140,
                condition = as.character(
                    condition_v4140
                ),
                subcluster = as.character(
                    anti_subcluster_v4140
                )
            ),
        by = "cell"
    )

write.csv(
    score_df,
    file.path(
        TAB_DIR,
        "06_cell_level_functional_scores_v4.14.0.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 17. Subcluster-level functional heatmap
# ------------------------------------------------------------------------------

subcluster_program <- score_df %>%
    pivot_longer(
        cols = all_of(
            names(PROGRAM_USE)
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
    subcluster_program,
    file.path(
        TAB_DIR,
        "07_subcluster_functional_program_means_v4.14.0.csv"
    ),
    row.names = FALSE
)

heat_df <- subcluster_program %>%
    pivot_wider(
        names_from = subcluster,
        values_from = mean_score
    )

heat_mat <- as.matrix(
    heat_df[
        ,
        -1,
        drop = FALSE
    ]
)

rownames(heat_mat) <- heat_df$program

heat_z <- t(
    scale(
        t(
            heat_mat
        )
    )
)

heat_z[
    !is.finite(
        heat_z
    )
] <- 0

HEAT_LIMIT <- 2.0

heat_plot <- pmax(
    pmin(
        heat_z,
        HEAT_LIMIT
    ),
    -HEAT_LIMIT
)

write.csv(
    data.frame(
        program = rownames(
            heat_z
        ),
        heat_z,
        check.names = FALSE
    ),
    file.path(
        TAB_DIR,
        "08_subcluster_functional_program_row_zscore_v4.14.0.csv"
    ),
    row.names = FALSE
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "05_Anti_inflammatory_subcluster_functional_heatmap_v4.14.0.pdf"
    ),
    width = 9,
    height = 6.5
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
        -HEAT_LIMIT,
        HEAT_LIMIT,
        length.out = 102
    ),
    border_color = "white",
    fontsize_row = 9,
    fontsize_col = 8,
    angle_col = 45,
    main = paste0(
        "Anti-inflammatory-MΦ internal functional states\n",
        "subcluster mean score | row z-score"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 18. Program score UMAPs
# ------------------------------------------------------------------------------

umap_emb <- Embeddings(
    anti_int,
    reduction = "anti.umap.rpca"
)

umap_df <- tibble(
    cell = rownames(
        umap_emb
    ),
    UMAP_1 = umap_emb[, 1],
    UMAP_2 = umap_emb[, 2]
) %>%
    left_join(
        score_df,
        by = "cell"
    )

program_umap_plots <- list()

for (program_name in names(PROGRAM_USE)) {

    df <- umap_df %>%
        transmute(
            cell = cell,
            UMAP_1 = UMAP_1,
            UMAP_2 = UMAP_2,
            score = .data[[program_name]]
        )

    limits_now <- quantile(
        df$score,
        probs = c(
            0.02,
            0.98
        ),
        na.rm = TRUE
    )

    midpoint_now <- median(
        df$score,
        na.rm = TRUE
    )

    p <- ggplot(
        df,
        aes(
            UMAP_1,
            UMAP_2,
            color = score
        )
    ) +
        geom_point(
            size = 0.55,
            alpha = 0.90
        ) +
        scale_color_gradient2(
            low = "#0033FF",
            mid = "#FFFFFF",
            high = "#FF1A1A",
            midpoint = midpoint_now,
            limits = limits_now,
            oob = scales::squish
        ) +
        coord_equal() +
        labs(
            title = program_name,
            x = NULL,
            y = NULL,
            color = "Score"
        ) +
        theme_classic(
            base_size = 9
        ) +
        theme(
            plot.title = element_text(
                face = "bold",
                hjust = 0.5
            )
        )

    program_umap_plots[[program_name]] <- p
}

p_program_umaps <- wrap_plots(
    program_umap_plots,
    ncol = 3
) +
    plot_annotation(
        title =
            "Anti-inflammatory-MΦ functional programs on RPCA UMAP"
    )

save_pdf(
    "06_Anti_inflammatory_functional_program_UMAPs_v4.14.0.pdf",
    p_program_umaps,
    14,
    12
)

# ------------------------------------------------------------------------------
# 19. Representative-gene custom DotPlot by subcluster
# ------------------------------------------------------------------------------

REP_USE <- intersect(
    REPRESENTATIVE_GENES,
    rownames(
        rna_data
    )
)

dot_rows <- list()

for (cl in levels(
    anti_int$anti_subcluster_v4140
)) {

    cells <- colnames(
        anti_int
    )[
        as.character(
            anti_int$anti_subcluster_v4140
        ) ==
            cl
    ]

    if (length(cells) == 0L) {
        next
    }

    mat <- rna_data[
        REP_USE,
        cells,
        drop = FALSE
    ]

    dot_rows[[cl]] <- tibble(
        subcluster = cl,
        gene = REP_USE,
        avg_expr = as.numeric(
            Matrix::rowMeans(
                mat
            )
        ),
        pct_expr = 100 *
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

            z[
                !is.finite(
                    z
                )
            ] <- 0

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
        "09_subcluster_representative_gene_DotPlot_numeric_v4.14.0.csv"
    ),
    row.names = FALSE
)

p_dot <- ggplot(
    dot_df,
    aes(
        x = gene,
        y = subcluster
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
            "Anti-inflammatory-MΦ internal-state gene validation",
        x =
            NULL,
        y =
            "Subcluster"
    ) +
    theme_classic(
        base_size = 9
    ) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        axis.text.x = element_text(
            angle = 60,
            hjust = 1
        )
    )

save_pdf(
    "07_Anti_inflammatory_subcluster_representative_gene_DotPlot_v4.14.0.pdf",
    p_dot,
    14,
    6.5
)

# ------------------------------------------------------------------------------
# 20. Sample × subcluster functional heatmap
# ------------------------------------------------------------------------------

sample_subcluster_program <- score_df %>%
    pivot_longer(
        cols = all_of(
            names(PROGRAM_USE)
        ),
        names_to = "program",
        values_to = "score"
    ) %>%
    group_by(
        sample,
        condition,
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
    sample_subcluster_program,
    file.path(
        TAB_DIR,
        "10_sample_subcluster_functional_program_means_v4.14.0.csv"
    ),
    row.names = FALSE
)

sample_heat <- sample_subcluster_program %>%
    mutate(
        column_id = paste(
            subcluster,
            sample,
            sep = " | "
        )
    ) %>%
    select(
        program,
        column_id,
        mean_score
    ) %>%
    pivot_wider(
        names_from = column_id,
        values_from = mean_score
    )

sample_heat_mat <- as.matrix(
    sample_heat[
        ,
        -1,
        drop = FALSE
    ]
)

rownames(
    sample_heat_mat
) <- sample_heat$program

sample_heat_z <- t(
    scale(
        t(
            sample_heat_mat
        )
    )
)

sample_heat_z[
    !is.finite(
        sample_heat_z
    )
] <- 0

sample_heat_plot <- pmax(
    pmin(
        sample_heat_z,
        2
    ),
    -2
)

grDevices::cairo_pdf(
    file.path(
        FIG_DIR,
        "08_Anti_inflammatory_sample_subcluster_functional_heatmap_v4.14.0.pdf"
    ),
    width = 13,
    height = 6.5
)

pheatmap::pheatmap(
    sample_heat_plot,
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
    fontsize_row = 9,
    fontsize_col = 6,
    angle_col = 45,
    main = paste0(
        "Anti-inflammatory-MΦ internal states by biological sample\n",
        "functional-program means | row z-score"
    )
)

grDevices::dev.off()

# ------------------------------------------------------------------------------
# 21. Sham -> Tx subcluster effect table
# ------------------------------------------------------------------------------

subcluster_treatment_effect <- abundance %>%
    filter(
        condition %in%
            c(
                "Sham",
                "Tx"
            )
    ) %>%
    group_by(
        subcluster,
        condition
    ) %>%
    summarise(
        mean_percent =
            mean(
                percent_within_Anti_inflammatory,
                na.rm = TRUE
            ),
        min_percent =
            min(
                percent_within_Anti_inflammatory,
                na.rm = TRUE
            ),
        max_percent =
            max(
                percent_within_Anti_inflammatory,
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
    subcluster_treatment_effect,
    file.path(
        TAB_DIR,
        "11_subcluster_abundance_effect_Tx_vs_Sham_v4.14.0.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 22. Save annotated object
# ------------------------------------------------------------------------------

OUTPUT_RDS <- file.path(
    RDS_DIR,
    "Mouse_Mphi_AntiInflammatory_heterogeneity_CleanB_v4.14.0.rds"
)

saveRDS(
    anti_int,
    OUTPUT_RDS
)

# ------------------------------------------------------------------------------
# 23. README
# ------------------------------------------------------------------------------

readme <- c(

    "Mouse MASH Anti-inflammatory-MΦ heterogeneity Clean-B v4.14.0",
    "",

    paste0(
        "Input: ",
        INPUT_RDS
    ),
    "",

    "Parent classification:",
    "  v4.8.4 macrophage class is fixed.",
    "  This analysis does NOT redefine the five major MΦ classes.",
    "",

    "Primary subset:",
    "  Anti-inflammatory-MΦ only",
    "",

    "Integration:",
    "  biological-sample split",
    "  LogNormalize",
    paste0(
        "  variable features = ",
        NFEATURES
    ),
    "  RPCA anchors",
    paste0(
        "  integration dims = ",
        min(DIMS_USE),
        "-",
        max(DIMS_USE)
    ),
    paste0(
        "  IntegrateData k.weight = ",
        K_WEIGHT
    ),
    "  UMAP method = uwot",
    "  UMAP metric = cosine",
    "",

    "Resolution sweep:",
    paste(
        RESOLUTIONS,
        collapse = ", "
    ),
    paste0(
        "Provisional default resolution = ",
        DEFAULT_RESOLUTION
    ),
    "  Compare through Res2.5 before fixing the final internal-state resolution.",
    "",

    "Outputs:",
    "  01 resolution sweep UMAP",
    "  02 final Res1.2 subcluster UMAP",
    "  03 condition UMAP",
    "  04 subcluster abundance",
    "  05 subcluster functional heatmap",
    "  06 functional-program UMAPs",
    "  07 representative-gene DotPlot",
    "  08 sample x subcluster functional heatmap",
    "",

    "Interpretation principles:",
    "  Use biological samples for abundance interpretation.",
    "  STD vs CDAHFD is descriptive n=1 vs n=1.",
    "  Sham vs Tx is n=2 vs n=2.",
    "  Subcluster labels should be assigned only after marker/program review.",
    "",

    "Next:",
    "  use validated internal states to inform output-ligand analysis v4.14.1"
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_v4.14.0.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.14.0.txt"
    )
)

# ------------------------------------------------------------------------------
# 24. Final
# ------------------------------------------------------------------------------

msg("DONE.")
msg("Output: ", OUTPUT_DIR)
msg("Annotated RDS: ", OUTPUT_RDS)

print(
    resolution_qc_df
)

print(
    subcluster_treatment_effect
)
