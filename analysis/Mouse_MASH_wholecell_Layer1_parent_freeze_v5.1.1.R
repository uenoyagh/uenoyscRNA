#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH whole-cell scRNA-seq
# Conservative Layer-1 parent freeze + target subset export
# v5.1.1
#
# INPUT:
#   Mouse_MASH_wholecell_Res2_Layer1_annotated_v5.1.0.rds
#
# PURPOSE:
#   1) Keep Res 2.0 as the frozen whole-cell parent resolution.
#   2) Correct clear Layer-1 failures identified in the v5.1.0 audit.
#   3) Define BROAD, CONSERVATIVE parent populations for:
#        - Hepatocyte
#        - HSC/Mesenchymal
#        - Endothelial (LSEC + vascular + boundary clusters)
#   4) Export independent subset RDS files for lineage-specific reclustering.
#
# IMPORTANT:
#   - No RPCA reintegration is performed.
#   - No lineage-specific final annotation is imposed here.
#   - Borderline clusters are retained conservatively when biologically relevant.
#   - Macrophage Clean-B v4.16.5 remains an independent frozen analysis.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(5110)

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(patchwork)
})

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

DATA_ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
    DATA_ROOT,
    "Mouse_MASH_RDS",
    "WholeCell_Layer1_Annotation_v5.1.0",
    "RDS",
    "Mouse_MASH_wholecell_Res2_Layer1_annotated_v5.1.0.rds"
)

OUT_ROOT <- file.path(
    DATA_ROOT,
    "Mouse_MASH_RDS",
    "WholeCell_Layer1_ParentFreeze_v5.1.1"
)

DIR_RDS <- file.path(OUT_ROOT, "RDS")
DIR_TAB <- file.path(OUT_ROOT, "Tables")
DIR_FIG <- file.path(OUT_ROOT, "Figures")
DIR_LOG <- file.path(OUT_ROOT, "Logs")

for (d in c(OUT_ROOT, DIR_RDS, DIR_TAB, DIR_FIG, DIR_LOG)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------------------------
# 1. Fixed columns
# ------------------------------------------------------------------------------

RES2_COL <- "wholecell_res.2_v502"
LAYER1_V510 <- "wholecell_layer1_FINAL_v510"
LAYER1_V511 <- "wholecell_layer1_FINAL_v511"
UMAP_USE <- "umapRPCA"

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

save_pdf <- function(p, filename, width, height) {
    ggsave(
        filename = file.path(DIR_FIG, filename),
        plot = p,
        device = cairo_pdf,
        width = width,
        height = height,
        units = "in",
        limitsize = FALSE
    )
}

safe_count <- function(x) {
    sort(table(x), decreasing = TRUE)
}

# ------------------------------------------------------------------------------
# 3. Load
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_RDS)) {
    stop("Input v5.1.0 RDS not found:\n", INPUT_RDS)
}

msg("Loading v5.1.0 RDS: ", INPUT_RDS)
obj <- readRDS(INPUT_RDS)

stopifnot(inherits(obj, "Seurat"))
stopifnot(RES2_COL %in% colnames(obj@meta.data))
stopifnot(LAYER1_V510 %in% colnames(obj@meta.data))
stopifnot(UMAP_USE %in% Reductions(obj))

cluster <- as.character(obj@meta.data[[RES2_COL]])
old <- as.character(obj@meta.data[[LAYER1_V510]])

# ------------------------------------------------------------------------------
# 4. Conservative corrections
# ------------------------------------------------------------------------------
#
# These corrections are deliberately minimal.
#
# v5.1.0 audit interpretation:
#
# cluster 31:
#   Strong B-cell evidence but v5.1.0 consensus became Other.
#   -> Correct to B.
#
# cluster 32:
#   NK/T boundary. Do NOT force a pure subtype.
#   -> Keep existing v5.1.0 lineage unless it was Other; then classify as NK/T.
#
# cluster 58:
#   Monocyte/neutrophil boundary.
#   -> Preserve as Myeloid-boundary rather than forcing lineage.
#
# cluster 38:
#   HSC legacy vs macrophage marker discordance.
#   -> Keep out of STRICT HSC parent; retain in BROAD HSC candidate pool.
#
# cluster 45:
#   LSEC legacy vs macrophage marker discordance.
#   -> Keep in BROAD endothelial candidate pool but not strict endothelial.
#
# cluster 64:
#   HSC legacy vs NK evidence.
#   -> Keep out of strict HSC; retain only as HSC candidate for audit.
#
# cluster 68:
#   HSC legacy vs cholangiocyte/mesothelial evidence.
#   -> Keep out of strict HSC; retain only as HSC candidate for audit.
#
# clusters 36,57,69:
#   LSEC / vascular endothelial boundary.
#   -> Explicitly retain in broad endothelial parent.
# ------------------------------------------------------------------------------

new <- old

# Clear correction
new[cluster == "31"] <- "B"

# Boundary labels only when v5.1.0 was Other
new[cluster == "32" & old == "Other"] <- "NK/T boundary"
new[cluster == "58"] <- "Monocyte/Neutrophil boundary"

obj@meta.data[[LAYER1_V511]] <- new

# ------------------------------------------------------------------------------
# 5. Parent-population definitions
# ------------------------------------------------------------------------------

# Hepatocyte:
# v5.1.0 hepatocyte annotation was comparatively clean.
hep_strict <- (
    obj@meta.data[[LAYER1_V511]] == "Hepatocyte"
)

# HSC:
# strict = current HSC/Mesenchymal excluding discordant clusters
HSC_DISCORDANT_CLUSTERS <- c("38", "64", "68")

hsc_strict <- (
    obj@meta.data[[LAYER1_V511]] == "HSC/Mesenchymal" &
        !(cluster %in% HSC_DISCORDANT_CLUSTERS)
)

# broad candidate = strict + discordant HSC legacy clusters
hsc_broad <- (
    hsc_strict |
        cluster %in% HSC_DISCORDANT_CLUSTERS
)

# Endothelial:
# strict = LSEC or vascular endothelial, excluding cluster 45
ENDOTHELIAL_BOUNDARY_CLUSTERS <- c("36", "57", "69")
ENDOTHELIAL_DISCORDANT_CLUSTERS <- c("45")

endo_strict <- (
    obj@meta.data[[LAYER1_V511]] %in%
        c("LSEC", "Vascular endothelial") &
        !(cluster %in% ENDOTHELIAL_DISCORDANT_CLUSTERS)
)

# broad = strict + LSEC/vascular boundaries + discordant cluster 45
endo_broad <- (
    endo_strict |
        cluster %in%
            c(
                ENDOTHELIAL_BOUNDARY_CLUSTERS,
                ENDOTHELIAL_DISCORDANT_CLUSTERS
            )
)

obj$parent_Hepatocyte_v511 <- hep_strict
obj$parent_HSC_strict_v511 <- hsc_strict
obj$parent_HSC_broad_v511 <- hsc_broad
obj$parent_Endothelial_strict_v511 <- endo_strict
obj$parent_Endothelial_broad_v511 <- endo_broad

# ------------------------------------------------------------------------------
# 6. Parent audit table
# ------------------------------------------------------------------------------

parent_audit <- tibble(
    cell = colnames(obj),
    cluster_res2 = cluster,
    layer1_v510 = old,
    layer1_v511 = new,
    Hepatocyte_parent = hep_strict,
    HSC_strict = hsc_strict,
    HSC_broad = hsc_broad,
    Endothelial_strict = endo_strict,
    Endothelial_broad = endo_broad
)

write.csv(
    parent_audit,
    file.path(
        DIR_TAB,
        "01_cell_level_parent_membership_v5.1.1.csv"
    ),
    row.names = FALSE
)

cluster_parent_audit <- parent_audit %>%
    group_by(
        cluster_res2,
        layer1_v510,
        layer1_v511
    ) %>%
    summarise(
        n_cells = n(),
        Hepatocyte_n = sum(Hepatocyte_parent),
        HSC_strict_n = sum(HSC_strict),
        HSC_broad_n = sum(HSC_broad),
        Endothelial_strict_n = sum(Endothelial_strict),
        Endothelial_broad_n = sum(Endothelial_broad),
        .groups = "drop"
    ) %>%
    arrange(as.numeric(cluster_res2))

write.csv(
    cluster_parent_audit,
    file.path(
        DIR_TAB,
        "02_Res2_cluster_parent_membership_v5.1.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Correction table
# ------------------------------------------------------------------------------

correction_table <- tibble(
    cluster = c(
        "31", "32", "58", "38", "45", "64", "68",
        "36", "57", "69"
    ),
    issue = c(
        "Strong B-cell marker evidence; v5.1.0 Other",
        "NK/T boundary",
        "Monocyte/Neutrophil boundary",
        "HSC legacy vs macrophage evidence",
        "LSEC legacy vs macrophage evidence",
        "HSC legacy vs NK evidence",
        "HSC legacy vs cholangiocyte/mesothelial evidence",
        "LSEC/vascular boundary",
        "LSEC/vascular boundary",
        "LSEC/vascular boundary"
    ),
    v5_1_1_action = c(
        "Correct Layer-1 to B",
        "Do not force subtype; flag boundary if needed",
        "Set boundary label",
        "Exclude strict HSC; retain broad HSC candidate",
        "Exclude strict endothelial; retain broad endothelial candidate",
        "Exclude strict HSC; retain broad HSC candidate",
        "Exclude strict HSC; retain broad HSC candidate",
        "Retain broad endothelial",
        "Retain broad endothelial",
        "Retain broad endothelial"
    )
)

write.csv(
    correction_table,
    file.path(
        DIR_TAB,
        "03_manual_corrections_and_boundary_policy_v5.1.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Parent summary
# ------------------------------------------------------------------------------

parent_summary <- tibble(
    population = c(
        "Hepatocyte",
        "HSC strict",
        "HSC broad",
        "Endothelial strict",
        "Endothelial broad"
    ),
    n_cells = c(
        sum(hep_strict),
        sum(hsc_strict),
        sum(hsc_broad),
        sum(endo_strict),
        sum(endo_broad)
    )
) %>%
    mutate(
        percent_whole = 100 * n_cells / ncol(obj)
    )

write.csv(
    parent_summary,
    file.path(
        DIR_TAB,
        "04_parent_population_summary_v5.1.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Whole-cell UMAP after conservative correction
# ------------------------------------------------------------------------------

layer_order <- c(
    "Hepatocyte",
    "HSC/Mesenchymal",
    "LSEC",
    "Vascular endothelial",
    "Kupffer/Macrophage",
    "Monocyte",
    "Monocyte/Neutrophil boundary",
    "Dendritic",
    "Neutrophil",
    "NK",
    "NK/T boundary",
    "T",
    "B",
    "Plasma",
    "Cholangiocyte",
    "Mesothelial",
    "Cycling",
    "Platelet",
    "Other"
)

obj@meta.data[[LAYER1_V511]] <- factor(
    obj@meta.data[[LAYER1_V511]],
    levels = layer_order
)

palette <- c(
    "Hepatocyte" = "#00A6A6",
    "HSC/Mesenchymal" = "#FF4FA3",
    "LSEC" = "#00CFFF",
    "Vascular endothelial" = "#2D7FF9",
    "Kupffer/Macrophage" = "#E53935",
    "Monocyte" = "#FF8C00",
    "Monocyte/Neutrophil boundary" = "#D49A00",
    "Dendritic" = "#8E44AD",
    "Neutrophil" = "#B8860B",
    "NK" = "#6A5ACD",
    "NK/T boundary" = "#536DFE",
    "T" = "#4169E1",
    "B" = "#00A86B",
    "Plasma" = "#7B1FA2",
    "Cholangiocyte" = "#00C853",
    "Mesothelial" = "#A1887F",
    "Cycling" = "#F06292",
    "Platelet" = "#795548",
    "Other" = "#9E9E9E"
)

p_all <- DimPlot(
    obj,
    reduction = UMAP_USE,
    group.by = LAYER1_V511,
    cols = palette,
    raster = FALSE,
    pt.size = 0.20
) +
    ggtitle(
        "Mouse MASH whole-cell Layer-1 parent freeze — Res 2.0"
    ) +
    theme_classic(base_size = 11) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 14,
            hjust = 0.5
        ),
        legend.text = element_text(size = 9)
    ) +
    guides(
        colour = guide_legend(
            override.aes = list(size = 3.5)
        )
    )

save_pdf(
    p_all,
    "01_wholecell_Layer1_parent_freeze_Res2_v5.1.1.pdf",
    11,
    8.5
)

# ------------------------------------------------------------------------------
# 10. Strict vs broad target UMAP
# ------------------------------------------------------------------------------

make_highlight <- function(
    cells_use,
    title,
    highlight_color
) {
    DimPlot(
        obj,
        reduction = UMAP_USE,
        cells.highlight = cells_use,
        cols = "grey88",
        cols.highlight = highlight_color,
        pt.size = 0.15,
        sizes.highlight = 0.32,
        raster = FALSE
    ) +
        ggtitle(title) +
        theme_classic(base_size = 10) +
        theme(
            plot.title = element_text(
                face = "bold",
                hjust = 0.5
            ),
            legend.position = "none"
        )
}

p_hep <- make_highlight(
    colnames(obj)[hep_strict],
    "Hepatocyte parent",
    "#00A6A6"
)

p_hsc_s <- make_highlight(
    colnames(obj)[hsc_strict],
    "HSC strict",
    "#FF4FA3"
)

p_hsc_b <- make_highlight(
    colnames(obj)[hsc_broad],
    "HSC broad candidate",
    "#FF4FA3"
)

p_endo_s <- make_highlight(
    colnames(obj)[endo_strict],
    "Endothelial strict",
    "#00CFFF"
)

p_endo_b <- make_highlight(
    colnames(obj)[endo_broad],
    "Endothelial broad candidate",
    "#00CFFF"
)

p_targets <- (
    p_hep |
        p_hsc_s |
        p_hsc_b
) /
    (
        p_endo_s |
            p_endo_b |
            plot_spacer()
    ) +
    plot_annotation(
        title =
            "Whole-cell parent populations for lineage-specific reclustering"
    )

save_pdf(
    p_targets,
    "02_target_parent_strict_vs_broad_UMAP_v5.1.1.pdf",
    16,
    10
)

# ------------------------------------------------------------------------------
# 11. Export subset RDS files
# ------------------------------------------------------------------------------

msg("Creating Hepatocyte subset...")
hep <- subset(
    obj,
    cells = colnames(obj)[hep_strict]
)

msg("Creating HSC BROAD subset...")
hsc <- subset(
    obj,
    cells = colnames(obj)[hsc_broad]
)

msg("Creating Endothelial BROAD subset...")
endo <- subset(
    obj,
    cells = colnames(obj)[endo_broad]
)

# Preserve provenance explicitly
hep$wholecell_parent_source_v511 <- "Hepatocyte"
hsc$wholecell_parent_source_v511 <- ifelse(
    hsc$parent_HSC_strict_v511,
    "HSC_strict",
    "HSC_boundary_candidate"
)
endo$wholecell_parent_source_v511 <- ifelse(
    endo$parent_Endothelial_strict_v511,
    "Endothelial_strict",
    "Endothelial_boundary_candidate"
)

HEP_RDS <- file.path(
    DIR_RDS,
    "Mouse_MASH_Hepatocyte_parent_v5.1.1.rds"
)

HSC_RDS <- file.path(
    DIR_RDS,
    "Mouse_MASH_HSC_Mesenchymal_parent_BROAD_v5.1.1.rds"
)

ENDO_RDS <- file.path(
    DIR_RDS,
    "Mouse_MASH_Endothelial_parent_BROAD_v5.1.1.rds"
)

WHOLE_RDS <- file.path(
    DIR_RDS,
    "Mouse_MASH_wholecell_Res2_Layer1_parent_frozen_v5.1.1.rds"
)

msg("Saving whole-cell frozen parent RDS...")
saveRDS(obj, WHOLE_RDS, compress = FALSE)

msg("Saving Hepatocyte parent RDS...")
saveRDS(hep, HEP_RDS, compress = FALSE)

msg("Saving HSC parent RDS...")
saveRDS(hsc, HSC_RDS, compress = FALSE)

msg("Saving Endothelial parent RDS...")
saveRDS(endo, ENDO_RDS, compress = FALSE)

# ------------------------------------------------------------------------------
# 12. Subset provenance table
# ------------------------------------------------------------------------------

subset_manifest <- tibble(
    subset = c(
        "Whole-cell frozen",
        "Hepatocyte parent",
        "HSC/Mesenchymal BROAD parent",
        "Endothelial BROAD parent"
    ),
    n_cells = c(
        ncol(obj),
        ncol(hep),
        ncol(hsc),
        ncol(endo)
    ),
    source = c(
        "Whole-cell Res2",
        "Layer1 Hepatocyte",
        "HSC strict + clusters 38/64/68",
        "LSEC + vascular + clusters 36/45/57/69"
    ),
    file = c(
        WHOLE_RDS,
        HEP_RDS,
        HSC_RDS,
        ENDO_RDS
    )
)

write.csv(
    subset_manifest,
    file.path(
        DIR_TAB,
        "05_subset_RDS_manifest_v5.1.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 13. Recommended next-stage parameters
# ------------------------------------------------------------------------------

next_stage <- tibble(
    lineage = c(
        "Hepatocyte",
        "HSC/Mesenchymal",
        "Endothelial"
    ),
    parent_RDS = c(
        HEP_RDS,
        HSC_RDS,
        ENDO_RDS
    ),
    initial_PCs_to_test = c(
        "1:30",
        "1:30",
        "1:30"
    ),
    resolutions_to_test = c(
        "0.4,0.6,0.8,1.0,1.2,1.5",
        "0.4,0.6,0.8,1.0,1.2,1.5",
        "0.4,0.6,0.8,1.0,1.2,1.5"
    ),
    biological_goal = c(
        "zonation + injury/stress states",
        "qHSC + activated/ECM + contractile states",
        "LSEC states + vascular endothelial separation"
    )
)

write.csv(
    next_stage,
    file.path(
        DIR_TAB,
        "06_recommended_lineage_reclustering_plan_v5.1.1.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# 14. README
# ------------------------------------------------------------------------------

readme <- c(
    "Mouse MASH whole-cell Layer-1 parent freeze v5.1.1",
    "===================================================",
    "",
    "Frozen parent resolution:",
    "Res 2.0",
    "",
    "Whole-cell parent cluster column:",
    RES2_COL,
    "",
    "Key policy:",
    "Do not re-integrate the whole-cell object at this stage.",
    "Use Res 2.0 only to define broad parent lineages.",
    "Perform lineage-specific PCA/neighbors/clustering after subsetting.",
    "",
    "Hepatocyte:",
    "Use v5.1.0 Hepatocyte Layer-1 population directly.",
    "",
    "HSC/Mesenchymal:",
    "Strict HSC excludes clusters 38, 64 and 68.",
    "BROAD HSC parent retains these clusters for lineage-specific audit.",
    "",
    "Endothelial:",
    "Strict endothelial = LSEC + vascular endothelial excluding cluster 45.",
    "BROAD endothelial retains boundary clusters 36, 45, 57 and 69.",
    "",
    "Clear Layer-1 correction:",
    "Cluster 31 -> B.",
    "",
    "Boundary policy:",
    "Cluster 58 is explicitly retained as Monocyte/Neutrophil boundary.",
    "Borderline HSC/endothelial clusters are not discarded before",
    "lineage-specific reclustering.",
    "",
    "Macrophage:",
    "The frozen macrophage Clean-B analysis through v4.16.5 is independent",
    "and is not modified by this script.",
    "",
    "Next step:",
    "Run independent lineage-resolution audits for Hepatocyte,",
    "HSC/Mesenchymal, and Endothelial parent RDS files."
)

writeLines(
    readme,
    file.path(
        OUT_ROOT,
        "README_WholeCell_Layer1_ParentFreeze_v5.1.1.txt"
    )
)

# ------------------------------------------------------------------------------
# 15. Session info / final report
# ------------------------------------------------------------------------------

capture.output(
    sessionInfo(),
    file = file.path(
        DIR_LOG,
        "sessionInfo_v5.1.1.txt"
    )
)

msg("DONE.")
msg("Output root: ", OUT_ROOT)
msg("Whole-cell frozen RDS: ", WHOLE_RDS)
msg("Hepatocyte RDS: ", HEP_RDS)
msg("HSC BROAD RDS: ", HSC_RDS)
msg("Endothelial BROAD RDS: ", ENDO_RDS)

print(parent_summary)
print(subset_manifest)
