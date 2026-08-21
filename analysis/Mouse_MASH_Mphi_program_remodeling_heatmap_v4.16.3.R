#!/usr/bin/env Rscript

# ==============================================================================
# Mouse MASH MΦ
# Program-level remodeling heatmap
# v4.16.3
#
# PURPOSE
#   Manuscript-oriented summary of Sham -> Tx functional-program remodeling.
#
# DESIGN
#   Rows    = FINAL MΦ subtypes
#   Columns = functional programs
#   Fill    = mean Tx-vs-Sham effect across Tx17 and Tx5
#   Symbol  = replicate direction consistency
#
# INPUT
#   Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds
#
# FINAL ANNOTATION
#   macrophage_class_Res2_FINAL_v4145
#
# IMPORTANT
#   - Parent clustering remains Res2.0.
#   - No reannotation.
#   - Sham vs Tx uses biological sample-level pseudobulk.
#   - n=2 Sham vs n=2 Tx: emphasize effect size and replicate concordance.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(4163)

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

INPUT_RDS <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "RDS",
    "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

OUTPUT_DIR <- file.path(
    ROOT,
    "Mouse_MASH_Mphi_RDS",
    "Mphi_Res2_CleanB_FINAL_v4.14.5",
    "FINAL_Program_Remodeling_v4.16.3"
)

FIG_DIR <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR <- file.path(OUTPUT_DIR, "Tables")
LOG_DIR <- file.path(OUTPUT_DIR, "Logs")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

ASSAY_USE <- "RNA"
FINAL_CLASS_COL <- "macrophage_class_Res2_FINAL_v4145"

CLASS_ORDER <- c(
    "Inflammatory-Mphi",
    "Anti-inflammatory-Mphi",
    "ECM-associated inflammatory-Mphi",
    "Repair/Resolution-Mphi",
    "Lipid-associated/TREM2-Mphi"
)

CLASS_LABELS <- c(
    "Inflammatory-Mphi" = "Inflammatory-MΦ",
    "Anti-inflammatory-Mphi" = "Anti-inflammatory-MΦ",
    "ECM-associated inflammatory-Mphi" = "ECM-associated inflammatory-MΦ",
    "Repair/Resolution-Mphi" = "Repair/Resolution-MΦ",
    "Lipid-associated/TREM2-Mphi" = "Lipid-associated/TREM2-MΦ"
)

SAMPLE_ORDER <- c("Sham1", "Sham20", "Tx17", "Tx5")

PROGRAMS <- list(
    Inflammatory = c(
        "Il1b","Tnf","Cxcl10","Ccl2","Stat1"
    ),
    Anti_inflammatory = c(
        "Cd163","Mrc1","Il1rn","Mertk","Igf1","Hmox1"
    ),
    IL10_STAT3 = c(
        "Il10ra","Il10rb","Jak1","Tyk2","Stat3","Socs3","Bcl3","Il1rn"
    ),
    ECM_associated = c(
        "Thbs1","Fn1","Tgfb1","Col1a1","Col1a2","Col3a1"
    ),
    Repair_Resolution = c(
        "Mfge8","Gas6","Mmp13","Mmp14","Plau"
    ),
    Lipid_TREM2 = c(
        "Trem2","Gpnmb","Cd9","Lpl","Apoe"
    )
)

PROGRAM_ORDER <- c(
    "Inflammatory",
    "Anti_inflammatory",
    "IL10_STAT3",
    "ECM_associated",
    "Repair_Resolution",
    "Lipid_TREM2"
)

PROGRAM_LABELS <- c(
    "Inflammatory" = "Inflammatory",
    "Anti_inflammatory" = "Anti-inflammatory",
    "IL10_STAT3" = "IL10 / STAT3",
    "ECM_associated" = "ECM-associated",
    "Repair_Resolution" = "Repair / Resolution",
    "Lipid_TREM2" = "Lipid / TREM2"
)

required_packages <- c(
    "Seurat","SeuratObject","Matrix","dplyr","tidyr","tibble",
    "ggplot2","scales"
)

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
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
    library(scales)
})

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
        LayerData(object, assay = assay, layer = layer),
        error = function(e) NULL
    )
    if (is.null(x)) {
        x <- tryCatch(
            GetAssayData(object, assay = assay, slot = layer),
            error = function(e) NULL
        )
    }
    x
}

canonical_condition <- function(sample_name) {
    x <- as.character(sample_name)
    out <- rep(NA_character_, length(x))
    out[grepl("^Sham", x, ignore.case = TRUE)] <- "Sham"
    out[grepl("^Tx", x, ignore.case = TRUE)] <- "Tx"
    factor(out, levels = c("Sham", "Tx"))
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

if (!file.exists(INPUT_RDS)) {
    stop("FINAL RDS not found:\n", INPUT_RDS)
}

msg("Loading FINAL RDS: ", INPUT_RDS)

mphi <- readRDS(INPUT_RDS)

if (!inherits(mphi, "Seurat")) {
    stop("Input object is not a Seurat object.")
}

if (!FINAL_CLASS_COL %in% colnames(mphi@meta.data)) {
    stop("FINAL annotation column missing: ", FINAL_CLASS_COL)
}

DefaultAssay(mphi) <- ASSAY_USE

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
    stop("RNA counts layer not available.")
}

meta_cols <- colnames(mphi@meta.data)

SAMPLE_COL <- first_existing(
    meta_cols,
    c("sample_4group","sample","sample_id","orig.ident")
)

if (is.na(SAMPLE_COL)) {
    stop("Sample metadata column not found.")
}

meta <- mphi@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
        cell = cell,
        sample = as.character(.data[[SAMPLE_COL]]),
        macrophage_class = as.character(.data[[FINAL_CLASS_COL]])
    ) %>%
    filter(
        sample %in% SAMPLE_ORDER,
        macrophage_class %in% CLASS_ORDER
    ) %>%
    mutate(
        condition = canonical_condition(sample)
    )

# ------------------------------------------------------------------------------
# Program-gene audit
# ------------------------------------------------------------------------------

program_gene_audit <- bind_rows(
    lapply(
        names(PROGRAMS),
        function(program_now) {
            tibble(
                program = program_now,
                gene = PROGRAMS[[program_now]],
                detected = PROGRAMS[[program_now]] %in% rownames(counts)
            )
        }
    )
)

write.csv(
    program_gene_audit,
    file.path(TAB_DIR, "01_program_gene_audit_v4.16.3.csv"),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# Sample-level subtype pseudobulk
# ------------------------------------------------------------------------------

msg("Building sample-level subtype pseudobulk...")

pb_rows <- list()

for (class_now in CLASS_ORDER) {
    for (sample_now in SAMPLE_ORDER) {

        cells_now <- meta$cell[
            meta$macrophage_class == class_now &
            meta$sample == sample_now
        ]

        if (length(cells_now) == 0L) next

        pb_counts <- Matrix::rowSums(
            counts[, cells_now, drop = FALSE]
        )

        lib_size <- sum(pb_counts)

        pb_rows[[paste(class_now, sample_now, sep = "__")]] <- tibble(
            macrophage_class = class_now,
            sample = sample_now,
            condition = as.character(canonical_condition(sample_now)),
            gene = rownames(counts),
            raw_count = as.numeric(pb_counts),
            CPM = as.numeric(pb_counts) / lib_size * 1e6,
            logCPM = log2(as.numeric(pb_counts) / lib_size * 1e6 + 1)
        )
    }
}

pb_long <- bind_rows(pb_rows)

write.csv(
    pb_long,
    file.path(TAB_DIR, "02_sample_subtype_pseudobulk_all_genes_v4.16.3.csv"),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# Program score per subtype x sample
# ------------------------------------------------------------------------------

msg("Calculating program scores by subtype and biological sample...")

program_score_rows <- list()

for (class_now in CLASS_ORDER) {
    for (sample_now in SAMPLE_ORDER) {

        dat_now <- pb_long %>%
            filter(
                macrophage_class == class_now,
                sample == sample_now
            )

        if (nrow(dat_now) == 0L) next

        for (program_now in PROGRAM_ORDER) {

            genes_now <- intersect(
                PROGRAMS[[program_now]],
                dat_now$gene
            )

            if (length(genes_now) == 0L) next

            score_now <- dat_now %>%
                filter(gene %in% genes_now)

            program_score_rows[[
                paste(class_now, sample_now, program_now, sep = "__")
            ]] <- tibble(
                macrophage_class = class_now,
                sample = sample_now,
                condition = as.character(canonical_condition(sample_now)),
                program = program_now,
                n_genes = length(genes_now),
                program_score = mean(score_now$logCPM, na.rm = TRUE)
            )
        }
    }
}

program_scores <- bind_rows(program_score_rows) %>%
    mutate(
        macrophage_class = factor(
            macrophage_class,
            levels = CLASS_ORDER
        ),
        sample = factor(
            sample,
            levels = SAMPLE_ORDER
        ),
        condition = factor(
            condition,
            levels = c("Sham", "Tx")
        ),
        program = factor(
            program,
            levels = PROGRAM_ORDER
        )
    )

write.csv(
    program_scores,
    file.path(TAB_DIR, "03_program_scores_by_sample_subtype_v4.16.3.csv"),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# Tx17 / Tx5 effects relative to mean Sham
# ------------------------------------------------------------------------------

effect_wide <- program_scores %>%
    select(
        macrophage_class,
        program,
        sample,
        program_score
    ) %>%
    pivot_wider(
        names_from = sample,
        values_from = program_score
    )

required_cols <- c("Sham1","Sham20","Tx17","Tx5")

missing_cols <- setdiff(
    required_cols,
    colnames(effect_wide)
)

if (length(missing_cols) > 0L) {
    stop(
        "Missing sample columns in program score table: ",
        paste(missing_cols, collapse = ", ")
    )
}

program_effects <- effect_wide %>%
    mutate(
        Sham_mean = (Sham1 + Sham20) / 2,
        Tx_mean = (Tx17 + Tx5) / 2,

        Tx17_minus_Sham =
            Tx17 - Sham_mean,

        Tx5_minus_Sham =
            Tx5 - Sham_mean,

        mean_Tx_minus_Sham =
            Tx_mean - Sham_mean,

        replicate_concordant =
            sign(Tx17_minus_Sham) ==
            sign(Tx5_minus_Sham),

        replicate_effect_range =
            abs(Tx17_minus_Sham - Tx5_minus_Sham)
    )

write.csv(
    program_effects,
    file.path(TAB_DIR, "04_program_effects_Tx_vs_Sham_v4.16.3.csv"),
    row.names = FALSE
)

# ------------------------------------------------------------------------------
# Heatmap data
# ------------------------------------------------------------------------------

heat_df <- program_effects %>%
    transmute(
        macrophage_class = factor(
            macrophage_class,
            levels = rev(CLASS_ORDER)
        ),
        program = factor(
            program,
            levels = PROGRAM_ORDER
        ),
        mean_effect = mean_Tx_minus_Sham,
        concordant = replicate_concordant,
        Tx17_effect = Tx17_minus_Sham,
        Tx5_effect = Tx5_minus_Sham,
        replicate_range = replicate_effect_range
    )

write.csv(
    heat_df,
    file.path(TAB_DIR, "05_heatmap_plot_data_v4.16.3.csv"),
    row.names = FALSE
)

DISPLAY_LIMIT <- max(
    0.5,
    as.numeric(
        quantile(
            abs(heat_df$mean_effect),
            probs = 0.95,
            na.rm = TRUE
        )
    )
)

msg(
    "Heatmap display limit +/- ",
    round(DISPLAY_LIMIT, 3)
)

# ------------------------------------------------------------------------------
# Figure 1: Main manuscript heatmap
# ------------------------------------------------------------------------------

p_heat <- ggplot(
    heat_df,
    aes(
        x = program,
        y = macrophage_class,
        fill = mean_effect
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.8
    ) +
    geom_point(
        data = heat_df %>%
            filter(concordant),
        aes(
            x = program,
            y = macrophage_class
        ),
        inherit.aes = FALSE,
        shape = 21,
        fill = "black",
        color = "black",
        size = 2.8,
        stroke = 0.4
    ) +
    scale_fill_gradient2(
        low = "#0033FF",
        mid = "#FFFFFF",
        high = "#FF1A1A",
        midpoint = 0,
        limits = c(
            -DISPLAY_LIMIT,
            DISPLAY_LIMIT
        ),
        oob = scales::squish,
        name = "Mean program effect\nTx − Sham\nlog2(CPM+1)"
    ) +
    scale_x_discrete(
        labels = PROGRAM_LABELS[PROGRAM_ORDER],
        position = "top"
    ) +
    scale_y_discrete(
        labels = function(x) CLASS_LABELS[x]
    ) +
    labs(
        title = "Functional-program remodeling within FINAL MΦ subtypes",
        subtitle = "Sham → Tx | color = mean effect | ● = same direction in Tx17 and Tx5",
        x = NULL,
        y = NULL
    ) +
    coord_fixed(ratio = 1) +
    theme_classic(base_size = 10) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 13
        ),
        plot.subtitle = element_text(
            size = 9
        ),
        axis.text.x = element_text(
            angle = 35,
            hjust = 0,
            vjust = 0,
            face = "bold",
            size = 9
        ),
        axis.text.y = element_text(
            face = "bold",
            size = 9
        ),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(
            size = 9
        )
    )

save_pdf(
    "01_FINAL_Mphi_program_remodeling_heatmap_v4.16.3.pdf",
    p_heat,
    10.5,
    7
)

# ------------------------------------------------------------------------------
# Figure 2: Labeled replicate-effects heatmap
# ------------------------------------------------------------------------------

heat_df_label <- heat_df %>%
    mutate(
        replicate_label = paste0(
            sprintf("%+.2f", Tx17_effect),
            " / ",
            sprintf("%+.2f", Tx5_effect)
        )
    )

p_heat_labeled <- ggplot(
    heat_df_label,
    aes(
        x = program,
        y = macrophage_class,
        fill = mean_effect
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.8
    ) +
    geom_text(
        aes(label = replicate_label),
        size = 2.7
    ) +
    scale_fill_gradient2(
        low = "#0033FF",
        mid = "#FFFFFF",
        high = "#FF1A1A",
        midpoint = 0,
        limits = c(
            -DISPLAY_LIMIT,
            DISPLAY_LIMIT
        ),
        oob = scales::squish,
        name = "Mean Tx − Sham"
    ) +
    scale_x_discrete(
        labels = PROGRAM_LABELS[PROGRAM_ORDER],
        position = "top"
    ) +
    scale_y_discrete(
        labels = function(x) CLASS_LABELS[x]
    ) +
    labs(
        title = "Functional-program remodeling within FINAL MΦ subtypes",
        subtitle = "Numbers in each tile = Tx17 effect / Tx5 effect",
        x = NULL,
        y = NULL
    ) +
    coord_fixed(ratio = 1) +
    theme_classic(base_size = 9) +
    theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(
            angle = 35,
            hjust = 0,
            vjust = 0,
            face = "bold"
        ),
        axis.text.y = element_text(
            face = "bold"
        ),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )

save_pdf(
    "02_FINAL_Mphi_program_remodeling_heatmap_replicate_effects_v4.16.3.pdf",
    p_heat_labeled,
    10.5,
    7
)

# ------------------------------------------------------------------------------
# Figure 3: Alternate program-effect dotplot
# ------------------------------------------------------------------------------

p_dot <- ggplot(
    program_effects,
    aes(
        x = mean_Tx_minus_Sham,
        y = factor(
            program,
            levels = rev(PROGRAM_ORDER)
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
            xend = mean_Tx_minus_Sham,
            yend = factor(
                program,
                levels = rev(PROGRAM_ORDER)
            )
        ),
        linewidth = 0.55
    ) +
    geom_point(
        aes(
            shape = replicate_concordant
        ),
        size = 2.8
    ) +
    facet_wrap(
        ~ macrophage_class,
        ncol = 3,
        scales = "free_x",
        labeller = as_labeller(
            CLASS_LABELS
        )
    ) +
    scale_y_discrete(
        labels = function(x) PROGRAM_LABELS[x]
    ) +
    scale_shape_manual(
        values = c(
            `TRUE` = 16,
            `FALSE` = 1
        ),
        labels = c(
            `TRUE` = "Tx17/Tx5 concordant",
            `FALSE` = "Discordant"
        ),
        name = "Replicate direction"
    ) +
    labs(
        title = "Program-level Tx vs Sham effects by FINAL MΦ subtype",
        subtitle = "Biological sample-level pseudobulk",
        x = "Mean program effect: Tx − Sham [log2(CPM+1)]",
        y = NULL
    ) +
    theme_classic(base_size = 9) +
    theme(
        plot.title = element_text(
            face = "bold"
        ),
        strip.text = element_text(
            face = "bold"
        )
    )

save_pdf(
    "03_FINAL_Mphi_program_remodeling_dotplot_v4.16.3.pdf",
    p_dot,
    11,
    8
)

# ------------------------------------------------------------------------------
# Replicate concordance summary
# ------------------------------------------------------------------------------

concordance_summary <- program_effects %>%
    group_by(
        macrophage_class
    ) %>%
    summarise(
        n_programs = n(),
        n_concordant = sum(
            replicate_concordant,
            na.rm = TRUE
        ),
        concordance_fraction = mean(
            replicate_concordant,
            na.rm = TRUE
        ),
        mean_abs_effect = mean(
            abs(mean_Tx_minus_Sham),
            na.rm = TRUE
        ),
        .groups = "drop"
    )

write.csv(
    concordance_summary,
    file.path(
        TAB_DIR,
        "06_program_replicate_concordance_summary_v4.16.3.csv"
    ),
    row.names = FALSE
)

program_ranking <- program_effects %>%
    arrange(
        macrophage_class,
        desc(abs(mean_Tx_minus_Sham))
    ) %>%
    mutate(
        direction = case_when(
            mean_Tx_minus_Sham > 0 ~ "Higher in Tx",
            mean_Tx_minus_Sham < 0 ~ "Lower in Tx",
            TRUE ~ "No change"
        )
    )

write.csv(
    program_ranking,
    file.path(
        TAB_DIR,
        "07_program_effect_ranking_v4.16.3.csv"
    ),
    row.names = FALSE
)

readme <- c(
    "Mouse MASH MΦ functional-program remodeling v4.16.3",
    "",
    paste0("Input FINAL RDS: ", INPUT_RDS),
    "",
    paste0("FINAL annotation column: ", FINAL_CLASS_COL),
    "",
    "Analysis design:",
    "  biological sample-level subtype pseudobulk",
    "  Sham1 / Sham20 / Tx17 / Tx5",
    "  functional program score = mean log2(CPM+1) across genes in each program",
    "",
    "Main figure:",
    "  rows = FINAL MΦ subtype",
    "  columns = functional program",
    "  fill = mean Tx - Sham program effect",
    "  black dot = Tx17 and Tx5 show the same direction",
    "",
    "Programs:",
    "  Inflammatory",
    "  Anti-inflammatory",
    "  IL10/STAT3",
    "  ECM-associated",
    "  Repair/Resolution",
    "  Lipid/TREM2",
    "",
    "Statistical interpretation:",
    "  n=2 Sham and n=2 Tx.",
    "  This analysis emphasizes biological-sample effect size and replicate",
    "  directional consistency rather than formal p-values."
)

writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README_program_remodeling_v4.16.3.txt"
    )
)

capture.output(
    sessionInfo(),
    file = file.path(
        LOG_DIR,
        "sessionInfo_v4.16.3.txt"
    )
)

msg("DONE.")
msg("Output: ", OUTPUT_DIR)

print(
    program_effects %>%
        arrange(
            macrophage_class,
            program
        )
)

print(
    concordance_summary
)
