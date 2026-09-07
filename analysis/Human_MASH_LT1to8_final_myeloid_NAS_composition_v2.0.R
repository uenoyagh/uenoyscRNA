#!/usr/bin/env Rscript

# Human_MASH_LT1to8_final_myeloid_NAS_composition_v2.0.R
#
# PURPOSE
# -------
# Quantify sample-level changes in the FINAL Human MΦ/Monocyte states
# across NAS0 -> NAS2 -> NAS4 -> NAS6.
#
# INPUTS
# ------
# 1) v1.9 final annotated myeloid object
# 2) v1.3 canonical whole-liver LT1-LT8 object, used only to obtain
#    per-sample whole-liver denominators.
#
# THIS SCRIPT DOES:
#   - calculate final myeloid abundance among captured whole-liver cells
#   - calculate each final state's fraction within the final myeloid compartment
#   - summarize per NAS stage using biological samples (n=2/stage)
#   - perform exploratory sample-level Spearman trend tests
#   - generate composition and trajectory figures
#   - run a sensitivity analysis excluding review-flagged clusters 13/14/16
#
# THIS SCRIPT DOES NOT:
#   - modify or save any RDS
#   - use cell-level tests
#   - perform DE / pseudobulk / pathway / cell-cell interaction analysis
#   - treat scRNA captured-cell fractions as absolute tissue cell abundance
#
# Statistical interpretation:
#   There are only 2 biological samples per NAS stage. Group summaries and
#   Spearman tests are exploratory/descriptive and should not be overinterpreted.

options(stringsAsFactors = FALSE, width = 220)

VERSION <- "v2.0"

MYELOID_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9/",
  "objects/Human_MASH_LT1to8_FinalMyeloidAnnotation_v1.9.rds"
)

WHOLE_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_Canonical_v1.3/",
  "objects/Human_MASH_LT1to8_Canonical_v1.3.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Human_MASH_RDS/Human_MASH_LT1to8_FinalMyeloid_NAS_Composition_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

for (p in c(MYELOID_RDS, WHOLE_RDS)) {
  if (!file.exists(p)) stop("Required RDS not found:\n", p)
}

suppressPackageStartupMessages({
  library(SeuratObject)
  library(ggplot2)
})

message("============================================================")
message("Human MASH LT1-LT8 final myeloid NAS composition ", VERSION)
message("============================================================")
message("READ-ONLY.")

my <- readRDS(MYELOID_RDS)
whole <- readRDS(WHOLE_RDS)

required_my <- c(
  "myeloid_state_FINAL_v19",
  "res0_8_cluster_frozen_v18",
  "sample_v17",
  "NAS_v17"
)

miss_my <- setdiff(required_my, colnames(my[[]]))
if (length(miss_my)) {
  stop("Missing myeloid metadata: ", paste(miss_my, collapse = ", "))
}

required_whole <- c(
  "canonical_sample_v13",
  "NAS_score_v13"
)

miss_whole <- setdiff(required_whole, colnames(whole[[]]))
if (length(miss_whole)) {
  stop("Missing whole-liver metadata: ", paste(miss_whole, collapse = ", "))
}

SAMPLE_ORDER <- c("LT7", "LT8", "LT5", "LT6", "LT3", "LT4", "LT1", "LT2")
EXPECTED_NAS <- c(
  LT7 = 0, LT8 = 0,
  LT5 = 2, LT6 = 2,
  LT3 = 4, LT4 = 4,
  LT1 = 6, LT2 = 6
)

STATE_ORDER <- c(
  "Resident/Kupffer-like_Mphi",
  "Repair/Resolution_Mphi",
  "Lipid-associated/TREM2_Mphi",
  "Inflammatory/Transitional_Mphi",
  "Inflammatory_Monocyte"
)

# ------------------------------------------------------------------
# 01 Whole-liver denominators
# ------------------------------------------------------------------

whole_sample <- as.character(whole$canonical_sample_v13)
whole_nas <- as.numeric(whole$NAS_score_v13)

whole_tab <- as.data.frame(
  table(sample = whole_sample),
  stringsAsFactors = FALSE
)
colnames(whole_tab)[2] <- "whole_liver_cells"

whole_tab$NAS <- unname(EXPECTED_NAS[whole_tab$sample])

if (anyNA(whole_tab$NAS)) {
  stop("Unexpected sample in whole-liver object.")
}

write.csv(
  whole_tab,
  file.path(TABDIR, "01_whole_liver_sample_denominators_v2.0.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 02 Final myeloid counts and captured-cell fraction
# ------------------------------------------------------------------

my_sample <- as.character(my$sample_v17)
my_nas <- as.numeric(my$NAS_v17)
my_state <- as.character(my$myeloid_state_FINAL_v19)
my_cluster <- as.character(my$res0_8_cluster_frozen_v18)

if (!setequal(unique(my_sample), names(EXPECTED_NAS))) {
  stop("Final myeloid object does not contain exactly LT1-LT8.")
}
if (!setequal(unique(my_state), STATE_ORDER)) {
  stop(
    "Unexpected final myeloid states:\n",
    paste(sort(unique(my_state)), collapse = " | ")
  )
}

my_total <- as.data.frame(
  table(sample = my_sample),
  stringsAsFactors = FALSE
)
colnames(my_total)[2] <- "final_myeloid_cells"

my_abundance <- merge(
  whole_tab,
  my_total,
  by = "sample",
  all.x = TRUE
)

my_abundance$final_myeloid_cells[is.na(my_abundance$final_myeloid_cells)] <- 0L
my_abundance$final_myeloid_fraction_of_whole <-
  my_abundance$final_myeloid_cells / my_abundance$whole_liver_cells

my_abundance <- my_abundance[
  match(SAMPLE_ORDER, my_abundance$sample),
  , drop = FALSE
]

write.csv(
  my_abundance,
  file.path(TABDIR, "02_final_myeloid_fraction_of_whole_by_sample_v2.0.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 03 State counts and fractions within final myeloid
# ------------------------------------------------------------------

state_counts <- as.data.frame(
  table(
    sample = my_sample,
    final_state = my_state
  ),
  stringsAsFactors = FALSE
)
colnames(state_counts)[3] <- "n_cells"

# Complete zero cells if any state absent from a sample.
grid <- expand.grid(
  sample = SAMPLE_ORDER,
  final_state = STATE_ORDER,
  stringsAsFactors = FALSE
)

state_counts <- merge(
  grid,
  state_counts,
  by = c("sample", "final_state"),
  all.x = TRUE
)
state_counts$n_cells[is.na(state_counts$n_cells)] <- 0L

state_counts <- merge(
  state_counts,
  my_total,
  by = "sample",
  all.x = TRUE
)

state_counts$NAS <- unname(EXPECTED_NAS[state_counts$sample])
state_counts$fraction_within_myeloid <-
  state_counts$n_cells / state_counts$final_myeloid_cells

state_counts$sample <- factor(
  state_counts$sample,
  levels = SAMPLE_ORDER
)
state_counts$final_state <- factor(
  state_counts$final_state,
  levels = STATE_ORDER
)

state_counts <- state_counts[
  order(state_counts$sample, state_counts$final_state),
  , drop = FALSE
]

write.csv(
  state_counts,
  file.path(TABDIR, "03_final_state_fraction_by_sample_v2.0.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 04 Biological-sample summaries by NAS
# ------------------------------------------------------------------

group_summary_list <- list()

for (st in STATE_ORDER) {
  for (nas_value in c(0, 2, 4, 6)) {
    z <- state_counts[
      as.character(state_counts$final_state) == st &
        state_counts$NAS == nas_value,
      "fraction_within_myeloid"
    ]

    group_summary_list[[length(group_summary_list) + 1L]] <- data.frame(
      final_state = st,
      NAS = nas_value,
      n_samples = length(z),
      mean_fraction = mean(z),
      sd_fraction = stats::sd(z),
      min_fraction = min(z),
      max_fraction = max(z),
      stringsAsFactors = FALSE
    )
  }
}

group_summary <- do.call(rbind, group_summary_list)

write.csv(
  group_summary,
  file.path(TABDIR, "04_final_state_fraction_by_NAS_summary_v2.0.csv"),
  row.names = FALSE
)

# Overall final-myeloid abundance by NAS.
abundance_summary <- do.call(
  rbind,
  lapply(c(0, 2, 4, 6), function(nas_value) {
    z <- my_abundance$final_myeloid_fraction_of_whole[
      my_abundance$NAS == nas_value
    ]

    data.frame(
      NAS = nas_value,
      n_samples = length(z),
      mean_fraction = mean(z),
      sd_fraction = stats::sd(z),
      min_fraction = min(z),
      max_fraction = max(z),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  abundance_summary,
  file.path(TABDIR, "05_total_myeloid_fraction_by_NAS_summary_v2.0.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 05 Exploratory sample-level Spearman trends
# ------------------------------------------------------------------

trend_list <- list()

for (st in STATE_ORDER) {
  z <- state_counts[
    as.character(state_counts$final_state) == st,
    , drop = FALSE
  ]

  z <- z[order(z$NAS, z$sample), , drop = FALSE]

  ct <- suppressWarnings(
    cor.test(
      z$NAS,
      z$fraction_within_myeloid,
      method = "spearman",
      exact = FALSE
    )
  )

  trend_list[[st]] <- data.frame(
    outcome = st,
    denominator = "final_myeloid",
    n_samples = nrow(z),
    spearman_rho = unname(ct$estimate),
    p_value_exploratory = ct$p.value,
    stringsAsFactors = FALSE
  )
}

# Total myeloid fraction among whole-liver captured cells.
ct_total <- suppressWarnings(
  cor.test(
    my_abundance$NAS,
    my_abundance$final_myeloid_fraction_of_whole,
    method = "spearman",
    exact = FALSE
  )
)

trend_list[["TOTAL_MYELOID"]] <- data.frame(
  outcome = "TOTAL_final_myeloid_fraction_of_whole",
  denominator = "whole_liver_captured_cells",
  n_samples = nrow(my_abundance),
  spearman_rho = unname(ct_total$estimate),
  p_value_exploratory = ct_total$p.value,
  stringsAsFactors = FALSE
)

trend_table <- do.call(rbind, trend_list)
rownames(trend_table) <- NULL

write.csv(
  trend_table,
  file.path(TABDIR, "06_exploratory_sample_level_Spearman_v2.0.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 07 Sensitivity analysis excluding review-flagged res0.8 clusters
#     13,14,16
# ------------------------------------------------------------------

REVIEW_CLUSTERS <- c("13", "14", "16")

sens_keep <- !(my_cluster %in% REVIEW_CLUSTERS)

sens_sample <- my_sample[sens_keep]
sens_state <- my_state[sens_keep]

sens_total <- as.data.frame(
  table(sample = sens_sample),
  stringsAsFactors = FALSE
)
colnames(sens_total)[2] <- "sensitivity_myeloid_cells"

sens_counts <- as.data.frame(
  table(
    sample = sens_sample,
    final_state = sens_state
  ),
  stringsAsFactors = FALSE
)
colnames(sens_counts)[3] <- "n_cells"

sens_grid <- expand.grid(
  sample = SAMPLE_ORDER,
  final_state = STATE_ORDER,
  stringsAsFactors = FALSE
)

sens_counts <- merge(
  sens_grid,
  sens_counts,
  by = c("sample", "final_state"),
  all.x = TRUE
)
sens_counts$n_cells[is.na(sens_counts$n_cells)] <- 0L

sens_counts <- merge(
  sens_counts,
  sens_total,
  by = "sample",
  all.x = TRUE
)

sens_counts$fraction_within_myeloid <-
  sens_counts$n_cells / sens_counts$sensitivity_myeloid_cells
sens_counts$NAS <- unname(EXPECTED_NAS[sens_counts$sample])

write.csv(
  sens_counts,
  file.path(TABDIR, "07_sensitivity_excluding_clusters13_14_16_v2.0.csv"),
  row.names = FALSE
)

# Compare primary vs sensitivity fractions.
primary_compare <- state_counts[
  , c("sample", "final_state", "fraction_within_myeloid"),
  drop = FALSE
]
colnames(primary_compare)[3] <- "primary_fraction"

sens_compare <- sens_counts[
  , c("sample", "final_state", "fraction_within_myeloid"),
  drop = FALSE
]
colnames(sens_compare)[3] <- "sensitivity_fraction"

sensitivity_delta <- merge(
  primary_compare,
  sens_compare,
  by = c("sample", "final_state"),
  all = TRUE
)
sensitivity_delta$absolute_delta <-
  abs(
    sensitivity_delta$primary_fraction -
      sensitivity_delta$sensitivity_fraction
  )

write.csv(
  sensitivity_delta,
  file.path(TABDIR, "08_primary_vs_sensitivity_fraction_delta_v2.0.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------

palette_values <- c(
  "Resident/Kupffer-like_Mphi"      = "#00A6D6",
  "Repair/Resolution_Mphi"          = "#00C853",
  "Lipid-associated/TREM2_Mphi"     = "#F0007A",
  "Inflammatory/Transitional_Mphi"  = "#FF7A00",
  "Inflammatory_Monocyte"           = "#7B1FA2"
)

# 1) 100% stacked composition by biological sample
plot_comp <- state_counts
plot_comp$sample <- factor(plot_comp$sample, levels = SAMPLE_ORDER)
plot_comp$final_state <- factor(plot_comp$final_state, levels = STATE_ORDER)

p_comp <- ggplot(
  plot_comp,
  aes(
    x = sample,
    y = fraction_within_myeloid,
    fill = final_state
  )
) +
  geom_col(width = 0.82) +
  scale_fill_manual(values = palette_values, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle(
    "Human MASH LT1-LT8 | final myeloid-state composition",
    "Biological samples ordered NAS0 -> NAS2 -> NAS4 -> NAS6"
  )

ggsave(
  file.path(FIGDIR, "01_final_state_composition_by_sample_v2.0.pdf"),
  p_comp,
  width = 10.5,
  height = 6.8
)

# 2) State trajectories: biological sample points + NAS mean
sample_plot <- state_counts
sample_plot$final_state <- factor(sample_plot$final_state, levels = STATE_ORDER)

mean_plot <- group_summary
mean_plot$final_state <- factor(mean_plot$final_state, levels = STATE_ORDER)

p_traj <- ggplot() +
  geom_point(
    data = sample_plot,
    aes(
      x = NAS,
      y = fraction_within_myeloid,
      color = final_state
    ),
    size = 2.4,
    alpha = 0.78,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  geom_line(
    data = mean_plot,
    aes(
      x = NAS,
      y = mean_fraction,
      color = final_state,
      group = final_state
    ),
    linewidth = 0.9
  ) +
  geom_point(
    data = mean_plot,
    aes(
      x = NAS,
      y = mean_fraction,
      color = final_state
    ),
    size = 3.2
  ) +
  scale_color_manual(values = palette_values, drop = FALSE) +
  scale_x_continuous(breaks = c(0, 2, 4, 6)) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.title = element_blank()
  ) +
  labs(
    x = "NAS",
    y = "Fraction within final myeloid cells"
  ) +
  ggtitle(
    "Human MASH | final myeloid-state trajectories across NAS",
    "Points = biological samples; larger connected points = stage mean"
  )

ggsave(
  file.path(FIGDIR, "02_final_state_fraction_trajectory_by_NAS_v2.0.pdf"),
  p_traj,
  width = 10.5,
  height = 7.2
)

# 3) Faceted state trajectories to avoid visual overlap
p_facet <- ggplot(
  sample_plot,
  aes(
    x = NAS,
    y = fraction_within_myeloid
  )
) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_line(
    data = mean_plot,
    aes(
      x = NAS,
      y = mean_fraction,
      group = 1
    ),
    linewidth = 0.9
  ) +
  geom_point(
    data = mean_plot,
    aes(
      x = NAS,
      y = mean_fraction
    ),
    size = 3.0
  ) +
  facet_wrap(~ final_state, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = c(0, 2, 4, 6)) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  theme_classic(base_size = 10) +
  labs(
    x = "NAS",
    y = "Fraction within final myeloid cells"
  ) +
  ggtitle(
    "Human MASH | sample-level myeloid-state fractions",
    "n=2 biological samples per NAS stage; descriptive"
  )

ggsave(
  file.path(FIGDIR, "03_final_state_fraction_facets_by_NAS_v2.0.pdf"),
  p_facet,
  width = 10.5,
  height = 10
)

# 4) Total final-myeloid captured-cell fraction
p_total <- ggplot(
  my_abundance,
  aes(
    x = NAS,
    y = final_myeloid_fraction_of_whole
  )
) +
  geom_point(size = 2.8) +
  geom_line(
    data = abundance_summary,
    aes(
      x = NAS,
      y = mean_fraction,
      group = 1
    ),
    linewidth = 0.9
  ) +
  geom_point(
    data = abundance_summary,
    aes(
      x = NAS,
      y = mean_fraction
    ),
    size = 3.3
  ) +
  scale_x_continuous(breaks = c(0, 2, 4, 6)) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x, 1), "%")
  ) +
  theme_classic(base_size = 11) +
  labs(
    x = "NAS",
    y = "Final myeloid cells / captured whole-liver cells"
  ) +
  ggtitle(
    "Human MASH | captured final-myeloid fraction",
    "Descriptive scRNA-seq cell composition; not absolute tissue abundance"
  )

ggsave(
  file.path(FIGDIR, "04_total_final_myeloid_fraction_of_whole_v2.0.pdf"),
  p_total,
  width = 8.5,
  height = 6.5
)

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

max_delta <- max(sensitivity_delta$absolute_delta, na.rm = TRUE)

summary_lines <- c(
  paste0("Human MASH LT1-LT8 final myeloid NAS composition ", VERSION),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "STATUS: SUCCESS",
  "",
  paste0("Final myeloid RDS: ", MYELOID_RDS),
  paste0("Canonical whole-liver RDS: ", WHOLE_RDS),
  "",
  "BIOLOGICAL DESIGN:",
  "NAS0: LT7, LT8",
  "NAS2: LT5, LT6",
  "NAS4: LT3, LT4",
  "NAS6: LT1, LT2",
  "n = 2 biological samples per NAS stage",
  "",
  "ANALYSES:",
  "- final myeloid fraction among captured whole-liver cells",
  "- final-state fraction within final myeloid cells",
  "- biological-sample summaries by NAS",
  "- exploratory sample-level Spearman trends",
  "- sensitivity excluding review clusters 13/14/16",
  "",
  "EXPLORATORY SPEARMAN:",
  paste(
    capture.output(print(trend_table, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  paste0(
    "Maximum absolute primary-vs-sensitivity state-fraction delta: ",
    sprintf("%.4f", max_delta)
  ),
  "",
  "INTERPRETATION:",
  "- Use biological samples, not cells, as the unit of inference.",
  "- With n=2/stage, figures and effect direction are more informative than p-values.",
  "- Whole-liver captured-cell fractions are compositional scRNA-seq measurements,",
  "  not absolute histologic cell abundance.",
  "- Do not infer mechanisms from composition alone.",
  "",
  "NEXT STEP:",
  "Review the NAS trajectories. Then select biologically relevant final states for",
  "sample-level pseudobulk DE/pathway analysis."
)

writeLines(
  summary_lines,
  file.path(TABDIR, "09_run_summary_v2.0.txt"),
  useBytes = TRUE
)

message("============================================================")
message("v2.0 COMPLETE")
message("Output: ", OUTDIR)
message("No RDS modified or saved.")
message("============================================================")
