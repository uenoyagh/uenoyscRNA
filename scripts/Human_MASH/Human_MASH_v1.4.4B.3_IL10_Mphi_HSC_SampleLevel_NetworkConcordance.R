#!/usr/bin/env Rscript

# ==============================================================================
# Human_MASH_v1.4.4B.3_IL10_Mphi_HSC_SampleLevel_NetworkConcordance.R
#
# Purpose
# -------
# Final sample-level descriptive integration of the mouse host mechanism:
#
#   engineered human graft / IL10 anchor
#       -> IL10R-competent primary macrophage pool
#       -> macrophage-state remodeling
#       -> Repair/Resolution-Mphi Pdgfb
#       -> Mphi->ECM-HSC PDGF communication / HSC PDGF program
#       -> ECM-activated / cycling HSC burden
#
# This script DOES NOT test or claim:
#   - IL10 protein secretion
#   - receptor occupancy
#   - pSTAT3 activation
#   - direct cell-autonomous IL10R -> Pdgfb repression
#   - PDGF as the sole mediator of fibrosis improvement
#
# Statistics:
#   descriptive only, n=2 Sham vs n=2 Tx.
#
# External delivery HDDs are not required.
# ==============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

VERSION <- "v1.4.4B.3"
ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------

B2_DIR <- file.path(
  ROOT,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4B.2",
  "analysis/v1.4.4B.2_graft_receiver_weighted_Mphi_HSC_network_integration"
)

B2_POOLED <- file.path(
  B2_DIR,
  "tables/v1.4.4B.2_03_pooled_receiver_competence_by_sample.csv"
)

B2_IL10_STATE <- file.path(
  B2_DIR,
  "tables/v1.4.4B.2_02_receiver_weighted_state_by_sample.csv"
)

C21_ABUND <- file.path(
  ROOT,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.3C.2.1",
  "analysis/v1.4.3C.2.1_mouse_Mphi_IL10_responder_mechanistic_audit",
  "tables/v1.4.3C.2.1_03_class_abundance_by_sample.tsv"
)

C22_PDGFB <- file.path(
  ROOT,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.3C.2.2",
  "analysis/v1.4.3C.2.2_mouse_Mphi_IL10R_PDGFB_bridge_audit",
  "tables/v1.4.3C.2.2_03_RepairResolution_unsplit_Pdgfb_by_sample.tsv"
)

C22_BOUNDARY <- file.path(
  ROOT,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.3C.2.2",
  "analysis/v1.4.3C.2.2_mouse_Mphi_IL10R_PDGFB_bridge_audit",
  "tables/v1.4.3C.2.2_11_interpretation.tsv"
)

V636_EVIDENCE <- file.path(
  ROOT,
  "Mouse_MASH_Interaction/PDGF_FN1_interaction_closure_v6.3.6",
  "Tables/05_PDGF_FN1_evidence_matrix_v6.3.6.csv"
)

V636_HSC <- file.path(
  ROOT,
  "Mouse_MASH_Interaction/PDGF_FN1_interaction_closure_v6.3.6",
  "Tables/03_ECM_HSC_burden_replicate_summary_v6.3.6.csv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4B.3",
  "analysis/v1.4.4B.3_IL10_Mphi_HSC_sample_level_network_concordance"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
LOG <- file.path(OUT, "logs")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)

sink(file.path(LOG, paste0(VERSION, "_run_log.txt")), split = TRUE)
on.exit(sink(), add = TRUE)

cat("============================================================================\n")
cat("Human MASH ", VERSION, "\n", sep = "")
cat("IL10 -> Mphi -> HSC sample-level network concordance\n")
cat("============================================================================\n\n")

for (f in c(
  B2_POOLED, B2_IL10_STATE, C21_ABUND, C22_PDGFB,
  C22_BOUNDARY, V636_EVIDENCE, V636_HSC
)) {
  if (!file.exists(f)) stop("Required input missing: ", f)
}

SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")
COND <- c(Sham1 = "Sham", Sham20 = "Sham", Tx17 = "Tx", Tx5 = "Tx")

PRIMARY_STATES <- c(
  "Inflammatory-Mphi",
  "Anti-inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

relation4 <- function(s1, s2, t1, t2) {
  z <- c(s1, s2, t1, t2)
  if (any(!is.finite(z))) return("Incomplete_four_sample_data")
  if (min(c(t1, t2)) > max(c(s1, s2))) {
    return("Tx_higher_complete_replicate_range_separation")
  }
  if (max(c(t1, t2)) < min(c(s1, s2))) {
    return("Tx_lower_complete_replicate_range_separation")
  }
  "Replicate_ranges_overlap"
}

direction_call <- function(delta, tol = 1e-12) {
  if (!is.finite(delta)) return("NA")
  if (delta > tol) return("Tx_up")
  if (delta < -tol) return("Tx_down")
  "No_change"
}

extract_sample_values <- function(df, row_idx) {
  setNames(
    as.numeric(df[row_idx, SAMPLES, drop = TRUE]),
    SAMPLES
  )
}

add_metric <- function(store, idx, layer, metric, values,
                       expected_direction, source, boundary) {
  sham_mean <- mean(values[c("Sham1", "Sham20")])
  tx_mean <- mean(values[c("Tx17", "Tx5")])
  delta <- tx_mean - sham_mean
  obs <- direction_call(delta)

  concord <- if (expected_direction == "Tx_up") {
    obs == "Tx_up"
  } else if (expected_direction == "Tx_down") {
    obs == "Tx_down"
  } else {
    NA
  }

  store[[idx]] <- data.frame(
    layer = layer,
    metric = metric,
    Sham1 = values["Sham1"],
    Sham20 = values["Sham20"],
    Tx17 = values["Tx17"],
    Tx5 = values["Tx5"],
    Sham_mean = sham_mean,
    Tx_mean = tx_mean,
    Tx_minus_Sham = delta,
    observed_direction = obs,
    expected_model_direction = expected_direction,
    model_direction_concordant = concord,
    replicate_relation = relation4(
      values["Sham1"], values["Sham20"],
      values["Tx17"], values["Tx5"]
    ),
    source = source,
    boundary = boundary,
    stringsAsFactors = FALSE
  )

  store
}

# ------------------------------------------------------------------------------
# 1. Build sample-level macrophage matrix
# ------------------------------------------------------------------------------

pooled <- read.csv(B2_POOLED, stringsAsFactors = FALSE, check.names = FALSE)
il10_pool <- pooled[pooled$human_ligand == "IL10", , drop = FALSE]

if (nrow(il10_pool) != 4L) {
  stop("Expected four IL10 pooled-receiver sample rows.")
}

il10_pool <- il10_pool[match(SAMPLES, il10_pool$sample), , drop = FALSE]
if (any(is.na(il10_pool$sample))) stop("IL10 sample matching failed.")

il10_state <- read.csv(B2_IL10_STATE, stringsAsFactors = FALSE, check.names = FALSE)
il10_state <- il10_state[
  il10_state$human_ligand == "IL10" &
    il10_state$sample %in% SAMPLES &
    il10_state$macrophage_state %in% PRIMARY_STATES,
  ,
  drop = FALSE
]

if (nrow(il10_state) != 20L) {
  stop("Expected 20 IL10 sample x state rows.")
}

c21 <- read.delim(C21_ABUND, stringsAsFactors = FALSE, check.names = FALSE)
c21 <- c21[
  c21$sample %in% SAMPLES &
    c21$macrophage_class %in% PRIMARY_STATES,
  ,
  drop = FALSE
]

if (nrow(c21) != 20L) stop("Expected 20 C.2.1 primary-state abundance rows.")

c22 <- read.delim(C22_PDGFB, stringsAsFactors = FALSE, check.names = FALSE)
c22 <- c22[c22$sample %in% SAMPLES, , drop = FALSE]
if (nrow(c22) != 4L) stop("Expected four Repair/Resolution Pdgfb sample rows.")

sample_matrix <- data.frame(
  sample = SAMPLES,
  condition = unname(COND[SAMPLES]),
  stringsAsFactors = FALSE
)

sample_matrix$IL10R_dual_positive_fraction_primary5_Mphi <-
  il10_pool$pooled_receptor_positive_fraction

# state fractions
for (st in PRIMARY_STATES) {
  nm <- gsub("[^A-Za-z0-9]+", "_", st)
  z <- c21[c21$macrophage_class == st, , drop = FALSE]
  z <- z[match(SAMPLES, z$sample), , drop = FALSE]

  sample_matrix[[paste0("Mphi_fraction_", nm)]] <- z$fraction_total_Mphi
}

# receiver-weighted IL10R pool contribution by state
for (st in PRIMARY_STATES) {
  nm <- gsub("[^A-Za-z0-9]+", "_", st)
  z <- il10_state[il10_state$macrophage_state == st, , drop = FALSE]
  z <- z[match(SAMPLES, z$sample), , drop = FALSE]

  sample_matrix[[paste0("IL10R_weighted_fraction_", nm)]] <-
    z$receiver_weighted_fraction_of_primary5_Mphi
}

# Repair/Resolution Pdgfb
c22 <- c22[match(SAMPLES, c22$sample), , drop = FALSE]
sample_matrix$RepairResolution_Pdgfb_mean_per_cell_CP10k <-
  c22$Pdgfb_mean_per_cell_CP10k
sample_matrix$RepairResolution_Pdgfb_pseudobulk_CP10k <-
  c22$Pdgfb_pseudobulk_CP10k

write.csv(
  sample_matrix,
  file.path(TAB, paste0(VERSION, "_01_sample_level_macrophage_matrix.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 2. Add v6.3.6 sample-level PDGF/HSC metrics
# ------------------------------------------------------------------------------

ev <- read.csv(V636_EVIDENCE, stringsAsFactors = FALSE, check.names = FALSE)
hsc <- read.csv(V636_HSC, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(SAMPLES %in% colnames(ev))) {
  stop("v6.3.6 evidence matrix lacks sample columns.")
}
if (!all(SAMPLES %in% colnames(hsc))) {
  stop("v6.3.6 HSC table lacks sample columns.")
}

# Identify exact rows without silently choosing multiple matches.
pick_one <- function(df, pattern, label) {
  idx <- grep(pattern, df$metric, ignore.case = TRUE)
  if (length(idx) != 1L) {
    stop(
      "Could not uniquely identify ", label,
      ". Matching rows = ", length(idx),
      "; pattern = ", pattern
    )
  }
  idx
}

idx_cellchat <- pick_one(
  ev,
  "CellChat PDGF supported strength",
  "CellChat PDGF Repair/Resolution-Mphi -> ECM-HSC"
)

# The PDGF mitogenic-module row is expected in the evidence matrix.
idx_mitogenic_candidates <- grep(
  "mitogenic|PDGF.*module|module.*PDGF",
  ev$metric,
  ignore.case = TRUE
)

if (length(idx_mitogenic_candidates) == 1L) {
  idx_mitogenic <- idx_mitogenic_candidates
} else {
  # Print candidate metrics to make failure auditable rather than guessing.
  cat("Available PDGF metrics:\n")
  print(ev[ev$axis == "PDGF", c("metric"), drop = FALSE])
  stop(
    "Could not uniquely identify ECM-HSC PDGF mitogenic-module row. ",
    "Matching rows = ", length(idx_mitogenic_candidates)
  )
}

idx_ecm_frac <- pick_one(
  hsc,
  "^ECM-activated HSC fraction within HSC3$",
  "ECM-activated HSC fraction"
)

idx_cycling_frac <- pick_one(
  hsc,
  "^Cycling-like fraction within ECM-activated HSC$",
  "cycling fraction within ECM-HSC"
)

idx_cycling_burden_candidates <- grep(
  "burden",
  hsc$metric,
  ignore.case = TRUE
)
if (length(idx_cycling_burden_candidates) == 1L) {
  idx_cycling_burden <- idx_cycling_burden_candidates
} else {
  cat("Available HSC burden metrics:\n")
  print(hsc[, "metric", drop = FALSE])
  stop(
    "Could not uniquely identify cycling ECM-HSC burden row. ",
    "Matching rows = ", length(idx_cycling_burden_candidates)
  )
}

sample_matrix$CellChat_PDGF_RepairResolution_to_ECMHSC <-
  unname(extract_sample_values(ev, idx_cellchat)[SAMPLES])

sample_matrix$ECMHSC_PDGF_mitogenic_module <-
  unname(extract_sample_values(ev, idx_mitogenic)[SAMPLES])

sample_matrix$ECM_HSC_fraction_within_HSC3 <-
  unname(extract_sample_values(hsc, idx_ecm_frac)[SAMPLES])

sample_matrix$Cycling_fraction_within_ECM_HSC <-
  unname(extract_sample_values(hsc, idx_cycling_frac)[SAMPLES])

sample_matrix$Cycling_ECM_HSC_burden <-
  unname(extract_sample_values(hsc, idx_cycling_burden)[SAMPLES])

write.csv(
  sample_matrix,
  file.path(TAB, paste0(VERSION, "_02_sample_level_full_network_matrix.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 3. Directional evidence ladder
# ------------------------------------------------------------------------------

metrics <- list()
mi <- 1L

vals <- setNames(
  sample_matrix$IL10R_dual_positive_fraction_primary5_Mphi,
  sample_matrix$sample
)
metrics <- add_metric(
  metrics, mi,
  "Host_receiver",
  "IL10R-dual-positive fraction across primary five Mphi states",
  vals,
  "Tx_up",
  "v1.4.4B.2",
  "Receptor-transcript competence only; not receptor activation or pSTAT3."
)
mi <- mi + 1L

# Anti-inflammatory Mphi
col_ai <- "Mphi_fraction_Anti_inflammatory_Mphi"
vals <- setNames(sample_matrix[[col_ai]], sample_matrix$sample)
metrics <- add_metric(
  metrics, mi,
  "Mphi_state_remodeling",
  "Anti-inflammatory-Mphi fraction",
  vals,
  "Tx_up",
  "v1.4.3C.2.1 / v1.4.4B.3",
  "State abundance association; not a direct IL10 effect."
)
mi <- mi + 1L

# Repair/Resolution Mphi
col_rr <- "Mphi_fraction_Repair_Resolution_Mphi"
vals <- setNames(sample_matrix[[col_rr]], sample_matrix$sample)
metrics <- add_metric(
  metrics, mi,
  "Mphi_state_remodeling",
  "Repair/Resolution-Mphi fraction",
  vals,
  "Tx_up",
  "v1.4.3C.2.1 / v1.4.4B.3",
  "State abundance association; not a direct IL10 effect."
)
mi <- mi + 1L

# Inflammatory Mphi
col_inf <- "Mphi_fraction_Inflammatory_Mphi"
vals <- setNames(sample_matrix[[col_inf]], sample_matrix$sample)
metrics <- add_metric(
  metrics, mi,
  "Mphi_state_remodeling",
  "Inflammatory-Mphi fraction",
  vals,
  "Tx_down",
  "v1.4.3C.2.1 / v1.4.4B.3",
  "State abundance association."
)
mi <- mi + 1L

# Repair Pdgfb per cell
vals <- setNames(
  sample_matrix$RepairResolution_Pdgfb_mean_per_cell_CP10k,
  sample_matrix$sample
)
metrics <- add_metric(
  metrics, mi,
  "Mphi_effector",
  "Repair/Resolution-Mphi mean Pdgfb per-cell CP10k",
  vals,
  "Tx_down",
  "v1.4.3C.2.2 / v6.3.6",
  "Direct IL10R-status -> Pdgfb repression was not demonstrated."
)
mi <- mi + 1L

# CellChat PDGF
vals <- setNames(
  sample_matrix$CellChat_PDGF_RepairResolution_to_ECMHSC,
  sample_matrix$sample
)
metrics <- add_metric(
  metrics, mi,
  "Mphi_to_HSC",
  "CellChat PDGF: Repair/Resolution-Mphi -> ECM-HSC",
  vals,
  "Tx_down",
  "v6.3.6",
  "Selective communication estimate; not proof of ligand flux."
)
mi <- mi + 1L

# HSC PDGF module
vals <- setNames(
  sample_matrix$ECMHSC_PDGF_mitogenic_module,
  sample_matrix$sample
)
metrics <- add_metric(
  metrics, mi,
  "HSC_receiver_program",
  "ECM-HSC PDGF mitogenic module",
  vals,
  "Tx_down",
  "v6.3.6",
  "Transcriptomic receiver program; not a direct biochemical activity assay."
)
mi <- mi + 1L

# ECM HSC fraction
vals <- setNames(
  sample_matrix$ECM_HSC_fraction_within_HSC3,
  sample_matrix$sample
)
metrics <- add_metric(
  metrics, mi,
  "HSC_state",
  "ECM-activated HSC fraction within HSC3",
  vals,
  "Tx_down",
  "v6.3.6",
  "Association with treatment; PDGF is not assumed to be the sole cause."
)
mi <- mi + 1L

# Cycling fraction
vals <- setNames(
  sample_matrix$Cycling_fraction_within_ECM_HSC,
  sample_matrix$sample
)
metrics <- add_metric(
  metrics, mi,
  "HSC_state",
  "Cycling-like fraction within ECM-activated HSC",
  vals,
  "Tx_down",
  "v6.3.6",
  "Association with treatment."
)
mi <- mi + 1L

# Cycling burden
vals <- setNames(
  sample_matrix$Cycling_ECM_HSC_burden,
  sample_matrix$sample
)
metrics <- add_metric(
  metrics, mi,
  "HSC_state",
  "Cycling ECM-HSC burden within HSC3",
  vals,
  "Tx_down",
  "v6.3.6",
  "Association with treatment."
)

ladder <- do.call(rbind, metrics)

write.csv(
  ladder,
  file.path(TAB, paste0(VERSION, "_03_directional_network_evidence_ladder.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Network-level concordance summary
# ------------------------------------------------------------------------------

summary <- data.frame(
  n_metrics = nrow(ladder),
  n_model_direction_concordant = sum(
    ladder$model_direction_concordant %in% TRUE,
    na.rm = TRUE
  ),
  fraction_model_direction_concordant = mean(
    ladder$model_direction_concordant,
    na.rm = TRUE
  ),
  n_complete_replicate_range_separation = sum(
    ladder$replicate_relation %in% c(
      "Tx_higher_complete_replicate_range_separation",
      "Tx_lower_complete_replicate_range_separation"
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary,
  file.path(TAB, paste0(VERSION, "_04_network_concordance_summary.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Explicit direct-chain negative boundary
# ------------------------------------------------------------------------------

c22_boundary <- read.delim(
  C22_BOUNDARY,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.csv(
  c22_boundary,
  file.path(TAB, paste0(VERSION, "_05_direct_IL10R_Pdgfb_negative_boundary.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Descriptive sample-level correlation matrix
#
# This is deliberately secondary/exploratory. n=4 is too small for inference.
# We report Spearman rho only as a descriptive consistency check.
# ------------------------------------------------------------------------------

corr_vars <- c(
  "IL10R_dual_positive_fraction_primary5_Mphi",
  "Mphi_fraction_Anti_inflammatory_Mphi",
  "Mphi_fraction_Repair_Resolution_Mphi",
  "Mphi_fraction_Inflammatory_Mphi",
  "RepairResolution_Pdgfb_mean_per_cell_CP10k",
  "CellChat_PDGF_RepairResolution_to_ECMHSC",
  "ECMHSC_PDGF_mitogenic_module",
  "ECM_HSC_fraction_within_HSC3",
  "Cycling_ECM_HSC_burden"
)

corr <- cor(
  sample_matrix[, corr_vars, drop = FALSE],
  method = "spearman",
  use = "pairwise.complete.obs"
)

write.csv(
  corr,
  file.path(TAB, paste0(VERSION, "_06_exploratory_Spearman_rho_n4.csv"))
)

corr_note <- data.frame(
  item = c(
    "sample_n",
    "inference",
    "use"
  ),
  value = c(
    "4 total biological samples",
    "Descriptive only; no p-values and no inferential interpretation",
    "Secondary consistency check only; condition separation can dominate correlations"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  corr_note,
  file.path(TAB, paste0(VERSION, "_07_exploratory_correlation_boundary.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Publication-oriented evidence-ladder figure
# ------------------------------------------------------------------------------

plot_ladder <- ladder
plot_ladder$metric <- factor(
  plot_ladder$metric,
  levels = rev(plot_ladder$metric)
)

plot_ladder$direction_score <- ifelse(
  plot_ladder$observed_direction == "Tx_up",
  1,
  ifelse(plot_ladder$observed_direction == "Tx_down", -1, 0)
)

plot_ladder$replicate_complete <- grepl(
  "complete_replicate_range_separation",
  plot_ladder$replicate_relation,
  fixed = TRUE
)

p1 <- ggplot(
  plot_ladder,
  aes(
    x = direction_score,
    y = metric,
    shape = replicate_complete
  )
) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_point(size = 4) +
  scale_x_continuous(
    breaks = c(-1, 0, 1),
    labels = c("Tx lower", "No change", "Tx higher"),
    limits = c(-1.25, 1.25)
  ) +
  labs(
    title = "IL10-responsive macrophage to HSC network: directional evidence ladder",
    subtitle = "Observed Tx-vs-Sham direction across linked mechanistic layers",
    x = NULL,
    y = NULL,
    shape = "Complete replicate-\nrange separation",
    caption = paste(
      "Descriptive n=2/group. Directional concordance supports a network-level model;",
      "it does not establish a direct intracellular IL10R->Pdgfb chain."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG, paste0(VERSION, "_01_directional_network_evidence_ladder.pdf")),
  p1,
  width = 10.5,
  height = 7.5
)

ggsave(
  file.path(FIG, paste0(VERSION, "_01_directional_network_evidence_ladder.png")),
  p1,
  width = 10.5,
  height = 7.5,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 8. Four-sample standardized trajectory figure
#
# Per metric z-scaling is used only for visualization across unlike units.
# For expected Tx-down metrics, sign is reversed so upward = model-concordant.
# ------------------------------------------------------------------------------

traj <- ladder
rows <- list()
ri <- 1L

for (i in seq_len(nrow(traj))) {
  vals <- as.numeric(traj[i, SAMPLES, drop = TRUE])

  if (sd(vals) == 0 || !is.finite(sd(vals))) {
    z <- rep(0, length(vals))
  } else {
    z <- as.numeric(scale(vals))
  }

  if (traj$expected_model_direction[i] == "Tx_down") {
    z <- -z
  }

  for (j in seq_along(SAMPLES)) {
    rows[[ri]] <- data.frame(
      sample = SAMPLES[j],
      condition = unname(COND[SAMPLES[j]]),
      layer = traj$layer[i],
      metric = traj$metric[i],
      direction_aligned_z = z[j],
      stringsAsFactors = FALSE
    )
    ri <- ri + 1L
  }
}

traj_long <- do.call(rbind, rows)

write.csv(
  traj_long,
  file.path(TAB, paste0(VERSION, "_08_direction_aligned_sample_trajectory_long.csv")),
  row.names = FALSE
)

# Layer-average is visualization only.
layer_avg <- aggregate(
  direction_aligned_z ~ sample + condition + layer,
  data = traj_long,
  FUN = mean
)

layer_order <- unique(ladder$layer)
layer_avg$layer <- factor(layer_avg$layer, levels = layer_order)
layer_avg$sample <- factor(layer_avg$sample, levels = SAMPLES)

p2 <- ggplot(
  layer_avg,
  aes(
    x = layer,
    y = direction_aligned_z,
    group = sample,
    linetype = condition
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ sample, ncol = 2) +
  labs(
    title = "Four-sample direction-aligned network trajectory",
    subtitle = "Per-metric z-scaled visualization; upward is aligned with the proposed treatment model",
    x = NULL,
    y = "Direction-aligned z (visualization only)",
    caption = "Not an inferential score; metrics have different biological units."
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG, paste0(VERSION, "_02_four_sample_direction_aligned_trajectory.pdf")),
  p2,
  width = 11,
  height = 8
)

ggsave(
  file.path(FIG, paste0(VERSION, "_02_four_sample_direction_aligned_trajectory.png")),
  p2,
  width = 11,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 9. Interpretation / provenance
# ------------------------------------------------------------------------------

boundary <- data.frame(
  item = c(
    "primary_claim",
    "IL10R_boundary",
    "pSTAT3_boundary",
    "direct_Pdgfb_boundary",
    "PDGF_boundary",
    "statistics_boundary",
    "correlation_boundary",
    "fibrosis_boundary"
  ),
  value = c(
    "Network-level convergence: Tx is associated with expansion/remodeling of an IL10R-competent macrophage pool, selective attenuation of macrophage-to-HSC PDGF signaling, and reduced activated/cycling HSC burden",
    "Il10ra+Il10rb transcript co-detection indicates receptor-transcript competence only",
    "Stat3 transcript measurements do not establish STAT3 phosphorylation/activation",
    "v1.4.3C.2.2 did not demonstrate a Tx-specific replicate-consistent direct IL10R-status -> Pdgfb relationship",
    "PDGF is interpreted as a selective mitogenic/proliferative arm, not global HSC suppression and not the sole mediator",
    "n=2 biological samples per condition; descriptive only",
    "n=4 Spearman correlations are exploratory consistency checks only and must not be used as inferential evidence",
    "Reduced HSC fibrogenic activity is mechanistically coherent with fibrosis improvement, but scRNA-seq alone does not prove the complete causal route to histologic fibrosis regression"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  boundary,
  file.path(TAB, paste0(VERSION, "_09_interpretation_boundaries.csv")),
  row.names = FALSE
)

prov <- data.frame(
  key = c(
    "version",
    "B2_pooled_receiver",
    "B2_state_receiver",
    "C21_state_abundance",
    "C22_Pdgfb",
    "C22_negative_boundary",
    "v636_evidence",
    "v636_HSC_burden"
  ),
  value = c(
    VERSION,
    B2_POOLED,
    B2_IL10_STATE,
    C21_ABUND,
    C22_PDGFB,
    C22_BOUNDARY,
    V636_EVIDENCE,
    V636_HSC
  ),
  stringsAsFactors = FALSE
)

write.csv(
  prov,
  file.path(TAB, paste0(VERSION, "_10_provenance.csv")),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(LOG, paste0(VERSION, "_sessionInfo.txt"))
)

# ------------------------------------------------------------------------------
# 10. Console summary
# ------------------------------------------------------------------------------

cat("============================================================================\n")
cat("NETWORK DIRECTIONAL EVIDENCE LADDER\n")
cat("============================================================================\n")

print(
  ladder[, c(
    "layer",
    "metric",
    "Sham_mean",
    "Tx_mean",
    "Tx_minus_Sham",
    "observed_direction",
    "expected_model_direction",
    "model_direction_concordant",
    "replicate_relation"
  )],
  row.names = FALSE
)

cat("\n============================================================================\n")
cat("NETWORK CONCORDANCE SUMMARY\n")
cat("============================================================================\n")
print(summary, row.names = FALSE)

cat("\n============================================================================\n")
cat(VERSION, " COMPLETE\n")
cat("============================================================================\n")
cat("Output:", OUT, "\n")
cat("Sample-level IL10 -> Mphi -> HSC matrix: COMPLETE\n")
cat("Directional evidence ladder: COMPLETE\n")
cat("Exploratory n=4 correlation matrix: COMPLETE (descriptive only)\n")
cat("Direct IL10R -> Pdgfb negative boundary preserved: YES\n")
cat("No pSTAT3 claim made: YES\n")
cat("External delivery HDDs required: NO\n")
