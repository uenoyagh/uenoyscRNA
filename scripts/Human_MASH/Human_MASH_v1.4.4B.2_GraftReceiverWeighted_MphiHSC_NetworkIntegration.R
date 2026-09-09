#!/usr/bin/env Rscript

# ==============================================================================
# Human_MASH_v1.4.4B.2_GraftReceiverWeighted_MphiHSC_NetworkIntegration.R
#
# Purpose
# -------
# Build the next network-level bridge:
#
#   latest human graft ligands
#       -> mouse Mphi receptor-transcript competence
#       -> receptor-competent Mphi pool after state-composition remodeling
#       -> existing Mphi-to-HSC PDGF evidence
#       -> activated/cycling HSC burden
#
# This analysis does NOT infer receptor activation and does NOT claim a direct
# cell-autonomous IL10R -> Pdgfb causal chain.
#
# Primary macrophage universe:
#   five final v4.14.5 states, excluding "Other", matching v1.4.3C.2.1.
#
# Statistics:
#   descriptive only, n=2 Sham and n=2 Tx.
# ==============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

VERSION <- "v1.4.4B.2"
ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

B1 <- file.path(
  ROOT,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4B.1",
  "analysis/v1.4.4B.1_graft_ligand_mouse_Mphi5_receiver_bridge",
  "tables/v1.4.4B.1_07_integrated_graft_ligand_mouse_Mphi5_bridge.csv"
)

C21_ABUND <- file.path(
  ROOT,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.3C.2.1",
  "analysis/v1.4.3C.2.1_mouse_Mphi_IL10_responder_mechanistic_audit",
  "tables/v1.4.3C.2.1_03_class_abundance_by_sample.tsv"
)

C22_BOUNDARY <- file.path(
  ROOT,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.3C.2.2",
  "analysis/v1.4.3C.2.2_mouse_Mphi_IL10R_PDGFB_bridge_audit",
  "tables/v1.4.3C.2.2_11_interpretation.tsv"
)

HSC_EVIDENCE <- file.path(
  ROOT,
  "Mouse_MASH_Interaction/PDGF_FN1_interaction_closure_v6.3.6",
  "Tables/05_PDGF_FN1_evidence_matrix_v6.3.6.csv"
)

HSC_BURDEN <- file.path(
  ROOT,
  "Mouse_MASH_Interaction/PDGF_FN1_interaction_closure_v6.3.6",
  "Tables/03_ECM_HSC_burden_replicate_summary_v6.3.6.csv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS/Human_MASH_12sample_Context_v1.4.4B.2",
  "analysis/v1.4.4B.2_graft_receiver_weighted_Mphi_HSC_network_integration"
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
cat("Graft receiver-weighted Mphi -> HSC network integration\n")
cat("============================================================================\n\n")

for (f in c(B1, C21_ABUND, C22_BOUNDARY, HSC_EVIDENCE, HSC_BURDEN)) {
  if (!file.exists(f)) stop("Required input missing: ", f)
}

PRIMARY_STATES <- c(
  "Inflammatory-Mphi",
  "Anti-inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")
CONDITION <- c(Sham1 = "Sham", Sham20 = "Sham", Tx17 = "Tx", Tx5 = "Tx")

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

# ------------------------------------------------------------------------------
# 1. Load and validate B.1
# ------------------------------------------------------------------------------

b1 <- read.csv(B1, stringsAsFactors = FALSE, check.names = FALSE)

req <- c(
  "human_ligand",
  "receiver_complex",
  "macrophage_state",
  "Sham1_n_cells",
  "Sham20_n_cells",
  "Tx17_n_cells",
  "Tx5_n_cells",
  "Sham1_receptor_positive_cells",
  "Sham20_receptor_positive_cells",
  "Tx17_receptor_positive_cells",
  "Tx5_receptor_positive_cells",
  "Sham1_receptor_positive_fraction",
  "Sham20_receptor_positive_fraction",
  "Tx17_receptor_positive_fraction",
  "Tx5_receptor_positive_fraction",
  "human_transition_call",
  "human_Tx_vs_SfLB_log2FC_per_cell_descriptive"
)

missing_req <- setdiff(req, colnames(b1))
if (length(missing_req)) {
  stop("B.1 schema mismatch: ", paste(missing_req, collapse = ", "))
}

if (nrow(b1) != 50L) stop("Expected 50 B.1 rows, observed ", nrow(b1))
if (length(unique(b1$human_ligand)) != 10L) stop("Expected 10 ligands.")
if (!setequal(unique(b1$macrophage_state), PRIMARY_STATES)) {
  stop("B.1 state universe does not match frozen primary five states.")
}

# Cell denominators must be invariant across ligands for a sample/state.
denom_long <- list()
di <- 1L

for (ss in SAMPLES) {
  ncolx <- paste0(ss, "_n_cells")
  for (st in PRIMARY_STATES) {
    z <- unique(b1[b1$macrophage_state == st, ncolx])
    if (length(z) != 1L) {
      stop("B.1 denominator is not invariant: ", ss, " / ", st)
    }
    denom_long[[di]] <- data.frame(
      sample = ss,
      condition = unname(CONDITION[ss]),
      macrophage_state = st,
      n_cells = as.numeric(z),
      stringsAsFactors = FALSE
    )
    di <- di + 1L
  }
}

denom <- do.call(rbind, denom_long)

sample_totals <- aggregate(
  n_cells ~ sample + condition,
  data = denom,
  FUN = sum
)
names(sample_totals)[names(sample_totals) == "n_cells"] <- "primary5_Mphi_total"

# ------------------------------------------------------------------------------
# 2. Validate against C.2.1 frozen class abundance
# ------------------------------------------------------------------------------

c21 <- read.delim(C21_ABUND, stringsAsFactors = FALSE, check.names = FALSE)

c21 <- c21[
  c21$sample %in% SAMPLES &
    c21$macrophage_class %in% PRIMARY_STATES,
  ,
  drop = FALSE
]

val <- merge(
  denom,
  c21[, c(
    "sample",
    "condition",
    "macrophage_class",
    "n_cells",
    "sample_total_Mphi",
    "fraction_total_Mphi"
  )],
  by.x = c("sample", "condition", "macrophage_state"),
  by.y = c("sample", "condition", "macrophage_class"),
  all.x = TRUE,
  suffixes = c("_B1", "_C21"),
  sort = FALSE
)

val$n_cell_abs_difference <- abs(val$n_cells_B1 - val$n_cells_C21)

if (any(!is.finite(val$n_cell_abs_difference)) ||
    any(val$n_cell_abs_difference != 0)) {
  stop("B.1 vs C.2.1 state denominator validation FAILED.")
}

totcheck <- merge(
  sample_totals,
  unique(c21[, c("sample", "sample_total_Mphi")]),
  by = "sample",
  all.x = TRUE
)
totcheck$total_abs_difference <- abs(
  totcheck$primary5_Mphi_total - totcheck$sample_total_Mphi
)

if (any(totcheck$total_abs_difference != 0)) {
  stop("Primary five-state total differs from C.2.1 sample_total_Mphi.")
}

write.csv(
  val,
  file.path(TAB, paste0(VERSION, "_01_B1_vs_C21_denominator_validation.csv")),
  row.names = FALSE
)

cat("B.1 vs C.2.1 primary 5-state denominator validation: PASS\n")
print(totcheck, row.names = FALSE)
cat("\n")

# ------------------------------------------------------------------------------
# 3. Receiver-weighted abundance by ligand x state x sample
#
# receiver_weighted_fraction_of_primary5_Mphi
# = receptor-positive cells in state / all primary-five Mphi in sample
#
# This combines:
#   state abundance * within-state receptor-positive fraction.
# ------------------------------------------------------------------------------

rw_rows <- list()
ri <- 1L

for (i in seq_len(nrow(b1))) {
  r <- b1[i, ]

  for (ss in SAMPLES) {
    n_state <- as.numeric(r[[paste0(ss, "_n_cells")]])
    n_pos <- as.numeric(r[[paste0(ss, "_receptor_positive_cells")]])
    pf <- as.numeric(r[[paste0(ss, "_receptor_positive_fraction")]])

    total <- sample_totals$primary5_Mphi_total[
      sample_totals$sample == ss
    ]

    if (length(total) != 1L) stop("Sample total lookup failed: ", ss)

    state_fraction <- n_state / total
    receiver_weighted <- n_pos / total

    identity_error <- abs(receiver_weighted - state_fraction * pf)

    rw_rows[[ri]] <- data.frame(
      human_ligand = r$human_ligand,
      receiver_complex = r$receiver_complex,
      human_transition_call = r$human_transition_call,
      human_Tx_vs_SfLB_log2FC_per_cell_descriptive =
        r$human_Tx_vs_SfLB_log2FC_per_cell_descriptive,
      sample = ss,
      condition = unname(CONDITION[ss]),
      macrophage_state = r$macrophage_state,
      sample_primary5_Mphi_total = total,
      state_n_cells = n_state,
      state_fraction_of_primary5_Mphi = state_fraction,
      receptor_positive_cells = n_pos,
      receptor_positive_fraction_within_state = pf,
      receiver_weighted_fraction_of_primary5_Mphi = receiver_weighted,
      receiver_weighted_cells_per_1000_primary5_Mphi = receiver_weighted * 1000,
      identity_error = identity_error,
      stringsAsFactors = FALSE
    )
    ri <- ri + 1L
  }
}

rw <- do.call(rbind, rw_rows)

if (max(rw$identity_error, na.rm = TRUE) > 1e-12) {
  stop("Receiver-weighted identity check failed.")
}

write.csv(
  rw,
  file.path(TAB, paste0(VERSION, "_02_receiver_weighted_state_by_sample.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Pooled receptor-competent fraction across the primary five-state Mphi pool
# ------------------------------------------------------------------------------

pooled_sample <- aggregate(
  cbind(receptor_positive_cells, state_n_cells) ~
    human_ligand + receiver_complex + sample + condition,
  data = rw,
  FUN = sum
)

pooled_sample$pooled_receptor_positive_fraction <-
  pooled_sample$receptor_positive_cells / pooled_sample$state_n_cells

pooled_sample$pooled_receptor_positive_cells_per_1000_Mphi <-
  pooled_sample$pooled_receptor_positive_fraction * 1000

write.csv(
  pooled_sample,
  file.path(TAB, paste0(VERSION, "_03_pooled_receiver_competence_by_sample.csv")),
  row.names = FALSE
)

pooled_rows <- list()
pi <- 1L

for (lig in unique(b1$human_ligand)) {
  z <- pooled_sample[pooled_sample$human_ligand == lig, , drop = FALSE]

  getv <- function(ss) {
    x <- z$pooled_receptor_positive_fraction[z$sample == ss]
    if (length(x) != 1L) return(NA_real_)
    x
  }

  s1 <- getv("Sham1")
  s2 <- getv("Sham20")
  t1 <- getv("Tx17")
  t2 <- getv("Tx5")

  h <- b1[b1$human_ligand == lig, , drop = FALSE][1, ]

  pooled_rows[[pi]] <- data.frame(
    human_ligand = lig,
    receiver_complex = h$receiver_complex,
    human_transition_call = h$human_transition_call,
    human_Tx_vs_SfLB_log2FC_per_cell_descriptive =
      h$human_Tx_vs_SfLB_log2FC_per_cell_descriptive,
    Sham1_pooled_fraction = s1,
    Sham20_pooled_fraction = s2,
    Tx17_pooled_fraction = t1,
    Tx5_pooled_fraction = t2,
    Sham_mean_pooled_fraction = mean(c(s1, s2)),
    Tx_mean_pooled_fraction = mean(c(t1, t2)),
    Tx_minus_Sham_pooled_fraction =
      mean(c(t1, t2)) - mean(c(s1, s2)),
    pooled_replicate_relation = relation4(s1, s2, t1, t2),
    inference =
      "descriptive_receiver_competence_of_primary5_Mphi_pool_n2_per_condition",
    stringsAsFactors = FALSE
  )
  pi <- pi + 1L
}

pooled_summary <- do.call(rbind, pooled_rows)

write.csv(
  pooled_summary,
  file.path(TAB, paste0(VERSION, "_04_pooled_receiver_Sham_Tx_summary.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. IL10-specific decomposition
# ------------------------------------------------------------------------------

il10_rw <- rw[rw$human_ligand == "IL10", , drop = FALSE]

il10_condition <- aggregate(
  cbind(
    state_fraction_of_primary5_Mphi,
    receptor_positive_fraction_within_state,
    receiver_weighted_fraction_of_primary5_Mphi
  ) ~ condition + macrophage_state,
  data = il10_rw,
  FUN = mean
)

sh <- il10_condition[il10_condition$condition == "Sham", , drop = FALSE]
tx <- il10_condition[il10_condition$condition == "Tx", , drop = FALSE]

il10_delta <- merge(
  sh,
  tx,
  by = "macrophage_state",
  suffixes = c("_Sham_mean", "_Tx_mean"),
  sort = FALSE
)

il10_delta$Tx_minus_Sham_state_fraction <-
  il10_delta$state_fraction_of_primary5_Mphi_Tx_mean -
  il10_delta$state_fraction_of_primary5_Mphi_Sham_mean

il10_delta$Tx_minus_Sham_within_state_IL10R_fraction <-
  il10_delta$receptor_positive_fraction_within_state_Tx_mean -
  il10_delta$receptor_positive_fraction_within_state_Sham_mean

il10_delta$Tx_minus_Sham_receiver_weighted_fraction <-
  il10_delta$receiver_weighted_fraction_of_primary5_Mphi_Tx_mean -
  il10_delta$receiver_weighted_fraction_of_primary5_Mphi_Sham_mean

il10_delta <- il10_delta[
  order(-il10_delta$Tx_minus_Sham_receiver_weighted_fraction),
  ,
  drop = FALSE
]

write.csv(
  il10_delta,
  file.path(TAB, paste0(VERSION, "_05_IL10_receiver_pool_state_decomposition.csv")),
  row.names = FALSE
)

# Contribution of each state to all IL10R-dual-positive primary-five Mphi
il10_contrib <- il10_rw
il10_contrib$total_IL10R_positive_in_sample <- ave(
  il10_contrib$receptor_positive_cells,
  il10_contrib$sample,
  FUN = sum
)
il10_contrib$state_fraction_of_IL10R_positive_pool <-
  il10_contrib$receptor_positive_cells /
  il10_contrib$total_IL10R_positive_in_sample

write.csv(
  il10_contrib,
  file.path(TAB, paste0(VERSION, "_06_IL10R_positive_pool_contribution_by_sample_state.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Existing direct-boundary and Mphi->HSC evidence
# ------------------------------------------------------------------------------

c22 <- read.delim(
  C22_BOUNDARY,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

hsc_ev <- read.csv(
  HSC_EVIDENCE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

hsc_burden <- read.csv(
  HSC_BURDEN,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.csv(
  c22,
  file.path(TAB, paste0(VERSION, "_07_direct_IL10R_Pdgfb_boundary_from_C22.csv")),
  row.names = FALSE
)

write.csv(
  hsc_ev,
  file.path(TAB, paste0(VERSION, "_08_existing_PDGF_FN1_evidence_v636.csv")),
  row.names = FALSE
)

write.csv(
  hsc_burden,
  file.path(TAB, paste0(VERSION, "_09_existing_ECM_HSC_burden_v636.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Compact network evidence chain
# ------------------------------------------------------------------------------

chain <- list()
ci <- 1L

il10_pool <- pooled_summary[
  pooled_summary$human_ligand == "IL10",
  ,
  drop = FALSE
]

chain[[ci]] <- data.frame(
  layer = "Graft_to_host_receiver",
  metric = "IL10R-dual-positive fraction across primary five Mphi states",
  Sham_mean = il10_pool$Sham_mean_pooled_fraction,
  Tx_mean = il10_pool$Tx_mean_pooled_fraction,
  Tx_minus_Sham = il10_pool$Tx_minus_Sham_pooled_fraction,
  evidence_relation = il10_pool$pooled_replicate_relation,
  interpretation =
    "Tx is associated with expansion of the receptor-transcript-competent primary Mphi pool.",
  boundary =
    "Transcript competence only; not IL10 binding, receptor occupancy, JAK activation, or pSTAT3.",
  stringsAsFactors = FALSE
)
ci <- ci + 1L

# Carry all PDGF evidence rows forward without changing the v6.3.6 metrics.
if (all(c("axis", "metric", "Sham_mean", "Tx_mean", "Tx_minus_Sham",
          "evidence_grade") %in% colnames(hsc_ev))) {
  pdgf <- hsc_ev[hsc_ev$axis == "PDGF", , drop = FALSE]
  for (j in seq_len(nrow(pdgf))) {
    chain[[ci]] <- data.frame(
      layer = paste0("Mphi_to_HSC_", pdgf$evidence_layer[j]),
      metric = pdgf$metric[j],
      Sham_mean = pdgf$Sham_mean[j],
      Tx_mean = pdgf$Tx_mean[j],
      Tx_minus_Sham = pdgf$Tx_minus_Sham[j],
      evidence_relation = pdgf$evidence_grade[j],
      interpretation =
        "Existing v6.3.6 PDGF-axis evidence; carried forward without recalculation.",
      boundary =
        "Selective pathway-level association; not proof that PDGF is the sole mediator.",
      stringsAsFactors = FALSE
    )
    ci <- ci + 1L
  }
}

if (all(c("metric", "Sham_mean", "Tx_mean", "Tx_minus_Sham",
          "evidence_grade") %in% colnames(hsc_burden))) {
  for (j in seq_len(nrow(hsc_burden))) {
    chain[[ci]] <- data.frame(
      layer = "HSC_burden",
      metric = hsc_burden$metric[j],
      Sham_mean = hsc_burden$Sham_mean[j],
      Tx_mean = hsc_burden$Tx_mean[j],
      Tx_minus_Sham = hsc_burden$Tx_minus_Sham[j],
      evidence_relation = hsc_burden$evidence_grade[j],
      interpretation =
        "Existing v6.3.6 HSC-state evidence; carried forward without recalculation.",
      boundary =
        "Association is coherent with network remodeling but does not identify a sole causal pathway.",
      stringsAsFactors = FALSE
    )
    ci <- ci + 1L
  }
}

chain_df <- do.call(rbind, chain)

write.csv(
  chain_df,
  file.path(TAB, paste0(VERSION, "_10_integrated_IL10_Mphi_HSC_network_evidence_chain.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Figures
# ------------------------------------------------------------------------------

# Figure 1: pooled receptor competence across ligands, biological samples shown
pdat <- pooled_sample
lig_order <- pooled_summary$human_ligand[
  order(-pooled_summary$Tx_minus_Sham_pooled_fraction)
]
pdat$human_ligand <- factor(pdat$human_ligand, levels = lig_order)

p1 <- ggplot(
  pdat,
  aes(
    x = human_ligand,
    y = pooled_receptor_positive_fraction,
    shape = condition
  )
) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 3
  ) +
  stat_summary(
    aes(group = condition),
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    position = position_dodge(width = 0.5)
  ) +
  labs(
    title = "Receptor-competent fraction of the primary five-state macrophage pool",
    subtitle = "Each point is one biological sample; descriptive Sham vs Tx comparison",
    x = "Human graft ligand",
    y = "Receptor-complex positive fraction",
    caption = "Transcript co-detection/availability only; not pathway activation."
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG, paste0(VERSION, "_01_pooled_receiver_competence_by_ligand.pdf")),
  p1,
  width = 11,
  height = 6
)
ggsave(
  file.path(FIG, paste0(VERSION, "_01_pooled_receiver_competence_by_ligand.png")),
  p1,
  width = 11,
  height = 6,
  dpi = 300
)

# Figure 2: IL10R-positive cell contribution per 1000 primary Mphi
p2dat <- il10_rw
p2dat$macrophage_state <- factor(
  p2dat$macrophage_state,
  levels = PRIMARY_STATES
)

p2 <- ggplot(
  p2dat,
  aes(
    x = sample,
    y = receiver_weighted_cells_per_1000_primary5_Mphi,
    fill = macrophage_state
  )
) +
  geom_col() +
  labs(
    title = "Composition of the IL10R-dual-positive primary macrophage pool",
    subtitle = "Receiver-weighted abundance = IL10R-dual-positive cells per 1000 primary five-state Mphi",
    x = NULL,
    y = "IL10R-dual-positive cells / 1000 Mphi",
    fill = "Mphi state",
    caption = "This integrates state abundance and within-state Il10ra+Il10rb transcript co-detection."
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG, paste0(VERSION, "_02_IL10R_positive_pool_state_composition.pdf")),
  p2,
  width = 9,
  height = 6
)
ggsave(
  file.path(FIG, paste0(VERSION, "_02_IL10R_positive_pool_state_composition.png")),
  p2,
  width = 9,
  height = 6,
  dpi = 300
)

# Figure 3: change decomposition for IL10 receiver-weighted pool
p3dat <- il10_delta
p3dat$macrophage_state <- factor(
  p3dat$macrophage_state,
  levels = rev(PRIMARY_STATES)
)

p3 <- ggplot(
  p3dat,
  aes(
    x = Tx_minus_Sham_receiver_weighted_fraction,
    y = macrophage_state
  )
) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  labs(
    title = "State-level contribution to the Tx-associated change in the IL10R-competent Mphi pool",
    subtitle = "Positive values indicate a larger IL10R-dual-positive contribution after transplantation",
    x = "Tx - Sham receiver-weighted fraction of primary Mphi",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG, paste0(VERSION, "_03_IL10_receiver_pool_change_decomposition.pdf")),
  p3,
  width = 9,
  height = 5.5
)
ggsave(
  file.path(FIG, paste0(VERSION, "_03_IL10_receiver_pool_change_decomposition.png")),
  p3,
  width = 9,
  height = 5.5,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 9. Boundaries and provenance
# ------------------------------------------------------------------------------

boundaries <- data.frame(
  item = c(
    "primary_macrophage_universe",
    "receiver_weighted_definition",
    "statistics_boundary",
    "IL10_activation_boundary",
    "direct_IL10R_Pdgfb_boundary",
    "PDGF_boundary",
    "human_graft_boundary",
    "nonIL10_direction_boundary"
  ),
  value = c(
    "Five final v4.14.5 macrophage states; Other excluded, matching v1.4.3C.2.1 primary denominator",
    "state fraction among primary five Mphi multiplied by within-state receptor-positive fraction; algebraically receptor-positive cells divided by all primary five Mphi in that sample",
    "Descriptive only; n=2 Sham and n=2 Tx",
    "Il10ra+Il10rb transcript co-detection supports response competence only; it is not pSTAT3 evidence",
    "v1.4.3C.2.2 remains negative: no Tx-specific replicate-consistent direct cell-autonomous IL10R-status to Pdgfb suppression",
    "v6.3.6 supports selective attenuation of the Mphi-to-HSC PDGF mitogenic arm, not global HSC suppression and not PDGF as the sole mediator",
    "SfLB to Tx is an engraftment-associated comparison, not an IL10-only perturbation",
    "Non-IL10 graft ligands are candidate signals only; receptor availability does not establish whether their net biological action is beneficial, neutral, or adverse"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  boundaries,
  file.path(TAB, paste0(VERSION, "_11_interpretation_boundaries.csv")),
  row.names = FALSE
)

prov <- data.frame(
  key = c(
    "version",
    "B1_input",
    "C21_abundance",
    "C22_direct_boundary",
    "HSC_evidence",
    "HSC_burden"
  ),
  value = c(
    VERSION,
    B1,
    C21_ABUND,
    C22_BOUNDARY,
    HSC_EVIDENCE,
    HSC_BURDEN
  ),
  stringsAsFactors = FALSE
)

write.csv(
  prov,
  file.path(TAB, paste0(VERSION, "_12_provenance.csv")),
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
cat("POOLED PRIMARY-5 Mphi RECEIVER COMPETENCE\n")
cat("============================================================================\n")

show <- pooled_summary[
  order(-pooled_summary$Tx_minus_Sham_pooled_fraction),
  c(
    "human_ligand",
    "Sham_mean_pooled_fraction",
    "Tx_mean_pooled_fraction",
    "Tx_minus_Sham_pooled_fraction",
    "pooled_replicate_relation"
  ),
  drop = FALSE
]
print(show, row.names = FALSE)

cat("\n============================================================================\n")
cat("IL10 RECEIVER-POOL STATE DECOMPOSITION\n")
cat("============================================================================\n")
print(
  il10_delta[, c(
    "macrophage_state",
    "state_fraction_of_primary5_Mphi_Sham_mean",
    "state_fraction_of_primary5_Mphi_Tx_mean",
    "Tx_minus_Sham_state_fraction",
    "receptor_positive_fraction_within_state_Sham_mean",
    "receptor_positive_fraction_within_state_Tx_mean",
    "Tx_minus_Sham_receiver_weighted_fraction"
  )],
  row.names = FALSE
)

cat("\n============================================================================\n")
cat(VERSION, " COMPLETE\n")
cat("============================================================================\n")
cat("Output:", OUT, "\n")
cat("B.1 vs C.2.1 denominator validation: PASS\n")
cat("Receiver-weighted primary Mphi integration: COMPLETE\n")
cat("IL10 Mphi-pool decomposition: COMPLETE\n")
cat("Existing Mphi->HSC / HSC-burden evidence carried forward: COMPLETE\n")
cat("No inferential test performed (n=2/group).\n")
cat("No pSTAT3 or direct IL10R->Pdgfb causal claim made.\n")
cat("External delivery HDDs required: NO\n")
