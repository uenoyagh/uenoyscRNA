#!/usr/bin/env Rscript

# ==============================================================================
# Human_MASH_v1.4.4B.1_GraftLigand_MouseMphi5_ReceiverBridge.R
#
# Purpose
# -------
# Integrate the latest human graft-ligand transition/source results with the
# final mouse macrophage (Mphi) state object and quantify receptor-complex
# transcript competence in Sham1/Sham20 vs Tx17/Tx5.
#
# Primary analysis:
#   - 5 final mouse Mphi states
#   - descriptive n=2 Sham vs n=2 Tx
#   - no inferential test
#
# Context:
#   - "Other" Mphi-like cells are retained in context tables but excluded from
#     primary 5-state ranking/figures.
#
# Receptor-complex logic is frozen from v1.4.3C.1.1:
#   IL10    -> Il10ra AND Il10rb
#   FGF1    -> ANY(Fgfr1,Fgfr2,Fgfr3,Fgfr4)
#   TNFSF12 -> Tnfrsf12a
#   VEGFA   -> Flt1 OR Kdr
#   CXCL12  -> Cxcr4 OR Ackr3
#   TGFB3   -> Tgfbr1 AND Tgfbr2
#   BMP5/2  -> Bmpr2 AND ANY(Bmpr1a,Acvr1)
#   HBEGF/BTC -> Egfr OR Erbb4
#
# Important boundaries
# --------------------
# - Receptor transcript co-detection = response competence only.
# - It is NOT ligand binding, receptor occupancy, pathway activation, or pSTAT3.
# - SfLB -> Tx is not an IL10-only perturbation.
# - Bmpr1b is outside the relevant measured feature space; BMPR logic follows
#   the validated v1.4.3C.1.1 definition.
# - This script does not claim direct IL10R -> Pdgfb repression.
# ==============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

VERSION <- "v1.4.4B.1"

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

MOUSE_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS",
  "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

HUMAN_A2 <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.4A.2",
  "analysis",
  "v1.4.4A.2_probe_aware_ligand_transition_audit",
  "tables",
  "03_ligand_transition_classification_v1.4.4A.2.csv"
)

HUMAN_A31 <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.4A.3.1",
  "analysis",
  "v1.4.4A.3.1_ligand_source_lineage_atlas_full_denominator_fix",
  "tables",
  "04_nonIL10_dominant_source_and_source_shift_summary_v1.4.4A.3.1.csv"
)

OUT <- file.path(
  ROOT,
  "Human_MASH_RDS",
  "Human_MASH_12sample_Context_v1.4.4B.1",
  "analysis",
  "v1.4.4B.1_graft_ligand_mouse_Mphi5_receiver_bridge"
)

TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
LOG <- file.path(OUT, "logs")

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOG, paste0(VERSION, "_run_log.txt"))
sink(LOGFILE, split = TRUE)
on.exit(sink(), add = TRUE)

cat("============================================================================\n")
cat("Human MASH ", VERSION, "\n", sep = "")
cat("Human graft ligand -> mouse Mphi 5-state receptor bridge\n")
cat("============================================================================\n\n")

required_files <- c(MOUSE_RDS, HUMAN_A2, HUMAN_A31)
for (f in required_files) {
  if (!file.exists(f)) stop("Required input missing: ", f)
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

read_csv_safe <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(NA_real_)
  mean(x[is.finite(x)])
}

safe_min <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(NA_real_)
  min(x[is.finite(x)])
}

safe_max <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(NA_real_)
  max(x[is.finite(x)])
}

replicate_relation <- function(sham1, sham20, tx17, tx5) {
  z <- c(sham1, sham20, tx17, tx5)
  if (any(!is.finite(z))) return("Incomplete_four_sample_data")
  sham <- c(sham1, sham20)
  tx <- c(tx17, tx5)
  if (min(tx) > max(sham)) {
    return("Tx_higher_complete_replicate_range_separation")
  }
  if (max(tx) < min(sham)) {
    return("Tx_lower_complete_replicate_range_separation")
  }
  "Replicate_ranges_overlap"
}

get_counts <- function(obj) {
  out <- tryCatch(
    LayerData(obj, assay = "RNA", layer = "counts"),
    error = function(e1) {
      tryCatch(
        GetAssayData(obj, assay = "RNA", layer = "counts"),
        error = function(e2) {
          GetAssayData(obj, assay = "RNA", slot = "counts")
        }
      )
    }
  )
  out
}

# ------------------------------------------------------------------------------
# 1. Frozen ligand -> receptor-complex map
# ------------------------------------------------------------------------------

ligand_map <- data.frame(
  human_ligand = c(
    "IL10",
    "FGF1",
    "TNFSF12",
    "VEGFA",
    "CXCL12",
    "TGFB3",
    "BMP5",
    "BMP2",
    "HBEGF",
    "BTC"
  ),
  receiver_complex = c(
    "IL10R_Il10ra_AND_Il10rb",
    "FGFR_ANY_Fgfr1_4",
    "TWEAKR_Tnfrsf12a",
    "VEGFR_core_Flt1_OR_Kdr",
    "CXCL12R_Cxcr4_OR_Ackr3",
    "TGFBR_Tgfbr1_AND_Tgfbr2",
    "BMPR_Bmpr2_AND_typeI",
    "BMPR_Bmpr2_AND_typeI",
    "EGFR_Egfr_OR_Erbb4",
    "EGFR_Egfr_OR_Erbb4"
  ),
  receptor_logic = c(
    "Il10ra AND Il10rb",
    "ANY(Fgfr1,Fgfr2,Fgfr3,Fgfr4)",
    "Tnfrsf12a",
    "Flt1 OR Kdr",
    "Cxcr4 OR Ackr3",
    "Tgfbr1 AND Tgfbr2",
    "Bmpr2 AND ANY(Bmpr1a,Acvr1)",
    "Bmpr2 AND ANY(Bmpr1a,Acvr1)",
    "Egfr OR Erbb4",
    "Egfr OR Erbb4"
  ),
  stringsAsFactors = FALSE
)

complex_genes <- list(
  IL10R_Il10ra_AND_Il10rb = c("Il10ra", "Il10rb"),
  FGFR_ANY_Fgfr1_4 = c("Fgfr1", "Fgfr2", "Fgfr3", "Fgfr4"),
  TWEAKR_Tnfrsf12a = c("Tnfrsf12a"),
  VEGFR_core_Flt1_OR_Kdr = c("Flt1", "Kdr"),
  CXCL12R_Cxcr4_OR_Ackr3 = c("Cxcr4", "Ackr3"),
  TGFBR_Tgfbr1_AND_Tgfbr2 = c("Tgfbr1", "Tgfbr2"),
  BMPR_Bmpr2_AND_typeI = c("Bmpr2", "Bmpr1a", "Acvr1"),
  EGFR_Egfr_OR_Erbb4 = c("Egfr", "Erbb4")
)

complex_positive <- function(complex_name, count_mat) {
  if (complex_name == "IL10R_Il10ra_AND_Il10rb") {
    return((count_mat["Il10ra", ] > 0) & (count_mat["Il10rb", ] > 0))
  }
  if (complex_name == "FGFR_ANY_Fgfr1_4") {
    return(Matrix::colSums(count_mat[c("Fgfr1","Fgfr2","Fgfr3","Fgfr4"), , drop = FALSE] > 0) > 0)
  }
  if (complex_name == "TWEAKR_Tnfrsf12a") {
    return(count_mat["Tnfrsf12a", ] > 0)
  }
  if (complex_name == "VEGFR_core_Flt1_OR_Kdr") {
    return((count_mat["Flt1", ] > 0) | (count_mat["Kdr", ] > 0))
  }
  if (complex_name == "CXCL12R_Cxcr4_OR_Ackr3") {
    return((count_mat["Cxcr4", ] > 0) | (count_mat["Ackr3", ] > 0))
  }
  if (complex_name == "TGFBR_Tgfbr1_AND_Tgfbr2") {
    return((count_mat["Tgfbr1", ] > 0) & (count_mat["Tgfbr2", ] > 0))
  }
  if (complex_name == "BMPR_Bmpr2_AND_typeI") {
    return(
      (count_mat["Bmpr2", ] > 0) &
      ((count_mat["Bmpr1a", ] > 0) | (count_mat["Acvr1", ] > 0))
    )
  }
  if (complex_name == "EGFR_Egfr_OR_Erbb4") {
    return((count_mat["Egfr", ] > 0) | (count_mat["Erbb4", ] > 0))
  }
  stop("Unknown receiver complex: ", complex_name)
}

write.csv(
  ligand_map,
  file.path(TAB, paste0(VERSION, "_01_frozen_ligand_receptor_map.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 2. Read latest human graft ligand evidence
# ------------------------------------------------------------------------------

a2 <- read_csv_safe(HUMAN_A2)
a31 <- read_csv_safe(HUMAN_A31)

needed_a2 <- c(
  "gene",
  "SfLB1_shared_TRUE_UMI_per_1000",
  "SfLB2_shared_TRUE_UMI_per_1000",
  "Tx17_shared_TRUE_UMI_per_1000",
  "Tx5_shared_TRUE_UMI_per_1000",
  "Tx_vs_SfLB_log2FC_per_cell_descriptive",
  "preexisting_in_both_SfLB_strict",
  "detected_in_both_Tx_strict",
  "replicate_relation",
  "transition_call",
  "IL10_special_boundary"
)

miss_a2 <- setdiff(needed_a2, colnames(a2))
if (length(miss_a2) > 0) {
  stop("A.2 missing expected columns: ", paste(miss_a2, collapse = ", "))
}

needed_a31 <- c(
  "gene",
  "SfLB_dominant_source",
  "SfLB_dominant_source_fraction",
  "Tx_dominant_source",
  "Tx_dominant_source_fraction",
  "dominant_source_switch",
  "Tx_mean_source_fraction_from_Tx_emergent_lineages",
  "SfLB_mean_source_fraction_from_lineages_lost_after_engraftment"
)

miss_a31 <- setdiff(needed_a31, colnames(a31))
if (length(miss_a31) > 0) {
  stop("A.3.1 missing expected columns: ", paste(miss_a31, collapse = ", "))
}

a2$gene <- toupper(trimws(a2$gene))
a31$gene <- toupper(trimws(a31$gene))

human <- merge(
  ligand_map,
  a2,
  by.x = "human_ligand",
  by.y = "gene",
  all.x = TRUE,
  sort = FALSE
)

human <- merge(
  human,
  a31,
  by.x = "human_ligand",
  by.y = "gene",
  all.x = TRUE,
  sort = FALSE
)

human <- human[
  match(ligand_map$human_ligand, human$human_ligand),
  ,
  drop = FALSE
]

human$source_note <- ifelse(
  human$human_ligand == "IL10",
  "Engineered anchor; IL10 source-lineage detail handled in v1.4.4A.3.1",
  paste0(
    "SfLB dominant=", human$SfLB_dominant_source,
    " (", sprintf("%.3f", human$SfLB_dominant_source_fraction), "); ",
    "Tx dominant=", human$Tx_dominant_source,
    " (", sprintf("%.3f", human$Tx_dominant_source_fraction), ")"
  )
)

write.csv(
  human,
  file.path(TAB, paste0(VERSION, "_02_human_graft_ligand_evidence.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 3. Read final mouse Mphi object and validate metadata / feature space
# ------------------------------------------------------------------------------

obj <- readRDS(MOUSE_RDS)

if (!inherits(obj, "Seurat")) stop("Mouse input is not a Seurat object.")
if (ncol(obj) != 18029L) {
  stop("Unexpected final mouse Mphi cell count: ", ncol(obj), " (expected 18029)")
}
if (!"RNA" %in% Assays(obj)) stop("RNA assay missing.")

required_meta <- c(
  "sample_4group",
  "condition_4group",
  "macrophage_class_Res2_FINAL_v4145"
)

miss_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(miss_meta) > 0) {
  stop("Missing final Mphi metadata columns: ", paste(miss_meta, collapse = ", "))
}

counts <- get_counts(obj)

all_receptor_genes <- unique(unlist(complex_genes, use.names = FALSE))
missing_receptor_genes <- setdiff(all_receptor_genes, rownames(counts))

feature_audit <- data.frame(
  gene = all_receptor_genes,
  present_in_final_mouse_Mphi_RNA = all_receptor_genes %in% rownames(counts),
  stringsAsFactors = FALSE
)

write.csv(
  feature_audit,
  file.path(TAB, paste0(VERSION, "_03_mouse_receptor_feature_audit.csv")),
  row.names = FALSE
)

if (length(missing_receptor_genes) > 0) {
  stop(
    "Required receptor genes missing from final mouse Mphi RNA assay: ",
    paste(missing_receptor_genes, collapse = ", ")
  )
}

cat("Mouse final Mphi RDS: PASS\n")
cat("ncell =", ncol(obj), "\n")
cat("RNA features =", nrow(counts), "\n")
cat("Required receptor genes present =", length(all_receptor_genes), "/", length(all_receptor_genes), "\n\n")

# ------------------------------------------------------------------------------
# 4. Define primary 5-state and context state universe
# ------------------------------------------------------------------------------

PRIMARY_STATES <- c(
  "Inflammatory-Mphi",
  "Anti-inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi"
)

CONTEXT_STATES <- c(PRIMARY_STATES, "Other")
PRIMARY_SAMPLES <- c("Sham1", "Sham20", "Tx17", "Tx5")

md <- obj@meta.data
md$sample_bridge <- as.character(md$sample_4group)
md$condition_bridge <- as.character(md$condition_4group)
md$state_bridge <- as.character(md$macrophage_class_Res2_FINAL_v4145)

state_levels_found <- sort(unique(md$state_bridge))
unexpected_states <- setdiff(state_levels_found, CONTEXT_STATES)
missing_primary_states <- setdiff(PRIMARY_STATES, state_levels_found)

if (length(unexpected_states) > 0) {
  stop("Unexpected final Mphi state(s): ", paste(unexpected_states, collapse = ", "))
}
if (length(missing_primary_states) > 0) {
  stop("Missing primary Mphi state(s): ", paste(missing_primary_states, collapse = ", "))
}

sample_state_counts <- as.data.frame(
  table(
    sample = md$sample_bridge,
    macrophage_state = md$state_bridge,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)
names(sample_state_counts)[3] <- "n_cells"

write.csv(
  sample_state_counts,
  file.path(TAB, paste0(VERSION, "_04_final_Mphi_sample_state_counts.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Sample x state x receptor-complex competence
# ------------------------------------------------------------------------------

complex_rows <- list()
gene_rows <- list()
ci <- 1L
gi <- 1L

for (ss in PRIMARY_SAMPLES) {
  for (st in CONTEXT_STATES) {

    idx <- which(md$sample_bridge == ss & md$state_bridge == st)
    n <- length(idx)

    if (n == 0L) next

    cond <- unique(md$condition_bridge[idx])
    if (length(cond) != 1L) {
      stop("Condition is not unique for sample/state: ", ss, " / ", st)
    }

    needed <- counts[all_receptor_genes, idx, drop = FALSE]

    # per-gene support
    for (gg in all_receptor_genes) {
      x <- needed[gg, , drop = TRUE]

      gene_rows[[gi]] <- data.frame(
        sample = ss,
        condition = cond,
        macrophage_state = st,
        primary_5state = st %in% PRIMARY_STATES,
        gene = gg,
        n_cells = n,
        positive_cells = sum(x > 0),
        positive_fraction = sum(x > 0) / n,
        total_UMI = sum(x),
        UMI_per_1000_cells = sum(x) / n * 1000,
        stringsAsFactors = FALSE
      )
      gi <- gi + 1L
    }

    # receptor complex support
    for (cc in names(complex_genes)) {
      cp <- complex_positive(cc, needed)

      complex_rows[[ci]] <- data.frame(
        sample = ss,
        condition = cond,
        macrophage_state = st,
        primary_5state = st %in% PRIMARY_STATES,
        receiver_complex = cc,
        receptor_logic = ligand_map$receptor_logic[
          match(cc, ligand_map$receiver_complex)
        ],
        n_cells = n,
        complex_positive_cells = sum(cp),
        complex_positive_fraction = sum(cp) / n,
        stringsAsFactors = FALSE
      )
      ci <- ci + 1L
    }
  }
}

gene_support <- do.call(rbind, gene_rows)
complex_support <- do.call(rbind, complex_rows)

# receptor_logic can map to duplicate ligand rows but same logic; repair if needed
logic_lookup <- unique(ligand_map[, c("receiver_complex", "receptor_logic")])
complex_support$receptor_logic <- logic_lookup$receptor_logic[
  match(complex_support$receiver_complex, logic_lookup$receiver_complex)
]

write.csv(
  gene_support,
  file.path(TAB, paste0(VERSION, "_05_receptor_gene_support_by_sample_state.csv")),
  row.names = FALSE
)

write.csv(
  complex_support,
  file.path(TAB, paste0(VERSION, "_06_receptor_complex_by_sample_state.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Ligand x primary Mphi-state integrated Sham-vs-Tx summary
# ------------------------------------------------------------------------------

integrated_rows <- list()
ii <- 1L

for (lig in ligand_map$human_ligand) {

  cc <- ligand_map$receiver_complex[match(lig, ligand_map$human_ligand)]
  hm <- human[human$human_ligand == lig, , drop = FALSE]

  for (st in PRIMARY_STATES) {

    z <- complex_support[
      complex_support$receiver_complex == cc &
        complex_support$macrophage_state == st,
      ,
      drop = FALSE
    ]

    get_pf <- function(sample_name) {
      v <- z$complex_positive_fraction[z$sample == sample_name]
      if (length(v) != 1L) return(NA_real_)
      as.numeric(v)
    }

    get_n <- function(sample_name) {
      v <- z$n_cells[z$sample == sample_name]
      if (length(v) != 1L) return(NA_real_)
      as.numeric(v)
    }

    get_posn <- function(sample_name) {
      v <- z$complex_positive_cells[z$sample == sample_name]
      if (length(v) != 1L) return(NA_real_)
      as.numeric(v)
    }

    sham1 <- get_pf("Sham1")
    sham20 <- get_pf("Sham20")
    tx17 <- get_pf("Tx17")
    tx5 <- get_pf("Tx5")

    sham_mean <- safe_mean(c(sham1, sham20))
    tx_mean <- safe_mean(c(tx17, tx5))

    integrated_rows[[ii]] <- data.frame(
      human_ligand = lig,
      receiver_complex = cc,
      receptor_logic = ligand_map$receptor_logic[
        match(lig, ligand_map$human_ligand)
      ],

      human_transition_call = hm$transition_call,
      human_ligand_replicate_relation = hm$replicate_relation,
      human_Tx_vs_SfLB_log2FC_per_cell_descriptive =
        hm$Tx_vs_SfLB_log2FC_per_cell_descriptive,
      human_preexisting_in_both_SfLB_strict =
        hm$preexisting_in_both_SfLB_strict,
      human_detected_in_both_Tx_strict =
        hm$detected_in_both_Tx_strict,

      SfLB_dominant_source = hm$SfLB_dominant_source,
      SfLB_dominant_source_fraction = hm$SfLB_dominant_source_fraction,
      Tx_dominant_source = hm$Tx_dominant_source,
      Tx_dominant_source_fraction = hm$Tx_dominant_source_fraction,
      dominant_source_switch = hm$dominant_source_switch,
      source_note = hm$source_note,

      macrophage_state = st,

      Sham1_n_cells = get_n("Sham1"),
      Sham20_n_cells = get_n("Sham20"),
      Tx17_n_cells = get_n("Tx17"),
      Tx5_n_cells = get_n("Tx5"),

      Sham1_receptor_positive_cells = get_posn("Sham1"),
      Sham20_receptor_positive_cells = get_posn("Sham20"),
      Tx17_receptor_positive_cells = get_posn("Tx17"),
      Tx5_receptor_positive_cells = get_posn("Tx5"),

      Sham1_receptor_positive_fraction = sham1,
      Sham20_receptor_positive_fraction = sham20,
      Tx17_receptor_positive_fraction = tx17,
      Tx5_receptor_positive_fraction = tx5,

      Sham_mean_receptor_positive_fraction = sham_mean,
      Tx_mean_receptor_positive_fraction = tx_mean,
      Tx_min_receptor_positive_fraction = safe_min(c(tx17, tx5)),
      Tx_max_receptor_positive_fraction = safe_max(c(tx17, tx5)),
      Tx_minus_Sham_receptor_positive_fraction = tx_mean - sham_mean,
      receptor_replicate_relation =
        replicate_relation(sham1, sham20, tx17, tx5),

      inference =
        "descriptive_receptor_transcript_competence_n2_per_condition_not_pathway_activation",

      stringsAsFactors = FALSE
    )

    ii <- ii + 1L
  }
}

integrated <- do.call(rbind, integrated_rows)

# rank primary Mphi states within each ligand by Tx receptor competence
integrated$Tx_receiver_rank_within_ligand <- ave(
  -integrated$Tx_mean_receptor_positive_fraction,
  integrated$human_ligand,
  FUN = function(x) rank(x, ties.method = "min", na.last = "keep")
)

integrated <- integrated[
  order(
    match(integrated$human_ligand, ligand_map$human_ligand),
    integrated$Tx_receiver_rank_within_ligand,
    integrated$macrophage_state
  ),
  ,
  drop = FALSE
]

write.csv(
  integrated,
  file.path(TAB, paste0(VERSION, "_07_integrated_graft_ligand_mouse_Mphi5_bridge.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Top receiver states per ligand
# ------------------------------------------------------------------------------

top_rows <- list()
ti <- 1L

for (lig in ligand_map$human_ligand) {
  z <- integrated[integrated$human_ligand == lig, , drop = FALSE]
  z <- z[
    order(
      -z$Tx_mean_receptor_positive_fraction,
      -z$Tx_min_receptor_positive_fraction,
      z$macrophage_state
    ),
    ,
    drop = FALSE
  ]

  top_rows[[ti]] <- z[
    seq_len(min(3L, nrow(z))),
    c(
      "human_ligand",
      "human_transition_call",
      "Tx_dominant_source",
      "receiver_complex",
      "receptor_logic",
      "macrophage_state",
      "Tx_receiver_rank_within_ligand",
      "Sham_mean_receptor_positive_fraction",
      "Tx17_receptor_positive_fraction",
      "Tx5_receptor_positive_fraction",
      "Tx_mean_receptor_positive_fraction",
      "Tx_min_receptor_positive_fraction",
      "Tx_minus_Sham_receptor_positive_fraction",
      "receptor_replicate_relation",
      "human_Tx_vs_SfLB_log2FC_per_cell_descriptive"
    ),
    drop = FALSE
  ]

  ti <- ti + 1L
}

top3 <- do.call(rbind, top_rows)

write.csv(
  top3,
  file.path(TAB, paste0(VERSION, "_08_top3_mouse_Mphi_receivers_per_human_ligand.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. IL10-focused state table
# ------------------------------------------------------------------------------

il10 <- integrated[
  integrated$human_ligand == "IL10",
  ,
  drop = FALSE
]

il10 <- il10[
  order(-il10$Tx_mean_receptor_positive_fraction),
  ,
  drop = FALSE
]

write.csv(
  il10,
  file.path(TAB, paste0(VERSION, "_09_IL10_mouse_Mphi5_receiver_state_summary.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9. Human-ligand priority subset
# ------------------------------------------------------------------------------

priority <- integrated[
  integrated$human_ligand == "IL10" |
    grepl(
      "Tx_higher_complete_replicate_range_separation",
      integrated$human_transition_call,
      fixed = TRUE
    ),
  ,
  drop = FALSE
]

write.csv(
  priority,
  file.path(TAB, paste0(VERSION, "_10_priority_graft_ligand_mouse_Mphi5_bridge.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 10. Figures
# ------------------------------------------------------------------------------

state_order <- c(
  "Inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Lipid-associated/TREM2-Mphi",
  "Anti-inflammatory-Mphi",
  "Repair/Resolution-Mphi"
)

plot_df <- integrated
plot_df$macrophage_state <- factor(
  plot_df$macrophage_state,
  levels = rev(state_order)
)
plot_df$human_ligand <- factor(
  plot_df$human_ligand,
  levels = ligand_map$human_ligand
)

p1 <- ggplot(
  plot_df,
  aes(
    x = human_ligand,
    y = macrophage_state,
    fill = Tx_mean_receptor_positive_fraction
  )
) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(
    aes(label = sprintf("%.2f", Tx_mean_receptor_positive_fraction)),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#FFFFFF",
    mid = "#7FC8F8",
    high = "#0033FF",
    midpoint = 0.25,
    limits = c(0, max(plot_df$Tx_mean_receptor_positive_fraction, na.rm = TRUE)),
    name = "Tx mean\npositive fraction"
  ) +
  labs(
    title = "Human graft ligand -> mouse macrophage receptor competence",
    subtitle = "Final 5 macrophage states; Tx17/Tx5 mean receptor-complex positive fraction",
    x = "Human graft ligand",
    y = NULL,
    caption = "Receptor transcript competence only; not ligand binding or pathway activation."
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG, paste0(VERSION, "_01_Tx_receiver_competence_heatmap.pdf")),
  p1,
  width = 11.5,
  height = 5.8,
  units = "in"
)

ggsave(
  file.path(FIG, paste0(VERSION, "_01_Tx_receiver_competence_heatmap.png")),
  p1,
  width = 11.5,
  height = 5.8,
  units = "in",
  dpi = 300
)

max_abs_delta <- max(
  abs(plot_df$Tx_minus_Sham_receptor_positive_fraction),
  na.rm = TRUE
)
if (!is.finite(max_abs_delta) || max_abs_delta == 0) max_abs_delta <- 1

p2 <- ggplot(
  plot_df,
  aes(
    x = human_ligand,
    y = macrophage_state,
    fill = Tx_minus_Sham_receptor_positive_fraction
  )
) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(
    aes(label = sprintf("%+.2f", Tx_minus_Sham_receptor_positive_fraction)),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(-max_abs_delta, max_abs_delta),
    name = "Tx - Sham\npositive fraction"
  ) +
  labs(
    title = "Treatment-associated change in mouse macrophage receptor competence",
    subtitle = "Descriptive Sham1/Sham20 vs Tx17/Tx5 comparison; n=2 biological samples per condition",
    x = "Human graft ligand",
    y = NULL,
    caption = "Positive values indicate higher receptor-complex transcript co-detection in Tx."
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG, paste0(VERSION, "_02_Tx_minus_Sham_receiver_change_heatmap.pdf")),
  p2,
  width = 11.5,
  height = 5.8,
  units = "in"
)

ggsave(
  file.path(FIG, paste0(VERSION, "_02_Tx_minus_Sham_receiver_change_heatmap.png")),
  p2,
  width = 11.5,
  height = 5.8,
  units = "in",
  dpi = 300
)

# Integrated bridge bubble plot:
# x = human ligand change, point size = Tx receptor competence,
# fill = Tx-Sham receptor change.
bubble <- integrated
bubble$macrophage_state <- factor(
  bubble$macrophage_state,
  levels = rev(state_order)
)
bubble$human_ligand <- factor(
  bubble$human_ligand,
  levels = ligand_map$human_ligand
)

p3 <- ggplot(
  bubble,
  aes(
    x = human_ligand,
    y = macrophage_state,
    size = Tx_mean_receptor_positive_fraction,
    fill = Tx_minus_Sham_receptor_positive_fraction
  )
) +
  geom_point(shape = 21, color = "black", stroke = 0.35) +
  scale_size_continuous(
    range = c(2.5, 10),
    name = "Tx mean receptor\npositive fraction"
  ) +
  scale_fill_gradient2(
    low = "#0033FF",
    mid = "#FFFFFF",
    high = "#FF1A1A",
    midpoint = 0,
    limits = c(-max_abs_delta, max_abs_delta),
    name = "Tx - Sham\nreceptor fraction"
  ) +
  labs(
    title = "Integrated graft-ligand -> mouse macrophage receiver landscape",
    subtitle = "Human graft ligand evidence is integrated in the companion table; bubble size shows Tx receiver competence",
    x = "Human graft ligand",
    y = NULL,
    caption = "Descriptive receptor-transcript landscape; no causal receptor activation inference."
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(FIG, paste0(VERSION, "_03_integrated_graft_ligand_Mphi_receiver_bubble.pdf")),
  p3,
  width = 11.5,
  height = 5.8,
  units = "in"
)

ggsave(
  file.path(FIG, paste0(VERSION, "_03_integrated_graft_ligand_Mphi_receiver_bubble.png")),
  p3,
  width = 11.5,
  height = 5.8,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# 11. Interpretation boundaries / provenance
# ------------------------------------------------------------------------------

boundaries <- data.frame(
  item = c(
    "analysis_unit",
    "statistics_boundary",
    "receptor_boundary",
    "IL10_activation_boundary",
    "direct_Pdgfb_boundary",
    "SfLB_Tx_boundary",
    "BMPR_logic",
    "primary_state_boundary",
    "unmapped_latest_graft_ligands"
  ),
  value = c(
    "Mouse final macrophage states from v4.14.5; primary Sham1/Sham20 versus Tx17/Tx5",
    "Descriptive only; n=2 biological samples per condition; no inferential test",
    "Receptor transcript co-detection supports competence only and is not ligand binding, receptor occupancy, or signaling activation",
    "Il10ra/Il10rb transcript competence does not establish JAK activation or pSTAT3",
    "This analysis does not claim direct cell-autonomous IL10R activation suppresses Pdgfb; v1.4.3C.2.2 remains the negative boundary",
    "SfLB to Tx combines transgene status, transplantation/engraftment, host microenvironment, selection/state transition, and possible residual xenogeneic innate pressure",
    "BMPR = Bmpr2 AND ANY(Bmpr1a,Acvr1); Bmpr1b is not used, matching validated v1.4.3C.1.1 logic",
    "Primary figures/ranking use five final Mphi states; Other is retained only in context tables",
    "FST, IGFBP2, IGFBP4, INHBA, VEGFB are not assigned receptor complexes in frozen v1.4.3C.1.1 and are intentionally not added here without a separate receptor-definition audit"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  boundaries,
  file.path(TAB, paste0(VERSION, "_11_interpretation_boundaries.csv")),
  row.names = FALSE
)

provenance <- data.frame(
  key = c(
    "version",
    "mouse_Mphi_RDS",
    "human_ligand_transition",
    "human_source_shift",
    "mouse_annotation_column",
    "sample_column",
    "condition_column",
    "receptor_definition_source"
  ),
  value = c(
    VERSION,
    MOUSE_RDS,
    HUMAN_A2,
    HUMAN_A31,
    "macrophage_class_Res2_FINAL_v4145",
    "sample_4group",
    "condition_4group",
    "v1.4.3C.1.1 frozen receptor-complex logic"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  provenance,
  file.path(TAB, paste0(VERSION, "_12_provenance.csv")),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(LOG, paste0(VERSION, "_sessionInfo.txt"))
)

# ------------------------------------------------------------------------------
# 12. Console summary
# ------------------------------------------------------------------------------

cat("============================================================================\n")
cat("PRIMARY TOP RECEIVER STATE PER HUMAN LIGAND\n")
cat("============================================================================\n")

for (lig in ligand_map$human_ligand) {
  z <- integrated[integrated$human_ligand == lig, , drop = FALSE]
  z <- z[order(-z$Tx_mean_receptor_positive_fraction), , drop = FALSE]
  top <- z[1, ]

  cat(
    sprintf(
      "%-8s -> %-34s | Tx mean=%.4f | Sham mean=%.4f | delta=%+.4f | %s\n",
      lig,
      top$macrophage_state,
      top$Tx_mean_receptor_positive_fraction,
      top$Sham_mean_receptor_positive_fraction,
      top$Tx_minus_Sham_receptor_positive_fraction,
      top$receptor_replicate_relation
    )
  )
}

cat("\n============================================================================\n")
cat("IL10 RECEIVER LANDSCAPE ACROSS FINAL 5 Mphi STATES\n")
cat("============================================================================\n")

print(
  il10[, c(
    "macrophage_state",
    "Sham1_receptor_positive_fraction",
    "Sham20_receptor_positive_fraction",
    "Tx17_receptor_positive_fraction",
    "Tx5_receptor_positive_fraction",
    "Sham_mean_receptor_positive_fraction",
    "Tx_mean_receptor_positive_fraction",
    "Tx_minus_Sham_receptor_positive_fraction",
    "receptor_replicate_relation"
  )],
  row.names = FALSE
)

cat("\n============================================================================\n")
cat(VERSION, " COMPLETE\n")
cat("============================================================================\n")
cat("Output:", OUT, "\n")
cat("Primary 5-state bridge: COMPLETE\n")
cat("Human ligand transition/source integration: COMPLETE\n")
cat("Receptor-complex transcript competence: COMPLETE\n")
cat("No inferential test performed (n=2/group).\n")
cat("No pSTAT3 or direct IL10R->Pdgfb causal claim made.\n")
cat("External delivery HDDs required: NO\n")
