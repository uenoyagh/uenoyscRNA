#!/usr/bin/env Rscript

# ============================================================
# Mouse MASH Mphi unbiased pseudobulk Tx-response ranking v6.4.9.1
#
# PURPOSE
#   After curated IL10-response adjudication, perform an unbiased
#   sample-level pseudobulk effect-size screen within each primary
#   macrophage state.
#
#   Main questions:
#     1) How strongly is Cdcp1 suppressed within TREM2-Mphi
#        relative to the transcriptome-wide treatment response?
#     2) Is Cdcp1 suppression shared across multiple Mphi states?
#     3) How do Spp1 and Tgfb1 rank?
#     4) Which genes show the most consistent Tx-associated
#        up/down shifts across all five states?
#
# STATISTICS
#   Biological sample is the unit.
#   Sham n=2, Tx n=2.
#   Effect-size ranking only; no inferential p-values.
#
# NORMALIZATION
#   Aggregate raw RNA counts by sample x Mphi state.
#   Normalize to pseudobulk CPM.
#   log2CPM = log2(CPM + 0.5)
#
# FILTER
#   Within each state, retain genes with:
#     total pseudobulk counts >= 20 across the four samples
#     AND CPM > 1 in >= 2 of 4 samples
#
# No reclustering.
# No source-RDS overwrite.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(Matrix)
  library(ggplot2)
})

set.seed(20261118)

ROOT <- "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk"

MPHI_RDS <- file.path(
  ROOT,
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_CleanB_FINAL_v4.14.5",
  "RDS",
  "Mouse_Mphi_Res2_CleanB_FINAL_annotated_v4.14.5.rds"
)

OUT <- file.path(
  ROOT,
  "Mouse_MASH_RDS",
  "Mouse_MASH_Mphi_Unbiased_Pseudobulk_TxResponse_v6.4.9.1"
)

DIR_TABLE <- file.path(OUT, "tables")
DIR_FIG <- file.path(OUT, "figures")
DIR_OBJ <- file.path(OUT, "objects")
DIR_LOG <- file.path(OUT, "logs")
DIR_PROV <- file.path(OUT, "provenance")

for (d in c(
  DIR_TABLE,
  DIR_FIG,
  DIR_OBJ,
  DIR_LOG,
  DIR_PROV
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

LOG_FILE <- file.path(
  DIR_LOG,
  "Mouse_MASH_Mphi_Unbiased_Pseudobulk_TxResponse_v6.4.9.1.log"
)

zz <- file(LOG_FILE, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")

cat("============================================================\n")
cat("Mouse MASH Mphi unbiased pseudobulk Tx-response v6.4.9.1\n")
cat("Start:", format(Sys.time()), "\n")
cat("============================================================\n\n")
cat("Patch note: v6.4.9.1 fixes tibble gene-name shadowing in focus-gene export; ranking definitions unchanged.\n\n")

stopifnot(file.exists(MPHI_RDS))

SAMPLES <- c(
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

sample_to_group <- c(
  "Sham1" = "Sham",
  "Sham20" = "Sham",
  "Tx17" = "Tx",
  "Tx5" = "Tx"
)

PRIMARY_MPHI_STATES <- c(
  "Anti-inflammatory-Mphi",
  "ECM-associated inflammatory-Mphi",
  "Inflammatory-Mphi",
  "Lipid-associated/TREM2-Mphi",
  "Repair/Resolution-Mphi"
)

TREM2_STATE <- "Lipid-associated/TREM2-Mphi"

FOCUS_GENES <- c(
  "Cdcp1",
  "Spp1",
  "Tgfb1",
  "Socs3",
  "Bcl3",
  "Il1rn",
  "Sbno2",
  "Nfil3",
  "Cebpb",
  "Dusp1",
  "Tnfaip3",
  "Nfkbia",
  "Zfp36",
  "Tnf",
  "Il1b",
  "Ccl2",
  "Nlrp3"
)

STATE_COL_PRIMARY <- "macrophage_class_Res2_FINAL_v4145_char"
STATE_COL_FALLBACK <- "macrophage_class_Res2_FINAL_v4145"

SAMPLE_CANDIDATES <- c(
  "sample_4group",
  "sample_for_annotation",
  "sample_FIXED2",
  "sample"
)

# ------------------------------------------------------------
# 1. Load and resolve metadata
# ------------------------------------------------------------

mphi <- readRDS(MPHI_RDS)
DefaultAssay(mphi) <- "RNA"

meta_cols <- colnames(mphi[[]])

state_col <- if (
  STATE_COL_PRIMARY %in% meta_cols
) {
  STATE_COL_PRIMARY
} else if (
  STATE_COL_FALLBACK %in% meta_cols
) {
  STATE_COL_FALLBACK
} else {
  stop("Final Mphi state column not found.")
}

sample_col <- SAMPLE_CANDIDATES[
  SAMPLE_CANDIDATES %in% meta_cols
][1]

if (
  length(sample_col) == 0 ||
  is.na(sample_col)
) {
  stop("Supported sample column not found.")
}

md <- mphi[[]] %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    sample =
      as.character(
        .data[[sample_col]]
      ),
    sender_state =
      as.character(
        .data[[state_col]]
      )
  ) %>%
  filter(
    sample %in% SAMPLES,
    sender_state %in%
      PRIMARY_MPHI_STATES
  ) %>%
  mutate(
    group =
      unname(
        sample_to_group[sample]
      )
  )

if (!all(SAMPLES %in% unique(md$sample))) {
  stop("Missing biological sample.")
}

if (!all(PRIMARY_MPHI_STATES %in% unique(md$sender_state))) {
  stop("Missing primary Mphi state.")
}

# ------------------------------------------------------------
# 2. Aggregate raw counts sample x state
# ------------------------------------------------------------

counts <- GetAssayData(
  mphi,
  assay = "RNA",
  layer = "counts"
)

counts <- counts[
  ,
  md$cell,
  drop = FALSE
]

sample_state_id <- paste(
  md$sample,
  md$sender_state,
  sep = "|||"
)

ids <- unique(sample_state_id)

pb_list <- lapply(
  ids,
  function(id) {

    cells_id <- md$cell[
      sample_state_id ==
        id
    ]

    Matrix::rowSums(
      counts[
        ,
        cells_id,
        drop = FALSE
      ]
    )
  }
)

pb_counts <- do.call(
  cbind,
  pb_list
)

rownames(pb_counts) <- rownames(counts)
colnames(pb_counts) <- ids

pb_meta <- tibble(
  sample_state_id = ids
) %>%
  separate(
    sample_state_id,
    into = c(
      "sample",
      "sender_state"
    ),
    sep = "\\|\\|\\|",
    remove = FALSE
  ) %>%
  mutate(
    group =
      unname(
        sample_to_group[sample]
      ),
    cells =
      vapply(
        sample_state_id,
        function(id) {
          sum(
            sample_state_id ==
              id
          )
        },
        numeric(1)
      )
  )

# correct cell counts using md rather than self-referential vector
pb_meta$cells <- vapply(
  pb_meta$sample_state_id,
  function(id) {
    sum(
      sample_state_id ==
        id
    )
  },
  numeric(1)
)

write_tsv(
  pb_meta,
  file.path(
    DIR_TABLE,
    "01_pseudobulk_sample_state_inventory.tsv"
  )
)

# ------------------------------------------------------------
# 3. CPM normalization
# ------------------------------------------------------------

pb_libsize <- colSums(pb_counts)

if (any(
  !is.finite(pb_libsize) |
    pb_libsize <= 0
)) {
  stop("Invalid pseudobulk library size.")
}

pb_cpm <- t(
  t(pb_counts) /
    pb_libsize
) * 1e6

pb_log2cpm <- log2(
  pb_cpm + 0.5
)

# ------------------------------------------------------------
# 4. Per-state effect-size ranking
# ------------------------------------------------------------

rank_one_state <- function(
  state_name
) {

  cols <- pb_meta %>%
    filter(
      sender_state ==
        state_name
    )

  if (nrow(cols) != 4) {
    stop(
      "Expected 4 pseudobulk samples for state: ",
      state_name
    )
  }

  mat_count <- pb_counts[
    ,
    cols$sample_state_id,
    drop = FALSE
  ]

  mat_cpm <- pb_cpm[
    ,
    cols$sample_state_id,
    drop = FALSE
  ]

  mat_log <- pb_log2cpm[
    ,
    cols$sample_state_id,
    drop = FALSE
  ]

  sham_cols <- cols$sample_state_id[
    cols$group ==
      "Sham"
  ]

  tx_cols <- cols$sample_state_id[
    cols$group ==
      "Tx"
  ]

  total_counts <-
    rowSums(
      mat_count
    )

  cpm_gt1_n <-
    rowSums(
      mat_cpm > 1
    )

  keep <-
    total_counts >= 20 &
    cpm_gt1_n >= 2

  mean_sham_log2cpm <-
    rowMeans(
      mat_log[
        ,
        sham_cols,
        drop = FALSE
      ]
    )

  mean_tx_log2cpm <-
    rowMeans(
      mat_log[
        ,
        tx_cols,
        drop = FALSE
      ]
    )

  mean_sham_cpm <-
    rowMeans(
      mat_cpm[
        ,
        sham_cols,
        drop = FALSE
      ]
    )

  mean_tx_cpm <-
    rowMeans(
      mat_cpm[
        ,
        tx_cols,
        drop = FALSE
      ]
    )

  out <- tibble(
    gene =
      rownames(mat_log),
    sender_state =
      state_name,
    total_counts =
      total_counts,
    cpm_gt1_n =
      cpm_gt1_n,
    retained =
      keep,
    mean_CPM_Sham =
      mean_sham_cpm,
    mean_CPM_Tx =
      mean_tx_cpm,
    mean_log2CPM_Sham =
      mean_sham_log2cpm,
    mean_log2CPM_Tx =
      mean_tx_log2cpm,
    delta_log2CPM =
      mean_tx_log2cpm -
      mean_sham_log2cpm
  ) %>%
    filter(retained) %>%
    arrange(delta_log2CPM) %>%
    mutate(
      rank_down =
        row_number(),
      rank_up =
        n() -
        row_number() +
        1,
      percentile_down =
        rank_down /
        n(),
      percentile_up =
        rank_up /
        n()
    )

  out
}

ranked <- bind_rows(
  lapply(
    PRIMARY_MPHI_STATES,
    rank_one_state
  )
)

write_tsv(
  ranked,
  file.path(
    DIR_TABLE,
    "02_ALL_state_gene_effect_size_rankings.tsv.gz"
  )
)

# ------------------------------------------------------------
# 5. TREM2-focused ranking
# ------------------------------------------------------------

trem2_rank <- ranked %>%
  filter(
    sender_state ==
      TREM2_STATE
  ) %>%
  arrange(
    delta_log2CPM
  )

write_tsv(
  trem2_rank,
  file.path(
    DIR_TABLE,
    "03_TREM2_Mphi_gene_effect_size_ranking.tsv"
  )
)

# ------------------------------------------------------------
# 6. Focus-gene ranks
# ------------------------------------------------------------

focus_rank <- ranked %>%
  filter(
    gene %in%
      FOCUS_GENES
  ) %>%
  arrange(
    gene,
    sender_state
  )

write_tsv(
  focus_rank,
  file.path(
    DIR_TABLE,
    "04_focus_gene_rank_across_states.tsv"
  )
)

# ------------------------------------------------------------
# 7. Cross-state consistency
# ------------------------------------------------------------

shared <- ranked %>%
  select(
    gene,
    sender_state,
    delta_log2CPM
  ) %>%
  pivot_wider(
    names_from =
      sender_state,
    values_from =
      delta_log2CPM
  )

state_cols <- intersect(
  PRIMARY_MPHI_STATES,
  colnames(shared)
)

shared <- shared %>%
  rowwise() %>%
  mutate(
    states_available =
      sum(
        !is.na(
          c_across(
            all_of(
              state_cols
            )
          )
        )
      ),

    states_down =
      sum(
        c_across(
          all_of(
            state_cols
          )
        ) < 0,
        na.rm = TRUE
      ),

    states_up =
      sum(
        c_across(
          all_of(
            state_cols
          )
        ) > 0,
        na.rm = TRUE
      ),

    median_delta_log2CPM =
      median(
        c_across(
          all_of(
            state_cols
          )
        ),
        na.rm = TRUE
      ),

    mean_delta_log2CPM =
      mean(
        c_across(
          all_of(
            state_cols
          )
        ),
        na.rm = TRUE
      )
  ) %>%
  ungroup() %>%
  filter(
    states_available >= 4
  ) %>%
  arrange(
    median_delta_log2CPM
  ) %>%
  mutate(
    cross_state_down_rank =
      row_number(),
    cross_state_up_rank =
      n() -
      row_number() +
      1
  )

write_tsv(
  shared,
  file.path(
    DIR_TABLE,
    "05_cross_state_consistency_ranking.tsv"
  )
)

# ------------------------------------------------------------
# 8. Top unbiased genes
# ------------------------------------------------------------

top_by_state <- ranked %>%
  group_by(
    sender_state
  ) %>%
  arrange(
    delta_log2CPM,
    .by_group = TRUE
  ) %>%
  summarise(
    top_down_20 =
      paste(
        head(
          gene,
          20
        ),
        collapse = ";"
      ),
    top_up_20 =
      paste(
        tail(
          gene,
          20
        ),
        collapse = ";"
      ),
    genes_ranked = n(),
    .groups = "drop"
  )

write_tsv(
  top_by_state,
  file.path(
    DIR_TABLE,
    "06_top20_up_down_by_state.tsv"
  )
)

shared_top <- bind_rows(
  shared %>%
    slice_head(
      n = 50
    ) %>%
    mutate(
      direction =
        "shared_down"
    ),
  shared %>%
    slice_tail(
      n = 50
    ) %>%
    arrange(
      desc(
        median_delta_log2CPM
      )
    ) %>%
    mutate(
      direction =
        "shared_up"
    )
)

write_tsv(
  shared_top,
  file.path(
    DIR_TABLE,
    "07_top50_cross_state_shared_up_down.tsv"
  )
)

# ------------------------------------------------------------
# 9. Cdcp1 summary
# ------------------------------------------------------------

cdcp1_summary <- focus_rank %>%
  filter(
    gene ==
      "Cdcp1"
  ) %>%
  left_join(
    shared %>%
      filter(
        gene ==
          "Cdcp1"
      ) %>%
      select(
        gene,
        states_available,
        states_down,
        states_up,
        median_delta_log2CPM,
        mean_delta_log2CPM,
        cross_state_down_rank,
        cross_state_up_rank
      ),
    by = "gene"
  ) %>%
  arrange(
    sender_state
  )

write_tsv(
  cdcp1_summary,
  file.path(
    DIR_TABLE,
    "08_PRIMARY_Cdcp1_transcriptome_rank_summary.tsv"
  )
)

# ------------------------------------------------------------
# 10. Per-sample focus-gene pseudobulk values
# ------------------------------------------------------------

focus_available <- intersect(
  FOCUS_GENES,
  rownames(pb_log2cpm)
)

focus_sample <- bind_rows(
  lapply(
    focus_available,
    function(gene_name) {

      cpm_values <-
        as.numeric(
          pb_cpm[
            gene_name,
            ,
            drop = TRUE
          ]
        )

      log2cpm_values <-
        as.numeric(
          pb_log2cpm[
            gene_name,
            ,
            drop = TRUE
          ]
        )

      tibble(
        sample_state_id =
          colnames(pb_log2cpm),
        gene =
          rep(
            gene_name,
            ncol(pb_log2cpm)
          ),
        CPM =
          cpm_values,
        log2CPM =
          log2cpm_values
      )
    }
  )
) %>%
  left_join(
    pb_meta,
    by = "sample_state_id"
  ) %>%
  select(
    gene,
    sample,
    group,
    sender_state,
    cells,
    CPM,
    log2CPM
  ) %>%
  arrange(
    gene,
    sender_state,
    sample
  )

write_tsv(
  focus_sample,
  file.path(
    DIR_TABLE,
    "09_focus_gene_pseudobulk_by_sample_state.tsv"
  )
)

# ------------------------------------------------------------
# 10. Focus-gene filter-status audit
#
# Keep the unbiased ranking filter unchanged, but report why
# a focus gene may be absent from a given state's ranked table.
# ------------------------------------------------------------

focus_filter_status <- bind_rows(
  lapply(
    PRIMARY_MPHI_STATES,
    function(state_name) {

      cols <- pb_meta %>%
        filter(
          sender_state ==
            state_name
        )

      mat_count <- pb_counts[
        ,
        cols$sample_state_id,
        drop = FALSE
      ]

      mat_cpm <- pb_cpm[
        ,
        cols$sample_state_id,
        drop = FALSE
      ]

      genes_here <- intersect(
        FOCUS_GENES,
        rownames(mat_count)
      )

      tibble(
        gene =
          genes_here,
        sender_state =
          state_name,
        total_counts =
          rowSums(
            mat_count[
              genes_here,
              ,
              drop = FALSE
            ]
          ),
        cpm_gt1_n =
          rowSums(
            mat_cpm[
              genes_here,
              ,
              drop = FALSE
            ] > 1
          )
      ) %>%
        mutate(
          retained_for_unbiased_rank =
            total_counts >= 20 &
            cpm_gt1_n >= 2
        )
    }
  )
) %>%
  arrange(
    gene,
    sender_state
  )

write_tsv(
  focus_filter_status,
  file.path(
    DIR_TABLE,
    "10_focus_gene_filter_status_across_states.tsv"
  )
)

# ------------------------------------------------------------
# 11. Figure: TREM2 rank plot
# ------------------------------------------------------------

plot_trem2 <- trem2_rank %>%
  mutate(
    highlight =
      gene %in%
        FOCUS_GENES
  )

p1 <- ggplot(
  plot_trem2,
  aes(
    x = rank_down,
    y = delta_log2CPM
  )
) +
  geom_point(
    size = 0.7,
    alpha = 0.5
  ) +
  geom_point(
    data =
      plot_trem2 %>%
      filter(highlight),
    size = 2
  ) +
  geom_text(
    data =
      plot_trem2 %>%
      filter(highlight),
    aes(
      label = gene
    ),
    size = 3,
    nudge_y = 0.15,
    check_overlap = TRUE
  ) +
  geom_hline(
    yintercept = 0
  ) +
  theme_classic(
    base_size = 10
  ) +
  labs(
    title =
      "TREM2-Mphi unbiased pseudobulk Tx-response ranking",
    x =
      "Down-regulation rank",
    y =
      "Tx - Sham mean log2CPM"
  )

ggsave(
  file.path(
    DIR_FIG,
    "01_TREM2_unbiased_gene_rank.pdf"
  ),
  p1,
  width = 9,
  height = 6
)

# ------------------------------------------------------------
# 12. Save compact results
# ------------------------------------------------------------

saveRDS(
  list(
    pseudobulk_meta =
      pb_meta,
    focus_rank =
      focus_rank,
    cross_state =
      shared,
    cdcp1_summary =
      cdcp1_summary,
    focus_sample =
      focus_sample,
    focus_filter_status =
      focus_filter_status
  ),
  file.path(
    DIR_OBJ,
    "Mouse_MASH_Mphi_Unbiased_Pseudobulk_TxResponse_v6.4.9.1_results.rds"
  )
)

# ------------------------------------------------------------
# 13. Integrity / provenance
# ------------------------------------------------------------

integrity <- tibble(
  item = c(
    "Mphi_source",
    "samples_n",
    "primary_Mphi_states_n",
    "pseudobulk_profiles_expected",
    "pseudobulk_profiles_observed",
    "normalization",
    "filter_total_counts",
    "filter_CPM",
    "inferential_p_values_computed",
    "source_RDS_overwritten",
    "reclustering_performed"
  ),
  value = c(
    MPHI_RDS,
    as.character(
      length(SAMPLES)
    ),
    as.character(
      length(
        PRIMARY_MPHI_STATES
      )
    ),
    as.character(
      length(SAMPLES) *
        length(
          PRIMARY_MPHI_STATES
        )
    ),
    as.character(
      nrow(pb_meta)
    ),
    "raw-count pseudobulk CPM; log2(CPM+0.5)",
    ">=20 across four samples within state",
    ">1 CPM in >=2/4 samples within state",
    "FALSE",
    "FALSE",
    "FALSE"
  )
)

write_tsv(
  integrity,
  file.path(
    DIR_PROV,
    "integrity_summary_v6.4.9.1.tsv"
  )
)

cat("\n")
cat("============================================================\n")
cat("PSEUDOBULK INVENTORY\n")
cat("============================================================\n")
print(pb_meta)

cat("\n")
cat("============================================================\n")
cat("FOCUS GENE RANKS\n")
cat("============================================================\n")
print(focus_rank)

cat("\n")
cat("============================================================\n")
cat("CDCP1 TRANSCRIPTOME RANK SUMMARY\n")
cat("============================================================\n")
print(cdcp1_summary)

cat("\n")
cat("============================================================\n")
cat("TOP CROSS-STATE SHARED DOWN\n")
cat("============================================================\n")
print(
  shared %>%
    slice_head(
      n = 30
    ) %>%
    select(
      gene,
      states_available,
      states_down,
      states_up,
      median_delta_log2CPM,
      mean_delta_log2CPM,
      cross_state_down_rank
    )
)

cat("\n")
cat("============================================================\n")
cat("TOP TREM2 DOWN\n")
cat("============================================================\n")
print(
  trem2_rank %>%
    slice_head(
      n = 30
    ) %>%
    select(
      gene,
      mean_CPM_Sham,
      mean_CPM_Tx,
      delta_log2CPM,
      rank_down,
      percentile_down
    )
)

cat("\n")
cat("Output directory:\n")
cat(OUT, "\n\n")

cat("STATUS: SUCCESS\n")
cat("ANALYSIS COMPLETED\n")
cat("End:", format(Sys.time()), "\n")

sink(type = "message")
sink(type = "output")
close(zz)
