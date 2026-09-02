suppressPackageStartupMessages({
  library(ggplot2)
})

VERSION <- "v6.7.8"

BASE <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/"
)

V675 <- file.path(
  BASE,
  "Mouse_MASH_LSEC_v6.7.5"
)

V677 <- file.path(
  BASE,
  "Mouse_MASH_LSEC_v6.7.7"
)

OUTDIR <- file.path(
  BASE,
  paste0("Mouse_MASH_LSEC_", VERSION)
)

FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH LSEC mechanistic closure\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

# =========================================================
# 1. Sample-level module scores
# =========================================================

module_file <- file.path(
  V675,
  "tables",
  "LSEC_module_scores_by_sample_v6.7.5.csv"
)

if (!file.exists(module_file)) {
  stop("Missing module-score file: ", module_file)
}

ms <- read.csv(
  module_file,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

sample_order <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

ms$sample <- factor(
  ms$sample,
  levels=sample_order
)

median_cols <- grep(
  "^MS_.*_v675_median$",
  colnames(ms),
  value=TRUE
)

if (length(median_cols) == 0) {
  stop("No median module-score columns found.")
}

module_name <- function(x) {

  x <- sub(
    "^MS_",
    "",
    x
  )

  x <- sub(
    "_v675_median$",
    "",
    x
  )

  x
}

modules <- vapply(
  median_cols,
  module_name,
  character(1)
)

# =========================================================
# 2. Long table + within-module Z score
# =========================================================

long_list <- list()

for (i in seq_along(median_cols)) {

  col <- median_cols[i]
  mod <- modules[i]

  x <- data.frame(
    sample=as.character(ms$sample),
    condition=ms$condition,
    module=mod,
    score=ms[[col]],
    stringsAsFactors=FALSE
  )

  s <- sd(
    x$score,
    na.rm=TRUE
  )

  if (
    is.na(s) ||
    s == 0
  ) {

    x$zscore <- 0

  } else {

    x$zscore <-
      (
        x$score -
        mean(
          x$score,
          na.rm=TRUE
        )
      ) / s
  }

  long_list[[mod]] <- x
}

long <- do.call(
  rbind,
  long_list
)

rownames(long) <- NULL

long$sample <- factor(
  long$sample,
  levels=sample_order
)

write.csv(
  long,
  file.path(
    TABDIR,
    "LSEC_module_scores_long_all6samples_v6.7.8.csv"
  ),
  row.names=FALSE
)

# =========================================================
# 3. Six-sample module heatmap
# =========================================================

p_heat <- ggplot(
  long,
  aes(
    x=sample,
    y=module,
    fill=zscore
  )
) +
  geom_tile(
    color="grey85",
    linewidth=0.3
  ) +
  geom_text(
    aes(
      label=sprintf(
        "%.2f",
        score
      )
    ),
    size=2.7
  ) +
  scale_fill_gradient2(
    low="#0033FF",
    mid="#FFFFFF",
    high="#FF1A1A",
    midpoint=0
  ) +
  theme_classic(base_size=11) +
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1
    ),
    axis.title=element_blank()
  ) +
  labs(
    title=
      "Mouse MASH LSEC programs across six samples",
    subtitle=
      "Fill = within-module Z score; text = raw median module score",
    fill="Z score"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_module_heatmap_all6samples_v6.7.8.pdf"
  ),
  p_heat,
  width=9,
  height=7
)

# =========================================================
# 4. Disease-axis vs treatment-axis summary
#
# CDHFD vs STD:
# descriptive only because n=1 vs n=1
#
# Tx vs Sham:
# mean of two biological samples/group
# =========================================================

contrast <- data.frame()

for (i in seq_along(median_cols)) {

  col <- median_cols[i]
  mod <- modules[i]

  get_value <- function(sample_name) {

    z <- ms[
      as.character(ms$sample) ==
        sample_name,
      col
    ]

    if (length(z) != 1) {
      stop(
        "Unexpected sample count for ",
        sample_name
      )
    }

    as.numeric(z)
  }

  STD <- get_value("STD_rep1")
  CDHFD <- get_value("CDHFD_rep1")

  sham_values <- c(
    get_value("Sham1"),
    get_value("Sham20")
  )

  tx_values <- c(
    get_value("Tx17"),
    get_value("Tx5")
  )

  disease_delta <- CDHFD - STD

  sham_mean <- mean(
    sham_values,
    na.rm=TRUE
  )

  tx_mean <- mean(
    tx_values,
    na.rm=TRUE
  )

  tx_delta <- tx_mean -
    sham_mean

  direction_class <- if (
    disease_delta == 0 ||
    tx_delta == 0
  ) {

    "Near_zero_axis"

  } else if (
    sign(disease_delta) !=
      sign(tx_delta)
  ) {

    "Opposite_direction"

  } else {

    "Same_direction"
  }

  contrast <- rbind(
    contrast,
    data.frame(
      module=mod,

      STD=STD,
      CDHFD=CDHFD,
      disease_delta_CDHFD_minus_STD=
        disease_delta,

      Sham1=sham_values[1],
      Sham20=sham_values[2],
      Sham_mean=sham_mean,

      Tx17=tx_values[1],
      Tx5=tx_values[2],
      Tx_mean=tx_mean,

      tx_delta_Tx_minus_Sham=
        tx_delta,

      direction_class=
        direction_class,

      stringsAsFactors=FALSE
    )
  )
}

write.csv(
  contrast,
  file.path(
    TABDIR,
    "LSEC_disease_vs_treatment_module_direction_v6.7.8.csv"
  ),
  row.names=FALSE
)

# =========================================================
# 5. Disease-vs-treatment direction plot
# =========================================================

p_axis <- ggplot(
  contrast,
  aes(
    x=disease_delta_CDHFD_minus_STD,
    y=tx_delta_Tx_minus_Sham,
    label=module
  )
) +
  geom_hline(
    yintercept=0,
    linetype=2
  ) +
  geom_vline(
    xintercept=0,
    linetype=2
  ) +
  geom_point(
    size=3
  ) +
  geom_text(
    size=3,
    nudge_y=0.015,
    check_overlap=TRUE
  ) +
  theme_classic(base_size=12) +
  labs(
    title=
      "LSEC disease-axis vs transplantation-axis programs",
    subtitle=paste0(
      "X: CDHFD−STD is descriptive (n=1 each); ",
      "Y: mean Tx−mean Sham (n=2 each)"
    ),
    x=
      "Disease axis: CDHFD − STD",
    y=
      "Transplantation axis: Tx − Sham"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_disease_vs_treatment_program_direction_v6.7.8.pdf"
  ),
  p_axis,
  width=8.5,
  height=7
)

# =========================================================
# 6. Ambient-aware Hallmark consensus
# =========================================================

hallmark_file <- file.path(
  V677,
  "tables",
  "Hallmark_consensus_exploratory_padj025_v6.7.7.csv"
)

if (!file.exists(hallmark_file)) {
  stop(
    "Missing Hallmark consensus: ",
    hallmark_file
  )
}

hall <- read.csv(
  hallmark_file,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

if (
  all(
    c(
      "NES_Primary",
      "NES_Shared"
    ) %in% colnames(hall)
  )
) {

  hall$mean_NES_recalc <-
    (
      hall$NES_Primary +
      hall$NES_Shared
    ) / 2

  hall <- hall[
    order(
      -abs(
        hall$mean_NES_recalc
      )
    ),
    ,
    drop=FALSE
  ]
}

write.csv(
  hall,
  file.path(
    TABDIR,
    "LSEC_Hallmark_consensus_closure_v6.7.8.csv"
  ),
  row.names=FALSE
)

# =========================================================
# 7. GO BP consensus
# =========================================================

gobp_file <- file.path(
  V677,
  "tables",
  "GOBP_consensus_exploratory_padj025_v6.7.7.csv"
)

if (
  file.exists(gobp_file)
) {

  gobp <- read.csv(
    gobp_file,
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  if (
    all(
      c(
        "NES_Primary",
        "NES_Shared"
      ) %in% colnames(gobp)
    )
  ) {

    gobp$mean_NES_recalc <-
      (
        gobp$NES_Primary +
        gobp$NES_Shared
      ) / 2

    gobp <- gobp[
      order(
        -abs(
          gobp$mean_NES_recalc
        )
      ),
      ,
      drop=FALSE
    ]
  }

  write.csv(
    gobp,
    file.path(
      TABDIR,
      "LSEC_GOBP_consensus_closure_v6.7.8.csv"
    ),
    row.names=FALSE
  )
}

# =========================================================
# 8. LSEC-supported cross-analysis genes
# =========================================================

candidate_file <- file.path(
  V677,
  "tables",
  "LSEC_supported_cross_analysis_candidates_v6.7.7.csv"
)

if (!file.exists(candidate_file)) {
  stop(
    "Missing candidate table: ",
    candidate_file
  )
}

cand <- read.csv(
  candidate_file,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

cand <- cand[
  cand$same_direction %in%
    c(TRUE, "TRUE", 1),
  ,
  drop=FALSE
]

cand <- cand[
  order(
    cand$max_FDR,
    -abs(cand$mean_logFC)
  ),
  ,
  drop=FALSE
]

cand$Tx_direction <- ifelse(
  cand$mean_logFC > 0,
  "Tx_up",
  "Tx_down"
)

write.csv(
  cand,
  file.path(
    TABDIR,
    "LSEC_supported_cross_analysis_candidates_closure_v6.7.8.csv"
  ),
  row.names=FALSE
)

# =========================================================
# 9. Candidate concordance plot
# =========================================================

top_candidate <- head(
  cand,
  30
)

p_gene <- ggplot(
  cand,
  aes(
    x=logFC_Primary,
    y=logFC_Shared
  )
) +
  geom_hline(
    yintercept=0,
    linetype=2
  ) +
  geom_vline(
    xintercept=0,
    linetype=2
  ) +
  geom_abline(
    slope=1,
    intercept=0,
    linetype=3
  ) +
  geom_point(
    size=1.8,
    alpha=0.7
  ) +
  geom_text(
    data=top_candidate,
    aes(label=gene),
    size=2.8,
    nudge_y=0.05,
    check_overlap=TRUE
  ) +
  theme_classic(base_size=12) +
  labs(
    title=
      "LSEC-supported Tx-response candidates",
    subtitle=
      "Primary_no_QC versus Shared_core",
    x=
      "Primary logFC (Tx vs Sham)",
    y=
      "Shared-core logFC (Tx vs Sham)"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_supported_candidate_concordance_v6.7.8.pdf"
  ),
  p_gene,
  width=8,
  height=7
)

# =========================================================
# 10. Concise closure summary
# =========================================================

top_up <- cand[
  cand$mean_logFC > 0,
  ,
  drop=FALSE
]

top_down <- cand[
  cand$mean_logFC < 0,
  ,
  drop=FALSE
]

top_up <- head(
  top_up,
  10
)

top_down <- head(
  top_down,
  10
)

summary_lines <- c(

  "# Mouse MASH LSEC mechanistic closure v6.7.8",
  "",

  "## Analysis framework",
  "- Final annotation scaffold: clean LSEC RPCA resolution 0.3",
  "- Clean LSEC cells: 9,190",
  "- Primary Tx vs Sham statistical unit: biological sample",
  "- Sham n=2; Tx n=2",
  "- STD vs CDHFD module comparison is descriptive only (n=1 each)",
  "- Hepatocyte-dominant transcripts were excluded from primary pathway interpretation",
  "",

  "## Core interpretation",
  paste0(
    "- LSEC shows substantial sample-to-sample heterogeneity; ",
    "Sham and Tx do not form strongly separated pseudobulk groups."
  ),
  paste0(
    "- Tx-associated LSEC changes therefore should not be interpreted ",
    "as a uniform global normalization response."
  ),
  paste0(
    "- Ambient-aware pathway analysis identifies reproducible ",
    "LSEC remodeling programs shared by Primary and Shared-core analyses."
  ),
  paste0(
    "- Hallmark results include Tx-associated TNFA/NFKB, hypoxia, ",
    "inflammatory response and TGF-beta programs, rather than a simple ",
    "anti-inflammatory shift."
  ),
  paste0(
    "- Coagulation shows an opposite/downward Tx-associated direction ",
    "in the concordant Hallmark analysis."
  ),
  paste0(
    "- Individual gene-level evidence is exploratory because no ",
    "non-hepatocyte candidate reaches conventional strong FDR significance."
  ),
  "",

  "## Top LSEC-supported Tx-up candidates",
  paste(
    top_up$gene,
    collapse=", "
  ),
  "",

  "## Top LSEC-supported Tx-down candidates",
  paste(
    top_down$gene,
    collapse=", "
  ),
  "",

  "## Conclusion",
  paste0(
    "LSEC is best interpreted as a heterogeneous secondary remodeling ",
    "compartment after transplantation rather than a compartment showing ",
    "a strong uniform therapeutic normalization."
  ),
  paste0(
    "The strongest mechanistic conclusions of the mouse dataset remain ",
    "centered on macrophage-HSC and macrophage-hepatocyte axes; LSEC ",
    "findings are supportive/exploratory."
  )
)

writeLines(
  summary_lines,
  file.path(
    OUTDIR,
    "LSEC_mechanistic_closure_summary_v6.7.8.md"
  )
)

# =========================================================
# 11. Terminal summary
# =========================================================

cat("\n=== MODULE DIRECTION TABLE ===\n")

print(
  contrast[
    ,
    c(
      "module",
      "disease_delta_CDHFD_minus_STD",
      "tx_delta_Tx_minus_Sham",
      "direction_class"
    )
  ]
)

cat("\n=== TOP HALLMARK CONSENSUS ===\n")

print(
  head(
    hall,
    20
  )
)

cat("\n=== TOP LSEC-SUPPORTED Tx-UP ===\n")

print(
  top_up[
    ,
    c(
      "gene",
      "logFC_Primary",
      "logFC_Shared",
      "mean_logFC",
      "max_FDR"
    )
  ]
)

cat("\n=== TOP LSEC-SUPPORTED Tx-DOWN ===\n")

print(
  top_down[
    ,
    c(
      "gene",
      "logFC_Primary",
      "logFC_Shared",
      "mean_logFC",
      "max_FDR"
    )
  ]
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.7.8.txt"
  )
)

cat("\n====================================================\n")
cat("v6.7.8 COMPLETE\n")
cat("LSEC mechanistic closure summary generated\n")
cat("No source RDS modified\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
