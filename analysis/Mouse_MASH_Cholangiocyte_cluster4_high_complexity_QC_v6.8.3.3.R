suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

VERSION <- "v6.8.3.3"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.2/objects/",
  "Mouse_MASH_Cholangiocyte_res0.3_audit_v6.8.2.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_", VERSION
)

FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Cholangiocyte cluster 4 high-complexity QC audit\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing input RDS: ", INPUT_RDS)
}

obj <- readRDS(INPUT_RDS)
md <- obj@meta.data

required <- c(
  "sample",
  "condition",
  "Chol_res03_v682",
  "nCount_RNA",
  "nFeature_RNA"
)

missing_cols <- setdiff(
  required,
  colnames(md)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata: ",
    paste(missing_cols, collapse=", ")
  )
}

md$cell <- rownames(md)
md$cluster <- as.character(md$Chol_res03_v682)

# Confirmed contamination clusters from v6.8.3.2
remove_clusters <- c("3", "9", "10", "11", "12")

md$retained_v6832 <-
  !md$cluster %in% remove_clusters

md$cluster4 <-
  md$cluster == "4"

sample_order <- c(
  "STD_rep1",
  "CDHFD_rep1",
  "Sham1",
  "Sham20",
  "Tx17",
  "Tx5"
)

# ---------------------------------------------------------
# Within-sample reference thresholds
#
# Reference:
# retained Cholangiocytes excluding cluster 4.
#
# Extreme threshold:
# Q3 + 3 x IQR
# Diagnostic only; NOT an exclusion threshold yet.
# ---------------------------------------------------------

threshold_list <- list()
summary_list <- list()

for (s in sample_order) {

  ref <- md[
    md$sample == s &
      md$retained_v6832 &
      !md$cluster4,
    ,
    drop=FALSE
  ]

  c4 <- md[
    md$sample == s &
      md$cluster4,
    ,
    drop=FALSE
  ]

  if (nrow(ref) == 0) {
    stop("No reference cells for sample: ", s)
  }

  q_count <- quantile(
    ref$nCount_RNA,
    probs=c(0.25, 0.50, 0.75),
    na.rm=TRUE
  )

  q_feature <- quantile(
    ref$nFeature_RNA,
    probs=c(0.25, 0.50, 0.75),
    na.rm=TRUE
  )

  iqr_count <- q_count[3] - q_count[1]
  iqr_feature <- q_feature[3] - q_feature[1]

  count_threshold <-
    as.numeric(q_count[3] + 3 * iqr_count)

  feature_threshold <-
    as.numeric(q_feature[3] + 3 * iqr_feature)

  threshold_list[[s]] <- data.frame(
    sample=s,
    reference_n=nrow(ref),
    ref_nCount_Q1=as.numeric(q_count[1]),
    ref_nCount_median=as.numeric(q_count[2]),
    ref_nCount_Q3=as.numeric(q_count[3]),
    nCount_extreme_threshold=count_threshold,
    ref_nFeature_Q1=as.numeric(q_feature[1]),
    ref_nFeature_median=as.numeric(q_feature[2]),
    ref_nFeature_Q3=as.numeric(q_feature[3]),
    nFeature_extreme_threshold=feature_threshold,
    stringsAsFactors=FALSE
  )

  if (nrow(c4) > 0) {

    extreme_count <-
      c4$nCount_RNA > count_threshold

    extreme_feature <-
      c4$nFeature_RNA > feature_threshold

    summary_list[[s]] <- data.frame(
      sample=s,
      cluster4_n=nrow(c4),

      cluster4_median_nCount=
        median(c4$nCount_RNA),

      reference_median_nCount=
        median(ref$nCount_RNA),

      nCount_ratio=
        median(c4$nCount_RNA) /
        median(ref$nCount_RNA),

      cluster4_median_nFeature=
        median(c4$nFeature_RNA),

      reference_median_nFeature=
        median(ref$nFeature_RNA),

      nFeature_ratio=
        median(c4$nFeature_RNA) /
        median(ref$nFeature_RNA),

      extreme_nCount_n=
        sum(extreme_count),

      extreme_nCount_pct=
        mean(extreme_count) * 100,

      extreme_nFeature_n=
        sum(extreme_feature),

      extreme_nFeature_pct=
        mean(extreme_feature) * 100,

      extreme_both_n=
        sum(
          extreme_count &
            extreme_feature
        ),

      extreme_both_pct=
        mean(
          extreme_count &
            extreme_feature
        ) * 100,

      stringsAsFactors=FALSE
    )
  }
}

threshold_df <- do.call(
  rbind,
  threshold_list
)

summary_df <- do.call(
  rbind,
  summary_list
)

rownames(threshold_df) <- NULL
rownames(summary_df) <- NULL

write.csv(
  threshold_df,
  file.path(
    TABDIR,
    "Cholangiocyte_cluster4_within_sample_QC_thresholds_v6.8.3.3.csv"
  ),
  row.names=FALSE
)

write.csv(
  summary_df,
  file.path(
    TABDIR,
    "Cholangiocyte_cluster4_high_complexity_summary_v6.8.3.3.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Per-cell cluster 4 QC flags
# ---------------------------------------------------------

c4 <- md[
  md$cluster4,
  ,
  drop=FALSE
]

c4$nCount_extreme_threshold <- NA_real_
c4$nFeature_extreme_threshold <- NA_real_

for (s in sample_order) {

  th <- threshold_df[
    threshold_df$sample == s,
    ,
    drop=FALSE
  ]

  idx <- c4$sample == s

  c4$nCount_extreme_threshold[idx] <-
    th$nCount_extreme_threshold

  c4$nFeature_extreme_threshold[idx] <-
    th$nFeature_extreme_threshold
}

c4$extreme_nCount <-
  c4$nCount_RNA >
  c4$nCount_extreme_threshold

c4$extreme_nFeature <-
  c4$nFeature_RNA >
  c4$nFeature_extreme_threshold

c4$extreme_both <-
  c4$extreme_nCount &
  c4$extreme_nFeature

write.csv(
  c4,
  file.path(
    TABDIR,
    "Cholangiocyte_cluster4_per_cell_high_complexity_audit_v6.8.3.3.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Plot dataset
# ---------------------------------------------------------

plot_md <- md[
  md$retained_v6832,
  ,
  drop=FALSE
]

plot_md$QC_class <- ifelse(
  plot_md$cluster4,
  "Cluster4_QC_watch",
  "Other_retained_Cholangiocyte"
)

plot_md$sample <- factor(
  plot_md$sample,
  levels=sample_order
)

# ---------------------------------------------------------
# nCount boxplot
# ---------------------------------------------------------

p_count <- ggplot(
  plot_md,
  aes(
    x=QC_class,
    y=nCount_RNA
  )
) +
  geom_boxplot(
    outlier.shape=NA
  ) +
  scale_y_log10() +
  facet_wrap(
    ~ sample,
    scales="free_y"
  ) +
  theme_classic(base_size=10) +
  theme(
    axis.text.x=
      element_text(
        angle=45,
        hjust=1
      )
  ) +
  labs(
    title=
      "Cluster 4 RNA-count complexity vs other retained Cholangiocytes",
    x=NULL,
    y="nCount_RNA (log10 scale)"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_cluster4_nCount_by_sample_v6.8.3.3.pdf"
  ),
  p_count,
  width=13,
  height=8
)

# ---------------------------------------------------------
# nFeature boxplot
# ---------------------------------------------------------

p_feature <- ggplot(
  plot_md,
  aes(
    x=QC_class,
    y=nFeature_RNA
  )
) +
  geom_boxplot(
    outlier.shape=NA
  ) +
  scale_y_log10() +
  facet_wrap(
    ~ sample,
    scales="free_y"
  ) +
  theme_classic(base_size=10) +
  theme(
    axis.text.x=
      element_text(
        angle=45,
        hjust=1
      )
  ) +
  labs(
    title=
      "Cluster 4 feature complexity vs other retained Cholangiocytes",
    x=NULL,
    y="nFeature_RNA (log10 scale)"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_cluster4_nFeature_by_sample_v6.8.3.3.pdf"
  ),
  p_feature,
  width=13,
  height=8
)

# ---------------------------------------------------------
# nCount vs nFeature
# ---------------------------------------------------------

p_scatter <- ggplot(
  plot_md,
  aes(
    x=nCount_RNA,
    y=nFeature_RNA,
    shape=QC_class
  )
) +
  geom_point(
    alpha=0.35,
    size=0.7
  ) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(
    ~ sample
  ) +
  theme_classic(base_size=10) +
  labs(
    title=
      "Cholangiocyte RNA complexity: cluster 4 QC watch",
    x="nCount_RNA (log10)",
    y="nFeature_RNA (log10)",
    shape=NULL
  )

ggsave(
  file.path(
    FIGDIR,
    "Cholangiocyte_cluster4_nCount_vs_nFeature_v6.8.3.3.pdf"
  ),
  p_scatter,
  width=13,
  height=8
)

# ---------------------------------------------------------
# Terminal report
# ---------------------------------------------------------

cat("\n=== WITHIN-SAMPLE QC THRESHOLDS ===\n")
print(
  threshold_df,
  row.names=FALSE
)

cat("\n=== CLUSTER 4 HIGH-COMPLEXITY SUMMARY ===\n")
print(
  summary_df,
  row.names=FALSE
)

cat("\n=== CLUSTER 4 EXTREME-BOTH COUNTS ===\n")

print(
  table(
    sample=c4$sample,
    extreme_both=c4$extreme_both
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.3.3.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.3.3 COMPLETE\n")
cat("Diagnostic QC only\n")
cat("No cells removed\n")
cat("No object rewritten\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
