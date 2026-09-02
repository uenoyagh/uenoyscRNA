suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20260902)

VERSION <- "v6.7.5"
RES_COL <- "LSECclean_res0.3"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_v6.7.4.2/objects/",
  "Mouse_MASH_LSEC_sample_specific_audit_v6.7.4.2.rds"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_", VERSION
)

OBJDIR <- file.path(OUTDIR, "objects")
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")

dir.create(OBJDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(TABDIR, recursive=TRUE, showWarnings=FALSE)

cat("====================================================\n")
cat("Mouse MASH LSEC annotation freeze\n")
cat("Version:", VERSION, "\n")
cat("Final annotation scaffold: res0.3\n")
cat("====================================================\n\n")

obj <- readRDS(INPUT_RDS)

if (!RES_COL %in% colnames(obj@meta.data)) {
  stop("Missing ", RES_COL)
}

if (!"umapLSECclean" %in% Reductions(obj)) {
  stop("Missing umapLSECclean")
}

cl <- as.character(obj@meta.data[[RES_COL]])

annotation_map <- c(
  "0" = "Inflammatory_stress_high_LSEC",
  "1" = "Homeostatic_like_LSEC",
  "2" = "Tx5_enriched_LSEC_state",
  "3" = "Wnt_angiocrine_high_LSEC",
  "4" = "Low_quality_ambient_enriched_LSEC",
  "5" = "Tx17_enriched_Cd209b_Ctsj_LSEC",
  "6" = "Cycling_LSEC"
)

evidence_map <- c(
  "0" = "cross_sample_primary",
  "1" = "cross_sample_primary",
  "2" = "sample_specific_exploratory",
  "3" = "cross_sample_primary",
  "4" = "QC_flagged",
  "5" = "sample_specific_exploratory",
  "6" = "rare_cross_sample"
)

if (!all(unique(cl) %in% names(annotation_map))) {
  stop("Unexpected res0.3 cluster.")
}

obj$LSEC_state_v675 <-
  unname(annotation_map[cl])

obj$LSEC_state_evidence_v675 <-
  unname(evidence_map[cl])

obj$LSEC_QC_flag_v675 <-
  ifelse(
    cl == "4",
    "QC_flagged",
    "Pass"
  )

obj$LSEC_sample_specific_flag_v675 <-
  ifelse(
    cl %in% c("2","5"),
    "Sample_specific_exploratory",
    "No"
  )

cat("=== ANNOTATION COUNTS ===\n")
print(
  table(
    obj$LSEC_state_v675
  )
)

cat("\n=== EVIDENCE CLASS ===\n")
print(
  table(
    obj$LSEC_state_evidence_v675
  )
)

# ---------------------------------------------------------
# Annotation table
# ---------------------------------------------------------

decision_table <- data.frame(
  cluster=names(annotation_map),
  state=unname(annotation_map),
  evidence=unname(evidence_map),
  stringsAsFactors=FALSE
)

decision_table$n_cells <- as.integer(
  table(
    factor(
      cl,
      levels=names(annotation_map)
    )
  )
)

write.csv(
  decision_table,
  file.path(
    TABDIR,
    "LSEC_res0.3_annotation_decisions_v6.7.5.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# State composition
# ---------------------------------------------------------

state_by_sample <- as.data.frame(
  table(
    sample=obj$sample,
    condition=obj$condition,
    state=obj$LSEC_state_v675
  )
)

sample_totals <- aggregate(
  Freq ~ sample,
  data=state_by_sample,
  FUN=sum
)

colnames(sample_totals)[2] <- "sample_total"

state_by_sample <- merge(
  state_by_sample,
  sample_totals,
  by="sample"
)

state_by_sample$fraction <-
  state_by_sample$Freq /
  state_by_sample$sample_total

write.csv(
  state_by_sample,
  file.path(
    TABDIR,
    "LSEC_state_fraction_by_sample_v6.7.5.csv"
  ),
  row.names=FALSE
)

p_fraction <- ggplot(
  state_by_sample,
  aes(
    x=sample,
    y=fraction,
    fill=state
  )
) +
  geom_col() +
  theme_classic(base_size=11) +
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1
    )
  ) +
  labs(
    title="LSEC state composition by sample",
    x=NULL,
    y="Fraction of LSEC cells",
    fill="LSEC state"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_state_fraction_by_sample_v6.7.5.pdf"
  ),
  p_fraction,
  width=11,
  height=6.5
)

# ---------------------------------------------------------
# Annotated UMAP
# ---------------------------------------------------------

p_umap <- DimPlot(
  obj,
  reduction="umapLSECclean",
  group.by="LSEC_state_v675",
  label=TRUE,
  repel=TRUE,
  pt.size=0.45
) +
  theme_classic(base_size=12) +
  theme(
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    legend.title=element_blank()
  ) +
  ggtitle(
    "Mouse MASH LSEC states - v6.7.5"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_annotated_UMAP_v6.7.5.pdf"
  ),
  p_umap,
  width=11,
  height=7
)

# ---------------------------------------------------------
# Module definitions
# ---------------------------------------------------------

DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(
  obj,
  assay="RNA",
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

modules <- list(

  LSEC_identity = c(
    "Clec4g","Stab1","Stab2","Lyve1",
    "Fcgr2b","Mrc1","Oit3","Dnase1l3"
  ),

  Wnt_angiocrine = c(
    "Wnt2","Wnt9b","Rspo3","Bmp2"
  ),

  Capillarization = c(
    "Cd34","Vwf","Plvap","Emcn","Esm1"
  ),

  Inflammatory_adhesion = c(
    "Icam1","Vcam1","Sele","Selp",
    "Ccl2","Cxcl9","Cxcl10"
  ),

  IFN_response = c(
    "Isg15","Ifit1","Ifit2","Ifit3",
    "Irf7","Stat1"
  ),

  Angiogenic = c(
    "Apln","Aplnr","Kdr",
    "Pgf","Angpt2","Esm1"
  ),

  Stress_IEG = c(
    "Fos","Fosb","Jun","Junb",
    "Jund","Atf3","Egr1",
    "Hspa1a","Hspa1b"
  ),

  MHCII = c(
    "H2-Aa","H2-Ab1","H2-Eb1","Cd74"
  ),

  Cycling = c(
    "Mki67","Pcna","Mcm5","Cdt1",
    "Rad51","Brca1"
  )
)

module_cols <- character()

for (nm in names(modules)) {

  genes <- intersect(
    modules[[nm]],
    rownames(obj)
  )

  if (length(genes) < 2) {
    next
  }

  temp_name <- paste0(
    "TEMP_", nm, "_"
  )

  obj <- AddModuleScore(
    obj,
    features=list(genes),
    assay="RNA",
    name=temp_name,
    seed=20260902
  )

  old_col <- paste0(
    temp_name,
    "1"
  )

  new_col <- paste0(
    "MS_",
    nm,
    "_v675"
  )

  obj@meta.data[[new_col]] <-
    obj@meta.data[[old_col]]

  obj@meta.data[[old_col]] <- NULL

  module_cols <- c(
    module_cols,
    new_col
  )
}

# ---------------------------------------------------------
# Sample-level module summaries
# ---------------------------------------------------------

md <- obj@meta.data

sample_summary <- do.call(
  rbind,
  lapply(
    unique(as.character(md$sample)),
    function(s) {

      x <- md[
        as.character(md$sample) == s,
        ,
        drop=FALSE
      ]

      out <- data.frame(
        sample=s,
        condition=as.character(x$condition[1]),
        n_cells=nrow(x),
        stringsAsFactors=FALSE
      )

      for (m in module_cols) {

        out[[paste0(m,"_median")]] <-
          median(
            x[[m]],
            na.rm=TRUE
          )

        out[[paste0(m,"_mean")]] <-
          mean(
            x[[m]],
            na.rm=TRUE
          )
      }

      out
    }
  )
)

write.csv(
  sample_summary,
  file.path(
    TABDIR,
    "LSEC_module_scores_by_sample_v6.7.5.csv"
  ),
  row.names=FALSE
)

# ---------------------------------------------------------
# Sham vs Tx sample-level figure
# ---------------------------------------------------------

shamtx <- sample_summary[
  sample_summary$condition %in%
    c("Sham","Tx"),
  ,
  drop=FALSE
]

median_cols <- grep(
  "_median$",
  colnames(shamtx),
  value=TRUE
)

long <- do.call(
  rbind,
  lapply(
    median_cols,
    function(m) {

      data.frame(
        sample=shamtx$sample,
        condition=shamtx$condition,
        module=sub(
          "^MS_",
          "",
          sub(
            "_v675_median$",
            "",
            m
          )
        ),
        score=shamtx[[m]],
        stringsAsFactors=FALSE
      )
    }
  )
)

p_modules <- ggplot(
  long,
  aes(
    x=condition,
    y=score,
    group=sample
  )
) +
  geom_point(
    size=2.8
  ) +
  facet_wrap(
    ~module,
    scales="free_y",
    ncol=3
  ) +
  theme_classic(base_size=11) +
  labs(
    title="LSEC programs: Sham vs Tx",
    x=NULL,
    y="Sample median module score"
  )

ggsave(
  file.path(
    FIGDIR,
    "LSEC_module_scores_Shams_vs_Tx_sample_level_v6.7.5.pdf"
  ),
  p_modules,
  width=11,
  height=9
)

# ---------------------------------------------------------
# Preserve annotated object
# ---------------------------------------------------------

saveRDS(
  obj,
  file.path(
    OBJDIR,
    "Mouse_MASH_LSEC_annotated_v6.7.5.rds"
  ),
  compress=FALSE
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.7.5.txt"
  )
)

cat("\n====================================================\n")
cat("v6.7.5 COMPLETE\n")
cat("Cells retained:", ncol(obj), "\n")
cat("Final annotation scaffold: res0.3\n")
cat("No cells removed in v6.7.5\n")
cat("Cluster 4 remains QC-flagged\n")
cat("Clusters 2 and 5 remain sample-specific exploratory\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
