suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

set.seed(20260902)

VERSION <- "v6.8.0.1"

INPUT_RDS <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/",
  "RDS3_annotation_visualization_v4.1.1/objects/",
  "RDS3_with_visualization_metadata_v4.1.1.rds"
)

AUDIT_CSV <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_Cholangiocyte_v6.8.0/tables/",
  "Cycling_biliary_rescue_audit_per_cell_v6.8.0.csv"
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
cat("Cycling -> Cholangiocyte rescue validation\n")
cat("Version:", VERSION, "\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_RDS)) {
  stop("Missing whole-cell RDS.")
}

if (!file.exists(AUDIT_CSV)) {
  stop("Missing v6.8.0 rescue audit.")
}

whole <- readRDS(INPUT_RDS)

audit <- read.csv(
  AUDIT_CSV,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

audit$rescue_candidate <- as.logical(
  audit$rescue_candidate
)

candidate_cells <- audit$cell[
  audit$rescue_candidate
]

candidate_cells <- intersect(
  candidate_cells,
  colnames(whole)
)

cat("v6.8.0 rescue candidates:", length(candidate_cells), "\n")

if (length(candidate_cells) == 0) {
  stop("No rescue candidates found.")
}

cyc <- subset(
  whole,
  cells=candidate_cells
)

DefaultAssay(cyc) <- "RNA"

cyc <- NormalizeData(
  cyc,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

panels <- list(

  Biliary_core = c(
    "Krt19","Krt7","Epcam","Sox9",
    "Muc1","Hnf1b","Cftr","Slc4a2"
  ),

  Biliary_support = c(
    "Krt8","Krt18","Klf5",
    "Tacstd2","Prom1","Spp1","Mmp7"
  ),

  Cycling = c(
    "Mki67","Top2a","Birc5",
    "Ube2c","Cenpf","Pcna"
  ),

  Hepatocyte = c(
    "Alb","Ttr","Apoa1","Fabp1",
    "Cps1","Ass1","Hnf4a"
  ),

  Myeloid = c(
    "Ptprc","Lyz2","Adgre1",
    "Tyrobp","Aif1"
  ),

  Endothelial = c(
    "Pecam1","Cdh5","Kdr",
    "Klf2","Stab2"
  ),

  Mesenchymal = c(
    "Col1a1","Col3a1",
    "Pdgfra","Pdgfrb",
    "Lrat","Rbp1","Rgs5"
  ),

  Lymphoid_neutrophil = c(
    "Cd3d","Cd3e","Nkg7",
    "S100a8","S100a9","Retnlg"
  )
)

panels <- lapply(
  panels,
  function(x) intersect(x, rownames(cyc))
)

for (nm in names(panels)) {

  genes <- panels[[nm]]

  if (length(genes) >= 2) {

    cyc <- AddModuleScore(
      cyc,
      features=list(genes),
      name=paste0("MS_", nm, "_"),
      assay="RNA",
      seed=20260902
    )
  }
}

for (nm in names(panels)) {

  old <- paste0("MS_", nm, "_1")

  if (old %in% colnames(cyc@meta.data)) {

    new <- paste0("MS_", nm, "_v6801")

    cyc@meta.data[[new]] <-
      cyc@meta.data[[old]]

    cyc@meta.data[[old]] <- NULL
  }
}

counts <- GetAssayData(
  cyc,
  assay="RNA",
  layer="counts"
)

count_hits <- function(genes) {

  genes <- intersect(
    genes,
    rownames(counts)
  )

  if (length(genes) == 0) {
    return(rep(0, ncol(counts)))
  }

  Matrix::colSums(
    counts[genes, , drop=FALSE] > 0
  )
}

biliary_core_genes <- c(
  "Krt19","Krt7","Epcam","Sox9",
  "Muc1","Hnf1b","Cftr","Slc4a2"
)

cycling_genes <- c(
  "Mki67","Top2a","Birc5",
  "Ube2c","Cenpf","Pcna"
)

hep_genes <- c(
  "Alb","Ttr","Apoa1","Fabp1",
  "Cps1","Ass1","Hnf4a"
)

myeloid_genes <- c(
  "Ptprc","Lyz2","Adgre1",
  "Tyrobp","Aif1"
)

endo_genes <- c(
  "Pecam1","Cdh5","Kdr","Stab2"
)

mes_genes <- c(
  "Col1a1","Col3a1",
  "Pdgfra","Pdgfrb",
  "Lrat","Rbp1","Rgs5"
)

lymph_neut_genes <- c(
  "Cd3d","Cd3e","Nkg7",
  "S100a8","S100a9","Retnlg"
)

cyc$biliary_core_hits_v6801 <-
  count_hits(biliary_core_genes)

cyc$cycling_hits_v6801 <-
  count_hits(cycling_genes)

cyc$hepatocyte_hits_v6801 <-
  count_hits(hep_genes)

cyc$myeloid_hits_v6801 <-
  count_hits(myeloid_genes)

cyc$endothelial_hits_v6801 <-
  count_hits(endo_genes)

cyc$mesenchymal_hits_v6801 <-
  count_hits(mes_genes)

cyc$lymph_neut_hits_v6801 <-
  count_hits(lymph_neut_genes)

cyc$competing_lineage_v6801 <-
  cyc$hepatocyte_hits_v6801 >= 2 |
  cyc$myeloid_hits_v6801 >= 2 |
  cyc$endothelial_hits_v6801 >= 2 |
  cyc$mesenchymal_hits_v6801 >= 2 |
  cyc$lymph_neut_hits_v6801 >= 2

cyc$stringent_biliary_rescue_v6801 <-
  cyc$biliary_core_hits_v6801 >= 2 &
  cyc$cycling_hits_v6801 >= 1 &
  !cyc$competing_lineage_v6801

md <- cyc@meta.data

out <- data.frame(
  cell=rownames(md),
  sample=as.character(md$sample),
  condition=as.character(md$condition),

  biliary_core_hits=
    md$biliary_core_hits_v6801,

  cycling_hits=
    md$cycling_hits_v6801,

  hepatocyte_hits=
    md$hepatocyte_hits_v6801,

  myeloid_hits=
    md$myeloid_hits_v6801,

  endothelial_hits=
    md$endothelial_hits_v6801,

  mesenchymal_hits=
    md$mesenchymal_hits_v6801,

  lymph_neut_hits=
    md$lymph_neut_hits_v6801,

  competing_lineage=
    md$competing_lineage_v6801,

  stringent_biliary_rescue=
    md$stringent_biliary_rescue_v6801,

  stringsAsFactors=FALSE
)

write.csv(
  out,
  file.path(
    TABDIR,
    "Cycling_biliary_stringent_rescue_per_cell_v6.8.0.1.csv"
  ),
  row.names=FALSE
)

sample_summary <- do.call(
  rbind,
  lapply(
    split(out, out$sample),
    function(x) {

      data.frame(
        sample=x$sample[1],
        permissive_candidates=nrow(x),
        stringent_rescue=sum(
          x$stringent_biliary_rescue
        ),
        competing_lineage=sum(
          x$competing_lineage
        ),
        stringent_fraction=
          mean(
            x$stringent_biliary_rescue
          ),
        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(sample_summary) <- NULL

write.csv(
  sample_summary,
  file.path(
    TABDIR,
    "Cycling_biliary_stringent_rescue_by_sample_v6.8.0.1.csv"
  ),
  row.names=FALSE
)

Idents(cyc) <- factor(
  cyc$sample,
  levels=c(
    "STD_rep1",
    "CDHFD_rep1",
    "Sham1",
    "Sham20",
    "Tx17",
    "Tx5"
  )
)

all_features <- unique(
  unlist(
    panels,
    use.names=FALSE
  )
)

p1 <- DotPlot(
  cyc,
  features=all_features,
  assay="RNA",
  dot.scale=6
) +
  RotatedAxis() +
  theme_classic(base_size=9) +
  theme(
    axis.title=element_blank()
  ) +
  ggtitle(
    "Cycling biliary rescue candidates: marker audit by sample"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cycling_biliary_candidates_marker_DotPlot_by_sample_v6.8.0.1.pdf"
  ),
  p1,
  width=21,
  height=6.5
)

cyc$rescue_class_v6801 <- ifelse(
  cyc$stringent_biliary_rescue_v6801,
  "Stringent_biliary_rescue",
  "Rejected_or_ambiguous"
)

Idents(cyc) <- factor(
  cyc$rescue_class_v6801,
  levels=c(
    "Stringent_biliary_rescue",
    "Rejected_or_ambiguous"
  )
)

p2 <- DotPlot(
  cyc,
  features=all_features,
  assay="RNA",
  dot.scale=7
) +
  RotatedAxis() +
  theme_classic(base_size=9) +
  theme(
    axis.title=element_blank()
  ) +
  ggtitle(
    "Cycling biliary rescue: stringent rescued vs rejected"
  )

ggsave(
  file.path(
    FIGDIR,
    "Cycling_biliary_stringent_vs_rejected_DotPlot_v6.8.0.1.pdf"
  ),
  p2,
  width=21,
  height=4.5
)

cat("\n=== PERMISSIVE CANDIDATES ===\n")
cat(ncol(cyc), "\n")

cat("\n=== STRINGENT RESCUE ===\n")
cat(
  sum(
    cyc$stringent_biliary_rescue_v6801
  ),
  "\n"
)

cat("\n=== COMPETING-LINEAGE POSITIVE ===\n")
cat(
  sum(
    cyc$competing_lineage_v6801
  ),
  "\n"
)

cat("\n=== STRINGENT RESCUE BY SAMPLE ===\n")
print(sample_summary)

cat("\n=== BILIARY CORE HIT DISTRIBUTION ===\n")
print(
  table(
    cyc$biliary_core_hits_v6801
  )
)

cat("\n=== CYCLING HIT DISTRIBUTION ===\n")
print(
  table(
    cyc$cycling_hits_v6801
  )
)

saveRDS(
  cyc,
  file.path(
    OUTDIR,
    "Cycling_biliary_rescue_validation_v6.8.0.1.rds"
  ),
  compress=FALSE
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTDIR,
    "sessionInfo_v6.8.0.1.txt"
  )
)

cat("\n====================================================\n")
cat("v6.8.0.1 COMPLETE\n")
cat("No cells added to Cholangiocyte parent yet\n")
cat("Output:", OUTDIR, "\n")
cat("====================================================\n")
