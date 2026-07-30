############################################################
## 02_validate_RDS3_markers_v2.R
## RDS3 全クラスタ Marker抽出（presto対応・再開可能）
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(openxlsx)
  library(presto)
})

options(future.globals.maxSize = 32 * 1024^3)

############################################################
## 設定
############################################################

rds_file <-
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/results/Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"

outdir <- "Mouse_MASH_RDS3_validation/Phase2_Markers_v2"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "PerCluster"), showWarnings = FALSE)
dir.create(file.path(outdir, "Logs"), showWarnings = FALSE)

############################################################
## Load
############################################################

cat("Loading RDS...\n")
obj <- readRDS(rds_file)

DefaultAssay(obj) <- "RNA"

if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj <- JoinLayers(obj)
}

############################################################
## Cluster column
############################################################

cluster_candidates <- c(
  "cluster_for_R8plot_FIXED2",
  "cluster_for_R8plot",
  "seurat_clusters"
)

cluster_col <- cluster_candidates[
  cluster_candidates %in% colnames(obj@meta.data)
][1]

if (is.na(cluster_col))
  stop("Cluster column not found.")

Idents(obj) <- obj[[cluster_col]][,1]

cat("Cluster column :", cluster_col, "\n")

############################################################
## Resume
############################################################

clusters <- levels(Idents(obj))

completed <- list.files(
  file.path(outdir,"PerCluster"),
  pattern="Cluster_.*csv"
)

completed <- sub("Cluster_","",completed)
completed <- sub(".csv","",completed,fixed=TRUE)

todo <- setdiff(clusters, completed)

cat(length(todo),"clusters remaining\n")

############################################################
## Marker extraction
############################################################

for(cl in todo){

  cat("----------------------------------------\n")
  cat("Cluster",cl,"\n")

  mk <- FindMarkers(
    object=obj,
    ident.1=cl,
    only.pos=TRUE,
    test.use="wilcox",
    logfc.threshold=0.25,
    min.pct=0.10
  )

  mk$gene <- rownames(mk)
  mk$cluster <- cl

  fwrite(
    mk,
    file.path(
      outdir,
      "PerCluster",
      paste0("Cluster_",cl,".csv")
    )
  )

  write(
    paste(Sys.time(),"Finished cluster",cl),
    file=file.path(outdir,"Logs","FindMarkers_progress.log"),
    append=TRUE
  )

}

############################################################
## Merge
############################################################

cat("Merging...\n")

files <- list.files(
  file.path(outdir,"PerCluster"),
  full.names=TRUE,
  pattern="csv$"
)

allmarkers <- rbindlist(
  lapply(files,fread),
  fill=TRUE
)

fwrite(
  allmarkers,
  file.path(outdir,"AllClusterMarkers.csv")
)

############################################################
## Top10/20/30
############################################################

top10 <-
  allmarkers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC,n=10)

top20 <-
  allmarkers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC,n=20)

top30 <-
  allmarkers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC,n=30)

fwrite(top10,file.path(outdir,"Top10Markers.csv"))
fwrite(top20,file.path(outdir,"Top20Markers.csv"))
fwrite(top30,file.path(outdir,"Top30Markers.csv"))

############################################################
## Excel
############################################################

wb <- createWorkbook()

addWorksheet(wb,"AllMarkers")
writeData(wb,"AllMarkers",allmarkers)

addWorksheet(wb,"Top10")
writeData(wb,"Top10",top10)

addWorksheet(wb,"Top20")
writeData(wb,"Top20",top20)

addWorksheet(wb,"Top30")
writeData(wb,"Top30",top30)

saveWorkbook(
  wb,
  file.path(outdir,"RDS3_ClusterMarkers_v2.xlsx"),
  overwrite=TRUE
)

############################################################

cat("\nDone.\n")
