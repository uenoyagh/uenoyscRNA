#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(4840)

INPUT_RDS <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk",
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_dominant_program_v4.8.0",
  "RDS",
  "Mouse_Mphi_Res2_dominant_program_annotated_v4.8.0.rds"
)

OUTPUT_DIR <- file.path(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk",
  "Mouse_MASH_Mphi_RDS",
  "Mphi_Res2_manual_annotation_v4.8.4"
)

FIG_OUT_DIR <- file.path(OUTPUT_DIR, "Figures")
CSV_OUT_DIR <- file.path(OUTPUT_DIR, "Tables")
RDS_OUT_DIR <- file.path(OUTPUT_DIR, "RDS")

dir.create(FIG_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CSV_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

UMAP_REDUCTION <- "mphi.umap.rpca"
RES2_CLUSTER_COL <- "cluster_res2"
SAMPLE_COL <- "sample_4group"
CONDITION_COL <- "condition_4group"

CONDITION_ORDER <- c("STD","CDAHFD","Sham","Tx")

CLASS_ORDER <- c(
  "Inflammatory-Mphi",
  "Anti-inflammatory-Mphi",
  "Fibrogenic-Mphi",
  "Repair/Resolution-Mphi",
  "Lipid-associated/TREM2-Mphi",
  "Other"
)

CONDITION_COLORS <- c(
  "STD"="#2F65FF",
  "CDAHFD"="#F04444",
  "Sham"="#777777",
  "Tx"="#F28C18"
)

CLASS_COLORS <- c(
  "Inflammatory-Mphi"="#E31A1C",
  "Anti-inflammatory-Mphi"="#1478FF",
  "Fibrogenic-Mphi"="#B218B2",
  "Repair/Resolution-Mphi"="#00A65A",
  "Lipid-associated/TREM2-Mphi"="#F28C18",
  "Other"="#B5B5B5"
)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

manual_map <- tibble::tribble(
  ~cluster, ~macrophage_class,
  "0","Anti-inflammatory-Mphi",
  "1","Inflammatory-Mphi",
  "2","Anti-inflammatory-Mphi",
  "3","Fibrogenic-Mphi",
  "4","Other",
  "5","Inflammatory-Mphi",
  "6","Lipid-associated/TREM2-Mphi",
  "7","Other",
  "8","Other",
  "9","Lipid-associated/TREM2-Mphi",
  "10","Inflammatory-Mphi",
  "11","Lipid-associated/TREM2-Mphi",
  "12","Anti-inflammatory-Mphi",
  "13","Lipid-associated/TREM2-Mphi",
  "14","Anti-inflammatory-Mphi",
  "15","Other",
  "16","Lipid-associated/TREM2-Mphi",
  "17","Repair/Resolution-Mphi",
  "18","Lipid-associated/TREM2-Mphi",
  "19","Inflammatory-Mphi",
  "20","Other",
  "21","Fibrogenic-Mphi",
  "22","Lipid-associated/TREM2-Mphi",
  "23","Anti-inflammatory-Mphi",
  "24","Other",
  "25","Inflammatory-Mphi",
  "26","Inflammatory-Mphi",
  "27","Fibrogenic-Mphi"
)

numeric_cluster_levels <- function(x) {
  z <- unique(as.character(x))
  suppressWarnings(n <- as.numeric(z))
  if (all(!is.na(n))) z[order(n)] else sort(z)
}

save_pdf <- function(filename, plot, width, height) {
  ggsave(
    filename=file.path(FIG_OUT_DIR, filename),
    plot=plot,
    device=cairo_pdf,
    width=width,
    height=height,
    units="in",
    limitsize=FALSE
  )
}

save_png_600 <- function(filename, plot, width, height) {
  ggsave(
    filename=file.path(FIG_OUT_DIR, filename),
    plot=plot,
    device="png",
    width=width,
    height=height,
    units="in",
    dpi=600,
    limitsize=FALSE
  )
}

R8TONE_BASE <- c(
  "#F04444","#F28C18","#D8B400","#55A600",
  "#00B85A","#00BFC4","#2F65FF","#B84BE8"
)

make_r8tone_palette <- function(levels_vec) {
  levels_vec <- as.character(levels_vec)
  n <- length(levels_vec)

  tone2 <- c(
    "#D62828","#D96D00","#B49700","#3F8F00",
    "#00994B","#009DA6","#174ED1","#9634C9"
  )
  tone3 <- c(
    "#FF6B6B","#FFA340","#E8C83C","#78C928",
    "#29C978","#28D3D8","#5B82FF","#CD6AF0"
  )
  tone4 <- c(
    "#A91E1E","#B75400","#8C7600","#307000",
    "#00793B","#007A82","#103AA8","#74269F"
  )

  pool <- c(R8TONE_BASE,tone2,tone3,tone4)
  cols <- pool[seq_len(n)]
  names(cols) <- levels_vec
  cols
}

umap_theme <- function(base_size=12) {
  theme_classic(base_size=base_size) +
    theme(
      plot.title=element_text(face="bold"),
      strip.text=element_text(face="bold"),
      panel.border=element_rect(color="black",fill=NA,linewidth=0.6)
    )
}

mphi <- readRDS(INPUT_RDS)

if (!inherits(mphi,"Seurat")) stop("Input object is not a Seurat object.")

required_cols <- c(RES2_CLUSTER_COL,SAMPLE_COL,CONDITION_COL)
missing_cols <- setdiff(required_cols,colnames(mphi@meta.data))
if (length(missing_cols)>0L) {
  stop("Missing metadata column(s): ",paste(missing_cols,collapse=", "))
}

if (!(UMAP_REDUCTION %in% Reductions(mphi))) {
  stop("UMAP reduction not found: ",UMAP_REDUCTION)
}

cluster_levels <- numeric_cluster_levels(mphi@meta.data[[RES2_CLUSTER_COL]])

missing_clusters <- setdiff(cluster_levels,manual_map$cluster)
if (length(missing_clusters)>0L) {
  stop("Manual map is missing cluster(s): ",paste(missing_clusters,collapse=", "))
}

mphi@meta.data[[RES2_CLUSTER_COL]] <- factor(
  as.character(mphi@meta.data[[RES2_CLUSTER_COL]]),
  levels=cluster_levels
)

mphi@meta.data[[CONDITION_COL]] <- factor(
  as.character(mphi@meta.data[[CONDITION_COL]]),
  levels=CONDITION_ORDER
)

class_map <- setNames(manual_map$macrophage_class,manual_map$cluster)

mphi$macrophage_class_Res2_v484 <- factor(
  unname(class_map[as.character(mphi@meta.data[[RES2_CLUSTER_COL]])]),
  levels=CLASS_ORDER
)

if (anyNA(mphi$macrophage_class_Res2_v484)) {
  stop("NA class detected after manual mapping.")
}

write.csv(
  manual_map,
  file.path(CSV_OUT_DIR,"00_Res2_cluster_manual_class_map_v4.8.4.csv"),
  row.names=FALSE
)

saveRDS(
  mphi,
  file.path(RDS_OUT_DIR,"Mouse_Mphi_Res2_manual_class_annotated_v4.8.4.rds"),
  compress=FALSE
)

cluster_palette <- make_r8tone_palette(cluster_levels)

p_cluster_umap <- DimPlot(
  mphi,
  reduction=UMAP_REDUCTION,
  group.by=RES2_CLUSTER_COL,
  label=TRUE,
  repel=TRUE,
  label.size=6.0,
  raster=FALSE,
  pt.size=0.95,
  shuffle=TRUE,
  seed=4840,
  cols=cluster_palette
) +
  labs(
    title="MΦ-only RPCA Res 2.0 clusters",
    subtitle="R8tone | manual annotation v4.8.4",
    color="Res2 cluster"
  ) +
  umap_theme(12)

save_pdf(
  "01_Mphi_Res2_cluster_UMAP_R8tone_v4.8.4.pdf",
  p_cluster_umap,11.5,8.5
)

p_class_umap <- DimPlot(
  mphi,
  reduction=UMAP_REDUCTION,
  group.by="macrophage_class_Res2_v484",
  label=TRUE,
  repel=TRUE,
  label.size=5.0,
  raster=FALSE,
  pt.size=1.00,
  shuffle=TRUE,
  seed=4840,
  cols=CLASS_COLORS
) +
  labs(
    title="MΦ-only RPCA Res 2.0: manual class annotation v4.8.4",
    subtitle="Inflammatory / Anti-inflammatory / Fibrogenic / Repair-Resolution / Lipid-TREM2 / Other",
    color="MΦ class"
  ) +
  umap_theme(12)

save_pdf(
  "02_Mphi_Res2_manual_class_UMAP_v4.8.4.pdf",
  p_class_umap,11,8.5
)

save_png_600(
  "02_Mphi_Res2_manual_class_UMAP_v4.8.4_600dpi.png",
  p_class_umap,11,8.5
)

p_class_split <- DimPlot(
  mphi,
  reduction=UMAP_REDUCTION,
  group.by="macrophage_class_Res2_v484",
  split.by=CONDITION_COL,
  raster=FALSE,
  pt.size=1.10,
  shuffle=TRUE,
  seed=4840,
  cols=CLASS_COLORS,
  ncol=2
) +
  labs(
    title="MΦ manual classes by condition",
    subtitle="STD / CDAHFD / Sham / Tx | shared Res2 RPCA/UMAP coordinates",
    color="MΦ class"
  ) +
  umap_theme(11)

save_pdf(
  "03_Mphi_Res2_manual_class_UMAP_by_condition_v4.8.4.pdf",
  p_class_split,14,11
)

sample_class <- mphi@meta.data %>%
  transmute(
    sample=as.character(.data[[SAMPLE_COL]]),
    condition=factor(
      as.character(.data[[CONDITION_COL]]),
      levels=CONDITION_ORDER
    ),
    macrophage_class=factor(
      as.character(macrophage_class_Res2_v484),
      levels=CLASS_ORDER
    )
  )

sample_condition <- sample_class %>% distinct(sample,condition)

count_raw <- sample_class %>%
  count(sample,condition,macrophage_class,name="n_cells")

sample_abundance <- tidyr::expand_grid(
  sample=unique(sample_class$sample),
  macrophage_class=CLASS_ORDER
) %>%
  left_join(sample_condition,by="sample") %>%
  left_join(
    count_raw,
    by=c("sample","condition","macrophage_class")
  ) %>%
  mutate(n_cells=tidyr::replace_na(n_cells,0L)) %>%
  group_by(sample) %>%
  mutate(
    total_mphi=sum(n_cells),
    fraction=ifelse(total_mphi>0,n_cells/total_mphi,0),
    percent=100*fraction
  ) %>%
  ungroup() %>%
  mutate(
    condition=factor(condition,levels=CONDITION_ORDER),
    macrophage_class=factor(macrophage_class,levels=CLASS_ORDER)
  )

write.csv(
  sample_abundance,
  file.path(CSV_OUT_DIR,"04_Mphi_manual_class_abundance_by_sample_v4.8.4.csv"),
  row.names=FALSE
)

condition_summary <- sample_abundance %>%
  group_by(condition,macrophage_class) %>%
  summarise(
    n_samples=n_distinct(sample),
    mean_cells=mean(n_cells),
    mean_fraction=mean(fraction),
    mean_percent=mean(percent),
    min_percent=min(percent),
    max_percent=max(percent),
    .groups="drop"
  )

write.csv(
  condition_summary,
  file.path(CSV_OUT_DIR,"05_Mphi_manual_class_condition_summary_v4.8.4.csv"),
  row.names=FALSE
)

p_fraction <- ggplot(
  condition_summary,
  aes(x=condition,y=mean_percent,fill=condition)
) +
  geom_col(width=0.68,alpha=0.90) +
  geom_point(
    data=sample_abundance,
    aes(x=condition,y=percent),
    inherit.aes=FALSE,
    shape=21,
    size=2.7,
    fill="white",
    stroke=0.75,
    position=position_jitter(width=0.07,height=0)
  ) +
  facet_wrap(~macrophage_class,scales="free_y",ncol=3) +
  scale_fill_manual(values=CONDITION_COLORS) +
  labs(
    title="MΦ manual classes across four conditions",
    subtitle="Bars = condition mean; white points = biological samples",
    x=NULL,
    y="Fraction of total MΦ (%)",
    fill=NULL
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title=element_text(face="bold"),
    strip.text=element_text(face="bold"),
    axis.text.x=element_text(angle=25,hjust=1),
    legend.position="top"
  )

save_pdf(
  "06_Mphi_manual_class_fraction_4conditions_v4.8.4.pdf",
  p_fraction,14,9
)

p_composition <- ggplot(
  condition_summary,
  aes(x=condition,y=mean_fraction,fill=macrophage_class)
) +
  geom_col(width=0.76) +
  scale_fill_manual(values=CLASS_COLORS,drop=FALSE) +
  scale_y_continuous(
    labels=scales::percent_format(accuracy=1),
    expand=c(0,0)
  ) +
  labs(
    title="Mean MΦ composition: manual annotation v4.8.4",
    subtitle="Mean biological-sample fractions within each condition",
    x=NULL,
    y="Mean fraction of MΦ",
    fill="MΦ class"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title=element_text(face="bold"),
    axis.text.x=element_text(angle=25,hjust=1),
    legend.position="right"
  )

save_pdf(
  "07_Mphi_manual_class_mean_composition_v4.8.4.pdf",
  p_composition,11.5,7
)

change_table <- condition_summary %>%
  select(condition,macrophage_class,mean_fraction,mean_percent) %>%
  pivot_wider(
    names_from=condition,
    values_from=c(mean_fraction,mean_percent)
  ) %>%
  mutate(
    delta_Disease_percent=mean_percent_CDAHFD-mean_percent_STD,
    delta_Treatment_percent=mean_percent_Tx-mean_percent_Sham
  )

write.csv(
  change_table,
  file.path(CSV_OUT_DIR,"08_Mphi_manual_class_Disease_Treatment_change_v4.8.4.csv"),
  row.names=FALSE
)

p_cluster_class <- manual_map %>%
  mutate(
    cluster=factor(cluster,levels=cluster_levels),
    macrophage_class=factor(macrophage_class,levels=CLASS_ORDER)
  ) %>%
  ggplot(
    aes(x=cluster,y=1,fill=macrophage_class)
  ) +
  geom_tile(color="white",linewidth=0.50) +
  geom_text(aes(label=cluster),size=3.8,fontface="bold") +
  scale_fill_manual(values=CLASS_COLORS,drop=FALSE) +
  scale_y_continuous(breaks=NULL) +
  labs(
    title="Res2 cluster -> manual macrophage class v4.8.4",
    x="Res2 cluster",
    y=NULL,
    fill="MΦ class"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title=element_text(face="bold"),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    legend.position="bottom"
  )

save_pdf(
  "09_Res2_cluster_to_manual_class_map_v4.8.4.pdf",
  p_cluster_class,14,4.5
)

p_summary <- (
  p_class_umap /
  p_class_split /
  p_fraction /
  p_composition /
  p_cluster_class
) +
  patchwork::plot_layout(
    heights=c(1.1,1.4,1.1,0.9,0.55)
  ) +
  patchwork::plot_annotation(
    title="Mouse MASH MΦ Res2 manual annotation v4.8.4",
    subtitle="Manual cluster classification fixed for downstream analysis",
    theme=theme(
      plot.title=element_text(face="bold",size=20),
      plot.subtitle=element_text(size=12)
    )
  )

save_pdf(
  "10_Mphi_Res2_manual_annotation_summary_v4.8.4.pdf",
  p_summary,17,27
)

summary_lines <- c(
  "Mouse MASH MΦ Res2 manual annotation v4.8.4",
  "",
  "Important:",
  "  RPCA not recalculated.",
  "  UMAP not recalculated.",
  "  Res2 clustering not recalculated.",
  "  Only cluster-to-class annotation changed.",
  "",
  "First inspect:",
  "  00 cluster-class map CSV",
  "  01 cluster UMAP",
  "  02 manual class UMAP",
  "  03 condition-split UMAP",
  "  06 class fractions",
  "  07 mean composition",
  "  08 disease/treatment changes",
  "  09 cluster-class map"
)

writeLines(
  summary_lines,
  file.path(
    OUTPUT_DIR,
    "README_Mphi_Res2_manual_annotation_v4.8.4.txt"
  )
)

capture.output(
  sessionInfo(),
  file=file.path(
    OUTPUT_DIR,
    "sessionInfo_v4.8.4.txt"
  )
)

message("DONE: ",OUTPUT_DIR)
