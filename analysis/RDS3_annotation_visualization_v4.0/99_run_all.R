#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath("analysis/RDS3_annotation_visualization_v4.0/99_run_all.R", mustWork = FALSE)
}
analysis_dir <- dirname(script_path)
project_root_guess <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = FALSE)

source(file.path(project_root_guess, "config", "config.R"))
source(file.path(project_root_guess, "R", "01_utils_v3.2.R"))
source(file.path(project_root_guess, "R", "04_plot_export_v3.2.R"))
source(file.path(project_root_guess, "R", "04_plot_export_v4.0.R"))

required <- c(
  "Seurat", "SeuratObject", "ggplot2", "ggrepel",
  "patchwork", "openxlsx"
)
assert_packages(required)

external_root <- detect_external_root(CFG$external_ssd_candidates)
base_dir <- file.path(
  external_root,
  dirname(CFG$output_relative_dir)
)

v34_dir <- file.path(base_dir, "RDS3_annotation_validation_v3.4")
rds_path <- file.path(
  v34_dir,
  "objects",
  "RDS3_with_annotation_validation_v3.4.rds"
)
audit_path <- file.path(
  v34_dir,
  "tables",
  "annotation_audit_v3.4.csv"
)

if (!file.exists(rds_path)) stop("v3.4 RDS not found: ", rds_path)
if (!file.exists(audit_path)) stop("v3.4 audit table not found: ", audit_path)

output_dir <- file.path(base_dir, "RDS3_annotation_visualization_v4.0")
fig_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
object_dir <- file.path(output_dir, "objects")
ensure_dirs(c(output_dir, fig_dir, table_dir, object_dir))

object <- readRDS(rds_path)
audit <- utils::read.csv(audit_path, check.names = FALSE)

reduction <- CFG$reduction
cluster_col <- CFG$cluster_col
if (!cluster_col %in% colnames(object[[]])) {
  cluster_col <- detect_metadata_column(
    object[[]], CFG$cluster_col, CFG$cluster_candidates, "cluster"
  )
}

object <- derive_summary_annotation_v40(object)
object <- derive_annotation_difference_v40(
  object,
  current_col = CFG$annotation_col,
  recommended_col = "vote_ueno_summary_v40"
)

# Interpretation UMAPs
plots <- list(
  lineage = publish_umap_r8_v40(
    object, "vote_ueno_lineage_v34", reduction,
    ueno_lineage_palette_v40,
    "RDS3: Ueno lineage",
    pt_size = 0.24, label_size = 3.7
  ),
  celltype = publish_umap_r8_v40(
    object, "vote_ueno_celltype_v34", reduction,
    ueno_celltype_palette_v40,
    "RDS3: Ueno cell type",
    pt_size = 0.24, label_size = 3.3
  ),
  subtype = publish_umap_r8_v40(
    object, "vote_ueno_subtype_v34", reduction,
    ueno_subtype_palette_v40,
    "RDS3: Ueno subtype",
    pt_size = 0.24, label_size = 3.0
  ),
  summary = publish_umap_r8_v40(
    object, "vote_ueno_summary_v40", reduction,
    ueno_subtype_palette_v40,
    "RDS3: comprehensive liver-cell annotation",
    pt_size = 0.26, label_size = 3.2, legend_ncol = 1
  ),
  difference = publish_umap_r8_v40(
    object, "annotation_difference_v40", reduction,
    annotation_difference_palette_v40,
    "RDS3: current vs recommended annotation",
    pt_size = 0.24, label_size = 3.8
  ),
  confidence_lineage = publish_confidence_umap_v40(
    object, "vote_ueno_lineage_v34_confidence", reduction,
    "RDS3: lineage confidence"
  ),
  confidence_celltype = publish_confidence_umap_v40(
    object, "vote_ueno_celltype_v34_confidence", reduction,
    "RDS3: cell-type confidence"
  ),
  confidence_subtype = publish_confidence_umap_v40(
    object, "vote_ueno_subtype_v34_confidence", reduction,
    "RDS3: subtype confidence"
  )
)

for (nm in names(plots)) {
  save_pdf(
    plots[[nm]],
    file.path(fig_dir, paste0("UMAP_", nm, "_R8_square.pdf")),
    width = 11,
    height = 11
  )
}

# A wider explanatory version with larger labels and two-column legend.
p_summary_wide <- publish_umap_r8_v40(
  object,
  "vote_ueno_summary_v40",
  reduction,
  ueno_subtype_palette_v40,
  "RDS3: comprehensive liver-cell annotation",
  pt_size = 0.28,
  label_size = 3.6,
  legend_ncol = 2
)
save_pdf(
  p_summary_wide,
  file.path(fig_dir, "UMAP_summary_R8_explanatory_wide.pdf"),
  width = 15,
  height = 11
)

# Cluster audit PDF
pages <- make_cluster_audit_pages_v40(
  object = object,
  audit_table = audit,
  cluster_col = cluster_col,
  reduction = reduction,
  current_col = CFG$annotation_col
)

audit_pdf <- file.path(fig_dir, "Cluster_annotation_audit_v4.0.pdf")
grDevices::cairo_pdf(audit_pdf, width = 14, height = 8.5, onefile = TRUE)
for (p in pages) print(p)
grDevices::dev.off()

# Summary tables
meta <- object[[]]
summary_table <- as.data.frame(table(
  summary_annotation = meta$vote_ueno_summary_v40,
  useNA = "ifany"
))
summary_table <- summary_table[order(-summary_table$Freq), , drop = FALSE]

difference_table <- as.data.frame(table(
  status = meta$annotation_difference_v40,
  useNA = "ifany"
))

safe_write_csv(
  summary_table,
  file.path(table_dir, "summary_annotation_cell_counts_v4.0.csv")
)
safe_write_csv(
  difference_table,
  file.path(table_dir, "annotation_difference_counts_v4.0.csv")
)

saveRDS(
  object,
  file.path(object_dir, "RDS3_with_visualization_metadata_v4.0.rds"),
  compress = FALSE
)

write_excel_report(
  file.path(table_dir, "RDS3_annotation_visualization_v4.0.xlsx"),
  list(
    Summary_counts = summary_table,
    Difference_counts = difference_table,
    Annotation_audit = audit
  )
)

message("v4.0 visualization outputs completed: ", output_dir)
