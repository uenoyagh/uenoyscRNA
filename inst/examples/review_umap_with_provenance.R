# Replace with the actual RDS path before running this example.
rds_path <- "/actual/path/to/your_object.rds"
seurat_object <- readRDS(rds_path)

result <- review_umap(
  object = seurat_object,
  group_by = "celltype_for_R8plot_FIXED2",
  split_by = "sample_display_FIXED2",
  rds_file = rds_path,
  analysis_date = Sys.Date(),
  output_dir = "results",
  filename = "review_umap.pdf"
)
