# uenoyscRNA Annotation Validation v3.2
# User-editable configuration
# Framework/scripts: internal SSD
# RDS and outputs: external SSD

CFG <- list(
  framework_version = "uenoyscRNA Framework v1.0",
  pipeline_version = "RDS3 annotation validation v3.2",

  # Internal SSD project root
  project_root = "/Users/uenoya/Projects/uenoyscRNA",

  # External SSD candidates are checked in this order.
  external_ssd_candidates = c(
    "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk",
    "/Volumes/SSD990_uenoy"
  ),

  rds_relative_path = file.path(
    "Mouse_MASH_RDS",
    "Mouse_object_with_FIXED2_R8tone_sample_celltype_metadata.rds"
  ),

  output_relative_dir = file.path(
    "Mouse_MASH_RDS",
    "RDS3_annotation_validation_v3.2"
  ),

  assay = "RNA",
  data_layer = "data",
  counts_layer = "counts",
  reduction = "umapRPCA",

  # Metadata columns. NULL enables automatic detection.
  cluster_col = "cluster_for_R8plot_FIXED2",
  annotation_col = "celltype_for_R8plot_FIXED2",
  sample_col = "sample_display_FIXED2",
  condition_col = "condition_FIXED2",

  cluster_candidates = c(
    "seurat_clusters", "cluster", "clusters",
    "RNA_snn_res.3", "RNA_snn_res.2.5",
    "integrated_snn_res.3", "integrated_snn_res.2.5"
  ),
  annotation_candidates = c(
    "celltype", "cell_type", "annotation", "celltype_annotation",
    "layer1", "Layer1", "major_celltype", "CellType"
  ),
  sample_candidates = c(
    "sample", "sample_id", "orig.ident", "Sample", "sample_name"
  ),
  condition_candidates = c(
    "condition", "group", "treatment", "Condition", "Group"
  ),

  # Parallel settings
  use_future = TRUE,
  workers = max(1L, min(6L, parallel::detectCores(logical = TRUE) - 1L)),
  future_max_size_gb = 32,

  # Marker detection
  marker_engine = "presto",
  min_cells_per_cluster = 20L,
  min_pct = 0.10,
  logfc_threshold = 0.10,
  only_pos = FALSE,
  top_n_markers = 100L,
  top_n_plot = 10L,

  # Marker voting
  score_method = "mean_expression",
  min_present_genes = 2L,
  positive_weight = 1.0,
  negative_weight = 1.0,
  ambiguity_delta = 0.15,
  low_confidence_score = 0.05,

  # Plot settings
  umap_point_size = 0.35,
  umap_label = TRUE,
  dotplot_scale_min = -2.5,
  dotplot_scale_max = 2.5,
  violin_ncol = 4L,
  pdf_width = 14,
  pdf_height = 10,

  # Resume/checkpoint behavior
  resume = TRUE,
  overwrite = FALSE,
  save_intermediate_rds = TRUE,
  seed = 1234L
)
