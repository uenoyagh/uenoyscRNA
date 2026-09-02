suppressPackageStartupMessages({
  library(utils)
})

VERSION <- "v6.7.2.1"

INDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_v6.7.2/tables"
)

OUTDIR <- paste0(
  "/Volumes/SSD990_uenoy/scRNA_MASLD_MASH_pseudobulk/",
  "Mouse_MASH_RDS/Mouse_MASH_LSEC_", VERSION
)

TABDIR <- file.path(OUTDIR, "tables")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

INFILE <- file.path(
  INDIR,
  "FindAllMarkers_res0.6_all_positive_v6.7.2.csv"
)

cat("============================================\n")
cat("Mouse MASH LSEC marker-table correction\n")
cat("Version:", VERSION, "\n")
cat("============================================\n\n")

if (!file.exists(INFILE)) {
  stop("Missing input file: ", INFILE)
}

x <- read.csv(
  INFILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c(
  "gene",
  "cluster",
  "p_val_adj"
)

missing <- setdiff(
  required,
  colnames(x)
)

if (length(missing) > 0) {
  stop(
    "Missing column(s): ",
    paste(missing, collapse = ", ")
  )
}

fc_candidates <- c(
  "avg_log2FC",
  "avg_logFC"
)

fc_col <- intersect(
  fc_candidates,
  colnames(x)
)[1]

if (is.na(fc_col)) {
  stop("No log-fold-change column found.")
}

sig <- x[
  is.finite(x$p_val_adj) &
    x$p_val_adj < 0.05,
  ,
  drop = FALSE
]

sig <- sig[
  order(
    as.numeric(as.character(sig$cluster)),
    -sig[[fc_col]],
    sig$p_val_adj
  ),
  ,
  drop = FALSE
]

write.csv(
  sig,
  file.path(
    TABDIR,
    "FindAllMarkers_res0.6_SIGNIFICANT_v6.7.2.1.csv"
  ),
  row.names = FALSE
)

top20 <- do.call(
  rbind,
  lapply(
    split(sig, sig$cluster),
    function(y) {

      y <- y[
        order(
          -y[[fc_col]],
          y$p_val_adj
        ),
        ,
        drop = FALSE
      ]

      head(y, 20)
    }
  )
)

rownames(top20) <- NULL

write.csv(
  top20,
  file.path(
    TABDIR,
    "FindAllMarkers_res0.6_SIGNIFICANT_TOP20_per_cluster_v6.7.2.1.csv"
  ),
  row.names = FALSE
)

cat("Significant markers:", nrow(sig), "\n")

cat("\n=== TOP 10 SIGNIFICANT MARKERS ===\n")

clusters <- sort(
  unique(
    as.numeric(
      as.character(top20$cluster)
    )
  )
)

for (cl in clusters) {

  y <- top20[
    as.character(top20$cluster) ==
      as.character(cl),
    ,
    drop = FALSE
  ]

  cat(
    "Cluster ", cl, ": ",
    paste(
      head(y$gene, 10),
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

cat("\n============================================\n")
cat("v6.7.2.1 COMPLETE\n")
cat("============================================\n")
