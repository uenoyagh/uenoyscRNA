\name{detect_metadata}
\alias{detect_metadata}
\title{Detect shared metadata roles in a Seurat object}
\usage{
detect_metadata(
  object,
  overrides = list(),
  required_roles = c("celltype", "condition"),
  strict = TRUE
)
}
\arguments{
\item{object}{A Seurat object.}
\item{overrides}{Named list containing explicit columns for any of \code{celltype}, \code{cluster}, \code{condition}, and \code{sample}.}
\item{required_roles}{Metadata roles that must be detected.}
\item{strict}{Logical; stop when a required role cannot be detected.}
}
\value{An object of class \code{uenoy_metadata_map} containing selected columns and a detection report.}
\description{Identifies likely cell-type, cluster, condition, and sample columns. Explicit overrides have highest priority.}
