\name{write_metadata_detection_report}
\alias{write_metadata_detection_report}
\title{Write a metadata-detection report}
\usage{
write_metadata_detection_report(metadata_map, path)
}
\arguments{
\item{metadata_map}{Result returned by \code{detect_metadata()}.}
\item{path}{Output CSV path.}
}
\value{The normalized output path, invisibly.}
\description{Writes the detection report from a shared metadata map to CSV.}
