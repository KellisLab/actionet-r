% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{annotate.archetypes.using.labels}
\alias{annotate.archetypes.using.labels}
\title{Annotate archetypes using prior cell annotations
(It uses t-test on the archetype footprint matrix (H))}
\usage{
annotate.archetypes.using.labels(
  ace,
  labels,
  archetype.slot = "H_merged",
  algorithm = "ttest"
)
}
\arguments{
\item{ace}{ACTIONet output object}

\item{labels}{Annotation of interest (clusters, celltypes, etc.) to test enrichment}
}
\value{
A named list: \itemize{
\item Label: Inferred archetype labels
\item Confidence: Confidence of inferred labels
\item Enrichment: Full enrichment matrix
}
}
\description{
Annotate archetypes using prior cell annotations
(It uses t-test on the archetype footprint matrix (H))
}
\examples{
arch.annot <- annotate.archetypes.using.labels(ace, sce$celltypes)
}
