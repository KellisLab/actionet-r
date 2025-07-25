% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputation.R
\name{imputeFeatures}
\alias{imputeFeatures}
\title{Gene expression imputation using network diffusion.}
\usage{
imputeFeatures(
  ace,
  features,
  algorithm = c("actionet", "pca"),
  features_use = NULL,
  alpha = 0.85,
  norm_method = "pagerank_sym",
  thread_no = 0,
  max_it = 5,
  assay_name = "logcounts",
  reduction_slot = "action",
  net_slot = "actionet"
)
}
\arguments{
\item{ace}{ACTIONetExperiment object containing output of 'run.ACTIONet()'.}

\item{features_use}{A vector of features of length NROW(ace) or the name of a column of rowData(ace) containing the genes given in 'genes'.}

\item{alpha}{Depth of diffusion between (0, 1).
The larger it is, the deeper the diffusion, which results in less nonzeros (default = 0.85).}

\item{thread_no}{Number of parallel threads}

\item{assay_name}{Slot in the ace object with normalized counts.}

\item{genes}{The list of genes to perform imputation on.}

\item{diffusion_it}{Number of diffusion iterations (default = 5)}
}
\value{
Imputed gene expression matrix. Column names are set with imputed genes names and rows are cells.
}
\description{
Gene expression imputation using network diffusion.
}
\examples{
imputed.genes <- impute.genes.using.ACTIONet(ace, c("CD14", "CD19", "CD3G"))
plot.ACTIONet.gradient(ace, imputed.genes[, 1])
}
