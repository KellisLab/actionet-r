
.normalize.default <- function(
    ace,
    assay_name = "counts",
    assay_out = "logcounts",
    scale_param = NULL,
    trans_func = NULL
) {
  S <- SummarizedExperiment::assays(ace)[[assay_name]]
  S <- normalize.matrix(S, dim = 2, scale_param = scale_param, trans_func = trans_func)
  rownames(S) <- rownames(ace)
  colnames(S) <- colnames(ace)
  SummarizedExperiment::assays(ace)[[assay_out]] <- S
  return(ace)
}

.scale.matrix <- function(X, dim, scale_fac = NULL) {
  if(is.null(scale_fac)) {
    return(X)
  }

  if (length(scale_fac) == 1) {
    X <- X * scale_fac
  } else {
    if (is.matrix(X)) {
      X = switch (dim,
                  "1" = scaleMatrixDense(X, scale_fac, 1), # Better memory.
                  "2" = scaleMatrixDense(X, scale_fac, 0)) # Absolutely better.
    } else {
      if (!is(X, "CsparseMatrix")) {
        X = as(X, "CsparseMatrix")
      }
      X = switch (
        dim,
        "1" = scaleMatrixSparse(X, scale_fac, 1), # Better memory.
        "2" = X %*% Matrix::Diagonal(n = length(scale_fac), x = scale_fac)
      )
    }
  }
  return(X)
}
