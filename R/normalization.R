#' @export
normalize.ace <- function(
    ace,
    norm_method = c("default", "multiBatchNorm"),
    batch_attr = NULL,
    assay_name = "counts",
    assay_out = "logcounts",
    scale_param = stats::median,
    trans_func = base::log1p,
    BPPARAM = SerialParam()) {
  norm_method <- match.arg(norm_method)

  if (is.null(norm_method)) {

  } else if (norm_method == "multiBatchNorm") {
    ace <- normalize.multiBatchNorm(
      ace = ace,
      batch_attr = batch_attr,
      assay_name = assay_name,
      assay_out = assay_out,
      BPPARAM = BPPARAM
    )
  } else {
    ace <- .normalize.default(
      ace = ace,
      assay_name = assay_name,
      assay_out = assay_out,
      scale_param = scale_param,
      trans_func = trans_func
    )
    norm_method <- "default"
  }

  metadata(ace)$norm_method <- norm_method

  return(ace)
}

#' @importFrom batchelor multiBatchNorm
#' @export
normalize.multiBatchNorm <- function(ace,
                                     batch_attr,
                                     assay_name = "counts",
                                     assay_out = "logcounts",
                                     BPPARAM = SerialParam(),
                                     norm.args = list(),
                                     min.mean = 1,
                                     subset.row = NULL,
                                     normalize.all = FALSE,
                                     preserve.single = TRUE) {
  batch_attr <- ACTIONetExperiment::get.data.or.split(ace, attr = batch_attr, to_return = "data")
  sce_temp <- as(ace, "SingleCellExperiment")

  sce_temp <- batchelor::multiBatchNorm(
    sce_temp,
    batch = batch_attr,
    norm.args = norm.args,
    assay.type = assay_name,
    min.mean = min.mean,
    subset.row = subset.row,
    normalize.all = normalize.all,
    preserve.single = preserve.single,
    BPPARAM = BPPARAM
  )

  SummarizedExperiment::assays(ace)[[assay_out]] <- SummarizedExperiment::assays(sce_temp)[["logcounts"]]

  return(ace)
}


#' @export
normalize.matrix <- function(S,
                             dim = 2,
                             scale_param = NULL,
                             trans_func = NULL) {
  if (!is.matrix(S) && !ACTIONetExperiment:::is.sparseMatrix(S)) {
    err <- sprintf("`S` must be `matrix` or `sparseMatrix`.\n")
    stop(err)
  }

  if (!dim %in% c(1, 2)) {
    err <- sprintf("Invalid `dim`.\n")
    stop(err)
  }

  if (!is.null(scale_param)) {
    if (!is.function(scale_param) && !is.numeric(scale_param)) {
      err <- sprintf("`scale_param` must be `function` or `numeric`.\n")
      stop(err)
    } else if (!(length(scale_param) == 1) ||
      !(length(scale_param) == dim(S)[dim])) {
      err <- sprintf("`scale_param` must be of length 1 or dim(S) (%d).\n", dim(S)[dim])
    }
  }

  if (!is.null(trans_func) && !is.function(trans_func)) {
    err <- sprintf("`trans_func` must be an element-wise function.\n")
    stop(err)
  }

  # Library sizes are ABSOLUTE sums!!! See p-norm formula.
  lib_sizes <- switch(dim,
    "1" = Matrix::rowSums(abs(S)),
    "2" = Matrix::colSums(abs(S))
  )
  norm_factor <- lib_sizes # Require untransformed lib_sizes later
  norm_factor[norm_factor == 0] <- 1
  norm_factor <- 1.0 / norm_factor

  S <- .scale.matrix(X = S, dim = dim, scale_fac = norm_factor)

  if (is.function(scale_param)) {
    scale_param <- scale_param(lib_sizes)
  }
  S <- .scale.matrix(S, dim = dim, scale_fac = scale_param)

  if (!is.null(trans_func)) {
    S <- trans_func(S)
  }

  return(S)
}


.normalize.default <- function(
    ace,
    assay_name = "counts",
    assay_out = "logcounts",
    scale_param = NULL,
    trans_func = NULL) {
  S <- SummarizedExperiment::assays(ace)[[assay_name]]
  S <- normalize.matrix(S, dim = 2, scale_param = scale_param, trans_func = trans_func)
  rownames(S) <- rownames(ace)
  colnames(S) <- colnames(ace)
  SummarizedExperiment::assays(ace)[[assay_out]] <- S
  return(ace)
}


.scale.matrix <- function(X, dim, scale_fac = NULL) {
  dn <- dimnames(X)
  if (is.null(scale_fac)) {
    return(X)
  }

  if (length(scale_fac) == 1) {
    X <- X * scale_fac
  } else {
    if (is.matrix(X)) {
      X <- switch(dim,
        "1" = C_scaleMatrixDense(X, scale_fac, 1), # Better memory.
        "2" = C_scaleMatrixDense(X, scale_fac, 0) # Absolutely better.
      )
    } else {
      if (!is(X, "CsparseMatrix")) {
        X <- as(X, "CsparseMatrix")
      }
      X <- switch(dim,
        "1" = C_scaleMatrixSparse(X, scale_fac, 1), # Better memory.
        "2" = X %*% Matrix::Diagonal(n = length(scale_fac), x = scale_fac)
      )
    }
  }
  dimnames(X) <- dn
  return(X)
}
