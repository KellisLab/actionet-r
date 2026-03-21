#' @export
normalize.ace <- function(
    adata = NULL,
    layer = "counts",
    layer_out = "logcounts",
    scale_param = stats::median,
    trans_func = base::log2,
    pseudocount = 1,
    assay_name = NULL,
    assay_out = NULL,
    ace = NULL) {
  adata <- .resolve_container_arg(adata = adata, ace = ace)
  layer <- .resolve_layer_arg(
    layer = layer,
    assay_name = assay_name,
    default = "counts",
    layer_missing = missing(layer),
    assay_name_missing = missing(assay_name)
  )
  if (!missing(assay_out) && !is.null(assay_out)) {
    .Deprecated(msg = "'assay_out' is deprecated; use 'layer_out' instead.")
    if (missing(layer_out) || is.null(layer_out)) {
      layer_out <- assay_out
    }
  }
  adata <- .validate_ace(adata, allow_se_like = TRUE, as_ace = TRUE, return_elem = TRUE)
  S <- .validate_assay(adata, assay_name = layer, return_elem = TRUE)
  S <- normalize.matrix(
    S,
    dim = 2,
    scale_param = scale_param,
    trans_func = trans_func,
    pseudocount = pseudocount)
  rownames(S) <- .actionet_rownames(adata)
  colnames(S) <- .actionet_colnames(adata)
  adata <- .set_layer_matrix(adata, layer_out, S)
  if (is.null(adata$X)) {
    adata <- .set_layer_matrix(adata, NULL, S)
  }
  return(adata)
}

#' @importFrom batchelor multiBatchNorm
#' @export
normalize.multiBatchNorm <- function(adata = NULL,
                                     batch_attr,
                                     layer = "counts",
                                     layer_out = "logcounts",
                                     BPPARAM = SerialParam(),
                                     norm.args = list(),
                                     min.mean = 1,
                                     subset.row = NULL,
                                     normalize.all = FALSE,
                                     preserve.single = TRUE,
                                     assay_name = NULL,
                                     assay_out = NULL,
                                     ace = NULL) {
  adata <- .resolve_container_arg(adata = adata, ace = ace)
  layer <- .resolve_layer_arg(
    layer = layer,
    assay_name = assay_name,
    default = "counts",
    layer_missing = missing(layer),
    assay_name_missing = missing(assay_name)
  )
  if (!missing(assay_out) && !is.null(assay_out)) {
    .Deprecated(msg = "'assay_out' is deprecated; use 'layer_out' instead.")
    if (missing(layer_out) || is.null(layer_out)) {
      layer_out <- assay_out
    }
  }
  adata <- .validate_ace(adata, allow_se_like = TRUE, as_ace = TRUE, return_elem = TRUE)
  batch_attr <- .validate_vector_attr(adata, attr = batch_attr, return_type = "data")
  sce_temp <- adata$as_SingleCellExperiment()

  sce_temp <- batchelor::multiBatchNorm(
    sce_temp,
    batch = batch_attr,
    norm.args = norm.args,
    assay.type = layer,
    min.mean = min.mean,
    subset.row = subset.row,
    normalize.all = normalize.all,
    preserve.single = preserve.single,
    BPPARAM = BPPARAM
  )

  adata <- .set_layer_matrix(adata, layer_out, SummarizedExperiment::assays(sce_temp)[["logcounts"]])
  uns <- .get_uns(adata)
  uns$norm_method <- "multiBatchNorm"
  adata <- .set_uns(adata, uns)

  return(adata)
}


#' @export
normalize.matrix <- function(S,
                             dim = 2,
                             scale_param = NULL,
                             trans_func = NULL,
                             pseudocount = 0) {
  if (!is.matrix(S) && !.is_sparse_matrix(S)) {
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
    if (.is_sparse_matrix(S)) {
      S@x <- trans_func(S@x + pseudocount)
    } else {
      S[S != 0] <- trans_func(S[S != 0] + pseudocount)
    }
  }

  return(S)
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
