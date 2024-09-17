#' Takes a `ACTIONetExperiment` object and adds the reduced kernel matrix
#'
#' @param ace ACTIONetExperiment object.
#' @param reduced_dim Dimension of reduced kernel matrix. Number of singular vectors to estimate. Passed to <code>runSVD()</code>.
#' @param assay_name Name of assay to reduce.
#' @param reduction_out Name of slot to store reduced output.
#' @param algorithm Singular value decomposition algorithm. Passed to <code>runSVD()</code>.
#' @param max_it Number of SVD iterations. If `NULL`: 1000 for "ilrb", 5 otherwise.
#' @param seed Random seed.
#'
#' @return ACTIONetExperiment object with reduction in colMaps(ace).
#'
#' @examples
#' ace <- import.ace.from.10X(input_path)
#' ace <- reduce.ace(ace)
#' @export
reduce.ace <- function(ace,
                       reduced_dim = 50,
                       assay_name = "logcounts",
                       reduction_out = "action",
                       algorithm = c("irlb", "halko", "feng"),
                       max_it = NULL,
                       seed = 0, verbose = TRUE) {

  ace = .validate_ace(ace, as_ace = TRUE, allow_se_like = TRUE, return_elem = TRUE)
  if (is.null(rownames(ace))) {
    rownames(ace) <- ACTIONetExperiment:::.default_rownames(NROW(ace))
  } else {
    rownames(ace) <- make.unique(rownames(ace), sep = "_")
  }

  if (is.null(colnames(ace))) {
    colnames(ace) <- ACTIONetExperiment:::.default_colnames(NCOL(ace))
  } else {
    colnames(ace) <- make.unique(colnames(ace), sep = "_")
  }

  if (is.null(assay_name)) {
    if ("default_assay" %in% names(metadata(ace))) {
      message(sprintf("Input assay_name is NULL. Setting assay_name to the metadata(ace)[['default_assay']]"))
      assay_name <- metadata(ace)[["default_assay"]]
    } else {
      message(sprintf("Input assay_name is NULL. Setting assay_name to logcounts"))
      assay_name <- "logcounts"
    }
  }
  S <- .validate_assay(ace, assay_name = assay_name, sparse_type = "CsparseMatrix", return_elem = TRUE)

  algorithm = match.arg(algorithm)
  algorithm = switch(algorithm, "irlb" = 0, "halko" = 1, "feng" = 2)

  if(is.null(max_it)) {
    max_it = ifelse(algorithm == 0, 1000, 5)
  }

  if (is.matrix(S)) {
    reduction.out <- C_reduceKernelDense(S, k = reduced_dim, svd_alg = algorithm, max_it = max_it, seed = seed, verbose = verbose)
  } else {
    reduction.out <- C_reduceKernelSparse(S, k = reduced_dim, svd_alg = algorithm, max_it = max_it, seed = seed, verbose = verbose)
  }

  S_r <- reduction.out$S_r
  colnames(S_r) <- colnames(ace)
  rownames(S_r) <- paste0("dim_", 1:NROW(S_r))
  colMaps(ace)[[reduction_out]] <- Matrix::t(S_r)
  colMapTypes(ace)[[reduction_out]] <- "reduction"

  metadata(ace)[["default_reduction"]] <- reduction_out
  metadata(ace)[["default_assay"]] <- assay_name

  V <- reduction.out$V
  colnames(V) <- paste0("V", 1:NCOL(V))
  rowMaps(ace)[[sprintf("%s_V", reduction_out)]] <- V
  rowMapTypes(ace)[[sprintf("%s_V", reduction_out)]] <- "internal"

  A <- reduction.out$A
  colnames(A) <- paste0("A", 1:NCOL(A))
  rowMaps(ace)[[sprintf("%s_A", reduction_out)]] <- A
  rowMapTypes(ace)[[sprintf("%s_A", reduction_out)]] <- "internal"

  B <- reduction.out$B
  colnames(B) <- paste0("B", 1:NCOL(B))
  colMaps(ace)[[sprintf("%s_B", reduction_out)]] <- B
  colMapTypes(ace)[[sprintf("%s_B", reduction_out)]] <- "internal"

  metadata(ace)[[sprintf("%s_sigma", reduction_out)]] <- reduction.out$sigma

  return(ace)
}
