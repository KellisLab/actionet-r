#' Compute the reduced kernel matrix and decomposition
#'
#' @param adata AnnData, `SummarizedExperiment`, `SingleCellExperiment`, or matrix-like input.
#' @param k Dimension of reduced kernel matrix. Number of singular vectors to estimate. Passed to <code>runSVD()</code>.
#' @param algorithm Singular value decomposition algorithm. Passed to <code>runSVD()</code>.
#' @param max_it Number of SVD iterations. If `NULL`: 1000 for "ilrb", 5 otherwise.
#' @param seed Random seed.
#' @param verbose Print status messages.
#' @param layer Layer to reduce. `NULL` uses `adata$X`.
#' @param reduction_slot Entry in `adata$obsm` in which to store the reduced kernel matrix, and the prefix for related decomposition outputs.
#' @param return_raw Return raw output regardless of container type.
#' @param assay_name Deprecated alias for `layer`.
#' @param obj Deprecated alias for `adata`.
#'
#' @return AnnData object with reduction stored in canonical ACTIONet slots, or the raw decomposition output when `return_raw = TRUE`.
#'
#' @export
reduceKernel <- function(
    adata = NULL,
    k = 30,
    algorithm = c("irlb", "halko", "feng"),
    max_it = NULL,
    seed = 0,
    verbose = TRUE,
    layer = "logcounts",
    reduction_slot = "action",
    return_raw = FALSE,
    assay_name = NULL,
    obj = NULL) {
  adata <- .resolve_container_arg(adata = adata, obj = obj)
  layer <- .resolve_layer_arg(
    layer = layer,
    assay_name = assay_name,
    default = "logcounts",
    layer_missing = missing(layer),
    assay_name_missing = missing(assay_name)
  )

  is_ace <- .is_se_like(adata)

  if (is_ace && !return_raw) {
    adata <- .validate_ace(adata, as_ace = TRUE, allow_se_like = TRUE, fix_dimnames = TRUE, return_elem = TRUE, error_on_fail = TRUE)
  }

  X <- .ace_or_assay(
    adata,
    assay_name = layer,
    sparse_type = "CsparseMatrix",
    allow_se_like = TRUE,
    return_elem = TRUE
  )

  # Bare matrix inputs are assumed to be legacy genes x cells orientation.
  # Transpose to cells x genes before calling C++ (which now expects cells x genes).
  # AnnData containers are already cells x genes and do not need transposition.
  if (!is_ace) {
    X <- if (.is_sparse_matrix(X)) Matrix::t(X) else t(as.matrix(X))
  }

  algorithm <- match.arg(algorithm, several.ok = TRUE)[1]
  algorithm <- switch(algorithm,
    "irlb" = 0,
    "halko" = 1,
    "feng" = 2
  )

  if (is.null(max_it)) {
    max_it <- ifelse(algorithm == 0, 1000, 5)
  }

  if (is.matrix(X)) {
    out <- C_reduceKernelDense(X, k = k, svd_alg = algorithm, max_it = max_it, seed = seed, verbose = verbose)
  } else {
    out <- C_reduceKernelSparse(X, k = k, svd_alg = algorithm, max_it = max_it, seed = seed, verbose = verbose)
  }

  if (is_ace && !return_raw) {
    S_r <- out$S_r
    rownames(S_r) <- .actionet_colnames(adata)   # cells are rows in obsm
    colnames(S_r) <- paste0("dim_", seq_len(NCOL(S_r)))
    colMaps(adata)[[reduction_slot]] <- S_r
    colMapTypes(adata)[[reduction_slot]] <- "reduction"

    V <- out$U
    colnames(V) <- paste0("U", seq_len(NCOL(V)))
    rowMaps(adata)[[sprintf("%s_U", reduction_slot)]] <- V
    rowMapTypes(adata)[[sprintf("%s_U", reduction_slot)]] <- "internal"

    A <- out$A
    colnames(A) <- paste0("A", seq_len(NCOL(A)))
    rowMaps(adata)[[sprintf("%s_A", reduction_slot)]] <- A
    rowMapTypes(adata)[[sprintf("%s_A", reduction_slot)]] <- "internal"

    B <- out$B
    colnames(B) <- paste0("B", seq_len(NCOL(B)))
    colMaps(adata)[[sprintf("%s_B", reduction_slot)]] <- B
    colMapTypes(adata)[[sprintf("%s_B", reduction_slot)]] <- "internal"

    uns <- .get_uns(adata)
    uns[[sprintf("%s_params", reduction_slot)]] <- c(.as_plain_list(uns[[sprintf("%s_params", reduction_slot)]]), list(sigma = out$sigma))
    adata <- .set_uns(adata, uns)

    return(adata)
  }

  return(out)
}


#' @export
decompACTION <- function(
    X,
    k_min,
    k_max,
    max_it = 100,
    tol = 1e-16,
    thread_no = 0) {
  if (k_min < 2) {
    stop("'k_min' must be >=2")
  }

  if (k_max < k_min) {
    stop("'k_max' must be >='k_min'")
  }

  .validate_matrix(X, var_name = "X", matrix_type = "dense", force_type = FALSE, return_elem = FALSE)

  out <- C_decompACTION(
    S_r = X,
    k_min = k_min,
    k_max = k_max,
    max_it = max_it,
    tol = tol,
    thread_no = thread_no
  )

  return(out)
}


#' @export
collectArchetypes <- function(C_trace,
                              H_trace,
                              specificity_th = -3,
                              min_cells_per_arch = 2) {
  out <- C_collectArchetypes(
    C_trace = C_trace,
    H_trace = H_trace,
    spec_th = specificity_th,
    min_obs = min_cells_per_arch
  )

  return(out)
}


#' @export
mergeArchetypes <- function(
    S_r,
    C_stacked,
    H_stacked,
    thread_no = 0) {
  if (!all(dim(C_stacked) == rev(dim(H_stacked)))) {
    err <- sprintf("Dimensions of `C_stacked` (%s) and `H_stacked` (%s) are incompatible: C must be archetypes x cells_or_k and H must be cells_or_k x archetypes.\n",
                   paste(dim(C_stacked), collapse = "x"),
                   paste(dim(H_stacked), collapse = "x"))
    stop(err)
  }

  out <- C_mergeArchetypes(
    S_r = S_r,
    C_stacked = C_stacked,
    H_stacked = H_stacked,
    thread_no = thread_no
  )

  return(out)
}
