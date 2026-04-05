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
smoothKernel <- function(
    adata = NULL,
    norm_method = "pagerank",
    alpha = 0.85,
    max_it = 5,
    reduction_slot = "action",
    net_slot = "actionet",
    thread_no = 0,
    return_raw = FALSE,
    ace = NULL) {
  adata <- .resolve_container_arg(adata = adata, ace = ace)
  adata <- .validate_ace(adata, allow_se_like = TRUE, return_elem = TRUE, error_on_fail = TRUE)

  vars <- list(
    U = rowMaps(adata)[[sprintf("%s_U", reduction_slot)]],
    A = rowMaps(adata)[[sprintf("%s_A", reduction_slot)]],
    B = colMaps(adata)[[sprintf("%s_B", reduction_slot)]],
    sigma = .get_uns(adata)[[sprintf("%s_params", reduction_slot)]][["sigma"]]
  )

  if (any(sapply(vars, is.null))) {
    nullvars <- sprintf("%s_%s", reduction_slot, names(vars)[which(sapply(vars, is.null))])
    err <- sprintf("'%s' missing from 'ace'. Did you run 'reduceKernel()'?", paste(nullvars, collapse = ","))
    stop(err)
  }

  S_r <- .validate_map(
    adata,
    map_slot = reduction_slot,
    matrix_type = "dense",
    force_type = TRUE,
  )

  G <- .validate_net(
    adata,
    net_slot = net_slot,
    matrix_type = "sparse",
    force_type = TRUE,
  )

  U <- vars$U
  A <- vars$A
  B <- vars$B
  sigma <- vars$sigma

  U_right <- as.matrix(S_r %*% Matrix::Diagonal(length(sigma), 1.0 / sigma))
  SVD.out <- C_perturbedSVD(U, sigma, U_right, -A, B)
  V.smooth <- networkDiffusion(
    obj = G,
    scores = SVD.out$v,
    norm_method = norm_method,
    alpha = alpha,
    thread_no = thread_no,
    max_it = max_it
  )

  H <- V.smooth %*% diag(SVD.out$d)

  if (return_raw == TRUE) {
    out <- list(U = U, SVD.out = SVD.out, V.smooth = V.smooth, H = H)
    return(out)
  } else {
    W <- SVD.out$u
    rownames(W) <- .actionet_rownames(adata)
    smooth_red_name <- sprintf("%s_smooth", reduction_slot)
    smooth_U_name <- sprintf("%s_U", reduction_slot)
    rowMaps(adata)[[smooth_U_name]] <- W
    colMaps(adata)[[smooth_red_name]] <- H
    rowMapTypes(adata)[[smooth_U_name]] <- colMapTypes(adata)[[smooth_red_name]] <- "internal"
    return(adata)
  }
}


#' @export
runSVD <- function(
    X,
    k = 30,
    algorithm = c("irlb", "halko", "feng"),
    max_it = NULL,
    seed = 0,
    verbose = TRUE) {
  algorithm <- match.arg(algorithm)
  algorithm <- switch(algorithm,
    "irlb" = 0,
    "halko" = 1,
    "feng" = 2
  )

  if (is.null(max_it)) {
    max_it <- ifelse(algorithm == 0, 1000, 5)
  }

  X <- .validate_matrix(X)
  if (is.matrix(X)) {
    out <- C_runSVDDense(A = X, k = k, max_it = max_it, seed = seed, algorithm = algorithm, verbose = verbose)
  } else {
    out <- C_runSVDSparse(A = X, k = k, max_it = max_it, seed = seed, algorithm = algorithm, verbose = verbose)
  }

  return(out)
}
