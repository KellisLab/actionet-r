#' @export
run.SVD <- function(
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


decomp.simplex <- function(X,
                           B,
                           return_raw = FALSE) {
  if (!is.matrix(X) || !is.matrix(B)) {
    err <- sprintf("'X' and 'W0' must be of type 'matrix'.\n")
    warning(err)
  }

  H <- run_simplex_regression(A = as.matrix(X), B = as.matrix(B))

  if (return_raw == TRUE) {
    out <- H
  } else {
    out <- list(W = X, H = H)
  }
  return(out)
}


decomp.SPA <- function(X,
                       k = 30,
                       return_raw = FALSE) {
  if (!is.matrix(X)) {
    err <- sprintf("'X' must be of type 'matrix'.\n")
    warning(err)
    X <- as.matrix(X)
  }

  SPA.out <- run_SPA(X, k)

  if (return_raw == TRUE) {
    out <- SPA.out
  } else {
    W <- X[, SPA.out$selected_columns]
    H <- run_simplex_regression(W, X)
    misc <- list(selected_cols = SPA.out$selected_columns, norms = SPA.out$norms)

    out <- list(W = W, H = H, misc = misc)
  }
  return(out)
}

decomp.AA <- function(X,
                      k = NULL,
                      W0 = NULL,
                      max_iter = 50,
                      min_delta = 1e-300,
                      return_raw = FALSE) {
  if (!is.matrix(X)) {
    err <- sprintf("'X' must be of type 'matrix'.\n")
    warning(err)
    X <- as.matrix(X)
  }

  if (!is.null(W0) && !is.matrix(W0)) {
    err <- sprintf("'W0' must be of type 'matrix' if given.\n")
    warning(err)
    W0 <- as.matrix(W0)
  }

  if (is.null(W0)) {
    if (is.null(k)) {
      err <- sprintf("'k' cannot be 'NULL' if 'W0' is not given.\n")
      stop(err)
    }
    W0 <- matrix(0, NROW(X), k)
  }

  AA.out <- run_AA(X, W0 = W0, max_it = max_iter, min_delta = min_delta)

  if (return_raw == TRUE) {
    out <- AA.out
  } else {
    out <- list(
      W = AA.out$W,
      H = AA.out$H,
      misc = list(C = AA.out$C)
    )
  }
  return(out)
}

#' Run ACTION_decomposition
decomp.ACTION <- function(X,
                          k_min,
                          k_max,
                          specificity_th = -3,
                          min_cells_per_arch = 2,
                          max_it = 100,
                          tol = 1e-16,
                          thread_no = 0,
                          norm = 1) {
  if (!is.matrix(X)) {
    err <- sprintf("'X' must be of type 'matrix'.\n")
    warning(err)
    X <- as.matrix(X)
  }

  action.out <- decompACTION(
    S_r = X,
    k_min = k_min,
    k_max = k_max,
    max_it = max_it,
    tol = tol,
    thread_no = thread_no
  )

  # Prune nonspecific and/or unreliable archetypes
  pruning.out <- .collectArchetypes(
    C_trace = action.out$C,
    H_trace = action.out$H,
    specificity_th = specificity_th,
    min_cells_per_arch = min_cells_per_arch
  )
  message("")
  # Identiy equivalent classes of archetypes and group them together
  C_stacked <- pruning.out$C_stacked
  H_stacked <- pruning.out$H_stacked
  unification.out <- .mergeArchetypes(
    S_r = X,
    C_stacked = C_stacked,
    H_stacked = H_stacked,
    normalization = 0,
    thread_no = thread_no
  )

  H <- unification.out$H_merged
  W <- X %*% unification.out$C_merged
  misc <- list(
    H = action.out$H,
    C = action.out$C,
    H_stacked = pruning.out$H_stacked,
    C_stacked = pruning.out$C_stacked,
    H_merged = unification.out$H_merged,
    C_merged = unification.out$C_merged,
    assigned_archetype = unification.out$assigned_archetype
  )
  out <- list(W = W, H = H, misc = misc)

  return(out)
}


.collectArchetypes <- function(C_trace,
                               H_trace,
                               specificity_th = -3,
                               min_cells_per_arch = 2) {
  out <- collectArchetypes(
    C_trace = C_trace,
    H_trace = H_trace,
    spec_th = specificity_th,
    min_obs = min_cells_per_arch
  )

  return(out)
}


.mergeArchetypes <- function(S_r,
                             C_stacked,
                             H_stacked,
                             normalization = 0,
                             thread_no = 0) {
  if (!all(dim(C_stacked) == rev(dim(H_stacked)))) {
    err <- sprintf("Dimensions of `C_stacked` and transposed `H_stacked` do not match.\n")
    stop(err)
  }

  out <- mergeArchetypes(
    S_r = S_r,
    C_stacked = C_stacked,
    H_stacked = H_stacked,
    norm = normalization,
    thread_no = thread_no
  )

  return(out)
}
