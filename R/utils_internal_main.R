
#' Prune nonspecific and/or unreliable archetypes
.run.collectArchetypes <- function(ace,
                                   C_trace,
                                   H_trace,
                                   specificity_th = -3,
                                   min_cells_per_arch = 2) {
  .validate_ace(ace, allow_null = FALSE, return_elem = FALSE)

  pruning.out <- .collectArchetypes(
    C_trace = C_trace,
    H_trace = H_trace,
    specificity_th = specificity_th,
    min_cells_per_arch = min_cells_per_arch
  )

  colMaps(ace)[["H_stacked"]] <- Matrix::t(as(pruning.out$H_stacked, "sparseMatrix"))
  colMapTypes(ace)[["H_stacked"]] <- "internal"

  colMaps(ace)[["C_stacked"]] <- as(pruning.out$C_stacked, "sparseMatrix")
  colMapTypes(ace)[["C_stacked"]] <- "internal"

  return(ace)
}


#' Identiy equivalent classes of archetypes and group them together
.run.mergeArchetypes <- function(ace,
                                 reduction_slot = "action",
                                 C_stacked_slot = "C_stacked",
                                 H_stacked_slot = "H_stacked",
                                 normalization = 0,
                                 merged_suffix = "merged",
                                 footprint_slot_name = "assigned_archetype",
                                 thread_no = 0,
                                 return_raw = FALSE) {
  .validate_ace(ace, allow_null = FALSE, return_elem = FALSE)

  S_r <- .validate_map(
    ace = ace,
    map_slot = reduction_slot,
    matrix_type = "dense",
    force_type = TRUE
  )
  S_r <- Matrix::t(S_r)

  C_stacked <- .validate_map(
    ace = ace,
    map_slot = C_stacked_slot,
    matrix_type = "dense",
    force_type = TRUE
  )

  H_stacked <- .validate_map(
    ace = ace,
    map_slot = H_stacked_slot,
    matrix_type = "dense",
    force_type = TRUE
  )
  H_stacked <- Matrix::t(H_stacked)

  unification.out <- .mergeArchetypes(
    S_r = S_r,
    C_stacked = C_stacked,
    H_stacked = H_stacked,
    normalization = normalization,
    thread_no = thread_no
  )

  if (return_raw == TRUE) {
    return(unification.out)
  } else {
    Ht_merged <- as(Matrix::t(unification.out$H_merged), "sparseMatrix")
    colMaps(ace)[[sprintf("H_%s", merged_suffix)]] <- Ht_merged
    colMapTypes(ace)[[sprintf("H_%s", merged_suffix)]] <- "internal"

    colMaps(ace)[[sprintf("C_%s", merged_suffix)]] <- as(unification.out$C_merged, "sparseMatrix")
    colMapTypes(ace)[[sprintf("C_%s", merged_suffix)]] <- "internal"

    obs <- .get_obs_data(ace)
    obs[[footprint_slot_name]] <- c(unification.out$assigned_archetype)
    ace <- .set_obs_data(ace, obs)

    return(ace)
  }
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
