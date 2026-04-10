#' @export
runACTION <- function(
    adata = NULL,
    k_min = 2,
    k_max = 30,
    reduction_slot = "action",
    prenormalize = TRUE,
    min_obs = 2,
    max_it = 50,
    tol = 1e-100,
    spec_th = -3,
    thread_no = 0,
    return_c_matrices = TRUE,
    merged_suffix = "merged",
    archetype_slot_out = "assigned_archetype",
    ace = NULL) {
  adata <- .resolve_container_arg(adata = adata, ace = ace)
  adata <- .validate_ace(
    adata,
    as_ace = TRUE,
    allow_se_like = TRUE,
    return_elem = TRUE
  )

  S_r <- .validate_map(
    ace = adata,
    map_slot = reduction_slot,
    matrix_type = "dense",
    force_type = TRUE,
    return_elem = TRUE
  )

  if (prenormalize) {
    S_r <- normalize.matrix(
      S_r,
      dim = 1, # cells are rows in obsm (cells x k); normalize each cell vector
      scale_param = NULL,
      trans_func = NULL
    )
  }

  out <- C_runACTION(
    S_r = S_r,
    k_min = k_min,
    k_max = k_max,
    max_it = max_it,
    tol = tol,
    spec_th = spec_th,
    min_obs = min_obs,
    thread_no = thread_no,
    return_c_matrices = return_c_matrices
  )

  colMaps(adata)[["H_stacked"]] <- as(out$H_stacked, "sparseMatrix")
  colMapTypes(adata)[["H_stacked"]] <- "internal"

  colMaps(adata)[[sprintf("H_%s", merged_suffix)]] <- as(out$H_merged, "sparseMatrix")
  colMapTypes(adata)[[sprintf("H_%s", merged_suffix)]] <- "internal"

  if (return_c_matrices) {
    colMaps(adata)[["C_stacked"]] <- as(out$C_stacked, "sparseMatrix")
    colMapTypes(adata)[["C_stacked"]] <- "internal"

    colMaps(adata)[[sprintf("C_%s", merged_suffix)]] <- as(out$C_merged, "sparseMatrix")
    colMapTypes(adata)[[sprintf("C_%s", merged_suffix)]] <- "internal"
  }

  obs <- .get_obs_data(adata)
  obs[[archetype_slot_out]] <- c(out$assigned_archetypes)
  adata <- .set_obs_data(adata, obs)
  return(adata)
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
