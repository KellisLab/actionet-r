#' Perform batch correction on AnnData or `SingleCellExperiment`-like inputs.

#' @export
correctBatchEffectFastMNN <- function(
    adata = NULL,
    batch_attr = NULL,
    layer = NULL,
    reduced_dim = 50,
    MNN_k = 20,
    reduction_out = "MNN",
    BPPARAM = SerialParam(),
    assay_name = NULL,
    ace = NULL) {
  .check_and_load_package(c("scran", "SingleCellExperiment", "batchelor", "BiocParallel"))
  adata <- .resolve_container_arg(adata = adata, ace = ace)
  layer <- .resolve_layer_arg(
    layer = layer,
    assay_name = assay_name,
    default = NULL,
    layer_missing = missing(layer),
    assay_name_missing = missing(assay_name)
  )
  adata <- .validate_ace(adata, as_ace = TRUE, allow_se_like = TRUE, fix_dimnames = TRUE, return_elem = TRUE)
  S <- .validate_assay(adata, assay_name = layer, return_elem = TRUE)
  batch_vec <- .validate_vector_attr(
    adata,
    attr = batch_attr,
    return_type = "data",
    attr_name = "batch_attr",
    obj_name = "adata",
    return_elem = TRUE
  )

  IDX <- split(seq_len(length(batch_vec)), f = batch_vec)
  merge_order <- order(sapply(IDX, function(idx) length(idx)), decreasing = TRUE)

  set.seed(0)
  mnn.out <- batchelor::fastMNN(
    S,
    batch = batch_vec,
    k = MNN_k,
    d = reduced_dim,
    auto.merge = FALSE,
    merge.order = merge_order,
    cos.norm = FALSE,
    BPPARAM = BPPARAM
  )

  S_r <- SingleCellExperiment::reducedDims(mnn.out)[["corrected"]]
  rownames(S_r) <- .actionet_colnames(adata)
  colnames(S_r) <- sapply(seq_len(dim(S_r)[2]), function(i) sprintf("PC%d", i))

  colMaps(adata)[[reduction_out]] <- S_r
  colMapTypes(adata)[[reduction_out]] <- "reduction"

  V <- rowData(mnn.out)[["rotation"]]
  colnames(V) <- paste0("V", seq_len(NCOL(V)))
  rowMaps(adata)[[sprintf("%s_V", reduction_out)]] <- V
  rowMapTypes(adata)[[sprintf("%s_V", reduction_out)]] <- "internal"

  invisible(gc())

  return(adata)
}

# TODO: Allow design_mat to be multiple input types.
#' @export
correctBatchEffect <- function(
    adata = NULL,
    batches = NULL,
    design = NULL,
    reduction_slot = "action",
    corrected_suffix = "orth",
    layer = "logcounts",
    assay_name = NULL,
    ace = NULL) {
  adata <- .resolve_container_arg(adata = adata, ace = ace)
  layer <- .resolve_layer_arg(
    layer = layer,
    assay_name = assay_name,
    default = "logcounts",
    layer_missing = missing(layer),
    assay_name_missing = missing(assay_name)
  )
  if (!is.null(batches) && !is.null(design)) {
    err <- sprintf("Only one of 'batches' or 'design' must be specified")
    stop(err)
  } else {
    if (!is.null(batches)) {
      batches <- .validate_vector_attr(
        adata,
        attr = batches,
        return_type = "data",
        keep_factor = TRUE,
        attr_name = "batches",
        obj_name = "adata",
        return_elem = TRUE
      )
      batches <- as.factor(batches)
      if (length(levels(batches)) < 2) {
        err <- sprintf("'batches' must have >=2 levels")
        stop(err)
      }
      design_mat <- model.matrix(~ 0 + batches)
    } else {
      design_mat <- .make_design_mat(design = design, data = .get_obs_data(adata), remove_intercept = TRUE)
      if (NROW(design_mat) != .actionet_ncol(adata)) {
        err <- sprintf("Size of 'design' does not match the number of observations in 'adata'")
      }
    }
  }

  adata <- .validate_ace(adata, allow_se_like = FALSE, fix_dimnames = TRUE, return_elem = TRUE, error_on_fail = TRUE)
  S <- .validate_assay(adata, assay_name = layer, error_on_fail = TRUE, return_elem = TRUE)
  S_r <- .validate_map(adata, map_slot = reduction_slot, matrix_type = "dense", force_type = TRUE, return_elem = TRUE)

  B <- .validate_map(adata, map_slot = sprintf("%s_B", reduction_slot), matrix_type = "dense", force_type = TRUE, return_elem = TRUE, row = FALSE)
  U <- .validate_map(adata, map_slot = sprintf("%s_U", reduction_slot), matrix_type = "dense", force_type = TRUE, return_elem = TRUE, row = TRUE)
  A <- .validate_map(adata, map_slot = sprintf("%s_A", reduction_slot), matrix_type = "dense", force_type = TRUE, return_elem = TRUE, row = TRUE)
  sigma <- .get_uns(adata)[[sprintf("%s_params", reduction_slot)]][["sigma"]]

  if (length(sigma) != NCOL(S_r)) {
    err <- sprintf("Size of 'sigma' in action params does not match reduction\nRecommend re-running 'reduceKernel()'")
    stop(err)
  }

  if (is.matrix(S)) {
    out <- C_orthogonalizeBatchEffect_full(
      S = S,
      old_S_r = S_r,
      old_U = U,
      old_A = A,
      old_B = B,
      old_sigma = sigma,
      design = design_mat
    )
  } else {
    out <- C_orthogonalizeBatchEffect(
      S = S,
      old_S_r = S_r,
      old_U = U,
      old_A = A,
      old_B = B,
      old_sigma = sigma,
      design = design_mat
    )
  }
  S_r <- out$S_r
  colnames(S_r) <- .actionet_colnames(adata)
  rownames(S_r) <- sapply(seq_len(NROW(S_r)), function(i) sprintf("Dim%d", i))
  name_Sr <- sprintf("%s_%s", reduction_slot, corrected_suffix)
  colMaps(adata)[[name_Sr]] <- Matrix::t(S_r)
  colMapTypes(adata)[[name_Sr]] <- "reduction"


  U <- out$U
  colnames(U) <- sapply(seq_len(dim(U)[2]), function(i) sprintf("U%d", i))
  name_U <- sprintf("%s_U_%s", reduction_slot, corrected_suffix)
  rowMaps(adata)[[name_U]] <- U
  rowMapTypes(adata)[[name_U]] <- "internal"


  A <- out$A
  colnames(A) <- sapply(seq_len(dim(A)[2]), function(i) sprintf("A%d", i))
  name_A <- sprintf("%s_A_%s", reduction_slot, corrected_suffix)
  rowMaps(adata)[[name_A]] <- A
  rowMapTypes(adata)[[name_A]] <- "internal"


  B <- out$B
  colnames(B) <- sapply(seq_len(dim(B)[2]), function(i) sprintf("B%d", i))
  name_B <- sprintf("%s_B_%s", reduction_slot, corrected_suffix)
  colMaps(adata)[[name_B]] <- B
  colMapTypes(adata)[[name_B]] <- "internal"

  name_sigma <- sprintf("%s_%s_params", reduction_slot, corrected_suffix)
  uns <- .get_uns(adata)
  uns[[name_sigma]] <- list(sigma = out$sigma)
  adata <- .set_uns(adata, uns)

  return(adata)
}
