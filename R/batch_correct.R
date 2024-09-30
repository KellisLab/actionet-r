#' Perform batch correction on `ACTIONetExperiment` and `SingleCellExperiment` objects.

#' @export
correctBatchEffectFastMNN <- function(
    ace,
    batch_attr = NULL,
    assay_name = NULL,
    reduced_dim = 50,
    MNN_k = 20,
    reduction_out = "MNN",
    BPPARAM = SerialParam()) {

  ACTIONetExperiment:::.check_and_load_package(c("scran", "SingleCellExperiment", "batchelor", "BiocParallel"))
  ace <- .validate_ace(ace, as_ace = TRUE, allow_se_like = TRUE, fix_dimnames = TRUE, return_elem = TRUE)
  S <- .validate_assay(ace, assay_name = assay_name, return_elem = TRUE)
  batch_vec <- .validate_vector_attr(
    ace,
    attr = batch_attr,
    return_type = "data",
    attr_name = "batch_attr",
    obj_name = "ace",
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
  rownames(S_r) <- colnames(ace)
  colnames(S_r) <- sapply(seq_len(dim(S_r)[2]), function(i) sprintf("PC%d", i))

  colMaps(ace)[[reduction_out]] <- S_r
  colMapTypes(ace)[[reduction_out]] <- "reduction"

  V <- rowData(mnn.out)[["rotation"]]
  colnames(V) <- paste0("V", seq_len(NCOL(V)))
  rowMaps(ace)[[sprintf("%s_V", reduction_out)]] <- V
  rowMapTypes(ace)[[sprintf("%s_V", reduction_out)]] <- "internal"

  invisible(gc())

  return(ace)
}

# TODO: Allow design_mat to be multiple input types.
#' @export
correctBatchEffect <- function(
    ace,
    design_mat,
    reduction_slot = "action",
    corrected_suffix = "orth",
    assay_name = NULL) {
  ace <- .validate_ace(ace, allow_se_like = FALSE, fix_dimnames = TRUE, return_elem = TRUE, error_on_fail = TRUE)
  S <- .validate_assay(ace, assay_name = assay_name, error_on_fail = TRUE, return_elem = TRUE)
  S_r <- .validate_map(ace, map_slot = reduction_slot, matrix_type = "dense", force_type = TRUE, return_elem = TRUE)

  # V <- rowMaps(ace)[[sprintf("%s_V", reduction_slot)]]
  # A <- rowMaps(ace)[[sprintf("%s_A", reduction_slot)]]
  # B <- colMaps(ace)[[sprintf("%s_B", reduction_slot)]]
  # sigma <- S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]

  B <- .validate_map(ace, map_slot = sprintf("%s_B", reduction_slot), matrix_type = "dense", force_type = TRUE, return_elem = TRUE, row = FALSE)
  V <- .validate_map(ace, map_slot = sprintf("%s_V", reduction_slot), matrix_type = "dense", force_type = TRUE, return_elem = TRUE, row = TRUE)
  A <- .validate_map(ace, map_slot = sprintf("%s_A", reduction_slot), matrix_type = "dense", force_type = TRUE, return_elem = TRUE, row = TRUE)
  sigma <- S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]

  if (length(sigma) != NCOL(S_r)) {
    err <- sprintf("Size of 'sigma' in 'metadata(ace)' does not match reduction\nRecommend re-running 'reduceKernel()'")
    stop(err)
  }

  if (is.matrix(S)) {
    out <- C_orthogonalizeBatchEffect_full(
      S = S,
      old_S_r = S_r,
      old_V = V,
      old_A = A,
      old_B = B,
      old_sigma = sigma,
      design = design_mat
    )
  } else {
    out <- C_orthogonalizeBatchEffect(
      S = S,
      old_S_r = S_r,
      old_V = V,
      old_A = A,
      old_B = B,
      old_sigma = sigma,
      design = design_mat
    )
  }
  S_r <- out$S_r
  colnames(S_r) <- colnames(ace)
  rownames(S_r) <- sapply(seq_len(NROW(S_r)), function(i) sprintf("Dim%d", i))
  name_Sr <- sprintf("%s_%s", reduction_slot, corrected_suffix)
  colMaps(ace)[[name_Sr]] <- Matrix::t(S_r)
  colMapTypes(ace)[[name_Sr]] <- "reduction"


  V <- out$V
  colnames(V) <- sapply(seq_len(dim(V)[2]), function(i) sprintf("V%d", i))
  name_V <- sprintf("%s_V_%s", reduction_slot, corrected_suffix)
  rowMaps(ace)[[name_V]] <- V
  rowMapTypes(ace)[[name_V]] <- "internal"


  A <- out$A
  colnames(A) <- sapply(seq_len(dim(A)[2]), function(i) sprintf("A%d", i))
  name_A <- sprintf("%s_A_%s", reduction_slot, corrected_suffix)
  rowMaps(ace)[[name_A]] <- A
  rowMapTypes(ace)[[name_A]] <- "internal"


  B <- out$B
  colnames(B) <- sapply(seq_len(dim(B)[2]), function(i) sprintf("B%d", i))
  name_B <- sprintf("%s_B_%s", reduction_slot, corrected_suffix)
  colMaps(ace)[[name_B]] <- B
  colMapTypes(ace)[[name_B]] <- "internal"

  name_sigma <- sprintf("%s_sigma_%s", reduction_slot, corrected_suffix)
  S4Vectors::metadata(ace)[[name_sigma]] <- out$sigma

  return(ace)
}

# orthogonalize.ace.batch.simple <- function(ace,
#                                            batch_attr,
#                                            reduction_slot = "action",
#                                            corrected_out = "ACTION_ortho",
#                                            assay_name = NULL) {
#   ace <- as(ace, "ACTIONetExperiment")
#
#   if (is.null(assay_name)) {
#     if ("default_assay" %in% names(metadata(ace))) {
#       message(sprintf("Input assay_name is NULL. Setting assay_name to the metadata(ace)[['default_assay']]"))
#       assay_name <- metadata(ace)[["default_assay"]]
#     } else {
#       message(sprintf("Input assay_name is NULL. Setting assay_name to logcounts"))
#       assay_name <- "logcounts"
#     }
#   }
#   .validate_assay(ace, assay_name = assay_name, return_elem = FALSE)
#
#
#   batch_attr <- ACTIONetExperiment::get.data.or.split(ace, attr = batch_attr, to_return = "data")
#   batch_attr <- as.factor(batch_attr)
#   design_mat <- stats::model.matrix(~batch_attr)
#
#   ace <- correctBatchEffect(
#     ace,
#     design_mat,
#     reduction_slot = reduction_slot,
#     corrected_out = corrected_out,
#     assay_name = assay_name
#   )
#
#   return(ace)
# }

#' reduce.and.batch.orthogonalize.ace <- function(ace,
#'                                                design_mat,
#'                                                reduced_dim = 50,
#'                                                max_iter = 1000,
#'                                                assay_name = NULL,
#'                                                reduction_out = "action",
#'                                                corrected_out = "ACTION_ortho",
#'                                                seed = 0,
#'                                                SVD_algorithm = 0) {
#'   ace <- as(ace, "ACTIONetExperiment")
#'
#'   if (is.null(assay_name)) {
#'     if ("default_assay" %in% names(metadata(ace))) {
#'       message(sprintf("Input assay_name is NULL. Setting assay_name to the metadata(ace)[['default_assay']]"))
#'       assay_name <- metadata(ace)[["default_assay"]]
#'     } else {
#'       message(sprintf("Input assay_name is NULL. Setting assay_name to logcounts"))
#'       assay_name <- "logcounts"
#'     }
#'   }
#'   .validate_assay(ace, assay_name = assay_name, return_elem = FALSE)
#'
#'   if (!is.matrix(design_mat)) {
#'     err <- sprintf("'design_mat' must be a matrix.\n")
#'     stop(err)
#'   }
#'
#'   ace <- reduce.ace(
#'     ace = ace,
#'     reduced_dim = reduced_dim,
#'     max_iter = max_iter,
#'     assay_name = assay_name,
#'     reduction_out = reduction_out,
#'     seed = seed,
#'     SVD_algorithm = SVD_algorithm
#'   )
#'
#'   ace <- correctBatchEffect(
#'     ace = ace,
#'     design_mat = design_mat,
#'     reduction_slot = reduction_out,
#'     corrected_out = corrected_out,
#'     assay_name = assay_name
#'   )
#'
#'   return(ace)
#' }
