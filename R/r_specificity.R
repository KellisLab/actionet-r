#' Computes feature specificity scores for each cluster.
#' @export
computeFeatureSpecificity <- function(
    adata = NULL,
    labels,
    labels_use = NULL,
    layer = "logcounts",
    map_out_prefix = "cluster",
    return_lower = FALSE,
    thread_no = 0,
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
  is_ace <- .validate_ace(adata, allow_se_like = TRUE, error_on_fail = FALSE, return_elem = FALSE)
  X <- .ace_or_assay(obj = adata, assay_name = layer, allow_se_like = TRUE, obj_name = "obj", return_elem = TRUE)

  labels <- .validate_vector_attr(
    adata,
    attr = labels,
    groups_use = labels_use,
    return_type = "data",
    attr_name = "labels",
    return_elem = TRUE
  )

  na_mask <- is.na(labels)
  if (any(na_mask)) {
    labels <- labels[!na_mask]
    X <- X[, !na_mask, drop = FALSE]
  }

  obs_factor <- factor(labels)
  obs_labels <- as.numeric(obs_factor)
  obs_keys <- levels(obs_factor)

  # Compute gene specificity for each cluster
  if (is.matrix(X)) {
    out <- C_computeFeatureSpecificityDense(X, labels = obs_labels, thread_no = thread_no)
  } else {
    out <- C_computeFeatureSpecificitySparse(X, labels = obs_labels, thread_no = thread_no)
  }

  out <- lapply(out, function(scores) {
    colnames(scores) <- obs_keys
    rownames(scores) <- .actionet_rownames(adata)
    return(scores)
  })

  if (is_ace && !return_raw) {
    spec_slot_out <- sprintf("%s_upper", map_out_prefix)
    rowMaps(adata)[[spec_slot_out]] <- out[["upper_significance"]]
    rowMapTypes(adata)[[spec_slot_out]] <- "reduction"

    if (return_lower && !is.null(out[["lower_significance"]])) {
      spec_slot_lower <- sprintf("%s_lower", map_out_prefix)
      rowMaps(adata)[[spec_slot_lower]] <- out[["lower_significance"]]
      rowMapTypes(adata)[[spec_slot_lower]] <- "reduction"
    }

    return(adata)
  }

  return(out)
}


#' Computes feature specificity scores for each archetype.
#' @export
archetypeFeatureSpecificity <- function(
    adata = NULL,
    layer = "logcounts",
    map_slot = "archetype_footprint",
    map_out_prefix = "archetype",
    thread_no = 0,
    return_raw = FALSE,
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
  adata <- .validate_ace(adata, allow_se_like = TRUE, error_on_fail = TRUE, return_elem = TRUE)
  X <- .validate_assay(ace = adata, assay_name = layer)

  H <- .validate_map(
    ace = adata,
    map_slot = map_slot,
    matrix_type = "dense",
    force_type = TRUE
  )
  H <- Matrix::t(H)

  if (is.matrix(X)) {
    out <- C_archetypeFeatureSpecificityDense(X, H = H, thread_no = thread_no)
  } else {
    out <- C_archetypeFeatureSpecificitySparse(X, H = H, thread_no = thread_no)
  }

  out <- lapply(out, function(scores) {
    colnames(scores) <- paste("A", seq_len(NCOL(scores)), sep = "")
    rownames(scores) <- .actionet_rownames(adata)
    return(scores)
  })


  if (!return_raw) {
    prof_slot_out <- sprintf("%s_feat_profile", map_out_prefix)
    rowMaps(adata)[[prof_slot_out]] <- out[["archetypes"]]
    rowMapTypes(adata)[[prof_slot_out]] <- "internal"

    spec_slot_out <- sprintf("%s_feat_specificity_upper", map_out_prefix)
    rowMaps(adata)[[spec_slot_out]] <- out[["upper_significance"]]
    rowMapTypes(adata)[[spec_slot_out]] <- "reduction"

    if (!is.null(out[["lower_significance"]])) {
      lower_slot_out <- sprintf("%s_feat_specificity_lower", map_out_prefix)
      rowMaps(adata)[[lower_slot_out]] <- out[["lower_significance"]]
      rowMapTypes(adata)[[lower_slot_out]] <- "reduction"
    }

    return(adata)
  }

  return(out)
}
