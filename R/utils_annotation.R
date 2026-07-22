.encode_markers <- function(
    obj,
    markers,
    features_use = NULL,
    obj_name = "obj",
    features_name = "features_use") {
  if (any(class(markers) == "list")) {
    mn <- names(markers)
  } else {
    mn <- colnames(markers)
  }

  if (is.null(mn) || any(is.na(mn))) {
    err <- sprintf("'markers' contains unnamed entries")
    stop(err)
  }

  if (any(duplicated(mn))) {
    err <- sprintf("'markers' contains duplicated labels")
    stop(err)
  }

  features_use <- .get_features(obj, features_use = features_use, allow_empty = FALSE)
  if (any(class(markers) %in% c("list", "data.frame"))) {
    X <- sapply(markers, function(x) {
      as.numeric(features_use %in% x)
    })
    X <- as(X, "CsparseMatrix")
  } else if (is.matrix(markers) || .is_sparse_matrix(markers)) {
    if (any(!is.finite(markers))) {
      err <- sprintf("'markers' contains non-numeric values")
      stop(err)
    }

    n_features <- if (.is_anndata(obj)) .n_vars(obj) else nrow(obj)
    if (NROW(markers) != n_features) {
      err <- sprintf("NROW(%s) does not match NROW(markers)", obj_name)
      stop(err)
    }
    X <- as(markers, "dMatrix")
    X[X != 0] <- 1
  } else {
    err <- sprintf("'markers' must be one of: 'list', 'data.frame', 'matrix', 'sparseMatrix'")
    stop(err)
  }

  rownames(X) <- features_use

  cs_zero <- (Matrix::colSums(X) == 0)
  if (any(cs_zero)) {
    if (all(cs_zero)) {
      err <- sprintf("No markers in '%s'", features_name)
      stop(err)
    }

    dropped <- names(which(cs_zero))
    for (k in dropped) {
      wrn <- sprintf("Label '%s' has no markers", k)
      message(wrn)
    }
  }

  return(X)
}


# Shared marker-branch pipeline for `annotateClusters` and `annotateArchetypes`.
#
# Mirrors Python `_annotate_from_markers` in
# `actionet-python/src/actionet/annotation/annotation.py`. Given upper/lower
# feature-specificity matrices (rows = features, cols = groups), a marker
# specification, and the group axis row names, it forms
# `pmax(upper - lower, 0)`, encodes markers via `.encode_markers`, calls
# `C_assess_enrichment`, and returns the standard annotation list.
#
# Args:
#   adata: AnnData/ACE-like container passed through to `.encode_markers`.
#   markers: marker specification (list/data.frame/matrix); see
#     `.encode_markers`.
#   features_use: feature-labels vector or `var` column name; forwarded to
#     `.encode_markers`.
#   upper_sig, lower_sig: feature x group specificity matrices. `lower_sig`
#     may be `NULL`, in which case `feat_spec = upper_sig`.
#   row_names: character vector of length `ncol(upper_sig)` used to name
#     rows of the returned enrichment matrix.
#   thread_no: threads for `C_assess_enrichment`.
.annotate_from_markers <- function(adata,
                                   markers,
                                   features_use,
                                   upper_sig,
                                   lower_sig,
                                   row_names,
                                   thread_no) {
  upper_sig <- as.matrix(upper_sig)
  if (!is.null(lower_sig)) {
    lower_sig <- as.matrix(lower_sig)
    feat_spec <- upper_sig - lower_sig
  } else {
    feat_spec <- upper_sig
  }
  feat_spec[feat_spec < 0] <- 0
  colnames(feat_spec) <- row_names

  marker_mat <- .encode_markers(
    adata,
    markers = markers,
    features_use = features_use,
    obj_name = "adata"
  )
  marker_mat <- as(marker_mat, "CsparseMatrix")

  # Row-count parity check between marker matrix and specificity matrix.
  # Both are constructed row-by-row over the full feature axis of `adata`
  # (see caller comments in R/annotation.R), so they are positionally
  # aligned even when their rownames differ.
  if (nrow(marker_mat) != nrow(feat_spec)) {
    stop(sprintf(
      "Feature axis mismatch: marker matrix has %d rows but specificity matrix has %d rows.",
      nrow(marker_mat), nrow(feat_spec)
    ))
  }
  if (!is.null(rownames(marker_mat))) {
    rownames(feat_spec) <- rownames(marker_mat)
  }

  enrich <- C_assess_enrichment(
    scores = feat_spec,
    associations = marker_mat,
    thread_no = thread_no
  )
  enrichment <- Matrix::t(enrich$logPvals)
  rownames(enrichment) <- colnames(feat_spec)
  colnames(enrichment) <- colnames(marker_mat)
  enrichment[!is.finite(enrichment)] <- 0

  annots <- colnames(enrichment)[apply(enrichment, 1, which.max)]
  conf <- apply(enrichment, 1, max)
  names(annots) <- rownames(enrichment)
  names(conf) <- rownames(enrichment)

  list(
    labels = annots,
    confidence = conf,
    enrichment = enrichment
  )
}
