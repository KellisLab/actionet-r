#' Filter rows and columns of an ACTIONet-compatible container.
#'
#' @param adata AnnData or compatible container.
#' @param layer Layer to filter on. `NULL` uses `adata$X`.
#' @param assay_name Deprecated alias for `layer`.
#' @param ace Deprecated alias for `adata`.
#' @export
filterActionet <- function(
    adata = NULL,
    layer = NULL,
    min_cells_per_feat = NULL,
    min_feats_per_cell = NULL,
    min_umis_per_cell = NULL,
    max_umis_per_cell = NULL,
    return_fil_ace = TRUE,
    assay_name = NULL,
    ace = NULL) {
  adata <- .resolve_container_arg(adata = adata, ace = ace)
  layer <- .resolve_layer_arg(
    layer = layer,
    assay_name = assay_name,
    default = NULL,
    layer_missing = missing(layer),
    assay_name_missing = missing(assay_name)
  )
  adata <- .validate_ace(adata, as_ace = TRUE, allow_se_like = TRUE, fix_dimnames = TRUE, return_elem = TRUE)
  init_dim <- .actionet_dim(adata)
  init_dnames <- list(.actionet_rownames(adata), .actionet_colnames(adata))

  X <- .validate_assay(adata, assay_name = layer, sparse_type = "CsparseMatrix", return_elem = TRUE)

  dimnames(X) <- list(1:NROW(X), 1:NCOL(X))

  i <- 0
  repeat {
    prev_dim <- dim(X)
    rows_mask <- rep(TRUE, NROW(X))
    cols_mask <- rep(TRUE, NCOL(X))
    if (!is.null(min_umis_per_cell)) {
      umi_mask <- Matrix::colSums(X) >= min_umis_per_cell
      cols_mask <- cols_mask & umi_mask
    }

    if (!is.null(max_umis_per_cell)) {
      umi_mask <- Matrix::colSums(X) <= max_umis_per_cell
      cols_mask <- cols_mask & umi_mask
    }

    if (!is.null(min_feats_per_cell)) {
      feature_mask <- Matrix::colSums(.validate_matrix(X > 0)) >= min_feats_per_cell
      cols_mask <- cols_mask & feature_mask
    }

    if (!is.null(min_cells_per_feat)) {
      if ((min_cells_per_feat < 1) & (min_cells_per_feat > 0)) {
        min_fc <- ceiling(min_cells_per_feat * prev_dim[2])
      } else {
        min_fc <- min_cells_per_feat
      }
      cell_count_mask <- Matrix::rowSums(.validate_matrix(X > 0)) >= min_fc
      rows_mask <- rows_mask & cell_count_mask
    }

    X <- X[rows_mask, cols_mask]
    invisible(gc())
    i <- i + 1
    if (all(dim(X) == prev_dim)) {
      break
    }
  }
  adata <- .subset_actionet_container(
    adata,
    features = as.numeric(rownames(X)),
    cells = as.numeric(colnames(X))
  )
  invisible(gc())

  if (return_fil_ace) {
    return(adata)
  } else {
    fil_cols_mask <- !(init_dnames[[2]] %in% .actionet_colnames(adata))
    fil_rows_mask <- !(init_dnames[[1]] %in% .actionet_rownames(adata))

    fil_cols_list <- data.frame(
      name = init_dnames[[2]][fil_cols_mask],
      idx = which(fil_cols_mask)
    )

    fil_rows_list <- data.frame(
      name = init_dnames[[1]][fil_rows_mask],
      idx = which(fil_rows_mask)
    )

    fil_list <- list(
      cols_filtered = fil_cols_list,
      rows_filtered = fil_rows_list
    )

    return(fil_list)
  }
}


#' Filter rows and columns of an ACTIONet-compatible container by observation attribute.
#'
#' @param adata AnnData or compatible container.
#' @param by Observation attribute used to define groups before filtering.
#' @param layer Layer to filter on. `NULL` uses `adata$X`.
#' @param assay_name Deprecated alias for `layer`.
#' @param ace Deprecated alias for `adata`.
#' @export
filterActionetByAttr <- function(
    adata = NULL,
    by,
    layer = NULL,
    min_cells_per_feat = NULL,
    min_feats_per_cell = NULL,
    min_umis_per_cell = NULL,
    max_umis_per_cell = NULL,
    assay_name = NULL,
    ace = NULL) {
  adata <- .resolve_container_arg(adata = adata, ace = ace)
  layer <- .resolve_layer_arg(
    layer = layer,
    assay_name = assay_name,
    default = NULL,
    layer_missing = missing(layer),
    assay_name_missing = missing(assay_name)
  )
  adata <- .validate_ace(adata, as_ace = TRUE, allow_se_like = TRUE, fix_dimnames = TRUE, return_elem = TRUE)
  IDX <- .validate_vector_attr(adata, attr = by, return_type = "split", dim = 2)

  if (any(duplicated(.actionet_rownames(adata)))) {
    msg <- sprintf("Adding suffix to duplicate rownames.\n")
    warning(msg)
    adata <- .set_actionet_rownames(adata, make.unique(.actionet_rownames(adata)))
  }
  if (any(duplicated(.actionet_colnames(adata)))) {
    msg <- sprintf("Adding suffix to duplicate colnames.\n")
    warning(msg)
    adata <- .set_actionet_colnames(adata, make.unique(.actionet_colnames(adata)))
  }

  fil_names <- lapply(IDX, function(idx) {
    fil_list <- filterActionet(
      adata = .subset_actionet_container(adata, cells = idx),
      layer = layer,
      min_cells_per_feat = min_cells_per_feat,
      min_umis_per_cell = min_umis_per_cell,
      max_umis_per_cell = max_umis_per_cell,
      min_feats_per_cell = min_feats_per_cell,
      return_fil_ace = FALSE
    )

    return(fil_list)
  })

  fil_col <- lapply(fil_names, function(i) i[["cols_filtered"]]$name)
  fil_col <- Reduce(union, fil_col)


  fil_row <- lapply(fil_names, function(i) i[["rows_filtered"]]$name)
  fil_row <- Reduce(union, fil_row)

  keep_row <- which(!(.actionet_rownames(adata) %in% fil_row))
  keep_col <- which(!(.actionet_colnames(adata) %in% fil_col))

  adata <- .subset_actionet_container(adata, features = keep_row, cells = keep_col)
  adata <- .set_obs_data(adata, droplevels(.get_obs_data(adata)))
  adata <- .set_feature_data(adata, droplevels(.get_feature_data(adata)))

  invisible(gc())
  return(adata)
}

#' @export
filter.ace <- function(...) {
  .Deprecated("filterActionet")
  filterActionet(...)
}

#' @export
filter.ace.by.attr <- function(...) {
  .Deprecated("filterActionetByAttr")
  filterActionetByAttr(...)
}
