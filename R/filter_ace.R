#' Filter columns and rows of `ACTIONetExperiment` or `SummarizedExperiment`-like object.
#' @export
filter.ace <- function(
    ace,
    assay_name = "counts",
    min_cells_per_feat = NULL,
    min_feats_per_cell = NULL,
    min_umis_per_cell = NULL,
    max_umis_per_cell = NULL,
    return_fil_ace = TRUE) {
  init_dim <- dim(ace)
  init_dnames <- dimnames(ace)

  X <- .validate_assay(ace, assay_name = assay_name, sparse_type = "CsparseMatrix", return_elem = TRUE)

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
  ace <- ace[as.numeric(rownames(X)), as.numeric(colnames(X))]
  invisible(gc())

  if (return_fil_ace) {
    return(ace)
  } else {
    fil_cols_mask <- !(init_dnames[[2]] %in% colnames(ace))
    fil_rows_mask <- !(init_dnames[[1]] %in% rownames(ace))

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


#' Filter columns and rows of `ACTIONetExperiment` or `SummarizedExperiment` object by column attribute.
#' @export
filter.ace.by.attr <- function(
    ace,
    by,
    assay_name = "counts",
    min_cells_per_feat = NULL,
    min_feats_per_cell = NULL,
    min_umis_per_cell = NULL,
    max_umis_per_cell = NULL) {
  IDX <- ACTIONetExperiment::get.data.or.split(ace, attr = by, to_return = "split")

  if (any(duplicated(rownames(ace)))) {
    msg <- sprintf("Adding suffix to duplicate rownames.\n")
    warning(msg)
    rownames(ace) <- make.unique(rownames(ace))
  }
  if (any(duplicated(colnames(ace)))) {
    msg <- sprintf("Adding suffix to duplicate colnames.\n")
    warning(msg)
    colnames(ace) <- make.unique(colnames(ace))
  }

  fil_names <- lapply(IDX, function(idx) {
    fil_list <- filter.ace(
      ace = ace[, idx],
      assay_name = assay_name,
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

  keep_row <- which(!(rownames(ace) %in% fil_row))
  keep_col <- which(!(colnames(ace) %in% fil_col))

  ace <- ace[keep_row, keep_col]
  colData(ace) <- droplevels(colData(ace))
  rowData(ace) <- droplevels(rowData(ace))

  invisible(gc())
  return(ace)
}
