.is_se_like <- function(obj) {
  allowed_classes <- c(
    "AbstractAnnData",
    "ACTIONetExperiment",
    "SummarizedExperiment",
    "RangedSummarizedExperiment",
    "SingleCellExperiment"
  )

  if (any(class(obj) %in% allowed_classes)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

.validate_ace <- function(
    obj,
    as_ace = FALSE,
    allow_se_like = FALSE,
    allow_null = FALSE,
    obj_name = "ace",
    fix_dimnames = FALSE,
    return_elem = FALSE,
    error_on_fail = TRUE) {
  if (is.null(obj) && !allow_null) {
    err <- sprintf("'%s' cannot be 'NULL'.\n", obj_name)
    stop(err)
  }

  if (!.is_anndata(obj)) {
    if (!.is_se_like(obj) || (!allow_se_like && !as_ace)) {
      if (error_on_fail) {
        err <- sprintf(
          "'%s' must be an AnnData object%s.\n",
          obj_name,
          ifelse(allow_se_like || as_ace, " or coercible container", "")
        )
        stop(err)
      }
      return(FALSE)
    }
    obj <- toAnnData(obj)
  }

  if (fix_dimnames) {
    obj <- .ensure_unique_dimnames(obj)
  }

  if (return_elem == TRUE) {
    return(obj)
  } else {
    return(TRUE)
  }
}

.validate_assay <- function(
    ace,
    assay_name = NULL,
    obj_name = "ace",
    matrix_type = "either",
    sparse_type = "CsparseMatrix",
    force_type = FALSE,
    error_on_fail = FALSE,
    return_elem = TRUE) {
  ace <- .validate_ace(
    ace,
    as_ace = FALSE,
    allow_se_like = TRUE,
    allow_null = FALSE,
    obj_name = obj_name,
    return_elem = TRUE,
    error_on_fail = TRUE
  )

  if (!is.null(assay_name) && !(assay_name %in% .get_layer_names(ace))) {
    if (error_on_fail) {
      stop(sprintf("'%s' is not a layer of '%s'.\n", assay_name, obj_name))
    }
    return(FALSE)
  }
  x <- .get_layer_matrix(ace, layer = assay_name, transpose = TRUE)

  if (force_type == TRUE) {
    x <- .validate_matrix(
      x = x,
      var_name = ifelse(
        is.null(assay_name),
        sprintf("%s$X", obj_name),
        sprintf("%s$layers[['%s']]", obj_name, assay_name)
      ),
      matrix_type = matrix_type,
      sparse_type = sparse_type,
      force_type = force_type,
      return_elem = TRUE
    )
  }

  if (return_elem == TRUE) {
    return(x)
  } else {
    return(TRUE)
  }
}

.validate_map <- function(
    ace,
    map_slot = NULL,
    obj_name = "ace",
    matrix_type = "either",
    force_type = FALSE,
    row = FALSE,
    return_elem = TRUE) {
  ace <- .validate_ace(
    ace,
    allow_null = FALSE,
    obj_name = obj_name,
    return_elem = TRUE,
    error_on_fail = TRUE
  )

  if (row == TRUE) {
    if (!(map_slot %in% names(rowMaps(ace)))) {
      stop(sprintf("'%s' is not an attribute of 'rowMaps(%s)'.\n", map_slot, obj_name))
    }
    x <- rowMaps(ace)[[map_slot]]
  } else {
    if (!(map_slot %in% names(colMaps(ace)))) {
      stop(sprintf("'%s' is not an attribute of 'colMaps(%s)'.\n", map_slot, obj_name))
    }
    x <- colMaps(ace)[[map_slot]]
  }

  if (force_type == TRUE) {
    x <- .validate_matrix(
      x = x,
      var_name = ifelse(
        row,
        sprintf("rowMaps(%s)$%s", obj_name, map_slot),
        sprintf("colMaps(%s)$%s", obj_name, map_slot)
      ),
      matrix_type = matrix_type,
      force_type = force_type,
      return_elem = TRUE
    )
  }

  if (return_elem == TRUE) {
    return(x)
  } else {
    return(NULL)
  }
}


.validate_net <- function(
    ace,
    net_slot = NULL,
    obj_name = "ace",
    matrix_type = "sparse",
    sparse_type = c("dMatrix", "CsparseMatrix"),
    force_type = FALSE,
    row = FALSE,
    return_elem = TRUE) {
  ace <- .validate_ace(
    ace,
    allow_null = FALSE,
    obj_name = obj_name,
    return_elem = TRUE,
    error_on_fail = TRUE
  )

  if (row == TRUE) {
    if (!(net_slot %in% names(rowNets(ace)))) {
      stop(sprintf("'%s' is not an attribute of 'rowNets(%s)'.\n", net_slot, obj_name))
    }
    x <- rowNets(ace)[[net_slot]]
  } else {
    if (!(net_slot %in% names(colNets(ace)))) {
      stop(sprintf("'%s' is not an attribute of 'colNets(%s)'.\n", net_slot, obj_name))
    }
    x <- colNets(ace)[[net_slot]]
  }

  if (force_type == TRUE) {
    x <- .validate_matrix(
      x = x,
      var_name = ifelse(
        row,
        sprintf("rowNets(%s)$%s", obj_name, net_slot),
        sprintf("colNets(%s)$%s", obj_name, net_slot)
      ),
      matrix_type = matrix_type,
      sparse_type = sparse_type,
      force_type = force_type,
      return_elem = TRUE
    )
  }

  if (return_elem == TRUE) {
    return(x)
  } else {
    return(NULL)
  }
}

# To be deprecated in favor of `.validate_vector_attr()`
.validate_attr <- function(
    obj,
    attr,
    return_type = "data",
    obj_name = "obj",
    attr_name = "attr",
    dim = 2,
    match_row = FALSE,
    return_elem = TRUE) {
  if (is.null(obj)) {
    err <- sprintf("'%s' cannot be 'NULL'.\n", obj_name)
    stop(err)
  }

  if (is.null(attr)) {
    err <- sprintf("'%s' cannot be 'NULL'.\n", attr_name)
    stop(err)
  }

  if (.is_se_like(obj)) {
    attr <- .validate_vector_attr(
      obj = obj,
      attr = attr,
      return_type = return_type,
      dim = dim,
      attr_name = attr_name,
      obj_name = obj_name,
      return_elem = TRUE
    )
  } else {
    attr <- c(attr)

    if (match_row == TRUE) {
      if (!(length(attr) == nrow(obj))) {
        err <- sprintf("'length(%s)' must equal 'NROW(%s)'.\n", attr_name, obj_name)
        stop(err)
      }
    } else {
      if (!(length(attr) == ncol(obj))) {
        err <- sprintf("'length(%s)' must equal 'NCOL(%s)'.\n", attr_name, obj_name)
        stop(err)
      }
    }
  }

  if (return_elem == TRUE) {
    return(attr)
  } else {
    return(NULL)
  }
}

# Replaces `.validate_attr()` and `ACTIONetExperiment:::.get.data.or.split()`
.validate_vector_attr <- function(
    obj,
    attr,
    groups_use = NULL,
    return_type = c("data", "split"),
    dim = 2,
    keep_factor = FALSE,
    attr_name = "attr",
    obj_name = "obj",
    return_elem = TRUE) {
  return_type <- match.arg(return_type, several.ok = TRUE)[1]


  if (length(attr) == 1) {
    if (!.is_se_like(obj)) {
      stop(sprintf("length(%s) must be %d for given object type", attr_name, .actionet_dim(obj)[dim]))
    }
    data_vec <- switch(dim,
      .get_feature_data(obj)[[attr]],
      .get_obs_data(obj)[[attr]]
    )
  } else {
    if (length(attr) != .actionet_dim(obj)[dim]) {
      err <- sprintf(
        "length(%s) (%d) does not match %s(%s) (%d)",
        attr_name,
        length(attr),
        ifelse(dim == 1, "NROW", "NCOL"),
        obj_name,
        .actionet_dim(obj)[dim]
      )
      stop(err)
    }
    data_vec <- attr
  }
  if (is.null(data_vec)) {
    err <- sprintf("Invalid attribute (%s)", attr_name)
    stop(err)
  } else {
    if (is.factor(data_vec)) {
      if (keep_factor) {
        data_vec <- droplevels(data_vec)
      } else {
        data_vec <- as.character(data_vec)
      }
    }
  }
  idx <- seq_len(.actionet_dim(obj)[dim])
  if (!is.null(groups_use)) {
    na_mask <- (!data_vec %in% groups_use)
    data_vec[na_mask] <- NA
    split_vec <- data_vec[!na_mask]
    if (length(split_vec) == 0) {
      stop("Invalid 'groups_use'")
    }

    idx_list <- split(idx[!na_mask], f = split_vec, drop = TRUE)
  } else {
    idx_list <- split(idx, f = data_vec)
  }

  if (return_type == "data") {
    out <- data_vec
  } else if (return_type == "split") {
    out <- idx_list
  }

  if (return_elem == TRUE) {
    return(out)
  } else {
    return(TRUE)
  }
}


.validate_matrix <- function(
    x,
    var_name = "x",
    matrix_type = c("either", "sparse", "dense"),
    sparse_type = c("dMatrix", "CsparseMatrix"),
    force_type = FALSE,
    return_elem = TRUE) {
  # Select first default if default vector passed by calling function.
  matrix_type <- match.arg(matrix_type, several.ok = TRUE)[1]
  sparse_type <- match.arg(sparse_type, several.ok = TRUE)[1]

  if (matrix_type == "sparse") {
    if (.is_sparse_matrix(x)) {
      if (!is(x, sparse_type)) {
        x <- as(x, sparse_type)
      }
    } else if (force_type == TRUE) {
      x <- as(x, sparse_type)
    } else {
      err <- sprintf("'%s' must be 'sparseMatrix'.\n", var_name)
      stop(err)
    }
  } else if (matrix_type == "dense") {
    if (!is.matrix(x)) {
      if (force_type == TRUE) {
        x <- as.matrix(x)
      } else {
        err <- sprintf("'%s' must be 'matrix'.\n", var_name)
        stop(err)
      }
    }
  } else {
    if (.is_sparse_matrix(x)) {
      if (!is(x, sparse_type)) {
        x <- as(x, sparse_type)
      }
    } else if (!is.matrix(x)) {
      err <- sprintf("'%s' must be 'matrix' or 'sparseMatrix'.\n", var_name)
      stop(err)
    }
  }

  if (return_elem == TRUE) {
    return(x)
  } else {
    return(NULL)
  }
}


.ace_or_assay <- function(
    obj,
    assay_name = NULL,
    allow_se_like = FALSE,
    matrix_type = c("either", "sparse", "dense"),
    force_type = FALSE,
    sparse_type = c("dMatrix", "CsparseMatrix"),
    obj_name = NULL,
    return_elem = TRUE) {
  if (.is_se_like(obj)) {
    x <- .validate_assay(
      ace = obj,
      assay_name = assay_name,
      obj_name = ifelse(is.null(obj_name), "ace", obj_name),
      matrix_type = matrix_type,
      sparse_type = sparse_type,
      force_type = force_type,
      error_on_fail = TRUE,
      return_elem = return_elem
    )
  } else {
    x <- .validate_matrix(
      x = obj,
      var_name = ifelse(is.null(obj_name), "obj", obj_name),
      matrix_type = matrix_type,
      sparse_type = sparse_type,
      force_type = force_type,
      return_elem = return_elem
    )
  }

  return(x)
}


.ace_or_map <- function(
    obj,
    map_slot = NULL,
    matrix_type = c("either", "sparse", "dense"),
    force_type = FALSE,
    obj_name = NULL,
    row = FALSE,
    transpose_map = FALSE,
    return_elem = TRUE) {
  if (.is_se_like(obj)) {
    x <- .validate_map(
      ace = obj,
      map_slot = map_slot,
      obj_name = ifelse(is.null(obj_name), "ace", obj_name),
      matrix_type = matrix_type,
      force_type = force_type,
      row = row,
      return_elem = return_elem
    )
    if (transpose_map) {
      x <- Matrix::t(x)
    }
  } else {
    x <- .validate_matrix(
      x = obj,
      var_name = ifelse(is.null(obj_name), "obj", obj_name),
      matrix_type = matrix_type,
      force_type = force_type,
      return_elem = return_elem
    )
  }

  return(x)
}


.ace_or_net <- function(
    obj,
    net_slot = NULL,
    matrix_type = c("either", "sparse", "dense"),
    sparse_type = c("CsparseMatrix", "dMatrix"),
    force_type = FALSE,
    obj_name = NULL,
    row = FALSE,
    return_elem = TRUE) {
  if (.is_se_like(obj)) {
    x <- .validate_net(
      ace = obj,
      net_slot = net_slot,
      obj_name = ifelse(is.null(obj_name), "ace", obj_name),
      matrix_type = matrix_type,
      sparse_type = sparse_type,
      force_type = force_type,
      row = row,
      return_elem = return_elem
    )
  } else {
    x <- .validate_matrix(
      x = obj,
      var_name = ifelse(is.null(obj_name), "obj", obj_name),
      matrix_type = matrix_type,
      sparse_type = sparse_type,
      force_type = force_type,
      return_elem = return_elem
    )
  }

  return(x)
}

.warn_dgt <- function(X) {
  if (any(class(X) %in% "dgTMatrix")) {
    msg <- sprintf("Computing on sparse matrix type 'dgTMatrix' will be very slow\nRecommend casting to 'dgCMatrix'")
    warning(msg)
  }
}
