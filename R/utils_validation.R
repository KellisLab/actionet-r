.validate_ace <- function(
    obj,
    as_ace = FALSE,
    allow_se_like = FALSE,
    allow_null = FALSE,
    obj_name = "ace",
    return_elem = FALSE,
    error_on_fail = TRUE) {
  if (is.null(obj) && !allow_null) {
    err <- sprintf("'%s' cannot be 'NULL'.\n", obj_name)
    stop(err)
  }

  if (any(class(obj) != "ACTIONetExperiment")) {
    if (allow_se_like) {
      if (!class(obj) %in% c("SummarizedExperiment", "RangedSummarizedExperiment", "SingleCellExperiment")) {
        if (error_on_fail) {
          err <- sprintf("'%s' must be 'ACTIONetExperiment' or inherit from 'SummarizedExperiment'.\n", obj_name)
          stop(err)
        } else {
          return(FALSE)
        }
      }
    } else {
      if (error_on_fail) {
        err <- sprintf("'%s' must be 'ACTIONetExperiment'.\n", obj_name)
        stop(err)
      } else {
        return(FALSE)
      }
    }
    if (as_ace) {
      obj <- as(obj, "ACTIONetExperiment")
    }
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
    return_elem = TRUE) {
  .validate_ace(
    ace,
    as_ace = FALSE,
    allow_se_like = TRUE,
    allow_null = FALSE,
    obj_name = obj_name,
    return_elem = FALSE,
    error_on_fail = TRUE
  )

  if (!(assay_name %in% names(assays(ace)))) {
    err <- sprintf("'%s' is not an assay of '%s'.\n", assay_name, obj_name)
    stop(err)
  }
  x <- SummarizedExperiment::assays(ace)[[assay_name]]

  if (force_type == TRUE) {
    x <- .validate_matrix(
      x = x,
      var_name = sprintf("assays(%s)$%s", obj_name, assay_name),
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

.validate_map <- function(
    ace,
    map_slot = NULL,
    obj_name = "ace",
    matrix_type = "either",
    force_type = FALSE,
    row = FALSE,
    return_elem = TRUE) {
  .validate_ace(
    ace,
    allow_null = FALSE,
    obj_name = obj_name,
    return_elem = FALSE,
    error_on_fail = TRUE
  )

  if (row == TRUE) {
    if (!(map_slot %in% names(rowMaps(ace)))) {
      err <- sprintf("'%s' is not an attribute of 'rowMaps(%s)'.\n", map_slot, obj_name)
      stop(err)
    }
    x <- rowMaps(ace)[[map_slot]]
  } else {
    if (!(map_slot %in% names(colMaps(ace)))) {
      err <- sprintf("'%s' is not an attribute of 'colMaps(%s)'.\n", map_slot, obj_name)
      stop(err)
    }
    x <- colMaps(ace)[[map_slot]]
  }

  if (force_type == TRUE) {
    x <- .validate_matrix(
      x = x,
      var_name = elseif(row, sprintf("rowMaps(%s)$%s", obj_name, map_slot), sprintf("colMaps(%s)$%s", obj_name, map_slot)),
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
  .validate_ace(
    ace,
    allow_null = FALSE,
    obj_name = obj_name,
    return_elem = FALSE,
    error_on_fail = TRUE
  )

  if (row == TRUE) {
    if (!(net_slot %in% names(rowNets(ace)))) {
      err <- sprintf("'%s' is not an attribute of 'rowNets(%s)'.\n", net_slot, obj_name)
      stop(err)
    }
    x <- rowNets(ace)[[net_slot]]
  } else {
    if (!(net_slot %in% names(colNets(ace)))) {
      err <- sprintf("'%s' is not an attribute of 'colNets(%s)'.\n", net_slot, obj_name)
      stop(err)
    }
    x <- colNets(ace)[[net_slot]]
  }

  if (force_type == TRUE) {
    x <- .validate_matrix(
      x = x,
      var_name = elseif(row, sprintf("rowNets(%s)$%s", obj_name, map_slot), sprintf("colNets(%s)$%s", obj_name, map_slot)),
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

.validate_attr <- function(
    obj,
    attr,
    return_type = "data",
    obj_name = "obj",
    attr_name = "attr",
    split_d = 2,
    match_row = FALSE,
    return_elem = TRUE) {
  if (is.null(obj)) {
    err <- sprintf("'%s' cannot be 'NULL'.\n", var_name)
    stop(err)
  }

  if (is.null(attr)) {
    err <- sprintf("'%s' cannot be 'NULL'.\n", attr_name)
    stop(err)
  }

  if (is(obj, "ACTIONetExperiment")) {
    attr <- ACTIONetExperiment::get.data.or.split(obj, attr = attr, to_return = return_type, d = split_d)
  } else {
    attr <- c(attr)

    if (match_row == TRUE) {
      if (!(length(attr) == NROW(obj))) {
        err <- sprintf("'length(%s)' must equal 'NROW(%s)'.\n", attr_name, obj_name)
        stop(err)
      }
    } else {
      if (!(length(attr) == NCOL(obj))) {
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
    if (ACTIONetExperiment:::is.sparseMatrix(x)) {
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
    if (ACTIONetExperiment:::is.sparseMatrix(x)) {
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
    obj_name = NULL,
    return_elem = TRUE) {
  if (is(obj, "ACTIONetExperiment")) {
    x <- .validate_assay(
      ace = obj,
      assay_name = assay_name,
      obj_name = ifelse(is.null(obj_name), "ace", obj_name),
      matrix_type = matrix_type,
      force_type = force_type,
      return_elem = return_elem
    )
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


.ace_or_map <- function(
    obj,
    map_slot = NULL,
    matrix_type = c("either", "sparse", "dense"),
    force_type = FALSE,
    obj_name = NULL,
    row = FALSE,
    transpose_map = FALSE,
    return_elem = TRUE) {
  if (is(obj, "ACTIONetExperiment")) {
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
  if (is(obj, "ACTIONetExperiment")) {
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
