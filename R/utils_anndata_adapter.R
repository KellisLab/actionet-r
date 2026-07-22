.is_anndata <- function(x) {
  inherits(x, "AbstractAnnData")
}

.is_backed_anndata <- function(x) {
  .is_anndata(x) && any(grepl("HDF5AnnData", class(x), fixed = TRUE))
}

.abort_if_backed_anndata <- function(x, arg = "adata") {
  if (.is_backed_anndata(x)) {
    stop(sprintf(
      "'%s' is a backed AnnData object. Backed object support is not implemented in actionet-r yet; materialize it in memory first.",
      arg
    ))
  }
}

.default_feature_names <- function(d) {
  paste0("feat_", seq_len(d))
}

.default_cell_names <- function(d) {
  paste0("sam_", seq_len(d))
}

.is_sparse_matrix <- function(x) {
  any(is(x) == "sparseMatrix")
}

.fast_col_means <- function(mat) {
  Matrix::colSums(mat) / nrow(mat)
}

.fast_row_means <- function(mat) {
  Matrix::rowSums(mat) / ncol(mat)
}

.tscalet <- function(x, center = TRUE, scale = TRUE) {
  t(scale(t(x), center = center, scale = scale))
}

.check_and_load_package <- function(pkg_names) {
  for (pk in pkg_names) {
    if (!requireNamespace(pk, quietly = TRUE)) {
      stop(sprintf("Package '%s' is not installed.\n", pk))
    }
  }
  invisible(TRUE)
}

.as_plain_list <- function(x) {
  if (is.null(x)) {
    return(list())
  }
  if (is.list(x) && !methods::is(x, "List")) {
    return(x)
  }
  if (methods::is(x, "List")) {
    return(as.list(x))
  }
  as.list(x)
}

.canonical_key <- function(key, alias_map) {
  if (length(key) == 0 || is.na(key) || is.null(key)) {
    return(key)
  }
  if (key %in% names(alias_map)) {
    alias_map[[key]]
  } else {
    key
  }
}

.canonicalize_named_list <- function(x, alias_map) {
  x <- .as_plain_list(x)
  if (length(x) == 0) {
    return(x)
  }

  out <- list()
  x_names <- names(x)
  if (is.null(x_names)) {
    return(x)
  }

  for (i in seq_along(x)) {
    key <- .canonical_key(x_names[[i]], alias_map)
    out[[key]] <- x[[i]]
  }

  out
}

.expand_aliases <- function(x, alias_map) {
  x <- .canonicalize_named_list(x, alias_map)
  if (length(x) == 0) {
    return(x)
  }

  out <- x
  reverse_map <- split(names(alias_map), alias_map)
  for (canonical in names(reverse_map)) {
    if (canonical %in% names(x)) {
      for (legacy in reverse_map[[canonical]]) {
        if (!(legacy %in% names(out))) {
          out[[legacy]] <- x[[canonical]]
        }
      }
    }
  }

  out
}

.ACTIONET_OBSM_ALIASES <- c(
  ACTION = "action",
  ACTIONet2D = "umap_2d_actionet",
  ACTIONet3D = "umap_3d_actionet",
  denovo_color = "colors_actionet"
)

.ACTIONET_VARM_ALIASES <- c(
  arch_feat_spec = "archetype_feat_specificity_upper"
)

.ACTIONET_OBSP_ALIASES <- c(
  ACTIONet = "actionet"
)

.ACTIONET_VARP_ALIASES <- character(0)

.ACTIONET_UNS_ALIASES <- c(
  ACTION_sigma = "action_params",
  action_sigma = "action_params"
)

.obs_names <- function(adata) {
  rownames(adata)
}

.var_names <- function(adata) {
  colnames(adata)
}

.set_obs_names <- function(adata, value) {
  rownames(adata) <- value
  adata
}

.set_var_names <- function(adata, value) {
  colnames(adata) <- value
  adata
}

.n_obs <- function(adata) {
  nrow(adata)
}

.n_vars <- function(adata) {
  ncol(adata)
}

.actionet_nrow <- function(obj) {
  if (.is_anndata(obj)) {
    return(.n_obs(obj))  # cells
  }
  nrow(obj)
}

.actionet_ncol <- function(obj) {
  if (.is_anndata(obj)) {
    return(.n_vars(obj))  # genes
  }
  ncol(obj)
}

.actionet_dim <- function(obj) {
  c(.actionet_nrow(obj), .actionet_ncol(obj))
}

.actionet_rownames <- function(obj) {
  if (.is_anndata(obj)) {
    return(.var_names(obj))  # genes (features) — stays in varm row-space
  }
  rownames(obj)
}

.actionet_colnames <- function(obj) {
  if (.is_anndata(obj)) {
    return(.obs_names(obj))  # cells (observations)
  }
  colnames(obj)
}

.set_actionet_rownames <- function(obj, value) {
  if (.is_anndata(obj)) {
    return(.set_var_names(obj, value))
  }
  rownames(obj) <- value
  obj
}

.set_actionet_colnames <- function(obj, value) {
  if (.is_anndata(obj)) {
    return(.set_obs_names(obj, value))
  }
  colnames(obj) <- value
  obj
}

#' Deduplicate obs/var names.
#'
#' @param x Character vector of names.
#' @param kind One of "obs" or "var"; used only for the warning message.
#' @param warn Whether to emit a single count-only warning when duplicates are
#'   present. Never enumerates the affected names since obs_names collisions on
#'   scRNA-seq data routinely exceed thousands.
#'
#' @return `make.unique(as.character(x), sep = "_")`.
#' @keywords internal
#' @noRd
.dedup_names <- function(x, kind = c("obs", "var"), warn = TRUE) {
  kind <- match.arg(kind)
  x <- as.character(x)
  dupes <- duplicated(x)
  n_dupes <- sum(dupes)
  if (n_dupes > 0L && warn) {
    warning(sprintf(
      "Duplicated %s_names detected (%d duplicates); deduplicated via make.unique(sep = '_').",
      kind,
      n_dupes
    ), call. = FALSE)
  }
  make.unique(x, sep = "_")
}

.ensure_unique_dimnames <- function(adata) {
  adata <- .as_inmemory_anndata(adata, arg = "adata")

  obs_names <- .obs_names(adata)
  var_names <- .var_names(adata)

  if (is.null(obs_names) || anyNA(obs_names) || any(obs_names == "")) {
    obs_names <- .default_cell_names(.n_obs(adata))
  } else {
    obs_names <- .dedup_names(obs_names, kind = "obs")
  }

  if (is.null(var_names) || anyNA(var_names) || any(var_names == "")) {
    var_names <- .default_feature_names(.n_vars(adata))
  } else {
    var_names <- .dedup_names(var_names, kind = "var")
  }

  adata <- .set_obs_names(adata, obs_names)
  adata <- .set_var_names(adata, var_names)

  adata
}

.as_inmemory_anndata <- function(x, arg = "adata") {
  if (!.is_anndata(x)) {
    stop(sprintf("'%s' must be an AnnData object.", arg))
  }

  .abort_if_backed_anndata(x, arg = arg)

  if (inherits(x, "AnnDataView")) {
    x <- x$as_InMemoryAnnData()
  }

  x
}

.get_layer_names <- function(adata) {
  names(.as_plain_list(adata$layers))
}

.get_layer_matrix <- function(adata, layer = NULL, transpose = FALSE, allow_null = FALSE) {
  adata <- .as_inmemory_anndata(adata)

  mat <- if (is.null(layer)) {
    adata$X
  } else {
    adata$layers[[layer]]
  }

  if (is.null(mat)) {
    if (allow_null) {
      return(NULL)
    }

    if (is.null(layer)) {
      stop("AnnData.X is NULL.")
    }
    stop(sprintf("Layer '%s' does not exist.", layer))
  }

  if (transpose) {
    mat <- Matrix::t(mat)
  }

  mat
}

.set_layer_matrix <- function(adata, layer = NULL, value, transpose = FALSE) {
  adata <- .as_inmemory_anndata(adata)
  if (transpose) {
    value <- if (.is_sparse_matrix(value)) Matrix::t(value) else t(as.matrix(value))
  }

  if (is.null(layer)) {
    adata$X <- value
  } else {
    adata$layers[[layer]] <- value
  }

  adata
}

.get_row_data_df <- function(adata) {
  as.data.frame(adata$var, stringsAsFactors = FALSE)
}

.set_row_data_df <- function(adata, value) {
  adata <- .as_inmemory_anndata(adata)
  adata$var <- as.data.frame(value, stringsAsFactors = FALSE)
  adata
}

.get_feature_data <- function(obj) {
  if (.is_anndata(obj)) {
    return(.get_row_data_df(obj))
  }
  as.data.frame(SummarizedExperiment::rowData(obj), stringsAsFactors = FALSE)
}

.set_feature_data <- function(obj, value) {
  if (.is_anndata(obj)) {
    return(.set_row_data_df(obj, value))
  }
  SummarizedExperiment::rowData(obj) <- value
  obj
}

.get_col_data_df <- function(adata) {
  as.data.frame(adata$obs, stringsAsFactors = FALSE)
}

.set_col_data_df <- function(adata, value) {
  adata <- .as_inmemory_anndata(adata)
  adata$obs <- as.data.frame(value, stringsAsFactors = FALSE)
  adata
}

.get_obs_data <- function(obj) {
  if (.is_anndata(obj)) {
    return(.get_col_data_df(obj))
  }
  as.data.frame(SummarizedExperiment::colData(obj), stringsAsFactors = FALSE)
}

.set_obs_data <- function(obj, value) {
  if (.is_anndata(obj)) {
    return(.set_col_data_df(obj, value))
  }
  SummarizedExperiment::colData(obj) <- value
  obj
}

.get_uns <- function(adata) {
  .as_plain_list(adata$uns)
}

.set_uns <- function(adata, value) {
  adata <- .as_inmemory_anndata(adata)
  adata$uns <- .as_plain_list(value)
  adata
}

.canonical_obsm_key <- function(key) {
  .canonical_key(key, .ACTIONET_OBSM_ALIASES)
}

.canonical_varm_key <- function(key) {
  .canonical_key(key, .ACTIONET_VARM_ALIASES)
}

.canonical_obsp_key <- function(key) {
  .canonical_key(key, .ACTIONET_OBSP_ALIASES)
}

.canonical_varp_key <- function(key) {
  .canonical_key(key, .ACTIONET_VARP_ALIASES)
}

.get_collection <- function(adata, slot_name, alias_map) {
  adata <- .as_inmemory_anndata(adata)
  .expand_aliases(.as_plain_list(adata[[slot_name]]), alias_map)
}

.set_collection <- function(adata, slot_name, value, alias_map) {
  adata <- .as_inmemory_anndata(adata)
  adata[[slot_name]] <- .canonicalize_named_list(value, alias_map)
  adata
}

.get_map_types <- function(adata, key_name, alias_map) {
  uns <- .get_uns(adata)
  types <- uns[[key_name]]
  .expand_aliases(.as_plain_list(types), alias_map)
}

.set_map_types <- function(adata, key_name, value, alias_map) {
  uns <- .get_uns(adata)
  uns[[key_name]] <- .canonicalize_named_list(value, alias_map)
  .set_uns(adata, uns)
}

colMaps <- function(object, all = TRUE) {
  if (.is_anndata(object)) {
    return(.get_collection(object, "obsm", .ACTIONET_OBSM_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    return(ACTIONetExperiment::colMaps(object, all = all))
  }
  stop("'object' does not support 'colMaps()'.")
}

`colMaps<-` <- function(object, value) {
  if (.is_anndata(object)) {
    return(.set_collection(object, "obsm", value, .ACTIONET_OBSM_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    ACTIONetExperiment::colMaps(object) <- value
    return(object)
  }
  stop("'object' does not support 'colMaps<-'.")
}

rowMaps <- function(object, all = TRUE) {
  if (.is_anndata(object)) {
    return(.get_collection(object, "varm", .ACTIONET_VARM_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    return(ACTIONetExperiment::rowMaps(object, all = all))
  }
  stop("'object' does not support 'rowMaps()'.")
}

`rowMaps<-` <- function(object, value) {
  if (.is_anndata(object)) {
    return(.set_collection(object, "varm", value, .ACTIONET_VARM_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    ACTIONetExperiment::rowMaps(object) <- value
    return(object)
  }
  stop("'object' does not support 'rowMaps<-'.")
}

colNets <- function(object, ...) {
  if (.is_anndata(object)) {
    return(.get_collection(object, "obsp", .ACTIONET_OBSP_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    return(ACTIONetExperiment::colNets(object, ...))
  }
  stop("'object' does not support 'colNets()'.")
}

`colNets<-` <- function(object, value) {
  if (.is_anndata(object)) {
    return(.set_collection(object, "obsp", value, .ACTIONET_OBSP_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    ACTIONetExperiment::colNets(object) <- value
    return(object)
  }
  stop("'object' does not support 'colNets<-'.")
}

rowNets <- function(object, ...) {
  if (.is_anndata(object)) {
    return(.get_collection(object, "varp", .ACTIONET_VARP_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    return(ACTIONetExperiment::rowNets(object, ...))
  }
  stop("'object' does not support 'rowNets()'.")
}

`rowNets<-` <- function(object, value) {
  if (.is_anndata(object)) {
    return(.set_collection(object, "varp", value, .ACTIONET_VARP_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    ACTIONetExperiment::rowNets(object) <- value
    return(object)
  }
  stop("'object' does not support 'rowNets<-'.")
}

colMapTypes <- function(object, all = TRUE) {
  if (.is_anndata(object)) {
    return(.get_map_types(object, "actionet_colMapTypes", .ACTIONET_OBSM_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    return(ACTIONetExperiment::colMapTypes(object, all = all))
  }
  stop("'object' does not support 'colMapTypes()'.")
}

`colMapTypes<-` <- function(object, value) {
  if (.is_anndata(object)) {
    return(.set_map_types(object, "actionet_colMapTypes", value, .ACTIONET_OBSM_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    ACTIONetExperiment::colMapTypes(object) <- value
    return(object)
  }
  stop("'object' does not support 'colMapTypes<-'.")
}

rowMapTypes <- function(object, all = TRUE) {
  if (.is_anndata(object)) {
    return(.get_map_types(object, "actionet_rowMapTypes", .ACTIONET_VARM_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    return(ACTIONetExperiment::rowMapTypes(object, all = all))
  }
  stop("'object' does not support 'rowMapTypes()'.")
}

`rowMapTypes<-` <- function(object, value) {
  if (.is_anndata(object)) {
    return(.set_map_types(object, "actionet_rowMapTypes", value, .ACTIONET_VARM_ALIASES))
  }
  if (inherits(object, "ACTIONetExperiment") && requireNamespace("ACTIONetExperiment", quietly = TRUE)) {
    ACTIONetExperiment::rowMapTypes(object) <- value
    return(object)
  }
  stop("'object' does not support 'rowMapTypes<-'.")
}

.transpose_matrix <- function(mat) {
  if (is.null(mat)) {
    return(NULL)
  }
  if (.is_sparse_matrix(mat)) {
    Matrix::t(mat)
  } else if (is.matrix(mat)) {
    t(mat)
  } else {
    Matrix::t(mat)
  }
}

.resolve_x_layer <- function(assay_names, x_layer) {
  if (length(assay_names) == 0) {
    return(NULL)
  }
  if (is.null(x_layer)) {
    return(assay_names[[1]])
  }
  if (!(x_layer %in% assay_names)) {
    stop(sprintf(
      "'x_layer' = '%s' not found among available assays: %s",
      x_layer,
      paste(assay_names, collapse = ", ")
    ), call. = FALSE)
  }
  x_layer
}

.fold_sigma_into_uns <- function(uns, ace_meta) {
  if (length(ace_meta) == 0) {
    return(uns)
  }
  for (nm in names(ace_meta)) {
    val <- ace_meta[[nm]]
    if (grepl("_sigma$", nm)) {
      params_nm <- sub("_sigma$", "_params", nm)
      existing <- uns[[params_nm]]
      if (is.null(existing) || !is.list(existing)) {
        existing <- list()
      }
      existing[["sigma"]] <- val
      uns[[params_nm]] <- existing
      # Retain the flat key too for consumers reading metadata directly.
      uns[[nm]] <- val
    } else {
      uns[[nm]] <- val
    }
  }
  uns
}

.unfold_sigma_from_uns <- function(uns) {
  # Symmetric reverse of .fold_sigma_into_uns: for any 'foo_params' with a
  # $sigma child, emit 'foo_sigma' at the top level (in addition to
  # 'foo_params', which is retained). Internal book-keeping keys are dropped.
  if (length(uns) == 0) {
    return(list())
  }
  internal_keys <- c("X_name", "actionet_colMapTypes", "actionet_rowMapTypes")
  meta <- list()
  for (nm in names(uns)) {
    if (nm %in% internal_keys) next
    val <- uns[[nm]]
    meta[[nm]] <- val
    if (grepl("_params$", nm) && is.list(val) && !is.null(val[["sigma"]])) {
      sigma_key <- sub("_params$", "_sigma", nm)
      if (is.null(meta[[sigma_key]])) {
        meta[[sigma_key]] <- val[["sigma"]]
      }
    }
  }
  meta
}

.se_to_anndata_direct <- function(x, x_layer = NULL, ace = NULL) {
  assays_list <- as.list(SummarizedExperiment::assays(x))
  assay_names <- names(assays_list)
  if (is.null(assay_names) && length(assays_list) > 0) {
    assay_names <- paste0("assay", seq_along(assays_list))
    names(assays_list) <- assay_names
  }

  x_layer <- .resolve_x_layer(assay_names, x_layer)

  X <- if (!is.null(x_layer)) .transpose_matrix(assays_list[[x_layer]]) else NULL

  layers <- list()
  for (nm in assay_names) {
    if (identical(nm, x_layer)) next
    layers[[nm]] <- .transpose_matrix(assays_list[[nm]])
  }

  obs_df <- as.data.frame(
    SummarizedExperiment::colData(x),
    stringsAsFactors = FALSE,
    optional = TRUE
  )
  var_df <- as.data.frame(
    SummarizedExperiment::rowData(x),
    stringsAsFactors = FALSE,
    optional = TRUE
  )

  obs_names <- if (is.null(colnames(x))) .default_cell_names(ncol(x)) else .dedup_names(colnames(x), "obs")
  var_names <- if (is.null(rownames(x))) .default_feature_names(nrow(x)) else .dedup_names(rownames(x), "var")

  .apply_dimnames <- function(mat) {
    if (is.null(mat)) return(NULL)
    dimnames(mat) <- list(obs_names, var_names)
    mat
  }
  X <- .apply_dimnames(X)
  layers <- lapply(layers, .apply_dimnames)

  # anndataR::AnnData() rejects zero-column dataframes with row indexing; rebuild
  # obs/var with row.names anchored to dimnames.
  if (ncol(obs_df) == 0) {
    obs_df <- data.frame(row.names = obs_names)
  } else {
    rownames(obs_df) <- obs_names
  }
  if (ncol(var_df) == 0) {
    var_df <- data.frame(row.names = var_names)
  } else {
    rownames(var_df) <- var_names
  }

  adata <- anndataR::AnnData(
    X = X,
    layers = if (length(layers) > 0) layers else NULL,
    obs = obs_df,
    var = var_df
  )
  adata <- .as_inmemory_anndata(adata, arg = "x")

  # Metadata: for ACE, defer to .fold_sigma_into_uns below (which handles both
  # _sigma folding and plain passthrough). For bare SE/SCE, copy metadata()
  # directly.
  se_meta <- .as_plain_list(S4Vectors::metadata(x))
  uns <- .get_uns(adata)
  if (!inherits(x, "ACTIONetExperiment") && length(se_meta) > 0) {
    for (nm in names(se_meta)) {
      uns[[.canonical_key(nm, .ACTIONET_UNS_ALIASES)]] <- se_meta[[nm]]
    }
  }

  # SingleCellExperiment: preserve reducedDims -> obsm.
  if (inherits(x, "SingleCellExperiment") && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    rd <- SingleCellExperiment::reducedDims(x)
    if (length(rd) > 0) {
      obsm <- .as_plain_list(adata$obsm)
      for (nm in names(rd)) {
        obsm[[.canonical_obsm_key(nm)]] <- as.matrix(rd[[nm]])
      }
      adata$obsm <- obsm
    }
  }

  # ACE-specific slots (only when a live ACE was supplied).
  if (!is.null(ace) && inherits(ace, "ACTIONetExperiment")) {
    col_maps <- .canonicalize_named_list(
      ACTIONetExperiment::colMaps(ace, all = TRUE),
      .ACTIONET_OBSM_ALIASES
    )
    row_maps <- .canonicalize_named_list(
      ACTIONetExperiment::rowMaps(ace, all = TRUE),
      .ACTIONET_VARM_ALIASES
    )
    col_nets <- .canonicalize_named_list(
      ACTIONetExperiment::colNets(ace),
      .ACTIONET_OBSP_ALIASES
    )
    row_nets <- .canonicalize_named_list(
      ACTIONetExperiment::rowNets(ace),
      .ACTIONET_VARP_ALIASES
    )

    if (length(col_maps) > 0) {
      obsm <- .as_plain_list(adata$obsm)
      # ACE colMaps entries may be SummarizedExperiment-wrapped; unwrap to raw matrix.
      for (nm in names(col_maps)) {
        val <- col_maps[[nm]]
        if (inherits(val, "SummarizedExperiment")) {
          val <- SummarizedExperiment::assay(val, 1L)
        }
        obsm[[nm]] <- val
      }
      adata$obsm <- obsm
    }
    if (length(row_maps) > 0) {
      varm <- .as_plain_list(adata$varm)
      for (nm in names(row_maps)) {
        val <- row_maps[[nm]]
        if (inherits(val, "SummarizedExperiment")) {
          val <- SummarizedExperiment::assay(val, 1L)
        }
        varm[[nm]] <- val
      }
      adata$varm <- varm
    }
    if (length(col_nets) > 0) {
      adata$obsp <- c(.as_plain_list(adata$obsp), col_nets)
    }
    if (length(row_nets) > 0) {
      adata$varp <- c(.as_plain_list(adata$varp), row_nets)
    }

    uns <- .fold_sigma_into_uns(uns, .as_plain_list(S4Vectors::metadata(ace)))

    col_types <- .as_plain_list(ACTIONetExperiment::colMapTypes(ace, all = TRUE))
    row_types <- .as_plain_list(ACTIONetExperiment::rowMapTypes(ace, all = TRUE))
    if (length(col_types) > 0) {
      uns[["actionet_colMapTypes"]] <- .canonicalize_named_list(col_types, .ACTIONET_OBSM_ALIASES)
    }
    if (length(row_types) > 0) {
      uns[["actionet_rowMapTypes"]] <- .canonicalize_named_list(row_types, .ACTIONET_VARM_ALIASES)
    }
  }

  # Record which assay landed in .X so the reverse round-trip can restore the name.
  if (!is.null(x_layer)) {
    uns[["X_name"]] <- x_layer
  }

  adata$uns <- uns
  adata
}

#' Convert supported ACTIONet inputs to in-memory AnnData.
#'
#' Delegates transposition to explicit slot-by-slot copies (no
#' `anndataR::as_AnnData()` intermediate), so behavior does not depend on
#' `anndataR`'s S3 method table.
#'
#' @param x AnnData, matrix, sparse matrix, `SummarizedExperiment`,
#'   `SingleCellExperiment`, or `ACTIONetExperiment`. Bare matrix / sparse
#'   matrix inputs are treated as genes x cells (SE convention) and transposed
#'   to cells x genes for AnnData.
#' @param x_layer Optional character. When `x` carries multiple assays
#'   (`SummarizedExperiment` / `SingleCellExperiment` / `ACTIONetExperiment`),
#'   the assay named here is placed in `.X` and the remaining assays go into
#'   `.layers`. Defaults to the first assay. Ignored for AnnData and bare
#'   matrix inputs.
#'
#' @return An in-memory `anndataR::AnnData` object.
#' @export
toAnnData <- function(x, x_layer = NULL) {
  if (.is_anndata(x)) {
    x <- .as_inmemory_anndata(x, arg = "x")
    return(.ensure_unique_dimnames(x))
  }

  if (.is_sparse_matrix(x) || is.matrix(x)) {
    # Dedup dimnames BEFORE constructing AnnData; anndataR reads obs_names /
    # var_names off obs/var row.names, and rejects X whose dimnames disagree.
    obs_names <- if (is.null(colnames(x))) .default_cell_names(ncol(x)) else .dedup_names(colnames(x), "obs")
    var_names <- if (is.null(rownames(x))) .default_feature_names(nrow(x)) else .dedup_names(rownames(x), "var")

    x_mat <- x
    dimnames(x_mat) <- list(var_names, obs_names)
    X_t <- if (.is_sparse_matrix(x_mat)) Matrix::t(x_mat) else t(as.matrix(x_mat))
    obs_df <- data.frame(row.names = obs_names)
    var_df <- data.frame(row.names = var_names)

    adata <- anndataR::AnnData(
      X = X_t,
      obs = obs_df,
      var = var_df
    )
    return(.as_inmemory_anndata(adata, arg = "x"))
  }

  if (inherits(x, "ACTIONetExperiment")) {
    .check_and_load_package("ACTIONetExperiment")
    # ACTIONetExperiment::colMaps() / rowMaps() / colMapTypes() / rowMapTypes()
    # internally call unqualified `assays()` and `metadata()`; attach the
    # relevant Bioconductor packages for the duration of this call so those
    # generics resolve.
    for (pk in c("SummarizedExperiment", "S4Vectors")) {
      if (!paste0("package:", pk) %in% search()) {
        suppressPackageStartupMessages(attachNamespace(pk))
      }
    }
    return(.se_to_anndata_direct(x, x_layer = x_layer, ace = x))
  }

  if (inherits(x, "SingleCellExperiment") || inherits(x, "SummarizedExperiment")) {
    return(.se_to_anndata_direct(x, x_layer = x_layer, ace = NULL))
  }

  stop("'x' must be an AnnData, matrix, sparseMatrix, SummarizedExperiment, SingleCellExperiment, or ACTIONetExperiment object.")
}

.reverse_key <- function(key, alias_map) {
  reverse_map <- split(names(alias_map), alias_map)
  if (key %in% names(reverse_map)) {
    reverse_map[[key]][[1]]
  } else {
    key
  }
}

#' Convert AnnData to `ACTIONetExperiment`.
#'
#' Builds an `ACTIONetExperiment` directly from the AnnData slots (no
#' `as_SingleCellExperiment()` / `as.ACTIONetExperiment()` intermediate), so
#' `colMaps` / `rowMaps` / `colNets` / `rowNets` are populated exactly once
#' with a consistent element type.
#'
#' The assay name for `.X` is recovered from `adata$uns$X_name` if present,
#' otherwise defaults to `"X"`. Layers keep their original names.
#'
#' Metadata round-trips symmetrically: any `uns[[<foo>_params]]$sigma` is
#' restored as `metadata(ace)[[<foo>_sigma]]` (mirror of the forward
#' `_sigma`-into-`_params` folding done by [toAnnData()]).
#'
#' Requires the `ACTIONetExperiment` package (in `Suggests:`); errors with an
#' actionable message when it is not installed.
#'
#' @param adata AnnData or compatible input (see [toAnnData()]).
#'
#' @return An `ACTIONetExperiment` object with `rows = features (genes)`,
#'   `cols = samples (cells)`.
#' @export
toACTIONetExperiment <- function(adata) {
  .check_and_load_package("ACTIONetExperiment")

  adata <- toAnnData(adata)

  # Build assays: cells x genes -> genes x cells.
  uns <- .get_uns(adata)
  x_name <- if (!is.null(uns[["X_name"]]) && nzchar(uns[["X_name"]])) uns[["X_name"]] else "X"

  assays_list <- list()
  X <- adata$X
  if (!is.null(X)) {
    assays_list[[x_name]] <- .transpose_matrix(X)
  }
  for (nm in names(.as_plain_list(adata$layers))) {
    if (identical(nm, x_name)) next
    assays_list[[nm]] <- .transpose_matrix(adata$layers[[nm]])
  }

  obs_df <- .get_col_data_df(adata)
  var_df <- .get_row_data_df(adata)

  rownames_var <- .var_names(adata)
  rownames_obs <- .obs_names(adata)

  ace <- ACTIONetExperiment::ACTIONetExperiment(
    assays = assays_list,
    rowData = if (ncol(var_df) > 0) S4Vectors::DataFrame(var_df) else S4Vectors::DataFrame(row.names = rownames_var),
    colData = if (ncol(obs_df) > 0) S4Vectors::DataFrame(obs_df) else S4Vectors::DataFrame(row.names = rownames_obs)
  )

  rownames(ace) <- rownames_var
  colnames(ace) <- rownames_obs

  # obsm -> colMaps (raw matrix; canonical key reversed to legacy alias if needed)
  for (nm in names(.as_plain_list(adata$obsm))) {
    ACTIONetExperiment::colMaps(ace)[[.reverse_key(nm, .ACTIONET_OBSM_ALIASES)]] <- adata$obsm[[nm]]
  }
  for (nm in names(.as_plain_list(adata$varm))) {
    ACTIONetExperiment::rowMaps(ace)[[.reverse_key(nm, .ACTIONET_VARM_ALIASES)]] <- adata$varm[[nm]]
  }
  for (nm in names(.as_plain_list(adata$obsp))) {
    ACTIONetExperiment::colNets(ace)[[.reverse_key(nm, .ACTIONET_OBSP_ALIASES)]] <- adata$obsp[[nm]]
  }
  for (nm in names(.as_plain_list(adata$varp))) {
    ACTIONetExperiment::rowNets(ace)[[.reverse_key(nm, .ACTIONET_VARP_ALIASES)]] <- adata$varp[[nm]]
  }

  # Optional map-type sidecars.
  if (!is.null(uns[["actionet_colMapTypes"]])) {
    types <- .as_plain_list(uns[["actionet_colMapTypes"]])
    if (length(types) > 0) {
      ACTIONetExperiment::colMapTypes(ace) <- types
    }
  }
  if (!is.null(uns[["actionet_rowMapTypes"]])) {
    types <- .as_plain_list(uns[["actionet_rowMapTypes"]])
    if (length(types) > 0) {
      ACTIONetExperiment::rowMapTypes(ace) <- types
    }
  }

  # Metadata: symmetric _sigma reverse plus passthrough.
  meta <- .unfold_sigma_from_uns(uns)
  if (length(meta) > 0) {
    S4Vectors::metadata(ace) <- meta
  }

  ace
}

.subset_actionet_container <- function(adata, features = NULL, cells = NULL) {
  adata <- toAnnData(adata)

  if (is.null(features)) {
    features <- seq_len(.n_vars(adata))
  }
  if (is.null(cells)) {
    cells <- seq_len(.n_obs(adata))
  }

  out <- adata[cells, features]
  if (inherits(out, "AnnDataView")) {
    out <- out$as_InMemoryAnnData()
  }
  .ensure_unique_dimnames(out)
}

.resolve_layer_arg <- function(
    layer = NULL,
    assay_name = NULL,
    default = NULL,
    layer_missing = FALSE,
    assay_name_missing = FALSE) {
  if (!assay_name_missing && !is.null(assay_name)) {
    .Deprecated(msg = "'assay_name' is deprecated; use 'layer' instead.")
    if (is.null(layer)) {
      layer <- assay_name
    }
  }

  if (layer_missing && (assay_name_missing || is.null(assay_name))) {
    layer <- default
  }

  layer
}

.resolve_container_arg <- function(adata = NULL, ace = NULL, obj = NULL, default_obj = NULL) {
  out <- default_obj

  if (!is.null(obj)) {
    .Deprecated(msg = "'obj' is deprecated; use 'adata' instead.")
    out <- obj
  }
  if (!is.null(ace)) {
    .Deprecated(msg = "'ace' is deprecated; use 'adata' instead.")
    out <- ace
  }
  if (!is.null(adata)) {
    out <- adata
  }

  out
}
