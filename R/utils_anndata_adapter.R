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

.ensure_unique_dimnames <- function(adata) {
  adata <- .as_inmemory_anndata(adata, arg = "adata")

  obs_names <- .obs_names(adata)
  var_names <- .var_names(adata)

  if (is.null(obs_names) || anyNA(obs_names) || any(obs_names == "")) {
    obs_names <- .default_cell_names(.n_obs(adata))
  } else {
    obs_names <- make.unique(as.character(obs_names), sep = "_")
  }

  if (is.null(var_names) || anyNA(var_names) || any(var_names == "")) {
    var_names <- .default_feature_names(.n_vars(adata))
  } else {
    var_names <- make.unique(as.character(var_names), sep = "_")
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

.copy_ace_payload_to_anndata <- function(adata, ace) {
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
    adata$obsm <- c(.as_plain_list(adata$obsm), col_maps)
  }
  if (length(row_maps) > 0) {
    adata$varm <- c(.as_plain_list(adata$varm), row_maps)
  }
  if (length(col_nets) > 0) {
    adata$obsp <- c(.as_plain_list(adata$obsp), col_nets)
  }
  if (length(row_nets) > 0) {
    adata$varp <- c(.as_plain_list(adata$varp), row_nets)
  }

  uns <- .get_uns(adata)
  ace_meta <- .as_plain_list(S4Vectors::metadata(ace))
  if (length(ace_meta) > 0) {
    for (nm in names(ace_meta)) {
      if (grepl("_sigma$", nm)) {
        canonical_nm <- sub("_sigma$", "_params", .canonical_key(sub("_sigma$", "", nm), .ACTIONET_OBSM_ALIASES))
        uns[[canonical_nm]] <- c(.as_plain_list(uns[[canonical_nm]]), list(sigma = ace_meta[[nm]]))
      } else {
        uns[[.canonical_key(nm, .ACTIONET_UNS_ALIASES)]] <- ace_meta[[nm]]
      }
    }
  }
  adata$uns <- uns

  col_types <- .as_plain_list(ACTIONetExperiment::colMapTypes(ace, all = TRUE))
  row_types <- .as_plain_list(ACTIONetExperiment::rowMapTypes(ace, all = TRUE))
  if (length(col_types) > 0) {
    adata <- .set_map_types(adata, "actionet_colMapTypes", col_types, .ACTIONET_OBSM_ALIASES)
  }
  if (length(row_types) > 0) {
    adata <- .set_map_types(adata, "actionet_rowMapTypes", row_types, .ACTIONET_VARM_ALIASES)
  }

  adata
}

#' Convert supported ACTIONet inputs to in-memory AnnData.
#'
#' @param x AnnData, matrix, sparse matrix, `SummarizedExperiment`,
#'   `SingleCellExperiment`, or `ACTIONetExperiment`.
#'
#' @return An in-memory `anndataR::AnnData` object.
#' @export
toAnnData <- function(x) {
  if (.is_anndata(x)) {
    x <- .as_inmemory_anndata(x, arg = "x")
    return(.ensure_unique_dimnames(x))
  }

  if (.is_sparse_matrix(x) || is.matrix(x)) {
    adata <- anndataR::AnnData(
      # Assume legacy genes x cells orientation; transpose to cells x genes for AnnData
      X = if (.is_sparse_matrix(x)) Matrix::t(x) else t(as.matrix(x))
    )
    adata <- .ensure_unique_dimnames(adata)
    adata <- .set_obs_names(adata, if (is.null(colnames(x))) .default_cell_names(ncol(x)) else make.unique(colnames(x), sep = "_"))
    adata <- .set_var_names(adata, if (is.null(rownames(x))) .default_feature_names(nrow(x)) else make.unique(rownames(x), sep = "_"))
    return(adata)
  }

  if (inherits(x, "ACTIONetExperiment")) {
    .check_and_load_package("ACTIONetExperiment")
    sce <- as(x, "SingleCellExperiment")
    adata <- anndataR::as_AnnData(sce)
    adata <- .as_inmemory_anndata(adata, arg = "x")
    assay_names <- names(SummarizedExperiment::assays(x))
    if (length(assay_names) > 0) {
      adata$X <- adata$layers[[assay_names[[1]]]]
    }
    adata <- .copy_ace_payload_to_anndata(adata, x)
    return(.ensure_unique_dimnames(adata))
  }

  if (inherits(x, "SummarizedExperiment") || inherits(x, "SingleCellExperiment")) {
    adata <- anndataR::as_AnnData(x)
    adata <- .as_inmemory_anndata(adata, arg = "x")
    assay_names <- names(SummarizedExperiment::assays(x))
    if (length(assay_names) > 0) {
      adata$X <- adata$layers[[assay_names[[1]]]]
    }
    return(.ensure_unique_dimnames(adata))
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
#' @param adata AnnData or compatible input.
#'
#' @return An `ACTIONetExperiment` object.
#' @export
toACTIONetExperiment <- function(adata) {
  .check_and_load_package("ACTIONetExperiment")
  suppressPackageStartupMessages(
    require("ACTIONetExperiment", character.only = TRUE)
  )
  adata <- toAnnData(adata)

  sce <- adata$as_SingleCellExperiment()
  ace <- ACTIONetExperiment::as.ACTIONetExperiment(sce)

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

  uns <- .get_uns(adata)
  if ("action_params" %in% names(uns) && !is.null(uns[["action_params"]][["sigma"]])) {
    S4Vectors::metadata(ace)[["action_sigma"]] <- uns[["action_params"]][["sigma"]]
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
