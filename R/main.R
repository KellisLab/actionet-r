#' Run main ACTIONet pipeline
#'
#' @param adata AnnData or compatible container containing an appropriate reduction in `reduction_slot`.
#' @param k_min Minimum depth of decompositions. (default=2)
#' @param k_max Maximum depth of decompositions. (default=30)
#' @param layer Layer to be used. `NULL` uses `adata$X`. (default='logcounts')
#' @param reduction_slot Entry in `adata$obsm` containing the reduced kernel. (default='action')
#' @param net_slot_out Name of entry in `adata$obsp` to store the ACTIONet adjacency matrix. (default='actionet')
#' @param min_obs Minimum number of observations required to construct an archetype. (default=2)
#' @param max_it Maximum number of iterations for ACTION algorithm. (default=50)
#' @param spec_th Defines the stringency of pruning nonspecific archetypes.
#' The larger the value, the more archetypes will be filtered out. (default=-3)
#' @param network_metric Distance metric with which to compute cell-to-cell similarity during network construction. Options are 'jsd' (Jensen-Shannon divergence), L2-norm ('l2'), and inner product ('ip'). (default='jsd')
#' @param network_algorithm Algorithm to use for network construction. Options are k-nearest neighbors ('knn') and k*-nearest neighbors ('k*nn'). (default='k*nn')
#' @param network_density Density factor of ACTIONet graph. (default=1)
#' @param network_M HNSW graph connectivity parameter for network construction. (default=16)
#' @param network_ef_construction HNSW construction search breadth for network construction. For `network_algorithm="k*nn"`, the effective value is `max(network_ef_construction, kNN)`. (default=200)
#' @param network_ef HNSW query search breadth for network construction. For `network_algorithm="k*nn"`, the effective value is `max(network_ef, kNN)`. (default=200)
#' @param network_k Number of nearest neighbors for `network_algorithm="knn"`. (default=100)
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors. (default=TRUE)
#' @param layout_method Algorithm for computing plot layout. Options are UMAP ("umap") or t-UMAP ("tumap"). (default="umap")
#' @param layout_epochs Number of epochs for SGD algorithm. (default=250)
#' @param layout_parallel Run layout construction using multiple cores. May result in marginally different outputs across runs due to parallelization-induced randomization. (default=TRUE)
#' @param compute_specificity_parallel Run feature specificity enrichment using multiple cores. Setting this to `TRUE` on large datasets may cause an out of memory crash. (default=FALSE)
#' @param thread_no Number of parallel threads. (default=0)
#' @param seed Seed for random initialization. (default=0)
#' @param assay_name Deprecated alias for `layer`.
#' @param ace Deprecated alias for `adata`.
#'
#' @return An AnnData object populated with canonical ACTIONet outputs.
#'
#' @examples
#' adata <- runACTIONet(adata)
#' @export
runACTIONet <- function(adata = NULL,
                        k_min = 2,
                        k_max = 30,
                        layer = "logcounts",
                        reduction_slot = "action",
                        net_slot_out = "actionet",
                        min_obs = 2,
                        max_it = 50,
                        spec_th = -3,
                        network_metric = "jsd",
                        network_algorithm = "k*nn",
                        network_density = 1,
                        network_M = 16,
                        network_ef_construction = 200,
                        network_ef = 200,
                        network_k = 100,
                        mutual_edges_only = TRUE,
                        layout_method = c("umap", "tumap", "largevis"),
                        layout_epochs = 100,
                        layout_spread = 1.0,
                        layout_min_dist = 1.0,
                        layout_3d = TRUE,
                        layout_parallel = TRUE,
                        compute_specificity_parallel = FALSE,
                        thread_no = 0,
                        seed = 0,
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
  adata <- .validate_ace(adata, as_ace = TRUE, allow_se_like = TRUE, return_elem = TRUE)

  layout_method <- tolower(layout_method)
  layout_method <- match.arg(layout_method, several.ok = FALSE)

  .validate_assay(adata, assay_name = layer, return_elem = FALSE)
  .validate_map(adata, map_slot = reduction_slot, return_elem = FALSE)

  adata <- runACTION(
    adata = adata,
    k_min = k_min,
    k_max = k_max,
    reduction_slot = reduction_slot,
    prenormalize = TRUE,
    min_obs = min_obs,
    max_it = max_it,
    spec_th = spec_th,
    thread_no = thread_no
  )

  # Build ACTIONet
  adata <- buildNetwork(
    adata = adata,
    algorithm = network_algorithm,
    distance_metric = network_metric,
    density = network_density,
    thread_no = thread_no,
    M = network_M,
    ef_construction = network_ef_construction,
    ef = network_ef,
    mutual_edges_only = mutual_edges_only,
    k = network_k,
    map_slot = "H_stacked",
    net_slot_out = net_slot_out
  )

  adata <- networkCentrality(
    adata = adata,
    label_attr = "assigned_archetype",
    algorithm = "local_coreness",
    thread_no = thread_no,
    net_slot = net_slot_out,
    attr_out = "node_centrality"
  )

  # Smooth archetype footprints
  adata <- networkDiffusion(
    adata = adata,
    scores = colMaps(adata)[["H_merged"]],
    norm_method = "pagerank",
    thread_no = thread_no,
    net_slot = net_slot_out,
    map_slot_out = "archetype_footprint"
  )

  # Use archetypal reduction as initial coordinates for uwot
  # Need to reduce it to 3D coordinate space.
  # red.out <- runSVD(
  #   X = scale(colMaps(ace)[["archetype_footprint"]]),
  #   k = 3,
  #   seed = seed,
  #   verbose = FALSE
  # )
  # initial_coordinates <- scale(red.out$u)

  slot_layout <- sprintf("%s_%s", layout_method, net_slot_out)
  initial_coordinates <- colMaps(adata)[["archetype_footprint"]]

  layout_args <- list(
    method = layout_method,
    n_components = 2,
    min_dist = layout_min_dist,
    spread = layout_spread,
    n_epochs = layout_epochs,
    net_slot = net_slot_out,
    seed = seed,
    thread_no = ifelse(layout_parallel, thread_no, 1)
  )
  layout_args$map_slot_out <- sprintf("%s_%dd_%s", layout_args$method, layout_args$n_components, layout_args$net_slot)

  adata <- do.call(
    layoutNetwork,
    c(
      list(
        adata = adata,
        initial_coordinates = scale(initial_coordinates)
      ),
      layout_args
    )
  )

  if (layout_3d) {
    # Warm-start 3D layout using 2D embedding
    initial_coordinates <- cbind(
      colMaps(adata)[[layout_args$map_slot_out]],
      initial_coordinates[, 3]
    )

    layout_args$n_components <- 3
    layout_args$n_epochs <- layout_epochs / 2
    layout_args$map_slot_out <- sprintf("%s_%dd_%s", layout_args$method, layout_args$n_components, layout_args$net_slot)

    adata <- do.call(
      layoutNetwork,
      c(
        list(
          adata = adata,
          initial_coordinates = scale(initial_coordinates)
        ),
        layout_args
      )
    )

    adata <- computeNodeColors(
      adata = adata,
      embedding_slot = layout_args$map_slot_out,
      color_slot_out = sprintf("colors_%s", layout_args$net_slot),
      thread_no = thread_no
    )
  }

  # Compute gene specificity for each archetype
  adata <- archetypeFeatureSpecificity(
    adata = adata,
    layer = layer,
    map_slot = "archetype_footprint",
    thread_no = ifelse(compute_specificity_parallel, thread_no, 1),
    return_raw = FALSE
  )

  return(adata)
}

#' @export
runACTION <- function(
    adata = NULL,
    k_min = 2,
    k_max = 30,
    reduction_slot = "action",
    prenormalize = TRUE,
    min_obs = 2,
    max_it = 50,
    tol = 1e-100,
    spec_th = -3,
    thread_no = 0,
    merged_suffix = "merged",
    archetype_slot_out = "assigned_archetype",
    ace = NULL) {
  adata <- .resolve_container_arg(adata = adata, ace = ace)
  adata <- .validate_ace(
    adata,
    as_ace = TRUE,
    allow_se_like = TRUE,
    return_elem = TRUE
  )

  S_r <- .validate_map(
    ace = adata,
    map_slot = reduction_slot,
    matrix_type = "dense",
    force_type = TRUE,
    return_elem = TRUE
  )

  if (prenormalize) {
    S_r <- normalize.matrix(
      S_r,
      dim = 1, # cells are rows in obsm (cells x k); normalize each cell vector
      scale_param = NULL,
      trans_func = NULL
    )
  }

  out <- C_runACTION(
    S_r = S_r,
    k_min = k_min,
    k_max = k_max,
    max_it = max_it,
    tol = tol,
    spec_th = spec_th,
    min_obs = min_obs,
    thread_no = thread_no
  )

  colMaps(adata)[["H_stacked"]] <- as(out$H_stacked, "sparseMatrix")
  colMapTypes(adata)[["H_stacked"]] <- "internal"

  colMaps(adata)[["C_stacked"]] <- as(out$C_stacked, "sparseMatrix")
  colMapTypes(adata)[["C_stacked"]] <- "internal"

  colMaps(adata)[[sprintf("H_%s", merged_suffix)]] <- as(out$H_merged, "sparseMatrix")
  colMapTypes(adata)[[sprintf("H_%s", merged_suffix)]] <- "internal"

  colMaps(adata)[[sprintf("C_%s", merged_suffix)]] <- as(out$C_merged, "sparseMatrix")
  colMapTypes(adata)[[sprintf("C_%s", merged_suffix)]] <- "internal"

  obs <- .get_obs_data(adata)
  obs[[archetype_slot_out]] <- c(out$assigned_archetypes)
  adata <- .set_obs_data(adata, obs)
  return(adata)
}
