#' Run main ACTIONet pipeline
#'
#' @param ace `ACTIONetExperiment` containing an appropriate reduction in 'reduction_slot'.
#' @param k_min Minimum depth of decompositions. (default=2)
#' @param k_max Maximum depth of decompositions. (default=30)
#' @param assay_name Name of assay to be used. (default='logcounts')
#' @param reduction_slot Slot in colMaps(ace) containing reduced kernel. (default='ACTION')
#' @param net_slot_out  Name of slot in colMaps(ace) to store ACTIONet adjacency matrix. (default='ACTIONet')
#' @param min_obs Minimum number of observations required to construct an archetype. (default=2)
#' @param max_it Maximum number of iterations for ACTION algorithm. (default=50)
#' @param spec_th Defines the stringency of pruning nonspecific archetypes.
#' The larger the value, the more archetypes will be filtered out. (default=-3)
#' @param network_metric Distance metric with which to compute cell-to-cell similarity during network construction. Options are 'jsd' (Jensen-Shannon divergence), L2-norm ('l2'), and inner product ('ip'). (default='jsd')
#' @param network_algorithm Algorithm to use for network construction. Options are k-nearest neighbors ('knn') and k*-nearest neighbors ('k*nn'). (default='k*nn')
#' @param network_density Density factor of ACTIONet graph. (default=1)
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors. (default=TRUE)
#' @param layout_method Algorithm for computing plot layout. Options are UMAP ("umap") or t-UMAP ("tumap"). (default="umap")
#' @param layout_epochs Number of epochs for SGD algorithm. (default=250)
#' @param layout_parallel Run layout construction using multiple cores. May result in marginally different outputs across runs due to parallelization-induced randomization. (default=TRUE)
#' @param compute_specificity_parallel Run feature specificity enrichment using multiple cores. Setting this to `TRUE` on large datasets may cause an out of memory crash. (default=FALSE)
#' @param thread_no Number of parallel threads. (default=0)
#' @param seed Seed for random initialization. (default=0)
#'
#' @return An ACTIONetExperiment object
#'
#' @examples
#' ace <- runACTIONet(ace)
#' @export
runACTIONet <- function(ace,
                        k_min = 2,
                        k_max = 30,
                        assay_name = NULL,
                        reduction_slot = "action",
                        net_slot_out = "actionet",
                        min_obs = 2,
                        max_it = 50,
                        spec_th = -3,
                        network_metric = "jsd",
                        network_algorithm = "k*nn",
                        network_density = 1,
                        mutual_edges_only = TRUE,
                        layout_method = c("umap", "tumap", "largevis"),
                        layout_epochs = 100,
                        layout_spread = 1.0,
                        layout_min_dist = 1.0,
                        layout_3d = TRUE,
                        layout_parallel = TRUE,
                        compute_specificity_parallel = FALSE,
                        thread_no = 0,
                        seed = 0) {
  ace <- .validate_ace(ace, as_ace = TRUE, allow_se_like = TRUE, return_elem = TRUE)

  layout_method <- tolower(layout_method)
  layout_method <- match.arg(layout_method, several.ok = FALSE)

  if (is.null(assay_name)) {
    if ("default_assay" %in% names(metadata(ace))) {
      assay_name <- metadata(ace)[["default_assay"]]
      message(sprintf("Input assay_name is NULL. Setting assay_name to the metadata(ace)[['default_assay']] => %s", assay_name))
    } else {
      message(sprintf("Input assay_name is NULL. Setting assay_name to logcounts"))
      assay_name <- "logcounts"
    }
  }
  .validate_assay(ace, assay_name = assay_name, return_elem = FALSE)

  if (is.null(reduction_slot)) {
    if ("default_reduction" %in% names(metadata(ace))) {
      reduction_slot <- metadata(ace)[["default_reduction"]]
      message(sprintf("Input reduction_slot is NULL. Setting reduction_slot to the metadata(ace)[['default_reduction']] => %s", reduction_slot))
    } else {
      message(sprintf("Input reduction_slot is NULL. Setting reduction_slot to ACTION"))
      reduction_slot <- "ACTION"
    }
  }
  .validate_map(ace = ace, map_slot = reduction_slot, return_elem = FALSE)

  ace <- runACTION(
    ace,
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
  ace <- buildNetwork(
    ace,
    algorithm = network_algorithm,
    distance_metric = network_metric,
    density = network_density,
    thread_no = thread_no,
    mutual_edges_only = mutual_edges_only,
    map_slot = "H_stacked",
    net_slot_out = net_slot_out
  )

  ace <- networkCentrality(
    ace,
    label_attr = "assigned_archetype",
    algorithm = "local_coreness",
    thread_no = thread_no,
    net_slot = net_slot_out,
    attr_out = "node_centrality"
  )

  # Smooth archetype footprints
  ace <- networkDiffusion(
    ace,
    scores = colMaps(ace)[["H_merged"]],
    algorithm = "pagerank",
    thread_no = thread_no,
    net_slot = net_slot_out,
    map_slot_out = "archetype_footprint"
  )

  # Use archetypal reduction as initial coordinates for uwot
  # Need to reduce it to 3D coordinate space.
  red.out <- runSVD(
    X = scale(colMaps(ace)[["archetype_footprint"]]),
    k = 3,
    seed = seed,
    verbose = FALSE
  )
  initial_coordinates <- scale(red.out$u)

  slot_layout <- sprintf("%s_%s", layout_method, net_slot_out)
  layout_args <- list(
    method = layout_method,
    n_components = 2,
    min_dist = layout_min_dist,
    spread = layout_spread,
    n_epochs = layout_epochs,
    net_slot = net_slot_out,
    seed = seed,
    thread_no = ifelse(layout_parallel, thread_no, 1),
    map_slot_out = sprintf("%s_2d", slot_layout)
  )

  ace <- do.call(
    layoutNetwork,
    c(
      list(
        obj = ace,
        initial_coordinates = initial_coordinates
      ),
      layout_args
    )
  )

  if (layout_3d) {
    # Warm-start 3D layout using 2D embedding
    initial_coordinates <- cbind(
      colMaps(ace)[[layout_args$map_slot_out]],
      initial_coordinates[, 3]
    )

    layout_args$n_components <- 3
    layout_args$n_epochs <- layout_epochs / 2
    layout_args$map_slot_out <- sprintf("%s_3d", slot_layout)

    ace <- do.call(
      layoutNetwork,
      c(
        list(
          obj = ace,
          initial_coordinates = scale(initial_coordinates)
        ),
        layout_args
      )
    )

    ace <- computeNodeColors(
      ace,
      embedding_slot = layout_args$map_slot_out,
      color_slot_out = sprintf("%s_colors", slot_layout),
      thread_no = thread_no
    )
  }

  # Compute gene specificity for each archetype
  ace <- archetypeFeatureSpecificity(
    ace = ace,
    assay_name = assay_name,
    map_slot = "archetype_footprint",
    thread_no = ifelse(compute_specificity_parallel, thread_no, 1),
    return_raw = FALSE
  )

  return(ace)
}

#' Reconstructs the ACTIONet graph with the new parameters (uses prior decomposition)
#'
#' @param ace ACTIONetExperiment object.
#' @param network_density Density factor of ACTIONet graph. (default=1)
#' @param network_metric Distance metric with which to compute cell-to-cell similarity during network construction. Options are 'jsd' (Jensen-Shannon divergence), L2-norm ('l2'), and inner product ('ip'). (default='jsd')
#' @param algorithm Algorithm to use for network construction. Options are k-nearest neighbors ('knn') and k*-nearest neighbors ('k*nn'). (default='k*nn')
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors. (default=TRUE)
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout (default=50).
#' @param layout_epochs Number of epochs for SGD algorithm. (default=250)
#' @param layout_algorithm Algorithm for computing plot layout. Options are UMAP ("umap") or t-UMAP ("tumap"). (default="umap")
#' @param layout_parallel Run layout construction using multiple cores. May result in marginally different outputs across runs due to parallelization-induced randomization. (default=TRUE)
#' @param thread_no Number of parallel threads. (default=0)
#' @param reduction_slot Slot in colMaps(ace) containing reduced kernel. (default='ACTION')
#' @param net_slot_out Name of slot in colMaps(ace) to store ACTIONet adjacency matrix. (default='ACTIONet')
#' @param seed Seed for random initialization. (default=0)
#'
#' @return ace Updated ace object
#'
#' @examples
#' plot.ACTIONet(ace)
#' ace.updated <- rebuildACTIONet(ace, network_density = 0.1)
#' plot.ACTIONet(ace.updated)
#' @export
rebuildACTIONet <- function(ace,
                            network_density = 1,
                            network_metric = "jsd",
                            algorithm = "k*nn",
                            mutual_edges_only = TRUE,
                            layout_spread = 1.0,
                            layout_min_dist = 1.0,
                            layout_gamma = 1.0,
                            layout_epochs = 100,
                            layout_algorithm = c("umap", "tumap"),
                            layout_parallel = TRUE,
                            thread_no = 0,
                            reduction_slot = "ACTION",
                            net_slot_out = "ACTIONet",
                            H_stacked_slot = "H_stacked",
                            reduction_normalization = 1,
                            seed = 0) {
  set.seed(seed)
  layout_algorithm <- tolower(layout_algorithm)
  layout_algorithm <- match.arg(layout_algorithm, several.ok = FALSE)

  .validate_ace(ace, return_elem = FALSE)
  .validate_map(ace, map_slot = reduction_slot, return_elem = FALSE)

  if (!(sprintf("%s_normalized", reduction_slot) %in% names(colMaps(ace)))) {
    ace <- normalize.reduction(ace, reduction_slot = reduction_slot, reduction_normalization = reduction_normalization)
  }
  reduction_slot <- sprintf("%s_normalized", reduction_slot)

  # re-Build ACTIONet
  H_stacked <- .validate_map(
    ace = ace,
    map_slot = H_stacked_slot,
    matrix_type = "dense",
    force_type = TRUE
  )
  H_stacked <- Matrix::t(H_stacked)

  G <- buildNetwork(
    H = H_stacked,
    algorithm = algorithm,
    distance_metric = network_metric,
    density = network_density,
    thread_no = thread_no,
    mutual_edges_only = mutual_edges_only
  )
  colNets(ace)[[net_slot_out]] <- G

  # Layout ACTIONet
  ace <- rerunLayout(
    ace = ace,
    algorithm = layout_algorithm,
    n_epochs = layout_epochs,
    spread = layout_spread,
    min_dist = layout_min_dist,
    gamma = layout_gamma,
    net_slot = net_slot_out,
    thread_no = ifelse(layout_parallel, thread_no, 1),
    seed = seed
  )

  return(ace)
}


#' Rerun layout on the ACTIONet graph with new parameters
#'
#' @param ace ACTIONetExperiment object.
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout (default=50).
#' @param layout_epochs Number of epochs for SGD algorithm (default=250).
#' @param layout_algorithm Algorithm for computing plot layout. Options are UMAP ("umap") or t-UMAP ("tumap"). (default="umap")
#' @param network_density Density factor of ACTIONet graph (default=1).
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors (default=TRUE).
#' @param thread_no Number of parallel threads (default=0).
#' @param reduction_slot Slot in colMaps(ace) containing reduced kernel (default='ACTION').
#' @param net_slot Slot in colMaps(ace) containing ACTIONet adjacency matrix (default='ACTIONet').
#' @param seed Seed for random initialization (default=0).
#'
#' @return ace Updated ace object
#'
#' @examples
#' plot.ACTIONet(ace)
#' ace.updated <- rerunLayout(ace, compactness = 20)
#' plot.ACTIONet(ace.updated)
#' @export
rerunLayout <- function(ace,
                        initial_coordinates = NULL,
                        algorithm = c("umap", "tumap"),
                        min_dist = 1.0,
                        spread = 1.0,
                        gamma = 1.0,
                        n_epochs = 100,
                        learning_rate = 1.0,
                        net_slot = "ACTIONet",
                        seed = 0,
                        thread_no = 0) {
  algorithm <- tolower(algorithm)
  algorithm <- match.arg(algorithm, several.ok = FALSE)

  if (is.null(initial_coordinates)) {
    red.out <- runSVD_full(scale(ace$archetype_footprint), k = 3, algorithm = 0, max_it = 1000, verbose = 0, seed = 0)
    initial_coordinates <- scale(red.out$u[, 1:3])
  }

  ace <- .run.layoutNetwork(
    ace = ace,
    initial_coordinates = initial_coordinates,
    algorithm = algorithm,
    min_dist = min_dist,
    spread = spread,
    gamma = gamma,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    net_slot = net_slot,
    thread_no = thread_no,
    seed = seed
  )

  return(ace)
}


#' @export
rerun.archetype.unification <- function(ace,
                                        backbone_density = 0.5,
                                        resolution = 1.0,
                                        min_cluster_size = 3,
                                        footprint_alpha = 0.85,
                                        compute_specificity_parallel = FALSE,
                                        assay_name = "logcounts",
                                        reduction_slot = "ACTION",
                                        C_stacked_slot = "C_stacked",
                                        H_stacked_slot = "H_stacked",
                                        net_slot = "ACTIONet",
                                        merged_suffix = "merged",
                                        reduction_normalization = 1,
                                        thread_no = 0) {
  .validate_ace(ace, return_elem = FALSE)

  G <- .validate_net(
    ace = ace,
    net_slot = net_slot,
    matrix_type = "sparse",
    force_type = FALSE
  )

  if (!(sprintf("%s_normalized", reduction_slot) %in% names(colMaps(ace)))) {
    ace <- normalize.reduction(ace, reduction_slot = reduction_slot, reduction_normalization = reduction_normalization)
  }
  reduction_slot <- sprintf("%s_normalized", reduction_slot)

  ace <- .run.mergeArchetypes(
    ace = ace,
    reduction_slot = reduction_slot,
    C_stacked_slot = C_stacked_slot,
    H_stacked_slot = H_stacked_slot,
    normalization = 0,
    merged_suffix = merged_suffix,
    thread_no = thread_no
  )

  colData(ace)[["node_centrality"]] <- networkCentrality(
    obj = G,
    label_attr = colData(ace)[["assigned_archetype"]],
    algorithm = "local_coreness"
  )

  archetype_footprint <- networkDiffusion(
    obj = G,
    scores = colMaps(ace)[[sprintf("H_%s", merged_suffix)]],
    algorithm = "pagerank",
    alpha = footprint_alpha,
    thread_no = thread_no,
    max_it = 5,
    tol = 1e-8,
    net_slot = NULL
  )
  colMaps(ace)$archetype_footprint <- archetype_footprint

  # Compute gene specificity for each archetype
  ace <- archetypeFeatureSpecificity(
    ace = ace,
    assay_name = assay_name,
    map_slot = "archetype_footprint",
    thread_no = ifelse(compute_specificity_parallel, thread_no, 1),
    return_raw = FALSE
  )

  ace <- constructBackbone(ace, backbone_density = backbone_density)

  return(ace)
}

#' @export
runACTION <- function(
    ace,
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
    archetype_slot_out = "assigned_archetype") {
  ace <- .validate_ace(
    ace,
    as_ace = TRUE,
    allow_se_like = TRUE,
    return_elem = TRUE
  )

  S_r <- .validate_map(
    ace = ace,
    map_slot = reduction_slot,
    matrix_type = "dense",
    force_type = TRUE,
    return_elem = TRUE
  )

  if (prenormalize) {
    S_r <- normalize.matrix(
      S_r,
      dim = 1, # obs are stored rows in colMaps()
      scale_param = NULL,
      trans_func = NULL
    )
  }

  out <- C_runACTION(
    S_r = Matrix::t(S_r),
    k_min = k_min,
    k_max = k_max,
    max_it = max_it,
    tol = tol,
    spec_th = spec_th,
    min_obs = min_obs,
    thread_no = thread_no
  )

  colMaps(ace)[["H_stacked"]] <- Matrix::t(as(out$H_stacked, "sparseMatrix"))
  colMapTypes(ace)[["H_stacked"]] <- "internal"

  colMaps(ace)[["C_stacked"]] <- as(out$C_stacked, "sparseMatrix")
  colMapTypes(ace)[["C_stacked"]] <- "internal"

  colMaps(ace)[[sprintf("H_%s", merged_suffix)]] <- as(Matrix::t(out$H_merged), "sparseMatrix")
  colMapTypes(ace)[[sprintf("H_%s", merged_suffix)]] <- "internal"

  colMaps(ace)[[sprintf("C_%s", merged_suffix)]] <- as(out$C_merged, "sparseMatrix")
  colMapTypes(ace)[[sprintf("C_%s", merged_suffix)]] <- "internal"

  colData(ace)[[archetype_slot_out]] <- c(out$assigned_archetype)
  return(ace)
}
