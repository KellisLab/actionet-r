#' @export
clusterNetwork <- function(
    obj,
    objective_function = c("modularity", "CPM"),
    resolution_parameter = 1.0,
    initial_membership = NULL,
    n_iterations = 3,
    min_size = 3,
    net_slot = "actionet",
    attr_out = NULL,
    return_raw = FALSE) {
  is_installed <- require(igraph)

  if (!is_installed) {
    stop("Package 'igraph' is not installed")
  }

  # algorithm <- match.arg(algorithm)
  algorithm <- "leiden"
  objective_function <- match.arg(objective_function)
  is_ace <- .validate_ace(obj, error_on_fail = FALSE, return_elem = FALSE)

  G <- .ace_or_net(
    obj = obj,
    net_slot = net_slot,
    matrix_type = "sparse",
    force_type = TRUE,
    obj_name = "obj"
  )

  if (!is.null(initial_membership)) {
    initial_membership <- .validate_vector_attr(
      obj,
      attr = initial_membership,
      return_type = "data",
      attr_name = "initial_membership",
      return_elem = TRUE
    )
    initial_membership <- as.numeric(as.factor(initial_membership))
  }

  ig <- igraph::graph_from_adjacency_matrix(
    adjmatrix = G,
    mode = "undirected",
    weighted = TRUE,
    diag = TRUE,
    add.colnames = NULL,
    add.rownames = NA
  )

  # Put this into separate hidden functions once more modes are supported.
  if (algorithm == "leiden") {
    comm <- igraph::cluster_leiden(
      graph = ig,
      objective_function = objective_function,
      weights = NULL,
      resolution_parameter = resolution_parameter,
      beta = 0.01,
      initial_membership = initial_membership,
      n_iterations = n_iterations,
      vertex_weights = NULL
    )
  } else {
    stop("Invalid algorithm")
  }

  clusters <- igraph::membership(comm)
  c0 <- names(which(table(clusters) < min_size))
  if (length(c0) > 0) {
    clusters[clusters %in% as.numeric(c0)] <- 0
  }


  if (is_ace && !return_raw) {
    if (is.null(attr_out)) {
      attr_out <- sprintf("%s_%s", algorithm, net_slot)
    }
    colData(obj)[[attr_out]] <- clusters
    return(obj)
  }

  return(clusters)
}
