#' Annotate clusters using prior cell annotations
#' (It uses Fisher's exact test for computing overlaps -- approximate HGT is used)
#'
#' @param ace ACTIONet output object
#' @param clusters Cluster
#' @param labels Annotation of interest (clusters, celltypes, etc.) to test enrichment
#'
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#' }
#'
#' @examples
#' arch.annot <- annotate.clusters.using.labels(ace, ace$clusters, sce$celltypes)
#' @export
annotate.clusters.using.labels <- function(ace,
                                           clusters,
                                           labels) {
  clusters <- .preprocess_annotation_labels(clusters, ace)
  Labels <- .preprocess_annotation_labels(labels, ace)

  pop.size <- length(Labels)
  pos.size <- table(Labels)

  logPvals <- sapply(sort(unique(clusters)), function(i) {
    idx <- which(clusters == i)
    sample.size <- length(idx)
    success.size <- sapply(sort(unique(Labels)), function(i) {
      sum(Labels[idx] == i)
    })

    logPval <- HGT_tail(
      population.size = pop.size,
      success.count = pos.size,
      sample.size = sample.size,
      observed.success = success.size
    )

    return(logPval)
  })


  cl.Annot <- names(clusters)[match(sort(unique(clusters)), clusters)]
  Annot <- names(Labels)[match(sort(unique(Labels)), Labels)]

  colnames(logPvals) <- cl.Annot
  rownames(logPvals) <- Annot

  clusterLabels <- Annot[apply(logPvals, 2, which.max)]

  cellLabels <- match(clusterLabels[clusters], Annot)
  names(cellLabels) <- clusterLabels[clusters]

  res <- list(
    Labels = clusterLabels,
    cellLabels = cellLabels,
    Enrichment = logPvals
  )

  return(res)
}


#' Annotate clusters using known marker genes
#' (It uses permutation test on cluster specificity scores)
#'
#' @param ace ACTIONet output object
#' @param marker.genes A list of lists (each a set of markers for a given cell type)
#' @param specificity.slot.name An entry in the rowMaps(ace), precomputed using clusterFeatureSpecificity() function
#' @param rand.sample.no Number of random permutations (default=1000)
#'
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#' }
#'
#' @examples
#' data("curatedMarkers_human") # pre-packaged in ACTIONet
#' marker.genes <- curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' ace <- clusterFeatureSpecificity(ace, ace$clusters, "cluster_specificity_scores")
#' arch.annot <- annotate.clusters.using.markers(ace, marker.genes = marker.genes, specificity.slot.name = "cluster_specificity_scores")
#' @export
annotate.clusters.using.markers <- function(ace,
                                            marker.genes,
                                            specificity.slot.name,
                                            rand.sample.no = 1000) {
  if (!(specificity.slot.name %in% names(rowMaps(ace)))) {
    message(sprintf("%s does not exist in rowMaps(ace)", specificity.slot.name))
  }

  if (is.matrix(marker.genes) | ACTIONetExperiment:::is.sparseMatrix(marker.genes)) {
    marker.genes <- apply(marker.genes, 2, function(x) {
      rownames(marker.genes)[x >
        0]
    })
  }

  specificity.panel <- Matrix::t(as.matrix(log1p(rowMaps(ace)[[specificity.slot.name]])))

  GS.names <- names(marker.genes)
  if (is.null(GS.names)) {
    GS.names <- sapply(1:length(GS.names), function(i) sprintf("Celltype %s", i))
  }

  markers.table <- do.call(rbind, lapply(names(marker.genes), function(celltype) {
    genes <- marker.genes[[celltype]]

    if (length(genes) == 0) {
      err <- sprintf("No markers left.\n")
      stop(err, call. = FALSE)
    }

    signed.count <- sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
    is.signed <- signed.count > 0

    if (!is.signed) {
      df <- data.frame(
        Gene = genes,
        Direction = +1,
        Celltype = celltype,
        stringsAsFactors = FALSE
      )
    } else {
      pos.genes <- (as.character(sapply(
        genes[grepl("+", genes, fixed = TRUE)],
        function(gene) stringr::str_replace(gene, stringr::fixed("+"), "")
      )))
      neg.genes <- (as.character(sapply(
        genes[grepl("-", genes, fixed = TRUE)],
        function(gene) stringr::str_replace(gene, stringr::fixed("-"), "")
      )))

      df <- data.frame(
        Gene = c(pos.genes, neg.genes),
        Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))),
        Celltype = celltype,
        stringsAsFactors = FALSE
      )
    }
  }))

  markers.table <- markers.table[markers.table$Gene %in% colnames(specificity.panel), ]

  if (dim(markers.table)[1] == 0) {
    err <- sprintf("No markers left.\n")
    stop(err, call. = FALSE)
  }

  specificity.panel <- specificity.panel[, markers.table$Gene]

  IDX <- split(1:dim(markers.table)[1], markers.table$Celltype)

  print("Computing significance scores")
  set.seed(0)
  Z <- sapply(IDX, function(idx) {
    markers <- (as.character(markers.table$Gene[idx]))
    directions <- markers.table$Direction[idx]
    mask <- markers %in% colnames(specificity.panel)

    A <- as.matrix(specificity.panel[, markers[mask]])
    sgn <- as.numeric(directions[mask])
    stat <- A %*% sgn

    rand.stats <- sapply(1:rand.sample.no, function(i) {
      rand.samples <- sample.int(dim(specificity.panel)[2], sum(mask))
      rand.A <- as.matrix(specificity.panel[, rand.samples])
      rand.stat <- rand.A %*% sgn
    })

    cell.zscores <- as.numeric((stat - apply(rand.stats, 1, mean)) / apply(
      rand.stats,
      1, sd
    ))

    return(cell.zscores)
  })

  Z[is.na(Z)] <- 0
  Labels <- colnames(Z)[apply(Z, 1, which.max)]

  Labels.conf <- apply(Z, 1, max)

  names(Labels) <- rownames(specificity.panel)
  names(Labels.conf) <- rownames(specificity.panel)
  rownames(Z) <- rownames(specificity.panel)

  out <- list(
    Label = Labels,
    Confidence = Labels.conf,
    Enrichment = Z
  )

  return(out)
}


#' Annotate arbitary feature score matrix using known marker genes
#' (It uses permutation test on cluster specificity scores)
#'
#' @param profile An arbitrary matrix with rows corresponding to features and columns to any given annotation/grouping of cells
#' @param markers A list of lists (each a set of markers for a given cell type)
#'
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#' }
#'
#' @export
annotate.profile.using.markers <- function(profile, markers) {
  marker_mat <- .preprocess_annotation_markers(markers, rownames(profile))

  enrichment.out <- assess_enrichment(profile, marker_mat)

  X <- enrichment.out$logPvals
  rownames(X) <- colnames(marker_mat)
  colnames(X) <- colnames(profile)

  return(X)
}

#' @export
clusterNetwork <- function(
    obj,
    algorithm = c("leiden"),
    resolution_parameter = 1.0,
    initial_membership = NULL,
    n_iterations = 3,
    net_slot = "actionet",
    attr_out = NULL,
    return_raw = TRUE) {
  is_installed <- require(igraph)

  if (!is_installed) {
    stop("Package 'igraph' is not installed")
  }

  algorithm <- tolower(algorithm)
  is_ace <- .validate_ace(obj, error_on_fail = FALSE, return_elem = FALSE)

  G <- .ace_or_net(
    obj = obj,
    net_slot = net_slot,
    matrix_type = "sparse",
    force_type = TRUE,
    obj_name = "obj"
  )

  if (!is.null(initial_membership)) {
    initial_membership <- .validate_attr(
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
      objective_function = "modularity",
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
  if (is_ace && !return_raw) {
    if (is.null(attr_out)) {
      attr_out <- sprintf("%s_%s", algorithm, net_slot)
    }
    colData(obj)[[attr_out]] <- clusters
    return(obj)
  }

  return(clusters)
}


#' @export
clusterCells <- function(ace, algorithm = "leiden",
                         cluster_name = "leiden",
                         resolution_parameter = 1.0,
                         initial_clusters = NULL,
                         seed = 0,
                         net_slot = "ACTIONet") {
  if (!is.null(initial_clusters)) {
    if (is.character(initial_clusters)) {
      initial_clusters <- as.factor(initial_clusters)
    }
  }
  if (algorithm == "fix") {
    cl <- initial_clusters
  } else {
    if (is.null(initial_clusters)) {
      initial_clusters <- ace$assigned_archetype
    }
    cl <- clusterNetwork(ace,
      algorithm = algorithm,
      resolution_parameter = resolution_parameter,
      initial_clusters = initial_clusters,
      seed = seed,
      net_slot = net_slot
    )
  }

  colData(ace)[[cluster_name]] <- cl
  ace <- computeGeneSpecifity.ace(ace, cl, out_name = cluster_name)

  return(ace)
}
