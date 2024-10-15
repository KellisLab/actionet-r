.groupedRowSums <- function(S, group_vec) {
  if (ACTIONetExperiment:::is.sparseMatrix(S)) {
    mat <- C_computeGroupedRowSumsSparse(S, sample_assignments = group_vec)
  } else {
    mat <- C_computeGroupedRowSumsDense(S, sample_assignments = group_vec)
  }
  return(mat)
}


.groupedRowMeans <- function(S, group_vec) {
  if (ACTIONetExperiment:::is.sparseMatrix(S)) {
    mat <- C_computeGroupedRowMeansSparse(S, sample_assignments = group_vec)
  } else {
    mat <- C_computeGroupedRowMeansDense(S, sample_assignments = group_vec)
  }
  return(mat)
}


.groupedRowVars <- function(S, group_vec) {
  if (ACTIONetExperiment:::is.sparseMatrix(S)) {
    mat <- C_computeGroupedRowVarsSparse(S, sample_assignments = group_vec)
  } else {
    mat <- C_computeGroupedRowVarsDense(S, sample_assignments = group_vec)
  }
  return(mat)
}


#' @export
aggregateMatrix <- function(S,
                            group_vec,
                            method = c("sum", "mean", "var")) {
  method <- match.arg(method, several.ok = FALSE)

  if (any(is.na(group_vec))) {
    err <- sprintf("'NA' values in 'group_vec'.\n")
    stop(err)
  }

  lf <- factor(.validate_attr(
    S,
    attr = group_vec,
    obj_name = "S",
    attr_name = "group_vec"
  ))
  labels <- as.numeric(lf)
  keys <- levels(lf)

  if (ACTIONetExperiment:::is.sparseMatrix(S) &&
    !is(S, "dMatrix")) {
    S <- as(S, "dMatrix")
  }

  if (method == "sum") {
    mat <- .groupedRowSums(S, labels)
  } else if (method == "mean") {
    mat <- .groupedRowMeans(S, labels)
  } else if (method == "var") {
    mat <- .groupedRowVars(S, labels)
  }

  colnames(mat) <- keys
  rownames(mat) <- rownames(S)
  return(mat)
}

warnifnot <- function(cond) {
  if (cond == TRUE) {
    warning("PASSED")
  } else {
    warning("FAILED")
  }
}

verify_aces <- function(ace1, ace2) {
  ###############################################################
  ###############################################################
  ################ Check ACTION reduction #######################
  ###############################################################
  ###############################################################
  A1 <- ace1$ACTION
  A2 <- ace2$ACTION

  deltaA <- sum(abs(A1 - A2)) / length(A1)
  warnifnot(deltaA < 1e-5)
  print(sprintf("Delta ACTION (reduction) = %.2e", deltaA))

  ###############################################################
  ###############################################################
  ############## Check ACTION decomposition #####################
  ###############################################################
  ###############################################################

  # check multi-level decomposition
  C1 <- colMaps(ace1)$C_stacked
  H1 <- colMaps(ace1)$H_stacked
  C2 <- colMaps(ace2)$C_stacked
  H2 <- colMaps(ace2)$H_stacked

  n1 <- ncol(H1)
  n2 <- ncol(H2)
  warnifnot(n1 == n2)
  print(sprintf("Same number of multi-level archetypes [%d]", n1))

  deltaH <- sum(abs(H1 - H2)) / length(H1)
  warnifnot(deltaH < 1e-5)
  print(sprintf("Delta H (multi-level) = %.2e", deltaH))

  deltaC <- sum(abs(C1 - C2)) / length(H1)
  warnifnot(deltaC < 1e-5)
  print(sprintf("Delta C (multi-level) = %.2e", deltaC))

  # check multi-level decomposition
  C1 <- colMaps(ace1)$C_merged
  H1 <- colMaps(ace1)$H_merged
  C2 <- colMaps(ace2)$C_merged
  H2 <- colMaps(ace2)$H_merged

  n1 <- ncol(H1)
  n2 <- ncol(H2)
  warnifnot(n1 == n2)
  print(sprintf("Same number of multi-resolution archetypes [%d]", n1))

  deltaH <- sum(abs(H1 - H2)) / length(H1)
  warnifnot(deltaH < 1e-5)
  print(sprintf("Delta H (multi-resolution) = %.2e", deltaH))

  deltaC <- sum(abs(C1 - C2)) / length(H2)
  warnifnot(deltaC < 1e-5)
  print(sprintf("Delta C (multi-resolution) = %.2e", deltaC))

  archs1 <- ace1$assigned_archetype
  archs2 <- ace2$assigned_archetype
  mismatch_perc <- 100 * sum(archs1 != archs2) / ncol(ace1)
  warnifnot(mismatch_perc < 0.5)
  print(sprintf("%.02f %% archetype assignment mismatch", mismatch_perc))


  ###############################################################
  ###############################################################
  ########## Check archetype feature specificity  ###############
  ###############################################################
  ###############################################################
  spec1 <- round(ace1$arch_feat_spec, 3)
  spec2 <- round(ace2$arch_feat_spec, 3)

  deltaSpec <- sum(abs(spec1 - spec2)) / length(spec1)
  warnifnot(deltaSpec < 1e-3)
  print(sprintf("Delta archetype feature specificity = %.2e", deltaSpec))


  ###############################################################
  ###############################################################
  ############## Check network construction #####################
  ###############################################################
  ###############################################################
  net1 <- round(ace1$ACTIONet, 3)
  net2 <- round(ace2$ACTIONet, 3)

  mismatch.edges <- 100 * sum(net1 != net2) / length((net1@i))
  warnifnot(mismatch.edges < 0.5)
  print(sprintf("%.02f %% ACTIONet edges mismatch", mismatch.edges))


  ###############################################################
  ###############################################################
  ############# Check network visualization  ####################
  ###############################################################
  ###############################################################
  ## 2D
  coor2D1 <- round(ace1$ACTIONet2D, 1)
  coor2D2 <- round(ace2$ACTIONet2D, 1)

  mismatch.2D <- 100 * sum(coor2D1 != coor2D2) / length(coor2D1)
  warnifnot(mismatch.2D < 0.5)
  print(sprintf("%.02f %% 2D mismatch", mismatch.2D))


  ## 3D
  coor3D1 <- round(ace1$ACTIONet3D, 1)
  coor3D2 <- round(ace2$ACTIONet3D, 1)

  mismatch.3D <- 100 * sum(coor3D1 != coor3D2) / length(coor3D1)
  warnifnot(mismatch.3D < 0.5)
  print(sprintf("%.02f %% 3D mismatch", mismatch.3D))

  ## Colors
  colors1 <- round(ace1$colors_actionet, 1)
  colors2 <- round(ace2$colors_actionet, 1)

  mismatch.colors <- 100 * sum(colors1 != colors2) / length(colors1)
  warnifnot(mismatch.colors < 0.5)
  print(sprintf("%.02f %% colors mismatch", mismatch.colors))
}

export_minimal_sce <- function(ace, export_logcounts = FALSE) {
  if (export_logcounts) {
    sce <- SingleCellExperiment(assays = list(counts = counts(ace), logcounts = logcounts(ace)))
  } else {
    sce <- SingleCellExperiment(assays = list(counts = counts(ace)))
  }
  colData(sce) <- colData(ace)
  rowData(sce) <- rowData(ace)

  return(sce)
}
