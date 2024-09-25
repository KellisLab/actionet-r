#' @export
findMarkers.ACTIONet <- function(ace,
                                 cluster_attr,
                                 top_genes = 10,
                                 most_specific = FALSE,
                                 features_use = NULL,
                                 feat_subset = NULL,
                                 assay_name = "logcounts",
                                 thread_no = 0,
                                 to_return = c("data.frame", "df", "list")) {
  to_return <- match.arg(to_return)

  sa <- ACTIONetExperiment::get.data.or.split(ace, attr = cluster_attr, to_return = "levels")
  features_use <- .get_feature_vec(ace, features_use = features_use)

  specificity <- computeFeatureSpecificity(
    obj = ace,
    cluster_attr = sa$index,
    assay_name = assay_name,
    thread_no = thread_no,
    return_raw = TRUE
  )

  feat_spec <- specificity[["upper_significance"]] - specificity[["lower_significance"]]
  feat_spec[feat_spec < 0] <- 0
  rownames(feat_spec) <- features_use
  colnames(feat_spec) <- sa$keys

  if (!is.null(feat_subset)) {
    feat_spec <- feat_spec[rownames(feat_spec) %in% feat_subset, , drop = FALSE]
  }

  if (most_specific == TRUE) {
    W <- select.top.k.features(
      feat_spec,
      top_features = top_genes,
      normalize = FALSE,
      reorder_columns = FALSE
    )
    feat_spec_top <- apply(W, 2, function(v) rownames(W)[order(v, decreasing = TRUE)][1:top_genes])
  } else {
    feat_spec_top <- sapply(colnames(feat_spec), function(type) {
      c <- feat_spec[, type]
      names(head(sort(c, decreasing = TRUE), top_genes))
    })
  }

  df <- data.frame(feat_spec_top)

  if (to_return == "list") {
    return(as.list(df))
  } else {
    return(df)
  }
}

computeGeneSpecifity.ACTIONet <- function(ace, f, out_name = "cond", pos_only = T, blacklist_pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT") {
  require(ACTIONet)
  print("Running computeGeneSpecifity.ACTIONet()")

  if (class(f) != "factor") {
    warning("f must be a factor")
    return(ace)
  }
  f <- droplevels(f)

  S <- logcounts(ace)

  if (is.matrix(S)) {
    out <- compute_cluster_feature_specificity_full(S, as.numeric(f))
  } else {
    out <- compute_cluster_feature_specificity(S, as.numeric(f))
  }
  metadata(ace)[[sprintf("%s_feature_specificity_ACTIONet", out_name)]] <- out

  scores <- as.matrix(out$upper_significance - out$lower_significance)
  colnames(scores) <- levels(f)
  rownames(scores) <- rownames(ace)

  if (pos_only == T) {
    scores[scores < 0] <- 0
  } # Only "positive markers" [negative would be markers in other levels]

  blacklisted.rows <- grep(blacklist_pattern, rownames(ace), ignore.case = T)
  scores[blacklisted.rows, ] <- 0

  rowMaps(ace)[[sprintf("%s_feature_specificity", out_name)]] <- scores
  rowMapTypes(ace)[[sprintf("%s_feature_specificity", out_name)]] <- "reduction"

  return(ace)
}
