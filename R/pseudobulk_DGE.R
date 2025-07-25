#' @export
get.pseudobulk.SE <- function(
    ace,
    sample_attr,
    ensemble = FALSE,
    bins = 20,
    assay_name = "counts",
    col_data = NULL,
    pseudocount = 0,
    with_S = FALSE,
    with_E = FALSE,
    with_V = FALSE,
    min_cells_per_batch = 3,
    BPPARAM = BiocParallel::SerialParam()) {
  group_vec <- ACTIONetExperiment::get.data.or.split(ace, attr = sample_attr, to_return = "data")
  sample_counts <- table(group_vec)
  good_samples <- sample_counts >= min_cells_per_batch

  if (!any(good_samples)) {
    msg <- sprintf("No samples remaining.")
    warning(msg)
    return(NULL)
  } else if (!all(good_samples)) {
    old_samples <- names(sample_counts)
    sample_counts <- sample_counts[good_samples]
    ace <- ace[, ace[[sample_attr]] %in% names(sample_counts)]
    group_vec <- ACTIONetExperiment::get.data.or.split(ace, attr = sample_attr, to_return = "data")
    bad_sample_names <- setdiff(old_samples, names(sample_counts))
    msg <- sprintf("Samples Dropped: %s\n", paste0(bad_sample_names, collapse = ", "))
    message(msg)
  }

  counts_mat <- SummarizedExperiment::assays(ace)[[assay_name]]
  sample_names <- levels(factor(group_vec))

  se_assays <- list()

  se_assays$counts <- aggregateMatrix(counts_mat, group_vec = group_vec, method = "sum") + pseudocount

  if (with_E == TRUE) {
    se_assays$mean <- aggregateMatrix(counts_mat, group_vec = group_vec, method = "mean")
  }

  if (with_V == TRUE) {
    se_assays$var <- aggregateMatrix(counts_mat, group_vec = group_vec, method = "var")
  }

  if (ensemble == TRUE) {
    if (!any(with_S, with_E, with_V)) {
      err <- sprintf("No ensemble assays to make.\n")
      stop(err)
    }

    IDX <- ACTIONetExperiment::get.data.or.split(ace, attr = sample_attr, to_return = "split")
    counts_list <- bplapply(IDX, function(idx) counts_mat[, idx, drop = FALSE], BPPARAM = BPPARAM)

    mr_assays <- .make_ensemble_assays(
      counts_list = counts_list,
      bins = bins,
      pseudocount = pseudocount,
      with_S = with_S,
      with_E = with_E,
      with_V = with_V,
      BPPARAM = BPPARAM
    )
    se_assays <- c(se_assays, mr_assays$assays)
  }

  cd <- data.frame(
    n_cells = c(sample_counts[sample_names]),
    sample = sample_names
  )

  if (!is.null(col_data)) {
    md <- col_data[col_data[[sample_attr]] %in% sample_names, , drop = FALSE]
    md <- md[match(sample_names, md[[sample_attr]]), , drop = FALSE]
    cd <- data.frame(md, cd)
  }
  rownames(cd) <- sample_names
  cd <- droplevels(cd)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = se_assays,
    colData = cd,
    rowData = SummarizedExperiment::rowData(ace)
  )

  if (ensemble == TRUE) {
    S4Vectors::metadata(se)[["bins"]] <- bins
  }

  invisible(gc())
  return(se)
}


.make_ensemble_assays <- function(
    counts_list,
    bins,
    pseudocount,
    with_S = FALSE,
    with_E = FALSE,
    with_V = FALSE,
    BPPARAM = BiocParallel::SerialParam()) {
  mr_lists <- bplapply(names(counts_list), function(n) {
    S <- Matrix::t(counts_list[[n]])
    bin_IDX <- round(seq(1, (1 - bins^-1) * nrow(S), length.out = bins))

    out <- list()

    if (with_S == TRUE) {
      S.sorted <- apply(S, 2, function(s) cumsum(sort(s)))
      cs <- S.sorted[nrow(S.sorted), ]
      S.prefix_sum <- sapply(bin_IDX, function(idx) cs - S.sorted[idx, , drop = FALSE])
      out <- list(Sp = S.prefix_sum)
    }

    if (with_E == TRUE) {
      E.sorted <- apply(S, 2, function(s) rev(cumsum(sort(s, decreasing = TRUE)) / seq.int(length(s))))
      E.prefix_sum <- sapply(bin_IDX, function(idx) E.sorted[idx, , drop = FALSE])
      out$Ep <- E.prefix_sum
    }
    if (with_V == TRUE) {
      V.sorted <- apply(S, 2, function(s) rev(roll_var(sort(s, decreasing = T))))
      V.sorted[is.na(V.sorted)] <- 0
      V.prefix_sum <- sapply(bin_IDX, function(idx) V.sorted[idx, , drop = FALSE])
      out$Vp <- V.prefix_sum
    }
    return(out)
  }, BPPARAM = BPPARAM)

  mr_assays <- list()

  if (with_S == TRUE) {
    S_list <- lapply(1:bins, function(i) {
      sapply(mr_lists, function(L) L$Sp[, i, drop = FALSE]) + pseudocount
    })
    names(S_list) <- paste0("S", 1:bins)
    mr_assays <- c(mr_assays, S_list)
  }

  if (with_E == TRUE) {
    E_list <- lapply(1:bins, function(i) {
      sapply(mr_lists, function(L) L$Ep[, i, drop = FALSE])
    })
    names(E_list) <- paste0("E", 1:bins)
    mr_assays <- c(mr_assays, E_list)
  }

  if (with_V == TRUE) {
    V_list <- lapply(1:bins, function(i) {
      sapply(mr_lists, function(L) L$Vp[, i, drop = FALSE])
    })
    names(V_list) <- paste0("V", 1:bins)
    mr_assays <- c(mr_assays, V_list)
  }

  mr_out <- list(assays = mr_assays)

  invisible(gc())
  return(mr_out)
}

#' @export
run.ensemble.pseudobulk.DESeq <- function(
    se,
    design,
    bins = NULL,
    bins_use = NULL,
    slot_prefix = "S",
    p_adj_method = "fdr",
    BPPARAM = BiocParallel::SerialParam()) {
  ACTIONetExperiment:::.check_and_load_package("DESeq2")

  if (is.null(bins)) {
    bins <- S4Vectors::metadata(se)$bins
  }

  if (is.null(bins_use)) {
    bins_use <- 1:bins
  }

  dds_out <- bplapply(paste0(slot_prefix, 1:bins), function(S) {
    cts <- SummarizedExperiment::assays(se)[[S]]
    invisible({
      dds <- DESeq2::DESeqDataSetFromMatrix(
        cts,
        design = design,
        colData = SummarizedExperiment::colData(se),
        rowData = SummarizedExperiment::rowData(se)
      )

      dds <- DESeq2::DESeq(dds, betaPrior = FALSE)
    })
    return(dds)
  }, BPPARAM = BPPARAM)

  dds_res <- lapply(dds_out, function(dds) {
    dds_res <- results(dds)
    out <- data.frame(rowData(dds), dds_res)
    out <- out[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "dispOutlier")]
    return(out)
  })

  outlier_count <- Matrix::rowSums(sapply(dds_res, function(dds) dds$dispOutlier))

  mat_Bstatm <- sapply(dds_res, function(dds) dds$baseMean)
  bm_mean <- Matrix::rowMeans(mat_Bstatm, na.rm = TRUE)

  mat_lfc <- sapply(dds_res, function(dds) dds$log2FoldChange)
  lfc_mean <- apply(mat_lfc, 1, function(r) mean(r, na.rm = TRUE, trim = 0.25))

  mat_se <- sapply(dds_res, function(dds) dds$lfcSE)
  se_mean <- apply(mat_se, 1, function(r) mean(r, na.rm = TRUE, trim = 0.25))

  mat.stat <- sapply(dds_res, function(dds) dds$stat)
  stat_mean <- Matrix::rowMeans(mat.stat, na.rm = TRUE)

  mat_pval <- sapply(dds_res, function(dds) dds$pvalue)
  mat_pval[is.na(mat_pval)] <- 1
  mat_pval <- Matrix::t(-log10(mat_pval))
  corr_pvals <- combine.logPvals(mat_pval, base = 10)
  corr_pvals <- 10^-corr_pvals

  res <- data.frame(
    SummarizedExperiment::rowData(se),
    baseMean = bm_mean,
    log2FCMean = lfc_mean,
    lfcSEMean = se_mean,
    statMean = stat_mean,
    pvalue = corr_pvals,
    padj = stats::p.adjust(corr_pvals, method = p_adj_method),
    dispOutlierFrac = outlier_count / bins
  )

  invisible(gc())
  return(res)
}

#' @export
run.ensemble.pseudobulk.Limma <- function(
    se,
    design,
    bins = NULL,
    bins_use = NULL,
    variable_name = NULL,
    min_covered_samples = 2,
    p_adj_method = "fdr",
    BPPARAM = BiocParallel::SerialParam()) {
  ACTIONetExperiment:::.check_and_load_package(c("SummarizedExperiment", "limma"))

  if (class(se) != "SummarizedExperiment") {
    stop("'se' must be an object of type 'SummarizedExperiment'.")
  }

  if (is.null(bins)) {
    bins <- S4Vectors::metadata(se)[["bins"]]
  }

  if (is.null(bins_use)) {
    bins_use <- 1:bins
  }

  design_mat <- .make_design_mat(design, data = SummarizedExperiment::colData(se))

  design_list <- .preprocess_design_matrix_and_var_names(design_mat, variable_name)
  design_mat <- design_list$design_mat
  variable_name <- design_list$variable_name

  out <- bplapply(bins_use, function(i) {
    E <- SummarizedExperiment::assays(se)[[paste0("E", i)]]
    V <- SummarizedExperiment::assays(se)[[paste0("V", i)]]
    W <- 1 / V

    slot_E <- paste0("E", i)
    slot_V <- paste0("V", i)

    tbl <- variance.adjusted.limma(
      se = se,
      slot_E = slot_E,
      slot_V = slot_V,
      design = design_mat,
      variable_name = variable_name,
      min_covered_samples = min_covered_samples
    )
    return(tbl)
  }, BPPARAM = BPPARAM)

  common <- lapply(out, function(tbl) rownames(tbl))
  common <- Reduce(intersect, common)
  out <- lapply(out, function(tbl) tbl[rownames(tbl) %in% common, ])

  mat_Bstatm <- sapply(out, function(tbl) tbl$AveExpr)
  bm_mean <- Matrix::rowMeans(mat_Bstatm, na.rm = TRUE)

  mat_lfc <- sapply(out, function(tbl) tbl$logFC)
  lfc_mean <- apply(mat_lfc, 1, function(r) mean(r, na.rm = TRUE, trim = 0.25))

  mat_tstat <- sapply(out, function(tbl) tbl$t)
  t_mean <- Matrix::rowMeans(mat_tstat, na.rm = TRUE)

  mat_Bstat <- sapply(out, function(tbl) tbl$B)
  B_mean <- Matrix::rowMeans(mat_Bstat, na.rm = TRUE)

  mat_pval <- sapply(out, function(tbl) tbl$P.Value)
  mat_pval[is.na(mat_pval)] <- 1
  mat_pval <- Matrix::t(-log10(mat_pval))
  corr_pvals <- combine.logPvals(mat_pval, base = 10)
  corr_pvals <- 10^-corr_pvals

  rdat <- SummarizedExperiment::rowData(se)[rownames(se) %in% common, ]
  res <- data.frame(
    rdat,
    aveExpr = bm_mean,
    log2FCMean = lfc_mean,
    tStatMean = t_mean,
    BStatMean = B_mean,
    pvalue = corr_pvals,
    padj = stats::p.adjust(corr_pvals, method = p_adj_method)
  )

  invisible(gc())
  return(res)
}

#' @export
variance.adjusted.limma <- function(
    se,
    design,
    slot_E = "counts",
    slot_V = NULL,
    W_mat = NULL,
    variable_name = NULL,
    min_covered_samples = 3) {
  ACTIONetExperiment:::.check_and_load_package(c("SummarizedExperiment", "limma"))

  if (class(se) != "SummarizedExperiment") {
    stop("se must be an object of type 'SummarizedExperiment'.")
  }

  design_mat <- .make_design_mat(design, data = SummarizedExperiment::colData(se))

  design_list <- .preprocess_design_matrix_and_var_names(design_mat, variable_name)
  design_mat <- design_list$design_mat
  variable_name <- design_list$variable_name
  selected_vars <- Matrix::colSums(design_mat > 0) >= min_covered_samples

  E <- SummarizedExperiment::assays(se)[[slot_E]]

  if (!is.null(W_mat)) {
    W <- W_mat
  } else if (!is.null(slot_V)) {
    V <- SummarizedExperiment::assays(se)[[slot_V]]
    W <- 1 / V
  } else {
    W <- NULL
  }

  if (!is.null(W)) {
    W_masked <- Matrix::t(apply(W, 1, function(w) {
      mask <- (!is.na(w) & is.finite(w) & (w != 0))

      upper <- stats::median(w[mask] + 3 * stats::mad(w[mask]))
      w[w > upper] <- 0

      lower <- stats::median(w[mask] - 3 * stats::mad(w[mask]))
      w[w < lower] <- 0

      w[!mask] <- 0
      return(w)
    }))
    W_masked[W_masked == 0] <- 1e-16

    fil_dm <- apply(design_mat[, selected_vars, drop = FALSE], 2, function(x) {
      mm <- (x > 0)
      v <- as.numeric((Matrix::rowSums(W_masked[, mm, drop = FALSE]) == 0) > 0)
      return(v)
    })

    fil_idx <- which(Matrix::rowSums(fil_dm) > 0)

    selected_feats <- setdiff(1:NROW(W_masked), fil_idx)
    lm_weights <- W_masked[selected_feats, ]
  } else {
    selected_feats <- 1:NROW(se)
    lm_weights <- NULL
  }

  suppressWarnings({
    fit <- limma::lmFit(
      object = E[selected_feats, ],
      design = design_mat[, selected_vars, drop = FALSE],
      weights = lm_weights
    )
  })

  suppressWarnings({
    contrast_mat <- limma::makeContrasts(contrasts = variable_name, levels = design_mat)
  })

  cfit <- limma::contrasts.fit(fit, contrast_mat[selected_vars])
  suppressWarnings(efit <- limma::eBayes(cfit, trend = FALSE, robust = TRUE))
  tbl <- limma::topTable(efit, number = Inf, sort.by = "none")
  rownames(tbl) <- rownames(se)[selected_feats]

  return(tbl)
}

# #' @export
# get.pseudobulk.SE.old <- function(
#     ace,
#     sample_attr,
#     ensemble = FALSE,
#     bins = 20,
#     assay = "counts",
#     col_data = NULL,
#     pseudocount = 0,
#     with_S = FALSE,
#     with_E = FALSE,
#     with_V = FALSE,
#     min_cells_per_batch = 3,
#     BPPARAM = BiocParallel::SerialParam()) {
#   IDX <- ACTIONetExperiment::get.data.or.split(ace, attr = sample_attr, to_return = "split")
#   good_batches <- sapply(IDX, length) >= min_cells_per_batch

#   if (!any(good_batches)) {
#     msg <- sprintf("No samples remaining.")
#     warning(msg)
#     return(NULL)
#   } else if (!all(good_batches)) {
#     old_batches <- names(IDX)
#     ace <- ace[, ace[[sample_attr]] %in% names(good_batches[good_batches])]
#     IDX <- ACTIONetExperiment::get.data.or.split(ace, attr = sample_attr, to_return = "split")
#     bad_batch_names <- setdiff(old_batches, names(IDX))
#     msg <- sprintf("Samples Dropped: %s\n", paste0(bad_batch_names, collapse = ", "))
#     message(msg)
#   }

#   counts_mat <- SummarizedExperiment::assays(ace)[[assay]]
#   sample_names <- names(IDX)
#   counts_list <- bplapply(IDX, function(idx) counts_mat[, idx, drop = FALSE], BPPARAM = BPPARAM)

#   se_assays <- list()

#   S0 <- do.call(cbind, bplapply(counts_list, Matrix::rowSums, BPPARAM = BPPARAM)) + pseudocount
#   se_assays$counts <- S0

#   if (with_E == TRUE) {
#     E0 <- do.call(cbind, bplapply(counts_list, ACTIONetExperiment:::fastRowMeans, BPPARAM = BPPARAM))
#     se_assays$mean <- E0
#   }

#   if (with_V == TRUE) {
#     V0 <- do.call(cbind, bplapply(counts_list, MatrixGenerics::rowVars, BPPARAM = BPPARAM))
#     se_assays$var <- V0
#   }

#   if (ensemble == TRUE) {
#     if (!any(with_S, with_E, with_V)) {
#       err <- sprintf("No ensemble assays to make.\n")
#       stop(err)
#     }

#     mr_assays <- .make_ensemble_assays(
#       counts_list = counts_list,
#       bins = bins,
#       pseudocount = pseudocount,
#       with_S = with_S,
#       with_E = with_E,
#       with_V = with_V,
#       BPPARAM = BPPARAM
#     )
#     se_assays <- c(se_assays, mr_assays$assays)
#   }

#   n_cells <- sapply(counts_list, NCOL)
#   nnz_feat_mean <- sapply(counts_list, function(X) mean(Matrix::colSums(X > 0)))
#   cd <- data.frame(
#     n_cells = n_cells,
#     nnz_feat_mean = nnz_feat_mean,
#     sample = factor(sample_names)
#   )

#   if (!is.null(col_data)) {
#     md <- col_data[col_data[[sample_attr]] %in% sample_names, , drop = FALSE]
#     md <- md[match(sample_names, md[[sample_attr]]), , drop = FALSE]
#     cd <- data.frame(md, cd)
#   }
#   rownames(cd) <- sample_names
#   cd <- droplevels(cd)
#   se <- SummarizedExperiment::SummarizedExperiment(
#     assays = se_assays,
#     colData = cd,
#     rowData = SummarizedExperiment::rowData(ace)
#   )

#   if (ensemble == TRUE) {
#     S4Vectors::metadata(se)[["bins"]] <- bins
#   }

#   invisible(gc())
#   return(se)
# }
