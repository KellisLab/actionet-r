library(actionet)

make_test_counts <- function() {
  set.seed(1)
  x <- matrix(abs(rnorm(60)), nrow = 6, ncol = 10)
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  x
}

test_that("toAnnData preserves orientation and assays", {
  counts <- make_test_counts()

  adata <- toAnnData(counts)
  expect_s3_class(adata, "AbstractAnnData")
  expect_equal(dim(adata$X), c(ncol(counts), nrow(counts)))
  expect_equal(rownames(adata), colnames(counts))
  expect_equal(colnames(adata), rownames(counts))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = counts + 1
    )
  )
  adata_sce <- toAnnData(sce)
  # First assay lands in .X; the rest in .layers (default x_layer behavior).
  expect_equal(as.matrix(adata_sce$X), t(counts))
  expect_equal(adata_sce$uns[["X_name"]], "counts")
  expect_true("logcounts" %in% names(adata_sce$layers))
  expect_false("counts" %in% names(adata_sce$layers))
  expect_equal(as.matrix(adata_sce$layers[["logcounts"]]), t(counts + 1))
})

test_that("AnnData container path matches matrix/raw compute paths", {
  counts <- make_test_counts()
  adata <- toAnnData(counts)

  red <- reduceKernel(adata, k = 3, layer = NULL, verbose = FALSE)
  raw_red <- reduceKernel(counts, k = 3, verbose = FALSE, return_raw = TRUE)

  expect_true("action" %in% names(red$obsm))
  expect_true("action_U" %in% names(red$varm))
  expect_true("action_A" %in% names(red$varm))
  expect_true("action_B" %in% names(red$obsm))
  expect_true("action_params" %in% names(red$uns))
  # Both AnnData and bare-matrix paths now route through cells x genes C++ contract.
  # raw_red$S_r is cells x k; obsm[["action"]] is also cells x k.
  expect_equal(unname(as.matrix(red$obsm[["action"]])), unname(raw_red$S_r), tolerance = 1e-8)

  act <- runACTION(adata = red, k_min = 2, k_max = 3, thread_no = 1)
  expect_true(all(c("H_stacked", "H_merged", "C_stacked", "C_merged") %in% names(act$obsm)))
  expect_true("assigned_archetype" %in% names(act$obs))

  net <- buildNetwork(adata = act, thread_no = 1)
  expect_true("actionet" %in% names(net$obsp))
})

test_that("deprecated compatibility wrappers still forward", {
  counts <- make_test_counts()
  adata <- toAnnData(counts)
  adata$layers[["logcounts"]] <- adata$X
  adata$obs[["batch"]] <- rep(c("a", "b"), each = 5)

  expect_warning(
    reduceKernel(adata, k = 3, assay_name = "logcounts", verbose = FALSE),
    "deprecated"
  )

  expect_warning(
    filter.ace(adata, layer = NULL, min_feats_per_cell = 1),
    "deprecated"
  )
  filtered <- suppressWarnings(
    filter.ace(adata, layer = NULL, min_feats_per_cell = 1)
  )
  expect_s3_class(filtered, "AbstractAnnData")

  expect_warning(
    filter.ace.by.attr(adata, by = "batch", layer = NULL, min_feats_per_cell = 1),
    "deprecated"
  )
})

test_that("AnnData plotting dispatch works", {
  skip_if_not_installed("ggplot2")

  counts <- make_test_counts()
  adata <- toAnnData(counts)
  adata$obsm[["umap_2d_actionet"]] <- matrix(seq_len(ncol(counts) * 2), ncol(counts), 2)
  adata$obs[["label"]] <- rep(c("a", "b"), each = 5)

  p <- plot(adata, label_attr = "label", coordinate_attr = "umap_2d_actionet")
  expect_s3_class(p, "ggplot")
})

test_that("backed AnnData objects fail fast", {
  tmp <- tempfile(fileext = ".h5ad")
  counts <- make_test_counts()
  adata <- toAnnData(counts)
  anndataR::write_h5ad(adata, tmp)

  backed <- anndataR::read_h5ad(tmp, as = "HDF5AnnData", mode = "r")
  expect_error(
    toAnnData(backed),
    "Backed object support is not implemented"
  )
})

test_that("optional ACE conversion round-trips canonical keys", {
  skip_if_not_installed("ACTIONetExperiment")

  counts <- make_test_counts()
  adata <- toAnnData(counts)
  adata$obsm[["action"]] <- matrix(seq_len(ncol(counts) * 2), ncol(counts), 2)
  adata$varm[["archetype_feat_specificity_upper"]] <- matrix(seq_len(nrow(counts) * 2), nrow(counts), 2)
  adata$obsp[["actionet"]] <- Matrix::Diagonal(ncol(counts))
  adata$uns[["action_params"]] <- list(sigma = c(1, 2))

  ace <- toACTIONetExperiment(adata)
  roundtrip <- toAnnData(ace)

  expect_true("action" %in% names(roundtrip$obsm))
  expect_true("archetype_feat_specificity_upper" %in% names(roundtrip$varm))
  expect_true("actionet" %in% names(roundtrip$obsp))
  expect_equal(as.matrix(roundtrip$X), as.matrix(adata$X))
})

test_that("toAnnData from live ACE preserves orientation, assays, and metadata", {
  skip_if_not_installed("ACTIONetExperiment")

  counts <- make_test_counts()
  ace <- ACTIONetExperiment::ACTIONetExperiment(
    assays = list(counts = counts, logcounts = log1p(counts))
  )
  rownames(ace) <- rownames(counts)
  colnames(ace) <- colnames(counts)
  ACTIONetExperiment::colMaps(ace)[["action"]] <- matrix(seq_len(ncol(counts) * 2), ncol(counts), 2)
  ACTIONetExperiment::rowMaps(ace)[["arch_feat_spec"]] <- matrix(seq_len(nrow(counts) * 2), nrow(counts), 2)
  ACTIONetExperiment::colNets(ace)[["ACTIONet"]] <- Matrix::Diagonal(ncol(counts))
  S4Vectors::metadata(ace)[["action_sigma"]] <- c(1, 2)
  S4Vectors::metadata(ace)[["foo_sigma"]] <- 3.14

  adata <- toAnnData(ace)
  expect_s3_class(adata, "AbstractAnnData")
  # Orientation: obs = cells, var = genes.
  expect_equal(dim(adata), c(ncol(counts), nrow(counts)))
  expect_equal(rownames(adata), colnames(counts))
  expect_equal(colnames(adata), rownames(counts))
  # First assay ("counts") lands in .X; the rest in layers.
  expect_equal(as.matrix(adata$X), t(counts))
  expect_true("logcounts" %in% names(adata$layers))
  expect_equal(as.matrix(adata$layers[["logcounts"]]), t(log1p(counts)))
  # X_name recorded for round-trip.
  expect_equal(adata$uns[["X_name"]], "counts")
  # Slots are canonicalized (alias -> canonical key).
  expect_true("action" %in% names(adata$obsm))
  expect_true("archetype_feat_specificity_upper" %in% names(adata$varm))
  expect_true("actionet" %in% names(adata$obsp))
  # _sigma folding: forward-side.
  expect_true("action_params" %in% names(adata$uns))
  expect_equal(adata$uns[["action_params"]][["sigma"]], c(1, 2))
  expect_true("foo_params" %in% names(adata$uns))
  expect_equal(adata$uns[["foo_params"]][["sigma"]], 3.14)
})

test_that("toAnnData x_layer selects which assay lands in .X", {
  skip_if_not_installed("ACTIONetExperiment")

  counts <- make_test_counts()
  ace <- ACTIONetExperiment::ACTIONetExperiment(
    assays = list(counts = counts, logcounts = log1p(counts))
  )
  rownames(ace) <- rownames(counts)
  colnames(ace) <- colnames(counts)

  adata <- toAnnData(ace, x_layer = "logcounts")
  expect_equal(as.matrix(adata$X), t(log1p(counts)))
  expect_true("counts" %in% names(adata$layers))
  expect_false("logcounts" %in% names(adata$layers))
  expect_equal(adata$uns[["X_name"]], "logcounts")

  expect_error(toAnnData(ace, x_layer = "nope"), "'x_layer' = 'nope'")
})

test_that("colMaps element type is stable across round-trip", {
  skip_if_not_installed("ACTIONetExperiment")

  counts <- make_test_counts()
  ace <- ACTIONetExperiment::ACTIONetExperiment(assays = list(counts = counts))
  rownames(ace) <- rownames(counts)
  colnames(ace) <- colnames(counts)
  m <- matrix(seq_len(ncol(counts) * 3), ncol(counts), 3)
  # Use the legacy ACE key ("ACTION"); toAnnData() canonicalizes to "action".
  ACTIONetExperiment::colMaps(ace)[["ACTION"]] <- m

  adata <- toAnnData(ace)
  # No SE-wrapping drift: obsm entry is a plain matrix, not a SummarizedExperiment.
  expect_true("action" %in% names(adata$obsm))
  expect_false(inherits(adata$obsm[["action"]], "SummarizedExperiment"))
  expect_equal(unname(as.matrix(adata$obsm[["action"]])), unname(m))

  ace2 <- toACTIONetExperiment(adata)
  # Reverse-mapped key back to the legacy alias.
  expect_true("ACTION" %in% names(ACTIONetExperiment::colMaps(ace2)))
  cm2 <- ACTIONetExperiment::colMaps(ace2)[["ACTION"]]
  # Regression: previous double-population yielded SE-wrapped entries here.
  expect_false(inherits(cm2, "SummarizedExperiment"))
  expect_equal(unname(as.matrix(cm2)), unname(m))
})

test_that("_sigma metadata round-trips symmetrically", {
  skip_if_not_installed("ACTIONetExperiment")

  counts <- make_test_counts()
  ace <- ACTIONetExperiment::ACTIONetExperiment(assays = list(counts = counts))
  rownames(ace) <- rownames(counts)
  colnames(ace) <- colnames(counts)
  S4Vectors::metadata(ace)[["action_sigma"]] <- c(1, 2)
  S4Vectors::metadata(ace)[["foo_sigma"]] <- 3.14
  S4Vectors::metadata(ace)[["plain_key"]] <- "hello"

  adata <- toAnnData(ace)
  ace2 <- toACTIONetExperiment(adata)
  meta <- S4Vectors::metadata(ace2)

  expect_equal(meta[["action_sigma"]], c(1, 2))
  expect_equal(meta[["foo_sigma"]], 3.14)
  expect_equal(meta[["plain_key"]], "hello")
})

test_that("bare SummarizedExperiment converts to AnnData", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = make_test_counts())
  )
  adata <- toAnnData(se)
  expect_s3_class(adata, "AbstractAnnData")
  expect_equal(dim(adata), c(ncol(se), nrow(se)))
  expect_equal(as.matrix(adata$X), t(SummarizedExperiment::assay(se, "counts")))
  expect_equal(adata$uns[["X_name"]], "counts")
})

test_that("matrix path with duplicate rownames warns count-only and dedupes", {
  x <- matrix(1, nrow = 3, ncol = 4)
  rownames(x) <- c("g1", "g1", "g1")
  colnames(x) <- c("c1", "c2", "c3", "c4")

  expect_warning(
    adata <- toAnnData(x),
    regexp = "Duplicated var_names detected \\(2 duplicates\\)"
  )
  # Message must not enumerate the offending names.
  w <- tryCatch(toAnnData(x), warning = function(w) conditionMessage(w))
  expect_false(any(grepl("g1(_1|_2|,)", w)))
  # Result should have deduped var names.
  expect_true(all(!duplicated(colnames(adata))))
})
