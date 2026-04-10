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
  expect_equal(as.matrix(adata_sce$X), t(counts))
  expect_true(all(c("counts", "logcounts") %in% names(adata_sce$layers)))
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
