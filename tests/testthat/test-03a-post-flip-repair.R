library(actionet)

# ── Helpers ────────────────────────────────────────────────────────────────────

make_adata <- function(n_cells = 12, n_genes = 8, seed = 42) {
  set.seed(seed)
  mat <- matrix(abs(rnorm(n_cells * n_genes)), nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("g", seq_len(n_genes))
  colnames(mat) <- paste0("c", seq_len(n_cells))
  adata <- toAnnData(mat)
  adata$layers[["logcounts"]] <- adata$X
  adata
}

make_full_adata <- function(n_cells = 15, n_genes = 10, seed = 7) {
  adata <- make_adata(n_cells = n_cells, n_genes = n_genes, seed = seed)
  adata <- reduceKernel(adata, k = 4, layer = "logcounts", verbose = FALSE)
  adata <- runACTION(adata = adata, k_min = 2, k_max = 3, thread_no = 1)
  adata <- buildNetwork(adata = adata, thread_no = 1)
  adata
}

# ── 1. layoutNetwork() regression: must use svd$u (cell-space) ────────────────

test_that("layoutNetwork uses cell-space SVD output (svd$u) for default init", {
  adata <- make_full_adata()

  # Must not error; previous bug would pass svd$v (genes x k) which fails the
  # NROW(initial_coordinates) == n_obs(adata) shape guard.
  result <- layoutNetwork(adata, n_components = 2, verbose = FALSE)

  expect_s3_class(result, "AbstractAnnData")
  embedding_key <- grep("umap_2d", names(result$obsm), value = TRUE)[[1]]
  emb <- result$obsm[[embedding_key]]
  expect_equal(nrow(emb), nrow(adata))
  expect_equal(ncol(emb), 2L)
})

test_that("layoutNetwork accepts explicit initial_coordinates as before", {
  adata <- make_full_adata()
  init <- matrix(rnorm(nrow(adata) * 2), nrow = nrow(adata), ncol = 2)

  result <- layoutNetwork(adata, initial_coordinates = init,
                          n_components = 2, verbose = FALSE)
  expect_s3_class(result, "AbstractAnnData")
})

# ── 2. filterActionet() regression: cells-on-rows orientation ─────────────────

test_that("filterActionet removes cells below UMI threshold (rows in cells x genes)", {
  set.seed(1)
  n_cells <- 10
  n_genes <- 6
  # Build a matrix where cell 1 has very low total counts
  mat <- matrix(abs(rnorm(n_cells * n_genes)) + 1, nrow = n_genes, ncol = n_cells)
  mat[, 1] <- 0.01  # cell 1 has near-zero counts in all genes
  rownames(mat) <- paste0("g", seq_len(n_genes))
  colnames(mat) <- paste0("c", seq_len(n_cells))
  adata <- toAnnData(mat)
  adata$layers[["logcounts"]] <- adata$X

  # With threshold just above 0.01*n_genes, cell 1 (c1) should be removed.
  threshold <- 0.05 * n_genes
  filtered <- filterActionet(adata, layer = "logcounts",
                             min_umis_per_cell = threshold)

  expect_s3_class(filtered, "AbstractAnnData")
  expect_lt(nrow(filtered), nrow(adata))      # fewer cells
  expect_false("c1" %in% rownames(filtered))  # cell c1 removed
  expect_equal(ncol(filtered), ncol(adata))   # all genes retained
})

test_that("filterActionet removes features below cell-count threshold (cols in cells x genes)", {
  set.seed(2)
  n_cells <- 10
  n_genes <- 8
  mat <- matrix(abs(rnorm(n_cells * n_genes)) + 1, nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("g", seq_len(n_genes))
  colnames(mat) <- paste0("c", seq_len(n_cells))
  # gene g1 is expressed in only 1 cell
  mat["g1", ] <- 0
  mat["g1", "c1"] <- 5
  adata <- toAnnData(mat)
  adata$layers[["logcounts"]] <- adata$X

  # Require a gene to appear in >= 3 cells; g1 appears in only 1 cell.
  filtered <- filterActionet(adata, layer = "logcounts",
                             min_cells_per_feat = 3)

  expect_s3_class(filtered, "AbstractAnnData")
  expect_lt(ncol(filtered), ncol(adata))      # fewer genes
  expect_false("g1" %in% colnames(filtered))  # gene g1 removed
  expect_equal(nrow(filtered), nrow(adata))   # all cells retained
})

test_that("filter.ace deprecated wrapper still works", {
  adata <- make_adata()
  expect_warning(
    result <- filter.ace(adata, layer = NULL, min_feats_per_cell = 1),
    "deprecated"
  )
  expect_s3_class(result, "AbstractAnnData")
})

# ── 3. filterActionetByAttr() regression ─────────────────────────────────────

test_that("filterActionetByAttr groups and filters correctly", {
  set.seed(3)
  n_cells <- 12
  n_genes <- 6
  mat <- matrix(abs(rnorm(n_cells * n_genes)) + 1, nrow = n_genes, ncol = n_cells)
  # Make cell c1 (batch A) and cell c7 (batch B) have very low counts.
  mat[, 1] <- 0.01
  mat[, 7] <- 0.01
  rownames(mat) <- paste0("g", seq_len(n_genes))
  colnames(mat) <- paste0("c", seq_len(n_cells))
  adata <- toAnnData(mat)
  adata$layers[["logcounts"]] <- adata$X
  adata$obs[["batch"]] <- c(rep("A", 6), rep("B", 6))

  threshold <- 0.05 * n_genes
  filtered <- filterActionetByAttr(adata, by = "batch", layer = "logcounts",
                                   min_umis_per_cell = threshold)

  expect_s3_class(filtered, "AbstractAnnData")
  expect_lt(nrow(filtered), nrow(adata))
  # Both low-count cells from each batch should be removed.
  expect_false("c1" %in% rownames(filtered))
  expect_false("c7" %in% rownames(filtered))
})

# ── 4. Bare-matrix contract: legacy genes x cells ─────────────────────────────

test_that("bare matrix (genes x cells) and AnnData path produce identical reduction", {
  set.seed(5)
  counts <- matrix(abs(rnorm(80)), nrow = 8, ncol = 10)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  adata <- toAnnData(counts)

  # AnnData path
  red_adata <- reduceKernel(adata, k = 3, layer = NULL, verbose = FALSE)

  # Bare matrix path: counts is genes x cells (legacy).
  # After the fix, reduceKernel transposes internally, so S_r is cells x k.
  raw_red <- reduceKernel(counts, k = 3, verbose = FALSE, return_raw = TRUE)

  expect_equal(dim(raw_red$S_r), c(ncol(counts), 3L))
  expect_equal(dim(red_adata$obsm[["action"]]), c(nrow(adata), 3L))
  expect_equal(
    unname(as.matrix(red_adata$obsm[["action"]])),
    unname(raw_red$S_r),
    tolerance = 1e-8
  )
})
