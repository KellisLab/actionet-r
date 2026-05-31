library(actionet)

make_layout_adata <- function(n_cells = 12L, n_genes = 6L, seed = 11L) {
  set.seed(seed)
  counts <- matrix(abs(rnorm(n_genes * n_cells)), nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("g", seq_len(n_genes))
  colnames(counts) <- paste0("c", seq_len(n_cells))
  adata <- toAnnData(counts)
  adata$layers[["logcounts"]] <- adata$X

  i <- c(seq_len(n_cells), seq_len(n_cells))
  j <- c((seq_len(n_cells) %% n_cells) + 1L, seq_len(n_cells))
  x <- rep(1, length(i))
  adata$obsp[["actionet"]] <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n_cells, n_cells))
  adata
}

make_init <- function(n_cells, n_dims = 3L, seed = 99L) {
  set.seed(seed)
  matrix(rnorm(n_cells * n_dims), nrow = n_cells, ncol = n_dims)
}

test_that("layoutNetwork deterministic rng_type is reproducible for the same seed", {
  adata <- make_layout_adata()
  init <- make_init(nrow(adata))

  emb_a <- layoutNetwork(
    adata = adata,
    initial_coordinates = init,
    method = "umap",
    n_components = 2,
    n_epochs = 20,
    rng_type = "deterministic",
    seed = 1,
    thread_no = 1,
    verbose = FALSE,
    return_raw = TRUE
  )

  emb_b <- layoutNetwork(
    adata = adata,
    initial_coordinates = init,
    method = "umap",
    n_components = 2,
    n_epochs = 20,
    rng_type = "deterministic",
    seed = 1,
    thread_no = 1,
    verbose = FALSE,
    return_raw = TRUE
  )

  expect_equal(emb_a, emb_b)
})

test_that("layoutNetwork supports method='largevis'", {
  adata <- make_layout_adata()
  init <- make_init(nrow(adata))

  emb <- layoutNetwork(
    adata = adata,
    initial_coordinates = init,
    method = "largevis",
    n_components = 2,
    n_epochs = 20,
    seed = 0,
    thread_no = 1,
    verbose = FALSE,
    return_raw = TRUE
  )

  expect_equal(dim(emb), c(nrow(adata), 2L))
})

test_that("layoutNetwork enforces ai/aj requirements for leopold methods", {
  adata <- make_layout_adata()
  init <- make_init(nrow(adata))
  ai <- seq(0.1, 1.0, length.out = nrow(adata))

  expect_error(
    layoutNetwork(
      adata = adata,
      initial_coordinates = init,
      method = "leopold",
      n_components = 2,
      n_epochs = 10,
      thread_no = 1,
      verbose = FALSE,
      return_raw = TRUE
    ),
    "ai"
  )

  expect_error(
    layoutNetwork(
      adata = adata,
      initial_coordinates = init,
      method = "leopold2",
      ai = ai,
      n_components = 2,
      n_epochs = 10,
      thread_no = 1,
      verbose = FALSE,
      return_raw = TRUE
    ),
    "ai.*aj"
  )
})

test_that("layoutNetwork supports method='leopold' and method='leopold2'", {
  adata <- make_layout_adata()
  init <- make_init(nrow(adata))
  ai <- seq(0.1, 1.0, length.out = nrow(adata))
  aj <- seq(0.2, 1.1, length.out = nrow(adata))

  emb_l1 <- layoutNetwork(
    adata = adata,
    initial_coordinates = init,
    method = "leopold",
    ai = ai,
    n_components = 2,
    n_epochs = 20,
    seed = 0,
    thread_no = 1,
    verbose = FALSE,
    return_raw = TRUE
  )
  expect_equal(dim(emb_l1), c(nrow(adata), 2L))

  emb_l2 <- layoutNetwork(
    adata = adata,
    initial_coordinates = init,
    method = "leopold2",
    ai = ai,
    aj = aj,
    n_components = 2,
    n_epochs = 20,
    seed = 0,
    thread_no = 1,
    verbose = FALSE,
    return_raw = TRUE
  )
  expect_equal(dim(emb_l2), c(nrow(adata), 2L))
})

# ---------------------------------------------------------------------------
# Regression tests for the canonical disconnected-vertex repair
# (mirror of `actionet-python/tests/test_layout_no_axis_artifact.py`).
#
# Without the repair, fully isolated vertices receive zero updates from the
# batch optimizer and remain frozen at their seed coordinate; with the repair
# (canonical umap-learn `simplicial_set_embedding` behavior) they get a
# per-component offset and per-coordinate jitter that breaks any axis-aligned
# seed structure.
# ---------------------------------------------------------------------------

make_disconnected_adata <- function(n_main = 80L, n_orphans = 30L, seed = 11L) {
  n_obs <- n_main + n_orphans
  set.seed(seed)
  counts <- matrix(abs(rnorm(4L * n_obs)), nrow = 4L, ncol = n_obs)
  rownames(counts) <- paste0("g", seq_len(4L))
  colnames(counts) <- paste0("c", seq_len(n_obs))
  adata <- toAnnData(counts)
  adata$layers[["logcounts"]] <- adata$X

  # Ring on the first n_main vertices; the remaining n_orphans vertices are
  # fully isolated (no edges).
  i <- c(seq_len(n_main), seq_len(n_main))
  j <- c((seq_len(n_main) %% n_main) + 1L, seq_len(n_main))
  x <- rep(1, length(i))
  adata$obsp[["actionet"]] <- Matrix::sparseMatrix(
    i = i, j = j, x = x, dims = c(n_obs, n_obs)
  )
  list(adata = adata, n_main = n_main, n_orphans = n_orphans)
}

axis_aligned_init <- function(n_obs, n_components = 3L) {
  init <- matrix(0, nrow = n_obs, ncol = n_components)
  init[, 1] <- seq(-1, 1, length.out = n_obs)
  init
}

test_that("layoutNetwork repair_disconnected moves orphans off-axis", {
  built <- make_disconnected_adata(n_main = 80L, n_orphans = 30L)
  adata <- built$adata
  init <- axis_aligned_init(nrow(adata), n_components = 3L)

  emb <- layoutNetwork(
    adata = adata,
    initial_coordinates = init,
    method = "umap",
    n_components = 2,
    n_epochs = 50,
    seed = 0,
    thread_no = 1,
    verbose = FALSE,
    repair_disconnected = TRUE,
    return_raw = TRUE
  )

  orphan_idx <- (built$n_main + 1L):nrow(emb)
  orphan_y <- emb[orphan_idx, 2]

  # No orphan should remain on the seed axis (Y == 0) after repair.
  on_axis <- abs(orphan_y) < 1e-3
  expect_false(any(on_axis))

  # And the overall axis-aligned fraction across all cells should be small.
  on_x_axis <- abs(emb[, 2]) < 1e-3
  on_y_axis <- abs(emb[, 1]) < 1e-3
  expect_lt(mean(on_x_axis | on_y_axis), 0.05)
})

test_that("layoutNetwork repair_disconnected = FALSE keeps orphans frozen at seed", {
  built <- make_disconnected_adata(n_main = 80L, n_orphans = 30L)
  adata <- built$adata
  init <- axis_aligned_init(nrow(adata), n_components = 3L)

  emb <- layoutNetwork(
    adata = adata,
    initial_coordinates = init,
    method = "umap",
    n_components = 2,
    n_epochs = 50,
    seed = 0,
    thread_no = 1,
    verbose = FALSE,
    repair_disconnected = FALSE,
    return_raw = TRUE
  )

  orphan_idx <- (built$n_main + 1L):nrow(emb)
  orphan_coords <- emb[orphan_idx, , drop = FALSE]
  seed_orphans <- init[orphan_idx, 1:2, drop = FALSE]

  # Orphans must match their seed up to float32 precision (the optimizer
  # rounds initial coords through `fmat`).
  expect_equal(orphan_coords, seed_orphans, tolerance = 1e-6)

  # Consequently, every orphan still lies on the seed X axis.
  expect_true(all(abs(orphan_coords[, 2]) < 1e-6))
})

test_that("layoutNetwork repair_disconnected jitters duplicate-seed single component", {
  # Fully connected ring (single component), but every vertex starts at the
  # same coordinate. Without jitter the optimizer keeps them collapsed; with
  # jitter the symmetry is broken and the points spread.
  n_obs <- 40L
  set.seed(7L)
  counts <- matrix(abs(rnorm(4L * n_obs)), nrow = 4L, ncol = n_obs)
  rownames(counts) <- paste0("g", seq_len(4L))
  colnames(counts) <- paste0("c", seq_len(n_obs))
  adata <- toAnnData(counts)
  adata$layers[["logcounts"]] <- adata$X

  i <- c(seq_len(n_obs), seq_len(n_obs))
  j <- c((seq_len(n_obs) %% n_obs) + 1L, seq_len(n_obs))
  x <- rep(1, length(i))
  adata$obsp[["actionet"]] <- Matrix::sparseMatrix(
    i = i, j = j, x = x, dims = c(n_obs, n_obs)
  )

  init <- matrix(rep(c(0.5, 0.0, 0.0), each = n_obs), nrow = n_obs, ncol = 3L)

  emb_repair <- layoutNetwork(
    adata = adata,
    initial_coordinates = init,
    method = "umap",
    n_components = 2,
    n_epochs = 50,
    seed = 0,
    thread_no = 1,
    verbose = FALSE,
    repair_disconnected = TRUE,
    return_raw = TRUE
  )

  emb_no_repair <- layoutNetwork(
    adata = adata,
    initial_coordinates = init,
    method = "umap",
    n_components = 2,
    n_epochs = 50,
    seed = 0,
    thread_no = 1,
    verbose = FALSE,
    repair_disconnected = FALSE,
    return_raw = TRUE
  )

  centered_repair <- sweep(emb_repair, 2, colMeans(emb_repair))
  centered_no_repair <- sweep(emb_no_repair, 2, colMeans(emb_no_repair))
  spread_repair <- sd(sqrt(rowSums(centered_repair^2)))
  spread_no_repair <- sd(sqrt(rowSums(centered_no_repair^2)))

  expect_gt(spread_repair, spread_no_repair)
})
