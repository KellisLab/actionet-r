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

test_that("layoutNetwork deterministic rng_type is reproducible across seeds", {
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
    seed = 999,
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
