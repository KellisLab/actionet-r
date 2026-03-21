test_that("runACTIONet forwards knn network M and k", {
  adata <- list(
    obsm = list(action = matrix(seq_len(6), nrow = 3, ncol = 2)),
    obs = list(),
    obsp = list(),
    varm = list()
  )
  forwarded <- new.env(parent = emptyenv())

  local_mocked_bindings(
    .resolve_container_arg = function(adata = NULL, ace = NULL) adata,
    .resolve_layer_arg = function(...) "logcounts",
    .validate_ace = function(adata, ...) adata,
    .validate_assay = function(...) invisible(NULL),
    .validate_map = function(...) invisible(NULL),
    colMaps = function(object, all = TRUE) object$obsm,
    runACTION = function(adata, ...) {
      H <- matrix(
        c(
          1.0, 0.0,
          0.8, 0.2,
          0.2, 0.8
        ),
        nrow = 3,
        byrow = TRUE
      )
      adata$obsm[["H_stacked"]] <- H
      adata$obsm[["H_merged"]] <- H
      adata$obs[["assigned_archetype"]] <- c("a", "b", "c")
      adata
    },
    buildNetwork = function(adata, ..., M, k, net_slot_out) {
      forwarded$M <- M
      forwarded$k <- k
      adata$obsp[[net_slot_out]] <- Matrix::Diagonal(3)
      adata
    },
    networkCentrality = function(adata, ..., attr_out) {
      adata$obs[[attr_out]] <- rep(1, 3)
      adata
    },
    networkDiffusion = function(adata, ..., map_slot_out) {
      adata$obsm[[map_slot_out]] <- matrix(
        c(
          1.0, 0.0,
          0.5, 0.5,
          0.0, 1.0
        ),
        nrow = 3,
        byrow = TRUE
      )
      adata
    },
    layoutNetwork = function(adata, ..., map_slot_out, n_components) {
      adata$obsm[[map_slot_out]] <- matrix(0, nrow = 3, ncol = n_components)
      adata
    },
    archetypeFeatureSpecificity = function(adata, ...) adata,
    .package = "actionet"
  )

  out <- runACTIONet(
    adata = adata,
    network_algorithm = "knn",
    network_M = 24,
    network_k = 7,
    layout_3d = FALSE,
    thread_no = 1
  )

  expect_equal(forwarded$M, 24)
  expect_equal(forwarded$k, 7)
  expect_true("actionet" %in% names(out$obsp))
})
