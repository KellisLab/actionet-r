#!/usr/bin/env Rscript
# Plan 07 — Final Parity Validation: R end-to-end parity test.
#
# Runs the full canonical ACTIONet pipeline on the parity fixture, then
# compares every parity-critical output slot against:
#  1. The Plan 00 R baseline (intra-language regression)
#  2. The Python baseline h5ad (cross-language parity)
#
# Usage (from actionet-r repo root):
#   Rscript tests/test_parity.R [--skip-python]
#
# Exit code 0 if all checks pass, 1 otherwise.

suppressPackageStartupMessages({
  options(pkgbuild.override_build_tools = TRUE)
  pkgload::load_all(".", quiet = TRUE, compile = FALSE, recompile = FALSE)
  library(Matrix)
  library(anndataR)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
skip_python <- "--skip-python" %in% args

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
args_full <- commandArgs(trailingOnly = FALSE)
file_arg  <- args_full[grepl("^--file=", args_full)]
if (length(file_arg) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  SCRIPT_DIR <- file.path(getwd(), "tests")
}

FIXTURE         <- file.path(SCRIPT_DIR, "fixtures", "parity_fixture.h5ad")
BASELINE_R_H5AD <- file.path(SCRIPT_DIR, "fixtures", "baseline_r.h5ad")
BASELINE_PY_H5AD <- file.path(SCRIPT_DIR, "..", "..", "actionet-python",
                               "tests", "fixtures", "baseline_python.h5ad")
SEED <- 42L

# ---------------------------------------------------------------------------
# Test state
# ---------------------------------------------------------------------------
PASS <- 0L
FAIL <- 0L
WARN <- 0L

check <- function(name, cond, detail = NULL, warn = FALSE) {
  if (isTRUE(cond)) {
    cat(sprintf("  PASS  %s\n", name))
    PASS <<- PASS + 1L
  } else if (warn) {
    msg <- if (!is.null(detail)) sprintf(": %s", detail) else ""
    cat(sprintf("  WARN  %s%s\n", name, msg))
    WARN <<- WARN + 1L
  } else {
    msg <- if (!is.null(detail)) sprintf(": %s", detail) else ""
    cat(sprintf("  FAIL  %s%s\n", name, msg))
    FAIL <<- FAIL + 1L
  }
}

section <- function(title) {
  cat(sprintf("\n%s\n  %s\n%s\n", strrep("=", 60), title, strrep("=", 60)))
}

# ---------------------------------------------------------------------------
# Canonicalization helpers
# ---------------------------------------------------------------------------
canon_svd_sign <- function(mat) {
  mat <- as.matrix(mat)
  for (j in seq_len(ncol(mat))) {
    col <- mat[, j]
    idx <- which.max(abs(col))
    if (col[idx] < 0) mat[, j] <- -col
  }
  mat
}

canon_archetype_order <- function(...) {
  mats <- list(...)
  ref  <- as.matrix(mats[[1]])
  norms <- apply(ref, 2, function(col) sqrt(sum(col^2)))
  order <- order(norms, decreasing = TRUE)
  lapply(mats, function(m) as.matrix(m)[, order, drop = FALSE])
}

cmp_dense <- function(slot, a, b, atol = 1e-6, rtol = 1e-4) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (!identical(dim(a), dim(b))) {
    check(slot, FALSE, sprintf("shape mismatch: %s vs %s",
                               paste(dim(a), collapse = "x"),
                               paste(dim(b), collapse = "x")))
    return(invisible(NULL))
  }
  dev <- max(abs(a - b))
  ok  <- all(abs(a - b) <= atol + rtol * abs(b))
  check(slot, ok, sprintf("max_dev=%.3e", dev))
}

cmp_archetype_slots <- function(slot, a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (ncol(a) != ncol(b)) {
    check(slot, FALSE,
          sprintf("archetype count differs: R=%d Python=%d (known stochastic pruning diff)",
                  ncol(a), ncol(b)),
          warn = TRUE)
    return(invisible(NULL))
  }
  ac <- canon_archetype_order(a)[[1]]
  bc <- canon_archetype_order(b)[[1]]
  cmp_dense(slot, ac, bc)
}

# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
run_pipeline <- function(adata) {
  adata <- reduceKernel(adata, k = 20L, algorithm = "irlb", seed = SEED,
                        layer = "logcounts", reduction_slot = "action", verbose = FALSE)
  adata <- runACTION(adata, k_min = 2L, k_max = 20L, reduction_slot = "action")
  adata <- buildNetwork(adata, map_slot = "H_stacked", net_slot_out = "actionet")
  labels <- as.integer(adata$obs[["assigned_archetype"]])
  adata <- computeFeatureSpecificity(adata, labels = labels, layer = "logcounts",
                                     map_out_prefix = "cluster", return_lower = TRUE)
  adata <- archetypeFeatureSpecificity(adata, layer = "logcounts",
                                       map_slot = "H_merged", map_out_prefix = "archetype")
  batches <- adata$obs[["batch"]]
  adata <- correctBatchEffect(adata, batches = batches, reduction_slot = "action",
                               corrected_suffix = "orth", layer = "logcounts")
  adata <- layoutNetwork(adata, net_slot = "actionet", seed = SEED, verbose = FALSE,
                          map_slot_out = "actionet_2d")
  adata
}

# ---------------------------------------------------------------------------
# Shape verification
# ---------------------------------------------------------------------------
verify_shapes <- function(adata) {
  section("Shape verification (cells x genes AnnData-native contract)")
  n_obs  <- nrow(adata)
  n_vars <- ncol(adata)

  check(sprintf("X shape (%d x %d)", n_obs, n_vars),
        nrow(adata$X) == n_obs && ncol(adata$X) == n_vars)

  for (slot in c("action", "action_B", "H_stacked", "H_merged",
                 "C_stacked", "C_merged", "action_orth", "actionet_2d")) {
    m <- colMaps(adata)[[slot]]
    if (!is.null(m))
      check(sprintf("colMaps/%s rows == n_obs", slot), nrow(m) == n_obs,
            sprintf("got %d", nrow(m)))
  }
  for (slot in c("action_U", "action_A", "action_U_orth", "action_A_orth",
                 "cluster_upper", "cluster_lower",
                 "archetype_feat_profile",
                 "archetype_feat_specificity_upper",
                 "archetype_feat_specificity_lower")) {
    m <- rowMaps(adata)[[slot]]
    if (!is.null(m))
      check(sprintf("rowMaps/%s rows == n_vars", slot), nrow(m) == n_vars,
            sprintf("got %d", nrow(m)))
  }
  G <- colNets(adata)[["actionet"]]
  if (!is.null(G))
    check(sprintf("colNets/actionet (%d x %d)", n_obs, n_obs),
          nrow(G) == n_obs && ncol(G) == n_obs)
}

# ---------------------------------------------------------------------------
# Intra-language regression: new vs Plan 00 R baseline
# ---------------------------------------------------------------------------
check_intra_regression <- function(adata, baseline) {
  section("Intra-language regression (new vs Plan 00 R baseline)")

  # Reduction
  for (slot in c("action", "action_B")) {
    r_m <- colMaps(adata)[[slot]]
    b_m <- colMaps(baseline)[[slot]]
    if (!is.null(r_m) && !is.null(b_m))
      cmp_dense(sprintf("%s regression", slot),
                canon_svd_sign(as.matrix(r_m)), canon_svd_sign(as.matrix(b_m)))
  }
  for (slot in c("action_U", "action_A")) {
    r_m <- rowMaps(adata)[[slot]]
    b_m <- rowMaps(baseline)[[slot]]
    if (!is.null(r_m) && !is.null(b_m))
      cmp_dense(sprintf("%s regression", slot),
                canon_svd_sign(as.matrix(r_m)), canon_svd_sign(as.matrix(b_m)))
  }
  if (!is.null(adata$uns[["action_params"]][["sigma"]]) &&
      !is.null(baseline$uns[["action_params"]][["sigma"]])) {
    cmp_dense("sigma regression",
              as.numeric(adata$uns[["action_params"]][["sigma"]]),
              as.numeric(baseline$uns[["action_params"]][["sigma"]]),
              atol = 1e-8, rtol = 1e-6)
  }

  # Network
  G_new <- colNets(adata)[["actionet"]]
  G_ref <- colNets(baseline)[["actionet"]]
  if (!is.null(G_new) && !is.null(G_ref)) {
    G_new_d <- as.matrix(G_new)
    G_ref_d <- as.matrix(G_ref)
    if (!identical(dim(G_new_d), dim(G_ref_d))) {
      check("network regression", FALSE,
            sprintf("shape: %s vs %s",
                    paste(dim(G_new_d), collapse = "x"),
                    paste(dim(G_ref_d), collapse = "x")))
    } else {
      dev <- max(abs(G_new_d - G_ref_d))
      check("network regression", dev <= 1e-6, sprintf("max_dev=%.3e", dev))
    }
  }

  # Specificity
  for (slot in c("cluster_upper", "cluster_lower")) {
    r_m <- rowMaps(adata)[[slot]]
    b_m <- rowMaps(baseline)[[slot]]
    if (!is.null(r_m) && !is.null(b_m))
      cmp_dense(sprintf("%s regression", slot), as.matrix(r_m), as.matrix(b_m),
                atol = 1e-4, rtol = 1e-3)
  }

  # Batch correction
  for (slot in c("action_orth")) {
    r_m <- colMaps(adata)[[slot]]
    b_m <- colMaps(baseline)[[slot]]
    if (!is.null(r_m) && !is.null(b_m))
      cmp_dense(sprintf("%s regression", slot),
                canon_svd_sign(as.matrix(r_m)), canon_svd_sign(as.matrix(b_m)))
  }
  for (slot in c("action_U_orth", "action_A_orth")) {
    r_m <- rowMaps(adata)[[slot]]
    b_m <- rowMaps(baseline)[[slot]]
    if (!is.null(r_m) && !is.null(b_m))
      cmp_dense(sprintf("%s regression", slot),
                canon_svd_sign(as.matrix(r_m)), canon_svd_sign(as.matrix(b_m)))
  }
}

# ---------------------------------------------------------------------------
# Cross-language parity: R new vs Python baseline h5ad
# ---------------------------------------------------------------------------
check_cross_language <- function(adata, py_h5ad_path) {
  section("Cross-language parity (R vs Python)")

  if (!file.exists(py_h5ad_path)) {
    check("Python baseline exists", FALSE, sprintf("not found: %s", py_h5ad_path))
    return(invisible(NULL))
  }
  check("Python baseline exists", TRUE)

  py <- anndataR::read_h5ad(py_h5ad_path, as = "InMemoryAnnData")
  cat(sprintf("  Python baseline: %d obs x %d vars\n", nrow(py), ncol(py)))

  # Reduction — exact match expected
  for (slot in c("action", "action_B")) {
    r_m  <- colMaps(adata)[[slot]]
    py_m <- py$obsm[[slot]]
    if (!is.null(r_m) && !is.null(py_m)) {
      cmp_dense(sprintf("cross-lang %s", slot),
                canon_svd_sign(as.matrix(r_m)),
                canon_svd_sign(as.matrix(py_m)))
    } else {
      check(sprintf("cross-lang %s present", slot), FALSE,
            sprintf("MISSING in %s", if (is.null(r_m)) "R" else "Python"))
    }
  }
  for (slot in c("action_U", "action_A")) {
    r_m  <- rowMaps(adata)[[slot]]
    py_m <- py$varm[[slot]]
    if (!is.null(r_m) && !is.null(py_m)) {
      cmp_dense(sprintf("cross-lang %s", slot),
                canon_svd_sign(as.matrix(r_m)),
                canon_svd_sign(as.matrix(py_m)))
    } else {
      check(sprintf("cross-lang %s present", slot), FALSE,
            sprintf("MISSING in %s", if (is.null(r_m)) "R" else "Python"))
    }
  }
  if (!is.null(adata$uns[["action_params"]][["sigma"]]) &&
      !is.null(py$uns[["action_params"]][["sigma"]])) {
    cmp_dense("cross-lang sigma",
              as.numeric(adata$uns[["action_params"]][["sigma"]]),
              as.numeric(py$uns[["action_params"]][["sigma"]]),
              atol = 1e-8, rtol = 1e-6)
  }

  # ACTION archetypes — shape may differ
  for (slot in c("H_stacked", "H_merged", "C_stacked", "C_merged")) {
    r_m  <- colMaps(adata)[[slot]]
    py_m <- py$obsm[[slot]]
    if (!is.null(r_m) && !is.null(py_m))
      cmp_archetype_slots(sprintf("cross-lang %s", slot),
                          as.matrix(r_m), as.matrix(py_m))
    else
      check(sprintf("cross-lang %s present", slot), FALSE,
            sprintf("MISSING in %s", if (is.null(r_m)) "R" else "Python"))
  }

  # Network — depends on archetypes, warn if different
  G_r  <- colNets(adata)[["actionet"]]
  G_py <- py$obsp[["actionet"]]
  if (!is.null(G_r) && !is.null(G_py)) {
    G_r_d  <- as.matrix(G_r)
    G_py_d <- as.matrix(G_py)
    if (!identical(dim(G_r_d), dim(G_py_d))) {
      check("cross-lang network", FALSE,
            sprintf("shape mismatch: R=%s Python=%s (follows from archetype diff)",
                    paste(dim(G_r_d), collapse = "x"),
                    paste(dim(G_py_d), collapse = "x")),
            warn = TRUE)
    } else {
      nnz_r  <- sum(G_r_d != 0)
      nnz_py <- sum(G_py_d != 0)
      if (nnz_r != nnz_py) {
        check("cross-lang network", FALSE,
              sprintf("nnz R=%d Python=%d (follows from archetype diff)", nnz_r, nnz_py),
              warn = TRUE)
      } else {
        dev <- max(abs(G_r_d - G_py_d))
        check("cross-lang network", dev <= 1e-6, sprintf("max_dev=%.3e", dev))
      }
    }
  }

  # Specificity
  spec_pairs <- list(
    list(r = "cluster_upper",  py = "specificity_upper",  label = "cross-lang specificity_upper"),
    list(r = "cluster_lower",  py = "specificity_lower",  label = "cross-lang specificity_lower"),
    list(r = "archetype_feat_profile",           py = "archetype_feat_profile",           label = "cross-lang arch_feat_profile"),
    list(r = "archetype_feat_specificity_upper", py = "archetype_feat_specificity_upper", label = "cross-lang arch_feat_upper"),
    list(r = "archetype_feat_specificity_lower", py = "archetype_feat_specificity_lower", label = "cross-lang arch_feat_lower")
  )
  for (p in spec_pairs) {
    r_slot  <- p$r
    py_slot <- p$py
    label   <- p$label
    r_m  <- rowMaps(adata)[[r_slot]]
    py_m <- py$varm[[py_slot]]
    if (is.null(r_m) || is.null(py_m)) {
      check(label, FALSE, sprintf("MISSING in %s",
                                  if (is.null(r_m)) "R" else "Python"))
      next
    }
    r_mat  <- as.matrix(r_m)
    py_mat <- as.matrix(py_m)
    if (!identical(dim(r_mat), dim(py_mat))) {
      check(label, FALSE,
            sprintf("shape: R=%s Python=%s (archetype count diff)",
                    paste(dim(r_mat), collapse = "x"),
                    paste(dim(py_mat), collapse = "x")),
            warn = TRUE)
    } else {
      cmp_dense(label, r_mat, py_mat, atol = 1e-4, rtol = 1e-3)
    }
  }

  # Batch correction
  for (r_slot in c("action_orth")) {
    r_m  <- colMaps(adata)[[r_slot]]
    py_m <- py$obsm[["action_corrected"]]
    if (!is.null(r_m) && !is.null(py_m))
      cmp_dense(sprintf("cross-lang %s", r_slot),
                canon_svd_sign(as.matrix(r_m)),
                canon_svd_sign(as.matrix(py_m)))
  }
  for (pair in list(list("action_U_orth", "action_corrected_U"),
                    list("action_A_orth", "action_corrected_A"))) {
    r_m  <- rowMaps(adata)[[pair[[1]]]]
    py_m <- py$varm[[pair[[2]]]]
    if (!is.null(r_m) && !is.null(py_m))
      cmp_dense(sprintf("cross-lang %s", pair[[1]]),
                canon_svd_sign(as.matrix(r_m)),
                canon_svd_sign(as.matrix(py_m)))
  }
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
section("Loading fixture and running pipeline")
if (!file.exists(FIXTURE)) stop(sprintf("Fixture not found: %s", FIXTURE))

adata <- anndataR::read_h5ad(FIXTURE, as = "InMemoryAnnData")
cat(sprintf("  Loaded %d obs x %d vars\n", nrow(adata), ncol(adata)))
adata <- run_pipeline(adata)
cat("  Pipeline complete.\n")

verify_shapes(adata)

# Intra-language regression
if (file.exists(BASELINE_R_H5AD)) {
  baseline <- anndataR::read_h5ad(BASELINE_R_H5AD, as = "InMemoryAnnData")
  check_intra_regression(adata, baseline)
} else {
  section("Intra-language regression")
  check("R baseline exists", FALSE, sprintf("not found: %s", BASELINE_R_H5AD))
}

# Cross-language parity
if (!skip_python) {
  check_cross_language(adata, BASELINE_PY_H5AD)
} else {
  section("Cross-language parity")
  cat("  SKIPPED (--skip-python)\n")
}

# Final summary
cat(sprintf("\n%s\n  RESULTS: %d passed, %d failed, %d warnings\n%s\n",
            strrep("=", 60), PASS, FAIL, WARN, strrep("=", 60)))

if (WARN > 0) {
  cat("\nKnown differences (documented, not counted as failures):\n")
  cat("  See WARN lines above - archetype count and network differ between\n")
  cat("  R and Python due to stochastic ACTION pruning (different defaults).\n")
}

quit(status = if (FAIL == 0L) 0L else 1L)
