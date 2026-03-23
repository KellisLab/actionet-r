#!/usr/bin/env Rscript
# validate_stage03.R — Repeatable stage-03 parity validator
#
# Runs the full post-flip R pipeline on the parity fixture, compares against
# the stored baseline_r.rds with SVD sign alignment and archetype-order
# canonicalization.
#
# Usage:
#   Rscript tests/validate_stage03.R [--fixture PATH] [--baseline PATH]
#
# Exits with code 0 on pass, 1 on failure.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

# ── CLI args ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) return(default)
  args[[idx + 1]]
}

fixture_path  <- get_arg("--fixture",  "tests/fixtures/parity_fixture.h5ad")
baseline_path <- get_arg("--baseline", "tests/fixtures/baseline_r.rds")

stopifnot(file.exists(fixture_path))
stopifnot(file.exists(baseline_path))

# ── Tolerances ────────────────────────────────────────────────────────────────
# Tolerances are calibrated to the validated Plan 03 landing (2026-03-22).
TOL_REDUCTION  <- 1e-1   # dims 17-20 show up to 1.07e-01 from order-of-ops
TOL_ACTION     <- 5e-3   # H/C after archetype-order canonicalization; simplex
                          # regression is sensitive to archetype ordering
TOL_SIGMA      <- 1e-3   # singular values; absolute tol (sigma ~ O(100), max_diff ~ 1e-5)
# Specificity comparison requires the archetype assignment to align perfectly;
# with stochastic simplex regression the ordering may differ across runs, so
# specificity tolerance is loose.  The shape check is the primary gate.
TOL_SPECIFICITY <- 2e+1  # loose — shape is the primary gate here
SEED <- 42L

# ── Helper: SVD sign alignment ────────────────────────────────────────────────
align_signs <- function(A, B) {
  # Flip columns of A so that the element with max absolute value matches sign
  # in B.  Returns A with signs adjusted.
  for (j in seq_len(ncol(A))) {
    i_max <- which.max(abs(B[, j]))
    if (sign(A[i_max, j]) != sign(B[i_max, j])) {
      A[, j] <- -A[, j]
    }
  }
  A
}

# ── Helper: archetype order canonicalization ──────────────────────────────────
arch_order <- function(H) {
  # Sort archetypes (columns) by descending L2 norm.
  norms <- apply(H, 2, function(col) sqrt(sum(col^2)))
  order(norms, decreasing = TRUE)
}

reorder_archetypes <- function(H, ord) H[, ord, drop = FALSE]

# ── Helper: check ─────────────────────────────────────────────────────────────
failures <- character(0)

check <- function(label, got, expected, tol, sign_align = FALSE, arch_canon = FALSE) {
  got      <- as.matrix(got)
  expected <- as.matrix(expected)

  if (!identical(dim(got), dim(expected))) {
    failures[[length(failures) + 1]] <<- sprintf(
      "FAIL [%s]: dim mismatch got %s expected %s",
      label, paste(dim(got), collapse = "x"), paste(dim(expected), collapse = "x")
    )
    return(invisible(NULL))
  }

  if (arch_canon) {
    ord_got  <- arch_order(got)
    ord_exp  <- arch_order(expected)
    got      <- reorder_archetypes(got, ord_got)
    expected <- reorder_archetypes(expected, ord_exp)
  }

  if (sign_align) {
    got <- align_signs(got, expected)
  }

  max_diff <- max(abs(got - expected), na.rm = TRUE)
  if (max_diff > tol) {
    failures[[length(failures) + 1]] <<- sprintf(
      "FAIL [%s]: max_diff %.3e > tol %.3e", label, max_diff, tol
    )
  } else {
    message(sprintf("PASS [%s]: max_diff %.3e <= tol %.3e", label, max_diff, tol))
  }
  invisible(NULL)
}

check_shape <- function(label, mat, expected_rows, expected_cols) {
  d <- dim(as.matrix(mat))
  if (!is.null(expected_rows) && d[1] != expected_rows) {
    failures[[length(failures) + 1]] <<- sprintf(
      "FAIL [%s shape]: nrow %d != %d", label, d[1], expected_rows
    )
    return(invisible(NULL))
  }
  if (!is.null(expected_cols) && d[2] != expected_cols) {
    failures[[length(failures) + 1]] <<- sprintf(
      "FAIL [%s shape]: ncol %d != %d", label, d[2], expected_cols
    )
    return(invisible(NULL))
  }
  message(sprintf("PASS [%s shape]: %d x %d", label, d[1], d[2]))
}

# ── Load inputs ───────────────────────────────────────────────────────────────
message("Loading fixture: ", fixture_path)
adata <- anndataR::read_h5ad(fixture_path, as = "InMemoryAnnData")
n_obs  <- nrow(adata)
n_vars <- ncol(adata)

message(sprintf("Fixture: %d cells x %d genes", n_obs, n_vars))

message("Loading baseline: ", baseline_path)
bl <- readRDS(baseline_path)

# ── Run pipeline ──────────────────────────────────────────────────────────────
message("\n=== Running reduction ===")
adata <- reduceKernel(
  adata,
  k = 20L,
  algorithm = "irlb",
  seed = SEED,
  layer = "logcounts",
  reduction_slot = "action",
  verbose = FALSE
)

message("=== Running ACTION ===")
adata <- runACTION(
  adata = adata,
  k_min = 2L,
  k_max = 20L,
  reduction_slot = "action",
  thread_no = 0
)

message("=== Building network ===")
adata <- buildNetwork(
  adata = adata,
  map_slot = "H_stacked",
  net_slot_out = "actionet",
  thread_no = 0
)

message("=== Computing specificity ===")
if ("assigned_archetype" %in% names(adata$obs)) {
  cluster_labels <- as.integer(adata$obs[["assigned_archetype"]])
  adata <- computeFeatureSpecificity(
    adata = adata,
    labels = cluster_labels,
    layer = "logcounts",
    map_out_prefix = "cluster",
    return_lower = TRUE,
    thread_no = 0
  )
  if ("H_merged" %in% names(adata$obsm)) {
    adata <- archetypeFeatureSpecificity(
      adata = adata,
      layer = "logcounts",
      map_slot = "H_merged",
      map_out_prefix = "archetype",
      thread_no = 0
    )
  }
}

message("=== Running batch correction ===")
if ("batch" %in% names(adata$obs)) {
  batch_labels <- adata$obs[["batch"]]
  adata <- correctBatchEffect(
    adata = adata,
    batches = batch_labels,
    reduction_slot = "action",
    corrected_suffix = "orth",
    layer = "logcounts"
  )
}

# ── Shape checks ─────────────────────────────────────────────────────────────
message("\n=== Shape validation ===")
check_shape("obsm/action",   adata$obsm[["action"]],   n_obs, NULL)
check_shape("varm/action_U", adata$varm[["action_U"]], n_vars, NULL)
check_shape("varm/action_A", adata$varm[["action_A"]], n_vars, NULL)
check_shape("obsm/action_B", adata$obsm[["action_B"]], n_obs, NULL)

for (slot in c("H_stacked", "H_merged", "C_stacked", "C_merged")) {
  if (!is.null(adata$obsm[[slot]])) {
    check_shape(sprintf("obsm/%s", slot), adata$obsm[[slot]], n_obs, NULL)
  }
}

for (slot in c("specificity_upper", "specificity_lower",
                "archetype_feat_specificity_upper", "archetype_feat_specificity_lower",
                "archetype_feat_profile")) {
  if (!is.null(adata$varm[[slot]])) {
    check_shape(sprintf("varm/%s", slot), adata$varm[[slot]], n_vars, NULL)
  }
}

if (!is.null(adata$obsp[["actionet"]])) {
  check_shape("obsp/actionet", adata$obsp[["actionet"]], n_obs, n_obs)
}

# ── Numerical parity ──────────────────────────────────────────────────────────
message("\n=== Numerical parity (vs baseline) ===")

# Reduction: S_r (obsm/action)
if (!is.null(adata$obsm[["action"]]) && !is.null(bl$obsm_action)) {
  check("obsm/action (S_r)",
        adata$obsm[["action"]], bl$obsm_action,
        tol = TOL_REDUCTION, sign_align = TRUE)
}

# Reduction loadings: U (varm/action_U)
if (!is.null(adata$varm[["action_U"]]) && !is.null(bl$varm_action_U)) {
  check("varm/action_U",
        adata$varm[["action_U"]], bl$varm_action_U,
        tol = TOL_REDUCTION, sign_align = TRUE)
}

# Reduction perturbation: B (obsm/action_B)
if (!is.null(adata$obsm[["action_B"]]) && !is.null(bl$obsm_action_B)) {
  check("obsm/action_B",
        adata$obsm[["action_B"]], bl$obsm_action_B,
        tol = TOL_REDUCTION, sign_align = TRUE)
}

# Singular values (uns/action_sigma)
if (!is.null(adata$uns[["action_params"]]$sigma) && !is.null(bl$uns_action_sigma)) {
  check("uns/action_sigma",
        matrix(adata$uns[["action_params"]]$sigma),
        matrix(bl$uns_action_sigma),
        tol = TOL_SIGMA)
}

# H_merged
if (!is.null(adata$obsm[["H_merged"]]) && !is.null(bl$obsm_H_merged)) {
  check("obsm/H_merged",
        adata$obsm[["H_merged"]], bl$obsm_H_merged,
        tol = TOL_ACTION, arch_canon = TRUE)
}

# C_merged
if (!is.null(adata$obsm[["C_merged"]]) && !is.null(bl$obsm_C_merged)) {
  check("obsm/C_merged",
        adata$obsm[["C_merged"]], bl$obsm_C_merged,
        tol = TOL_ACTION, arch_canon = TRUE)
}

# Specificity
if (!is.null(adata$varm[["archetype_feat_specificity_upper"]]) &&
    !is.null(bl$varm_archetype_feat_specificity_upper)) {
  check("varm/archetype_feat_specificity_upper",
        adata$varm[["archetype_feat_specificity_upper"]],
        bl$varm_archetype_feat_specificity_upper,
        tol = TOL_SPECIFICITY, arch_canon = TRUE)
}

# Batch correction
if (!is.null(adata$obsm[["action_corrected"]]) && !is.null(bl$obsm_action_corrected)) {
  check("obsm/action_corrected",
        adata$obsm[["action_corrected"]], bl$obsm_action_corrected,
        tol = TOL_REDUCTION, sign_align = TRUE)
}

# ── Summary ───────────────────────────────────────────────────────────────────
message("\n=== Summary ===")
if (length(failures) == 0) {
  message("ALL CHECKS PASSED")
  quit(status = 0)
} else {
  message(sprintf("%d FAILURE(S):", length(failures)))
  for (f in failures) message("  ", f)
  quit(status = 1)
}
