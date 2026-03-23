#!/usr/bin/env Rscript
# validate_stage06.R — Stage 06: Unified Specificity validator
#
# Verifies:
#  1. computeFeatureSpecificity is non-mutating (S unchanged after call)
#  2. Specificity output shapes are (n_vars, n_clusters)
#  3. annotateCells runs without dimension errors with cells x genes S
#  4. annotateCells output shapes are correct
#
# Usage:
#   cd actionet-r
#   Rscript tests/validate_stage06.R

suppressPackageStartupMessages({
  options(pkgbuild.override_build_tools = TRUE)
  pkgload::load_all(".", quiet = TRUE, compile = FALSE, recompile = FALSE)
})

PASS <- 0L
FAIL <- 0L

check <- function(name, cond, detail = NULL) {
  if (isTRUE(cond)) {
    cat(sprintf("  PASS  %s\n", name))
    PASS <<- PASS + 1L
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
# Load fixture
# ---------------------------------------------------------------------------
section("Loading fixture")

fixture_path <- file.path("tests", "fixtures", "parity_fixture.h5ad")
if (!file.exists(fixture_path)) {
  stop(sprintf("Fixture not found: %s", fixture_path))
}

adata <- anndataR::read_h5ad(fixture_path, as = "InMemoryAnnData")
n_obs  <- nrow(adata)   # cells
n_vars <- ncol(adata)   # genes
cat(sprintf("  Loaded %d obs x %d vars\n", n_obs, n_vars))

gene_names <- adata$var_names

# Cluster labels for specificity
set.seed(42)
n_clusters <- 5L
labels_vec <- sample(paste0("C", seq_len(n_clusters)), n_obs, replace = TRUE)

# ---------------------------------------------------------------------------
# Section 1: Non-mutating check
# ---------------------------------------------------------------------------
section("1. Non-mutating specificity check")

# Snapshot expression matrix before call
X_layer <- adata$layers[["logcounts"]]
if (is.null(X_layer)) X_layer <- adata$X
X_before_copy <- X_layer + 0  # force deep copy

out_raw <- computeFeatureSpecificity(
  adata,
  labels   = labels_vec,
  layer    = "logcounts",
  return_raw = TRUE,
  thread_no  = 1L
)

X_after <- if (!is.null(adata$layers[["logcounts"]])) adata$layers[["logcounts"]] else adata$X

if (inherits(X_before_copy, "sparseMatrix")) {
  check("S unchanged after specificity (data)",
        isTRUE(all.equal(X_before_copy@x, X_after@x)))
  check("S unchanged after specificity (dims)",
        all(dim(X_before_copy) == dim(X_after)))
} else {
  check("S unchanged after specificity",
        isTRUE(all.equal(as.matrix(X_before_copy), as.matrix(X_after))))
}

# ---------------------------------------------------------------------------
# Section 2: Specificity output shapes
# ---------------------------------------------------------------------------
section("2. Specificity output shapes (n_vars x n_clusters)")

check("upper_significance is a matrix",
      is.matrix(out_raw$upper_significance))
check("upper_significance: nrow == n_vars",
      nrow(out_raw$upper_significance) == n_vars,
      sprintf("got %d, expected %d", nrow(out_raw$upper_significance), n_vars))
check("upper_significance: ncol == n_clusters",
      ncol(out_raw$upper_significance) == n_clusters,
      sprintf("got %d, expected %d", ncol(out_raw$upper_significance), n_clusters))
check("lower_significance shape matches upper",
      all(dim(out_raw$lower_significance) == dim(out_raw$upper_significance)))
check("average_profile shape: (n_vars, n_clusters)",
      !is.null(out_raw$average_profile) &&
      nrow(out_raw$average_profile) == n_vars &&
      ncol(out_raw$average_profile) == n_clusters,
      sprintf("got %s", if (is.null(out_raw$average_profile)) "NULL"
              else paste(dim(out_raw$average_profile), collapse = " x ")))
check("upper_significance all non-negative",
      all(out_raw$upper_significance >= 0))
check("lower_significance all non-negative",
      all(out_raw$lower_significance >= 0))

# ---------------------------------------------------------------------------
# Section 2b: Sparse vs Dense in-memory parity
# ---------------------------------------------------------------------------
section("2b. Sparse vs Dense in-memory parity")

X_layer_dense <- as.matrix(if (!is.null(adata$layers[["logcounts"]])) adata$layers[["logcounts"]] else adata$X)

adata_dense <- adata$clone()
adata_dense$X <- X_layer_dense

out_dense <- computeFeatureSpecificity(
  adata_dense,
  labels   = labels_vec,
  layer    = NULL,
  return_raw = TRUE,
  thread_no  = 1L
)

check("Sparse == Dense: upper shape matches",
      all(dim(out_raw$upper_significance) == dim(out_dense$upper_significance)))
check("Sparse == Dense: upper values match",
      isTRUE(all.equal(out_raw$upper_significance, out_dense$upper_significance,
                        tolerance = 1e-10)),
      sprintf("max_diff=%.2e",
              max(abs(out_raw$upper_significance - out_dense$upper_significance))))
check("Sparse == Dense: lower values match",
      isTRUE(all.equal(out_raw$lower_significance, out_dense$lower_significance,
                        tolerance = 1e-10)),
      sprintf("max_diff=%.2e",
              max(abs(out_raw$lower_significance - out_dense$lower_significance))))
check("Sparse == Dense: average_profile values match",
      isTRUE(all.equal(out_raw$average_profile, out_dense$average_profile,
                        tolerance = 1e-10)),
      sprintf("max_diff=%.2e",
              max(abs(out_raw$average_profile - out_dense$average_profile))))

# ---------------------------------------------------------------------------
# Section 3: annotateCells with cells x genes S
# ---------------------------------------------------------------------------
section("3. annotateCells accepts cells x genes S")

# Build network if not present
if (!"actionet" %in% names(adata$obsp)) {
  set.seed(0)
  nz <- 500L
  rows <- sample.int(n_obs, nz, replace = TRUE)
  cols <- sample.int(n_obs, nz, replace = TRUE)
  G <- Matrix::sparseMatrix(i = rows, j = cols, x = rep(1.0, nz),
                            dims = c(n_obs, n_obs))
  G <- G + Matrix::t(G)
  G@x <- rep(1.0, length(G@x))
  adata$obsp[["actionet"]] <- G
}

n_marker <- min(5L, n_vars)
markers <- list(
  type_A = gene_names[seq_len(n_marker)],
  type_B = gene_names[seq(n_marker + 1L, 2L * n_marker)]
)

tryCatch({
  res <- annotateCells(
    ace          = adata,
    markers      = markers,
    method       = "vision",
    assay_name   = "logcounts",
    thread_no    = 1L
  )
  check("annotateCells returns labels",
        !is.null(res$labels))
  check("annotateCells labels length == n_obs",
        length(res$labels) == n_obs,
        sprintf("got %d", length(res$labels)))
  check("annotateCells enrichment rows == n_obs",
        nrow(res$enrichment) == n_obs,
        sprintf("got %d", nrow(res$enrichment)))
  check("annotateCells enrichment cols == n_celltypes",
        ncol(res$enrichment) == 2L,
        sprintf("got %d", ncol(res$enrichment)))
}, error = function(e) {
  check("annotateCells runs without error", FALSE, conditionMessage(e))
})

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat(sprintf("\n%s\n  RESULTS: %d passed, %d failed\n%s\n",
            strrep("=", 60), PASS, FAIL, strrep("=", 60)))

quit(status = if (FAIL == 0L) 0L else 1L)
