#!/usr/bin/env Rscript
#
# Baseline capture script for Plan 00 — Parity Baseline Infrastructure.
#
# Runs the full canonical ACTIONet pipeline on the shared parity fixture and
# saves all parity-critical slots to:
#   - tests/fixtures/baseline_r.h5ad  (full AnnData via anndataR)
#   - tests/fixtures/baseline_r.rds   (named list of extracted matrices)
#
# Usage (from the actionet-r repo root):
#   Rscript tests/generate_baseline.R
#
# Requirements: actionet (installed or devtools::load_all()), anndataR
#

suppressPackageStartupMessages({
  library(Matrix)
  library(anndataR)
})

# Always load actionet from the source to get the current dev-backed version
# (the installed package may be outdated).
if (file.exists(file.path(getwd(), "DESCRIPTION"))) {
  message("Loading actionet via pkgload::load_all() from current directory ...")
  suppressPackageStartupMessages(
    {
      options(pkgbuild.override_build_tools = TRUE)
      pkgload::load_all(".", quiet = TRUE, reset = TRUE, compile = FALSE, recompile = FALSE)
    }
  )
} else if (requireNamespace("actionet", quietly = TRUE)) {
  message("Loading installed actionet package ...")
  suppressPackageStartupMessages(library(actionet))
} else {
  stop(
    "Package 'actionet' is not installed and this script is not being run ",
    "from the actionet-r repo root. Please install actionet or cd to the repo root."
  )
}

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
args_full <- commandArgs(trailingOnly = FALSE)
file_arg  <- args_full[grepl("^--file=", args_full)]
if (length(file_arg) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  # Fallback when sourced interactively or path detection fails
  SCRIPT_DIR <- file.path(getwd(), "tests")
}

FIXTURE_PATH <- file.path(SCRIPT_DIR, "fixtures", "parity_fixture.h5ad")
OUT_H5AD     <- file.path(SCRIPT_DIR, "fixtures", "baseline_r.h5ad")
OUT_RDS      <- file.path(SCRIPT_DIR, "fixtures", "baseline_r.rds")
SEED         <- 42L

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.check_fixture <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Fixture not found: ", path, "\n",
      "Run libactionet/test/fixtures/generate_fixture.py first."
    )
  }
}

.extract_arrays <- function(adata) {
  arrays <- list()

  .get <- function(collection, key, label) {
    val <- tryCatch(collection[[key]], error = function(e) NULL)
    if (!is.null(val)) {
      if (is(val, "sparseMatrix")) {
        arrays[[label]] <<- as.matrix(val)
      } else {
        arrays[[label]] <<- as.matrix(val)
      }
      message(sprintf("  [OK]      %-40s  %s", label, paste(dim(arrays[[label]]), collapse = "x")))
    } else {
      message(sprintf("  [MISSING] %s", label))
    }
  }

  # Reduction
  .get(colMaps(adata),  "action",   "obsm_action")
  .get(colMaps(adata),  "action_B", "obsm_action_B")
  .get(rowMaps(adata),  "action_U", "varm_action_U")
  .get(rowMaps(adata),  "action_A", "varm_action_A")
  uns <- tryCatch(adata$uns, error = function(e) list())
  if ("action_params" %in% names(uns) && "sigma" %in% names(uns[["action_params"]])) {
    arrays[["uns_action_sigma"]] <- as.numeric(uns[["action_params"]][["sigma"]])
    message(sprintf("  [OK]      %-40s  length=%d", "uns_action_sigma", length(arrays[["uns_action_sigma"]])))
  } else {
    message("  [MISSING] uns_action_params/sigma")
  }

  # ACTION
  .get(colMaps(adata), "H_stacked", "obsm_H_stacked")
  .get(colMaps(adata), "H_merged",  "obsm_H_merged")
  .get(colMaps(adata), "C_stacked", "obsm_C_stacked")
  .get(colMaps(adata), "C_merged",  "obsm_C_merged")
  obs_data <- tryCatch(adata$obs, error = function(e) NULL)
  if (!is.null(obs_data) && "assigned_archetype" %in% colnames(obs_data)) {
    arrays[["obs_assigned_archetype"]] <- as.integer(obs_data[["assigned_archetype"]])
    message(sprintf("  [OK]      %-40s  length=%d", "obs_assigned_archetype", length(arrays[["obs_assigned_archetype"]])))
  } else {
    message("  [MISSING] obs/assigned_archetype")
  }

  # Network
  nets <- tryCatch(colNets(adata), error = function(e) list())
  if ("actionet" %in% names(nets)) {
    G <- nets[["actionet"]]
    arrays[["obsp_actionet"]] <- as.matrix(G)
    message(sprintf("  [OK]      %-40s  %s", "obsp_actionet", paste(dim(arrays[["obsp_actionet"]]), collapse = "x")))
  } else {
    message("  [MISSING] obsp/actionet")
  }

  # Specificity (cluster) — R stores only upper by default; request lower too
  .get(rowMaps(adata), "cluster_upper",  "varm_specificity_upper")
  .get(rowMaps(adata), "cluster_lower",  "varm_specificity_lower")

  # Specificity (archetype)
  .get(rowMaps(adata), "archetype_feat_profile",             "varm_archetype_feat_profile")
  .get(rowMaps(adata), "archetype_feat_specificity_upper",   "varm_archetype_feat_specificity_upper")
  .get(rowMaps(adata), "archetype_feat_specificity_lower",   "varm_archetype_feat_specificity_lower")

  # Batch correction — slot names: action_{corrected_suffix}, action_U_{corrected_suffix}, action_A_{corrected_suffix}
  .get(colMaps(adata), "action_orth",   "obsm_action_corrected")
  .get(rowMaps(adata), "action_U_orth", "varm_action_corrected_U")
  .get(rowMaps(adata), "action_A_orth", "varm_action_corrected_A")

  # Layout
  .get(colMaps(adata), "actionet_2d", "obsm_actionet_2d")

  arrays
}

# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

run_pipeline <- function(adata) {
  message("Step 1: reduceKernel ...")
  adata <- reduceKernel(
    adata,
    k            = 20L,
    algorithm    = "irlb",
    seed         = SEED,
    layer        = "logcounts",
    reduction_slot = "action",
    verbose      = FALSE
  )

  message("Step 2: runACTION ...")
  adata <- runACTION(
    adata,
    k_min          = 2L,
    k_max          = 20L,
    reduction_slot = "action"
  )

  message("Step 3: buildNetwork ...")
  adata <- buildNetwork(
    adata,
    map_slot     = "H_stacked",
    net_slot_out = "actionet"
  )

  # Cluster labels from archetype assignment
  obs_data <- adata$obs
  labels   <- as.integer(obs_data[["assigned_archetype"]])

  message("Step 4: computeFeatureSpecificity ...")
  adata <- computeFeatureSpecificity(
    adata,
    labels         = labels,
    layer          = "logcounts",
    map_out_prefix = "cluster",
    return_lower   = TRUE
  )

  message("Step 5: archetypeFeatureSpecificity ...")
  adata <- archetypeFeatureSpecificity(
    adata,
    layer          = "logcounts",
    map_slot       = "H_merged",
    map_out_prefix = "archetype"
  )

  message("Step 6: correctBatchEffect ...")
  obs_data <- adata$obs
  batches  <- obs_data[["batch"]]
  adata <- correctBatchEffect(
    adata,
    batches           = batches,
    reduction_slot    = "action",
    corrected_suffix  = "orth",
    layer             = "logcounts"
  )

  message("Step 7: layoutNetwork ...")
  adata <- layoutNetwork(
    adata,
    net_slot     = "actionet",
    seed         = SEED,
    verbose      = FALSE,
    map_slot_out = "actionet_2d"
  )

  adata
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

.check_fixture(FIXTURE_PATH)

message("Loading fixture: ", FIXTURE_PATH)
adata <- anndataR::read_h5ad(FIXTURE_PATH, as = "InMemoryAnnData")
message(sprintf("  Shape: %d cells x %d genes", adata$n_obs(), adata$n_vars()))

message("\nRunning pipeline (seed=", SEED, ")...")
adata <- run_pipeline(adata)

message("\nSaving full AnnData -> ", OUT_H5AD)
if (file.exists(OUT_H5AD)) file.remove(OUT_H5AD)
adata$write_h5ad(OUT_H5AD)

message("Extracting parity arrays -> ", OUT_RDS)
arrays <- .extract_arrays(adata)
saveRDS(arrays, OUT_RDS)

message("\nSaved ", length(arrays), " arrays:")
for (nm in sort(names(arrays))) {
  v <- arrays[[nm]]
  message(sprintf("  %-50s  dim=%s", nm, paste(dim(v), collapse = "x")))
}

message("\nR baseline complete.")
