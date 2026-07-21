# AnnData Support Audit — actionet-r

Date: 2026-07-21

Scope: Inventory of AnnData vs `ACTIONetExperiment` (ACE) support across
exported R functions in `actionet-r`, the state of concatenation in
`anndataR`, and applicability of the `actionet-python/io` submodule as a
template for extending AnnData functionality in R.

Related context:

- [context/PROJECT_CONTEXT.md](../context/PROJECT_CONTEXT.md)
- [context/DECISIONS.md](../context/DECISIONS.md) — AnnData is the canonical container; ACE remains as an optional compatibility path
- [context/AGENT_PLAYBOOK.md](../context/AGENT_PLAYBOOK.md)

---

## Executive summary

- AnnData support is broad and near-complete. The adapter shim in
  [`R/utils_anndata_adapter.R`](../R/utils_anndata_adapter.R) plus
  `.validate_ace()` in [`R/utils_validation.R`](../R/utils_validation.R)
  provides the "accepts both containers" contract used by most exported
  functions.
- Three ACE-flavored gaps remain where an AnnData input would fail today:
  `get.pseudobulk.SE`, `plot.individual.gene`, and (independently broken)
  `annotate.archetypes.using.markers`.
- Two functions moved too far in the opposite direction and now reject ACE:
  `annotateCells` and `correctBatchEffect`. Trivial one-line fixes.
- `anndataR` 1.0.0 has **no** native concatenation support. Upstream tracks
  it in [scverse/anndataR#325](https://github.com/scverse/anndataR/issues/325),
  unmerged as of Nov 28, 2025.
- The `actionet-python/io` submodule is **not** a source for concatenation
  (it contains none). It *is* a strong template for the separate problem of
  backed HDF5 AnnData lifecycle (subset/checkpoint/matrix streaming), which
  `anndataR` also lacks and which the R adapter currently short-circuits by
  aborting on `HDF5AnnData` inputs.

---

## Q1. Does every function accepting ACE have full AnnData support?

**Almost. Three real gaps and one paper cut.**

### Canonical "accepts both" pattern

```711:727:R/utils_anndata_adapter.R
.resolve_container_arg <- function(adata = NULL, ace = NULL, obj = NULL, default_obj = NULL) {
  out <- default_obj

  if (!is.null(obj)) {
    .Deprecated(msg = "'obj' is deprecated; use 'adata' instead.")
    out <- obj
  }
  if (!is.null(ace)) {
    .Deprecated(msg = "'ace' is deprecated; use 'adata' instead.")
    out <- ace
  }
  if (!is.null(adata)) {
    out <- adata
  }

  out
}
```

```17:44:R/utils_validation.R
.validate_ace <- function(
    obj,
    as_ace = FALSE,
    allow_se_like = FALSE,
    ...
    ) {
  ...
  if (!.is_anndata(obj)) {
    if (!.is_se_like(obj) || (!allow_se_like && !as_ace)) {
      if (error_on_fail) {
        stop(...)
      }
      return(FALSE)
    }
    obj <- toAnnData(obj)   # ACE → AnnData conversion happens here
  }
```

Functions that call `.validate_ace(..., allow_se_like = TRUE, as_ace = TRUE, return_elem = TRUE)`
(or use the direct-dispatch getters `colMaps` / `rowMaps` / `colNets` / `rowNets` / `.get_obs_data` / `.get_feature_data`)
accept both containers.

### CATEGORY A — Accepts both (majority of the API)

Via `.validate_ace(..., as_ace = TRUE, allow_se_like = TRUE)` (ACE path coerces to AnnData through `toAnnData()`):

- `runACTIONet` — [R/main.R:35-193](../R/main.R)
- `runACTION` — [R/r_action.R:1-72](../R/r_action.R)
- `reduceKernel` — [R/r_reduction.R:18-107](../R/r_reduction.R)
- `smoothKernel` — [R/r_reduction.R:110-182](../R/r_reduction.R)
- `computeFeatureSpecificity` — [R/r_specificity.R:3-72](../R/r_specificity.R)
- `archetypeFeatureSpecificity` — [R/r_specificity.R:77-136](../R/r_specificity.R)
- `layoutNetwork` — [R/r_visualization.R:2-194](../R/r_visualization.R)
- `computeNodeColors` — [R/r_visualization.R:197-226](../R/r_visualization.R)
- `filterActionet` / `filterActionetByAttr` / `filter.ace` / `filter.ace.by.attr` — [R/filter_ace.R](../R/filter_ace.R)
- `correctBatchEffectFastMNN` — [R/batch_correct.R:4-64](../R/batch_correct.R)
- `normalize.ace` — [R/normalization.R:2-41](../R/normalization.R)
- `normalize.multiBatchNorm` — [R/normalization.R:45-94](../R/normalization.R)
- `buildNetwork` / `networkDiffusion` / `networkCentrality` / `propagateLabels` — [R/network_tools.R](../R/network_tools.R)
- `clusterNetwork` — [R/clustering.R:2-99](../R/clustering.R)
- `imputeFeatures` / `impute.genes.using.archetypes` / `impute.specific.genes.using.archetypes` — [R/imputation.R](../R/imputation.R)
- `assess.TF.activities.from.archetypes` — [R/enrichment.R:49-56](../R/enrichment.R)
- `assess.peakset.enrichment.from.archetypes` — [R/enrichment.R:105-146](../R/enrichment.R)

Via shim getters only (no `.validate_ace` call; ACE path is direct, no `toAnnData` conversion):

- `getFeatureAbundance` / `getMitoAbundance` — [R/feature_stats.R](../R/feature_stats.R)
- `findMarkers` — [R/marker_detection.R:2-73](../R/marker_detection.R)
- `plotFeatureDist` / `plotMitoDist` — [R/feature_stats_plot.R](../R/feature_stats_plot.R)
- `annotate.archetypes.using.labels` — [R/annotation.R:174-231](../R/annotation.R)
- `annotate.clusters.using.labels` — [R/annotation.R:454-499](../R/annotation.R)
- `annotate.clusters.using.markers` — [R/annotation.R:522-633](../R/annotation.R)
- `plot.ACTIONet` / `plot.ACTIONet.gradient` / `plot.ACTIONet.interactive` / `plot.AbstractAnnData` / `plot.InMemoryAnnData` / `plot.ACTIONetExperiment` / `plot.ACTIONet.feature.view` / `plot.ACTIONet.gene.view` / `plot.top.k.genes` / `plot.archetype.selected.genes` / `plot.ACTIONet.archetype.footprint` — [R/plots.R](../R/plots.R)

### CATEGORY B — AnnData-only (rejects ACE)

Both use `.validate_ace(..., allow_se_like = FALSE)` with no `as_ace = TRUE`, so ACE inputs error out. Given the migration goal in `context/DECISIONS.md` this is inconsistent with the rest of the API.

- **`annotateCells`** — [R/annotation.R:22-101](../R/annotation.R)
  ```37:44:R/annotation.R
    method <- match.arg(method)
    norm_method <- match.arg(norm_method)
    norm_method <- ifelse(norm_method == "pagerank_sym", 2, 0)

    .validate_ace(ace, allow_se_like = FALSE, return_elem = FALSE, error_on_fail = TRUE)
    X <- .encode_markers(ace, markers = markers, features_use = features_use, obj_name = "ace")
    S <- .validate_assay(ace, assay_name = assay_name, matrix_type = "sparse", force_type = TRUE, error_on_fail = TRUE, return_elem = TRUE)
    G <- .validate_net(ace, net_slot = net_slot, matrix_type = "sparse", force_type = TRUE, return_elem = TRUE)
  ```
- **`correctBatchEffect`** — [R/batch_correct.R:68-182](../R/batch_correct.R). Same pattern at line 113.

Fix in both cases: swap `allow_se_like = FALSE` for `allow_se_like = TRUE, as_ace = TRUE`.

### CATEGORY C — ACE-only (would fail on AnnData)

- **`get.pseudobulk.SE`** — [R/pseudobulk_DGE.R:2-99](../R/pseudobulk_DGE.R). Calls `SummarizedExperiment::assays(ace)[[assay_name]]`, `SummarizedExperiment::rowData(ace)`, and subsets `ace[, mask]` directly. Never runs `toAnnData()`.
  ```26:34:R/pseudobulk_DGE.R
      ace <- ace[, ace[[sample_attr]] %in% names(sample_counts)]
      group_vec <- .validate_vector_attr(ace, attr = sample_attr, return_type = "data")
      ...
    counts_mat <- SummarizedExperiment::assays(ace)[[assay_name]]
  ```
  Migration path: swap direct SE calls for `.get_layer_matrix()`, `.get_feature_data()`, and `.subset_actionet_container()` ([R/utils_anndata_adapter.R:674-689](../R/utils_anndata_adapter.R)).

- **`plot.individual.gene`** — [R/plots.R:666-723](../R/plots.R). Uses `SummarizedExperiment::assays(ace)[[assay_name]]` at line 705. Also has a likely argument-order bug at line 672: `.preprocess_annotation_labels(ace, labels)` but the helper signature is `(labels, ace)` ([R/utils_internal_misc.R:38](../R/utils_internal_misc.R)).

### Broken independent of container type

- **`annotate.archetypes.using.markers`** — [R/annotation.R:252-273](../R/annotation.R) calls `assess.geneset.enrichment.from.archetypes(...)` which is not defined in the current `R/` (only in `R/_old_R/`). Likely a stale reference.

### Paper cut

- **`clusterNetwork`** — [R/clustering.R:2-99](../R/clustering.R) never calls `.resolve_container_arg`, so it silently rejects the deprecated `adata = ` / `ace = ` aliases. Cosmetic, but breaks the uniformity of the API.

### CATEGORY D — Matrix-only exports (out of scope for this audit)

`mergeArchetypes`, `runSVD`, `normalize.matrix`, `aggregateMatrix`, `assess.TF.activities.from.scores`, `assess.geneset.enrichment.from.scores`, `assess.genesets`, `annotate.profile.using.markers`, `geneset.enrichment.gProfiler`, `assess.geneset.enrichment.gProfiler`, `select.top.k.features`, `plot.top.k.features`, `variance.adjusted.limma`, `run.ensemble.pseudobulk.DESeq`, `run.ensemble.pseudobulk.Limma`.

---

## Q2. Is ACE still fully supported as a backwards-compatible alternative?

**Yes, for practical purposes.** ACE remains a first-class input across the modern API:

- Most exported functions coerce ACE to AnnData at the boundary via `toAnnData()` and return AnnData. `toACTIONetExperiment()` (defined at [R/utils_anndata_adapter.R:643](../R/utils_anndata_adapter.R)) is available for the reverse trip.
- A subset (`.ace_or_assay`, `.ace_or_map`, `.ace_or_net` code paths) preserves the container type end-to-end: ACE in → ACE out.
- The exceptions are the two AnnData-only functions in Category B (trivially fixable) plus `get.pseudobulk.SE` in Category C (inherently ACE-flavored since it produces a `SummarizedExperiment`).

This matches the intent in [context/DECISIONS.md](../context/DECISIONS.md):
"ACE remains supported only through optional compatibility converters and deprecated wrappers."

---

## Q3. AnnData concatenation in R — can it be added?

### State of `anndataR`

- **Installed version:** 1.0.0 (Bioconductor 3.22, 2025-10-30).
- **No native concatenation.** Verified three ways:
  1. NAMESPACE exports only 7 functions: `AnnData`, `AnnDataView`, `as_AnnData`, `generate_dataset`, `get_generator_types`, `read_h5ad`, `write_h5ad`.
  2. `grep("concat|rbind|cbind|merge|append|combine|bind", ls(asNamespace("anndataR")), ignore.case=TRUE)` → `character(0)`.
  3. GitHub issue [scverse/anndataR#325](https://github.com/scverse/anndataR/issues/325) tracks concat; open since Aug 2025. Maintainer comment Nov 28, 2025: *"for the `anndata::concat` function I currently don't have a replacement yet, so this would need to be implemented before I'm fully able to claim anndata is superseded by anndataR."*

### What `anndataR` gives us to build on

- `[i, j]` subsetting → returns lazy `AnnDataView` (composable). See method signature:
  ```r
  `[.AbstractAnnData` <- function(x, i, j, drop = TRUE, ...) {
    if (inherits(x, "AnnDataView")) {
      return(x$subset(i, j))
    }
    AnnDataView$new(x, i, j)
  }
  ```
- Direct read/write on all slots (`X`, `layers`, `obs`, `var`, `obs_names`, `var_names`, `obsm`, `varm`, `obsp`, `varp`, `uns`) via active bindings.
- `write_h5ad` / `read_h5ad` for persistence.
- Backed HDF5 support via `HDF5AnnData` (modes `r`, `r+`, `a`, `w`, `w-`, `x`); Zarr on `main`.
- Conversions: `as_SingleCellExperiment`, `as_Seurat`, `as_InMemoryAnnData`, `as_ReticulateAnnData`, `as_HDF5AnnData`.

### Feasibility

- **In-memory concat:** straightforward. `~300-500 LOC`, no C++. Standard `Matrix::` + `dplyr::bind_rows` work: union/intersect obs+var indices, stack `.X` and `.layers`, block-diagonal for `.obsp`/`.varp`, `.uns` merge policy.
- **Backed concat:** harder. Requires streaming writes into a fresh HDF5 file, plus compression handling. Blocked by the same problem as `subset_backed_inplace` — see Q3-adjacent below.

### Recommended concat scope (parking-lot decision)

Match the reference (`SummarizedExperiment::rbind`/`cbind` and `anndata.concat`) on *logically-concatable-axes only*, plus reasonable extras:

- `concat_anndata(list, axis = "obs")` — cells across shared-feature datasets. Union (outer) or intersection (inner) of `var_names`. Stack `.X` and `.layers` row-wise. Row-bind `.obs` (fill missing columns with `NA`). Row-bind `.obsm` slots present in all inputs. Block-diagonal for `.obsp`. Drop `.varm`/`.varp` unless every input carries the same slot (feature-space; must agree). `.uns` merge policy: `"first"` / `"drop"` / `"unique"`.
- `axis = "var"` — mirror of the above.
- ACE parity check: verify `rbind`/`cbind` semantics from
  [../ACTIONetExperiment](../../ACTIONetExperiment/) match. The `SummarizedExperiment`
  `rbind`/`cbind` methods it inherits require identical rowData/colData for the
  perpendicular axis, which is close to option B.

---

## Q3-adjacent. `actionet-python/io` submodule as a template

**The Python `io` submodule contains no concatenation code.** `np.concatenate` appears exactly once in the entire `actionet-python` codebase, in `visualization/qc.py` for plotting.

It *does* solve four other problems that `anndataR` also lacks. The R adapter currently short-circuits by rejecting backed AnnData:

```9:16:R/utils_anndata_adapter.R
.abort_if_backed_anndata <- function(x, arg = "adata") {
  if (.is_backed_anndata(x)) {
    stop(sprintf(
      "'%s' is a backed AnnData object. Backed object support is not implemented in actionet-r yet; materialize it in memory first.",
      arg
    ))
  }
}
```

If backed support becomes a goal, the Python `io` submodule is a strong template.

### Port candidates (decreasing priority)

- **A. `subset_backed_inplace` + `materialize_backed`** — atomic on-disk row/col rewrite of an h5ad, only safe way to shrink a backed object without doubling file size. Highest value; `anndataR` has no equivalent.
- **B. `checkpoint_backed` + `_DirtyTracker` + optional HDF5 repack** — flush partial annotation updates without a full file rewrite. Solves the "`write_h5ad` doubles the file every save" problem.
- **C. `MatrixSource`** — unified chunked accessor over dense/sparse × in-memory/backed with row/col aggregations (`row_sums`, `col_sums`, `nnz_row_counts`, `nnz_col_counts`, `row_sum_of_squares`, `feature_subset`, `apply_rowwise`).
- **D. `LazyTransform` / `create_lazy_transform`** — deferred library-size + log1p on backed `.X`. Only useful once R wraps the same `libactionet` backed operator via Rcpp.

### Required supporting infrastructure

- `anndata_io.append_to_anndata` + private HDF5 writers (`_write_matrix`, `_write_dataframe_to_h5`, `_write_uns_value_to_file`, etc.). Would need reimplementation in R against `rhdf5` / `hdf5r`, respecting the AnnData `encoding-type` / `encoding-version` schema.
- `compression.CompressionPolicy` + storage-metadata helpers.
- `persist._refresh_backed_handle` / `_DirtyTracker` — small, self-contained.

### Do not port

- `operator.py` — pybind11-specific.
- `lazy_transform.py` — depends on `operator.py`.
- AnnData 0.12/0.13-specific workarounds — R will have different quirks.
- Python `StringDtype` coercions — no-op in R.

### Suggested port order (smallest, most-independent first)

1. `compression.py` (pure metadata, no deps).
2. `anndata_io.py` HDF5 writer + validation (depends on 1).
3. `persist.py` in-memory + eager write path + `_DirtyTracker` (depends on 2).
4. `checkpoint.py` (depends on 2, 3).
5. `subset.py` `subset_backed_inplace` + `materialize_backed` (depends on 1, 2, 3).
6. `matrix_source.py` `MatrixSource` — can happen in parallel with 2-5.
7. (`operator.py` + `lazy_transform.py`) — deferred until R wraps `libactionet` backed operator.

---

## Follow-up workstreams (not scheduled)

1. Close the three ACE→AnnData gaps: rewrite `get.pseudobulk.SE`, fix `plot.individual.gene`, restore or replace `assess.geneset.enrichment.from.archetypes`.
2. Make `annotateCells` and `correctBatchEffect` accept ACE like the rest of the API (one-line changes each).
3. Implement in-memory `concat_anndata()` — logically-concatable-axes-only, ACE-parity scope.
4. (Optional, larger) Port the `actionet-python/io` backed-AnnData lifecycle: `subset_backed_inplace`, `materialize_backed`, `checkpoint_backed`, `MatrixSource`.

Nothing is scheduled from this audit. Cross-reference this document when
picking up any of the above.

---

## Appendix: transcript

- Audit R functions AnnData vs ACE support: [3face953-89b1-4f56-b56a-b5ea5e8822d1](3face953-89b1-4f56-b56a-b5ea5e8822d1)
- Inspect anndataR concatenation support: [1aac4b64-0c69-45c0-b667-9639e481f38d](1aac4b64-0c69-45c0-b667-9639e481f38d)
- Inspect actionet-python io module: [5b0b4ada-6291-47da-8f2d-95cf019920f1](5b0b4ada-6291-47da-8f2d-95cf019920f1)
