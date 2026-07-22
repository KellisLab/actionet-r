## Primary
* Fix `networkDiffusion` to keep dimnames.
* Deprecate and remove old R functions
* Edit the package 'DESCRIPTION'.
* Edit the exports in 'NAMESPACE', and add necessary imports.
* Finish rerun archetype merging
* Separation ACTION and ACTIONet
* Fix plotting
 * interactive 3D coords
 * Generalize "correctBatchEffect" for matrix input after fixing back-end
 * Incorporate bact ortho into reduceKernel (requires generalized reduceKernel)
* Replace `plot.marker.enrichment.by.archetype.heatmap` and `plot.marker.enrichment.by.cluster.heatmap`
* Change coordinate ionitialization to "archetype_footprint" in `runACTIONet()`
* In `plot.ACTIONet()`, plotting order should have `NA` points on bottom.
* Fix bug in `plot.ACTIONet()` and `.layout_plot_labels()` when `label_attr` contains empty string (`""`)
* Missing R wrappers `assess_enrichment()` and `XICOR()` (only `C_assess_enrichment` and `C_XICOR` Rcpp entry points exist). Broken call sites still using the unprefixed names:
  * `assess_enrichment()`: `R/annotation.R:794` (`annotate.profile.using.markers`); `R/enrichment.R:23` (`assess.TF.activities.from.scores`), `R/enrichment.R:83` (`assess.geneset.enrichment.from.scores`), `R/enrichment.R:139` (`assess.peakset.enrichment.from.archetypes`).
  * `R/enrichment.R:86` and `R/enrichment.R:142` also read `$thresholds` from the enrichment result; the C++ binding now returns `$peak_rank_idx` (positional 0-based rank, not a score threshold). Rename the key at the call site when adding the R wrapper.
  * Fix path: either add thin R wrappers (`assess_enrichment <- function(...) C_assess_enrichment(...)` and same for `XICOR`) or route each caller to the `C_*` name directly.
* Force all zero rows removed if `min_cells_per_feat` > 0 in `filter.ace`
* Fix impute with single gene

## Secondary
* Deprecate `compute_specificity_parallel`
* Change ggtheme to bw
* Normalize only subset of features

## Done
* New marker detection
* Working and improved ledien via igraph
* New qc functions
* Removed `C_fastSpMatViewSum` due to inaccuracy when multithreading.
* Add arbitrary pseudo count
* Finish `annotateClusters` (full parity with Python `annotate_clusters`; AnnData-first signature; fixed specificity slot lookup; added de-novo `computeFeatureSpecificity` path; added error when pre-computed specificity is missing; switched to `C_assess_enrichment`/`C_XICOR` and `.encode_markers`).
* Modernize `annotateArchetypes` (2026-07-21). AnnData-first signature (`adata`, `ace=` deprecated, `archetype_slot`, `specificity_key`, `features_use`, `layer`, `thread_no`, `archetype_specificity_slot=` deprecated with back-compat suffix stripping). Marker branch reads `{specificity_key}_upper` / `_lower` and forms `pmax(upper - lower, 0)` (previously only upper); falls back to `archetypeFeatureSpecificity(return_raw = TRUE)` when `specificity_key = NULL`. Labels/scores branches preserve the continuous-H XICOR path but now use `C_XICOR` and `.encode_markers`. Lowercase return keys (`labels`, `confidence`, `enrichment`) match `annotateClusters`. Roxygen `@export`; man page and `NAMESPACE` regenerated.
* Add Python parity `annotate_archetypes()` in `actionet-python/src/actionet/annotation/annotation.py` (exported at both `actionet.annotation` and `actionet` namespaces). Mirrors R semantics: markers via `_core.assess_enrichment` on the `pmax(upper - lower, 0)` specificity delta, labels/scores via `_core.xicor_matrix` on the continuous `adata.obsm[archetype_slot]`.
* Fix long-standing correctness bugs in the `libactionet` C++ core (2026-07-21; identical patch shipped in both actionet-r and actionet-python submodules):
  * `xicor` Z-score `ind` — switched to 1-based (`regspace(1, n)`) matching `XICOR::xicor` asymptotic. Prior 0-based version returned a systematically wrong Z on every input.
  * `rank_vec(method=1)` tie handling — now returns 1-based max-tie rank matching `R::rank(., ties.method="max")`.
  * `xicor` seed — now controls random tie-breaking on X via joint permutation + `stable_sort_index`; documented in the header.
  * `XICOR` matrix path — removed the `swap(X, Y)` + transpose optimization, which was silently wrong for the asymmetric xi statistic. Added per-column rank precomputation.
  * `assess_enrichment` — no longer mutates its `associations` argument; signature is now `const arma::sp_mat&`.
  * `assess_enrichment` output — the second returned matrix is now `peak_rank_idx` (0-based rank position, not a score threshold); the R and Python bindings were updated to expose the new name. Legacy `$thresholds` in `R/enrichment.R` call sites will need to be renamed when those wrappers are fixed (see wrapper-gap bullet above).
  * `assess_enrichment` — hoisted inner-loop scratch matrices, dropping 4× per-iteration allocations.
  * Golden regression tests added under `src/libactionet/tests/test_xicor_enrichment.cpp` (opt-in via `-DLIBACTIONET_BUILD_TESTS=ON`).
