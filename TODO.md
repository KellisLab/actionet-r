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
  * `assess_enrichment()`: `R/annotation.R:292` (`annotateArchetypes`), `R/annotation.R:794` (`annotate.profile.using.markers`); `R/enrichment.R:23` (`assess.TF.activities.from.scores`), `R/enrichment.R:83` (`assess.geneset.enrichment.from.scores`), `R/enrichment.R:139` (`assess.peakset.enrichment.from.archetypes`).
  * `XICOR()`: `R/annotation.R:307` and `R/annotation.R:325` (`annotateArchetypes` labels/scores branches).
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
