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
* Finish `annotateClusters`
  * Compare incomplete XICOR version to actionet-python
* Add error for `annotateClusters` for when feat_spec is not in object.
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
