## Primary
* Fix `networkDiffusion` to keep dimnames.
* Fix broken R functions
* Edit the package 'DESCRIPTION'.
* Edit the exports in 'NAMESPACE', and add necessary imports.
* Consolidate decompACTION and runACTION
* Add rerun archetype merging
* Separation ACTION and ACTIONet
* Creature "run pipeline" function
* Fix plotting
 * interactive 3D coords
 * Generalize "correctBatchEffect" for matrix input after fixing back-end
 * Incorporate bact ortho into reduceKernel (requires generalized reduceKernel)
* Replace `plot.marker.enrichment.by.archetype.heatmap` and `plot.marker.enrichment.by.cluster.heatmap`
* Change coordinate ionitialization to "archetype_footprint" in `runACTIONet()`
* In `plot.ACTIONet()`, plotting order should have `NA` points on bottom.
* Fix bug in `plot.ACTIONet()` and `.layout_plot_labels()` when `label_attr` contains empty string (`""`)

## Secondary


## Done
* New marker detection
* WOrking and improved ledien via igraph
* New qc functions
* Removed `C_fastSpMatViewSum` due to inaccuracy when multithreading.
