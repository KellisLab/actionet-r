---
name: Deprecated Unused Code Audit
overview: Read-only audit of `R/` in actionet-r identifying deprecated, dead, orphaned, and legacy code. Findings are grouped by severity so a human can decide what to delete, fix, or defer.
todos:
  - id: landmines
    content: "Fix remaining runtime landmines: bareword assess_enrichment in enrichment.R (3 sites) and annotate.profile.using.markers; roll_var in get.pseudobulk.SE; swapped-arg bug in plot.individual.gene; graveyard ref in annotate.archetypes.using.markers"
    status: pending
  - id: delete-dead-files
    content: Delete dead files R/alignment.R and R/projection.R (moved to R/_old_R/)
    status: completed
  - id: delete-dead-helpers
    content: "Delete verified-unused non-exported helpers: projectArchs (annotation.R), plotFeatureExpression / select.top.k.genes / gate.archetypes (plots.R), networkAutocorrelation + transitive chain (network_tools.R, enrichment.R), warnifnot / verify_aces / export_minimal_sce (utils_public.R), .create_design_formula (utils_stats.R)"
    status: pending
  - id: deprecate-legacy-annotation
    content: "Deprecate the legacy annotate.*.using.* family (5 exports) in favor of annotateArchetypes / annotateClusters / annotateCells; remove after a release. Note annotate.profile.using.markers has no direct replacement."
    status: pending
  - id: prune-old-R
    content: Delete R/_old_R/ contents after the two remaining live refs (annotate.archetypes.using.markers → assess.geneset.enrichment.from.archetypes; projectArchs → map.cell.scores.from.archetype.enrichment) are resolved
    status: pending
  - id: prune-rcpp-bindings
    content: Keep remaining 13 currently-unused C_* bindings in R/RcppExports.R (down from 15; annotation rewrite consumed C_assess_enrichment and C_XICOR)
    status: completed
  - id: shim-migration
    content: Route direct SummarizedExperiment access in pseudobulk_DGE.R and plots.R (plot.individual.gene) through the AnnData adapter shim
    status: pending
  - id: arg-routing
    content: "Add .resolve_container_arg to remaining offenders: plot.individual.gene and the pseudobulk_DGE.R family (get.pseudobulk.SE, run.ensemble.pseudobulk.DESeq/Limma, variance.adjusted.limma). Legacy annotate.*.using.* family also lacks routing but is being deprecated (see deprecate-legacy-annotation)."
    status: pending
  - id: consolidate-duplicates
    content: Pick canonical form for .validate_attr/.validate_vector_attr, .get_feature_vec/.get_features; migrate remaining callers
    status: pending
isProject: false
---


# Deprecated / unused code audit — actionet-r

Read-only audit of everything under [`R/`](R/) (top-level; `R/_old_R/` treated as a graveyard). Every entry below was verified by definition + call-site grep. No changes proposed — this is a triage document.

Cross-references: [context/DECISIONS.md](context/DECISIONS.md) (AnnData is canonical; ACE is compatibility-only), [plans/anndata_support_audit.md](plans/anndata_support_audit.md) (prior audit that flagged a subset of these issues), [TODO.md](TODO.md).

---

## 1. Landmines — will crash at runtime today

Exported functions that call symbols not defined anywhere in live `R/`.

- [`R/annotation.R`](R/annotation.R) line 262 — exported `annotate.archetypes.using.markers` calls `assess.geneset.enrichment.from.archetypes`, which lives only in [`R/_old_R/enrichment_ext.R`](R/_old_R/enrichment_ext.R) line 14. The new `annotateArchetypes` supersedes this function; see §5.
- [`R/enrichment.R`](R/enrichment.R) line 23, 83, 139 — exported `assess.TF.activities.from.archetypes` / `.from.scores`, `assess.geneset.enrichment.from.scores`, `assess.peakset.enrichment.from.archetypes` all call bareword `assess_enrichment`; only `C_assess_enrichment` exists at [`R/RcppExports.R`](R/RcppExports.R) line 318.
- [`R/annotation.R`](R/annotation.R) line 928 — exported `annotate.profile.using.markers` calls bareword `assess_enrichment` (same missing symbol).
- [`R/pseudobulk_DGE.R`](R/pseudobulk_DGE.R) line 128 — inside exported `get.pseudobulk.SE` (`.make_ensemble_assays` branch), `roll_var` is undefined in `R/*.R`.
- [`R/plots.R`](R/plots.R) line 672 — exported `plot.individual.gene` calls `.preprocess_annotation_labels(ace, labels)` but the helper signature at [`R/utils_internal_misc.R`](R/utils_internal_misc.R) line 38 is `function(labels, ace = NULL)` — args swapped.

Non-exported callers that would crash if invoked (lower urgency but still broken):
- [`R/annotation.R`](R/annotation.R) line 691 — unexported `projectArchs` (still present at line 690) calls `map.cell.scores.from.archetype.enrichment` (only in [`R/_old_R/annotation_ext.R`](R/_old_R/annotation_ext.R) line 18).

Fixed since previous audit (2026-07-21):
- ~~`annotateArchetypes` calls bareword `assess_enrichment` / `XICOR`~~ — the rewrite now uses `C_assess_enrichment` (via `.annotate_from_markers`) and `C_XICOR` directly.
- ~~`R/alignment.R` / `R/projection.R` broken refs~~ — files moved to `R/_old_R/` (see §2).

---

## 2. Dead files (nothing exported, nothing called externally) — DONE

Both moved to `R/_old_R/` on 2026-07-21:

- ~~`R/alignment.R`~~ — now [`R/_old_R/alignment.R`](R/_old_R/alignment.R).
- ~~`R/projection.R`~~ — now [`R/_old_R/projection.R`](R/_old_R/projection.R).

Follow-up: these are still referenced from §1 (landmines in [`R/annotation.R`](R/annotation.R)) but only via unrelated symbols. Full deletion of `R/_old_R/` is tracked in §9.

---

## 3. Dead non-exported helpers (verified: only referenced at definition site)

- [`R/annotation.R`](R/annotation.R) — `projectArchs` (line 690). Still unexported, still calls the broken graveyard ref `map.cell.scores.from.archetype.enrichment` (see §1/§9).
- [`R/plots.R`](R/plots.R) — `plotFeatureExpression` (line 743), `select.top.k.genes` (line 856), `gate.archetypes` (line 1072).
- [`R/network_tools.R`](R/network_tools.R) — `networkAutocorrelation` (line 279). Transitively kills `assess.categorical.autocorrelation`, `compute.phi`, `assess.continuous.autocorrelation`, `compute.Geary.C` in [`R/enrichment.R`](R/enrichment.R) lines 315-430.
- [`R/utils_public.R`](R/utils_public.R) — `warnifnot` (line 98), `verify_aces` (line 106), `export_minimal_sce` (line 229).
- [`R/utils_stats.R`](R/utils_stats.R) — `.create_design_formula` (line 1).

Resolved since previous audit (2026-07-21):
- ~~`annotateArchetypes`~~, ~~`annotateClusters`~~ — both now exported ([`NAMESPACE`](NAMESPACE) lines 25, 27) and functional. `annotateArchetypes` at [`R/annotation.R`](R/annotation.R) line 356; `annotateClusters` at line 549. Both route through `.resolve_container_arg` + `.validate_ace(..., as_ace = TRUE, allow_se_like = TRUE)`, use `.annotate_from_markers` at [`R/utils_annotation.R`](R/utils_annotation.R) line 87 (which calls `C_assess_enrichment`), and use `C_XICOR` directly in the labels/scores branches. This also incidentally consumes two of the previously-unused `C_*` bindings from §6 (`C_XICOR`, `C_assess_enrichment`).

---

## 4. Explicit `.Deprecated()` shims (alive by design)

Kept for backwards compat; candidates for eventual removal.

- [`R/plots.R`](R/plots.R) line 3 — `plot.ACTIONetExperiment` → advises `toAnnData()` + `plot()`. Still `S3method(plot,ACTIONetExperiment)` in [`NAMESPACE`](NAMESPACE) line 11.
- [`R/filter_ace.R`](R/filter_ace.R) lines 181, 187 — `filter.ace`, `filter.ace.by.attr` → forward to `filterActionet` / `filterActionetByAttr`. S3 methods still in NAMESPACE.
- [`R/normalization.R`](R/normalization.R) lines 21, 67 — `assay_out` arg deprecated in favor of `layer_out` in `normalize.ace` / `normalize.multiBatchNorm`.
- [`R/utils_anndata_adapter.R`](R/utils_anndata_adapter.R) lines 698, 715, 719 — arg-alias deprecations (`assay_name`→`layer`, `obj`/`ace`→`adata`).

Comment-marked "to be deprecated" (no `.Deprecated()` call yet):
- [`R/utils_validation.R`](R/utils_validation.R) line 207 — `.validate_attr` marked "To be deprecated in favor of `.validate_vector_attr()`".
- [`R/utils_internal_misc.R`](R/utils_internal_misc.R) line 23 — `.get_feature_vec` marked "Deprecated". Still called from `annotation.R:259,286,355`, `imputation.R:139,168`, `plots.R:673`.

---

## 5. Duplicate / superseded function pairs

- `.validate_attr` ([`R/utils_validation.R`](R/utils_validation.R) line 208) vs `.validate_vector_attr` (same file) — the latter is used from [`R/pseudobulk_DGE.R`](R/pseudobulk_DGE.R).
- `.get_feature_vec` ([`R/utils_internal_misc.R`](R/utils_internal_misc.R) line 24) vs `.get_features` (same file line 1) — deprecated form is still the more-called one.
- `filter.ace` / `filter.ace.by.attr` vs `filterActionet` / `filterActionetByAttr` — see §4.
- `plot.ACTIONetExperiment` vs `plot.AbstractAnnData` / `plot.InMemoryAnnData` ([`R/plots.R`](R/plots.R) lines 35, 40) — see §4.
- **New (2026-07-21):** modern `annotateArchetypes` / `annotateClusters` supersede the old snake-case S3-style family (`annotate.archetypes.using.markers`, `annotate.archetypes.using.labels`, `annotate.clusters.using.markers`, `annotate.clusters.using.labels`, `annotate.profile.using.markers`). Both are still exported ([`NAMESPACE`](NAMESPACE) lines 20-24 vs 25/27). The old family also carries the §1 landmines and the §8 arg-routing gaps. Recommended path: mark the old family with `.Deprecated()` pointing at the new camelCase functions, then remove after a release. `annotate.profile.using.markers` and `annotate.clusters.using.markers` are the tricky ones — `annotateClusters` covers cluster-mode + markers via `specificity_key`, but the matrix-only `annotate.profile.using.markers` signature has no direct replacement.
- Related: `projectArchs` (unexported, line 690) is superseded by nothing exported and is currently broken (§1). Candidate for deletion rather than migration.

---

## 6. Unused Rcpp bindings surfaced to R (kept — decision recorded)

**Decision (2026-07-21): keep all remaining.** Retain as R-side entry points; targets for wiring up the §1 landmines.

Since the annotation rewrite (2026-07-21), `C_assess_enrichment` and `C_XICOR` are now called from [`R/utils_annotation.R`](R/utils_annotation.R) line 126 and [`R/annotation.R`](R/annotation.R) lines 445, 469, 641, 666 — so they are no longer "unused".

Thirteen `C_*` symbols in [`R/RcppExports.R`](R/RcppExports.R) still have no `R/*.R` caller:

- `C_runAA` (line 20), `C_runSimplexRegression` (89), `C_runSPA` (101), `C_orthogonalizeBasal` (167), `C_orthogonalizeBasal_full` (171), `C_autocorrelation_Moran_parametric` (283), `C_autocorrelation_Moran` (287), `C_autocorrelation_Geary` (291), `C_normalizeMatrixSparse` (379), `C_normalizeMatrixDense` (383), `C_MWM_hungarian` (408), `C_MWM_rank1` (412), `C_xicor` (416).

Note: `C_xicor` (lowercase, line 416) is distinct from `C_XICOR` (uppercase, line 420); the annotation rewrite uses the uppercase one.

---

## 7. Legacy container access smells

Live code that reaches directly into `SummarizedExperiment` on a user-supplied `ace`, bypassing the adapter shim at [`R/utils_anndata_adapter.R`](R/utils_anndata_adapter.R). Per [context/DECISIONS.md](context/DECISIONS.md), ACE is compatibility-only.

- [`R/pseudobulk_DGE.R`](R/pseudobulk_DGE.R) line 33 — `SummarizedExperiment::assays(ace)[[assay_name]]` in exported `get.pseudobulk.SE`.
- [`R/pseudobulk_DGE.R`](R/pseudobulk_DGE.R) line 90 — `SummarizedExperiment::rowData(ace)` in `get.pseudobulk.SE`.
- [`R/plots.R`](R/plots.R) line 705 — `SummarizedExperiment::assays(ace)[[assay_name]]` in exported `plot.individual.gene`; also breaks on AnnData.
- Additional occurrences at lines 188, 193-194, 230, 268, 275-276, 315, 345, 352, 357 of [`R/pseudobulk_DGE.R`](R/pseudobulk_DGE.R) operate on locally-constructed SEs (legitimate but flagged for review).

Unqualified `rowData()` / `assays()` calls (fragile — depend on `SummarizedExperiment` being attached at call time):
- [`R/batch_correct.R`](R/batch_correct.R) line 56 — `rowData(mnn.out)`.
- [`R/pseudobulk_DGE.R`](R/pseudobulk_DGE.R) line 204 — `rowData(dds)`.
- [`R/projection.R`](R/projection.R) line 83 — `assays(bulk)` (file is dead anyway).

---

## 8. Broken deprecated-arg routing

Exported functions that take `ace` positionally but do NOT call `.resolve_container_arg`, so the deprecated-alias contract silently doesn't apply.

- [`R/annotation.R`](R/annotation.R) — `annotate.archetypes.using.labels` (line 177), `annotate.archetypes.using.markers` (line 255), `annotate.clusters.using.labels` (line 730), `annotate.clusters.using.markers` (line 798), `annotate.profile.using.markers` (line 925). All belong to the superseded snake-case family in §5.
- [`R/plots.R`](R/plots.R) line 666 — `plot.individual.gene`.
- [`R/pseudobulk_DGE.R`](R/pseudobulk_DGE.R) — `get.pseudobulk.SE`, `run.ensemble.pseudobulk.DESeq`, `run.ensemble.pseudobulk.Limma`, `variance.adjusted.limma`.

Resolved since previous audit (2026-07-21):
- ~~`annotateCells`~~ — now routes through `.resolve_container_arg(adata = adata, ace = ace)` at [`R/annotation.R`](R/annotation.R) line 43 and `.validate_ace(..., as_ace = TRUE, allow_se_like = TRUE)` at line 44. This also closes the "AnnData-only, rejects ACE" bug flagged in [plans/anndata_support_audit.md](plans/anndata_support_audit.md) Category B.
- ~~`annotateArchetypes`~~, ~~`annotateClusters`~~ — routed via `.resolve_container_arg` at lines 368, 561 (new implementations).

---

## 9. Live references into the graveyard

`R/_old_R/` should have zero live callers. Two remain, both in [`R/annotation.R`](R/annotation.R) and both in code superseded by the new `annotateArchetypes` (§5):

- Line 262 → `assess.geneset.enrichment.from.archetypes` in [`R/_old_R/enrichment_ext.R`](R/_old_R/enrichment_ext.R). Called from the legacy `annotate.archetypes.using.markers`.
- Line 691 → `map.cell.scores.from.archetype.enrichment` in [`R/_old_R/annotation_ext.R`](R/_old_R/annotation_ext.R). Called from the unused/broken `projectArchs`.

Newly-graveyarded (2026-07-21): [`R/_old_R/alignment.R`](R/_old_R/alignment.R), [`R/_old_R/projection.R`](R/_old_R/projection.R) — no live callers.

Everything else in `R/_old_R/` (`old_main.R`, `SCINET_interface.R`, `clusters_ext.R`, remainder of `annotation_ext.R`/`enrichment_ext.R`) has no live callers and is safe to delete once the two references above are resolved (e.g. by deleting `annotate.archetypes.using.markers` and `projectArchs`).

---

## 10. Possibly-dead shim helpers (uncertain)

Not conclusive because they may be intended API surface, but currently unused:

- [`R/utils_anndata_adapter.R`](R/utils_anndata_adapter.R) — `.tscalet` (line 38, only used from dead [`R/alignment.R`](R/alignment.R)), `.fast_col_means` (line 30, only used by dead `gate.archetypes`), `.set_feature_data` (line 313), `.set_obs_data` (line 338), `.set_row_data_df` (line 300), `.set_col_data_df` (line 325), `.actionet_nrow` (line 164), `.actionet_ncol` (line 171).
- `plot.ACTIONet` at [`R/plots.R`](R/plots.R) line 75 — an internal branch references `plotFeatureExpression`; the branch is unreachable given current dispatch and `plotFeatureExpression` is dead.

---

## Suggested triage buckets (for a follow-up cleanup PR)

1. **Delete-safe now** — §3 (dead helpers, incl. `projectArchs`). Full `R/_old_R/` deletion once the two §9 refs land.
2. **Fix, don't delete** — §1 (remaining landmines: wire `assess_enrichment` sites to `C_assess_enrichment`; add `roll_var`; fix `plot.individual.gene` arg order), §7 (route through adapter), §8 (add `.resolve_container_arg` to the non-legacy offenders).
3. **Deprecate then remove** — §5 legacy `annotate.*.using.*` family (5 exports) in favor of the new `annotateArchetypes` / `annotateClusters` / `annotateCells`. Simultaneously resolves the annotation-side entries in §1, §8, and §9.
4. **Defer** — §4 (deprecation shims kept intentionally), §10 (uncertain shim helpers).
5. **Policy call** — remaining §5 duplicate pairs (`.validate_attr` vs `.validate_vector_attr`, `.get_feature_vec` vs `.get_features`).
6. **Settled** — §2 (dead alignment/projection files moved to `R/_old_R/`); §6 (unused `C_*` bindings retained; count now 13, down from 15); `annotateArchetypes` + `annotateClusters` rewritten and exported (2026-07-21).

Full audit transcript with grep evidence: [Audit R deprecated/unused code](c2cac021-0839-48ce-a6ca-4f2856398da7).
