# Decisions (ADR-lite)

This document records **deliberate architectural and operational decisions** for the ACTIONet ecosystem. These decisions are considered settled unless explicitly revised.

---

## Software architecture

### Multi-repo structure

**Decision:** Maintain separate repositories for:

- `libactionet` (C++ core)
- `actionet-r` (R front-end)
- `actionet-python` (Python front-end)
- `ACTIONetExperiment` (optional R compatibility container for migration)

**Rationale:**

- Clear separation of concerns
- Independent packaging and release cycles (C++ / CRAN-style / PyPI-style)
- Avoids monorepo friction while preserving coordination via shared specs

### Container contract

**Decision:**

- AnnData is the canonical data container contract across front-ends
- `actionet-r` uses `anndataR` as its primary container implementation
- `ACTIONetExperiment` remains supported only through optional compatibility converters and deprecated wrappers

**Rationale:**

- Reduces cross-language drift
- Simplifies interoperability between R and Python front-ends
- Removes duplicated container logic from `actionet-r`

---

## Language bindings

### C++ core + wrappers

**Decision:**

- C++ core built with **CMake**
- R bindings via **Rcpp**
- Python bindings via **pybind11**

**Rationale:**

- Mature, stable tooling
- Explicit control over ABI and performance
- Good compatibility with HPC environments

---

## Front-end prioritization

### Python vs R

**Decision:**

- Python front-end is the **performance-first and pipeline-critical interface**
- R front-end remains supported and serves as a **reference for semantics and outputs**

**Rationale:**

- R performance and ecosystem limitations at scale
- Python integration with pipeline and HPC workflows
- Preserve backward compatibility for existing R users

---

<!-- ## Reproducibility and stability

### Output contracts
**Decision:**
- Output formats, directory structures, and file naming are treated as **contracts**
- Changes require explicit documentation and migration plans

### Versioning
**Decision:**
- Critical dependencies (especially `actionet-python`) must be version-pinned and logged in pipeline runs

--- -->

## Change management

### Backward compatibility

**Decision:**

- Avoid breaking changes are allowed if justifed.
- Such changes must substantially improve:
  - Performance
  - Resource usage
  - User ease-of-use
  - Reproducibility

### Agent behavior

**Decision:**

- LLM/coding agents should not re-litigate decisions recorded in this document
- Deviations require explicit human approval

---

## UMAP / uwot disconnected-vertex repair (R surface)

**Decision:**

- `layoutNetwork()` exposes a `repair_disconnected = TRUE` argument, forwarded
  through `C_layoutNetwork` to `UwotArgs::repair_disconnected` in libactionet.
- The default is `TRUE`, matching the canonical umap-learn behavior. Existing
  R callers get the safeguard automatically without changing their code.
- The C++-side implementation is owned by libactionet and is documented in
  [`libactionet/context/DECISIONS.md`](../src/libactionet/context/DECISIONS.md)
  under "UMAP / uwot disconnected-vertex repair".

**Rationale:**

- The post-pruned graph in `optimize_layout_uwot` can leave fully isolated
  vertices that receive zero updates from the batch optimizer and remain
  frozen at their seed coordinate. With `runACTIONet`-style seeds (axes
  spanned by archetype footprints), such vertices land exactly on a
  coordinate axis, producing the cross-shaped artifact whose severity
  grows with N.
- The repair runs `connected_components_undirected` on the post-pruned graph
  and applies per-component recentering + offset (`N(0, 10)`) plus a small
  per-coordinate jitter (`N(0, 1e-4)`), seeded by `uwot_args.seed` so runs
  remain reproducible.
- Exposing the knob in R lets callers opt out for diagnostic or strict
  reproducibility scenarios while preserving the safeguard as the default.

**Reference:**

- libactionet ADR: "UMAP / uwot disconnected-vertex repair" (commit `0d17046`).
- Mirrored in `actionet-python` via `actionet.layout_network(repair_disconnected=...)`.

---
