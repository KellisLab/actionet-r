# Agent Playbook — actionet-r (R front-end via Rcpp)

## Purpose of this repository

This repository provides the **R user-facing interface** to the Actionet C++ core via **Rcpp**. It is the original front-end implementation and serves as a **reference for expected semantics and outputs** when maintaining parity across language bindings.

Upstream dependency:

- `libactionet` C++ core library
- `ACTIONetExperiment` R

Related front-end:

- `actionet-python` (pybind11)

---

## Repository layout (high level)

- `R/` — R functions and user-facing API
- `src/` — compiled code, Rcpp bridge, and `libactionet` submodule
- `man/` — Rd documentation
- `NAMESPACE`, `DESCRIPTION` — package metadata
- `inst/` — installed resources
- `data/` — bundled example or reference data
- `configure` — build configuration logic

---

## Role in the ecosystem

- Acts as the **canonical reference** for many behaviors and outputs
- Provides backward compatibility for existing R-based users
- Informs Python parity decisions unless an explicit divergence is documented

---

## What success looks like

- Stable, well-documented R API
- Clear mapping between R functions and C++ core behavior
- Documentation that accurately reflects current semantics
- Where applicable, alignment with Python front-end behavior

---

## Hard guardrails (must follow)

- Do not assume containerized or cloud-native execution environments.
- Maintain Rcpp binding hygiene (memory ownership, error handling).

---

## How to work safely in this repo

### Modifying existing functionality

1. Identify whether the change originates in the C++ core or the R layer.
2. Confirm whether Python front-end behavior should change as well.
3. Update R documentation (`man/`, roxygen comments) if semantics change.
4. Add or update tests if present (or add minimal regression checks).

### Adding new functionality

1. Define the R API (function signature, documentation) first.
2. Ensure naming and parameter conventions are consistent with existing API.
3. Implement Rcpp bindings with explicit ownership semantics.
4. Decide whether the feature should also be exposed in Python.

---

## Parity with Python front-end

- R is often the **reference implementation**.
- Python should mirror R outputs where feasible.
- If R diverges intentionally:
  - document the rationale

Recommended parity dimensions:

- function names
- parameter defaults
- return object structure
- column names and ordering
- determinism guarantees

---

## Build and environment assumptions

- Build system: CMake + Rcpp
- Typical environment: system R or conda-based R

Avoid:

- relying on system-wide nonstandard libraries
- embedding hard-coded paths or environment assumptions

---

## Common pitfalls

- Accidental semantic drift relative to Python
- Rcpp memory lifetime bugs
- Undocumented behavior changes

---

## When blocked, ask for

- Target `libactionet` version or commit
- Expected behavior or output examples
- Whether Python parity is required
- Build environment details (R version, compiler)
