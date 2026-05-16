# Implementation phases (piecewise rollout)

The original multi-site plan is split into **independent milestones**. Complete and validate each before starting the next.

## Phase A — Site registry and paths

- **Goal:** One YAML per study site; path helpers for `data/` and `figures/` without clobbering runs.
- **Deliverables:** [`sites/*.yaml`](../sites/), `SiteConfig` + `SitePaths` in code, [`.gitignore`](../.gitignore) updated for nested layouts.
- **Exit criteria:** Load McDermitt config; resolve correct directories for `layout: flat` (current repo) and `layout: nested` (future `data/sites/{id}/...`).

## Phase B — Shared Python package

- **Goal:** Remove duplicated `load_band`, `band_ratio`, `classify`, polygonization; CRS-safe MRDS prep.
- **Deliverables:** Editable install via [`pyproject.toml`](../pyproject.toml), package under [`src/critical_minerals_aster/`](../src/critical_minerals_aster/).
- **Exit criteria:** Notebooks 02–03 call the package; notebook 04 reprojects deposits with `to_crs(zones.crs)` instead of hardcoded UTM zone.

## Phase C — Granule selection policy

- **Goal:** Replace `results[7]` / `results[8]` with scoring (bbox coverage, band availability) + optional `granule_id` override in YAML.
- **Deliverables:** Module + notebook 01 refactor.
- **Exit criteria:** McDermitt still resolves to the same granule when override set; second site works with auto-pick.

## Phase D — CLI / batch runner

- **Goal:** `run --site mcdermitt` (and batch file) producing vectors, figures, and **tabular** summaries under `results/`.
- **Deliverables:** `python -m critical_minerals_aster ...` or `scripts/run_site.py`, provenance logging.
- **Exit criteria:** Fresh clone + env + one command reproduces McDermitt outputs (modulo EarthData fetch).

## Phase E — National synthesis

- **Goal:** Aggregate per-site CSV/Parquet; compare hit rates and commodity patterns across sites.
- **Deliverables:** Notebook or script under `notebooks/` or `scripts/`, README subsection on interpretation limits.
- **Exit criteria:** Two or more sites appear in a single summary table/figure.

## Phase F — Structural geology hooks (optional)

- **Goal:** Optional line layers; distance-to-fault and buffer flags for deposits and/or anomalies.
- **Deliverables:** Config keys for structure paths; `structure.py` with distance/buffer utilities; one pilot layer documented.

---

**Current status:** Phases **A** and **B** are implemented in-repo (site YAML, `src/critical_minerals_aster/`, notebooks 02–04 wired to the library where noted); **C–F** remain.
