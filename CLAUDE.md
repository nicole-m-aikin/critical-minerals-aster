# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project summary

Spectral alteration mapping of the McDermitt Caldera (NV/OR) using ASTER thermal infrared (TIR, B10–B14) band ratio analysis. Produces strong anomaly zone polygons and validates them against USGS MRDS mineral deposit data. Multi-site capable via per-site YAML configuration.

## Figure / plot quality check

After saving any figure, map, or plot (`.png`, `.jpg`, `.pdf`), **always open and visually inspect the output file** with the Read tool (which renders images inline) to verify:
- The hillshade / basemap covers the full data bbox (no white border or size mismatch)
- All intended layers are visible (zones, deposits, annotations, scale bar, north arrow)
- No obvious rendering artifacts (clipped labels, legend cut off, blank panels)

Reading the file path is not enough — you must **look at the rendered image** and describe what you see before marking the task complete.

## Environment setup

```bash
conda env create -f environment.yml
conda activate aster-minerals
pip install -e .
```

The editable install (`pip install -e .`) is required for notebooks and tests to import `critical_minerals_aster` from `src/`.

## Commands

**Run pipeline (single site):**
```bash
python -m critical_minerals_aster run --site mcdermitt
python -m critical_minerals_aster run --site mcdermitt --download   # fetch from EarthData first
python -m critical_minerals_aster run --site mcdermitt --skip-figures
```

**Batch and synthesis:**
```bash
python -m critical_minerals_aster run-batch --all-sites
python -m critical_minerals_aster synthesize
python -m critical_minerals_aster gee-parity          # local vs GEE TIR parity (8-site panel)
python -m critical_minerals_aster gee-parity --gee --project YOUR_GCP_PROJECT
```

**Tests:**
```bash
pytest tests/
pytest tests/test_classification.py   # single test file
pytest tests/test_classification.py::test_band_ratio_divide_by_zero  # single test
```

## Architecture

### Data flow

```
sites/{id}.yaml
    → SiteConfig (config.py)
    → SitePaths  (paths.py)       # resolves data/aster/, vectors/, figures/, results/
    → pipeline.run_site()         # orchestrates all steps
        ├── download_and_mosaic_aster()        → AST_05 emissivity QI/CI/MI mosaic + cls rasters
        ├── spectral.alteration_indices()      → QI, CI, MI (silica, carbonate, mafic slots)
        ├── classification.classify_percentiles() → 3-class maps per ratio
        ├── classification.combined_score()    → additive score 0–6
        ├── classification.vectorize_strong_zones() → GeoDataFrame (GeoJSON)
        ├── metrics.compute_site_summary()    → MRDS spatial join + CSV
        └── pipeline.write_provenance()       → results/{id}_provenance.json
```

### Key modules (`src/critical_minerals_aster/`)

| Module | Purpose |
|---|---|
| `config.py` | `SiteConfig` dataclass + YAML loader; `ClassificationParams`, `StructureLayer` |
| `paths.py` | `SitePaths` — all file/dir paths, layout-aware (`flat` vs `nested`) |
| `spectral.py` | ASTER TIR I/O, granule selection scoring, `alteration_ratios()` |
| `classification.py` | Percentile classification, `combined_score()`, polygon vectorization |
| `metrics.py` | MRDS spatial join, per-site summary CSV |
| `mrds.py` | MRDS CSV → GeoDataFrame with CRS reprojection |
| `structure.py` | Distance-to-fault and buffer annotation for deposits |
| `synthesis.py` | Aggregate `results/*_summary.csv` → national summary |
| `hyperspectral.py` | GEMx / AVIRIS CMR coverage lookup (no credentials; footprints only) |
| `gee_parity.py` | Local rasterio vs GEE TIR parity on an 8-site panel (recipe + probes; live `ee` optional) |
| `basemap.py` | Fetch Esri World Imagery via `contextily`, reproject to hillshade grid, cache to `data/basemap/{site_id}/rgb.tif` |
| `pipeline.py` | `run_site()` / `run_batch()` orchestration + figure generation; `FIG03_VERSION` controls cache invalidation |
| `__main__.py` | `argparse` CLI entry point |

### Site configuration (`sites/`)

Each site has a YAML with: `id`, `name`, `bbox_wgs84`, `granule_id` (or `null` for auto-select), `layout` (`flat` for current McDermitt, `nested` for multi-site), `classification` params, `temporal` range, and optional `structure_layers`.

`sites/index.yaml` lists site ids for batch runs.

**Layout difference:** `flat` writes to `data/aster/` and `data/vectors/`; `nested` writes to `data/sites/{id}/aster/` and `data/sites/{id}/vectors/`. All path logic lives in `SitePaths`.

### Spectral indices (Ninomiya 2003 / Ninomiya et al. 2005)

Computed on **AST_05 v004 daytime surface emissivity** (`spectral.alteration_indices`), stored in the legacy slot names:

- **Quartz index QI** (`silica` slot): B11² / (B10·B12)
- **Carbonate index CI** (`carbonate` slot): B13 / B14
- **Mafic index MI** (`mafic` slot): B12 / B13

Classification uses per-scene percentile thresholds (default 70th/90th). Combined score ≥ 3 defines "strong anomaly" zones. Thresholds are scene-relative; cross-site comparison of raw scores is not meaningful.

Legacy (`spectral.alteration_ratios`, still used by `gee_parity` and the bbox script): radiance ratios "silica" B13/B14 (actually CI), "carbonate" B13/B12 and "mafic" B12/B13 (reciprocals of each other). Radiance ratios track surface temperature, so zones moved between outcrop and valley fill with acquisition time — see `docs/data_verification.md` §2.4 and §6.

### Input product and mosaicking

`download_and_mosaic_aster(..., product="AST_05")` (default; `"AST_L1T"` = legacy radiance path):
- Selection (`select_mosaic_granules`): clear-sky (≤30% scene cloud) **day** scenes, May–Sep first, never mixing day and night, greedy bbox cover. A pinned `granule_id` selects the AST_05 granule of the same acquisition.
- Per-pixel QA (`ast05_qa_good`): QA_DataPlane quality = good, cloud = clear, cloud adjacency ≤ 1; all bands in [0.5, 1.0].
- Indices per granule, priority-stacked onto a 90 m bbox grid: the reference granule wins; later granules only fill gaps after mean/std matching on the overlap (no blending).
- AG100 v003 water and vegetation (NDVI > 0.3) masks removed before classification (`ag100.py`).
- Writes `{id}_mosaic_ratio_*` (QI/CI/MI), `_cls_*`, `_cls_valid`, `_TIR_B10` (B10 emissivity) and `{id}_mosaic_granules.json` (selection, QA rule, mask fractions).

Night AST_05 (`aster_night/`) is a replication check only (`scripts/day_night_replication.py`), never part of the score.

### Primary null

`scripts/site_specific_null_significance.py`: **bedrock-vs-cover stratified** null (SGMC "Unconsolidated" = cover; `lithology.cover_strata`). Each deposit's null probability is the zone coverage of its own stratum; exact Poisson-binomial p. The geometric null is kept as `*_geometric` columns.

### Outputs

- `data/vectors/strong_anomaly_zones.geojson` — vectorized strong anomaly polygons
- `figures/0*.png` — band ratio maps, classification maps, deposit overlay
- `results/{site_id}_summary.csv` — MRDS hit rates by commodity, Earth MRI category, and mineral system
- `results/{site_id}_provenance.json` — granule id, git commit, package versions

### Summary CSV row types

Each `results/*_summary.csv` contains four `row_type` values:

| `row_type` | Key column | Description |
|---|---|---|
| `site` | — | Site-level totals |
| `commodity` | `commodity_group` | Simplified commodity buckets (Gold/Silver, Uranium, etc.) |
| `earth_mri` | `earth_mri_category` | 8 Earth MRI priority groups (Energy, REE, Battery Metals…) |
| `mineral_system` | `mineral_system` | One of 24 OFR 2020-1042 mineral systems (Placer, Porphyry Cu-Mo-Au…) |

### Earth MRI / mineral-system classification (`mrds.py`)

Two independent classifiers operate on each MRDS row:

**`reclassify_mrds_earth_mri(df)`** — adds `earth_mri_category` from `commod1`.
Priority: Energy > REE > Battery Metals > PGM > Base Metals > Specialty/High-Tech > Gold/Silver > Industrial > Non-Critical.
Reflects the 2022 Final List of 50 Critical Minerals (includes scandium, cesium, rubidium, arsenic vs. the 2018 list).

**`reclassify_mrds_mineral_system(df)`** — adds `mineral_system` from a three-field cascade:
  1. `model` field (USGS deposit-model code string, e.g. `"53: Porphyry Cu"`) — most specific
  2. `dep_type` field (e.g. `"Placer"`, `"Stratabound"`)
  3. `commod1` field — broadest fallback

Coverage: ~84% of critical-mineral MRDS deposits are assigned a system; bulk non-critical (sand/gravel/stone) remain `"Unknown System"`.

**`is_critical_mineral(earth_mri_category)`** — returns `True` for all categories except `"Non-Critical"`.

### EarthData / ASTER data

Raw ASTER TIR rasters live under `data/aster/` / `data/sites/{id}/aster/` (not committed). AST_05 files are `{granule_id}_SRE_TIR_B{10-14}.tif` (uint16 × 0.001) plus `_SRE_TIR_QA_DataPlane.tif`; legacy L1T radiance is `{granule_id}_TIR_B{10-14}.tif` (archived in `aster_l1t_day/`, older mixed day/night set in `aster_archive_pre_clearday/`). The ASTER TIR subsystem has been off since 2026-01; V004 archive products remain downloadable. The `--download` flag triggers `earthaccess.login(strategy="interactive")` on first run; subsequent runs use cached credentials. ASTER SWIR bands (B04–B09) are not available in LP DAAC v004. GEMx / AVIRIS coverage of the 109 sites is inventoried by `scripts/gemx_coverage_report.py` (CMR metadata, no login); 71 sites have GEMx-era SWIR flight lines. That is coverage only — no reflectance ingest or mineral ID yet.

**GEE parity** (`gee_parity.py`, `docs/gee_parity.md`) ports the TIR recipe for 8 representative sites to Earth Engine (`ASTER/AST_L1T_003`) and benchmarks it against the local rasterio path. Default run is offline discrepancy probes (clip order, histogram percentiles, radiance UCC, EPSG:4326 resample). Live `ee` is optional (`pip install -e ".[gee]"`). This is a comparison module, not a 109-site rewrite.
